/*  Copyright (C) 2024-2026 P. David Buchan (pdbuchan@gmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// matrix.c - Derive the RGB/XYZ conversion matrices.

// gcc -Wall src/matrix.c -lm -o matrix

// Usage: ./matrix
// Output: reports to stdout

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>

// Function prototypes
int inputtext (char *);
int rgb_primaries (double **, int *);
int illum_white (double *, int);
int gaussjordan (int, double **);
int parse_int_string (const char *, int *);
double *allocate_doublemem (size_t);
double **allocate_doublememp (size_t);
char *allocate_strmem (size_t);

// Set some symbolic constants.
#define MAX_STRINGLEN 256  // Maximum number of characters per line

typedef struct {
  const char *name;
  double xr, yr;
  double xg, yg;
  double xb, yb;
  const char *source_quality;
  const char *source_note;
  const char *native_white_name;
  double xw, yw;
  int native_white_known;
} RGBSpace;

typedef struct {
  const char *name;
  double x, y;
  const char *source_note;
} WhitePoint;

// Primary coordinates are kept at the precision at which the defining standard,
// originating profile, or historical source specifies them.  Do not add digits
// merely by substituting a more precise chromaticity for a named illuminant: for
// standardized RGB spaces, the reference-white coordinates are part of the RGB
// space definition.
//
// "standard" means a normative standard/registry was available; "originator/profile"
// means an originating author or widely distributed defining profile was available;
// "historical/secondary" means no normative definition was located and the values
// are retained from the Kang/Pascale source set used by the original program.

static const RGBSpace rgb_spaces[] = {
  {"sRGB / ITU-R BT.709", 0.6400,0.3300, 0.3000, 0.6000, 0.1500, 0.0600,
   "standard", "IEC 61966-2-1 / ITU-R BT.709", "D65 (space definition)", 0.3127, 0.3290, 1},
  {"Adobe RGB (1998)", 0.6400, 0.3300, 0.2100, 0.7100, 0.1500,0.0600,
   "standard", "Adobe RGB (1998) Color Image Encoding", "D65 (space definition)", 0.3127, 0.3290, 1},
  {"Apple RGB", 0.6250, 0.3400, 0.2800, 0.5950, 0.1550, 0.0700,
   "originator/profile", "legacy Apple/ColorSync profile definition; cross-checked against Adobe profile set", "D65 (legacy profile)", 0.3127, 0.3290, 1},
  {"Best RGB", 0.7347, 0.2653, 0.2150, 0.7750, 0.1300, 0.0350,
   "originator/profile", "Bruce Lindbloom working-space definition", "D50 (working-space definition)", 0.3457, 0.3585, 1},
  {"Beta RGB", 0.6888, 0.3112, 0.1986, 0.7551, 0.1265, 0.0352,
   "originator/profile", "Bruce Lindbloom working-space definition", "D50 (working-space definition)", 0.3457, 0.3585, 1},
  {"Bruce RGB", 0.6400, 0.3300, 0.2800, 0.6500, 0.1500, 0.0600,
   "originator/profile", "Bruce Lindbloom working-space definition", "D65 (working-space definition)", 0.3127, 0.3290, 1},
  {"CIE 1931 RGB (2 degree observer)",
   0.7346900, 0.2653100, 0.2736827, 0.7174214, 0.1665347, 0.0088840,
   "historical/derived", "700.0, 546.1 and 435.8 nm primaries evaluated from official CIE 1931 2 degree CMFs", "equal-energy E", 1.0 / 3.0, 1.0 / 3.0, 1},
  {"CIE 1964/RGB (legacy secondary table)", 0.7232, 0.2768, 0.1248, 0.8216, 0.1616, 0.0134,
   "historical/secondary", "retained from Kang; CIE standardized a 10 degree XYZ observer, not this RGB working space", "equal-energy E (historical convention)", 1.0 / 3.0, 1.0 / 3.0, 1},
  {"ColorMatch RGB", 0.6300, 0.3400, 0.2950, 0.6050, 0.1500, 0.0750,
   "originator/profile", "legacy ColorMatch profile definition", "D50 (profile definition)", 0.3457, 0.3585, 1},
  {"Don RGB 4", 0.6960, 0.3000, 0.2150, 0.7650, 0.1300, 0.0350,
   "originator/profile", "Don Hutcheson / HutchColor working-space definition", "D50 (working-space definition)", 0.3457, 0.3585, 1},
  {"EBU Tech. 3213-E", 0.6400, 0.3300, 0.2900, 0.6000, 0.1500, 0.0600,
   "standard", "EBU Tech. 3213-E / ITU-R BT.470 625-line systems", "D65", 0.3127, 0.3290, 1},
  {"eciRGB v2", 0.6700, 0.3300, 0.2100, 0.7100, 0.1400, 0.0800,
   "standard", "ECI / ISO 22028-2 registry definition", "D50 (space definition)", 0.3457, 0.3585, 1},
  {"Ekta Space PS5", 0.6950, 0.3050, 0.2600, 0.7000, 0.1100, 0.0050,
   "originator/profile", "Joseph Holmes working-space/profile definition", "D50 (profile definition)", 0.3457, 0.3585, 1},
  {"Eureka RGB", 0.6915, 0.3083, 0.0000, 1.0000, 0.1440, 0.0296,
   "historical/secondary", "retained from Kang/Pascale", NULL, 0.0, 0.0, 0},
  {"Extended RGB", 0.7010, 0.2990, 0.1700, 0.7960, 0.1310, 0.0460,
   "historical/secondary", "retained from Kang/Pascale", NULL, 0.0, 0.0, 0},
  {"Guild RGB", 0.7000, 0.3000, 0.2550, 0.7200, 0.1500, 0.0500,
   "historical/secondary", "retained from Kang/Pascale", NULL, 0.0, 0.0, 0},
  {"Ink-jet RGB", 0.7000, 0.3000, 0.2500, 0.7200, 0.1300, 0.0500,
   "historical/secondary", "retained from Kang/Pascale", NULL, 0.0, 0.0, 0},
  {"Judd-Wyszecki RGB", 0.7347, 0.2653, 0.0743, 0.8338, 0.1741, 0.0050,
   "historical/secondary", "retained from Kang/Pascale; not a modern RGB working-space standard", NULL, 0.0, 0.0, 0},
  {"Kress RGB", 0.6915, 0.3083, 0.1547, 0.8059, 0.1440, 0.0297,
   "historical/secondary", "retained from Kang; primaries correspond approximately to spectral choices", NULL, 0.0, 0.0, 0},
  {"Laser RGB (Starkweather)", 0.7117241, 0.2882321, 0.0328204, 0.8029257, 0.1632099, 0.0119374,
   "historical/derived", "633, 514 and 442 nm laser primaries evaluated from official CIE 1931 2 degree CMFs", NULL, 0.0, 0.0, 0},
  {"NTSC (1953)", 0.6700, 0.3300, 0.2100, 0.7100, 0.1400, 0.0800,
   "standard/historical", "NTSC 1953 / ITU-R BT.470 System M", "C (system definition)", 0.310, 0.316, 1},
  {"PAL / SECAM (BT.470 625-line)", 0.6400, 0.3300, 0.2900, 0.6000, 0.1500, 0.0600,
   "standard", "ITU-R BT.470 systems B/G and EBU", "D65", 0.3127, 0.3290, 1},
  {"ProPhoto RGB (ROMM primaries)", 0.7347, 0.2653, 0.1596, 0.8404, 0.0366, 0.0001,
   "originator/profile", "ProPhoto profile uses the ROMM RGB primary set", "D50", 0.3457, 0.3585, 1},
  {"ROMM RGB / RIMM RGB", 0.7347, 0.2653, 0.1596, 0.8404, 0.0366, 0.0001,
   "standard", "ISO 22028-2 ROMM RGB / I3A 7466 RIMM RGB", "D50 (encoding definition)", 0.3457, 0.3585, 1},
  {"ROM RGB", 0.8730, 0.1440, 0.1750, 0.9270, 0.0850, 0.0001,
   "historical/secondary", "retained from Kang/Pascale; no normative definition located", NULL, 0.0, 0.0, 0},
  {"SGI RGB", 0.6250, 0.3400, 0.2800, 0.5950, 0.1550, 0.0700,
   "historical/secondary", "retained from Kang/Pascale; same primary set as legacy Apple RGB", NULL, 0.0, 0.0, 0},
  {"SMPTE-C / SMPTE 170M", 0.6300, 0.3400, 0.3100, 0.5950, 0.1550, 0.0700,
   "standard", "SMPTE-C / SMPTE 170M; cross-checked in ITU-R BT.2380", "D65", 0.3127, 0.3290, 1},
  {"SMPTE 240M", 0.6300, 0.3400, 0.3100, 0.5950, 0.1550, 0.0700,
   "standard", "SMPTE 240M; cross-checked in ITU-R BT.2380", "D65", 0.3127, 0.3290, 1},
  {"Sony P-22 phosphors", 0.6250, 0.3400, 0.2800, 0.5950, 0.1550, 0.0700,
   "historical/secondary", "legacy device-phosphor coordinates retained from Kang/Pascale", NULL, 0.0, 0.0, 0},
  {"Adobe Wide Gamut RGB", 0.7347, 0.2653, 0.1152, 0.8264, 0.1566, 0.0177,
   "originator/profile", "legacy Adobe Wide Gamut RGB profile coordinates; retained rather than forced to a secondary table", "D50 (profile definition)", 0.3457, 0.3585, 1},
  {"Wright RGB", 0.7260, 0.2740, 0.1547, 0.8059, 0.1440, 0.0297,
   "historical/secondary", "retained from Kang/Pascale; historical experimental RGB system", NULL, 0.0, 0.0, 0},
  {"Usami RGB", 0.7347, 0.2653, -0.0860, 1.0860, 0.0957, -0.0314,
   "historical/secondary", "retained from Kang; imaginary green/blue primaries intentionally lie outside the spectral locus", NULL, 0.0, 0.0, 0}
};

#define N_RGB_SPACES ((int) (sizeof (rgb_spaces) / sizeof (rgb_spaces[0])))

// General illuminant/white choices.  Where current CIE spectral data are
// available, x,y values below were recomputed from the official CIE SPD data
// at its published sampling interval with the official CIE 1931 or CIE 1964
// colour-matching functions.  These
// general illuminant values are intentionally separate from the nominal white
// coordinates that define standardized RGB spaces.

static const WhitePoint white_points[] = {
  {"CIE D65, 1931 2 degree", 0.31272687, 0.32902321, "official CIE D65 1 nm SPD + CIE 1931 CMFs"},
  {"CIE D65, 1964 10 degree", 0.31382365, 0.33099899, "official CIE D65 1 nm SPD + CIE 1964 CMFs"},
  {"CIE D50, 1931 2 degree", 0.34568422, 0.35850403, "official CIE D50 1 nm SPD + CIE 1931 CMFs"},
  {"CIE D50, 1964 10 degree", 0.34774768, 0.35953602, "official CIE D50 1 nm SPD + CIE 1964 CMFs"},
  {"CIE A, 1931 2 degree", 0.44757351, 0.40743944, "official CIE A 1 nm SPD + CIE 1931 CMFs"},
  {"CIE A, 1964 10 degree", 0.45117394, 0.40593660, "official CIE A 1 nm SPD + CIE 1964 CMFs"},
  {"CIE B, 1931 2 degree (legacy)", 0.34830,0.35160, "legacy published chromaticity; CIE B is obsolete and no current CIE SPD dataset was located"},
  {"CIE C, 1931 2 degree", 0.31005847, 0.31614971, "official CIE C SPD + CIE 1931 CMFs"},
  {"CIE C, 1964 10 degree", 0.31038866, 0.31905071, "official CIE C SPD + CIE 1964 CMFs"},
  {"CIE D55, 1931 2 degree", 0.33242410, 0.34742804, "official CIE D55 SPD + CIE 1931 CMFs"},
  {"CIE D55, 1964 10 degree", 0.33411634, 0.34876609, "official CIE D55 SPD + CIE 1964 CMFs"},
  {"ACES white point (D60-like), 1931 2 degree", 0.32168, 0.33767, "SMPTE ST 2065-1 / Academy ACES definition"},
  {"CIE D75, 1931 2 degree", 0.29902230, 0.31485274, "official CIE D75 SPD + CIE 1931 CMFs"},
  {"CIE D75, 1964 10 degree", 0.29967997, 0.31740324, "official CIE D75 SPD + CIE 1964 CMFs"},
  {"CIE daylight 9300 K (D93), 1931 2 degree", 0.28314501, 0.29711289, "CIE daylight-locus formula at 9300 K"},
  {"CIE daylight 9300 K (D93), 1964 10 degree", 0.28325, 0.30040, "CIE daylight components evaluated with CIE 1964 CMFs; rounded to avoid false precision"},
  {"Equal-energy E, 1931 2 degree", 1.0 / 3.0, 1.0 / 3.0, "mathematical equal-energy white"},
  {"Equal-energy E, 1964 10 degree", 1.0 / 3.0, 1.0 / 3.0, "mathematical equal-energy white"},
  {"CIE FL1, 1931 2 degree", 0.31306243, 0.33710648, "CIE 015:2018 FL1 SPD + CIE 1931 CMFs"},
  {"CIE FL1, 1964 10 degree", 0.31809880, 0.33548945, "CIE 015:2018 FL1 SPD + CIE 1964 CMFs"},
  {"CIE FL2, 1931 2 degree", 0.37206815, 0.37512256, "CIE 015:2018 FL2 SPD + CIE 1931 CMFs"},
  {"CIE FL2, 1964 10 degree", 0.37927483, 0.36722793, "CIE 015:2018 FL2 SPD + CIE 1964 CMFs"},
  {"CIE FL3, 1931 2 degree", 0.40909004, 0.39411713, "CIE 015:2018 FL3 SPD + CIE 1931 CMFs"},
  {"CIE FL3, 1964 10 degree", 0.41764468, 0.38312450, "CIE 015:2018 FL3 SPD + CIE 1964 CMFs"},
  {"CIE FL4, 1931 2 degree", 0.44018110, 0.40309069, "CIE 015:2018 FL4 SPD + CIE 1931 CMFs"},
  {"CIE FL4, 1964 10 degree", 0.44924770, 0.39060548, "CIE 015:2018 FL4 SPD + CIE 1964 CMFs"},
  {"CIE FL5, 1931 2 degree", 0.31375735, 0.34516065, "CIE 015:2018 FL5 SPD + CIE 1931 CMFs"},
  {"CIE FL5, 1964 10 degree", 0.31974054, 0.34236696, "CIE 015:2018 FL5 SPD + CIE 1964 CMFs"},
  {"CIE FL6, 1931 2 degree", 0.37787777, 0.38819415, "CIE 015:2018 FL6 SPD + CIE 1931 CMFs"},
  {"CIE FL6, 1964 10 degree", 0.38662283, 0.37837312, "CIE 015:2018 FL6 SPD + CIE 1964 CMFs"},
  {"CIE FL7, 1931 2 degree", 0.31285247, 0.32917418, "CIE 015:2018 FL7 SPD + CIE 1931 CMFs"},
  {"CIE FL7, 1964 10 degree", 0.31564564, 0.32950815, "CIE 015:2018 FL7 SPD + CIE 1964 CMFs"},
  {"CIE FL8, 1931 2 degree", 0.34580575, 0.35861758, "CIE 015:2018 FL8 SPD + CIE 1931 CMFs"},
  {"CIE FL8, 1964 10 degree", 0.34896556, 0.35931730, "CIE 015:2018 FL8 SPD + CIE 1964 CMFs"},
  {"CIE FL9, 1931 2 degree", 0.37409927, 0.37268420, "CIE 015:2018 FL9 SPD + CIE 1931 CMFs"},
  {"CIE FL9, 1964 10 degree", 0.37825426, 0.37038210, "CIE 015:2018 FL9 SPD + CIE 1964 CMFs"},
  {"CIE FL10, 1931 2 degree", 0.34578790, 0.35875793, "CIE 015:2018 FL10 SPD + CIE 1931 CMFs"},
  {"CIE FL10, 1964 10 degree", 0.35061017, 0.35430334, "CIE 015:2018 FL10 SPD + CIE 1964 CMFs"},
  {"CIE FL11, 1931 2 degree", 0.38053749, 0.37691531, "CIE 015:2018 FL11 SPD + CIE 1931 CMFs"},
  {"CIE FL11, 1964 10 degree", 0.38543539, 0.37109479, "CIE 015:2018 FL11 SPD + CIE 1964 CMFs"},
  {"CIE FL12, 1931 2 degree", 0.43702434, 0.40421500, "CIE 015:2018 FL12 SPD + CIE 1931 CMFs"},
  {"CIE FL12, 1964 10 degree", 0.44265468, 0.39706106, "CIE 015:2018 FL12 SPD + CIE 1964 CMFs"}
};

#define N_WHITE_POINTS ((int) (sizeof (white_points) / sizeof (white_points[0])))

int
main (void) {

  int i, j, rgb_index;
  double **p, *w, **pinv, *coeff, **c, **npm, ** npminv;

  fprintf (stdout, "\nDerivation of RGB/XYZ Conversion Matrices\n");

  // Allocate memory for various arrays.
  p = allocate_doublememp (3);
  pinv = allocate_doublememp (3);
  c = allocate_doublememp (3);
  npm = allocate_doublememp (3);
  npminv = allocate_doublememp (3);
  for (i=0; i<3; i++) {
    p[i] = allocate_doublemem (3);
    pinv[i] = allocate_doublemem (3);
    c[i] = allocate_doublemem (3);
    npm[i] = allocate_doublemem (3);
    npminv[i] = allocate_doublemem (3);
  }
  w = allocate_doublemem (3);
  coeff = allocate_doublemem (3);

  // XYZ CIE primary coordinates
  // Xx = 1, Xy = 0, Xz = 0
  // Yx = 0, Yy = 1, Yz = 0
  // Zx = 0, Zy = 0, Zz = 1

  // Choose RGB color primaries (i.e., RGB colorspace).
  for (;;) {
    if (rgb_primaries (p, &rgb_index) > -1) break;
  }

  // Choose white coordinates.
  for (;;) {
    if (illum_white (w, rgb_index) > -1) break;
  }

  // Show color primaries matrix p.
  fprintf (stdout, "\nColor primaries matrix p:\n");
  for (i=0; i<3; i++) {
    fprintf (stdout, "  ");
    for (j=0; j<3; j++) {
      fprintf (stdout, "%0.4lf ", p[i][j]);
      pinv[i][j] = p[i][j];  // Copy matrix p to matrix pinv for later in-place inversion.
    }
    fprintf (stdout, "\n");
  }

  // Show selected white chromaticity coordinates.
  fprintf (stdout, "\nWhite chromaticity coordinates (x y z):\n");
  fprintf (stdout, "%0.4lf %0.4lf %0.4lf\n", w[0], w[1], w[2]);

  // Populate white matrix.
  // Normalize to luminosity y.
  w[0] /= w[1];
  w[2] /= w[1];
  w[1] = 1.0;

  // Compute inverse of color primaries matrix p.
  if (gaussjordan (3, pinv) == EXIT_FAILURE) {
    fprintf (stderr, "ERROR: Unable to invert color primaries matrix.\n");
    exit (EXIT_FAILURE);
  }

  // Show inverse of color primaries matrix.
  fprintf (stdout, "\nInverse (pinv) of color primaries matrix:\n");
  for (i=0; i<3; i++) {
    fprintf (stdout, "  ");
    for (j=0; j<3; j++) {
      fprintf (stdout, "%0.4lf ", pinv[i][j]);
    }
    fprintf (stdout, "\n");
  }
  fprintf (stdout, "\n");

// Test to see if p * pinv = identity vector.
/*
double v;
  for (i=0; i<3; i++) {
    v = 0.0;
    for (j=0; j<3; j++) {
      v += p[i][j] * pinv[j][i];
    }
    fprintf (stdout, "%0.4lf\n", v);
  }
*/

  // Calculate RGB normalization coefficients.
  // coef = pinv * w
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      coeff[i] += pinv[i][j] * w[j];
    }
  }

  // Show normalization coefficients.
  fprintf (stdout, "Normalization coefficient vector:\n");
  for (i=0; i<3; i++) {
    fprintf (stdout, "  %0.4lf\n", coeff[i]);
  }
  fprintf (stdout, "\n");

  // Create diagonal coefficient matrix using normalization coefficients.
  c[0][0] = coeff[0];
  c[1][1] = coeff[1];
  c[2][2] = coeff[2];

  // Show diagonal coefficient matrix.
  fprintf (stdout, "Diagonal normalization coefficient matrix c:\n");
  for (i=0; i<3; i++) {
    fprintf (stdout, "  ");
    for (j=0; j<3; j++) {
      fprintf (stdout, "%0.4lf ", c[i][j]);
    }
    fprintf (stdout, "\n");
  }
  fprintf (stdout, "\n");

  // Compute NPM = P * C. Since C is diagonal, each primary column is
  // simply scaled by its corresponding normalization coefficient.
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      npm[i][j] = p[i][j] * coeff[j];
    }
  }

  // Copy npm to npminv so we can invert in-place.
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      npminv[i][j] = npm[i][j];
    }
  }
  if (gaussjordan (3, npminv) == EXIT_FAILURE) {
    fprintf (stderr, "ERROR: Unable to invert normalized primary matrix.\n");
    exit (EXIT_FAILURE);
  }

  fprintf (stdout, " --- Results ---\n\n");

  // Show inverse NPM matrix.
  fprintf (stdout, "Tristimulus conversion XYZ to RGB using inverse of Normalized Primary Matrix (npminv):\n");
  fprintf (stdout, "(R G B) = npminv * (X Y Z)\n");
  for (i=0; i<3; i++) {
    fprintf (stdout, "  ");
    for (j=0; j<3; j++) {
      fprintf (stdout, "%0.7lf ", npminv[i][j]);
    }
    fprintf (stdout, "\n");
  }
  fprintf (stdout, "\n");

  // Show NPM matrix.
  fprintf (stdout, "Tristimulus conversion RGB to XYZ using Normalized Primary Matrix (npm):\n");
  fprintf (stdout, "(X Y Z) = npm * (R G B)\n");
  for (i=0; i<3; i++) {
    fprintf (stdout, "  ");
    for (j=0; j<3; j++) {
      fprintf (stdout, "%0.7lf ", npm[i][j]);
    }
    fprintf (stdout, "\n");
  }
  fprintf (stdout, "\n");

  // Show rgb to xyz chromaticity conversion equations.
  fprintf (stdout, "rgb to xyz chromaticity conversion equations:\n");
  fprintf (stdout, "x = (%.7lfr + %.7lfg + %.7lfb) / (%.7lfr + %.7lfg + %.7lfb)\n", npm[0][0], npm[0][1], npm[0][2], npm[0][0] + npm[1][0] + npm[2][0], npm[0][1] + npm[1][1] + npm[2][1], npm[0][2] + npm[1][2] + npm[2][2]);
  fprintf (stdout, "y = (%.7lfr + %.7lfg + %.7lfb) / (%.7lfr + %.7lfg + %.7lfb)\n", npm[1][0], npm[1][1], npm[1][2], npm[0][0] + npm[1][0] + npm[2][0], npm[0][1] + npm[1][1] + npm[2][1], npm[0][2] + npm[1][2] + npm[2][2]);
  fprintf (stdout, "z = (%.7lfr + %.7lfg + %.7lfb) / (%.7lfr + %.7lfg + %.7lfb)\n", npm[2][0], npm[2][1], npm[2][2], npm[0][0] + npm[1][0] + npm[2][0], npm[0][1] + npm[1][1] + npm[2][1], npm[0][2] + npm[1][2] + npm[2][2]);

  // Show xyz to rgb chromaticity conversion equations.
  fprintf (stdout, "\nxyz to rgb chromaticity conversion equations:\n");
  fprintf (stdout, "r = (%.7lfx + %.7lfy + %.7lfz) / (%.7lfx + %.7lfy + %.7lfz)\n", npminv[0][0], npminv[0][1], npminv[0][2], npminv[0][0] + npminv[1][0] + npminv[2][0], npminv[0][1] + npminv[1][1] + npminv[2][1], npminv[0][2] + npminv[1][2] + npminv[2][2]);
  fprintf (stdout, "g = (%.7lfx + %.7lfy + %.7lfz) / (%.7lfx + %.7lfy + %.7lfz)\n", npminv[1][0], npminv[1][1], npminv[1][2], npminv[0][0] + npminv[1][0] + npminv[2][0], npminv[0][1] + npminv[1][1] + npminv[2][1], npminv[0][2] + npminv[1][2] + npminv[2][2]);
  fprintf (stdout, "b = (%.7lfx + %.7lfy + %.7lfz) / (%.7lfx + %.7lfy + %.7lfz)\n", npminv[2][0], npminv[2][1], npminv[2][2], npminv[0][0] + npminv[1][0] + npminv[2][0], npminv[0][1] + npminv[1][1] + npminv[2][1], npminv[0][2] + npminv[1][2] + npminv[2][2]);

  // Show luminance equation constants (second row of NPM matrix).
  fprintf (stdout, "\nLuminance equation constants YR, YG, and YB (second row of NPM matrix) (sometimes referred to as KR, KG, and KB)\n");
  fprintf (stdout,"Y = YR * R + YG * G + YB * B\n");
  fprintf (stdout, "  YR: %0.10lf\n", npm[1][0]);
  fprintf (stdout, "  YG: %0.10lf\n", npm[1][1]);
  fprintf (stdout, "  YB: %0.10lf\n", npm[1][2]);
  fprintf (stdout, "  SUM: YR + YG + YB = %0.10lf\n\n", npm[1][0] + npm[1][1] + npm[1][2]);

  // Free allocated memory.
  for (i=0; i<3; i++) {
    free (p[i]);
    free (pinv[i]);
    free (c[i]);
    free (npm[i]);
    free (npminv[i]);
  }
  free (p);
  free (pinv);
  free (npminv);
  free (w);
  free (coeff);
  free (c);
  free (npm);

  return (EXIT_SUCCESS);
}

// Gauss-Jordan in-place inversion of n*n matrix.
// Partial pivoting is used so an invertible matrix is not rejected merely
// because the current diagonal element is zero or very small.
int
gaussjordan (int n, double **matrix) {

  int i, j, k, pivot_row;
  double **augmented, pivot, factor, max_abs;
  double *rowtmp;

  if ((n <= 0) || (matrix == NULL)) return (EXIT_FAILURE);

  augmented = allocate_doublememp ((size_t) n);
  for (i=0; i<n; i++) {
    if (matrix[i] == NULL) {
      for (j=0; j<i; j++) free (augmented[j]);
      free (augmented);
      return (EXIT_FAILURE);
    }
    augmented[i] = allocate_doublemem ((size_t) (2 * n));
  }

  // Augment the matrix with the identity matrix.
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      augmented[i][j] = matrix[i][j];
      augmented[i][j + n] = (i == j) ? 1.0 : 0.0;
    }
  }

  // Perform Gauss-Jordan elimination with partial pivoting.
  for (i=0; i<n; i++) {
    pivot_row = i;
    max_abs = fabs (augmented[i][i]);
    for (j=i + 1; j<n; j++) {
      if (fabs (augmented[j][i]) > max_abs) {
        max_abs = fabs (augmented[j][i]);
        pivot_row = j;
      }
    }

    if (!isfinite (max_abs) || (max_abs <= DBL_EPSILON * 1024.0)) {
      for (j=0; j<n; j++) free (augmented[j]);
      free (augmented);
      return (EXIT_FAILURE);
    }

    if (pivot_row != i) {
      rowtmp = augmented[i];
      augmented[i] = augmented[pivot_row];
      augmented[pivot_row] = rowtmp;
    }

    pivot = augmented[i][i];
    for (j=0; j<(2 * n); j++) {
      augmented[i][j] /= pivot;
    }

    for (j=0; j<n; j++) {
      if (j == i) continue;
      factor = augmented[j][i];
      for (k=0; k<(2 * n); k++) {
        augmented[j][k] -= factor * augmented[i][k];
      }
    }
  }

  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      matrix[i][j] = augmented[i][j + n];
    }
  }

  for (i=0; i<n; i++) free (augmented[i]);
  free (augmented);

  return (EXIT_SUCCESS);
}

// Obtain a text string from standard input. It can include spaces.
int
inputtext (char *text) {

  int ch;
  size_t len;

  if (fgets (text, MAX_STRINGLEN, stdin) == NULL) {
    fprintf (stderr, "Unable to read text from standard input.\n");
    exit (EXIT_FAILURE);
  }

  len = strlen (text);

  // Remove trailing newline, and a preceding carriage return if present.
  if ((len > 0) && (text[len - 1] == '\n')) {
    text[--len] = '\0';
    if ((len > 0) && (text[len - 1] == '\r')) {
      text[--len] = '\0';
    }
    return (EXIT_SUCCESS);
  }

  // If the buffer is full, determine whether the input was exactly
  // MAX_STRINGLEN - 1 characters or was genuinely too long.
  if (len == MAX_STRINGLEN - 1) {

    ch = getchar ();

    // Exactly MAX_STRINGLEN - 1 characters followed by newline or EOF.
    if ((ch == '\n') || (ch == EOF)) {
      return (EXIT_SUCCESS);
    }

    // Handle CRLF after an exactly full input line.
    if (ch == '\r') {
      ch = getchar ();
      if ((ch == '\n') || (ch == EOF)) {
        return (EXIT_SUCCESS);
      }
    }

    // Discard the remainder of an overlong input line.
    while ((ch != '\n') && (ch != EOF)) {
      ch = getchar ();
    }

    fprintf (stderr, "Input text is too long; maximum is %d characters.\n",
             MAX_STRINGLEN - 1);
    exit (EXIT_FAILURE);
  }

  return (EXIT_SUCCESS);
}

// Convert a complete input string to int, allowing surrounding whitespace only.
int
parse_int_string (const char *text, int *value) {

  char *endptr;
  long parsed;

  if ((text == NULL) || (value == NULL)) return (EXIT_FAILURE);

  errno = 0;
  parsed = strtol (text, &endptr, 10);
  if ((errno == ERANGE) || (endptr == text)) return (EXIT_FAILURE);
  while (isspace ((unsigned char) *endptr)) endptr++;
  if ((*endptr != '\0') || (parsed < INT_MIN) || (parsed > INT_MAX)) {
    return (EXIT_FAILURE);
  }

  *value = (int) parsed;
  return (EXIT_SUCCESS);
}

// Color primary coordinates for RGB spaces and historical RGB primary sets.
// Returns: -1 if invalid selection, 0 if valid selection.
int
rgb_primaries (double **p, int *rgb_index) {

  int choice, i;
  double zr, zg, zb;
  char *temp;
  const RGBSpace *space;

  temp = allocate_strmem (MAX_STRINGLEN);

  fprintf (stdout, "\nChoose RGB color space / primary set:\n");
  for (i=0; i<N_RGB_SPACES; i++) {
    fprintf (stdout, " %2d - %s [%s]\n", i + 1, rgb_spaces[i].name, rgb_spaces[i].source_quality);
  }

  fprintf (stdout, "\nChoice? ");
  inputtext (temp);
  if (parse_int_string (temp, &choice) == EXIT_FAILURE) {
    fprintf (stderr, "ERROR: Cannot make integer of: %s\n", temp);
    free (temp);
    return (-1);
  }
  if ((choice < 1) || (choice > N_RGB_SPACES)) {
    fprintf (stderr, "Invalid choice.\n");
    free (temp);
    return (-1);
  }

  *rgb_index = choice - 1;
  space = &rgb_spaces[*rgb_index];

  zr = 1.0 - space->xr - space->yr;
  zg = 1.0 - space->xg - space->yg;
  zb = 1.0 - space->xb - space->yb;

  p[0][0] = space->xr;  p[0][1] = space->xg;  p[0][2] = space->xb;
  p[1][0] = space->yr;  p[1][1] = space->yg;  p[1][2] = space->yb;
  p[2][0] = zr;         p[2][1] = zg;         p[2][2] = zb;

  fprintf (stdout, "\nSelected primary set: %s\n", space->name);
  fprintf (stdout, "Data basis: %s\n", space->source_note);
  if (space->native_white_known) {
    fprintf (stdout, "Suggested native/reference white: %s, x = %.8lf, y = %.8lf\n",
             space->native_white_name, space->xw, space->yw);
  } else {
    fprintf (stdout, "Suggested native/reference white: none established by the source used.\n");
    fprintf (stdout, "Choose a white deliberately; changing the white changes the normalized primary matrix.\n");
  }

  free (temp);
  return (0);
}

// Illuminants - Choose white point.
// Choice 0, when available, uses the selected RGB space's own defining/native
// white coordinates.  The numbered general illuminant choices are independent
// colorimetric values and may contain more digits than a standardized RGB
// encoding's nominal reference white.
// Returns: -1 if invalid selection, 0 if valid selection.
int
illum_white (double *white_xyz, int rgb_index) {

  int choice, i;
  char *temp;
  const RGBSpace *space;
  const WhitePoint *wp;

  if ((rgb_index < 0) || (rgb_index >= N_RGB_SPACES)) {
    fprintf (stderr, "ERROR: Invalid RGB-space index in illum_white().\n");
    return (-1);
  }

  space = &rgb_spaces[rgb_index];
  temp = allocate_strmem (MAX_STRINGLEN);

  fprintf (stdout, "\nChoose white point:\n");
  if (space->native_white_known) {
    fprintf (stdout, "  0 - RECOMMENDED for %s: %s (x=%.8lf, y=%.8lf)\n",
             space->name, space->native_white_name, space->xw, space->yw);
  } else {
    fprintf (stdout, "  0 - unavailable: this primary set has no established native/reference white in the source used\n");
  }

  for (i=0; i<N_WHITE_POINTS; i++) {
    fprintf (stdout, " %2d - %s (x=%.8lf, y=%.8lf)\n",
             i + 1, white_points[i].name, white_points[i].x, white_points[i].y);
  }

  fprintf (stdout, "\nChoice? ");
  inputtext (temp);
  if (parse_int_string (temp, &choice) == EXIT_FAILURE) {
    fprintf (stderr, "ERROR: Cannot make integer of: %s\n", temp);
    free (temp);
    return (-1);
  }

  if (choice == 0) {
    if (!space->native_white_known) {
      fprintf (stderr, "No native/reference white is established for %s in the source used.\n", space->name);
      free (temp);
      return (-1);
    }
    white_xyz[0] = space->xw;
    white_xyz[1] = space->yw;
    fprintf (stdout, "Using defining/native white for %s: %s.\n", space->name, space->native_white_name);
  } else if ((choice >= 1) && (choice <= N_WHITE_POINTS)) {
    wp = &white_points[choice - 1];
    white_xyz[0] = wp->x;
    white_xyz[1] = wp->y;
    fprintf (stdout, "Using %s.\nData basis: %s\n", wp->name, wp->source_note);
  } else {
    fprintf (stderr, "Invalid selection.\n");
    free (temp);
    return (-1);
  }

  white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];

  free (temp);
  return (0);
}

// Allocate memory for an array of doubles.
double *
allocate_doublemem (size_t len) {

  void *tmp;

  if (len == 0) { 
    fprintf (stderr, "ERROR: Cannot allocate memory because len = %zu in allocate_doublemem().\n", len);
    exit (EXIT_FAILURE);
  }

  tmp = calloc (len, sizeof (double));
  if (tmp != NULL) {
    return (tmp);
  } else {
    fprintf (stderr, "ERROR: Cannot allocate memory for array in allocate_doublemem().\n");
    exit (EXIT_FAILURE);
  }
} 
  
// Allocate memory for an array of pointers to arrays of doubles.
double **
allocate_doublememp (size_t len) {

  void *tmp;

  if (len == 0) {
    fprintf (stderr, "ERROR: Cannot allocate memory because len = %zu in allocate_doublememp().\n", len);
    exit (EXIT_FAILURE);
  }

  tmp = calloc (len, sizeof (double *));
  if (tmp != NULL) {
    return (tmp);
  } else {
    fprintf (stderr, "ERROR: Cannot allocate memory for array in allocate_doublememp().\n");
    exit (EXIT_FAILURE);
  } 
}

// Allocate memory for an array of chars.
char *
allocate_strmem (size_t len) {

  void *tmp;

  if (len == 0) {
    fprintf (stderr, "ERROR: Cannot allocate memory because len = %zu in allocate_strmem().\n", len);
    exit (EXIT_FAILURE);
  }

  tmp = calloc (len, sizeof (char));
  if (tmp != NULL) {
    return (tmp);
  } else {
    fprintf (stderr, "ERROR: Cannot allocate memory for array in allocate_strmem().\n");
    exit (EXIT_FAILURE);
  }
}
