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

// locus.c - Produce the chromaticity locus in an sRGB bitmap using selected color-matching function (CMF).

// gcc -Wall src/locus.c -lm -o locus

// Usage: ./locus
// Inputs: CMF files:
//         data/CIE_xyz_1931_2deg.csv
//         data/CIE_xyz_1931_2deg_judd1951.csv
//         data/CIE_xyz_1931_2deg_judd1951_vos1978.csv
//         data/CIE_xyz_1964_10deg.csv
//         data/CIE_xyz_2006_10deg_lms_cones.csv
//         data/CIE_xyz_2006_2deg_lms_cones.csv
// Output: out.bmp

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>  // uint8_t
#include <string.h>
#include <math.h>  // fabs()
#include <errno.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>

typedef struct {
    uint8_t blue;   // Blue
    uint8_t green;  // Green
    uint8_t red;    // Red
} PIXEL;

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

/*
 * Primary coordinates are kept at the precision at which the defining standard,
 * originating profile, or historical source specifies them.  Do not add digits
 * merely by substituting a more precise chromaticity for a named illuminant: for
 * standardized RGB spaces, the reference-white coordinates are part of the RGB
 * space definition.
 *
 * "standard" means a normative standard/registry was available; "originator/profile"
 * means an originating author or widely distributed defining profile was available;
 * "historical/secondary" means no normative definition was located and the values
 * are retained from the Kang/Pascale source set used by the original program.
 */
static const RGBSpace rgb_spaces[] = {
  {"sRGB / ITU-R BT.709", 0.6400,0.3300, 0.3000,0.6000, 0.1500,0.0600,
   "standard", "IEC 61966-2-1 / ITU-R BT.709", "D65 (space definition)", 0.3127,0.3290, 1},
  {"Adobe RGB (1998)", 0.6400,0.3300, 0.2100,0.7100, 0.1500,0.0600,
   "standard", "Adobe RGB (1998) Color Image Encoding", "D65 (space definition)", 0.3127,0.3290, 1},
  {"Apple RGB", 0.6250,0.3400, 0.2800,0.5950, 0.1550,0.0700,
   "originator/profile", "legacy Apple/ColorSync profile definition; cross-checked against Adobe profile set", "D65 (legacy profile)", 0.3127,0.3290, 1},
  {"Best RGB", 0.7347,0.2653, 0.2150,0.7750, 0.1300,0.0350,
   "originator/profile", "Bruce Lindbloom working-space definition", "D50 (working-space definition)", 0.3457,0.3585, 1},
  {"Beta RGB", 0.6888,0.3112, 0.1986,0.7551, 0.1265,0.0352,
   "originator/profile", "Bruce Lindbloom working-space definition", "D50 (working-space definition)", 0.3457,0.3585, 1},
  {"Bruce RGB", 0.6400,0.3300, 0.2800,0.6500, 0.1500,0.0600,
   "originator/profile", "Bruce Lindbloom working-space definition", "D65 (working-space definition)", 0.3127,0.3290, 1},
  {"CIE 1931 RGB (2 degree observer)",
   0.7346900,0.2653100, 0.2736827,0.7174214, 0.1665347,0.0088840,
   "historical/derived", "700.0, 546.1 and 435.8 nm primaries evaluated from official CIE 1931 2 degree CMFs", "equal-energy E", 1.0/3.0,1.0/3.0, 1},
  {"CIE 1964/RGB (legacy secondary table)", 0.7232,0.2768, 0.1248,0.8216, 0.1616,0.0134,
   "historical/secondary", "retained from Kang; CIE standardized a 10 degree XYZ observer, not this RGB working space", "equal-energy E (historical convention)", 1.0/3.0,1.0/3.0, 1},
  {"ColorMatch RGB", 0.6300,0.3400, 0.2950,0.6050, 0.1500,0.0750,
   "originator/profile", "legacy ColorMatch profile definition", "D50 (profile definition)", 0.3457,0.3585, 1},
  {"Don RGB 4", 0.6960,0.3000, 0.2150,0.7650, 0.1300,0.0350,
   "originator/profile", "Don Hutcheson / HutchColor working-space definition", "D50 (working-space definition)", 0.3457,0.3585, 1},
  {"EBU Tech. 3213-E", 0.6400,0.3300, 0.2900,0.6000, 0.1500,0.0600,
   "standard", "EBU Tech. 3213-E / ITU-R BT.470 625-line systems", "D65", 0.3127,0.3290, 1},
  {"eciRGB v2", 0.6700,0.3300, 0.2100,0.7100, 0.1400,0.0800,
   "standard", "ECI / ISO 22028-2 registry definition", "D50 (space definition)", 0.3457,0.3585, 1},
  {"Ekta Space PS5", 0.6950,0.3050, 0.2600,0.7000, 0.1100,0.0050,
   "originator/profile", "Joseph Holmes working-space/profile definition", "D50 (profile definition)", 0.3457,0.3585, 1},
  {"Eureka RGB", 0.6915,0.3083, 0.0000,1.0000, 0.1440,0.0296,
   "historical/secondary", "retained from Kang/Pascale", NULL, 0.0,0.0, 0},
  {"Extended RGB", 0.7010,0.2990, 0.1700,0.7960, 0.1310,0.0460,
   "historical/secondary", "retained from Kang/Pascale", NULL, 0.0,0.0, 0},
  {"Guild RGB", 0.7000,0.3000, 0.2550,0.7200, 0.1500,0.0500,
   "historical/secondary", "retained from Kang/Pascale", NULL, 0.0,0.0, 0},
  {"Ink-jet RGB", 0.7000,0.3000, 0.2500,0.7200, 0.1300,0.0500,
   "historical/secondary", "retained from Kang/Pascale", NULL, 0.0,0.0, 0},
  {"Judd-Wyszecki RGB", 0.7347,0.2653, 0.0743,0.8338, 0.1741,0.0050,
   "historical/secondary", "retained from Kang/Pascale; not a modern RGB working-space standard", NULL, 0.0,0.0, 0},
  {"Kress RGB", 0.6915,0.3083, 0.1547,0.8059, 0.1440,0.0297,
   "historical/secondary", "retained from Kang; primaries correspond approximately to spectral choices", NULL, 0.0,0.0, 0},
  {"Laser RGB (Starkweather)", 0.7117241,0.2882321, 0.0328204,0.8029257, 0.1632099,0.0119374,
   "historical/derived", "633, 514 and 442 nm laser primaries evaluated from official CIE 1931 2 degree CMFs", NULL, 0.0,0.0, 0},
  {"NTSC (1953)", 0.6700,0.3300, 0.2100,0.7100, 0.1400,0.0800,
   "standard/historical", "NTSC 1953 / ITU-R BT.470 System M", "C (system definition)", 0.310,0.316, 1},
  {"PAL / SECAM (BT.470 625-line)", 0.6400,0.3300, 0.2900,0.6000, 0.1500,0.0600,
   "standard", "ITU-R BT.470 systems B/G and EBU", "D65", 0.3127,0.3290, 1},
  {"ProPhoto RGB (ROMM primaries)", 0.7347,0.2653, 0.1596,0.8404, 0.0366,0.0001,
   "originator/profile", "ProPhoto profile uses the ROMM RGB primary set", "D50", 0.3457,0.3585, 1},
  {"ROMM RGB / RIMM RGB", 0.7347,0.2653, 0.1596,0.8404, 0.0366,0.0001,
   "standard", "ISO 22028-2 ROMM RGB / I3A 7466 RIMM RGB", "D50 (encoding definition)", 0.3457,0.3585, 1},
  {"ROM RGB", 0.8730,0.1440, 0.1750,0.9270, 0.0850,0.0001,
   "historical/secondary", "retained from Kang/Pascale; no normative definition located", NULL, 0.0,0.0, 0},
  {"SGI RGB", 0.6250,0.3400, 0.2800,0.5950, 0.1550,0.0700,
   "historical/secondary", "retained from Kang/Pascale; same primary set as legacy Apple RGB", NULL, 0.0,0.0, 0},
  {"SMPTE-C / SMPTE 170M", 0.6300,0.3400, 0.3100,0.5950, 0.1550,0.0700,
   "standard", "SMPTE-C / SMPTE 170M; cross-checked in ITU-R BT.2380", "D65", 0.3127,0.3290, 1},
  {"SMPTE 240M", 0.6300,0.3400, 0.3100,0.5950, 0.1550,0.0700,
   "standard", "SMPTE 240M; cross-checked in ITU-R BT.2380", "D65", 0.3127,0.3290, 1},
  {"Sony P-22 phosphors", 0.6250,0.3400, 0.2800,0.5950, 0.1550,0.0700,
   "historical/secondary", "legacy device-phosphor coordinates retained from Kang/Pascale", NULL, 0.0,0.0, 0},
  {"Adobe Wide Gamut RGB", 0.7347,0.2653, 0.1152,0.8264, 0.1566,0.0177,
   "originator/profile", "legacy Adobe Wide Gamut RGB profile coordinates; retained rather than forced to a secondary table", "D50 (profile definition)", 0.3457,0.3585, 1},
  {"Wright RGB", 0.7260,0.2740, 0.1547,0.8059, 0.1440,0.0297,
   "historical/secondary", "retained from Kang/Pascale; historical experimental RGB system", NULL, 0.0,0.0, 0},
  {"Usami RGB", 0.7347,0.2653, -0.0860,1.0860, 0.0957,-0.0314,
   "historical/secondary", "retained from Kang; imaginary green/blue primaries intentionally lie outside the spectral locus", NULL, 0.0,0.0, 0}
};

#define N_RGB_SPACES ((int) (sizeof (rgb_spaces) / sizeof (rgb_spaces[0])))

/*
 * General illuminant/white choices.  Where current CIE spectral data are
 * available, x,y values below were recomputed from the official CIE SPD data
 * at its published sampling interval with the official CIE 1931 or CIE 1964
 * colour-matching functions.  These
 * general illuminant values are intentionally separate from the nominal white
 * coordinates that define standardized RGB spaces.
 */
static const WhitePoint white_points[] = {
  {"CIE D65, 1931 2 degree", 0.31272687,0.32902321, "official CIE D65 1 nm SPD + CIE 1931 CMFs"},
  {"CIE D65, 1964 10 degree", 0.31382365,0.33099899, "official CIE D65 1 nm SPD + CIE 1964 CMFs"},
  {"CIE D50, 1931 2 degree", 0.34568422,0.35850403, "official CIE D50 1 nm SPD + CIE 1931 CMFs"},
  {"CIE D50, 1964 10 degree", 0.34774768,0.35953602, "official CIE D50 1 nm SPD + CIE 1964 CMFs"},
  {"CIE A, 1931 2 degree", 0.44757351,0.40743944, "official CIE A 1 nm SPD + CIE 1931 CMFs"},
  {"CIE A, 1964 10 degree", 0.45117394,0.40593660, "official CIE A 1 nm SPD + CIE 1964 CMFs"},
  {"CIE B, 1931 2 degree (legacy)", 0.34830,0.35160, "legacy published chromaticity; CIE B is obsolete and no current CIE SPD dataset was located"},
  {"CIE C, 1931 2 degree", 0.31005847,0.31614971, "official CIE C SPD + CIE 1931 CMFs"},
  {"CIE C, 1964 10 degree", 0.31038866,0.31905071, "official CIE C SPD + CIE 1964 CMFs"},
  {"CIE D55, 1931 2 degree", 0.33242410,0.34742804, "official CIE D55 SPD + CIE 1931 CMFs"},
  {"CIE D55, 1964 10 degree", 0.33411634,0.34876609, "official CIE D55 SPD + CIE 1964 CMFs"},
  {"ACES white point (D60-like), 1931 2 degree", 0.32168,0.33767, "SMPTE ST 2065-1 / Academy ACES definition"},
  {"CIE D75, 1931 2 degree", 0.29902230,0.31485274, "official CIE D75 SPD + CIE 1931 CMFs"},
  {"CIE D75, 1964 10 degree", 0.29967997,0.31740324, "official CIE D75 SPD + CIE 1964 CMFs"},
  {"CIE daylight 9300 K (D93), 1931 2 degree", 0.28314501,0.29711289, "CIE daylight-locus formula at 9300 K"},
  {"CIE daylight 9300 K (D93), 1964 10 degree", 0.28325,0.30040, "CIE daylight components evaluated with CIE 1964 CMFs; rounded to avoid false precision"},
  {"Equal-energy E, 1931 2 degree", 1.0/3.0,1.0/3.0, "mathematical equal-energy white"},
  {"Equal-energy E, 1964 10 degree", 1.0/3.0,1.0/3.0, "mathematical equal-energy white"},
  {"CIE FL1, 1931 2 degree", 0.31306243,0.33710648, "CIE 015:2018 FL1 SPD + CIE 1931 CMFs"},
  {"CIE FL1, 1964 10 degree", 0.31809880,0.33548945, "CIE 015:2018 FL1 SPD + CIE 1964 CMFs"},
  {"CIE FL2, 1931 2 degree", 0.37206815,0.37512256, "CIE 015:2018 FL2 SPD + CIE 1931 CMFs"},
  {"CIE FL2, 1964 10 degree", 0.37927483,0.36722793, "CIE 015:2018 FL2 SPD + CIE 1964 CMFs"},
  {"CIE FL3, 1931 2 degree", 0.40909004,0.39411713, "CIE 015:2018 FL3 SPD + CIE 1931 CMFs"},
  {"CIE FL3, 1964 10 degree", 0.41764468,0.38312450, "CIE 015:2018 FL3 SPD + CIE 1964 CMFs"},
  {"CIE FL4, 1931 2 degree", 0.44018110,0.40309069, "CIE 015:2018 FL4 SPD + CIE 1931 CMFs"},
  {"CIE FL4, 1964 10 degree", 0.44924770,0.39060548, "CIE 015:2018 FL4 SPD + CIE 1964 CMFs"},
  {"CIE FL5, 1931 2 degree", 0.31375735,0.34516065, "CIE 015:2018 FL5 SPD + CIE 1931 CMFs"},
  {"CIE FL5, 1964 10 degree", 0.31974054,0.34236696, "CIE 015:2018 FL5 SPD + CIE 1964 CMFs"},
  {"CIE FL6, 1931 2 degree", 0.37787777,0.38819415, "CIE 015:2018 FL6 SPD + CIE 1931 CMFs"},
  {"CIE FL6, 1964 10 degree", 0.38662283,0.37837312, "CIE 015:2018 FL6 SPD + CIE 1964 CMFs"},
  {"CIE FL7, 1931 2 degree", 0.31285247,0.32917418, "CIE 015:2018 FL7 SPD + CIE 1931 CMFs"},
  {"CIE FL7, 1964 10 degree", 0.31564564,0.32950815, "CIE 015:2018 FL7 SPD + CIE 1964 CMFs"},
  {"CIE FL8, 1931 2 degree", 0.34580575,0.35861758, "CIE 015:2018 FL8 SPD + CIE 1931 CMFs"},
  {"CIE FL8, 1964 10 degree", 0.34896556,0.35931730, "CIE 015:2018 FL8 SPD + CIE 1964 CMFs"},
  {"CIE FL9, 1931 2 degree", 0.37409927,0.37268420, "CIE 015:2018 FL9 SPD + CIE 1931 CMFs"},
  {"CIE FL9, 1964 10 degree", 0.37825426,0.37038210, "CIE 015:2018 FL9 SPD + CIE 1964 CMFs"},
  {"CIE FL10, 1931 2 degree", 0.34578790,0.35875793, "CIE 015:2018 FL10 SPD + CIE 1931 CMFs"},
  {"CIE FL10, 1964 10 degree", 0.35061017,0.35430334, "CIE 015:2018 FL10 SPD + CIE 1964 CMFs"},
  {"CIE FL11, 1931 2 degree", 0.38053749,0.37691531, "CIE 015:2018 FL11 SPD + CIE 1931 CMFs"},
  {"CIE FL11, 1964 10 degree", 0.38543539,0.37109479, "CIE 015:2018 FL11 SPD + CIE 1964 CMFs"},
  {"CIE FL12, 1931 2 degree", 0.43702434,0.40421500, "CIE 015:2018 FL12 SPD + CIE 1931 CMFs"},
  {"CIE FL12, 1964 10 degree", 0.44265468,0.39706106, "CIE 015:2018 FL12 SPD + CIE 1964 CMFs"}
};

#define N_WHITE_POINTS ((int) (sizeof (white_points) / sizeof (white_points[0])))

// Function prototypes
int inputtext (char *);
int readline (FILE *, char *, int);
int parse_int_string (const char *, int *);
int parse_cmf_record (const char *, double *);
int choose_cmf (int *, char *, double *);
int load_cmf (int, char *, double **);
int cmf (double, int, double **, double *);
int xy2uv (double, double, int *, int *, int, int);
int plot (int, int, int *, uint8_t *, int, int);
int draw_num (int, int, char *, int *, uint8_t *, int, int);
int draw_line (int *, int *, int *, uint8_t *, int, int);
int within_polygon (double, double, double **, int);
int bmp (char *, uint8_t *, int, int);
void write_u16_le (FILE *, uint16_t);
void write_u32_le (FILE *, uint32_t);
void write_s32_le (FILE *, int32_t);
int rgb_primaries (double **);
int illum_white (double *);
int *allocate_intmem (size_t);
int **allocate_intmemp (size_t);
char *allocate_strmem (size_t);
double *allocate_doublemem (size_t);
double **allocate_doublememp (size_t);
uint8_t *allocate_ustrmem (size_t);
size_t image_buffer_size (int, int);

// Set some symbolic constants.
#define MAX_STRINGLEN 256  // Maximum number of characters per line

int
main (void) {

  int i, j, nlines, narray, *rgb, width, height, count, u, v;
  int **uv, add_axes, mark_white, uborder, vborder, mark_rgb;
  int *point1_uv, *point2_uv, fill_srgb, axes_width, axes_height;
  double INTERVAL, **cmxyz, lambda, **xyzbar, **xyz, range, sample_count, sum;
  double *white_xyz, val, **p, **polygon, *xyzvector, *rgb_double;
  uint8_t *buffer, *buffer2;
  char *filename, *temp;

  // Wavelength step size (interpolate if necessary)
  // INTERVAL is changed below to 10 nm for Judd 1951 data; even still, it's poorly behaved.
  INTERVAL = 0.1;

  // Side and bottom borders to add if axes are to be drawn.
  uborder = 44;  // Side border widths
  vborder = 40;  // Top & bottom border widths

  // Allocate memory for various arrays.
  filename = allocate_strmem (MAX_STRINGLEN);
  temp = allocate_strmem (MAX_STRINGLEN);
  rgb = allocate_intmem (3);
  white_xyz = allocate_doublemem (3);
  xyzvector = allocate_doublemem (3);
  rgb_double = allocate_doublemem (3);
  p = allocate_doublememp (3);
  polygon = allocate_doublememp (3);
  for (i=0; i<3; i++) {
    p[i] = allocate_doublemem (3);
    polygon[i] = allocate_doublemem (2);
  }
  point1_uv = allocate_intmem (2);
  point2_uv = allocate_intmem (2);
  buffer2 = NULL;

  // Ask for bitmap dimension.
  fprintf (stdout, "\nBitmap dimension (square) (px)? ");
  memset (temp, 0, MAX_STRINGLEN * sizeof (char));
  inputtext (temp);
  if ((parse_int_string (temp, &width) == EXIT_FAILURE) ||
      (width < 2) || (width > INT_MAX - (2 * uborder))) {
    fprintf (stderr, "ERROR: Invalid bitmap dimension: %s\n", temp);
    exit (EXIT_FAILURE);
  }
  height = width;  // Same as width so locus isn't distorted.

  // Allocate bitmap data without axes. Multiplication is checked before calloc().
  buffer = allocate_ustrmem (image_buffer_size (width, height));

  // Choose color-matching function (CMF).
  for (;;) {
    if (choose_cmf (&nlines, filename, &INTERVAL) > -1) break;
  }

  // Allocate memory for various arrays.
  cmxyz = allocate_doublememp ((size_t) nlines);
  for (i=0; i<(nlines); i++) {
    cmxyz[i] = allocate_doublemem (4);  // lambda, xbar, ybar, zbar
  }

  // Load color-matching function (CMF).
  for (;;) {
    if (load_cmf (nlines, filename, cmxyz) > -1) break;
  }

  // Ask whether to plot outline of rgb gamut on bitmap.
  fprintf (stdout, "\nPlot outline of RGB gamut on bitmap (y/n)? ");
  memset (temp, 0, MAX_STRINGLEN * sizeof (char));
  inputtext (temp);
  if ((temp[0] == 'y') || (temp[0] == 'Y')) {
    mark_rgb = 1;
    for (;;) {
      if (rgb_primaries (p) > -1) break;  // Choose RGB primaries.
    }
  } else {
    mark_rgb = 0;
  }

  // Ask whether to fill in the sRGB color gamut on bitmap.
  fprintf (stdout, "\nFill in the sRGB color gamut on bitmap (y/n)? ");
  memset (temp, 0, MAX_STRINGLEN * sizeof (char));
  inputtext (temp);
  if ((temp[0] == 'y') || (temp[0] == 'Y')) {
    fill_srgb = 1;
  } else {
    fill_srgb = 0;
  }

  // Ask whether to mark white point on bitmap.
  fprintf (stdout, "\nMark white point on bitmap (y/n)? ");
  memset (temp, 0, MAX_STRINGLEN * sizeof (char));
  inputtext (temp);
  if ((temp[0] == 'y') || (temp[0] == 'Y')) {
    mark_white = 1;
    for (;;) {
      if (illum_white (white_xyz) > -1) break;  // Choose white point.
    }
  } else {
    mark_white = 0;
  }

  // Allocate enough sample slots for the requested interpolation interval.
  range = cmxyz[nlines - 1][0] - cmxyz[0][0];
  sample_count = ceil (range / INTERVAL) + 1.0;
  if (!isfinite (sample_count) || (sample_count < 2.0) || (sample_count > (double) INT_MAX)) {
    fprintf (stderr, "ERROR: Invalid number of interpolated CMF samples.\n");
    exit (EXIT_FAILURE);
  }
  narray = (int) sample_count;
  xyzbar = allocate_doublememp ((size_t) narray);
  xyz = allocate_doublememp ((size_t) narray);
  uv = allocate_intmemp ((size_t) narray);
  for (i=0; i<narray; i++) {
    xyz[i] = allocate_doublemem (3);
    xyzbar[i] = allocate_doublemem (3);
    uv[i] = allocate_intmem (2);
  }

  // Ask whether to include axes on bitmap.
  fprintf (stdout, "\nInclude axes on bitmap? (borders will be added to bitmap size) (y/n)? ");
  memset (temp, 0, MAX_STRINGLEN * sizeof (char));
  inputtext (temp);
  if ((temp[0] == 'y') || (temp[0] == 'Y')) {
    add_axes = 1;
    axes_width = width + (2 * uborder);
    axes_height = height + (2 * vborder);
    buffer2 = allocate_ustrmem (image_buffer_size (axes_width, axes_height));
  } else {
    add_axes = 0;
    axes_width = 0;
    axes_height = 0;
  }

  // Loop through full range of wavelengths in nm steps defined by INTERVAL.
  count = 0;
  lambda = cmxyz[0][0];  // Start at first wavelength in CMF table.
  while ((count < narray) &&
         (lambda <= cmxyz[nlines - 1][0] + (fabs (cmxyz[nlines - 1][0]) * 1.0e-12 + 1.0e-12))) {

    if (lambda > cmxyz[nlines - 1][0]) lambda = cmxyz[nlines - 1][0];

    // Retrieve wavelength and color-matching coodinates xbar, ybar, zbar for requested wavelength.
    // Note: For single wavelength, the CMF X,Y,Z coordinates are the linear scene tristimulus values. i.e., no need to integrate
    // xbar = xyzbar[0], ybar = xyzbar[1], zbar = xyzbar[2]
    if (cmf (lambda, nlines, cmxyz, xyzbar[count]) == EXIT_FAILURE) {
      fprintf (stderr, "ERROR: lambda %.18lf is outside range of CMF data: %.18lf to %.18lf.\n", lambda, cmxyz[0][0], cmxyz[nlines-1][0]);
      exit (EXIT_FAILURE);
    }
//  fprintf (stdout, "XYZ CMF coordinates: %.12lf %.12lf %.12lf\n", xyzbar[count][0], xyzbar[count][1], xyzbar[count][2]);

    // Compute spectral chromaticity coordinates xyz.
    sum = xyzbar[count][0] + xyzbar[count][1] + xyzbar[count][2];
    if (!isfinite (sum) || (fabs (sum) <= DBL_MIN)) {
      fprintf (stderr, "ERROR: CMF tristimulus sum is zero or non-finite at %.12g nm.\n", lambda);
      exit (EXIT_FAILURE);
    }
    xyz[count][0] = xyzbar[count][0] / sum;
    xyz[count][1] = xyzbar[count][1] / sum;
    xyz[count][2] = xyzbar[count][2] / sum;
//  fprintf (stdout, "Spectral chromaticity coordinates xyz: %.12lf %.12lf %.12lf\n", xyz[count][0], xyz[count][1], xyz[count][2]);

    // Convert xy chromaticity coordinates to uv pixel coordinates.
    if (xy2uv (xyz[count][0], xyz[count][1], &u, &v, width, height) == EXIT_FAILURE) {
      fprintf (stderr, "ERROR: Invalid chromaticity coordinates at %.12g nm.\n", lambda);
      exit (EXIT_FAILURE);
    }
    uv[count][0] = u;
    uv[count][1] = v;

    count++;
    lambda = cmxyz[0][0] + (((double) count) * INTERVAL);
  }
  if (count < 2) {
    fprintf (stderr, "ERROR: Fewer than two chromaticity samples were generated.\n");
    exit (EXIT_FAILURE);
  }
//fprintf (stdout, "count: %i   narray: %i\n", count, narray);

  // Create a bitmap of the chromaticity locus.
  // Plot locus points as white pixels.
  rgb[0] = 255;
  rgb[1] = 255;
  rgb[2] = 255;
  for (i=0; i<count; i++) {
    plot (uv[i][0], uv[i][1], rgb, buffer, width, height);
  }

  // The line of purples directly joins the two wavelength endpoints of the spectral locus.
  draw_line (uv[0], uv[count - 1], rgb, buffer, width, height);

  // Fill in sRGB color gamut if requested.
  if (fill_srgb) {

    // Define a polygon by sRGB (BT.709) color primaries as vertices.
    polygon[0][0] = rgb_spaces[0].xr * (double) (width - 1); polygon[0][1] = rgb_spaces[0].yr * (double) (height - 1);  // Red
    polygon[1][0] = rgb_spaces[0].xg * (double) (width - 1); polygon[1][1] = rgb_spaces[0].yg * (double) (height - 1);  // Green
    polygon[2][0] = rgb_spaces[0].xb * (double) (width - 1); polygon[2][1] = rgb_spaces[0].yb * (double) (height - 1);  // Blue
   
    // Loop through all pixels in bitmap.
    for (u=0; u<width; u++) {
      for (v=0; v<height; v++) {

        // If pixel is within the sRGB color gamut, we plot it.
        if (within_polygon ((double) u, (double) v, polygon, 3)) {

          // Calculate chromaticity coordinates for current pixel.
          xyzvector[0] = (double) u / (double) (width - 1);
          xyzvector[1] = (double) v / (double) (height - 1);
          xyzvector[2] = 1.0 - xyzvector[0] - xyzvector[1];

          // Convert linear xyz chromaticity coordinates to linear rgb chromaticity coordinates.
          rgb_double[0] = ((3.2408357 * xyzvector[0]) - (1.5373195 * xyzvector[1]) - (0.4985901 * xyzvector[2])) /
                          ((2.3272512 * xyzvector[0]) + (0.1345891 * xyzvector[1]) + (0.6002182 * xyzvector[2]));
          rgb_double[1] = ((-0.9692294 * xyzvector[0]) + (1.8759400 * xyzvector[1]) + (0.0415544 * xyzvector[2])) /
                          ((2.3272512 * xyzvector[0]) + (0.1345891 * xyzvector[1]) + (0.6002182 * xyzvector[2]));
          rgb_double[2] = ((0.0556449 * xyzvector[0]) - (0.2040314 * xyzvector[1]) + (1.0572538 * xyzvector[2])) /
                          ((2.3272512 * xyzvector[0]) + (0.1345891 * xyzvector[1]) + (0.6002182 * xyzvector[2]));

          // Apply sRGB gamma-correction to linear rgb.
          for (i=0; i<3; i++) {
            if (rgb_double[i] > 0.0031308) {
              rgb_double[i] = (1.055 * pow (rgb_double[i], (1.0 / 2.4))) - 0.055;
            } else {
              rgb_double[i] *= 12.92;
            }
          }

          // Scale to achieve required range 0-255.
          for (i=0; i<3; i++) {
            rgb_double[i] *= 255.0;
          }

          // Convert to integer.
          for (i=0; i<3; i++) {
            rgb[i] = (int) (rgb_double[i] + 0.5);
          }

          // Clip any undershoot or overshoot resulting from the fact that
          // 8-bit RGB (0-255) has much smaller color gamut than xyz.
          for (i=0; i<3; i++) {
            if (rgb[i] < 0) rgb[i] = 0;
            if (rgb[i] > 255) rgb[i] = 255;
          }

//  fprintf (stdout, "\nsRGB coordinates (with sRGB gamma-correction): (%i, %i, %i)\n", rgb[0], rgb[1], rgb[2]);
//  fprintf (stdout, "\n");
          plot (u, v, rgb, buffer, width, height);

        }  // End if within locus
      }  // Next v
    }  // Next u
  }  // End if fill_srgb

  // Plot a marker at white point, if requested.
  if (mark_white) {
    if (fill_srgb) {  // Use black if coloring-in sRGB gamut.
      rgb[0] = 0;
      rgb[1] = 0;
      rgb[2] = 0;
    } else {  // Use white if not coloring-in sRGB gamut.
      rgb[0] = 255;
      rgb[1] = 255;
      rgb[2] = 255;
    }

    xy2uv (white_xyz[0], white_xyz[1], &u, &v, width, height);

    // Horizontal stroke
    for (i=-5; i<6; i++) {
      plot (u + i, v, rgb, buffer, width, height);
    }

    // Vertical stroke
    for (j=-5; j<6; j++) {
      plot (u, v + j, rgb, buffer, width, height);
    }
  }

  // Plot outline of selected RGB gamut if requested.
  if (mark_rgb) {
    rgb[0] = 255;
    rgb[1] = 255;
    rgb[2] = 255;
    xy2uv (p[0][0], p[1][0], &point1_uv[0], &point1_uv[1], width, height);  // Red primary uv coordinates
    xy2uv (p[0][1], p[1][1], &point2_uv[0], &point2_uv[1], width, height);  // Green primary uv coordinates
    draw_line (point1_uv, point2_uv, rgb, buffer, width, height);  // Line from red primary to green primary

    xy2uv (p[0][2], p[1][2], &point2_uv[0], &point2_uv[1], width, height);  // Blue primary uv coordinates
    draw_line (point1_uv, point2_uv, rgb, buffer, width, height);  // Line from red primary to blue primary
    xy2uv (p[0][1], p[1][1], &point1_uv[0], &point1_uv[1], width, height);  // Green primary uv coordinates
    draw_line (point1_uv, point2_uv, rgb, buffer, width, height);  // Line from green primary to blue primary
  }

  // Include axes on bitmap, if requested.
  if (add_axes) {

    // Use white pixels.
    rgb[0] = 255;
    rgb[1] = 255;
    rgb[2] = 255;

    // Copy locus to buffer2 which will contain axes.
    for (v=0; v<height; v++) {
      for (u=0; u<width; u++) {
        size_t src = ((size_t) u + (size_t) v * (size_t) width) * 3u;
        size_t dst = ((size_t) (u + uborder) +
                      (size_t) (v + vborder) * (size_t) axes_width) * 3u;
        buffer2[dst] = buffer[src];
        buffer2[dst + 1u] = buffer[src + 1u];
        buffer2[dst + 2u] = buffer[src + 2u];
      }
    }

    // Plot horizontal axis line from x = 0.0 through x = 1.0.
    point1_uv[0] = uborder;
    point1_uv[1] = vborder;
    point2_uv[0] = uborder + width - 1;
    point2_uv[1] = vborder;
    draw_line (point1_uv, point2_uv, rgb, buffer2, axes_width, axes_height);

    // Plot horizontal tick marks and numerical labels.
    for (i=0; i<=10; i++) {
      u = uborder + (int) lround ((double) i * (double) (width - 1) / 10.0);
      point1_uv[0] = u; point1_uv[1] = vborder - 10;
      point2_uv[0] = u; point2_uv[1] = vborder;
      draw_line (point1_uv, point2_uv, rgb, buffer2, axes_width, axes_height);

      val = (double) i / 10.0;
      snprintf (temp, MAX_STRINGLEN, "%0.1f", val);
      draw_num (u - 14, vborder - 24, temp, rgb, buffer2, axes_width, axes_height);
    }

    // Plot vertical axis line from y = 0.0 through y = 1.0.
    point1_uv[0] = uborder;
    point1_uv[1] = vborder;
    point2_uv[0] = uborder;
    point2_uv[1] = vborder + height - 1;
    draw_line (point1_uv, point2_uv, rgb, buffer2, axes_width, axes_height);

    // Plot vertical tick marks and numerical labels.
    for (i=0; i<=10; i++) {
      v = vborder + (int) lround ((double) i * (double) (height - 1) / 10.0);
      point1_uv[0] = uborder - 10; point1_uv[1] = v;
      point2_uv[0] = uborder; point2_uv[1] = v;
      draw_line (point1_uv, point2_uv, rgb, buffer2, axes_width, axes_height);

      val = (double) i / 10.0;
      snprintf (temp, MAX_STRINGLEN, "%0.1f", val);
      draw_num (0, v - 4, temp, rgb, buffer2, axes_width, axes_height);
    }

    if (bmp ("out.bmp", buffer2, axes_width, axes_height) == EXIT_FAILURE) {
      exit (EXIT_FAILURE);
    }

  } else {
    if (bmp ("out.bmp", buffer, width, height) == EXIT_FAILURE) {
      exit (EXIT_FAILURE);
    }
  }  // End if add_axes

  fprintf (stdout, "\n");

  // Free allocated memory.
  free (filename);
  free (temp);
  free (white_xyz);
  free (xyzvector);
  free (rgb_double);
  for (i=0; i<3; i++) {
    free (p[i]); 
    free (polygon[i]);
  } 
  free (p);
  free (polygon);
  free (point1_uv);
  free (point2_uv);
  free (rgb);
  for (i=0; i<nlines; i++) {
    free (cmxyz[i]);
  }
  free (cmxyz);
  for (i=0; i<narray; i++) {
    free (xyzbar[i]);
    free (xyz[i]);
    free (uv[i]);
  }
  free (xyzbar);
  free (xyz);
  free (uv);
  free (buffer);
  free (buffer2);

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

// Parse one numeric CMF record: wavelength, xbar, ybar, zbar.
// Blank lines, comments, and textual headers are ignored (return 0).
// A line that starts numerically but is malformed returns -1.
int
parse_cmf_record (const char *line, double *values) {

  const char *p;
  char extra;
  int n;

  if ((line == NULL) || (values == NULL)) return (-1);
  p = line;
  while (isspace ((unsigned char) *p)) p++;
  if ((*p == '\0') || (*p == '#') || (*p == ';')) return (0);

  if (!isdigit ((unsigned char) *p) && (*p != '+') && (*p != '-') && (*p != '.')) {
    return (0);
  }

  n = sscanf (p, "%lf %lf %lf %lf %c",
              &values[0], &values[1], &values[2], &values[3], &extra);
  if (n != 4) return (-1);
  if (!isfinite (values[0]) || !isfinite (values[1]) ||
      !isfinite (values[2]) || !isfinite (values[3])) return (-1);

  return (1);
}

// Read a single line of text from a csv text file.
// Convert commas to spaces.
// Returns -1 if EOF is encountered.
int
readline (FILE *fi, char *line, int limit) {

  int ch, i;

  if ((fi == NULL) || (line == NULL) || (limit < 2)) return (-2);

  i = 0;
  while (i < limit - 1) {
    ch = fgetc (fi);

    if (ch == EOF) {
      if (i == 0) return (-1);
      line[i] = '\0';
      return (0);
    }
    if (ch == '\r') continue;
    if (ch == ',') ch = ' ';
    if (ch == '\n') {
      line[i] = '\0';
      return (0);
    }

    line[i++] = (char) ch;
  }

  // A physical line that does not fit is rejected rather than silently truncated.
  line[i] = '\0';
  while ((ch = fgetc (fi)) != '\n' && ch != EOF) {
  }
  return (-2);
}

// Choose color-matching function (CMF).
// Returns: -1 if invalid selection, 0 if valid selection
int
choose_cmf (int *nlines, char *filename, double *INTERVAL) {

  int choice, parsed, status;
  double values[4];
  char *temp;
  FILE *fi;

  temp = allocate_strmem (MAX_STRINGLEN);

  fprintf (stdout, "\nChoose color-matching function (CMF):\n\n");
  fprintf (stdout, "  1 - 1964 10-deg XYZ CMFs (JIS Z 8701:1999)\n");
  fprintf (stdout, "  2 - 1931 2-deg XYZ CIE CMFs (CIE.15.2004)\n");
  fprintf (stdout, "  3 - 1931 2-deg XYZ CIE CMFs with Judd (1951) modifications [badly behaved data; won't interpolate here]\n");
  fprintf (stdout, "  4 - 1931 2-deg XYZ CIE CMFs with Judd (1951) and Vos (1978) modifications\n");
  fprintf (stdout, "  5 - 2006 2-deg XYZ CMFs transformed from the CIE (2006) 2-deg LMS cone fundamentals\n");
  fprintf (stdout, "  6 - 2006 10-deg XYZ CMFs transformed from the CIE (2006) 10-deg LMS cone fundamentals\n");
  fprintf (stdout, "  7 - Enter filename for CMFs (nm, xbar, ybar, zbar as .csv)\n");
  fprintf (stdout, "\nChoice? ");
  inputtext (temp);
  if (parse_int_string (temp, &choice) == EXIT_FAILURE) {
    fprintf (stderr, "ERROR: Cannot make integer of: %s\n", temp);
    free (temp);
    return (-1);
  }

  switch (choice) {
    case 1:
      snprintf (filename, MAX_STRINGLEN, "data/CIE_xyz_1964_10deg.csv");
      break;
    case 2:
      snprintf (filename, MAX_STRINGLEN, "data/CIE_xyz_1931_2deg.csv");
      break;
    case 3:
      snprintf (filename, MAX_STRINGLEN, "data/CIE_xyz_1931_2deg_judd1951.csv");
      *INTERVAL = 10.0;
      break;
    case 4:
      snprintf (filename, MAX_STRINGLEN, "data/CIE_xyz_1931_2deg_judd1951_vos1978.csv");
      break;
    case 5:
      snprintf (filename, MAX_STRINGLEN, "data/CIE_xyz_2006_2deg_lms_cones.csv");
      break;
    case 6:
      snprintf (filename, MAX_STRINGLEN, "data/CIE_xyz_2006_10deg_lms_cones.csv");
      break;
    case 7:
      fprintf (stdout, "Filename for csv CMFs? ");
      inputtext (filename);
      break;
    default:
      fprintf (stderr, "Invalid selection.\n");
      free (temp);
      return (-1);
  }

  fi = fopen (filename, "r");
  if (fi == NULL) {
    fprintf (stderr, "ERROR: Unable to open input file \"%s\".\n", filename);
    free (temp);
    exit (EXIT_FAILURE);
  }

  *nlines = 0;
  while ((status = readline (fi, temp, MAX_STRINGLEN)) != -1) {
    if (status == -2) {
      fprintf (stderr, "ERROR: Line in color-matching file exceeds %d characters.\n", MAX_STRINGLEN - 1);
      fclose (fi);
      free (temp);
      exit (EXIT_FAILURE);
    }
    parsed = parse_cmf_record (temp, values);
    if (parsed < 0) {
      fprintf (stderr, "ERROR: Malformed numeric CMF row: %s\n", temp);
      fclose (fi);
      free (temp);
      exit (EXIT_FAILURE);
    }
    if (parsed > 0) (*nlines)++;
  }

  if (*nlines < 2) {
    fprintf (stderr, "ERROR: Color-matching file must contain at least two numeric rows.\n");
    fclose (fi);
    free (temp);
    exit (EXIT_FAILURE);
  }

  fprintf (stdout, "%i lines in color-matching file: %s\n", *nlines, filename);
  fclose (fi);
  free (temp);
  return (0);
}

// Load color-matching function (CMF).
// Returns: -1 if invalid selection, 0 if valid selection
int
load_cmf (int nlines, char *filename, double **cmxyz) {

  int i, parsed, status;
  double values[4];
  char *temp;
  FILE *fi;

  if ((nlines < 2) || (filename == NULL) || (cmxyz == NULL)) return (EXIT_FAILURE);
  temp = allocate_strmem (MAX_STRINGLEN);

  fi = fopen (filename, "r");
  if (fi == NULL) {
    fprintf (stderr, "ERROR: Unable to open input file \"%s\".\n", filename);
    free (temp);
    exit (EXIT_FAILURE);
  }

  i = 0;
  while ((status = readline (fi, temp, MAX_STRINGLEN)) != -1) {
    if (status == -2) {
      fprintf (stderr, "ERROR: Line in color-matching file exceeds %d characters.\n", MAX_STRINGLEN - 1);
      fclose (fi);
      free (temp);
      exit (EXIT_FAILURE);
    }

    parsed = parse_cmf_record (temp, values);
    if (parsed < 0) {
      fprintf (stderr, "ERROR: Malformed numeric CMF row: %s\n", temp);
      fclose (fi);
      free (temp);
      exit (EXIT_FAILURE);
    }
    if (parsed == 0) continue;
    if (i >= nlines) {
      fprintf (stderr, "ERROR: CMF row count changed between count and load passes.\n");
      fclose (fi);
      free (temp);
      exit (EXIT_FAILURE);
    }

    cmxyz[i][0] = values[0];
    cmxyz[i][1] = values[1];
    cmxyz[i][2] = values[2];
    cmxyz[i][3] = values[3];

    if ((i > 0) && !(cmxyz[i][0] > cmxyz[i - 1][0])) {
      fprintf (stderr, "ERROR: CMF wavelengths must be strictly increasing; %.12g follows %.12g.\n",
               cmxyz[i][0], cmxyz[i - 1][0]);
      fclose (fi);
      free (temp);
      exit (EXIT_FAILURE);
    }
    i++;
  }

  fclose (fi);
  free (temp);

  if (i != nlines) {
    fprintf (stderr, "ERROR: Expected %d CMF rows but loaded %d.\n", nlines, i);
    exit (EXIT_FAILURE);
  }
  return (0);
}

// Return xyz color-matching function (CMF) coordinates by interpolating the data.
int
cmf (double lambda, int nlines, double **cmxyz, double *xyzbar) {

  int i, lb, ub;
  double frac;

  if ((nlines < 2) || (cmxyz == NULL) || (xyzbar == NULL) || !isfinite (lambda)) {
    return (EXIT_FAILURE);
  }

  // If lambda is above or below table data, return error.
  if ((lambda < cmxyz[0][0]) || (lambda > cmxyz[nlines-1][0])) {
    return (EXIT_FAILURE);

  // Otherwise, interpolate if necessary.
  } else {
 
    i = 0;  // Index of wavelengths in cmxyz. i.e., cmxyz[i][0] = wavelength i
    while (i < nlines) {
      if (cmxyz[i][0] < lambda) {
        i++;
      } else {
        break;
      }
    }

    // Lower bound
    if (i > 0) {
      lb = i - 1;
    } else {
      lb = i;
    }

    // Upper bound
    ub = lb + 1;

    if (!(cmxyz[ub][0] > cmxyz[lb][0])) return (EXIT_FAILURE);
    frac = (lambda - cmxyz[lb][0]) / (cmxyz[ub][0] - cmxyz[lb][0]);

    // Interpolate CMF coordinates.
    xyzbar[0] = (frac * (cmxyz[ub][1] - cmxyz[lb][1])) + cmxyz[lb][1];
    xyzbar[1] = (frac * (cmxyz[ub][2] - cmxyz[lb][2])) + cmxyz[lb][2];
    xyzbar[2] = (frac * (cmxyz[ub][3] - cmxyz[lb][3])) + cmxyz[lb][3];
  }

  return (EXIT_SUCCESS);
}

// Convert xy chromaticity coordinates to uv (bitmap px).
int
xy2uv (double x, double y, int *u, int *v, int width, int height) {

  double upos, vpos;

  if ((u == NULL) || (v == NULL) || (width < 1) || (height < 1) ||
      !isfinite (x) || !isfinite (y)) return (EXIT_FAILURE);

  // Pixel indices run from 0 through width - 1 and height - 1. Chromaticities
  // outside [0,1] are intentionally allowed here for imaginary RGB primaries;
  // plot() clips those portions to the bitmap.
  upos = x * (double) (width - 1);
  vpos = y * (double) (height - 1);
  if (!isfinite (upos) || !isfinite (vpos) ||
      (upos < (double) INT_MIN) || (upos > (double) INT_MAX) ||
      (vpos < (double) INT_MIN) || (vpos > (double) INT_MAX)) {
    return (EXIT_FAILURE);
  }

  *u = (int) lround (upos);
  *v = (int) lround (vpos);
  return (EXIT_SUCCESS);
}

// Plot one pixel, clipping coordinates that fall outside the bitmap.
int
plot (int u, int v, int *rgb, uint8_t *buffer, int width, int height) {

  size_t index;

  if ((rgb == NULL) || (buffer == NULL) || (width <= 0) || (height <= 0)) {
    return (EXIT_FAILURE);
  }
  if ((u < 0) || (u >= width) || (v < 0) || (v >= height)) {
    return (EXIT_FAILURE);
  }

  index = ((size_t) u + (size_t) v * (size_t) width) * 3u;
  buffer[index] = (uint8_t) rgb[2];  // B
  buffer[index + 1u] = (uint8_t) rgb[1];  // G
  buffer[index + 2u] = (uint8_t) rgb[0];  // R
  return (EXIT_SUCCESS);
}

// Draw decimal numbers to bitmap buffer.
// Input is terminated string text. e.g., "0.5" + '\0'
int
draw_num (int x0, int y0, char *text, int *rgb, uint8_t *buffer, int width, int height) {

  int x, y, digit, row, col;
  uint8_t line;

  // Font definitions of digits 0 - 9 (8 x 8 px)
  const unsigned char font[10][8] = {
    {0x7e, 0x81, 0x81, 0x81, 0x81, 0x81, 0x81, 0x7e},  // 0
    {0x08, 0x18, 0x28, 0x08, 0x08, 0x08, 0x08, 0x7e},  // 1
    {0xff, 0x01, 0x01, 0xff, 0x80, 0x80, 0x80, 0xff},  // 2
    {0xff, 0x01, 0x01, 0x7f, 0x01, 0x01, 0x01, 0xff},  // 3
    {0x84, 0x84, 0x84, 0xff, 0x04, 0x04, 0x04, 0x04},  // 4
    {0xff, 0x80, 0x80, 0xff, 0x01, 0x01, 0x01, 0xff},  // 5
    {0x7e, 0x80, 0x80, 0xff, 0x81, 0x81, 0x81, 0x7e},  // 6
    {0xff, 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40},  // 7
    {0x7e, 0x81, 0x81, 0x7e, 0x81, 0x81, 0x81, 0x7e},  // 8
    {0xff, 0x81, 0x81, 0xff, 0x01, 0x01, 0x01, 0x01}   // 9
  };

  // Starting coordinates.
  x = x0;
  y = y0;

  // Loop until string termination.
  while (*text != '\0') {

    // Draw digit.
    if (*text >= '0' && *text <= '9') {

      digit = *text - '0';  // Convert from ASCII to index 0 to 9.
//printf ("%i ", digit);
      for (row=0; row<8; row++) {
        line = font[digit][7 - row];
        for (col=0; col<8; col++) {
          if (line & (0x80 >> col)) {  // 0x80 - 10000000 in binary
            plot (x + col, y + row, rgb, buffer, width, height);
          }
        }
      }

    // Draw decimal point.
    } else if (*text == '.') {

      for (row=0; row<2; row++) {
        line = 0x18;
        for (col=0; col<8; col++) {
          if (line & (0x80 >> col)) {
            plot (x + col, y + row, rgb, buffer, width, height);
          }
        }
      }
    }
    x += 8 + 2;
    text++;
  }

  return (EXIT_SUCCESS);
}

// Draw a line between points p1 and p2 on bitmap.
// p1 and p2 are (u,v) vectors in units of pixels.
int
draw_line (int *p1, int *p2, int *rgb, uint8_t *buffer, int width, int height) {

  int64_t x0, y0, x1, y1, dx, dy, sx, sy, err, e2;

  if ((p1 == NULL) || (p2 == NULL) || (rgb == NULL) || (buffer == NULL)) {
    return (EXIT_FAILURE);
  }

  x0 = p1[0]; y0 = p1[1];
  x1 = p2[0]; y1 = p2[1];
  dx = llabs (x1 - x0);
  sx = (x0 < x1) ? 1 : -1;
  dy = -llabs (y1 - y0);
  sy = (y0 < y1) ? 1 : -1;
  err = dx + dy;

  for (;;) {
    plot ((int) x0, (int) y0, rgb, buffer, width, height);
    if ((x0 == x1) && (y0 == y1)) break;
    e2 = 2 * err;
    if (e2 >= dy) {
      err += dy;
      x0 += sx;
    }
    if (e2 <= dx) {
      err += dx;
      y0 += sy;
    }
  }

  return (EXIT_SUCCESS);
}

// Point-in-polygon test.
// Takes list uv of n points coordinates as polygon vertices.
// Returns 1 if test point (u,v) is inside polygon, 0 if outside.
int
within_polygon (double u, double v, double **uv, int n) {

  int i, j, inside;

  inside = 0;  // Default to coordinates are outside of polygon.
  for (i=0, j=n-1; i<n; j=i++) {
    if ((((uv[i][1] <= v) && (v < uv[j][1])) ||
       ((uv[j][1] <= v) && (v < uv[i][1]))) &&
       (u < (uv[j][0] - uv[i][0]) * (v - uv[i][1]) / (uv[j][1] - uv[i][1]) + uv[i][0]))
      inside = !inside;
  }

  return (inside);
}

// Create a bitmap output file.
int
bmp (char *filename, uint8_t *buffer, int width, int height) {

  int x, y;
  size_t row_size, image_size, file_size, c, padding_len;
  uint8_t padding[3] = {0, 0, 0};
  FILE *fo;

  if ((filename == NULL) || (buffer == NULL) || (width <= 0) || (height <= 0)) {
    return (EXIT_FAILURE);
  }

  row_size = ((size_t) width * 3u + 3u) & ~(size_t) 3u;
  if ((size_t) height > SIZE_MAX / row_size) {
    fprintf (stderr, "ERROR: BMP image size overflows size_t.\n");
    return (EXIT_FAILURE);
  }
  image_size = row_size * (size_t) height;
  if (image_size > UINT32_MAX - 54u) {
    fprintf (stderr, "ERROR: BMP is too large for the 32-bit BMP file-size fields.\n");
    return (EXIT_FAILURE);
  }
  file_size = 54u + image_size;
  padding_len = row_size - (size_t) width * 3u;

  fo = fopen (filename, "rb");
  if (fo != NULL) {
    fclose (fo);
    fprintf (stderr, "Output file %s already exists.\n", filename);
    return (EXIT_FAILURE);
  }
  fo = fopen (filename, "wb");
  if (fo == NULL) {
    fprintf (stderr, "Can't open output file %s.\n", filename);
    return (EXIT_FAILURE);
  }

  // BMP file header.
  write_u16_le (fo, 0x4d42);
  write_u32_le (fo, (uint32_t) file_size);
  write_u16_le (fo, 0);
  write_u16_le (fo, 0);
  write_u32_le (fo, 54);

  // BMP information header.
  write_u32_le (fo, 40);
  write_s32_le (fo, (int32_t) width);
  write_s32_le (fo, (int32_t) height);
  write_u16_le (fo, 1);
  write_u16_le (fo, 24);
  write_u32_le (fo, 0);
  write_u32_le (fo, (uint32_t) image_size);
  write_s32_le (fo, 7874);
  write_s32_le (fo, 7874);
  write_u32_le (fo, 0);
  write_u32_le (fo, 0);

  c = 0;
  for (y=0; y<height; y++) {
    for (x=0; x<width; x++) {
      fputc ((int) buffer[c++], fo);
      fputc ((int) buffer[c++], fo);
      fputc ((int) buffer[c++], fo);
    }
    if ((padding_len > 0u) && (fwrite (padding, 1, padding_len, fo) != padding_len)) {
      fprintf (stderr, "ERROR: Unable to write BMP row padding.\n");
      fclose (fo);
      return (EXIT_FAILURE);
    }
  }

  if (ferror (fo)) {
    fprintf (stderr, "ERROR: Unable to write complete BMP file.\n");
    fclose (fo);
    return (EXIT_FAILURE);
  }
  if (fclose (fo) != 0) {
    fprintf (stderr, "ERROR: Unable to close BMP output file.\n");
    return (EXIT_FAILURE);
  }
  return (EXIT_SUCCESS);
}

// Write a unsigned little-endian 16-bit value to file.
void
write_u16_le (FILE *fo, uint16_t val) {

    fputc ((int) (val & 0xffu), fo);
    fputc ((int) ((val >> 8) & 0xffu), fo);

}

// Write a unsigned little-endian 32-bit value to file.
void
write_u32_le (FILE *fo, uint32_t val) {

    fputc ((int) (val & 0xffu), fo);
    fputc ((int) ((val >> 8) & 0xffu), fo);
    fputc ((int) ((val >> 16) & 0xffu), fo);
    fputc ((int) ((val >> 24) & 0xffu), fo);

}

// Write a signed little-endian 32-bit value to file.
void
write_s32_le (FILE *fo, int32_t val) {

    write_u32_le (fo, (uint32_t) val);  // Cast as unsigned in order to preserve bit pattern.

}

// Color primary coordinates for RGB spaces and historical RGB primary sets.
// Coordinates and names are kept in sync with the audited matrix.c data table.
int
rgb_primaries (double **p) {

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

  space = &rgb_spaces[choice - 1];
  zr = 1.0 - space->xr - space->yr;
  zg = 1.0 - space->xg - space->yg;
  zb = 1.0 - space->xb - space->yb;

  p[0][0] = space->xr;  p[0][1] = space->xg;  p[0][2] = space->xb;
  p[1][0] = space->yr;  p[1][1] = space->yg;  p[1][2] = space->yb;
  p[2][0] = zr;         p[2][1] = zg;         p[2][2] = zb;

  free (temp);
  return (0);
}

// Illuminants - Choose white point.
// Values are kept in sync with the audited matrix.c table; ACES D60-like is
// explicitly distinguished from a general CIE daylight illuminant.
int
illum_white (double *white_xyz) {

  int choice, i;
  char *temp;
  const WhitePoint *white;

  temp = allocate_strmem (MAX_STRINGLEN);
  fprintf (stdout, "\nChoose display white point for bitmap:\n");
  for (i=0; i<N_WHITE_POINTS; i++) {
    fprintf (stdout, " %2d - %s\n", i + 1, white_points[i].name);
  }
  fprintf (stdout, "\nChoice? ");
  inputtext (temp);
  if (parse_int_string (temp, &choice) == EXIT_FAILURE) {
    fprintf (stderr, "ERROR: Cannot make integer of: %s\n", temp);
    free (temp);
    return (-1);
  }
  if ((choice < 1) || (choice > N_WHITE_POINTS)) {
    fprintf (stderr, "Invalid choice.\n");
    free (temp);
    return (-1);
  }

  white = &white_points[choice - 1];
  white_xyz[0] = white->x;
  white_xyz[1] = white->y;
  white_xyz[2] = 1.0 - white->x - white->y;

  free (temp);
  return (0);
}

// Allocate memory for an array of ints.
int *
allocate_intmem (size_t len) {

  void *tmp;

  if (len == 0) {
    fprintf (stderr, "ERROR: Cannot allocate memory because len = %zu in allocate_intmem().\n", len);
    exit (EXIT_FAILURE);
  }

  tmp = calloc (len, sizeof (int));
  if (tmp != NULL) {
    return (tmp);
  } else {
    fprintf (stderr, "ERROR: Cannot allocate memory for array in allocate_intmem().\n");
    exit (EXIT_FAILURE);
  } 
}

// Allocate memory for an array of pointers to arrays of ints.
int **
allocate_intmemp (size_t len) {

  void *tmp;

  if (len == 0) {
    fprintf (stderr, "ERROR: Cannot allocate memory because len = %zu in allocate_intmemp().\n", len);
    exit (EXIT_FAILURE);
  }

  tmp = calloc (len, sizeof (int *));
  if (tmp != NULL) {
    return (tmp);
  } else {
    fprintf (stderr, "ERROR: Cannot allocate memory for array in allocate_intmemp().\n");
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

// Allocate memory for an array of unsigned chars.
uint8_t *
allocate_ustrmem (size_t len) {

  void *tmp;

  if (len == 0) {
    fprintf (stderr, "ERROR: Cannot allocate memory because len = %zu in allocate_ustrmem().\n", len);
    exit (EXIT_FAILURE);
  }

  tmp = calloc (len, sizeof (uint8_t));
  if (tmp != NULL) {
    return (tmp);
  } else {
    fprintf (stderr, "ERROR: Cannot allocate memory for array in allocate_ustrmem().\n");
    exit (EXIT_FAILURE);
  }
}


// Return byte size for a 24-bit RGB image after checking multiplication overflow.
size_t
image_buffer_size (int width, int height) {

  size_t pixels, row_size;

  if ((width <= 0) || (height <= 0)) {
    fprintf (stderr, "ERROR: Invalid image dimensions %d x %d.\n", width, height);
    exit (EXIT_FAILURE);
  }
  if ((size_t) width > SIZE_MAX / (size_t) height) {
    fprintf (stderr, "ERROR: Image dimensions overflow size_t.\n");
    exit (EXIT_FAILURE);
  }
  pixels = (size_t) width * (size_t) height;
  if (pixels > SIZE_MAX / 3u) {
    fprintf (stderr, "ERROR: RGB image byte count overflows size_t.\n");
    exit (EXIT_FAILURE);
  }

  // The program writes a classic BMP whose file-size and image-size fields are
  // 32-bit unsigned integers. Reject dimensions that can never be represented.
  row_size = ((size_t) width * 3u + 3u) & ~(size_t) 3u;
  if ((size_t) height > ((size_t) UINT32_MAX - 54u) / row_size) {
    fprintf (stderr, "ERROR: Image is too large for a classic 32-bit-size BMP file.\n");
    exit (EXIT_FAILURE);
  }

  return (pixels * 3u);
}
