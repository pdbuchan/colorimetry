/*  Copyright (C) 2024-2025 P. David Buchan (pdbuchan@gmail.com)

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

// Function prototypes
int inputtext (char *);
int rgb_primaries (double **);
int illum_white (double *);
int gaussjordan (int, double **);
int *allocate_intmem (int);
double *allocate_doublemem (int);
double **allocate_doublememp (int);
char *allocate_strmem (int);

// Set some symbolic constants.
#define MAXLEN 256  // Maximum number of characters per line

int
main (int argc, char **argv) {

  int i, j, k;
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
    if (rgb_primaries (p) > -1) break;
  }

  // Choose white coordinates.
  for (;;) {
    if (illum_white (w) > -1) break;
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

  // Show white color coordinates.
  fprintf (stdout, "\nWhite coordinates (x y z):\n");
  fprintf (stdout, "%0.4lf %0.4lf %0.4lf\n", w[0], w[1], w[2]);

  // Populate white matrix.
  // Normalize to luminosity y.
  w[0] /= w[1];
  w[2] /= w[1];
  w[1] = 1.0;

  // Compute inverse of color primaries matrix p.
  gaussjordan (3, pinv);

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

  // Compute NPM matrix, where NPM = P * C
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      for (k=0; k<3; k++) {
        npm[i][j] += p[i][k] * c[j][k];
      }
    }
  }

  // Copy npm to npminv so we can invert in-place.
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      npminv[i][j] = npm[i][j];
    }
  }
  gaussjordan (3, npminv);

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
int
gaussjordan (int n, double **matrix) {

  int i, j, k;
  double **augmented, pivot, factor;

  // Allocate memory for various arrays.
  augmented = allocate_doublememp (n);
  for (i=0; i<n; i++) {
    augmented[i] = allocate_doublemem (2 * n);  // Original matrix augmented by identity matrix
  }

  // Augment the matrix with the identity matrix.
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      augmented[i][j] = matrix[i][j];
      augmented[i][j + n] = (i == j) ? 1.0 : 0.0;
    }
  }

  // Perform Gauss-Jordan elimination.
  for (i=0; i<n; i++) {

    // Find the pivot element.
    pivot = augmented[i][i];
    if (fabs (pivot) < 1e-9) {  // If pivot is too small, the matrix is singular.
      fprintf (stdout, "ERROR: Singular matrix in gaussjordan().\n");
      exit (EXIT_FAILURE);
    }

    // Normalize the pivot row.
    for (j=0; j<(2 * n); j++) {
      augmented[i][j] /= pivot;
    }

    // Eliminate the current column in all rows except the current row.
    for (j=0; j<n; j++) {
      if (j != i) {
        factor = augmented[j][i];
        for (k=0; k<(2 * n); k++) {
          augmented[j][k] -= factor * augmented[i][k];
        }
      }
    }
  }

  // Copy the right half of the augmented matrix back to the original matrix.
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      matrix[i][j] = augmented[i][j + n];
    }
  }

  // Free allocated memory.
  for (i=0; i<n; i++) {
    free (augmented[i]);
  }
  free (augmented);

  return (EXIT_SUCCESS);
}

// Obtain a text string from standard input. It can include spaces.
int
inputtext (char *text) {

  // Request new text from standard input.
  fgets (text, MAXLEN, stdin);

  // Remove trailing newline, if there.
  if ((strnlen(text, MAXLEN) > 0) && (text[strnlen (text, MAXLEN) - 1] == '\n')) {
    text[strnlen (text, MAXLEN) - 1] = '\0';  // Replace newline with string termination.
  }

  return (EXIT_SUCCESS);
}

// Color primary coordinates for various RGB colorspaces
// References: Kang, Computational Color Technology (2006)
//             Pascale, A Review of RGB Color Spaces (2003)
int
rgb_primaries (double **p) {

  int choice;
  double xr, yr, zr, xg, yg, zg, xb, yb, zb;
  char *temp, *endptr;

  // Allocate memory for various arrays.
  temp = allocate_strmem (MAXLEN);

  fprintf (stdout, "\nChoose RGB colorspace:\n");
  fprintf (stdout, "  1 - sRGB (BT.709)\n");
  fprintf (stdout, "  2 - Adobe (1998)\n");
  fprintf (stdout, "  3 - Apple\n");
  fprintf (stdout, "  4 - Best\n");
  fprintf (stdout, "  5 - Beta\n");
  fprintf (stdout, "  6 - Bruse\n");
  fprintf (stdout, "  7 - CIE 2-deg observer\n");
  fprintf (stdout, "  8 - CIE 10-deg observer\n");
  fprintf (stdout, "  9 - ColorMatch\n");
  fprintf (stdout, " 10 - Don 4\n");
  fprintf (stdout, " 11 - EBU\n");
  fprintf (stdout, " 12 - ECI v2\n");
  fprintf (stdout, " 13 - Ekta Space PS5\n");
  fprintf (stdout, " 14 - Eureka\n");
  fprintf (stdout, " 15 - Extended\n");
  fprintf (stdout, " 16 - Guild\n");
  fprintf (stdout, " 17 - Ink-jet\n");
  fprintf (stdout, " 18 - Judd-Wyszecki\n");
  fprintf (stdout, " 19 - Kress\n");
  fprintf (stdout, " 20 - Laser (Starkweather)\n");
  fprintf (stdout, " 21 - NTSC (1953)\n");
  fprintf (stdout, " 22 - PAL / SECAM\n");
  fprintf (stdout, " 23 - ProPhoto\n");
  fprintf (stdout, " 24 - RIMM-ROMM\n");
  fprintf (stdout, " 25 - ROM\n");
  fprintf (stdout, " 26 - SGI\n");
  fprintf (stdout, " 27 - SMPTE-C (NTSC 1987)\n");
  fprintf (stdout, " 28 - SMPTE-240M\n");
  fprintf (stdout, " 29 - Sony P-22\n");
  fprintf (stdout, " 30 - Wide-Gamut\n");
  fprintf (stdout, " 31 - Wright\n");
  fprintf (stdout, " 32 - Usami\n");

  fprintf (stdout, "\nChoice? ");
  inputtext (temp);
  errno = 0;
  choice = (int) strtol (temp, &endptr, 10);
  if ((errno == ERANGE) || (errno == EINVAL) || (endptr == temp)) {
    fprintf (stderr, "ERROR: Cannot make integer of: %s\n", temp);
    exit (EXIT_FAILURE);
  }
  switch (choice) {

    // sRGB (BT.709)
    case 1:
      xr = 0.6400; yr = 0.3300;
      xg = 0.3000; yg = 0.6000;
      xb = 0.1500; yb = 0.0600;
      break;

    // Adobe (1998)
    case 2:
      xr = 0.6400; yr = 0.3300;
      xg = 0.2100; yg = 0.7100;
      xb = 0.1500; yb = 0.0600;
      break;

    // Apple
    case 3:
      xr = 0.6250; yr = 0.3400;
      xg = 0.2800; yg = 0.5950;
      xb = 0.1550; yb = 0.0700;
      break;

    // Best
    case 4:
      xr = 0.7347; yr = 0.2653;
      xg = 0.2150; yg = 0.7750;
      xb = 0.1300; yb = 0.0350;
      break;

    // Beta
    case 5:
      xr = 0.6888; yr = 0.3112;
      xg = 0.1986; yg = 0.7551;
      xb = 0.1265; yb = 0.0352;
      break;

    // Bruse
    case 6:
      xr = 0.6400; yr = 0.3300;
      xg = 0.2800; yg = 0.6500;
      xb = 0.1500; yb = 0.0600;
      break;

    // CIE 2-deg observer
    case 7:
      xr = 0.7347; yr = 0.2653;
      xg = 0.2737; yg = 0.7174;
      xb = 0.1665; yb = 0.0089;
      break;

    // CIE 10-deg observer
    case 8:
      xr = 0.7232; yr = 0.2768;
      xg = 0.1248; yg = 0.8216;
      xb = 0.1616; yb = 0.0134;
      break;

    // ColorMatch
    case 9:
      xr = 0.6300; yr = 0.3400;
      xg = 0.2950; yg = 0.6050;
      xb = 0.1500; yb = 0.0750;
      break;

    // Don
    case 10:
      xr = 0.6960; yr = 0.3000;
      xg = 0.2150; yg = 0.7650;
      xb = 0.1300; yb = 0.0350;
      break;

    // EBU
    case 11:
      xr = 0.6400; yr = 0.3300;
      xg = 0.2900; yg = 0.6000;
      xb = 0.1500; yb = 0.0600;
      break;

    // ECI v2
    case 12:
      xr = 0.6700; yr = 0.3300;
      xg = 0.2100; yg = 0.7100;
      xb = 0.1400; yb = 0.0800;
      break;

    // Ekta Space PS5
    case 13:
      xr = 0.6950; yr = 0.3050;
      xg = 0.2600; yg = 0.7000;
      xb = 0.1100; yb = 0.0050;
      break;

    // Eureka
    case 14:
      xr = 0.6915; yr = 0.3083;
      xg = 0.0000; yg = 1.0000;
      xb = 0.1440; yb = 0.0296;
      break;

    // Extended
    case 15:
      xr = 0.7010; yr = 0.2990;
      xg = 0.1700; yg = 0.7960;
      xb = 0.1310; yb = 0.0460;
      break;

    // Guild
    case 16:
      xr = 0.7000; yr = 0.3000;
      xg = 0.2550; yg = 0.7200;
      xb = 0.1500; yb = 0.0500;
      break;

    // Ink-jet
    case 17:
      xr = 0.7000; yr = 0.3000;
      xg = 0.2500; yg = 0.7200;
      xb = 0.1300; yb = 0.0500;
      break;

    // Judd-Wyszecki
    case 18:
      xr = 0.7347; yr = 0.2653;
      xg = 0.0743; yg = 0.8338;
      xb = 0.1741; yb = 0.0050;
      break;

    // Kress
    case 19:
      xr = 0.6915; yr = 0.3083;
      xg = 0.1547; yg = 0.8059;
      xb = 0.1440; yb = 0.0297;
      break;

    // Laser (Starkweather)
    case 20:
      xr = 0.7117; yr = 0.2882;
      xg = 0.0328; yg = 0.8029;
      xb = 0.1632; yb = 0.0119;
      break;

    // NTSC (1953)
    case 21:
      xr = 0.6700; yr = 0.3300;
      xg = 0.2100; yg = 0.7100;
      xb = 0.1400; yb = 0.0800;
      break;

    // PAL / SECAM
    case 22:
      xr = 0.6400; yr = 0.3300;
      xg = 0.2900; yg = 0.6000;
      xb = 0.1500; yb = 0.0600;
      break;

    // ProPhoto
    case 23:
      xr = 0.7347; yr = 0.2653;
      xg = 0.1596; yg = 0.8404;
      xb = 0.0366; yb = 0.0001;
      break;

    // RIMM-ROMM
    case 24:
      xr = 0.7347; yr = 0.2653;
      xg = 0.1596; yg = 0.8404;
      xb = 0.0366; yb = 0.0001;
      break;

    // ROM
    case 25:
      xr = 0.8730; yr = 0.1440;
      xg = 0.1750; yg = 0.9270;
      xb = 0.0850; yb = 0.0001;
      break;

    // SGI
    case 26:
      xr = 0.6250; yr = 0.3400;
      xg = 0.2800; yg = 0.5950;
      xb = 0.1550; yb = 0.0700;
      break;

    // SMPTE-C (NTSC 1987)
    case 27:
      xr = 0.6300; yr = 0.3400;
      xg = 0.3100; yg = 0.5950;
      xb = 0.1550; yb = 0.0700;
      break;

    // SMPTE-240M
    case 28:
      xr = 0.6700; yr = 0.3300;
      xg = 0.2100; yg = 0.7100;
      xb = 0.1500; yb = 0.0600;
      break;

    // Sony P-22
    case 29:
      xr = 0.6250; yr = 0.3400;
      xg = 0.2800; yg = 0.5950;
      xb = 0.1550; yb = 0.0700;
      break;

    // Wide-Gamut
    case 30:
      xr = 0.7347; yr = 0.2653;
      xg = 0.1152; yg = 0.8264;
      xb = 0.1566; yb = 0.0177;
      break;

    // Wright
    case 31:
      xr = 0.7260; yr = 0.2740;
      xg = 0.1547; yg = 0.8059;
      xb = 0.1440; yb = 0.0297;
      break;

    // Usami
    case 32:
      xr = 0.7347; yr = 0.2653;
      xg = -0.086; yg = 1.0860;
      xb = 0.0957; yb = -.0314;
      break;

    // Unknown
    default:
      fprintf (stderr, "Invalid choice.\n");
      return (-1);
  }

  // Compute z chromaticity coordinate of color primaries.
  zr = 1.0 - xr - yr;
  zg = 1.0 - xg - yg;
  zb = 1.0 - xb - yb;

  // Populate color primaries matrix p.
  p[0][0] = xr;  p[0][1] = xg;  p[0][2] = xb;
  p[1][0] = yr;  p[1][1] = yg;  p[1][2] = yb;
  p[2][0] = zr;  p[2][1] = zg;  p[2][2] = zb;

  // Free allocated memory.
  free (temp);

  return (0);  // Success
}

// Illuminants - Choose white point.
// Returns: -1 if invalid selection, 0 if valid selection

// References: 1. Danny Pascale, "A Review of RGB color spaces", Babel Color
//             2. Equivalent White Light Sources, and CIE Illuminants, HunterLab
//             3. CIE F-series Spectral Data, CIE 15.2:1986
//             4. Colorimetry, 4th Edition, CIE 015:2018, DOI: 10.25039/TR.015.2018
//             5. Tooms - Colour Reproduction in Electronic Imaging Stsrems (2015)
int
illum_white (double *white_xyz) {

  int choice;
  char *temp, *endptr;

  // Allocate memory for various arrays.
  temp = allocate_strmem (MAXLEN);

  fprintf (stdout, "\nChoose display white point for bitmap:\n");
  fprintf (stdout, "  1 - D65 1931 2-deg - 6504 K - Noon daylight (sRGB, BT.601, BT.709)\n");
  fprintf (stdout, "  2 - D65 1964 10-deg - 6504 K - Noon daylight\n");
  fprintf (stdout, "  3 - D50 1931 2-deg - 5003 K - Late afternoon\n");
  fprintf (stdout, "  4 - D50 1964 10-deg - 5003 K - Late afternoon\n");
  fprintf (stdout, "  5 - A 1931 2-deg - 2856 K - Tungsten filament\n");
  fprintf (stdout, "  6 - A 1964 10-deg - 2856 K - Tungsten filament\n");
  fprintf (stdout, "  7 - B 1931 2-deg - 4874 K - Noon daylight (obsolete)\n");
  fprintf (stdout, "  8 - C 1931 2-deg - 6774 K - Average daylight\n");
  fprintf (stdout, "  9 - C 1964 10-deg - 6774 K - Average daylight\n");
  fprintf (stdout, " 10 - D55 1931 2-deg - 5503 K - Mid-morning / Mid-afternoon\n");
  fprintf (stdout, " 11 - D55 1964 10-deg - 5503 K - Mid-morning / Mid-afternoon\n");
  fprintf (stdout, " 12 - D60 1931 2-deg - 5985 K - Daylight\n");
  fprintf (stdout, " 13 - D75 1931 2-deg - 7504 K - Northern daylight\n");
  fprintf (stdout, " 14 - D75 1964 10-deg - 7504 K - Northern daylight\n");
  fprintf (stdout, " 15 - D93 1931 2-deg - 9305 K - High-efficiency blue phosphor monitors, BT.2035, NTSC\n");
  fprintf (stdout, " 16 - D93 1964 10-deg - 9305 K - High-efficiency blue phosphor monitors, BT.2035, NTSC\n");
  fprintf (stdout, " 17 - E 1931 2-deg - 5454 K - Equal energy\n");
  fprintf (stdout, " 18 - E 1964 10-deg - 5454 K - Equal energy\n");
  fprintf (stdout, " 19 - F1 1931 2-deg - 6430 K - Florescent daylight\n");
  fprintf (stdout, " 20 - F1 1964 10-deg - 6430 K - Florescent daylight\n");
  fprintf (stdout, " 21 - F2 1931 2-deg - 4230 K - Cool white fluorescent\n");
  fprintf (stdout, " 22 - F2 1964 10-deg - 4230 K - Cool white fluorescent\n");
  fprintf (stdout, " 23 - F3 1931 2-deg - 3450 K - White fluorescent\n");
  fprintf (stdout, " 24 - F3 1964 10-deg - 3450 K - White fluorescent\n");
  fprintf (stdout, " 25 - F4 1931 2-deg - 2940 K - Warm white fluorescent\n");
  fprintf (stdout, " 26 - F4 1964 10-deg - 2940 K - Warm white fluorescent\n");
  fprintf (stdout, " 27 - F5 1931 2-deg - 6350 K - Daylight fluorescent\n");
  fprintf (stdout, " 28 - F5 1964 10-deg - 6350 K - Daylight fluorescent\n");
  fprintf (stdout, " 29 - F6 1931 2-deg - 4150 K - Light white fluorescent\n");
  fprintf (stdout, " 30 - F6 1964 10-deg - 4150 K - Light white fluorescent\n");
  fprintf (stdout, " 31 - F7 1931 2-deg - 6500 K - D65 daylight simulator\n");
  fprintf (stdout, " 32 - F7 1964 10-deg - 6500 K - D65 daylight simulator\n");
  fprintf (stdout, " 33 - F8 1931 2-deg - 5000 K - D50 simulator, Sylvania F40 Design 50\n");
  fprintf (stdout, " 34 - F8 1964 10-deg - 5000 K - D50 simulator, Sylvania F40 Design 50\n");
  fprintf (stdout, " 35 - F9 1931 2-deg - 4150 K - Cool white deluxe fluorescent\n");
  fprintf (stdout, " 36 - F9 1964 10-deg - 4150 K - Cool white deluxe fluorescent\n");
  fprintf (stdout, " 37 - F10 1931 2-deg - 5000 K - Philips TL85, Ultralume 50\n");
  fprintf (stdout, " 38 - F10 1964 10-deg - 5000 K - Philips TL85, Ultralume 50\n");
  fprintf (stdout, " 39 - F11 1931 2-deg - 4000 K - Philips TL84, Ultralume 40\n");
  fprintf (stdout, " 40 - F11 1964 10-deg - 4000 K - Philips TL84, Ultralume 40\n");
  fprintf (stdout, " 41 - F12 1931 2-deg - 3000 K - Philips TL83, Ultralume 30\n");
  fprintf (stdout, " 42 - F12 1964 10-deg - 3000 K - Philips TL83, Ultralume 30\n");
  fprintf (stdout, "\nChoice? ");
  inputtext (temp);
  errno = 0;
  choice = (int) strtol (temp, &endptr, 10);
  if ((errno == ERANGE) || (errno == EINVAL) || (endptr == temp)) {
    fprintf (stderr, "ERROR: Cannot make integer of: %s\n", temp);
    exit (EXIT_FAILURE);
  }
  switch (choice) {

    // D65 1931 2-deg - 6504 K - Noon daylight, BT.601, BT.709, sRGB
    // This is the sRGB standard white.
    case 1:
      white_xyz[0] = 0.31272;  // x
      white_xyz[1] = 0.32903;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // D65 1964 10-deg - 6504 K - Noon daylight
    case 2:
      white_xyz[0] = 0.31382;  // x
      white_xyz[1] = 0.33100;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // D50 1931 2-deg - 5003 K - Late afternoon
    case 3:
      white_xyz[0] = 0.34567;  // x
      white_xyz[1] = 0.35850;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // D50 1964 10-deg - 5003 K - Late afternoon
    case 4:
      white_xyz[0] = 0.34773;  // x
      white_xyz[1] = 0.35952;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // A 1931 2-deg - 2856 K - Tungsten filament
    case 5:
      white_xyz[0] = 0.44758;  // x
      white_xyz[1] = 0.40745;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // A 1964 10-deg - 2856 K - Tungsten filament
    case 6:
      white_xyz[0] = 0.45117;  // x
      white_xyz[1] = 0.40594;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // B 1931 2-deg - 4874 K - Noon daylight (obsolete)
    case 7:
      white_xyz[0] = 0.34830;  // x
      white_xyz[1] = 0.35160;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // C 1931 2-deg - 6774 K - Average daylight
    case 8:
      white_xyz[0] = 0.31006;  // x
      white_xyz[1] = 0.31615;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // C 1964 10-deg - 6774 K - Average daylight
    case 9:
      white_xyz[0] = 0.31039;  // x
      white_xyz[1] = 0.31905;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // D55 1931 2-deg - 5503 K - Mid-morning / Mid-afternoon
    case 10:
      white_xyz[0] = 0.33242;  // x
      white_xyz[1] = 0.34743;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // D55 1964 10-deg - 5503 K - Mid-morning / Mid-afternoon
    case 11:
      white_xyz[0] = 0.33411;  // x
      white_xyz[1] = 0.34877;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // D60 1931 2-deg - 5985 K - Daylight
    case 12:
      white_xyz[0] = 0.3217;
      white_xyz[1] = 0.3377;
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // D75 1931 2-deg - 7504 K - Northern daylight
    case 13:
      white_xyz[0] = 0.29902;  // x
      white_xyz[1] = 0.31485;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // D75 1964 10-deg - 7504 K - Northern daylight
    case 14:
      white_xyz[0] = 0.29968;  // x
      white_xyz[1] = 0.31740;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // D93 1931 2-deg - 9305 K - High-efficiency blue phosphor monitors, BT.2035, NTSC
    case 15:
      white_xyz[0] = 0.28315;  // x
      white_xyz[1] = 0.29711;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // D93 1964 10-deg - 9305 K - High-efficiency blue phosphor monitors, BT.2035, NTSC
    case 16:
      white_xyz[0] = 0.28327;  // x
      white_xyz[1] = 0.30043;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // E 1931 2-deg - 5454 K - Equal energy
    case 17:
      white_xyz[0] = 0.33333;  // x
      white_xyz[1] = 0.33333;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // E 1964 10-deg - 5454 K - Equal energy
    case 18:
      white_xyz[0] = 0.33333;  // x
      white_xyz[1] = 0.33333;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // F1 1931 2-deg - 6430 K - Florescent daylight
    case 19:
      white_xyz[0] = 0.31310;  // x
      white_xyz[1] = 0.33727;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // F1 1964 10-deg - 6430 K - Florescent daylight
    case 20:
      white_xyz[0] = 0.31811;  // x
      white_xyz[1] = 0.33559;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // F2 1931 2-deg - 4230 K - Cool white fluorescent
    case 21:
      white_xyz[0] = 0.37208;  // x
      white_xyz[1] = 0.37529;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // F2 1964 10-deg - 4230 K - Cool white fluorescent
    case 22:
      white_xyz[0] = 0.37925;  // x
      white_xyz[1] = 0.36733;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // F3 1931 2-deg - 3450 K - White fluorescent
    case 23:
      white_xyz[0] = 0.40910;  // x
      white_xyz[1] = 0.39430;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // F3 1964 10-deg - 3450 K - White fluorescent
    case 24:
      white_xyz[0] = 0.41761;  // x
      white_xyz[1] = 0.38324;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // F4 1931 2-deg - 2940 K - Warm white fluorescent
    case 25:
      white_xyz[0] = 0.44018;  // x
      white_xyz[1] = 0.40329;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // F4 1964 10-deg - 2940 K - Warm white fluorescent
    case 26:
      white_xyz[0] = 0.44920;  // x
      white_xyz[1] = 0.39074;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // F5 1931 2-deg - 6350 K - Daylight fluorescent
    case 27:
      white_xyz[0] = 0.31379;  // x
      white_xyz[1] = 0.34531;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // F5 1964 10-deg - 6350 K - Daylight fluorescent
    case 28:
      white_xyz[0] = 0.31975;  // x
      white_xyz[1] = 0.34246;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // F6 1931 2-deg - 4150 K - Light white fluorescent
    case 29:
      white_xyz[0] = 0.37790;  // x
      white_xyz[1] = 0.38835;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // F6 1964 10-deg - 4150 K - Light white fluorescent
    case 30:
      white_xyz[0] = 0.38660;  // x
      white_xyz[1] = 0.37847;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // F7 1931 2-deg - 6500 K - D65 daylight simulator
    case 31:
      white_xyz[0] = 0.31292;  // x
      white_xyz[1] = 0.32933;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // F7 1964 10-deg - 6500 K - D65 daylight simulator
    case 32:
      white_xyz[0] = 0.31569;  // x
      white_xyz[1] = 0.32960;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // F8 1931 2-deg - 5000 K - D50 simulator, Sylvania F40 Design 50
    case 33:
      white_xyz[0] = 0.34588;  // x
      white_xyz[1] = 0.35875;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // F8 1964 10-deg - 5000 K - D50 simulator, Sylvania F40 Design 50
    case 34:
      white_xyz[0] = 0.34902;  // x
      white_xyz[1] = 0.35939;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // F9 1931 2-deg - 4150 K - Cool white deluxe fluorescent
    case 35:
      white_xyz[0] = 0.37417;  // x
      white_xyz[1] = 0.37281;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // F9 1964 10-deg - 4150 K - Cool white deluxe fluorescent
    case 36:
      white_xyz[0] = 0.37829;  // x
      white_xyz[1] = 0.37045;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // F10 1931 2-deg - 5000 K - Philips TL85, Ultralume 50
    case 37:
      white_xyz[0] = 0.34609;  // x
      white_xyz[1] = 0.35986;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // F10 1964 10-deg - 5000 K - Philips TL85, Ultralume 50
    case 38:
      white_xyz[0] = 0.35090;  // x
      white_xyz[1] = 0.35444;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // F11 1931 2-deg - 4000 K - Philips TL84, Ultralume 40
    case 39:
      white_xyz[0] = 0.38052;  // x
      white_xyz[1] = 0.37713;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // F11 1964 10-deg - 4000 K - Philips TL84, Ultralume 40
    case 40:
      white_xyz[0] = 0.38541;  // x
      white_xyz[1] = 0.37123;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // F12 1931 2-deg - 3000 K - Philips TL83, Ultralume 30
    case 41:
      white_xyz[0] = 0.43695;  // x
      white_xyz[1] = 0.40441;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    // F12 1964 10-deg - 3000 K - Philips TL83, Ultralume 30
    case 42:
      white_xyz[0] = 0.44256;  // x
      white_xyz[1] = 0.39717;  // y
      white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];  // z
      break;

    default:
      fprintf (stderr, "Invalid selection.\n");
      return (-1);
  }

  // Free allocated memory.
  free (temp);

  return (0);  // Success
}

// Allocate memory for an array of ints.
int *
allocate_intmem (int len) {

  void *tmp;

  if (len <= 0) {
    fprintf (stderr, "ERROR: Cannot allocate memory because len = %i in allocate_intmem().\n", len);
    exit (EXIT_FAILURE);
  }

  tmp = (int *) malloc (len * sizeof (int));
  if (tmp != NULL) {
    memset (tmp, 0, len * sizeof (int));
    return (tmp);
  } else {
    fprintf (stderr, "ERROR: Cannot allocate memory for array in allocate_intmem().\n");
    exit (EXIT_FAILURE);
  }
}

// Allocate memory for an array of doubles.
double *
allocate_doublemem (int len)
{ 
  void *tmp;
  
  if (len <= 0) { 
    fprintf (stderr, "ERROR: Cannot allocate memory because len = %i in allocate_doublemem().\n", len);
    exit (EXIT_FAILURE);
  }
  
  tmp = (double *) malloc (len * sizeof (double));
  if (tmp != NULL) {
    memset (tmp, 0, len * sizeof (double));
    return (tmp);
  } else {
    fprintf (stderr, "ERROR: Cannot allocate memory for array in allocate_doublemem().\n");
    exit (EXIT_FAILURE);
  }
} 
  
// Allocate memory for an array of pointers to arrays of doubles.
double **
allocate_doublememp (int len)
{ 
  void *tmp;
  
  if (len <= 0) {
    fprintf (stderr, "ERROR: Cannot allocate memory because len = %i in allocate_doublememp().\n", len);
    exit (EXIT_FAILURE);
  }
  
  tmp = (double **) malloc (len * sizeof (double *));
  if (tmp != NULL) {
    memset (tmp, 0, len * sizeof (double *));
    return (tmp);
  } else {
    fprintf (stderr, "ERROR: Cannot allocate memory for array in allocate_doublememp().\n");
    exit (EXIT_FAILURE);
  } 
}

// Allocate memory for an array of chars.
char *
allocate_strmem (int len) {

  void *tmp;
    
  if (len <= 0) {
    fprintf (stderr, "ERROR: Cannot allocate memory because len = %i in allocate_strmem().\n", len);
    exit (EXIT_FAILURE);
  }
  
  tmp = (char *) malloc (len * sizeof (char));
  if (tmp != NULL) {
    memset (tmp, 0, len * sizeof (char));
    return (tmp);
  } else {
    fprintf (stderr, "ERROR: Cannot allocate memory for array in allocate_strmem().\n");
    exit (EXIT_FAILURE);
  }
}
