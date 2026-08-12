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

// wl.c - Provide color-matching function (CMF) coordinates, corresponding chromaticities, and sRGB coordinates for a given wavelength.
//        Note that none of the pure wavelength stimuli can actually be reproduced by sRGB; the sRGB color gamut triangle lies within
//        the chromaticity locus. The results here are therefore not equivalent to the actual monochromatic stimuli.

// gcc -Wall src/wl.c -lm -o wl

// Usage: ./wl
// Inputs: CMF files:
//         data/CIE_xyz_1931_2deg.csv
//         data/CIE_xyz_1931_2deg_judd1951.csv
//         data/CIE_xyz_1931_2deg_judd1951_vos1978.csv
//         data/CIE_xyz_1964_10deg.csv
//         data/CIE_xyz_2006_10deg_lms_cones.csv
//         data/CIE_xyz_2006_2deg_lms_cones.csv
// Output: reports to stdout

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>  // uint8_t
#include <string.h>
#include <math.h>  // for fabs(), pow()
#include <errno.h>

// Function prototypes
int inputtext (char *);
int readline (FILE*, char*, int);
int gaussjordan (int, double **);
int cmf (double, int, double **, double *);
int illum_white (double *);
int *allocate_intmem (int);
char *allocate_strmem (int);
double *allocate_doublemem (int);
double **allocate_doublememp (int);
uint8_t *allocate_ustrmem (int);

// Set some symbolic constants.
#define MAXLEN 256  // Maximum number of characters per line

int
main (int argc, char **argv) {

  int choice, i, j, k, nlines, *rgb, apply_gamma;
  double *white_xyz, **cmxyz, lambda, *xyzbar, *xyz;
  double zr, zg, zb, *rgb_double, **p, *w, **pinv, **c, *coeff, **npm, **npminv;
  char *filename, *temp, *endptr;
  FILE *fi;

  // Allocate memory for various arrays.
  filename = allocate_strmem (MAXLEN);
  temp = allocate_strmem (MAXLEN);
  xyzbar = allocate_doublemem (3);
  xyz = allocate_doublemem (3);
  rgb = allocate_intmem (3);
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
  rgb_double = allocate_doublemem (3);
  white_xyz = allocate_doublemem (3);

  // Choose color-matching function (CMF).
  fprintf (stdout, "\nChoose color-matching function (CMF):\n\n");
  fprintf (stdout, "  1 - 1964 10-deg XYZ CMFs (JIS Z 8701:1999)\n");
  fprintf (stdout, "  2 - 1931 2-deg XYZ CIE CMFs (CIE.15.2004)\n");
  fprintf (stdout, "  3 - 1931 2-deg XYZ CIE CMFs with Judd (1951) modifications\n");
  fprintf (stdout, "  4 - 1931 2-deg XYZ CIE CMFs with Judd (1951) and Vos (1978) modifications\n");
  fprintf (stdout, "  5 - 2006 2-deg XYZ CMFs transformed from the CIE (2006) 2-deg LMS cone fundamentals\n");
  fprintf (stdout, "  6 - 2006 10-deg XYZ CMFs transformed from the CIE (2006) 10-deg LMS cone fundamentals\n");
  fprintf (stdout, "  7 - Enter filename for CMFs (nm, xbar, ybar, zbar as .csv)\n");
  fprintf (stdout, "\nChoice? ");
  memset (temp, 0, MAXLEN * sizeof (char));
  inputtext (temp);
  errno = 0;
  choice = (int) strtol (temp, &endptr, 10);
  if ((errno == ERANGE) || (errno == EINVAL) || (endptr == temp)) {
    fprintf (stderr, "ERROR: Cannot make integer of: %s\n", temp);
    exit (EXIT_FAILURE);
  }
  switch (choice) {

    case 1:
      strncpy (filename, "data/CIE_xyz_1964_10deg.csv", MAXLEN);  // 1964 10-deg XYZ CMFs (JIS Z 8701:1999)
      break;
    case 2:
      strncpy (filename, "data/CIE_xyz_1931_2deg.csv", MAXLEN);  // 1931 2-deg XYZ CIE CMFs (CIE.15.2004)
      break;
    case 3:
      strncpy (filename, "data/CIE_xyz_1931_2deg_judd1951.csv", MAXLEN);  // 1931 2-deg XYZ CIE CMFs with Judd (1951) modifications
      break;
    case 4:
      strncpy (filename, "data/CIE_xyz_1931_2deg_judd1951_vos1978.csv", MAXLEN);  // 1931 2-deg XYZ CIE CMFs with Judd (1951) and Vos (1978) modifications
      break;
    case 5:
      strncpy (filename, "data/CIE_xyz_2006_2deg_lms_cones.csv", MAXLEN);  // 2006 2-deg XYZ CMFs transformed from the CIE (2006) 2-deg LMS cone fundamentals
      break;
    case 6:
      strncpy (filename, "data/CIE_xyz_2006_10deg_lms_cones.csv", MAXLEN);  // 2006 10-deg XYZ CMFs transformed from the CIE (2006) 10-deg LMS cone fundamentals
      break;
    case 7:
      fprintf (stdout, "Filename for csv CMFs? ");
      memset (temp, 0, MAXLEN * sizeof (char));
      inputtext (filename);
      break;
    default:
      fprintf (stderr, "Unknown CMF choice.\n");
      exit (EXIT_FAILURE);
  }

  // Open color-matching functions csv file.
  fi = fopen (filename, "r");
  if (fi == NULL) {
    fprintf (stderr, "ERROR: Unable to open input file \"%s\".\n", filename);
    exit (EXIT_FAILURE);
  }

  // Count lines in input file.
  nlines = 0;
  for (;;) {

    if (readline (fi, temp, MAXLEN) == -1) {
      break;  // Reached end of file.
    }

    if ((temp[0] >= '0') && (temp[0]<= '9')) nlines++;

  }  // Next line of input file.
  fprintf (stdout, "%i lines in color-matching file: %s\n", nlines, filename);
  rewind (fi);

  // Allocate memory for various arrays.
  cmxyz = allocate_doublememp (nlines);
  for (i=0; i<nlines; i++) {
    cmxyz[i] = allocate_doublemem (4);  // lambda, xbar, ybar, zbar
  }

  // Read color-matching functions into array cmxyz.
  // cmxyz array: lambda, xbar, ybar, zbar.
  for (i=0; i<nlines; i++) {

    // Read line from input file.
    memset (temp, 0, MAXLEN * sizeof (char));
    if (readline (fi, temp, MAXLEN) == -1) {
      fprintf (stderr, "ERROR: Cannot read color-matching file \"%s\".\n", filename);
      exit (EXIT_FAILURE);
    }

    // Extract values from line of text.
    sscanf (temp, "%lf %lf %lf %lf", &cmxyz[i][0], &cmxyz[i][1], &cmxyz[i][2], &cmxyz[i][3]);
//printf ("%.0lf %.12lf %.12lf %.12lf\n", cmxyz[i][0], cmxyz[i][1], cmxyz[i][2], cmxyz[i][3]);

  }  // Next line

  // Close input file.
  fclose (fi);

  // sRGB uses the same color primaries as BT.709.

  // Red BT.709 color primary coordinates in xyz.
  const double xr = 0.640;
  const double yr = 0.330;
  zr = 1.0 - (xr + yr);

  // Green BT.709 color primary coordinates in xyz.
  const double xg = 0.300;
  const double yg = 0.600;
  zg = 1.0 - (xg + yg);

  // Blue BT.709 color primary coordinates in xyz.
  const double xb = 0.150;
  const double yb = 0.060;
  zb = 1.0 - (xb + yb);

  // Choose white point for display.
  // NOTE: sRGB standard uses D65, so this is generally the most appropriate. Other options are for interest.
  illum_white (white_xyz);

  // Populate color primaries matrix p with BT.709 rgb primary coordinates.
  p[0][0] = xr;  p[0][1] = xg;  p[0][2] = xb;
  p[1][0] = yr;  p[1][1] = yg;  p[1][2] = yb;
  p[2][0] = zr;  p[2][1] = zg;  p[2][2] = zb;

  // Populate white matrix with XYZ coordinates.
  // White matrix is normalized so that Y = 1.0.
  w[0] = white_xyz[0] / white_xyz[1];
  w[1] = 1.0;
  w[2] = white_xyz[2] / white_xyz[1];

  // Calculate BT.709 color primaries matrix p.
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      pinv[i][j] = p[i][j];  // Copy matrix p to matrix pinv for later in-place inversion.
    }
  }

  // Compute inverse of color primaries matrix p.
  gaussjordan (3, pinv);

  // Calculate RGB normalization coefficients.
  // coef = pinv * w
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      coeff[i] += pinv[i][j] * w[j];
    }
  }

  // Create diagonal coefficient matrix using normalization coefficients.
  c[0][0] = coeff[0];
  c[1][1] = coeff[1];
  c[2][2] = coeff[2];

  // Compute normalized primary matrix (NPM) matrix, where NPM = P * C.
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      for (k=0; k<3; k++) {
        npm[i][j] += p[i][k] * c[j][k];
      }
    }
  }

  // Find inverse of NPM matrix in order to convert XYZ to sRGB tristimulus values.
  // Copy npm to npminv first.
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      npminv[i][j] = npm[i][j];
    }
  }
  gaussjordan (3, npminv);

  // Show inverse NPM matrix to convert XYZ to RGB.
  // RGB = npminv * XYZ
/*
  fprintf (stdout, "Inverse NPM matrix to convert XYZ to RGB:\n");
  for (i=0; i<3; i++) {
    fprintf (stdout, "  ");
    for (j=0; j<3; j++) {
      fprintf (stdout, "%0.7lf ", npminv[i][j]);
    }
    fprintf (stdout, "\n");
  }
  fprintf (stdout, "\n");
*/

  // Ask for wavelength.
  fprintf (stdout, "\nWavelength (nm)? ");
  memset (temp, 0, MAXLEN * sizeof (char));
  inputtext (temp);
  errno = 0;
  lambda = strtol (temp, &endptr, 10);
  if ((errno == ERANGE) || (errno == EINVAL) || (endptr == temp)) {
    fprintf (stderr, "ERROR: Cannot make double of: %s\n", temp);
    exit (EXIT_FAILURE);
  }

  // Check range of requested wavelength.
  if ((lambda < cmxyz[0][0]) || (lambda > cmxyz[nlines-1][0])) {
    fprintf (stderr, "ERROR: Requested wavelength %.3lf is outside of range of color-matching functions.\n", lambda);
    fprintf (stderr, "       Available range is from %.3lf to %.3lf.\nz", cmxyz[0][0], cmxyz[nlines-1][0]);
    exit (EXIT_FAILURE);
  }

  // Ask whether to apply sRGB gamma-correction.
  fprintf (stdout, "\nApply sRGB gamma-correction (should be applied for standard sRGB image display)? ");
  memset (temp, 0, MAXLEN * sizeof (char));
  inputtext (temp);
  if ((temp[0] == 'y') || (temp[0] == 'Y')) {
    apply_gamma = 1;
  } else {
    apply_gamma = 0;
  }

  // Retrieve wavelength and color-matching coodinates xbar, ybar, zbar for requested wavelength.
  // Note: For single wavelength, the CMF X,Y,Z coordinates are the linear scene tristimulus values. i.e., no need to integrate
  // xbar = xyzbar[0], ybar = xyzbar[1], zbar = xyzbar[2]
  if (cmf (lambda, nlines, cmxyz, xyzbar) == EXIT_FAILURE) {
    fprintf (stderr, "ERROR: lambda %.18lf is outside range of CMF data: %.18lf to %.18lf.\n", lambda, cmxyz[0][0], cmxyz[nlines-1][0]);
    exit (EXIT_FAILURE);
  }

  fprintf (stdout, "\nXYZ CMF coordinates: %.12lf %.12lf %.12lf\n", xyzbar[0], xyzbar[1], xyzbar[2]);

  // Compute spectral chromaticity coordinates xyz.
  xyz[0] = xyzbar[0] / (xyzbar[0] + xyzbar[1] + xyzbar[2]);
  xyz[1] = xyzbar[1] / (xyzbar[0] + xyzbar[1] + xyzbar[2]);
  xyz[2] = xyzbar[2] / (xyzbar[0] + xyzbar[1] + xyzbar[2]);

  fprintf (stdout, "\nSpectral chromaticity coordinates xyz: %.12lf %.12lf %.12lf\n", xyz[0], xyz[1], xyz[2]);

  // Calculate linear srgb values from xyz coordinates.
  // rgb = npminv * XYZ
  memset (rgb_double, 0.0, 3 * sizeof (double));
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      rgb_double[i] += npminv[i][j] * xyzbar[j];
    }
  }
  fprintf (stdout, "\nRGB coordinates (no gamma correction, no scaling, no clipping): %.12lf %.12lf %.12lf\n", rgb_double[0], rgb_double[1], rgb_double[2]);
  fprintf (stdout, "\nRGB coordinates (no gamma correction, with scaling, no clipping): %.12lf %.12lf %.12lf\n", rgb_double[0] * 255.0, rgb_double[1] * 255.0, rgb_double[2] * 255.0);

  // Apply sRGB gamma-correction to linear rgb, if requested.
  if (apply_gamma) {
    for (i=0; i<3; i++) {
      if (rgb_double[i] > 0.0031308) {
        rgb_double[i] = (1.055 * pow (rgb_double[i], (1.0 / 2.4))) - 0.055;
      } else {
        rgb_double[i] *= 12.92;
      }
    }
  }
  fprintf (stdout, "\nRGB coordinates (with gamma correction, no scaling, no clipping): %.12lf %.12lf %.12lf\n", rgb_double[0], rgb_double[1], rgb_double[2]);

  // Scale to achieve required range 0-255.
  for (i=0; i<3; i++) {
    rgb_double[i] *= 255.0;
  }

  // Convert to integer.
  for (i=0; i<3; i++) {
    rgb[i] = (int) (rgb_double[i] + 0.5);
  }

  // Clip any undershoot or overshoot resulting from the fact that
  // 8-bit RGB (0-255) has a somewhat different color gamut than XYZ.
  for (i=0; i<3; i++) {
    if (rgb[i] < 0) rgb[i] = 0;
    if (rgb[i] > 255) rgb[i] = 255;
  }

  fprintf (stdout, "\nRGB coordinates (with sRGB gamma-correction): (%i, %i, %i)\n", rgb[0], rgb[1], rgb[2]);
  fprintf (stdout, "\n");

  // Free allocated memory.
  free (filename);
  free (temp);
  free (white_xyz);
  free (xyzbar);
  free (xyz);
  free (rgb);
  for (i=0; i<nlines; i++) {
    free (cmxyz[i]);
  }
  free (cmxyz);
  for (i=0; i<3; i++) {
    free (p[i]);
    free (pinv[i]);
    free (c[i]);
    free (npm[i]);
    free (npminv[i]);
  }
  free (p);
  free (pinv);
  free (w);
  free (coeff);
  free (c);
  free (npm);
  free (npminv);
  free (rgb_double);

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

// Read a single line of text from a csv text file.
// Convert commas to spaces.
// Returns -1 if EOF is encountered.
int
readline (FILE *fi, char *line, int limit) {

  int i, n;

  i = 0;  // i is pointer to byte in line.
  while (i < limit) {

    // Grab next byte from file.
    n = fgetc (fi);

    // End of file reached.
    // Tell calling function, by returning -1, that we're at end of file, so it won't call readline() again.
    if (n == EOF) {

      // If there's no end of line at the end of the file, ensure string termination.
      if (i > 0) {
        line[i] = 0;
        return (0);
      }
      return (-1);
    }

    // Found a carriage return. Ignore it.
    if (n == '\r') {
      continue;
    }

    // Found a comma; convert to space.
    if (n == ',') {
      n = ' ';
    }

    // Found a newline. Change to 0 for string termination.
    // Break out of loop since this is the end of the current line.
    if (n == '\n') {
      line[i] = 0;  // Replace with 0 for string termination.
      return (0);
    }

    // Seems to be a valid character. Keep it.
    line[i] = n;
    i++;
  }

  // Advance to next line.
  n = 0;
  while ((n != '\n') && (n != EOF)) {
    n = fgetc (fi);
  }

  return (0);
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

// Return xyz color-matching function (CMF) coordinates by interpolating the data.
int
cmf (double lambda, int nlines, double **cmxyz, double *xyzbar) {

  int i, lb, ub;
  double frac;

  // If lambda is above or below table data, return error.
  if ((lambda < cmxyz[0][0]) || (lambda > cmxyz[nlines-1][0])) {
    return (EXIT_FAILURE);

  // Otherwise, interpolate if necessary.
  } else {
 
    i = 0;  // Index of wavelengths in cmxyz. i.e., cmxyz[i][0] = wavelength i
    while ((cmxyz[i][0] < lambda) && (i < nlines)) i++;

    // Lower bound
    if (i > 0) {
      lb = i - 1;
    } else {
      lb = i;
    }

    // Upper bound
    ub = lb + 1;

    frac = (lambda - cmxyz[lb][0]) / (cmxyz[ub][0] - cmxyz[lb][0]);

    // Interpolate CMF coordinates.
    xyzbar[0] = (frac * (cmxyz[ub][1] - cmxyz[lb][1])) + cmxyz[lb][1];
    xyzbar[1] = (frac * (cmxyz[ub][2] - cmxyz[lb][2])) + cmxyz[lb][2];
    xyzbar[2] = (frac * (cmxyz[ub][3] - cmxyz[lb][3])) + cmxyz[lb][3];
  }

  return (EXIT_SUCCESS);
}

// Illuminants - Choose white point.
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
  fprintf (stdout, "NOTE: sRGB standard uses 1931 2-deg D65, so this is most appropriate for sRGB bitmap output. Other options are for interest.\n");
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
      exit (EXIT_FAILURE);
  }

  // Free allocated memory.
  free (temp);

  return (EXIT_SUCCESS);
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

// Allocate memory for an array of unsigned chars.
uint8_t *
allocate_ustrmem (int len) {
    
  void *tmp;
  
  if (len <= 0) {
    fprintf (stderr, "ERROR: Cannot allocate memory because len = %i in allocate_ustrmem().\n", len);
    exit (EXIT_FAILURE);
  }

  tmp = (uint8_t *) malloc (len * sizeof (uint8_t));
  if (tmp != NULL) {
    memset (tmp, 0, len * sizeof (uint8_t));
    return (tmp);
  } else {
    fprintf (stderr, "ERROR: Cannot allocate memory for array in allocate_ustrmem().\n");
    exit (EXIT_FAILURE);
  }
}
