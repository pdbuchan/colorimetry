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
#include <ctype.h>
#include <limits.h>
#include <float.h>

// Function prototypes
int inputtext (char *);
int parse_int_string (const char *, int *);
int parse_double_string (const char *, double *);
int parse_cmf_record (const char *, double *);
int readline (FILE *, char *, int);
int count_cmf_rows (const char *);
int load_cmf (int, const char *, double **);
int gaussjordan (int, double **);
int cmf (double, int, double **, double *);
int illum_white (double *);
int *allocate_intmem (size_t);
char *allocate_strmem (size_t);
double *allocate_doublemem (size_t);
double **allocate_doublememp (size_t);

// Set some symbolic constants.
#define MAX_STRINGLEN 256  // Maximum number of characters per line

typedef struct {
  const char *name;
  double x, y;
  const char *source_note;
} WhitePoint;

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

  int choice, i, j, nlines, *rgb, apply_gamma;
  double *white_xyz, **cmxyz, lambda, *xyzbar, *xyz;
  double zr, zg, zb, *rgb_double, **p, *w, **pinv, *coeff, **npm, **npminv;
  char *filename, *temp;

  // Allocate memory for various arrays.
  filename = allocate_strmem (MAX_STRINGLEN);
  temp = allocate_strmem (MAX_STRINGLEN);
  xyzbar = allocate_doublemem (3);
  xyz = allocate_doublemem (3);
  rgb = allocate_intmem (3);
  p = allocate_doublememp (3);
  pinv = allocate_doublememp (3);
  npm = allocate_doublememp (3);
  npminv = allocate_doublememp (3);
  for (i=0; i<3; i++) {
    p[i] = allocate_doublemem (3);
    pinv[i] = allocate_doublemem (3);
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
  memset (temp, 0, MAX_STRINGLEN * sizeof (char));
  inputtext (temp);
  if (parse_int_string (temp, &choice) == EXIT_FAILURE) {
    fprintf (stderr, "ERROR: Cannot make integer of: %s\n", temp);
    exit (EXIT_FAILURE);
  }
  switch (choice) {

    case 1:
      snprintf (filename, MAX_STRINGLEN, "data/CIE_xyz_1964_10deg.csv");  // 1964 10-deg XYZ CMFs (JIS Z 8701:1999)
      break;
    case 2:
      snprintf (filename, MAX_STRINGLEN, "data/CIE_xyz_1931_2deg.csv");  // 1931 2-deg XYZ CIE CMFs (CIE.15.2004)
      break;
    case 3:
      snprintf (filename, MAX_STRINGLEN, "data/CIE_xyz_1931_2deg_judd1951.csv");  // 1931 2-deg XYZ CIE CMFs with Judd (1951) modifications
      break;
    case 4:
      snprintf (filename, MAX_STRINGLEN, "data/CIE_xyz_1931_2deg_judd1951_vos1978.csv");  // 1931 2-deg XYZ CIE CMFs with Judd (1951) and Vos (1978) modifications
      break;
    case 5:
      snprintf (filename, MAX_STRINGLEN, "data/CIE_xyz_2006_2deg_lms_cones.csv");  // 2006 2-deg XYZ CMFs transformed from the CIE (2006) 2-deg LMS cone fundamentals
      break;
    case 6:
      snprintf (filename, MAX_STRINGLEN, "data/CIE_xyz_2006_10deg_lms_cones.csv");  // 2006 10-deg XYZ CMFs transformed from the CIE (2006) 10-deg LMS cone fundamentals
      break;
    case 7:
      fprintf (stdout, "Filename for csv CMFs? ");
      memset (temp, 0, MAX_STRINGLEN * sizeof (char));
      inputtext (filename);
      break;
    default:
      fprintf (stderr, "Unknown CMF choice.\n");
      exit (EXIT_FAILURE);
  }

  // Count and load validated CMF records. Text headers, comments, and blank
  // lines are ignored consistently in both passes.
  nlines = count_cmf_rows (filename);
  if (nlines < 0) exit (EXIT_FAILURE);
  if (nlines < 2) {
    fprintf (stderr, "ERROR: Color-matching file must contain at least two numeric rows.\n");
    exit (EXIT_FAILURE);
  }
  fprintf (stdout, "%i lines in color-matching file: %s\n", nlines, filename);

  cmxyz = allocate_doublememp ((size_t) nlines);
  for (i=0; i<nlines; i++) {
    cmxyz[i] = allocate_doublemem (4u);  // lambda, xbar, ybar, zbar
  }
  if (load_cmf (nlines, filename, cmxyz) == EXIT_FAILURE) {
    fprintf (stderr, "ERROR: Unable to load color-matching file %s.\n", filename);
    exit (EXIT_FAILURE);
  }

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
  for (;;) {
    if (illum_white (white_xyz) > -1) break;
  }

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
  if (gaussjordan (3, pinv) == EXIT_FAILURE) {
    fprintf (stderr, "ERROR: Unable to invert sRGB/BT.709 primary matrix.\n");
    exit (EXIT_FAILURE);
  }

  // Calculate RGB normalization coefficients.
  // coef = pinv * w
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      coeff[i] += pinv[i][j] * w[j];
    }
  }

  // Compute NPM = P * C. Since the normalization matrix C is diagonal,
  // scale each primary column directly by its coefficient.
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      npm[i][j] = p[i][j] * coeff[j];
    }
  }

  // Find inverse of NPM matrix in order to convert XYZ to sRGB tristimulus values.
  // Copy npm to npminv first.
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      npminv[i][j] = npm[i][j];
    }
  }
  if (gaussjordan (3, npminv) == EXIT_FAILURE) {
    fprintf (stderr, "ERROR: Unable to invert normalized primary matrix.\n");
    exit (EXIT_FAILURE);
  }

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
  memset (temp, 0, MAX_STRINGLEN * sizeof (char));
  inputtext (temp);
  if (parse_double_string (temp, &lambda) == EXIT_FAILURE) {
    fprintf (stderr, "ERROR: Cannot make finite wavelength of: %s\n", temp);
    exit (EXIT_FAILURE);
  }

  // Check range of requested wavelength.
  if ((lambda < cmxyz[0][0]) || (lambda > cmxyz[nlines-1][0])) {
    fprintf (stderr, "ERROR: Requested wavelength %.3lf is outside of range of color-matching functions.\n", lambda);
    fprintf (stderr, "       Available range is from %.3lf to %.3lf.\n", cmxyz[0][0], cmxyz[nlines-1][0]);
    exit (EXIT_FAILURE);
  }

  // Ask whether to apply the sRGB transfer function.
  fprintf (stdout, "\nApply the sRGB transfer function (required for standard sRGB image display)? ");
  memset (temp, 0, MAX_STRINGLEN * sizeof (char));
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
  {
    double xyz_sum = xyzbar[0] + xyzbar[1] + xyzbar[2];
    if (!isfinite (xyz_sum) || (fabs (xyz_sum) <= DBL_MIN)) {
      fprintf (stderr, "ERROR: CMF tristimulus sum is zero or non-finite at %.6lf nm.\n", lambda);
      exit (EXIT_FAILURE);
    }
    xyz[0] = xyzbar[0] / xyz_sum;
    xyz[1] = xyzbar[1] / xyz_sum;
    xyz[2] = xyzbar[2] / xyz_sum;
  }

  fprintf (stdout, "\nSpectral chromaticity coordinates xyz: %.12lf %.12lf %.12lf\n", xyz[0], xyz[1], xyz[2]);

  // Calculate linear srgb values from xyz coordinates.
  // rgb = npminv * XYZ
  memset (rgb_double, 0.0, 3 * sizeof (double));
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      rgb_double[i] += npminv[i][j] * xyzbar[j];
    }
    if (!isfinite (rgb_double[i])) {
      fprintf (stderr, "ERROR: Non-finite RGB value at wavelength %.6lf nm.\n", lambda);
      exit (EXIT_FAILURE);
    }
  }
  fprintf (stdout, "\nLinear RGB coordinates (no transfer function, no scaling, no clipping): %.12lf %.12lf %.12lf\n",
           rgb_double[0], rgb_double[1], rgb_double[2]);
  fprintf (stdout, "\nLinear RGB coordinates scaled by 255 (no clipping): %.12lf %.12lf %.12lf\n",
           rgb_double[0] * 255.0, rgb_double[1] * 255.0, rgb_double[2] * 255.0);

  // Apply the sRGB transfer function to linear RGB, if requested.
  if (apply_gamma) {
    for (i=0; i<3; i++) {
      if (rgb_double[i] > 0.0031308) {
        rgb_double[i] = (1.055 * pow (rgb_double[i], (1.0 / 2.4))) - 0.055;
      } else {
        rgb_double[i] *= 12.92;
      }
      if (!isfinite (rgb_double[i])) {
        fprintf (stderr, "ERROR: Non-finite transfer-function result at %.6lf nm.\n", lambda);
        exit (EXIT_FAILURE);
      }
    }
  }
  if (apply_gamma) {
    fprintf (stdout, "\nRGB coordinates (with sRGB transfer function, no scaling, no clipping): %.12lf %.12lf %.12lf\n",
             rgb_double[0], rgb_double[1], rgb_double[2]);
  } else {
    fprintf (stdout, "\nLinear RGB coordinates retained without a transfer function: %.12lf %.12lf %.12lf\n",
             rgb_double[0], rgb_double[1], rgb_double[2]);
  }

  // Clip to [0,1] and quantize to 8-bit RGB without converting an
  // out-of-range floating-point value to int.
  for (i=0; i<3; i++) {
    if (rgb_double[i] <= 0.0) {
      rgb[i] = 0;
    } else if (rgb_double[i] >= 1.0) {
      rgb[i] = 255;
    } else {
      rgb[i] = (int) lround (rgb_double[i] * 255.0);
    }
  }

  if (apply_gamma) {
    fprintf (stdout, "\n8-bit clipped RGB coordinates (sRGB transfer function applied): (%i, %i, %i)\n",
             rgb[0], rgb[1], rgb[2]);
  } else {
    fprintf (stdout, "\n8-bit clipped linear RGB coordinates (no transfer function): (%i, %i, %i)\n",
             rgb[0], rgb[1], rgb[2]);
  }
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
    free (npm[i]);
    free (npminv[i]);
  }
  free (p);
  free (pinv);
  free (w);
  free (coeff);
  free (npm);
  free (npminv);
  free (rgb_double);

  return (EXIT_SUCCESS);
}

// Obtain a text string from standard input. It can include spaces.
int
inputtext (char *text) {

  size_t len;

  if (fgets (text, MAX_STRINGLEN, stdin) == NULL) {
    fprintf (stderr, "Unable to read text from standard input.\n");
    exit (EXIT_FAILURE);
  }

  len = strlen (text);
  if ((len > 0u) && (text[len - 1u] == '\n')) {
    text[len - 1u] = '\0';
  } else if (len == MAX_STRINGLEN - 1u) {
    int ch;
    while (((ch = getchar ()) != '\n') && (ch != EOF)) {
      // Discard the remainder of an overlong input line.
    }
    fprintf (stderr, "Input text is too long; maximum is %d characters.\n", MAX_STRINGLEN - 2);
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
  if ((*endptr != '\0') || (parsed < INT_MIN) || (parsed > INT_MAX)) return (EXIT_FAILURE);
  *value = (int) parsed;
  return (EXIT_SUCCESS);
}

// Convert a complete input string to a finite double, allowing surrounding whitespace only.
int
parse_double_string (const char *text, double *value) {

  char *endptr;
  double parsed;

  if ((text == NULL) || (value == NULL)) return (EXIT_FAILURE);
  errno = 0;
  parsed = strtod (text, &endptr);
  if ((errno == ERANGE) || (endptr == text) || !isfinite (parsed)) return (EXIT_FAILURE);
  while (isspace ((unsigned char) *endptr)) endptr++;
  if (*endptr != '\0') return (EXIT_FAILURE);
  *value = parsed;
  return (EXIT_SUCCESS);
}

// Parse one numeric CMF record: wavelength, xbar, ybar, zbar.
// Blank lines, comments, and textual headers are ignored (return 0).
// Once a line begins with a numeric field, all four finite values are required.
int
parse_cmf_record (const char *line, double *values) {

  const char *p;
  char *endptr;
  int i;

  if ((line == NULL) || (values == NULL)) return (-1);
  p = line;
  while (isspace ((unsigned char) *p)) p++;
  if ((*p == '\0') || (*p == '#') || (*p == ';')) return (0);

  for (i=0; i<4; i++) {
    errno = 0;
    values[i] = strtod (p, &endptr);
    if (endptr == p) return ((i == 0) ? 0 : -1);
    if ((errno == ERANGE) || !isfinite (values[i])) return (-1);
    p = endptr;
    while (isspace ((unsigned char) *p)) p++;
  }

  if (*p != '\0') return (-1);
  return (1);
}

// Read one physical line from a CSV/text file, converting commas to spaces.
// Returns 0 for a line, -1 for EOF, -2 for an overlong line, and -3 for I/O error.
int
readline (FILE *fi, char *line, int limit) {

  size_t i, len;
  int ch;

  if ((fi == NULL) || (line == NULL) || (limit < 2)) return (-3);

  if (fgets (line, limit, fi) == NULL) {
    return (feof (fi) ? -1 : -3);
  }

  len = strlen (line);
  if ((len > 0u) && (line[len - 1u] == '\n')) {
    line[--len] = '\0';
  } else {
    // fgets() filled the buffer or reached EOF without a newline. Peek one byte
    // to distinguish a valid final/exact-fit line from a genuinely overlong one.
    ch = fgetc (fi);
    if (ch != EOF) {
      if (ch != '\n') {
        while (((ch = fgetc (fi)) != '\n') && (ch != EOF)) {
        }
        if (ferror (fi)) return (-3);
        return (-2);
      }
    } else if (ferror (fi)) {
      return (-3);
    }
  }

  if ((len > 0u) && (line[len - 1u] == '\r')) {
    line[--len] = '\0';
  }
  for (i=0u; i<len; i++) {
    if (line[i] == ',') line[i] = ' ';
  }

  return (0);
}

int
count_cmf_rows (const char *filename) {

  int count, parsed, status;
  double values[4];
  char line[MAX_STRINGLEN] = {0};
  FILE *fi;

  if (filename == NULL) return (-1);
  fi = fopen (filename, "r");
  if (fi == NULL) {
    fprintf (stderr, "ERROR: Unable to open input file %s.\n", filename);
    return (-1);
  }

  count = 0;
  while ((status = readline (fi, line, MAX_STRINGLEN)) != -1) {
    if (status == -2) {
      fprintf (stderr, "ERROR: Line in color-matching file exceeds %d characters.\n", MAX_STRINGLEN - 1);
      fclose (fi);
      return (-1);
    }
    if (status == -3) {
      fprintf (stderr, "ERROR: I/O error while reading color-matching file.\n");
      fclose (fi);
      return (-1);
    }
    parsed = parse_cmf_record (line, values);
    if (parsed < 0) {
      fprintf (stderr, "ERROR: Malformed numeric CMF row: %s\n", line);
      fclose (fi);
      return (-1);
    }
    if (parsed > 0) {
      if (count == INT_MAX) {
        fclose (fi);
        return (-1);
      }
      count++;
    }
  }

  if (ferror (fi)) {
    fclose (fi);
    return (-1);
  }
  fclose (fi);
  return (count);
}

int
load_cmf (int nlines, const char *filename, double **cmxyz) {

  int i, parsed, status;
  double values[4];
  char line[MAX_STRINGLEN] = {0};
  FILE *fi;

  if ((nlines < 2) || (filename == NULL) || (cmxyz == NULL)) return (EXIT_FAILURE);
  fi = fopen (filename, "r");
  if (fi == NULL) return (EXIT_FAILURE);

  i = 0;
  while ((status = readline (fi, line, MAX_STRINGLEN)) != -1) {
    if (status == -2) {
      fprintf (stderr, "ERROR: Line in color-matching file exceeds %d characters.\n", MAX_STRINGLEN - 1);
      fclose (fi);
      return (EXIT_FAILURE);
    }
    if (status == -3) {
      fprintf (stderr, "ERROR: I/O error while reading color-matching file.\n");
      fclose (fi);
      return (EXIT_FAILURE);
    }
    parsed = parse_cmf_record (line, values);
    if (parsed < 0) {
      fprintf (stderr, "ERROR: Malformed numeric CMF row: %s\n", line);
      fclose (fi);
      return (EXIT_FAILURE);
    }
    if (parsed == 0) continue;
    if (i >= nlines) {
      fclose (fi);
      return (EXIT_FAILURE);
    }

    cmxyz[i][0] = values[0];
    cmxyz[i][1] = values[1];
    cmxyz[i][2] = values[2];
    cmxyz[i][3] = values[3];

    if ((i > 0) && !(cmxyz[i][0] > cmxyz[i - 1][0])) {
      fprintf (stderr, "ERROR: CMF wavelengths must be strictly increasing; %.12g follows %.12g.\n",
               cmxyz[i][0], cmxyz[i - 1][0]);
      fclose (fi);
      return (EXIT_FAILURE);
    }
    i++;
  }

  if (ferror (fi)) {
    fclose (fi);
    return (EXIT_FAILURE);
  }
  fclose (fi);
  if (i != nlines) return (EXIT_FAILURE);
  return (EXIT_SUCCESS);
}

// Gauss-Jordan in-place inversion of n*n matrix with partial pivoting.
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

  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      augmented[i][j] = matrix[i][j];
      augmented[i][j + n] = (i == j) ? 1.0 : 0.0;
    }
  }

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
    for (j=0; j<(2 * n); j++) augmented[i][j] /= pivot;
    for (j=0; j<n; j++) {
      if (j == i) continue;
      factor = augmented[j][i];
      for (k=0; k<(2 * n); k++) augmented[j][k] -= factor * augmented[i][k];
    }
  }

  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) matrix[i][j] = augmented[i][j + n];
  }
  for (i=0; i<n; i++) free (augmented[i]);
  free (augmented);
  return (EXIT_SUCCESS);
}

// Return XYZ color-matching-function coordinates by linear interpolation.
int
cmf (double lambda, int nlines, double **cmxyz, double *xyzbar) {

  int i, lb, ub;
  double frac;

  if ((nlines < 2) || (cmxyz == NULL) || (xyzbar == NULL) || !isfinite (lambda)) {
    return (EXIT_FAILURE);
  }
  if ((lambda < cmxyz[0][0]) || (lambda > cmxyz[nlines - 1][0])) return (EXIT_FAILURE);

  i = 0;
  while ((i < nlines) && (cmxyz[i][0] < lambda)) i++;
  if (i == 0) {
    lb = 0;
    ub = 1;
  } else if (i >= nlines) {
    return (EXIT_FAILURE);
  } else {
    lb = i - 1;
    ub = i;
  }

  if (!(cmxyz[ub][0] > cmxyz[lb][0])) return (EXIT_FAILURE);
  frac = (lambda - cmxyz[lb][0]) / (cmxyz[ub][0] - cmxyz[lb][0]);
  xyzbar[0] = (frac * (cmxyz[ub][1] - cmxyz[lb][1])) + cmxyz[lb][1];
  xyzbar[1] = (frac * (cmxyz[ub][2] - cmxyz[lb][2])) + cmxyz[lb][2];
  xyzbar[2] = (frac * (cmxyz[ub][3] - cmxyz[lb][3])) + cmxyz[lb][3];
  return (EXIT_SUCCESS);
}

// Choose the white used to derive the RGB conversion matrix.
// Choice 0 is the nominal D65 chromaticity that defines standard sRGB.
int
illum_white (double *white_xyz) {

  int choice, i;
  char *temp;
  const WhitePoint *white;

  if (white_xyz == NULL) return (-1);
  temp = allocate_strmem (MAX_STRINGLEN);

  fprintf (stdout, "\nChoose display/reference white point:\n");
  fprintf (stdout, "  0 - RECOMMENDED for standard sRGB: D65 space definition (x=0.31270000, y=0.32900000)\n");
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
    white_xyz[0] = 0.3127;
    white_xyz[1] = 0.3290;
    fprintf (stdout, "Using nominal D65 from the sRGB/BT.709 space definition.\n");
  } else if ((choice >= 1) && (choice <= N_WHITE_POINTS)) {
    white = &white_points[choice - 1];
    white_xyz[0] = white->x;
    white_xyz[1] = white->y;
    fprintf (stdout, "Using %s.\nData basis: %s\n", white->name, white->source_note);
    fprintf (stdout, "NOTE: this selected-white RGB system is not standard sRGB.\n");
  } else {
    fprintf (stderr, "Invalid selection.\n");
    free (temp);
    return (-1);
  }

  white_xyz[2] = 1.0 - white_xyz[0] - white_xyz[1];
  free (temp);
  return (0);
}

// Allocate memory for an array of ints.
int *
allocate_intmem (size_t len) {

  int *tmp;

  if (len == 0u) {
    fprintf (stderr, "ERROR: Cannot allocate zero-length int array.\n");
    exit (EXIT_FAILURE);
  }
  tmp = calloc (len, sizeof (*tmp));
  if (tmp == NULL) {
    fprintf (stderr, "ERROR: Cannot allocate memory for int array.\n");
    exit (EXIT_FAILURE);
  }
  return (tmp);
}

// Allocate memory for an array of chars.
char *
allocate_strmem (size_t len) {

  char *tmp;

  if (len == 0u) {
    fprintf (stderr, "ERROR: Cannot allocate zero-length string.\n");
    exit (EXIT_FAILURE);
  }
  tmp = calloc (len, sizeof (*tmp));
  if (tmp == NULL) {
    fprintf (stderr, "ERROR: Cannot allocate memory for string.\n");
    exit (EXIT_FAILURE);
  }
  return (tmp);
}

// Allocate memory for an array of doubles.
double *
allocate_doublemem (size_t len) {

  double *tmp;

  if (len == 0u) {
    fprintf (stderr, "ERROR: Cannot allocate zero-length double array.\n");
    exit (EXIT_FAILURE);
  }
  tmp = calloc (len, sizeof (*tmp));
  if (tmp == NULL) {
    fprintf (stderr, "ERROR: Cannot allocate memory for double array.\n");
    exit (EXIT_FAILURE);
  }
  return (tmp);
}

// Allocate memory for an array of pointers to arrays of doubles.
double **
allocate_doublememp (size_t len) {

  double **tmp;

  if (len == 0u) {
    fprintf (stderr, "ERROR: Cannot allocate zero-length pointer array.\n");
    exit (EXIT_FAILURE);
  }
  tmp = calloc (len, sizeof (*tmp));
  if (tmp == NULL) {
    fprintf (stderr, "ERROR: Cannot allocate memory for pointer array.\n");
    exit (EXIT_FAILURE);
  }
  return (tmp);
}

