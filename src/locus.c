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

typedef struct {
    uint8_t blue;   // Blue
    uint8_t green;  // Green
    uint8_t red;    // Red
} PIXEL;

// Function prototypes
int inputtext (char *);
int readline (FILE*, char*, int);
int choose_cmf (int *, char *, double *);
int load_cmf (int, char *, double **);
int cmf (double, int, double **, double *);
int xy2uv (double, double, int *, int *, int, int);
int plot (int, int, int *, uint8_t *, int, int);
int draw_num (int, int, char *, int *, unsigned char *, int, int);
int draw_line (int *, int *, int *, uint8_t *, int, int);
int within_polygon (double, double, double **, int);
int bmp (char *, uint8_t *, int, int);
void write_u16_le (FILE *, uint16_t);
void write_u32_le (FILE *, uint32_t);
void write_s32_le (FILE *, int32_t);
int rgb_primaries (double **);
int illum_white (double *);
int *allocate_intmem (int);
int **allocate_intmemp (int);
char *allocate_strmem (int);
double *allocate_doublemem (int);
double **allocate_doublememp (int);
uint8_t *allocate_ustrmem (int);

// Set some symbolic constants.
#define MAXLEN 256  // Maximum number of characters per line

int
main (int argc, char **argv) {

  int i, j, nlines, narray, *rgb, width, height, count, u, v, *left, *right;
  int **uv, u_centroid, v_centroid, min, max, add_axes, mark_white, du, dv, uborder, vborder, mark_rgb;
  int *point1_uv, *point2_uv, fill_srgb;
  double INTERVAL, **cmxyz, lambda, **xyzbar, **xyz;
  double *white_xyz, val, **p, **polygon, *xyzvector, *rgb_double;
  uint8_t *buffer, *buffer2;
  char *filename, *temp, *endptr;

  // Wavelength step size (interpolate if necessary)
  // INTERVAL is changed below to 10 nm for Judd 1951 data; even still, it's poorly behaved.
  INTERVAL = 0.1;

  // Side and bottom borders to add if axes are to be drawn.
  uborder = 44;  // Side border widths
  vborder = 40;  // Top & bottom border widths

  // Allocate memory for various arrays.
  filename = allocate_strmem (MAXLEN);
  temp = allocate_strmem (MAXLEN);
  rgb = allocate_intmem (3);
  white_xyz = allocate_doublemem (3);
  xyzvector = allocate_doublemem (3);
  rgb_double = allocate_doublemem (3);
  left = allocate_intmem (2);
  right = allocate_intmem (2);
  p = allocate_doublememp (3);
  polygon = allocate_doublememp (3);
  for (i=0; i<3; i++) {
    p[i] = allocate_doublemem (3);
    polygon[i] = allocate_doublemem (2);
  }
  point1_uv = allocate_intmem (2);
  point2_uv = allocate_intmem (2);

  // Ask for bitmap dimension.
  fprintf (stdout, "\nBitmap dimension (square) (px)? ");
  memset (temp, 0, MAXLEN * sizeof (char));
  inputtext (temp);
  errno = 0;
  width = (int) strtol (temp, &endptr, 10);
  if ((errno == ERANGE) || (errno == EINVAL) || (endptr == temp)) {
    fprintf (stderr, "ERROR: Cannot make integer of: %s\n", temp);
    exit (EXIT_FAILURE);
  }
  height = width;  // Same as width so locus isn't distorted.

  // Allocate memory for various arrays.
  buffer = allocate_ustrmem (width * height * 3);  // For bitmap data without axes.
  buffer2 = allocate_ustrmem ((width + (2 * uborder)) * (height + (2 * vborder)) * 3);  // For bitmap data with axes.

  // Choose color-matching function (CMF).
  for (;;) {
    if (choose_cmf (&nlines, filename, &INTERVAL) > -1) break;
  }

  // Allocate memory for various arrays.
  cmxyz = allocate_doublememp (nlines);
  for (i=0; i<(nlines); i++) {
    cmxyz[i] = allocate_doublemem (4);  // lambda, xbar, ybar, zbar
  }

  // Load color-matching function (CMF).
  for (;;) {
    if (load_cmf (nlines, filename, cmxyz) > -1) break;
  }

  // Ask whether to plot outline of rgb gamut on bitmap.
  fprintf (stdout, "\nPlot outline of RGB gamut on bitmap (y/n)? ");
  memset (temp, 0, MAXLEN * sizeof (char));
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
  memset (temp, 0, MAXLEN * sizeof (char));
  inputtext (temp);
  if ((temp[0] == 'y') || (temp[0] == 'Y')) {
    fill_srgb = 1;
  } else {
    fill_srgb = 0;
  }

  // Ask whether to mark white point on bitmap.
  fprintf (stdout, "\nMark white point on bitmap (y/n)? ");
  memset (temp, 0, MAXLEN * sizeof (char));
  inputtext (temp);
  if ((temp[0] == 'y') || (temp[0] == 'Y')) {
    mark_white = 1;
    for (;;) {
      if (illum_white (white_xyz) > -1) break;  // Choose white point.
    }
  } else {
    mark_white = 0;
  }

  // Allocate memory for various arrays.
  narray = (int) (((cmxyz[nlines-1][0] - cmxyz[0][0]) / INTERVAL) + 1.5);
  xyzbar = allocate_doublememp (narray);
  xyz = allocate_doublememp (narray);
  uv = allocate_intmemp (narray);
  for (i=0; i<narray; i++) {
    xyz[i] = allocate_doublemem (3);
    xyzbar[i] = allocate_doublemem (3);
    uv[i] = allocate_intmem (2);
  }

  // Ask whether to include axes on bitmap.
  fprintf (stdout, "\nInclude axes on bitmap? (borders will be added to bitmap size) (y/n)? ");
  memset (temp, 0, MAXLEN * sizeof (char));
  inputtext (temp);
  if ((temp[0] == 'y') || (temp[0] == 'Y')) {
    add_axes = 1;
  } else {
    add_axes = 0;
  }

  // Loop through full range of wavelengths in nm steps defined by INTERVAL.
  count = 0;
  lambda = cmxyz[0][0];  // Start at first wavelength in CMF table.
  while (lambda <= cmxyz[nlines-1][0]) {

    // Retrieve wavelength and color-matching coodinates xbar, ybar, zbar for requested wavelength.
    // Note: For single wavelength, the CMF X,Y,Z coordinates are the linear scene tristimulus values. i.e., no need to integrate
    // xbar = xyzbar[0], ybar = xyzbar[1], zbar = xyzbar[2]
    if (cmf (lambda, nlines, cmxyz, xyzbar[count]) == EXIT_FAILURE) {
      fprintf (stderr, "ERROR: lambda %.18lf is outside range of CMF data: %.18lf to %.18lf.\n", lambda, cmxyz[0][0], cmxyz[nlines-1][0]);
      exit (EXIT_FAILURE);
    }
//  fprintf (stdout, "XYZ CMF coordinates: %.12lf %.12lf %.12lf\n", xyzbar[count][0], xyzbar[count][1], xyzbar[count][2]);

    // Compute spectral chromaticity coordinates xyz.
    xyz[count][0] = xyzbar[count][0] / (xyzbar[count][0] + xyzbar[count][1] + xyzbar[count][2]);
    xyz[count][1] = xyzbar[count][1] / (xyzbar[count][0] + xyzbar[count][1] + xyzbar[count][2]);
    xyz[count][2] = xyzbar[count][2] / (xyzbar[count][0] + xyzbar[count][1] + xyzbar[count][2]);
//  fprintf (stdout, "Spectral chromaticity coordinates xyz: %.12lf %.12lf %.12lf\n", xyz[count][0], xyz[count][1], xyz[count][2]);

    // Convert xy chromaticity coordinates to uv pixel coordinates.
    xy2uv (xyz[count][0], xyz[count][1], &u, &v, width, height);
    uv[count][0] = u;
    uv[count][1] = v;

    count++;
    lambda = cmxyz[0][0] + (((double) count) * INTERVAL);
  }
//fprintf (stdout, "count: %i   narray: %i\n", count, narray);

  // Create a bitmap of the chromaticity locus.
  // Plot locus points as white pixels.
  rgb[0] = 255;
  rgb[1] = 255;
  rgb[2] = 255;
  for (i=0; i<narray; i++) {
    plot (uv[i][0], uv[i][1], rgb, buffer, width, height);
  }

  // Calculate coordinates of the centroid of polygon defined by pixel coordinates of chromaticities.
  u_centroid = 0;
  v_centroid = 0;
  for (i=0; i<narray; i++) {
    u_centroid += uv[i][0];
    v_centroid += uv[i][1];
  }
  u_centroid = (int) ((double) u_centroid / (double) narray);
  v_centroid = (int) ((double) v_centroid / (double) narray);
//printf ("Centroid (x,y): %i %i\n", u_centroid, v_centroid);

  // Find lowest point to right of centroid. i.e., rightmost end of locus
  min = height;
  for (i=0; i<narray; i++) {
    if (uv[i][0] > u_centroid) {
      if (uv[i][1] < min) {
        min = uv[i][1];
        right[0] = uv[i][0];
        right[1] = uv[i][1];
      }
    }
  }
//  printf ("rightmost: %i %i\n", right[0], right[1]);

  // Find rightmost point to the left and below centroid. i.e., other end of locus
  max = 0;
  for (i=0; i<narray; i++) {
    if ((uv[i][1] < v_centroid) && (uv[i][0] < u_centroid)) {
      if (uv[i][0] > max) {
        max = uv[i][0];
        left[0] = uv[i][0];
        left[1] = uv[i][1];
      }
    }
  }
//  printf ("left: %i %i\n", left[0], left[1]);

  // Plot purple line.
  draw_line (left, right, rgb, buffer, width, height);

  // Fill in sRGB color gamut if requested.
  if (fill_srgb) {

    // Define a polygon by sRGB (BT.709) color primaries as vertices.
    polygon[0][0] = 0.640 * (double) width; polygon[0][1] = 0.330 * (double) height;  // Red
    polygon[1][0] = 0.300 * (double) width; polygon[1][1] = 0.600 * (double) height;  // Green
    polygon[2][0] = 0.150 * (double) width; polygon[2][1] = 0.060 * (double) height;  // Blue
   
    // Loop through all pixels in bitmap.
    for (u=0; u<width; u++) {
      for (v=0; v<height; v++) {

        // If pixel is within the sRGB color gamut, we plot it.
        if (within_polygon ((double) u, (double) v, polygon, 3)) {

          // Calculate chromaticity coordinates for current pixel.
          xyzvector[0] = (double) u / (double) width;
          xyzvector[1] = (double) v / (double) height;
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
    i = 0;  // Index of buffer
    for (v=0; v<height; v++) {
      for (u=0; u<width; u++) {
        buffer2[((u + uborder) * 3) + ((v + vborder) * (width + (2 * uborder)) * 3)] = buffer[i];
        buffer2[((u + uborder) * 3) + ((v + vborder) * (width + (2 * uborder)) * 3) + 1] = buffer[i + 1];
        buffer2[((u + uborder) * 3) + ((v + vborder) * (width + (2 * uborder)) * 3) + 2] = buffer[i + 2];
        i += 3;
      }
    }

    // Plot horizontal axis line.
    du = (int) (width / 10);
    point1_uv[0] = uborder; point1_uv[1] = vborder;
    point2_uv[0] = uborder + (10 * du); point2_uv[1] = vborder;
    draw_line (point1_uv, point2_uv, rgb, buffer2, width + (2 * uborder), height + (2 * vborder));

    // Plot horizontal axis tick marks.
    u = 0;
    for (i=0; i<=10; i++) {
      point1_uv[0] = u + uborder; point1_uv[1] = vborder - 10;
      point2_uv[0] = u + uborder; point2_uv[1] = vborder;
      draw_line (point1_uv, point2_uv, rgb, buffer2, width + (2 * uborder), height + (2 * vborder));
      u += du;
    }

    // Plot numerical labels for horizontal axis.
    val = 0.0;
    u = uborder - 14;  // Starting horizontal position (px)
    v = vborder - 24;  // Vertical position (px)
    for (i=0; i<=10; i++) {
      memset (temp, 0, MAXLEN * sizeof (char));
      sprintf (temp, "%0.1lf", val);
      draw_num (u, v, temp, rgb, buffer2, width + (2 * uborder), height + (2 * vborder));
      val += 0.1;
      u += du;
    }

    // Plot vertical axis line.
    dv = (int) (height / 10);
    point1_uv[0] = uborder; point1_uv[1] = vborder;
    point2_uv[0] = uborder; point2_uv[1] = vborder + (10 * dv);
    draw_line (point1_uv, point2_uv, rgb, buffer2, width + (2 * uborder), height + (2 * vborder));

    // Plot vertical axis tick marks.
    v = 0;
    for (i=0; i<=10; i++) {
      point1_uv[0] = uborder - 10; point1_uv[1] = v + vborder;
      point2_uv[0] = uborder; point2_uv[1] = v + vborder;
      draw_line (point1_uv, point2_uv, rgb, buffer2, width + (2 * uborder), height + (2 * vborder));
      v += dv;
    }

    // Plot numerical labels for vertical axis.
    val = 0.0;
    u = 0;  // Horizontal position (px)
    v = vborder - 4;  // Starting vertical position (px)
    for (i=0; i<=10; i++) {
      memset (temp, 0, MAXLEN * sizeof (char));
      sprintf (temp, "%0.1lf", val);
      draw_num (u, v, temp, rgb, buffer2, width + (2 * uborder), height + (2 * vborder));
      val += 0.1;
      v += dv;
    }

    // Create bitmap output file.
    bmp ("out.bmp", buffer2, width + (2 * uborder), height + (2 * vborder));

  // No axes to plot.
  } else {

    // Create bitmap output file.
    bmp ("out.bmp", buffer, width, height);

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
  free (left);
  free (right);
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

// Choose color-matching function (CMF).
// Returns: -1 if invalid selection, 0 if valid selection
int
choose_cmf (int *nlines, char *filename, double *INTERVAL) {

  int choice;
  char *temp, *endptr;
  FILE *fi;

  // Allocate memory for various arrays.
  temp = allocate_strmem (MAXLEN);

  fprintf (stdout, "\nChoose color-matching function (CMF):\n\n");
  fprintf (stdout, "  1 - 1964 10-deg XYZ CMFs (JIS Z 8701:1999)\n");
  fprintf (stdout, "  2 - 1931 2-deg XYZ CIE CMFs (CIE.15.2004)\n");
  fprintf (stdout, "  3 - 1931 2-deg XYZ CIE CMFs with Judd (1951) modifications [badly behaved data; won't interpolate here]\n");
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
      (*INTERVAL) = 10.0; 
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
      fprintf (stderr, "Invalid selection.\n");
      return (-1);
  }

  // Open color-matching functions csv file.
  fi = fopen (filename, "r");
  if (fi == NULL) {
    fprintf (stderr, "ERROR: Unable to open input file \"%s\".\n", filename);
    exit (EXIT_FAILURE);
  }

  // Count lines in input file.
  (*nlines) = 0;
  for (;;) {

    if (readline (fi, temp, MAXLEN) == -1) {
      break;  // Reached end of file.
    }

    if ((temp[0] >= '0') && (temp[0]<= '9')) (*nlines)++;

  }  // Next line of input file.
  fprintf (stdout, "%i lines in color-matching file: %s\n", (*nlines), filename);

  // Close input file.
  fclose (fi);

  // Free allocated memory.
  free (temp);

  return (0);  // Success
}

// Load color-matching function (CMF).
// Returns: -1 if invalid selection, 0 if valid selection
int
load_cmf (int nlines, char *filename, double **cmxyz) {

  int i;
  char *temp;
  FILE *fi;

  // Allocate memory for various arrays.
  temp = allocate_strmem (MAXLEN);

  // Open color-matching functions csv file.
  fi = fopen (filename, "r");
  if (fi == NULL) {
    fprintf (stderr, "ERROR: Unable to open input file \"%s\".\n", filename);
    exit (EXIT_FAILURE);
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

  // Free allocated memory.
  free (temp);

  return (0);  // Success
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

  // Chromaticities xy range from 0 to 1.
  // Pixel positions range from 0 to width, 0 to height.
  *u = (int) (x * width);
  *v = (int) (y * height);

  return (EXIT_SUCCESS);
}

// Plot a pixel in bitmap buffer as white.
// Ranges: 0 <= u <= width, 0 <= v <= height
int
plot (int u, int v, int *rgb, uint8_t *buffer, int width, int height) {

  int index;

  index = (u * 3) + (v * width * 3);

  buffer[index] = (uint8_t) rgb[2];  // B
  buffer[index + 1] = (uint8_t) rgb[1];  // G
  buffer[index + 2] = (uint8_t) rgb[0];  // R

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

  int i, j;
  double slope, intercept;

  // Non-vertical line
  if (p1[0] != p2[0]) {

    // Calculate slope and intercept of line.
    slope = ((double) p2[1] - (double) p1[1]) / ((double) p2[0] - (double) p1[0]);
    intercept = (double) p2[1] - (slope * (double) p2[0]);

    // Plot line.
    // Choose to iterate on horizontal or vertical px depending upon which has greatest difference.
    // This ensures densely populated pixel lines.
    if (fabs (p2[0] - p1[0]) > fabs (p2[1] - p1[1])) {
      if (p2[0] > p1[0]) {
        for (i=p1[0]; i<p2[0]; i++) {
          j = (int) ((((double) i * slope) + intercept) + 0.5);
          plot (i, j, rgb, buffer, width, height);
        }
      } else {
        for (i=p2[0]; i<p1[0]; i++) {
          j = (int) ((((double) i * slope) + intercept) + 0.5);
          plot (i, j, rgb, buffer, width, height);
        }
      }
    } else {
      if (p2[1] > p1[1]) {
        for (j=p1[1]; j<p2[1]; j++) {
          i = (int) (((double) j - intercept) / slope);
          plot (i, j, rgb, buffer, width, height);
        }
      } else {
        for (j=p2[1]; j<p1[1]; j++) {
          i = (int) (((double) j - intercept) / slope);
          plot (i, j, rgb, buffer, width, height);
        }
      }
    }

  // Vertical line
  } else {
    if (p2[1] > p1[1]) {
      for (j=p1[1]; j<p2[1]; j++) {
        plot (p1[0], j, rgb, buffer, width, height);
      }
    } else {
      for (j=p2[1]; j<p1[1]; j++) {
        plot (p1[0], j, rgb, buffer, width, height);
      }
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

  int row_size, image_size, x, y, c;
  uint8_t padding[3] = {0, 0, 0};  // Padding to make each row 4 bytes aligned
  FILE *fo;

  // Open output file.
  fo = fopen (filename, "r");
  if (fo != NULL) {
    fprintf (stderr, "Output file %s already exists.\n", filename);
    exit (EXIT_FAILURE);
  }
  fo = fopen (filename, "wb");
  if (fo == NULL) {
    printf ("Can't open output file %s.\n", filename);
    exit (EXIT_FAILURE);
  }

  // Calculate the padding required for each row.
  row_size = (width * 3 + 3) & (~3);  // Each row must be a multiple of 4 bytes
  image_size = row_size * height;

  // BMP file header
  write_u16_le (fo, 0x4d42);              // File type, should be "BM"
  write_u32_le (fo, 54 + image_size);     // Size of the file (bytes)
  write_u16_le (fo, 0);                   // Reserved (set to 0)
  write_u16_le (fo, 0);                   // Reserved (set to 0)
  write_u32_le (fo, 54);                  // Offset (bytes) to the start of the pixel data

  // BMP information header
  write_u32_le (fo, 40);                  // Size of this header (40 bytes)
  write_s32_le (fo, width);               // Width of the image (px)
  write_s32_le (fo, height);              // Height of the image (px)
  write_u16_le (fo, 1);                   // Number of color planes (always 1)
  write_u16_le (fo, 24);                  // Bits per pixel (24 for RGB)
  write_u32_le (fo, 0);                   // Compression method (0 for none)
  write_u32_le (fo, image_size);          // Size of the image data (bytes)
  write_s32_le (fo, 7874);                // Horizontal resolution (in pixels per meter) (200 DPI)
  write_s32_le (fo, 7874);                // Vertical resolution (in pixels per meter) (200 DPI)
  write_u32_le (fo, 0);                   // Number of colors used (0 for 2^24)
  write_u32_le (fo, 0);                   // Important colors (0 for all)

  // Loop through each row and write the pixels.
  c = 0;  // Index of buffer array
  for (y=0; y<height; y++) {
    for (x=0; x<width; x++) {

      // Write the blue, green, and red values (24-bit color).
      fputc (buffer[c], fo);  // B
      c++;
      fputc (buffer[c], fo);  // G
      c++;
      fputc (buffer[c], fo);  // R
      c++;
    }

    // Write padding, if necessary.
    fwrite (padding, 1, row_size - width * 3, fo);
  }

  // Close output file.
  fclose (fo);

  return (EXIT_SUCCESS);
}

// Write a unsigned little-endian 16-bit value to file.
void
write_u16_le (FILE *fo, uint16_t val) {

    fputc (val & 0xff, fo);
    fputc ((val >> 8) & 0xff, fo);

}

// Write a unsigned little-endian 32-bit value to file.
void
write_u32_le (FILE *fo, uint32_t val) {

    fputc (val & 0xff, fo);
    fputc ((val >> 8) & 0xff, fo);
    fputc ((val >> 16) & 0xff, fo);
    fputc ((val >> 24) & 0xff, fo);

}

// Write a signed little-endian 32-bit value to file.
void
write_s32_le (FILE *fo, int32_t val) {

    write_u32_le (fo, (uint32_t) val);  // Cast as unsigned in order to preserve bit pattern.

}

// Color primary coordinates for various RGB colorspaces
// Returns: -1 if invalid selection, 0 if valid selection

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

// Allocate memory for an array of pointers to arrays of ints.
int **
allocate_intmemp (int len) {

  void *tmp;

  if (len <= 0) {
    fprintf (stderr, "ERROR: Cannot allocate memory because len = %i in allocate_intmemp().\n", len);
    exit (EXIT_FAILURE);
  }

  tmp = (int **) malloc (len * sizeof (int *));
  if (tmp != NULL) {
    memset (tmp, 0, len * sizeof (int *));
    return (tmp);
  } else {
    fprintf (stderr, "ERROR: Cannot allocate memory for array in allocate_intmemp().\n");
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
