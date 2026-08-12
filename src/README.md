# Colorimetry programs

This directory contains four interactive C programs for exploring colorimetry from color-matching-function data.

| Program | Purpose | Output |
|---|---|---|
| [`locus.c`](locus.c) | Plot the chromaticity locus for a selected color-matching function, with optional RGB-gamut outline, sRGB-gamut fill, white point, and axes. | `out.bmp` |
| [`matrix.c`](matrix.c) | Derive RGB-to-XYZ and XYZ-to-RGB conversion matrices for a selected RGB color space and white point. | Standard output |
| [`vis-spec.c`](vis-spec.c) | Produce an sRGB bitmap representation of the visible spectrum for a selected color-matching function and white point. | `out.bmp` |
| [`wl.c`](wl.c) | Report CMF tristimulus values, spectral chromaticity coordinates, and corresponding sRGB values for a selected wavelength. | Standard output |

## Building

From the repository root:

```sh
make
```

This builds `locus`, `matrix`, `vis-spec`, and `wl` in the repository root. The three programs that use CMF tables expect the bundled files under `data/`, so run the programs from the repository root.

Individual programs can also be built directly, for example:

```sh
gcc -Wall src/locus.c -lm -o locus
```

## Running

Each program is interactive:

```sh
./locus
./matrix
./vis-spec
./wl
```

`locus` and `vis-spec` refuse to overwrite an existing `out.bmp`; rename or remove that file before generating another bitmap.

## Notes

The programs deliberately expose intermediate colorimetric quantities and choices rather than hiding them behind a library API. This reflects the project's goal of working from first principles.

The source code is licensed under the GNU General Public License, version 3 or later, as stated in each source file.
