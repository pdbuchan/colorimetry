# Colorimetry

C programs and CIE color-matching-function data for exploring chromaticity, RGB/XYZ conversion matrices, monochromatic wavelengths, and visible-spectrum rendering from first principles.

The illustrated project description is available on the [Colorimetry GitHub Pages website](https://pdbuchan.github.io/colorimetry/).

## Overview

After working on [PGS subtitle extraction](https://github.com/pdbuchan/subtitle-tools/tree/main/pgs), which involves conversion of BT.709 YCbCr to the sRGB color space, I continued investigating other aspects of colorimetry. This repository contains the C-language tools developed during that work.

The intent is to work from first principles to the greatest extent practical. The programs therefore work directly with CIE color-matching functions (CMFs). The bundled CMF data were obtained from the [Colour & Vision Research Laboratory (CVRL)](http://www.cvrl.org/).

For display on a typical computer monitor, sRGB uses the D65 white point and the same RGB chromaticity primaries as BT.709. The programs can apply the sRGB transfer function to linear RGB values where appropriate. Other RGB primaries and white points are also available in some programs for comparison and experimentation.

The original CIE 1931 standard observer was based on a 2° field of view. Because cone distribution across the retina is not uniform, the CIE later standardized a 10° observer in 1964 for applications involving a larger field of view. This repository includes both 2° and 10° CMF datasets, as well as modified and later datasets.

> **Important gamut limitation:** monochromatic stimuli on the spectral locus cannot, in general, be reproduced exactly in sRGB because their chromaticities lie outside the sRGB gamut triangle. An sRGB image of the visible spectrum is therefore only a display approximation. The same limitation applies when a chromaticity diagram fills regions outside the sRGB triangle with sRGB pixel values.

## Programs

| Program | Description |
|---|---|
| [`locus.c`](src/locus.c) | Produces a bitmap of the chromaticity locus for a selected CMF. Optional features include an RGB-gamut outline, sRGB-gamut fill, white-point marker, and axes. |
| [`matrix.c`](src/matrix.c) | Derives RGB-to-XYZ and XYZ-to-RGB conversion matrices for a selected RGB color space and white point. |
| [`wl.c`](src/wl.c) | Reports XYZ CMF values, spectral chromaticity coordinates, and corresponding sRGB values for a selected wavelength. |
| [`vis-spec.c`](src/vis-spec.c) | Produces an sRGB bitmap approximation of the visible spectrum using a selected CMF and white point. |

See [`src/README.md`](src/README.md) for program-specific notes.

## Color-matching-function data

| File | Data set | Range | Sampling |
|---|---|---:|---:|
| [`CIE_xyz_1931_2deg.csv`](data/CIE_xyz_1931_2deg.csv) | CIE 1931 2° XYZ CMFs | 360–830 nm | 1 nm |
| [`CIE_xyz_1931_2deg_judd1951.csv`](data/CIE_xyz_1931_2deg_judd1951.csv) | CIE 1931 2° XYZ CMFs with Judd (1951) modification | 370–770 nm | 10 nm |
| [`CIE_xyz_1931_2deg_judd1951_vos1978.csv`](data/CIE_xyz_1931_2deg_judd1951_vos1978.csv) | CIE 1931 2° XYZ CMFs with Judd (1951) and Vos (1978) modifications | 380–825 nm | 5 nm |
| [`CIE_xyz_1964_10deg.csv`](data/CIE_xyz_1964_10deg.csv) | CIE 1964 10° XYZ CMFs | 360–830 nm | 1 nm |
| [`CIE_xyz_2006_2deg_lms_cones.csv`](data/CIE_xyz_2006_2deg_lms_cones.csv) | 2006 2° XYZ CMFs transformed from CIE 2006 2° LMS cone fundamentals | 390–830 nm | 0.1 nm |
| [`CIE_xyz_2006_10deg_lms_cones.csv`](data/CIE_xyz_2006_10deg_lms_cones.csv) | 2006 10° XYZ CMFs transformed from CIE 2006 10° LMS cone fundamentals | 390–830 nm | 0.1 nm |

The files contain four comma-separated fields per row with no header: wavelength in nanometres, followed by x̄, ȳ, and z̄. See [`data/README.md`](data/README.md) for details.

## Building

A simple Makefile is included. From the repository root:

```sh
make
```

This builds four executables in the repository root:

```text
locus
matrix
vis-spec
wl
```

The programs are interactive. For example:

```sh
./wl
```

The CMF-reading programs use the datasets in `data/` for their built-in selections, so they should normally be run from the repository root. `locus` and `vis-spec` write `out.bmp` and intentionally refuse to overwrite an existing file of that name.

To remove generated executables and `out.bmp`:

```sh
make clean
```

## Repository layout

```text
colorimetry/
├── README.md
├── Makefile
├── LICENSE
├── .gitignore
├── data/
│   ├── README.md
│   └── CIE_xyz_*.csv
├── src/
│   ├── README.md
│   ├── locus.c
│   ├── matrix.c
│   ├── vis-spec.c
│   └── wl.c
└── docs/
    ├── README.md
    ├── index.html
    ├── .nojekyll
    └── assets/
        ├── colorimetry.css
        └── images/
            └── README.md
```

The `docs/` directory is ready to be selected as the GitHub Pages publishing source for the `main` branch.

## License and data provenance

The C source files are licensed under the GNU General Public License, version 3 or later; see [`LICENSE`](LICENSE).

The bundled CMF datasets were obtained from CVRL and originate from CIE and related colorimetric work. The software license does not by itself determine redistribution rights for third-party data; consult the original data sources for any applicable terms.
