# Color-matching-function data

The CSV files in this directory contain CIE XYZ color-matching-function data obtained from the [Colour & Vision Research Laboratory (CVRL)](http://www.cvrl.org/).

Each row has four comma-separated fields and no header:

```text
wavelength_nm,xbar,ybar,zbar
```

| File | Data set | Wavelength range | Sampling |
|---|---|---:|---:|
| `CIE_xyz_1931_2deg.csv` | CIE 1931 2° XYZ color-matching functions | 360–830 nm | 1 nm |
| `CIE_xyz_1931_2deg_judd1951.csv` | CIE 1931 2° XYZ CMFs with Judd (1951) modification | 370–770 nm | 10 nm |
| `CIE_xyz_1931_2deg_judd1951_vos1978.csv` | CIE 1931 2° XYZ CMFs with Judd (1951) and Vos (1978) modifications | 380–825 nm | 5 nm |
| `CIE_xyz_1964_10deg.csv` | CIE 1964 10° XYZ color-matching functions | 360–830 nm | 1 nm |
| `CIE_xyz_2006_2deg_lms_cones.csv` | 2006 2° XYZ CMFs transformed from the CIE 2006 2° LMS cone fundamentals | 390–830 nm | 0.1 nm |
| `CIE_xyz_2006_10deg_lms_cones.csv` | 2006 10° XYZ CMFs transformed from the CIE 2006 10° LMS cone fundamentals | 390–830 nm | 0.1 nm |

The C programs use these filenames for their built-in CMF choices. A user-supplied CSV may also be selected by filename in programs that offer that option.

The source-code license does not automatically establish redistribution terms for third-party data. Consult CVRL and the underlying CIE data sources for any applicable terms when redistributing the datasets.
