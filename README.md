# ISIMIP non-climate input conversion
This repository provides code to convert various types of non-climate inputs
provided by the Inter-Sectoral Impact Model Intercomparison Project
(https://www.isimip.org/) into inputs that can be used with LPJmL.

The purpose is to document how source data as provided by ISIMIP during phase
3 a & b are converted into inputs that can be used with LPJmL.

It may also help other users to generate inputs for their versions of LPJmL.

See https://github.com/PIK-LPJmL/ISIMIP_climate_conversion for conversion
scripts for ISIMIP climate data.

## Requirements

- R
- R packages `ncdf4`, `raster`, `udunits2` and
  [`lpjmlkit`](https://github.com/PIK-LPJmL/lpjmlkit)

## Files

- `create_deposition_inputs_ISIMIP3a.R`, `create_deposition_inputs_ISIMIP3b.R`:
  R scripts to convert N deposition inputs
- `LICENSE`: copy of GNU AFFERO GENERAL PUBLIC LICENSE
- `README.md`: this file

## Notes

See the [README for lpjmlkit](https://github.com/PIK-LPJmL/lpjmlkit/blob/master/README.md)
for instructions how to install this R package because it is not available
through CRAN.
