## EnviroPolate

This repo contains a series of scripts for interpolating data from a sparse grid onto a dense grid. Included are scripts for fitting several different GAM 2d-spline models, and performing cross-validation to choose the best model for interpolation. An example is provided using the follwing environmental variables: PPFD (µmol m-2 s-1), Air temperature (°C), Air humidity (%), Air VPD (kPa). The example data consists of these variables, collected every 15 minutes for one day in the PhenoArcch phenotyping platform. The sparse grid consists of data collected from 10 sensors distributed in the platform, and the dense grid is a physical metric coordinate system to interpolate PPFD down to the individual pot-level. 

## Example

These scripts should be run in sequential order:

* 1_env_interpolation_model_comparison.R
* 2_env_interpolation.R
* 3_env_interpolation_plots.R

The 