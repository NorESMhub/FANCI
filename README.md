# FANCI
The sectional aerosol model for NorESM (Flexible Aerosol Number and Composition Interactions)

## config
Contains configuration files for FANCI that are used/executed during the build process.
More information in the [config documentation](config/CONFIG.md)

## src_cam
Contains files from OsloAero and CAM that are currently necessary to make the FANCI code run. These files will be phased out in the future.

## src
The FANCI source code.

## ../pp_fanci
* chemistry.F90
* chem_mech.in

## Getting this thing running
A quick guide is in the [Setup documentation](SETUP.md) - wiki coming soon.

## How to contribute
PR's should be made to the NorESMhub repositories for CAM and FANCI respectively.

NorESMhub/CAM:
* PR's from contributors should go into the branch sectional_develop
* Contributors are responsible to keep their clones and forks up to date with the newest CAM version from NorESMhub/CAM/sectional_develop
* @Trudigard is responsible to keep NorESMhub/CAM/sectional_develop up to date with NorESMhub/CAM/noresm_develop

NorESMhub/FANCI:
* PR's from contributors should go into sectional_develop
* Contributors are responsible to keep their clones and forks up to date with NorESMhub/FANCI/sectional_develop
