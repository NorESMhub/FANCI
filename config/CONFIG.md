# FANCI configuration
User input and fanci_config.py. These configuration files are used to specify the bin structure, as well as the aerosol species and their properties. fanci_config.py is called by cime_config/buildnml during the build process when a FANCI compset is being used.
The general workflow is:

### Input
- **configuration_file.ini** : Meant for user edits to adjust bin configuration or aerosol species and their properties
- **chem_mech.in** : Lists all gas species and their chemical reactions but *without* aerosol tracers

### Processing
- **fanci_config.py** : Executed in buildnml, creates files read by the model and puts them in the CASEDIR

### Output
- **fanci_aerosol_properties_nl** : Namelist for the general aerosol properties, included in atm_in (CASEDIR)
- **fanci_aerosol_species_properties * nspecies** : Species specific namelists, one per aerosol species in atm_in (CASEDIR)
- **chem_mech.in** : Edited chem_mech.in including all aerosol tracers (CASEDIR)

## User input options *.ini
E.g. fanci.ini.
Specifies the bin configuration, aerosol species and their properties.

### Bin specs
- **nbin**: Number of bins
- **radius_1**: *Center* radius of the smallest bin [nm]
- **radius_N**: *Center* radius of the largest bin [nm]

### Range specs
- **ranges**: Boolean indicating whether size resolution for chemical species should be lower than total aerosol number concentration. If `ranges == True`: The bins (aerosol number concentration) and ranges (chemical composition) have different size resolutions. The composition in a number of adjacent bins (one range) is averaged. If `ranges == False`: The scheme works like a classical sectional scheme where each species has one tracer per bin.
- **range_bounds**: The *boundary* radii for the ranges [nm]. These range bounds get adjusted to the closest bin bound, once they are calculated in fanci_config.py

### Species
All species have their own section. fanci_config.py will assume that every section (denoted by the [] parentheses) after [RANGE SPECS] is a new species.
Each species needs all the properties specified below:

- **[NAME]**                      : The name for the species, this is only relevant in fanci_config.py
- **active** = \<True>/\<False>       : Is this species an active component or not. If set to true, the corresponding tracers are added to chem_mech.in and the fanci_aerosol_properties_nl namelist in atm_in
- **short_name** = \<XX>           : A short string identifier for each species, e.g. *DU* or *SO4*
- **long_name** = \<some string>   : A long name for the species, e.g. *Dust aerosol*
- **range_bounds** = \<lower>, \<upper> : The *boundary* radii for this species. These should be the same as two values in range_bounds, in order to avoid species ending up in unexpected bins.
- **composition** = \<X2Yy3>       : Chemical formula for the compound, this is used for the tracer specification in chem_mech.in
- **density** = 0.0             : Density of the aerosol species [kg/m3]
- **molecular_weight** = 0.0    : Molecular weight of the aerosol species in [kg/kmol]
- **mixed** = \<True>/\<False>     : When True, the aerosol is internally mixed with all other aerosol of this kind. False means only externally mixed.
- **kappa** = 0.0               : Hygroscopicity parameter kappa, current values from SALSA (Oslo_aero: hygroscopicity according to eq. 4 in Abdul-Razzak and Ghan (doi: 10.1029/1999JD901161/), calculated as `nr Ions * osmotic coefficient * soluble mass fraction / molecular weight`)


## fanci_config.py

### Features
- Reads and processes the configuration settings in the *.ini file described above
- Calculates bin specifications for the aerosol size distribution
- Generates species specifig settings for the aerosol tracers
- Outputs necessary files for the sectional aerosol model in NorESM (namelists and chem_mech)

The script contains two main functions: *fanci_config* and *add_FANCI_nl* that called by *buildnml*.
*fanci_config* is called by buildnml before the chemical pre-processor to set the aerosol configuration.
*add_FANCI_nl* is called by buildnml after the namelist file atm_in has been created and attaches the new sectional aerosol namelists to it.

### Input files
- **Aerosol Configuration file (*.ini)**: The aerosol configuration file described above.
- **Chemistry mechanism file path (chem_mech.in)**: The original chemistry mechanism file for NorESM, which will be updated with aerosol tracers. Usually a chem_mech.in from an existing pp_folder, or a -usr_mech_infile from the runscript.
- **Atmospheric namelist file (atm_in**): The original CAM namelist required for NorESM, which will be modified to include the sectional aerosol namelists generated in fanci_config

### Output files
- **Namelist file (atm_in) + edits**: fanci_config.py creates additional namelist entries that are inserted into the CAM atm_in namelist file via the *add_FANCI_nl* call in buildnml following alphabetical order for namelists.

     #### fanci_properties_nl
    - **fanci_nspecies_tot**: Number of species in the model. This corresponds to the number of sections in the *.ini file describing an aerosol species
    - **fanci_nbins**: Number of bins from the *.ini file
    - **fanci_nranges**: Number of chemical ranges
    - **fanci_nspecies**: Number of species in each range
    - **fanci_bin_bounds**: The boundary radii for each bin [nm]
    - **fanci_bin_centers**: The center radius for each bin [nm]
    - **fanci_range_bounds**: A list of the indices of smallest:largest bin in a range

     #### fanci_properties_aerosol_nl (one per species)
    - **fanci_aerosol_name**: Short name of the aerosol species (e.g. "DU")
    - **fanci_aerosol_range**: Indices of the smallest:largest range including the species
    - **fanci_aerosol_mixed**: True if internally mixed
    - **fanci_aerosol_density**: Density of the aerosol species ()
    - **fanci_aerosol_weight**: Molecular weight of the aerosol species ()

- **Updated chemistry mechanism file (my_chem_mech.in)**: The modified version of chem_mech.in, now including aerosol tracers in the form of e.g. *num_1, num_2, ..., DU_R3, DU_R4*
Modified chem_mech.in file

## Usage instructions
*buildnml* will automatically call fanci_config.py, if a sectional aerosol configuration file is defined either via user input or by using a FANCI compset.

### Running from command line
**Arguments**
- --aerconf: Path to the aerosol configuration file (*.ini)
- --chem_mech: Path to the chemistry mechanism file (chem_mech.in)
- --chem_mech_new: Path where the modified chemistry mechanism_file should be saved
- --atm_in: Path to the original namelist file
- --atm_in_new: Path where the new atm_in file should be saved

**Example**

`python fanci_config.py --aerconf config.ini --chem_mech chem_mech.in --chem_mech_new my_chem_mech.in --atm_in atm_in --atm_in_new atm_in_new`

## Integration in NorESM build process
- `bld/definition.xml` : Defines name of the chemistry $chem for the perl files (e.g. fanci). Also includes CAM_AEROSOL_CONFIG_HASH in order for it to be cached (?)
- `bld/build-namelist` : Sets namelist values using the $chem definition -> some of these may need to be extracted to port the aerosol model! E.g. scavenging coefficients, etc.
- `bld/configure` : add pp_fanci chemistry, caches cam_aerosol_config_hash to config_cache, give priority to FANCI src paths
- `bld/namelist_files/use_cases` : Contains use cases for the FANCI compsets, currently only defaults
- `cime_config/config_component.xml` : Detailed specifications for FANCI compsets, what chemistry to use etc., also xml variable definitions for cam_aerosol_config and cam_aerosol_config_hash
- `cime_config/config_compsets.xml` : Name and definition for FANCI compsets
- `cime_config/buildcpp` : added cam_aerosol_config, and cam_aerosol_config_hash to add hash to cache
- `cime_config/buildnml` : Checks if aerosol config file xml variable is set. If so, run fanci_config.fanci_config, buildcpp and cache the hash for aerosol config. Later call fanci_config.add_FANCI_nl to add new namelists to atm_in
- `cime_config/QueryBuildCache.py`: script added to read xml files
- `src/chemistry/pp_fancil` : contains initial chem_mech.in and chemistry.F90 files (+ other initial chem files, not needed?)