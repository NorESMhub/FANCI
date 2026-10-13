# ==============================================================================
# Python script that reads config.ini with parameter settings
# for the FANCI sectional aerosol model in NorESM
# Output:
# atm_in: edited initial atm_in namelist file with settings for FANCI aerosols
# my_chem_mech.in: edited chem_mech.in to contain aerosol tracers
# ==============================================================================

import os
import logging
import configparser
import argparse
import math
from datetime import datetime

#logging.basicConfig(level=logging.INFO) # basic level would be "warning"
logger = logging.getLogger("fanci_config")

class AeroConfigError(Exception):
    pass

def _parse_range(config, section, variable):
    ''' Parse a range from the config file and make a list
    Parameters:
        config : Instance of ConfigParser class (aerosol config file)
        section (str) : Section of the config.ini file in [] that should be read
        variable (str) : variable in the config file to be read
    Returns:
        parsed_range (list(float)) : The values of the "variable"
        '''
    range_str = config.get(section, variable).split(',')
    parsed_range = [float(i) for i in range_str]
    return parsed_range

# ==============================================================================
# Types for bin, range and species information
# ==============================================================================
class _AerosolSpecies:
    ''' Class representing an aerosol species.

    Attributes:
        active (bool) : whether the species should be included in the simulation
        short_name (str) : short name of the aerosol species
        species_type (str) : type of the species, e.g. 'dust', 's-organic', 'p-organic'
        long_name (str) : long name of the aerosol species
        composition (str) : composition to trace aerosol mass
        mixed (bool) : True for internally mixed, false for externally mixed aerosol
        range_bnds (list(float)) : range bounds within which the species exists
        range_idx (list(int)) : indices of ranges within which the species exists
        kappa (float) : species specific hygroscopicity parameter
        molecular_weight (float) : molecular weight of the aerosol species
    '''

    def __init__(self, config, species):
        ''' Initializes an aerosol species object with values from the aerosol
            configuration file.

        Parameters:
            config : Instance of ConfigParser class (aerosol config file)
            species (str) : Section with aerosol name of the config.ini file in [] that should be read
        '''
        try:
            self.active = config.getboolean(species, 'active', fallback=False)
            self.short_name = config.get(species, 'short_name')
            self.species_type = config.get(species, 'species_type')
            self.long_name = config.get(species, 'long_name', fallback=self.short_name)
            self.composition = config.get(species, 'composition')
            self.density = config.getfloat(species, 'density')
            self.molecular_weight = config.getfloat(species, 'molecular_weight')
            self.mixed = config.getboolean(species, 'mixed', fallback=True)
            self.range_bnds = _parse_range(config, species, 'range_bounds')
            self.kappa = config.getfloat(species, 'kappa')
            logger.info(f"Successfully parsed species attributes for '{species}'")
        except Exception as e:
            logger.error(f"Error parsing species attributes for '{species}': {e}")
            raise AeroConfigError("ERROR: Species attributes in configuration file missing")
        self.range_idx = []


    def check_range_bnds(self, range_bnds):
        ''' Check whether range bounds for the individual species
        in the configuration file are valid (that is equal to the specified range bounds)
        and raise an exception if not.

        Parameters:
            range_bnds (list(float)) : list of range bounds read from the config file
        '''
        if not all(i in range_bounds for i in self.range_bnds):
            logger.error(f"Invalid range bounds: {self.range_bnds}")
            raise AeroConfigError(f"ERROR: Species range bound not equal to range bounds")

    def get_range_idx(self, range_bnds):
        ''' Adjust the species bounds to the new range bounds
        (adjusted in the _RangeSpecs adjust_range_bnds function).

        Parameters:
            range_bnds list(float) : The range bounds that the species range bounds should be adjusted to
        Attributes:
            range_bnds list(float) : The species specific range bounds of the class _AerosolSpecies to be adjusted
            range_idx list(int) : The indices of the ranges the species live in
        '''
        abs_diff_lo = [abs(rb - self.range_bnds[0]) for rb in range_bnds] # calculate absolute difference to each range_bound value
        idx0 = abs_diff_lo.index(min(abs_diff_lo))                        # get the index from this range bound value
        abs_diff_hi = [abs(rb - self.range_bnds[1]) for rb in range_bnds]
        idx1 = abs_diff_hi.index(min(abs_diff_hi))
        self.range_bnds = [range_bnds[idx0], range_bnds[idx1]]
        # get range indices for each species
        self.range_idx.extend([i for i in range(1, len(range_bnds))
                                if range_bnds[i] > self.range_bnds[0]
                                and range_bnds[i] <= self.range_bnds[1]])

# TODO: bin settings, maybe add 'method' to allow for other than volume ratio. e.g. 'custom' -> user defined bin bounds
class _BinSpecs:
    ''' Class representing all attributes associated with the size distribution in the FANCI sectional aerosol model.

    Attributes:
        N (int) : Total number of bins read from the sectional aerosol configuration file
        r_1 (float) : Center radius of the smallest bin (nm)
        r_N (float) : Center radius of the largest bin (nm)
        r (list(float)) : list of the center radii of all bins (nm)
        r_bnds (list(float)) : list of the radii at bin boundaries (nm)
    '''

    def __init__(self, config):
        ''' Initializes an instance of the _BinSpecs class

        Parameters:
            config : Instance of ConfigParser class (aerosol config file)
        '''
        self.N = config.getint('BIN SPECS', 'nbin')
        self.r_1 = config.getfloat('BIN SPECS', 'radius_1')
        self.r_N = config.getfloat('BIN SPECS', 'radius_N')
        self.r = None
        self.r_bnds = None

    def calc_bins(self):
        ''' Calculates the actual bin specifications with the given number of bins and radius
        information from the sectional aerosol configuration file using the
        volume ratio approach described in Jacobson, M. Z. (2005). Fundamentals of Atmospheric Modeling.
        Cambridge University Press., p. 451 ff.

        Attributes:
            N (int) : Total number of bins read from the sectional aerosol configuration file
            r_1 (float) : Center radius of the smallest bin (nm)
            r_N (float) : Center radius of the largest bin (nm)
            r list(float) : List of center radii of all bins (nm)
            r_bnds list(float) : list of radii at bin boundaries (nm)
        '''
        V_rat = (self.r_N / self.r_1) ** (3 / (self.N-1))             # volume ratio - Formula (13.3) in Jacobson
        v_0 = 4/3 * math.pi * (self.r_1) ** 3                         # smallest volume - sperical volume for smallest center bin
        v = [v_0 * V_rat ** i for i in range(self.N)]                 # All volumes - Formula (13.2) in Jacobson
        v_lo = [(2*v[i]) / (1+V_rat) for i in range(self.N)]          # lower radius bounds - Formula (13.7) in Jacobson
        v_hi = V_rat*v_lo[-1]                                         # upper bnd for largest bin - Formula (13.6) in Jacobson
        v_bnds = v_lo + [v_hi]                                        # all bounds to one array
        self.r = [(v[i]*3/(4*math.pi))**(1/3) for i in range(self.N)] # center radii - from volumes
        self.r_bnds = [(v_bnds[i]*4/(3*math.pi))**(1/3) for i in range(self.N+1)] # radius bounds - from volumes

class _RangeSpecs:
    ''' Class representing all attributes associated with a chemical composition range in the
    FANCI sectional aerosol model.

    Attributes:
        ranges (bool) : True if the model should average the chemistry for a range of bins
        range_bnds (list(float)) : The radii at the range boundaries (nm)
        range_bnd_bin_idx (list(int)) : Indices of bins within a range
        nspecies (list(int)) : number of species in each range
    '''
    def __init__(self, config):
        ''' Initializes an instance of the _RangeSpecs class.

        Parameters:
            config : Instance of ConfigParser class (aerosol config file)
        '''
        self.ranges = config.getboolean('RANGE SPECS', 'ranges', fallback=True)
        self.range_bnds = _parse_range(config, 'RANGE SPECS', 'range_bounds')
        self.range_bnd_bin_idx = []
        self.nspecies = []
    def adjust_range_bnds(self, r_bnds):
        ''' Function to adjust the soft range bounds given in the configuration file to the
        radius bounds calculated in the _BinSpecs calc_bins routine.

        Parameters:
            r_bnds (list(float)) : Radius bounds (nm) calculated in the _BinSpecs calc_bin function
        Attributes:
            ranges (bool) : True if the model should average the chemistry for a range of bins
            range_bnds (list(float)) : The radii at the range boundaries (nm)
            range_bnd_bin_idx (list(tuple(int,int))) : Indices of bins within a range (lower bin index, upper bin index)
        '''
        if self.ranges:
            self.range_bnds[-1] = r_bnds[-1] # set lowest range bound equal to lowest bin bound
            for i in range(0, len(self.range_bnds)-1):
                # find closest bin bound to each range bound
                abs_diff_lo = [abs(radb - self.range_bnds[i]) for radb in r_bnds]
                abs_diff_hi = [abs(radb - self.range_bnds[i+1]) for radb in r_bnds]
                idx_lo = abs_diff_lo.index(min(abs_diff_lo))
                idx_hi = abs_diff_hi.index(min(abs_diff_hi))
                self.range_bnds[i] = r_bnds[idx_lo] # set the range_bnd to the radius_bnd that is closest
                self.range_bnd_bin_idx.append((idx_lo+1, idx_hi)) # idx_lo +1 -> first one will be 1, since fortran arrays start at 1
                                                                  # idx_hi -> not plus one, since this is the index of the r_bnds and technically we would need to calculate -1 to get fortran indices
        else:
            # If no ranges (classic sectional scheme), set range bounds equal to bin bounds TODO: is this necessary?
            self.range_bnds = r_bnds
    def get_nspecies(self, species_obj_list):
            n = len(self.range_bnds)-1
            self.nspecies = [0] * n
            for r in range(n):
                self.nspecies[r] = sum([1 for obj in species_obj_list if r+1 in obj.range_idx])

def fanci_config(aerconf_file, chemconf, chem_infile, FANCI_in):
    ''' Main function called from buildnml if a compset with the FANCI sectional aerosol
    model is used. This function is used to initialize instances of the classes above using
    information from the sectional aerosol configuration file.
    It will create a new chemistry mechanism file with added tracers for each aerosol species,
    as well as a temporary namelist file with the necessary parameters for the model.

    Parameters:
        aerconf_file (str) : The full path to the aerosol configuration file.
                             The name of the file can be changed with the xml variable
                             CAM_AEROSOL_CONFIG. Currently the path is set to
                             srcroot/src/chemistry/FANCI/config/CAM_AEROSOL_CONFIG
                             in the buildnml script.
        chemconf (str) :     The full path to the chemistry mechanism file. This is the chem_mech.in
                             file in the pp_ chemistry folder associated with the compset
        chem_infile (str) :  The full path to the modified chemconf file with added aerosol tracers.
                             This file is then added to the casefolder.
        FANCI_in (str) : Full path where the temporary namelist for the sectional aerosol model
                             is written out. This file is deleted in buildnml after the contents are added
                             to atm_in
    '''

    # ==============================================================================
    # Read input from config file
    # ==============================================================================
    # INPUT
    config = configparser.ConfigParser()
    try:
        config.read(aerconf_file)
    except:
        raise AeroConfigError('ERROR: Config file does not exist or bad file format')

    # ==============================================================================
    # Initialize
    # ==============================================================================´

    # initialize bins
    bin_specs = _BinSpecs(config)
    bin_specs.calc_bins()

    # initialize ranges
    range_specs = _RangeSpecs(config)

    # initialize aerosol species
    species_obj_list = []
    for section in config.sections():
        if section not in ['BIN SPECS', 'RANGE SPECS']:
            species_obj_list.append(_AerosolSpecies(config, section))

    [species_obj.check_range_bnds for species_obj in species_obj_list]  # check species range bnds

    # adjust range range_bnds
    range_specs.adjust_range_bnds(bin_specs.r_bnds)

    # get species range indices
    [species_obj.get_range_idx(range_specs.range_bnds) for species_obj in species_obj_list]

    active_species_obj_list = [species for species in species_obj_list if species.active]
    nspecies_tot = len(active_species_obj_list)
    range_specs.get_nspecies(active_species_obj_list)
    # =====================================================================
    # Write to temporary FANCI namelist file
    # =====================================================================

    with open(FANCI_in, "w") as f:
        f.write("&fanci_properties_nl\n")
        f.write(" fanci_nspecies_tot       =  ")
        f.write(f"{nspecies_tot} \n")
        f.write(" fanci_nbins       =  ")
        f.write(f"{bin_specs.N} \n")
        f.write(" fanci_nranges       =  ")
        f.write(f"{len(range_specs.range_bnds)-1} \n")
        f.write(" fanci_nspecies       =  ")
        for i in range(0, len(range_specs.range_bnds)-1):
            f.write(f"{range_specs.nspecies[i]}")
            if i != len(range_specs.range_bnds)-2:
                f.write(',')
        f.write("\n")
        f.write(" fanci_bin_bounds     =  ")
        for i in range(bin_specs.N):
            f.write(f"'{bin_specs.r_bnds[i]:.3f}D0:{bin_specs.r_bnds[i+1]:.3f}D0'")
            if i != bin_specs.N-1:
                f.write(', ')
        f.write("\n")

        f.write(" fanci_bin_centers        =  ")
        for i in range(bin_specs.N):
            f.write(f"'{bin_specs.r[i]:.3f}D0'")
            if i != bin_specs.N-1:
                f.write(', ')
        f.write("\n")

        f.write(" fanci_range_bounds       =  ")
        for i in range(len(range_specs.range_bnds)-1):
            f.write(f"'{range_specs.range_bnd_bin_idx[i][0]}:{range_specs.range_bnd_bin_idx[i][1]}'")
            if i != len(range_specs.range_bnds)-2:
                f.write(', ')
        f.write("\n")
        f.write("/\n")

        for species in species_obj_list:
            if species.active:
                f.write("&fanci_properties_aerosol_nl\n")
                f.write(f" fanci_aerosol_name      =  '{species.short_name}' \n")
                f.write(f" fanci_aerosol_type = '{species.species_type}' \n")
                f.write(" fanci_aerosol_range       =  ")
                f.write(f"'{species.range_idx[0]}:{species.range_idx[-1]}' \n")
                f.write(f" fanci_aerosol_mixed     =  .{species.mixed}. \n")
                f.write(f" fanci_aerosol_density   =  {species.density:.3f} \n")
                f.write(f" fanci_aerosol_weight    =  {species.molecular_weight:.3f} \n")
                f.write(f" fanci_aerosol_kappa =    {species.kappa:.3f} \n")
                f.write("/\n")


    # ==============================================================================
    # Prepare output for chem_mech.in file
    # The code below reads the chem_mech.in file line by line and looks for keywords
    # "Solution" and "Implicit". Below these, the composition of the tracers and the
    # names of the tracers are added. The file is then written out to the name and
    # path specified in buildnml
    # ==============================================================================

    composition_list = []
    implicit_list = []

    for species in species_obj_list:
        if species.active:
            for range_idx in species.range_idx:
                composition_list.append(
                    f"{species.short_name}_R{range_idx} -> {species.composition}"
                    ) # test
                implicit_list.append(
                    f"{species.short_name}_R{range_idx}"
                    )

    for i in range(1,bin_specs.N+1):
        composition_list.append(
            f"num_{i} -> H"
            )
        implicit_list.append(
            f"num_{i}"
    )

    with open(chemconf, 'r') as chem_file:
        lines = chem_file.readlines()

    modified_chem = []

    for line in lines:
    # write lines from the old file for the new file
        modified_chem.append(line)
    # add species composition
        if 'Solution' in line and not 'End' in line and not 'Classes' in line:
            for composition in composition_list:
                if any(chemline == f"{composition}\n" for chemline in lines): # check if the lines have already been added to chem_mech
                    pass
                else:
                    modified_chem.append(f"{composition}\n")
    # add species for advection
        if 'Implicit' in line and not 'End' in line:
            for implicit in implicit_list:
                if any(chemline == f"{implicit}\n" for chemline in lines):
                    pass
                else:
                    modified_chem.append(f"{implicit}\n")
    # write out to my_chem_mech.in
    with open(chem_infile, 'w') as file:
        file.writelines(modified_chem)

def add_FANCI_nl(tmp_atm_nlfile, FANCI_in, atm_nlfile):
    ''' Function to modify the atm_in namelist file and add the FANCI namelists
    to it. Called by buildnml. The namelists in the original atm_in file are sorted alphabetically,
    the new namelists are inserted in alphabetical order.
    Each species gets its own namelist, all are named & FANCI_properties_aerosol_nl and
    iterated through by the nl reader.

    Parameters:
        tmp_atm_nlfile : full path to the original atm_in file, temporarily moved to tmp_atm_in
        FANCI_in : full path to the FANCI sectional namelists created in fanci_config.fanci_config
        atm_nlfile : full path to the final atm_in, combined tmp_atm_nlfile and FANCI_in
    '''
    modified_atm_in = []

    with open(tmp_atm_nlfile, 'r') as f1:                  # read atm_in
        lines = f1.readlines()
    with open(FANCI_in, 'r') as f2:                # read sectional nl file
        lines_fanci = f2.readlines()
    idx = 0
    for line in lines:
        if "&" in line and line > lines_fanci[0]:        # write entries before FANCI nl
            break
        else:
            modified_atm_in.append(line)
            idx += 1
    for line_oslo in lines_fanci:                        # write FANCI nl
        modified_atm_in.append(f"{line_oslo}")
    for line in lines[idx:]:                                # write entries after FANCI nl
        modified_atm_in.append(f"{line}")

    with open(atm_nlfile, 'w') as atm_infile:               # write out modified nlfile to "atm_in"
        atm_infile.writelines(modified_atm_in)

def _main_func():
    parser = argparse.ArgumentParser(description=("Process aerosol configuration for the"
    "sectional aerosol scheme in NorESM, write the namelist and add the tracers to chemistry."))
    parser.add_argument('--aerconf', required=True, help='Path to the aerosol configuration file')
    parser.add_argument('--chem_mech', required=True, help='Path to the initial chem_mech.in file')
    parser.add_argument('--chem_mech_new', required=True, help='Path to the modified chem_mech.in file')
    parser.add_argument('--atm_in', required=True, help='Path to the original atm_in file')
    parser.add_argument('--atm_in_new', required=True, help='Path to new atm_in file with sectional aerosol info')
    args = parser.parse_args()

    FANCI_in = f'FANCI_tmp_{datetime.today().strftime('%Y%m%d%H%M%S')}' # temporary nl file with sectional info

    try:
        if not os.path.isfile(args.aerconf):
            raise AeroConfigError('ERROR: Specified aerosol configuration file does not exist')
        if not os.path.isfile(args.chem_mech):
            raise AeroConfigError('ERROR: Specified chem_mech.in file does not exist')
        if not os.path.isfile(args.atm_in):
            raise AeroConfigError('ERROR: Specified atm_in file does not exist')
    except AeroConfigError as errmsg:
        logger.error(errmsg)


    fanci_config(args.aerconf, args.chem_mech, args.chem_mech_new, FANCI_in)
    add_FANCI_nl(args.atm_in, FANCI_in, args.atm_in_new)

    os.remove(FANCI_in)

if __name__ == "__main__":
    _main_func()