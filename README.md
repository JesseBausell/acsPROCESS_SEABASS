# acsPROCESS_SEABASS

This Matlab script processes raw spectral absorption (a) and attenuation (c) measurements, as sampled in natural water bodies using Wet Labs Spectral Absorption and Attenuation Sensor (ac-s). acsPROCESS_SEABASS uses up to date processing protocols (as of June 2019) to process ac-s data, which it outputs as both individual and depth-binned spectra. All output data products are formatted as to be consistent with NASA's SeaWiFS Bio-Optical Archive and Storage System (SeaBASS).

Inputs:
metadata_HeaderFile_acs.txt - ascii file (.txt) containing metadata required to process raw ac-s data (see below)
purewater_absorption.mat - MAT file containing pure-water a calibration spectrum (see Purewater_SpecBuilder)
purewater_attenuation.mat - MAT file containing pure-water c calibration spectrum (see Purewater_SpecBuilder)

Outputs:
Station_#_acs.txt - Seabass-formatted ascii file containing processed ac spectra
Station_#_a_bin#.txt* - Seabass-formatted ascii file(s) containing processed & depth-binned a spectra
Station_#_c_bin#.txt* - Seabass-formatted ascii file(s) containing processed & depth-binned c spectra
Station_#ac_offsets.txt - Seabass-formatted ascii file containing the spectral offsets that were used to correct for the reliance of ac-s  on two holographic gratings for sampling
Station_#_ACS.fig - figure depicting vertical position of ac-s across time (of cast). Used to time-synchronize ac-s and hs6 cast data.
* - multiple files can be made, depending on how many depth-bin sizes user specifies in metadata_HeaderFile_acs.txt

Required Matlab Scripts and Functions:
ACS_fileREADER_MS.m
BIN_acDATA.m
HoloGrater_ACS_4.m
lambda_INTERPOLATE.m
metaData_Reader.m
PureWater_Extract3.m
salinityCOND.m
sullivan_VALS.m

Required data files:
Seabass_header_ACS5.mat
Sullivan_VALUES.mat

Program Description:
acs_PROCESS_SEABASS processes raw field-collected ac-s measurements following a series of steps. It is outfitted to process raw data contained in ac-s ascii files regardless of ac-s channel number, wavelengths, or orientation of column field headers. Steps are outlined below:
  1. Reads ascii data into Matlab
  2. Calculates water column salinity using measurements conductivity (CTD)
  3. Performs correction for spectral "jumps" caused by ac-s sampling using two holographic gratings
  4. Subtracts a/c pure-water calibration spectra from field-measured ac-s data**
  5. Corrects for the optical effects temperature and salinity using Sullivan et al. (2006)
  6. Corrects for scattering using Rottgers et al. (2013)
  7. QA/QC ac-s data. Paired ac spectra are flagged and removed if:
    a. c spectrum contains value less than zero or greater than 4 /m (channel centered at wavelength < 715 nm)
    b. a spectrum contains value less than zero or greater than c value measured by the same channel (channel centered at wavelength < 715 nm)
