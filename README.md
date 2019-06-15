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

