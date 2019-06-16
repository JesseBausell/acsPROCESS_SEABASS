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
    a. c spectrum contains value less than zero or greater than 4 /m (400-700 nm)
    b. a spectrum contains value less than zero or greater than c value measured by the same channel (400-700 nm)
  8. Processed ac-s data are timestamped
  9. Produces Seabass-formatted ascii (.txt) file containing spectral offsets (a & c) calculated from holographic grating correction 
  10. Produces SeaBASS-formatted ascii (.txt) file containing time-stamped a/c spectra with depths at which they were sampled
  11. Produces SeaBASS-formatted ascii (.txt) file(s) containing depth-binned a/c average spectra and standard deviations. 
  12. Produces Matlab plot (.fig) detailing ac-s water column position (depth) over time (spectrum index). Plot will contain a user-selected reference point and a time-stamp (which is listed above the plot). Assuming ac-s was deployed simultaneously with hs6, hs6PROCESS_SEABASS.m will use this information to synchronize ac-s and hs6 data. In the event that ac-s and hs6 were deployed independently, the .fig file can be ignored.
  
User Instructions:
  1. Fill out metadata_HeaderFile_acs.txt (as specified below)
  2. Run acsPROCESS_SEABASS using Matlab command window.
  3. Select appropriate metadata_HeaderFile_acs.txt file when prompted. 
  4. Select appropriate pure-water absorption (MAT) file when prompted. (file is created using Purewater_SpecBuilder.m)
  5. Select appropriate pure-water attenuation (MAT) file when prompted. (file is created using Purewater_SpecBuilder.m)
  6. Time-stamp ac-s data for potential syncrhonization with hs6
    a. Matlab will produce a "time series" plot indicating ac-s depth over "time" (spectum index).
    b. User selects a reference point on this plot by entering the index of the desired point into the command window.***
    c. User is asked to confirm his/her selection on command window with y/n keys. If user rejects his/her selection, he/she will be           prompted to try again.
    
*** User-selected reference point (step 6b) can be used to synchroize ac-s cast with other instruments (provided they were deployed simultaneoulsy as a single unit). I recommend chosing the point at which ac-s is on the sea surface immediately before its first descent, as this point can be easily referenced for other instruments, regardless of differences in their sampling frequency. Matlab's data cursor (see figure toolbar) can assist user in determining the index of his/her desired point. 

Filling out metadata_HeaderFile_acs.txt:
acsPROCESS_SEABASS relies on metadata_HeaderFile_acs.txt to process ac-s data. All information (excluding pure-water MAT files) should be included in this header. A header template (metadata_HeaderFile_acs.txt) indicating important fields is provided in GitHub acsPROCESS_SEABASS repository. When filling out this header file, the first three headers (indicating user instructions) should be left alone. Required information fields contain = signs. USER SHOULD ONLY ALTER TEXT APPEARING ON THE RIGHT HAND SIDE OF =. User should indicate unavailability of desired information with "NA". DO NOT DELETE ROWS! Below are fields contained in metadata_HeaderFile_acs.txt and instructions on how to fill them out. Spaces should never be used in header fields; use underscore instead (_).

data_file_name= indicate name of ascii file containing unprocessed ac-s data. This file is generated using WET Labs Archive File Processing (WAP) software program, which merges a/c data with CTD data using nearest neighbor approach. Prior to running acsPROCESS_SEABASS, user must open  WAP-generated ascii file and maually indicate Conductivity, Temperature, and Depth column headers inside as follows: COND*, TEMP* and DEPTH*. All other column headers should remain untouched.

data_file_name=pathway for aforementioned WAP-generated ac-s ascii file (data_file_name). This pathway should include the folder in which sits, and should be ended using "/" or "\" for mac and pc respectively. 

affiliations=name of company or research institution with which investigators are affiliated. 

investigators=lists of investigators. Multiple names should be separated by commas and _ should be used in place of spaces.

contact=email of principle investigator

experiment=name of experiment or field campaign 

station=field station number 

latitude=latitude of field station. This should be indicated in decimal form. DO NOT format in minutes or seconds. Do not include Roman letters. South should be indicated with a negative sign.

longitude=longitude of field station. This should be indicated in decimal form. DO NOT format in minutes or seconds. Do not include Roman letters. West should be indicated with a negative sign.

documents=additional documents user wishes to submit to SeaBASS. DO NOT INDICATE kudelalab_ACS_readme.pdf. This is printed automatically in output files.

water_depth=bottom depth of the field station in meters. Numerals only. Do not include units.

calibration_files=names of original ac-s pure-water (DAT) calibration files from which MAT files were generated. Separate with a comma, and do not include spaces. Pure-water absorption file should come first. 

date(yyyymmdd)=indicate date on which ac-s was deployed.

start_time(military_time:HH:MM:SS)=military time at which ac-s cast was initiated. This is indicated in ac-s cast summary file.

time_lag(seconds)=elapsed time between cast initiation and data acquisition. This is indicated in ac-s cast summary file.

bin_size=desired depth bin-sizes for binning. User can include as many as he/she wishes.

first_grating=highest chanel number found on ac-s first holographic grating. This can also be thought of as the wavelength index directly before (left side) the "jump" in unprocessed a/c spectra.
