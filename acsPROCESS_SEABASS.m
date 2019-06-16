% acsPROCESS_SEABASS
% Jesse Bausell
% September 13, 2017
%
% This matlab script processes WET Labs ac-s data for submission of NASA's 
% SEABASS data repository. The program requires a metadata header file 
% and pure-water absorption and attenuation specra (see readme). It 
% formats data into SEABASS-submittable .txt files. These contain single
% specra as well as data binned to user specifications (as indicated in
% metadata header file.
%
% Required scripts and functions:
% ACS_fileREADER_MS.m
% BIN_acDATA.m
% HoloGrater_ACS_4.m
% lambda_INTERPOLATE.m
% metaData_Reader.m
% PureWater_Extract3.m
% salinityCOND.m
% sullivan_VALS.m
% 
% Required data files:
% Seabass_header_ACS5.mat
% Sullivan_VALUES.mat

clear; close all; clc;

ACS_fileREADER_MS;
% This program will read a WetLabs ac-s data file into our program. It
% produces raw absorption and attenuation matrices, as well as temperature,
% salinity, and depth arrays.

%%   1. Performs holographic grating correction
% WetLabs ac-s takes measurements using two holographic gratings. Its 87-89
% channels are partitioned onto these gradients, which take measurements
% asynchronously. This can cause an artifact in raw ac-s spectra, exhibited
% as a "jump" in individual spectrum. This section of code corrects for
% this jump. It also collects the offsets that are interpolated in order to
% correct spectra.

diFF_C = nan(l_AC,1); % Nan array to deposit attenuation offsets 
diFF_A = nan(l_AC,1); % Nan array to deposit absorption offsets 

for hh = 1:l_AC 
    % Corrects holographic grating artifact for each individual
    % absorption/attenuation spectrum
    [C_CORR(hh,:), diFF_C(hh)] = HoloGrater_ACS_4(c_wv,C_CORR(hh,:),grate_IND); %attenuation
    [A_CORR(hh,:), diFF_A(hh)] = HoloGrater_ACS_4(a_wv,A_CORR(hh,:),grate_IND); %absorption
end
    
%%   2. Subtracts pure water values 
% This section corrects ac-s spectra by substracting pure water absorption
% and attenuation spectra from measured spectra.


[A_CORR, C_CORR, wcal_a_infile, wcal_c_infile] = PureWater_Extract3(A_CORR,C_CORR,nanmean(T_insitu),grate_IND);
% This line of code subtracts pure water a and c from ACS spectra. It
% prompts the user to choose two .mat files (absorption & attenuation spectra 
% averaged from ACS pure water instrument calibration files (.dat). The program than normalizes them to a
% temperature of 15 degrees C and subtracts them (original temperature is
% included in the .mat files).

%%   3. Performs Sullivan et al. temperture and salinity corrections
% Absorption and attenuation are impacted by temperature and salinity. 
% Sullivan et al. (2006) (see readme) provides temperature adjustment coefficients
% (used in the function PureWater_Extract3), as well as salinity adjustment 
% coefficients for absorption and attenuation spectra. Salinity adjustment
% coefficients assume an ambient water temperature of 15 degrees (C).
% Using Sullivan et al. (2006), ac-s data are adjsuted to a standard
% temperature of 15 degrees, salinity corrected, and then re-adjusted to
% the average ambient water temperature in which they were sampled.

Tnorm = 15; %standard water temperature temperature for salnity adjustement.
load('Sullivan_VALUES');  %load sullivan parameters

for ii = 1:length(c_wv) 
    % This loop normalizes absorption and attenuation data one wavelength
    % at a time. 
       
    [temp_DEP_a, salt_DEP_a] = sullivan_VALS(a_wv(ii),Sullivan_VALUES);
    [temp_DEP_c, salt_DEP_c, salt_DEP_c] = sullivan_VALS(c_wv(ii),Sullivan_VALUES);
    % The above function retrieves the Sullivan et al. (2006) temperature  
    % and salinity adjustment coefficients for absorption and attenuation 
    % for a given wavelength. Because absorption and attenuation
    % wavelengths are different, the program is run for each one
    % separately.
    
    C_tNORM = C_CORR(:,ii)-(T_insitu - Tnorm)*temp_DEP_c;
    A_tNORM = A_CORR(:,ii)-(T_insitu - Tnorm)*temp_DEP_a;
    % Adjust absorption and attenuation to a temperature of 15 degrees C
    % using Sullivan et al. (2006) adjustment coefficients.
    
    C_tNORM1 = C_tNORM - S_insitu*salt_DEP_c;
    A_tNORM1 = A_tNORM - S_insitu*salt_DEP_a;
    % Salinity-adjust absorption and attenuation using salinity measured 
    % by the ac-s CTD.
    
    C_tNORM2 = C_tNORM1-(Tnorm-T_insitu)*temp_DEP_c;
    A_tNORM2 = A_tNORM1-(Tnorm-T_insitu)*temp_DEP_a;
    % Re-adjust absorption and attenuation back to the original temperature
    % in which ac-s was deployed.
    
    C_CORR(:,ii) = C_tNORM2;
    A_CORR(:,ii) = A_tNORM2;
    % Replace absorption and attenuation matrices with salinity-corrected
    % data.
    
end

%%   4. Performs Rottgers et al. scattering correction. 
% ac-s measurements are subject to scattering-induced errors which can
% inflate field-sampled absorption. Absorption data are corrected for these
% errors using methods described in Rottgers et al. (2013) (see readme)

a_m715_INT = nan(length(deptH),1); 
c_m715_INT = nan(length(deptH),1);
% Creates nan array to place absorption and attenuation interpolated at
% a wavelength of 715 nm.

    for ii = 1:length(deptH)
    % Interpolate absorption and attenuation to wavelength of 715 nm (a715).
       a_m715_INT(ii) = lambda_INTERPOLATE(a_wv,A_CORR(ii,:),715); %absorption
       c_m715_INT(ii) = lambda_INTERPOLATE(c_wv,C_CORR(ii,:),715); %attenuation
    end

% Remove negative a715 values in order to prevent them from later becoming
% imaginary. Here they are replaced by NaNs.
a_715_IND = find(a_m715_INT < 0); %Find index of negative a715
a_m715_INT(a_715_IND) = NaN; % Replace negative values with NaN.

% Create variables for Rottgers et al. (2013) equations
e_c_INVERSE = 1; %constant in the correction equation (1/ec)
a_715 = 0.212*(a_m715_INT.^1.135); % Adjusts asorption at 715nm to account for scatter

% denominator component for Rottgers scatter correction equation
deNOM = (e_c_INVERSE*c_m715_INT-a_m715_INT); 

for ii = 1:length(a_wv) 
    % Correct absorption for scatter-related errors by wavelength usign
    % Rogttgers et al. (2013)'s scatter correction method, Equation 2b.
    A_CORR(:,ii) = A_CORR(:,ii) - (a_m715_INT-a_715).*((e_c_INVERSE*C_CORR(:,ii)-A_CORR(:,ii))./deNOM);
end

% Create referece variables for time and depth before they are altered.
deptH2 = deptH; % Create an alternate depth variable. Important in time stamp ID
time_insitu2 = time_insitu; % Create an alternate time index. Important in time stamp ID

% As explained above, absorption spectra for which a715<0 were eliminated
% before Rottgers correction was applied. Because each absorption spectrum
% has a corresponding depth and attenuation, depths and attenuations which
% correspond to eliminated absorption spectra must also be flagged and
% removed.
C_CORR(a_715_IND,:) = NaN; %Remove corresponding attenuation spectra 
deptH(a_715_IND) = NaN; %Remove corresponding depth
time_insitu(a_715_IND) = NaN; %Remove corresponding time

%%   5. QA/QC Data (automated)
% Automoatic quality control for ac-s absorption and attenuation data.
% Spectra with attenuation values greater than 4 and aborption values greater 
% than their corresponding attenuation (wavelength dependent) are flagged, 
% replacing them with NaN. Spectra with negative absorption/attenuation values 
% are also flagged for removal. NOTE: If an absorption spectrum is flagged,
% it's corresponding attenuation spectrum is also flagged and vice versa,
% however this is NOT applied to wavelengths greater than 715 nm (due to
% the Rottgers correction).
c_min = min(C_CORR(:,1:79),[],2); % find minimum attenuation value for each spectrum
c_max = max(C_CORR(:,1:79),[],2); % find maximum attenuation value for each spectrum
[cERR_y, cERR_x] = find(c_min < 0 | c_max > 4); 
% Find row indices of all rows with attenuations less than zero OR greater than 4.

% Flag all attenuation spectra (and corresponding absorption, depth, time stamps, 
% and holographic grating offsets) that contain negative values or values
% greater than 4 m^-1.
C_CORR(cERR_y,:) = NaN;
A_CORR(cERR_y,:) = NaN;
diFF_C(cERR_y) = NaN;
diFF_A(cERR_y) = NaN;
deptH(cERR_y) = NaN;
time_insitu(cERR_y) = NaN;

a_min = min(A_CORR(:,1:79),[],2); % Find minimum absorption value in each spectrum
% Find instances where absorption > attenuation at wavelengths < 715nm
booL = A_CORR(:,1:79) > C_CORR(:,1:79); % Find where absorption > attenuation
booL_max = max(booL,[],2); % Find spectra where absorption > attenuation
% Find indices at which absorption < 0 or greater absorption > attenuation
[aERR_y, aERR_x] = find(a_min < 0 | booL_max > 0); 

% Flag all absorption spectra (and corresponding attenuation, depth, time stamps, 
% and holographic grating offsets) that contain negative values or values
% greater than corresponding attenuation.
C_CORR(aERR_y,:) = NaN;
A_CORR(aERR_y,:) = NaN;
diFF_C(aERR_y) = NaN;
diFF_A(aERR_y) = NaN;
deptH(aERR_y) = NaN;% = abs_MATRIX{end-1}; %gives us the depth array of RR_MATRIX
time_insitu(aERR_y) = NaN;% = time_insitu{1}; %gives us the time in milliseconds
% NaN all rows that don't fit our selection criteria

nan_IND = find(isnan(deptH)); % Ensure that all NaN flagged spectra are accounted for

% Eliminate NaN flagged spectra
diFF_C(nan_IND) = [];
diFF_A(nan_IND) = [];
C_CORR(nan_IND,:) = [];
A_CORR(nan_IND,:) = [];
deptH(nan_IND) = [];
time_insitu(nan_IND) = [];

%%   6. Create Time Variable
% Seabass requires a timestamp for each processed ac spectrum. This
% requires converting seconds into military time using the time at which
% ac-s was deployed (start time) as well as the time lag, which represents
% the elapsed time between ac-s being turned on and data collection being
% initiated.

cruise = experiment; % Rename the experiement variable for use in header
load('Seabass_header_ACS5.mat'); % Generic header field data
% Add time lags to time arrays 
time_insitu = time_insitu + str2double(lag); % time of Qa/QC'd spectra (sec)
time_insitu2 = time_insitu2 + str2double(lag); % original unaltered time (sec)
% Convert time_insitu into decimal day
time_insitu = round(time_insitu)/(24*3600);
time_insitu2 = time_insitu2/(24*3600);
% Add starting time to time_insitu. All in mat time.
time_insitu = time_insitu+datenum(startTIME);
time_insitu2 = time_insitu2+datenum(startTIME);
% Craete cell arrays to hold time strings
tiME = cell(length(time_insitu),1);
tiME2 = cell(length(time_insitu2),1);

for ii = 1:length(time_insitu)
    % Create time strings for each QA/QC'd ac spectrum by converting times
    % (in seconds) into military time
    tiME{ii} = datestr(time_insitu(ii),'HH:MM:SS');
end

for ii = 1:length(time_insitu2)
    % Create time strings for reference time variable by converting it from 
    % its current form (seconds) into military time
    tiME2{ii} = datestr(time_insitu2(ii),'HH:MM:SS.FFF');
end


%%   7. Select and plot a reference time
% In the event that Hydroscat-6 backscattering data sampled concurrently
% with ac-s data are to be submitted along side of them, it is paramount
% that ac-s and hs6 time stamps are synchronoized. Seabass_hs6PROCESS
% automates this process using a visual reference point produced in this
% section. 

keY = 1; % Used to cycle through the while loop
while 1
    % This while loop prompts user to select a reference point using the
    % cursor
    if isequal(keY,1)
        disp('Select a point as reference time'); %displays text
        plot(deptH2,'b','LineWidth',2); %plots depth over index
        xlabel('Index','FontSize',20); ylabel('Depth (m)','FontSize',20); %labels
        title('ACS Cast Depth','FontSize',20); % title of figure
        set(gca,'ydir','reverse'); %invert y-axis
        time_IND = input('Select a time point: '); % prompts user to select an index
        close all; clc; keY = keY + 1; % advances to the next section
    
    elseif isequal(keY,2)
        plot(deptH2,'b','LineWidth',2); hold on; %plots depth over index (again)
        scatter(time_IND,deptH2(time_IND),'or','MarkerFaceColor','r'); %plots the user-selected point
        xlabel('Index','FontSize',20); ylabel('Depth (m)','FontSize',20); %x and y labels
        title('ACS Cast Depth','FontSize',20); % title
        set(gca,'ydir','reverse'); %inverts the y-axis
    
        while 1
            % asks user if he/she is satisfied
            question = input('Accept this point as reference time? (y/n): ','s');
            if strcmpi(question,'y') % if user is satsified with reference time, break the loop
                keY = keY+1; break
            elseif strcmpi(question,'n') % if user is unsatisfied, start over again
                keY = 1; close all; break
            else % if the selection is nonsensical
                disp('Invalid selection. Must select y or n.')
            end
        end
    
    else % if user originally indicated satsifaction...
        title(tiME{time_IND});  % put the reference time on the figure as the title
        savefig([in_DIR cruise '_' station '_ACS']); % save the figure for future reference
        break
    end 
end

%% 8. Create an associated holographic offset file
% Calculated offsets used to adjust spectral "jumps" between the two ac-s
% holographic gratings are placed into a text file. The file name is also
% included onto the header of the processed ac-s data as an 'associated
% file'.

offset_NAME = [cruise '_' station 'ac_offsets' '.txt']; % File name
fid_off = fopen([in_DIR offset_NAME],'w'); % Create a file with the above name
% Headers
fprintf(fid_off,'%s\n','Time,Depth,Offset_agp,Offset_cgp');
fprintf(fid_off,'%s\n','HH:MM:SS,m,1/m,1/m');
fprintf(fid_off,'%s\n',['! This file contains the ac offsets used to offset the ac spectral "jumps" in ' cruise '_' station ]);
fprintf(fid_off,'%s\n','! This file contains the ac offsets used to offset the ac spectral "jumps"');

for ii = 1:length(deptH)
    % Add offset data on spectrum at a time
    fprintf(fid_off,'%s\t',tiME{ii});    
    fprintf(fid_off,'%f\t',deptH(ii));
    fprintf(fid_off,'%f\t',diFF_A(ii));   
    fprintf(fid_off,'%f\n',diFF_C(ii));    
end

fclose(fid_off); % Close the file to complete it

%%   9. Create the ac-s data file formatted for Seabass.
% This section of code produces a .txt file of processed ac-s data with a
% header that is formatted for Seabass.

units = []; % empty repository to place units of measurement for variables
printFODDER = []; % empty repository to place file identifiers for variables
fields_a = []; % empty repository for absorption wavelenths
fields_c = []; % empty repository for attenuation wavelengths
for ii = 1:length(a_wv)
    % For every ac-s channel create the following
    fields_a = [fields_a ',agp' num2str(a_wv(ii))]; % absorption header
    fields_c = [fields_c ',cgp' num2str(c_wv(ii))]; % attenuation header
    units = [units ',1/m']; % unit of measurement
    printFODDER = [printFODDER '%1.6f\t']; %format specifier
end    

% Combines depth array with absorption and attenuation matrices
acs_MAT = [deptH A_CORR C_CORR]; 
% Create text file for Seabass-ready ac-s data
fid_acs = fopen([in_DIR cruise '_' station '_acs' '.txt'],'w'); 
% Create Header in the file:
fprintf(fid_acs,'%s\n',headER(1).c); %Line 1
fprintf(fid_acs,'%s\n',[headER(2).c datestr(datenum(date),'yyyymmdd')]); %Line 2
fprintf(fid_acs,'%s\n',[headER(3).c affiliations]); %Line 3
fprintf(fid_acs,'%s\n',[headER(4).c investigators]); %Line 4
fprintf(fid_acs,'%s\n',[headER(5).c contact]); %Line 5
fprintf(fid_acs,'%s\n',[headER(6).c experiment]); %Line 6
fprintf(fid_acs,'%s\n',[headER(7).c cruise '_' dat]); %Line 7
fprintf(fid_acs,'%s\n',[headER(8).c station]); %Line 8
fprintf(fid_acs,'%s\n',headER(9).c); %Line 9
fprintf(fid_acs,'%s\n',[headER(10).c lon '[deg]']); %Line 10
fprintf(fid_acs,'%s\n',[headER(11).c lon '[deg]']); %Line 11
fprintf(fid_acs,'%s\n',[headER(12).c lat '[deg]']); %Line 12
fprintf(fid_acs,'%s\n',[headER(13).c lat '[deg]']); %Line 13
fprintf(fid_acs,'%s\n',[headER(15).c dat]); %Line 14
fprintf(fid_acs,'%s\n',[headER(16).c dat]); %Line 15
fprintf(fid_acs,'%s\n',[headER(17).c startTIME '[GMT]']); %Line 16
fprintf(fid_acs,'%s\n',[headER(18).c tiME{end} '[GMT]']); %Line 17
fprintf(fid_acs,'%s\n',[headER(19).c in_FILE]); %Line 18
fprintf(fid_acs,'%s\n',['/associated_files=' offset_NAME]); %Line 19
fprintf(fid_acs,'%s\n',[headER(20).c 'NA']); %Line 20
fprintf(fid_acs,'%s\n',headER(21).c); %Line 21
if ~strcmpi(doC,'NA');
    fprintf(fid_acs,'%s\n',[headER(22).c ',' doC]); %Line 22
else
    fprintf(fid_acs,'%s\n',headER(22).c); %Line 22
end
fprintf(fid_acs,'%s\n',headER(23).c); %Line 23
fprintf(fid_acs,'%s\n',headER(24).c); %Line 24
fprintf(fid_acs,'%s\n',[headER(25).c D]); %Line 25
fprintf(fid_acs,'%s\n',[headER(26).c cal_FILE_ac]); %Line 26
fprintf(fid_acs,'%s\n',[headER(27).c fields_a fields_c]); %Line 27
fprintf(fid_acs,'%s\n',[headER(28).c units units]); %Line 28
%fprintf(fid_acs,'%s\n',headER(29).c); %Line 29

fprintf(fid_acs,'%s\n','! Summary of ac-s Processing Steps:'); %Line 29

wvl1 = min([floor(a_wv(grate_IND)) floor(c_wv(grate_IND))]);
wvl2 = max([ceil(a_wv(grate_IND+2)) ceil(c_wv(grate_IND+2))]);
fprintf(fid_acs,'%s\n',['!   1. Correct for the spectral "jump" in the center of ac spectra (between ' num2str(wvl1) '-' num2str(wvl2) 'nm) caused by the ac-s using separate detectors for upper and lower wavelengths']); %Line 30
fprintf(fid_acs,'%s\n','!   2. Adjust ac spectra using pure water calibration file'); %Line 31
fprintf(fid_acs,'%s\n','!   3. Adjust ac spectra for temperature and salinity biases using Sullivan et al. (2006)'); %Line 32
fprintf(fid_acs,'%s\n','!   4. Correct a spectra for scattering biases using Rottgers et al. (2013)'); %Line 33
fprintf(fid_acs,'%s\n','!   5. Exclude contaminated ac spectra'); %Line 34

%fprintf(fid_acs,'%s\n',['! The spectral "jump", (referenced in step 1) is corrected by subtracting an offset from ac channels ' num2str(grate_IND+1)  '-' num2str(length(a_wv)) '. Offset is calculated by interpolating ac for channel ' num2str(grate_IND+1) '. Interpolations were performed with matlabs "spline" function - channels ' num2str(grate_IND-2) '-' num2str(grate_IND) ' were used as inputs.']); %Line 35
fprintf(fid_acs,'%s\n',['! The spectral "jump", (referenced in step 1) is corrected by' ...
' subtracting an offset from ac channels ' num2str(grate_IND+1)  '-' ...
num2str(length(a_wv)) '. Offset is calculated by interpolating ac for channel ' ...
num2str(grate_IND+1) '. Interpolations were performed with matlabs "spline" function - channels ' ...
num2str(grate_IND-2) '-' num2str(grate_IND) ' were used as inputs.']); %Line 35
fprintf(fid_acs,'%s\n','! Above-mentioned ac offset values are listed under "/associated_files"'); %Line 36
fprintf(fid_acs,'%s\n','! Contaminated c spectra are defined as spectra for which c is greater than 4 at any point, c is negative between 400-700nm, or corresponding a spectra are contaminated'); %Line 37
fprintf(fid_acs,'%s\n','! Contaminated a spectra are defined as spectra for which a is greater than the maximum c value at any point, a is negative between 400-700nm, or corresponding c spectra are contaminated'); %Line 38
fprintf(fid_acs,'%s\n','! Data are likely to contain multiple downcasts and upcasts. All downcasts and upcasts are included in this submission.'); %Line 39
fprintf(fid_acs,'%s\n','! Data have been processed using code written and made avaiable by Jesse Bausell (email: jbausell@ucsc.edu, GitHub: JesseBausell).'); %Line 40
fprintf(fid_acs,'%s\n','/end_header'); %Line 41

for ii = 1:length(deptH)
    % Print processed ac-s data line by line
    fprintf(fid_acs,'%s\t',tiME{ii}); % time stamp
    fprintf(fid_acs,'%f\t',deptH(ii)); % depth
    fprintf(fid_acs,printFODDER,A_CORR(ii,:)); % absorption spectrum
    fprintf(fid_acs,[printFODDER(1:end-1) 'n'],C_CORR(ii,:)); % attenuation spectrum
end

fclose(fid_acs); % Complete the ac-s .txt file by closing it.


%% 10. Create binned ac-s files (absorption and attenuation are separate).

% Reopen newly made ac-s .txt file to extract the header to extract header
% data
fid_acs2 = fopen([in_DIR cruise '_' station '_acs' '.txt']); % Open combined ac-s file
for ii = 1:41
    % Read the entire 41 line header into a matlab structure
    headER(ii).acs = fgetl(fid); 
end
fclose(fid_acs2); % close the newly made ac-s .txt file

BIN_acDATA; % Creates binned absorption and attenuation files using user-specified bin sizes.