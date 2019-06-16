% BIN_acDATA
% Jesse Bausell
% July 7, 2018
%
% This program takes processed seabass data and depth-bins it into
% user-selected bin sizes. It then builds Seabass submitable files
% containing binned average absorption and attenuation spectra, along with
% their standard deviations. 

        
%% 1. Depth-bin ac-s data according to user-specified bin sizes.
% In this section of code, ac-s data are binned using the different
% binsizes that the user selects.

% Ready the binning by providing the necessary arrays to hold binned data.
% Each of the two cell arrays contains a cell for each binsize. Absorption
% and attenuation data will be split up and placed inside two different
% arrays, with each cell within an holding the binned data plus depth bin
% medians.
binned_CELLS_c = cell(1,length(depth_BINSIZE)); %Holds binned attenuation data
binned_CELLS_a = cell(1,length(depth_BINSIZE)); %Holds binned absorption data
max_DEPTH = ceil(max(deptH)); % Maximum depth of the ac-s 
[v,h] = size(C_CORR); % Dimensions of absorption/attenuation matrix (they are the same size)

for ii = depth_BINSIZE
    % For each user-selected depth binsize create matrices for absorption
    % and attenuation. First column is the depth bins, followed by
    % bin-averaged absorption/attenuation spectra, followed by
    % absorption/attenuation standard deviations.

    % For each cell, create nan matrices for all of the binned mean and
    % standard deviations cell.
    BINNED_MAT_c = nan(ceil(max_DEPTH/ii),2*h+1);
    BINNED_MAT_a = nan(ceil(max_DEPTH/ii),2*h+1); 

    for jj = 0:ii:max_DEPTH-ii
        % This for loop fills the empty nan matrix
        vert_IND = length(0:ii:jj); % vertical index for BINNED_MAT
        BINNED_MAT_c(vert_IND,1) = (jj+jj+ii)/2; %Insert binned depth
        BINNED_MAT_a(vert_IND,1) = (jj+jj+ii)/2; %Insert binned depth
        % The following four lines bin absorption and attenuation
        BIN_IND = find(deptH >=jj & deptH <jj+ii); % Find appropriate depth bin
        BINNED_MAT_c(vert_IND,2:h+1) = nanmean(C_CORR(BIN_IND,:),1); %binned attenuation avg
        BINNED_MAT_c(vert_IND,h+2:2*h+1) = nanstd(C_CORR(BIN_IND,:),0); %binned attenuation std
        BINNED_MAT_a(vert_IND,2:h+1) = nanmean(A_CORR(BIN_IND,:),1); %binned absorption avg
        BINNED_MAT_a(vert_IND,h+2:2*h+1) = nanstd(A_CORR(BIN_IND,:),0); %binned absorption std
    end
    nNAN_IND = find(~isnan(BINNED_MAT_a(:,2))); % Find any NaN rows in matrices
    binned_CELLS_c{ii == depth_BINSIZE} = BINNED_MAT_c(nNAN_IND,:); % Extract NaNs from attenuation
    binned_CELLS_a{ii == depth_BINSIZE} = BINNED_MAT_a(nNAN_IND,:); % Extract NaNs from absorption
end

%% 2. Create headers for binned absorption and attenuation files

for i = 1:25
    % Copy first lines 1-25 of the header from the combined acs file
    headER(i).c = headER(i).acs; % binned attenuation header
    headER(i).a = headER(i).acs; % binned absorption header
end

% line 26 in header lists pure-water calibration files for absorption and
% attenuation.
linE_26 = headER(26).acs; %pull line 26 of combined acs file
equaL = regexpi(linE_26,'='); %find index of = sign
commA = regexpi(cal_FILE_ac,','); % find index of comma separating the two files
headER(26).c = [linE_26(1:equaL) cal_FILE_ac(commA(1)+1:end)]; % Attenuation pure-water file here
headER(26).a = [linE_26(1:equaL) cal_FILE_ac(1:commA(1)-1)]; % Absorption pure-water file here

% line 27 in header lists fields of ac-s data. This must be broken up for
% binned data because absorption and attenuation are formatted separately.
linE_27 = headER(27).acs; % pull line 27 from combined acs file
ENDER_COMMAc = regexpi(linE_27,'cgp'); % Find all attenuation fields in line 27
ENDER_COMMAa = regexpi(linE_27,'agp'); % Find all absorption fields in line 27

% Construct line 27 for absorption and attenuation fields.
linE_27c = ['/fields=depth,' linE_27(ENDER_COMMAc(1):end)]; % Include ONLY attenuation fields
linE_27a = ['/fields=depth,' linE_27(ENDER_COMMAa(1):ENDER_COMMAc(1)-2)]; % Include ONLY absorption fields
commASc = regexpi([linE_27c ','],','); % Find comma indices in line 27 (attenuation) & add comma to the end.
commASa = regexpi([linE_27a ','],','); % Find comma indices in line 27 (absorption) & add comma to the end.
% Fields for absorption and attenuation standard deviations must be added
% onto line 27 for both variables.
SD_c = []; SD_a = []; % Empty repositories to hold strings absorption and attenutation standard deviations
for ii = 1:length(commASa)-1
    % Append absorption and attenuation standard deviation fields onto line
    % 27 one channel (wavelength) at a time. Take each field one at a time,
    % append '_SD' to the end of it, and then place it into the repository
    % variables.
    SD_c = [SD_c ',' linE_27c(commASc(ii)+1:commASc(ii+1)-1) '_SD'];% standard deviations - attenuation channels
    SD_a = [SD_a ',' linE_27a(commASa(ii)+1:commASa(ii+1)-1) '_SD'];% standard deviations - absorption channels
end
headER(27).c = [linE_27c SD_c]; % Combine line 27 with standard deviation fields attenuation
headER(27).a = [linE_27a SD_a]; % Combine line 27 with standard deviation fields absorption

%line 28 in header lists units of measurement for the fields
linE_28 = headER(28).c; %modify this line to fit with binned data files 
headER(28).c = ['/units=' headER(28).acs(17:end)]; % units of measurement for attenuation 
headER(28).a = ['/units=' headER(28).acs(17:end)]; % units of measurement for absorption 

% lines 29-42. These are taken from combined acs file header, except for
% line 35, which is created here. Line 35 changes with bin size so only
% the first part of the line is created here.
for ii = 29:34 
    % Create header lines 29-34
    headER(ii).c = headER(ii).acs;
    headER(ii).a = headER(ii).acs;
end
% Create the first part of header line 35
headER(35).c = '!   6. Calculate spectral means and standard deviations for ';
headER(35).a = '!   6. Calculate spectral means and standard deviations for ';
for ii = 36:42
    % Create header lines 36-42 out of combined acs file header lines 35-41
    headER(ii).c = headER(ii-1).acs;
    headER(ii).a = headER(ii-1).acs;
end

%% 3. Write binned absorption and attenuation .txt files for Seabass submission
% At this point, BIN_acDATA has formulated headers for binned absorption
% and attenuation data, as well as the binned data itself. Now it is time
% to combine headers with binned ac-s data and write the data files in
% Seabass-compatible format.
 
for kk = depth_BINSIZE
    % Create binned data file for each user-specified bin size.
    fid_c = fopen([in_DIR cruise '_' station '_c_bin' num2str(kk)  '.txt'],'w'); %create binned attenuation .txt file
    fid_a = fopen([in_DIR cruise '_' station '_a_bin' num2str(kk)  '.txt'],'w'); %create binned absorption .txt file
    for ll = 1:length(headER) 
    % Create header in binned absorption and attenuation files line by line
        if isequal(ll,35)  
            % line 35 needs to be customized according to bin size
            fprintf(fid_c,'%s\n',[headER(ll).c ' ' num2str(kk) ' m depth bins.']);
            fprintf(fid_a,'%s\n',[headER(ll).a ' ' num2str(kk) ' m depth bins.']);
        else
            % all other lines in the 42 line header do not need to be
            % customized. They can be printed into binned data files as is.
            fprintf(fid_c,'%s\n',headER(ll).c); % attenuation header
            fprintf(fid_a,'%s\n',headER(ll).a); % absorption header
        end
    end
    BINNED_MAT = binned_CELLS_c{kk==depth_BINSIZE}; % Create variable with binned attenuation data matrix
    [l_bin,w_bin] = size(BINNED_MAT); % Get size dimensions of this binned attenuation data matrix
    for mm = 1:l_bin
        % print data into .txt files one line at a time
        fprintf(fid_c,['%1.3f\t' printFODDER printFODDER(1:end-1) 'n'],binned_CELLS_c{kk==depth_BINSIZE}(mm,:)); %attenuation
        fprintf(fid_a,['%1.3f\t' printFODDER printFODDER(1:end-1) 'n'],binned_CELLS_a{kk==depth_BINSIZE}(mm,:)); %absorption
    end
    fclose(fid_c); % close binned attenuation file
    fclose(fid_a); % close binned absorption file
end