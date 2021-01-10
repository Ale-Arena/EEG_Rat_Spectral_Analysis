%%
clear all
close all
clc

%% LOADING FOR TIME-FREQUENCY ANALYSIS
% Loading an EEG or LFP signal that must be organized in a matrix as: 
% samples x channels x trials

load('j15_rec1_Wake_PREPROCESSED_HBP');

%
% DEFINE PARAMETERS (if not already specified)
%

eeg = eeg;      % matrix: samples x channels x trials
time = time;    % time vector of trial (1 x samples)
fs = fs;        % sampling frequency in [Hz]

clear eeg_diff

%% WAVELET CONVOLUTION  - STEP 1: CREATION OF KERNEL (Family of Wavelets)
% A wavelet convolution is perfomed on each EEG/LFP trial in order to
% extract powers and phases for all time-frequency-channel points.
% 1) The kernel of the convolution is created (family of wavelets)
% 2) Convolution is performed
% 3) Power and Phases are extracted for each time-frequency-channel-trial points

%
% CREATION OF A FAMILY OF WAVELETS:
% (40 Morlet wavelets of 3 cycles and 4 seconds, linearly spanning from 1 to 40 Hz)
%

% Definition of parameters:
wavelet_time = -2:1/fs:2; % wavelet lenght
n = 3;                    % number of cycles of wavelets

lowest_frequency  =   1;  % lowest pk [Hz]
highest_frequency =  40;  % highest pk [Hz]
num_wavelets      =  (highest_frequency-lowest_frequency)+1;  % number of pk (frequency bands) 

% create a vector of peak frequencies:
frequencies = linspace(lowest_frequency,highest_frequency,num_wavelets);

% initialize the matrix of wavelet family:
wavelet_family = zeros(num_wavelets,length(wavelet_time));
 
% Loop through frequencies and make a family of wavelets.
for fi=1:num_wavelets
    % create a sine wave at this frequency
    sinewave = exp(2*1i*pi*frequencies(fi).*wavelet_time); % the "1i" makes it a complex wavelet
    % create a Gaussian window
    gaus_win = exp(-wavelet_time.^2./(2*(n/(2*pi*frequencies(fi)))^2));
    % create wavelet via element-by-element multiplication of the sinewave and gaussian window
    wavelet_family(fi,:) = sinewave.*gaus_win;
end
clear fi

%
% PLOT WAVELETS FOR VISUALIZATION AND QUALITY CHECK
% (All the wavelet need to be centered in 0s and 
%  their amplitude have to converge to zero at the time extremities)
%

figure (2)

% Plot wavelets
subplot(2,1,1)
plot(wavelet_time,real(wavelet_family)') 
title('Family of wavelet 1')

% finally, image the wavelet family.
subplot(2,1,2)
imagesc(wavelet_time,frequencies,real(wavelet_family))
axis xy % equivalent to "set(gca,'ydir','normal')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Family of wavelet 2')

axis xy 
xlabel('Time (s)')

%% WAVELET CONVOLUTION  - STEP 2: PERFORM CONVOLUTION
% A wavelet convolution is perfomed on each EEG/LFP trial in order to
% extract powers and phases for all time-frequency-channel points.
% 1) The kernel of the convolution is created (family of wavelets)
% 2) Convolution is performed
% 3) Power and Phases are extracted for each time-frequency-channel-trial points


%
% DEFINE THE DATA TO ANALYSE AND RUN CONVOLUTION
%

prestim = abs(time(1)); % in [s], define the max abs time pre stimulus
X1= -1; % in [s], define the lower time for TIME FREQUENCY ANALYSIS
X2= 2; % in [s], define the upper time for TIME FREQUENCY ANALYSIS

% define subset data for CONVOLUTION
data = eeg (((fs.*X1)+(fs.*prestim)+1):((fs.*X2)+(fs.*prestim)+1),:,:);

% time_data
time_data = X1:1/fs:X2;

% CONVOLUTION: loop through frequencies and compute synchronization
CONV = zeros (size(wavelet_family,1),size(data,1),size(data,3),size(data,2));

for fr = 1:size(wavelet_family,1);
for sw = 1:size(data,3);
for ch = 1:size(data,2);
    
    CONV(fr,:,sw,ch) = conv (squeeze(data(:,ch,sw)),squeeze(wavelet_family(fr,:)), 'same'); 

end 
end
end

clear fr sw ch
clear gaus_win highest_frequency lowest_frequency n  num_wavelets prestim
clear sinewave wavelet_family wavelet_time


%% WAVELET CONVOLUTION  - STEP 3: EXTRACTION OF POWER AND PHASES
% A wavelet convolution is perfomed on each EEG/LFP trial in order to
% extract powers and phases for all time-frequency-channel points.
% 1) The kernel of the convolution is created (family of wavelets)
% 2) Convolution is performed
% 3) Power and Phases are extracted for each time-frequency-channel-trial points


%
% EXTRACT POWERS and NORMALIZATION:
%

% Raw power extraction
powers = (abs(CONV(:,:,:,:)).^2); 
% Matrix: frequency x samples x trials x channels

% Power Normalization: from uV^2 to % of baseline to dB:
% Obtain Mean Baseline for each single frequency band
b1 = -.5;    % [s] time range of baseline
b2 = -.2;
% Mean baseline across time for each trial, channel and frequency band:
power_baseline = squeeze(mean (powers(:,((fs.*b1)+(fs.*abs(X1))+1):((fs.*b2)+(fs.*abs(X1))+1),:,:),2));
% Mean baseline across time and trials for each frequency band and channel:
power_baseline = squeeze(mean (power_baseline,2));

% Normalization on the mean baseline across time and trials  
powers_NORM = zeros(size(powers,1),size(powers,2),size(powers,3),size(powers,4));

for fr = 1:size(powers,1); % frequencies
for sw = 1:size(powers,3); % sweeps
for ch = 1:size(powers,4); % channels
    
  % Normalization of each trial on the mean baseline across trials and time for each freq. band and channel
    powers_NORM (fr,:,sw,ch) =  powers(fr,:,sw,ch)./power_baseline(fr,ch);

end
end
end
clear fr sw ch 
clear power_baseline

% Mean Normalized Power across trials
powers_M_NORM = squeeze(mean (powers_NORM, 3));

% Power Conversion to dB
decibels =  10.* (log10(powers_NORM));
decibels_M =  10.* (log10(powers_M_NORM));



%
% PHASE EXTRACTION and CLUSTERING COMPUTATION:
%

% Raw Phase Extraction
phases = angle(CONV(:,:,:,:));    % exctract phases from convolution

% Compute the average vector of phases of each trials at all time points and peak frequencies
% by Euler's formula (vector lengths set to 1)
% NB: av_vector_phase = average phase-angle among trials [radiants]
%     av_vector_length= magnitude of phase-locking across trials [between 0 and 1]

% av_vector_phase = squeeze (angle(mean(exp(1i.*(phases(:,:,:,:))),3))); 
av_vector_length = squeeze (abs(mean(exp(1i.*(phases(:,:,:,:))),3))); % Inter Trial Phase Clustering



%%  BOOTSTRAP STATISTIC 1 - COMPUTE STATISTICAL THRESHOLD FOR dB
% It finds 1 positive and 1 negative thresholds based on the surrogate distribution of the maximum and
% minimum normalized powers of the ensamble average across trials (resampled) 
% for each frequency and channel

%
% DEFINE PARAMETERS:
%

db1 = b1;    %[s] time range of baseline
db2 = b2;
baseline = powers_NORM (:,((fs.*db1)+(fs.*abs(X1))+1):((fs.*db2)+(fs.*abs(X1))+1),:,:); % baseline matrix from which it calculates the distribution (Matrix: frequency x samples x trials x channels)
n_straps = 500; % total number of resamplings (bootstraps; usually between 500 and 1000 is enough)
alpha = 0.05;   % significant level (probability of rejecting the null hypothesis when it's true) can be 0,05 or less
n_samples = size (baseline,2); % number of samples
n_sweeps = size (baseline,3);  % number of sweeps
n_freq = size(baseline,1);     % number of frequencies
n_channel = size(baseline,4);  % number of channels

%
% LOOP OVER BOOTSTRAP SAMPLES AND CREATE SURROGATE DISTRIBUTIONS 
%

bootstrap_dB_max = zeros(n_freq,n_channel,n_straps);
bootstrap_dB_min = zeros(n_freq,n_channel,n_straps);

for ii = 1:n_straps
    
        % Take randomly (n = n_sweeps) samples from data (indexes of sweeps, 1:nsweeps), 
        % with replacement:
        resampled_sweeps = datasample(1:n_sweeps,n_sweeps); 
        
        % Take the actual trials corresponding to the randomized indexes:
        resampled_baseline_dB = baseline(:,:,resampled_sweeps,:);
         
        % Make the ensamble average of the new resempled baseline matrix and dB conversion:
        resampled_averageBaseline_dB = squeeze(mean(resampled_baseline_dB,3));
        resampled_averageBaseline_dB = 10.* (log10(resampled_averageBaseline_dB));
   
        bootstrap_dB_max(:,:,ii) = max(resampled_averageBaseline_dB,[],2);
        bootstrap_dB_min(:,:,ii) = min(resampled_averageBaseline_dB,[],2);

end

%
% FIND THRESHOLDS
%

T_max_percentile_dB = zeros(n_freq,n_channel); % Positive threshold
T_min_percentile_dB = zeros(n_freq,n_channel); % Negative threshold

for ch = 1:n_channel
    for fr = 1:n_freq
    T_max_percentile_dB(fr,ch) = prctile(bootstrap_dB_max(fr,ch,:),100.*(1-alpha));
    T_min_percentile_dB(fr,ch) = prctile(bootstrap_dB_min(fr,ch,:),100.*alpha);
    end
end
clear ii ch fr

clear baseline n_straps alpha n_samples n_sweeps n_freq n_channel
clear resampled_averageBaseline_dB  resampled_baseline_dB resampled_sweeps
clear bootstrap_dB_max bootstrap_dB_min  db1 db2


%%  BOOTSTRAP STATISTIC 2 - COMPUTE STATISTICAL THRESHOLD FOR ITPC
% It finds 1 positive threshold based on the surrogate distribution 
% of the maximum Inter Trial Phase Clustering across trials (resampled) 
% for each frequency and channel

%
% DEFINE PARAMETERS:
%

p1 = b1;    % [s] time range of baseline
p2 = b2;

% Take the phases of baseline
Phases_baseline = phases(:,((fs.*p1)+(fs.*abs(X1))+1):((fs.*p2)+(fs.*abs(X1))+1),:,:); 

% Euler's Formula:
Euler_baseline = (exp(1i*(Phases_baseline(:,:,:,:)))); % Partial Inter Trial Phase Clustering formula

% Define Boostrap Parameters:
baseline = Euler_baseline; % baseline matrix from which it calculates the surrogate distribution (Matrix: frequency x samples x trials x channels)
n_straps = 500; % total number of resamplings (bootstraps; usually between 500 and 1000 is enough)
alpha = 0.01;   % significant level (probability of rejecting the null hypothesis when it's true) can be 0,05 or less
n_sweeps = size (baseline,3);  % number of sweeps
n_freq = size(baseline,1);     % number of frequencies
n_chii = size(baseline,4);  % number of channels

%
% LOOP OVER BOOTSTRAP SAMPLES AND CREATE SURROGATE DISTRIBUTIONS 
%

bootstrap_ITPC_max = zeros(n_freq,n_chii,n_straps);

for bb =   1:n_straps
    
        % Take randomly (n = n_sweeps) samples from data (indexes of sweeps, 1:nsweeps), 
        % with replacement:
        resampled_sweeps = datasample(1:n_sweeps,n_sweeps); 
        
        % Take the actual sweeps corresponding to the randomized indexes:
        resampled_baseline = baseline(:,:,resampled_sweeps,:);
        
        % Make the Inter Trial Phase Clustering (ITPC) of the new resempled baseline matrix:
        resampled_ITPC_baseline = squeeze(abs(mean(resampled_baseline,3)));
   
       
        % Find the maximum value of the new ITPC of
        % running baseline and create a matrix with MAX for each
        % frequecy and channel:
        bootstrap_ITPC_max(:,:,bb) = max(resampled_ITPC_baseline,[],2);
       
end
clear bb

%
% FIND THRESHOLDS
%

T_max_percentile_ITPC = zeros(n_freq, n_chii);
for chii = 1:n_chii
for fr = 1:n_freq
    
    T_max_percentile_ITPC(fr,chii) = prctile(bootstrap_ITPC_max(fr,chii,:),100.*(1-alpha));

end
end
clear ii chii chjj fr


clear baseline n_straps alpha n_samples n_sweeps n_freq n_channel
clear resampled_averageBaseline_ITPC  resampled_baseline_ITPC resampled_sweeps
clear bootstrap_ITPC_max_LF  bootstrap_ITPC_max_HF bootstrap_ITPC_max
clear resampled_ITPC_baseline_HF resampled_ITPC_baseline_LF
clear fH_start fH_x1 fH_end fH_x2 fL_start fL_x1 fL_end fL_x2 
clear Euler_baseline n_chii Phases_baseline resampled_ITPC_baseline resampled_baseline
clear p1 p2

%%  BOOTSTRAP STATISTIC 3 - THRESHOLDING ON dB and ITPC
% It performs thresholding on both powers (dB signal) and phase-locking (ITPC signal)
% in order to conserve only significant variations from baseline

%
% TAKE THE SIGNIFICANT MEAN RELATIVE POWERS (dB)
%

% Set not significant values to zero
% For decibels:
decibels_M_sig = zeros(size(decibels_M,1), size(decibels_M,2), size(decibels_M,3));

for fr = 1:size(decibels_M,1);
for sp = 1:size(decibels_M,2);
for ch = 1:size(decibels_M,3);
   
    if decibels_M (fr, sp, ch) > T_max_percentile_dB(fr,ch);  
       decibels_M_sig (fr, sp, ch) = decibels_M (fr, sp, ch); 
    else
        if decibels_M (fr, sp, ch) < T_min_percentile_dB(fr,ch);  
           decibels_M_sig (fr, sp, ch) = decibels_M (fr, sp, ch); 
        else
            decibels_M_sig (fr, sp, ch) = 0;
        end
    end  

end
end
end
clear fr sp ch


%
% TAKE THE SIGNIFICANT ITPC 
%

    av_vector_length_sig = zeros(size(av_vector_length,1), size(av_vector_length,2), size(av_vector_length,3));

    for fr = 1:size(av_vector_length,1);
    for sp = 1:size(av_vector_length,2);
    for ch = 1:size(av_vector_length,3);
   
     if av_vector_length (fr, sp, ch) > T_max_percentile_ITPC(fr,ch) ;  
        av_vector_length_sig (fr, sp, ch) = av_vector_length (fr, sp, ch);
     else
        av_vector_length_sig (fr, sp, ch) = 0; 
     end  

    end
    end
    end

    clear fr sp ch
    clear phases powers powers_NORM decibels CONV
    
%% DETECTION OF THE INTERUPTION OF PHASE-LOCKING (ITPC_drop time)
% It finds the last significant ITPC in a certain time and frequency range

%
% PARAMETERS 
%

% Post stimulus time window to use:
prestim = abs(time_data(1)); %[s]
ITCP_post_x1 = 0; %[s]
ITCP_post_x2 = 0.8; %[s]

% Frequency Range to use:
fH_start = 8; % [Hz] bottom limit of High Frequency range
fH_x1 = find(frequencies == fH_start);
fH_end = 40;   % [Hz] upper limit
fH_x2 = find(frequencies == fH_end);
   

% Extract the poststimulus matrix of the significant ITCP (average vector
% lengths) for all frequencies and channels: 
    ITCP_postim_H_sig = av_vector_length_sig (fH_x1:fH_x2,((fs.*ITCP_post_x1)+(fs.*prestim)+1):((fs.*ITCP_post_x2)+(fs.*prestim)+1),:);
% Average across frequencies:    
    ITCP_postim_av_H = squeeze(mean (ITCP_postim_H_sig,1));
% Create its time vector:
    time_postim = ITCP_post_x1:1/fs:ITCP_post_x2;    
        
%
% ITPC_drop Time FOR EACH CHANNEL
%

% Find the sample point corrisponding to the last synchronous time point
% whithin the defined frequency range for each channel:
    ITCP_index_lastsig = zeros(1, size (ITCP_postim_av_H,2));
    
    for ch = 1:size(ITCP_postim_av_H,2)
 
    % IF no phase clustering between 0 and .8s (end) -> No syncrony at all -> set to first sample (0s)
    if mean(ITCP_postim_av_H(:,ch), 1) == 0 
           ITCP_index_lastsig(:,ch) = find(time_postim==ITCP_post_x1);
    else
    % FIND the last significant ITPC in the range    
    ITCP_index_lastsig(:,ch) = find(ITCP_postim_av_H(:,ch),1,'last');
    
    end
    end
    clear ch
   
% Convert indexes in seconds:
ITPC_drop = time_postim(ITCP_index_lastsig);
ITPC_drop_M = squeeze(mean(ITPC_drop,2)); % mean across channels

ITPC_freq = [fH_start, fH_end];


clear fH_x1 fH_x2 fH_end fH_start ITCP_index_lastsig ITCP_post_x1 ITCP_post_x2
clear ITCP_postim_av_H ITCP_postim_H_sig av_vector_length time_postim

   
   
%% QUANTIFICATION OF EARLY HIGH FREQUENCY POWER (OFF PERIOD, if any)
% It quantifies the High Frequency power dynamic early after the stimulation
% and detects the OFF period (suppression of high frequencies) if present

%
% MEAN HF POWER AFTER STIMULATION (OFF period, if HF power < 0dB)
% NB: Obtained from each channel and then averaged
%

% Define Time windows [ms]
x1 = 0.08;  
x2 = 0.18;

% Define the frequency band of interest
fH_start = 20; % [Hz] bottom limit
fH_x1 = find(frequencies == fH_start);
fH_end = 40;  % [Hz] upper limit
fH_x2 = find(frequencies == fH_end);


% NB: conversion into % of baseline
dB_HF_allCH = squeeze (mean(10.^((decibels_M_sig((fH_x1:fH_x2),:,:)./10)),1));

OFF_dB_p = mean(dB_HF_allCH (((fs.* x1)+(fs.*prestim)+1):((fs.* x2)+(fs.*prestim)+1), :),1);
OFF_dB = 10.*(log10(OFF_dB_p)); % NB: reconverted in dB
OFF_dB_M = 10.*(log10(mean (OFF_dB_p))); % NB: reconverted in dB

OFF_frange = [fH_start, fH_end]; 

clear x1 x2 fH_start fH_x1 fH_end fH_x2 OFF_dB_p HF_Wind1_dB_allCH

%
% QUANTIFICATION OF:
% SUPPRESSION MIN (in High Frequency Range) - Minimum Peak
% SUPPRESSION end and start (in High Frequency Range)  - Starting and ending point of OFF period
%

% NB: Obtained from each channel and then averaged

% Define time window for High Frequency Range:
x1_H = 0.0;  % [s]
x2_H = 0.3;

% HIG FREQ range;

% Extract matrix
T = time_data(:,((fs.* x1_H)+(fs.*prestim)+1):((fs.* x2_H)+(fs.*prestim)+1));
% NB: reconverted in dB
P = 10.*(log10(dB_HF_allCH (((fs.* x1_H)+(fs.*prestim)+1):((fs.* x2_H)+(fs.*prestim)+1),:)));
% Percentage 
Perc = dB_HF_allCH (((fs.* x1_H)+(fs.*prestim)+1):((fs.* x2_H)+(fs.*prestim)+1),:);

% Obtain minimum power for all channels and its time:
[MIN_power,MIN_index] = min(P,[],1);
MIN_time = T(MIN_index);

% Obtain the mean minimum power and time among channels:
MIN_power_mean = 10.*(log10(mean (10.^((MIN_power)./10))));
MIN_time_mean = mean (MIN_time);

% Take only the suppression points (negative voltage)
% Find only the beginning and the end points of the suppression period
P_logic = P < 0;
P_logic_cross = zeros(size (P,1) ,size (P,2));
   
    for kk = 1:size (P,1)-1;
    for ch = 1:size (P,2);

     % 1st_ set +1 = upward crossings; set -1 = downward crossings
    P_logic_cross(kk,ch) = abs(P_logic(kk,ch) - P_logic(kk+1,ch));
    
    % if suppression last longer then the considered window (0,3s), set the
    % end of suppression at the lenght of the window
    if P_logic(end,ch) == 1
       P_logic_cross(end,ch) = 1;
    end
    
    end
    end
 clear kk ch 
 
% Find the indexes of the beginning and the ending points
SUPP_index_start = zeros (1,size (P,2));
SUPP_index_end = zeros (1,size (P,2));

for ch = 1:size(P,2);
   
    % find the first crossing after the suppression peak 'SUPP_index_end'
    % and the absolute first 'SUPP_index_start'
    % if no suppression, set index of start and end to 1 (time = 0)
    
    E = P_logic_cross(MIN_index(ch):size(P_logic_cross,1),ch);
    if mean (E) == 0;
      SUPP_index_end(ch) = 1;
    else
    SUPP_index_end(ch)  = find (E,1, 'first');
    SUPP_index_end(ch) = SUPP_index_end(ch)+MIN_index(ch)-1;
    end
    
    S = P_logic_cross(1:MIN_index(ch),ch);
    if mean (S) == 0;
      SUPP_index_start(ch) = 1;
    else
    SUPP_index_start(ch)  = find (S,1, 'first') + 1;
    end
    
end
clear ch

% Obtain the time points of the START and END of Suppression 
SUPP_end = zeros (1,size (P,2));
SUPP_start= zeros (1,size (P,2));

for ch = 1:size(P,2);
SUPP_end(ch) = T(SUPP_index_end(ch));
SUPP_start(ch) = T(SUPP_index_start(ch));
end
clear ch

% Obtain the mean start and end time of suppression (only across channels that present suppression)
SUPP_endm = mean(SUPP_end(SUPP_end>0),2);
SUPP_startm = mean(SUPP_start(SUPP_start>0),2);

% Obtain the suppression duration and the mean duration (only across channels that present suppression):
SUPP_duration = SUPP_end - SUPP_start;
SUPP_mduration = mean(SUPP_duration(SUPP_duration>0),2);

% Find the channels without suppression:
B = SUPP_end == 0;
SUPP_free_ch = find (B);
clear B 

%
% QUANTIFICATION OF AMPLITUDE OF EARLY EVOKED RESPONSE (root mean squared)
%

% Obtain the RMS of the early ERP:
x_rms_1 = 0.006; % [s]
x_rms_2 = 0.05;

early_ERP = squeeze(mean(data (((fs.* x_rms_1)+(fs.*prestim)+1):((fs.* x_rms_2)+(fs.*prestim)+1), :,:),3));
RMS_early_ERP = rms(early_ERP,1);
RMSm_early_ERP = mean(RMS_early_ERP,2);

clear x_rms_1 x_rms_2 early_ERP SUPP_index_end SUPP_index_start
clear S E P_logic_cross P_logic dB_HF_allCH MIN_index P Perc T x1_H x2_H


%% LATE HF POWER QUANTIFICATION
%
% FIND THE TIME POINT OF FIRST HFpower AFTER SUPPRESSION for each channel
% FIND ALSO THE TIME POINT OF THE CORRESPONDING ITPC
% MAKE AVERAGE ACROSS THOSE CHANNELS
% CALCULATE PERSENTAGE OF CHANNELS THAT HAVE HFpower AFTER SUPPRESSION
%

% Choose windows:
t1= 0.1; % time interval [s]
t2= 0.8;
fdb1 = 20; % frequency window [Hz]
fdb2 = 40; 

T = time_data(:,((fs.* t1)+(fs.*prestim)+1):((fs.* t2)+(fs.*prestim)+1));
power_wind = decibels_M_sig (:,((fs.* t1)+(fs.*prestim)+1):((fs.* t2)+(fs.*prestim)+1),:);
% NB: conversion into % of baseline
power_wind_NORM = squeeze (mean(10.^((power_wind((fdb1:fdb2),:,:)./10)),1));
% NB: reconverted in dB
power_wind_dB = 10.*(log10(power_wind_NORM));


power_wind_logic = power_wind_dB > 0;

 power_wind_index = zeros(1, size(power_wind_logic,2));
 for ii=1:size(power_wind_logic,2)
     if mean (power_wind_logic(:,ii),1) == 0
         power_wind_index(ii) = 0;
     else
    power_wind_index(ii) = find (power_wind_logic(:,ii), 1, 'first');
     end
end
clear ii

power_channels_aaa = power_wind_index > 0;
power_channels = find(power_channels_aaa); % channels with significant HFpower

power_wind_time_aaa = power_wind_index(power_channels);
power_wind_time = T(power_wind_time_aaa); % onset time of significant HFpower

power_channel_ITPC = ITPC_drop(power_channels); % ITPCdrop of channels with significant HFpower

% Results:
Latepower_time_mean = mean(power_wind_time,2);% mean time of late HFpower onset
Latepower_time = power_wind_time;

Latepower_ITPC_mean = mean(power_channel_ITPC,2); % mean ITPC drop of channel with late HFpower 
Latepower_ITPC = power_channel_ITPC;

Latepower_chpersent = size(power_channels,2)/size(data,2); % persentage of channels with  late HF power

clear power_wind_time_aaa power_channels power_channels_aaa
clear power_wind_index power_wind_logic power_wind_dB power_wind_NORM T power_wind
clear t1 t2 fH1 fH2 fdb1 fdb2 power_channel_ITPC power_wind_time 


%%  FIGURE
% It plots: 
% _ [TOP LEFT,RIGHT] The mean ERP from 1 channel with 5 single trials overimposed, 
% plus the same channel in bold, in the butterfly plot with the mean ERPs from all other channels.
% _ [SECOND LINE LEFT,RIGHT] The mean power spectrum [dB] from one channel and the average across a
% specified frequency range for all the channels
% _ [THIRD LINE LEFT,RIGHT] The phase-locking across trials [ITPC] from one channel and the average across a
% specified frequency range for all the channels
% _ [BOTTOM LINE] The phase-locking across trials [ITPC] from one single channel, averaged across a
% specified frequency range 
%

clf
figure (1)

%
% Single channel
% PARAMETERS:
%

x1 = -.2; % [s]  % Time
x2 =  .6; % [s]
y1 = -200;% [uV] % ERP amplitude
y2 =  200;% [uV]
s1 = -8;  % [dB] % dB in surface plot
s2 = +8;  % [dB]

ch =  12; % number of channel

sw1 = 5;  % number of sweep

f1 = 1;   % frequency to plot (start)
f2 = 40;  % frequency to plot (end)

% bottom frequency limit for average ITPC (start) [Hz]
f_start = ITPC_freq (1); 
% f_start = freqency_range_LOW (1); 
fi_logic_start = frequencies==f_start;
fi_index_start = find(fi_logic_start);
clear fi_logic_start

% upper frequency limit for average ITPC (end)  [Hz]
f_end =  ITPC_freq(2);  
% f_end = freqency_range_LOW (2); 
fi_logic_end = frequencies==f_end;
fi_index_end = find(fi_logic_end);
clear fi_logic_end


%
% PLOT SINGLE CHANNEL:
%

% PLOT ERPs, ENSAMBLE AVERAGE AND SINGLE SWEEPS:
subplot(4,2,1)
% Wakefulness:
plot (time_data, data (:, ch, sw1), 'LineWidth',0.5,'Color',[0.5 0.5 0.5])
hold on
plot (time_data, data (:, ch, sw1+1), 'LineWidth',0.5,'Color',[0.5 0.5 0.5])
hold on
plot (time_data, data (:, ch, sw1+2), 'LineWidth',0.5,'Color',[0.5 0.5 0.5])
hold on
plot (time_data, data (:, ch, sw1+3), 'LineWidth',0.5,'Color',[0.5 0.5 0.5])
hold on
plot (time_data, data (:, ch, sw1+4), 'LineWidth',0.5,'Color',[0.5 0.5 0.5])
hold on
plot(time_data,mean(data(:,ch,:),3), 'LineWidth',1,'Color',[0 0 0])

colorbar
xlabel('Time [s]'), ylabel('Amplitude [uV]')
axis([x1,x2, y1,y2])
title( ['Mean ERP and 5 single ERPs - ch' num2str(ch) '' ])
hold on % Last ITCP_drop
      plot (ITPC_drop(ch),linspace(y1,y2), '. k')

      
% PLOT SIGNIFICANT dB AT EACH FREQUENCY
subplot(4,2,3)
contourf(time_data,frequencies,decibels_M_sig(:,:,ch),100,'linecolor','none')
colorbar
colormap (parula)
axis([x1,x2, f1,f2])
set(gca,'clim',[s1 s2])
xlabel('Time [s]'), ylabel('Frequency [Hz]')
title(['dB Power - ch' num2str(ch) ''])
hold on % Last ITCP_drop
      plot (ITPC_drop(ch),linspace(1,80), '. k')
hold on % Line at 20 Hz
      plot (linspace(x1,x2),20, '. k')

     
% PLOT ITPC OF ALL FREQUENCY BANDS IN TIME:
subplot(4,2,5)
contourf(time_data,frequencies,av_vector_length_sig(:,:,ch),100,'linecolor','none')
colorbar
colormap (parula)
axis([x1,x2, f1,f2])
set(gca,'clim',[0 .8])
xlabel('Time [s]'), ylabel('Frequencies [Hz]')
title(['ITPC - ch' num2str(ch) ''])
hold on % Last ITCP_drop
      plot (ITPC_drop(ch),linspace(1,80), '. k')

    
% PLOT mean ITPC:
subplot(4,2,7)
plot (time_data,squeeze(mean(av_vector_length_sig(fi_index_start:fi_index_end,:,ch),1)),'LineWidth',1.5,'Color',[0 0 0])
colorbar
axis([x1,x2, 0,1])
xlabel('Time [s]')
ylabel('ITPC')
title(['ITPC beetween ' num2str(fi_index_start) ' & ' num2str(fi_index_end) ' Hz - ch' num2str(ch) ''])
hold on % Last ITCP_drop
      plot (ITPC_drop(ch),linspace(0,1), '. k')   

%      
% All channels   
% PARAMETERS:
%

z1 = -5; % [dB]
z2 = +5; % [dB]
itpc1 = 0; % ITPC
itpc2 = .5;

MEDIA1 = squeeze(mean(data,3));

% Define start and end of the range:
    f_start_ITCP = ITPC_freq(1);
    f_end_ITCP = ITPC_freq(2);
    
    f_start_dB = OFF_frange(1);
    f_end_dB = OFF_frange(2);
    
% bottom frequency limit for average ITPC (start) [Hz]
fi_logic_start_ITCP = frequencies==f_start_ITCP;
fi_index_start_ITCP = find(fi_logic_start_ITCP);

% upper frequency limit for average ITPC (end)  [Hz]
fi_logic_end_ITCP = frequencies==f_end_ITCP;
fi_index_end_ITCP = find(fi_logic_end_ITCP);

% bottom frequency limit for average dB (start) [Hz]
fi_logic_start_dB = frequencies==f_start_dB;
fi_index_start_dB = find(fi_logic_start_dB);

% upper frequency limit for average dB (end)  [Hz]
fi_logic_end_dB = frequencies==f_end_dB;
fi_index_end_dB = find(fi_logic_end_dB);


% Obtain the average dB across frequencies:
av_dB_w = squeeze (10.*(log10(mean(10.^((decibels_M_sig(fi_index_start_dB:fi_index_end_dB,:,:))./10),1))));

% Obtain the average ITPC across frequencies:
av_vector_length_M_0w = squeeze(mean(av_vector_length_sig(fi_index_start_ITCP:fi_index_end_ITCP,:,:),1));

%
% PLOT ALL CHANNELS
%

% Plot BUTTERFLY PLOT
subplot (4,2,2);
    plot (time_data, MEDIA1 (:, :), 'LineWidth',0.5,'Color',[0.5 0.5 0.5])
    colorbar
    axis([x1,x2, y1,y2])
    xlabel ('Time [s]')
    ylabel ('Amplitude [uV]')
    title (['Mean ERPs - all channels (ch' num2str(ch) ' in bold)'])
hold on
    plot (time_data, MEDIA1 (:, ch), 'LineWidth',1,'Color',[0 0 0])
    axis([x1,x2, y1,y2])


% PLOT dB FOR EACH CHANNEL IN TIME:
subplot (4,2,4)
contourf(time_data,1:size(av_dB_w,2),av_dB_w(:,:)',80,'linecolor','none')
colorbar
axis([x1,x2, 1,size(decibels_M_sig,3)])
set(gca,'clim',[z1 z2])
xlabel('Time [s]'), ylabel('Channel')
title(['dB Power beetween ' num2str(f_start_dB) ' & ' num2str(f_end_dB) '  Hz - all channels'])

% PLOT ITPC FOR EACH CHANNEL IN TIME:
subplot(4,2,6)
contourf(time_data,1:size(av_vector_length_M_0w,2),av_vector_length_M_0w(:,:)',80,'linecolor','none')
colorbar
axis([x1,x2, 1,size(av_vector_length_sig,3)])
colormap (parula)
set(gca,'clim',[itpc1 itpc2])
xlabel('Time [s]'), ylabel('Channel')
title(['ITPC beetween ' num2str(f_start_ITCP) ' & ' num2str(f_end_ITCP) '  Hz - all channels'])


clear z1 z2 itpc1 itpc2 MEDIA1 f_start_ITCP  f_end_ITCP  f_start_dB   f_end_dB 
clear fi_logic_start_ITCP fi_index_start_ITCP fi_logic_end_ITCP fi_index_end_ITCP 
clear fi_logic_start_dB fi_index_start_dB fi_logic_end_dB fi_index_end_dB 
clear av_dB_w av_vector_length_M_0w x1 x2 y1 y2 s1 s2 ch sw1 f1 f2 


%% SAVING

save ('j07_r1_Light_Sevo1-itpc005_ERPCAUSALITY');


