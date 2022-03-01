%Code for matlab
%https://github.com/mikexcohen/AnalyzingNeuralTimeSeries/blob/main/chapter26.m

% Global variables
global EEG_srate
global EEG_pnts
global EEG_times
EEG_srate = 1000
EEG_pnts = 4001
EEG_times = linspace(-2000,2000,EEG_pnts)
% clear all to clear variables global

% Select All subjects or Individual
subjects_N = true; subjects_1 = false
% N subjects
if subjects_N == true
    disp('Loading mat file N Subjects \\')
    myfullname = '/Users/borri/github/EEG_Cadiz_2021/LFP_EEG.mat'
    EEG = load(myfullname)
    print('EEG.LFP_EEG.data.covert_cond1.[STN,FRONTAL,OCCIPITAL] || EEG.LFP_EEG.data.overt_cond2.[STN,FRONTAL,OCCIPITAL]')
    EEG = EEG.LFP_EEG

    msg = 'N_subjects'
    result = analysis_Nsubjects(EEG, msg)
else
    % Select one subject
    disp('Loading mat file 1 Subject \\')
    myfullname = '/Users/borri/github/EEG_Cadiz_2021/subj5.mat'
    EEG = load(myfullname)
    EEG = EEG.subj5
    msg = 'subject5'
    result = analysis1subjects(EEG, msg)
end 

function result = analysis_1subjects(EEG, msg)
%This function computes ISPC over trial for 1 subject
disp(['File ' myfullname ' loaded!' ])
% Channel origin and target
ntrial = 2
condition = 'covert_cond1'
region1 = 'STN' 
region2 = 'FRONTAL'
trialid = ntrial
cha_tot1 = size(EEG.(condition).(region1){trialid}(:,:))
cha_tot2 = size(EEG.(condition).(region2){trialid}(:,:))
chaid1 = cha_tot1(1) ; chaid2 = cha_tot2(1) 
list_cha1 = 1:(chaid1);list_cha2 = 1:(chaid2)
%cha_data1 = EEG.(condition).(region1){trialid}(chaid,:)
%cha_data2 = EEG.(condition).(region){trialid}(chaid+1,:)
cha_data1 = mean(EEG.(condition).(region1){trialid}(list_cha1,:))
cha_data2 = mean(EEG.(condition).(region2){trialid}(list_cha2,:))

EEG_trials1 = size(EEG.(condition).(region1),2)
EEG_trials2 = size(EEG.(condition).(region2),2)
EEG_trials = min(EEG_trials1,EEG_trials2)


%%%%
% Figure 26.1 phase angle differences between 2 electrodes over time
%%%%
center_freq = 5;
time        = -1:1/EEG_srate:1; % time for wavelet
wavelet     = exp(2*1i*pi*center_freq.*time) .* exp(-time.^2./(2*(4/(2*pi*center_freq))^2))/center_freq;
% FFT parameters
n_wavelet     = length(time);
n_data        = EEG_pnts;
n_convolution = n_wavelet+n_data-1;

% FFT of wavelet
fft_wavelet = fft(wavelet,n_convolution);

% initialize output time-frequency data
phase_data = zeros(2,EEG_pnts);
real_data  = zeros(2,EEG_pnts);


end 


function  result = analysisNsubjects(EEG, msg)
%This function computes ISPC over trial for N subjects
disp(['File ' myfullname ' loaded!' ])

end 










%%% Figs 26.1 and 26.2
%%% 1 trial 1 subject

%EEG.covert_cond1.STN{n_trial}
n_trials_cond = size(EEG.covert_cond1.STN,2)
n_channels = size(EEG.covert_cond1.STN{2},1)
n_points = size(EEG.covert_cond1.STN{2},2); %==4001
% Set a trial we are interested in
ntrial = 2
% Select channeks for one trial and same condition same region
% EEG.<condition>.<region>{<trialid>}(<chanekid>,:)
% condition= covert_cond1|overt_cond2
% regions = OCCIPITAL|FRONTAL|STN
% trialid = ~115 or size(EEG.overt_cond2.OCCIPITAL,2)
% access structure with variables
condition = 'covert_cond1'
region = 'FRONTAL'
trialid = ntrial
chaid = 1
cha_data1 = EEG.(condition).(region){trialid}(chaid,:)
cha_data2 = EEG.(condition).(region){trialid}(chaid+1,:)
EEG_trials = size(EEG.(condition).(region),2)
%cha_data1 = EEG.covert_cond1.STN{ntrial}(1,:)
%cha_data2 = EEG.covert_cond1.STN{ntrial}(2,:)


%ndata = 4000
%n_wavelet = length(time);
%n_conv    = n_wavelet+ndata-1;
%EEGdata 0..4000 float
%min_freq = 2; max_freq = 64; num_frex = 30
%fft_data = fft(squeeze(es{1}(1,:)),n_conv);

%%%%
% Figure 26.1 phase angle differences between 2 electrodes over time
%%%%
center_freq = 5;
time        = -1:1/EEG_srate:1; % time for wavelet
wavelet     = exp(2*1i*pi*center_freq.*time) .* exp(-time.^2./(2*(4/(2*pi*center_freq))^2))/center_freq;
% FFT parameters
n_wavelet     = length(time);
n_data        = EEG_pnts;
n_convolution = n_wavelet+n_data-1;

% FFT of wavelet
fft_wavelet = fft(wavelet,n_convolution);

% initialize output time-frequency data
phase_data = zeros(2,EEG_pnts);
real_data  = zeros(2,EEG_pnts);

for chani=1:2
    %fft_data = fft(squeeze(EEG.data(chanidx(chani),:,1)),n_convolution);
	fft_data = fft(squeeze(cha_data+chani),n_convolution)
	%if chani==1
	%	fft_data = fft(squeeze(cha_data1),n_convolution)
	%else
	%	fft_data = fft(squeeze(cha_data2),n_convolution)
    convolution_result_fft = ifft(fft_wavelet.*fft_data,n_convolution) * sqrt(4/(2*pi*center_freq));
    convolution_result_fft = convolution_result_fft(half_of_wavelet_size+1:end-half_of_wavelet_size);
 
    % collect real and phase data
    phase_data(chani,:) = angle(convolution_result_fft);
    real_data(chani,:)  = real(convolution_result_fft);
end
 % open and name figure
%figure, set(gcf,'Name','Movie magic minimizes the mystery.','Number','off');
timemin=0 ; timemax=4000;
figure, set(gca,'xlim',[timemin, timemax],'ylim',[min(real_data(:)) max(real_data(:))])
% draw the filtered signals
subplot(321)
filterplotH1 = plot(EEG.times(1),real_data(1,:),'b');
hold on
filterplotH2 = plot(EEG.times(1),real_data(2,:),'m');
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
title([ 'Filtered signal at ' num2str(center_freq) ' Hz' ])

% draw the phase angle time series
subplot(322)
phaseanglesH1 = plot(phase_data(1,:),'b');
hold on
phaseanglesH2 = plot(phase_data(2,:),'m');
set(gca,'xlim',[timemin timemax],'ylim',[-pi pi]*1.1,'ytick',-pi:pi/2:pi)
xlabel('Time (ms)')
ylabel('Phase angle (radian)')
title('Phase angle time series')

% draw phase angle differences in cartesian space
subplot(323)
filterplotDiffH1 = plot(real_data(1,:)-real_data(2,:),'b');
set(gca,'xlim',[timemin timemax],'ylim',[-10 10])
xlabel('Time (ms)')
ylabel('Voltage (\muV)')
title([ 'Filtered signal at ' num2str(center_freq) ' Hz' ])

% draw the phase angle time series
subplot(324)
phaseanglesDiffH1 = plot(phase_data(1,:)-phase_data(2,:),'b');
set(gca,'xlim',[timemin timemax],'ylim',[-pi pi]*2.2,'ytick',-2*pi:pi/2:pi*2)
xlabel('Time (ms)')
ylabel('Phase angle (radian)')
title('Phase angle time series')

% draw phase angles in polar space
subplot(325)
polar2chanH1 = polar([phase_data(1,:) phase_data(1,:)],'b');
hold on
polar2chanH2 = polar([phase_data(2,:) phase_data(2,:)],'m');
polar2chanH2 = polar([phase_data(1,:) phase_data(2,)]',repmat([0 1],1,1)','m');
title('Phase angles from two channels')
 
% draw phase angle differences in polar space
subplot(326)
polarAngleDiffH = polar([zeros(1,EEG_pnts) phase_data(2,:)-phase_data(1,:)]',repmat([0 1],1,1)','k');
polar(phase_data(2,:)-phase_data(1,:))
title('Phase angle differences from two channels')


%%% Figs 26.3
%%% Over trials

freqs2use  = logspace(log10(4),log10(30),15); % 4-30 Hz in 15 steps
times2save = -400:20:800;
timewindow = linspace(1.5,3,length(freqs2use)); % number of cycles on either end of the center point (1.5 means a total of 3 cycles))
baselinetm = [-400 -200];
% wavelet and FFT parameters
time          = -1:1/EEG_srate:1;
half_wavelet  = (length(time)-1)/2;
num_cycles    = logspace(log10(4),log10(8),length(freqs2use));
n_wavelet     = length(time);
n_data        = EEG_pnts*EEG_trials;
n_convolution = n_wavelet+n_data-1;

% time in indices
times2saveidx = dsearchn(EEG_times',times2save');
baselineidx   = dsearchn(times2save',baselinetm');

% initialize
ispc = zeros(length(freqs2use),length(times2save));
ps   = zeros(length(freqs2use),length(times2save));


% YS get the data of trials as concatenated
concattrial1 = EEG.covert_cond1.STN{1}(1,:)
concattrial2 = EEG.covert_cond1.STN{1}(2,:)
for trial=2:EEG_trials
	concattrial1= cat(2,concattrial1,EEG.covert_cond1.STN{trial}(1,:));
	concattrial2= cat(2,concattrial2,EEG.covert_cond1.STN{trial}(2,:));
end


% data FFTs
data_fft1 = fft(concattrial1,n_convolution);
data_fft2 = fft(concattrial2,n_convolution); %fft(reshape(EEG.data(chanidx(2),:,:),1,n_data),n_convolution);


for fi=1:length(freqs2use)
    
    % create wavelet and take FFT
    s = num_cycles(fi)/(2*pi*freqs2use(fi));
    wavelet_fft = fft( exp(2*1i*pi*freqs2use(fi).*time) .* exp(-time.^2./(2*(s^2))) ,n_convolution);
    
    % phase angles from channel 1 via convolution
    convolution_result_fft = ifft(wavelet_fft.*data_fft1,n_convolution);
    convolution_result_fft = convolution_result_fft(half_wavelet+1:end-half_wavelet);
    phase_sig1 = angle(reshape(convolution_result_fft,EEG_pnts,EEG_trials));
    
    % phase angles from channel 2 via convolution
    convolution_result_fft = ifft(wavelet_fft.*data_fft2,n_convolution);
    convolution_result_fft = convolution_result_fft(half_wavelet+1:end-half_wavelet);
    phase_sig2 = angle(reshape(convolution_result_fft,EEG_pnts,EEG_trials));
    
    % phase angle differences
    phase_diffs = phase_sig1-phase_sig2;
    
    % compute ICPS over trials
    ps(fi,:) = abs(mean(exp(1i*phase_diffs(times2saveidx,:)),2));
    
    % compute time window in indices for this frequency
    time_window_idx = round((1000/freqs2use(fi))*timewindow(fi)/(1000/EEG_srate));
%     time_window_idx = round(300/(1000/EEG_srate)); % set 300 to 100 for figure 3c/d
    
    for ti=1:length(times2save)
        
        % compute phase synchronization
        phasesynch = abs(mean(exp(1i*phase_diffs(times2saveidx(ti)-time_window_idx:times2saveidx(ti)+time_window_idx,:)),1));
        
        % average over trials
        ispc(fi,ti) = mean(phasesynch);
    end
end % end frequency loop

contourf(times2save,freqs2use,ispc-repmat(mean(ispc(:,baselineidx(1):baselineidx(2)),2),1,size(ispc,2)),20,'linecolor','none')
set(gca,'clim',[-.08 .08],'yscale','log','ytick',round(logspace(log10(freqs2use(1)),log10(freqs2use(end)),8)))
xlabel('Time (ms)'), ylabel('Frequency (Hz)')

figure
plot(freqs2use,(1000./freqs2use).*timewindow*2,'o-','markerface','k')
hold on
plot(freqs2use,(1000./freqs2use).*timewindow(1)*2,'ro-','markerface','m')
ylabel('Window width (ms)'), xlabel('Frequency (Hz)')
legend({'variable windows';'fixed 3*f window'})

