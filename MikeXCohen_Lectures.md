Mike X Cohen

# Lecture 1 
https://www.youtube.com/watch?v=ukjuFUghieg&list=PLn0OLiymPak0t1moK3sn4Sl1seXlEOPHT
Most if not all of neuroscience, is a source separation problem.
There is a source that is projecting to the sensors. For example, we cant measure attention directly but we can measure the activity of the neurons, in the real world we dont have direct access to the sources or latent constructs, but we can meaure the observables, neurons, time points etc.
The problem is that we can easily identify which sources feed the sensor, but we can fidn ways, mathematically to reconstruct the source.
Note that just by looking at the manifest variables themselves is insufficient to isolate the latent constructs.
How can sources be separated?
1. Anatomical source separation: based on anatomy i just look at a brain structure
2. Cognitive source separation: Clever experiment ot isolate one cog. process eg short term memory clearly, if the experiment is good, separatwd from other processes like long term memory, attention etc.
3. Temporal ss 
4. Spatial ss (see 5)

5. Statistical ss: use statistical characteristics of data, using descriptive statistics (covariane, spectral features etc.) no inferentila stats like p values.
Thus, statistical ss leads to time and spatial ss.
Sourec separation via Temporal or spatial filtering.

Temporal filtering is pretty much the same as spectral separation.
Amount of energy at a range of frequencies.
When 2 signals are mixed in the spectrum (signals have power in same range of frequencies) is impossible to separate them using purely spectral methods. But there are othe methods different from spectral separation (overlapping energy of signal and noise).

# Lecture 2
https://www.youtube.com/watch?v=Bmt89hHyxuM&list=PLn0OLiymPak0t1moK3sn4Sl1seXlEOPHT&index=2

# Lecture 3
https://www.youtube.com/watch?v=JMB9nZNGVyk&list=PLn0OLiymPak0t1moK3sn4Sl1seXlEOPHT&index=3
Preprocessing data:
Steps: 
* Import into MatL
* High pass filter (eg 0.5 Hz) avoud slow drifts artifacts
* import channel locations
* rereference eog, ekg, emg
--till here data is 2D : time x channels
* epoch data around stimuli
--now is 3D time x channels x epochs or trials
* substract pre stimuli baseline: avg voltage potential from each electrode fro eg -100ms to 0, and subbstract that average to the rest of the time course (similar to high pass filter) now all the trials (trial == a stimulus repetition) are in the same scale
* manual/visual trial rejection (or with algorithm)
* mark bad electrodes
* average reference eeg channels
* run ica to clean data (see components that want to reject)

# Lecture 4
Artifact visual removal/interpretation
https://www.youtube.com/watch?v=VDqwfP0mlfU&list=PLn0OLiymPak0t1moK3sn4Sl1seXlEOPHT&index=4

# Lecture 5
Topographical maps
https://www.youtube.com/watch?v=2edwDBSPLFs&list=PLn0OLiymPak0t1moK3sn4Sl1seXlEOPHT&index=5

# Lecture 6
Time Domain analysis (ERP) methods
https://www.youtube.com/watch?v=iFWrVzLYop0&list=PLn0OLiymPak0t1moK3sn4Sl1seXlEOPHT&index=6

ERP is the time domain averaging method.
Time-Frequency analysis negates the necessity of doing ERP analysis, because TF analysis encompasses ERP.
T-F analysis is a general framework that allows you to extract more information from the EEG, like spectral and connectivity information.

ERP are time domain analysis: 
Plot all the trials (n) for one channel.
the ERP from all trials  is the avergae of all voltage levels from each trial at each time point.
The ERP is about one order of magnitude smaller than the single trial variability. We loose meaningful information. 
ERP is so computationallly simple to compute that is why is being used since the beginning, but in now days we can compute ERPs with a mobile phone(people were computing ERPs in the 50s, 60s, 70s before digital computers)

In simple time domain ERP information is lost: because of non phase lock dynamics. If there is energy in a narrow frequency band that is phase lock (non randomized) you lose that info.
Averaging on the time domian we loose info, we need to average on the T-F domain


# Lecture 7
https://www.youtube.com/watch?v=3hk4z3yrMzk&list=PLn0OLiymPak0t1moK3sn4Sl1seXlEOPHT&index=7
If trials with phase-locked then trial average ERP is OK. But when non-phased locked ERP losses info.
T-F is crucial non only because phase-locked but because the brain works at different frequencies with specific rythmic activity

# Lecture 8
https://www.youtube.com/watch?v=s2MfmIx8wv4&list=PLn0OLiymPak0t1moK3sn4Sl1seXlEOPHT&index=8
 

ANTS book

PIV Ch 23. PCA
PV Ch 25. Intro to Connectivity analysis

Ch 26. Phase based connectivity


*PV Ch 25. Intro to Connectivity analysis*

Connectivity refers to 2 or more signals considered at a time. This most commonly refers to signals from two electrondes. Connectivity can be power and phase based, linear and non linear.
Multivariate or all to all or one to all are in reality mass bivariate, because we only take two at a time.

Why are most connectivity measures? Easier to interpret conceptually.
Trully or almost multivariate is graphs and ivariate correlations can inflate or misrepresent estimates of relationships if the network structure is actually multivariate. This is particularly relevant for brain connectivity because the brain is a highly multivariate system. For task-related connectivity this potential inflation is mitigated somewhat by condition comparisons because inflated connectivity estimates should affect all conditions.

Phase and power likely reflect different neurophysiological dynamocs, with phase likely reflecting the timing of activity within a neural population and power likely reflecting the number of neurons or spatial extent of the neural. 
Many connectivity results can be confounded by volume conduction.  

PHASE based connectivity
Phase-based connectivity analyses (described in greater detail in chapter 26) rely on the dis- tribution of phase angle differences between two electrodes, with the idea that when neural populations are functionally coupled, the timing of their oscillatory processes, as measured through phase, becomes synchronized.
 Phase-based measures rely on precise temporal relationships, usually in the identical frequency band, and are therefore susceptible to temporal jitter or uncertainty in the precise timing of experiment events. These temporal uncertainties can have more significant effects at higher frequencies. Phase-based measures do not provide compelling evidence for directionality

POWER based connectivity (27)
Correlating time-frequency power between two electrodes across time or over trials. These correlations can be computed between activity in the same or different frequencies and at the same or different time points. Power-based connectivity measures are arguably the most similar to connectivity measures often used in fMRI. the correlated fluctuations in activity are relatively slower, compared to phase-based connec- tivity measures. Power-based connectivity measures are also relatively insensitive to temporal jitter and uncertainty

Granger (28)
MI (29)
Cross frequency connectivity (30)
There is a potentially huge search space (frequencies × frequencies × electrodes × electrodes × conditions × time), 
Graph theory (31)

Confound of Volume Conduction
Volume conduction is a potential confound that can lead to spurious connectivity results. sources in the brain generate large electromagnetic fields that are measured by more than one EEG electrode or MEG sensor, thus introducing


27. Power based connectivity
Functional connectivity analysis based on fluctuations in power.

Phase-based connectivity analyses assume that the connectivity is instantaneous and at the same frequency. 
Power-based connectivity analyses do not have this constraint, which make them more flexible for hypothesis-driven as well as exploratory analyses. (However power based assumes zero phase-lag)

Within-subject power correlations should be done with the Spearman and not the Pearson correlation: power data are nonnormally distributed, and they contain outliers.






###############################################################################################
Who Is Fourier?

Chapter 1 shows brilliantly the intuition behind Fourier idea that adding up many simple waves you can get a wave as complex as you like. 
Next step, Chapter 2, is how to break a complex wave into simple waves.

The Fourier series allows us to treat any wave as complicated as we loke as the sum of siple waves, as long as the complicated wave is periodic.



Chapter 10.
Euler's formal: a way to describe sin and cos usong i and e numbers.



#######################
Cadiz paper notes. POWER BASED CONNECTIVITY
------------------
Phase based connectivity assumes that the connectivity is instantaneous and at the same frequency, power based connectivity does not have this constraint.

Spearman correlation (rank correlation) eliminates the influence of outliers.
Power data are nonnormally distributed.
Real EEG data, time-frequency power data are nonnormally distributed and contain outliers. A logarithm-base-ten transform helps the distribution look more normal, although in this case, the distribution is still technically nonnormal, according to the Kolmogorov-Smirnov test.

Note that trial-averaged baseline-corrected power (decibel or percentage change) is often normally distributed and therefore can be tested using Pearson correlations.

Typically, power data are averaged across all trials, and then a baseline normalization is applied. 

 The baseline is a period of time, typically a few hundred mil- liseconds before the start of the trial, when little or no task-related processing is expected.
2396851886

power_corr_trialsov_occ-cov_occ.png
27.6 A
Power correlations over trials. The correlation is calculated using predetermined time-frequency windows (see x- and y-axis labels; each dot corresponds to one trial). In this case there is no relationship. 


TF Power and Baseline normalization 217
The frequency spectrum of data tends to show decreasing power at increasing frequencies. This is not charactwristic of EEG but also characterizes the relationship between power and frequency of many signals (radio, radiation from the big bang). This decrease in power with frequecny follows a 1/f shape. the more general form is c/f^x which is a power law, one variable is a power function of another variable (frequency), and this is why it is difficult to visualize power from a large range of frequency bands simultaneously, typically, theta will have much more power than beta .

Limitations:
1. difficult to visualize power across a large range of fre- quency bands
2. raw power values change in scale as a function of frequency, lower frequencies will usually show a seemingly larger effect than higher frequencies in terms of the overall magnitude. 
3. aggregating effects across subjects can be difficult with raw power value, because raw power is affected by individual anatomical properties (skull thickness, sulcal anatomy, recording environment, electrodes impedance...)
4. task- related changes in power can be difficult to disentangle from background activity, this is especially true for frequencies with high power in baseline periods eg alpha over occipital electrondes.
5. raw power values are not normally distributed because they cannot be negative and they are strongly positively skewed

To address this limitations: 1/f shape can be attenuated by taking the logarithm of the power values. But this alone doesnt fix the 5 limitations. Fortunately, there is one way to address them all: use a baseline normalization.
There are at least 3 posssible baseline normalization methods:

i. Decibel Conversion
The decibel (dB) is a ratio between the strength of one signal (frequency-band-specific power) and the strength of another signal (a baseline level of power in that same frequency band).
The decibel scale is convenient because it is robust to many limitations of raw power, such as 1/f power scaling and subject-specific and electrode-specific idiosyncratic characteristics.



##### 
Prior to study connectivity, do TF analysis.
baselinetime = np.array([ -800, -500 ])
Condition 1.
Fro [200,800] and [1300,1700] O > C in alpha (8,10)
C1.1 [200,800]
C1.2 [1300,1700]



## PHASE BASED CONNECTIVITY
ISPC average of phase angle differences between electrodes over time.
ISPC_f = 1/n \sum e^{i\phi_{1t}-\phi_{2t}}

n is the number of points, phi1 and phi2 are the two channels, and ISPC is always computed at a frequency f. We are substracting phase angle between to electrodes

i) 26.1 26.2 ISPC from one trial
ii) 26.3 ISPC for trials
----

ISPC is symmetric(non directional)
#Fig 26.1 (A) filtered signal (bandpass) and phase angles of the 2 channels
          (B) filtered signal and phase angles differences betweeen the 2 channels
          (C) individual phase of each electrode, see how they are distributed in polar space. (D) substraction of the difference, see if they cluster(nonuniformilly distrb)

# Fig 26.2 to prove that ISPC does not depend on the phase angle themselves, only the clusterng of phase angle differences


--

