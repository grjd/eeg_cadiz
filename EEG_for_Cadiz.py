#EEG_for_Cadiz.py
#https://mark-kramer.github.io/Case-Studies-Python/02.html#case-study-data
#https://github.com/lyndond/Analyzing_Neural_Time_Series
from scipy.io import loadmat 
from numpy.fft import fft, ifft
from numpy.polynomial.polynomial import polyfit
from scipy.stats import spearmanr, pearsonr, kstest, rankdata
import matplotlib.pyplot as plt
import os
import pandas as pd
import numpy as np
import pdb
from matplotlib import cm

figures_dir = '/Users/borri/github/EEG_Cadiz_2021/figures/Frontal-STN/thr1std' #normaliz_db normaliz_div, normaliz_pc, normaliz_z
figures_dir = '/Users/borri/github/EEG_Cadiz_2021/figures/18-11-2021'
# fig_name = os.path.join(figures_dir, 'occ_ov_avg.png')
# plt.savefig(fig_name)

# fig, ax = plt.subplots(2, figsize=(12,5))
# fig.suptitle('EEG signal Occipital cov vs ov')
# ax[0].plot(t_ms, occ_cov_avg)
# ax[0].plot([0, 0], [-40,40], 'k--', lw=1)
# ax[0].set_ylabel('Voltage ')
# ax[1].set_ylabel('Voltage ')
# ax[1].set_xlabel('time ')
# ax[1].plot(t_ms, occ_ov_avg)
# ax[1].plot([0, 0], [-40,40], 'k--', lw=1)
# fig_name = os.path.join(figures_dir, 'occ_ov_12_avg.png')
# plt.savefig(fig_name)

# Visualize one trial at a time
# get dataframe trials x time (for each channel)
#plt.imshow()

# Chapter 27. Power Based connectivity

def power_partialcorr_all_freqs(data_seed,data_target,data_control, seed_chan,target_chan, control_chan,t_ms):
	"""27.7 with control channel
	"""

	EEGtrials = 115
	s_rate = 1000
	min_freq = 2
	max_freq = 40
	num_frex = 15
	wavelet_cycles = 4.
	frequencies = np.logspace(np.log10(min_freq), np.log10(max_freq),num_frex)

	times2save = t_ms
	times2save_idx = [np.argmin(np.abs(t_ms - t)) for t in times2save]
	tf_corr_data = np.zeros((len(frequencies), len(times2save), 2))

	timesize = min(data_seed.shape[0],data_target.shape[0],data_control.shape[0])
	trialsize = min(data_seed.shape[1],data_target.shape[1],data_control.shape[1])
	data_seed= data_seed.iloc[:,0:trialsize]
	data_target= data_target.iloc[:,0:trialsize]
	data_control= data_control.iloc[:,0:trialsize]
	data_seed=data_seed.to_numpy().flatten()
	data_target=data_target.to_numpy().flatten()
	data_control=data_control.to_numpy().flatten()
	#force same size all signals

	def bivar_corr(x, y):
		return 1-6*np.sum((x-y)**2)/(EEGtrials*(EEGtrials**2-1))

	for fi in range(len(frequencies)):
		conv_result_seed = rankdata(filt(data_seed, frequencies[fi], wavelet_cycles, s_rate)[times2save_idx, :], axis=1)
		conv_result_target = rankdata(filt(data_target, frequencies[fi], wavelet_cycles, s_rate)[times2save_idx,:], axis=1)
		conv_result_control = rankdata(filt(data_control, frequencies[fi], wavelet_cycles,s_rate)[times2save_idx, :], axis=1)

		for ti in range(len(times2save)):
			r_st = bivar_corr(conv_result_seed[ti], conv_result_target[ti])
			r_sc = bivar_corr(conv_result_seed[ti], conv_result_control[ti])
			r_tc = bivar_corr(conv_result_control[ti], conv_result_target[ti])
			tf_corr_data[fi, ti, 0 ] = r_st
			# partial correlation
			tf_corr_data[fi, ti,1] = (r_st-r_sc*r_tc) / (np.sqrt(1-r_sc**2)*np.sqrt(1-r_tc**2))

	msg = 'Partial Rank Corr_' + seed_chan + '-' + target_chan+'|' + control_chan
	fig, ax = plt.subplots(1, 2)
	vmin = 1.1*tf_corr_data.min()
	vmax = 1.1*tf_corr_data.max()

	ax[0].contourf(times2save, frequencies, tf_corr_data[..., 1], vmin=vmin, vmax=vmax)
	#ax[0].set(yscale='log')
	ax[0].set_title(msg)
	ax[0].set_xlabel(r'Time(ms)')
	ax[0].set_ylabel("Frequency(Hz)")
	#ax[1].contourf(times2save, frequencies, tf_corr_data[..., 1], vmin=0, vmax=.6)
	#ax[1].set(yscale='log')
	msg = '(log scale)'
	cs = ax[1].contourf(times2save, frequencies, tf_corr_data[..., 1], vmin=vmin, vmax=vmax)
	ax[1].set(yscale='log')
	ax[1].set_title(msg)
	ax[1].set_xlabel(r'Time(ms)')
	ax[1].set_ylabel("log Frequency(Hz)")
	fig.colorbar(cs, ax=ax[1], shrink=0.9)

	fig_name = os.path.join(figures_dir, 'Partial Rank Corr_' + seed_chan + '-' + target_chan +'|' + control_chan +'.png')
	fig.tight_layout()
	plt.savefig(fig_name,bbox_inches='tight',dpi=400)
	# 
	pdb.set_trace()

def power_corr_all_freqs(data_seed,data_target,seed_chan,target_chan,t_ms):
	"""27.7 no partial corr
	"""

	EEGtrials = 115
	s_rate = 1000
	EEGsrate = 1000
	min_freq = 2
	max_freq = 40
	num_frex = 15
	wavelet_cycles = 4.
	frequencies = np.logspace(np.log10(min_freq), np.log10(max_freq),num_frex)

	times2save = t_ms
	times2save_idx = [np.argmin(np.abs(t_ms - t)) for t in times2save]
	tf_corr_data = np.zeros((len(frequencies), len(times2save), 2))

	data_seed=data_seed.to_numpy().flatten()
	data_target=data_target.to_numpy().flatten()
	

	def bivar_corr(x, y):
		return 1-6*np.sum((x-y)**2)/(EEGtrials*(EEGtrials**2-1))

	for fi in range(len(frequencies)):
		conv_result_seed = rankdata(filt(data_seed, frequencies[fi], wavelet_cycles, s_rate)[times2save_idx, :], axis=1)
		conv_result_target = rankdata(filt(data_target, frequencies[fi], wavelet_cycles, s_rate)[times2save_idx,:], axis=1)
		#conv_result_control = rankdata(filt(data_control, frequencies[fi], wavelet_cycles)[t_ms, :], axis=1)

		for ti in range(len(times2save)):
			r_st = bivar_corr(conv_result_seed[ti], conv_result_target[ti])
			#r_sc = bivar_corr(conv_result_seed[ti], conv_result_control[ti])
			#r_tc = bivar_corr(conv_result_control[ti], conv_result_target[ti])

			tf_corr_data[fi, ti, 0 ] = r_st

			# partial correlation
			#tf_corr_data[fi, ti,1] = (r_st-r_sc*r_tc) / (np.sqrt(1-r_sc**2)*np.sqrt(1-r_tc**2)
	msg = 'Rank Corr_' + seed_chan + '-' + target_chan	
	fig, ax = plt.subplots(1, 2)
	vmin = 1.1*tf_corr_data.min()
	vmax = 1.1*tf_corr_data.max()

	ax[0].contourf(times2save, frequencies, tf_corr_data[..., 0], vmin=vmin, vmax=vmax)
	#ax[0].set(yscale='log')
	ax[0].set_title(msg)
	ax[0].set_xlabel(r'Time(ms)')
	ax[0].set_ylabel("Frequency(Hz)")
	#ax[1].contourf(times2save, frequencies, tf_corr_data[..., 1], vmin=0, vmax=.6)
	#ax[1].set(yscale='log')
	msg = 'Rank Corr (log scale)'
	cs = ax[1].contourf(times2save, frequencies, tf_corr_data[..., 0], vmin=vmin, vmax=vmax)
	ax[1].set(yscale='log')
	ax[1].set_title(msg)
	ax[1].set_xlabel(r'Time(ms)')
	ax[1].set_ylabel("log Frequency(Hz)")
	fig.colorbar(cs, ax=ax[1], shrink=0.9)

	fig_name = os.path.join(figures_dir, 'vmM' + msg +'.png')
	fig.tight_layout()
	plt.savefig(fig_name,bbox_inches='tight',dpi=400)


def power_correlation_trials_aux(data1,data2,channel1,channel2,time_win1,time_win2,center_freq1,center_freq2,EEGtrials,EEGtimes):
	""" 27.6  time-frequency-electrode windows to use a priori for each electrode
	"""
	timepoints = data1.shape[1]
	wavelet_cycles = 4.5
	EEGsrate = 1000 #==s_rate = 1000
	time_idx1 = np.zeros_like(time_win1)
	time_idx2 = np.zeros_like(time_win2)
	#channel1= data1.name; channel2= data2.name
	msg=channel1 + '' + str(center_freq1) + 'Hz' + '-' + channel2  + '' + str(center_freq2) + 'Hz'
	for i in range(2):
		time_idx1[i] = np.argmin(np.abs(EEGtimes-time_win1[i]))
		time_idx2[i] = np.argmin(np.abs(EEGtimes-time_win2[i]))
	time = np.arange(-1, 1+1/EEGsrate, 1/EEGsrate)
	half_of_wavelet_size = len(time)//2
	wavelet_cycles = 4.5
	# data times x trials
	#data1_rs = np.reshape(data1,values, (data1.shape[0]*data1.shape[1]), order="f")

	data1 = data1.to_numpy().flatten()
	data2 = data2.to_numpy().flatten()
	
	analytic_signal1 = filt(data1, center_freq1, wavelet_cycles, EEGsrate, EEGtrials, timepoints)
	analytic_signal2 = filt(data2, center_freq2, wavelet_cycles, EEGsrate, EEGtrials, timepoints)
	tf_window_data1 = analytic_signal1[time_idx1[0]:time_idx1[1], :].mean(0)
	tf_window_data2 = analytic_signal2[time_idx2[0]:time_idx2[1], :].mean(0)

	#Plot Power TF window correlation (27.6 A)
	fig, ax = plt.subplots(1,1)
	#ax[0].plot(tf_window_data1, tf_window_data2, '.')
	#ax[0].set(title=f'tf corr',xlabel=f'{channel1}: {time_win1[0]}-{time_win1[1]}; {center_freq1}Hz',ylabel=f'{channel2}: {time_win2[0]}-{time_win2[1]};{center_freq2}Hz')
	ax.plot(rankdata(tf_window_data1), rankdata(tf_window_data2), '.')
	ax.set(title='rank corr',xlabel=f'{channel1}: {time_win1[0]}-{time_win1[1]}; {center_freq1}Hz',ylabel=f'{channel2}: {time_win2[0]}-{time_win2[1]};{center_freq2}Hz')
	fig_name = os.path.join(figures_dir, 'Rank_power_corr_trials' + msg +'.png')
	plt.savefig(fig_name,bbox_inches='tight',dpi=400)
	
	# 27.6 B power correlation at each time point over trials
	corr_ts = np.zeros_like(EEGtimes)
	for i in range(timepoints):
		corr_ts[i]  = spearmanr(analytic_signal1[i,:], analytic_signal2[i,:])[0]

	fig, ax = plt.subplots(1,1)
	ax.plot(corr_ts)
	ax.axvline(x=2000, color='k', linestyle='--')
	ax.set(title='Power corr at each time across trials: '+ msg,xlabel='Time(ms)',ylabel=r'Spearman $\rho$')
	fig_name = os.path.join(figures_dir, 'Rank_ppower_corr_ts' + msg +'.png')
	plt.savefig(fig_name,bbox_inches='tight',dpi=400)
	
	# 27.6 C
	# Spectrum over trials and freq
	times2save = np.arange(-200, 2000, 20)
	frex = np.logspace(np.log10(2), np.log10(40), 20)
	times2save_idx = [np.argmin(np.abs(EEGtimes - t)) for t in times2save]
	seeddata_rank = rankdata(tf_window_data2)
	expl_corrs = np.zeros((len(frex), len(times2save)))
	for fi in range(len(frex)):
		analytic_signal1 = filt(data1, frex[fi], wavelet_cycles, EEGsrate, EEGtrials, timepoints)
		for ti in range(len(times2save)):
			y = analytic_signal1[times2save_idx[ti],:]
			expl_corrs[fi, ti] = 1 - 6 * sum((seeddata_rank -rankdata(y))**2) / (EEGtrials*(EEGtrials**2-1))

	fig, ax = plt.subplots(1,1)
	cs = ax.contourf(times2save, frex, expl_corrs,  origin='lower')
	fig.colorbar(cs, ax=ax, shrink=0.9)
	ax.set(title=f'target {channel1} corr over trials from seed {channel2}, {center_freq2}Hz'f'{time_win2[0]}-{time_win2[1]}ms', xlabel='time (ms)', ylabel='frequency (Hz)')
	fig_name = os.path.join(figures_dir, 'Spectrum_overtrials' + msg +'.png')
	plt.savefig(fig_name,bbox_inches='tight',dpi=400)
	# contourf only above threshold

	fig, ax = plt.subplots(1,1)
	threshold = 0.3
	threshold = expl_corrs.mean() + 1.5*expl_corrs.std()

	#mask_th = np.ma.masked_outside(expl_corrs, -threshold, threshold,0)
	expl_corrs_abs = np.abs(expl_corrs)
	#mask_th = np.where(np.logical_and(expl_corrs<= -.2, expl_corrs>= .2) ,expl_corrs,0)
	mask_th = np.where(expl_corrs>= threshold ,expl_corrs,0)
	cs2 = ax.contourf(times2save, frex, mask_th,  origin='lower',cmap=plt.cm.Greys_r)
	fig.colorbar(cs2, ax=ax, shrink=0.9)
	ax.set(title=f'Target {channel1} corr > {round(threshold,2)} From Seed {channel2}, {center_freq2}Hz'f'{time_win2[0]}-{time_win2[1]}ms', xlabel='time (ms)', ylabel='frequency (Hz)')
	fig_name = os.path.join(figures_dir, 'Spectrum_overtrials '+ 'Mask_' + str(round(threshold,2)) + msg +'.png')
	plt.savefig(fig_name,bbox_inches='tight',dpi=400)

	return tf_window_data1, tf_window_data2, expl_corrs


def power_corr_trials(data1,data2,sensor1,sensor2,time_win1,time_win2,center_freq1,center_freq2,EEGtimes):
	""" 27.6  time-frequency-electrode windows to use a priori for each electrode
	"""
	EEGtrials = 115
	EEGsrate = 1000
	time_idx1 = np.zeros_like(time_win1)
	time_idx2 = np.zeros_like(time_win2)
	msg=sensor1 + '-' + sensor2
	for i in range(2):
		time_idx1[i] = np.argmin(np.abs(EEGtimes-time_win1[i]))
		time_idx2[i] = np.argmin(np.abs(EEGtimes-time_win2[i]))
	time = np.arange(-1, 1+1/EEGsrate, 1/EEGsrate)
	half_of_wavelet_size = len(time)//2
	wavelet_cycles = 4.5
	# data times x trials
	#data1_rs = np.reshape(data1,values, (data1.shape[0]*data1.shape[1]), order="f")
	print('analytic_signal1')

	#data1 = np.reshape(data1, (data1.shape[0]*data1.shape[1]), order="f")
	data1=data1.to_numpy().flatten()
	data2=data2.to_numpy().flatten()
	#data2 = np.reshape(EEGdata[EEGchanlocslabels==sensor2], (EEGpnts*EEGtrials), order="f"
	analytic_signal1 = filt(data1, center_freq1, wavelet_cycles,EEGsrate)
	analytic_signal2 = filt(data2, center_freq1, wavelet_cycles,EEGsrate)
	tf_window_data1 = analytic_signal1[time_idx1[0]:time_idx1[1], :].mean(0)
	tf_window_data2 = analytic_signal2[time_idx2[0]:time_idx2[1], :].mean(0)

	fig, ax = plt.subplots(1,2)
	ax[0].plot(tf_window_data1, tf_window_data2, '.')
	ax[0].set(title=f'tf corr',xlabel=f'{sensor1}: {time_win1[0]}-{time_win1[1]}; {center_freq1}Hz',ylabel=f'{sensor2}: {time_win2[0]}-{time_win2[1]};{center_freq2}Hz')

	ax[1].plot(rankdata(tf_window_data1), rankdata(tf_window_data2), '.')
	ax[1].set(title='rank corr')
	fig_name = os.path.join(figures_dir, 'power_corr_trials' + msg +'.png')
	plt.savefig(fig_name)
	
	# power correlation at each tiemn point over trials
	corr_ts = np.zeros_like(EEGtimes)
	for i in range(4001):
		corr_ts[i]  = spearmanr(analytic_signal1[i,:], analytic_signal2[i,:])[0]
	
	fig, ax = plt.subplots(1,1)
	ax.plot(corr_ts)
	fig_name = os.path.join(figures_dir, 'power_corr_overtrials' + msg +'.png')
	plt.savefig(fig_name)

	# Spectrum over trials and freq
	times2save = np.arange(-200, 1225, 25)
	frex = np.logspace(np.log10(2), np.log10(40), 20)
	times2save_idx = [np.argmin(np.abs(EEGtimes - t)) for t in times2save]
	seeddata_rank = rankdata(tf_window_data2)
	expl_corrs = np.zeros((len(frex), len(times2save)))
	for fi in range(len(frex)):
		analytic_signal1 = filt(data1, frex[fi],wavelet_cycles, EEGsrate)
		for ti in range(len(times2save)):
			y = analytic_signal1[times2save_idx[ti],:]
			expl_corrs[fi, ti] = 1 - 6 * sum((seeddata_rank -rankdata(y))**2) / (EEGtrials*(EEGtrials**2-1))

	fig, ax = plt.subplots(1,1)

	ax.contourf(times2save, frex, expl_corrs)
	ax.set(title=f'corr over trials from seed {sensor2}, {center_freq2}Hz and'f'{time_win2[0]}-{time_win2[1]}ms', xlabel='time (ms)', ylabel='frequency (Hz)')
	fig_name = os.path.join(figures_dir, 'Spectrum_overtrials' + msg +'.png')
	plt.savefig(fig_name)


def connectivity_2channels(data1, data2, trials1, trials2, timepoints, center_freq, wavelet_cycles, msg=None):
	"""
	"""
	s_rate = 1000
	#center_freq = 10
	#trial2plot = 9
	wavelet_cycles = 4.5
	EEGtimes = np.linspace(-2000,2000, num=timepoints)

	channel1= msg.split('-')[0]
	channel2= msg.split('-')[1].split('_freq')[0]#msg.split('-')[1][1]
	timei, timee = 1000, 1800
	times2plot = np.argwhere(np.logical_and(EEGtimes> timei, EEGtimes<timee)).squeeze()

	filtered1 = filt(data1, center_freq[0], wavelet_cycles, s_rate, trials1, timepoints)
	filtered2 = filt(data2, center_freq[1], wavelet_cycles, s_rate, trials2, timepoints)
	min_tt = min(filtered1.shape[1], filtered2.shape[1])
	# force having same number of trials
	filtered1  = filtered1[:,0:min_tt]
	filtered2 = filtered2[:,0:min_tt]	
	if filtered1[:].shape[1] > 0:
		filtered1 = filtered1[:].mean(axis=1)
		filtered2 = filtered2[:].mean(axis=1)

	#  Plot series and correlation coefficients
	fig, ax = plt.subplots(1,1)

	#ax.plot(EEGtimes[times2plot], filtered1[times2plot])
	#ax.plot(EEGtimes[times2plot], filtered2[times2plot])	
	#ax.plot(filtered1[:]) 
	#ax.plot(filtered2[:])	
	ax.plot(filtered1[times2plot])  
	ax.plot(filtered2[times2plot])
	ax.legend()
	fig_name = os.path.join(figures_dir, 'Signal:' + msg +'.png')
	plt.savefig(fig_name)

	fig, ax = plt.subplots(1,1)
	ax.plot(filtered1[times2plot], filtered2[times2plot], '.')
	ax.set(xlabel=channel1, ylabel=channel2,title=f"Pearson r{pearsonr(filtered1[times2plot], filtered2[times2plot])[0]:.2f}")
	ax.legend()
	fig_name = os.path.join(figures_dir, 'PwPw_corr: ' + msg +'.png')
	plt.savefig(fig_name)

	fig, ax = plt.subplots(1,1)
	ax.plot(rankdata(filtered1[times2plot]), rankdata(filtered2[times2plot]), '.')
	ax.set(xlabel=channel1, ylabel=channel2,title=f"spearmanr r{spearmanr(filtered1[times2plot], filtered2[times2plot])[0]:.2f}")
	ax.legend()
	fig_name = os.path.join(figures_dir, 'PwSp_corr: ' + msg +'.png')
	plt.savefig(fig_name)

	return filtered1, filtered2


def normality_test(ts):
	"""
	"""
	print(kstest(ts.flatten(), 'norm'))
	print(kstest(np.log10(ts.flatten()), 'norm'))
	print(kstest(np.random.randn(10000), 'norm'))

def filt(data, center_freq, wavelet_cycles,s_rate, trials=None, timepoints=None):
	"""
	"""

	time = np.arange(-1, 1 + 1/s_rate, 1/s_rate)  # time for wavelet from -1 to 1 in secs
	sg = wavelet_cycles/(2*np.pi*center_freq)  # stdev of Gaussian
	wavelet = np.exp(2*np.pi*1j*center_freq*time) * np.exp(-time**2 / (2*sg**2))

	# fft params
	n_wavelet = len(wavelet)
	n_data = len(data)
	#n_data = data.shape[1]
	n_convolution = n_data+n_wavelet-1
	half_of_wavelet_size = len(wavelet)//2

	# fft of wavelet and eeg data
	convolution_result_fft = ifft(fft(wavelet, n_convolution) * fft(data, n_convolution))*np.sqrt(sg)
	filtered = convolution_result_fft[half_of_wavelet_size:-half_of_wavelet_size]
	filtered_dim = np.abs(np.reshape(filtered, (timepoints,trials), order="f"))**2
	#filtered = np.abs(filtered)**2
	#filtered = np.abs(np.reshape(filtered, (EEGpnts, EEGtrials), order="f"))**2
	return filtered_dim


def plot_power_dist(data, timepoints=None,trials=None, msg=None):
	"""
	"""


	s_rate = 1000
	center_freq = 10
	wavelet_cycles = 4.5

	convolution_result_fft = filt(data, center_freq, wavelet_cycles,s_rate, trials, timepoints)
	convolution_result_fft =  convolution_result_fft[100:-100]
	# plot distr power data
	fig, ax = plt.subplots(1,2)
	ax[0].hist(convolution_result_fft.flatten())
	ax[0].set_title("Dist power")
	ax[0].set_xlabel(r'$muV^2$ Power')
	ax[0].set_ylabel("")
	ax[1].hist(np.log10(convolution_result_fft.flatten()))
	ax[1].set_title(r"Dist $log_{10}$ power")
	#ax[1].set_xlabel(r'$log_10$ ppwer')

	fig_name = os.path.join(figures_dir, msg + '.png')
	plt.savefig(fig_name)
	return convolution_result_fft



def load_subj5_data(sub5_mat_fpath):
	"""
	"""
	data = loadmat('subj5.mat')
	data.keys()
	#The keys that start and end with two underscores ( __ ) are private and contain information 
	#about the MATLAB file; we will not need those keys here. 
	#dict_keys(['__header__', '__version__', '__globals__', 'subj5'])
	# Extract the two conditions: covert, overt
	data = data['subj5']
	covert = data[0]['covert_cond1']
	overt = data[0]['overt_cond2']

	# Electrodes
	occ_ov = overt[0]['OCCIPITAL']
	fro_ov = overt[0]['FRONTAL']
	stn_ov = overt[0]['STN']
	occ_cov = covert[0]['OCCIPITAL']
	fro_cov = covert[0]['FRONTAL']
	stn_cov = covert[0]['STN']
	# Print shape of signal
	time_pts = occ_cov[0][0][0][0].shape[1]
	t_ms = np.linspace(-2000,2000, num=time_pts)
	cha_occ = occ_cov[0][0][0][0].shape[0]
	cha_fro = fro_cov[0][0][0][0].shape[0]
	cha_stn = stn_cov[0][0][0][0].shape[0]
	print('OCC covert signal shape=' + str(occ_cov[0][0][0][0].shape) + ' trials=' + str(occ_cov[0][0][0].shape[0]))
	print('OCC overt signal shape=' + str(occ_ov[0][0][0][0].shape) + ' trials=' + str(occ_ov[0][0][0].shape[0]))
	print('FTL covert signal shape=' + str(fro_cov[0][0][0][0].shape) + ' trials=' + str(fro_cov[0][0][0].shape[0]))
	print('FTL overt signal shape=' + str(fro_ov[0][0][0][0].shape) + ' trials=' + str(fro_ov[0][0][0].shape[0]))
	print('STN covert signal shape=' + str(stn_cov[0][0][0][0].shape) + ' trials=' + str(stn_cov[0][0][0].shape[0]))
	print('STN overt signal shape=' + str(stn_ov[0][0][0][0].shape) + ' trials=' + str(stn_ov[0][0][0].shape[0]))

	occ_cov_arr = occ_cov[0][0][0]
	fro_cov_arr = fro_cov[0][0][0]
	stn_cov_arr = stn_cov[0][0][0]

	occ_ov_arr = occ_ov[0][0][0]
	fro_ov_arr = fro_ov[0][0][0]
	stn_ov_arr = stn_ov[0][0][0]
	
	# Get dataframe trials X (channels x time)
	trial = occ_cov_arr[5]

	occ_trials = occ_ov_arr.shape[0]
	fro_trials = fro_ov_arr.shape[0]
	stn_trials = stn_ov_arr.shape[0]
	
	# mean of all channels for the trial
	channel_mu = trial.mean(0)

	# time x trial
	df_occ_ov = pd.DataFrame([]);df_occ_cov = pd.DataFrame([])
	df_fro_ov = pd.DataFrame([]);df_fro_cov = pd.DataFrame([])
	df_stn_ov = pd.DataFrame([]);df_stn_cov = pd.DataFrame([])
	for i in np.arange(0, occ_trials):
		df_occ_ov['T'+str(i)] = occ_ov_arr[i].mean(axis=0)
		df_occ_cov['T'+str(i)] = occ_cov_arr[i].mean(axis=0)
	for i in np.arange(0, fro_trials):
		df_fro_ov['T'+str(i)] = fro_ov_arr[i].mean(axis=0)
		df_fro_cov['T'+str(i)] = fro_cov_arr[i].mean(axis=0)
	for i in np.arange(0, stn_trials):
		df_stn_ov['T'+str(i)] = stn_ov_arr[i].mean(axis=0)
		df_stn_cov['T'+str(i)] = stn_cov_arr[i].mean(axis=0)

	# Transpose (trials x time)
	df_occ_ov = df_occ_ov.T
	# Plot df
	ntrials = occ_cov_arr.shape[0]
	fig = plt.figure()
	plt.imshow(df_occ_ov,cmap='BuPu',extent=[t_ms[0], t_ms[-1], 0, ntrials],aspect='auto',origin='lower' )
	plt.xlabel('Time[s]') # Label the axes
	plt.ylabel('OCC OVERT Trial #')
	plt.colorbar() # Show voltage to color mapping
	plt.vlines(0.25, 1, 1000, 'k', lw=2)

	fig_name = os.path.join(figures_dir, 'IMocc_ov_avg.png')
	plt.savefig(fig_name)

	# Mean accross trials
	ntrials = occ_cov[0][0][0][:].shape[0]
	occ_cov[0][0][0][:].mean(axis=0)

	# Do the average of all channels
	# occ_cov[0][0][0][110][8] trial channel
	occ_cov_avg = occ_cov[0][0][0][0][:].mean(axis=0)
	# specific trial
	#occ_cov_avg = occ_cov[0][0][0][1][:].mean(axis=0)
	#occ_cov_avg = occ_cov[0][0][0][110][:].mean(axis=0)
	###
	occ_ov_avg = occ_ov[0][0][0][0][:].mean(axis=0)
	fro_cov_avg = fro_cov[0][0][0][0][:].mean(axis=0)
	fro_ov_avg = fro_ov[0][0][0][0][:].mean(axis=0)
	stn_cov_avg = stn_cov[0][0][0][0][:].mean(axis=0)
	stn_ov_avg = stn_ov[0][0][0][0][:].mean(axis=0)
	averages = [occ_cov_avg, occ_ov_avg, fro_cov_avg, fro_ov_avg, stn_cov_avg, stn_ov_avg]
	timextrials = [df_occ_cov, df_occ_ov, df_fro_cov,df_fro_ov,df_stn_cov,df_stn_ov]
	return averages, timextrials

def plot_trials_time(df, t_ms, msg=None):
	""" df must be in Trials x time
	"""
	# Plot df
	ntrials = df.shape[0]
	fig = plt.figure()
	plt.imshow(df,cmap='BuPu',extent=[t_ms[0], t_ms[-1], 0, ntrials],aspect='auto',origin='lower' )
	plt.xlabel('Time[s]') # Label the axes
	plt.ylabel(msg + 'Trial #')
	plt.colorbar() # Show voltage to color mapping
	#plt.vlines(0.25, 1, 1000, 'k', lw=2)
	fig_name = os.path.join(figures_dir, msg + '.png')
	plt.savefig(fig_name)

	plt.savefig(fig_name)

def aggregate_subjects(subjects):
	"""IN list f dataframe one for each subject. Out the mean across subjects
	"""
	#pdb.set_trace()
	df_cond_roi = pd.concat(subjects, axis=0)
	cond_roi_mean = df_cond_roi.mean(axis=0)
	return cond_roi_mean

def aggregate_trials_subjects(subjects):
	"""
	"""
	
	print('Aggreagte subjects (trials x time')
	n_subjects = len(subjects)
	min_t = 1000
	for i in range(0, n_subjects):
		trials = subjects[0].shape[0]
		if trials < min_t:min_t = trials
	df_cond_roi = pd.concat(subjects, axis=0)
	#pdb.set_trace()
	df_tt = pd.DataFrame([]);
	for tt in range(0, min_t):
		df_tt['T'+str(tt)] = df_cond_roi.loc['T'+str(tt)].mean(axis=0)
	
	df_tt= df_tt.T
	return df_tt


def group_analysis():
	"""
	"""
	fpath = '/Users/borri/github/EEG_Cadiz_2021/LFP_EEG.mat'
	fdata = loadmat(fpath)
	print(fdata['__header__'])
	eeg = fdata['LFP_EEG']
	eeg_labels = eeg['labels']
	eeg_data = eeg['data'][0][0]
	covert = eeg_data['covert_cond1']
	overt = eeg_data['overt_cond2']
	##
	covert_stn = covert[0][0]['STN']
	covert_frontal = covert[0][0]['FRONTAL']
	covert_occipital = covert[0][0]['OCCIPITAL']

	
	
	overt_stn = overt[0][0]['STN']
	overt_frontal = overt[0][0]['FRONTAL']
	overt_occipital = overt[0][0]['OCCIPITAL']

	s1_ov_occ = overt_occipital[0][0]['subj1'][0][0];s1_cov_occ = covert_occipital[0][0]['subj1'][0][0]
	s2_ov_occ = overt_occipital[0][0]['subj2'][0][0];s2_cov_occ = covert_occipital[0][0]['subj2'][0][0]
	s3_ov_occ = overt_occipital[0][0]['subj3'][0][0];s3_cov_occ = covert_occipital[0][0]['subj3'][0][0]
	s4_ov_occ = overt_occipital[0][0]['subj4'][0][0];s4_cov_occ = covert_occipital[0][0]['subj4'][0][0]
	s5_ov_occ = overt_occipital[0][0]['subj5'][0][0];s5_cov_occ = covert_occipital[0][0]['subj5'][0][0]
	s6_ov_occ = overt_occipital[0][0]['subj6'][0][0];s6_cov_occ = covert_occipital[0][0]['subj6'][0][0]
	s7_ov_occ = overt_occipital[0][0]['subj7'][0][0];s7_cov_occ = covert_occipital[0][0]['subj7'][0][0]

	s1_ov_fro = overt_frontal[0][0]['subj1'][0][0];s1_cov_fro = covert_frontal[0][0]['subj1'][0][0]
	s2_ov_fro = overt_frontal[0][0]['subj2'][0][0];s2_cov_fro = covert_frontal[0][0]['subj2'][0][0]
	s3_ov_fro = overt_frontal[0][0]['subj3'][0][0];s3_cov_fro = covert_frontal[0][0]['subj3'][0][0]
	s4_ov_fro = overt_frontal[0][0]['subj4'][0][0];s4_cov_fro = covert_frontal[0][0]['subj4'][0][0]
	s5_ov_fro = overt_frontal[0][0]['subj5'][0][0];s5_cov_fro = covert_frontal[0][0]['subj5'][0][0]
	s6_ov_fro = overt_frontal[0][0]['subj6'][0][0];s6_cov_fro = covert_frontal[0][0]['subj6'][0][0]
	s7_ov_fro = overt_frontal[0][0]['subj7'][0][0];s7_cov_fro = covert_frontal[0][0]['subj7'][0][0]

	s1_ov_stn = overt_stn[0][0]['subj1'][0][0];s1_cov_stn = covert_stn[0][0]['subj1'][0][0]
	s2_ov_stn = overt_stn[0][0]['subj2'][0][0];s2_cov_stn = covert_stn[0][0]['subj2'][0][0]
	s3_ov_stn = overt_stn[0][0]['subj3'][0][0];s3_cov_stn = covert_stn[0][0]['subj3'][0][0]
	s4_ov_stn = overt_stn[0][0]['subj4'][0][0];s4_cov_stn = covert_stn[0][0]['subj4'][0][0]
	s5_ov_stn = overt_stn[0][0]['subj5'][0][0];s5_cov_stn = covert_stn[0][0]['subj5'][0][0]
	s6_ov_stn = overt_stn[0][0]['subj6'][0][0];s6_cov_stn = covert_stn[0][0]['subj6'][0][0]
	s7_ov_stn = overt_stn[0][0]['subj7'][0][0];s7_cov_stn = covert_stn[0][0]['subj7'][0][0]
	
	# mean of all channels for the trial
	num_trials = s1_ov_occ.shape[1]
	num_channels = s1_ov_occ[0][0].shape[0]

	# time x trial
	df_s1_cov_occ = pd.DataFrame([]);df_s2_cov_occ = pd.DataFrame([])
	df_s3_cov_occ = pd.DataFrame([]);df_s4_cov_occ = pd.DataFrame([])
	df_s5_cov_occ = pd.DataFrame([]);df_s6_cov_occ = pd.DataFrame([])
	df_s7_cov_occ = pd.DataFrame([]);
	df_s1_ov_occ = pd.DataFrame([]);df_s2_ov_occ = pd.DataFrame([])
	df_s3_ov_occ = pd.DataFrame([]);df_s4_ov_occ = pd.DataFrame([])
	df_s5_ov_occ = pd.DataFrame([]);df_s6_ov_occ = pd.DataFrame([])
	df_s7_ov_occ = pd.DataFrame([]);
	df_s1_cov_fro = pd.DataFrame([]);df_s2_cov_fro = pd.DataFrame([])
	df_s3_cov_fro = pd.DataFrame([]);df_s4_cov_fro = pd.DataFrame([])
	df_s5_cov_fro = pd.DataFrame([]);df_s6_cov_fro = pd.DataFrame([])
	df_s7_cov_fro = pd.DataFrame([]);
	df_s1_ov_fro = pd.DataFrame([]);df_s2_ov_fro = pd.DataFrame([])
	df_s3_ov_fro = pd.DataFrame([]);df_s4_ov_fro = pd.DataFrame([])
	df_s5_ov_fro = pd.DataFrame([]);df_s6_ov_fro = pd.DataFrame([])
	df_s7_ov_fro = pd.DataFrame([]);
	df_s1_cov_stn = pd.DataFrame([]);df_s2_cov_stn = pd.DataFrame([])
	df_s3_cov_stn = pd.DataFrame([]);df_s4_cov_stn = pd.DataFrame([])
	df_s5_cov_stn = pd.DataFrame([]);df_s6_cov_stn = pd.DataFrame([])
	df_s7_cov_stn = pd.DataFrame([]);
	df_s1_ov_stn = pd.DataFrame([]);df_s2_ov_stn = pd.DataFrame([])
	df_s3_ov_stn = pd.DataFrame([]);df_s4_ov_stn = pd.DataFrame([])
	df_s5_ov_stn = pd.DataFrame([]);df_s6_ov_stn = pd.DataFrame([])
	df_s7_ov_stn = pd.DataFrame([]);

	# Overt occipital  subject x condition x area
	for i in np.arange(0, s1_ov_occ.shape[1]):
		df_s1_ov_occ['T'+str(i)] = s1_ov_occ[0,i].mean(axis=0)
	for i in np.arange(0, s2_ov_occ.shape[1]):
		df_s2_ov_occ['T'+str(i)] = s2_ov_occ[0,i].mean(axis=0)
	for i in np.arange(0, s3_ov_occ.shape[1]):	
		df_s3_ov_occ['T'+str(i)] = s3_ov_occ[0,i].mean(axis=0)
	for i in np.arange(0, s4_ov_occ.shape[1]):	
		df_s4_ov_occ['T'+str(i)] = s4_ov_occ[0,i].mean(axis=0)
	for i in np.arange(0, s5_ov_occ.shape[1]):
		df_s5_ov_occ['T'+str(i)] = s5_ov_occ[0,i].mean(axis=0)
	for i in np.arange(0, s6_ov_occ.shape[1]):
		df_s6_ov_occ['T'+str(i)] = s6_ov_occ[0,i].mean(axis=0)
	for i in np.arange(0, s7_ov_occ.shape[1]):
		df_s7_ov_occ['T'+str(i)] = s7_ov_occ[0,i].mean(axis=0)
	#Covert occipital
	for i in np.arange(0, s1_cov_occ.shape[1]):
		df_s1_cov_occ['T'+str(i)] = s1_cov_occ[0,i].mean(axis=0)
	for i in np.arange(0, s2_cov_occ.shape[1]):
		df_s2_cov_occ['T'+str(i)] = s2_cov_occ[0,i].mean(axis=0)
	for i in np.arange(0, s3_cov_occ.shape[1]):	
		df_s3_cov_occ['T'+str(i)] = s3_cov_occ[0,i].mean(axis=0)
	for i in np.arange(0, s4_cov_occ.shape[1]):	
		df_s4_cov_occ['T'+str(i)] = s4_cov_occ[0,i].mean(axis=0)
	for i in np.arange(0, s5_cov_occ.shape[1]):
		df_s5_cov_occ['T'+str(i)] = s5_cov_occ[0,i].mean(axis=0)
	for i in np.arange(0, s6_cov_occ.shape[1]):
		df_s6_cov_occ['T'+str(i)] = s6_cov_occ[0,i].mean(axis=0)
	for i in np.arange(0, s7_cov_occ.shape[1]):
		df_s7_cov_occ['T'+str(i)] = s7_cov_occ[0,i].mean(axis=0)	

	# Frontal
	# Overt Frontal  subject x condition x area
	for i in np.arange(0, s1_ov_fro.shape[1]):
		df_s1_ov_fro['T'+str(i)] = s1_ov_fro[0,i].mean(axis=0)
	for i in np.arange(0, s2_ov_fro.shape[1]):
		df_s2_ov_fro['T'+str(i)] = s2_ov_fro[0,i].mean(axis=0)
	for i in np.arange(0, s3_ov_fro.shape[1]):	
		df_s3_ov_fro['T'+str(i)] = s3_ov_fro[0,i].mean(axis=0)
	for i in np.arange(0, s4_ov_fro.shape[1]):	
		df_s4_ov_fro['T'+str(i)] = s4_ov_fro[0,i].mean(axis=0)
	for i in np.arange(0, s5_ov_fro.shape[1]):
		df_s5_ov_fro['T'+str(i)] = s5_ov_fro[0,i].mean(axis=0)
	for i in np.arange(0, s6_ov_fro.shape[1]):
		df_s6_ov_fro['T'+str(i)] = s6_ov_fro[0,i].mean(axis=0)
	for i in np.arange(0, s7_ov_fro.shape[1]):
		df_s7_ov_fro['T'+str(i)] = s7_ov_fro[0,i].mean(axis=0)
	#Covert occipital
	for i in np.arange(0, s1_cov_fro.shape[1]):
		df_s1_cov_fro['T'+str(i)] = s1_cov_fro[0,i].mean(axis=0)
	for i in np.arange(0, s2_cov_fro.shape[1]):
		df_s2_cov_fro['T'+str(i)] = s2_cov_fro[0,i].mean(axis=0)
	for i in np.arange(0, s3_cov_fro.shape[1]):	
		df_s3_cov_fro['T'+str(i)] = s3_cov_fro[0,i].mean(axis=0)
	for i in np.arange(0, s4_cov_fro.shape[1]):	
		df_s4_cov_fro['T'+str(i)] = s4_cov_fro[0,i].mean(axis=0)
	for i in np.arange(0, s5_cov_fro.shape[1]):
		df_s5_cov_fro['T'+str(i)] = s5_cov_fro[0,i].mean(axis=0)
	for i in np.arange(0, s6_cov_fro.shape[1]):
		df_s6_cov_fro['T'+str(i)] = s6_cov_fro[0,i].mean(axis=0)
	for i in np.arange(0, s7_cov_fro.shape[1]):
		df_s7_cov_fro['T'+str(i)] = s7_cov_fro[0,i].mean(axis=0)

	# STN
	# Overt STN  subject x condition x area
	for i in np.arange(0, s1_ov_stn.shape[1]):
		df_s1_ov_stn['T'+str(i)] = s1_ov_stn[0,i].mean(axis=0)
	for i in np.arange(0, s2_ov_stn.shape[1]):
		df_s2_ov_stn['T'+str(i)] = s2_ov_stn[0,i].mean(axis=0)
	for i in np.arange(0, s3_ov_stn.shape[1]):	
		df_s3_ov_stn['T'+str(i)] = s3_ov_stn[0,i].mean(axis=0)
	for i in np.arange(0, s4_ov_stn.shape[1]):	
		df_s4_ov_stn['T'+str(i)] = s4_ov_stn[0,i].mean(axis=0)
	for i in np.arange(0, s5_ov_stn.shape[1]):
		df_s5_ov_stn['T'+str(i)] = s5_ov_stn[0,i].mean(axis=0)
	for i in np.arange(0, s6_ov_stn.shape[1]):
		df_s6_ov_stn['T'+str(i)] = s6_ov_stn[0,i].mean(axis=0)
	for i in np.arange(0, s7_ov_stn.shape[1]):
		df_s7_ov_stn['T'+str(i)] = s7_ov_stn[0,i].mean(axis=0)
	#Covert occipital
	for i in np.arange(0, s1_cov_stn.shape[1]):
		df_s1_cov_stn['T'+str(i)] = s1_cov_stn[0,i].mean(axis=0)
	for i in np.arange(0, s2_cov_stn.shape[1]):
		df_s2_cov_stn['T'+str(i)] = s2_cov_stn[0,i].mean(axis=0)
	for i in np.arange(0, s3_cov_stn.shape[1]):	
		df_s3_cov_stn['T'+str(i)] = s3_cov_stn[0,i].mean(axis=0)
	for i in np.arange(0, s4_cov_stn.shape[1]):	
		df_s4_cov_stn['T'+str(i)] = s4_cov_stn[0,i].mean(axis=0)
	for i in np.arange(0, s5_cov_stn.shape[1]):
		df_s5_cov_stn['T'+str(i)] = s5_cov_stn[0,i].mean(axis=0)
	for i in np.arange(0, s6_cov_stn.shape[1]):
		df_s6_cov_stn['T'+str(i)] = s6_cov_stn[0,i].mean(axis=0)
	for i in np.arange(0, s7_cov_stn.shape[1]):
		df_s7_cov_stn['T'+str(i)] = s7_cov_stn[0,i].mean(axis=0)

	
	#Transpose df_occ_ov = df_occ_ov.T 
	df_s1_cov_stn = df_s1_cov_stn.T; df_s2_cov_stn = df_s2_cov_stn.T; df_s3_cov_stn = df_s3_cov_stn.T; 
	df_s4_cov_stn = df_s4_cov_stn.T; df_s5_cov_stn = df_s5_cov_stn.T; df_s6_cov_stn = df_s6_cov_stn.T; 
	df_s7_cov_stn = df_s7_cov_stn.T;

	df_s1_ov_stn = df_s1_ov_stn.T; df_s2_ov_stn = df_s2_ov_stn.T; df_s3_ov_stn = df_s3_ov_stn.T;
	df_s4_ov_stn = df_s4_ov_stn.T; df_s5_ov_stn = df_s5_ov_stn.T;df_s6_ov_stn = df_s6_ov_stn.T;
	df_s7_ov_stn = df_s7_ov_stn.T; 

	df_s1_cov_fro = df_s1_cov_fro.T; df_s2_cov_fro = df_s2_cov_fro.T; df_s3_cov_fro = df_s3_cov_fro.T; 
	df_s4_cov_fro = df_s4_cov_fro.T; df_s5_cov_fro = df_s5_cov_fro.T; df_s6_cov_fro = df_s6_cov_fro.T; 
	df_s7_cov_fro = df_s7_cov_fro.T;

	df_s1_ov_fro = df_s1_ov_fro.T; df_s2_ov_fro = df_s2_ov_fro.T; df_s3_ov_fro = df_s3_ov_fro.T;
	df_s4_ov_fro = df_s4_ov_fro.T; df_s5_ov_fro = df_s5_ov_fro.T;df_s6_ov_fro = df_s6_ov_fro.T;
	df_s7_ov_fro = df_s7_ov_fro.T; 

	df_s1_cov_occ = df_s1_cov_occ.T; df_s2_cov_occ = df_s2_cov_occ.T; df_s3_cov_occ = df_s3_cov_occ.T; 
	df_s4_cov_occ = df_s4_cov_occ.T; df_s5_cov_occ = df_s5_cov_occ.T; df_s6_cov_occ = df_s6_cov_occ.T; 
	df_s7_cov_occ = df_s7_cov_occ.T;

	df_s1_ov_occ = df_s1_ov_occ.T; df_s2_ov_occ = df_s2_ov_occ.T; df_s3_ov_occ = df_s3_ov_occ.T;
	df_s4_ov_occ = df_s4_ov_occ.T; df_s5_ov_occ = df_s5_ov_occ.T; df_s6_ov_occ = df_s6_ov_occ.T;
	df_s7_ov_occ = df_s7_ov_occ.T; 

	# Get aggregate of all subjects: condition_area
	ov_stn = [df_s1_ov_stn, df_s2_ov_stn, df_s3_ov_stn, df_s4_ov_stn, df_s5_ov_stn, df_s6_ov_stn, df_s7_ov_stn]
	cov_stn = [df_s1_cov_stn, df_s2_cov_stn, df_s3_cov_stn, df_s4_cov_stn, df_s5_cov_stn, df_s6_cov_stn, df_s7_cov_stn]
	ov_fro = [df_s1_ov_fro, df_s2_ov_fro, df_s3_ov_fro, df_s4_ov_fro, df_s5_ov_fro, df_s6_ov_fro, df_s7_ov_fro]
	cov_fro = [df_s1_cov_fro, df_s2_cov_fro, df_s3_cov_fro, df_s4_cov_fro, df_s5_cov_fro, df_s6_cov_fro, df_s7_cov_fro]
	ov_occ = [df_s1_ov_occ, df_s2_ov_occ, df_s3_ov_occ, df_s4_ov_occ, df_s5_ov_occ, df_s6_ov_occ, df_s7_ov_occ]
	cov_occ = [df_s1_cov_occ, df_s2_cov_occ, df_s3_cov_occ, df_s4_cov_occ, df_s5_cov_occ, df_s6_cov_occ, df_s7_cov_occ]
	# (1 x time)
	ov_stn_arr = aggregate_subjects(ov_stn)
	cov_stn_arr = aggregate_subjects(cov_stn)
	# only 3 trials remove
	del ov_fro[3]
	ov_fro_arr = aggregate_subjects(ov_fro)
	cov_fro_arr = aggregate_subjects(cov_fro)
	# only 3 trials remove
	del ov_occ[3]
	ov_occ_arr = aggregate_subjects(ov_occ)
	cov_occ_arr = aggregate_subjects(cov_occ)
	# (trials x time)

	ov_stn_df = aggregate_trials_subjects(ov_stn)
	cov_stn_df = aggregate_trials_subjects(cov_stn)
	ov_fro_df = aggregate_trials_subjects(ov_fro)
	cov_fro_df = aggregate_trials_subjects(cov_fro)
	ov_occ_df = aggregate_trials_subjects(ov_occ)
	cov_occ_df = aggregate_trials_subjects(cov_occ)

	trials_mean = [ov_stn_arr,cov_stn_arr,ov_fro_arr,cov_fro_arr,ov_occ_arr,cov_occ_arr]
	trials_per_time = [ov_stn_df,cov_stn_df,ov_fro_df,cov_fro_df,ov_occ_df,cov_occ_df]

	return trials_mean, trials_per_time
	pdb.set_trace()

	time_pts = df_s7_cov_stn.shape[1];t_ms = np.linspace(-2000,2000, num=time_pts)
	msg = 'S7_covert_STN'
	plot_trials_time(df_s7_cov_stn, t_ms, msg)
	df_s7_ov_stn = df_s7_ov_stn.T;time_pts = df_s7_cov_stn.shape[1];t_ms = np.linspace(-2000,2000, num=time_pts)
	msg = 'S7_overt_STN'
	plot_trials_time(df_s7_ov_stn, t_ms, msg)
	
	time_pts = df_s6_cov_stn.shape[1];t_ms = np.linspace(-2000,2000, num=time_pts)
	msg = 'S6_covert_STN'
	plot_trials_time(df_s7_cov_stn, t_ms, msg)
	df_s6_ov_stn = df_s6_ov_stn.T;time_pts = df_s6_cov_stn.shape[1];t_ms = np.linspace(-2000,2000, num=time_pts)
	msg = 'S6_overt_STN'
	plot_trials_time(df_s6_ov_stn, t_ms, msg)

	time_pts = df_s5_cov_stn.shape[1];t_ms = np.linspace(-2000,2000, num=time_pts)
	msg = 'S5_covert_STN'
	plot_trials_time(df_s5_cov_stn, t_ms, msg)
	df_s5_ov_stn = df_s5_ov_stn.T;time_pts = df_s5_cov_stn.shape[1];t_ms = np.linspace(-2000,2000, num=time_pts)
	msg = 'S5_overt_STN'
	plot_trials_time(df_s5_ov_stn, t_ms, msg)
	
	time_pts = df_s4_cov_stn.shape[1];t_ms = np.linspace(-2000,2000, num=time_pts)
	msg = 'S4_covert_STN'
	plot_trials_time(df_s4_cov_stn, t_ms, msg)
	df_s4_ov_stn = df_s4_ov_stn.T;time_pts = df_s4_cov_stn.shape[1];t_ms = np.linspace(-2000,2000, num=time_pts)
	msg = 'S4_overt_STN'
	plot_trials_time(df_s4_ov_stn, t_ms, msg)

	time_pts = df_s3_cov_stn.shape[1];t_ms = np.linspace(-2000,2000, num=time_pts)
	msg = 'S3_covert_STN'
	plot_trials_time(df_s3_cov_stn, t_ms, msg)
	df_s3_ov_stn = df_s3_ov_stn.T;time_pts = df_s3_cov_stn.shape[1];t_ms = np.linspace(-2000,2000, num=time_pts)
	msg = 'S3_overt_STN'
	plot_trials_time(df_s3_ov_stn, t_ms, msg)

	time_pts = df_s2_cov_stn.shape[1];t_ms = np.linspace(-2000,2000, num=time_pts)
	msg = 'S2_covert_STN'
	plot_trials_time(df_s2_cov_stn, t_ms, msg)
	df_s2_ov_stn = df_s2_ov_stn.T;time_pts = df_s2_cov_stn.shape[1];t_ms = np.linspace(-2000,2000, num=time_pts)
	msg = 'S2_overt_STN'
	plot_trials_time(df_s2_ov_stn, t_ms, msg)

	time_pts = df_s1_cov_stn.shape[1];t_ms = np.linspace(-2000,2000, num=time_pts)
	msg = 'S1_covert_STN'
	plot_trials_time(df_s1_cov_stn, t_ms, msg)
	df_s1_ov_stn = df_s1_ov_stn.T;time_pts = df_s1_cov_stn.shape[1];t_ms = np.linspace(-2000,2000, num=time_pts)
	msg = 'S1_overt_STN'
	plot_trials_time(df_s1_ov_stn, t_ms, msg)

 
def time_frequency_data(channel, msg=None):
	"""18.2 no color scaling works optimally for all frequency bands."""

	EEGsrate = 1000
	EEGpnts = channel.shape[1]

	# average per trial
	EEGdata = channel.mean(axis=0)
	# wavelet parameters
	min_freq = 2; max_freq = 64; num_frex = 30 #ax_freq = 128
	frequencies = np.logspace(np.log10(min_freq),np.log10(max_freq),num_frex)
	time = np.arange(-1,1+1/EEGsrate,1/EEGsrate)
	half_of_wavelet_size = len(time)//2
	pdb.set_trace()
	def nextpow2(n):
		m_f = np.log2(n)
		m_i = np.ceil(m_f)
		return 2**m_i

	#FFT parameters (use next-power-of-2)
	n_wavelet     = len(time)
	n_data        = EEGpnts;	
	n_convolution = n_wavelet+n_data-1;
	n_conv_pow2   = int(nextpow2(n_convolution))
	wavelet_cycles= 4
	fft_data = fft(np.squeeze(EEGdata),n_conv_pow2)
	#initialize output time-frequency data
	tf_data = np.zeros([len(frequencies),EEGpnts])
	for fi in range(len(frequencies)):
		#create wavelet and get its FFT
		wavelet = (np.pi*frequencies[fi]*np.sqrt(np.pi))**(-.5) * np.exp(2*1j*np.pi*frequencies[fi]*time)* np.exp(-time**2/(2*( wavelet_cycles /(2*np.pi*frequencies[fi]))**2))/frequencies[fi]
		fft_wavelet = fft(wavelet,n_conv_pow2)
		#run convolution
		convolution_result_fft = ifft(fft_wavelet*fft_data,n_conv_pow2)
		convolution_result_fft = convolution_result_fft[:n_convolution] # note: here we remove the extra points from the power-of-2 FFT
		convolution_result_fft = convolution_result_fft[half_of_wavelet_size:-half_of_wavelet_size]
		#put power data into time-frequency matrix
		tf_data[fi,:] = np.abs(convolution_result_fft)**2

	ytickskip = np.arange(2,num_frex+4,4) 
	
	fig, ax = plt.subplots(1, 1)
	ax.imshow(tf_data,origin="lower",aspect="auto",cmap=cm.jet, vmin=0, vmax=25,extent=[-500,1500,0,30])
	plt.setp(plt.gca(),'yticks',ytickskip,'yticklabels',np.round(frequencies[ytickskip-1]),'xlim',[-500, 1500])
	ax.set_title('Color limit of 0 to 25 ' + msg) 
	ax.set_xlabel('Time (ms)')
	ax.set_ylabel('Frequency (Hz)')
	fig_name = os.path.join(figures_dir, 'tf_25' + msg + '.png')
	plt.savefig(fig_name, bbox_inches='tight',dpi=400)

	fig, ax = plt.subplots(1, 1)
	ax.imshow(np.log10(tf_data),origin="lower",aspect="auto",cmap=cm.jet, vmin=-4, vmax=4,extent=[-500,1500,0,30])
	plt.setp(plt.gca(),'yticks',ytickskip,'yticklabels',np.round(frequencies[ytickskip-1]),'xlim',[-500, 1500])
	ax.set_title('Color limit of log10 -4 to 4 ' + msg)
	fig_name = os.path.join(figures_dir, 'tf_log' + msg + '.png')
	plt.savefig(fig_name, bbox_inches='tight',dpi=400)
	return tf_data


def plot_normalized_power(signalconverted, msg=None):
	"""signalconverted is db percent or z
	"""
	fig, ax = plt.subplots(1, 2)

def plot_2signals(tf_data,dbconverted,msg=None):
	"""
	"""
	EEGpnts = tf_data.shape[1]
	EEGtimes = np.linspace(-2000,2000, num=EEGpnts) 
	# wavelet parameters
	min_freq = 2; max_freq = 128; num_frex = 30
	frequencies = np.logspace(np.log10(min_freq),np.log10(max_freq),num_frex)
	time2plot = 100 # in ms
	timeidx = np.argmin(np.abs(EEGtimes-time2plot))
	fig, ax = plt.subplots(1, 2)
	ax[0].plot(frequencies,tf_data[:,timeidx])
	ax[0].set_title( 'Power spectrum at ' + str(EEGtimes[timeidx]) + ' ms' )
	ax[0].set_ylabel(r'Raw power (\muV^2)')
	ax[0].set_xlabel('Frequency (Hz)')

	ax[1].plot(frequencies,dbconverted[:,timeidx])
	ax[1].set_ylabel('Baseline-normalized power (dB)')
	ax[1].set_xlabel('Frequency (Hz)')

	fig_name = os.path.join(figures_dir, 'power_dec_at-'+ str(time2plot) + ' ms ' + msg + '.png')
	plt.savefig(fig_name, bbox_inches='tight',dpi=400)



def baseline_normalization_percent(tf_data, msg=None):
	"""
	"""
	EEGpnts = tf_data.shape[1]
	# wavelet parameters
	min_freq = 2; max_freq = 128; num_frex = 30
	frequencies = np.logspace(np.log10(min_freq),np.log10(max_freq),num_frex)

	EEGtimes = np.linspace(-2000,2000, num=EEGpnts) 
	baselinetime = np.array([ -800, -500 ])
	baselineidx =np.zeros(2,dtype=int)
	# convert baseline window time to indices
	baselineidx[0]=np.argmin(np.abs(EEGtimes-baselinetime[0]))
	baselineidx[1]=np.argmin(np.abs(EEGtimes-baselinetime[1]))
	baseline_power = np.mean(tf_data[:,baselineidx[0]:baselineidx[1]],axis=1)
	pctchange = 100 * (tf_data - np.reshape(baseline_power,[len(baseline_power),1]))/np.reshape(baseline_power,[len(baseline_power),1])
	return pctchange

def baseline_normalization_baselinediv(tf_data, msg=None):
	"""
	"""
	EEGpnts = tf_data.shape[1]
	# wavelet parameters
	min_freq = 2; max_freq = 128; num_frex = 30
	frequencies = np.logspace(np.log10(min_freq),np.log10(max_freq),num_frex)

	EEGtimes = np.linspace(-2000,2000, num=EEGpnts) 
	baselinetime = np.array([ -800, -500 ])
	baselineidx =np.zeros(2,dtype=int)
	# convert baseline window time to indices
	baselineidx[0]=np.argmin(np.abs(EEGtimes-baselinetime[0]))
	baselineidx[1]=np.argmin(np.abs(EEGtimes-baselinetime[1]))
	baseline_power = np.mean(tf_data[:,baselineidx[0]:baselineidx[1]],axis=1)
	baselinediv = tf_data / np.reshape(baseline_power,[len(baseline_power),1])
	return baselinediv

def baseline_normalization_zeta(tf_data, msg=None):
	"""
	"""
	EEGpnts = tf_data.shape[1]
	# wavelet parameters
	min_freq = 2; max_freq = 128; num_frex = 30
	frequencies = np.logspace(np.log10(min_freq),np.log10(max_freq),num_frex)

	EEGtimes = np.linspace(-2000,2000, num=EEGpnts) 
	baselinetime = np.array([ -800, -500 ])
	baselineidx =np.zeros(2,dtype=int)
	# convert baseline window time to indices
	baselineidx[0]=np.argmin(np.abs(EEGtimes-baselinetime[0]))
	baselineidx[1]=np.argmin(np.abs(EEGtimes-baselinetime[1]))
	baseline_power = tf_data[:,baselineidx[0]:baselineidx[1]]
	baselineZ = (tf_data - np.reshape(np.mean(baseline_power,axis=1),[len(baseline_power),1]))/np.reshape(np.std(baseline_power,axis=1),[len(baseline_power),1])
	return baselineZ

def baseline_normalization_decibels(tf_data, msg=None):
	"""
	"""
	
	EEGpnts = tf_data.shape[1]
	# wavelet parameters
	min_freq = 2; max_freq = 128; num_frex = 30
	frequencies = np.logspace(np.log10(min_freq),np.log10(max_freq),num_frex)

	EEGtimes = np.linspace(-2000,2000, num=EEGpnts) 
	baselinetime = np.array([ -800, -500 ])
	baselineidx =np.zeros(2,dtype=int)
	# convert baseline window time to indices
	baselineidx[0]=np.argmin(np.abs(EEGtimes-baselinetime[0]))
	baselineidx[1]=np.argmin(np.abs(EEGtimes-baselinetime[1]))
	#dB-correct
	baseline_power = np.mean(tf_data[:,baselineidx[0]:baselineidx[1]],axis=1)
	dbconverted = 10*np.log10(tf_data/np.reshape(baseline_power,[len(baseline_power),1]))
	
	fig, ax = plt.subplots(1, 1)
	cs = ax.contourf(EEGtimes,frequencies,dbconverted,40,vmin=-12,vmax=12,cmap=cm.jet)
	plt.setp(plt.gca(),'yticks',np.round(np.logspace(np.log10(frequencies[0]),np.log10(frequencies[-1])*100)/100))
	ax.set_yscale('log')
	ax.set_title('Decibel conversion -12 +12 dB ' + msg)
	ax.set_xlabel('Time (ms)')
	ax.set_ylabel('Frequency (Hz)')

	fig.colorbar(cs, ax=ax, shrink=0.9)

	fig_name = os.path.join(figures_dir, 'baseline_tf' + msg + '.png')
	plt.savefig(fig_name, bbox_inches='tight',dpi=400)
	return dbconverted

def baseline_normalize_channel(channel, msg):
	""" channel, return 4 type of normalization
	"""
	# Convert data to TF data
	channel_tf = time_frequency_data(channel, msg)
	channel_db = baseline_normalization_decibels(channel_tf, msg)
	channel_z = baseline_normalization_zeta(channel_tf, msg)
	channel_pc = baseline_normalization_percent(channel_tf, msg)
	channel_div = baseline_normalization_baselinediv(channel_tf, msg)
	return channel_db, channel_z, channel_pc, channel_div, channel_tf

def compare_power_strectra(tf1,tf2, msg=None):
	"""
	"""
	print('Comparing ' + msg)
	
	tf_diff = tf1 - tf2

	min_freq = 2; max_freq = 64; num_frex = 30
	ytickskip = np.arange(2,num_frex+4,4) 
	frequencies = np.logspace(np.log10(min_freq),np.log10(max_freq),num_frex)
	#time = np.arange(-1,1+1/EEGsrate,1/EEGsrate)

	# fig = plt.figure()
	# ax1 = fig.add_subplot(221)
	# cm1 = plt.cm.Reds
	# cm1.set_under('white')
	# ax1.imshow(tf1,origin="lower",aspect="auto",vmin=0, vmax=25,extent=[-100,2000,0,30], cmap = cm1)
	# plt.setp(plt.gca(),'yticks',ytickskip,'yticklabels',np.round(frequencies[ytickskip-1]),'xlim',[-100, 2000])

	# ax2 = fig.add_subplot(223)
	# cm2 = plt.cm.Blues
	# cm2.set_under('white')
	# ax2.imshow(tf2, origin="lower",aspect="auto",vmin=0, vmax=25,extent=[-100,2000,0,30],cmap = cm1)
	# plt.setp(plt.gca(),'yticks',ytickskip,'yticklabels',np.round(frequencies[ytickskip-1]),'xlim',[-100, 2000])

	# ax3 = fig.add_subplot(222)
	# ax3.imshow(tf1, origin="lower",aspect="auto",vmin=0, vmax=25,extent=[-100,2000,0,30], alpha=0.5,cmap = cm1)
	# ax3.imshow(tf2, origin="lower",aspect="auto",vmin=0, vmax=25,extent=[-100,2000,0,30],alpha=0.5,cmap = cm2)
	# plt.setp(plt.gca(),'yticks',ytickskip,'yticklabels',np.round(frequencies[ytickskip-1]),'xlim',[-100, 2000])
	# #fig.colorbar(im3)
	# fig_name = os.path.join(figures_dir, 'all3C11_PowerSpectra_' + msg + '.png')
	# plt.savefig(fig_name, bbox_inches='tight',dpi=400)
	# pdb.set_trace()


	fig, ax = plt.subplots(1, 1)
	im = ax.imshow(tf1,origin="lower",aspect="auto",cmap=cm.Blues, vmin=0, vmax=25,extent=[-100,2000,0,30],alpha=.5)
	im2 = ax.imshow(tf2,origin="lower",aspect="auto",cmap=cm.Reds, vmin=0, vmax=25,extent=[-100,2000,0,30],alpha=.5)
	plt.axvline(x=500, color='grey', linestyle='--');plt.axvline(x=800, color='grey', linestyle='--')
	plt.axvline(x=1200, color='grey', linestyle='--');#plt.axvline(x=1600, color='grey', linestyle='--')
	plt.setp(plt.gca(),'yticks',ytickskip,'yticklabels',np.round(frequencies[ytickskip-1]),'xlim',[-100, 2000])
	ax.set_title('Power spectra Diff ' + msg.split()[0] + ' (B) -' + msg.split()[-1] + ' (R)' ) 
	ax.set_xlabel('Time (ms)')
	ax.set_ylabel('Frequency (Hz)')
	cbar1 = fig.colorbar(im, shrink=0.7); #cbar1.ax.set_xlabel(msg.split()[0], rotation=270)
	cbar2 = fig.colorbar(im2, shrink=0.7);#cbar2.ax.set_xlabel(msg.split()[1], rotation=270)
	fig_name = os.path.join(figures_dir,  'C1_PowerSpectra_Z_' + msg + '.png')
	plt.savefig(fig_name, bbox_inches='tight',dpi=400)

def connectivity_pair(channel1, channel2, center_freq):
	""" plot power correlation between 2 channels at a specific frequency
	""" 
	EEGsrate = 1000
	# center_freq = 4
	wavelet_cycles = 4.5
	msg = channel1.name + '-' +  channel2.name + '_freq_' + str(center_freq[0]) + '_' + str(center_freq[1])
	data1  = channel1.to_numpy();
	EEGpnts = data1.shape[1]
	EEGtrials1 = data1.shape[0]
	data1 = np.reshape(data1, (EEGpnts*EEGtrials1), order="f")
	data2  = channel2.to_numpy();
	EEGtrials2 = data2.shape[0]
	data2 = np.reshape(data2, (EEGpnts*EEGtrials2), order="f")
	center_freq = [10, 20] #alpha, beta
	fil1, fil2 = connectivity_2channels(data1, data2, EEGtrials1, EEGtrials2, EEGpnts, center_freq, wavelet_cycles,msg)
	return fil1, fil2

def power_correlation_trials(df1, df2, center_freq1, center_freq2):
	"""
	"""
	#df1 = ov_fro; df2 = ov_stn;
	#center_freq2=4;center_freq1=9
	sensor1=df1.name; sensor2=df2.name; #sensor3= 'stn_cov'
	# pass dataframes with same number of trials
	min_tt = min(df1.shape[0], df2.shape[0])
	df1 = df1[0:min_tt]; df2 = df2[0:min_tt]; 
	#time_win1=(-500,-100); #time_win2=(200,1000);
	#time_win1=(500,800); time_win2=(500,800);
	tf_window_data1, tf_window_data2 = power_corr_trials_new(df1, df2, sensor1, sensor2, time_win1, time_win2, center_freq1,center_freq2,min_tt,t_ms)
	return tf_window_data1, tf_window_data2 



def plot_difference_spectra(expl_corrs1, expl_corrs2,msg):
	"""
	"""
	expl_corrs = np.abs(expl_corrs1) - np.abs(expl_corrs2)
	times2save = np.arange(-200, 2000, 20)
	frex = np.logspace(np.log10(2), np.log10(40), 20)
	fig, ax = plt.subplots(1,1)
	cs = ax.contourf(times2save, frex, expl_corrs,  origin='lower')
	fig.colorbar(cs, ax=ax, shrink=0.9)
	#ax.set(title=f'target {channel1} corr over trials from seed {channel2}, {center_freq2}Hz'f'{time_win2[0]}-{time_win2[1]}ms', xlabel='time (ms)', ylabel='frequency (Hz)')
	ax.set(title=msg, xlabel='time (ms)', ylabel='frequency (Hz)')
	fig_name = os.path.join(figures_dir, 'Diff_of_Spectrum_overtrials' + msg +'.png')
	plt.savefig(fig_name,bbox_inches='tight',dpi=400)
	# contourf only above threshold
	

	fig, ax = plt.subplots(1,1)
	threshold = 0.3
	threshold = expl_corrs.mean() + 1.5*expl_corrs.std()
	#mask_th = np.ma.masked_outside(expl_corrs, -threshold, threshold,0)
	#expl_corrs_abs = np.abs(expl_corrs)
	expl_corrs_abs = np.abs(expl_corrs1) - np.abs(expl_corrs2) 
	#mask_th = np.where(np.logical_and(expl_corrs<= -.2, expl_corrs>= .2) ,expl_corrs,0)
	mask_th = np.where(expl_corrs>= threshold ,expl_corrs,0)
	cs2 = ax.contourf(times2save, frex, mask_th,  origin='lower',cmap=plt.cm.Greys_r)
	fig.colorbar(cs2, ax=ax, shrink=0.9)
	ax.set(title=msg, xlabel='time (ms)', ylabel='frequency (Hz)')
	fig_name = os.path.join(figures_dir, 'Diff_of_Spectrum_overtrials '+ 'Mask_' + str(round(threshold,2)) + msg +'.png')
	plt.savefig(fig_name,bbox_inches='tight',dpi=400)
	pdb.set_trace()

def connectivity_analysis(ov_stn, cov_stn, ov_fro, cov_fro, ov_occ, cov_occ):
	"""
	"""
	## Covert alpha frontal beta stn
	center_freq = [10,20]
	print('Correlation between channels Covert cond....\n')
	# Correlation Covert Cortex Frontal vs STN Covert
	fil1_c, fil2_c = connectivity_pair(cov_fro, cov_stn, center_freq)
	# Correlation Overt Cortex Frontal vs STN Overt
	fil1_o, fil2_o = connectivity_pair(ov_fro, ov_stn, center_freq)
	pdb.set_trace()		



	print('Seed Target Power based connectivity ....\n')
	# Seed Target Power based connectivity 
	
	# alpha Frontal vs beta STN Overt
	center_freq_seed=10; center_freq_target = 20
	time_win_seed=(1000,1800); time_win_target=(1000,1800)
	# Target/ Seed
	data1, data2 =  ov_fro,ov_stn# ov_fro
	EEGpnts1, EEGpnts2 = data1.shape[1],data2.shape[1]
	EEGtrials1, EEGtrials2 = data1.shape[0], data2.shape[0]
	EEGpnts, EEGtrials = np.min([EEGpnts1, EEGpnts2]), np.min([EEGtrials1, EEGtrials2])
	EEGtimes = np.linspace(-2000,2000, num=EEGpnts) #t_ms
	# make sure same number of trials both df
	data1_tr, data2_tr = data1[:EEGtrials], data2[:EEGtrials]
	tf_window_data1o, tf_window_data2o,expl_corrs_o1 = power_correlation_trials_aux(data1_tr, data2_tr, data1.name, data2.name,time_win_target,time_win_seed, center_freq_target,center_freq_seed, EEGtrials,EEGtimes)
	pdb.set_trace()	
	# STN theta
	center_freq_target = 5
	tf_window_data1o, tf_window_data2o,expl_corrs_o2 = power_correlation_trials_aux(data1_tr, data2_tr, data1.name, data2.name,time_win_target,time_win_seed, center_freq_target,center_freq_seed,EEGtrials,EEGtimes)
	
	# alpha Frontal vs beta STN COvert
	data1, data2 = cov_stn, cov_fro
	center_freq_seed=9; center_freq_target = 19
	EEGpnts1, EEGpnts2 = data1.shape[1],data2.shape[1]
	EEGtrials1, EEGtrials2 = data1.shape[0], data2.shape[0]
	EEGpnts, EEGtrials = np.min([EEGpnts1, EEGpnts2]), np.min([EEGtrials1, EEGtrials2])
	data1_tr, data2_tr = data1[:EEGtrials], data2[:EEGtrials]
	tf_window_data1c, tf_window_data2c, expl_corrs_c1 = power_correlation_trials_aux(data1_tr, data2_tr, data1.name, data2.name,time_win_target,time_win_seed, center_freq_target,center_freq_seed,EEGtrials,EEGtimes)

	# STN theta
	center_freq_target = 5
	tf_window_data1c, tf_window_data2c,expl_corrs_c2 = power_correlation_trials_aux(data1_tr, data2_tr, data1.name, data2.name,time_win_target,time_win_seed, center_freq_target,center_freq_seed,EEGtrials,EEGtimes)

	msg = 'Ov-Cv_Target_' + data1.name + '_f_' + str(center_freq_target) + ' Seed_' + data2.name + '_f_' + str(center_freq_seed)
	plot_difference_spectra(expl_corrs_o1, expl_corrs_c1, msg)
	pdb.set_trace()
	

def Phase_based_connectivity():
	"""
	"""
	#https://github.com/lyndond/Analyzing_Neural_Time_Series/blob/master/notebooks/chapter19.ipynb





	
###################################################################
def main():

	trials_mean, trials_per_time = group_analysis()

	#ov_stn, cov_stn, ov_fro, cov_fro, ov_occ, cov_occ = trials_mean[0],trials_mean[1],trials_mean[2],trials_mean[3],trials_mean[4],trials_mean[5]
	ov_stn, cov_stn, ov_fro, cov_fro, ov_occ, cov_occ = trials_per_time[0],trials_per_time[1],trials_per_time[2],trials_per_time[3],trials_per_time[4],trials_per_time[5]
	# baseline normalization type = [db, z, pc, div]
	
	# ov_stn_db, ov_stn_z, ov_stn_pc, ov_stn_div = baseline_normalize_channel(ov_stn, 'ov_stn'); cov_stn_db, cov_stn_z, cov_stn_pc, cov_stn_div = baseline_normalize_channel(cov_stn, 'cov_stn')
	# ov_fro_db, ov_fro_z, ov_fro_pc, ov_fro_div = baseline_normalize_channel(ov_fro, 'ov_fro'); cov_fro_db, cov_fro_z, cov_fro_pc, cov_fro_div = baseline_normalize_channel(cov_fro, 'cov_fro')
	# ov_occ_db, ov_occ_z, ov_occ_pc, ov_occ_div = baseline_normalize_channel(ov_occ, 'ov_occ'); cov_occ_db, cov_occ_z, cov_occ_pc, cov_occ_div = baseline_normalize_channel(cov_occ, 'cov_occ')
	
	# Compare conditions 
	# C1 Frontal [8,10] O > C
	tf_ov_fro = baseline_normalize_channel(ov_fro, 'ov_fro')[1]; tf_cov_fro = baseline_normalize_channel(cov_fro, 'cov_fro')[1]# C1 Frontal [8,10] O > C
	compare_power_strectra(tf_ov_fro,tf_cov_fro, 'ov_fro - cov_fro')
	compare_power_strectra(tf_cov_fro,tf_ov_fro, 'cov_fro - ov_fro')
	# C2 STN [17,21,10] [800,1200] O > C
	#tf_ov_stn = baseline_normalize_channel(ov_stn, 'ov_stn')[1]; tf_cov_stn = baseline_normalize_channel(cov_stn, 'cov_stn')[1]# C2 STN  O > C
	#compare_power_strectra(tf_ov_stn,tf_cov_stn, 'ov_stn - cov_stn')
	#compare_power_strectra(tf_cov_stn,tf_ov_stn, 'cov_stn - ov_stn')
	#pdb.set_trace()

	# C3 OCC  C > O
	tf_ov_occ = baseline_normalize_channel(ov_occ, 'ov_occ')[1]; tf_cov_occ = baseline_normalize_channel(cov_occ, 'cov_occ')[1]# C2 STN  O > C
	compare_power_strectra(tf_ov_occ,tf_cov_occ, 'ov_occ - cov_occ')
	compare_power_strectra(tf_cov_occ,tf_ov_occ, 'cov_occ - ov_occ')

	# pdb.set_trace()
	# ov_stn_db.name = 'ov_stn_db'; cov_stn_db.name = 'cov_stn_db'; ov_fro_db.name = 'ov_fro_db'; cov_fro_db.name = 'cov_fro_db'; ov_occ_db.name = 'ov_occ_db';cov_occ_db.name = 'cov_occ_db'
	# ov_stn_z.name = 'ov_stn_z'; cov_stn_z.name = 'cov_stn_z'; ov_fro_z.name = 'ov_fro_z'; cov_fro_z.name = 'cov_fro_z'; ov_occ_z.name = 'ov_occ_z';cov_occ_z.name = 'cov_occ_z'
	# ov_stn_pc.name = 'ov_stn_pc'; cov_stn_pc.name = 'cov_stn_pc'; ov_fro_pc.name = 'ov_fro_pc'; cov_fro_pc.name = 'cov_fro_pc'; ov_occ_pc.name = 'ov_occ_pc';cov_occ_pc.name = 'cov_occ_pc'
	# ov_stn_div.name = 'ov_stn_div'; cov_stn_div.name = 'cov_stn_div'; ov_fro_div.name = 'ov_fro_div'; cov_fro_div.name = 'cov_fro_div'; ov_occ_div.name = 'ov_occ_div';cov_occ_div.name = 'cov_occ_div'
	
	#tf_ov_stn = time_frequency_data(ov_stn, 'ov_stn')
	#tf_ov_stn_dec = baseline_normalization_decibels(tf_ov_stn, 'ov_stn')
	#tf_ov_stn_zeta = baseline_normalization_zeta(tf_ov_stn, 'ov_stn')
	#tf_ov_stn_percent = baseline_normalization_percent(tf_ov_stn, 'ov_stn')
	#tf_ov_stn_baselinediv = baseline_normalization_baselinediv(tf_ov_stn, 'ov_stn')
	#plot_2signals(tf_ov_stn,tf_ov_stn_dec,'ov_stn-Db')
	#plot_2signals(tf_ov_stn,tf_ov_stn_zeta,'ov_stn-Zeta')
	#plot_2signals(tf_ov_stn,tf_ov_stn_percent,'ov_stn-Percent')
	#plot_2signals(tf_ov_stn,tf_ov_stn_baselinediv,'ov_stn-Divide')

	ov_stn.name = 'ov_stn'; cov_stn.name = 'cov_stn'; ov_fro.name = 'ov_fro'; cov_fro.name = 'cov_fro'; ov_occ.name = 'ov_occ';cov_occ.name = 'cov_occ'
	time_pts = 4001
	t_ms = np.linspace(-2000,2000, num=time_pts)

	for channel in [ov_stn, cov_stn, ov_fro, cov_fro, ov_occ, cov_occ]:
		EEGpnts = channel.shape[1]; EEGtrials = channel.shape[0]
		data = channel.to_numpy()
		data = np.reshape(data, (EEGpnts*EEGtrials), order="f")
		convolution_result_fft = plot_power_dist(data,EEGpnts,EEGtrials, channel.name)
		normality_test(convolution_result_fft)

	# Select data raw or normalize
	normalized = ['db', 'z', 'pc', 'div', 'raw']; normalized = normalized[-1]
	if	normalized == 'db':
		ov_stn, cov_stn, ov_fro, cov_fro, ov_occ, cov_occ = ov_stn_db, cov_stn_db, ov_fro_db, cov_fro_db, ov_occ_db, cov_occ_db  
	if	normalized == 'z':
		ov_stn, cov_stn, ov_fro, cov_fro, ov_occ, cov_occ = ov_stn_z, cov_stn_z, ov_fro_z, cov_fro_z, ov_occ_z, cov_occ_z 
	if	normalized == 'pc':
		ov_stn, cov_stn, ov_fro, cov_fro, ov_occ, cov_occ = ov_stn_pc, cov_stn_pc, ov_fro_pc, cov_fro_pc, ov_occ_pc, cov_occ_pc  
	if	normalized == 'div':
		ov_stn, cov_stn, ov_fro, cov_fro, ov_occ, cov_occ = ov_stn_div, cov_stn_div, ov_fro_div, cov_fro_div, ov_occ_div, cov_occ_div  
	else:
		print('Raw Power No normalized')
	print('Calling to connectivity_analysis....\n')

	connectivity_analysis(ov_stn, cov_stn, ov_fro, cov_fro, ov_occ, cov_occ)	




	#######################
	#######################
	sub5_mat_fpath = '/Users/borri/github/EEG_Cadiz_2021/subj5.mat'[df_s1_ov_stn, df_s2_ov_stn, df_s3_ov_stn, df_s4_ov_stn, df_s5_ov_stn, df_s6_ov_stn, df_s7_ov_stn]
	averages, timextrials = load_subj5_data(sub5_mat_fpath) 
	
	occ_cov=averages[0]; occ_ov= averages[1];fro_cov=averages[2];fro_ov=averages[3];stn_cov=averages[4];stn_ov=averages[5];
	df_occ_cov=timextrials[0];df_occ_ov=timextrials[1];df_fro_cov=timextrials[2];df_fro_ov=timextrials[3];df_stn_cov=timextrials[4];df_stn_ov=timextrials[5];
	time_pts = 4001
	t_ms = np.linspace(-2000,2000, num=time_pts) #EEGtimes
	#convolution_result_fft = plot_pwer_dist(occ_cov)
	#normality_test(convolution_result_fft)

	#connectivity_2channels(occ_cov, fro_cov, t_ms, 'occ_cov-fro_cov')
	#connectivity_2channels(stn_cov, fro_cov, t_ms, 'stn_cov-fro_cov')
	#connectivity_2channels(stn_ov, fro_ov, t_ms, 'stn_ov-fro_ov')
	sensor1='occ_cov'; sensor2='fro_cov'; sensor3= 'stn_cov'
	time_win1=(-500,-100); time_win2=(200,1000);
	center_freq1=6;center_freq2=10
	# time freq with seed
	#power_corr_trials(df_occ_cov,df_fro_cov,sensor1,sensor2,time_win1,time_win2,center_freq1,center_freq2,t_ms)
	# time freq all frequencies (data_seed,data_target,seed_chan,target_chan,t_ms)
	power_corr_all_freqs(df_occ_cov,df_fro_cov,sensor1,sensor2,t_ms)
	power_partialcorr_all_freqs(df_occ_cov,df_fro_cov,df_stn_cov,sensor1,sensor2,sensor3, t_ms)

	
if __name__ == "__name__":
	
	main()
