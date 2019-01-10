from __future__ import division
from numpy.fft import rfft
from numpy import argmax, mean, diff, log
from matplotlib.mlab import find
#from scipy.signal import blackmanharris, fftconvolve
from time import time
import sys
import wave
import struct
import numpy as np
import matplotlib.pyplot as plt
import copy
#try:
#    import soundfile as sf
#except ImportError:
#    from scikits.audiolab import flacread


#from parabolic import parabolic

def blackmanharris(n):
	#return np.array( [1.0 for i in xrange(n)] )
	c = [0.35875, 0.48829, 0.14128, 0.01168]
	duepi_div = (2.0*np.pi)/(n-1)
	w = []
	for i in xrange(n):
		ph = duepi_div*i
		x = c[0] - c[1]*np.cos(ph) + c[2]*np.cos(2.0*ph) - c[3]*np.cos(3.0*ph)
		w.append(x)
	#return np.ones(n)
	return w

def freq_from_crossings(sig, fs):
    """
    Estimate frequency by counting zero crossings
    """
    # Find all indices right before a rising-edge zero crossing
    indices = find((sig[1:] >= 0) & (sig[:-1] < 0))

    # Naive (Measures 1000.185 Hz for 1000 Hz, for instance)
    # crossings = indices

    # More accurate, using linear interpolation to find intersample
    # zero-crossings (Measures 1000.000129 Hz for 1000 Hz, for instance)
    crossings = [i - sig[i] / (sig[i+1] - sig[i]) for i in indices]

    # Some other interpolation based on neighboring points might be better.
    # Spline, cubic, whatever

    return fs / mean(diff(crossings))


def freq_from_fft(sig, fs):
    """
    Estimate frequency from peak of FFT
    """
    # Compute Fourier transform of windowed signal
    windowed = sig * blackmanharris(len(sig))
    f = rfft(windowed)

    # Find the peak and interpolate to get a more accurate peak
    i = argmax(abs(f))  # Just use this for less-accurate, naive version
    #true_i = parabolic(log(abs(f)), i)[0]
    true_i = i

    # Convert to equivalent frequency
    return fs * true_i / len(windowed)


def freq_from_autocorr(sig, fs):
    """
    Estimate frequency using autocorrelation
    """
    # Calculate autocorrelation (same thing as convolution, but with
    # one input reversed in time), and throw away the negative lags
    #corr = fftconvolve(sig, sig[::-1], mode='full')
    corr = corr[len(corr)//2:]

    # Find the first low point
    d = diff(corr)
    start = find(d > 0)[0]

    # Find the next peak after the low point (other than 0 lag).  This bit is
    # not reliable for long signals, due to the desired peak occurring between
    # samples, and other peaks appearing higher.
    # Should use a weighting function to de-emphasize the peaks at longer lags.
    peak = argmax(corr[start:]) + start
    px, py = parabolic(corr, peak)

    return fs / px


def freq_from_HPS(sig, fs):
    """
    Estimate frequency using harmonic product spectrum (HPS)

    """
    windowed = sig * blackmanharris(len(sig))

    #from pylab import subplot, plot, log, copy, show

    # harmonic product spectrum:
    c = abs(rfft(windowed))
    maxharms = 8
    plt.subplot(maxharms, 1, 1)
    plt.plot(log(c))
    for x in range(2, maxharms):
        a = copy.copy(c[::x])  # Should average or maximum instead of decimating
        # max(c[::x],c[1::x],c[2::x],...)
        c = c[:len(a)]
        i = argmax(abs(c))
        #true_i = parabolic(abs(c), i)[0]
        true_i = i
        print 'Pass %d: %f Hz' % (x, fs * true_i / len(windowed))
        c *= a
        plt.subplot(maxharms, 1, x)
        plt.plot(log(c))
    plt.show()

def music_frequencies():
	r = pow(2, 1./12.)
	s = f = 27.5
	t = [f]
	l = ['A', 'A#', 'B', 'C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#']
	i = 0
	n = 9
	lab = [l[0] + "0"]
	while f < 8000.:
		f *= r
		t.append(f)
		i += 1
		n += 1
		lab.append( l[i%12] + str(int(n/12)) )
	return t, lab
	
filename = sys.argv[1]

print 'Reading file "%s"\n' % filename
f = wave.open(filename,'r')
#signal = f.readframes(f.getnframes())
bps = f.getsampwidth()
nfr = max(f.getnframes(), 8192)
nch = f.getnchannels()
fs = f.getframerate()
print nfr, bps, nch, fs
tmps = f.readframes(nfr)
framesize = bps*nch
tmp = []
tmp2 = []
offset = 0
if bps == 1:
	for k in xrange(nfr):
		tmp.append( struct.unpack('b',tmps[offset:offset+bps])[0] )
		if nch == 2:
			tmp2.append( struct.unpack('b',tmps[offset+bps:offset+2*bps])[0] )
		offset += framesize
elif bps == 2:
	for k in xrange(nfr):
		tmp.append( struct.unpack('h',tmps[offset:offset+bps])[0] )
		if nch == 2:
			tmp2.append( struct.unpack('h',tmps[offset+bps:offset+2*bps])[0] )
		offset += framesize
totsignal = np.array(tmp) / 32767.0
totsignal2 = np.array(tmp2) / 32767.0
print len(totsignal), len(totsignal2)

f.close()

#plot
#plt.subplot(211)
#plt.plot(totsignal)
#plt.subplot(212)
#plt.plot(totsignal2)
#plt.show()
#sys.exit(0)

mfrq,mlab=music_frequencies()
#print mfrq

#spectrogram
#winsize = 16384
#ovlp=winsize-8820
mycm='gist_heat' # COLORMAP
# alcune colormap disponibili 'gray', 'gist_gray', 'gist_yarg', 'Spectral', 'seismic', 'cool'
# 'bone', 'BuGn', 'gist_heat'
# None per usare la colormap di default
winsize = 8192
ovlp=winsize-128
if nch == 2:
	sub1 = plt.subplot(211)
else:
	sub1 = plt.subplot()	
sub1.specgram(totsignal, NFFT=winsize, noverlap=ovlp, Fs=fs, scale='dB', mode='psd', cmap=mycm)
sub1.set_yscale('symlog')
sub1.set_yticks(mfrq)
sub1.set_yticklabels(mlab)
sub1.grid(True, axis='y')
sub1.set_ylim([25., 8000.])
if nch == 2:
	sub2 = plt.subplot(212, sharey=sub1, sharex=sub1)
	sub2.specgram(totsignal2, NFFT=winsize, noverlap=ovlp, Fs=fs, scale='dB', mode='psd', cmap=mycm)
	sub2.set_yscale('symlog')
	sub2.set_yticks(mfrq)
	sub2.set_yticklabels(mlab)
	sub2.grid(True, axis='y')
	sub2.set_ylim([25., 8000.])
plt.subplots_adjust(top=0.95, bottom=0.05, hspace=0.10)	
plt.show()
sys.exit(0)

print 'Calculating frequency from FFT:',
start_time = time()
print '%f Hz' % freq_from_fft(signal, fs)
print 'Time elapsed: %.3f s\n' % (time() - start_time)

print 'Calculating frequency from zero crossings:',
start_time = time()
print '%f Hz' % freq_from_crossings(signal, fs)
print 'Time elapsed: %.3f s\n' % (time() - start_time)

#print 'Calculating frequency from autocorrelation:',
#start_time = time()
#print '%f Hz' % freq_from_autocorr(signal, fs)
#print 'Time elapsed: %.3f s\n' % (time() - start_time)

print 'Calculating frequency from harmonic product spectrum:'
start_time = time()
freq_from_HPS(signal, fs)
print 'Time elapsed: %.3f s\n' % (time() - start_time)
