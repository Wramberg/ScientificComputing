from __future__ import division
import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/tausen/p8/git/code/')

params = {'legend.fontsize': 15,
          'legend.linewidth': 2}
plt.rcParams.update(params)

N = int(3.2e3)
fif = 800e3
f0 = 1.9e6
kco = 500e3
fs = 0


def getinput():
    global N, fif, fs

    I = np.genfromtxt("gmsk_I.dat")[:N]
    Q = np.genfromtxt("gmsk_Q.dat")[:N]
    os = np.genfromtxt("gmsk_os.dat")
    fsymb = np.genfromtxt("gmsk_fsymb.dat")

    fs = fsymb * os
    tn = np.arange(len(I)) / fs
    x = I * np.cos(2 * np.pi * fif * tn) - Q * np.sin(2 * np.pi * fif * tn)

    phix = np.unwrap(np.arctan2(Q, I))

    return tn, x, phix


def gety(tn, phiy):
    global f0
    y = np.cos(2 * np.pi * f0 * tn + phiy)
    return y


def loop(x, tn):
    global fs, f0, fif, kco, N

    # Values of s[-1], u[-1] and v[-1]
    sm1 = 1
    um1 = 0
    vm1 = 0

    # Get coefficients of 1st order Butterworth lowpass filter
    b, a = signal.butter(1, 200e3 / (fs / 2), 'low')

    # Initialize arrays for signals
    u    = np.empty(N)
    s    = np.empty(N)
    v    = np.empty(N)
    phiy = np.empty(N)

    # Calculate signals for n = 0
    u[0]    = x[0] * sm1
    v[0]    = u[0] * b[0] + um1 * b[1] - vm1 * a[1]
    phiy[0] = 2 * np.pi * kco * 1 / fs * np.sum(v[:1])
    s[0]    = np.cos(2 * np.pi * fif * tn[0] + phiy[0])

    # Main loop
    for n in range(1, N):
        u[n]    = x[n] * s[n - 1]
        v[n]    = u[n] * b[0] + u[n - 1] * b[1] - v[n - 1] * a[1]
        phiy[n] = 2 * np.pi * kco * 1 / fs * np.sum(v[:n + 1])
        s[n]    = np.cos(2 * np.pi * fif * tn[n] + phiy[n])

    return u, v, phiy, s


def plotphases(phix, phiy, phiycomp=None):
    plt.figure(figsize=(11,2.75), dpi=200)   
    plt.plot(phix, label='Phase of input signal')
    plt.plot(phiy, label='Phase of output signal')
    if phiycomp is not None:
        plt.plot(phiycomp, label='Compensated phase of output signal')
    plt.grid()
    plt.xlim(0, N)
    plt.xlabel('Samples')
    plt.legend(loc=3, prop={'size': 12}, labelspacing=0.25)
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.175)
    if phiycomp is not None:
        plt.savefig('phasesignalscomp.eps') 
    else:
        plt.savefig('phasesignals.eps') 


def getphaseerror(phix, phiy):
    try:
        diff = phix[:phiy.size] - phiy
    except:
        raise Exception('phiy.size must be < than phix.size')
    err = np.sqrt(np.mean(diff[250:] ** 2))
    return diff, err


def phasecomp(phix, phiy):
    # Find phase offset
    diff = phix[250:] - phiy[250:]
    poff = np.mean(diff)

    # Find time offset
    maxoffset = 50
    rmserr = np.empty(maxoffset)
    for i in np.arange(0, maxoffset):
        dummy, rmserr[i] = getphaseerror(phix[200:], phiy[200+maxoffset-i:]+poff)

    toff = maxoffset - np.argmin(rmserr)
    if toff < 0:
        raise Exception('Got negative time offset, time delay should be positive')
    return phiy[toff:] + poff


def getperiodogram(x):
    global fs
    freqs, psd = signal.welch(x, fs=fs, window='boxcar', nperseg=300, noverlap=None, detrend='constant', return_onesided=True)
    return freqs, psd


def plotpsd(x, y):
    freqsx, psdx = getperiodogram(x)
    freqsy, psdy = getperiodogram(y)
    plt.figure(figsize=(11,2.75), dpi=200)  
    plt.plot(freqsx/1e3, psdx, label='PSD of input signal')
    plt.plot(freqsy/1e3, psdy, label='PSD of output signal')
    plt.grid()
    plt.xlim(0, 3.5e3)
    plt.xlabel('Frequency [kHz]')
    plt.legend(loc=1, prop={'size': 12}, labelspacing=0.25)
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.175)
    plt.savefig('psds.eps') 


def plotphaseerr(x):
    plt.figure(figsize=(11,2.75), dpi=200)  
    plt.plot(x)
    plt.grid()
    plt.xlabel('Samples')
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.175)
    plt.savefig('phaseerr.eps')  
  

tn, x, phix = getinput()
u, v, phiy, s = loop(x, tn)
y = gety(tn, phiy)
phiycomp = phasecomp(phix, phiy)

plotphases(phix, phiy, phiycomp)
plotpsd(x, y)

diff, rms = getphaseerror(phix, phiycomp)
plotphaseerr(diff)

#outfreq = f0 + kco * v
#plt.figure()
#plt.figure(figsize=(11,2.75), dpi=200)  
#plt.plot(outfreq/1e6)
#plt.grid()
#plt.ylabel('Frequency [MHz]')
#plt.xlabel('Samples')
#plt.subplots_adjust(left=0.1, right=0.9, bottom=0.175)
#plt.savefig('outfreq.eps')  

#plt.show()
