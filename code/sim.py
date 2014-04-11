from __future__ import division
import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/tausen/p8/git/code/')
import basics


N = int(3e3)
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

if fs == 0:
    tn, x, phix = getinput()

u, v, phiy, s = loop(x, tn)
