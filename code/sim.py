from __future__ import division
import numpy as np
# import matplotlib.pyplot as plt
import scipy.signal as signal


N = int(3e3)
fif = 800e3
f0 = 1.9e6
kco = 500e3
fs = 0

def nco(x, fs):
    global f0, kco

    N = x.size
    y = np.zeros(N)
    for n in range(N):
        y[n] = np.sin(2 * np.pi * f0 * n / fs + 2 * np.pi * kco * 1 / fs * np.sum(x[0:n]))

    t = np.arange(0, N) / fs
    return t, y


def loadfiles():
    global N, fif, fs

    I = np.genfromtxt("gmsk_I.dat")[:N]
    Q = np.genfromtxt("gmsk_Q.dat")[:N]
    os = np.genfromtxt("gmsk_os.dat")
    fsymb = np.genfromtxt("gmsk_fsymb.dat")

    fs = fsymb * os

    tn = np.arange(len(I)) / fs
    x = I * np.cos(2 * np.pi * fif * tn) - Q * np.sin(2 * np.pi * fif * tn)

    phix = np.unwrap(np.arctan2(I, Q))

    return tn, x, phix


def getz(t):
    global f0, fif

    z = 2 * np.cos(2 * np.pi * (f0 - fif) * t)

    return z


def gety(x, z, tn):
    global fs, f0, fif, kco
    
    # n = -1
    sm1 = 0
    um1 = 0
    vm1 = 0

    b, a = signal.butter(1, 200e3 / (fs / 2), 'low')

    N = len(x)
    u = np.empty(N)
    s = np.empty(N)
    v = np.empty(N)
    y = np.empty(N)
    phiy = np.empty(N)
    for n in range(N):
        if n == 0:
            u[n] = x[n] * sm1
            v[n] = u[n] * b[0] + um1 * b[1] + vm1 * a[1]
            y[n] = 0
        else:
            u[n] = x[n] * s[n - 1]
            v[n] = u[n] * b[0] + u[n - 1] * b[1] + v[n - 1] * a[1]
            phiy[n] = 2 * np.pi * kco * 1 / fs * sum(v[:n])
            y[n] = np.cos(2 * np.pi * f0 * tn[n] + phiy[n])

    return u, v


tn, x, phix = loadfiles()
z = getz(tn)
u, v = gety(x, z)
