import numpy as np

def Pulse(P0,u0,chirp):
    return np.sqrt(P0)*u0*np.exp(1j*chirp)

def Intensity(A0):
    return np.abs(A0)**2

def Spectrum(A0):
    return np.abs(np.fft.fftshift(np.fft.ifft(A0)))**2

def chirp (A,t):
     phi_inst = np.unwrap(np.angle(A))
     chirp = np.gradient(phi_inst, t)
     return chirp

def Energy(A,t):
    return np.trapz(Intensity(A),t)


def FWHM(x, y):
    y = np.asarray(y)
    half_max = 0.5 * np.max(y)

    idx = np.where(y >= half_max)[0]
    if len(idx) < 2:
        return 0.0  # no se puede definir FWHM

    return x[idx[-1]] - x[idx[0]]
