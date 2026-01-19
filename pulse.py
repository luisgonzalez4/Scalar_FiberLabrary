import numpy as np
def Pulse(P0,u0,chirp):
    return np.sqrt(P0)*u0*np.exp(-1j*chirp)

def Intensity(A0):
    return np.abs(A0)**2

def Spectrum(A0):
        return np.abs(np.fft.fftshift(np.fft.ifft(A0)))**2

def chirp (t,A):
     phi_inst = np.unwrap(np.angle(A))
     chirp= np.gradient(phi_inst, t)
     return chirp

def Energy(A,t):
    return np.trapz(Intensity(A),t)


def FHWM(x, y):
    half_max = max(y) / 2
    index = np.where(y >= half_max)[0]
    xi = x[index[0]]
    xf = x[index[-1]]
    FWHM = xf - xi 
    return FWHM
