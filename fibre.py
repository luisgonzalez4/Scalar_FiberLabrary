import numpy as np
import pulse as pls
c = 299792
def SMF(A,w,beta2,beta3,gamma,L,step):
     h=L*1e9/step
     D2=0.5j * beta2 * w**2
     D3=1j/6 * beta3 * w**3
     D=np.exp((D2+D3)*h/2)
     for i in range(0,step):
        NL= 1j*gamma*(np.abs(A)**2)
        A=np.fft.fft(np.fft.ifft(A)*D)
        A=A*np.exp(NL*h)
        A=np.fft.fft(np.fft.ifft(A)*D)
     return A

def Coupler(A,alfa):
    A_ro=np.sqrt(alfa)*A
    A_bo=1j*np.sqrt(1-alfa)*A
    return A_ro,A_bo


def DopedFibre(A,t,w,lambda0,BW,beta2,gamma,g0,Esat,L,step):
    h=L*1e9/step
    D2=0.5j * beta2 * w**2
    Omega_g=-((2*np.pi*c*BW)/(lambda0**2) )
    GainShape=1/(1+(w/Omega_g)**2)
    for i in range(step):

        NL = 1j * gamma * np.abs(A)**2

        E = pls.Energy(A, t)
        g_sat = g0 / (1 + E / Esat)  
        Linear = np.exp((D2 + 0.5 * g_sat * GainShape) * h/2)

        A = np.fft.fft(np.fft.ifft(A) * Linear)
        A = A * np.exp(NL * h)


        E = pls.Energy(A, t)
        g_sat = g0 / (1 + E / Esat)  
        Linear = np.exp((D2 + 0.5 * g_sat * GainShape) * h/2)

        A = np.fft.fft(np.fft.ifft(A) * Linear)
    return A

def SA(A,q0,psat):
    T=1-(q0 / (1+ (pls.Intensity(A)/psat)) )
    Aout=T*A
    return Aout