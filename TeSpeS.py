import numpy as np
c = 2.998e8
def window(nt,Tmax,lambda0):   
    dt = 2 * Tmax / nt 
    omega0 = (2 * np.pi * c) / lambda0 
    print(omega0)
    t = (np.arange(-nt/2, nt/2)) * dt
    omega = (np.pi / Tmax) * np.fft.fftshift(np.arange(-nt/2, nt/2))
    lambda_values = (2 * np.pi * c) / (omega + omega0) 
    lambda_nm =np.flip(np.fft.fftshift(lambda_values * 1e9 ))
    print("-" * 50)
    print('Data window:')
    print('dt = ', dt,'[sec]')
    print('Wavelength array:',int(lambda_nm[-1]),'to',int(lambda_nm[0]),'[nm]') 
    print("-" * 50)
    return dt,t,omega,lambda_nm

def filter_gaussian(A,w,lambda0,lambda_c,BW_nm,m_order):
    wc = (2 * np.pi * c) / lambda_c 
    w0 = (2 * np.pi * c) / lambda0 
    Omega=w+w0
    Num=Omega-wc
    Den=-(2*np.pi*c*BW_nm)/(lambda0**2)
    filter=np.exp(-0.5*(Num/Den)**(2*m_order))
    Acplx=(np.fft.ifft(A)*filter)
    return np.fft.fft(Acplx)