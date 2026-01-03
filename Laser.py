import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import TeSpeS as ts # Temporal-Spectral Space
import pulse as ps # Pulse functions
import fibre as fb # Fibre propagation, SMF, Cuplers, AF, atc...


plt.rc('font', family='')  
TamFuente = 10
                        ###    Time and Spectral array ###
nt = 2**13      
Tmax = 60e-12       # [-tmax/2 tmax/2]
lambda0=1550e-9
dt,t,omega,lambda_nm = ts.window(nt,Tmax,lambda0)
                        ###    Pulse propities   ###
t0=1e-12
phase=np.random.uniform(0,2*np.pi,nt)
P0=np.random.uniform(0,150,nt)
u0=np.exp(-(t**2)/(2*t0**2))
A0=np.sqrt(P0)*u0*np.exp(1j*phase)
#df = pd.read_csv("/Users/luisgonzalez/Desktop/Intracavity-evolution-Mabed/Ain.txt", header=None)
#A0 = np.array(df.iloc[:, 0], dtype=complex)
RT=10
dfspec=np.zeros((RT,nt))
dfTemporal=np.zeros((RT,nt))
I=np.zeros((RT,nt))
E=np.zeros(RT)
print('Laser')
plt.plot(t/t0,ps.Intensity(A0))
plt.show()
plt.ion()
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6, 2.5))
A=A0
for i in range(RT):

    A=fb.DopedFibre(A,t,omega,lambda0,10e-9,4e-26,6e-3,0.73,1e-9,10,2000)
    A=fb.SMF(A,omega,-2.17e-26,8.6e-41,1.1e-3,3,1000)
    A=fb.SMF(A,omega,-2.17e-26,8.6e-41,1.1e-3,4,1000)
    Aout,A=fb.Coupler(A,0.1)
    A=fb.SA(A,0.9,150)
    A=ts.filter_gaussian(A,omega,lambda0,lambda0,10e-9,2)
    A=fb.SMF(A,omega,-2.17e-26,8.6e-41,1.1e-3,8,1000)
    E[i]=ps.Energy(Aout,t)
    print('Round trip: ', i+1)

    ax1.cla()
    dB = 10 * np.log10(ps.Spectrum(Aout) / np.max(ps.Spectrum(Aout)))
    ax1.plot(lambda_nm, dB  , color='black', linewidth=0.5)
    ax1.set_xlabel('Wavelength, nm')
    ax1.set_ylabel('Spectral Intensity, dB')
    ax1.set_title(f'Round Trip {i+1}')
    #ax1.grid(True)
    #ax1.set_xlim([1350, 1840])
    #ax1.set_ylim([-50, 0])

    ax2.cla()
    ax2.plot(t / t0, ps.Intensity(Aout), color='black', linewidth=0.5)
    #ax2.set_xlim([-20,20])
    ax2.set_ylim([0,max(ps.Intensity(Aout))])
    ax2.set_xlabel('Time, ps')
    ax2.set_ylabel('Intensity, W')
    #ax2.grid(True)
    plt.tight_layout()
    plt.pause(0.1)
    plt.savefig(f"Laser-Estable_{i+1}.png", dpi=300, bbox_inches='tight') 
np.savetxt("/Users/luisgonzalez/Desktop/Ain.txt",A, delimiter=",")

plt.ioff()
plt.show()

plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1) 

dB=10*np.log10(ps.Spectrum(Aout)/max(ps.Spectrum(Aout)))
plt.plot(lambda_nm, dB , color='black', linewidth=0.5)
plt.xlim([1400,1700])
plt.ylim([-150,0])
plt.xlabel('Wavelength, nm')
plt.ylabel('Spectral Intensity, dB')

plt.subplot(1, 2, 2) 
plt.plot(t / t0, ps.Intensity(Aout), color='black', linewidth=0.5)
plt.xlim([-20,20])
plt.ylim([-1,max(ps.Intensity(Aout))])
plt.xlabel('Time, ps')
plt.ylabel('Intensity, W')
#plt.grid(True)
plt.tight_layout()
plt.show()

plt.plot(np.arange(1,RT+1,1),E*1e9,'ko-')
plt.xlabel('Rount-Trip')
plt.ylabel('Energy, nJ')
plt.show()