import matplotlib.pyplot as plt
import numpy as np
import cmath
from pandas import read_csv
from scipy import signal

from matplotlib.ticker import FormatStrFormatter

import warnings

def LCL1(Li, Cf, Lg, freq):
    
    omega = np.zeros(len(freq), dtype=complex)
    
    det_s1 = np.zeros(len(freq), dtype=complex)
    det_s0 = np.zeros(len(freq), dtype=complex)
    
    num_s3 = np.zeros(len(freq), dtype=complex)
    num_s2 = np.zeros(len(freq), dtype=complex)
    num_s1 = np.zeros(len(freq), dtype=complex)
    num_s0 = np.zeros(len(freq), dtype=complex)

    y = np.zeros(len(freq), dtype=complex)
    y_phase = np.zeros(len(freq), dtype=float)
    y_gain = np.zeros(len(freq), dtype=float)
    
    
    omega = freq*2.0*np.pi*1.0j


    det_s0 = 1

    num_s3 = ((Lg*Li*Cf)*np.power(omega, 3))
    num_s2 = 0
    num_s1 = ((Li + Lg)*omega)
    num_s0 = 0
    

    y = ((det_s0) / (num_s3 + num_s2 + num_s1 + num_s0))
    for i in range(0, len(y)):
        y_phase[i] = cmath.phase(y[i])
        
    y_phase = np.unwrap(y_phase) * (180/np.pi)
    y_phase[y_phase > -90]=-270
         
    y_gain = 20*np.log10((np.abs(y)))
    
    return [y_gain, y_phase]
def LCL2(Li, Cf, Rc, Lg, freq):
    
    omega = np.zeros(len(freq), dtype=complex)
    
    det_s1 = np.zeros(len(freq), dtype=complex)
    det_s0 = np.zeros(len(freq), dtype=complex)
    
    num_s3 = np.zeros(len(freq), dtype=complex)
    num_s2 = np.zeros(len(freq), dtype=complex)
    num_s1 = np.zeros(len(freq), dtype=complex)
    num_s0 = np.zeros(len(freq), dtype=complex)

    y = np.zeros(len(freq), dtype=complex)
    y_phase = np.zeros(len(freq), dtype=float)
    y_gain = np.zeros(len(freq), dtype=float)
    
    
    omega = freq*2.0*np.pi*1.0j

    det_s1 = Rc*Cf*omega
    det_s0 = 1

    num_s3 = ((Lg*Li*Cf)*np.power(omega, 3))
    num_s2 = (((Li + Lg)*Cf*Rc)*np.power(omega, 2))
    num_s1 = ((Li + Lg)*omega)
    num_s0 = 0
    

    y = (det_s1 + det_s0) / (num_s3 + num_s2 + num_s1 + num_s0)
    for i in range(0, len(y)):
        y_phase[i] = cmath.phase(y[i])
        
    y_phase = np.unwrap(y_phase) * (180/np.pi)
         
    y_gain = 20*np.log10((np.abs(y)))
    
    return [y_gain, y_phase]
def LCL3(Li, Ri, Cf, Rc, Lg, Rg, freq):
    
    omega = np.zeros(len(freq), dtype=complex)
    
    det_s1 = np.zeros(len(freq), dtype=complex)
    det_s0 = np.zeros(len(freq), dtype=complex)
    
    num_s3 = np.zeros(len(freq), dtype=complex)
    num_s2 = np.zeros(len(freq), dtype=complex)
    num_s1 = np.zeros(len(freq), dtype=complex)
    num_s0 = np.zeros(len(freq), dtype=complex)

    y = np.zeros(len(freq), dtype=complex)
    y_phase = np.zeros(len(freq), dtype=float)
    y_gain = np.zeros(len(freq), dtype=float)
    
    

    omega = freq*2.0*np.pi*1.0j

    det_s1 = Rc*Cf*omega
    det_s0 = 1

    num_s3 = ((Lg*Li*Cf)*np.power(omega, 3))
    num_s2 = ((Cf*(Lg*(Rc + Ri) + Li*(Rc + Rg)))*np.power(omega, 2))
    num_s1 = ((Lg + Li + Cf*(Rc*Rg + Rc*Ri + Rg*Ri))*omega)
    num_s0 = (Rg + Ri)
    

    y = (det_s1 + det_s0) / (num_s3 + num_s2 + num_s1 + num_s0)
    for i in range(0, len(y)):
        y_phase[i] = cmath.phase(y[i])
        
    y_phase = np.unwrap(y_phase) * (180/np.pi)
         
    y_gain = 20*np.log10((np.abs(y)))
    
    return [y_gain, y_phase]

def LCL3_load(Li, Ri, Cf, Rc, Lg, Rg, Rl, Ll, freq):
    
    omega = np.zeros(len(freq), dtype=complex)
    
    det_s1 = np.zeros(len(freq), dtype=complex)
    det_s0 = np.zeros(len(freq), dtype=complex)
    
    num_s3 = np.zeros(len(freq), dtype=complex)
    num_s2 = np.zeros(len(freq), dtype=complex)
    num_s1 = np.zeros(len(freq), dtype=complex)
    num_s0 = np.zeros(len(freq), dtype=complex)
    
    load = np.zeros(len(freq), dtype=complex)

    y = np.zeros(len(freq), dtype=complex)
    y_phase = np.zeros(len(freq), dtype=float)
    y_gain = np.zeros(len(freq), dtype=float)
    
    

    omega = freq*2.0*np.pi*1.0j

    det_s1 = Rc*Cf*omega
    det_s0 = 1

    num_s3 = (((Lg+Ll)*Li*Cf)*np.power(omega, 3))
    num_s2 = ((Cf*((Lg+Ll)*(Rc + Ri) + Li*(Rc + (Rg+Rl))))*np.power(omega, 2))
    num_s1 = (((Lg+Ll) + Li + Cf*(Rc*(Rg+Rl) + Rc*Ri + (Rg+Rl)*Ri))*omega)
    num_s0 = ((Rg+Rl) + Ri)
    
    load = ((Rl + Ll * omega))

    y = ((det_s1 + det_s0) / (num_s3 + num_s2 + num_s1 + num_s0))
    for i in range(0, len(y)):
        y_phase[i] = cmath.phase(y[i])
        
    y_phase = np.unwrap(y_phase) * (180/np.pi)
         
    y_gain = 20*np.log10((np.abs(y)))
    
    return [y_gain, y_phase]

#df = read_csv('LCL1Ohm_1.CSV')
#df = read_csv('LCL_2coil_1Ohm_1.CSV')

n = 10

Cf = 30e-6
Rc = Rf = 220e-3

Li = 1e-3
Ri = 300e-3

Lg = 47e-6
Rg = 45.6e-3

Ll = 4.33e-3
Rl = 1 + 0.2


freq = np.linspace(10, 15e3, int(10e3))
    
y_phase_LCL1 = np.zeros(len(freq), dtype=float)
y_gain_LCL1 = np.zeros(len(freq), dtype=float)
y_phase_LCL2 = np.zeros(len(freq), dtype=float)
y_gain_LCL2 = np.zeros(len(freq), dtype=float)
y_phase_LCL3 = np.zeros(len(freq), dtype=float)
y_gain_LCL3 = np.zeros(len(freq), dtype=float)
y_phase_LCL3_load = np.zeros(len(freq), dtype=float)
y_gain_LCL3_load = np.zeros(len(freq), dtype=float)

[y_gain_LCL1, y_phase_LCL1] = LCL1(Li, Cf, Lg, freq)
[y_gain_LCL2, y_phase_LCL2] = LCL2(Li, Cf, Rc, Lg, freq)
[y_gain_LCL3, y_phase_LCL3] = LCL3(Li, Ri, Cf, Rc, Lg, Rg, freq)
[y_gain_LCL3_load, y_phase_LCL3_load] = LCL3_load(Li, Ri, Cf, Rc, Lg, Rg, Rl, Ll, freq)

fig1, axs1 = plt.subplots(2)

axs1[0].semilogx(freq, y_gain_LCL1, label=r'$H\left(s\right)_{LCL_{1}}$')     # Bode magnitude plot
axs1[0].semilogx(freq, y_gain_LCL2, label=r'$H\left(s\right)_{LCL_{2}}$') 
axs1[0].semilogx(freq, y_gain_LCL3, label=r'$H\left(s\right)_{LCL_{3}}$')
#axs1[0].semilogx(freq, y_gain_LCL3_load, label=r'$H\left(s\right)_{LCL_{3}}$ with load')
#axs1[0].semilogx(df['Freq'].rolling(window =n).mean(), 20*np.log10(df['Ia'].rolling(window =n).mean()/df['Ean'].rolling(window =n).mean()), color='magenta', label='Measured')
axs1[0].set_ylabel('Gain [dB]')
axs1[0].set_xlim(200, 15e3)
#axs1[0].set_ylim(-100, 20)
axs1[0].grid(b=None, which='major', axis='both')
axs1[0].grid(b=None, which='minor', axis='both')
axs1[0].legend()

axs1[1].semilogx(freq, y_phase_LCL1, label=r'$H\left(s\right)_{LCL_{1}}$')    # Bode magnitude plot
axs1[1].semilogx(freq, y_phase_LCL2, label=r'$H\left(s\right)_{LCL_{2}}$') 
axs1[1].semilogx(freq, y_phase_LCL3, label=r'$H\left(s\right)_{LCL_{3}}$')
#axs1[1].semilogx(freq, y_phase_LCL3_load, label=r'$H\left(s\right)_{LCL_{3}}$ with load')
axs1[1].set_xlabel('Frequency [Hz]')
axs1[1].set_ylabel('Phase [Degree]')
axs1[1].set_xlim(200, 15e3)
#axs1[1].set_ylim(-251, -79)
axs1[1].grid(b=None, which='major', axis='both')
axs1[1].grid(b=None, which='minor', axis='both')
axs1[1].legend()
