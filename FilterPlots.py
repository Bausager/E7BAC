# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 10:22:27 2022

@author: Bausa
"""


import matplotlib.pyplot as plt
import numpy as np
import cmath
from pandas import read_csv

from matplotlib.ticker import FormatStrFormatter


df1 = read_csv('1 Order.csv')
df2 = read_csv('2 Order.csv')
df3 = read_csv('3 Order.csv')



def firstOrder(R, C, freq):
    omega = freq*2.0*np.pi
    
    omega_0 = 1/(R*C)
    out = np.zeros(1, dtype=complex)
    out = 1.0 / (1 + 1j*(omega/omega_0))
    
    return out

def secondOrder(K, R1, C1, R2, C2, freq):
    
    omega = freq*2.0*np.pi
    
    x1 = np.zeros(1, dtype=complex)
    x1 = (omega**2)*R1*C1*R2*C2
    
    x2 = np.zeros(1, dtype=complex)
    x2 = 1.0j*omega*((1.0 - K)*R1*C1 + R1*C2 + R2*C2)
    
    out = np.zeros(1, dtype=complex)
    out = (1.0/(1.0 - x1 + x2)) * K
    
    return out

def thirdOrder(K, R, C, R1, C1, R2, C2, freq):
    omega = freq*2.0*np.pi
    omega_0 = 1/(R*C)
    
    firstOrder = np.zeros(1, dtype=complex)
    firstOrder = 1.0 / (1 + 1j*(omega/omega_0))
    
    x1 = np.zeros(1, dtype=complex)
    x1 = (omega**2)*R1*C1*R2*C2
    
    x2 = np.zeros(1, dtype=complex)
    x2 = 1.0j*omega*((1.0 - K)*R1*C1 + R1*C2 + R2*C2)
    
    out = np.zeros(1, dtype=complex)
    out = ((1.0/(1.0 - x1 + x2)) * K) * firstOrder
    
    return out


# Tolerances
Cap_tolerance = 0.1
Res_tolerance = 0.01

# First Order Filter
C = 68e-9
R = 470

# Second Order Filter
C1 = 200e-9
C2 = 100e-9
R1 = 220
R2 = 220
K = 1

m = 16
M = 9


N = int(1e3)
freq = np.linspace(10, 10e3, N)


y = np.zeros(N, dtype=complex)
y_phase = np.zeros(N, dtype=float)
y_gain = np.zeros(N, dtype=float)
y = firstOrder(R, C, freq)
for i in range(0, len(y)):
    y_phase[i] = cmath.phase(y[i]) * (180/np.pi)
y_gain = 20*np.log10((np.abs(y)))


y_max = np.zeros(N, dtype=complex)
y_phase_max = np.zeros(N, dtype=float)
y_gain_max = np.zeros(N, dtype=float)
y_max = firstOrder(R*(1 - (Res_tolerance*0.5)), C*(1 - (Cap_tolerance*0.5)), freq)
for i in range(0, len(y)):
    y_phase_max[i] = cmath.phase(y_max[i]) * (180/np.pi)
y_gain_max = 20*np.log10((np.abs(y_max)))


y_min = np.zeros(N, dtype=complex)
y_phase_min = np.zeros(N, dtype=float)
y_gain_min = np.zeros(N, dtype=float)
y_min = firstOrder(R*(1 + (Res_tolerance*0.5)), C*(1 + (Cap_tolerance*0.5)), freq)
for i in range(0, len(y)):
    y_phase_min[i] = cmath.phase(y_min[i]) * (180/np.pi)
y_gain_min = 20*np.log10((np.abs(y_min)))



fig1, axs1 = plt.subplots(2, 2, figsize=(m, M))

axs1[0, 0].semilogx(freq, y_gain, color='black')
axs1[0, 0].semilogx(freq, y_gain_max, color='black', linestyle='--')
axs1[0, 0].semilogx(freq, y_gain_min, color='black', linestyle='-.')
axs1[0, 0].semilogx(df1['#Frequency (Hz)'], df1['Channel 2 Magnitude (dB)'], color='red')
axs1[0, 0].set_ylabel('Gain [dB]')
axs1[0, 0].title.set_text('Bodeplot Gain')
axs1[0, 0].set_ylim([-5, 0.1])
axs1[0, 0].grid()

axs1[0, 1].plot(freq, y_gain, color='black', label='No tolerance')
axs1[0, 1].plot(freq, y_gain_max, color='black', linestyle='--', label='Max Gain')
axs1[0, 1].plot(freq, y_gain_min, color='black', linestyle='-.', label='Min Gain')
axs1[0, 1].plot(df1['#Frequency (Hz)'], df1['Channel 2 Magnitude (dB)'], color='red', label='Measured')
axs1[0, 1].title.set_text('Bodeplot Gain (Zoomed)')
axs1[0, 1].legend(loc=1)
axs1[0, 1].set_xlim([2.5e3, 7.5e3])
axs1[0, 1].set_ylim([-5, 0.1])
axs1[0, 1].grid()



axs1[1, 0].semilogx(freq, y_phase, color='black')
axs1[1, 0].semilogx(freq, y_phase_max, color='black', linestyle='--')
axs1[1, 0].semilogx(freq, y_phase_min, color='black', linestyle='-.')
axs1[1, 0].semilogx(df1['#Frequency (Hz)'], df1['Channel 2 Phase (deg)'], color='red')
axs1[1, 0].set_ylabel('Phase [Degree]')
axs1[1, 0].set_xlabel('Frequency [Hz]')
axs1[1, 0].title.set_text('Bodeplot Phase-Shift')
axs1[1, 0].grid()

axs1[1, 1].plot(freq, y_phase, color='black', label='No tolerance')
axs1[1, 1].plot(freq, y_phase_max, color='black', linestyle='--', label='Max Phase-Shift')
axs1[1, 1].plot(freq, y_phase_min, color='black', linestyle='-.', label='Min Phase-Shift')
axs1[1, 1].plot(df1['#Frequency (Hz)'], df1['Channel 2 Phase (deg)'], color='red', label='Measured')
axs1[1, 1].title.set_text('Bodeplot Phase-Shift (Zoomed)')
axs1[1, 1].set_xlabel('Frequency [Hz]')
axs1[1, 1].legend(loc=1)
axs1[1, 1].set_xlim([2.5e3, 7.5e3])
axs1[1, 1].set_ylim([-60, -30])
axs1[1, 1].grid()
fig1.tight_layout()




y = np.zeros(N, dtype=complex)
y_ph = np.zeros(N, dtype=float)
y_gain = np.zeros(N, dtype=float)
y = secondOrder(K, R1, C1, R2, C2, freq)
for i in range(0, len(y)):
    y_ph[i] = cmath.phase(y[i]) * (180/np.pi)
y_gain = 20*np.log10((np.abs(y)))


y_max = np.zeros(N, dtype=complex)
y_ph_max = np.zeros(N, dtype=float)
y_gain_max = np.zeros(N, dtype=float)
y_max = secondOrder(K,
                         R1*(1 + (Res_tolerance*0.5)), 
                             C1*(1 + (Cap_tolerance*0.5)), 
                                 R2*(1 - (Res_tolerance*0.5)), 
                                     C2*(1 - (Cap_tolerance*0.5)), freq)
for i in range(0, len(y_max)):
    y_ph_max[i] = cmath.phase(y_max[i]) * (180/np.pi)
y_gain_max = 20*np.log10((np.abs(y_max)))  


y_min = np.zeros(N, dtype=complex)
y_ph_min = np.zeros(N, dtype=float)
y_gain_min = np.zeros(N, dtype=float)
y_min = secondOrder(K,
                         R1*(1 - (Res_tolerance*0.5)), 
                             C1*(1 - (Cap_tolerance*0.5)), 
                                 R2*(1 + (Res_tolerance*0.5)), 
                                     C2*(1 + (Cap_tolerance*0.5)), freq)
for i in range(0, len(y_max)):
    y_ph_min[i] = cmath.phase(y_min[i])  * (180/np.pi)
y_gain_min = 20*np.log10((np.abs(y_min))) 


fig2, axs2 = plt.subplots(2, 2, figsize=(m, M))

axs2[0, 0].semilogx(freq, y_gain, color='black')
axs2[0, 0].semilogx(freq, y_gain_max, color='black', linestyle='--')
axs2[0, 0].semilogx(freq, y_gain_min, color='black', linestyle='-.')
axs2[0, 0].semilogx(df2['#Frequency (Hz)'], df2['Channel 2 Magnitude (dB)'], color='red')
axs2[0, 0].set_ylabel('Gain [dB]')
axs2[0, 0].title.set_text('Bodeplot Gain')
axs2[0, 0].set_ylim([-5, 0.1])
axs2[0, 0].grid()

axs2[0, 1].plot(freq, y_gain, color='black', label='No tolerance')
axs2[0, 1].plot(freq, y_gain_max, color='black', linestyle='--', label='Max Gain')
axs2[0, 1].plot(freq, y_gain_min, color='black', linestyle='-.', label='Min Gain')
axs2[0, 1].plot(df2['#Frequency (Hz)'], df2['Channel 2 Magnitude (dB)'], color='red', label='Measured')
axs2[0, 1].title.set_text('Bodeplot Gain (Zoomed)')
axs2[0, 1].legend(loc=1)
axs2[0, 1].set_xlim([2.5e3, 7.5e3])
axs2[0, 1].set_ylim([-5, 0.1])
axs2[0, 1].grid()




axs2[1, 0].semilogx(freq, y_ph, color='black')
axs2[1, 0].semilogx(freq, y_ph_max, color='black', linestyle='--')
axs2[1, 0].semilogx(freq, y_ph_min, color='black', linestyle='-.')
axs2[1, 0].semilogx(df2['#Frequency (Hz)'], df2['Channel 2 Phase (deg)'], color='red')
axs2[1, 0].title.set_text('Bodeplot Phase-Shift')
axs2[1, 0].set_xlabel('Frequency [Hz]')
axs2[1, 0].set_ylabel('Phase [Degree]')
axs2[1, 0].grid()


axs2[1, 1].plot(freq, y_ph, color='black', label='No tolerance')
axs2[1, 1].plot(freq, y_ph_max, color='black', linestyle='--', label='Max Phase-Shift')
axs2[1, 1].plot(freq, y_ph_min, color='black', linestyle='-.', label='Min Phase-Shift')
axs2[1, 1].plot(df2['#Frequency (Hz)'], df2['Channel 2 Phase (deg)'], color='red', label='Measured')
axs2[1, 1].title.set_text('Bodeplot Phase-Shift (Zoomed)')
axs2[1, 1].set_xlabel('Frequency [Hz]')
axs2[1, 1].legend(loc=1)
axs2[1, 1].set_xlim([2.5e3, 7.5e3])
axs2[1, 1].set_ylim([-105, -75])
axs2[1, 1].grid()
fig2.tight_layout()




y = np.zeros(N, dtype=complex)
y_ph = np.zeros(N, dtype=float)
y_gain = np.zeros(N, dtype=float)
y = thirdOrder(K, R, C, R1, C1, R2, C2, freq)
for i in range(0, len(y)):
    y_ph[i] = cmath.phase(y[i]) * (180/np.pi)
y_gain = 20*np.log10((np.abs(y)))


y_max = np.zeros(N, dtype=complex)
y_ph_max = np.zeros(N, dtype=float)
y_gain_max = np.zeros(N, dtype=float)
y_max = thirdOrder(K, 
                   R*(1 - (Res_tolerance*0.5)),
                      C*(1 - (Cap_tolerance*0.5)),
                         R1*(1 + (Res_tolerance*0.5)), 
                             C1*(1 + (Cap_tolerance*0.5)), 
                                 R2*(1 - (Res_tolerance*0.5)), 
                                     C2*(1 - (Cap_tolerance*0.5)), freq)
for i in range(0, len(y_max)):
    y_ph_max[i] = cmath.phase(y_max[i]) * (180/np.pi)
y_gain_max = 20*np.log10((np.abs(y_max)))  


y_min = np.zeros(N, dtype=complex)
y_ph_min = np.zeros(N, dtype=float)
y_gain_min = np.zeros(N, dtype=float)
y_min = thirdOrder(K, 
                   R*(1 + (Res_tolerance*0.5)),
                      C*(1 + (Cap_tolerance*0.5)),
                         R1*(1 - (Res_tolerance*0.5)), 
                             C1*(1 - (Cap_tolerance*0.5)), 
                                 R2*(1 + (Res_tolerance*0.5)), 
                                     C2*(1 + (Cap_tolerance*0.5)), freq)
for i in range(0, len(y_max)):
    y_ph_min[i] = cmath.phase(y_min[i])  * (180/np.pi)
y_gain_min = 20*np.log10((np.abs(y_min))) 

fig3, axs3 = plt.subplots(2, 2, figsize=(m, M))

axs3[0, 0].semilogx(freq, y_gain, color='black')
axs3[0, 0].semilogx(freq, y_gain_max, color='black', linestyle='--')
axs3[0, 0].semilogx(freq, y_gain_min, color='black', linestyle='-.')
axs3[0, 0].semilogx(df3['#Frequency (Hz)'], df3['Channel 2 Magnitude (dB)'], color='red')
axs3[0, 0].set_ylim([-5, 0.1])
axs3[0, 0].set_ylabel('Gain [dB]')
axs3[0, 0].title.set_text('Bodeplot Gain')
axs3[0, 0].grid()

axs3[0, 1].plot(freq, y_gain, color='black', label='No tolerance')
axs3[0, 1].plot(freq, y_gain_max, color='black', linestyle='--', label='Max Gain')
axs3[0, 1].plot(freq, y_gain_min, color='black', linestyle='-.', label='Min Gain')
axs3[0, 1].plot(df3['#Frequency (Hz)'], df3['Channel 2 Magnitude (dB)'], color='red', label='Measured')
axs3[0, 1].title.set_text('Bodeplot Gain (Zoomed)')
axs3[0, 1].legend(loc=1)
axs3[0, 1].set_xlim([2.5e3, 7.5e3])
axs3[0, 1].set_ylim([-10, 0.1])
axs3[0, 1].grid()



axs3[1, 0].semilogx(freq, y_ph, color='black')
axs3[1, 0].semilogx(freq, y_ph_max, color='black', linestyle='--')
axs3[1, 0].semilogx(freq, y_ph_min, color='black', linestyle='-.')
axs3[1, 0].semilogx(df3['#Frequency (Hz)'], df3['Channel 2 Phase (deg)'], color='red')
axs3[1, 0].title.set_text('Bodeplot Phase-Shift')
axs3[1, 0].set_xlabel('Frequency [Hz]')
axs3[1, 0].set_ylabel('Phase [Degree]')
axs3[1, 0].grid()


axs3[1, 1].plot(freq, y_ph, color='black', label='No tolerance')
axs3[1, 1].plot(freq, y_ph_max, color='black', linestyle='--', label='Max Phase-Shift')
axs3[1, 1].plot(freq, y_ph_min, color='black', linestyle='-.', label='Min Phase-Shift')
axs3[1, 1].plot(df3['#Frequency (Hz)'], df3['Channel 2 Phase (deg)'], color='red', label='Measured')
axs3[1, 1].title.set_text('Bodeplot Phase-Shift (Zoomed)')
axs3[1, 1].set_xlabel('Frequency [Hz]')
axs3[1, 1].legend(loc=1)
axs3[1, 1].set_xlim([2.5e3, 7.5e3])
axs3[1, 1].set_ylim([-105-45, -75-45])
axs3[1, 1].grid()
fig3.tight_layout()

fig1.savefig('FirstOrderFilter.png')
fig2.savefig('SecondOrderFilter.png')
fig3.savefig('ThirdOrderFilter.png')

plt.show()




