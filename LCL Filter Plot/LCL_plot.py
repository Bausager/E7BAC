import numpy as np

from scipy import signal

import matplotlib.pyplot as plt

Li_tol = 0.3
Lg_tol = 0.15
C_tol = 0.1
R_tol = 0.05

Cf = 30 *1e-6
Rf = 220*1e-3
Li = 20 * 1e-3
Lg = 47 *1e-6

m = 16
M = 9


Cf_min = Cf * (1 - C_tol) 
Cf_max = Cf * (1 + C_tol)

Rf_min = Rf * (1 - R_tol) 
Rf_max = Rf * (1 + R_tol) 

Li_min = Li * (1 - Li_tol) 
Li_max = Li * (1 + Li_tol) 

Lg_min = Lg * (1 - Lg_tol) 
Lg_max = Lg * (1 + Lg_tol) 


fig1, axs1 = plt.subplots(2, 1, figsize=(m, M))

sys = signal.TransferFunction([Cf*Rf, 1], [Li*Lg*Cf, (Li + Lg)*Cf*Rf, (Li + Lg), 0])

sys_min = signal.TransferFunction([Cf_max*Rf_max, 1], [Li_max*Lg_max*Cf_max, (Li_max + Lg_max)*Cf_max*Rf_max, (Li_max + Lg_max), 0])

sys_max = signal.TransferFunction([Cf_min*Rf_min, 1], [Li_min*Lg_min*Cf_min, (Li_min + Lg_min)*Cf_min*Rf_min, (Li_min + Lg_min), 0])

w, mag, phase = signal.bode(sys)
w_min, mag_min, phase_min = signal.bode(sys_min)
w_max, mag_max, phase_max = signal.bode(sys_max)




axs1[0].semilogx(w/(2*np.pi), mag, color='r', label='No tolerance')    # Bode magnitude plot
axs1[0].semilogx(w_min/(2*np.pi), mag_min, color='black', linestyle='-.', label='Min Gain')    # Bode magnitude plot
axs1[0].semilogx(w_max/(2*np.pi), mag_max, color='black', linestyle='--', label='Max Gain')    # Bode magnitude plot
#axs1[0].set_xlabel('Frequency [Hz]')
axs1[0].set_ylabel('Gain [dB]')
axs1[0].set_xlim(10, 15e3)
axs1[0].set_ylim(-100, 20)
axs1[0].grid()

axs1[1].semilogx(w/(2*np.pi), phase, color='r', label='No tolerance')    # Bode magnitude plot
axs1[1].semilogx(w_min/(2*np.pi), phase_min, color='black', linestyle='-.', label='Min Gain')    # Bode magnitude plot
axs1[1].semilogx(w_max/(2*np.pi), phase_max, color='black', linestyle='--', label='Max Gain')    # Bode magnitude plot
axs1[1].set_xlabel('Frequency [Hz]')
axs1[1].set_ylabel('Phase [Degree]')
axs1[1].set_xlim(10, 15e3)
axs1[1].set_ylim(-251, -79)
axs1[1].grid()
axs1[1].legend()


