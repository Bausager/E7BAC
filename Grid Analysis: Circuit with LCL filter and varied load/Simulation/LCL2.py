# -*- coding: utf-8 -*-
"""
Created on Sat Mar 19 08:37:32 2022

@author: Bausa
"""

import matplotlib.pyplot as plt
import numpy as np
from pandas import read_csv

from tqdm import tqdm
import time
import warnings


plt.close('all')

df = read_csv('LCLFilters.txt', sep='\t')

E_ab_0mH_15V = df['V(15V_0mH_Ea,15V_0mH_Eb)']
E_ac_0mH_15V = df['V(15V_0mH_Ea,15V_0mH_Ec)']
E_bc_0mH_15V = df['V(15V_0mH_Eb,15V_0mH_Ec)']
E_a_0mH_15V = (E_ab_0mH_15V - (-E_ac_0mH_15V))/3
E_b_0mH_15V = (E_bc_0mH_15V - E_ab_0mH_15V)/3
E_c_0mH_15V = ((-E_ac_0mH_15V) - E_bc_0mH_15V)/3
Vrms_0mH_15V = np.sqrt(np.mean(E_a_0mH_15V**2))
Ia_0mH_15V = df['I(15v_0mh_l1)']
Ib_0mH_15V = df['I(15v_0mh_l2)']
Ic_0mH_15V = df['I(15v_0mh_l3)']
Q_0mH_15V = -((1/np.sqrt(3)) * ((Ia_0mH_15V *(E_c_0mH_15V - E_b_0mH_15V)) + 
                        (Ib_0mH_15V *(E_a_0mH_15V - E_c_0mH_15V)) + 
                        (Ic_0mH_15V * (E_b_0mH_15V - E_a_0mH_15V))))
P_0mH_15V = -((E_a_0mH_15V * Ia_0mH_15V) + 
      (E_b_0mH_15V * Ib_0mH_15V) + 
      (E_c_0mH_15V * Ic_0mH_15V))

E_ab_0mH_20V = df['V(20V_0mH_Ea,20V_0mH_Eb)']
E_ac_0mH_20V = df['V(20V_0mH_Ea,20V_0mH_Ec)']
E_bc_0mH_20V = df['V(20V_0mH_Eb,20V_0mH_Ec)']
E_a_0mH_20V = (E_ab_0mH_20V - (-E_ac_0mH_20V))/3
E_b_0mH_20V = (E_bc_0mH_20V - E_ab_0mH_20V)/3
E_c_0mH_20V = ((-E_ac_0mH_20V) - E_bc_0mH_20V)/3
Vrms_0mH_20V = np.sqrt(np.mean(E_a_0mH_20V**2))
Ia_0mH_20V = df['I(20v_0mh_l1)']
Ib_0mH_20V = df['I(20v_0mh_l2)']
Ic_0mH_20V = df['I(20v_0mh_l3)']
Q_0mH_20V = -((1/np.sqrt(3)) * ((Ia_0mH_20V *(E_c_0mH_20V - E_b_0mH_20V)) + 
                        (Ib_0mH_20V *(E_a_0mH_20V - E_c_0mH_20V)) + 
                        (Ic_0mH_20V * (E_b_0mH_20V - E_a_0mH_20V))))
P_0mH_20V = -((E_a_0mH_20V * Ia_0mH_20V) + 
      (E_b_0mH_20V * Ib_0mH_20V) + 
      (E_c_0mH_20V * Ic_0mH_20V))

E_ab_0mH_25V = df['V(25V_0mH_Ea,25V_0mH_Eb)']
E_ac_0mH_25V = df['V(25V_0mH_Ea,25V_0mH_Ec)']
E_bc_0mH_25V = df['V(25V_0mH_Eb,25V_0mH_Ec)']
E_a_0mH_25V = (E_ab_0mH_25V - (-E_ac_0mH_25V))/3
E_b_0mH_25V = (E_bc_0mH_25V - E_ab_0mH_25V)/3
E_c_0mH_25V = ((-E_ac_0mH_25V) - E_bc_0mH_25V)/3
Vrms_0mH_25V = np.sqrt(np.mean(E_a_0mH_25V**2))
Ia_0mH_25V = df['I(25v_0mh_l1)']
Ib_0mH_25V = df['I(25v_0mh_l2)']
Ic_0mH_25V = df['I(25v_0mh_l3)']
Q_0mH_25V = -((1/np.sqrt(3)) * ((Ia_0mH_25V *(E_c_0mH_25V - E_b_0mH_25V)) + 
                        (Ib_0mH_25V *(E_a_0mH_25V - E_c_0mH_25V)) + 
                        (Ic_0mH_25V * (E_b_0mH_25V - E_a_0mH_25V))))
P_0mH_25V = -((E_a_0mH_25V * Ia_0mH_25V) + 
      (E_b_0mH_25V * Ib_0mH_25V) + 
      (E_c_0mH_25V * Ic_0mH_25V))

E_ab_0mH_30V = df['V(30V_Ea,30V_Eb)']
E_ac_0mH_30V = df['V(30V_Ea,30V_Ec)']
E_bc_0mH_30V = df['V(30V_Eb,30V_Ec)']
E_a_0mH_30V = (E_ab_0mH_30V - (-E_ac_0mH_30V))/3
E_b_0mH_30V = (E_bc_0mH_30V - E_ab_0mH_30V)/3
E_c_0mH_30V = ((-E_ac_0mH_30V) - E_bc_0mH_30V)/3
Vrms_0mH_30V = np.sqrt(np.mean(E_a_0mH_30V**2))
Ia_0mH_30V = df['I(30v_0mh_l1)']
Ib_0mH_30V = df['I(30v_0mh_l2)']
Ic_0mH_30V = df['I(30v_0mh_l3)']
Q_0mH_30V = -((1/np.sqrt(3)) * ((Ia_0mH_30V *(E_c_0mH_30V - E_b_0mH_30V)) + 
                        (Ib_0mH_30V *(E_a_0mH_30V - E_c_0mH_30V)) + 
                        (Ic_0mH_30V * (E_b_0mH_30V - E_a_0mH_30V))))
P_0mH_30V = -((E_a_0mH_30V * Ia_0mH_30V) + 
      (E_b_0mH_30V * Ib_0mH_30V) + 
      (E_c_0mH_30V * Ic_0mH_30V))


fig, axs = plt.subplots(3,4)

axs[0,0].set_ylabel('Voltage [V]')
axs[0,0].set_title(r'$V_{RMS}$: ' + f'{round(Vrms_0mH_15V,3)}' + '; 0mH & 15 Ohm')
axs[0,0].plot(df['time'], E_ab_0mH_15V, color='r', alpha=0.2, label=r'$E_{ab}$')
axs[0,0].plot(df['time'],E_ac_0mH_15V, color='b', alpha=0.2, label=r'$E_{ac}$')
axs[0,0].plot(df['time'],E_bc_0mH_15V, color='g', alpha=0.2, label=r'$E_{bc}$')
axs[0,0].plot(df['time'],E_a_0mH_15V, color='r', label=r'$E_{a}$')
axs[0,0].plot(df['time'],E_b_0mH_15V, color='b', label=r'$E_{b}$')
axs[0,0].plot(df['time'],E_c_0mH_15V, color='g', label=r'$E_{c}$')
axs[0,0].grid()
axs[1,0].set_ylabel('Ampere [A]')
axs[1,0].plot(df['time'],Ia_0mH_15V, color='r', label=r'$I_{a}$')
axs[1,0].plot(df['time'],Ib_0mH_15V, color='b', label=r'$I_{a}$')
axs[1,0].plot(df['time'],Ic_0mH_15V, color='g', label=r'$I_{a}$')
axs[1,0].grid()
axs[2,0].set_ylabel('Power [W/VAr]')
axs[2,0].set_xlabel('Time [s]')
axs[2,0].plot(df['time'],Q_0mH_15V, color='r', label=r'$Q$')
axs[2,0].plot(df['time'],P_0mH_15V, color='b', label=r'$P$')
axs[2,0].legend(loc=5)
axs[2,0].grid()

axs[0,1].set_title(r'$V_{RMS}$: ' + f'{round(Vrms_0mH_20V,3)}' + '; 0mH & 15 Ohm')
axs[0,1].plot(df['time'], E_ab_0mH_20V, color='r', alpha=0.2, label=r'$E_{ab}$')
axs[0,1].plot(df['time'],E_ac_0mH_20V, color='b', alpha=0.2, label=r'$E_{ac}$')
axs[0,1].plot(df['time'],E_bc_0mH_20V, color='g', alpha=0.2, label=r'$E_{bc}$')
axs[0,1].plot(df['time'],E_a_0mH_20V, color='r', label=r'$E_{a}$')
axs[0,1].plot(df['time'],E_b_0mH_20V, color='b', label=r'$E_{b}$')
axs[0,1].plot(df['time'],E_c_0mH_20V, color='g', label=r'$E_{c}$')
axs[0,1].grid()
axs[1,1].plot(df['time'],Ia_0mH_20V, color='r', label=r'$I_{a}$')
axs[1,1].plot(df['time'],Ib_0mH_20V, color='b', label=r'$I_{a}$')
axs[1,1].plot(df['time'],Ic_0mH_20V, color='g', label=r'$I_{a}$')
axs[1,1].grid()
axs[2,1].set_xlabel('Time [s]')
axs[2,1].plot(df['time'],Q_0mH_20V, color='r', label=r'$Q$')
axs[2,1].plot(df['time'],P_0mH_20V, color='b', label=r'$P$')
axs[2,1].legend(loc=5)
axs[2,1].grid()

axs[0,2].set_title(r'$V_{RMS}$: ' + f'{round(Vrms_0mH_25V,3)}' + '; 0mH & 15 Ohm')
axs[0,2].plot(df['time'],E_ab_0mH_25V, color='r', alpha=0.2, label=r'$E_{ab}$')
axs[0,2].plot(df['time'],E_ac_0mH_25V, color='b', alpha=0.2, label=r'$E_{ac}$')
axs[0,2].plot(df['time'],E_bc_0mH_25V, color='g', alpha=0.2, label=r'$E_{bc}$')
axs[0,2].plot(df['time'],E_a_0mH_25V, color='r', label=r'$E_{a}$')
axs[0,2].plot(df['time'],E_b_0mH_25V, color='b', label=r'$E_{b}$')
axs[0,2].plot(df['time'],E_c_0mH_25V, color='g', label=r'$E_{c}$')
axs[0,2].grid()
axs[1,2].plot(df['time'],Ia_0mH_25V, color='r', label=r'$I_{a}$')
axs[1,2].plot(df['time'],Ib_0mH_25V, color='b', label=r'$I_{a}$')
axs[1,2].plot(df['time'],Ic_0mH_25V, color='g', label=r'$I_{a}$')
axs[1,2].grid()
axs[2,2].set_xlabel('Time [s]')
axs[2,2].plot(df['time'],Q_0mH_25V, color='r', label=r'$Q$')
axs[2,2].plot(df['time'],P_0mH_25V, color='b', label=r'$P$')
axs[2,2].legend(loc=5)
axs[2,2].grid()

axs[0,3].set_title(r'$V_{RMS}$: ' + f'{round(Vrms_0mH_30V,3)}' + '; 0mH & 15 Ohm')
axs[0,3].plot(df['time'],E_ab_0mH_30V, color='r', alpha=0.2, label=r'$E_{ab}$')
axs[0,3].plot(df['time'],E_ac_0mH_30V, color='b', alpha=0.2, label=r'$E_{ac}$')
axs[0,3].plot(df['time'],E_bc_0mH_30V, color='g', alpha=0.2, label=r'$E_{bc}$')
axs[0,3].plot(df['time'],E_a_0mH_30V, color='r', label=r'$E_{a}$')
axs[0,3].plot(df['time'],E_b_0mH_30V, color='b', label=r'$E_{b}$')
axs[0,3].plot(df['time'],E_c_0mH_30V, color='g', label=r'$E_{c}$')
axs[0,3].grid()
axs[1,3].plot(df['time'],Ia_0mH_30V, color='r', label=r'$I_{a}$')
axs[1,3].plot(df['time'],Ib_0mH_30V, color='b', label=r'$I_{a}$')
axs[1,3].plot(df['time'],Ic_0mH_30V, color='g', label=r'$I_{a}$')
axs[1,3].grid()
axs[2,3].set_xlabel('Time [s]')
axs[2,3].plot(df['time'],Q_0mH_30V, color='r', label=r'$Q$')
axs[2,3].plot(df['time'],P_0mH_30V, color='b', label=r'$P$')
axs[2,3].legend(loc=5)
axs[2,3].grid()

E_ab_433mH_15V = df['V(15V_4.33mH_Ea,15V_4.33mH_Eb)']
E_ac_433mH_15V = df['V(15V_4.33mH_Ea,15V_4.33mH_Ec)']
E_bc_433mH_15V = df['V(15V_0mH_Eb,15V_0mH_Ec)']
E_a_433mH_15V = (E_ab_433mH_15V - (-E_ac_433mH_15V))/3
E_b_433mH_15V = (E_bc_433mH_15V - E_ab_433mH_15V)/3
E_c_433mH_15V = ((-E_ac_433mH_15V) - E_bc_433mH_15V)/3
Vrms_433mH_15V = np.sqrt(np.mean(E_a_433mH_15V**2))
Ia_433mH_15V = df['I(15v_4.33mh_l1)']
Ib_433mH_15V = df['I(15v_4.33mh_l2)']
Ic_433mH_15V = df['I(15v_4.33mh_l3)']
Q_433mH_15V = -((1/np.sqrt(3)) * ((Ia_433mH_15V *(E_c_433mH_15V - E_b_433mH_15V)) + 
                        (Ib_433mH_15V *(E_a_433mH_15V - E_c_433mH_15V)) + 
                        (Ic_433mH_15V * (E_b_433mH_15V - E_a_433mH_15V))))
P_433mH_15V = -((E_a_433mH_15V * Ia_433mH_15V) + 
      (E_b_433mH_15V * Ib_433mH_15V) + 
      (E_c_433mH_15V * Ic_433mH_15V))


E_ab_433mH_20V = df['V(20V_4.33mH_Ea,20V_4.33mH_Eb)']
E_ac_433mH_20V = df['V(20V_4.33mH_Ea,20V_4.33mH_Ec)']
E_bc_433mH_20V = df['V(20V_4.33mH_Eb,20V_4.33mH_Ec)']
E_a_433mH_20V = (E_ab_433mH_20V - (-E_ac_433mH_20V))/3
E_b_433mH_20V = (E_bc_433mH_20V - E_ab_433mH_20V)/3
E_c_433mH_20V = ((-E_ac_433mH_20V) - E_bc_433mH_20V)/3
Vrms_433mH_20V = np.sqrt(np.mean(E_a_433mH_20V**2))
Ia_433mH_20V = df['I(20v_4.33mh_l1)']
Ib_433mH_20V = df['I(20v_4.33mh_l2)']
Ic_433mH_20V = df['I(20v_4.33mh_l3)']
Q_433mH_20V = -((1/np.sqrt(3)) * ((Ia_433mH_20V *(E_c_433mH_20V - E_b_433mH_20V)) + 
                        (Ib_433mH_20V *(E_a_433mH_20V - E_c_433mH_20V)) + 
                        (Ic_433mH_20V * (E_b_433mH_20V - E_a_433mH_20V))))
P_433mH_20V = -((E_a_433mH_20V * Ia_433mH_20V) + 
      (E_b_433mH_20V * Ib_433mH_20V) + 
      (E_c_433mH_20V * Ic_433mH_20V))

E_ab_433mH_25V = df['V(25V_4.33mH_Ea,25V_4.33mH_Eb)']
E_ac_433mH_25V = df['V(25V_4.33mH_Ea,25V_4.33mH_Ec)']
E_bc_433mH_25V = df['V(25V_4.33mH_Eb,25V_4.33mH_Ec)']
E_a_433mH_25V = (E_ab_433mH_25V - (-E_ac_433mH_25V))/3
E_b_433mH_25V = (E_bc_433mH_25V - E_ab_433mH_25V)/3
E_c_433mH_25V = ((-E_ac_433mH_25V) - E_bc_433mH_25V)/3
Vrms_433mH_25V = np.sqrt(np.mean(E_a_433mH_25V**2))
Ia_433mH_25V = df['I(25v_4.33mh_l1)']
Ib_433mH_25V = df['I(25v_4.33mh_l2)']
Ic_433mH_25V = df['I(25v_4.33mh_l3)']
Q_433mH_25V = -((1/np.sqrt(3)) * ((Ia_433mH_25V *(E_c_433mH_25V - E_b_433mH_25V)) + 
                        (Ib_433mH_25V *(E_a_433mH_25V - E_c_433mH_25V)) + 
                        (Ic_433mH_25V * (E_b_433mH_25V - E_a_433mH_25V))))
P_433mH_25V = -((E_a_433mH_25V * Ia_433mH_25V) + 
      (E_b_433mH_25V * Ib_433mH_25V) + 
      (E_c_433mH_25V * Ic_433mH_25V))

E_ab_433mH_30V = df['V(30V_4.33mH_Ea,30V_4.33mH_Eb)']
E_ac_433mH_30V = df['V(30V_4.33mH_Ea,30V_4.33mH_Ec)']
E_bc_433mH_30V = df['V(30V_4.33mH_Eb,30V_4.33mH_Ec)']
E_a_433mH_30V = (E_ab_433mH_30V - (-E_ac_433mH_30V))/3
E_b_433mH_30V = (E_bc_433mH_30V - E_ab_433mH_30V)/3
E_c_433mH_30V = ((-E_ac_433mH_30V) - E_bc_433mH_30V)/3
Vrms_433mH_30V = np.sqrt(np.mean(E_a_433mH_30V**2))
Ia_433mH_30V = df['I(30v_4.33mh_l1)']
Ib_433mH_30V = df['I(30v_4.33mh_l2)']
Ic_433mH_30V = df['I(30v_4.33mh_l3)']
Q_433mH_30V = -((1/np.sqrt(3)) * ((Ia_433mH_30V *(E_c_433mH_30V - E_b_433mH_30V)) + 
                        (Ib_433mH_30V *(E_a_433mH_30V - E_c_433mH_30V)) + 
                        (Ic_433mH_30V * (E_b_433mH_30V - E_a_433mH_30V))))
P_433mH_30V = -((E_a_433mH_30V * Ia_433mH_30V) + 
      (E_b_433mH_30V * Ib_433mH_30V) + 
      (E_c_433mH_30V * Ic_433mH_30V))

fig1, axs1 = plt.subplots(3,4)

axs1[0,0].set_title(r'$V_{RMS}$: ' + f'{round(Vrms_433mH_15V,3)}' + '; 4.33mH & 15 Ohm')
axs1[0,0].set_ylabel('Voltage [V]')
axs1[0,0].plot(df['time'], E_ab_433mH_15V, color='r', alpha=0.2, label=r'$E_{ab}$')
axs1[0,0].plot(df['time'],E_ac_433mH_15V, color='b', alpha=0.2, label=r'$E_{ac}$')
axs1[0,0].plot(df['time'],E_bc_433mH_15V, color='g', alpha=0.2, label=r'$E_{bc}$')
axs1[0,0].plot(df['time'],E_a_433mH_15V, color='r', label=r'$E_{a}$')
axs1[0,0].plot(df['time'],E_b_433mH_15V, color='b', label=r'$E_{b}$')
axs1[0,0].plot(df['time'],E_c_433mH_15V, color='g', label=r'$E_{c}$')
axs1[0,0].grid()
axs1[1,0].set_ylabel('Ampere [A]')
axs1[1,0].plot(df['time'],Ia_433mH_15V, color='r', label=r'$I_{a}$')
axs1[1,0].plot(df['time'],Ib_433mH_15V, color='b', label=r'$I_{a}$')
axs1[1,0].plot(df['time'],Ic_433mH_15V, color='g', label=r'$I_{a}$')
axs1[1,0].grid()
axs1[2,0].set_ylabel('Power [W/VAr]')
axs1[2,0].set_xlabel('Time [s]')
axs1[2,0].plot(df['time'],Q_433mH_15V, color='r', label=r'$Q$')
axs1[2,0].plot(df['time'],P_433mH_15V, color='b', label=r'$P$')
axs1[2,0].legend(loc=5)
axs1[2,0].grid()


axs1[0,1].set_title(r'$V_{RMS}$: ' + f'{round(Vrms_433mH_20V,3)}' + '; 4.33mH & 15 Ohm')
axs1[0,1].plot(df['time'], E_ab_433mH_20V, color='r', alpha=0.2, label=r'$E_{ab}$')
axs1[0,1].plot(df['time'],E_ac_433mH_20V, color='b', alpha=0.2, label=r'$E_{ac}$')
axs1[0,1].plot(df['time'],E_bc_433mH_20V, color='g', alpha=0.2, label=r'$E_{bc}$')
axs1[0,1].plot(df['time'],E_a_433mH_20V, color='r', label=r'$E_{a}$')
axs1[0,1].plot(df['time'],E_b_433mH_20V, color='b', label=r'$E_{b}$')
axs1[0,1].plot(df['time'],E_c_433mH_20V, color='g', label=r'$E_{c}$')
axs1[0,1].grid()
axs1[1,1].plot(df['time'],Ia_433mH_20V, color='r', label=r'$I_{a}$')
axs1[1,1].plot(df['time'],Ib_433mH_20V, color='b', label=r'$I_{a}$')
axs1[1,1].plot(df['time'],Ic_433mH_20V, color='g', label=r'$I_{a}$')
axs1[1,1].grid()
axs1[2,1].set_xlabel('Time [s]')
axs1[2,1].plot(df['time'],Q_433mH_20V, color='r', label=r'$Q$')
axs1[2,1].plot(df['time'],P_433mH_20V, color='b', label=r'$P$')
axs1[2,1].legend(loc=5)
axs1[2,1].grid()

axs1[0,2].set_title(r'$V_{RMS}$: ' + f'{round(Vrms_433mH_25V,3)}' + '; 433mH & 15 Ohm')
axs1[0,2].plot(df['time'],E_ab_433mH_25V, color='r', alpha=0.2, label=r'$E_{ab}$')
axs1[0,2].plot(df['time'],E_ac_433mH_25V, color='b', alpha=0.2, label=r'$E_{ac}$')
axs1[0,2].plot(df['time'],E_bc_433mH_25V, color='g', alpha=0.2, label=r'$E_{bc}$')
axs1[0,2].plot(df['time'],E_a_433mH_25V, color='r', label=r'$E_{a}$')
axs1[0,2].plot(df['time'],E_b_433mH_25V, color='b', label=r'$E_{b}$')
axs1[0,2].plot(df['time'],E_c_433mH_25V, color='g', label=r'$E_{c}$')
axs1[0,2].grid()
axs1[1,2].plot(df['time'],Ia_433mH_25V, color='r', label=r'$I_{a}$')
axs1[1,2].plot(df['time'],Ib_433mH_25V, color='b', label=r'$I_{a}$')
axs1[1,2].plot(df['time'],Ic_433mH_25V, color='g', label=r'$I_{a}$')
axs1[1,2].grid()
axs1[2,2].set_xlabel('Time [s]')
axs1[2,2].plot(df['time'],Q_433mH_25V, color='r', label=r'$Q$')
axs1[2,2].plot(df['time'],P_433mH_25V, color='b', label=r'$P$')
axs1[2,2].legend(loc=5)
axs1[2,2].grid()

axs1[0,3].set_title(r'$V_{RMS}$: ' + f'{round(Vrms_0mH_30V,3)}' + '; 433mH & 15 Ohm')
axs1[0,3].plot(df['time'],E_ab_433mH_30V, color='r', alpha=0.2, label=r'$E_{ab}$')
axs1[0,3].plot(df['time'],E_ac_433mH_30V, color='b', alpha=0.2, label=r'$E_{ac}$')
axs1[0,3].plot(df['time'],E_bc_433mH_30V, color='g', alpha=0.2, label=r'$E_{bc}$')
axs1[0,3].plot(df['time'],E_a_433mH_30V, color='r', label=r'$E_{a}$')
axs1[0,3].plot(df['time'],E_b_433mH_30V, color='b', label=r'$E_{b}$')
axs1[0,3].plot(df['time'],E_c_433mH_30V, color='g', label=r'$E_{c}$')
axs1[0,3].grid()
axs1[1,3].plot(df['time'],Ia_433mH_30V, color='r', label=r'$I_{a}$')
axs1[1,3].plot(df['time'],Ib_433mH_30V, color='b', label=r'$I_{a}$')
axs1[1,3].plot(df['time'],Ic_433mH_30V, color='g', label=r'$I_{a}$')
axs1[1,3].grid()
axs1[2,3].set_xlabel('Time [s]')
axs1[2,3].plot(df['time'],Q_433mH_30V, color='r', label=r'$Q$')
axs1[2,3].plot(df['time'],P_433mH_30V, color='b', label=r'$P$')
axs1[2,3].legend(loc=5)
axs1[2,3].grid()


E_ab_866mH_15V = df['V(15V_8.66mH_Ea,15V_8.66mH_Eb)']
E_ac_866mH_15V = df['V(15V_8.66mH_Ea,15V_8.66mH_Ec)']
E_bc_866mH_15V = df['V(15V_8.66mH_Eb,15V_8.66mH_Ec)']
E_a_866mH_15V = (E_ab_866mH_15V - (-E_ac_866mH_15V))/3
E_b_866mH_15V = (E_bc_866mH_15V - E_ab_866mH_15V)/3
E_c_866mH_15V = ((-E_ac_866mH_15V) - E_bc_866mH_15V)/3
Vrms_866mH_15V = np.sqrt(np.mean(E_a_866mH_15V**2))
Ia_866mH_15V = df['I(15v_8.66mh_l1)']
Ib_866mH_15V = df['I(15v_8.66mh_l2)']
Ic_866mH_15V = df['I(15v_8.66mh_l3)']
Q_866mH_15V = -((1/np.sqrt(3)) * ((Ia_866mH_15V *(E_c_866mH_15V - E_b_866mH_15V)) + 
                        (Ib_866mH_15V *(E_a_866mH_15V - E_c_866mH_15V)) + 
                        (Ic_866mH_15V * (E_b_866mH_15V - E_a_866mH_15V))))
P_866mH_15V = -((E_a_866mH_15V * Ia_866mH_15V) + 
      (E_b_866mH_15V * Ib_866mH_15V) + 
      (E_c_866mH_15V * Ic_866mH_15V))


E_ab_866mH_20V = df['V(20V_8.66mH_Ea,20V_8.66mH_Eb)']
E_ac_866mH_20V = df['V(20V_8.66mH_Ea,20V_8.66mH_Ec)']
E_bc_866mH_20V = df['V(20V_8.66mH_Eb,20V_8.66mH_Ec)']
E_a_866mH_20V = (E_ab_866mH_20V - (-E_ac_866mH_20V))/3
E_b_866mH_20V = (E_bc_866mH_20V - E_ab_866mH_20V)/3
E_c_866mH_20V = ((-E_ac_866mH_20V) - E_bc_866mH_20V)/3
Vrms_866mH_20V = np.sqrt(np.mean(E_a_866mH_20V**2))
Ia_866mH_20V = df['I(20v_8.66mh_l1)']
Ib_866mH_20V = df['I(20v_8.66mh_l2)']
Ic_866mH_20V = df['I(20v_8.66mh_l3)']
Q_866mH_20V = -((1/np.sqrt(3)) * ((Ia_866mH_20V *(E_c_866mH_20V - E_b_866mH_20V)) + 
                        (Ib_866mH_20V *(E_a_866mH_20V - E_c_866mH_20V)) + 
                        (Ic_866mH_20V * (E_b_866mH_20V - E_a_866mH_20V))))
P_866mH_20V = -((E_a_866mH_20V * Ia_866mH_20V) + 
      (E_b_866mH_20V * Ib_866mH_20V) + 
      (E_c_866mH_20V * Ic_866mH_20V))

E_ab_866mH_25V = df['V(25V_8.66mH_Ea,25V_8.66mH_Eb)']
E_ac_866mH_25V = df['V(25V_8.66mH_Ea,25V_8.66mH_Ec)']
E_bc_866mH_25V = df['V(25V_8.66mH_Eb,25V_8.66mH_Ec)']
E_a_866mH_25V = (E_ab_866mH_25V - (-E_ac_866mH_25V))/3
E_b_866mH_25V = (E_bc_866mH_25V - E_ab_866mH_25V)/3
E_c_866mH_25V = ((-E_ac_866mH_25V) - E_bc_866mH_25V)/3
Vrms_866mH_25V = np.sqrt(np.mean(E_a_866mH_25V**2))
Ia_866mH_25V = df['I(25v_8.66mh_l1)']
Ib_866mH_25V = df['I(25v_8.66mh_l2)']
Ic_866mH_25V = df['I(25v_8.66mh_l3)']
Q_866mH_25V = -((1/np.sqrt(3)) * ((Ia_866mH_25V *(E_c_866mH_25V - E_b_866mH_25V)) + 
                        (Ib_866mH_25V *(E_a_866mH_25V - E_c_866mH_25V)) + 
                        (Ic_866mH_25V * (E_b_866mH_25V - E_a_866mH_25V))))
P_866mH_25V = -((E_a_866mH_25V * Ia_866mH_25V) + 
      (E_b_866mH_25V * Ib_866mH_25V) + 
      (E_c_866mH_25V * Ic_866mH_25V))

E_ab_866mH_30V = df['V(30V_8.66mH_Ea,30V_8.66mH_Eb)']
E_ac_866mH_30V = df['V(30V_8.66mH_Ea,30V_8.66mH_Ec)']
E_bc_866mH_30V = df['V(30V_8.66mH_Eb,30V_8.66mH_Ec)']
E_a_866mH_30V = (E_ab_866mH_30V - (-E_ac_866mH_30V))/3
E_b_866mH_30V = (E_bc_866mH_30V - E_ab_866mH_30V)/3
E_c_866mH_30V = ((-E_ac_866mH_30V) - E_bc_866mH_30V)/3
Vrms_866mH_30V = np.sqrt(np.mean(E_a_866mH_30V**2))
Ia_866mH_30V = df['I(30v_8.66mh_l1)']
Ib_866mH_30V = df['I(30v_8.66mh_l2)']
Ic_866mH_30V = df['I(30v_8.66mh_l3)']
Q_866mH_30V = -((1/np.sqrt(3)) * ((Ia_866mH_30V *(E_c_866mH_30V - E_b_866mH_30V)) + 
                        (Ib_866mH_30V *(E_a_866mH_30V - E_c_866mH_30V)) + 
                        (Ic_866mH_30V * (E_b_866mH_30V - E_a_866mH_30V))))
P_866mH_30V = -((E_a_866mH_30V * Ia_866mH_30V) + 
      (E_b_866mH_30V * Ib_866mH_30V) + 
      (E_c_866mH_30V * Ic_866mH_30V))

fig2, axs2 = plt.subplots(3,4)

axs2[0,0].set_ylabel('Voltage [V]')
axs2[0,0].set_title(r'$V_{RMS}$: ' + f'{round(Vrms_866mH_15V,3)}' + '; 0mH & 15 Ohm')
axs2[0,0].plot(df['time'], E_ab_866mH_15V, color='r', alpha=0.2, label=r'$E_{ab}$')
axs2[0,0].plot(df['time'],E_ac_866mH_15V, color='b', alpha=0.2, label=r'$E_{ac}$')
axs2[0,0].plot(df['time'],E_bc_866mH_15V, color='g', alpha=0.2, label=r'$E_{bc}$')
axs2[0,0].plot(df['time'],E_a_866mH_15V, color='r', label=r'$E_{a}$')
axs2[0,0].plot(df['time'],E_b_866mH_15V, color='b', label=r'$E_{b}$')
axs2[0,0].plot(df['time'],E_c_866mH_15V, color='g', label=r'$E_{c}$')
axs2[0,0].grid()
axs2[1,0].set_ylabel('Ampere [A]')
axs2[1,0].plot(df['time'],Ia_866mH_15V, color='r', label=r'$I_{a}$')
axs2[1,0].plot(df['time'],Ib_866mH_15V, color='b', label=r'$I_{a}$')
axs2[1,0].plot(df['time'],Ic_866mH_15V, color='g', label=r'$I_{a}$')
axs2[1,0].grid()
axs2[2,0].set_ylabel('Power [W/VAr]')
axs2[2,0].set_xlabel('Time [s]')
axs2[2,0].plot(df['time'],Q_866mH_15V, color='r', label=r'$Q$')
axs2[2,0].plot(df['time'],P_866mH_15V, color='b', label=r'$P$')
axs2[2,0].legend(loc=5)
axs2[2,0].grid()


axs2[0,1].set_title(r'$V_{RMS}$: ' + f'{round(Vrms_866mH_20V,3)}' + '; 866mH & 15 Ohm')
axs2[0,1].plot(df['time'], E_ab_866mH_20V, color='r', alpha=0.2, label=r'$E_{ab}$')
axs2[0,1].plot(df['time'],E_ac_866mH_20V, color='b', alpha=0.2, label=r'$E_{ac}$')
axs2[0,1].plot(df['time'],E_bc_866mH_20V, color='g', alpha=0.2, label=r'$E_{bc}$')
axs2[0,1].plot(df['time'],E_a_866mH_20V, color='r', label=r'$E_{a}$')
axs2[0,1].plot(df['time'],E_b_866mH_20V, color='b', label=r'$E_{b}$')
axs2[0,1].plot(df['time'],E_c_866mH_20V, color='g', label=r'$E_{c}$')
axs2[0,1].grid()
axs2[1,1].plot(df['time'],Ia_866mH_20V, color='r', label=r'$I_{a}$')
axs2[1,1].plot(df['time'],Ib_866mH_20V, color='b', label=r'$I_{a}$')
axs2[1,1].plot(df['time'],Ic_866mH_20V, color='g', label=r'$I_{a}$')
axs2[1,1].grid()
axs2[2,1].set_xlabel('Time [s]')
axs2[2,1].plot(df['time'],Q_866mH_20V, color='r', label=r'$Q$')
axs2[2,1].plot(df['time'],P_866mH_20V, color='b', label=r'$P$')
axs2[2,1].legend(loc=5)
axs2[2,1].grid()

axs2[0,2].set_title(r'$V_{RMS}$: ' + f'{round(Vrms_866mH_25V,3)}' + '; 866mH & 15 Ohm')
axs2[0,2].plot(df['time'],E_ab_866mH_25V, color='r', alpha=0.2, label=r'$E_{ab}$')
axs2[0,2].plot(df['time'],E_ac_866mH_25V, color='b', alpha=0.2, label=r'$E_{ac}$')
axs2[0,2].plot(df['time'],E_bc_866mH_25V, color='g', alpha=0.2, label=r'$E_{bc}$')
axs2[0,2].plot(df['time'],E_a_866mH_25V, color='r', label=r'$E_{a}$')
axs2[0,2].plot(df['time'],E_b_866mH_25V, color='b', label=r'$E_{b}$')
axs2[0,2].plot(df['time'],E_c_866mH_25V, color='g', label=r'$E_{c}$')
axs2[0,2].grid()
axs2[1,2].plot(df['time'],Ia_866mH_25V, color='r', label=r'$I_{a}$')
axs2[1,2].plot(df['time'],Ib_866mH_25V, color='b', label=r'$I_{a}$')
axs2[1,2].plot(df['time'],Ic_866mH_25V, color='g', label=r'$I_{a}$')
axs2[1,2].grid()
axs2[2,2].set_xlabel('Time [s]')
axs2[2,2].plot(df['time'],Q_866mH_25V, color='r', label=r'$Q$')
axs2[2,2].plot(df['time'],P_866mH_25V, color='b', label=r'$P$')
axs2[2,2].legend(loc=5)
axs2[2,2].grid()

axs2[0,3].set_title(r'$V_{RMS}$: ' + f'{round(Vrms_866mH_30V,3)}' + '; 866mH & 15 Ohm')
axs2[0,3].plot(df['time'],E_ab_866mH_30V, color='r', alpha=0.2, label=r'$E_{ab}$')
axs2[0,3].plot(df['time'],E_ac_866mH_30V, color='b', alpha=0.2, label=r'$E_{ac}$')
axs2[0,3].plot(df['time'],E_bc_866mH_30V, color='g', alpha=0.2, label=r'$E_{bc}$')
axs2[0,3].plot(df['time'],E_a_866mH_30V, color='r', label=r'$E_{a}$')
axs2[0,3].plot(df['time'],E_b_866mH_30V, color='b', label=r'$E_{b}$')
axs2[0,3].plot(df['time'],E_c_866mH_30V, color='g', label=r'$E_{c}$')
axs2[0,3].grid()
axs2[1,3].plot(df['time'],Ia_866mH_30V, color='r', label=r'$I_{a}$')
axs2[1,3].plot(df['time'],Ib_866mH_30V, color='b', label=r'$I_{a}$')
axs2[1,3].plot(df['time'],Ic_866mH_30V, color='g', label=r'$I_{a}$')
axs2[1,3].grid()
axs2[2,3].set_xlabel('Time [s]')
axs2[2,3].plot(df['time'],Q_866mH_30V, color='r', label=r'$Q$')
axs2[2,3].plot(df['time'],P_866mH_30V, color='b', label=r'$P$')
axs2[2,3].legend(loc=5)
axs2[2,3].grid()




N = 4

V1 = np.zeros(N, dtype='float')
V2 = np.zeros(N, dtype='float')
V3 = np.zeros(N, dtype='float')
V1 = [Vrms_0mH_15V, Vrms_0mH_20V, Vrms_0mH_25V, Vrms_0mH_30V]
V2 = [Vrms_433mH_15V, Vrms_433mH_20V, Vrms_433mH_25V, Vrms_433mH_30V]
V3 = [Vrms_866mH_15V, Vrms_866mH_20V, Vrms_866mH_25V, Vrms_866mH_30V]

P1 = np.zeros(N, dtype='float')
P2 = np.zeros(N, dtype='float')
P3 = np.zeros(N, dtype='float')
P1 = [P_0mH_15V[0], P_0mH_20V[0], P_0mH_25V[0], P_0mH_30V[0]]
P2 = [P_433mH_15V[0], P_433mH_20V[0], P_433mH_25V[0], P_433mH_30V[0]]
P3 = [P_866mH_15V[0], P_866mH_20V[0], P_866mH_25V[0], P_866mH_30V[0]]

Q1 = np.zeros(N, dtype='float')
Q2 = np.zeros(N, dtype='float')
Q3 = np.zeros(N, dtype='float')
Q1 = [Q_0mH_15V[0], Q_0mH_20V[0], Q_0mH_25V[0], Q_0mH_30V[0]]
Q2 = [Q_433mH_15V[0], Q_433mH_20V[0], Q_433mH_25V[0], Q_433mH_30V[0]]
Q3 = [Q_866mH_15V[0], Q_866mH_20V[0], Q_866mH_25V[0], Q_866mH_30V[0]]


Y_meas1 = np.zeros([N, 3], dtype='float')
Y_meas2 = np.zeros([N, 3], dtype='float')
Y_meas3 = np.zeros([N, 3], dtype='float')


for i in range(0, N):
    Y_meas1[i, 0] = V1[i]
    Y_meas1[i, 1] = P1[i]
    Y_meas1[i, 2] = Q1[i]
    
    Y_meas2[i, 0] = V2[i]
    Y_meas2[i, 1] = P2[i]
    Y_meas2[i, 2] = Q2[i]
    
    Y_meas3[i, 0] = V3[i]
    Y_meas3[i, 1] = P3[i]
    Y_meas3[i, 2] = Q3[i]


plt.close('all')

n = 44
Y1 = np.zeros([n, 4], dtype='float')
Y2 = np.zeros([n, 4], dtype='float')
Y3 = np.zeros([n, 4], dtype='float')



def J(Eg, R, L, V, P, Q):
    

    div1 = ( ((R**2) + (L**2))*((P**2) + (Q**2)) )/(9 * (V**2))
    div2 = ( (R*P) + (L*Q) )/(3)
    
    y = (V**2) - (Eg**2) + div1 - (2*div2)
    return np.power(y,2)


def Stepping(N, Y, Y_meas):
    
    Y = Y[Y[:, 0].argsort()]
    
    c = np.ceil(0.25*N)
    for i in range(0, N):
        if(c < i and i <= c*2):
            Y[i,0] = -1
            Y[i,1] = Y[int(i-c),1] * np.random.uniform(0.99,1.01)
            Y[i,2] = Y[int(i-c),2] * np.random.uniform(0.99,1.01)
            Y[i,3] = Y[int(i-c),3] * np.random.uniform(0.99,1.01)
        if(i >= c*2):
            Y[i, 1] = np.random.uniform(-60,60)
            Y[i, 2] = np.random.uniform(-25,25)
            Y[i, 3] = np.random.uniform(-1,1)

    for i in range(0, len(Y[:,0])):
        Y[i,0] = 0
        for ii in range(0, len(Y_meas[:,0])):
            Y[i,0] += J(Y[i,1], Y[i,2], Y[i,3], Y_meas[ii,0], Y_meas[ii,1], Y_meas[ii,2])

    Y = Y[Y[:, 0].argsort()]

    return Y


M = int(10e3)
c = 1000

E1 = np.zeros([int(M/c),4], dtype='float')
E2 = np.zeros([int(M/c),4], dtype='float')
E3 = np.zeros([int(M/c),4], dtype='float')

k = 0
for i in tqdm(range(0, M)):
    Y1 = Stepping(N, Y1, Y_meas1)
    Y2 = Stepping(N, Y2, Y_meas2)
    Y3 = Stepping(N, Y3, Y_meas3)
    
    if(i%c == 0):
        E1[k,:] = Y1[0,:]
        E2[k,:] = Y2[0,:]
        E3[k,:] = Y3[0,:]
        k += 1



E1[:,3] = (E1[:,3]/(np.pi*2*50))*1000
E2[:,3] = (E2[:,3]/(np.pi*2*50))*1000
E3[:,3] = (E3[:,3]/(np.pi*2*50))*1000

fig3, axs3 = plt.subplots(5,3)
x = np.linspace(0, int(M), int(M/c))
x1 = np.linspace(1, int(len(Y_meas1[:,0])), int(len(Y_meas1[:,0])))
axs3[0,0].set_title(f'0mH, 15 Ohm, {N*3} samples & {n} variables')
axs3[0,0].plot(x[1:], E1[1:,0], label=r'$Error$.' + f' (Best:{round(E1[-1,0],3)})')
axs3[0,0].set_ylabel('Error')
axs3[0,0].grid()
axs3[0,0].legend(loc=9)

axs3[1,0].plot(x, E1[:,1], label=r'$e_{g}.' + f' (Best:{round(E1[-1,1],3)})')
axs3[1,0].set_ylabel('Voltage [V]')
axs3[1,0].grid()
axs3[1,0].legend(loc=9)

axs3[2,0].plot(x, E1[:,2], label=r'$R.' + f' (Best:{round(E1[-1,2],3)})')
axs3[2,0].set_ylabel('Resistance [Ohm]')
axs3[2,0].grid()
axs3[2,0].legend(loc=9)

axs3[3,0].plot(x, E1[:,3], label=r'$L.' + f' (Best:{round(E1[-1,3],3)})')
axs3[3,0].grid()
axs3[3,0].legend(loc=9)
axs3[3,0].set_ylabel('Inductance [H]')
axs3[3,0].set_xlabel('Iterations')
 
axs3[4,0].plot(x1, Y_meas1[:,0], color='r', label=r'$V$')
#axs3[4,0].plot(x, Y_meas1[:,1], color='b', label=r'$P$')
axs3[4,0].plot(x1, Y_meas1[:,2], color='g', label=r'$Q$')
axs3[4,0].legend(loc=2)
axs3[4,0].set_ylim([-1,30])
axs3[4,0].grid()

#

axs3[0,1].set_title(f'4.33mH, 15 Ohm, {N*3} samples & {n} variables')
axs3[0,1].plot(x[1:], E2[1:,0], label=r'$Error$.' + f' (Best:{round(E2[-1,0],3)})')
axs3[0,1].grid()
axs3[0,1].legend(loc=9)

axs3[1,1].plot(x, E2[:,1], label=r'$e_{g}.' + f' (Best:{round(E2[-1,1],3)})')
axs3[1,1].grid()
axs3[1,1].legend(loc=9)

axs3[2,1].plot(x, E2[:,2], label=r'$R.' + f' (Best:{round(E2[-1,2],3)})')
axs3[2,1].grid()
axs3[2,1].legend(loc=9)

axs3[3,1].plot(x, E2[:,3], label=r'$L.' + f' (Best:{round(E2[-1,3],3)})')
axs3[3,1].grid()
axs3[3,1].legend(loc=9)
axs3[3,1].set_xlabel('Iterations')

axs3[4,1].plot(x1, Y_meas2[:,0], color='r', label=r'$V$')
#axs3[4,1].plot(x, Y_meas2[:,1], color='b', label=r'$P$')
axs3[4,1].plot(x1, Y_meas2[:,2], color='g', label=r'$Q$')
axs3[4,1].legend(loc=2)
axs3[4,1].set_ylim([-1,30])
axs3[4,1].grid()

#

axs3[0,2].set_title(f'8.66mH, 15 Ohm, {N*3} samples & {n} variables')
axs3[0,2].plot(x[1:], E3[1:,0], label=r'$Error$.' + f' (Best:{round(E3[-1,0],3)})')
axs3[0,2].grid()
axs3[0,2].legend(loc=9)

axs3[1,2].plot(x, E3[:,1], label=r'$e_{g}.' + f' (Best:{round(E3[-1,1],3)})')
axs3[1,2].grid()
axs3[1,2].legend(loc=9)

axs3[2,2].plot(x, E3[:,2], label=r'$R.' + f' (Best:{round(E3[-1,2],3)})')
axs3[2,2].grid()
axs3[2,2].legend(loc=9)

axs3[3,2].plot(x, E3[:,3], label=r'$L.' + f' (Best:{round(E3[-1,3],3)})')
axs3[3,2].grid()
axs3[3,2].legend(loc=9)
axs3[3,2].set_xlabel('Iterations')

axs3[4,2].plot(x1, Y_meas3[:,0], color='r', label=r'$V$')
#axs3[4,2].plot(x, Y_meas3[:,1], color='b', label=r'$P$')
axs3[4,2].plot(x1, Y_meas3[:,2], color='g', label=r'$Q$')
axs3[4,2].legend(loc=2)
axs3[4,2].set_ylim([-1,30])
axs3[4,2].grid()















#%%
M = int(100e3)
c = 1000.0

E1 = np.zeros([int(M),4], dtype='float')
E2 = np.zeros([int(M),4], dtype='float')
E3 = np.zeros([int(M),4], dtype='float')

k1 = 0
k2 = 0
k3 = 0

end = 0.1

cnt = 0
while(1):
    Y1 = Stepping(N, Y1, Y_meas1)
    if(cnt%c == 0):
        print(f'Y1: {cnt}, Error: {Y1[0,0]}')
        E1[k1,:] = Y1[0,:]
        k1 += 1
    if(E1[k1-1,0] <= end):
        k1 = k1 - 1
        break
    cnt += 1
x1 = np.linspace(0, cnt, k1)
cnt = 0
while(1):
    Y2 = Stepping(N, Y2, Y_meas2)
    if(cnt%c == 0):
        print(f'Y2: {cnt}, Error: {Y2[0,0]}')
        E2[k2,:] = Y2[0,:]
        k2 += 1
    if(E2[k2-1,0] <= end):
        k2 = k2 - 1
    cnt += 1
x2 = np.linspace(0, cnt, k2)
cnt = 0   
while(1):
    Y3 = Stepping(N, Y3, Y_meas3)
    if(cnt%c == 0):
        print(f'Y3: {cnt}, Error: {Y3[0,0]}')
        E3[k3,:] = Y3[0,:]
        k3 += 1
    if(E3[k3-1,0] <= end):
        k3 = k3 - 1
        break
    cnt += 1
x3 = np.linspace(0, cnt, k3)

# =============================================================================
# for i in tqdm(range(0, M)):
#     Y1 = Stepping(N, Y1, Y_meas1)
#     Y2 = Stepping(N, Y2, Y_meas2)
#     Y3 = Stepping(N, Y3, Y_meas3)
#     if(i%c == 0):
#         E1[k,:] = Y1[0,:]
#         E2[k,:] = Y2[0,:]
#         E3[k,:] = Y3[0,:]
#         k += 1
# =============================================================================




fig3, axs3 = plt.subplots(5,3)
x = np.linspace(0, N, N)
axs3[0,0].set_title(f'0mH, 15 Ohm, {N*3} samples & {n} variables')
axs3[0,0].plot(x1[1:], E1[1:k1,0], label=r'$Error$.' + f' (Best:{round(E1[k1,0],3)})')
axs3[0,0].set_ylabel('Error')
axs3[0,0].grid()
axs3[0,0].legend(loc=9)

axs3[1,0].plot(x1, E1[:k1,1], label=r'$e_{g}.' + f' (Best:{round(E1[k1,1],3)})')
axs3[1,0].set_ylabel('Voltage [V]')
axs3[1,0].grid()
axs3[1,0].legend(loc=9)

axs3[2,0].plot(x1, E1[:k1,2], label=r'$R.' + f' (Best:{round(E1[k1,2],3)})')
axs3[2,0].set_ylabel('Resistance [Ohm]')
axs3[2,0].grid()
axs3[2,0].legend(loc=9)

axs3[3,0].plot(x1, E1[:k1,3], label=r'$L.' + f' (Best:{round(E1[k1,3],3)})')
axs3[3,0].grid()
axs3[3,0].legend(loc=9)
axs3[3,0].set_ylabel('Inductance [H]')
axs3[3,0].set_xlabel('Iterations')
 
axs3[4,0].plot(x, Y_meas1[:,0], color='r', label=r'$V$')
#axs3[4,0].plot(x, Y_meas1[:,1], color='b', label=r'$P$')
axs3[4,0].plot(x, Y_meas1[:,2], color='g', label=r'$Q$')
axs3[4,0].legend(loc=2)
axs3[4,0].set_ylim([-1,30])
axs3[4,0].grid()

#

axs3[0,1].set_title(f'4.33mH, 15 Ohm, {n} samples & {N*3} variables')
axs3[0,1].plot(x2, E2[1:k2,0], label=r'$Error$.' + f' (Best:{round(E2[k2,0],3)})')
axs3[0,1].grid()
axs3[0,1].legend(loc=9)

axs3[1,1].plot(x2, E2[:k2,1], label=r'$e_{g}.' + f' (Best:{round(E2[k2,1],3)})')
axs3[1,1].grid()
axs3[1,1].legend(loc=9)

axs3[2,1].plot(x2, E2[:k2,2], label=r'$R.' + f' (Best:{round(E2[k2,2],3)})')
axs3[2,1].grid()
axs3[2,1].legend(loc=9)

axs3[3,1].plot(x2, E2[:k2,3], label=r'$L.' + f' (Best:{round(E2[k2,3],3)})')
axs3[3,1].grid()
axs3[3,1].legend(loc=9)
axs3[3,1].set_xlabel('Iterations')

axs3[4,1].plot(x, Y_meas2[:,0], color='r', label=r'$V$')
#axs3[4,1].plot(x, Y_meas2[:,1], color='b', label=r'$P$')
axs3[4,1].plot(x, Y_meas2[:,2], color='g', label=r'$Q$')
axs3[4,1].legend(loc=2)
axs3[4,1].set_ylim([-1,30])
axs3[4,1].grid()

#

axs3[0,2].set_title(f'8.66mH, 15 Ohm, {n} samples & {N*3} variables')
axs3[0,2].plot(x3, E3[0:k3,0], label=r'$Error$.' + f' (Best:{round(E3[k3,0],3)})')
axs3[0,2].grid()
axs3[0,2].legend(loc=9)

axs3[1,2].plot(x3, E3[:k3,1], label=r'$e_{g}.' + f' (Best:{round(E3[k3,1],3)})')
axs3[1,2].grid()
axs3[1,2].legend(loc=9)

axs3[2,2].plot(x3, E3[:k3,2], label=r'$R.' + f' (Best:{round(E3[k3,2],3)})')
axs3[2,2].grid()
axs3[2,2].legend(loc=9)

axs3[3,2].plot(x3, E3[:k3,3], label=r'$L.' + f' (Best:{round(E3[k3,3],3)})')
axs3[3,2].grid()
axs3[3,2].legend(loc=9)
axs3[3,2].set_xlabel('Iterations')

axs3[4,2].plot(x, Y_meas3[:,0], color='r', label=r'$V$')
#axs3[4,2].plot(x, Y_meas3[:,1], color='b', label=r'$P$')
axs3[4,2].plot(x, Y_meas3[:,2], color='g', label=r'$Q$')
axs3[4,2].legend(loc=2)
axs3[4,2].set_ylim([-1,30])
axs3[4,2].grid()