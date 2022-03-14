import matplotlib.pyplot as plt
import numpy as np
import cmath
from pandas import read_csv
from scipy import signal

from matplotlib.ticker import FormatStrFormatter

import warnings

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

df1_1 = read_csv('Phase1_1.CSV')
df2_1 = read_csv('Phase2_1.CSV')
df3_1 = read_csv('Phase3_1.CSV')

df1_2 = read_csv('Phase1_2.CSV')
df2_2 = read_csv('Phase2_2.CSV')
df3_2 = read_csv('Phase3_2.CSV')

df1_3 = read_csv('Phase1_3.CSV')
df2_3 = read_csv('Phase2_3.CSV')
df3_3 = read_csv('Phase3_3.CSV')

n = 10

E_ab_1 = df1_1['Channel 1 (V)'].rolling(n).mean()
E_ac_1 = df2_1['Channel 1 (V)'].rolling(n).mean()
E_bc_1 = df3_1['Channel 1 (V)'].rolling(n).mean()

E_ab_2 = df1_2['Channel 1 (V)'].rolling(n).mean()
E_ac_2 = df2_2['Channel 1 (V)'].rolling(n).mean()
E_bc_2 = df3_2['Channel 1 (V)'].rolling(n).mean()

E_ab_3 = df1_3['Channel 1 (V)'].rolling(n).mean()
E_ac_3 = df2_3['Channel 1 (V)'].rolling(n).mean()
E_bc_3 = df3_3['Channel 1 (V)'].rolling(n).mean()

E_a_1 = (E_ab_1 - (-E_ac_1))/3
E_b_1 = (E_bc_1 - E_ab_1)/3
E_c_1 = ((-E_ac_1) - E_bc_1)/3

E_a_2 = (E_ab_2 - (-E_ac_2))/3
E_b_2 = (E_bc_2 - E_ab_2)/3
E_c_2 = ((-E_ac_2) - E_bc_2)/3

E_a_3 = (E_ab_3 - (-E_ac_3))/3
E_b_3 = (E_bc_3 - E_ab_3)/3
E_c_3 = ((-E_ac_3) - E_bc_3)/3


I_a_1 = df1_1['Channel 2 (V)'].rolling(n).mean()
I_b_1 = df2_1['Channel 2 (V)'].rolling(n).mean()
I_c_1 = df3_1['Channel 2 (V)'].rolling(n).mean()

I_a_2 = df1_2['Channel 2 (V)'].rolling(n).mean()
I_b_2 = df2_2['Channel 2 (V)'].rolling(n).mean()
I_c_2 = df3_2['Channel 2 (V)'].rolling(n).mean()

I_a_3 = df1_3['Channel 2 (V)'].rolling(n).mean()
I_b_3 = df2_3['Channel 2 (V)'].rolling(n).mean()
I_c_3 = df3_3['Channel 2 (V)'].rolling(n).mean()

Q_1 = ((1/np.sqrt(3)) * ((I_a_1 *(E_c_1 - E_b_1)) + (I_b_1 *(E_a_1 - E_c_1)) + (I_c_1 * (E_b_1 - E_a_1))))
Q_2 = ((1/np.sqrt(3)) * ((I_a_2 *(E_c_2 - E_b_2)) + (I_b_2 *(E_a_2 - E_c_2)) + (I_c_2 * (E_b_2 - E_a_2))))
Q_3 = ((1/np.sqrt(3)) * ((I_a_3 *(E_c_3 - E_b_3)) + (I_b_3 *(E_a_3 - E_c_3)) + (I_c_3 * (E_b_3 - E_a_3))))

P_1 = ((E_a_1 * I_a_1) + (E_b_1 * I_b_1) + (E_c_1 * I_c_1))
P_2 = ((E_a_2 * I_a_2) + (E_b_2 * I_b_2) + (E_c_2 * I_c_2))
P_3 = ((E_a_3 * I_a_3) + (E_b_3 * I_b_3) + (E_c_3 * I_c_3))

fig1, axs1 = plt.subplots(3,3)
fig1.tight_layout()

axs1[0,0].plot(df1_1['Time (s)'].rolling(n).mean(), E_ab_1, color='r', label='U_ab')
axs1[0,0].plot(df2_1['Time (s)'].rolling(n).mean(), -E_ac_1, color='g', label='U_ca')
axs1[0,0].plot(df3_1['Time (s)'].rolling(n).mean(), E_bc_1, color='b', label='U_bc')

axs1[0,0].plot(df1_1['Time (s)'].rolling(n).mean(), E_a_1, alpha=0.2, color='r', label='U_a')
axs1[0,0].plot(df2_1['Time (s)'].rolling(n).mean(), E_b_1, alpha=0.2, color='g', label='U_b')
axs1[0,0].plot(df3_1['Time (s)'].rolling(n).mean(), E_c_1, alpha=0.2, color='b', label='U_c')

axs1[0,0].set_ylabel('Amplitude [V]')
axs1[0,0].grid()
axs1[0,0].set_title(r'Load: $15 \Omega$')
axs1[0,0].legend(loc=5)

axs1[1,0].plot(df1_1['Time (s)'].rolling(n).mean(), I_a_1, color='r', label='I_a')
axs1[1,0].plot(df2_1['Time (s)'].rolling(n).mean(), I_b_1, color='g', label='I_b')
axs1[1,0].plot(df3_1['Time (s)'].rolling(n).mean(), I_c_1, color='b', label='I_c')
axs1[1,0].set_ylabel('Amplitude [A]')
axs1[1,0].grid()
axs1[1,0].legend(loc=5)

axs1[2,0].plot(df1_1['Time (s)'].rolling(n).mean(), Q_1, label='Q')
axs1[2,0].plot(df1_1['Time (s)'].rolling(n).mean(), P_1, label='P')
axs1[2,0].set_ylabel('Power [W/VAr]')
axs1[2,0].set_xlabel('Time [s]')
axs1[2,0].grid()
axs1[2,0].legend(loc=5)



axs1[0,1].plot(df1_2['Time (s)'].rolling(n).mean(), E_ab_2, color='r', label='U_ab')
axs1[0,1].plot(df2_2['Time (s)'].rolling(n).mean(), -E_ac_2, color='g', label='U_ca')
axs1[0,1].plot(df3_2['Time (s)'].rolling(n).mean(), E_bc_2, color='b', label='U_bc')

axs1[0,1].plot(df1_2['Time (s)'].rolling(n).mean(), E_a_2, alpha=0.2, color='r', label='U_a')
axs1[0,1].plot(df2_2['Time (s)'].rolling(n).mean(), E_b_2, alpha=0.2, color='g', label='U_b')
axs1[0,1].plot(df3_2['Time (s)'].rolling(n).mean(), E_c_2, alpha=0.2, color='b', label='U_c')
axs1[0,1].set_title(r'Load: $15 \Omega + 4.33\,mH$')
axs1[0,1].grid()
axs1[0,1].legend(loc=5)

axs1[1,1].plot(df1_2['Time (s)'].rolling(n).mean(), I_a_2, color='r', label='I_a')
axs1[1,1].plot(df2_2['Time (s)'].rolling(n).mean(), I_b_2, color='g', label='I_b')
axs1[1,1].plot(df3_2['Time (s)'].rolling(n).mean(), I_c_2, color='b', label='I_c')
axs1[1,1].grid()
axs1[1,1].legend(loc=5)

axs1[2,1].plot(df1_2['Time (s)'].rolling(n).mean(), Q_2, label='Q')
axs1[2,1].plot(df1_2['Time (s)'].rolling(n).mean(), P_2, label='P')
axs1[2,1].set_xlabel('Time [s]')
axs1[2,1].grid()
axs1[2,1].legend(loc=5)


axs1[0,2].plot(df1_3['Time (s)'].rolling(n).mean(), E_ab_3, color='r', label='U_ab')
axs1[0,2].plot(df2_3['Time (s)'].rolling(n).mean(), -E_ac_3, color='g', label='U_ca')
axs1[0,2].plot(df3_3['Time (s)'].rolling(n).mean(), E_bc_3, color='b', label='U_bc')

axs1[0,2].plot(df1_3['Time (s)'].rolling(n).mean(), E_a_3, alpha=0.2, color='r', label='U_a')
axs1[0,2].plot(df2_3['Time (s)'].rolling(n).mean(), E_b_3, alpha=0.2, color='g', label='U_b')
axs1[0,2].plot(df3_3['Time (s)'].rolling(n).mean(), E_c_3, alpha=0.2, color='b', label='U_c')
axs1[0,2].set_title(r'Load: $15 \Omega + 8.66\,mH$')
axs1[0,2].grid()
axs1[0,2].legend(loc=5)

axs1[1,2].plot(df1_3['Time (s)'].rolling(n).mean(), I_a_3, color='r', label='I_a')
axs1[1,2].plot(df2_3['Time (s)'].rolling(n).mean(), I_b_3, color='g', label='I_b')
axs1[1,2].plot(df3_3['Time (s)'].rolling(n).mean(), I_c_3, color='b', label='I_c')
axs1[1,2].grid()
axs1[1,2].legend(loc=5)

axs1[2,2].plot(df1_3['Time (s)'].rolling(n).mean(), Q_3, label='Q')
axs1[2,2].plot(df1_3['Time (s)'].rolling(n).mean(), P_3, label='P')
axs1[2,2].set_xlabel('Time [s]')
axs1[2,2].grid()
axs1[2,2].legend(loc=5)


