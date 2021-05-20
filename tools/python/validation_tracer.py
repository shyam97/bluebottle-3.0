import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

case1 = np.array(pd.read_csv('case1.csv',delimiter=','))
case2 = np.array(pd.read_csv('case2.csv',delimiter=','))
case3 = np.array(pd.read_csv('case3.csv',delimiter=','))
case4 = np.array(pd.read_csv('case4.csv',delimiter=','))
case5 = np.array(pd.read_csv('case5.csv',delimiter=','))
case6 = np.array(pd.read_csv('case6.csv',delimiter=','))
case7 = np.array(pd.read_csv('case7.csv',delimiter=','))
case8 = np.array(pd.read_csv('case8.csv',delimiter=','))
case9 = np.array(pd.read_csv('case9.csv',delimiter=','))

from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='r', lw=2),
                Line2D([0], [0], color='g', lw=2),
                Line2D([0], [0], color='b', lw=2),
                Line2D([0], [0], color='k', linestyle='dotted', lw=2),
                Line2D([0], [0], color='k', linestyle = '-.', lw=2),
                Line2D([0], [0], color='k', linestyle='-', lw=2)]

plt.figure(num=1,figsize=(8,8))
plt.plot(case3[:,0],case3[:,1],'b-')
plt.plot(case6[:,0],case6[:,1],'g-')
plt.plot(case9[:,0],case9[:,1],'r-')
plt.plot(case2[:,0],case2[:,1],'b-.')
plt.plot(case5[:,0],case5[:,1],'g-.')
plt.plot(case8[:,0],case8[:,1],'r-.')
plt.plot(case1[:,0],case1[:,1],'b.')
plt.plot(case4[:,0],case4[:,1],'g.')
plt.plot(case7[:,0],case7[:,1],'r.')
plt.xlim([0,2])
plt.ylim([0,1000])
plt.legend(custom_lines,['D = 5e-3 $\mathrm{m^2/s}$',
            'D = 5e-2 $\mathrm{m^2/s}$','D = 5e-1 $\mathrm{m^2/s}$',
            '$\Delta$t = 1e-2 s','$\Delta$t = 1e-3 s','$\Delta$t = 1e-4 s'])
plt.xlabel('Time [s]')
plt.ylabel('Tracers in domain [-]')
plt.title('Variation in D and $\Delta$t')
plt.savefig('variation.png',bbox_inches='tight',format='png')
