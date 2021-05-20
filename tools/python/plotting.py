import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

tracercount = np.array(pd.read_csv("test/tracercount_1.00e-03.csv",delimiter=','))

# tracercount = np.array(tracercount)
plt.figure(num=2,figsize=(6,6),dpi=150)
plt.plot(tracercount[:,0],tracercount[:,1],'r-')
plt.plot(tracercount[:,0],tracercount[:,2],'b-.')
plt.xlim([0,1])
#plt.ylim([0,1000])
plt.xlabel('Time [s]')
plt.ylabel('MSD [$\mathrm{m^2}$]')
#plt.title('D=%.0e, $\Delta$t=%.0e' %(D,delt))
plt.savefig("test/wowie.png" ,bbox_inches='tight',format='png')
plt.clf()
