import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
  
dir_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
dir_path = dir_path + "/sim/output"
# print(dir_path)

img_dir = dir_path + '/images'
if not os.path.exists(img_dir):
  os.makedirs(img_dir)
  
myPath = dir_path+"/tracer"
n = len([f for f in os.listdir(myPath) 
     if f.endswith('.csv') and os.path.isfile(os.path.join(myPath, f))])
#print(n)

timedata = np.array(pd.read_csv(dir_path+"/tracerinfo.csv",header=None))

def plotter(time,extents):
    part = np.array(pd.read_csv(dir_path+"/particleinfo/particle-%d.csv" %time, header=None))
    nparts = len(part)

    fig = plt.figure(num=1,figsize=(8,8))
    ax = fig.add_subplot(111,projection='3d')

    for k in range(nparts):
      [cx,cy,cz,r]=part[k,1:]
      u = np.linspace(0, 2 * np.pi, 100)
      v = np.linspace(0, np.pi, 100)
      x = ( r * np.outer(np.cos(u), np.sin(v)) + cx)
      y = ( r * np.outer(np.sin(u), np.sin(v)) + cy)
      z = (r * np.outer(np.ones(np.size(u)), np.cos(v)) + cz)
      ax.plot_surface(x, y, z, color='b')

    ax.set_xlim3d(left=extents[0],right=extents[1])
    ax.set_ylim3d(bottom=extents[2],top=extents[3])
    ax.set_zlim3d(bottom=extents[4],top=extents[5])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.savefig(img_dir+"/parts.png", bbox_inches='tight', format='png')
    plt.clf()
    
extents = np.array(pd.read_csv(dir_path+'/domaindata.csv',header=None))
extents = extents[0]

plotter(0,extents)