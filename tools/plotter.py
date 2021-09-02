import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

def set_axes_equal(ax):

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

dir_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
dir_path = dir_path + "/sim/output/"

img_dir = dir_path + 'images'
if not os.path.exists(img_dir):
  os.makedirs(img_dir)
  
tracerpath = dir_path + "tracer/"
n = len([f for f in os.listdir(tracerpath) 
     if f.endswith('.csv') and os.path.isfile(os.path.join(tracerpath, f))])


timedata = np.array(pd.read_csv(dir_path + "tracerinfo.csv",header=None))

def plotter(time,loc_array,extents,n,color_array):
    part = np.array(pd.read_csv(dir_path + "particleinfo/particle-%d.csv" %time, header=None))
    nparts = len(part)

    fig = plt.figure(num=1,figsize=(14,9))
    ax = fig.add_subplot(111,projection='3d')

    for k in range(nparts):
      [cx,cy,cz,r]=part[k,1:]
      u = np.linspace(0, 2 * np.pi, 100)
      v = np.linspace(0, np.pi, 100)
      x = ( r * np.outer(np.cos(u), np.sin(v)) + cx)
      y = ( r * np.outer(np.sin(u), np.sin(v)) + cy)
      z = (r * np.outer(np.ones(np.size(u)), np.cos(v)) + cz)
      ax.plot_surface(x, y, z, color='b')

    X = np.array([-1,1])
    Y = np.array([-0.5,0.5])
    Z = np.array([-0.5,0.5])

    max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max()
    Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X.max()+X.min())
    Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y.max()+Y.min())
    Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Z.max()+Z.min())

    for xb, yb, zb in zip(Xb, Yb, Zb):
        ax.plot([xb], [yb], [zb], 'w')

    ax.scatter(loc_array[:,0],loc_array[:,1],loc_array[:,2],alpha=1,c=color_array,s=2)

    ax.set_xlim3d(left=extents[0],right=extents[1])
    ax.set_ylim3d(bottom=extents[2],top=extents[3])
    ax.set_zlim3d(bottom=extents[4],top=extents[5])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.title('t=%.3f' %timedata[int(time/100),1])
    plt.savefig(img_dir+"/%d.png" %time, bbox_inches='tight', format='png')
    plt.clf()
    
extents = np.array(pd.read_csv(dir_path+'domaindata.csv',header=None))
extents = extents[0]

for j in range(n-1,n):
    i = j
    if j%10 == 0 : print(j,"/",n)
    csvname = dir_path + "tracer/tracer-%d.csv" %i
    stuff = np.array(pd.read_csv(csvname))
    type_array = stuff[:,4]

    color_array = []
    for x in type_array:
        if x==0:
            color_array.append('r')
        else:
            color_array.append('b')

    loc_array = stuff[:,1:4]

    plotter(i,loc_array,extents,n,color_array)

