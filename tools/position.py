import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

dir_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
dir_path = dir_path + "/sim/output/"

img_dir = dir_path + 'images'
if not os.path.exists(img_dir):
  os.makedirs(img_dir)
  
tracerpath = dir_path + "tracer/"
# n = len([f for f in os.listdir(tracerpath) 
#      if f.endswith('.csv') and os.path.isfile(os.path.join(tracerpath, f))])

timedata = np.array(pd.read_csv(dir_path + "tracerinfo.csv",header=None))
n = len(timedata)

extents = np.array(pd.read_csv(dir_path+'domaindata.csv',header=None))
extents = extents[0]

bins = np.linspace(extents[0],extents[1],num=100, endpoint=True)

def plotter(i,type0array,type1array):

    fig = plt.figure(num=1,figsize=(6,12))
    ax = fig.add_subplot(211)
    plt.hist(type0array,bins=bins)
    ax = fig.add_subplot(212)
    plt.hist(type1array,bins=bins)
    plt.savefig(img_dir+"/hist_%d.png" %i, bbox_inches='tight', format='png')
    plt.clf()
    
type0array = []
type1array = []

for j in range(n-50000,n):
    i = j
    if j%100 == 0 : print(j,"/",n)
    csvname = dir_path + "tracer/tracer-%d.csv" %i
    stuff = np.array(pd.read_csv(csvname))
    loc_array = stuff[:,1:4]
    type_array = stuff[:,4]

    for x in range(len(stuff[:,4])):
        if stuff[x,4]:
            type1array.append(loc_array[x,0])
        else:
            type0array.append(loc_array[x,0])

type0array = np.array(type0array) 
type1array = np.array(type1array)

plotter(i,type0array,type1array)

