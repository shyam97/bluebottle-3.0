import numpy as np
from numpy.core.fromnumeric import trace
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os, csv

dir_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
dir_path = dir_path + "/sim/output/"

img_dir = dir_path + 'images'
if not os.path.exists(img_dir):
  os.makedirs(img_dir)
  
tracerpath = dir_path + "tracer/"
particlepath = dir_path + "particleinfo/"
# n = len([f for f in os.listdir(tracerpath) 
#      if f.endswith('.csv') and os.path.isfile(os.path.join(tracerpath, f))])

nparts = len(np.array(pd.read_csv(particlepath + "particle-0.csv",header=None)))
print(nparts)

ntracers = len(np.array(pd.read_csv(tracerpath + "tracer-0.csv",header=None)))
print(ntracers)

timedata = np.array(pd.read_csv(dir_path + "tracerinfo.csv",header=None))
n = len(timedata)

extents = np.array(pd.read_csv(dir_path+'domaindata.csv',header=None))
extents = extents[0]

wallarea = (extents[1]-extents[0])*(extents[5]-extents[4])
domvol = (extents[1]-extents[0])*(extents[3]-extents[2])*(extents[5]-extents[4])
concforc = ntracers/(domvol)
H = (extents[3]-extents[2])
D = 1e-4
sherwoodlog=[]
countab = []
countba = []
xcount = []
ycount1 = []
ycount2 = []
countcu = 0

for i in range(10000,n-1):

    xcount.append(timedata[i,1])

    count01 = 0
    count10 = 0
    
    data1 = np.array(pd.read_csv(tracerpath + "tracer-%d.csv" %i,header=None))
    data2 = np.array(pd.read_csv(tracerpath + "tracer-%d.csv" %(i+1), header=None))

    for j in range(ntracers):
        type_now = data1[j,4]
        type_nex = data2[j,4]

        if type_now==0 and type_nex==1:
            ycount1.append(timedata[i,1])
            count01+=1
            countcu += 0.5
        
        if type_now==1 and type_nex==0:
            ycount2.append(timedata[i,1])
            count10+=1
            countcu += 0.5

    countab.append(count01)
    countba.append(count10)
    # print(countab,countba)

    massflux = countcu/(xcount[-1] - xcount[0])

    h = massflux/(wallarea*concforc)
    sherwood = h*H/D
    sherwoodlog.append(sherwood)

    if i%100 == 0: 
        print(i,"/",n, "\tTime=%.4f" %timedata[i,1],"\tSh=%.4f" %sherwood, "\tCount=%d" %(countcu*2))

fig = plt.figure(num=1,figsize=(14,9))
ax = fig.add_subplot(211)
ax.scatter(xcount,countab)
ax = fig.add_subplot(212)
ax.scatter(xcount,countba)
plt.savefig(img_dir+"/counts.png", bbox_inches='tight', format='png')
plt.clf()

bins = np.linspace(xcount[0],xcount[-1],num=100, endpoint=True)

fig = plt.figure(num=1,figsize=(6,12))
ax = fig.add_subplot(211)
plt.hist(ycount1,bins=bins)
ax = fig.add_subplot(212)
plt.hist(ycount2,bins=bins)
plt.savefig(img_dir+"/hist_sherwood.png", bbox_inches='tight', format='png')
plt.clf()

fig = plt.figure(num=1,figsize=(6,6))
ax = fig.add_subplot(211)
plt.plot(xcount,sherwoodlog)
ax.set_xlabel("Time")
ax.set_ylabel("Sherwood Number")
ax.set_ylim([0,20])
plt.savefig(img_dir+"/sherwood.png",bbox_inches='tight',format='png')
plt.clf()

header = ["Cumulative Count", "Time of Count", "Mass Flux", "Wall Area", "Number of Tracers", "Domain Volume",
           "Driving Force", "h", "H", "D", "Sherwood"]
data = [countcu, xcount[-1] - xcount[0], massflux, wallarea, ntracers, domvol, concforc, h, H, D, sherwood]

with open(img_dir+'/sherwood.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(header)
    writer.writerow(data)

print("\n\nCumulative count =",countcu)
print("Total calculation time =",xcount[-1] - xcount[0])
print("Mass flux =",massflux)

print("Wall area =", wallarea)

print("Number of tracers =",ntracers)
print("Domain volume =",domvol)
print("Concentration driving force =",concforc)

print('h =',h)
print('H =',H)
print('D =',D)

print('Sh =',sherwood)

