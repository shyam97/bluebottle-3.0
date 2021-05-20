# -*- coding: utf-8 -*-
"""tracer.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1O7BmWN8xUeJgvZt-PigJpvQiTyMav4OW
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import bluebottle_flow_reader as bbflow
import bluebottle_particle_reader as bbparts
import sys,os

# data_dir = sys.argv[1]

img_dir = 'images_sphere'
if not os.path.exists(img_dir):
  os.makedirs(img_dir)

# times = bbflow.init(data_dir)
# grid = bbflow.read_flow_extents(data_dir)
# [coordx,coordy,coordz] = bbflow.read_flow_position(data_dir)

# D = 2.0204e-9
D = 5e-2
delt = 1e-3

u = np.zeros((201,201,201))
v = np.zeros((201,201,201))
w = np.zeros((201,201,201))

locp = [0,0,0]
r = 0.05

theta = np.linspace(0,2*np.pi,num=100)
cx = r*np.cos(theta)
cy = r*np.sin(theta)

coordx = np.linspace(-0.1,0.1,num=201)
coordy = np.linspace(-0.1,0.1,num=201)
coordz = np.linspace(-0.1,0.1,num=201)
# print(coordx.shape)
# print(coordy.shape)
# print(coordz.shape)
# print(grid)
# extents = [grid[3],grid[4],grid[6],grid[7],grid[9],grid[10]]
extents = [coordx[0],coordx[-1],coordy[0],coordy[-1],coordz[0],coordz[-1]]
# print(extents)

"""#Tracer Class"""

class tracers:
  def __init__(self,extents):
    xmin,xmax,ymin,ymax,zmin,zmax=extents
    x = xmin
    # y = ymin + (ymax-ymin)*np.random.random()
    # y = (ymin+ymax)/2 + (ymax-ymin)*(np.random.random()-0.5)/2
    y = 0
    # z = zmin + (zmax-zmin)*np.random.random()
    # z = (zmin+zmax)/2 + (zmax-zmin)*(np.random.random()-0.5)/2
    z = (zmin+zmax)/2
    # print([x,y,z])
    self.loc = [x,y,z]

  def updater(self,delt,vel0,rand0):
    pos1 = self.loc + delt*vel0 + rand0
    return pos1

  def finalize(self,pos):
      self.loc=pos

"""#Useful Functions"""

def wall_checker(loc,extents):
  x,y,z = loc
  xmin,xmax,ymin,ymax,zmin,zmax=extents

  if x<xmin:
    [x,y,z] = [xmin+(xmin-x),y,z]
  elif x>xmax:
    [x,y,z] = [xmax-(x-xmax),y,z]
    # return True
  if y<ymin:
    [x,y,z] = [x,ymin+(ymin-y),z]
  elif y>ymax:
    [x,y,z] = [x,ymax-(y-ymax),z]
  if z<zmin:
    [x,y,z] = [x,y,zmin+(zmin-z)]
  elif z>zmax:
    [x,y,z] = [x,y,zmax-(z-zmax)]
  return [x,y,z]

#-------------------------------------------------------------------------------

def sphere_checker(loc,px,py,pz,r):
  for l in range(0,len(px)):
    if distance(loc,px[l],py[l],pz[l]) <= r:
      return True
    else:
      return False

def sphere_checker_single(loc,locp,r):
    if distance(loc,locp[0],locp[1],locp[2]) <= r:
        return True
    else:
        return False

#-------------------------------------------------------------------------------

def distance(loc,x,y,z):
  return ((loc[0]-x)**2 + (loc[1]-y)**2 +(loc[2]-z)**2)**0.5

#-------------------------------------------------------------------------------

def index_finder(loc,coordx,coordy,coordz):

  for xnear in coordx:
    if xnear>loc[0]:
      xnear = np.squeeze(np.where(coordx == xnear)) - 1
      break
  for ynear in coordy:
    if ynear>loc[1]:
      ynear = np.squeeze(np.where(coordy == ynear)) - 1
      break
  for znear in coordz:
    if znear>loc[2]:
      znear = np.squeeze(np.where(coordz == znear)) - 1
      break
  return xnear, ynear, znear

#-------------------------------------------------------------------------------

def interpolate(loc,u,v,w,coordx,coordy,coordz):
  xnear, ynear, znear = index_finder(loc,coordx,coordy,coordz)
  px,py,pz = loc

  hx = coordx[xnear+1]-coordx[xnear]
  hy = coordy[ynear+1]-coordy[ynear]
  hz = coordz[znear+1]-coordz[znear]

  cx = (px - coordx[xnear])/(coordx[xnear+1]-coordx[xnear])
  cy = (py - coordy[ynear])/(coordy[ynear+1]-coordy[ynear])
  cz = (pz - coordz[znear])/(coordz[znear+1]-coordz[znear])

  tempx1 = (1-cx)*u[xnear,ynear,znear] + cx*u[xnear+1,ynear,znear]
  tempx2 = (1-cx)*u[xnear,ynear+1,znear] + cx*u[xnear+1,ynear+1,znear]
  tempx3 = (1-cx)*u[xnear,ynear,znear+1] + cx*u[xnear+1,ynear,znear+1]
  tempx4 = (1-cx)*u[xnear,ynear+1,znear+1] + cx*u[xnear+1,ynear+1,znear+1]
  tempy1 = (1-cy)*tempx1 + cy*tempx2
  tempy2 = (1-cy)*tempx3 + cy*tempx4
  resultu = (1-cz)*tempy1 + cz*tempy2

  tempx1 = (1-cx)*v[xnear,ynear,znear] + cx*v[xnear+1,ynear,znear]
  tempx2 = (1-cx)*v[xnear,ynear+1,znear] + cx*v[xnear+1,ynear+1,znear]
  tempx3 = (1-cx)*v[xnear,ynear,znear+1] + cx*v[xnear+1,ynear,znear+1]
  tempx3 = (1-cx)*v[xnear,ynear+1,znear+1] + cx*v[xnear+1,ynear+1,znear+1]
  tempy1 = (1-cy)*tempx1 + cy*tempx2
  tempy2 = (1-cy)*tempx3 + cy*tempx4
  resultv = (1-cz)*tempy1 + cz*tempy2

  tempx1 = (1-cx)*w[xnear,ynear,znear] + cx*w[xnear+1,ynear,znear]
  tempx2 = (1-cx)*w[xnear,ynear+1,znear] + cx*w[xnear+1,ynear+1,znear]
  tempx3 = (1-cx)*w[xnear,ynear,znear+1] + cx*w[xnear+1,ynear,znear+1]
  tempx3 = (1-cx)*w[xnear,ynear+1,znear+1] + cx*w[xnear+1,ynear+1,znear+1]
  tempy1 = (1-cy)*tempx1 + cy*tempx2
  tempy2 = (1-cy)*tempx3 + cy*tempx4
  resultw = (1-cz)*tempy1 + cz*tempy2

  return np.array([resultu, resultv, resultw])

#-------------------------------------------------------------------------------

def reflector(loc1,loc2,locp,r):
  if distance(locp,loc2[0],loc2[1],loc2[2]) >=r:
      return loc2
  else:
      deltax = loc2[0]-loc1[0]
      deltay = loc2[1]-loc1[1]
      deltaz = loc2[2]-loc1[2]

      locx = loc1
      delta = 10

      while delta<1e9:

        locx_copy = locx.copy()
        locx = locx + np.array([deltax,deltay,deltaz])/delta

        if distance(locp,locx[0],locx[1],locx[2]) <=r:
          locx = locx_copy.copy()
          delta *= 10
      # print(locx)
      dsmag = np.linalg.norm(loc2-locx)
      di = (loc2 - locx)/dsmag
      dn = (locx - locp)/np.linalg.norm(locx-locp)

      ds = di - 2 * np.dot(dn,di) * dn
      return ds*dsmag + locx

#-------------------------------------------------------------------------------

def randomizer(D,delt):
  magnitude = (6*D*delt)**0.5
  argument1 = np.pi*(1-2*np.random.random())
  argument2 = np.pi*(1-2*np.random.random())
  xcomponent = magnitude*np.cos(argument1)*1#np.sin(argument2)
  ycomponent = magnitude*np.sin(argument1)*1#np.sin(argument2)
  zcomponent = magnitude*np.cos(argument2)*0
  return np.array([xcomponent,ycomponent,zcomponent])

#-------------------------------------------------------------------------------

def plotter(time,loc_array,extents):
    fig = plt.figure(num=1,figsize=(8,8))
    ax = fig.add_subplot(111,projection='3d')
    ax.scatter(loc_array[:,0],loc_array[:,1],loc_array[:,2],c='red',s=1,alpha=1,
        zorder=-1)
    ax.plot(cx,cy,'b-')
    ax.set_xlim3d(left=extents[0],right=extents[1])
    ax.set_ylim3d(bottom=extents[2],top=extents[3])
    ax.set_zlim3d(bottom=extents[4],top=extents[5])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.title("Time = %.3f" %time)
    plt.savefig(img_dir+"/t=%.3f.png" %time, bbox_inches='tight', format='png')
    fig.clf()

"""#Main"""

"""
1. Interpolate velocities
2. Update position
    2.a. Check wall reflection
    2.b. Check sphere reflection
3. Finalize position

"""

tracer_array = [tracers(extents) for i in range(1000)]

loc_array = []
for t in tracer_array:
    loc_array.append(t.loc)
loc_array = np.array(loc_array)
plotter(0,loc_array,extents)

tracercount=[]

time = 0
count = 0
while time<1:
    time+=delt
    i = 0
    while i < len(tracer_array):
        t = tracer_array[i]
        vel0 = interpolate(t.loc,u,v,w,coordx,coordy,coordz)
        rand0 = randomizer(D,delt)
        newloc = t.updater(delt,vel0,rand0)

        if wall_checker(newloc,extents) == True:
            del tracer_array[i]
            i-=1
        else:
            newloc = wall_checker(newloc,extents)
            newloc = reflector(t.loc,newloc,locp,r)
            t.finalize(newloc)

        i+=1

    if len(tracer_array)==0:
        print("All particles have exited the domain.")
        sys.exit()

    if time>=(count+1)*0.005:
        loc_array = []
        for t in tracer_array:
            loc_array.append(t.loc)
        loc_array = np.array(loc_array)
        plotter(time,loc_array,extents)
        count+=1

    print('Time = %.3f, tracers = %d.' %(time,len(tracer_array)), end="\r")
    tracercount.append([time,len(tracer_array)])

tracercount = np.array(tracercount)
plt.figure(num=1)
plt.plot(tracercount[:,0],tracercount[:,1])
plt.savefig(img_dir+"/count_%.2e.png" %delt ,bbox_inches='tight',format='png')
plt.clf()
np.savetxt(img_dir+'/tracercount_%.2e.csv' %delt,tracercount,delimiter=',')
print("Simulation has reached its end time.")
