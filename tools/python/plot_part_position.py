################################################################################
################################## BLUEBOTTLE ##################################
################################################################################
#
#  Copyright 2012 - 2018 Adam Sierakowski and Daniel Willen,
#                         The Johns Hopkins University
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
#  Please contact the Johns Hopkins University to use Bluebottle for
#  commercial and/or for-profit applications.
################################################################################

#!/usr/bin/env python

# PURPOSE:
#   Pull the particle position time series from Bluebottle CGNS output files and
#   plot starting at the given <start_time>. If no <start_time> is given, use all
#   available data.
#
# USAGE:
#   ./plot_part_position.py <./path/to/sim/output> <start_time>
#
# OUTPUT
#   <./path/to/sim/output>/img/part_pos_time.png

# Imports:
import sys, os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import bluebottle_particle_reader as bbparts
import bluebottle_flow_reader as bbflow

#########

# Parse output directory from commandline
if len(sys.argv) >= 2:    # output directory given
  data_dir = sys.argv[1]

  if len(sys.argv) >= 3:  # start time given
    t_start= sys.argv[2]

else:                     # nothing given
  print("plot_part_position error: Invalid commandline arguments.")
  print("Usage: ")
  print("   ./plot_part_position.py <./path/to/sim/output> <start_time>")
  print(" or")
  print("   ./plot_part_position.py <./path/to/sim/output>")
  sys.exit()

# Init the reader
times = bbparts.init(data_dir)

grid = bbflow.read_flow_extents()

# Get nparts
bbparts.open(times[0])
nparts = bbparts.read_nparts()
bbparts.close()

# Init data arays
x = np.zeros(nparts)
y = np.zeros(nparts)
z = np.zeros(nparts)
u = np.zeros((len(times), nparts))
v = np.zeros((len(times), nparts))
w = np.zeros((len(times), nparts))
t = np.zeros(len(times))
# Loop over time and pull data
for tt,time in enumerate(times):
  bbparts.open(time)

  t[tt] = bbparts.read_time()

  (x,y,z) = bbparts.read_part_position()
  u[tt,:]=x
  v[tt,:]=y
  w[tt,:]=z

  bbparts.close()

  # Plot
  fig = plt.figure(num=1)

  ax1 = fig.add_subplot(111,projection='3d')
  ax1.scatter(x, y, z)
  ax1.set_xlabel("$x$")
  ax1.set_ylabel("$y$")
  ax1.set_zlabel("$z$")
  ax1.set_xlim3d(left=grid[3],right=grid[4])
  ax1.set_ylim3d(bottom=grid[6],top=grid[7])
  ax1.set_zlim3d(bottom=grid[9],top=grid[10])
  plt.title("t=%.2f" %t[tt])

  # Save figure
  img_dir = data_dir + "/img/"
  if not os.path.exists(img_dir):
    os.makedirs(img_dir)
  plt.savefig(img_dir + "part_pos_%.2f.png" %t[tt], bbox_inches='tight', format='png')

  plt.clf()

np.savetxt(img_dir+"part_pos_x.csv",u,delimiter=",")
np.savetxt(img_dir+"part_pos_y.csv",v,delimiter=",")
np.savetxt(img_dir+"part_pos_z.csv",w,delimiter=",")
np.savetxt(img_dir+"t.csv",t,delimiter=",")
