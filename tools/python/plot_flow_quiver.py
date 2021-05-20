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

# bluebottle_reader python module example code

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import bluebottle_flow_reader as bbflow

#########

# Parse output directory from commandline
if len(sys.argv) >= 2:    # output directory given
  data_dir = sys.argv[1]

  if len(sys.argv) >= 3:  # start time given
    t_start= sys.argv[2]

else:                     # nothing given
  print("plot_flow_streamlines error: Invalid commandline arguments.")
  print("Usage: ")
  print("   ./plot_flow_streamlines <./path/to/sim/output> <start_time>")
  print(" or")
  print("   ./plot_flow_streamlines.py <./path/to/sim/output>")
  sys.exit()

# initialize the reader
times = bbflow.init(data_dir)
grid = bbflow.read_flow_extents()

# Pull grid positions
(x,y,z) = bbflow.read_flow_position()

# visit all outputted time values
for time in times:
    # open the CGNS file for this particular output time
    print(time)
    bbflow.open(time)

    # read the CGNS file
    t = bbflow.read_time()
    (u,v,w) = bbflow.read_flow_velocity()

    fig = plt.figure(num=1)
    ax1 = fig.add_subplot(111,projection='3d')
    ax1.quiver(x, y, z, u, v, w, normalize=True)
    ax1.set_xlabel("$x$")
    ax1.set_ylabel("$y$")
    ax1.set_zlabel("$z$")
    ax1.set_xlim3d(left=grid[3],right=grid[4])
    ax1.set_ylim3d(bottom=grid[6],top=grid[7])
    ax1.set_zlim3d(bottom=grid[9],top=grid[10])
    plt.title("t=%.2f" %t[tt])

    # close the CGNS file
    bbflow.close()
