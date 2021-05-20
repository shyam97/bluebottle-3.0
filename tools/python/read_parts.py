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

# bluebottle_particle_reader python module example code

import sys, getopt
import numpy as np
import bluebottle_particle_reader as bbparts

# Parse output directory from commandline
if len(sys.argv) >= 2:    # output directory given
  data_dir = sys.argv[1]

  if len(sys.argv) >= 3:  # start time given
    t_start= sys.argv[2]

else:                     # nothing given
  print("read_parts error: Invalid commandline arguments.")
  print("Usage: ")
  print("   ./read_parts.py <./path/to/sim/output> <start_time>")
  print(" or")
  print("   ./read_parts.py <./path/to/sim/output>")
  sys.exit()

# initialize the reader
times = bbparts.init(data_dir)

# visit all outputted time values
for time in times:
  # open the CGNS file for this particular output time
  bbparts.open(time)

  # read the CGNS file
  t = bbparts.read_time()
  n = bbparts.read_nparts()
  (x,y,z) = bbparts.read_part_position()
  (u,v,w) = bbparts.read_part_velocity()

  print("time = ", time, "t =", t, "n =", n)

  print(u)
  sys.exit()

  # close the CGNS file
  bbparts.close()
