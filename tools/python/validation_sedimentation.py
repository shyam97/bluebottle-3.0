import sys, os, getopt
import numpy as np
import matplotlib.pyplot as plt
import bluebottle_flow_reader as bbflow
import bluebottle_particle_reader as bbparts

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

times = bbparts.init(data_dir)
timesf = bbflow.init(data_dir)

bbflow.open(timesf[0])
rho,nu = bbflow.read_flow_params()
bbflow.close()

bbparts.open(times[0])
nparts = bbparts.read_nparts()
bbparts.close()

vel = np.zeros((len(times), nparts))
zvel = np.zeros((len(times), nparts))
t = np.zeros(len(times))
Rep = np.zeros((len(times), nparts))
r = np.zeros(nparts)

for tt,time in enumerate(times):
  bbparts.open(time)

  t[tt] = bbparts.read_time()
  r = bbparts.read_part_radius()

  (u,v,w) = bbparts.read_part_velocity()

  zvel[tt,:] = w

  vel[tt,:] = ((u**2) + (v**2) + (w**2))**0.5

  Rep[tt,:] = vel[tt,:]*r/nu

  bbparts.close()

# Plot
fig = plt.figure(figsize=(6,6),dpi=120)

ax1 = fig.add_subplot(111)
plt.plot(t, zvel)
plt.xlabel("$t$")
plt.ylabel("$w$")
plt.ylim([-0.14,0])
plt.xlim([0,4.5])

# ax1 = fig.add_subplot(212)
# plt.plot(t, vel)
# plt.xlabel("$t$")
# plt.ylabel("$u$")

# Save figure
img_dir = data_dir + "/img/"
if not os.path.exists(img_dir):
  os.makedirs(img_dir)
plt.savefig(img_dir + "vel.png", bbox_inches='tight', format='png')

np.savetxt(img_dir+'t.csv',t,delimiter=',')
np.savetxt(img_dir+'Re.csv',Rep,delimiter=',')
np.savetxt(img_dir+'zvel.csv',zvel,delimiter=',')
