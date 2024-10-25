import rebound
from rebound import hash as h
import numpy as np
import matplotlib.pyplot as plt

# "getting a closer look"
# ERROR: they never collide (asteroids escape instead)

N_pl = 6
sa = rebound.Simulationarchive("archive.bin")
sim = sa[186] # 186,000 years

archive = "archive_new.bin"

tstart = sim.t # current time
tout = 1.0 # new output resolution (more detail, starting from the start time tstart)
tend = 1000 + tstart # end at the next timestep

sim.save_to_file(archive, interval = tout, delete_file = True)
sim.collision_resolve = "halt"
times = np.arange(tstart, tend, tout)
Nsteps = len(times)

for i in range(Nsteps):
    sim.integrate(times[i], exact_finish_time = 0)
    print(f"Time {sim.t/1e6} Myr --Fraction Done {(sim.t - tstart) / (tend - tstart)} -- # of clones {sim.N - N_pl}")

sa_new = rebound.Simulationarchive("archive_new.bin")

pid = 110
d_Jupiter = np.zeros(len(sa_new))
d_Saturn = np.zeros(len(sa_new))
for i, sim in enumerate(sa_new):
    d_Jupiter[i] = ((sim.particles[h(pid)].x - sim.particles[h(4)].x)**2 + (sim.particles[h(pid)].y - sim.particles[h(4)].y)**2 + (sim.particles[h(pid)].z - sim.particles[h(4)].z)**2)**0.5
    d_Saturn[i] = ((sim.particles[h(pid)].x - sim.particles[h(4)].x)**2 + (sim.particles[h(pid)].y - sim.particles[h(4)].y)**2 + (sim.particles[h(pid)].z - sim.particles[h(4)].z)**2)**0.5

plt.plot(sim.t, d_Jupiter, label = "Jupiter")
plt.plot(sim.t, d_Saturn, label = "Saturn")
plt.xlabel("t Myr")
plt.ylabel("Distance to Planet (AU)")
plt.legend()

plt.show()