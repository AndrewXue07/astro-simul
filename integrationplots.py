import rebound 
from rebound import hash as h # same thing as h = rebound.hash
import matplotlib.pyplot as plt
import numpy as np

sa = rebound.Simulationarchive("archive.bin")

# print(f"Number of snapshots: {len(sa)}")
# print(f"Time of first and last snapshot: {sa.tmin}, {sa.tmax}")

N_pl = 6 # number of planets + Sun
N_tp = 20 # number of "massless" test particles

# # initial conditions: 
# sim = sa[0] # 0th element in simulation archive list
# rebound.OrbitPlot(sim, unitlabel = "[AU]", color = (N_pl-1)*["black"] + N_tp*["red"], xlim = [-5, 5], ylim = [-5, 5])
# plt.title("Initial Conditions")
# # plt.show()


# # latest conditions/snapshot:
# sim = sa[-1]
# rebound.OrbitPlot(sim, unitlabel = "[AU]", color = (N_pl-1)*["black"] + N_tp*["red"], xlim = [-5, 5], ylim = [-5, 5])
# plt.title("Latest Conditions (right before last asteroid leaves Solar System)")

# sim.status()
# orbits = sim.orbits()
# for orbit in orbits:
#     print(orbit)

# plt.show()



# # another plot: extract information about a single particle
# t = np.zeros(len(sa)) # times
# a = np.zeros(len(sa)) # semi-major axis
# e = np.zeros(len(sa)) # eccentricity


# for i, sim in enumerate(sa):
#     pid = 100
#     t[i] = sim.t/1e6

#     try:
#         a[i] = sim.particles[h(pid)].a
#         e[i] = sim.particles[h(pid)].e
#     except rebound.ParticleNotFound:
#         a = a[:i]
#         e = e[:i]
#         t = t[:i]
#         break

# plt.plot(t, a, label = "semi-major axis")
# q = a*(1-e)
# Q = a*(1+e)
# plt.plot(t, q, label = "q")
# plt.plot(t, Q, label = "Q")

# plt.xlabel("time (Million Years)")
# plt.ylabel("a (AU)")
# plt.title("Particle {0}".format(pid))
# plt.legend()
# plt.yscale("log")

# plt.show()



# # movie of particle's orbit
# pid = 100
# count = 0 # start on first frame
# fig, ax = plt.subplots() 
# op1 = rebound.OrbitPlot(sa[0], orbit_style = "solid", lw = 1, particles = [1, 2, 3, 4, 5, 6])
# op1.particles.set_sizes([0])
# op2 = rebound.OrbitPlot(sa[0], fig = op1.fig, ax = op1.ax, orbit_style = "solid", lw = 1, particles = [h(pid)], color = "red")
# op2.particles.set_sizes([0])


# for i, sim in enumerate(sa): 
#     if i%100: # plot every 100th frame
#         continue 
#     try: 
#         op1 = rebound.OrbitPlot(sim, fig = fig, ax = ax, orbit_style = "solid", lw = 1, particles = [1, 2, 3, 4, 5, 6]) # first particle
#         op1.particles.set_sizes([0]) # don't plot the particles (only plot orbits)
#         op2 = rebound.OrbitPlot(sim, fig = fig, ax = ax, orbit_style = "solid", lw = 1, particles = [h(pid)], color = "red") # graph on same figre/axes as op1
#         op2.particles.set_sizes([0])

#         ax.set_aspect("equal")
#         ax.set_xlim(-8, 8)
#         ax.set_ylim(-8, 8)
#         ax.set_xlabel("x [AU]")
#         ax.set_ylabel("y [AU]")
#         ax.set_title("Particle {0}".format(pid))

#         ax.text(-4, -7, "t = {:4.1f} Myr".format(sim.t/1e6))
#         fig.savefig("frame_{:04d}.png".format(count))

#         count += 1 # next frame
#         ax.clear() # clear axis to reuse it in the next frame

#     except rebound.ParticleNotFound:
#         break



# # "residence plot" showing change in semi-major axis and eccentricity of all clones over time
# a_all = []
# e_all = []

# for i, sim in enumerate(sa): # every frame
#     for pid in range(100, 100+N_tp):
#         try:
#             a_all.append(sim.particles[h(pid)].a)
#             e_all.append(sim.particles[h(pid)].e)
#         except rebound.ParticleNotFound:
#             continue 

# amin, amax = 0, 12
# emin, emax = 0, 1
# h2d, xedge, yedge, im = plt.hist2d(a_all, e_all, range = [[amin, amax], [emin, emax]], bins = (30, 30))
# plt.xlabel("a (AU)")
# plt.ylabel("e")

# plt.colorbar() # units = ADU counts 

# a_init = 3.260026162454334
# e_init = 0.6157226608577462
# plt.plot(a_init, e_init, "*", color = "r")

# # plt.show()

# # where particles are at last timestep:
# sim = sa[-1]
# a_tp = []
# e_tp = []

# for pid in range(100, 100+N_tp):
#     try:
#         a_tp.append(sim.particles[h(pid)].a)
#         e_tp.append(sim.particles[h(pid)].e)
#     except:
#         continue

# plt.plot(a_tp, e_tp, "co")

# a_pl = []
# e_pl = []
# for id in range(1, N_pl):
#     a_pl.append(sim.particles[id].a)
#     e_pl.append(sim.particles[id].e)
# plt.plot(a_pl, e_pl, "yo")



# fig, ax = plt.subplots()
# cbar = fig.colorbar(im, ax = ax)
# cax = cbar.ax

# count = 0
# for i, sim in enumerate(sa):
#     if i%100: # only do every 100 snapshots
#         continue

# amin, amax = 0, 12
# emin, emax = 0, 1
# h2d, xedge, yedge, im = plt.hist2d(a_all, e_all, range = [[amin, amax], [emin, emax]], bins = (30, 30))
# plt.xlabel("a (AU)")
# plt.ylabel("e")

    # fig.savefig("frame_{:04d}.png".format(count))
    # ax.clear() # clear axis to reuse it 
    # count += 1



plt.show()