import rebound
import matplotlib.pyplot as plt
import numpy as np
from math import *
import os

sim = rebound.Simulation()
sim.integrator = "mercurius"
sim.units = ("yr", "AU", "Msun")
date = "JD2460499.898553"
discard_file_name = "discards.txt"

sim.add("Sun", date = date, hash = 0)
h = rebound.hash 

try:
    os.remove(discard_file_name) # clear the discard file if it exists
except:
    pass 

discard_file_name = "discards.txt"
def collision_discard_log(sim_pointer, collision, discard_file_name = discard_file_name): # sim_pointer = pointer to collision
    sim = sim_pointer.contents # simulation object from pointer
    id_p1 = sim.particles[collision.p1].hash.value
    id_p2 = sim.particles[collision.p2].hash.value

    discard_file = open(discard_file_name, "a")

    if id_p1 > id_p2:
        print(f"Particle{id_p1} collidied with {id_p2} at {sim.t} years")
        print(f"Particle{id_p1} collidied with {id_p2} at {sim.t} years", file = discard_file)
        print(f"Removing particle{id_p1}")
        ToRemove = 1 
    else: 
        print(f"Particle{id_p2} collidied with {id_p1} at {sim.t} years")
        print(f"Particle{id_p2} collidied with {id_p1} at {sim.t} years", file = discard_file)
        print(f"Removing particle{id_p2}")
        ToRemove = 2

    discard_file.close()
    return ToRemove # id of particle to be removed 

# print(sim.particles[h(0)].m)
# print("gravitational constant * solar mass:", sim.G*sim.particles[h(0)].m) # 0 = list index but h(0) is the hash value
sim.particles[h(0)].r = 0.00465047

sim.add("Venus", hash = 1)
sim.particles[h(1)].r = 4.04537843*10**-5
sim.add("Earth", hash = 2)
sim.particles[h(2)].r = 4.26352e-5
sim.add("Mars", hash = 3)
sim.particles[h(3)].r = 2.26608e-5
sim.add("Jupiter", hash = 4)
sim.particles[h(4)].r = 0.00046732617
sim.add("Saturn", hash = 5)
sim.particles[h(5)].r = 0.000389256877

N_pl = 6 # number of planets + Sun
N_tp = 20 # number of "massless" test particles
sim.N_active = N_pl # number of active "massive" objects (for example neglect asteroid's mass)
sim.add(a = 3.260026162454334, e = 0.6157226608577462, inc = radians(9.99806183380943), Omega = radians(148.00565303763332), omega = radians(142.3811217193129), M = radians(1.975142962970665), hash = 100)
sim.particles[h(100)].r = 6.6846e-9

x = sim.particles[h(100)].x
y = sim.particles[h(100)].y
z = sim.particles[h(100)].z
vx = sim.particles[h(100)].vx
vy = sim.particles[h(100)].vy
vz = sim.particles[h(100)].vz

for i in range(1, N_tp):
    sim.add(x = x + np.random.uniform(-10**-4, 10**-4), y = y + np.random.uniform(-10**-4, 10**-4), 
            z = z + np.random.uniform(-10**-4, 10**-4), vx = vx + np.random.uniform(-10**-4, 10**-4), 
            vy = vy + np.random.uniform(-10**-4, 10**-4), vz = vz + np.random.uniform(-10**-4, 10**-4), 
            r = 6.6846e-9, hash = 100 + i)

# for p in sim.particles:
#     print(p.hash.value)
# print(sim.N)

sim.exit_max_distance = 1000 # after 1000 AU, asteroid = ejected from Solar System
sim.collision = "direct" # radii of two particles overlap
sim.collision_resolve = collision_discard_log # resolving collisions
sim.move_to_com() # move center of image to center of mass (near Sun barycenter)
tend = 50e6 # time end -> run the simulation for 50 million years
tout = 1000 # time output -> output every 1000 years 
sim.dt = sim.particles[h(1)].P / 25 # use 1/25 of the smallest period (Venus) as the time step 

archive = "archive.bin"
sim.save_to_file(archive, interval = tout, delete_file = True) # save snapshots of every particle's position/velocity to archive file
times = np.arange(0, tend, tout) # array of all the times for outputting data
Nsteps = len(times)

for i in range(Nsteps): # number of times to output stuff
    try: 
        sim.integrate(times[i], exact_finish_time = 0) # exact_finish_time = 0 means that it won't always be exactly 1000 but rather something more convenient depending on situation to be faster
    except rebound.Escape as error: # unless if a particle escapes (sim.exit_max_distance = 1000)
        
        for j in range(sim.N): # for each of the particles
            p = sim.particles[j]
            d2 = p.x**2 + p.y**2 + p.z**2 # distance to Sun squared

            if d2 > (sim.exit_max_distance**2): 
                index = j # save index but don't remove it because otherwise it messes up for loop's indexing
        
        pid = sim.particles[index].hash.value 
        print(f"Particle {pid} was too far from the Sun at {sim.t} years")

        discard_file = open(discard_file_name, "a")
        print(f"Particle {pid} was too far from the Sun at {sim.t} years", file = discard_file)
        discard_file.close()

        sim.remove(index = index) # finally remove the ejected asteroid's index
    
    print(f"Time {sim.t / 1e6} Myr-- Fraction Done {sim.t / tend}-- # of clones {sim.N - N_pl}")

    if sim.N <= N_pl: # if number of particles is less than or equal to # of planets + Sun (N_pl), end integration
        print("No more test particles, ending simulation")
        break



rebound.OrbitPlot(sim, unitlabel = "[AU]", color = (N_pl-1)*["black"] + (N_tp)*["red"]) # only top view
rebound.OrbitPlotSet(sim, unitlabel = "[AU]", xlim = [-5, 5], ylim = [-5, 5], color = (N_pl-1)*["black"] + (N_tp)*["red"]) # x, y, and z views
plt.show()