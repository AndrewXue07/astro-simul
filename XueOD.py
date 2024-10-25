import numpy as np
from math import *
from odlib import *
from XueOD5 import mog
from XueOD4 import *
import matplotlib.pyplot as plt

k = 0.0172020989484 # Gaussian gravitational constant; 1 day = 0.0172... Gd
cAU = 173.144643267 # speed of light in AU/mean solar day
eps = math.radians(23.4384668053) # Earth's obliquity

def retrieveData():
    # file = open("XueInput.txt")
    # file = open("testinput_7_23_2024.txt")
    # file = open("gitelsonInputTestCase2.txt")
    file = open("finalorbitalelements_input.txt")

    times = []
    RAs = []
    decs = []
    earth_sun_vectors = []

    for line in file.readlines():
        tempArr = line.split(" ")

        times.append(float(tempArr[0])) # time of middle of observation in Julian days
        RAs.append(HMStoDeg(tempArr[1], False))
        decs.append(DMStoDeg(tempArr[2], False))
        earth_sun_vectors.append((float(tempArr[3]), float(tempArr[4]), float(tempArr[5])))
    
    file.close()
    return times, RAs, decs, earth_sun_vectors

def calculateRhoUnitVectors(RA, dec): # takes in decimal RA/dec in degrees; calculates p_x hat, p_y hat, p_z hat
    rho_x = cos(radians(RA)) * cos(radians(dec))
    rho_y = sin(radians(RA)) * cos(radians(dec))
    rho_z = sin(radians(dec))
    
    return rho_x, rho_y, rho_z

def calcPosVector(e_s_arr_0, rho_unit_arr1, rho_unit_arr2, rho_unit_arr3, c1, c2, c3): # p1, p2, p3
    e_s_arr = np.array(e_s_arr_0)

    D0 = np.dot(rho_unit_arr1, np.cross(rho_unit_arr2, rho_unit_arr3))

    D1_arr = [0, 0, 0]
    D2_arr = [0, 0, 0]
    D3_arr = [0, 0, 0]

    # print("earth sun 3", e_s_arr[2])
    for i in range(3):
        D1_arr[i] = np.dot(np.cross(e_s_arr[i], rho_unit_arr2), rho_unit_arr3)
        D2_arr[i] = np.dot(np.cross(rho_unit_arr1, e_s_arr[i]), rho_unit_arr3)
        D3_arr[i] = np.dot(rho_unit_arr1, np.cross(rho_unit_arr2, e_s_arr[i]))

        # print(f"Iteration {i}, D1: {D1_arr}, D2: {D2_arr}, D3: {D3_arr}")

    # print(f"D0: {D0}, D1: {D1_arr}, D2: {D2_arr}, D3: {D3_arr}")
    # print(f"c1: {c1}, c2 {c2}, c3 {c3}")

    rho_1 = (c1*D1_arr[0] + c2*D1_arr[1] + c3*D1_arr[2]) / (c1*D0)
    rho_2 = (c1*D2_arr[0] + c2*D2_arr[1] + c3*D2_arr[2]) / (c2*D0)
    rho_3 = (c1*D3_arr[0] + c2*D3_arr[1] + c3*D3_arr[2]) / (c3*D0)

    # print(f"rho1: {rho_1}, rho2: {rho_2}, rho3: {rho_3}")

    r1 = float(rho_1) * rho_unit_arr1 - e_s_arr[0]
    r2 = float(rho_2) * rho_unit_arr2 - e_s_arr[1]
    r3 = float(rho_3) * rho_unit_arr3 - e_s_arr[2]

    # print(f"r1: {r1}, r2: {r2}, r3: {r3}")
    return r1, r2, r3, rho_1, rho_2, rho_3

def calcVelocity(r1, r2, r3, tau0, tau1, tau3):
    r12_dot = (r2 - r1) / (-tau1) # calculate average velocities and assume them as instantenous velocities
    r23_dot = (r3 - r2) / tau3
    # print(f"r12_dot {r12_dot}, r23_dot {r23_dot}")

    r2_dot = (tau3 / tau0)*r12_dot - (tau1 / tau0)*r23_dot
    # print(r2_dot)
    return r2_dot

def lightTravelCorrection(timeArr, rhoArr):
    newTimes = [0, 0, 0]
    for i in range(3):
        newTimes[i] = timeArr[i] - (k*rhoArr[i]/(cAU))
    
    return newTimes

def get_f_and_g(tau1, tau3, r2, r2_dot):
    f1, f3, g1, g3 = mog(tau1, tau3, r2, r2_dot, 4)

    return f1, f3, g1, g3

def get_c_and_d(f1, f3, g1, g3):
    c1 = g3 / (f1*g3 - g1*f3)
    c3 = -g1 / (f1*g3 - g1*f3)

    d1 = -f3 / (f1*g3 - f3*g1)
    d3 = f1 / (f1*g3 - f3*g1)

    return c1, c3, d1, d3

def equatorialToEcliptical(r2):
    return np.linalg.inv([[1, 0, 0], [0, cos(eps), -sin(eps)], [0, sin(eps), cos(eps)]]) @ r2


def main(useMonteCarlo, RA_list, dec_list, todaysDate): 
    
    times, RAs, decs, earth_sun_vectors = retrieveData()
    if useMonteCarlo: # if we would like to manually input RAs and decs instead (not from a separate .txt file, but within the .py file)
        RAs, decs = RA_list, dec_list

    rho_x1, rho_y1, rho_z1 = calculateRhoUnitVectors(RAs[0], decs[0]) # first RA, dec pair
    rho_x2, rho_y2, rho_z2 = calculateRhoUnitVectors(RAs[1], decs[1]) # second RA, dec pair
    rho_x3, rho_y3, rho_z3 = calculateRhoUnitVectors(RAs[2], decs[2]) # third RA, dec pair

    # print("earth sun vectors", earth_sun_vectors)
    tau0 = k*(times[2] - times[0]) # t3 - t1
    tau1 = k*(times[0] - times[1]) # t1 - t2
    tau3 = k*(times[2] - times[1]) # t3 - t2
    c1 = tau3 / tau0 # initial c1, c2, c3 guesses
    c2 = -1
    c3 = -tau1 / tau0 
    r1, r2, r3, rho_1scalar, rho_2scalar, rho_3scalar = calcPosVector(earth_sun_vectors, np.array([rho_x1, rho_y1, rho_z1]), np.array([rho_x2, rho_y2, rho_z2]), np.array([rho_x3, rho_y3, rho_z3]), c1, c2, c3)
    r2_dot = calcVelocity(r1, r2, r3, tau0, tau1, tau3)


    rho2 = rho_2scalar * np.array([rho_x2, rho_y2, rho_z2])
    # print("first rho2:", rho2)
    rho2_current = [-1000000, 1000, 1000]
    rho2_arr = []
    while(abs(magnitude(rho2) - magnitude(rho2_current)) > 10**-10):
        rho2_arr.append(magnitude(rho2))
        # print("iteration")
        # print(f"Current difference: {abs(magnitude(rho2) - magnitude(rho2_current))}")
        # print("iteration", r2, r2_dot)
        # print(f"rho2: {rho2}")
        newTimes = lightTravelCorrection(times, [rho_1scalar, rho_2scalar, rho_3scalar])
        # print(newTimes)
        f1, f3, g1, g3 = get_f_and_g(tau1, tau3, r2, r2_dot)
        # print(f"f1: {f1}, f3: {f3}, g1: {g1}, g3: {g3}")
        c1, c3, d1, d3 = get_c_and_d(f1, f3, g1, g3)
        # print(c1, c3)

        rho2 = rho_2scalar * np.array([rho_x2, rho_y2, rho_z2])
        rho_x1, rho_y1, rho_z1 = calculateRhoUnitVectors(RAs[0], decs[0]) # first RA, dec pair
        rho_x2, rho_y2, rho_z2 = calculateRhoUnitVectors(RAs[1], decs[1]) # second RA, dec pair
        rho_x3, rho_y3, rho_z3 = calculateRhoUnitVectors(RAs[2], decs[2]) # third RA, dec pair

        tau0 = k*(newTimes[2] - newTimes[0]) # t3 - t1
        tau1 = k*(newTimes[0] - newTimes[1]) # t1 - t2
        tau3 = k*(newTimes[2] - newTimes[1]) # t3 - t2
        # print("earth sun vectors", earth_sun_vectors)
        r1, r2, r3, rho_1scalar, rho_2scalar, rho_3scalar = calcPosVector(earth_sun_vectors, np.array([rho_x1, rho_y1, rho_z1]), np.array([rho_x2, rho_y2, rho_z2]), np.array([rho_x3, rho_y3, rho_z3]), c1, c2, c3)
        r2_dot = d1*r1 + d3*r3

        newTimes = lightTravelCorrection(newTimes, [rho_1scalar, rho_2scalar, rho_3scalar])
        # print(newTimes)
        f1, f3, g1, g3 = get_f_and_g(tau1, tau3, r2, r2_dot)
        # print(f"f1: {f1}, f3: {f3}, g1: {g1}, g3: {g3}")
        c1, c3, d1, d3 = get_c_and_d(f1, f3, g1, g3)
        c2 = -1
        # print(c1, c3)

        rho2_current = rho_2scalar * np.array([rho_x2, rho_y2, rho_z2])
        # print("rho2_current:", rho2_current)
        # print("rho2", rho2)
        # print("iteration")

    # print(f"r2: {r2}, r2_dot: {r2_dot}")
    pos_ecliptical = equatorialToEcliptical(r2)
    velocity_ecliptical = equatorialToEcliptical(r2_dot)

    h = getAngMomentum(pos_ecliptical, velocity_ecliptical)
    a = getSemiMajorAxis(pos_ecliptical, velocity_ecliptical)
    e = getEccentricity(h, pos_ecliptical, velocity_ecliptical, a)
    i = getInclination(h)
    Omega = getOmega(h, i) # longitude of ascending node
    omega, V = getomega(a, e, i, Omega, h, pos_ecliptical, velocity_ecliptical) # argument of perihelion
    # print(f"cosE: {(1-magnitude(pos_ecliptical) / a) / e}")
    meanAnomaly = getMeanAnomaly(V, a, e, pos_ecliptical)
    meanAnomaly = getNewMeanAnomaly(newTimes[1], todaysDate, meanAnomaly, a)

    # print(f"difference in days: {times[1]}, and difference is {times[1] - todaysDate}. In gaussian days it's {(times[1] - todaysDate) * k}")
    # print("unit vectors", rho_x2, rho_y2, rho_z2)
    # print("rho2 scalar", rho_2scalar)

    if not useMonteCarlo:
        print(f"Asteroid's position: {pos_ecliptical} AU (r2)")
        print(f"Asteroid's Velocity at t2: {(k)*velocity_ecliptical} AU/day (r2_dot)")
        # print(r2_dot)
        # print(magnitude(rho2_current))
        print(f"Range Vector at t2: {rho2_current} (rho2)")
        
        print(f"Semi-Major Axis: {a} AU")
        print(f"Eccentricity: {e}")
        print(f"Inclination: {i} Degrees")
        print(f"Longitude of Ascending Node: {Omega} Degrees")
        print(f"Argument of Perihelion: {omega} Degrees")
        # print("V", V)
        print(f"Mean Anomaly: {meanAnomaly} Degrees")

        plt.plot(np.arange(0, len(rho2_arr)), rho2_arr)
        plt.title("Magnitude of ρ_2 over time with the method of Gauss")
        plt.xlabel("Iteration #")
        plt.ylabel("Magnitude of ρ_2")
        plt.show()
    return a, e, i, Omega, omega, meanAnomaly
 
if __name__ == "__main__":
    main(False, 0, 0, 2460518.750000)