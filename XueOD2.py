# runs fine

"""Determine the six orbital elements (a, e, i, Omega, omega, M)"""
# given: position vector; velocity vector 
from odlib import *
import numpy as np
import math

mew = 1

def getAngMomentum(position, velocity): # angular momentum = r x p
    return np.cross(position, velocity)

def getSemiMajorAxis(position, velocity): # calculate a
    return ((2/magnitude(position)) - ((magnitude(velocity)**2) / mew)) **(-1)

def getEccentricity(h, position, velocity, a):
    numerator = np.cross(position, velocity)
    # print(numerator)
    return math.sqrt(1 - ((numerator[0]) **2 + (numerator[1]) **2 + (numerator[2]) **2)/(mew * a))

def getInclination(h):
    return math.degrees(math.acos((h[2]) / (math.sqrt((h[0]) **2 + (h[1]) **2 + (h[2]) **2))))

def getOmega(h, i):
    cosVal = -h[1] / (magnitude(h) * math.sin(math.radians(i)))
    sinVal = h[0] / (magnitude(h) * math.sin(math.radians(i)))

    return quadrantCheck(cosVal, sinVal)

def getomega(a, e, i, Omega, h, position, velocity): 
    cosV = (1/e) * ((a*(1-e**2)) / magnitude(position) - 1)
    sinV = (a / magnitude(h)) * ((1-e**2) / e) * ((np.dot(position, velocity)) / magnitude(position))
    # print(f"cosV: {cosV} and sinV: {sinV}")

    V = quadrantCheck(cosV, sinV)

    cosU = (position[0] * math.cos(math.radians(Omega)) + position[1] * math.sin(math.radians(Omega))) / magnitude(position)
    sinU = (position[2]) / (magnitude(position) * math.sin(math.radians(i)))
    # print(f"cosU: {cosU} and sinU: {sinU}")

    U = quadrantCheck(cosU, sinU)

    # print("U", U)
    # print("V", V)
    if U - V > 0:
        return U - V
    else: 
        return 360 + (U - V)


def main():
    input = open("Xue_Input.txt")
    lineArr = input.readlines()

    position = np.array(lineArr[0].strip().split(" ")).astype(float)
    # print(position)
    velocity = (1/0.0172020989484) * np.array(lineArr[1].strip().split(" ")).astype(float)
    
    h = getAngMomentum(position, velocity)


    # determine orbital elements
    a = getSemiMajorAxis(position, velocity)
    print(f"Semi-Major Axis: \n Expected = 1.05671892483881, Calculated = {a}, % error = {100*(abs(a-1.05671892483881))/1.05671892483881}%")

    e = getEccentricity(h, position, velocity, a) 
    print(f"Eccentricity: \n Expected = 0.3442798363212599, Calculated = {e}, % error = {100*(abs(e-.3442798363212599))/.3442798363212599}%")

    i = getInclination(h) 
    print(f"Inclination: \n Expected = 25.15375982051783, Calculated = {i}, % error = {100*(abs(i-25.15375982051783))/25.15375982051783}%")

    Omega = getOmega(h, i) 
    print(f"Longitude of Ascending Node: \n Expected = 236.2502480179119, Calculated = {Omega}, % error = {100*(abs(Omega - 236.2502480179119))/236.2502480179119}%")

    omega = getomega(a, e, i, Omega, h, position, velocity)
    print(f"Argument of Perihelion: \n Expected = 255.5316184725226, Calculated = {omega}, % error = {100*(abs(omega - 255.5316184725226))/255.5316184725226}%")
    print("Units = AU or degrees")

    input.close()


if __name__ == "__main__":
    main()

# output angular momentum vector rounded to 6 decimal places 
# units: AU (distance), Gaussian days (time)

# write "runs fine" when submitting 
# create zip files with all necessary files 