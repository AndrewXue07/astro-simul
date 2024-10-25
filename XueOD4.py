"""Determine the six orbital elements (a, e, i, Omega, omega, M)"""
# given: position vector; velocity vector 
from odlib import *
import numpy as np
import math

# use orbital elements at a given time to produce RA/dec for any time
# earth-sun vector on July 13, 2018, 00:00:00 UT: [-3.530809229063387E-01, 8.746426747111206E-01, 3.791364871830791E-01]

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

def getOmega(h, i): # longitude of ascending node
    cosVal = -h[1] / (magnitude(h) * math.sin(math.radians(i)))
    sinVal = h[0] / (magnitude(h) * math.sin(math.radians(i)))

    return quadrantCheck(cosVal, sinVal)

def getomega(a, e, i, Omega, h, position, velocity): # argument of perihelion
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
        return U - V, V
    else: 
        return 360 + (U - V), V

def getMeanAnomaly(v, a, e, position): # output mean anomaly and time of last perihelion passage (Julian days) and Julian date
    cosE = (1-magnitude(position) / a) / e
    E = 0 
    # print("v", v)

    if (v < 180 and v >= 0):
        E = math.acos(cosE)
    elif (v >= 180 and v <= 360):
        E = 2*math.pi - math.acos(cosE)
    else:
        print("Look at E value (negative or above 2pi)")
    
    return math.degrees(E - e*math.sin(E))

def getNewMeanAnomaly(lastTime, targetTime, M_0, a):
    period = math.sqrt((4*(math.pi **2)) * a **3) # units = Gaussian days
    # print("period", period)
    M = M_0 + 360 * ((targetTime - lastTime) / 58.13244087) / period

    return M % 360
    # if M > 360: # or M < -2*math.pi:
    #     # print("M", M)
    #     return M % 360
    # else:
    #     # print("M", M)

    #     return M

def getEclipticAnomaly(M, e):
    E = M # inital guess for E
    # print(f"M, {math.radians(M)}")
    M_guess = E - e*math.sin(E)

    while (abs(M_guess - M) > 10**-12):
        M_guess = E - e*math.sin(E)
        E = E - (M - (E - e*math.sin(E))) / (e*math.cos(E) - 1)
    
    # print(f"E: {E}")
    
    return E

def getPosVector(a, e, E):
    # print("a", a)
    # print("e", e)
    # print("E", E)
    return np.array([a*math.cos(E) - a*e, a*(math.sqrt(1-e**2))*math.sin(E), 0])



def orbital_to_ecliptic(position, omega, i, Omega):

    # print(f"omega {omega}, i {i}, Omega {Omega}")
    return np.array([[math.cos(math.radians(Omega)), -math.sin(math.radians(Omega)), 0], 
                     [math.sin(math.radians(Omega)), math.cos(math.radians(Omega)), 0], 
                     [0, 0, 1]]) @ np.array([
                     [1, 0, 0], 
                     [0, math.cos(math.radians(i)), -math.sin(math.radians(i))], 
                     [0, math.sin(math.radians(i)), math.cos(math.radians(i))]]) @ np.array([
                    [math.cos(math.radians(omega)), -math.sin(math.radians(omega)), 0], 
                     [math.sin(math.radians(omega)), math.cos(math.radians(omega)), 0], 
                     [0, 0, 1]]) @ position
# orbital to ecliptic: 
# 1: 
def rotate_w(position, omega): # rotate argument of perihelion to align with ascending node
    # omega = -omega
    return np.array([[math.cos(math.radians(omega)), -math.sin(math.radians(omega)), 0], 
                     [math.sin(math.radians(omega)), math.cos(math.radians(omega)), 0], 
                     [0, 0, 1]]) @ position

#2: 
def rotate_i(position, i): # rotate inclination to math the North Ecliptic Pole
    # i = -i
    return np.array([[1, 0, 0], 
                     [0, math.cos(math.radians(i)), -math.sin(math.radians(i))], 
                     [0, math.sin(math.radians(i)), math.cos(math.radians(i))]]) @ position

#3: 
def rotate_Omega(position, Omega):
    # Omega = -Omega
    return np.array([[math.cos(math.radians(Omega)), -math.sin(math.radians(Omega)), 0], [math.sin(math.radians(Omega)), math.cos(math.radians(Omega)), 0], [0, 0, 1]]) @ position

def ecliptical_to_equatorial(position, tilt):

    
    # tilt = -tilt
    return np.array([[1, 0, 0], 
                     [0, math.cos(math.radians(tilt)), -math.sin(math.radians(tilt))], 
                     [0, math.sin(math.radians(tilt)), math.cos(math.radians(tilt))]]) @ position

def equatorialRect_to_equatorialPolar(earth_dist_to_sun, position):
    # print("earth sun vector", earth_dist_to_sun)
    # print(f"position: {position}")

    rho = earth_dist_to_sun + position
    print("rho", rho)
    # print("rho[2] -> z component", rho[2])

    # print("position", position)
    dec = math.degrees(math.asin(rho[2] / magnitude(rho)))
    # print("dec", dec)
    # print(f"dec: {dec}")
    sinRA = rho[1] / (magnitude(rho) * math.cos(math.radians(dec)))
    cosRA = rho[0] / (magnitude(rho) * math.cos(math.radians(dec)))

    # print(f"sinRA {sinRA} and cosRA {cosRA}")
    RA = quadrantCheck(cosRA, sinRA)

    return RA, dec

def main():
    input = open("XueInput2.txt")
    lineArr = input.readlines()

    # position = np.array(lineArr[0].strip().split(" ")).astype(float)
    # # print(position)
    # velocity = (1/0.0172020989484) * np.array(lineArr[1].strip().split(" ")).astype(float)
    # h = getAngMomentum(position, velocity)

    lastTime = float(lineArr[6]) # Julian date
    # last_time_of_perihelion = float(lineArr[3])

    # determine orbital elements
    a = float(lineArr[0])
    e = float(lineArr[1])
    i = float(lineArr[2])
 
    Omega = float(lineArr[3])
    omega = float(lineArr[4])
    M = float(lineArr[5])

    targetTime = float(lineArr[9])
    M_current = getNewMeanAnomaly(lastTime, targetTime, M, a)

    earth_sun_vector = np.array(lineArr[8].strip().split(" ")).astype(float)
    # correct_RA_and_dec = np.array(lineArr[5].strip().split(" ")) # hours minutes seconds

    # print("new mean anomaly", math.degrees(M_current))
    # print(f"M_current: {M_current}")
    # print(f"i: {i}, omega {omega}, Omega {Omega}")
    E = getEclipticAnomaly(M_current, e)
    new_pos = getPosVector(a, e, E)
    # print("initial position", new_pos)
    new_pos = orbital_to_ecliptic(new_pos, omega, i, Omega)
    # print("new_pos before ecliptical to equatorial", new_pos)
    new_pos = ecliptical_to_equatorial(new_pos, 23.4354844788)
    RA, dec = equatorialRect_to_equatorialPolar(earth_sun_vector, new_pos)

    print(f"RA: Expected = 17:42:21.20, got {RAdecimalToHMS(RA)}, % error = {100*(RA - HMStoDeg("17:42:21.20", False)) / HMStoDeg("17:42:21.20", False)}")
    print(f"dec: Expected = 31:52:28.1, got {DECdecimalToDMS(dec)}, % error = {100*(dec - DMStoDeg("31:52:28.1", False)) / DMStoDeg("31:52:28.1", False)}")


    input.close()


if __name__ == "__main__":
    main()

# output angular momentum vector rounded to 6 decimal places 
# units: AU (distance), Gaussian days (time)

# write "runs fine" when submitting 
# create zip files with all necessary files 