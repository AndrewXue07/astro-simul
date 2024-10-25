import math
import numpy as np

def HMStoDeg(input, giveRad): # intake RA in hour/min/sec and output decimal degrees
    value = input.split(":") # list of hours, minutes, seconds
    if (giveRad):
        return math.radians((float(value[0]) * 15) + 15*(float(value[1]) / 60) + 15*(float(value[2]) / 3600))

    return (float(value[0]) * 15) + 15*(float(value[1]) / 60) + 15*(float(value[2]) / 3600)

def DMStoDeg(input, giveRad): # input declination in degree/arcmin/arcsec and output decimal degrees
    value = input.split(":") # list of degrees, minutes, seconds
    degree = float(value[0])
    arcmin = float(value[1])
    arcsec = float(value[2])

    if (giveRad):
        if (degree < 0):
            return math.radians(math.copysign((abs(degree) + (arcmin / 60) + (arcsec / 3600)),-1))

        return math.radians(abs(degree) + (arcmin / 60) + (arcsec / 3600))

    
    if (degree < 0):
        return math.copysign((abs(degree) + (arcmin / 60) + (arcsec / 3600)),-1)

    return abs(degree) + (arcmin / 60) + (arcsec / 3600)

def quadrantCheck(cosVal, sinVal): # take in degrees and return the version of it that matches a certain quadrant
    if (sinVal == 1):
        return 90
    elif (sinVal == -1):
        return 270
    elif (cosVal == 1):
        return 0
    elif (cosVal == -1):
        return 180
    
    if sinVal > 0:
        return math.degrees(math.acos(cosVal))
    elif sinVal < 0:
        return 360 - math.degrees(math.acos(cosVal))
    
def RAdecimalToHMS(RA): # print (not return) decimals in deg/min/sec
    decimal = (RA / 15) % 1 
    minutes = (60 * decimal)
    seconds = (minutes % 1) * 60
    minutes = math.floor(minutes)

    if RA > 0:
        integer = math.floor(RA / 15)
    else:
        integer = math.ceil(RA / 15)

    print(f"{integer}:{minutes}:{round(seconds, 2)}\"")

def DECdecimalToDMS(dec): # print (not return) decimals in deg/min/sec
    decimal = abs(dec) % 1
    minutes = (60 * decimal)
    seconds = 60 * (minutes % 1)
    minutes = math.floor(minutes)

    if dec > 0:
        integer = math.floor(dec)
    else:
        integer = math.ceil(dec)

    return str(integer) + ":" + str(minutes) + ":" + str(round(seconds, 2))
def rotationMatrix(vector, angle, axis):
    if axis == "x":
        rotMatrix = np.array([[1, 0, 0], [0, math.cos(math.radians(angle)), -math.sin(math.radians(angle))], [0, math.sin(math.radians(angle)), math.cos(math.radians(angle))]])
        return (rotMatrix @ vector)
    elif axis == "y":
        rotMatrix = np.array([[math.cos(math.radians(angle)), 0, math.sin(math.radians(angle))], [0, 1, 0], [-math.sin(math.radians(angle)), 0, math.cos(math.radians(angle))]])
        return (rotMatrix @ vector)

    else: # z axis
        rotMatrix = np.array([[math.cos(math.radians(angle)), -math.sin(math.radians(angle)), 0], [math.sin(math.radians(angle)), math.cos(math.radians(angle)), 0], [0, 0, 1]])
        return (rotMatrix @ vector)

def magnitude(h):
    return math.sqrt((h[0]) **2 + (h[1]) **2 + (h[2]) **2)


if (__name__ == "__main__"):
    # print(HMStoDeg(12, 3, 5.3, False))
    # print(DMStoDeg(-13, 45, 23.45, False))

    # print(quadrantCheck(math.cos(math.radians(0)), math.sin(math.radians(0))))
    # print(quadrantCheck(math.cos(math.radians(90)), math.sin(math.radians(90))))
    # print(quadrantCheck(math.cos(math.radians(130)), math.sin(math.radians(130))))
    # print(quadrantCheck(math.cos(math.radians(200)), math.sin(math.radians(200))))
    # print(quadrantCheck(math.cos(math.radians(300)), math.sin(math.radians(300))))
  
    # RAdecimalToHMS(-13.75651389)
    # DECdecimalToDMS(-13.75651389)

    # print(rotationMatrix(np.array([-0.23, 0.877, 0.34]), 120, "x"))
    # print(rotationMatrix(np.array([-0.23, 0.877, 0.34]), 45, "y"))
    # print(rotationMatrix(np.array([-0.23, 0.877, 0.34]), 70, "z"))

    # print(HMStoDeg("19:10:27.78", False))
    print("RA")
    print(DECdecimalToDMS(18.5809044728542)) # dec not RA function because of how the data is formatted 
    print("dec")
    print(DECdecimalToDMS(7.08998781599412))