import numpy as np
from odlib import *

# input: file containing (x, y) centroids and (RA, dec) J2000 coordinates of 10 reference stars
# (x, y) centroid of unknown object (centroid_AX.py)

# output: 6 plate constants (degrees and degrees/pixels)
# (RA, dec) J2000 coordinates of unknown object
# uncertainty of (RA, dec)

# inverse: np.linalg.inv()
# dot product: np.dot()

def LSPR(file, x_centroid, y_centroid): 
    with open(file) as data: # x centroid, y centroid, RA (hours), declination (degrees)
        lines = data.readlines()
        stars = []

        for i in range(len(lines)):
            stars.append(lines[i].split())

        # print(stars)
        RA, dec, RA_uncert, dec_uncert, RA_constants, dec_constants = calcCoordinates(stars, x_centroid, y_centroid)
        
        # print(RA_uncert, dec_uncert)
        print(
    f"""***************
plate constants
***************
b1: {RA_constants[0]} deg
b2: {dec_constants[0]} deg
a11: {RA_constants[1]} deg/pix
a12: {RA_constants[2]} deg/pix
a21: {dec_constants[1]} deg/pix
a22: {dec_constants[2]} deg/pix
***************
uncertainty
***************
 RA: {RA_uncert} arcsec
Dec: {dec_uncert} arcsec
*********************************************
astrometry for 
(x, y) = ({x_centroid}, {y_centroid})
*********************************************
 RA = {RA}
Dec = {dec}
""")


def calcCoordinates(stars, x_centroid, y_centroid): # calculate RA or declination with the two linear equations
    RA_constants, dec_constants = calculatePlateConstants(stars)

    RA = RA_constants[0] + (RA_constants[1] * x_centroid) + (RA_constants[2] * y_centroid)
    dec = dec_constants[0] + (dec_constants[1] * x_centroid) + (dec_constants[2] * y_centroid)

    RA_uncert, dec_uncert = calculateUncertainty(stars, RA_constants, dec_constants)
    # print(type(RA), type(dec))
    RA = RAdecimalToHMS(RA)
    dec = DECdecimalToDMS(dec)

    # print("RA", RA)
    # print("dec", dec)
    # print("RA_uncert", RA_uncert)
    # print("dec_uncert", dec_uncert)

    return RA, dec, RA_uncert, dec_uncert, RA_constants, dec_constants

def calculatePlateConstants(stars): # with matrices
    firstArr = np.array([[len(stars), 0., 0.],
                        [0., 0., 0.],
                        [0., 0., 0.]]
                        )
    # take inverse later 

    # calculate RA constants
    secondArr = np.array([0., 0., 0.])
    
    for i in range(len(stars)):
        # first array in matrix multiplication
        firstArr[0, 1] += float(stars[i][0])
        firstArr[0, 2] += float(stars[i][1]) 
        # first row done
        
        firstArr[1, 0] += float(stars[i][0])
        firstArr[1, 1] += float(stars[i][0]) **2
        firstArr[1, 2] += float(stars[i][0]) * float(stars[i][1])
        # second row done

        firstArr[2, 0] += float(stars[i][1])
        firstArr[2, 1] += float(stars[i][0]) * float(stars[i][1])
        firstArr[2, 2] += float(stars[i][1]) **2
        # third row done

        # second matrix in matrix multiplication 
        secondArr[0] += HMStoDeg(stars[i][2], False)
        secondArr[1] += HMStoDeg(stars[i][2], False) * float(stars[i][0])
        secondArr[2] += HMStoDeg(stars[i][2], False) * float(stars[i][1])

    RA_constants = np.linalg.inv(firstArr) @ secondArr
    # print(RA_constants)

    # calculate declination constants
    thirdArr = np.array([0., 0., 0.])
    for i in range(len(stars)):
        thirdArr[0] += DMStoDeg(stars[i][3], False)
        thirdArr[1] += DMStoDeg(stars[i][3], False) * float(stars[i][0])
        thirdArr[2] += DMStoDeg(stars[i][3], False) * float(stars[i][1])

    dec_constants = np.linalg.inv(firstArr) @ thirdArr
    # print(dec_constants)

    return RA_constants, dec_constants

def calculateUncertainty(stars, RA_constants, dec_constants):
    # calculate RA uncertainty

    RA_uncert = 0.
    RA_uncert_numerator = 0.

    for i in range(len(stars)):
        RA_uncert_numerator += (HMStoDeg(stars[i][2], False) - float(RA_constants[0]) - (float(RA_constants[1]) * float(stars[i][0])) - (float(RA_constants[2]) * float(stars[i][1]))) **2

    RA_uncert = math.sqrt(RA_uncert_numerator / (len(stars) - 3))
    RA_uncert = round(3600 * RA_uncert, 2)

    # calculate dec uncertainty

    dec_uncert = 0.
    dec_uncert_numerator = 0.

    for i in range(len(stars)):
        dec_uncert_numerator += (DMStoDeg(stars[i][3], False) - float(dec_constants[0]) - (float(dec_constants[1]) * float(stars[i][0])) - (float(dec_constants[2]) * float(stars[i][1]))) **2

    dec_uncert = math.sqrt(dec_uncert_numerator / (len(stars) - 3))
    dec_uncert = round(3600* dec_uncert, 2)



    return RA_uncert, dec_uncert

LSPR("LSPRtestinput1.txt", 484.35, 382.62)