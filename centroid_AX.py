""" Calculate the centroid of a star. 
    Checks: 
    inner_sky_radius > outer_sky_radius
    radius > inner_sky_radius
    radius, inner_sky_radius, or outer_sky_radius <= 0 
    bkg_method not one of the options (mode, median, mean)

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import math
import statistics


def findCentroid(fits_file, target_x, target_y, bkg_method, radius=3.5, inner_sky_radius=5.5, outer_sky_radius = 7.5):
    image = fits.getdata(fits_file)

    if (inner_sky_radius > outer_sky_radius or radius > inner_sky_radius):
        print("annulus or aperture sizes conflicting")
        return 
    elif (radius <= 0 or inner_sky_radius <= 0 or outer_sky_radius <= 0): 
        print("please input positive radii values")
        return 

    center_x = target_x
    center_y = target_y
    # inner_sky_radius = inner_sky_radius
    # outer_sky_radius = outer_sky_radius


    background_pixels = []
    background_overshoot = int(round(outer_sky_radius)) + 4
    for i in range(int(round(center_x)) - background_overshoot, int(round(center_x)) + background_overshoot):
        for j in range(int(round(center_y)) - background_overshoot, int(round(center_y)) + background_overshoot): 
            
            if ((i-center_x)**2 + (j-center_y)**2 < outer_sky_radius**2) and ((i-center_x)**2 + (j-center_y)**2 > inner_sky_radius**2):
                background_pixels.append(image[j, i])

    if (bkg_method == "median"):
        # print(background_pixels)
        background_subtract = statistics.median(background_pixels) # find background average/median
        # print(background_subtract)
    elif (bkg_method == "mean"):
        background_subtract = statistics.mean(background_pixels)
    elif (bkg_method == "mode"):
        background_subtract = statistics.mode(background_pixels)
    else: 
        print("Please input \"median,\" \"mean,\" or \"mode.\"")
        return

    for i in range(int(round(center_x)) - background_overshoot, int(round(center_x)) + background_overshoot):
        for j in range(int(round(center_y)) - background_overshoot, int(round(center_y)) + background_overshoot): 
            if (i-center_x)**2 + (j-center_y)**2  < radius **2:
                image[j, i] -= background_subtract


    # calculate y centroid
    yNum = 0 # numerator
    yDen = 0 # denominator
    for j in range(int(round(center_y)) - background_overshoot, int(round(center_y)) + background_overshoot):
        sum = 0
        for i in range(int(round(center_x)) - background_overshoot, int(round(center_x)) + background_overshoot): 
            # xNum += i * image[i]
            # xDen += i
            
            if (j < len(image) and j > 0 and i > 0 and i < len(image[j])):
                sum += image[j, i]
            
        
        yNum += j * sum
        yDen += sum
        # if yDen == 0: 
        #     y_centroid = 0
        # else: 

        y_centroid = (yNum / yDen)

    # calculate x centroid
    xNum = 0 # numerator
    xDen = 0 # denominator
    for i in range(int(round(center_x)) - background_overshoot, int(round(center_x)) + background_overshoot):
        sum = 0
        for j in range(int(round(center_y)) - background_overshoot, int(round(center_y)) + background_overshoot): 
            # xNum += i * image[i]
            # xDen += i
            
            if (j < len(image) and j > 0 and i > 0 and i < len(image[j])):
                sum += image[j, i]
        
        xNum += i * sum
        xDen += sum
        # if xDen == 0: 
        #     x_centroid = 0
        # else: 
        x_centroid = (xNum / xDen)

    if (image[math.floor(y_centroid), math.floor(x_centroid)] < 1):
        print("No star found")



    # calculate uncertainty

    # y uncertainty
    N = 0
    yUncertaintyNumerator = 0
    # yUncertaintyDenominator = 0

    for j in range(int(round(center_y)) - background_overshoot, int(round(center_y)) + background_overshoot):
        for i in range(int(round(center_x)) - background_overshoot, int(round(center_x)) + background_overshoot): 
            if (j < len(image) and j > 0 and i > 0 and i < len(image[j])):
            
                N += image[j][i]
                yUncertaintyNumerator += image[j][i] * (i - y_centroid) **2

    # print(yUncertaintyNumerator)
    # print(N)
    # print(N)
    yUncertainty = math.sqrt(yUncertaintyNumerator / (N) / (N-1))

    # x uncertainty AND total signal
    # totalSignal = 0.
    N = 0
    xUncertaintyNumerator = 0

    for i in range(int(round(center_x)) - background_overshoot, int(round(center_x)) + background_overshoot):
        for j in range(int(round(center_y)) - background_overshoot, int(round(center_y)) + background_overshoot): 
            # totalSignal += j
            # print(j)
            if (j < len(image) and j > 0 and i > 0 and i < len(image[j])):
                N += image[j][i]
                xUncertaintyNumerator += image[j][i] * (j - x_centroid) **2

    # print(N + 1)
    xUncertainty = math.sqrt(xUncertaintyNumerator / (N) / (N-1))
    # y_centroid = 153.6

    return x_centroid, y_centroid, xUncertainty, yUncertainty
    # return 350.7806, 153.5709, xUncertainty, yUncertainty


# findCentroid("sampleimage.fits", 351, 154, "median", 3, 5, 10)
# centroid_x, centroid_y = findCentroid("sampleimage.fits", 351, 154, "median", 3, 5, 10)

# for i in range (5):
if __name__ == "__main__": 
    centroid_x, centroid_y, uncert_x, uncert_y = findCentroid("sampleimage.fits", 351, 154, "median", 3.5, 5.5, 7.5) # 3.5, 5.5, 7.5
    print(uncert_x, uncert_y)
    if abs(centroid_x - 350.7806) < 0.1 and abs(centroid_y - 153.5709) < 0.2:
        print("centroid calculation CORRECT")
    else:
        print(
            "centroid calculation INCORRECT, expected (350.7806, 153.5709), got ({}, {})".format(
            centroid_x, centroid_y))

    

# return x & y positions of center of target stars
# return uncertainties
# function should be able to select the size of the circular aperture and inner/out diamters of annulus
# user should be able to choose between median/mean/mode as a method to calculate background value 
# background pixel include: middle of pixel (x & y) must be within the aperture/annulus
