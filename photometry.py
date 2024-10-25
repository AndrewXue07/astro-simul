import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import math
import centroid_AX
import statistics

# aperture photometry
# output total signal S + uncertainty
# signal-to-noise ratio
# instrumental magnitude + uncertainty 

# total signal
def findTotalSignal(center_x, center_y, radius, inner_sky_radius, outer_sky_radius, image):
    # image = (image / 1000) * 255
    # plt.imshow(image)
    # plt.show()
    center_x = round(center_x)
    center_y = round(center_y)
    
    # first background subtract because my code for centroid_AX was bad (not defined in terms of functions)
    p_b = 0 # pixels in annulus; for noise/uncertainty
    backgroundSignal = 0 # for noise/uncertainty

    background_pixels = []
    background_overshoot = outer_sky_radius + 4
    for i in range(center_x - background_overshoot, center_x + background_overshoot): # determine background pixels
        for j in range(center_y - background_overshoot, center_y + background_overshoot): 
            
            if ((i-center_x)**2 + (j-center_y)**2 <= outer_sky_radius**2) and ((i-center_x)**2 + (j-center_y)**2 >= inner_sky_radius**2):
                background_pixels.append(image[j, i])
                p_b += 1
                # backgroundSignal += float(image[j, i])
                # print("background signal", backgroundSignal)

    background_subtract = statistics.median(background_pixels) # find median of background pixel values
    # print("background_subtract:", background_subtract)

    n_pixel = 0  # for noise/uncertainty
    totalSignal = 0 
    for i in range(center_x - background_overshoot, center_x + background_overshoot): # subtract background and find the total signal
        for j in range(center_y - background_overshoot, center_y + background_overshoot): 
            if (i - center_x)**2 + (j - center_y)**2  < radius **2:
                n_pixel += 1
                #image[j, i] -= float(background_subtract)
                #print(image[j,i].dtype)
                totalSignal += float(image[j, i])-background_subtract
                # print(totalSignal)
    
    for i in range(center_x - background_overshoot, center_x + background_overshoot): # determine background pixels
        for j in range(center_y - background_overshoot, center_y + background_overshoot): 
            
            if ((i-center_x)**2 + (j-center_y)**2 <= outer_sky_radius**2) and ((i-center_x)**2 + (j-center_y)**2 >= inner_sky_radius**2):
                backgroundSignal += float(image[j, i]) - background_subtract


    # p_b (pixels in annulus), n_pixel (pixels in aperture), backgroundSignal (signal of background annulus)
    noise = findNoise(totalSignal, backgroundSignal, p_b, n_pixel)

    # print("totalSignal:", totalSignal)
    return totalSignal, noise

# uncertainty of total signal
def findNoise(totalSignal, backgroundSignal, p_b, n_pixel): # noise = signal uncertainty
    readNoise = 11
    darkCurrent = 10
    gain = 1
    a_b = 1 + n_pixel / p_b

    # print("backgroundSignal", backgroundSignal)
    # print("p_b", p_b)

    # before major plus sign
    firstPart = (totalSignal + n_pixel * a_b * backgroundSignal + n_pixel * a_b * darkCurrent)
    # after major plus sign
    secondPart = n_pixel * a_b * (readNoise **2 + (gain **2) / 12)

    return math.sqrt(firstPart + secondPart)

# calculate signal-to-noise ratio
def calcSNR(totalSignal, noise):
    return totalSignal / noise

# calculate instrumental magnitude + uncertainty
def calcInstMag(totalSignal, SNR):
    instMag = -2.5 * math.log(totalSignal, 10)
    uncertainty = 1.0875 / SNR

    return instMag, uncertainty

def aperturePhotometry(img, centroid_x, centroid_y, mode, apertureRadius, innerRadius, outerRadius):
    centroid_x, centroid_y, uncert_x, uncert_y = centroid_AX.findCentroid(img, centroid_x, centroid_y, mode, apertureRadius, innerRadius, outerRadius)

    image = fits.getdata(img)
    # reposition aperture at centroid location (findTotalSignal has the new centroid_x/y locations)
    totalSignal, noise = findTotalSignal(centroid_x, centroid_y, apertureRadius, innerRadius, outerRadius, image)
    SNR = calcSNR(totalSignal, noise)
    instMag, instMagUncertainty= calcInstMag(totalSignal, SNR)

    # print(f"total signal {totalSignal}; noise {noise}")
    # print("SNR", SNR)
    # print(f"instrumental magnitude {instMag} and uncertainty {instMagUncertainty}")

    return instMag, uncert_x, uncert_y, instMagUncertainty, centroid_x, centroid_y, totalSignal, noise, SNR

def calculateTotalError(errorList):
    totalErrorSquared = 0
    for error in errorList:
        totalErrorSquared += error **2

    return math.sqrt(totalErrorSquared)

def main():
    part1_instMag, part1_uncert_x, part1_uncert_y, part1_instMagUncertainty, centroid_x, centroid_y, totalSignal, noise, SNR = aperturePhotometry("aptest.fit", 490, 293, "median", 5, 8, 13)
    print(f"Part 1 Aperture Photometry: \nCentroid at: {centroid_x, centroid_y}\nSignal: {totalSignal} +/- {noise} ADU\nSNR: {SNR}\nInstrumental Magnitude: {part1_instMag} +/- {part1_instMagUncertainty}")
    # ^ sample data from part 1 
    # diff_phot = fits.getdata("diff_phot.fit")

    with open("stars_test.txt") as stars:
        data_stars = stars.readlines()
        mag_difference_list = []

        for i in range(10): # create list of magnitude differences
            # print(i)
            star = data_stars[i].strip().split()
            star_instMag, _, __, ___, ____, _____, ______, _______, ________ = aperturePhotometry("diff_phot.fit", int(star[0]), int(star[1]), "median", 5, 8, 13)
            
            mag_difference_list.append(abs(star_instMag - float(star[2])))

        mean_difference = statistics.mean(mag_difference_list) # mean of magnitude differences
        # print(mean_difference)

        asteroid_instMag, uncert_x, uncert_y, instMagUncertainty, _, __, ___, ____, _____ = aperturePhotometry("diff_phot.fit", 546, 327, "median", 5, 8, 13)
        uncertainty = calculateTotalError([uncert_x, uncert_y, instMagUncertainty])
        m_catalog = mean_difference + asteroid_instMag
        
        # print(f"asteroid instMag {asteroid_instMag}")
        # print(f"uncertainty {uncertainty} m_catalog {m_catalog}")

        print(f"\nPart 2 Differential Photometry:\nC: {mean_difference}\ndm = {uncertainty}\nCatalog Magnitude = {m_catalog} +/- {instMagUncertainty}")

    # differential photometry:
    # determine zero-point offset of instrumental magnitude
    # apply aperture photometry on all stars of stars_test.txt
    # stars_test.txt: first two columns = (x, y); third column = catalog magnitude
    # take the average of differences between instrumental and catalog magnitudes (photometry lecture slides)
    # use this average to approximate the catalog magnitude of an asteroid at (x, y) = (546, 327)
    # propagate statistical/systematic uncertainties to the asteroid approximation 



if __name__ == "__main__": 
    main()