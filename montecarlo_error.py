import numpy as np
import XueOD
from odlib import *
import matplotlib.pyplot as plt

def renameThisFunction(RA_uncertainties, dec_uncertainties, fileName):
    a_jpl = 3.262856683191058
    e_jpl = 0.6170424842607419
    i_jpl = 10.01124250389698
    Omega_jpl = 147.9970125329559
    omega_jpl = 142.4457338116572
    MA_jpl = 1.966530717490744 # change this value later

    file = open(fileName)
    RAs = []
    decs = []

    for line in file.readlines():
        tempArr = line.split(" ")
        RAs.append(HMStoDeg(tempArr[1], False)) 
        decs.append(DMStoDeg(tempArr[2], False))

    # print(RAs, decs)
    
    a_vals = []
    e_vals = []
    i_vals = []
    Omega_vals = []
    omega_vals = []
    MA_vals = []
    
    for i in range(2000):
        # current_RA_offset = np.random.normal(0, 0.1)#uncertainty)
        # current_dec_offset = np.random.normal(0, 0.1)#uncertainty)

        # print("uncertainty:", current_offset)
        modified_RAs = RAs.copy()
        modified_decs = decs.copy()

        for j in range(len(RAs)): # remove this while loop and manually add (6 separate variables) to prevent the time stall? 
            modified_RAs[j] += np.random.normal(-RA_uncertainties[j], RA_uncertainties[j])#current_RA_offset#np.random.normal(0, 0.1)#RA_uncertainties[j])
            modified_decs[j] += np.random.normal(-dec_uncertainties[j], dec_uncertainties[j])#current_dec_offset#np.random.normal(0, 0.1)#dec_uncertainties[j])

        a, e, i, Omega, omega, meanAnomaly = XueOD.main(True, modified_RAs, modified_decs, 2460518.750000) # 2460518.750000 = due date
        a_vals.append(a)
        e_vals.append(e)
        i_vals.append(i)
        Omega_vals.append(Omega)
        omega_vals.append(omega)
        MA_vals.append(meanAnomaly)

    mean_a = np.mean(a_vals)
    std_a = np.std(a_vals)
    percent_error_a = (abs(mean_a - a_jpl) / a_jpl) * 100

    mean_e = np.mean(e_vals)
    std_e = np.std(e_vals)
    percent_error_e = (abs(mean_e - e_jpl) / e_jpl) * 100

    mean_i = np.mean(i_vals)
    std_i = np.std(i_vals)
    percent_error_i = (abs(mean_i - i_jpl) / i_jpl) * 100

    mean_Omega = np.mean(Omega_vals)
    std_Omega = np.std(Omega_vals)
    percent_error_Omega = (abs(mean_Omega - Omega_jpl) / Omega_jpl) * 100

    mean_omega = np.mean(omega_vals)
    std_omega = np.std(omega_vals)
    percent_error_omega = (abs(mean_omega - omega_jpl) / omega_jpl) *100

    mean_MA = np.mean(MA_vals)
    std_MA = np.std(MA_vals)
    percent_error_MA = (abs(mean_MA - MA_jpl) / MA_jpl) * 100

    print(f"Mean value of semi-major axis: {mean_a} with standard deviation of {std_a} and % error of {percent_error_a}%")
    print(f"Mean value of eccentricity: {mean_e} with standard deviation of {std_e} and % error of {percent_error_e}%")
    print(f"Mean value of inclination: {mean_i} with standard deviation of {std_i} and % error of {percent_error_i}%")
    print(f"Mean value of longitude of ascending node: {mean_Omega} with standard deviation of {std_Omega} and % error of {percent_error_Omega}%")
    print(f"Mean value of argument of perihelion: {mean_omega} with standard deviation of {std_omega} and % error of {percent_error_omega}%")
    print(f"Mean value of mean anomaly: {mean_MA} with standard deviation of {std_MA} and % error of {percent_error_MA}%")


    # chart fonts
    subTitleFont = {"size": 10}
    xAxisLabelFont = {"size": 8}
    yAxisLabelFont = {"size": 10}
    # semi-major axis
    plt.subplot(3, 2, 1)
    plt.hist(a_vals, bins = 100)
    plt.xlabel("Semi-Major Axis (a) Values", fontdict=xAxisLabelFont)
    plt.ylabel("Frequency", fontdict=yAxisLabelFont)
    plt.title("Semi-Major Axis Value Frequencies", fontdict=subTitleFont)
    plt.axvline(x = mean_a, ymin = 0, ymax = 1, c = "r")

    # eccentricity
    plt.subplot(3, 2, 2)
    plt.hist(e_vals, bins = 100)
    plt.xlabel("Eccentricity (e) Values", fontdict=xAxisLabelFont)
    plt.ylabel("Frequency", fontdict=yAxisLabelFont)
    plt.title("Eccentricity Value Frequences", fontdict=subTitleFont)
    plt.axvline(x = mean_e, ymin = 0, ymax = 1, c = "r")


    # inclination
    plt.subplot(3, 2, 3)
    plt.hist(i_vals, bins = 100)
    plt.xlabel("Inclination (i) Values", fontdict=xAxisLabelFont)
    plt.ylabel("Frequency", fontdict=yAxisLabelFont)
    plt.title("Inclination Value Frequences", fontdict=subTitleFont)
    plt.axvline(x = mean_i, ymin = 0, ymax = 1, c = "r")

    # longitude of ascending node
    plt.subplot(3, 2, 4)
    plt.hist(Omega_vals, bins = 100)
    plt.xlabel("Longitude of Ascending Node (Ω) Values", fontdict=xAxisLabelFont)
    plt.ylabel("Frequency", fontdict=yAxisLabelFont)
    plt.title("Longitude of Ascending Node Value Frequences", fontdict=subTitleFont)
    plt.axvline(x = mean_Omega, ymin = 0, ymax = 1, c = "r")

    # argument of perihelion
    plt.subplot(3, 2, 5)
    plt.hist(omega_vals, bins = 100)
    plt.xlabel("Argument of Perihelion (ω) Values", fontdict=xAxisLabelFont)
    plt.ylabel("Frequency", fontdict=yAxisLabelFont)
    plt.title("Argument of Perihelion Value Frequences", fontdict=subTitleFont)
    plt.axvline(x = mean_omega, ymin = 0, ymax = 1, c = "r")

    # mean anomaly
    plt.subplot(3, 2, 6)
    plt.hist(MA_vals, bins = 100)
    plt.xlabel("Mean Anomaly (M)) Values", fontdict=xAxisLabelFont)
    plt.ylabel("Frequency", fontdict=yAxisLabelFont)
    plt.title("Mean Anomaly Value Frequences", fontdict=subTitleFont)
    plt.axvline(x = mean_MA, ymin = 0, ymax = 1, c = "r")

    plt.suptitle("Monte Carlo Orbital Element Frequency Histograms")
    plt.tight_layout()
    plt.show()
    
def main():
    renameThisFunction([0.0007723697871, 0.000318, 0.000119], [0.0002239556032, 0.000154, 0.000134], "finalorbitalelements_input.txt")


if __name__ == "__main__":
    main()