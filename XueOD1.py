# runs fine

from odlib import *
import numpy as np

input = open("Xue_Input.txt")
lineArr = input.readlines()


def getAngMomentum(lines): # angular momentum = r x p
    position = np.array(lines[0].strip().split(" ")).astype(float)
    velocity = (1/0.0172020989484) * np.array(lines[1].strip().split(" ")).astype(float)
    # print(position)
    # print(velocity)

    return np.cross(position, velocity)

print(getAngMomentum(lineArr))

input.close()

# output angular momentum vector rounded to 6 decimal places 
# units: AU (distance), Gaussian days (time)

# write "runs fine" when submitting 
# create zip files with all necessary files 