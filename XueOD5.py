from odlib import *
import numpy as np
import math

def mog(tau1, tau3, r2, r2_dot, order): # order = 3rd or 4th
    r2_mag = magnitude(r2)
    # r2_dot_mag = magnitude(r2_dot)

    u = r2_mag **(-3)
    w = np.dot(r2_dot, r2_dot) / (r2_mag **2)
    z = np.dot(r2, r2_dot) / (r2_mag **2)

    if (order == 3): 
        f1 = 1 - (0.5*u*tau1 **2) + (0.5*u*z*tau1 **3)
        f3 = 1 - (0.5*u*tau3 **2) + (0.5*u*z*tau3 **3)

        g1 = tau1 - (1/6)*u*tau1 **3
        g3 = tau3 - (1/6)*u*tau3 **3

        return f1, f3, g1, g3

    
    if (order == 4):
        f1 = 1 - (0.5*u*tau1 **2) + (0.5*u*z*tau1 **3) + (1/24)*(-15*u*z**2 + 3*u*w - 2*u **2)*tau1 **4
        f3 = 1 - (0.5*u*tau3 **2) + (0.5*u*z*tau3 **3) + (1/24)*(-15*u*z**2 + 3*u*w - 2*u **2)*tau3 **4
        
        g1 = tau1 - (1/6)*u*tau1 **3 + 0.25*u*z*tau1 **4
        g3 = tau3 - (1/6)*u*tau3 **3 + 0.25*u*z*tau3 **4

        return f1, f3, g1, g3

if __name__ == "__main__":
    print("Test case with 3rd-order series")
    mog(-0.32618617484601165, 0.0508408854033231, [0.26799552002875776, -1.3726277901924608, -0.5026729612047128], [0.8456809141954584, -0.3838382184712308, 0.14215854191172816], 3)

    print("Test case with 4th-order series")
    mog(-0.3261857571141891, 0.05084081855693949, [0.26662393644794813, -1.381475976476564, -0.5048589337503169], [0.8442117090940343, -0.39728396707075087, 0.14202728258915864], 4)