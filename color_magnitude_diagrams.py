import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

tab = pd.read_csv("newgooddata.csv")
tab_isochrone = pd.read_table("tmp1721514650.txt", delim_whitespace = True, skiprows = 8)

x_scatterplot = np.array([tab["g"]]) - np.array([tab["r"]])
y_scatterplot = np.array(tab["r"]) - 9.8494
# M = np.array([r + 9.8494 for r in tab["r"]])
x_isochrone = np.array(tab_isochrone["sdss_g"]) - np.array(tab_isochrone["sdss_r"])
y_isochrone = np.array(tab_isochrone["sdss_r"])

plt.plot(x_isochrone, y_isochrone)
plt.scatter(x_scatterplot, y_scatterplot, s = 1, c = "purple")
plt.xlabel("g - r")
plt.ylabel("r")
plt.ylim(25, -5)
plt.savefig("xyplot.png")
# # center of main sequence = (1, 17)
# # 15 on isochrone: g-r=1.18; r=7.1506
# # difference (distance modulus): 17-7.1506 = 9.8494


# tmp2billionyears.txt
plt.clf()
tab = pd.read_csv("newgooddata.csv")
tab_isochrone = pd.read_table("tmp2billionyears.txt", delim_whitespace = True, skiprows = 8)

x_scatterplot = np.array([tab["g"]]) - np.array([tab["r"]])
y_scatterplot = np.array(tab["r"]) - (17 - 7.9086)
x_isochrone = np.array(tab_isochrone["sdss_g"]) - np.array(tab_isochrone["sdss_r"])
y_isochrone = np.array(tab_isochrone["sdss_r"])

plt.plot(x_isochrone, y_isochrone)
plt.scatter(x_scatterplot, y_scatterplot, s = 1, c = "purple")
plt.xlabel("g - r")
plt.ylabel("r")
plt.ylim(25, -5)
plt.savefig("xyplot2.png")


# 3 billion
plt.clf()
tab = pd.read_csv("newgooddata.csv")
tab_isochrone = pd.read_table("tmp3billion.txt", delim_whitespace = True, skiprows = 8)

x_scatterplot = np.array([tab["g"]]) - np.array([tab["r"]])
y_scatterplot = np.array(tab["r"]) - (17 - 9.6252)
x_isochrone = np.array(tab_isochrone["sdss_g"]) - np.array(tab_isochrone["sdss_r"])
y_isochrone = np.array(tab_isochrone["sdss_r"])

plt.plot(x_isochrone, y_isochrone)
plt.scatter(x_scatterplot, y_scatterplot, s = 1, c = "purple")
plt.xlabel("g - r")
plt.ylabel("r")
plt.ylim(25, -5)
plt.savefig("xyplot3.png")


# 4 billion
plt.clf()
tab = pd.read_csv("newgooddata.csv")
tab_isochrone = pd.read_table("tmp4billion.txt", delim_whitespace = True, skiprows = 8)

x_scatterplot = np.array([tab["g"]]) - np.array([tab["r"]])
y_scatterplot = np.array(tab["r"]) - (17 - 11.0495)
x_isochrone = np.array(tab_isochrone["sdss_g"]) - np.array(tab_isochrone["sdss_r"])
y_isochrone = np.array(tab_isochrone["sdss_r"])

plt.plot(x_isochrone, y_isochrone)
plt.scatter(x_scatterplot, y_scatterplot, s = 1, c = "purple")
plt.xlabel("g - r")
plt.ylabel("r")
plt.ylim(25, -5)
plt.savefig("xyplot4.png")

# # 5 billion
plt.clf()
tab = pd.read_csv("newgooddata.csv")
tab_isochrone = pd.read_table("tmp5billion.txt", delim_whitespace = True, skiprows = 8)

x_scatterplot = np.array([tab["g"]]) - np.array([tab["r"]])
y_scatterplot = np.array(tab["r"]) - (17 - 6.5946)
x_isochrone = np.array(tab_isochrone["sdss_g"]) - np.array(tab_isochrone["sdss_r"])
y_isochrone = np.array(tab_isochrone["sdss_r"])

plt.plot(x_isochrone, y_isochrone)
plt.scatter(x_scatterplot, y_scatterplot, s = 1, c = "purple")
plt.xlabel("g - r")
plt.ylabel("r")
plt.ylim(25, -5)
plt.savefig("xyplot5.png")


# 6 billion
plt.clf()
tab = pd.read_csv("newgooddata.csv")
tab_isochrone = pd.read_table("tmp6billion.txt", delim_whitespace = True, skiprows = 8)

x_scatterplot = np.array([tab["g"]]) - np.array([tab["r"]])
y_scatterplot = np.array(tab["r"]) - (17 - 6.4364)
x_isochrone = np.array(tab_isochrone["sdss_g"]) - np.array(tab_isochrone["sdss_r"])
y_isochrone = np.array(tab_isochrone["sdss_r"])

plt.plot(x_isochrone, y_isochrone)
plt.scatter(x_scatterplot, y_scatterplot, s = 1, c = "purple")
plt.xlabel("g - r")
plt.ylabel("r")
plt.ylim(25, -5)
plt.savefig("xyplot6.png")