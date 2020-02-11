# Anand viscoplastic model stress-strain curves

# MR Wei Gong
# Senior Design & Development Engineer

# Reference: Motalab, Mohammad, et al. "Determination of Anand constants for SAC solders using stress-strain or creep data." 13th InterSociety Conference on Thermal and Thermomechanical Phenomena in Electronic Systems. IEEE, 2012.

import numpy as np
import matplotlib.pyplot as plt

mat = "SAC305 - Herkommer"

# read constants from constants.txt
with open("src\\constants_{}.txt".format(mat), "r") as f:
    # Anand 9 constants
    Anand9      = np.array(f.readline().strip("\n").split(" "), dtype=float)
    # Temperature
    Temperature = np.array(f.readline().strip("\n").split(" "), dtype=float)
    # train rate
    dipsilon    = np.array(f.readline().strip("\n").split(" "), dtype=float)

# print
print("Anand 9 constants: {Anand9}".format(Anand9=Anand9))
print("Temperatures: {Temperature}".format(Temperature=Temperature))
print("Strain rates: {dipsilon}".format(dipsilon=dipsilon))

# Anand model 9 constants
s0 = Anand9[0]
QR = Anand9[1]
A  = Anand9[2]
xi = Anand9[3]
m  = Anand9[4]
h0 = Anand9[5]
ss = Anand9[6]
n  = Anand9[7]
a  = Anand9[8]

# units: MPa, K, s

for T in Temperature:
    for dip in dipsilon:

        # strain array
        ip  = np.arange(0, 0.02, 0.02/50)

        # stress calculation
        k    = dip/A*np.exp(QR/T)
        sig  = 1/xi*np.arcsinh(np.power(k, m))
        temp = np.power(ss*np.power(k, n) - s0, 1-a) + (a-1)*(h0*np.power(ss*np.power(k, n), -a))*ip
        sig *= ss*np.power(k, n)-np.power(temp, 1/(1-a))

        # plot
        plt.plot(ip, sig, label='T = {temperature}K, Strain Rate = {strain_rate}'.format(temperature=T, strain_rate=dip))

plt.xlabel("Strain")
plt.ylabel("Stress/MPa")
plt.legend(loc = 'lower right')
plt.title("Anand Viscoplastic Model Stress-Strain Curve of {}".format(mat))
plt.show()

