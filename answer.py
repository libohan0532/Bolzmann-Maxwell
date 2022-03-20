import numpy as np
from math import pi, e
import random
import matplotlib.pyplot as plt

ev= 1.602170565 * 10 ** (-19)
k = 1.380648 * 10 ** (-23)
class MBdist():
    global p
    # using this function to acquire posibillity
    def maxwell_boltzmann(self, temperature, mass, speed):
        mass = mass * 1.6726219 * 10 ** (-27)  # weight of H atom

        p = 4 * pi * speed ** 2 * (mass / (2 * pi * k * temperature)) ** (3.0 / 2.0) * \
            e ** (-mass * speed ** 2 / (2 * k * temperature))

        return p

    # prepared the cumulative weights
    def prepare(self):
        self.z[0] = self.maxwell_boltzmann(self.T, self.m, 0)
        for i in range(1, self.L):
            self.z[i] = self.z[i - 1] + self.maxwell_boltzmann(self.T, self.m, i)

    def __init__(self, T, m, L):
        self.T = T
        self.m = m
        self.L = L
        self.z = np.zeros(L)
        self.prepare()

    # takes bins, returns em full with n figures
    def sample_n(self, bins, n):
        # choose rand num from dist
        ran = random.choices(np.arange(self.L), cum_weights=self.z, k=n)#get enough number of sample
        # bin it
        for s in ran:
            ind = int(s * (len(bins) - 1) / self.L)
            bins[ind] = bins[ind] + 1

    def set_temperature(self, new_t):
        self.T = new_t
        self.prepare()

    def set_mass(self, new_m):
        self.m = new_m
        self.prepare()


r = 1  # how many subplots
l = 25  # x number of bars in histogram
xdiv = 1  # number of subplots in a row
T = 0.5*ev/k  # the temp given in the question
m = 1.008   # weight of H atom
L = 30000  # the largest scale of velocity in the system
o = MBdist(T, m, L)  # define the object to output
 # produce r subplots
bins = np.zeros(l, int)  # bins empty the array
o.set_temperature(T)
o.sample_n(bins, 1000)
plt.subplot(int(r / xdiv), xdiv, 1)
plt.xticks(np.arange(0, L, 2 * L / 10))
plt.bar(np.arange(0, L, L / l), bins, width=L / l, align='center', color='red', edgecolor='yellow')
plt.xlabel('speed (m/s)')
plt.ylabel('particle count')
plt.subplots_adjust(wspace=0.4, hspace=0.55)  # Set detailed figures for pic
mass = m * 1.6726219 * 10 ** (-27)
v = np.linspace(0,30000,10000,True)
fv = 4 * np.pi * v ** 2 * (mass / (2 * np.pi * k * T)) ** (3.0 / 2.0) * e ** (
                    -mass * v ** 2 / (2 * k * T))*1000000
plt.plot(v,fv,linewidth=5,color='cyan')
title = 'Distribution of the speeds of Hydrogen particles (m=1, KbT=0.5Ev)'
plt.title(title)

plt.savefig('answer.png', dpi=200)
del o
