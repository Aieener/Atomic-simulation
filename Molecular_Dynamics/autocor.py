# hw4.py
# HW4 Molecular Dynamics
# Author: Yuding Ai
# Date: 2017.03.21

import math
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import scipy.stats  as stats
import collections
import matplotlib as mpl
from matplotlib import rc
from itertools import groupby
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

def autocor():
    Vall = []
    T = [];
    with open("auto.txt","r") as file:
        count = 0
        for line in file:
            words = line.split()
            V_t = []
            for i in words:
                v = float(i)
                V_t.append(v)
            Vall.append(V_t)

    Auto = []
    T = []

    for k in range(5000):
        auto = 0
        for i in range(5000):
            for j in range(len(Vall[i])):
                vo = Vall[i][j]
                vt = Vall[i+k][j]
                auto = auto + vo*vt
        auto = auto/len(Vall[i])
        auto = auto/5000
        Auto.append(auto)
        T.append(k)
        print k


    fig = plt.figure()
    plt.plot(T,Auto,'c-',linewidth = 0.5,markersize = 5,label = r'Pressure')
    leg = plt.legend(prop={'size':10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'$N^{th}$ MD simulation'
    ylabel = r'Pressure ($eV/\AA^2$)'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'Pressure vs MD')
    fig.savefig("Auto",dpi = 300, bbox_inches ='tight')


def main():
    autocor()

main()
