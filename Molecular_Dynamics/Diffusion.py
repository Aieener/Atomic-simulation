# Diffusion.py
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

def diffuse():
    lnD = []
    T_rev = [];
    with open("selfdiffusion_v2.txt","r") as file:
        for line in file:
            words = line.split()
            trev = float(words[0])
            lnd = float(words[1])
            lnD.append(lnd)
            T_rev.append(trev)


    fig = plt.figure()
    plt.plot(T_rev,lnD,'ro-',linewidth = 0.5,markersize = 5,label = r'lnD')
    leg = plt.legend(prop={'size':10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'1/T $(K^{-1})$'
    ylabel = r'$ln(D/\AA^2fs^{-1})$'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim([0.012,0.021])
    plt.title(r'lnD vs 1/T')
    fig.savefig("selfDiffuse_v2",dpi = 300, bbox_inches ='tight')


def main():
    diffuse()

main()
