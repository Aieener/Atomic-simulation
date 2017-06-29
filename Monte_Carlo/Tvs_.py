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

def totE_VS_T():
    #data collect from each Monte Carlo simulations when reached equlibirum states
    KE_n = [0.17234,0.34468,0.51702,0.68936,0.8617,1.03404,1.20638,
	   1.37872,1.55106,1.7234,1.89574,2.06808,2.24042,2.41276,2.5851,2.75744]
    PE_n = [-11.7338,-11.6481,-11.554,-11.4399,-11.2739,-11.0459,-10.7113,-10.1422 ,-9.50342,-8.78682,
	    -7.29134,-6.65967,-4.79413,-1.25989,-0.689232,-0.339796]

    P = [-0.150643,-0.0882192,-0.190492,-0.179469,-0.0407762,-0.00422013,-0.134795,-0.407823,
           -0.101328,0.2125,-0.203906,-0.250221,-0.142897,-0.0401648,-0.0852417,-0.03337]
    E = [] # a list of E
    for i in range(len(KE_n)):
        e = KE_n[i] + PE_n[i]
        E.append(e)

    T = [5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80] # Temperature

    Cp = []
    for i in range(len(E)-1):
        cp = (E[i+1] - E[i])/5
        Cp.append(cp)

    Tcp = [7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5,57.5,62.5,67.5,72.5,77.5] # Temperature


    fig = plt.figure()	
    plt.plot(T,KE_n,'ro-',linewidth = 1,markersize = 5,label = r'KE')
    leg = plt.legend(prop={'size':10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'Temperature (K)'
    ylabel = r'Kinetic Energy (eV)'
    plt.ylim ([0,3])
    plt.xlim ([0,90])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'Kinetic Energy vs T')
    fig.savefig("KEvsT",dpi = 300, bbox_inches ='tight')

    fig = plt.figure()	
    plt.plot(T,P,'ro-',linewidth = 1,markersize = 5,label = r'Total pressure')
    leg = plt.legend(prop={'size':10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'Temperature (K)'
    ylabel = r'Total pressure ($eV/\AA^2$)'
    #  plt.ylim ([0,4])
    #  plt.xlim ([0,90])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'Total hydrostatic pressure vs T')
    fig.savefig("PvsT",dpi = 300, bbox_inches ='tight')



    fig = plt.figure()	
    plt.plot(T,PE_n,'bo-',linewidth = 1,markersize = 5,label = r'PE')
    leg = plt.legend(prop={'size':10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'Temperature (K)'
    ylabel = r'Potential Energy (eV)'
    plt.ylim ([-13,0])
    plt.xlim ([0,90])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'Potential Energy vs T')
    fig.savefig("PEvsT",dpi = 300, bbox_inches ='tight')

    fig = plt.figure()	
    plt.plot(T,E,'go-',linewidth = 1,markersize = 5,label = r'PE')
    leg = plt.legend(prop={'size':10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'Temperature (K)'
    ylabel = r'Total Energy (eV)'
    #  plt.ylim ([-13,0])
    plt.xlim ([0,90])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'Total Energy vs T')
    fig.savefig("EvsT",dpi = 300, bbox_inches ='tight')

    fig = plt.figure()	
    plt.plot(Tcp,Cp,'co-',linewidth = 1,markersize = 5,label = r'Cp')
    leg = plt.legend(prop={'size':10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'Temperature (K)'
    ylabel = r'Heat capacity (eV/K)'
    #  plt.ylim ([-13,-2])
    #  plt.xlim ([0,90])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'Heat capacity vs T')
    fig.savefig("CpvsT",dpi = 300, bbox_inches ='tight')


def main():
    totE_VS_T()


main()
