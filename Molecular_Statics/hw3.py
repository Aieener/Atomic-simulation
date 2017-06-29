# hw3.py
# HW3 Molecular statics relaxation
# Author: Yuding Ai
# Date: 2017.03.01

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

def plotconfig(txtname,filename):
    '''take the config.txt to plot the atoms'''
    X = [] # a list of Xpos
    Y = [] # a list of Ypos 
    Sxx = [] # a list of signa_xx
    Syy = [] # a list of signa_yy
    Sxy = [] # a list of signa_xy
    P = []
    with open(txtname,"r") as file:
        for line in file:
            words = line.split()
            x = float(words[0]) #take the value
            y = float(words[1]) #take the value
            sxx = float(words[2]) #take the value
            syy = float(words[3]) #take the value
            sxy = float(words[4]) #take the value
            p = -0.5*(syy+sxx)
            
            X.append(x); #append x value into X
            Y.append(y); #append y value into Y
            Sxx.append(sxx);
            Syy.append(syy);
            Sxy.append(sxy);
            P.append(p);

    fig = plt.figure()	
    plt.plot(X,Y,'ro',linewidth = 0.8)
    plt.ylim([-4,80])
    plt.xlim([-4,90])
    xlabel = r'$x_{coordinate}\ unit =\AA $'
    ylabel = r'$y_{coordinate}\ unit =\AA $'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title('Configuration of Atoms')
    fig.savefig(filename,dpi = 300, bbox_inches ='tight')
    
    #Plot the stresses as heatmap
    #---------------simgaxx--------------------------------
    fig, ax = plt.subplots(1)
    plt.scatter(X,Y,s = 150,marker = '8', c = Sxx,alpha=0.7)
    cbar = plt.colorbar()
    cbar.set_label(r"$eV/{\AA}^2$")
    xlabel = r'X_{pos}'
    ylabel = r'Y_{pos}'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.locator_params(axis='y', nticks=3)
    plt.locator_params(axis='x', nticks=8)
    fname = "sigma_xx" + filename 
    titlename = r"$\sigma_{xx}$"
    plt.title(titlename)
    fig.savefig(fname,dpi=300)

    #---------------simgayy--------------------------------
    fig, ax = plt.subplots(1)
    plt.scatter(X,Y,s = 150,marker = '8', c = Syy,alpha=0.7)
    cbar = plt.colorbar()
    cbar.set_label(r"$eV/{\AA}^2$")
    xlabel = r'X_{pos}'
    ylabel = r'Y_{pos}'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.locator_params(axis='y', nticks=3)
    plt.locator_params(axis='x', nticks=8)
    fname = "sigma_yy" + filename 
    titlename = r"$\sigma_{yy}$"
    plt.title(titlename)
    fig.savefig(fname,dpi=300)

    #---------------simgaxy--------------------------------
    fig, ax = plt.subplots(1)
    plt.scatter(X,Y,s = 150,marker = '8', c = Sxy,alpha=0.7)
    cbar = plt.colorbar()
    cbar.set_label(r"$eV/{\AA}^2$")
    xlabel = r'X_{pos}'
    ylabel = r'Y_{pos}'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.locator_params(axis='y', nticks=3)
    plt.locator_params(axis='x', nticks=8)
    fname = "sigma_xy" + filename 
    titlename = r"$\sigma_{xy}$"
    plt.title(titlename)
    fig.savefig(fname,dpi=300)

    #---------------  P -----------------------------------
    fig, ax = plt.subplots(1)
    plt.scatter(X,Y,s = 150,marker = '8', c = P,alpha=0.7)
    cbar = plt.colorbar()
    cbar.set_label(r"$eV/{\AA}^2$")
    xlabel = r'X_{pos}'
    ylabel = r'Y_{pos}'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.locator_params(axis='y', nticks=3)
    plt.locator_params(axis='x', nticks=8)
    fname = "hydrop" + filename 
    titlename = r"$\sigma_{xy}$"
    plt.title(titlename)
    fig.savefig(fname,dpi=300)




def plotrdf(txtname,filename):
    '''take the rdf.txt to plot the rdf function'''

    R = [] # a list of distance distribution
    G = [] # a list of distance distribution
    with open(txtname,"r") as file:
        for line in file:
            words = line.split()
            r = float(words[0]) #take the value
            g = float(words[1]) #take the value
            R.append(r); #append r value into r
            G.append(g); #append g(r) value into G
    fig = plt.figure()	

    if txtname == "rdf.txt":
        plt.stem(R, G, linefmt='b-', linewidth = 0.4,markerfmt='rx', basefmt='r-',label = r'RDF')
    else:
        plt.plot(R,G,'rx-',linewidth = 0.4,markersize = 1,label = r'RDF')
    xlabel = r'$r\ unit =\AA $'
    ylabel = r'$g(r)\ $'
    plt.title('Radial distribution function g(r)')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    fig.savefig(filename, dpi=180, bbox_inches='tight')

def plot_vsN(tx1,tx2,tx3,tx4,tx5,EN,SxxN,SyyN,SxyN,PN,FN,E8,F8):
    '''make plots for stress vs N; potential energy vs N
    maximum Force vs N and hydrostatic pressure vs N'''

    #-------------------------------------------------------------
    # step 1 load files
    #-------------------------------------------------------------
    #--------------------lambda = 0.5 ----------------------------
    N = [] # a list of interation
    E = [] # a list of PE 
    SXX = [] #sigma_xx
    SYY = [] #sigma_yy
    SXY = [] #sigma_xy
    P = [] #hydrostatic pressure
    F = [] #force
    with open(tx1,"r") as file:
        for line in file:
            words = line.split()
            n = float(words[0]) #take the value
            e = float(words[1]) #take the value
            sxx = float(words[2]) #take the value
            syy = float(words[3]) #take the value
            sxy = float(words[4]) #take the value
            p = float(words[5]) #take the value
            f = float(words[6]) #take the value

            N.append(n); #append x value into X
            E.append(e); #append y value into Y
            SXX.append(sxx); #append sxx value into SXX
            SYY.append(syy); #append syy value into SYY
            SXY.append(sxy); #append sxy value into SXY
            P.append(p);
            F.append(f);

    #--------------------lambda = 1 ------------------------------
    N1= [] # a list of interation
    E1= [] # a list of PE 
    SXX1= [] #sigma_xx
    SYY1= [] #sigma_yy
    SXY1= [] #sigma_xy
    P1= [] #hydrostatic pressure
    F1= [] #force
    with open(tx2,"r") as file:
        for line in file:
            words = line.split()
            n = float(words[0]) #take the value
            e = float(words[1]) #take the value
            sxx = float(words[2]) #take the value
            syy = float(words[3]) #take the value
            sxy = float(words[4]) #take the value
            p = float(words[5]) #take the value
            f = float(words[6]) #take the value

            N1.append(n); #append x value into X
            E1.append(e); #append y value into Y
            SXX1.append(sxx); #append sxx value into SXX
            SYY1.append(syy); #append syy value into SYY
            SXY1.append(sxy); #append sxy value into SXY
            P1.append(p);
            F1.append(f);

    #--------------------lambda = 1.5 ----------------------------
    N2= [] # a list of interation
    E2= [] # a list of PE 
    SXX2= [] #sigma_xx
    SYY2= [] #sigma_yy
    SXY2= [] #sigma_xy
    P2= [] #hydrostatic pressure
    F2= [] #force
    with open(tx3,"r") as file:
        for line in file:
            words = line.split()
            n = float(words[0]) #take the value
            e = float(words[1]) #take the value
            sxx = float(words[2]) #take the value
            syy = float(words[3]) #take the value
            sxy = float(words[4]) #take the value
            p = float(words[5]) #take the value
            f = float(words[6]) #take the value

            N2.append(n); #append x value into X
            E2.append(e); #append y value into Y
            SXX2.append(sxx); #append sxx value into SXX
            SYY2.append(syy); #append syy value into SYY
            SXY2.append(sxy); #append sxy value into SXY
            P2.append(p);
            F2.append(f);

    #--------------------lambda = 2 ------------------------------
    N3= [] # a list of interation
    E3= [] # a list of PE 
    SXX3= [] #sigma_xx
    SYY3= [] #sigma_yy
    SXY3= [] #sigma_xy
    P3= [] #hydrostatic pressure
    F3= [] #force
    with open(tx4,"r") as file:
        for line in file:
            words = line.split()
            n = float(words[0]) #take the value
            e = float(words[1]) #take the value
            sxx = float(words[2]) #take the value
            syy = float(words[3]) #take the value
            sxy = float(words[4]) #take the value
            p = float(words[5]) #take the value
            f = float(words[6]) #take the value

            N3.append(n); #append x value into X
            E3.append(e); #append y value into Y
            SXX3.append(sxx); #append sxx value into SXX
            SYY3.append(syy); #append syy value into SYY
            SXY3.append(sxy); #append sxy value into SXY
            P3.append(p);
            F3.append(f);

    #--------------------lambda = 8 ------------------------------
    N4= [] # a list of interation
    E4= [] # a list of PE 
    F4= [] #force
    with open(tx5,"r") as file:
        for line in file:
            words = line.split()
            n = float(words[0]) #take the value
            e = float(words[1]) #take the value
            p = float(words[5]) #take the value
            f = float(words[6]) #take the value

            N4.append(n); #append x value into X
            E4.append(e); #append y value into Y
            F4.append(f);

    #-------------------------------------------------------------
    #step 2
    #-------------------------------------------------------------

    #plot the plots

    #--------------------plot Potential energy vs N---------------

    fig = plt.figure()	

    plt.plot(N,E,'rx',linewidth = 0.3,markersize = 1,label = r'$\lambda = 0.5$')
    plt.plot(N,E1,'gx',linewidth = 0.3,markersize = 1,label = r'$\lambda = 1$')
    plt.plot(N,E2,'bx',linewidth = 0.3,markersize = 1,label = r'$\lambda = 1.5$')
    plt.plot(N,E3,'cx',linewidth = 0.3,markersize = 1,label = r'$\lambda = 2$')

    leg = plt.legend(prop={'size':10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'Iternations'
    ylabel = r'Potential Energy (eV)'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'Potential Energy vs N')
    fig.savefig(EN,dpi = 300, bbox_inches ='tight')


    #--------------------plot stresses vs N ----------------------
    # Sxx
    fig= plt.figure()
    plt.plot(N,SXX,'rx',linewidth = 0.3,markersize = 1,label = r'$\lambda = 0.5$')
    plt.plot(N,SXX1,'gx',linewidth = 0.3,markersize = 1,label = r'$\lambda = 1$')
    plt.plot(N,SXX2,'bx',linewidth = 0.3,markersize = 1,label = r'$\lambda = 1.5$')
    plt.plot(N,SXX3,'cx',linewidth = 0.3,markersize = 1,label = r'$\lambda = 2$')

    leg = plt.legend(prop={'size':10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'Iterations'
    ylabel = r'$\sigma_{xx}$ ($eV/\AA^{2}$)'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'$\sigma_{xx}$ vs N')
    fig.savefig(SxxN,dpi = 300, bbox_inches ='tight')

    # Syy
    fig= plt.figure()
    plt.plot(N,SYY,'rx',linewidth = 0.3,markersize = 1,label = r'$\lambda = 0.5$')
    plt.plot(N,SYY1,'gx',linewidth = 0.3,markersize = 1,label = r'$\lambda = 1$')
    plt.plot(N,SYY2,'bx',linewidth = 0.3,markersize = 1,label = r'$\lambda = 1.5$')
    plt.plot(N,SYY3,'cx',linewidth = 0.3,markersize = 1,label = r'$\lambda = 2$')

    leg = plt.legend(prop={'size':10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'Iterations'
    ylabel = r'$\sigma_{yy}$ ($eV/\AA^{2}$)'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'$\sigma_{yy}$ vs N')
    fig.savefig(SyyN,dpi = 300, bbox_inches ='tight')

    # Sxy
    fig= plt.figure()
    plt.plot(N,SXY,'rx',linewidth = 0.3,markersize = 1,label = r'$\lambda = 0.5$')
    plt.plot(N,SXY1,'gx',linewidth = 0.3,markersize = 1,label = r'$\lambda = 1$')
    plt.plot(N,SXY2,'bx',linewidth = 0.3,markersize = 1,label = r'$\lambda = 1.5$')
    plt.plot(N,SXY3,'cx',linewidth = 0.3,markersize = 1,label = r'$\lambda = 2$')

    leg = plt.legend(prop={'size':10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'Iterations'
    ylabel = r'$\sigma_{xy}$ ($eV/\AA^{2}$)'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'$\sigma_{xy}$ vs N')
    fig.savefig(SxyN,dpi = 300, bbox_inches ='tight')

    # P
    fig= plt.figure()
    plt.plot(N,P,'rx',linewidth = 0.3,markersize = 1,label = r'$\lambda = 0.5$')
    plt.plot(N,P1,'gx',linewidth = 0.3,markersize = 1,label = r'$\lambda = 1$')
    plt.plot(N,P2,'bx',linewidth = 0.3,markersize = 1,label = r'$\lambda = 1.5$')
    plt.plot(N,P3,'cx',linewidth = 0.3,markersize = 1,label = r'$\lambda = 2$')

    leg = plt.legend(prop={'size':10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'Iterations'
    ylabel = r'Hydrostatic pressure ($eV/\AA^{2}$)'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'Hydrostatic pressure vs N')
    fig.savefig(PN,dpi = 300, bbox_inches ='tight')


    #--------------------plot F_max vs N -------------------------
    fig= plt.figure()
    plt.plot(N,F,'r-',linewidth = 0.6,markersize = 1,label = r'$\lambda = 0.5$')
    plt.plot(N,F1,'g-',linewidth = 0.6,markersize = 1,label = r'$\lambda = 1$')
    plt.plot(N,F2,'b-',linewidth = 0.6,markersize = 1,label = r'$\lambda = 1.5$')
    plt.plot(N,F3,'c-',linewidth = 0.6,markersize = 1,label = r'$\lambda = 2$')

    leg = plt.legend(prop={'size':10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'Iterations'
    ylabel = r'Maximum Force ($eV/\AA$)'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'Maximum Force vs N')
    fig.savefig(FN,dpi = 300, bbox_inches ='tight')


    #--------------------plot F_8 vs N ---------------------------
    fig= plt.figure()
    plt.plot(N,F4,'r-',linewidth = 0.6,markersize = 1,label = r'$\lambda = 8$')
    leg = plt.legend(prop={'size':10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'Iterations'
    ylabel = r'Maximum Force ($eV/\AA$)'
    plt.ylim ([0,0.14])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'Maximum Force vs N')
    fig.savefig(F8,dpi = 300, bbox_inches ='tight')


    #--------------------plot E_8 vs N ---------------------------
    fig = plt.figure()	
    plt.plot(N,E4,'r-',linewidth = 0.6,markersize = 1,label = r'$\lambda = 8$')
    leg = plt.legend(prop={'size':10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'Iternations'
    ylabel = r'Potential Energy (eV)'
    plt.ylim ([-12,12])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'Potential Energy vs N')
    fig.savefig(E8,dpi = 300, bbox_inches ='tight')


def main():
    #-------------------------------------------------------------
    #plot the initial configuration and rdf 
    #-------------------------------------------------------------
    plotconfig("config.txt",'initial_config.png')
    plotrdf("rdf.txt","initial_rdf.png")

    #-------------------------------------------------------------
    #plot the final configuration and rdf of the 4 experiments
    #with different lambdas
    #-------------------------------------------------------------

    #----------------lambda = 0.5--------------------
    plotconfig("l=05_config.txt","l=05_config.png")
    plotrdf("l=05_rdf.txt","l=05_rdf.png")

    #----------------lambda = 1 ---------------------
    plotconfig("l=1_config.txt","l=1_config.png")
    plotrdf("l=1_rdf.txt","l=1_rdf.png")

    #----------------lambda = 1.5--------------------
    plotconfig("l=15_config.txt","l=15_config.png")
    plotrdf("l=15_rdf.txt","l=15_rdf.png")

    #----------------lambda = 2 ---------------------
    plotconfig("l=2_config.txt","l=2_config.png")
    plotrdf("l=2_rdf.txt","l=2_rdf.png")

    #----------------lambda = 8 ---------------------
    plotconfig("l=8_config.txt","l=8_config.png")
    plotrdf("l=8_rdf.txt","l=8_rdf.png")

    #------------------------------------------------
    #study the evolution lambda = 1.5
    #------------------------------------------------

    #iteration 250
    plotconfig("250config.txt","250config.png")
    plotrdf("250rdf.txt","250rdf.png")

    #iteration 500
    plotconfig("500config.txt","500config.png")
    plotrdf("500rdf.txt","500rdf.png")

    #iteration 1000
    plotconfig("1000config.txt","1000config.png")
    plotrdf("1000rdf.txt","1000rdf.png")

    #iteration 2000
    plotconfig("2000config.txt","2000config.png")
    plotrdf("2000rdf.txt","2000rdf.png")

    #------------------------------------------------
    #plot variables vs N
    #------------------------------------------------
    plot_vsN("l=0.5_vsN.txt","l=1_vsN.txt","l=1.5_vsN.txt",
            "l=2_vsN.txt","l=8_vsN.txt","E_vs_N.png",
            "Sxx_vs_N.png","Syy_vs_N.png","Sxy_vs_N.png",
            "P_vs_N.png","F_vs_N.png","P8_vs_N.png","F8_vs_N.png")

main()
