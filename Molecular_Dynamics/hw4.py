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
    plt.clim(-0.01,0.01)
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
    plt.clim(-0.01,0.01)
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
    plt.clim(-0.01,0.01)
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
    plt.clim(-0.01,0.01)
    xlabel = r'X_{pos}'
    ylabel = r'Y_{pos}'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.locator_params(axis='y', nticks=3)
    plt.locator_params(axis='x', nticks=8)
    fname = "hydrop" + filename 
    titlename = r"Hydrostatic pressure"
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

def E_vs_it():
    KE = [] # a list of KE 
    PE = [] # a list of PE 
    IT = [] # Iteration
    with open("MD_Energy.txt","r") as file:
        for line in file:
            words = line.split()
            ke = float(words[0]) #take the value
            pe = float(words[1]) #take the value
            it = float(words[2]) #take the value

            KE.append(ke); 
            PE.append(pe); 
            IT.append(it); 

    fig = plt.figure()	
    plt.plot(IT,KE,'r+',linewidth = 0.1,markersize = 1,label = r'KE')
    leg = plt.legend(prop={'size':10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'Iternations'
    ylabel = r'Kinetic Energy (eV)'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'Kinetic Energy vs N')
    fig.savefig("KEvsIT",dpi = 300, bbox_inches ='tight')


    fig = plt.figure()	
    plt.plot(IT,PE,'b+',linewidth = 0.1,markersize = 1,label = r'PE')
    leg = plt.legend(prop={'size':10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'Iternations'
    ylabel = r'Potential Energy (eV)'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'Potential Energy vs N')
    fig.savefig("PEvsIT",dpi = 300, bbox_inches ='tight')

def plottemp_pressure(filename):
    T = [] # a list of KE 
    P = [] # a list of PE 
    step = []
    s = 0
    with open(filename,"r") as file:
        for line in file:
            words = line.split()
            t = float(words[0]) #take the value
            p = float(words[4]) #take the value
            s= s+ 1

            T.append(t); 
            P.append(p); 
            step.append(s)
        

    fig = plt.figure()	
    plt.plot(step,T,'gx-',linewidth = 0.5,markersize = 8,label = r'Temperature')
    leg = plt.legend(prop={'size':10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'$N^{th}$ MD simulation'
    ylabel = r'Temperature (K)'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'Temperature vs MD')
    fig.savefig("temperature",dpi = 300, bbox_inches ='tight')

    fig = plt.figure()
    plt.plot(step,P,'cx-',linewidth = 0.5,markersize = 8,label = r'Pressure')
    leg = plt.legend(prop={'size':10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'$N^{th}$ MD simulation'
    ylabel = r'Pressure ($eV/\AA^2$)'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'Pressure vs MD')
    fig.savefig("Pressure",dpi = 300, bbox_inches ='tight')

def totE_VS_T():
    #data collect from each MD simulations
    #  KE = [0.181131,0.331347,0.548679,0.705396,0.924459,1.09682,1.18938,
           #  1.30866, 1.64581, 1.62688, 1.90752, 2.07408, 2.19489, 2.37099 ,
           #  2.58581,2.78542  ] # a list of KE
    #  PE = [-12.1695,-11.9623,-11.8072,-11.5919,-11.4187,-11.1689,-10.7799,
           #  -10.2842, -9.58792,-8.21997,-7.37905,-6.7233,-6.41896,-6.10916,
           #  -5.80244,-5.79847  ] # a list of PE
    KE_n = [0.170715,0.346807,0.528555,0.670939,0.867841,1.08324,1.12488,
	   1.43006,1.48821,1.60668,1.79389,2.05983,2.22374,2.52244,2.67557,2.73885]
    PE_n = [-12.1603,-11.9767,-11.7917,-11.5596,-11.3819,-11.1482,-10.6984,-10.3408 ,-9.49048,-8.22616,
	    -7.86945,-7.5226,-7.2282,-7.05292,-6.93157,-6.62429]

    P = [-0.0264,-0.05148,-0.07551,-0.10,-0.1300,-0.1519,-0.1814,-0.20616,
           -0.2266,-0.23037,-0.22961,-0.2160,-0.214,-0.212,-0.194,-0.185]
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
    plt.ylim ([0,4])
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
    plt.ylim ([-13,-5])
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
    plt.ylim ([-13,-2])
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

def autocor():
    T = [] 
    Auto = [] 
    with open("auto.txt","r") as file:
        for line in file:
            words = line.split()
            t = float(words[1]) #take the value
            a = float(words[0]) #take the value

            T.append(t); 
            Auto.append(a); 

    fig = plt.figure()	
    plt.plot(T,Auto,'c-',linewidth = 1,markersize = 5,label = r'Cp')
    leg = plt.legend(prop={'size':10})
    leg.get_frame().set_alpha(0.5)
    xlabel = r'Temperature (K)'
    ylabel = r'Heat capacity (eV/K)'
    #  plt.ylim ([-13,-2])
#    plt.xlim ([0,30])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(r'Heat capacity vs T')
    fig.savefig("Auto",dpi = 300, bbox_inches ='tight')



def main():
    totE_VS_T()

    #  autocor()
    #  E_vs_it()
    #  plotconfig("5K_MD_minimization.txt","config.png")
    #  plotrdf("5K_rdf_MDmini.txt","rdf.png")
    #  plottemp_pressure("part2_MD_temp_stress.txt")


main()
