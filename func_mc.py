# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 21:03:58 2019

@author: Nour Hameed

"""
import math
import random
import numpy as np
import matplotlib as plt

def angphi():  # Angle phi: 'angphi'
    rphi=random.random()
    phi=2*np.pi*rphi
    return phi

######################################################
def angtheta(E):
    G=E/0.511
    cont=1
    while cont:
        rx=random.random() 
        ry=random.random()
        rth=np.pi*rx 
        rKNc=2.75e-25*ry
        a1=np.sin(rth)/(1+G*(1-np.cos(rth)))**2
        a2=((G**2)*(1-np.cos(rth))**2)/(1+G*(1-np.cos(rth)))
        rKN=2.494672e-25*a1*(1+(np.cos(rth)**2)+a2)
        if rKN>rKNc:
            theta=rth
            cont=0
    return theta

######################################################
def angles(E,theta,theta1,phi,phi1):
    # Ángles SRL: 'angles'
    st=np.sin(theta)
    ct=np.cos(theta)
    st1=np.sin(theta1)
    ct1=np.cos(theta1)
    cp=np.cos(phi)
    sp=np.sin(phi)
    ct2=ct1*ct+st1*st*cp
    theta2=math.acos(ct2)
    st2=np.sin(theta2)
    cdphi=(ct-ct2*ct1)/(st1*st2)
    if cdphi>1:
        cdphi=1
    elif cdphi<-1:
        cdphi=-1
    sdphi=st*sp/st2
    if sdphi>0:
        phi2=math.acos(cdphi)+phi1
    else:
        phi2=2*np.pi-math.acos(cdphi)+phi1
    phi2=phi2%(2*np.pi)
    return phi2, theta2

######################################################
def coord(s,theta2,phi2,x,y,z):
    x1=s*np.sin(theta2)*np.cos(phi2)+x
    y1=s*np.sin(theta2)*np.sin(phi2)+y
    z1=s*np.cos(theta2)+z
    RA1=np.sqrt(x1**2+y1**2)
    return x1,y1,z1,RA1

#########################################################################
def dist_exp():
    N1=10000 #for the histogram
    np.histogram(-10*np.log(np.random.rand(N1)),20)
    np.figure()
    N2=100 #for trajectory
    for i in range(N2):
        x1=-10*np.log(random.random());
        x0=0;
        y0=(i+1)/N2; y1=y0;
        plt.plot([x0,x1],[y0,y1])

#####################################################
def Eabs(E0,Emin,Efot,ME,Edep):
    rango=E0-Emin
    j=math.ceil(ME*(Efot-Emin)/rango)-1
    if j>-1:
       if j==ME:
          j=ME-1
       Edep[0][j]=Edep[0][j]+1 
    return Edep

######################################################
def Ecapac(E,E1,Mz,z,H,capa,P,RA,R):
    # Energy - Cylinder: 'Ecapac'
    depox=E-E1;
    j=math.ceil(Mz*z/H); 
    r=math.ceil(P*RA/R);
    if r==0:
       r=1
    capa[j-1][r-1]=capa[j-1][r-1]+depox
    return capa

######################################################
def select_process(E,ZA,pmol,nam,mat,E1,faZ,b0,b1,b2):
    # Choice of process: 'select_process'
    G=E/0.511;
    sigmac0=((1+G)/(G**2)*(2*(1+G)/(1+2*G)-(np.log(1+2*G)/G))+np.log(1+2*G)/(2*G)-(1+3*G)/((1+2*G)**2))*4.989344e-25
    sigmac=np.multiply(sigmac0,ZA)
    sigmacm=sigmac.dot(np.transpose(nam))
    #
    if mat == 'NaI': # Change this 
        if E<=33.273e-3:
            sigmaefm=(pmol/6.022e23)*2.11831202e-4*(E1**(-2.892727))
        else:
            sigmaefm=(pmol/6.022e23)*2.35338947e-3*(E1**(-2.7778746))
    else:
        ZA2 = np.power(ZA,5)
        ZA3 = np.divide(faZ,G)
        ZA4 = np.multiply(ZA2,ZA3)
        ZA5 = b0+np.divide(b1,G)+np.divide(b2,np.power(G,2))
        ZA6 = np.multiply(ZA4,ZA5)
        sigmaef = ZA6*(2.829663e-33)
        sigmaefm=sigmaef.dot(np.transpose(nam))
    #
    C=sigmacm/sigmaefm
    comp=C/(C+1)
    Rp=random.random()
    if Rp<comp:
       control='Comp'
    else:
       control='Foto'
    return control

######################################################
def Erefle(E0,E1,ME,Er):
    # Energy distribution of the reflected photons : 'Erefle'
    Emax=E0 
    Emin=0.01 
    rango=Emax-Emin
    j=math.ceil(ME*(E1-Emin)/rango)-1
    Er[0][j]=Er[0][j]+1
    return Er
    
######################################################
def incid(d, Rmin, rf, geoh):  # interaction / impact
    costheta_max = d/np.sqrt(d**2+(Rmin+rf)**2);
    follow0=1
    while follow0:
        follow1=1
        while follow1:
            xf=2*rf*random.random()-rf; yf=2*rf*random.random()-rf;
            if xf**2+yf**2<=rf**2:
                follow1=0
    #
        x0=xf
        y0=yf 
        theta1=1e-10
        phi1=0
        follow2=1
        while follow2:
            follow3=1
            if geoh[0]=="i":
                 while follow3:
                     theta1=np.acos(1-random.random()*(1-costheta_max))
                     follow3=0
                 phi1=2*np.pi*random.random()
                 x0=d*np.tan(theta1)*np.cos(phi1)+xf
                 y0=d*np.tan(theta1)*np.sin(phi1)+yf
            RA0=np.sqrt(x0**2+y0**2)
            
            if RA0 <= Rmin:
                follow2=0
                follow0=0
            else:
                follow2=0
    z0=0
    return theta1, phi1, x0, y0, z0

######################################################
def Etrans(E0,E1,ME,Et):
    # Energy distribution of transmitted photons: 'Etrans'
    Emax=E0
    Emin=0.01
    rango=Emax-Emin
    j=np.ceil(ME*(E1-Emin)/rango)
    j = int(j)-1
    Et[0][j]=Et[0][j]+1
    return Et
    
######################################################
def energy_loss(E,theta):
    # Energy Loss: 'energy_loss'
    G=E/0.511;
    E1=E/(1+G*(1-np.cos(theta)))
    return E1
    
######################################################
def rlm(E,E1,ZA,N,nam,mat,pmol,faZ,b0,b1,b2):
    # Average - movement without interactions: 'rlm' 
    G=E1/0.511
    sigmac0=((1+G)/(G**2)*(2*(1+G)/(1+2*G)-(np.log(1+2*G)/G))+np.log(1+2*G)/(2*G)-(1+3*G)/((1+2*G)**2))*4.989344e-25
    sigmac=np.multiply(sigmac0,ZA)
    sigmacm=sigmac.dot(np.transpose(nam))
    #
    if mat == "NaI":
        if E<=33.273e-3:
            sigmaefm=(pmol/6.022e23)*2.11831202e-4*(E1**(-2.892727))
        else:
            sigmaefm=(pmol/6.022e23)*2.35338947e-3*(E1**(-2.7778746))
    else:
        ZA2 = np.power(ZA,5)
        ZA3 = np.divide(faZ,G)
        ZA4 = np.multiply(ZA2,ZA3)
        ZA5 = b0+np.divide(b1,G)+np.divide(b2,np.power(G,2))
        ZA6 = np.multiply(ZA4,ZA5)
        sigmaef = ZA6*(2.829663e-33)
        sigmaefm=sigmaef.dot(np.transpose(nam))
    seftot=sigmacm+sigmaefm
    chi=random.random()
    s=-np.log(chi)/(seftot*N)
    return s
    
####################################################
def threfle(B,theta2,numthr):
    # Angular distribution of the reflected photons: 'threfle'
    j=math.ceil(B*(theta2-np.pi/2)/(np.pi/2))-1
    numthr[0][j]=numthr[0][j]+1
    return numthr
    
####################################################
def thtr(B,theta2,numtht):
    # Angular distribution of transmitted photons: 'thtr'
    j=math.ceil(B*theta2/(np.pi/2))-1
    numtht[0][j]=numtht[0][j]+1
    return numtht
