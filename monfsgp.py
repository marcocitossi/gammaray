# -*- coding: utf-8 -*-
"""
Created on Sat Jun 22 20:43:24 2019

@author: Nour Hameed
"""
import numpy as np
import func_mc
from mpl_toolkits.mplot3d import Axes3D

def monfsgp(H,nfot,d,opt,op2, Rmin, rf, geoh, E0, ZA, N, nam, mat ,pmol, faZ, b0, b1, b2, NS, Emin,NT,NR,NM,NA,ME,Edep,Et,B,numtht, Er, numthr):
    print("monfsgp running")

    # Program - No figure
    nfot1 = int(nfot)
    for i in range(nfot1):
        theta1, phi1, x0, y0, z0 = func_mc.incid(d, Rmin, rf, geoh)
        E1=E0; E=E0; Efot=0; veces=0;
    
        # Graphic window   
        s=func_mc.rlm(E,E1,ZA,N,nam,mat,pmol,faZ,b0,b1,b2); 
        x=s*np.sin(theta1)*np.cos(phi1)+x0;
        y=s*np.sin(theta1)*np.sin(phi1)+y0;
        z=s*np.cos(theta1)+z0;
        if z>H:
            NS=NS+1; 
            choice=0;
        else:
            choice=1;
       
        while choice:
            # Choosing the process      
            control = func_mc.select_process(E,ZA,pmol,nam,mat,E1,faZ,b0,b1,b2)
            # Photoelectric effect      
            if (control=='Foto')|(E<Emin):
               E1=0;
            # Ecapap;
               Efot=Efot+(E-E1)
               Edep = func_mc.Eabs(E0,Emin,Efot,ME,Edep)
               NA=NA+1;
               break
          
            # Compton effect      
            if control=='Comp':
               veces=veces+1
               theta = func_mc.angtheta(E)
               E1 = func_mc.energy_loss(E,theta)
               Efot=Efot+(E-E1)
               phi = func_mc.angphi()
               s = func_mc.rlm(E,E1,ZA,N,nam,mat,pmol,faZ,b0,b1,b2)
               phi2, theta2 = func_mc.angles(E,theta,theta1,phi,phi1)
               x1,y1,z1,RA1 = func_mc.coord(s,theta2,phi2,x,y,z)

            if z1>H:
               Et = func_mc.Etrans(E0,E1,ME,Et)
               numtht = func_mc.thtr(B,theta2,numtht)
               Edep = func_mc.Eabs(E0,Emin,Efot,ME,Edep)
               NT=NT+1;
               break
             
            if z1<0:
               NR=NR+1;
               Er = func_mc.Erefle(E0,E1,ME,Er)
               numthr = func_mc.threfle(B,theta2,numthr)
               Edep = func_mc.Eabs(E0,Emin,Efot,ME,Edep)
               break
             
            # control
            theta1=theta2; phi1=phi2; x=x1; y=y1; z=z1; E=E1;     
         # choice
        if veces>1:
            NM=NM+1
    return NT,NS,NR,NA,nfot
