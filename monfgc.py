# -*- coding: utf-8 -*-
"""
Created on Sat Jun 22 18:32:44 2019

@author: Nour Hameed
"""

import matplotlib.pyplot as plt
import time as time
import func_mc
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

def monfgc(H,nfot,d,opt,op2, Rmin, rf, geoh, E0, ZA, N, nam, mat ,pmol, faZ, b0, b1, b2, NS, Emin,NT,NR,NM,NA,R):
    
	# monfgc #
    print("monfgc running")
	# Graphic window
    
    fig_fgp = plt.figure("Tracing the rays")
    ax = plt.axes(projection="3d")
    #ax = plt.axes()
    ax.view_init(0.0,0.0)
    ax.set_title("CYLINDRICAL GEOMETRY")
    ax.set_xlabel("Distance X (cm)")
    ax.set_ylabel("Distance Y (cm)")
    ax.set_zlabel("Depth (cm)")
    if op2=="t":
        time.sleep(2)
    nfot1 = int(nfot)
    for i in range(nfot1):
        theta1, phi1, x0, y0, z0 = func_mc.incid(d, Rmin, rf, geoh)
        E1=E0; E=E0; Efot=0; veces=0;
        # Point of the first interaction   
        s = func_mc.rlm(E,E1,ZA,N,nam,mat,pmol,faZ,b0,b1,b2)
        x=s*np.sin(theta1)*np.cos(phi1)+x0
        y=s*np.sin(theta1)*np.sin(phi1)+y0
        z=s*np.cos(theta1)+z0
        RA=np.sqrt(x**2+y**2)
        if op2=='t':
            ax.plot3D([x0,x],[y0,y],[z0,z],'-b');
        if (z>H)|(RA>R):
            choice4=0
        else:
            choice4=1
	   
        while choice4:
            # Choice of the process      
            control = func_mc.select_process(E,ZA,pmol,nam,mat,E1,faZ,b0,b1,b2)
        	# Photoelectric effect     
      
            if (control=='Foto'):
                Eper=1.0
                NA=NA+1
                break
		 
            if E<Emin:
                Eper=0.0
                NA=NA+1
                break
		 
	# Compton effect      
            if control=='Comp':
                veces=veces+1
                theta = func_mc.angtheta(E)
                E1 = func_mc.energy_loss(E, theta)
                Eper=1.0-(E0-E1)/E0
                phi = func_mc.angphi()
                s = func_mc.rlm(E,E1,ZA,N,nam,mat,pmol,faZ,b0,b1,b2)
                phi2, theta2 = func_mc.angles(E,theta,theta1,phi,phi1)
                x1,y1,z1,RA1 = func_mc.coord(s,theta2,phi2,x,y,z)
                if op2=='t':
                    ax.plot3D([x,x1],[y,y1],[z,z1],color='blue')
                else:
                    col=round(5.0*Eper)
                    if col==0:
                        ax.scatter3D(x,y,z,color='yellow')
                    elif col==1:
                        ax.scatter3D(x,y,z,color='yellow')
                    elif col==2:
                        ax.scatter3D(x,y,z,color='green')
                    elif col==3:
                        ax.scatter3D(x,y,z,color='green')
                    elif col==4:
                        ax.scatter3D(x,y,z,color='red')
                    else:
                        ax.scatter3D(x,y,z,color='red')
			  
                if (z1<0)|(z1>H)|(RA1>R):
                    break
			  
		    # control
            theta1=theta2; phi1=phi2; x=x1; y=y1; z=z1; RA=RA1; E=E1;     
	        # choice4
        if veces>1:
            NM=NM+1
	   
        if op2=='t':
    	    time.sleep(2) #
	   
    return NT,NS,NR,NA,nfot  
