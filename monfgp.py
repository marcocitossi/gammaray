# -*- coding: utf-8 -*-
"""
Created on Sat Jun 22 18:31:05 2019

@author: Nour Hameed
"""
import matplotlib.pyplot as plt
import numpy as np
import time as time
import func_mc
from mpl_toolkits.mplot3d import Axes3D

def monfgp(H,nfot,d,opt,op2, Rmin, rf, geoh, E0, ZA, N, nam, mat ,pmol, faZ, b0, b1, b2, NS, Emin,NT,NR,NM,NA):
    print("monfgp running")
    nfot1 = int(nfot)  
    fig_fgp = plt.figure("Tracing the rays",figsize=plt.figaspect(0.5)*1.5)

	#Graphic window
    ax = plt.axes(projection="3d")
    ax.set_xlim(-2*H,2*H)
    ax.set_ylim(-2*H,2*H)
    ax.set_zlim(-0.5*H,1.5*H)
    ax.plot3D([-2*H,-2*H,2*H,2*H,-2*H], [-2*H,2*H,2*H,-2*H,-2*H], [H,H,H,H,H], '-r')
    ax.plot3D([-2*H ,-2*H, 2*H, 2*H, -2*H], [-2*H, 2*H, 2*H, -2*H, -2*H], [0,0,0,0,0], '-r')
    ax.view_init(0.0,90.0)
    ax.set_title("PLANE GEOMETRY")
    ax.set_xlabel("Distance X (cm)")
    ax.set_ylabel("Distance Y (cm)")
    ax.set_zlabel("Depth (cm)")
    
    # Hide grid lines
    ax.grid(False)
    
    # Hide axes ticks
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
        
    if op2=="t":
        time.sleep(5)
    for i in range(nfot1):
        theta1, phi1, x0, y0, z0 = func_mc.incid(d, Rmin, rf, geoh)
        E1=E0; E=E0; Efot=0; veces=0;

    	# Point of the first interaction    
        s = func_mc.rlm(E,E1,ZA,N,nam,mat,pmol,faZ,b0,b1,b2)
        x=s*np.sin(theta1)*np.cos(phi1)+x0
        y=s*np.sin(theta1)*np.sin(phi1)+y0
        z=s*np.cos(theta1)+z0
    	   
        zt = min(d,0.4*H)
        x1 = x0 - (zt/np.cos(theta1))*np.sin(theta1)*np.cos(phi1)
        y1 = y0 - (zt/np.cos(theta1))*np.sin(theta1)*np.cos(phi1)
        z1 = z0 - zt
        #
        if op2=='t':
            ax.plot3D([x0,x],[y0,y],[z0,z],"blue")
            ax.plot3D([x1,x0],[y1,y0],[z1,z0],"blue") 
       
        if s>H:
            NS=NS+1
            choice=0
        else:
            choice=1
#####################################################
        while choice==1:
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
           
            # Efecto Compton
            if control=="Comp":
                veces=veces+1
                theta=func_mc.angtheta(E)
                E1 = func_mc.energy_loss(E,theta)
                Eper=1.0-(E0-E1)/E0
                phi = func_mc.angphi()
                s=func_mc.rlm(E,E1,ZA,N,nam,mat,pmol,faZ,b0,b1,b2)
                phi2,theta2=func_mc.angles(E,theta,theta1,phi,phi1)
                x1,y1,z1,RA1=func_mc.coord(s,theta2,phi2,x,y,z)
                if op2=='t':
                    ax.plot3D([x,x1],[y,y1],[z,z1],"blue")
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
    		   
                if z1>H:
                    NT=NT+1
                    x2=s*np.sin(theta2)*np.cos(phi2)+x1
                    y2=s*np.sin(theta2)*np.sin(phi2)+y1
                    z2=s*np.cos(theta2)+z1
                    if opt == 't':
                        ax.plot3D([x1, x2],[y1, y2],[z1, z2], "blue")
                        time.sleep(5)
                        break
                if z1<0:
                    NR=NR+1
                    break
                # control
                theta1=theta2; phi1=phi2; x=x1; y=y1; z=z1; E=E1;
        # choice
        if veces>1:
            NM=NM+1
        
        if op2=='t':
            time.sleep(2) 	   

    return NT,NS,NR,NA,nfot
