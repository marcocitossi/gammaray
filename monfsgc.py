import numpy as np
import func_mc
from mpl_toolkits.mplot3d import Axes3D

def monfsgp(H,nfot,d,opt,op2, Rmin, rf, geoh, E0, ZA, N, nam, mat ,pmol, faZ, b0, b1, b2, NS, Emin,NT,NR,NM,NA,R):
    print("monfsgp running")
    
    for i in range(nfot):
        theta1, phi1, x0, y0, z0 = func_mc.incid(d, Rmin, rf, geoh)
        E1=E0; E=E0; Efot=0; veces=0;
        
        # Point of first interaction    
        s = func_mc.rlm(E,E1,ZA,N,nam,mat,pmol,faZ,b0,b1,b2);
        x=s*np.sin(theta1)*np.cos(phi1)+x0;
        y=s*np.sin(theta1)*np.sin(phi1)+y0;
        z=s*np.cos(theta1)+z0;
        RA=np.sqrt(x^2+y^2);
        if (z>H)|(RA>R):
            choice4=0
        else:
            choice4=1
       
        while choice4:
            # Choosing the process     
            control = func_mc.select_process()
            # Photoelectric effect      
            if (control=="Foto")|(E<Emin):
                E1=0
                capa = func_mc.Ecapac(E,E1,Mz,z,H,capa,P,RA,R)
                Efot=Efot+(E-E1)
                Edep = func_mc.Eabs(E0,Emin,Efot,ME,Edep)
                NA=NA+1
                break
          
            # Compton effect   
            if control=="Comp":
                veces=veces+1;
                theta = func_mc.angtheta()
                E1 = func_mc.energy_loss(E,theta)
                Efot=Efot+(E-E1);
                phi = func_mc.angphi()
                s = func_mc.rlm(E,E1,ZA,N,nam,mat,pmol,faZ,b0,b1,b2)
                phi2, theta2 = func_mc.angles(E,theta,theta1,phi,phi1)
                x1,y1,z1,RA1 = func_mc.coord(s,theta2,phi2,x,y,z)
                capa = func_mc.Ecapac(E,E1,Mz,z,H,capa,P,RA,R)
                if (z1<0 or z1>H or RA1>R):
                    Edep = func_mc.Eabs(E0,Emin,Efot,ME,Edep)
                    break
                        
            # control
            theta1=theta2; phi1=phi2; x=x1; y=y1; z=z1; RA=RA1; E=E1;     
         # choice4
        if veces>1:
            NM=NM+1;
    return NT,NS,NR,NA,nfot
