import itertools
import math as m
import numpy as np
from scipy import special as s #special.erfi
import matplotlib.pyplot as plt
import cmath as c
import time

def S2(N):
    NN=N**2.0
    XX1=[]
    XX2=[]
    XX3=[]
    X1=range(-N,N+1) 
    for x1 in X1:
        for x2 in X1:
            for x3 in X1:
                xx1=x1**2.0
                xx2=x2**2.0
                xx3=x3**2.0
                if (xx1+xx2+xx3)<=NN:
                    XX1.append(x1)
                    XX2.append(x2)
                    XX3.append(x3)
##    (XX1,XX2,XX3)=(array(XX1),array(XX2),array(XX3))
    return (np.array(XX1),np.array(XX2),np.array(XX3))

t1 = time.time()
(x1,x2,x3)=S2(15)
t2 = time.time()

print("it took this long to generate",t2-t1)

print("number of terms in the sum = ",len(x1))



def Zcal(mLotwopi,EcmLotwopi, dvec, alpha, pres, J,mJ, N):
    omegaq = EcmLotwopi / 2.0
    
    xx=m.pow(omegaq,2)-m.pow(mLotwopi,2)

    rx,ry,rz=x1,x2,x3

    rr = rx*rx + ry*ry + rz*rz

    Etot=c.sqrt(m.pow(EcmLotwopi,2)+np.dot(dvec,dvec))
    
    omegk=np.sqrt(m.pow(mLotwopi,2) + rr)
    
    PP=np.dot(dvec,dvec)

    beta_vec = dvec/Etot
    
    gg = Etot/EcmLotwopi
    
    
    if PP not in [0]:
        amppar=(rx*dvec[0] + ry*dvec[1] + rz*dvec[2])/PP
    else:
        amppar=0

    rx_par,ry_par,rz_par=amppar*dvec[0],amppar*dvec[1],amppar*dvec[2] 
    rx_per,ry_per,rz_per=rx-rx_par,ry-ry_par,rz-rz_par
        

    rx_star = gg*(rx_par - beta_vec[0] * omegk) + rx_per
    ry_star = gg*(ry_par - beta_vec[1] * omegk) + ry_per
    rz_star = gg*(rz_par - beta_vec[2] * omegk) + rz_per

    rrstar = pow(rx_star,2)+pow(ry_star,2)+pow(rz_star,2)
    omegkstar=np.sqrt(m.pow(mLotwopi,2) + rrstar)

    ampS = omegkstar / omegk
        
    HCut=np.exp( - alpha*pow(rrstar-xx,N))

    ampS = ampS * HCut

    Ytensor = 1.0 #c.sqrt(4.0*pi) * Ylmr(rx_star,ry_star,rz_star,J,mJ)
        
    ampS = ampS * Ytensor
        
    denum =  pow(xx - rrstar ,N)
        
    S0 = sum(ampS / denum)
         
        
    #Int0 = XiN(xx,alpha, pres, N)
            

    tmp = S0 - 0.0 #Int0

    return tmp.real, tmp.imag



m1 = 4.0 / ( 2.0 * np.pi )
m2 = 4.0 / ( 2.0 * np.pi )
nP = np.array([0,1,1])
m1_sq = m1**2
m2_sq = m2**2

Ecm_o_m_start = 1.8
Ecm_o_m_stop = 4.0
Ecm_o_m_step = 0.05

alpha = 0.1
pres = 0

t1 = time.time()
for Ecm in np.arange (Ecm_o_m_start, Ecm_o_m_stop, Ecm_o_m_step):
   Ecm = Ecm * m1
   zeta = Zcal(m1,Ecm, nP, alpha, pres, 0,0, 15)
   print( Ecm / m1, " ", zeta )
t2 = time.time()

print("it took this long to calculate",t2-t1)
