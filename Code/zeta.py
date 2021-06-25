#! /usr/bin/env python

import itertools
import math as m
import numpy as np
from scipy import special as s #special.erfi
import matplotlib.pyplot as plt
import cmath as c
import time
from numpy.core.umath_tests import inner1d
from scipy import optimize


'''
 sqrt with branch cut aligned on positive real axis
'''
def mySqrt( z ):
   return 1j * c.sqrt( -z )

'''
 inner product for two list (native)
'''
def dot(A, B):
   return np.dot(A,B) #inner1d(A,B) #np.sum( A * B, axis=1) 

'''
 Function to generate 3-tuple "n = (nx,ny,nz)" list 
 subject to the constraint nx^2 + ny^2 + nz^2 <= nShell^2,
 where "nShell" indicates the number of shells
'''
def generateIntList( nShell ):
    intList = []
    rawList = np.array(np.arange(-nShell,nShell+1))
    for i in itertools.product(rawList,rawList,rawList):
        tmpList = np.array(i)
        if ( np.dot(tmpList,tmpList) <= nShell*nShell ):
            intList.append(tmpList)

    return np.asarray(intList)

#n = [ [a1,a2,a3], [b1,b2,b3], [c1,c2,c3], ... ]

## inner1d(n,n) = [ a.dot(a), b.dot(b), ... ]

## n.dot(unit(beta)) = [ a.dot(unit(beta)), b.dot(unit(beta)), ... ]

'''
 returns the unit vector of a vector "A"
'''
def unit_vector( A ):
    if( all( A==0 ) ):
        return 0.0
    else:
        return A / np.sqrt( np.dot(A,A) )

'''
 special multiply for array of arrays (A) and an array (B), used in the zeta function
'''
def multiply( A, B ):
   nA = len(A)
   return A.reshape(nA,1) * np.tile( B, (nA,1) )
   ### First, reshape A array to 2d (1xnA), and then copy b vector nA time to make (nA x nB) array
   ### alternatively, do (A.reshape(nA,1)).dot( B.reshape(1,3) ) 
   ### bad alternative:  np.full( (nB,nA), A ).T * np.full( (nA,nB) , B )


'''
 returns a Lorentz boosted four vector "(A0,A)" given some boost
 velocity "beta"
'''
def lorentzBoost( A0, A, beta ):
    gamma = 1.0 / np.sqrt( 1.0 - np.dot(beta,beta) )
    Apar  = multiply( A.dot( unit_vector(beta) ), unit_vector(beta) )
    Aper  = A - Apar
    A0_prime = gamma * ( A0 - A.dot( beta ) )
    A_prime  = gamma * ( Apar - multiply( A0, beta ) ) + Aper
    return A0_prime, A_prime

'''
 KSS cutoff function
'''
def cutoff( x_sq, n_sq, alpha ):
    return np.exp( -alpha * ( n_sq - x_sq ) )

'''
 KSS Zeta function for ell = m_ell = 0, given "x_sq", total momentum "nP"
 This function is in units of L / 2pi, e.g. m1 -> m1 * L / (2 * pi).
'''
def zetaFunction( x_sq, nP, alpha, nShell ):
    Ecm = np.sqrt( m1_sq + x_sq ) + np.sqrt( m2_sq + x_sq )
    beta = nP / np.sqrt( Ecm * Ecm + np.dot(nP,nP) )    

    n = intList
    n0 = np.sqrt( m1_sq + inner1d(n,n) )
    n0_cm, n_cm = lorentzBoost( n0, n, beta )
    tmp = cutoff( x_sq, inner1d(n_cm,n_cm), alpha ) * n0_cm / n0 
    tmp = tmp / ( inner1d(n_cm,n_cm ) - x_sq ) 
    sum0 = np.sum( tmp )

    eps = 1.0e-16 * ( 1.0 - m.copysign(1, x_sq) )
    x   = c.sqrt( x_sq + 1j * eps )
    z   = c.sqrt( alpha * x_sq )
    integral = 4.0 * m.pi * ( ( m.sqrt( 0.25 * m.pi / alpha ) ) * m.exp( alpha * x_sq )
                              - 0.5 * m.pi * x * s.erfi( z ) )
        
    return sum0 - integral - 4.0 * m.pi * ( 0.5 * 1j * m.pi * x ) # last term needed to compare to BBHO


'''
 convert Ecm to x_sq 
'''
def x_sq( Ecm ):
    return ( Ecm * Ecm - ( m1 + m2 )**2 ) * ( Ecm * Ecm - ( m1 - m2 )**2 ) / ( 4.0 * Ecm * Ecm )

'''
 running code
'''
alpha = .1 #( 1.0 / 3.0 )**4
nShell = 20

t1 = time.time()
intList = generateIntList( nShell )
t2 = time.time()
#print("it took this long to generate",t2-t1)

#print( lorentzBoost(v0,v,b) )

# m = m * L / ( 2*pi )
m1 = 4.0 / ( 2.0 * np.pi )
m2 = 4.0 / ( 2.0 * np.pi )
nP = np.array([0,0,1])
m1_sq = m1**2
m2_sq = m2**2

Ecm_o_m_start = 1.8
Ecm_o_m_stop = 4.0
Ecm_o_m_step = 0.01

#Ecm = 3.0 * m1
#zeta = zetaFunction( x_sq(Ecm), nP, alpha, nShell )
#print( Ecm / m1, " ", x_sq(Ecm), " ", zeta.real, " ", zeta.imag )

EcmRange = np.arange(Ecm_o_m_start, Ecm_o_m_stop, Ecm_o_m_step)
zetaData_real = []
zetaData_imag = []

t1 = time.time()
for Ecm in EcmRange:
   Ecm = Ecm * m1
   zeta = zetaFunction( x_sq(Ecm), nP, alpha, nShell )
   zetaData_real.append(zeta.real)
   zetaData_imag.append(zeta.imag)
   print( Ecm / m1, " ", x_sq(Ecm), " ", zeta.real, " ", zeta.imag )
t2 = time.time()
#print("it took this long to calculate",t2-t1)

plt.plot(EcmRange, zetaData_real, label="real")
plt.plot(EcmRange, zetaData_imag, label="imag")
plt.show()


print(" ")
#roots=optimize.brentq(lambda Ecm,nP,alpha,nShell: 1./zetaFunction(x_sq(Ecm), nP, alpha, nShell ) , 2,3,args=(nP,alpha,nShell))
#print(roots)
