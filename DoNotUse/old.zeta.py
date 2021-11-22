import itertools
import math as m
import numpy as np
from scipy import special as s #special.erfi
import matplotlib.pyplot as plt
import cmath as c
import time

'''
 sqrt with branch cut aligned on positive real axis
'''
def mySqrt( z ):
   return 1j * c.sqrt( -z )

#print( mySqrt( 1.0 + 0.01 * 1j), mySqrt( 1.0 - 0.01 * 1j) )

'''
 inner product for two list (native)
'''
def dot(A, B):
   if len(A) != len(B):
      return 0
   return sum(i[0] * i[1] for i in zip(A, B))

'''
 Function to generate 3-tuple "n = (nx,ny,nz)" list 
 subject to the constraint nx^2 + ny^2 + nz^2 <= nShell^2,
 where "nShell" indicates the number of shells
'''
def generateIntList( nShell ):
    intList = []
    rawList = [*range(-nShell, nShell+1, 1)]

    for i in itertools.product(rawList,rawList,rawList):
        tmpList = i
        if ( dot(tmpList,tmpList) <= nShell*nShell ):
            intList.append(tmpList)

    return intList

'''
 returns the unit vector of a vector "A"
'''
def unit_vector( A ):
    if( all( A==0 ) ):
        return 0.0
    else:
        return A / np.sqrt( dot(A,A) )

'''
 returns a Lorentz boosted four vector "(A0,A)" given some boost
 velocity "beta"
'''
def lorentzBoost( A0, A, beta ):
    gamma = 1.0 / np.sqrt( 1.0 - dot(beta,beta) )
    Apar  = dot( A, unit_vector(beta) ) * unit_vector(beta)
    Aper  = A - Apar
    A0_prime = gamma * ( A0 - dot(A,beta) )
    A_prime  = gamma * ( Apar - A0 * beta ) + Aper
    return A0_prime, A_prime

'''
 KSS cutoff function
'''
def cutoff( x_sq, n_sq, alpha ):
    return np.exp( -alpha * ( n_sq - x_sq ) )

'''
 KSS Zeta function for ell = m_ell = 0, given "x_sq", total momentum "nP"
'''
def zetaFunction( x_sq, nP, alpha, nShell ):
    Ecm = np.sqrt( m1_sq + x_sq ) + np.sqrt( m2_sq + x_sq )
    beta = nP / np.sqrt( Ecm * Ecm + dot(nP,nP) )
    sum = 0.0
    for n in intList:
       n0 = np.sqrt( m1_sq + dot(n,n) )
       n0_cm, n_cm = lorentzBoost( n0, n, beta )
       tmp = cutoff( x_sq, dot(n_cm,n_cm), alpha ) * n0_cm / n0
       sum = sum + tmp / ( dot(n_cm,n_cm ) - x_sq )


    eps = 1.0e-10 * ( 1.0 - m.copysign(1, x_sq) )
    x   = c.sqrt( x_sq + 1j * eps )
    z   = c.sqrt( alpha * x * x )
    integral = 4.0 * m.pi * ( ( m.sqrt( 0.25 * m.pi / alpha ) ) * m.exp( alpha * x_sq )
                              - 0.5 * m.pi * x * s.erfi( z ) )
        
    return sum - integral - 4.0 * m.pi * ( 0.5 * 1j * m.pi * x ) # last term needed to compare to BBHO



#def finteVolumeF( E, P ):
   

'''
 convert Ecm to x_sq 
'''
def x_sq( Ecm ):
    return ( Ecm * Ecm - ( m1 + m2 )**2 ) * ( Ecm * Ecm - ( m1 - m2 )**2 ) / ( 4.0 * Ecm * Ecm )

'''
 running code
'''
b = np.array([.25,0,0])
v = np.array([1,1,1])
v0 = np.sqrt( 1 + v.dot(v) )
alpha = 0.1 #( 1.0 / 3.0 )**4
nShell = 15

t1 = time.time()
intList = generateIntList( nShell )
t2 = time.time()
print("it took this long to generate",t2-t1)

#print( lorentzBoost(v0,v,b) )

# m = m * L / ( 2*pi )
m1 = 4.0 / ( 2.0 * np.pi )
m2 = 4.0 / ( 2.0 * np.pi )
nP = np.array([0,0,1])
m1_sq = m1**2
m2_sq = m2**2

Ecm_o_m_start = 1.8
Ecm_o_m_stop = 4.0
Ecm_o_m_step = 0.05


t1 = time.time()
for Ecm in np.arange (Ecm_o_m_start, Ecm_o_m_stop, Ecm_o_m_step):
   Ecm = Ecm * m1
   zeta = zetaFunction( x_sq(Ecm), nP, alpha, nShell )
   print( Ecm / m1, " ", x_sq(Ecm), " ", zeta.real, " ", zeta.imag )
t2 = time.time()
print("it took this long to calculate",t2-t1)

