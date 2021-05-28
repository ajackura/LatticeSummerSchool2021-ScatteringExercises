import itertools
import math as m
import numpy as np
import matplotlib.pyplot as plt
from numpy.core.umath_tests import inner1d

'''


'''


'''
 We work with an array of arrays, e.g. n = [ n0, n1, ..., nN ] where n0 = [ n0[0], n0[1], n0[2] ], etc.
'''

twoPi = 2.0 * m.pi

def free_spectrum( m1_sq, m2_sq, n, nP ):
    energy1 = np.sqrt( m1_sq + inner1d(n,n) )
    energy2 = np.sqrt( m2_sq + inner1d(nP-n,nP-n) )
    return energy1 + energy2

def generate_free_spectrum( mL, mRatio, nP, EcmMax ):
    m1 = mL / twoPi
    m2 = m1 * mRatio
    m1_sq = m1**2
    m2_sq = m2**2
    nMax = 3
    Energies = []
    intList = np.array(np.arange(-nMax,nMax+1))
    for i in itertools.product(intList,intList,intList):
        n = np.array(i)
        En = free_spectrum(m1_sq, m2_sq, n, nP)
        Ecm = np.sqrt( En**2 - np.dot(nP,nP) )
        Ecm = 2.0 * m.pi * Ecm / mL
        Ecm = np.round( Ecm, 6 )
        if Ecm not in Energies:
            if Ecm < EcmMax:
                Energies.append(Ecm)
    Energies.sort()
    return Energies

'''
 running code
'''

nP = np.array([0,0,2])
Ecm_o_m_Max = 5.0
mRatio = 1.0
for mL in np.arange (2.0, 8.0, 0.1):
   Ecm = generate_free_spectrum( mL, mRatio, nP, Ecm_o_m_Max )
   print( round(mL,6), " ", *Ecm )

