#! /usr/bin/env python

import itertools
import math as m
import numpy as np
import matplotlib.pyplot as plt
# from numpy.core.umath_tests import inner1d

'''


'''


'''
 We work with an array of arrays, e.g. n = [ n0, n1, ..., nN ] where n0 = [ n0[0], n0[1], n0[2] ], etc.
'''
twoPi = 2.0 * m.pi

def free_spectrum( m1_sq, m2_sq, n, nP ):
    energy1 = np.sqrt( m1_sq + np.inner(n,n) )
    energy2 = np.sqrt( m2_sq + np.inner(nP-n,nP-n) )
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



def get_spectrum_data(nP,mRatio,EcmMax,LRange):
    data=[]
    maxLevel = 0
    for mL in LRange:
        Ecm = generate_free_spectrum( mL, mRatio, nP, EcmMax )
        if len(Ecm) > maxLevel:
            maxLevel = len(Ecm)
        data.append(Ecm)

    level = 0
    L_start = 0
    levelData = []
    returnData = []
    for _ in range(maxLevel):
        for i in data:
            if len(i) >= level+1:
                levelData.append(i[level])
            else:
                L_start += 1
        returnData.append([L_list[L_start:],levelData.copy()])
        level += 1
        L_start = 0
        levelData.clear()

    return returnData

'''
 running code
'''
nPList = [[0,0,0],[0,0,1],[0,1,1],[1,1,1],[0,0,2]]
Lmin = 2.0
Lmax = 8.0
Lstep = 0.1
L_list = np.arange(Lmin,Lmax,Lstep)
mRatio = 1.0
Ecm_o_m_Max = 5.0

fig = plt.figure(figsize=(13,5))

for i,nP in enumerate(nPList):
    tmp = get_spectrum_data(nP,mRatio,Ecm_o_m_Max,L_list)
    plot = fig.add_subplot(1,len(nPList),i+1)
    plot.set_title(str(nP))
    plot.set_xlabel("L")
    if i == 0:
        plot.set_ylabel(r'$E_{cm}$')
    for j in tmp:
        plot.plot(j[0],j[1],color="tab:blue")

plt.show()

