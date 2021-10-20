#!/usr/bin/python3
import numpy as np
from numpy import tile
from scipy import linalg as LA
import math,cmath
pi = math.pi
cos = math.cos
sin = math.sin
exp=np.exp
sqrt = math.sqrt
log=np.log

def link(Nk, Nth, no_of_band_filled, evec):

    no_of_bonds=no_of_band_filled
    U1=np.zeros((2*Nth+1,2*Nk+1),dtype=complex)
    U2=np.zeros((2*Nth+1,2*Nk+1),dtype=complex)
    g1=np.zeros((no_of_bonds,no_of_bonds),dtype=complex)
    g2=np.zeros((no_of_bonds,no_of_bonds),dtype=complex)
    for nth in range(2*Nth+1):
        for nky in range(2*Nk+1):
            if (nth < 2*Nth):
                nthh=nth+1
            else:
                nthh=1
         
            for bond1 in range(no_of_bonds):
                for bond2 in range(no_of_bonds):
                    g1[bond1,bond2]=np.dot(evec[nth, nky,:,bond1],evec[nthh, nky,:,bond2])
                 
            S1=np.linalg.det(g1)
            U1[nth,nky]=S1/abs(S1)


    for nth in range(2*Nth+1):
        for nky in range(2*Nk+1):
            if (nky < 2*Nk):
                nkyy=nky+1
            else:
                nkyy=1
            for bond1 in range(no_of_bonds):
                for bond2 in range(no_of_bonds):
                    g2[bond1,bond2]=np.dot(evec[nth, nky,:,bond1],evec[nth, nkyy,:,bond2])
            S2=np.linalg.det(g2)
            U2[nth,nky]=S2/abs(S2)

    return U1, U2
