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

def field(Nky, Nth, U1,U2):

    
    Field=np.zeros((2*Nth+1,2*Nky+1),dtype=complex)
    for nth in range(2*Nth+1):
        for nky in range(2*Nky+1):
            if (nth < 2*Nth):
                nthh=nth+1
            else:
                nthh=1
   
            if (nky < 2*Nky):
                nkyy=nky+1
            else:
                nkyy=1
            Field[nth,nky]=log((U1[nth,nky]*U2[nthh,nky])/(U1[nth,nkyy]*U2[nth,nky]))
    mu=sum(sum(i) for i in Field)
    mu=mu/(4.0*pi*1j)                       
    
    return mu.real