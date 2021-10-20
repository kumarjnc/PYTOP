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

def bandstr(tx,ty,Nk,Ntheta,size_bzone,alpha,lmda,gama,square=True):
    if square:
    	assert (square==True), "Only square 2D lattice!"
    	print('Initialisation of tbg')
    	n=size_bzone
    	sigmaup=np.zeros((n,n),dtype=complex)
    	sigmadw=np.zeros((n,n),dtype=complex)
    	sigmaud=np.zeros((n,n),dtype=complex)
    	sigmadu=np.zeros((n,n),dtype=complex)
    	comp=np.zeros((2*n,2*n),dtype=complex)
    	energy=np.zeros((2*Ntheta+1,2*Nk+1,2*n))
    	energy_val=np.zeros((2*Ntheta+1,2*Nk+1,2*n,2*n),dtype=complex)
    	
#-----------------------------------------------------------------------------
    	A=0.0
    	nn=0.0
    	nth_m=[]
    	for nth in range(2*Ntheta+1):
    		theta1=pi/Ntheta*(nth-Ntheta)
    		nth_m.append(theta1)
    		nky_m=[]
    		for nky in range(2*Nk+1):
    			ky=pi/Nk*(nky-Nk)
    			nky_m.append(ky)
    	
    			for q in range(n):
    				for p in range(n):
    					x=p 
    					y=q
    					if(p==q):
    						sigmaup[q,p]=-ty*(exp(1j*ky)*exp(1j*2.0*(x)*pi/alpha)+exp(-1j*ky)*exp(-1j*2.0*(x)*pi/alpha))+(-1)**x*lmda+A-(-1)**y*nn



    			for q in range(n):
    				for p in range(n):
    					x=p 
    					y=q
    					if((p-q)==1):
    						sigmaup[q,p]=-tx*cos(2*pi*gama)

    			for q in range(n):
    				for p in range(n):
    					x=p 
    					y=q
    					if((q-p)==1):
    						sigmaup[q,p]=-tx*cos(2*pi*gama)

    			sigmaup[1,q]=-tx*cos(2*pi*gama)*exp(-1j*theta1)
    			sigmaup[q,1]=-tx*cos(2*pi*gama)*exp(1j*theta1)

    			for q in range(n):
    				for p in range(n):
    					x=p 
    					y=q
    					if(p==q):
    						sigmadw[q,p]=-ty*(exp(1j*ky)*exp(1j*2.0*(x)*pi/alpha)+exp(-1j*ky)*exp(-1j*2.0*(x)*pi/alpha))+(-1)**x*lmda+A-(-1)**y*nn



    			for q in range(n):
    				for p in range(n):
    					x=p 
    					y=q
    					if((p-q)==1):
    						sigmadw[q,p]=-tx*cos(2*pi*gama)

    			for q in range(n):
    				for p in range(n):
    					x=p 
    					y=q
    					if((q-p)==1):
    						sigmadw[q,p]=-tx*cos(2*pi*gama)

    			sigmadw[1,q]=-tx*cos(2*pi*gama)*exp(-1j*theta1)
    			sigmadw[q,1]=-tx*cos(2*pi*gama)*exp(1j*theta1)

    			for q in range(n):
    				for p in range(n):
    					x=p 
    					y=q
    					if((p-q)==1):
    						sigmadu[q,p]=-tx*(1j*sin(2*pi*gama))

    			for q in range(n):
    				for p in range(n):
    					x=p 
    					y=q
    					if((q-p)==1):
    						sigmadu[q,p]=-tx*(-1j*sin(2*pi*gama))

    			sigmadu[1,q]=-tx*(-1j*sin(2*pi*gama)*exp(-1j*theta1))
    			sigmadu[q,1]=-tx*(1j*sin(2*pi*gama)*exp(-1j*theta1))


    			for q in range(n):
    				for p in range(n):
    					x=p 
    					y=q
    					if((p-q)==1):
    						sigmaud[q,p]=-tx*(1j*sin(2*pi*gama))

    			for q in range(n):
    				for p in range(n):
    					x=p 
    					y=q
    					if((q-p)==1):
    						sigmaud[q,p]=-tx*(-1j*sin(2*pi*gama))

    			sigmaud[1,q]=-tx*(-1j*sin(2*pi*gama)*exp(1j*theta1))
    			sigmaud[q,1]=-tx*(1j*sin(2*pi*gama)*exp(1j*theta1))

    			
    			for q in range(n):
    				for p in range(n):
    					comp[2*q,2*p]=sigmaup[q,p]

    			for q in range(n):
    				for p in range(n):
    					comp[2*q-1,2*p-1]=sigmadw[q,p]
    			for q in range(n):
    				for p in range(n):
    					comp[2*q-1,2*p]=sigmadu[q,p]

    			for q in range(n):
    				for p in range(n):
    					comp[2*q,2*p-1]=sigmaud[q,p]
#--------------------------------------------Eigval and Eigvec----------------------------------
    			vk, wk=np.linalg.eigh(comp)
    			vkk=np.array(vk)
    			wkk=np.array(wk)
    			for i in range(2*n):
    				energy[nth,nky,i]=vkk[i].real
    			for i in range(2*n):
    				for j in range(2*n):
    					energy_val[nth,nky,i,j]=wkk[i,j]
    return energy, energy_val,nth_m , nky_m
