# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
import numpy as np
np.set_printoptions(linewidth=np.nan)
np.set_printoptions(threshold=np.inf)
np.set_printoptions(precision=2)

import matplotlib.pyplot as plt

from scipy.optimize import fsolve
import scipy.optimize as sopt


from scipy.linalg import block_diag
from scipy.linalg import solve
from scipy.linalg import lstsq
# -

pip install --upgrade pip


def symmetrize(A):
    return A+A.T


# +
nPts=5
nAs=5
eps=2

Dists=np.zeros((nPts,nPts))
Dists[np.arange(nPts-1),1+np.arange(nPts-1)]=2
Dists=symmetrize(Dists)
linked,linkedTo=np.nonzero(Dists)
sZ=nPts*nAs

Adj=np.zeros_like(Dists)
Adj[linked,linkedTo]=1
nLip=np.sum(Adj,dtype=np.int32)
sLip=nLip*nAs
    

P=np.ones((nAs,nPts))
P/=P.sum(axis=0)[None,:]

hatP=np.random.rand(nAs,nPts)
hatP/=hatP.sum()
hatPFlat=hatP.flatten()

#print(P)
#print(hatP)
#print(Dists)
#print(Adj)

# +
fig,axs=plt.subplots(nAs,figsize=(15,15))

for a in range(nAs):
    axs[a].scatter(np.arange(nPts),P[a,:],color='blue')
    axs[a].scatter(np.arange(nPts),hatP[a,:],color='red')
# -

# # Explicit Newton (pas top)

mu=1e-8


# +

def Constr(Z):
    cLip=np.zeros((nAs,nLip))
    cMass=np.zeros(nPts)
    for a in range(nAs):
        cLip[a,:]=np.exp(Dists[linked,linkedTo]/eps)*Z[a,linked]-Z[a,linkedTo]
    cMass=1-np.sum(Z*P, axis=0)
    return cLip,cMass

def Grad(Z):
    G=np.zeros_like(Z)
    G-=hatP/Z
    cLip,cMass=Constr(Z)
    for a in range(nAs):
        G[a,linked]-=mu*np.exp(Dists[linked,linkedTo]/eps)/cLip[a]
        G[a,linkedTo]+=mu/cLip[a]
    G+=mu*P/cMass.reshape(1,-1)
    return G.flatten()

def Hess(Z):
    cLip,cMass=Constr(Z)
    H=np.zeros((nAs*nPts,nAs*nPts))
    
    # Minimized energy term:
    H+=np.diag((hatP/Z**2).flatten())
    
    
    
    # Barriers:
    for a in range(nAs):
        # Mass constraint barrier:
        for b in range(nAs):
            H[a*nPts:(a+1)*nPts,b*nPts:(b+1)*nPts]+=mu*np.diag(P[a,:]*P[b,:]/cMass[:]**2)
            
        #Lip constraint barrier:
        H[a*nPts+linked,a*nPts+linkedTo]-=mu*np.exp(Dists[linked,linkedTo]/eps)/cLip[a]**2
        H[a*nPts+linkedTo,a*nPts+linked]-=mu*np.exp(Dists[linked,linkedTo]/eps)/cLip[a]**2
        H[a*nPts+linked,a*nPts+linked]+=mu*np.exp(2*Dists[linked,linkedTo]/eps)/cLip[a]**2
        H[a*nPts+linkedTo,a*nPts+linkedTo]+=mu/cLip[a]**2
    return H

def Energy(Z):
    cLip,cMass=Constr(Z)
    E=-np.sum(np.log(Z)*hatP)
    barLip=mu*np.sum(np.log(cLip))
    barMass=mu*np.sum(np.log(cMass))
    return E-barLip-barMass


# -

Z=np.ones((nAs,nPts))/(1+np.sum(P))
print(Constr(Z))
print(Grad(Z),Grad(Z).shape)
print(Hess(Z),np.linalg.det(Hess(Z)))
print(Energy(Z))

# +
lambdaMax=(3-np.sqrt(5))/2

def OneNewton(Z, verbose=False):
    H=Hess(Z)
    invH=np.linalg.inv(H)
    G=Grad(Z)
    try:
        p=solve(H, -G).reshape(nAs,nPts)
    except:
        print(np.linalg.det(Hess(Z)))
        print(Z)
    lambdaZ=np.sqrt(G.T@invH@G)
    if lambdaZ<lambdaMax:
        alpha=1
    else:
        print(lambdaZ)
        alpha=1/(1+lambdaZ)

    if verbose:
        print(alpha)
    
    return Z+alpha*p


# +
mu=1e-8
Z=np.ones((nAs,nPts))/(1+np.sum(P))

for i in range(10):
    print('E=',Energy(Z))
    Z=OneNewton(Z,True)
    cLip,cMass=Constr(Z)
    print('min Z=',np.min(Z))
    print('min Lip=',np.min(cLip))
    print('min Mass=',np.min(cMass))


# +
def fGrad(Z):
    return Grad(Z.reshape(nAs,nPts))

def fHess(Z):
    return Hess(Z.reshape(nAs,nPts))


# -

V=np.ones((nAs,nPts))/(1+np.sum(P,axis=0)).reshape(1,-1)
Z=V.flatten()
print(Constr(V))

root=fsolve(func=fGrad, x0=Z,fprime=fHess)
print(Constr(root.reshape(nAs,nPts)))

# # Using trust-cnstr scipy method

jacExp=np.zeros((nPts,sPhi))
PFlat=P.flatten()
for a in range(nAs):
    jacExp[:,a*nPts:(a+1)*nPts]=np.diag(PFlat[a*nPts:(a+1)*nPts])

constrExp=sopt.LinearConstraint(A=jacExp,ub=np.ones(nPts))

print(sLambdaLip)

constrLip=sopt.LinearConstraint(A=jacLip,lb=np.zeros(sLambdaLip))


def Objective(V):
    return -eps*np.sum(np.log(V)*hatPFlat)


# +
def GradObj(V):
    return -eps*hatPFlat/V

def HessObj(V):
    return eps*np.diag(hatPFlat/V**2)


# -

eps=2
V0=np.random.rand(sPhi)
print(Objective(V0))

Vopt=sopt.minimize(fun=Objective, x0=V0, method='trust-constr', jac=GradObj, hess=HessObj, constraints=[constrExp,constrLip],options={'disp':True})

RNorm=(Vopt.x*PFlat).reshape(nAs,nPts)
print(Objective(Vopt.x))

# +
fig,axs=plt.subplots(nAs,figsize=(15,15))

for a in range(nAs):
    axs[a].scatter(np.arange(nPts),RNorm[a,:],color='green')
    axs[a].scatter(np.arange(nPts),hatP[a,:],color='red')
# -

matR=np.zeros(((nPts+1)*nAs,nPts*nAs))
matR[:sPhi]+=np.identity(sPhi)
matR[sPhi:]+=block_diag(*[np.ones(nPts) for a in range(nAs)])
for a in range(nAs):
    for b in range(nAs):
        matR[a*nPts:(a+1)*nPts, b*nPts:(b+1)*nPts]-=np.diag(Vopt.x[a*nPts:(a+1)*nPts]*PFlat[a*nPts:(a+1)*nPts])

plt.imshow(matR)

np.linalg.det(matR[:sPhi])

massVec=np.zeros(nAs*(nPts+1))
massVec[-nAs:]=np.sum(hatP,axis=1)

Ropt, res, rnk, s=lstsq(matR, massVec)
Ropt=Ropt.reshape(nAs,nPts)

print(np.sum(res))

# +
fig,axs=plt.subplots(nAs,figsize=(15,15))

for a in range(nAs):
    axs[a].scatter(np.arange(nPts),Ropt[a,:],color='green')
    axs[a].scatter(np.arange(nPts),hatP[a,:],color='red')
# -


