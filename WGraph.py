from scipy.sparse import csr_matrix

class WGraph:
    def __init__(self, nPts, nAs, dMatrix):
        self.nPts=nPts
        self.nAs=nAs
        self.nHat=0
        self.hatP=np.zeros((nAs,nPts))
        #/!\ give dMatrix as sparse csr/csc/coo matrix.
        self.dMatrix=csr_matrix(dMatrix, shape=(nPts,nPts),copy=True)
        self.aMatrix=csr_matrix(dMatrix, shape=(nPts,nPts),dtype=int32,copy=True)
        self.aMatrix.data[:]=1
    
    def getDists(self, sparse=True):
        if sparse:
            return self.dMatrix
        else:
            return 
    def getAdjcy(self, sparse=True):
        return self.aMatrix

    def updateHatP(self,samples):
        """
        Updates value of P hat based on observations in samples. Samples is a 2 x nSamples np array of integers giving the
        observed class a and point x for each sample.
        """
        
        self.hatP*=self.nHat
        self.hatP[*samples]+=1
        self.nHat+=samples.shape[-1]
        self.hatP/=self.nHat
        
    def getHatP(self):
        return self.hatP


    def setP(self,P):
        self.P=P
    
    def getP(self):
        return self.P
    
    def computeConstr(self, Z):
        linked,linkedTo=self.aMatrix.nonzero()
        sumDeg=linked.shape[0]
        cLip=np.zeros((self.nAs,sumDeg))
        cMass=np.zeros(nPts)
        for a in range(nAs):
            cLip[a,:]=np.exp(self.dMatrix.data[:]/eps)*Z[a,linked]-Z[a,linkedTo]
        cMass=1-np.sum(Z*P, axis=0)
        return cLip,cMass
        
        
    def computeEnergy(self,Z):
        cLip, cMass=self.computeConstr(Z)

    def computeGrad(self):
        

    def computeHess(self):

    
