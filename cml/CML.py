from scipy import special as sp
import numpy as np
import math
import random

class CML():

    def __init__(self, initialMatrix,msize = -1):
        self.mat = self.__getGaussianMatrix(3)
        if(initialMatrix == 'gaussian') and (msize>0):
            self.mat = self.__getGaussianMatrix(msize)
        elif(initialMatrix == 'bessel') and (msize>0):
            self.mat = self.__getBesselMatrix(msize)
        elif(type(self.mat)==type(initialMatrix)) and (msize == -1):
            self.mat = initialMatrix
        else:
            raise Exception('Invalid CML input',initialMatrix,msize)
    
    def __getBesselMatrix(self,n):
        m = [[(10.0 * x / n - 5.0) ** 2.0 + (10.0 * y / n - 5.0) ** 2.0
            for x in range(0, n, 1)] for y in range(0, n, 1)]
        zBessel = sp.j0(m)
        minimo = np.min(zBessel)
        zBessel = zBessel - minimo
        maximo = np.max(zBessel)
        zBessel = zBessel / maximo
        return zBessel
        
    def __gaussian(self,x, y, meanX, meanY, stddev):
        return math.exp(-0.5 * (((x - meanX) / stddev) ** 2.0) -
                0.5 * (((y - meanY) / stddev) ** 2.0))
                
    def __getGaussianMatrix(self,n):
        #gaussian function
        m = [[self.__gaussian(10.0 * float(x) / n, 10.0 * float(y) / n, 5.0, 5.0, 1.5)
                for x in range(0, n, 1)] for y in range(0, n, 1)]
        return m

    def getCMLNeumann(self, mat,function, coupling,parameters=[]):
        neighborhood = [(1,0),(0,1),(-1,0),(0,-1)]
        return self.getCML(mat, neighborhood,function, coupling,parameters)

    def getCML(self, neighborhood,function, coupling,parameters=[]):
        outputMat = [row[:] for row in self.mat]
        rows = len(self.mat)
        cols = len(self.mat[0])
        for i in range(rows):
            for j in range(cols):
                outputMat[i][j] = coupling * function(self.mat[i][j],parameters)
                for n in neighborhood:
                    outputMat[i][j] += ((1.0 - coupling) / 4.0) * function(self.mat[(i+n[1]+rows) % rows][(j+n[0]+cols) % cols],parameters)
        self.mat = outputMat
        return outputMat

