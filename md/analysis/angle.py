'''
calPonitAngDist(iCoord, dVec, pointAtomIndexArr, row, numThetaArr)

'''

import numpy as np
import math
from scipy.optimize import curve_fit

def calPointAngDist(iCoord, dVec, pointAtomIndexArr, pointAtomNum, row, numThetaArr):
	bufThetaArr = np.zeros(numThetaArr)
	deltaTheta = 180/numThetaArr
	for j in range (1, pointAtomNum+1):
		colZ = 3*pointAtomIndexArr[j-1]+1 # colZ in dVec <- this 3 is different with cal2DPointRmsd!!
		dVecZ = dVec[row-1][colZ-1]	# dipole vector z value
		if dVecZ <= -1.00:
			theta = 179.9
		elif dVecZ >= 1.00:
			theta = 0.1
		else:
			theta = (180/math.pi) * math.acos(dVecZ)
		rowBufThetaArr = int(theta/deltaTheta) + 1
		bufThetaArr[rowBufThetaArr-1] += 1
	return bufThetaArr/pointAtomNum




