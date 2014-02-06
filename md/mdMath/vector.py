'''
vector.py

numpy is the very fast, easy module to the vector, array calculation.
But some complex vector calculation is not supported.
This module provides some complex vector functions based on the numpy module.

getVLength(v1) 
	: return the length of the v1.

getVAngle(v1, v2) 
	: return the angle between the v1 and v2.

getVFraction(v1):
	: Normalize the whole vector element as the sum of the vector elements is 1.
'''

import numpy as np
from math import sqrt, acos, pi 

def getVecLength(v1):
	return ((v1**2).sum())**0.5 



def getVecAngle(v1, v2):
	if len(v1) != len(v2):
		print "Same dimension vectors are needed."
		exit(1)

	cosAngle = np.dot(v1,v2) / (getVecLength(v1)*getVecLength(v2))
	angle = acos(cosAngle)*180/pi

	return angle




def getVecFraction(v1):
	dimension = len(v1)

	elementsSum = v1.sum()

	if elementsSum == 0:
		return v1
		
	else:
		fractionVec = v1/elementsSum
		return fractionVec

