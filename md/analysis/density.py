'''
density.py
	This module contains the analysis function used to calculte the density profile.


getZcoordSlabArray(zI, zF, zStep)
	: return the 
	1. array contains the reference z-coordinate when calculate the density.
	2. z axis slab Num.

getXYZcoordSlabArray(xI, xF, xStep, yI, yF, yStep, zI, zF, zStep)
	: return the 
	1. 4D list contains the reference xyz-coordinate when calculate the 2D density map
	as a function of the z coordinate.
	2. x, y, z axis slab Num respectively.

getSlabAtomIndexList(row, oCoord, zCoordSlabArray)
	: It makes the 
	1. assymetrical 2D array contains the atom index at each slab.
	2. 1D array contains the atom Num at each slab.

getXYZSlabAtomIndexList(row, oCoord, xyzCoordSlabArray)
	: return the 
	1. 4D list contains the index of the atoms corresponds to the each 3D cells.
	2. 3D array contains the atom Num at each cell.

getXYZSlabAtomNumArray(row, oCoord, xyzCoordSlabArray)
	: return the
	1. 3D array contains the atom Num at each cell.

getxyz3DDensity(oCoord, xI, xF, xStep, yI, yF, yStep, zI, zF, zStep)
	: return the 
	1. 4D(x*y*z*1) array contains the density of the target molecule.

'''

import numpy as np
import math

def getZcoordSlabArray(zI, zF, zStep):
	zcoordSlabArray = []
	slabNum = int(math.ceil((zF-zI)/zStep))

	for i in range (1, slabNum+1):
		zcoordSlabArray.append(zI + (i-1)*zStep + (0.5*zStep))

	return np.array(zcoordSlabArray), slabNum



def getXYZcoordSlabArray(xI, xF, xStep, yI, yF, yStep, zI, zF, zStep):
	slabXNum = int(math.ceil((xF-xI)/xStep))
	slabYNum = int(math.ceil((yF-yI)/yStep))
	slabZNum = int(math.ceil((zF-zI)/zStep))
	# last variable 3 means that each XYZ cell has XYZ-coordinate array.
	xyzcoordSlabArray = np.zeros([slabXNum, slabYNum, slabZNum, 3])

	for i in range (1, slabXNum+1):
		for j in range (1, slabYNum+1):
			for k in range (1, slabZNum+1):
				xyzcoordSlabArray[i-1][j-1][k-1][0] += (xI + (i-1)*xStep + (0.5*xStep))
				xyzcoordSlabArray[i-1][j-1][k-1][1] += (yI + (j-1)*yStep + (0.5*yStep))
				xyzcoordSlabArray[i-1][j-1][k-1][2] += (zI + (k-1)*zStep + (0.5*zStep))

	return xyzcoordSlabArray, slabXNum, slabYNum, slabZNum

	


def getSlabAtomIndexList(row, oCoord, zCoordSlabArray):
	slabNum = len(zCoordSlabArray)
	zStep = zCoordSlabArray[1]-zCoordSlabArray[0]
	zMin = zCoordSlabArray[0] - 0.5*zStep
	zMax = zCoordSlabArray[slabNum-1] + 0.5*zStep

	atomNum = int((np.shape(oCoord)[1]-1)/3)

	slabAtomIndexList = []
	for i in range (1, slabNum+1):
		slabAtomIndexList.append([])
	for i in range (1, atomNum+1):
		if oCoord[row-1][i*3] >= zMax or oCoord[row-1][i*3] < zMin:
			check = 1
		else:
			slabAtomIndexList[int((oCoord[row-1][i*3]-zMin) / zStep)].append(i)

	slabAtomNumArray = np.zeros(slabNum)
	for i in range (1, slabNum+1):
		slabAtomNumArray[i-1] += len(slabAtomIndexList[i-1])

	return slabAtomIndexList, slabAtomNumArray




def getXYZSlabAtomIndexList(row, oCoord, xyzCoordSlabArray):
	slabXNum = np.shape(xyzCoordSlabArray)[0]
	slabYNum = np.shape(xyzCoordSlabArray)[1]
	slabZNum = np.shape(xyzCoordSlabArray)[2]
	xStep = xyzCoordSlabArray[1][0][0][0] - xyzCoordSlabArray[0][0][0][0]
	yStep = xyzCoordSlabArray[0][1][0][1] - xyzCoordSlabArray[0][0][0][1]
	zStep = xyzCoordSlabArray[0][0][1][2] - xyzCoordSlabArray[0][0][0][2]
	xMin = xyzCoordSlabArray[0][0][0][0] - 0.5*xStep
	yMin = xyzCoordSlabArray[0][0][0][1] - 0.5*yStep
	zMin = xyzCoordSlabArray[0][0][0][2] - 0.5*zStep
	xMax = xyzCoordSlabArray[slabXNum-1][0][0][0] + 0.5*xStep
	yMax = xyzCoordSlabArray[0][slabYNum-1][0][1] + 0.5*yStep
	zMax = xyzCoordSlabArray[0][0][slabZNum-1][2] + 0.5*zStep

	atomNum = int((np.shape(oCoord)[1]-1)/3)

	slabAtomIndexList = np.empty([slabXNum, slabYNum, slabZNum, 0])
	slabAtomIndexList = slabAtomIndexList.tolist()
	for i in range (1, atomNum+1):
		if oCoord[row-1][i*3-2] >= xMax or oCoord[row-1][i*3-2] < xMin or \
		   oCoord[row-1][i*3-1] >= yMax or oCoord[row-1][i*3-1] < yMin or \
		   oCoord[row-1][i*3] >= zMax or oCoord[row-1][i*3] < zMin:
			check = 1
		else:
			slabAtomIndexList[int((oCoord[row-1][i*3-2]-xMin)/xStep)][int((oCoord[row-1][i*3-1]-yMin)/yStep)][int((oCoord[row-1][i*3]-zMin)/zStep)].append(i)

	slabAtomNumArray = np.zeros([slabXNum, slabYNum, slabZNum, 1])
	for i in range (1, slabXNum+1):
		for j in range (1, slabYNum+1):
			for k in range (1, slabZNum+1):
				slabAtomNumArray[i-1][j-1][k-1][0] += len(slabAtomIndexList[i-1][j-1][k-1])

	return slabAtomIndexList, slabAtomNumArray



def getXYZSlabAtomNumArray(row, oCoord, xyzCoordSlabArray):
	slabXNum = np.shape(xyzCoordSlabArray)[0]
	slabYNum = np.shape(xyzCoordSlabArray)[1]
	slabZNum = np.shape(xyzCoordSlabArray)[2]
	xStep = xyzCoordSlabArray[1][0][0][0] - xyzCoordSlabArray[0][0][0][0]
	yStep = xyzCoordSlabArray[0][1][0][1] - xyzCoordSlabArray[0][0][0][1]
	zStep = xyzCoordSlabArray[0][0][1][2] - xyzCoordSlabArray[0][0][0][2]
	xMin = xyzCoordSlabArray[0][0][0][0] - 0.5*xStep
	yMin = xyzCoordSlabArray[0][0][0][1] - 0.5*yStep
	zMin = xyzCoordSlabArray[0][0][0][2] - 0.5*zStep
	xMax = xyzCoordSlabArray[slabXNum-1][0][0][0] + 0.5*xStep
	yMax = xyzCoordSlabArray[0][slabYNum-1][0][1] + 0.5*yStep
	zMax = xyzCoordSlabArray[0][0][slabZNum-1][2] + 0.5*zStep

	atomNum = int((np.shape(oCoord)[1]-1)/3)

	slabAtomNumArray = np.zeros([slabXNum, slabYNum, slabZNum, 1])

	for i in range (1, atomNum+1):
		if oCoord[row-1][i*3-2] >= xMax or oCoord[row-1][i*3-2] < xMin or \
		   oCoord[row-1][i*3-1] >= yMax or oCoord[row-1][i*3-1] < yMin or \
		   oCoord[row-1][i*3] >= zMax or oCoord[row-1][i*3] < zMin:
			check = 1
		else:
			slabAtomNumArray[int((oCoord[row-1][i*3-2]-xMin)/xStep)][int((oCoord[row-1][i*3-1]-yMin)/yStep)][int((oCoord[row-1][i*3]-zMin)/zStep)] += 1

	return slabAtomNumArray



def getxy2DDensity(oCoord, xI, xF, xStep, yI, yF, yStep, zI, zF, zStep):
	rowNum = np.shape(oCoord)[0]
	molNum = int((np.shape(oCoord)[1]-1)/3)

	xyzSlab, slabXNum, slabYNum, slabZNum = getXYZcoordSlabArray(xI, xF, xStep, yI, yF, yStep, zI, zF, zStep)
	buf_slabAtomNumArray = np.zeros([slabXNum, slabYNum, slabZNum, 1])
	for row in range (1, rowNum+1):
		if row%100 == 0:
			print row

		slabAtomNumArray = getXYZSlabAtomNumArray(row, oCoord, xyzSlab)
		buf_slabAtomNumArray += slabAtomNumArray

	slabAtomDenArray = buf_slabAtomNumArray/rowNum

	return buf_slabAtomDenArray, slabZNum


