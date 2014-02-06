'''
slab.py
	this module contains the slab manipulating function
	When we analyze the MD trajectory, somtimes we need to 
	focus on the some specific fixed coordinates, not whole system.
	Before doing analysis, we should divide whole system into some specific
	slab wich contains coordinates we focus, count molecule in this slab
	and calculate some property within this slab, etc.
	 So I make one module seperately only dealing with slab-related job.

==============================================
checkPointAtomIndexArr(iCoord, pointAtomIndexArr, pointAtomNum, row, zCoord, zDelta)
	: when time goes t -> t+dt
	the atom can go out from certain slab which it resides at t.
	So we should delete these atoms in the slabAtomIndexArr
	Input
		List[n1, n2, .......]
								&
		float (atom number)
	Return
		List[n1, n2, ......]

checkSlabAtomIndexArr(ICoord, slabAtomIndexArr, slabAtomNumArr, row, zDim, zStep, slabNum)
	: when time goes t -> t+dt
	the atom can go out from certain slab which it resides at t.
	So we should delete these atoms in the slabAtomIndexArr
	Input
		List[[n11, n12, ...][n21, n22...] ... [nN1, nN2...]]
								&
		List[n1, n2, ... nN]
	Return
		List[[n11, n12, ...][n21, n22...] ... [nN1, nN2...]]

getPointAtomIndexArr(iCoord, zCoord, zDelta, atomNum, tI, dt)
	: this function 
	gets the index of atoms in the zCoord.
	Input
	coordinate file
		t   Atom1	Atom2	Atom3 .. 	AtomN
		.. 	x,y,z 	x,y,z 	x,y,z 		x,y,z
	Return
		List[n1, n2, ...]

getPointAtomNum(pointAtomIndexArr)
	: this function
	gets the index Arr and count the number of atoms 
	within each cell and make a Arr containg the number of atoms in each cell
	Input
		List[n1, n2, ...]
	Return
		int (num)

getPointAtomSurvPro(pointAtomNumRef, pointAtomNum)
	: this function
	gets the survival probability of the atoms in the zcoord between
	t and t+dt. Ref is t, noRef is t+dt
	Input
		int * 2
	Return
		int (probability)

getSlabZCoordArr(zDim, zStep, slabNum)
	: when we divide the vertial dimension into some slabs, 
	this function gets the 
	 Arr containing the z-coordinates which represents the each slabs.
	Return
		List[z1, z2, z3, ..., zN]

getSlabAtomIndexArr(iCoord, zDim, zStep, atomNum, slabNum, tI, dt)
	: this function
	distribute the atom into the each slab and 
	gets the index in the Arr
	Input
	coordinate file
		t   Atom1	Atom2	Atom3 .. 	AtomN
		.. 	x,y,z 	x,y,z 	x,y,z 		x,y,z
	Return 
		List[[n11, n12, ...][n21, n22...] ... [nN1, nN2...]]
			This List has slabNum inner bracket.

getSlabAtomNumArr(slabAtomIndexArr, slabNum)
	: this function
	gets the index Arr and count the number of atoms 
	within each cell and make a Arr containg the number of atoms in each cell
	Input
		List[[n11, n12, ...][n21, n22...] ... [nN1, nN2...]]
	Return
		List[n1, n2, ... nN]

getSlabAtomSurvArr(slabAtomNumArrRef, slabAtomNumArr, slabNum)
	: this function
	gets the survival probability of the atoms in the slab between
	t and t+dt. Ref is t, noRef is t+dt
	Input
		List[n1, n2, ... nN] * 2
	Return
		List[p1, p2, ... pN]
==============================================

All modules are related with the MDAnalysis module
'''
import numpy as np

def checkPointAtomIndexArr(iCoord, pointAtomIndexArr, pointAtomNum, row, zCoord, zDelta):
	deleteArr = []
	for j in range (1, pointAtomNum+1):
		colZ = 3*pointAtomIndexArr[j-1]+1 # column containing z coordinate of the atom
		if iCoord[row-1][colZ-1] >= zCoord+zDelta or iCoord[row-1][colZ-1] <= zCoord-zDelta:
			deleteArr.append(j-1)
	pointAtomIndexArr = np.delete(pointAtomIndexArr, deleteArr)
	return pointAtomIndexArr


def checkSlabAtomIndexArr(iCoord, slabAtomIndexArr, slabAtomNumArr, row, zDim, zStep, slabNum):
	deleteArr = []
	for i in range (1, slabNum+1):
		for j in range (1, slabAtomNumArr[i-1]+1):
			colZ = 3*slabAtomIndexArr[i-1][j-1]+1 # columnm containing z coordinate of the atom
			if int((iCoord[row-1][colZ-1]-zDim[0])/zStep) != i-1:
				deleteArr.append(j-1)
		slabAtomIndexArr[i-1] = np.delete(slabAtomIndexArr[i-1], deleteArr)
		deleteArr = []
	return slabAtomIndexArr


def getPointAtomIndexArr(iCoord, zCoord, zDelta, atomNum, tI, dt):
	pointAtomIndexArr = []
	### Atom is selected
	for i in range (1, atomNum+1):
		colZ = 3*i+1 	# column containing z coordinate of the atom
		rowI = int((tI-iCoord[0][0])/dt)+1 	# row represents t_i time data
		### Error checking : atoms outside the slab is neglected
		if iCoord[rowI-1][colZ-1] >= zCoord+zDelta or iCoord[rowI-1][colZ-1] <= zCoord-zDelta:
			a_a = 1 #Just Skip
		### If there is no Error, this atom is distributed
		else:
			pointAtomIndexArr.append(i)
	return pointAtomIndexArr


def getPointAtomNum(pointAtomIndexArr):
	pointAtomNum = len(pointAtomIndexArr)
	return pointAtomNum


def getPointAtomSurvPro(pointAtomNumRef, pointAtomNum):
	if pointAtomNum == 0:
		pointAtomSurvPro = 0
	else:
		pointAtomSurvPro = float(pointAtomNum) / float(pointAtomNumRef)
	return pointAtomSurvPro


def getSlabZCoordArr(zDim, zStep, slabNum):
	zCoordArr = []
	for i in range (1, slabNum+1):
		zCoord = zDim[0] + (i-0.5)*zStep
		zCoordArr.append(zCoord)
	return zCoordArr


def getSlabAtomIndexArr(iCoord, zDim, zStep, atomNum, slabNum, tI, dt):
	slabAtomIndexArr = []
	for i in range (1, slabNum+1):
		slabAtomIndexArr.append([]) # slabNum # inner bracket is included
	### Atom is distributed each inner bracket
	for i in range (1, atomNum+1):
		colZ = 3*i+1 	# column containing z coordinate of the atom
		rowI = int((tI-iCoord[0][0])/dt)+1 	# row represents t_i time data
		### Error checkgin : atoms outside the slab is neglected
		if iCoord[rowI-1][colZ-1] >= zDim[1] or iCoord[rowI-1][colZ-1] <= zDim[0]:
			a_a = 1 #Just Skip
		### If there is no Error, this atom is distributed
		else:
			slabAtomIndexArr[int((iCoord[rowI-1][colZ-1]-zDim[0])/zStep)].append(i)
	return slabAtomIndexArr


def getSlabAtomNumArr(slabAtomIndexArr, slabNum):
	slabAtomNumArr = []
	for i in range (1, slabNum+1):
		slabAtomNumArr.append(len(slabAtomIndexArr[i-1]))
	return slabAtomNumArr


def getSlabAtomSurvArr(slabAtomNumArrRef, slabAtomNumArr, slabNum):
	slabAtomSurvArr = []
	for i in range (1, slabNum+1):
		if slabAtomNumArr[i-1] == 0:
			slabAtomSurvArr.append(0)
		else:
			slabAtomSurvArr.append(float(slabAtomNumArr[i-1])/float(slabAtomNumArrRef[i-1]))
	return slabAtomSurvArr

