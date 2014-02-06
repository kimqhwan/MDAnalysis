'''
rmsd.py
	this module contains the rmsd related function

==============================================
cal2DLocalRmsd(iCoord, slabAtomIndexArr, slabAtomNumArr, slabAtomNumArrI, rowI, row, time, slabNum)
	: Calculate the local 2D RMSD
	Return
		List(time, Rmsd1, Rmsd2, Rmsd3, ... RmsdN)

cal2DPointRmsd(iCoord, pointAtomIndexArr, pointAtomNum, pointAtomNumI, rowI, row, time)
	: Calculate the local 2D RMSD
	Return
		List(time, Rmsd)

cal2DRmsd(iCoord0, iCoordt, atomNum)
	: Calculate the 2D RMSD
	Return 
		Float(Rmsd)

cal3DRmsd(iCoord0, iCoordt, atomNum)
	: Calculate the 3D RMSD
	Return 
		Float(Rmsd)

fit2DDiff(timeArr, rmsdArr)
	: Linear fit the RMSD and get the 2D diffusion
	Input
		List(t1, t2, ... , tN)
				&
		List(r1, r2, ... , rN)
	Return
		Float(Diffusion Coeff.)

fit3DDiff(timeArr, rmsdArr)
	: Linear fit the RMSD and get the 2D diffusion
	Input
		List(t1, t2, ... , tN)
				&
		List(r1, r2, ... , rN)
	Return
		Float(Diffusion Coeff.)

norm2DPointRmsd(rmsd, pointAtomSurvArr)
	: When RMSD is calculated, this results is from the N_i
	So we need to normalize this data using N_t/N_i, survival probaility.
	Input
		t RMSD1 
		.. ..   
						&
		List[P1 P2 P3 ... PN]

	Output
		t NRMSD1 
		.. ..    

norm2DLocalRmsd(rmsd, slabAtomSurvArr, slabNum)
	: When RMSD is calculated, this results is from the N_i
	So we need to normalize this data using N_t/N_i, survival probaility.
	Input
		t RMSD1 RMSD2 RMSD3 ... RMSDN 
		.. ..   ..    ..    ..  ..
						&
		P1 P2 P3 ... PN
		.. .. .. ..  ..
	Output
		t NRMSD1 NRMSD2 NRMSD3 ... NRMSDN
		.. ..    ..    ..       .. ..
==============================================

All modules are related with the MDAnalysis module
'''

import numpy as np
from scipy.optimize import curve_fit


def cal2DLocalRmsd(iCoord, slabAtomIndexArr, slabAtomNumArr, slabAtomNumArrI, rowI, row, time, slabNum):
	rmsd = []
	rmsd.append(time)
	local_rmsd = 0.
	for i in range (1, slabNum+1):
		for j in range (1, slabAtomNumArr[i-1]+1):
			colX = 3*slabAtomIndexArr[i-1][j-1]-1
			colY = 3*slabAtomIndexArr[i-1][j-1]
			local_rmsd += (iCoord[row-1][colX-1]-iCoord[rowI-1][colX-1])**2 \
						+ (iCoord[row-1][colY-1]-iCoord[rowI-1][colY-1])**2
		if slabAtomNumArr[i-1] == 0:
			rmsd.append(0)
		else:
			rmsd.append(local_rmsd/slabAtomNumArrI[i-1])
		local_rmsd = 0.
	return rmsd

def cal2DPointRmsd(iCoord, pointAtomIndexArr, pointAtomNum, pointAtomNumI, rowI, row, time):
	rmsd = []
	rmsd.append(time)
	point_rmsd = 0.
	for j in range (1, pointAtomNum+1):
		colX = 3*pointAtomIndexArr[j-1]-1
		colY = 3*pointAtomIndexArr[j-1]
		point_rmsd += (iCoord[row-1][colX-1]-iCoord[rowI-1][colX-1])**2 \
					+ (iCoord[row-1][colY-1]-iCoord[rowI-1][colY-1])**2
	if pointAtomNum == 0:
		rmsd.append(0)
	else:
		rmsd.append(point_rmsd/pointAtomNumI)
	local_rmsd = 0.
	return rmsd


def cal2DRmsd(iCoord0, iCoordt, atomNum):
	disp_xy = np.delete(iCoordt-iCoord0, 2, 1)

	summ = 0.
	for row in disp_xy:
		summ += (np.linalg.norm(row))**2

	rmsd = summ / atomNum

	return rmsd


def cal3DRmsd(iCoord0, iCoordt, atomNum):
	disp = iCoordt-iCoord0

	summ = 0.
	for row in disp:
		summ += (np.linalg.norm(row))**2

	rmsd = summ / atomNum

	return rmsd


def fit2DDiff(timeArr, rmsdArr):
	def func(t, D):
		return 4*D*t
	#### Delete the first and last 10% data before fitting###
	frameNum = len(timeArr)
	exclNum = int(frameNum/10.)	# first and last 10% are excluded
	newTimeArr = []
	newRmsdArr = []
	for i in range (1, frameNum):
		if i > exclNum and i < frameNum-exclNum:
			newTimeArr.append(timeArr[i-1])
			newRmsdArr.append(rmsdArr[i-1])
	#####

	D = float(curve_fit(func, np.asarray(newTimeArr), np.asarray(newRmsdArr))[0]*1000)

	return D


def fit3DDiff(timeArr, rmsdArr):
	def func(t, D):
		return 6*D*t
	#### Delete the first and last 10% data before fitting###
	frameNum = len(timeArr)
	exclNum = int(frameNum/10.)	# first and last 10% are excluded
	newTimeArr = []
	newRmsdArr = []
	for i in range (1, frameNum):
		if i > exclNum and i < frameNum-exclNum:
			newTimeArr.append(timeArr[i-1])
			newRmsdArr.append(rmsdArr[i-1])
	#####

	D = float(curve_fit(func, newTimeArr, newRmsdArr)[0]*100)

	return D


def norm2DPointRmsd(rmsd, pointAtomSurvArr):
	nRMSD = []
	nRMSDt = []
	rNum = np.shape(rmsd)[0]
	for row in range (1, rNum+1):
		nRMSDt.append(rmsd[row-1][0]) # append the time
		if pointAtomSurvArr[row-1] == 0.:
			nRMSDt.append(0.)
		else:
			nRMSDt.append(rmsd[row-1][1]/pointAtomSurvArr[row-1])
		nRMSD.append(nRMSDt)
		nRMSDt = []
	return nRMSD


def norm2DLocalRmsd(rmsd, slabAtomSurvArr, slabNum):
	nRMSD = []
	nRMSDt = []
	rNum = np.shape(rmsd)[0]
	for row in range (1, rNum+1):
		nRMSDt.append(rmsd[row-1][0]) # append the time
		for i in range (1, slabNum+1):
			if slabAtomSurvArr[row-1][i-1] == 0:
				nRMSDt.append(0.)
			else:
				nRMSDt.append(rmsd[row-1][i]/slabAtomSurvArr[row-1][i-1])
		nRMSD.append(nRMSDt)
		nRMSDt = []
	return nRMSD

