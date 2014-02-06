'''
correlation.py
	this module contains the correlation analysis function.


===========================================================
cal2DDipoleCorr(iArray, tI)
	: calculate 2D Dipole correlation

cal2DPointDipCorr(icorrd, dVec, pointAtomIndexArr, pointAtomNum, pointAtomNumI, rowI, row, time)
	: calculate 2D Dipole correlation within certain range of corrdinates
	Return
		List(time, Corr)

cal3DDipoleCorr(iArray, tI)
	: calculate 3D Dipole correlation

norm2DPointCorr(corr, pointAtomSurvArr)
	: When Correlation is calculated, this results is from the N_i
	So we need to normalize this data using N_t/N_i, survival probaility.
	Input
		t corr1 
		.. ..   
						&
		List[P1 P2 P3 ... PN]

	Output
		t NC1 
		.. ..    
===========================================================
'''

import numpy as np
import md.analysis.slab

def cal2DDipoleCorr(iArray, tI, lastTau):
	print "Calculate Dipole Correlation"

	oArray = []

	shape = np.shape(iArray)
	tNum = shape[0]
	atomNum = shape[1]-1
	wNum = int((shape[1]-1)/2)

	dt = iArray[1][0]-iArray[0][0]
	tF = iArray[tNum-1]


	buf_sum = 0.00
	buf_array = []

	
	for tau in range (0, lastTau+1):
		print tau
		buf_array.append(tau*dt)
		tauMax = (tNum-1)-tau # complex!!!!!!!!

		for tau0 in range (0, tauMax+1):

			# Summation of the molecule
			for col in range (2, atomNum+1, 2):

				p_x0=iArray[tau0][col-1]
				p_y0=iArray[tau0][col]
				p_xt=iArray[tau0+tau][col-1]
				p_yt=iArray[tau0+tau][col]
				buf_sum += p_x0*p_xt + p_y0*p_yt	

			averageSum = buf_sum/wNum
		averageSum = averageSum/tauMax
		buf_array.append(averageSum)

		oArray.append(buf_array)

		buf_sum = 0.00
		buf_array = []

	return oArray



def cal2DPointDipCorr(iCoord, dVec, zCoord, zDelta, tI, tF, lastTau):
	print "Calculate Dipole Correlation"

	corr = []
	rowCorr = []

	tNum = np.shape(iCoord)[0]
	atomNum = (np.shape(iCoord)[1]-1)/3

	dt = iCoord[1][0]-iCoord[0][0]

	
	for tau in range (0, int(lastTau)+1):
		print tau
		tauCorr = 0.
		rowCorr = []
		localCorr = 0.

		rowCorr.append(tau*dt)
		tauMax = (tNum-1)-tau # complex!!!!!!!!

		for tau0 in range (0, tauMax+1):

			t0t = tI + (tau0+tau)*dt

			pointAtomIndexArr = \
      	  		md.analysis.slab.getPointAtomIndexArr(iCoord, zCoord, zDelta, atomNum, t0t, dt)
			pointAtomNum = \
      	  		md.analysis.slab.getPointAtomNum(pointAtomIndexArr)

			rowI = tau0+1
			row = tau0+tau+1

			# Summation of the molecule
			for j in range (1, pointAtomNum+1):
				dVecColX = 2*pointAtomIndexArr[j-1]
				dVecColY = 2*pointAtomIndexArr[j-1]+1
				dotDip = dVec[row-1][dVecColX-1]*dVec[rowI-1][dVecColX-1]  + \
						dVec[row-1][dVecColY-1]*dVec[rowI-1][dVecColY-1]
				normDip = \
					np.sqrt(dVec[row-1][dVecColX-1]**2 + dVec[row-1][dVecColY-1]**2) * \
					np.sqrt(dVec[rowI-1][dVecColX-1]**2 + dVec[rowI-1][dVecColY-1]**2)
				if normDip == 0.:
					localCorr += 0.
				else:
					localCorr += dotDip / normDip
			if pointAtomNum == 0:
				localCorr += 0.
			else:
				localCorr /= pointAtomNum
			tauCorr += localCorr
		tauCorr /= tauMax
		rowCorr.append(tauCorr)
		print rowCorr
		corr.append(rowCorr)

	return corr



def cal3DDipoleCorr(iArray, tI):
	print "Calculate Dipole Correlation"

	oArray = np.array([])

	shape = np.shape(iArray)
	tNum = shape[0]
	atomNum = shape[1]-1
	wNum = int((shape[1]-1)/3)

	buf_sum = 0.00
	buf_array = np.array([])

	for row in range (1, tNum):
		if iArray[row][0] % 100 == 0:
			print iArray[row][0]

		buf_array.append(iArray[row][0] - tI)

		for col in range (2, atomNum+1, 3):
			p_x0=iArray[0][col-1]
			p_y0=iArray[0][col]
			p_z0=iArray[0][col+1]
			p_xt=iArray[row-1][col-1]
			p_yt=iArray[row-1][col]
			p_zt=iArray[row-1][col+1]
			buf_sum += p_x0*p_xt + p_y0*p_yt + p_z0*p_zt		

		averageSum = buf_sum/wNum
		buf_array.append(averageSum)

		oArray.append(buf_array)

		buf_sum = 0.00
		buf_array = np.array([])

	return oArray


def norm2DPointCorr(corr, pointAtomSurvArr):
	ncorr = []
	ncorrt = []
	ncorr.append([0.,1.]) # at t=0
	rNum = np.shape(corr)[0]
	for row in range (2, rNum+1):
		ncorrt.append(corr[row-1][0]) # append the time
		if pointAtomSurvArr[row-1] == 0.:
			ncorrt.append(0.)
		else:
			ncorrt.append(corr[row-1][1]/pointAtomSurvArr[row-1])
		ncorr.append(ncorrt)
		ncorrt = []
	return ncorr
