'''
hB.py
	calculate the total hydrogen bond number, dAProfile, (HB)n distribution.



checkDonorAcceptorDistance(dArray, aArray, length)
	: Check between donor-acceptor, return true-false.

checkDonorAcceptorDistance(dHArray, dAArray, angle)
	: Check angle between donor-hydrogen vector & donor-acceptor vector. return true-false.

getHB(iArray, time, length=0.35, angle=30.0) 
	: Get hydrogen bond profile from the coordinate file.
	Input -> coordinate file exported from the md.IOFile.waterCoordinateSpecificTime().
	Output -> hBNum, donor-acceptor profile.

getHBNumDistribution(dAProfile, wMolNum)
	: Get (HB)n distriution from the donor-acceptor profile.

getHBDenDistribution(hBNumEachMol, oCoord, zI, zF, zStep)
	: Get hydrogen bond density distribution as a function of z-coordinate.
	z-coordinate is set at this function with variables zI, zF, zStep.
'''


import numpy as np
import md.IOFile
import md.mdMath.vector
import md.analysis.density




def checkDonorAcceptorDistance(dArray, aArray, length):
	distVec = dArray - aArray
	if md.mdMath.vector.getVecLength(distVec) <= length:
		return 1
	else:
		return 0



def checkHydrogenDonorAcceptorAngle(dHArray, dAArray, angle):
	ang = md.mdMath.vector.getVecAngle(dHArray, dAArray)
	if ang <= angle:
		return 1
	else:
		return 0




def getHB(iArray, time, length=0.35, angle=30.0):
	wMolNum = np.shape(iArray)[0]
	hBNum = 0
	dAProfile = []

	dAProfile.append(int(time))

	for i in range(1, wMolNum+1):

		if i % 100 == 0:
			print i

		dArray = np.array([iArray[i-1][1], iArray[i-1][2], iArray[i-1][3]])
		h1Array = np.array([iArray[i-1][4], iArray[i-1][5], iArray[i-1][6]])
		h2Array = np.array([iArray[i-1][7], iArray[i-1][8], iArray[i-1][9]])
		dH1Array = h1Array - dArray
		dH2Array = h2Array - dArray
		for j in range(i+1, wMolNum+1):
			if i == j:
				check = 0
			else:
				aArray = np.array([iArray[j-1][1], iArray[j-1][2], iArray[j-1][3]])
				h3Array = np.array([iArray[j-1][4], iArray[j-1][5], iArray[j-1][6]])
				h4Array = np.array([iArray[j-1][7], iArray[j-1][8], iArray[j-1][9]])
				aH3Array = h3Array - aArray
				aH4Array = h4Array - aArray
				dAArray = aArray - dArray
				aDArray = dArray - aArray

				if checkDonorAcceptorDistance(dArray, aArray, length) == 0:
					check=0
				else:
					hBondDH1A = checkHydrogenDonorAcceptorAngle(dH1Array, dAArray, angle)
					hBondDH2A = checkHydrogenDonorAcceptorAngle(dH2Array, dAArray, angle)
					hBondAH3D = checkHydrogenDonorAcceptorAngle(aH3Array, aDArray, angle)
					hBondAH4D = checkHydrogenDonorAcceptorAngle(aH4Array, aDArray, angle)
					hBNum += (hBondDH1A+hBondDH2A+hBondAH3D+hBondAH4D)*2
					for k in range (1, hBondDH1A+1):
						dAProfile.append(int(i))
						dAProfile.append(int(j))
					for k in range (1, hBondDH2A+1):
						dAProfile.append(int(i))
						dAProfile.append(int(j))
					for k in range (1, hBondAH3D+1):
						dAProfile.append(int(j))
						dAProfile.append(int(i))
					for k in range (1, hBondAH4D+1):
						dAProfile.append(int(j))
						dAProfile.append(int(i))

	return hBNum, dAProfile, wMolNum




def getHBNumDistribution(dAProfile, wMolNum):
	hBNumDistribution = np.zeros(8)
	hBNumEachMol = np.zeros(wMolNum)
	hBNum = len(dAProfile)-1

	for i in range (2, hBNum+2):
		donorIndex = dAProfile[i-1]
		hBNumEachMol[donorIndex-1] += 1

	for i in range (1, wMolNum+1):
		if hBNumEachMol[i-1] >= 8:
			hBNumDistribution[len(hBNumDistribution)-1] += 1
		else:
			hBNumDistribution[hBNumEachMol[i-1]] += 1

	#np.insert(hBNumDistribution, 0, time)

	return hBNumDistribution, hBNumEachMol



def getHBDenDistribution(hBNumEachMol, oCoord, zI, zF, zStep):
	rowNum1 = np.shape(hBNumEachMol)[0]
	molNum1 = np.shape(hBNumEachMol)[1]
	rowNum2 = np.shape(oCoord)[0]
	molNum2 = int((np.shape(oCoord)[1]-1)/3)
	if rowNum1 != rowNum2 or molNum1 != molNum2:
		print "Hey, molNum or timeNumber is different for two file, hBing file and oCoord"
		exit(1)
	rowNum = rowNum1

	zCoordSlabArray, slabNum = md.analysis.density.getZcoordSlabArray(zI, zF, zStep)
	slabhBDenArray = np.zeros(slabNum)


	slabhBNumDistributionArray = []
	for i in range (1, slabNum+1):
		slabhBNumDistributionArray.append(np.zeros(8))
	slabhBNumDistributionArray = np.asarray(slabhBNumDistributionArray)



	for row in range (1, rowNum+1):
		buf_slabhBDenArray = []

		oSlabAtomIndexArray, oSlabAtomNumberArray = md.analysis.density.getSlabAtomIndexList(row, oCoord, zCoordSlabArray)

		slabhBNumArray = np.zeros(slabNum)

		for i in range (1, slabNum+1):
			for j in range (1, len(oSlabAtomIndexArray[i-1])+1):
				slabhBNumArray[i-1] += hBNumEachMol[row-1][oSlabAtomIndexArray[i-1][j-1]-1]
				slabhBNumDistributionArray[i-1][hBNumEachMol[row-1][oSlabAtomIndexArray[i-1][j-1]-1]] += 1

		for i in range (1, slabNum+1):
			if oSlabAtomNumberArray[i-1] == 0:
				buf_slabhBDenArray.append(0)
			else:
				buf_slabhBDenArray.append(slabhBNumArray[i-1]/oSlabAtomNumberArray[i-1])

		for i in range (1, slabNum+1):
			slabhBDenArray[i-1] += buf_slabhBDenArray[i-1]

	for i in range (1, slabNum+1):
		slabhBDenArray[i-1] = slabhBDenArray[i-1]/rowNum

	for i in range (1, slabNum+1):
		slabhBNumDistributionArray[i-1] = md.mdMath.vector.getVecFraction(slabhBNumDistributionArray[i-1])


	return slabhBDenArray, slabhBNumDistributionArray 



