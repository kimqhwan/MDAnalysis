import numpy as np 
import md.mdMath.vector

def getRowwiseAverage(inputArray):
	rowNumber = np.shape(inputArray)[0]
	columnNumber = np.shape(inputArray)[1]

	outputArray = np.zeros(columnNumber)

	for row in range (1, rowNumber+1):
		outputArray = md.mdMath.vector.getVectorSum(outputArray, inputArray[row-1])

	outputArray = md.mdMath.vector.getVectorScalarDivision(outputArray, rowNumber)

	return outputArray