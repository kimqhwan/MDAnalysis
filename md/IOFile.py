#!/usr/bin/env python
import sys
import numpy as np


def openFile(ifile):
	print "Open the i file"

	a = np.loadtxt(ifile)
	return a



def saveFile(ofile, oArray):
	print "Save the file"

	np.savetxt(ofile, oArray, fmt='%.4f')


def saveIntFile(ofile, oArray):
	print "Save the file"

	np.savetxt(ofile, oArray, fmt='%d')


def saveFileLineToLine(ofile, oLine):
	print "Save the file"



def wCoordSpecificTime(iArray, row):
	buf_iArray = []
	oArray = []

	for col in range (2, np.shape(iArray)[1]+1, 9):
		buf_iArray.append((col-1)/9+1) ## Append the molecule number
		buf_iArray.append(iArray[row-1][col-1]) ## Append the OW
		buf_iArray.append(iArray[row-1][col])
		buf_iArray.append(iArray[row-1][col+1])
		buf_iArray.append(iArray[row-1][col+2]) ## Append the HW1
		buf_iArray.append(iArray[row-1][col+3])
		buf_iArray.append(iArray[row-1][col+4])
		buf_iArray.append(iArray[row-1][col+5]) ## Append the HW2
		buf_iArray.append(iArray[row-1][col+6])
		buf_iArray.append(iArray[row-1][col+7])
		
		oArray.append(buf_iArray)
		buf_iArray = []

	return np.asarray(oArray)