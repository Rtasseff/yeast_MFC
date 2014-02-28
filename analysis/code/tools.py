# module for gerneal analysis tools
import numpy as np

	

def calcMSDStats(X):
	"""Calculate the statistic for use in MSD calcualtions.
	Specifically the sum of SD (squared displacment),
	the sum of squared SD, and the number of samples for the sum.
	for use on calculating MSD over multiple trajectories.
	MSD for only this trajectory would be sumSD/num.
	and variance would be sumSqSD/num - (sumSD/num)**2.
	Assumes a uniform time seris.
	Assumes rows of X corrispond with time and columns
	with corrdinates.
	"""
	n = len(X)
	sumSD = np.zeros(n)
	sumSqSD = np.zeros(n)
	num = np.zeros(n)
	for i in range(n):
		for j in range(n-i):
			SD = np.sum(((X[j] - X[j+i])**2))
			sumSD[i] += SD
			sumSqSD[i] += SD**2
			num[i] += 1
	return(sumSD,sumSqSD,num)


