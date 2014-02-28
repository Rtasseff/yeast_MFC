# python2.7 script to collect various metrics from output data
# assumes csv format, one complete file per step
# assumes format of radius,color,x dir, y dir, id, x, y, z


import numpy as np
# special module
import tools


inDir = '/Users/RTasseff/Projects/multiscale/yeast_MFC/analysis/data'
name = 'agent'
tStart = 36000
tStop = 150000
tStep = 300
outDir = '/Users/RTasseff/Projects/multiscale/yeast_MFC/analysis/output'
# time scaling factor change simulaiton units to hours, curnetly base line is 1 sec
tScale = 1./3600.



# first lets collect the time seris (this could get expensive)

# time seris dictionary
id2traj = {}
numAgents = np.array([])
# loop through data files
for i in range(tStart,tStop+tStep,tStep):
	finName = inDir+'/'+name+'_'+str(i)+'.csv'
	fin = open(finName)
	# first line is header:
	line = fin.next()
	count = 0
	# loop through all line
	for line in fin:
		count = count + 1
		data = line.strip().split(',')
		agentID = int(data[4])
		x = float(data[5])
		y = float(data[6])
		tmp = np.array([x,y])
		# check if this agent already exists
		if id2traj.has_key(agentID):
			# append to end, assuming here that this is in fact the next time step (no previous missing values)
			id2traj[agentID] = np.vstack((id2traj[agentID],np.array([x,y])))
		else:
			# create
			id2traj[agentID] = np.array([x,y])

	fin.close()
	numAgents = np.append(numAgents,count)

# number of time steps 
n = len(numAgents)

# get the mean squared displacment
sumSD = np.zeros(n)
num =  np.zeros(n)
sumSqSD =  np.zeros(n)


# loop over each agent
mMax = 0 
mMaxId = 0
for key, value in id2traj.iteritems():
	if len(value)>5:
		sumSD_tmp,sumSqSD_tmp, num_tmp = tools.calcMSDStats(value)
		m = len(sumSD_tmp)
		if m > mMax:
			mMaxId = key
			mMax = m
		num[:m] = num[:m] + num_tmp
		sumSD[:m] = sumSD[:m] + sumSD_tmp
		sumSqSD[:m] = sumSqSD[:m] + sumSqSD_tmp

msd = sumSD/num
vsd = sumSqSD/num - msd**2
errsd = np.sqrt(vsd/num)
dt = np.arange(n) * tStep * tScale




	


