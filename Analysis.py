
# Analysis Code for the Dendritic Growth, Marc Jacquart. 7 Dec. 2020

import numpy as np
import matplotlib.pyplot as plt
import time
import matplotlib.animation as animation

allFramesTab=[]
lastFrameTab=[] # {{x1,x2,x3,...,xn},{y1,y2,y3,...,yn}}
latticeSizeX=300;
latticeSizeY=300;




def lastFrameInTab(): 							# Fill the tab lastFrameTab with the elements of the last line of the text file

	txtFile=open("Simulation.txt",'r') 			# Open the data txt file in read mode
	lines=txtFile.read().splitlines()			# Split the different lines
	lastLine = lines[-1]						# Only the last is interesting here
	frameTab=[]
	tab=lastLine.split(";") 					# tab has elements of form ###,###
	for element in tab:
		frameTab.append(element.split(",")) 	# Split the x and y coordinate
	frameTab.pop() 								# Delete last element wich is the endline character
	lastFrameTab.append(frameTab)				# Result goes in lastFrameTab

def computeAllFramesTab():						# Fills the tab allFramesTab with all the frames of the simulation.txt result file (1 line = 1 print of the system, see the print frequency in the parameters text file)

	txtFile=open("Simulation.txt",'r') 			# Open the data txt file in read mode
	content=txtFile.readlines()
	for line in content:						# Read each line
		frameTab=[] 							# Create a result table for this frame, is set again to 0 each time
		tab=line.split(";") 					# tab has elements of form ###,###
		for element in tab:
			frameTab.append(element.split(",")) # Split the x and y coordinate
		frameTab.pop() 							# Delete last element wich is the endline character
		allFramesTab.append(frameTab); 			# Put the coordinates of the points of the frame in the tab





#Analysis:
def distance(r1,r2): 							# Compute the distance between atoms
	#print(float(r1[0]))
	result1 = np.sqrt((float(r1[0])-float(r2[0]))**2+(float(r1[1])-float(r2[1]))**2) #float to convert the string into double
	result2 = np.sqrt((float(r1[0])+latticeSizeX-float(r2[0]))**2+(float(r1[1])-float(r2[1]))**2) 	 	# Copy of r1 in the 8 lattice aroun the main, takes the minimum of the 9 distances			
	result3 = np.sqrt((float(r1[0])-float(r2[0]))**2+(float(r1[1])+latticeSizeY-float(r2[1]))**2) 		# Because with periodic bc, the atoms in the bottom let and top right corner are close to each other
	result4 = np.sqrt((float(r1[0])+latticeSizeX-float(r2[0]))**2+(float(r1[1])+latticeSizeY-float(r2[1]))**2) 
	result5 = np.sqrt((float(r1[0])+latticeSizeX-float(r2[0]))**2+(float(r1[1])-latticeSizeY-float(r2[1]))**2)
	result6 = np.sqrt((float(r1[0])-latticeSizeX-float(r2[0]))**2+(float(r1[1])-float(r2[1]))**2)
	result7 = np.sqrt((float(r1[0])-float(r2[0]))**2+(float(r1[1])-latticeSizeY-float(r2[1]))**2)
	result8 = np.sqrt((float(r1[0])-latticeSizeX-float(r2[0]))**2+(float(r1[1])+latticeSizeY-float(r2[1]))**2)
	result9 = np.sqrt((float(r1[0])-latticeSizeX-float(r2[0]))**2+(float(r1[1])-latticeSizeY-float(r2[1]))**2)
	result = np.amin([result1,result2,result3,result4,result5,result6,result7,result8,result9]) 	# Easy but slow way to take the minimal distance taking periodic condition into account
	#print(result)
	return 2.8*result 							# Convert in angstrom at the same time

def plotDistance():								# Compute the distance between each pair of atoms, not used in the end, lattice not big enough to have statistically significant results
	lastFrameInTab()
	lastTab=lastFrameTab[0]
	distanceTab=[]								# To fill in the results
	i=0
	N=len(lastTab)
	#N=100 #Smaller tests
	while i<N:
		j=i
		if (i%300==0):
			print(str(i)+" over "+str(N))		# To keep track of the histogram loading, str to convert int into string to concatenate
		
		while j< N: 							# Part 0 of enumerate contains the index
			#print(lastTab[i])
			d=distance(lastTab[i],lastTab[j])
			if d>0: 							# Possibility to chose d>0 to reduce the number of atoms to plot (reduce histogram filling time)
				distanceTab.append(d)
			j=j+1
		i=i+1
	print("Plotting in histogram the results...")
	print("Number of data to plot: "+str(len(distanceTab)))
	figDist=plt.figure()

	n, bins, patches=plt.hist(distanceTab, bins=212)		# Put everything in an histogram

	meanBins=[]												# Histogram gives values of bin border, manually computes the mean value of each bin
	k=0
	while k<(len(bins)-1):
		meanBins.append((bins[k+1]+bins[k])/2.0)
		k=k+1

	

	# Delete linear part:
	maxValue = max(n)										# Find the maximum value (to cut the linear part)
	maxIndex = np.argmax(n)									# And its index
	i=0
	while i<len(n):
		if i%10==0:
			print (i)
		n[i]=n[i]-(maxValue/meanBins[maxIndex])*meanBins[i]	# Substract the linear slope going trough (0,0) and the maximum to better see the peak
		i+=1

	plt.xlabel("distance [Å]")
	plt.ylabel("Counts [-]")
	plt.show()												# Plot before the cut of linear part
	
	plot2=plt.figure()
	plt.plot(meanBins, n,'k.')
	plt.xlabel("distance [Å]")
	plt.ylabel("Counts [-]")
	plt.show()												# Plot after the cut of linear part

	


def plotEvolution(): 							# Old way to plot the evolution of each print of the system with succesives plots, better with gif
	computeAllFramesTab()
	for frame in allFramesTab:
		xData= [(2.8*float(r[0])) for r in frame]
		yData= [(2.8*float(r[1])) for r in frame]

		plt.plot(xData,yData,"k.",markersize=80)
		plt.show()
		#time.sleep(0.01)
		plt.close()



def plotLast():									# Plot only the last configuration of the system
	lastFrameInTab()
	print("Number of particles:")
	print(len(lastFrameTab[0]))					# For info: Number of particles in the plot (#Atoms in the cluster if only 1 printed)


	print("Extracting data from lastFrameTab")
	xData= [(2.8*float(r[0])) for r in lastFrameTab[0]]
	yData= [(2.8*float(r[1])) for r in lastFrameTab[0]]
	print("Plotting last configuration")
	fig=plt.figure()
	fig.set_size_inches(10,8)
	plt.plot(xData,yData,"k.",markersize=2.1)
	plt.xlabel("X [Å]")
	plt.ylabel("Y [Å]")
	plt.show()




def plotGif():									# Plot a gif animation with all sucessive print: see the deposition in time
	computeAllFramesTab()
	#Animation
	fig=plt.figure()
	fig.set_size_inches(13,10)
	ax=plt.axes(xlim=(0,latticeSizeX),ylim=(0,latticeSizeY*0.866))

	graph, = ax.plot([], [], "k.",markersize=3) 

	# initialization function 
	def init(): 
		# creating an empty plot/frame 
		graph.set_data([], []) 
		return graph, 

	# lists to store x and y axis points 
	#xdata, ydata = [], [] 

	def animate(i):
		if(i%200==0):							# Keep track of the animation rendering time
			print(i)
		xData, yData = [], [] 
		xData= [float(r[0]) for r in allFramesTab[i]]
		yData= [float(r[1]) for r in allFramesTab[i]]
		#plt.plot(xData,yData,"k.")
		graph.set_data(xData, yData) 
		
		return graph,
	print("Analysis code:")
	anim = animation.FuncAnimation(fig, animate, init_func=init, 
								frames=len(allFramesTab), interval=1, blit=True) #This command makes the gif file

	#fig.show()
	anim.save('graph.gif', dpi=80, writer='imagemagick',fps=5)

def centerMass():								# Compute the center of mass of the cluster
	lastFrameInTab()
	N=len(lastFrameTab[0])
	print("N= ")
	print(N)
	xTot=0.0
	yTot=0.0
	for element in lastFrameTab[0]:				# Sum of the componant normalized with the number of atoms
		xTot+=float(element[0])
		yTot+=float(element[1])
	xCM=xTot/N
	yCM=yTot/N
	return xCM,yCM


def plotTime():									# For performance analysis: plot the totoal time and the time between two prints
	
	timeData=[]
	timeDiff=[]
	timeFile=open("Time.txt",'r')
	timeContent=timeFile.readlines()
	for line in timeContent:
		
		timeData.append(int(line.rstrip())) 	# rstrip() to erase eol character, int to transform string into int
		print(int(line.rstrip()))
				

	m=1
	while m<len(timeData):
		timeDiff.append(timeData[m]-timeData[m-1])
		m+=1
	xTime=np.linspace(0.0,100.0,len(timeData))
	xTimeDiff=np.linspace(0.0,100.0,len(timeDiff))
	figtime=plt.figure()
	plt.plot(xTime,timeData, 'k.')
	plt.xlabel("Physical time simulated [s]")
	plt.ylabel("Total computation time [s]")
	plt.show()
	figTimeDiff=plt.figure()
	plt.plot(xTimeDiff,timeDiff, 'k.')
	plt.xlabel("Physical time simulated [s]")
	plt.ylabel("Time to compute 5000 system evolutions [s]")
	plt.show()


def fractalDimension():							# Compute the fractal dimension D using N~r^D => linear fit on loglog plot
	xCM,yCM=centerMass()
	print(xCM)
	print(yCM)
	lastFrameInTab()
	#print(lastFrameTab)
	nSteps=60
	rMax=500
	print("Geometric center:")
	centerX=2.8*(latticeSizeX/3.5+0.5)+0.5*(latticeSizeX/2+0.5) # Center of tetrahedra =/= Center of mass (must be used)
	centerY=2.8*(latticeSizeX/2+0.5)*0.866
	print(centerX)
	print(centerY)

	rTab=np.logspace(0,np.log10(rMax), num=nSteps) 			# Power of 10, we want logspace to seem linear on loglog plot
	nTab=[]
	for r in rTab:
		n=0
		i=0
		while i<len(lastFrameTab[0]):
			if ((float(lastFrameTab[0][i][0])-xCM)**2+(float(lastFrameTab[0][i][1])-yCM)**2)<r**2: #Compute number of atoms in a circle of radius r around CM
				n=n+1
			i=i+1
		nTab.append(n)


	#for the "linear" fit on loglog plot:
	#1: change r and N into np.arrays:
	nTabNp = np.asarray(nTab, dtype=float)
	rTabNp = np.asarray(rTab, dtype=float)
	#2:compute log10:
	nTab10 = np.log10(nTabNp)
	rTab10 = np.log10(rTabNp)
	#3:DO the fit:
	rTab10Crop=rTab10[12:48]
	nTab10Crop=nTab10[12:48]
	
	coefFit,cov = np.polyfit(rTab10Crop, nTab10Crop, 1, cov=True)
	fitLine = np.poly1d(coefFit)
	print(fitLine)
	print (np.sqrt(np.diag(cov)))
	#Y=ax+b
	yN = fitLine(rTab10Crop)#y component of the fit
	#Plot the result
	
	figFractalDimension, leg = plt.subplots()
	plt.plot(rTab10,nTab10,'k.')
	plt.plot(rTab10Crop, yN, label='Linear fit: $log_{10}(N)=(1.71 \pm 0.01) log_{10}(r)+(0.33  \pm 0.01)$') 
	plt.xlabel("$log_{10}(r)$ [Å]")
	plt.ylabel("$log_{10}(N)$ [-]")
	legend=leg.legend(frameon=False)
	plt.show()

# Command to comment/uncomment to select the analysis

#fractalDimension()		# Compute the fractal dimension D using N~r^D => linear fit on loglog plot
#plotTime()				# For performance analysis: plot the totoal time and the time between two prints
#plotGif()				# Plot a gif animation with all sucessive print: see the deposition in time
#plotLast()				# Plot only the last configuration of the system
plotDistance()			# Compute the distance between each pair of atoms, not used in the end, lattice not big enough to have statistically significant results