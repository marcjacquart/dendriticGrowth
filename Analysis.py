import numpy as np
import matplotlib.pyplot as plt
import time
import matplotlib.animation as animation

allFramesTab=[]
lastFrameTab=[]
latticeSizeX=250;
latticeSizeY=250;




def lastFrameInTab():
	txtFile=open("Simulation.txt",'r') #Open the data txt file in read mode
	lines=txtFile.read().splitlines()
	lastLine = lines[-1]
	frameTab=[]
	tab=lastLine.split(";") #tab has elements of form ###,###
	for element in tab:
		frameTab.append(element.split(",")) #Split the x and y coordinate
	frameTab.pop() #delete last element wich is the endline character
	lastFrameTab.append(frameTab)

def computeAllFramesTab():

	txtFile=open("Simulation.txt",'r') #Open the data txt file in read mode
	content=txtFile.readlines()
	for line in content:
		frameTab=[] #Create a result table for this frame, is set again to 0 each time
		tab=line.split(";") #tab has elements of form ###,###
		for element in tab:
			frameTab.append(element.split(",")) #Split the x and y coordinate
		frameTab.pop() #delete last element wich is the endline character
		allFramesTab.append(frameTab); #put the coordinates of the points of the frame in the tab




#Analysis:
def distance(r1,r2):
	#print(float(r1[0]))
	result = np.sqrt((float(r1[0])-float(r2[0]))**2+(float(r1[1])-float(r2[1]))**2) #float to convert the string into double
	#print(result)
	return result

def plotDistance():
	lastFrameInTab()
	lastTab=lastFrameTab[0]
	distanceTab=[]
	i=0
	N=len(lastTab)
	#N=100 #Smaller tests
	while i<N:
		j=i
		if (i%300==0):
			print(str(i)+" over "+str(N))#To keep track of the histogram loading, str to convert int into string to concatenate
		
		while j< N: #part 0 of enumerate contains the index
			#print(lastTab[i])
			d=distance(lastTab[i],lastTab[j])
			if d>0: #Chose to reduce number to plot
				distanceTab.append(d)
			j=j+1
		i=i+1
	print("Plotting in histogram the results...")
	print("Number of data to plot: "+str(len(distanceTab)))
	figDist=plt.figure()
	plt.hist(distanceTab, bins=500)
	plt.show()

def plotEvolution(): #Old way with succesives plots, better with gif
	computeAllFramesTab()
	for frame in allFramesTab:
		xData= [float(r[0]) for r in frame]
		yData= [float(r[1]) for r in frame]

		plt.plot(xData,yData,"k.",markersize=80)
		plt.show()
		#time.sleep(0.01)
		plt.close()



def plotLast():
	lastFrameInTab()
	print("Extracting data from lastFrameTab")
	xData= [float(r[0]) for r in lastFrameTab[0]]
	yData= [float(r[1]) for r in lastFrameTab[0]]
	print("Plotting last configuration")
	fig=plt.figure()
	fig.set_size_inches(10,8)
	plt.plot(xData,yData,"k.",markersize=2)
	plt.show()




def plotGif():
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
		if(i%200==0):
			print(i)
		xData, yData = [], [] 
		xData= [float(r[0]) for r in allFramesTab[i]]
		yData= [float(r[1]) for r in allFramesTab[i]]
		#plt.plot(xData,yData,"k.")
		graph.set_data(xData, yData) 
		
		return graph,
	print("Analysis code:")
	anim = animation.FuncAnimation(fig, animate, init_func=init, 
								frames=len(allFramesTab), interval=1, blit=True)

	#fig.show()
	anim.save('graph.gif', dpi=80, writer='imagemagick',fps=5)

#plotGif()
plotLast()
plotDistance()