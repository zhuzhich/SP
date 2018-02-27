#Author: Zhicheng Zhu
#Email: zzhu3@lamar.edu
#Benchmark policy for rolling horizon comparison
#paper:
#A dynamic policy for grouping maintenance activities
#EJOR 1997, Wildeman et al.


import class_info as myClass


#decide component next MX time.
def arrange_mx(sysRunning):
	currentTime = sysRunning.clock
	endTime = sysRunning.endClock
	#step 1: find out available moving window. t1 = 0	
	for i in range(sysRunning.comRunningNum):
		#reset nextMxTime to scheduled time point
		#get the widow
		#set failed component's nextMxTime to currentTime
		sysRunning.comRunning[i].cal_movingWindow(currentTime,endTime)
	
	#sort runnning componenet based on their nextMxTime
	sysRunning.sort_running_component()
	
	#step 2: find the optimal group, maximize total saving
	sysRunning.find_group()
		
#
#
######################
#start from here
######################	
#

###
#Two class objects:
#1. sysRunning <- component_running <- component_parameters
#2. sysParams
###
###
#######
#Step 1: define sysParams: static system info
#######
I = 4				#number of components
TRoll = 20			#maximum current time
TInt = 20        #time interval from currentTime to 
T = 0+20			#time horizon.End_clock(abs_time) 
d = 5				#setup cost

sysParams = myClass.system_parameter(I, TInt, TRoll, d)
#define bounds
sysParams.set_wShapeBound([4,7])
sysParams.set_wScaleBound([4,20])
sysParams.set_cCrBound([6,16])
sysParams.set_cPrBound([1,1])
#init other parameters
sysParams.w_shape_init()	
sysParams.w_scale_init()	
sysParams.cCr_init()	
sysParams.cPr_init()	

########
#Step 2: define sysRunning: system grouping info
########
sysRunning = myClass.system_running(sysParams.TInt, sysParams.TRoll, sysParams.I)

########
#Step 3: 
#component parameters, component running info
#attach sysParams, sysRunning respectively
########
for i in range(sysParams.I):
	com_params = myClass.component_parameters(i, sysParams)
	com = myClass.component_running(i, com_params)
	sysRunning.add_com(com)
	sysRunning.add_comRunning(com)

######
#Step 4: begin rolling horizon
######	

arrange_mx(sysRunning)

for currentTime in range(1, sysParams.TRoll+1):
	TInt += -1						#####
	sysParams.set_TInt(TInt)
	sysRunning.set_time(currentTime, TInt)
	print (currentTime)
	#step 1: get residual life
	for i in range(sysRunning.comRunningNum):
		#use currentTime and index i to generate residual life
		sysRunning.comRunning[i].generate_residualLife(currentTime)	

	
	#step 2: check whether there is any component needs to be mx at current time
	sysRunning.mx()

	#Step 3: reschedule the mx plan	
	if sysRunning.needSchedule == True or \
		sysRunning.newOpp == True:
		arrange_mx(sysRunning)
		
	#step 4: re-mx if new opportunity emerges
	if sysRunning.newOpp == True:	
		sysRunning.mx()
	#step 5: after re-mx. check whether there is any component is replaced.
		if sysRunning.needSchedule == True:
			arrange_mx(sysRunning)
		#it shouldn't be any new opp any more
		if sysRunning.newOpp == False:
			print ("biggggg error!!!")

####		
sysRunning.cal_TotalCost()
sysRunning.print_data()
sysRunning.print_com()	
	
	
	
	
	
	
	
	
	
	