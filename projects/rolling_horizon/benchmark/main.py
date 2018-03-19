#Author: Zhicheng Zhu
#Email: zzhu3@lamar.edu
#Benchmark policy for rolling horizon comparison
#paper:
#A dynamic policy for grouping maintenance activities
#EJOR 1997, Wildeman et al.


import class_info as myClass
import time

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

def main(wScaleBound, setUpCost, crBound, ranSeed ):
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
	I = 10				#number of components
	TRoll = 20			#maximum current time
	TInt = 20        	#time interval from currentTime to 
	T = 0+20			#time horizon.End_clock(abs_time) 
	d = setUpCost		#setup cost#######################2#########

	sysParams = myClass.system_parameter(I, TInt, TRoll, d)
	#define bounds
	sysParams.set_wShapeBound([4,7])
	sysParams.set_wScaleBound(wScaleBound)###############1##########
	sysParams.set_cCrBound(crBound)			#############3##########
	sysParams.set_cPrBound([1,1])
	#init other parameters
	sysParams.w_shape_init()	
	sysParams.w_scale_init()	
	sysParams.cCr_init()	
	sysParams.cPr_init()
	sysParams.set_ranSeed(ranSeed)

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
		#print (currentTime)
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
	#sysRunning.print_data()
	#sysRunning.print_com()	
		
	return sysRunning.TotalCost



wScaleH = [9, 20]
wScaleL = [1, 8]
wScaleVector = [wScaleH]#[wScaleH, wScaleL]
dVector = [100, 5]
cCrH = [17, 27]
cCrL = [6, 16]
cCrVector = [cCrH]#[cCrH, cCrL]

counter = 0
print ("=============Benchmark=====================")

for wScale in wScaleVector:
	for d in dVector:
		for cCr in cCrVector:
			counter += 1
			print ("==============scale, d, ccr==================")
			print (wScale),
			print (d),
			print (cCr)
			counter1 = 0
			cost = []
			start_time = time.clock()
			for resLifeSeed in range(5):
				counter1 += 1
				ranSeed = resLifeSeed
				#
				#ranSeed = resLifeSeed + counter
				tmp = main(wScale, d, cCr, ranSeed)
				print ("cost[%d] = %d" %(counter1, tmp))
				cost.append(tmp)
			end_time = time.clock()
			time_ = end_time-start_time			
			print ("time             = %d" %time_)
			print ("cost")
			print (cost)
			print ("average cost")
			print (float(sum(cost))/len(cost))
			print ("=============================================")
				

	
	
	
	
	
	
	
	
	
	