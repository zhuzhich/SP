import math
import random
import numpy as np
from scipy.integrate import quad
class system_parameter():
	
	def set_TInt(self,t):
		self.TInt = t
		
	# weibull parameters
	def w_shape_init(self):
		I = self.I
		bound = self.wShapeBound
		for i in range(0, I):
			random.seed((i+1)*10) ###control the seed     
			temp = random.uniform(bound[0], bound[1])
			self.wShape.append(round(temp,1))

	def w_scale_init(self):
		I = self.I
		bound = self.wScaleBound
		for i in range(0, I):
			random.seed((i+1)*20) ###control the seed     
			temp = random.uniform(bound[0], bound[1])
			self.wScale.append(round(temp,1))	
	
	#cost parameters 			
	def cCr_init(self):
		I = self.I
		bound = self.cCrBound
		for i in range(0, I):
			random.seed((i+1)*30) ###control the seed     
			temp = random.uniform(bound[0], bound[1])
			self.cCr.append(round(temp,1))
			

	def cPr_init(self):
		I = self.I
		bound = self.cPrBound
		for i in range(0, I):
			random.seed((i+1)*40) ###control the seed     
			temp = random.uniform(bound[0], bound[1])
			self.cPr.append(round(temp,1))
	
	#some setting function
	def set_wShapeBound(self, bound):
		self.wShapeBound = bound
	def set_wScaleBound(self, bound):
		self.wScaleBound = bound
	def set_cCrBound(self, bound):
		self.cCrBound = bound
	def set_cPrBound(self, bound):
		self.cPrBound = bound
	def set_ranSeed(self, rs):
		self.ranSeed = rs
	def __init__(self,I, TInt, TRoll, d):
		#system parameters.
		#I, TInt, TRoll, W, d can be changed for "for" loop
		#Others are fixed
		
		#not fixed parameters
		self.I = I				#number of components
		self.TInt = TInt		#time interval. 
		self.T = TInt			#time horizon. currentTime + TInt = T
		self.TRoll = TRoll		#time points for rolling
		self.d = d				#setup cost

		#bounds
		self.wShapeBound = []
		self.wScaleBound = []
		self.cCrBound = []
		self.cPrBound = []
		
		#parameters
		self.wShape = []
		self.wScale = []
		self.cCr = []
		self.cPr = []
		
		#extra random seed for residual life
		self.ranSeed = 0

#########################################################################	
class component_parameters():
	def weibull_cdf(self, t):
		shape = float(self.wShape)
		scale = float(self.wScale)
		t = float(t)
		tmp = (t/scale)**shape
		part1 = math.exp(-tmp)
		return 1 - part1
	def weibull_pdf(self, t):
		shape = float(self.wShape)
		scale = float(self.wScale)
		t = float(t)
		part1 = shape/scale
		part2 = (t/scale)**(shape-1)
		tmp = (t/scale)**shape
		part3 = math.exp(-tmp)
		return part1*part2*part3
	def weibull_int(self, t):
		t = float(t)
		return t*self.weibull_pdf(t)
	def cost(self, t):
		fx = self.weibull_cdf(t)
		result = self.cCr*fx + self.cPr*(1-fx)
		return result							#return cost
	def cycle_length(self, t):
		x = t									#use x here, don't wanna make confusion.
		rx = 1 - self.weibull_cdf(x)
		#integrate from 0 to x, tf(t)dt.
		part1 = quad(self.weibull_int, 0, x)[0]
		part2 = rx*x
		result = part1 + part2
		return result	
	def cost_rate(self, t):
		cost = self.cost(t)
		cyc = self.cycle_length(t)
		result = float(cost)/cyc
		return result
	#h() function for component i	
	def func_h1(self, d_t):
		t = self.optInterval
		T = self.T
		if (t+d_t) > T or (t-d_t) < 0:
			return np.inf
		result = self.cost(t+d_t) + self.cost(t-d_t) - 2*self.cost(t)
		return result
	def find_optInterval(self):
		lb = 0.05
		ub = 30
		eps = 0.5
		a = lb
		last_res = self.cost_rate(a)
		while (1):
			b = a + eps
			new_res = self.cost_rate (b)
			#print (a,b)
			#print (last_res, new_res)
			if new_res <= last_res:	#improving
				last_res = new_res
				a = b
			else:					#worsing
				tmp = (a+b)/2
				break
			if a > ub:
				tmp = a
				break
		result = round(tmp)
		return result
	#calculate maximum moveable window for component i
	#h() is the extra cost by moving delta_t
	#d is the potential saving for moving
	#if extra cost is greater than d, means no need to move the component.
	def cal_deltaT(self):
		bound = 30
		for d_t in range(bound):
			if self.func_h1(d_t) - self.d >= 0:
				break
		result = d_t
		return result
#######	
	def __init__(self, i, sysInfo):
		self.index = i						#index
		self.wShape = sysInfo.wShape[i]		#shape parameter
		self.wScale = sysInfo.wScale[i]		#scale parameter
		self.cCr = sysInfo.cCr[i]			#cost of cr
		self.cPr = sysInfo.cPr[i]			#cost of pr
		self.d = sysInfo.d					#setup cost
		self.T = sysInfo.T
		self.ranSeed = sysInfo.ranSeed
		#optimal age-based mx interval
		self.optInterval = self.find_optInterval()
		#interval of moving window, which has a cost-saving > 0
		self.deltaT = self.cal_deltaT()
	
####
#class of component running information
##		
class component_running():
	
	def set_lastRepTime(self, t):
		self.lastRepTime = t
	def	add_repHistory(self, rep):
		self.repHistory.append(rep)
	def set_nextMxTime(self, t):
		self.nextMxTime = t
	def set_mxSch(self, t):
		self.mxSch = t	
	def reset_nextMxTime(self):
		self.set_mxSch(self.params.optInterval + max(self.lastRepTime, 0))
		self.set_nextMxTime(self.mxSch)
	def cal_movingWindow(self, currentTime, T):
		#reset nextMxTime to scheduled time
		self.reset_nextMxTime()	
		if self.residualLife[-1] == 0 and currentTime > self.lastRepTime:
			bound = [currentTime, currentTime]
			self.set_nextMxTime(currentTime)
		else:
			if self.mxSch > T:
				bound = [self.mxSch, self.mxSch]
			else:
				ub = self.nextMxTime + min(self.params.deltaT, self.params.optInterval - 1)
				ub = min(T, ub)
				lb = max(self.nextMxTime - self.params.deltaT, currentTime, self.lastRepTime + 1)
				if ub < lb:
					bound = [self.nextMxTime, self.nextMxTime]
				else:
					bound = [lb, ub]
		self.movingWindow = bound
	def generate_residualLife(self, currentTime):
		wShape = self.params.wShape
		wScale = self.params.wScale
		extraSeed = self.params.ranSeed
		ranSeed = currentTime+self.index + (extraSeed + 1)**2
		random.seed(ranSeed)
		ran_num = random.uniform(0,1)
		part1 = math.log(ran_num)
		s_inv = 1.0/wShape
		surv_time = currentTime - max(0,self.lastRepTime)
		part2 = (surv_time/wScale)**wShape
		part3 = part2 - part1
		LT = round((part3**s_inv)*wScale) - surv_time					
		res = max(0,LT)
		self.residualLife.append(res)
	def	replace(self):
		#replace itself
		if self.residualLife[-1] == 0:
			typeMx = 2
		else:
			typeMx = 1
		currentTime = self.nextMxTime
		repInfo = [currentTime, typeMx]
		self.add_repHistory(repInfo)	
		self.set_lastRepTime(currentTime)
		self.reset_nextMxTime()
	def print_data(self):
		print ("index         = "),
		print (self.index )
		print ("lastRepTime   = "),
		print (self.lastRepTime )
		print ("residualLife  = "),
		print (self.residualLife )
		print ("nextMxTime    = "),
		print (self.nextMxTime )
		print ("mxSch         = "),
		print (self.mxSch )
		print ("movingWindow  = "),
		print (self.movingWindow )
		print ("repHistory    = "),
		print (self.repHistory )
		print ("optInterval   = "),
		print (self.params.optInterval )
		print ("deltaT        = "),
		print (self.params.deltaT )		
	def __init__(self, i, params):
		self.index = i				#index
		self.params = params		#component parameters. class component_parameters
		
		self.lastRepTime = -1		#lastRepTime
		self.residualLife = [-1]	#residual life time
		self.nextMxTime = 0			#next replacement time. actual time
		self.reset_nextMxTime()	
		self.mxSch = self.params.optInterval	#next mx scheduled time
		self.movingWindow = []		#lb&ub of moving time. abs time
		self.repHistory = []		#[[time, typeOfMx],...,]: time:abs time, type of mx: 1 for PR, 2 for CR

###################################################################
	
class system_running():

	def add_com(self,com):
		self.com.append(com)
	def add_comRunning(self,com):
		self.comRunning.append(com)
		self.comRunningNum += 1
	def remove_comRunning(self, idx):
		self.comRunning.remove(idx)
		self.comRunningNum += -1
	def set_startClock(self,t):
		self.startClock = t
	def set_time(self, currentTime, TInt):
		self.clock = currentTime
		self.TInt = TInt
		self.endClock = currentTime + TInt
	def sort_running_component(self):
		#sort all running components based on nextMxTime. ascending.
		self.comRunning.sort(key=lambda component_running: component_running.nextMxTime)
	
	#try to group from component i to component j.
	#return [saving, time]
	def group(self, i, j):
		#step 1: find out the inter-set from i to j
		bound =	self.comRunning[j].movingWindow
		interSet = set(range(int(bound[0]),int(bound[1])+1))
		for k in range(i,j):
			bound = self.comRunning[k].movingWindow
			tmpSet = set(range(int(bound[0]),int(bound[1])+1))
			interSet = interSet & tmpSet
			if len(interSet) == 0:#empty
				result = [-np.inf, np.inf]
				return result				#[saving, time]
		optSaving = -np.inf
		optTime = np.inf
		interSet = list(interSet)
		for k in range(len(interSet)):
			saving = 0
			timeTmp = interSet[k]
			for kk in range(i, j+1):
				""" #optimal moving window is set, no need to check it. 
				if self.comRunning[kk].residualLife[-1]  == 0 and\
					interSet[k] != self.comRunning[kk].nextMxTime:
					saving = -np.inf
					break
				"""
				deltaT = abs(self.comRunning[kk].nextMxTime - interSet[k])
				saving += self.comRunning[kk].params.func_h1(deltaT)
			setUpCost = self.comRunning[kk].params.d
			saving = (j - i + 1 - 1)*setUpCost - saving	
			if saving > optSaving:
				optSaving = saving	
				optTime = timeTmp
		result = [optSaving, optTime]
		return result
			
	def find_group(self):
		#before calling this, assuming comRunning is sorted based on nextMxTime in ascending order
		TotalSaving = [[0]*self.comRunningNum for row in range(self.comRunningNum+1)]
		First = [-1]*self.comRunningNum
		FirstTime = [-1]*self.comRunningNum
		for j in range(self.comRunningNum):
			ij = range(j+1)
			ij.reverse()
			ji = ij
			First[j] = j
			FirstTime[j] = self.comRunning[j].nextMxTime
			for i in ji:
				result_ij = self.group(i,j)						#result_ij = [saving, time]
				TotalSaving[i][j] = result_ij[0]						#
				time_tmp = result_ij[1]
				if TotalSaving[i][j] + TotalSaving[0][max(i-1,0)] > TotalSaving[self.comRunningNum][j]:
					TotalSaving[self.comRunningNum][j] = TotalSaving[i][j] + TotalSaving[0][max(i-1,0)]
					First[j] = i
					FirstTime[j] = time_tmp
		self.TotalSaving[self.clock] = TotalSaving
		
		idxJ = self.comRunningNum - 1
		if idxJ < 0:
			print ("big error!!!!!!!!!!!No component")
			return
		while (idxJ >= 0):
			idxI = First[idxJ]
			timeJ = FirstTime[idxJ]
			for k in range(idxI, idxJ+1):
				self.comRunning[k].set_nextMxTime(timeJ)
			idxJ = idxI - 1
	
	#do the maintenance
	def mx(self):
		self.needSchedule = False	#need reschedule or not
		self.newOpp	= False			#is there any new opportunity or not		
		currentTime = self.clock	
		for i in range(self.comRunningNum):
			if self.comRunning[i].nextMxTime == currentTime:
				#scheduled mx
				self.comRunning[i].replace()
				self.needSchedule = True
				if self.comRunning[i].residualLife[-1] == 0:
					self.groupResult[i][currentTime] = 2
				else:
					self.groupResult[i][currentTime] = 1
			elif self.comRunning[i].residualLife == 0:
				self.newOpp = True
	
	def	cal_TotalCost(self):
		numCom = self.comRunningNum
		TRoll = self.TRoll
		cost = 0
		for j in range(TRoll):
			if sum(self.groupResult[u][j] for u in range(numCom)) >0:
				cost += self.comRunning[0].params.d
			for k in range(numCom):
				if self.groupResult[k][j] == 1:
					cost += self.comRunning[k].params.cPr
				if self.groupResult[k][j] == 2:
					cost += self.comRunning[k].params.cCr
		self.TotalCost = cost
	def print_com(self):
		#for i in range(len(self.com)):
		#	print ("=========info of com[%d]=========" %i)
		#	self.com[i].print_data()
		for i in range(len(self.comRunning)):
			print ("=========info of comRunning[%d]=========" %i)
			self.comRunning[i].print_data()	
	def print_data(self):
		print ("startClock  = " + str(self.startClock))
		print ("TimeInterval= " + str(self.TInt))
		print ("clock       = " + str(self.clock))
		print ("endClock    = " + str(self.endClock))
		print ("TotalSaving = " ),
		print (self.TotalSaving)
		print ("needSchedule= " ),
		print (self.needSchedule)
		print ("newOpp      = " ),
		print (self.newOpp)		
		print ("groupResult = " )		
		for row in range(len(self.com)):
			print (self.groupResult[row][:])
		print ("TotalCost   = " ),		
		print (self.TotalCost)		
	def __init__(self, TInt, TRoll, I):
		#absolute time 
		self.startClock = 0			#rolling starting time point
		self.clock = 0				#current running time
		self.TInt = TInt			#
		self.TRoll = TRoll			#
		self.endClock = TInt 		#maximum running time allowed
		
		#
		self.com = []				#component info
		self.comRunning = []		#component running set
		self.comRunningNum = 0
		#
		self.groupResult = [[0]*(TRoll+1) for row in range(I)]		#time set
		self.TotalSaving = {}		#total saving dic
		self.TotalCost = 0
		
		#status variable
		self.needSchedule = False	#need reschedule or not
		self.newOpp	= False			#is there any new opportunity or not