#Author: Zhicheng Zhu
#Email: zzhu3@lamar.edu
#Rolling horizon version of Algorithm 4: PH_Heuristic.
#Last Update: 03/21/2018

import os
import math
import random
import numpy as np
import time
import copy



#build up a linkedlist for accelerate the alg.
class Node:
	def __init__(self, data, pnext):
		#data:
		#	is replacement decisions [(a1,b1),(a2,b2)...]
		#	(a1,b1) means index (in B) a1 replaced at time b1
		#pnext:
		#	next j in B.
		self.data = data
		self.pnext = pnext


#create life time using the same seed
def create_LifeTime(LT_in, timeElapsed, residualLife):    
	global I
	global R
	global w
	global w_shape
	global w_scale
	global s

	#local variable
	I_Lo = I
	w_Lo = w
	R_Lo = R
	w_shape_Lo = w_shape
	w_scale_Lo = w_scale
	s_Lo = s
	
	if len(LT_in) > 0:
		LT_in_arr = np.array(LT_in)
		LT_arr = LT_in_arr[:,:,range(R_Lo)]
		for idx_w in range(0,w_Lo):
			for idx1 in range(0,I_Lo):
				i = idx1 + 1
				if residualLife[idx1] == 0:
					LT_arr[idx_w,idx1,0] = 0	#force to replace the first individual	
				else:
					random.seed(i+1+idx_w+1) ###control the seed  				
					ran_num = random.uniform(0,1)
					part1 = math.log(ran_num)
					s_inv = 1.0/w_shape_Lo[idx1]
					surv_time = timeElapsed[idx1]
					part2 = (surv_time/w_scale_Lo[idx1])**w_shape_Lo[idx1]
					part3 = part2 - part1
					tmp = round((part3**s_inv)*w_scale_Lo[idx1]) - surv_time					
					LT_arr[idx_w,idx1,0] = int(max(0,tmp))	
		LT = list(LT_arr)
		return LT
	
	######for the first time######################
	LT = []
	for idx_w in range(0,w_Lo):
		tmp = [[0]*R_Lo for row in range(I_Lo)]
		for idx1 in range(0,I_Lo):
			i = idx1 + 1
			for idx2 in range(0,R_Lo):
				r = idx2 + 1
				random.seed(i+r+idx_w+1) ###control the seed  
				ran_num = random.uniform(0,1)
				ran_num_log = -math.log(ran_num)
				s_inv = 1.0/w_shape_Lo[idx1]	
				LT1 = round((ran_num_log**s_inv)*w_scale_Lo[idx1])			
				##assign variable LT
				tmp[idx1][idx2] = int(max(1,LT1))
		#LT_tmp <- tmp
		LT.append(tmp)
	return LT

#cost of PR
def pCPR_init(pCPR):
	global I
	for i in range(0, I):
		pCPR.append(1)
#cost of CR: 
def pCCR_init(pCCR, crBound):
	global I
	for i in range(0, I):
		random.seed((i+1)*30) ###control the seed     
		temp = random.uniform(crBound[0],crBound[1])
		pCCR.append(round(temp,1))
	
def w_shape_init(w_shape):
	global I
	for i in range(0, I):
		random.seed((i+1)*10) ###control the seed     
		temp = random.uniform(4,7)
		w_shape.append(round(temp,1))
#weibull scale parameter: w_scale[i]=i, i=1...I
#u(1,8)

def w_scale_init(w_scale, wScaleBound):
	global I
	for i in range(0, I):
		random.seed((i+1)*20) ###control the seed     
		temp = random.uniform(wScaleBound[0],wScaleBound[1])
		w_scale.append(round(temp,1))
		
def cost_f(A,LT, w_pr, x_agr, rho):
	global I
	global T
	global pCPR
	global pCCR
	global d
	global s
	cost = 0
	x = []
	for i in range(I):
		ind = 0#individual index
		last_time = -1
		x.append(A[i][0])
		for j in range(T+1):
			if A[i][j] == 1:
				#if i == 15:
				#	print (j)
				if last_time == -1:
					elapse_time = j
				else:
					elapse_time = j - last_time
				last_time = j
				if elapse_time < LT[i][ind]: #smaller than life time, PR
					cost += pCPR[i]
				elif elapse_time == LT[i][ind]: #equals to life time, CR
					#print ("CR!!! elapse_time(%d,%d)=%d=LT(%d)" %(i,ind,elapse_time,LT[i][ind]))
					cost += pCCR[i]
				else:
					print ("error, elapse_time(%d,%d)=%d>LT(%d)" %(i,ind,elapse_time,LT[i][ind]))
				ind += 1

			#cost += sum(A[i][j] for j in range(T+1)) * pCPR[i]
		
	for j in range(T+1):			
		if sum(A[ii][j] for ii in range(I)) >0:
			cost += d
		
	#first-stage variables
	for i in range(len(x)):
		#if i == 1:
		#	cost += pCCR[i]*x[i] 
		#else:
		#	cost += pCPR[i]*x[i]
		if rho >=0: #rho = -1 means the first iteration
			cost += w_pr[i]*x[i]
	if rho >= 0: #rho = -1 means the first iteration
		x_array = np.array(x)
		tmp = (np.linalg.norm(x_array-x_agr))**2
		cost += rho/2*tmp
		
	return cost

def sort_B(B):
	B.sort(key=lambda x:x[0])
	return B

	
def A_to_X(A):
	x = []
	for i in range(len(A)):
		x.append(A[i][0])
	return x
	
def averge_cost(cost_opt_v):
	sum_cost = float(0)
	length_cost = len(cost_opt_v)
	for i in range(length_cost):
		sum_cost += cost_opt_v[i]
	return sum_cost/length_cost
#heuristic algorithm 
#input: 
#LT: life time of a scenario
#output:
#A_opt:optimal policy
def heurstic_alg(LT, w_pr, x_agr, rho):
	#initialization:
	global I
	global T
	t1_v = [0,1]
	t2_v = range(0,10)
	A_opt = []
	cost_opt = 10000000;
	for t1 in t1_v:
		for t2 in t2_v:
		#init for each iteration
			#print ("start for t1=%d,t2=%d" %(t1,t2))
			B = [[0]*2 for row in range(I)]
			A = [[0]*(T+1) for row in range(I)]
			t = [-1]*I
			D = [0]*I
			for idxI in range(I):
				B[idxI][0] = max(LT[idxI][0] - t1, 0)
				B[idxI][1] = idxI
			B = sort_B(B)	
		#main loop
			counter = 0
			while B[0][0] <= T:
				counter += 1
				D = [0]*I		# = 1 if component has been moved before.
				#shift procedure
				D_opt = copy.deepcopy(D)	
				t_opt = copy.deepcopy(t)	
				A_opt_L = copy.deepcopy(A)
				node_dict = {}
				cost_opt_L = 10000000	
				pole_i_opt = -1				
				for pole_i in range(I-1):			#the smallest opportunity.
					D_tenta = copy.deepcopy(D)
					t_tenta = copy.deepcopy(t)
					A_tenta = copy.deepcopy(A)
					cost_tenta = 0
					#init cost
					j = pole_i
					first_key = pole_i
					if pole_i in node_dict == False:
						node_data = []
						for i in range(pole_i+1, I):	#component that trying to move.
							#com_i = B[i][1]
							#com_j = B[j][1]
							if B[i][0] - B[j][0] <= t2 and \
								t_tenta[B[i][1]] < B[j][0] and \
								B[j][0] <= T:								
								node_data.append([i,B[j][0]])
							else:
								node_dict[j] = Node(node_data,i)
								if i in node_dict == True:
									break
								else:
									if i < I-1:	#i cannot be the last one
										j = i
										node_data = []
							if i == I - 1: #last element
								if j in node_dict == False:
									node_dict[j] = Node(node_data,-1)
								node_dict[j].pnext = -1				
					D_tenta[0] = 1
					t_tenta[B[0][1]] = B[0][0]
					#if t_tenta[B[0][1]] <= T:
					A_tenta[B[0][1]][B[0][0]] = 1
					k = first_key
					while (k in node_dict == True):
						if B[k][0] > T:
							break
						#replace k first
						if len(node_dict[k].data) > 0:
							D_tenta[k] = 1
							t_tenta[B[k][1]] = B[k][0]
							A_tenta[B[k][1]][B[k][0]] = 1
						#for other non-pole component
						data = node_dict[k].data
						k = node_dict[k].pnext
						for m in range(len(data)):
							idx = data[m][0]
							r_time = data[m][1]							
							com_m = B[idx][1]
							#assign value
							D_tenta[idx] = 1
							t_tenta[com_m] = r_time
							A_tenta[com_m][r_time] = 1							
					cost_tenta = cost_f(A_tenta,LT, w_pr, x_agr, rho)
					#update optimal
					if cost_tenta < cost_opt_L:
						pole_i_opt = pole_i
						D_opt = copy.deepcopy(D_tenta)
						t_opt = copy.deepcopy(t_tenta)
						A_opt_L = copy.deepcopy(A_tenta)
						cost_opt_L = cost_tenta
				#print (counter,pole_i_opt,D_opt,t_opt,cost_opt_L)
				#update from optimal
				D = copy.deepcopy(D_opt)
				t = copy.deepcopy(t_opt)
				A = copy.deepcopy(A_opt_L)
				#update some parameters
				for i in range(I):
					if D[i] == 1:
						com = B[i][1]
						next_ind = sum(A[com][ii] for ii in range(T+1))
						B[i][0] = max(LT[com][next_ind]-t1, 1) + t[com]
				B = sort_B(B)
			cost = cost_f(A, LT, w_pr, x_agr, rho)
			if cost < cost_opt:
				A_opt = copy.deepcopy(A)
				cost_opt = cost
				t1_opt = t1
				t2_opt = t2
	#print ("optimal (t1,t2,cost)=(%d,%d,%d)" %(t1_opt, t2_opt,cost_opt))
	return A_opt


	
	
def PH_alg(LT):
	global I
	global R
	global w
	
	rho = 50.0
	max_iter = 10
	w_pr = np.zeros((max_iter, w, I))#array type	
	x_agr = np.zeros(I)  #array type
	eps = 1.0e-4
	for iter_ in range(max_iter):
		x = []
		cost_opt_v = []
		for idxW1 in range(w):	
			#print ("iter_=%d,scen=%d" %(iter_,idxW1+1))
			if iter_ == 0:
				A = heurstic_alg(LT[idxW1],0,0,-1) #-1 means the first iteration
			else:
				A = heurstic_alg(LT[idxW1],w_pr[iter_-1][idxW1],x_agr,rho)
			x.append(A_to_X(A))
			#print ("x="),
			#print (A_to_X(A))
			cost_opt = cost_f(A,LT[idxW1], 0, 0, -1)
			cost_opt_v.append(cost_opt)
		cost_avg = averge_cost(cost_opt_v)
		x_array = np.array(x)  #convert to array
		#get a new aggregated x
		x_agr = 1.0/w*sum(x_array[i] for i in range(w))
		#print ("x_array = ")
		#print (x_array)
		#print ("x_agr = ")
		#print (x_agr)
		#print ("average cost= %f" %cost_avg)
		if iter_ > 0:
			tmp = 1.0/w*sum(np.linalg.norm(x_array[i]-x_agr) for i in range(w))
			#print ("converge distance = %f" %tmp)
			if tmp <= eps or iter_ == max_iter-1:
				print ("converge at iter_ = %d" %iter_)
				return list(x_agr)
				
		#get price w_pr
		for idxW1 in range(w):
			if iter_ == 0:
				w_pr[iter_][idxW1] = rho*(x_array[idxW1] - x_agr)
			else:
				w_pr[iter_][idxW1] = w_pr[iter_-1][idxW1] + rho*(x_array[idxW1] - x_agr)
	


def main(wScaleBound, setUpCost, crBound, ranSeed ):
	
	global I
	global T
	global w
	global s
	global R
	global d
	global pCPR
	global pCCR
	global w_shape
	global w_scale
	

	#independent variables
	I = 10			#number of components
	T = 20			#a fixed endClock
	TRoll = 20
	R = T + 2		#number of individuals
	w = 10			#number of scenarios
	#dependent variables
	d = setUpCost		#setup cost
	#print ("d=%d" %d)
	pCPR = []
	pCPR_init(pCPR)

	pCCR = []
	pCCR_init(pCCR, crBound)

	w_shape = []
	w_shape_init(w_shape)	

	w_scale = []
	w_scale_init(w_scale, wScaleBound)

	LT = []
	
	timeElapsed = [1]*I
	residualLife = [10]*I #just greater than 0
	resLifeV = []
	resLifeV.append(residualLife)
	s = 0
	
	LT = create_LifeTime(LT, timeElapsed, residualLife)					
	
	result = [[0]*(TRoll+1) for row in range(I)]
	for s in range(1, TRoll+1):
		print ("current time = %d"  %s)
		T = TRoll - s	#rtime length
		R = T + 2		#number of individuals

		#generate residual life
		for i in range(I):
			wShape = w_shape[i]
			wScale = w_scale[i]
			ranSeed1 = s+i + (ranSeed + 1)**2
			random.seed(ranSeed1)
			ran_num = random.uniform(0,1)
			part1 = math.log(ran_num)
			s_inv = 1.0/wShape
			surv_time = timeElapsed[i]
			part2 = (surv_time/wScale)**wShape
			part3 = part2 - part1
			res = round((part3**s_inv)*wScale) - surv_time					
			residualLife[i] = int(max(0,res))
			
		print ("residual life    = " ),
		print (residualLife)
		resLifeV.append(residualLife)			#update vector for check
		LT = create_LifeTime(LT, timeElapsed, residualLife)					
		x = PH_alg(LT)
		print ("x                = " ),
		print (x)
		#record result
		for i in range(I):
			result[i][s] = int(round(x[i]))
			if residualLife[i] == 0:
				result[i][s] = 2	#2 is CR
		#update next replacement time & LT
		for i in range(I):
			if result[i][s] > 0:
				timeElapsed[i] = 0
				LT_arr = np.array(LT)
				for idxW in range(w):
					tmp = list(LT_arr[idxW,i,:])
					tmp.append(0)
					tmp.pop(0)
					LT_arr[idxW,i,:] = np.array(tmp)
			else:
				timeElapsed[i] += 1
	optCost = 0	
	for j in range(TRoll+1):
		if sum(result[u][j] for u in range(I)) > 0:
			optCost += d
		for i in range(I):
			if result[i][j] == 1:
				optCost += pCPR[i]
			if result[i][j] == 2:
				optCost += pCCR[i] 
	print ("groupResult = " )		
	for row in range(I):
		print (result[row][:])
	print ("cost             = %d" %optCost)
	#tmp = np.array(resLifeV)
	#for i in range(I):
	#	tmp_i = tmp[:,i]
	#	print ("======for com[%d] residual life=======" %i)
	#	print (tmp_i)
	return optCost
				
 
################
## start
################
global I
global T
global w
global s
global R
global d
global pCPR
global pCCR
global w_shape
global w_scale


wScaleH = [9, 20]
wScaleL = [1, 8]
wScaleVector = [wScaleL]#[wScaleH, wScaleL]
dVector = [100]#[100, 5]
cCrH = [17, 27]
cCrL = [6, 16]
cCrVector = [cCrH]#[cCrH, cCrL]

counter = 0
print ("=============Heuristic=====================")
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
				

		
		
