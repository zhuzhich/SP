#Author: Zhicheng Zhu
#Email: zzhu3@lamar.edu
#A heuristic algorithm for multi-unit system maintenance problem.
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
def create_LifeTime(LT_in):    
	global I
	global R
	global w
	global w_shape
	global w_scale
	global pKesi
	global s

	#local variable
	I_Lo = I
	w_Lo = w
	R_Lo = R
	w_shape_Lo = w_shape
	w_scale_Lo = w_scale
	pKesi_Lo = pKesi
	s_Lo = s
	LT_tmp = []
	for idx_w in range(0,w_Lo):
		tmp = [[0]*R_Lo for row in range(I_Lo)]
		for idx1 in range(0,I_Lo):
			i = idx1
			for idx2 in range(0,R_Lo):
				r = idx2
				random.seed(i+r+idx_w) ###control the seed
				if r == 0:
					if pKesi_Lo[idx1] == 1:
						LT = 0
					else:    
						ran_num = random.uniform(0,1)
						ran_num_log = -math.log(ran_num)
						s_inv = 1.0/w_shape_Lo[idx1]	
						LT1 = round((ran_num_log**s_inv)*w_scale_Lo[idx1]) - s_Lo					
						LT = int(max(1,LT1))
				else:
					ran_num = random.uniform(0,1)
					ran_num_log = -math.log(ran_num)
					s_inv = 1.0/w_shape_Lo[idx1]
					LT1 = round((ran_num_log**s_inv)*w_scale_Lo[idx1])
					LT = int(max(1,LT1))
				##assign variable LT
				tmp[idx1][idx2] = LT
		#LT_tmp <- tmp
		LT_tmp.append(tmp)
	LT_in += LT_tmp
"""
	for idx_w in range(0,w_Lo):
		print ("scen=%d" %idx_w)
		for i in range(I_Lo):
			print ("["),
			for r in range(R_Lo):
				print (LT_in[idx_w][i][r]),
			print ("]")
"""

#cost of PR
def pCPR_init(pCPR):
	global I
	for i in range(0, I):
		pCPR.append(1)
#cost of CR: 
def pCCR_init(pCCR):
	global I
	for i in range(0, I):
		random.seed(i*30) ###control the seed     
		temp = random.uniform(6,16)
		pCCR.append(round(temp,1))
#kesi. Only the first component fails at current time
def pKesi_init(pKesi):
	global I
	for i in range(0, I):
		if i+1 == 1:
			pKesi.append(1)
		else:
			pKesi.append(0)		
def w_shape_init(w_shape):
	global I
	for i in range(0, I):
		random.seed(i*20) ###control the seed     
		temp = random.uniform(4,7)
		w_shape.append(round(temp,1))
#weibull scale parameter: w_scale[i]=i, i=1...I
#u(1,8)
def w_scale_init(w_scale):
	global I
	for i in range(0, I):
		random.seed(i*10) ###control the seed     
		temp = random.uniform(1,8)#(4,11)(1,8)
		w_scale.append(round(temp,1))
		
def cost_f(A,LT, w_pr, x_agr, rho):
	global I
	global T
	global pCPR
	global pCCR
	global d
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
	t2_v = range(10)
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
				B[idxI][0] = int(max(LT[idxI][0] - t1, 0))
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
					if (pole_i in node_dict) == False:
						node_data = []
						for i in range(j+1, I):	#component that trying to move.
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
								if (j in node_dict) == False:
									node_dict[j] = Node(node_data,-1)
								node_dict[j].pnext = -1				
					D_tenta[0] = 1
					t_tenta[B[0][1]] = B[0][0]
					#if t_tenta[B[0][1]] <= T:
					A_tenta[B[0][1]][B[0][0]] = 1
					k = first_key
					while ((k in node_dict) == True):
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
						B[i][0] = int(max(LT[com][next_ind]-t1, 1) + t[com])
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
	
	rho = 20.0
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
			#print ("A=", A);
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
			if tmp <= eps:
				print ("converge at iter_ = %d" %iter_)
				print ("x = ", x_agr)
				print ("average cost= %f" %cost_avg)
				return
				
		#get price w_pr
		for idxW1 in range(w):
			if iter_ == 0:
				w_pr[iter_][idxW1] = rho*(x_array[idxW1] - x_agr)
			else:
				w_pr[iter_][idxW1] = w_pr[iter_-1][idxW1] + rho*(x_array[idxW1] - x_agr)
	
		
global I
global T
global w
global T_ex
global s
global R
global d
global pCPR
global pCCR
global pKesi
global w_shape
global w_scale
global x
global directory

comp_list = [2]#[4,6,8]
time_list = [4]# t = 0, 1, 2, ..., T
scen_list = [8]
counter = 0
for idxI in comp_list:
	for idxT in time_list:
		for idxW in scen_list:
			#independent variables
			I = idxI			#number of components
			T = idxT			#time horizon
			w = idxW			#number of scenarios
			#dependent variables
			#T_ex = T + 40	#extended time horizon
			s = 2			#starting time
			R = T + 2		#number of individuals
			d = 5		#setup cost
			counter += 1
			print ("================(%d,%d,%d,%d)=============" %(counter,idxI,idxT,idxW))
			#print ("d=%d" %d)
			pCPR = []
			pCPR_init(pCPR)

			pCCR = []
			pCCR_init(pCCR)


			pKesi = []
			pKesi_init(pKesi)

			w_shape = []
			w_shape_init(w_shape)	

			w_scale = []
			w_scale_init(w_scale)
			
			LT = []
			create_LifeTime(LT)
			
			#print ("LT", LT);
			start_time = time.clock()
			
			PH_alg(LT)
			
			end_time = time.clock()
			time_ = end_time-start_time
			print ("time=%d" %time_)
 


		
		
