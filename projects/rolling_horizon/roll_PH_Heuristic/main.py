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
import class_info

def weibull_cdf(w_shape, w_scale, t):
	tmp = float(t)/w_scale;
	tmp1 = tmp ** w_shape;
	return 1 - math.exp(-tmp1);
		
def cond_fail_prob(w_shape, w_scale, age):
	tmp = weibull_cdf(w_shape, w_scale, age+1) -  weibull_cdf(w_shape, w_scale,age);
	tmp1 = 1 -  weibull_cdf(w_shape, w_scale,age);
	prob = float(tmp)/tmp1;	
	return prob;
	
def cost_f(A, sysInfo, w_pr, x_agr, rho):
	cost = 0;
	x = [];		#first stage solution, only shows replacement or not
	for i in range(sysInfo.nComponents):
		x.append(int(A[i][0]>0));
		for j in range(sysInfo.nStages):
			if A[i][j] == 1:
				cost += sysInfo.comInfoAll[i].cPR;
			elif A[i][j] == 2: 
				cost += sysInfo.comInfoAll[i].cCR;
	for j in range(sysInfo.nStages):			
		if sum(A[ii][j] for ii in range(sysInfo.nComponents)) >0:
			cost += sysInfo.cS;
	#if rho == -1:
	#	return cost;
		
	#first-stage variables
	for i in range(len(x)):
		if rho >=0: #rho = -1 means the first iteration
			cost += w_pr[i]*x[i]
	if rho >= 0: #rho = -1 means the first iteration
		x_array = np.array(x)
		tmp = (np.linalg.norm(x_array-x_agr))**2;
		cost += rho/2*tmp;
		
	return cost

def A_to_X(A):
	x = [];
	for i in range(len(A)):
		x.append(int(A[i][0]>0));
	return x
'''
def averge_cost(cost_opt_v):
	sum_cost = float(0)
	length_cost = len(cost_opt_v)
	for i in range(length_cost):
		sum_cost += cost_opt_v[i]
	return sum_cost/length_cost
'''
#heuristic algorithm 
#input: 
#LT: life time of a scenario
#output:
#A_opt:optimal policy
def heurstic_alg(sysInfo, scen, w_pr, x_agr, rho):
	#initialization:
	heuInfo = class_info.heuristic_info(sysInfo);
	
	I = sysInfo.nComponents;
	T = sysInfo.nStages;

	t1_v = [0,1,2,3]
	t2_v = range(10)
	A_opt = [];
	cost_opt = float("inf");
	for t1 in t1_v:
		for t2 in t2_v:
			#initialization
			tmpFlag = False
			tmpFlag1 = [0]*sysInfo.nComponents;
			if t1 >= 2:
				tmpFlag = True
				t1 += -2;
			heuInfo.__init__(sysInfo);
			for idxI in range(sysInfo.nComponents):
				heuInfo.currentWin[idxI][0] = \
							max(sysInfo.comInfoAll[idxI].LT[scen][0] - t1, 0);
				heuInfo.currentWin[idxI][1] = idxI;
				#if rho > 0 and w_pr[idxI] < -1*sysInfo.cS and tmpFlag == True:
				#	heuInfo.currentWin[idxI][0] = 0;
				if rho > 0 and w_pr[idxI] < -1*sysInfo.comInfoAll[idxI].cPR and tmpFlag == True:
					r = 0;
					t = 0;
					penaltyCost = 0;
					shiftT = max(sysInfo.comInfoAll[idxI].LT[scen][0] - t1, 1);
					while (t < T+shiftT):
						if t >= T:
							penaltyCost += sysInfo.comInfoAll[idxI].cPR;
							penaltyCost += sysInfo.cS;
						deltaT = max(sysInfo.comInfoAll[idxI].LT[scen][r] - t1, 1);
						t += deltaT;
						r += 1;
					if w_pr[idxI] + penaltyCost < 0:
						heuInfo.currentWin[idxI][0] = 0;
						tmpFlag1[idxI] = 1;
			heuInfo.sort_currentWin();
			heuInfo.shadowWin = copy.deepcopy(heuInfo.currentWin);
		#main loop
			counter = 0;
			while heuInfo.currentWin[0][0] <= T-1:
				counter += 1;
				cost_pole = float("inf");
				#in current window, replace the first one first
				repTime = heuInfo.currentWin[0][0];
				#if t2 == 1:
				#	aa = 1;
				heuInfo.update_i(sysInfo, scen,0, t1,repTime);
				heuInfoProto = copy.deepcopy(heuInfo);
				for pole_i in range(I-1):			#the smallest opportunity.
					heuInfoTenta = copy.deepcopy(heuInfoProto);
					#init cost
					j = pole_i;
					flag = True;
					groupCom = [0];
					penaltyCost = 0;
					for i in range(j + 1, I):
						if heuInfoTenta.currentWin[i][0] -\
							heuInfoTenta.currentWin[j][0] <=t2\
							and heuInfoTenta.lastIndTime[heuInfoTenta.currentWin[i][1]]\
								 < heuInfoTenta.currentWin[j][0]\
							and heuInfoTenta.currentWin[j][0] <= T - 1\
							and heuInfoTenta.currentWin[i][0] <= T - 1:
							#add penalty cost when an individual is shifted
							#step 1: find the original replacement time of the first individual replaced outside the planning horizon
							#        if the individual is not shifted
							#step 2: see the replacement time after shifting within the planning horizon or not
							r = heuInfoTenta.workingInd[heuInfoTenta.currentWin[i][1]];
							t = heuInfoTenta.currentWin[i][0];
							shiftT = heuInfoTenta.currentWin[i][0] - heuInfoTenta.currentWin[j][0];
							penaltyCostTmp = 0
							while (t < T+shiftT):
								if t >= T:
									if penaltyCostTmp > 0:
										penaltyCostTmp += sysInfo.cS;
									penaltyCostTmp += sysInfo.comInfoAll[heuInfoTenta.currentWin[i][1]].cPR;
								r += 1;
								deltaT = max(sysInfo.comInfoAll[heuInfoTenta.currentWin[i][1]].LT[scen][r] - t1, 1);
								t += deltaT;
							penaltyCost += penaltyCostTmp;
							repTime = heuInfoTenta.currentWin[j][0];
							heuInfoTenta.update_i(sysInfo, scen, i, t1, repTime);
							groupCom.append(i);
							# replace the pole first
							if flag == True and j != 0:
								flag = False;
								heuInfoTenta.update_i(sysInfo, scen, j, t1,repTime);
								groupCom.append(j);
						else:
							flag == True;
							j = i;				                  
					for i in range(I):
						if i not in groupCom and heuInfoTenta.currentWin[i][0] < T:	
							heuInfoTenta.update_i(sysInfo, scen, i, t1, heuInfoTenta.currentWin[i][0]);
					if 0:#tmpFlag == False:
						#add the handling of last individuals
						#fail all last individual see if they can be removed from the planning horizon
						for i in range(I):
							#if 	tmpFlag == True:
								#continue;
							lastR = heuInfoTenta.workingInd[i] - 1;
							if lastR < 0:
								continue;
							if lastR == 0 and tmpFlag1[i] == True:
								continue;
							tmpCount = 0;
							for t in range(T):
								lastT = T - t - 1;
								if heuInfoTenta.sol[i][lastT] > 0:
									if tmpCount == 0:
										tmpCount += 1;
									elif tmpCount == 1:
										break;
							lastLT = sysInfo.comInfoAll[i].LT[scen][lastR];
							if lastT + lastLT >= T:
								heuInfoTenta.sol[i][heuInfoTenta.lastIndTime[i]] = 0;
								heuInfoTenta.lastIndTime[i] = lastT;
								heuInfoTenta.workingInd[i] += -1;						
					
					heuInfoTenta.cost = cost_f(heuInfoTenta.sol,sysInfo, w_pr, x_agr, rho);
					if heuInfoTenta.cost + penaltyCost< cost_pole:
						cost_pole = heuInfoTenta.cost + penaltyCost;
						heuInfo = copy.deepcopy(heuInfoTenta);
						heuInfo.currentWin = \
						copy.deepcopy(heuInfo.shadowWin);
						heuInfo.sort_currentWin();
						heuInfo.shadowWin = \
						copy.deepcopy(heuInfo.currentWin);
				
			if  heuInfo.cost < cost_opt:
				opt_t1 = t1;
				opt_t2 = t2;
				cost_opt = heuInfo.cost;
				A_opt = heuInfo.sol;
	#print ("opt_t1,opt_t2");
	#print (opt_t1, opt_t2);
	#print (A_opt);
	return A_opt



	
	
def PH_alg(sysInfo):

	
	rho = 1.0;
	max_iter = 15;
	w_pr = np.zeros((sysInfo.nScenarios, sysInfo.nComponents));#array type	
	x_agr = np.zeros(sysInfo.nComponents);  #array type
	eps = 1.0e-1
	for iter_ in range(max_iter):
		x = []
		cost_opt_v = []
		for idxW1 in range(sysInfo.nScenarios):	
			#print ("iter_=%d,scen=%d" %(iter_,idxW1+1))
			if iter_ == 0:
				A = heurstic_alg(sysInfo,idxW1,0,0,-1) #-1 means the first iteration
			else:
				A = heurstic_alg(sysInfo, idxW1,w_pr[idxW1],x_agr,rho)
			#print ("A = ", A);
			x.append(A_to_X(A))
			#print ("x="),
			#print (A_to_X(A))
			#cost_opt = cost_f(A, 0, 0, -1)
			cost_opt = cost_f(A,sysInfo, 0, 0, -1);
			cost_opt_v.append(cost_opt)
		cost_avg = np.average(cost_opt_v)
		x_array = np.array(x)  #convert to array
		#get a new aggregated x
		x_agr = (1.0/sysInfo.nScenarios)*\
				sum(x_array[i] for i in range(sysInfo.nScenarios));
		
		#for i in range(sysInfo.nScenarios):
			#if x_agr[i] > 0.9:
				#x_agr[i] = 1.0;
			#if x_agr[i] < 0.1:
				#x_agr[i] = 0;
		#print ("iter_=", iter_)
		#print ("x_array = ", x_array)
		#print ("x_agr = ", x_agr)
		#print ("w_pr = ", w_pr)
		#print ("average cost= %f" %cost_avg)
		if iter_ > 0:
			tmp = (1.0/sysInfo.nScenarios)\
					*sum(np.linalg.norm(x_array[i]-x_agr) \
					for i in range(sysInfo.nScenarios));
			if tmp <= eps or iter_ == max_iter-1:
				print ("converge at iter_ = %d" %iter_)
				#print ("cost_avg = %f" %cost_avg)
				return list(x_agr)
				
				
		#get price w_pr
		for idxW1 in range(sysInfo.nScenarios):
			if iter_ == 0:
				w_pr[idxW1] = rho*(x_array[idxW1] - x_agr);
			else:
				w_pr[idxW1] = w_pr[idxW1] + rho*(x_array[idxW1] - x_agr);
	


def main(wShapeBound, wScaleBound, setUpCost, crBound, ranSeed ):

	nComponents = 8		#number of components
	nStages = 20			# t = 1, 2, ..., T
	T = nStages;			#remaining time horizon
	R = T + 20;				#number of individuals
	nScenarios = 250;		#number of scenarios
	cS = setUpCost			#setup cost	
	intvl = 1;
	
	sysInfo = class_info.system_info(nComponents, T, intvl, cS, nScenarios, ranSeed);
	
	kesi = [0]*nComponents;
	age = [0]*nComponents;
	cPR = [0]*nComponents;
	cCR = [0]*nComponents;
	w_shape = [0]*nComponents;
	w_scale = [0]*nComponents;

	for i in range(nComponents):
		#kesi
		if i==0:
			kesi[i] = 0;
		#age
		age[i] = 0;
		#cPR
		cPR[i] = 1;
		#cCR
		random.seed(i*30);
		temp = random.uniform(crBound[0], crBound[1]);
		cCR[i] = round(temp,1);			
		#shape
		random.seed(i*20)   
		temp = random.uniform(wShapeBound[0],wShapeBound[1])
		w_shape[i] = round(temp,1);			
		#scale
		random.seed(i*10)   
		temp = random.uniform(wScaleBound[0], wScaleBound[1])#(4,11)
		w_scale[i] = round(temp,1);		
		comInfo = class_info.component_info(i, w_shape[i], w_scale[i], age[i], kesi[i], \
											intvl, cCR[i], cPR[i], cS, nScenarios,R);
											# LT is initialized to be empty.
		sysInfo.add_com(comInfo);
		sysInfo.comInfoAll[i].create_LifeTime(sysInfo.ranSeed);
		#print (sysInfo.comInfoAll[i].LT);
	retCost = 0;

	for t in range(nStages):
		print ("t=",t);
		#sTime = time.clock();
		if t == nStages - 1:
			tmp = 0;
			for i in range(sysInfo.nComponents):
				if sysInfo.comInfoAll[i].initFail == 1:
					tmp += sysInfo.comInfoAll[i].cCR
					print ("fail",i);
				
			if tmp > 0:
				tmp += sysInfo.cS;	
			retCost += tmp;
		else:
			## solve current stage problem
			#flag = False;
			#if t == 6:
				#flag = True;
				#print("LT[0]", sysInfo.comInfoAll[0].LT);
				#print("LT[1]", sysInfo.comInfoAll[1].LT);
			x = PH_alg(sysInfo);	
			#print ("x=",x);
			res = copy.deepcopy(x);
			for i in range(sysInfo.nComponents):
				if sysInfo.comInfoAll[i].initFail == 1:
					res[i] = 2;
			print ("res=",res);
			x = [round(x[i]) for i in range(len(x))];	
			
			#calculate cost at current stage
			tmp = 0;
			for i in range(sysInfo.nComponents):
				tmp += sysInfo.comInfoAll[i].cPR*x[i];
				tmp += (sysInfo.comInfoAll[i].cCR - sysInfo.comInfoAll[i].cPR)*\
					sysInfo.comInfoAll[i].initFail;
				#if (x[i] == 1):
					#print (sysInfo.comInfoAll[i].initFail);
			if tmp > 0:
				tmp += sysInfo.cS;
			retCost += tmp;
			
			## prepare next stage:
			#generate the failure
			ageAfterMx = [sysInfo.comInfoAll[i].initAge*(1-x[i]) \
							for i in range(sysInfo.nComponents)];
			failProb = [sysInfo.comInfoAll[i].cond_fail_prob(ageAfterMx[i],ageAfterMx[i]+1) \
						for i in range(sysInfo.nComponents)];
			
			randProb = [];
			for i in range(sysInfo.nComponents):
				random.seed(ranSeed + t*100 + i);
				randProb.append(random.uniform(0,1));
			
			failState =  [int(randProb[i]<failProb[i]) for i in range(sysInfo.nComponents)];
			ageAfterMx = [ageAfterMx[i]+1 for i in range(sysInfo.nComponents)];
			
			#sysInfo.nStages += -1;
			sysInfo.nStages = max(sysInfo.nStages, 10)
			#R = sysInfo.nStages + 1;

			for i in range(sysInfo.nComponents):
				sysInfo.comInfoAll[i].initAge = ageAfterMx[i];
				sysInfo.comInfoAll[i].initFail = failState[i];
				#sysInfo.comInfoAll[i].nIndividuals = R;
				sysInfo.comInfoAll[i].create_LifeTime(sysInfo.ranSeed + t*100);
		#print("failProb", failProb);
		#print("randProb", randProb);
		#print("failState,", failState);
		#print("LT[0]", sysInfo.comInfoAll[0].LT);
		#print("LT[1]", sysInfo.comInfoAll[1].LT);
		#eTime = time.clock();
		#print ("time=", eTime-sTime);
							
	return retCost
				
 
################
## start
################
wShapeH = [4,7];
wShapeL = [1,3];
wShape = wShapeH;#wShapeH, wShapeL
wScaleH = [5, 10]
wScaleL = [1, 5]
wScaleNormal = [1,8]
wScaleVector = [wScaleL]#[wScaleH, wScaleL]
dVector = [5]#[100, 5]
cCrH = [17, 27]
cCrL = [6, 16]
cCrVector = [cCrL]#[cCrH, cCrL]

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
			aaa = 100;
			for resLifeSeed in range(aaa-2,aaa):
				#counter1 += 1
				ranSeed = resLifeSeed+counter1;
				#
				#ranSeed = resLifeSeed + counter
				tmp = main(wShape, wScale, d, cCr, ranSeed)
				print ("cost[%d] = %f" %(counter1, tmp))
				cost.append(tmp)
			end_time = time.clock()
			time_ = end_time-start_time			
			print ("time             = %d" %time_)
			print ("cost")
			print (cost)
			print ("average cost")
			print (float(sum(cost))/len(cost))
			print ("variance,", np.var(cost));
			print ("=============================================")
				

		
		
