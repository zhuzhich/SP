#Author: Zhicheng Zhu
#Email: zhicheng.zhu@ttu.edu, yisha.xiang@ttu.edu
#Benchmark policy of (n_i,N_i)
#paper:
#Cordinated condition-based repair strategies for components 
#multi-component maintenance system with discounts
#EJOR 1997, Wijnmalen et al.
#Last Update: Feb. 2020


import class_info as myClass
import numpy as np
import time
import random


#
#
######################
#start from here
######################	
#
###
#sysInfo <- component_info. These include running information
###
###
#######
#Step 0: define sysInfo
#######
##parameter needs to change
wShapeHigh = [4,7];
wShapeLow = [1,3];
wScaleHigh = [5,10];
wScaleLow = [1,5];
sysCsHigh = 100;
sysCsLow = 5;
cCrHigh = [17,27];
cCrLow = [6,16];
#
wShape = wShapeLow;#wShapeHigh, wShapeLow,
wScale = wScaleLow;#wScaleHigh, wScaleLow
sysCs = sysCsHigh;#sysCsHigh, sysCsLow
cCr = cCrLow;#cCrHigh,cCrLow
aaaa = 100;
#overall system parameters
numCom = 8;					#types of components
csLevel = 2;				#Don't change it!
#component information
numInd = [1]*numCom;			#number of individuals for each type
numStates = [0]*numCom;			#number of states for each type. Just initialization
cPr = [1]*numCom;				#pr cost 

#init system information
sysInfo = myClass.system_information();
sysInfo.I = numCom;
sysInfo.cs = sysCs;
sysInfo.csLevel = csLevel;		#only works for 2.

for i in range(sysInfo.I):
	comInfo = myClass.component_information();
	comInfo.index = i;
	comInfo.ind = numInd[i];
	comInfo.state = numStates[i];
	comInfo.cs = sysCs;
	comInfo.costPr = cPr[i];
	#comInfo.q = [[0.9, 0.1, 0,  0,  0,  0,  0],\
	#			 [0,   0.8, 0.2,0,  0,  0,  0],\
	#			 [0,   0,   0.8,0.1,0.1,0,  0],\
	#			 [0,   0,   0,  0.7,0.2,0.1,0],\
	#			 [0,   0,   0,  0,  0.6,0.2,0.2],\
	#			 [0,   0,   0,  0,  0,  0.5,0.5],\
	#			 [0,   0,   0,  0,  0,  0,  1]];
	#step 1. initialization of upper control limit & rule
	pi1 = 3;				#just a random number for testing
	comInfo.rule = [pi1]*sysInfo.csLevel;
	comInfo.pVec = [[0]*sysInfo.csLevel for j in range(2)];
	comInfo.pVec[0][0] = 1;
	comInfo.pVec[1][0] = 1;
	if sysInfo.csLevel == 2:	#if it is 2, do not over write it.
		comInfo.V = [0]*sysInfo.csLevel;
		comInfo.V[1] = sysInfo.cs;
	comInfo.wShapeBound = wShape;
	comInfo.wScaleBound = wScale;
	comInfo.cCrBound = cCr;
	comInfo.init_parameters();	#sample shape, scale, and cr cost parameters.
	comInfo.init_q();			#get # of states, and transition probability.
	comInfo.init_rule();		#get the initial rule
	comInfo.repair_cost();		#initialize repair cost
	comInfo.operational_cost();	#initialize operational cost
	sysInfo.comInfoAll.append(comInfo);
	
# calculate the aggregation info given current rule and pVector

counter = 0;
convergeFlag = False;					#loop it until converge 
while (convergeFlag == False and counter <20):
	counter += 1;
	#print("=======counter==========", counter);
	#sysInfo.print_data();
	#for i in range(sysInfo.I):
		#print("--------i=--------",i);
		#sysInfo.comInfoAll[i].print_data1();
		#sysInfo.comInfoAll[i].print_data2();
	sysInfo.calAggInfo();		#theta,phiRep,phiK 
	
	#step 2. probability vector 
	sysInfo.calpVec();		
	
	#step 3. policy update. Find the best policy given pVector
	sysInfo.policyIter();		
	
	#step 4.system cost calculation 
	sysInfo.calSysCost();
	
	#step 5. converge condition. Go to step 2 if not satisfied.
	convergeFlag = sysInfo.checkConverge();


print("=======counter==========", counter);
sysInfo.print_data();
#for i in range(sysInfo.I):
	#print("--------i=--------",i);
	#sysInfo.comInfoAll[i].print_data1();
	#sysInfo.comInfoAll[i].print_data2();
	


T = 20;
costAll = [];
for rep in range(aaaa,aaaa+5):
	age = [0]*sysInfo.I;
	failProb = [0]*sysInfo.I;
	cost = 0;
	for t in range(T):
		costTmp = 0;
		if t == T-1:
			for i in range(sysInfo.I):
				if age[i] == sysInfo.comInfoAll[i].state-1:
					costTmp += sysInfo.comInfoAll[i].costCr;
			if costTmp > 0:
				costTmp += sysInfo.cs;
			cost += costTmp
		else:
			mx = [0]*sysInfo.I;
			#1.take action
			for i in range(sysInfo.I):
				if age[i] >= sysInfo.comInfoAll[i].rule[0]:
					if age[i] == sysInfo.comInfoAll[i].state-1:
						costTmp += sysInfo.comInfoAll[i].costCr;
					else:
						costTmp += sysInfo.comInfoAll[i].costPr;
					mx[i] = 1;
			if costTmp > 0:
				costTmp += sysInfo.cs;
				for i in range(sysInfo.I):
					if age[i] >= sysInfo.comInfoAll[i].rule[1] and\
						age[i] < sysInfo.comInfoAll[i].rule[0]:
						costTmp += sysInfo.comInfoAll[i].costPr;
						mx[i] = 1;
			cost += costTmp;
			#2.cal the age at t+1:
			for i in range(sysInfo.I):
				if mx[i] == 1:
					age[i] = 0;
				failProb = sysInfo.comInfoAll[i].q[age[i]][-1];
				random.seed(rep+t*100+i);
				randProb = random.uniform(0,1);
				if randProb < failProb:
					#failed
					age[i] = sysInfo.comInfoAll[i].state-1;
				else:
					age[i] += 1;
	costAll.append(cost);
print("costAll", costAll);
print("mean", np.mean(costAll));
print("variance", np.var(costAll));
				





