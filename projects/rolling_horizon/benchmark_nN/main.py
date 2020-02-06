#Author: Zhicheng Zhu
#Email: zhicheng.zhu@ttu.edu, yisha.xiang@ttu.edu
#Benchmark policy of (n_i,N_i)
#paper:
#Cordinated condition-based repair strategies for components 
#multi-component maintenance system with discounts
#EJOR 1997, Wijnmalen et al.
#Last Update: Feb. 2020


import class_info as myClass
import time


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
#overall system parameters
numCom = 1;					#types of components
sysCs = 200;					#system level cs
csLevel = 2;				
#component information
numInd = [4];			#number of individuals for each type
numStates = [7];			#number of states for each type
comCs = [0];			#component level cs

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
	comInfo.cs = comCs[i];
	comInfo.q = [[0.9, 0.1, 0,  0,  0,  0,  0],\
				 [0,   0.8, 0.2,0,  0,  0,  0],\
				 [0,   0,   0.8,0.1,0.1,0,  0],\
				 [0,   0,   0,  0.7,0.2,0.1,0],\
				 [0,   0,   0,  0,  0.6,0.2,0.2],\
				 [0,   0,   0,  0,  0,  0.5,0.5],\
				 [0,   0,   0,  0,  0,  0,  1]];
	#step 1. initialization of upper control limit & rule
	pi1 = 3;				#just a random number for testing
	comInfo.rule = [pi1]*sysInfo.csLevel;
	comInfo.pVec = [[0]*sysInfo.csLevel for j in range(2)];
	comInfo.pVec[0][0] = 1;
	comInfo.pVec[1][0] = 1;
	if sysInfo.csLevel == 2:	#if it is 2, do not over write it.
		comInfo.V = [0]*sysInfo.csLevel;
		comInfo.V[1] = sysInfo.cs;
	
	comInfo.repair_cost();		#initialize repair cost
	comInfo.operational_cost();	#initialize operational cost
	sysInfo.comInfoAll.append(comInfo);
	
# calculate the aggregation info given current rule and pVector

counter = 0;
convergeFlag = False;					#loop it until converge 
while (convergeFlag == False and counter <20):
	counter += 1;
	print("counter", counter);
	sysInfo.print_data();
	sysInfo.comInfoAll[0].print_data1();
	sysInfo.comInfoAll[0].print_data2();
	sysInfo.calAggInfo();		#theta,phiRep,phiK 
	
	#step 2. probability vector 
	sysInfo.calpVec();		
	
	#step 3. policy update. Find the best policy given pVector
	sysInfo.policyIter();		
	
	#step 4.system cost calculation 
	sysInfo.calSysCost();
	
	#step 5. converge condition. Go to step 2 if not satisfied.
	convergeFlag = sysInfo.checkConverge();


print("counter", counter);
sysInfo.print_data();
sysInfo.comInfoAll[0].print_data1();
sysInfo.comInfoAll[0].print_data2();