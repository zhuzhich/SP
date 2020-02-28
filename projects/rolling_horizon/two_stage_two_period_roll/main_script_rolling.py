#Author: Zhicheng Zhu
#Email: zhicheng.zhu@ttu.edu, yisha.xiang@ttu.edu
#copyright @ 2019: Zhicheng Zhu. All right reserved.

#Info:
#script to control other main file for rolling horizion
#
#Last update: 03/09/2019
#!/usr/bin/python


import main_dynamic_solver

import class_info
import copy
import random
import time

import numpy as np
		
##########################
#start from here
##########################
# some fixed parameters
nComponents = 2;		#fix component numbers	
cS = 5;				#fix setup cost
intvl = 1;		#don't change this
nStages = 10;			

sysInfo = class_info.system_info(nComponents, nStages, intvl, cS);


kesi = [0]*nComponents;
age = [0]*nComponents;
cPR = [0]*nComponents;
cCR = [0]*nComponents;
w_shape = [0]*nComponents;
w_scale = [0]*nComponents;

for i in range(nComponents):
	#kesi
	if i==0:
		kesi[i] = 1;
	#age
	age[i] = 2;
	#cPR
	cPR[i] = 1;
	#cCR
	random.seed(i*30);
	temp = random.uniform(6,16);
	cCR[i] = round(temp,1);			
	#shape
	random.seed(i*20)   
	temp = random.uniform(4,7)
	w_shape[i] = round(temp,1);			
	#scale
	random.seed(i*10)   
	temp = random.uniform(1,8)#(4,11)
	w_scale[i] = round(temp,1);		
	comInfo = class_info.component_info(i, w_shape[i], w_scale[i], age[i], kesi[i], intvl, cCR[i], cPR[i], cS);
	sysInfo.add_com(comInfo);
	
rep = 10;
overAllCost = [];
sysInfo.nStages = 2;	#two-stage two-period rolling horizon
for repIdx in range(rep):	
	#initilization at the beginning of each replication
	age = [2]*sysInfo.nComponents;
	kesi = [1,0,0,0];
	for i in  range(sysInfo.nComponents):
		sysInfo.comInfoAll[i].initAge = age[i];
		sysInfo.comInfoAll[i].initFail = kesi[i];
	totalCost = 0;
	for t in range(nStages):
		if t == nStages - 1:
			tmp = 0;
			for i in range(sysInfo.nComponents):
				tmp += sysInfo.comInfoAll[i].cCR * sysInfo.comInfoAll[i].initFail;
			if tmp > 0:
				tmp += sysInfo.cS;
			totalCost += tmp;
		else:
			main_dynamic_solver.main(sysInfo);
			#print(t);
			
			#print(kesi);
			#print(sysInfo.solX);
			tmp = 0;
			for i in range(sysInfo.nComponents):
				tmp += sysInfo.comInfoAll[i].cPR * sysInfo.solX[i];
				tmp += (sysInfo.comInfoAll[i].cCR - sysInfo.comInfoAll[i].cPR) *sysInfo.comInfoAll[i].initFail 
			if tmp > 0:
				tmp += sysInfo.cS;
			totalCost += tmp;
			#move to next stage
			age = [sysInfo.comInfoAll[i].initAge*(1-sysInfo.solX[i]) for i in range(sysInfo.nComponents)];	#age after maintenance
			failProb = [sysInfo.comInfoAll[i].cond_fail_prob(age[i], age[i]+1) for i in range(sysInfo.nComponents)];
			#print (tmp);
			#print(age);
			#print (failProb);
			randProb = [random.uniform(0,1) for i in range(sysInfo.nComponents)];
			kesi = [int(randProb[i]<failProb[i]) for i in range(sysInfo.nComponents)];
			#print (kesi);
			for i in range(sysInfo.nComponents):
				sysInfo.comInfoAll[i].initFail = kesi[i];	
				sysInfo.comInfoAll[i].initAge = age[i] + 1;

	overAllCost.append(totalCost);
print ("all costs");
print (overAllCost);
print ("expected cost");	
print (np.mean(overAllCost));
print ("variance");	
print (np.var(overAllCost));
print ("max");
print (max(overAllCost));
print ("min");
print (min(overAllCost));
	
