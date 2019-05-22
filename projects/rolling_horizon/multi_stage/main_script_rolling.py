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
nComponents = 3;		#fix component numbers	
cS = 5;				#fix setup cost
intvl = 1;		#don't change this
nStages = 7;			

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
	random.seed(i*10)   
	temp = random.uniform(4,7)
	w_shape[i] = round(temp,1);			
	#scale
	random.seed(i*10)   
	temp = random.uniform(1,4)#(4,11)
	w_scale[i] = round(temp,1);		
	comInfo = class_info.component_info(i, w_shape[i], w_scale[i], age[i], kesi[i], intvl, cCR[i], cPR[i], cS);
	sysInfo.add_com(comInfo);


sysInfo1 = copy.deepcopy(sysInfo);
main_dynamic_solver.main(sysInfo1);
	
