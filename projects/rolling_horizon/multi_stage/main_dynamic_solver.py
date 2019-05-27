#Author: Zhicheng Zhu
#Email: zhicheng.zhu@ttu.edu, yisha.xiang@ttu.edu
#copyright @ 2019: Zhicheng Zhu. All right reserved.

#Info:
#main file to solve deterministic equivalent form of TBM model by back-tracking
#use enumeration
#Last update: 05/20/2019
#!/usr/bin/python

from __future__ import print_function
import sys
import cplex
import itertools
import time
from scipy.stats import gamma

def node_iter(sysInfo, curT, age, failState):
	#curT: current time
	#age: age after mx, vector in length of n
	bestObj = float("inf");
	bestSol = [-1]*sysInfo.nComponents;
	if curT == sysInfo.nStages-1:
		#the last stage
		tmpObj = 0;
		for i in range(sysInfo.nComponents):
			tmpObj += sysInfo.comInfoAll[i].cCR*failState[i];
		if tmpObj > 0:
			tmpObj += sysInfo.cS;
		bestObj = tmpObj;
	else:
		#not the last stage
		for solX in itertools.product([0,1], repeat = sysInfo.nComponents):
			solXL= list(solX);
			flag = False;
			ageAfter = [0]*sysInfo.nComponents;
			for i in range(sysInfo.nComponents):
				ageAfter[i] = age[i]*(1 - solXL[i]);
				if failState[i] > solXL[i]:
					flag = True;
			if flag == True:
				continue;
			tmpObj = 0;
			#calculate current node cost given solution solXL
			for i in range(sysInfo.nComponents):
				tmpObj += solXL[i]*sysInfo.comInfoAll[i].cPR;
				tmpObj += (sysInfo.comInfoAll[i].cCR - sysInfo.comInfoAll[i].cPR)*failState[i];
			if tmpObj > 0:
				tmpObj += sysInfo.cS;
			#print ("======t:===========");
			#print (curT);
			#print ("node:");
			#print (failState);
			#print ("solution:");
			#print (solX);
			
			#calculate expected cost at the next stage given solution solXL
			for kesi in itertools.product([0,1], repeat = sysInfo.nComponents):
				ageNext = [ageAfter[i]+1 for i in range(sysInfo.nComponents)];
				tmpCost = node_iter(sysInfo, curT+1, ageNext, kesi);
				tmpProb = 1;
				#print ("-------t:-------");
				#print (curT+1);
				#print ("node:");
				#print (kesi);
				#print ("cost:");
				#print (tmpCost);

				for i in range(sysInfo.nComponents):
					tmp = sysInfo.comInfoAll[i].cond_fail_prob(ageAfter[i],ageAfter[i]+1);
					tmpProb = tmpProb*(tmp*kesi[i] + (1 - tmp)*(1-kesi[i]));
				#print ("probability:");	
				#print (tmpProb);
				tmpObj += tmpProb*tmpCost;
				#print ("objective:");
				#print (tmpObj);				
			if tmpObj < bestObj:
				bestObj = tmpObj;
				bestSol = solXL;
	if curT == 0:
		res = [];
		res.append(bestObj);
		res.append(bestSol);
		return res;
	else:
		return bestObj;



def main(sysInfo):

	start_time = time.clock();

	age = [sysInfo.comInfoAll[i].initAge for i in range(sysInfo.nComponents)];
	failState = [sysInfo.comInfoAll[i].initFail for i in range(sysInfo.nComponents)];
	res = node_iter(sysInfo, 0, age, failState);
	bestObj0 = res[0];
	bestSol0 = res[1];
	'''
	## stage 1
	#generate first stage solution:
	bestObj0 = float("inf");
	bestSol0 = 0;
	for solX0 in itertools.product([0,1], repeat = sysInfo.nComponents):
		solXL_0 = list(solX0);
		flag1 = False;
		ageAfterMx_0 = [0]*sysInfo.nComponents;
		for i in range(sysInfo.nComponents):
			ageAfterMx_0[i] = sysInfo.comInfoAll[i].initAge*(1 - solXL_0[i]);
			if sysInfo.comInfoAll[i].initFail > solXL_0[i]:
				flag1 = True;
		if flag1 == True:
			continue;
		tmpObj0 = 0;
		for i in range(sysInfo.nComponents):
			tmpObj0 += solXL_0[i]*sysInfo.comInfoAll[i].cPR;
			tmpObj0 += (sysInfo.comInfoAll[i].cCR - sysInfo.comInfoAll[i].cPR)*sysInfo.comInfoAll[i].initFail;
		if tmpObj0 > 0:
			tmpObj0 += sysInfo.cS;		
	##stage 2
		#generate failure states at the second stage:
		for kesi1 in itertools.product([0,1], repeat = sysInfo.nComponents):
			#generate solution at the second stage:
			bestObj1 = float("inf");	#the best obj value for the current node
			for solX1 in itertools.product([0,1], repeat = sysInfo.nComponents):
				solXL_1 = list(solX1);
				flag1 = False;
				ageAfterMx_1 = [0]*sysInfo.nComponents;
				for i in range(sysInfo.nComponents):
					ageAfterMx_1[i] = ageAfterMx_0[i] + 1;
					ageAfterMx_1[i] = ageAfterMx_1[i]*(1 - solXL_1[i]);				
					if kesi1[i] > solXL_1[i]:
						flag1 = True;
				if flag1 == True:
					continue;			
				tmpObj1 = 0;	
				for i in range(sysInfo.nComponents):
					tmpObj1 += solXL_1[i]*sysInfo.comInfoAll[i].cPR;
					tmpObj1 += (sysInfo.comInfoAll[i].cCR - sysInfo.comInfoAll[i].cPR)*kesi1[i];
				if tmpObj1 > 0:
					tmpObj1 += sysInfo.cS;	
				#generate failure states at the third (last) stage:
	## stage 3
				for kesi2 in itertools.product([0,1], repeat = sysInfo.nComponents):
					tmpCost = 0;
					for i in range(sysInfo.nComponents):
						tmpCost += sysInfo.comInfoAll[i].cCR*kesi2[i];
					if tmpCost > 0:
						tmpCost += sysInfo.cS;
						tmpProb = 1;
						for i in range(sysInfo.nComponents):
							tmp = sysInfo.comInfoAll[i].cond_fail_prob(ageAfterMx_1[i],ageAfterMx_1[i]+1);
							tmpProb = tmpProb*(tmp*kesi2[i] + (1 - tmp)*(1-kesi2[i]));
						tmpObj1 += tmpProb*tmpCost;
				if tmpObj1 < bestObj1:
					bestObj1 = tmpObj1;
					
			tmpProb = 1;
			for i in range(sysInfo.nComponents):
				tmp = sysInfo.comInfoAll[i].cond_fail_prob(ageAfterMx_0[i],ageAfterMx_0[i]+1);
				tmpProb = tmpProb*(tmp*kesi1[i] + (1 - tmp)*(1-kesi1[i]));
			tmpObj0 += tmpProb * bestObj1;
		if tmpObj0 < bestObj0:
			bestObj0 = tmpObj0;
			bestSol0 = solXL_0;
	'''	

	end_time = time.clock();

	time_elapsed = end_time - start_time;
	sysInfo.time = time_elapsed;


	#print ("\n===============================main_2stage_solver, (m, n, t)=(%d,%d,%d)============" 
	#		%(sysInfo.comInfoAll[0].nStates, sysInfo.nComponents, sysInfo.nStages));

	print ("calculation time is %f"  %time_elapsed);


	sysInfo.objValue = bestObj0;
	print ("optimal objValue");
	print (bestObj0);
	print ("optimal solution");
	print (bestSol0);
		
		