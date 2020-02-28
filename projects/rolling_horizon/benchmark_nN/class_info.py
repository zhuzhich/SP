#Author: Zhicheng Zhu
#Email: zhicheng.zhu@ttu.edu, yisha.xiang@ttu.edu
#
#class information
#
#Last Update: Feb.2020
#

import math
import random
import numpy as np
import cplex
import copy
from scipy.integrate import quad
class system_information():
	
	def updateConvergeCondition(self):
		self.lastCostSys = copy.deepcopy(self.costSys);
		for i in range(self.I):
			self.comInfoAll[i].lastpVec = \
				copy.deepcopy(self.comInfoAll[i].pVec);
			self.comInfoAll[i].lastPhiRep = \
				copy.deepcopy(self.comInfoAll[i].phiRep);
			self.comInfoAll[i].lastRule = \
				copy.deepcopy(self.comInfoAll[i].rule);
		
	def checkConverge(self):
		if self.lastCostSys == -1:
			self.updateConvergeCondition();
			return False;
		ret = [0]*(self.I+1);
		if abs(self.lastCostSys - self.costSys) < 0.01:
			ret[-1] = 1;
		for i in range(self.I):
			if self.comInfoAll[i].checkConverge() == True:
				ret[i] = 1;
		self.updateConvergeCondition();
		if sum(ret) < len(ret):
			return False;
		else:
			return True;

	def calSysCost(self):
		cost = 0;
		for i in range(self.I):
			cost = cost + self.comInfoAll[i].ind*self.comInfoAll[i].g;
		cost1 = 1;
		tmp = 1;
		for i in range(self.I):
			tmp = tmp*((1-self.comInfoAll[i].phiRep)**self.comInfoAll[i].ind);
		cost1 += -1*tmp;
		tmp = 0;
		for i in range(self.I):
			tmp1 = 1;
			for ii in range(self.I):
				if i == ii:
					continue;
				tmp1 = tmp1*((1-self.comInfoAll[i].phiRep)**self.comInfoAll[i].ind)
			tmp = tmp+self.comInfoAll[i].ind*self.comInfoAll[i].phiRep\
					*((1-self.comInfoAll[i].phiRep)**(self.comInfoAll[i].ind-1))\
					*tmp1;
		cost1 += -1*tmp;
		cost1 = self.cs*cost1;
		for i in range(self.I):
			tmp1 = (1-self.comInfoAll[i].phiRep)**self.comInfoAll[i].ind;
			tmp2 = (1-self.comInfoAll[i].phiRep)**(self.comInfoAll[i].ind-1);
			tmp3 = 1-tmp1-self.comInfoAll[i].ind*self.comInfoAll[i].phiRep*tmp2;
			cost1 += self.comInfoAll[i].cs*tmp3;
		cost += cost1;
		self.costSys = cost;
		


	def policyIter(self):
		for i in range(self.I):
			self.comInfoAll[i].policyIter();

	def calpVec(self):
		#calculate pVector
		#only works for csLevel = 2
		for i in range(self.I):
			#p00
			tmp = (1-self.comInfoAll[i].phiK[0])**(self.comInfoAll[i].ind-1);
			for j in range(self.I):
				if j == i:
					continue;
				tmp = tmp*((1-self.comInfoAll[j].phiK[0])**(self.comInfoAll[j].ind));
			self.comInfoAll[i].pVec[0][0] = tmp;
			#p10
			tmp = (1-self.comInfoAll[i].phiK[1])**(self.comInfoAll[i].ind-1);
			#tmp = 1;
			for j in range(self.I):
				if j == i:
					continue;
				tmp = tmp*((1-self.comInfoAll[j].phiK[1])**(self.comInfoAll[j].ind));
			self.comInfoAll[i].pVec[1][0] = tmp;
			#p01
			self.comInfoAll[i].pVec[0][1] = 1-self.comInfoAll[i].pVec[0][0];
			#p11
			self.comInfoAll[i].pVec[1][1] = 1-self.comInfoAll[i].pVec[1][0];
			self.comInfoAll[i].pVec[1][1] = 1-self.comInfoAll[i].pVec[1][0];

	def calAggInfo(self):
		#calculate
		#theta, phiRep, phiK
		#this is component level, let each component runs itself.
		for i in range(self.I):
			self.comInfoAll[i].calAggInfo(self.csLevel);

	def print_data(self):
		print ("I", self.I);
		print ("cs", self.cs);
		print ("csLevel", self.csLevel);
		print ("costSys", self.costSys);
		print ("lastCostSys", self.lastCostSys);
		
	def __init__(self):
		
		#not fixed parameters
		self.I = 0;				#number of components
		self.cs = 0;				#setup cost
		self.csLevel = 0;
		#extra random seed for residual life
		self.comInfoAll = [];
		#convergence check
		self.costSys = -1;
		self.lastCostSys = -1;
		
#########################################################################	
class component_information():
	def checkConverge(self):
		if len(self.lastpVec) == 0 or len(self.lastRule) == 0 \
			or self.lastPhiRep == -1:
			return False;
		tmp = 0;
		ret = [0]*3;		#check 3 items.
		for i in range(len(self.pVec)):
			for j in range(len(self.pVec[0])):
				tmp = tmp + abs(self.lastpVec[i][j] - self.pVec[i][j]);
		if tmp < 0.01:
			ret[0] = 1;
		tmp = 0;
		for i in range(len(self.rule)):
			tmp = tmp + abs(self.lastRule[i] - self.rule[i]);
		if tmp < 0.01:
			ret[1] = 1;
		if abs(self.phiRep-self.lastPhiRep) < 0.01:
			ret[2] = 1;
		if sum(ret) < len(ret):
			return False;
		else:
			return True;
	
	def TQ(self, condition, discount, action):
		ret = 0;
		cost = 0;
		if action == 1 and condition >= 1:
			cost = self.costRep[condition] - self.V[discount];
		elif action == 0 and condition < self.state - 1:
			cost = self.costOp[condition];
		ret = ret + cost - self.Ts[action]*self.g;
		csLevel = 2;				#fix to 3;
		for j in range(self.state):
			for k in range(csLevel):
				tmp = 0;
				if action == 0 and condition <= self.state - 2\
					and j >= condition:
					tmp = self.pVec[int(j>=self.rule[0])][k]*self.q[condition][j];
				elif action == 1 and condition >= 1 and j == 0:
					tmp = self.pVec[int(j>=self.rule[0])][k];
				ret = ret + tmp*self.value[j][k];
		return ret;
					
		

	def policyIter(self):
		#policy iteration
		#it for csLevel = 2;
		#step 0. init self.value;
		csLevel = 2;			#fixed!
		self.value = [[-1]*csLevel for i in range(self.state)];
		#step 1
		
		#use cplex
		cpx = cplex.Cplex();
		cpx.objective.set_sense(cpx.objective.sense.minimize);
		
		#decision variable
		#v[i][0]
		varNameX = [];
		varX = []
		for i in range(self.rule[0]):
			name = "x"+str(i);
			varNameX.append(name);
			varX.append(cpx.variables.get_num());
			if i == 0:
				cpx.variables.add(obj = [0],\
						lb=[0.0], ub=[0.0],\
						#types=["C"],\
						names=[varNameX[-1]]);
			else:
				cpx.variables.add(obj = [0],\
						lb=[0.0], ub=[100000.0],\
						#types=["C"],\
						names=[varNameX[-1]]);
		#g
		varNameG = [];
		varG = []
		name = "g";
		varNameG.append(name);
		varG.append(cpx.variables.get_num());
		cpx.variables.add(obj = [0],\
			lb=[0.0], ub=[100000.0],\
			#types=["C"],\
			names=[varNameG[-1]]);
		
		#constraint.
		for i in range(self.rule[0]):
			nameVec = copy.deepcopy(varNameX);
			coeVec = [0]*self.rule[0];
			coeVec[i] += 1;					#add v[i][0]
			nameVec.append(varNameG[0])		#add g.
			coeVec.append(1);
			for k in range(csLevel):
				for j in range(i, self.rule[k]):
					coeVec[j] = coeVec[j] - self.q[i][j]*self.pVec[0][k];
			const = 0;
			for k in range(csLevel):
				for j in range(max(self.rule[k],i), self.state):
					const = const + self.q[i][j]*self.pVec[int(j>=self.rule[0])][k]*\
									(self.costRep[j]-self.V[k]);
			const += self.costOp[i];
			cpx.linear_constraints.add(
				lin_expr=[cplex.SparsePair(nameVec, coeVec)],
				senses=["E"], 
				range_values=[0.0],
				rhs=[const]);
		
		#solve it
		cpx.solve();
		#get result
		solution = cpx.solution;
		solutionAll = solution.get_values();
		#assign values
		for i in range(self.rule[0]):
			self.value[i][0] = solutionAll[i];
		self.g = solutionAll[-1];
		for k in range(csLevel):
			for i in range(self.state):
				if i >= self.rule[k]:
					self.value[i][k] = self.costRep[i]-self.V[k];
				else:
					self.value[i][k] = self.value[i][0];
		#step 2
		
		#flag = True;
		#count = 0;
		#while flag == True and count < 20:
			#count += 1;
		flag = False;
		#decrease:
		#k = 2;
		#for i in range(self.rule[k]):
		#	if self.TQ(i,k,1) < self.value[i][k]:
		#		flag = True;
		#		self.rule[k] = i;
		#		break;
			
		k = 1;
		for i in range(1, self.rule[k]):
			if round(self.TQ(i,k,1),4) <= round(self.value[i][k],4):
				flag = True;
				self.rule[k] = i;
				break;
		if flag == False:
			#increase:
			k = 1;
			for i in range(self.rule[k], self.rule[k-1]):
				ii = self.rule[k-1] - 1 - i + self.rule[k];
				if round(self.TQ(ii,k,0),4) <= round(self.value[i][k],4):
					self.rule[k] = ii + 1;
					flag = True;
					break;
	
	def calAggInfo(self, csLevel):
		#theta, phiRep, phiK
		#given current rule and pVector
		#use cplex
		cpx = cplex.Cplex();
		cpx.objective.set_sense(cpx.objective.sense.minimize);
		varNameX = [];	#it's theta
		varX = [];
		#decision variables
		#if types is specified, then it will be a MIP
		for j in range(self.state):
			varNameX.append("x"+str(j));
			varX.append(cpx.variables.get_num());
			cpx.variables.add(obj = [0],\
								lb=[0.0], ub=[1.0],\
								#types=["C"],\
								names=[varNameX[-1]]);
		#constraint 1
		for j in range(1, self.state):
			nameVec = [];
			coeVec = [];
			for i in range(self.state):
				nameVec.append(varNameX[i]);
				tmpCoe = 0;
				if i == j:
					tmpCoe = -1;
				#a version only works for csLevel = 2
				if i < self.rule[1]:
					tmpCoe += self.q[i][j];
				elif i >= self.rule[0]:
					tmpCoe += self.q[0][j];
				else:
					tmpCoe = tmpCoe + self.pVec[0][0]*self.q[i][j]\
							+ self.pVec[0][1]*self.q[0][j];
				#a general version.
				'''
				ii = csLevel;
				for k in range(csLevel):
					if k == 0:
						tmp = (i in range(self.rule[k], self.state));
					else:
						tmp = (i in range(self.rule[k], self.rule[k-1]));
					if tmp == True:
						ii = k;	#identify k
				if ii == 0:
					tmpCoe += self.q[0][j];
				elif ii == csLevel:
					tmpCoe += self.q[i][j];
				else:
					for k in range(csLevel):
						if k < ii:
							tmpCoe = tmpCoe + self.q[i][j]*self.pVec[0][k];
						else:
							tmpCoe = tmpCoe + self.q[0][j]*self.pVec[0][k];
				'''
				coeVec.append(tmpCoe);
			cpx.linear_constraints.add(
				lin_expr=[cplex.SparsePair(nameVec, coeVec)],
				senses=["E"], 
				range_values=[0.0],
				rhs=[0.0]);
		#constraint 2
		coeVec = [1]*self.state;
		cpx.linear_constraints.add(
				lin_expr=[cplex.SparsePair(varNameX, coeVec)],
				senses=["E"], 
				range_values=[0.0],
				rhs=[1.0]);
		#solve it
		cpx.solve();
		#get result
		solution = cpx.solution;
		solutionX = solution.get_values();
		self.theta = [];
		for i in range(self.state):
			self.theta.append(solutionX[i]);
		
		self.phiRep = 0;
		self.phiK = [0]*csLevel;
		for k in range(csLevel):
			for i in range(self.rule[k], self.state):
				self.phiK[k] += self.theta[i];
				if i >= self.rule[0]:
					self.phiRep = self.phiRep + self.pVec[1][k]*self.theta[i];
				else:
					self.phiRep = self.phiRep + self.pVec[0][k]*self.theta[i];
	def weibull_pdf(self, t):
		shape = float(self.wShape)
		scale = float(self.wScale)
		t = float(t)
		part1 = shape/scale
		part2 = (t/scale)**(shape-1)
		tmp = (t/scale)**shape
		part3 = math.exp(-tmp)
		return part1*part2*part3

		

	def weibull_cdf(self, t):
		tmp = float(t)/self.wScale;
		tmp1 = tmp ** self.wShape;
		return 1 - math.exp(-1*tmp1);
	def weibull_pdf_tr(self,t):
		return self.weibull_pdf(t)/self.weibull_cdf(self.state-1);
	def weibull_int(self, t):
		t = float(t)
		return t*self.weibull_pdf_tr(t)
	def weibull_cdf_tr(self, t):
		return self.weibull_cdf(t)/self.weibull_cdf(self.state - 1);
	def cost(self, t):
		fx = self.weibull_cdf_tr(t)
		result = self.costCr*fx + self.costPr*(1-fx)
		return result							#return cost
	def cycle_length(self, t):
		x = t									#use x here, don't wanna make confusion.
		rx = 1 - self.weibull_cdf_tr(x)
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

	def init_rule(self):
		optC = 10000000;
		optT = -1;
		for t in range(1,self.state-1):
			c = self.cost_rate(t);
			if c < optC:
				optT = t;
		self.rule = [t]*2;
	
	
		
	def init_q(self):
		#get transition probability matrix
		#1. get the 80% CDF life.
		designTime = 0.95;
		#if the reliability is less than 1-designTime, 
		#consider it failed.
		#assume the inspection interval is 1;
		prob = 0;
		state = -1;
		while prob < designTime:
			state += 1;
			prob = self.weibull_cdf(state);
		self.state = state + 1;
		self.q = [[0]*self.state for i in range(self.state)];
		for i in range(self.state):
			if i >= self.state - 2 :
				self.q[i][-1] = 1;
				continue;
			self.q[i][-1] = (self.weibull_cdf_tr(i+1) - self.weibull_cdf_tr(i))/(1-self.weibull_cdf_tr(i));
			self.q[i][i+1] = 1-self.q[i][-1];
	def init_parameters(self):
		#cr, shape, scale
		random.seed(self.index*30);
		self.costCr = round(random.uniform(self.cCrBound[0], self.cCrBound[1]), 1);
		random.seed(self.index*20);
		self.wShape = round(random.uniform(self.wShapeBound[0], self.wShapeBound[1]), 1);
		random.seed(self.index*10);
		self.wScale = round(random.uniform(self.wScaleBound[0], self.wScaleBound[1]), 1);
	
	def repair_cost(self):
		self.costRep = [self.costPr+self.cs]*self.state;
		self.costRep[-1] = self.costRep[-1] - self.costPr + self.costCr;
		self.costRep[0] = 1000000;
		'''
		self.costRep = [0]*self.state;
		for i in range(self.state):
			if i == 0:
				self.costRep[i] = 10000000;	#make it really large.
											#it is not supposed to be used
			else:
				self.costRep[i] = 280 + 10*(i+1);
		'''
		
	def operational_cost(self):
		self.costOp = [0]*self.state;
		self.costOp[-1] = 1000000;
		'''
		for i in range(self.state):
			if i == self.state - 1:
				self.costOp[i] = 10000000;
			else:
				self.costOp[i] = 5+10*(i+1);
		'''
		
#######	
	def print_data1(self):
		print ("index", self.index);
		print ("ind", self.ind);
		print ("state", self.state);
		print ("cs", self.cs);
		print ("costRep", self.costRep);
		print ("costOp", self.costOp);
		print ("wShape", self.wShape);
		print ("wScale", self.wScale);
		print ("costCr", self.costCr);
		
	def print_data2(self):
		print ("V", self.V);
		print ("transition matrix", self.q);
		print ("rule,", self.rule);
		print ("pVec", self.pVec);
		print ("theta", self.theta);
		print ("phiRep", self.phiRep);
		print ("phiK", self.phiK);
		print ("value", self.value);
		print ("g", self.g);
		print ("lastpVec", self.lastpVec);
		print ("lastRule", self.lastRule);
		print ("lastPhiRep", self.lastPhiRep);
	def __init__(self):
		self.index = 0		#index
		self.ind = 0;		#number of individuals
		self.state = 0;		#number of states
		#age parameters;
		self.wShapeBound = 0;
		self.wScaleBound = 0;
		self.cCrBound = 0;
		self.wShape = 0;
		self.wScale = 0;
		#cost
		self.cs = 0;		#system level setup cost
		self.costPr = 0;	#pr cost
		self.costCr = 0;	#cr cost;
		self.costRep = 0;	#Repairing cost. [1*state]
		self.costOp = 0;	#Operational cost. [1*state]		
		#running
		self.q = [];		#state transition probability matrix.[state*state]
		self.Ts = [1, 0];	#action = 0 -> 1, action=1 -> 0;
		self.rule = 0;		#current rule. [1*csLevel]
		self.pVec = 0;		#probability vector. [1*csLevel]*2
		self.theta = 0;		#prob. of staying at each state given current rule. [1*state]
		self.phiRep = 0;	#prob. of repairing given current rule.
		self.phiK = 0;		#prob. of staying above state rule[k]. [1*csLevel];
		self.value = 0;		#value function. [state*csLevel]
		self.V = [];		#discounts. [1*csLevel]
		self.g = 0;			#long-run average cost; A decision variable
		#convergence check;
		self.lastpVec = [];
		self.lastRule = [];
		self.lastPhiRep = -1;
		