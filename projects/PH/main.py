# Author Zhicheng Zhu
# Email: zzhu3@lamar.edu

# For files generation in folder nodedata, in order to run runph script
# Files generated
# 1. ScenarioStructure.dat: contains scenario tree info for runph
# 2. RootNode.dat: contains the fix info that doesn't change across scenarios
# 3. ScenNode1.dat - ScenNodew.dat, where w is the number of scenarios.
#	 These files contain the scenario specific data, which is the random life times
import os
import math
import random

#create or rewrite ScenarioStructure.dat
def create_scenarioStructure():
	global directory
	global w
	#use the local variable.
	dir_local = directory	
	w_Lo = w

	file_name = dir_local + "\\nodedata\\ScenarioStructure.dat"
	prob = 1.0/w
	#print (file_name)
	f = open(file_name,"w")
	nodes = []
	scenarios = []
	prefix_nodes = "ScenNode"
	prefix_scenarios = "Scen"
	for i in range(0,w_Lo):
		nodes.append(prefix_nodes+str(i+1))
		scenarios.append(prefix_scenarios+str(i+1))
	#print (nodes)
	#print (scenarios)
	stage1 = "FirstStage "
	stage2 = "SecondStage "
	root_node = "RootNode "
	text1 = "set Stages := " + stage1 + stage2 + " ;\n\n"
	#need appending
	text2 = "set Nodes := \n" + root_node + " \n"
	text3 = "param NodeStage :=  \n" + root_node + stage1 + " \n"
	text4 = "set Children[" + root_node.strip() + "]:= \n" 
	text5 = "param ConditionalProbability := \n" + root_node + "1.0 \n"
	text6 = "set Scenarios := \n"
	text7 = "param ScenarioLeafNode := \n"
	for i in range(0,w_Lo):
		text2 = text2 + nodes[i] + "\n"
		text3 = text3 + nodes[i]+ " " + stage2 + "\n"
		text4 = text4 + nodes[i]+ "\n"
		text5 = text5 + nodes[i] + " " + str(prob) + "\n"
		text6 = text6 + scenarios[i] + "\n" 
		text7 = text7 + scenarios[i] + " " + nodes[i] + "\n"

	text2 = text2 + ";\n\n"
	text3 = text3 + ";\n\n"
	text4 = text4 + ";\n\n"
	text5 = text5 + ";\n\n"
	text6 = text6 + ";\n\n"
	text7 = text7 + ";\n\n"

	text8 = "set StageVariables["+stage1.strip()+"] := \n" + "x[*]\n" +"z[*]\n;\n\n"
	text9 = "set StageVariables["+stage2.strip()+"] := \n" + "xs[*,*,*]\n" + "ws[*,*,*]\n" + \
												"zs[*]\n"+"u[*,*,*]\n"+"v[*,*,*]\n;\n\n"   
	#firstStageObj is the first stage cost variable in model, so is the secondStageObj	
	text10 = "#firstStageObj is the first stage cost variable in model, so is the secondStageObj\n"\
			+ "param StageCost := \n"\
			+ stage1 + " firstStageObj\n"\
			+ stage2 + " secondStageObj\n;\n\n"
	text11 = "param ScenarioBasedData := False ;"
	
	f.writelines(text1)
	f.writelines(text2)
	f.writelines(text3)
	f.writelines(text4)
	f.writelines(text5)
	f.writelines(text6)
	f.writelines(text7)
	f.writelines(text8)
	f.writelines(text9)
	f.writelines(text10)
	f.writelines(text11)

	f.close()

#create life time data file: scenNodex.dat
def create_lifeTime():
	global I
	global s
	global R
	global w
	global w_shape
	global w_scale
	global pKesi
	global directory

	#local variable
	I_Lo = I
	s_Lo = s
	R_Lo = R
	w_Lo = w
	w_shape_Lo = w_shape
	w_scale_Lo = w_scale
	pKesi_Lo = pKesi
	dir_local = directory	

	for idx in range(0,w_Lo):
		file_name = dir_local + "\\nodedata\\ScenNode"+str(idx+1)+".dat"
		f = open(file_name,"w")
		f.write("#pLT is two-dimensional [i,r]=[component, individual]\n#[i,*] means all data in row i.\n\n")
		f.write("param pLT := ")
		for idx1 in range(0,I_Lo):
			i = idx1 + 1
			f.write("\n["+str(i)+",*]\n")
			for idx2 in range(0,R_Lo):
				r = idx2 + 1
				random.seed(i+r+idx+1) ###control the seed
				if r == 1:
					if pKesi_Lo[idx1] == 1:
						LT = 0
					else:    
						ran_num = random.uniform(0,1)
						ran_num_log = -math.log(ran_num)
						s_inv = 1.0/w_shape_Lo[idx1]	
						LT1 = int((ran_num_log**s_inv)*w_scale_Lo[idx1]) - s_Lo					
						LT = max(1,LT1)
				else:
					ran_num = random.uniform(0,1)
					ran_num_log = -math.log(ran_num)
					s_inv = 1.0/w_shape_Lo[idx1]
					LT1 = int((ran_num_log**s_inv)*w_scale_Lo[idx1])
					LT = max(1,LT1)
				f.write(str(r)+ " " + " "*(2-len(str(r))) + str(LT)+ "    ")
				#f.write(str(r)+ " " + str(LT)+ "   ")
				if r%10 == 0:
					f.write("\n")
			f.write("\n")
		f.write(";")
		f.close()

#create root node data file: rootNode.dat
def create_rootNode():    
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
	global directory
	
	#local variable
	I_Lo = I
	T_Lo = T
	w_Lo = w
	T_ex_Lo = T_ex
	s_Lo = s
	R_Lo = R
	d_Lo = d
	pCPR_Lo = pCPR
	pCCR_Lo = pCCR	
	pKesi_Lo = pKesi
	w_shape_Lo = w_shape
	w_scale_Lo = w_scale
	dir_local = directory	
	
	file_name = dir_local + "\\nodedata\\RootNode.dat"
	f = open(file_name,"w")
	f.write("#Total component number\n")
	f.write("param pI := " + str(I_Lo)+";\n\n")

	f.write("#Time horizon\n")
	f.write("param pT := " + str(T_Lo)+";\n\n")

	f.write("#Number of Scenarios\n")
	f.write("param NUMSCEN := " + str(w_Lo)+";\n\n")

	f.write("#Probablity of each Scenario. Equally distributed\n")
	f.write("param prob := " + str(1.0/w_Lo)+";\n\n")	

	f.write("#Extended time horizon\n")
	f.write("param pT_ex := " + str(T_ex_Lo)+";\n\n")

	f.write("#Starting time\n")
	f.write("param ps := " + str(s_Lo)+";\n\n")

	f.write("#Number of individuals. R=T+2\n")
	f.write("param pR := " + str(R_Lo)+";\n\n")	

	f.write("#Setup cost\n")
	f.write("param pd := " + str(d_Lo)+";\n\n")

	text1 = "param pCPR := \n"
	text2 = "param pCCR := \n"
	text3 = "param pKesi := \n"
	text4 = "param w_shape := \n"
	text5 = "param w_scale := \n"
	for idx in range(0,I_Lo):
		i = idx + 1
		text1 = text1 + str(i) + " " + str(pCPR_Lo[idx]) + "\n"
		text2 = text2 + str(i) + " " + str(pCCR_Lo[idx]) + "\n"
		text3 = text3 + str(i) + " " + str(pKesi_Lo[idx]) + "\n"
		text4 = text4 + str(i) + " " + str(w_shape_Lo[idx])	 + "\n"
		text5 = text5 + str(i) + " " + str(w_scale_Lo[idx]) + "\n"

	text1 = text1 + ";\n\n"
	text2 = text2 + ";\n\n"
	text3 = text3 + ";\n\n"
	text4 = text4 + ";\n\n"
	text5 = text5 + ";\n\n"
	f.writelines(text1)
	f.writelines(text2)
	f.writelines(text3)
	f.writelines(text4)
	f.writelines(text5)
	f.close()

#cost of PR
def pCPR_init(pCPR):
	global I
	for i in range(0, I):
		pCPR.append(1)
#cost of CR: pCCR[i]=i*2, i=1...I
def pCCR_init(pCCR):
	global I
	for i in range(0, I):
		random.seed((i+1)*30) ###control the seed     
		temp = random.uniform(6,16)
		pCCR.append(round(temp,1))
		#pCCR.append((i+1)*2)
#kesi. Only the first component fails at current time
def pKesi_init(pKesi):
	global I
	for i in range(0, I):
		if i+1 == 1:
			pKesi.append(1)
		else:
			pKesi.append(0)		
#weibull shape parameter: w_shape[i]=6
#u(4,7)
def w_shape_init(w_shape):
	global I
	for i in range(0, I):
		random.seed((i+1)*10) ###control the seed     
		temp = random.uniform(4,7)
		w_shape.append(round(temp,1))
#weibull scale parameter: w_scale[i]=i, i=1...I
#u(1,8)
def w_scale_init(w_scale):
	global I
	for i in range(0, I):
		random.seed((i+1)*20) ###control the seed     
		temp = random.uniform(1,8)
		w_scale.append(round(temp,1))

#############################
#Step 1: set up parameters
#############################
directory = "C:\\Users\\zzhu3\\Documents\\codes\\SP\\projects\\PH"    
# ./main.ph
# ./models
# ./nodedata
comp_list = [4,6,8,10]
time_list = [10]
scen_list = [1000]
counter = 0
for I in comp_list:
	for T in time_list:
		for w in scen_list:
			counter += 1
			#if w == 100:
			#	continue
			#if T == 30:
			#	continue
			#I = 6 			#number of components
			#T = 10			#time horizon
			#w = 50			#number of scenarios

			T_ex = T + 40	#extended time horizon
			s = 2			#starting time
			R = T + 2		#number of individuals
			d = 5			#setup cost

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

			create_scenarioStructure() 
			create_rootNode()
			create_lifeTime()

			#############################
			#Use command line for runph
			############################# 
			#res_filename = str(counter) + "_" + str(I) + "_" + str(T) + "_" + str(w) + ".txt"
			res_filename = str(I) + "_" + str(T) + "_" + str(w) + ".txt"
			excmd = "runph --model-directory=models --instance-directory=nodedata --default-rho=1.0 --solver=cplex > " + res_filename
			os.system(excmd)
			#cd C:\Users\zzhu3\Documents\codes\SP\projects\PH
			#runef --model-directory=models --instance-directory=nodedata --solve --solver=cplex --traceback







