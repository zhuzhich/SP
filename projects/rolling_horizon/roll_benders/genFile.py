# Author Zhicheng Zhu
# Email: zzhu3@lamar.edu

# For files generation for Algorithm 1 or Algorithm 2
# Files generated
# 1. master.dat: generate master.dat file
# 2. ScenNodex.dat: subproblem x data file

import os
import math
import random


#create master.dat
def create_masterDataFile():    
	global I
	global w
	global d
	global directory
	
	#local variable
	I_Lo = I
	d_Lo = d
	w_Lo = w
	dir_local = directory	
	
	file_name = dir_local + "\\data\\master.dat"
	f = open(file_name,"w")
	f.write("#Total component number\n")
	f.write("param pI := " + str(I_Lo)+";\n\n")

	f.write("#Number of Scenarios\n")
	f.write("param NUMSCEN := " + str(w_Lo)+";\n\n")

	f.write("#Probablity of each Scenario. Equally distributed\n")
	f.write("param prob := " + str(1.0/w_Lo)+";\n\n")	

	f.write("#Setup cost\n")
	f.write("param pd := " + str(d_Lo)+";\n\n")

	f.close()

#create subproblem data file: subDataFile
def create_subDataFile():    
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
	global x_last
	

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
	#pKesi_Lo = pKesi
	w_shape_Lo = w_shape
	w_scale_Lo = w_scale
	#x_Lo = x
	dir_local = directory
	x_last_Lo = x_last
	is_first_time = False
	if x_last_Lo[0] == -1:
		is_first_time = True
	for idx_w in range(0,w_Lo):
		file_name = dir_local + "\\data\\ScenNode"+str(idx_w+1)+".dat"
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
		text6 = "param x := \n"
		#text3 and text6 will be done after lifetime creation.
		for idx in range(0,I_Lo):
			i = idx + 1
			text1 = text1 + str(i) + " " + str(pCPR_Lo[idx]) + "\n"
			text2 = text2 + str(i) + " " + str(pCCR_Lo[idx]) + "\n"
			#text3 = text3 + str(i) + " " + str(pKesi_Lo[idx]) + "\n"
			text4 = text4 + str(i) + " " + str(w_shape_Lo[idx])	 + "\n"
			text5 = text5 + str(i) + " " + str(w_scale_Lo[idx]) + "\n"
			#text6 = text6 + str(i) + " " + str(x_Lo[idx]) + "\n"

		text1 = text1 + ";\n\n"
		text2 = text2 + ";\n\n"
		#text3 = text3 + ";\n\n"
		text4 = text4 + ";\n\n"
		text5 = text5 + ";\n\n"
		#text6 = text6 + ";\n\n"
		f.writelines(text1)
		f.writelines(text2)
		#f.writelines(text3)
		f.writelines(text4)
		f.writelines(text5)
		#f.writelines(text6)
		#Begin to write life time into file
		f.write("#pLT is two-dimensional [i,r]=[component, individual]\n#[i,*] means all data in row i.\n\n")
		f.write("param pLT := ")
		for idx1 in range(0,I_Lo):
			i = idx1 + 1
			f.write("\n["+str(i)+",*]\n")
			for idx2 in range(0,R_Lo):
				r = idx2 + 1
				random.seed(i+r+idx_w+1) ###control the seed
				if r == 1:
					if is_first_time and i == 1: #first time. Decide initial state. Here, only first component fails.
						LT = 0
					else:    
						ran_num = random.uniform(0,1)
						s_inv = 1.0/w_shape_Lo[idx1]	
						if x[idx1] == 0:				#it's not replaced on last time epoch
							tmp = (float(s_Lo)/w_scale_Lo[idx1])**w_shape_Lo[idx1]
							ran_num_log = -math.log(ran_num*math.exp(-tmp)) 
							LT1 = int((ran_num_log**s_inv)*w_scale_Lo[idx1]) - s_Lo			
						else:
							ran_num_log = -math.log(ran_num)
							LT1 = int((ran_num_log**s_inv)*w_scale_Lo[idx1])
							if is_first_time:
								LT = max(1,LT1)	
							else:
								LT = max(0,LT1)	
					if LT == 0:
						pKesi[idx1] = 1
						x[idx1] = 1
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
		f.write("\n")		
		for idx in range(0,I_Lo):
			i = idx + 1
			text3 = text3 + str(i) + " " + str(pKesi[idx]) + "\n"
			text6 = text6 + str(i) + " " + str(x[idx]) + "\n"
		text3 = text3 + ";\n\n"
		text6 = text6 + ";\n\n"
		f.writelines(text3)
		f.writelines(text6)
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
def x_init(x):
	global I
	for i in range(0, I):
		idx = i + 1
		temp = 0
		if idx == 1: #or idx == 2 or idx == 3:
			temp = 1
		x.append(temp)

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
global x_last

def create_allFiles(I_in, T_in, w_in, x_last_in, t_roll_in, dir):		 	
#############################
#Step 1: set up parameters
#############################
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
	
	global x_last

	
	directory = dir    
	# ./main.ph
	# ./models
	# ./nodedata
	#print ("111111111111111111")
	#counter += 1
	#if w == 100:
	#	continue
	#if T == 30:
	#	continue
	I = I_in				#number of components

	w = w_in				#number of scenarios
	x_last = x_last_in		#replacement at last epoch
	T = T_in - t_roll_in	#time horizon

	T_ex = T + 40	#extended time horizon
	s = t_roll_in	#starting time
	R = T + 2		#number of individuals
	d = 5			#setup cost

	pCPR = []
	pCPR_init(pCPR)

	pCCR = []
	pCCR_init(pCCR)

	pKesi = [0]*I
	#pKesi_init(pKesi) #move to lifetime creation.

	w_shape = []
	w_shape_init(w_shape)	

	w_scale = []
	w_scale_init(w_scale)	

	x = [0]*I 			#should be integer feasible
	#x_init(x)			#move to lifetime creation
	
	create_masterDataFile() 
	create_subDataFile()

	#############################
	#Use command line for runph
	#############################
	#res_file_path = directory + "\\result"
	#res_filename = str(counter) + "_" + str(I) + "_" + str(T) + "_" + str(w) + ".txt"
	#res_file = res_file_path + "\\" + res_filename
	#excmd = "python run.py > " + res_file
	#excmd = "python run.py"
	#os.system(excmd)
		






