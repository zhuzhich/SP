Deterministic dxtensive form (DEF) to solve stochastic grouping problem in multi-component system.

For paper: xxxxxxxxx

Code author: Zhicheng Zhu
Email: zzhu3@lamar.edu, zhich.zhu@gmail.com
Last Update: 03/21/2018

Code structure:
===============
./run.py:
--------- 
	main script file
./ef.py:
-------- 
	extensive form model
./ef.dat: 
-------
	extensive form data.
./result/:
----------
	result files

Notice:
=======
1. No script to run multiple cases, i.e, different number of component and/or scenarios and/or time horizon.
2. For each experiment, one only need to amend file ./ef.dat for different parameters.
3. "REMEMBER": number of indidividuals changes along with time time horizon in ef.dat, i.e., R = T + 2

 