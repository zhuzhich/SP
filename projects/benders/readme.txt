Basic benders decomposition to solve stochastic grouping problem in multi-component system.
Which is the Algorithm 1 in paper:
xxxxxxxxx

Code author: Zhicheng Zhu
Email: zzhu3@lamar.edu, zhich.zhu@gmail.com
Last Update: 03/21/2018

Code structure:
===============
./data/:
------- 
	data file for each scenario. Named as ScenNode[x].dat, where [x] is the scenario index, starts from 1.
./linux_run/: 
-------------
	the job script for HPC (high performance computer). Ignore it.
./result/:
----------
	result files
./main.py:
-----------------
    main file for multi-cut implementation.
./genFile.py:
-------------
	generate scenario-specific data file in ./data
./sub.py:
--------------------
	relaxation version of sub-problem implementation.
./master.py:
--------------------
	relaxation version of master problem implementation.

Notice:
=======
1. We are using python-CPLEX to solve master problem by using branch-and-cut.
2. master problem is built by python-CPLEX.
3. B&C is embeded in CPLEX. So master problem is solved as integer problem.
4. Sub-problem is solved via Pyomo package, solver is CPLEX
 