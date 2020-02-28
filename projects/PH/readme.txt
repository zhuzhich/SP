This code trys to implement progressive hedging algorithm to solve stochastic grouping problem in multi-component system.
Which is the Algorithm 3 in paper:
xxxxxxxxx

Code author: Zhicheng Zhu
Email: zzhu3@lamar.edu, zhich.zhu@gmail.com
Last Update: 03/21/2018


Code structure:
===============
./backup/:
---------- 
	some backup files. Ignore.
./models/:
--------- 
	ReferenceModel.py: problem model. It has to be named as this because of "runph" script
./nodedata/: 
-------------
	ScnearioStructure.dat: register file. Register user-defined file to the hook of "runph" script.
	RootNode.dat: shared parameters across different scenarios.
	ScenNode[x].dat: scenario-specific parameters
./result/:
----------
	result files
./main.py:
-----------------
    main file.
	create files in ./nodedata/ .
	call runph script.

Notice:
=======
1. runph script from Pyomo/PySP package is used for solving standard progressive hedging problem.
2. runph requires:
	2.1 a model file: ./models/ReferenceModel.py
	2.2 a register file: ./nodedata/ScnearioStructure.dat
	2.3 some user-defined files for different scenarios
3. For more questions, one can refer http://www.pyomo.org/documentation/
 