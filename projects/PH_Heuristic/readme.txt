This code trys to implement progressive hedging based heuristc algorithm to solve stochastic grouping problem in multi-component system.
Which is the Algorithm 4 in paper:
xxxxxxxxx

Code author: Zhicheng Zhu
Email: zhicheng.zhu@ttu.edu, zhich.zhu@gmail.com
Last Update: 06/06/2019


Code structure:
===============
./old_version/:
---------- 
	1.zhu1_notbybook: init version. multi-group at each maintenance window. 
	  keep the components when the components ahead and behind are replaced. 
	  Short for: multi-group, keep the middle
	2.zhu2: multi-group, kill the middle
	3.xiang1: single-group, keep the middle
	4.xiang2: single-group, kill the middle
	
./result/:
----------
	result files
./main.py:
-----------------
    main file.
	improve, more efficient version of zhu1_notbybook.

Notice:
=======
1. see paper for more detail.
2. there is another version of PH heuristic in rolling horizon folder. That one is written 
better, easier to debug. So far, these two versions behave the same, and have the same results.
I'll leave this one as a benchmark. 