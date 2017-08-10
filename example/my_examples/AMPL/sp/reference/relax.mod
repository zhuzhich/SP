param I; 		#total number of components
param R; 		#max number of individuals
param W;		#scenarios
param delta; 	#interval
#param S; 		#end of time horizon, not using for now
param s; 		#start of time horizon
param d;		#set up cost

param c{1..I};  		#maintenance cost of component i
param w_shape{1..I}; 	#shape parameter of Weibull of component i
param w_scale{1..I};	#scale 

#
param T :=30;
param kesi{1..I};
param LT{i in 1..I, r in 1..R, w in 1..W} default round((-log(Uniform(0,1)))^(1/w_shape[i])*w_scale[i]);	#life time non-identical
param LTI{1..I};        #life time identical

#
#var x_1{1..I,0..T,1..R,1..W} binary >=0; #non-identical variable 
#var x_2{1..I,0..T,1..W} binary >=0;      #identical variable
#var z{0..T,1..W} binary >=0;				#system maintenance or not at t
#var x{1..I} binary >=0;

###relaxation
var x_1{1..I,0..T,1..R,1..W} >=0 , <=1; #non-identical variable 
var x_2{1..I,0..T,1..W} >=0 , <=1;      #identical variable
var z{0..T,1..W} >=0 , <=1;				#system maintenance or not at t
var x{1..I} >=0 , <=1;
var violate{w in 1..W} := 0;
var cost_w{w in 1..W} :=0;
#####
#objective
minimize cost:
sum{w in 1..W}(sum{i in 1..I}(sum{r in 1..R}c[i]*x_1[i,T,r,w] + sum{t in 0..T}c[i]*x_2[i,t,w]) + sum{t in 0..T}d*z[t,w])/W;

#constraints
s.t. cst_b {i in 1..I, t in 0..T-1, r in 1..R,w in 1..W}: 
	 	x_1[i,t,r,w]<=x_1[i,t+1,r,w];
s.t. cst_c {i in 1..I, t in 0..T-1, r in 1..R-1,w in 1..W}: 
		x_1[i,t+1,r+1,w]<=x_1[i,t,r,w];
s.t. cst_d {i in 1..I, t in 1..T,w in 1..W}: 				
		x_2[i,t,w] + sum{r in 1..R}(x_1[i,t,r,w]-x_1[i,t-1,r,w]) <=z[t,w];
s.t. cst_e {i in 1..I,w in 1..W}: 						  
		x_2[i,0,w] + x_1[i,0,1,w] <=z[0,w];
s.t. cst_f {i in 1..I, r in 1..R-1, w in 1..W, t in 0..T-LT[i,r+1,w]}: 	
		x_1[i,t,r,w]<=x_1[i,t+LT[i,r+1,w],r+1,w];
s.t. cst_g {i in 1..I, l in 0..T-LTI[i],w in 1..W}: 		  
		sum{t in l+1..l+LTI[i]}x_2[i,t,w]>=x_1[i,l,R,w];
s.t. cst_h {i in 1..I, t in 1..T,w in 1..W}:				  
		x_2[i,t,w]<=x_1[i,t-1,R,w];
s.t. cst_i {i in 1..I,w in 1..W}:							  
		x_2[i,0,w] = 0;
s.t. cst_j {w in 1..W, i in 1..I:LT[i,1,w]<=T}:				  
		x_1[i,LT[i,1,w],1,w]=1;
s.t. cst_k {i in 1..I, r in 2..R,w in 1..W}:
		x_1[i,0,r,w] = 0;
s.t. cst_l {i in 1..I,w in 1..W}:							  
		x[i]=x_1[i,0,1,w];							  
s.t. cst_m {i in 1..I}:							  
		x[i]>=kesi[i];

