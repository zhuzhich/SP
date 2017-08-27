param I; 		#total number of components
param T :=10 ;
param T_1:=20;
param R :=T+2; 		#max number of individuals
param W;		#scenarios
param delta; 	#interval
#param S; 		#end of time horizon, not using for now
param s; 		#start of time horizon
param d;		#set up cost

param c_pm{1..I};  		#maintenance cost of component i
param c_cm{1..I};

param w_shape{1..I}; 	#shape parameter of Weibull of component i
param w_scale{1..I};	#scale 

param x{1..I};
#

param kesi{1..I};
param LT{i in 1..I, r in 1..R, w in 1..W} default round((-log(Uniform(0,1)))^(1/w_shape[i])*w_scale[i]);	#life time non-identical

#non relaxation
#var x_1{1..I,0..T,1..R,1..W} binary >=0; #non-identical variable 
#var w_x{1..I,0..T_1,1..R,1..W} binary >=0;
#var u{1..I,0..T,1..R,1..W}  binary >=0;
#var v{1..I,0..T,1..R,1..W}  binary >=0;
#var z{0..T,1..W} binary >=0;				#system maintenance or not at t
#var x{1..I} binary >=0;

###relaxation
var x_1{1..I,0..T,1..R,1..W} >=0,<=1; #non-identical variable 
var w_x{1..I,0..T_1,1..R,1..W} >=0,<=1;
var u{1..I,0..T,1..R,1..W} >=0,<=1;
var v{1..I,0..T,1..R,1..W} >=0,<=1;
var z{1..T,1..W} >=0,<=1;				#system maintenance or not at t

#partial relaxation
#var x{1..I} binary >=0;
#
#var violate{w in 1..W} := 0;
#var cost_w{w in 1..W} :=0;


var cost_t;
#####
#objective

minimize cost:
sum{w in 1..W}
(
sum{i in 1..I}
(
c_cm[i]*w_x[i,LT[i,1,w],1,w]+c_pm[i]*(-w_x[i,LT[i,1,w],1,w])+sum{r in 2..R}
(
c_cm[i]*(x_1[i,T,r-1,w]+x_1[i,T,r,w]-sum{t in 0..T}(u[i,t,r,w]+v[i,t,r,w]))/2
+c_pm[i]*0.5*sum{t in 0..T}(u[i,t,r,w]+v[i,t,r,w])
)
)
+sum{t in 1..T}d*z[t,w]
)
/W;

#constraints
s.t. cst_b {i in 1..I, t in 0..T-1, r in 1..R,w in 1..W}: 
	 	x_1[i,t,r,w]<=x_1[i,t+1,r,w];
s.t. cst_c {i in 1..I, t in 0..T-1, r in 1..R-1,w in 1..W}: 
		x_1[i,t+1,r+1,w]<=x_1[i,t,r,w];
s.t. cst_d {i in 1..I, t in 1..T,w in 1..W}: 				
		sum{r in 1..R}(x_1[i,t,r,w]-x_1[i,t-1,r,w]) <=z[t,w];
#s.t. cst_e {i in 1..I,w in 1..W}: 						  
#		x_1[i,0,1,w] <=z[0,w];
s.t. cst_f {i in 1..I, r in 1..R-1, w in 1..W, t in 0..T-LT[i,r+1,w]}: 	
		x_1[i,t,r,w]<=x_1[i,t+LT[i,r+1,w],r+1,w];
s.t. cst_g {w in 1..W, i in 1..I:LT[i,1,w]<=T}:				  
		x_1[i,LT[i,1,w],1,w]=1;
s.t. cst_h {i in 1..I, r in 2..R,w in 1..W}:
		x_1[i,0,r,w] = 0;
s.t. cst_i {i in 1..I,w in 1..W}:							  
		x[i]=x_1[i,0,1,w];							  
#s.t. cst_j {i in 1..I}:							  
#		x[i]>=kesi[i];
		
s.t. cst_k {i in 1..I, r in 1..R, t in 1..T, w in 1..W}:
		w_x[i,t,r,w] = x_1[i,t,r,w]-x_1[i,t-1,r,w];
s.t. cst_l {i in 1..I,r in 1..R, w in 1..W}:
		w_x[i,0,r,w] = x_1[i,0,r,w];
s.t. cst_m {i in 1..I,r in 1..R,  w in 1..W, t in T+1..T+LT[i,r,w]-1}:
		w_x[i,t,r,w] = w_x[i,t-T,r,w];
s.t. cst_n {i in 1..I,r in 1..R, w in 1..W, t in T+LT[i,r,w]..T_1}:
		w_x[i,t,r,w] = 0;	
			
s.t. cst_o {i in 1..I, r in 2..R, w in 1..W,t in 0..T}:
		u[i,t,r,w]-v[i,t,r,w]=w_x[i,t+LT[i,r,w],r,w]-w_x[i,t,r-1,w];
s.t. cst_p {i in 1..I, w in 1..W, t in 0..T}:
		u[i,t,1,w]-v[i,t,1,w] = w_x[i,t+LT[i,1,w],1,w]-1;
s.t. cst_q {i in 1..I, r in 1..R, w in 1..W,t in 0..T}:
		u[i,t,r,w]+v[i,t,r,w] <= 1;

#s.t.  a1: x_1[1,0,1,1] = 1;
#s.t.  a2: x_1[1,3,2,1] = 1;
#s.t.  a3: x_1[1,6,3,1] = 1;
#s.t.  a4: x_1[1,10,4,1] = 1;
#s.t.  a5: x_1[1,15,5,1] = 1;
#s.t.  a6: x_1[1,20,6,1] = 1;
#s.t.  a7: x_1[1,24,7,1] = 1;
#s.t.  a8: x_1[1,28,8,1] = 1;
		


