#part 1: decision variable
var x1>=0; #first variable
var x2>=0; #second variable

#part 2: objective function
maximize z: 300*x1 + 200*x2

#part 3: constraints
subject M1: 2*x1 + x2 <=8;
subject M2:   x1 + 2*x2 <=8;