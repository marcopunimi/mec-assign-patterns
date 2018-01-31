
------------
SOURCE CODE STRUCTURE:
------------

- main.cpp (.hpp) = the main file of the algorithm
- InputDataApCloudletAssoc.cpp (.hpp) = it's the class containing the input and output data for the problem, togheter with some heuristics for the problem
- ap_cloudlet_assoc_model.cpp (.hpp) = contains the compact model of the problem and resolution method using the ILP general purpose solver of CPLEX
- ShortPathHeur.cpp (.hpp) = contains the code for the column generation algorithm
- ShortPathHeurBP.cpp (.hpp) = contains the code for the branch and price. it makes use of the ShortPathHeur code as CG
- supportFunctions.cpp (.hpp) = contains some support function to use throughout the code

------------
BUILDING
------------

To build the code, IBM ILOG CPLEX must be installed on the machine.
An example of the Makefile to build the code can be found in Makefile_ex file.

------------
COMMAND LINE PARAMETERS
------------

Command line parameters: 
* -dp filepath :	 complete path of the parameter file to read [mandatory]
* -dd filepath :	 complete path of the distance file to read [mandatory]
* -dt filepath :	 complete path of the traffic matrix file to read [mandatory]
* [-epgap double (0,1] ] : 	 percentage opt. gap for MIP
* [-c] : 	 execute continuous variant of the model
* [-sph] : 	 execute shortest path heuristic
* [-sphbp [exploration-strategy] [branching-rule]]: 	 execute shortest path heuristic with branch and price 
 	 - with exploration strategy: 
 	  -	pq = priority queue smaller LB 
 		 - df = depth first search (default) 
 	 - with branching rule: 
 		 - h highest fractional 
 		 - m most fractional (default) 
 		 - s split selection (recommended)
* [-t int] : time limit of execution (in seconds)
* [-h] : 	 help
 default execution is ILP compact formulation variant.

Executable will print both on standard output and on standard error (the latter for log purpose only)

 Example of command line uses:
* this command executes the ILP general solver of CPLEX over the compact model, with no time limit or percentage optimality gap limit
 - ./executable -dp /path/to/paramfile -dd /path/to/distancefile -dt /path/to/demandfile > output.out 2> output.log   
* this command executes the ILP general solver of CPLEX over the compact model, with one hour time limit (3600 seconds) and a limit of 1% of optimality gap
 - ./executable -dp /path/to/paramfile -dd /path/to/distancefile -dt /path/to/demandfile -t 3600 -epgap 0.01 > output.out 2> output.log  
* this command executes the CG algorithm over the model, stopping at the end of the CG (a time limit can be set with -t option)
 - ./executable -dp /path/to/paramfile -dd /path/to/distancefile -dt /path/to/demandfile -sph > output.out 2> output.log
* this command executes the branch-and-price algorithm over the model with the split selection branching rule (recommended), stopping when the tree search is empty (a time limit can be set with -t option). 
 - ./executable -dp /path/to/paramfile -dd /path/to/distancefile -dt /path/to/demandfile -spbh s > output.out 2> output.log

------------
INPUT FILES FORMAT
------------
* PARAMETER FILE (option -dp in command line)
     - 1 double = alpha parameter
     - 1 double = beta parameter
     - 1 double = U parameter
     - 1 double > 1 = facility percentage extra capacity w.r.t. peak demand in time
    
Example:
0.5   
0.5   
1
1.05
     
* DISTANCE MATRIX FILE (option -dd in command line)
     - 1 int = AP cardinality
     - 1 int = facility cardinality
     - matrix A x K double = distance AP to MEC facility (comma separated values)
     - matrix K x K double = distance facility to facility (comma separated values)

Example:
3
2
1, 1, 1
2, 2, 1
3, 4, 7
0, 2 
2, 0

* DEMAND FILE (option file -dt in command line)
     - 1 int = time slot cardinality
     - matrix N x T double = demand nodes in time (comma separeted values)

Example:
4
0, 5, 21, 5
1, 4, 66, 2
1, 44, 5, 7
