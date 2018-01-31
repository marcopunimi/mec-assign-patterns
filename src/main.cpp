//
//  main.cpp
//  ap_cloudlet_association
//
//  Created by Marco Premoli on 17/01/2017.
//  Copyright © 2017 Marco Premoli. All rights reserved.
//


#include "main.hpp"

int main(int argc, const char * argv[])
{
    int i;
    
    string filenameParam, filenameDistance, filenameDemand;
    ModelVariant modelVariant = Binary; //ModelVariant::Binary;
    
    // support variables for column generation total execution time
    time_t initExecutionTime;
    std::time(&initExecutionTime);
    int limitExecutionTime = -1;
    bool betaZero = false;
    double limitDistance = DBL_MAX;
    double epgap = 0.001;
    double beta_param = -1;
    double alpha_param = -1;
    bool single_assignment = false;
    
    double executionTime = 0.0;
    bool printHelp = false;
    ExplorationStrategy sphBPexplorationStrat = DepthFirst;
    BranchingRule sphBPbranchRule = MostFractional;
    
    string initSolDatafile = "";
    
    for (i = 1; i < argc; i++)
    { /* We will iterate over argv[] to get the parameters stored inside.
                                  * Note that we're starting on 1 because we don't need to know the
                                  * path of the program, which is stored in argv[0] */
        // Check that we haven't finished parsing already
        if (std::string(argv[i]) == "-dp") {
            // We know the next argument *should* be the filename:
            if (i + 1 < argc){
                filenameParam = argv[i + 1];
            }
            continue;
        }
        if (std::string(argv[i]) == "-dd") {
            // We know the next argument *should* be the filename:
            if (i + 1 < argc){
                filenameDistance = argv[i + 1];
            }
            continue;
        }
        if (std::string(argv[i]) == "-dt") {
            // We know the next argument *should* be the filename:
            if (i + 1 < argc){
                filenameDemand = argv[i + 1];
            }
            continue;
        }
        else if(std::string(argv[i]) == "-h"){
            printHelp = true;
            continue;
        }
        else if (std::string(argv[i]) == "-c") {
            // the countinous variant of the model will be executed
            modelVariant = Continuous; //ModelVariant::Continuous;
            continue;
        }
        else if (std::string(argv[i]) == "-bh") {
            // the countinous variant of the model will be executed
            modelVariant = BinaryHeuristic;
            continue;
        }
        else if (std::string(argv[i]) == "-sph") {
            // the countinous variant of the model will be executed
            modelVariant = ShortestPathHeuristic;
            continue;
        }
        else if (std::string(argv[i]) == "-sphc") {
            // the countinous variant of the model will be executed
            modelVariant = ShortestPathHeuristicCyclic;
            continue;
        }
        else if (std::string(argv[i]) == "-bc") {
            // the MILP of the compact model with cyclic variant will be executed
            modelVariant = BinaryCyclic;
            continue;
        }
        else if (std::string(argv[i]) == "-cc") {
            // the countinous relaxation of the compact model with cyclic variant will be executed
            modelVariant = ContinuousCyclic;
            continue;
        }
        else if (std::string(argv[i]) == "-sphbp") {
            // We know the next argument *should* be the filename:
            modelVariant = ShortestPathBranchAndPrice;
            for(int j = 1; j < 3; j++){
                if (i + j < argc)
                {
                    string branchPriceArg = std::string(argv[i + j]);
                    if(branchPriceArg == "pq"){
                        sphBPexplorationStrat = PriorityQueue;
                    }
                    else if(branchPriceArg == "df"){
                        sphBPexplorationStrat = DepthFirst;
                    }
                    else if(branchPriceArg== "h"){
                        sphBPbranchRule = HighestFractional;
                    }
                    else if(branchPriceArg== "m"){
                        sphBPbranchRule = MostFractional;
                    }
                    else if(branchPriceArg== "s"){
                        sphBPbranchRule = SplitAssignment;
                    }
                }
            }
            
            continue;
        }
        else if (std::string(argv[i]) == "-gah") {
            // the countinous variant of the model will be executed
            modelVariant = GeneralizedAssignHeuristic;
            continue;
        }
        else if (std::string(argv[i]) == "-epgap") {
            // We know the next argument *should* be the filename:
            if (i + 1 < argc){
                epgap = atof(argv[i + 1]);
            }
            continue;
        }
        else if (std::string(argv[i]) == "-t") {
            // We know the next argument *should* be the filename:
            if (i + 1 < argc){
                limitExecutionTime = atoi(argv[i + 1]);
            }
            continue;
        }
        else if (std::string(argv[i]) == "-bzero") {
            // We know the next argument *should* be the filename:
            if (i + 1 < argc){
                betaZero = true;
                limitDistance = atof(argv[i + 1]);
            }
            else{
                cout << "WARNING: input parameter bzero: no distance limit declared" << endl;
                cerr << "WARNING: input parameter bzero: no distance limit declared" << endl;
            }
            continue;
        }
        else if (std::string(argv[i]) == "-sol") {
            // We know the next argument *should* be the filename:
            if (i + 1 < argc){
                initSolDatafile = argv[i + 1];
            }
            continue;
        }
        else if (std::string(argv[i]) == "-beta") {
            // We know the next argument *should* be the filename:
            if (i + 1 < argc){
                beta_param = atof(argv[i + 1]);
            }
            continue;
        }
        else if (std::string(argv[i]) == "-alpha") {
            // We know the next argument *should* be the filename:
            if (i + 1 < argc){
                alpha_param = atof(argv[i + 1]);
            }
            continue;
        }
        else if (std::string(argv[i]) == "-single") {
            single_assignment = true;
            continue;
        }
    }
    
    if (filenameParam == "" || filenameDistance == "" || filenameDemand == "" || printHelp)
    {
        if((filenameParam == "" || filenameDistance == "" || filenameDemand == "")
           && !printHelp)
        {
            cout << "ERROR: no data filename!" << endl;
            cerr << "ERROR: no data filename!" << endl;
        }
        cerr << "Arguments: " << endl
            << " -dp filepath :\t complete path of the parameter file to read" << endl
            << " -dd filepath :\t complete path of the distance file to read" << endl
            << " -dt filepath :\t complete path of the traffic matrix file to read" << endl
            << " [-epgap double (0,1] ] : \t percentage opt. gap for MIP" << endl
            << " [-c] : \t execute countinous variant of the model" << endl
            << " [-bh] : \t execute binary greedy heuristic" << endl
            << " [-bc] : \t execute MILP of compact model with cyclic version" << endl
            << " [-cc] : \t execute continuous variant of compact model with cyclic version" << endl
            << " [-sph] : \t execute shortest path heuristic" << endl
            << " [-sphc] : \t execute shortest path heuristic with cyclic demand" << endl
            << " [-gah] : \t execute general assignment heuristic" << endl
            << " [-sphbp [exploration-strategy] [branching-rule]]: \t execute shortest path heuristic with branch and price \n \t with exploration strategy: " <<
                "\n \t \t  pq = priority queue smaller LB \n \t\t df = depth first search (default)"
                << " \n \t with branching rule: \n \t\t h highest fractional \n \t\t m most fractional (default) \n \t\t s split selection" << endl
            << " [-t int] : time limit of execution (in seconds)" << endl
            << " [-sol solFilePath : \t complete path to file containing assignment matrix of a feasible solution]" << endl
            << " [-h] : \t help" << endl
            << " default execution is ILP compact formulation variant" << endl;
        
        exit(-1);
    }
    
    InputDataApCloudletAssoc inputData = InputDataApCloudletAssoc(filenameParam,
                                                              filenameDistance,
                                                              filenameDemand,
                                                              limitExecutionTime,
                                                              betaZero,
                                                              limitDistance);
    if(!inputData.readOk)
    {
        cout << "error in reading data" << endl;
        cerr << "error in reading data" << endl;
        exit(-1);
    }
    else{
        // inputData.printAllData();
        inputData.setInitialExecutionTime(initExecutionTime);
        if(initSolDatafile != ""){
            inputData.readIntegerSolutionFromFile(initSolDatafile);
        }
    }
    
    if(alpha_param != -1){
        inputData.alpha = alpha_param;
        cout << "changed alpha to " << beta_param << endl;
    }
    if(beta_param != -1){
        inputData.beta = beta_param;
        cout << "changed beta to " << beta_param << endl;
    }
    
    if(single_assignment){
        inputData.single_assignment = single_assignment;
        cout << "single assignment through time" << beta_param << endl;
        modelVariant = Binary;
    }
    
    if(modelVariant== BinaryHeuristic){
        cout << " BINARY GREEDY HEURISTIC " << endl;
        cerr << " BINARY GREEDY HEURISTIC " << endl;
        inputData.greedy_heuristic();
    }
    else{
        IloEnv env;
        if(modelVariant == ShortestPathHeuristic){
            cout << " SHORTEST PATH COLUMN GENERATION " << endl;
            cerr << " SHORTEST PATH COLUMN GENERATION " << endl;
            ShortPathHeur shortestPathHeuristic(env, &inputData);
            if(shortestPathHeuristic.solve(XvarIndexes(), None, -1)){
                cout << " CG SOLUTION " << endl;
                cout << " final CG cost " << shortestPathHeuristic.getFinalLB() << endl;
                
                // print original variable values
                InputDataApCloudletAssoc::printVarValues(shortestPathHeuristic.getFinalVarValues());
                    
                // print best integer solution cost and variables
                inputData.printBestIntegerSolution();
                inputData.printAssignmentMatrixBestIntegerSolution();
            }
        }
//        else if(modelVariant == ShortestPathHeuristicCyclic)
//        {
//            cout << " SHORTEST PATH COLUMN GENERATION CYCLIC DEMAND " << endl;
//            cerr << " SHORTEST PATH COLUMN GENERATION CYCLIC DEMAND" << endl;
//            ShortPathHeurCycle shortestPathHeuristicCycle(env, &inputData);
//            if(shortestPathHeuristicCycle.solve())
//            {
//                cout << " CG SOLUTION " << endl;
//                cout << " final CG cost " << shortestPathHeuristicCycle.getFinalLB() << endl;
//
//                // print original variable values
//                InputDataApCloudletAssoc::printVarValues(shortestPathHeuristicCycle.getFinalVarValues());
//
//                // print best integer solution cost and variables
//                inputData.printBestIntegerSolution();
//                inputData.printAssignmentMatrixBestIntegerSolution();
//            }
//        }
        else if(modelVariant == Binary ||
                modelVariant == Continuous ||
                modelVariant == BinaryCyclic ||
                modelVariant == ContinuousCyclic)
        {
            if(modelVariant == Binary || modelVariant == BinaryCyclic)
            {
                cout << " BINARY COMPLETE MODEL " << (modelVariant == BinaryCyclic?"CYCLIC":"") << endl;
                cerr << " BINARY COMPLETE MODEL " << (modelVariant == BinaryCyclic?"CYCLIC":"") << endl;
            }
            else{
                cout << " CONTINUOUS COMPLETE MODEL " << (modelVariant == ContinuousCyclic?"CYCLIC":"") << endl;
                cerr << " CONTINUOUS COMPLETE MODEL " << (modelVariant == ContinuousCyclic?"CYCLIC":"")<< endl;
            }
            
            solveCompleteModel( env, &inputData, modelVariant, epgap, limitExecutionTime);
        }
//        else if(modelVariant == GeneralizedAssignHeuristic){
//            cout << " GENERALIZED ASSIGNMENT COLUMN GENERATION " << endl;
//            cerr << " GENERALIZED ASSIGNMENT COLUMN GENERATION " << endl;
//
////            cout << " FIRST SOLVE SHORTEST PATH HEURISTIC " << endl;
////            ShortPathHeur shortestPathHeuristic(env, &inputData);
////            if(shortestPathHeuristic.solve(XvarIndexes(), None, -1))
////            {
////                cout << " CG SOLUTION " << endl;
////                cout << " final CG cost " << shortestPathHeuristic.getFinalLB() << endl;
////
////                // print original variable values
////                InputDataApCloudletAssoc::printVarValues(shortestPathHeuristic.getFinalVarValues());
////
////                // print best integer solution cost and variables
////                inputData.printBestIntegerSolution();
////
////                cout << " try generalized assignment c.g. with best solution yet found " << endl;
////
////                GenerAssignHeur generalAssignHeuristic(env,&inputData);
////                generalAssignHeuristic.solve();
////            }
//
//            GenerAssignHeur generalAssignHeuristic(env,&inputData);
//            generalAssignHeuristic.solve();
//        }
        else if(modelVariant == ShortestPathBranchAndPrice)
        {
            std::stringstream modelStringToPrint;
            modelStringToPrint.str("");
            modelStringToPrint << " - SHORTEST PATH BRANCH-AND-PRICE - " << (sphBPexplorationStrat==PriorityQueue?"PQ":"DFS");
            if(sphBPbranchRule==HighestFractional){
                modelStringToPrint << " HIGHEST FRACTIONAL ";
            }
            else if(sphBPbranchRule == MostFractional){
                modelStringToPrint << " MOST FRACTIONAL ";
            }
            else if(sphBPbranchRule == SplitAssignment){
                modelStringToPrint << " SPLIT POSSIBLE ASSIGNMENT AP-TIME";
            }
            
            cout << modelStringToPrint.str() << endl;
            cerr << modelStringToPrint.str() << endl;
            
            ShortPathHeurBP shortestPathBP(env, &inputData, sphBPbranchRule);
            
            if(shortestPathBP.solve(sphBPexplorationStrat))
            {
                cout << " B-and-P ended " << endl;
                // print best integer solution cost and variables
                inputData.printBestIntegerSolution();
            }
            else{
                cout << " shortest path heuristic failed " << endl;
            }
        }
        env.end();
    }
    executionTime = computeTime(initExecutionTime);
    cout << "end of execution in " << executionTime << " seconds" << endl;
    cerr << "end of execution in " << executionTime << " seconds" << endl;
    
    return 0;
}

void solveCompleteModel(IloEnv env, InputDataApCloudletAssoc* inputData, ModelVariant modelVariant,
                         double epgap, int timelimit)
{
    
    bool cyclicVariant = (modelVariant == BinaryCyclic || modelVariant == ContinuousCyclic);
    Ap_cloudlet_assoc_model problem = Ap_cloudlet_assoc_model(env, inputData, modelVariant, cyclicVariant);
    
    IloCplex cplex(problem.model);
    
    if(modelVariant == Binary || modelVariant == BinaryCyclic)
    {
        cplex.setParam(IloCplex::EpGap, epgap);
    }
    else{
        cplex.setParam(IloCplex::Param::RootAlgorithm, CPX_ALG_PRIMAL);
    }
    
    // solve at root
    // cplex.setParam(IloCplex::NodeLim, 0);
    
    // set time limit to cplex execution if required
    if(timelimit !=-1)
    {
        cplex.setParam(IloCplex::TiLim, (double) timelimit);
    }
    
    cplex.exportModel("/Users/marcopremoli/Desktop/ti-bdc13/test/ap_cloudlet_model.lp");
    
    problem.solve(cplex);
    
    cplex.end();
}

