//
//  ShortPathHeur.cpp
//  ap_cloudlet_association
//
//  Created by Marco Premoli on 31/01/2017.
//  Copyright © 2017 Marco Premoli, Università Degli Studi di Milano. All rights reserved.
//
/*
 This file is part of mec-assign-patterns.
 
 mec-assign-patterns is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 mec-assign-patterns is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with mec-assign-patterns.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "ShortPathHeur.hpp"

ShortPathHeur::~ShortPathHeur()
{
    // destroy branching tree in heap
    destroyBranchingTree(branchingRoot);
}

// update dual variables from current cplex solution
void ShortPathHeur::updateDualValues(IloCplex* cplex, DualVariables & dualToUpdate)
{
    int i,j;
    
    for(i = 0; i < this->cloudletCapacity.getSize(); i++)
    {
        for(j =0 ; j < this->cloudletCapacity[i].getSize();j++)
        {
            dualToUpdate.lambda[i][j] = cplex->getDual(cloudletCapacity[i][j]);
        }
    }
    
    for(i = 0; i < this->valorizeZ.getSize();i++)
    {
        dualToUpdate.eta[i] = cplex->getDual(valorizeZ[i]);
    }
}

// add column variable to the model given an Association Path object, containing the AP and the assignment paths for each time-slot
void ShortPathHeur::addColumn(AssociationPath path)
{
    IloNumColumn newColumn(this->model.getEnv());
    
    // add variable to objective function
    // sum{t in T} sum{k in K} (beta*d[i,t]*m[i,k]*x[i,k,t] + sum{j in K} alpha*d[i,t+1]*l[k,j]*y[i,k,j,t+1] )
    double objFunContribution = 0.0;
    // for each hop of the assignment path (i.e. for each time-slot)
    for(int t =0; t < path.path.size(); t++)
    {
        // add contribution of the current cloudlet in the time slot
        int time = t;
        int cloudlet = path.path[t];
        
        objFunContribution += inputData->beta*inputData->nodeDemand[path.nodeId][time]*inputData->distanceApCloudlet[path.nodeId][cloudlet];
        
        if (t < path.path.size()-1)
        {
            // add contribution for the change of cloudlet in successive time slots
            int nextTime = t+1;
            int nextCloudlet = path.path[t+1];
            if(cloudlet != nextCloudlet)
            {
                objFunContribution += inputData->alpha*inputData->nodeDemand[path.nodeId][nextTime]*inputData->distanceCloudletCloudlet[cloudlet][nextCloudlet];
            }
        }
    }
    newColumn += this->objFun(objFunContribution);
    
    // add variable to cloudlet capacity constraints, for each time slot a single component for a single cloudlet
    for(int t =0; t < path.path.size(); t++)
    {
        int time = t;
        int cloudlet = path.path[t];
        newColumn += cloudletCapacity[time][cloudlet](-inputData->nodeDemand[path.nodeId][time]);
    }
    
    // add variable to z variable valorization constraint, one single component
    newColumn += valorizeZ[path.nodeId](1.0);
    
    // add column to z variables structure
    char nameVar[100];
    snprintf(nameVar, 100, "z_%d_%ld", path.nodeId, z[path.nodeId].getSize());
    z[path.nodeId].add(IloNumVar(newColumn, 0.0, 1.0, ILOFLOAT ,nameVar));
    // corresponding z value
    zval[path.nodeId].push_back(0.0);
}

// code to fix a variable in the model
void ShortPathHeur::fixVariable(XvarIndexes xvarIdx, FixingDirection fixDirection)
{
    if(fixDirection != None)
    {
        int i = xvarIdx.ap;
        int k = xvarIdx.cloudlet;
        int t = xvarIdx.time;
        
        if(fixDirection == AtZero)
        {
            // set to 0 all z variables of AP i using cloudlet k at time t
            for(int p = 0; p < addedPaths[i].size(); p++)
            {
                AssociationPath path = addedPaths[i][p];
                if(path.path[t] == k)
                {
                    z[i][p].setUB(0.0);
                }
            }
            // store the new forbidden assignment
            this->fixedVar_small[i][t][k] = 0;
        }
        else if(fixDirection == AtOne)
        {
            // set to 0 all z variables of AP i NOT USING cloudlet k at time t
            for(int p = 0; p < addedPaths[i].size(); p++)
            {
                AssociationPath path = addedPaths[i][p];
                if(path.path[t] != k)
                {
                    z[i][p].setUB(0.0);
                }
            }
            // store the new forced assignment
            this->fixedVar_small[i][t][k] = 1;
        }
    }
}

// given a new assignment path to add to the model, check if it already exist
bool ShortPathHeur::checkIfPathAlreadyExist(int ap, vector<int> path)
{
    // check all paths of the current AP, (try backwards...)
    for(int p = ((int) this->addedPaths[ap].size()-1); p >= 0 ;p--)
    {
        bool pathEqual = true;
        // check all hops of the current path
        for(int t = 0; t < this->addedPaths[ap][p].path.size();t++)
        {
            if(addedPaths[ap][p].path[t] != path[t])
            {
                // one of the hop does not coincide -> skip the check for the subsequent hops
                pathEqual = false;
                break;
            }
        }
        // an equal path has been found, return TRUE
        if(pathEqual)
            return true;
    }
    // all paths are not equal to the new path, return FALSE
    return false;
}

bool ShortPathHeur::solve()
{
    // solve without fixing variables --> root CG of B-&-P
    return solve(XvarIndexes(), None, -1.0);
}

// solve CG fixing a variable in a direction (i.e. to 0 or to 1) given a current best UB
bool ShortPathHeur::solve(XvarIndexes xvarIdx, FixingDirection fixDirection, double actualBestUB)
{
    IloCplex cplex(this->model);
    cplex.setParam(IloCplex::Param::RootAlgorithm, CPX_ALG_PRIMAL);
    bool result = solve(xvarIdx,fixDirection,actualBestUB, &cplex, fixDirection!=None);
    cplex.end();
    return result;
}

bool ShortPathHeur::solve(double actualBestUB, IloCplex*cplex){
    return solve(XvarIndexes(), None, actualBestUB, cplex, true);
}

// solve CG fixing a variable in a direction (i.e. to 0 or to 1) given a current best UB and a given IloCplex instance
bool ShortPathHeur::solve(XvarIndexes xvarIdx, FixingDirection fixDirection, double actualBestUB, IloCplex* cplex, bool useDummy)
{
    if(useDummy || fixDirection!=None)
    {
        // activate dummy variable when fixing variables
        if(!isDummyUsed)
        {
            dummy.setUB(1.0);
            isDummyUsed = true;
        }
        if(fixDirection!=None)
            fixVariable(xvarIdx, fixDirection);
    }
    
    // initialized cplex for the master problem
    //IloCplex cplex(model);
    
    //bool cgSuccess = true;
    isCgOK = true;
    
    double primalBound = 0.0, previousPrimalBound = 1e12, dualBound = 0.0,
        highestDualBound = -DBL_MAX, lowest_gap = DBL_MAX, lowest_gap_perc = 1.0,
        integerSolCost = inputData->getBestIntegerCost(), integerPrimalGapPerc = 1;
    double epsilon = 1e-4;
    double iterationGap = DBL_MAX;
    int cgNumberIteration = 0, i;
    
    // support variables for column generation total execution time
    time_t cgInitExecutionTime;
    std::time(&cgInitExecutionTime);
    double cgExecutionTime = 0;
    
    // support variable for log purpose
    std::stringstream cgIterationLogString;
    
    do
    {
        // support variables for column generation iteration execution time (and pricing execution time)
        time_t cgIterationInitTime;
        std::time(&cgIterationInitTime);
        double cgIterationExecTime = 0;
        double pricingIterationExecutionTime = 0;
        double masterTime = 0;
        
        // solve master
        bool isCplexSolved = false;
        isCplexSolved = cplex->solve();
        masterTime = computeTime(cgIterationInitTime);
        
        if(isCplexSolved)
        {
            if(cplex->getStatus() != IloAlgorithm::Infeasible && cplex->getCplexStatus() != IloCplex::Infeasible)
            {
                
                //if(!isDummyUsed)
                // {
                //     char lpFileName[300] = "";
                //     snprintf(lpFileName, 300, "/Users/marcopremoli/Desktop/ap_cloud_assoc_sph_%d.lp", cgNumberIteration );
                //     cplex->exportModel(lpFileName);
                // }
                
                // get primal bound
                primalBound = dualBound = cplex->getObjValue();
                cout << "cplex ok with o.f. value " << primalBound << endl;
                
                // save z variable values
                bool isIntegerSol = true;
                for(i = 0; i < z.getSize(); i++)
                {
                    for(int j = 0; j < z[i].getSize();j++)
                    {
                        double value = cplex->getValue(z[i][j]);
                        zval[i][j] = value;
                        if(IS_FRACTIONAL(value)) //  value > 1e-12 && value < 1-1e-12)
                        {
                            // check if solution is already integer
                            isIntegerSol = false;
                        }
                    }
                }
                this->isIntegerSolution = isIntegerSol;
                
                if(isDummyUsed)
                {
                    dummyVal = cplex->getValue(dummy);
                    cout << "dummy var value " << dummyVal << endl;
                    if(dummyVal > 1e-12)
                        this->isIntegerSolution = false;
                }
                
                if(isIntegerSolution){
                    cout << "CG INTEGER SOLUTION" << endl;
                }
                
                // get compact formulation variable values
                originalVarValues = this->getOriginalVarValuesFromZ();
                
                // update duals
                updateDualValues(cplex, this->dualVariables);
                
                // get new paths from the pricer with the current duals
                int countAddedPathIteration = 0;
                
                // if solution is already integer, no new paths have to be added
                time_t pricingInitTime;
                std::time(&pricingInitTime);
                vector<AssociationPath> newPaths = this->getNewPathsPricing(dualVariables);
                pricingIterationExecutionTime = computeTime(pricingInitTime);
                
                // add new columns and update dual bounds
                for(i = 0; i < newPaths.size();i++)
                {
                    if(newPaths[i].cost < -epsilon)
                    {
                        this->addColumn(newPaths[i]);
                        addedPaths[newPaths[i].nodeId].push_back(newPaths[i]);
                        
                        dualBound += newPaths[i].cost;
                        
                        cout << " path cost " << newPaths[i].cost << endl;
                        newPaths[i].printPath();
                        
                        countAddedPathIteration++;
                    }
                }
                
                // update iteration gap and overall higher dual and lower gap if needed
                bool lowestGapChanged = false;
                
                // primal dual gap in this iteration
                iterationGap = primalBound - dualBound;
                
                if(cgNumberIteration == 0)
                {
                    previousPrimalBound = primalBound;
                    highestDualBound = dualBound;
                    lowest_gap = iterationGap;
                }
                
                if(primalBound < previousPrimalBound)
                {
                    previousPrimalBound = primalBound;
                    lowestGapChanged = true;
                }
                
                if(dualBound > highestDualBound)
                {
                    highestDualBound = dualBound;
                    lowestGapChanged = true;
                }
                if((previousPrimalBound - highestDualBound) < lowest_gap)
                {
                    lowest_gap = previousPrimalBound - highestDualBound;
                    
                    if(lowest_gap !=0) lowest_gap_perc = (lowest_gap*100) / primalBound;
                    else lowest_gap_perc = 0;
                    
                    lowestGapChanged = true;
                }
                
                // try get integer solution
                if(cgNumberIteration > -1)
                {
                    if(!isDummyUsed)
                    {
                        computeIntegerSolution(integerSolCost, primalBound, false);
                    }
                }
                if(integerSolCost != -1)
                {
                    integerPrimalGapPerc = (integerSolCost - primalBound)/primalBound;
                }
                // get execution time for cg iteration
                cgIterationExecTime = computeTime(cgIterationInitTime);
                cgExecutionTime = computeTime(cgInitExecutionTime);
                // pricingIterationExecutionTime = cgIterationExecTime - masterTime;
                
                // print some info in final output file and in log file
                cgIterationLogString.str("");
                cgIterationLogString << "CG ITER "<< cgNumberIteration
                << " tExecTime " << cgExecutionTime
                << " #p: " << countAddedPathIteration
                << " iExecTime: " << cgIterationExecTime
                << " mTime: " << masterTime
                << " pricExecTime " << pricingIterationExecutionTime
                << " PB: " << primalBound
                << " hDB: " << highestDualBound
                << " IB: " << integerSolCost
                << " IPGAP: " << integerPrimalGapPerc
                << " lGAP: " << lowest_gap
                << " plGP " << lowest_gap_perc << " % "
                << " iDB: " << dualBound
                << " iGAP: " << iterationGap;
                
                if(inputData->getIntegralityGAP() > -1){
                    cgIterationLogString << " *GAP " << inputData->getMinActiveBranchNodeLB() << " *optGAP " << inputData->getIntegralityGAP();
                }
                
                cout << cgIterationLogString.str() << endl;
                cerr << cgIterationLogString.str() << endl;
            }
            else
            {
                // cplex execution had some problem, break the loop
                cout << "cplex status " << cplex->getStatus() << " " << cplex->getCplexStatus() << endl;
                cerr << "cplex status " << cplex->getStatus() << " " << cplex->getCplexStatus() << endl;
                isCgOK = false;
                break;
            }
        }
        else{
            // cplex execution had some problem, break the loop
            cout << "cplex status " << cplex->getStatus() << " " << cplex->getCplexStatus() << endl;
            cerr << "cplex status " << cplex->getStatus() << " " << cplex->getCplexStatus() << endl;
            isCgOK = false;
            break;
        }
        cgNumberIteration++;
        
        if(actualBestUB > -1 && highestDualBound != 0)
        {
            double gap_highLB_UB = (actualBestUB - highestDualBound)/highestDualBound;
            if(gap_highLB_UB <= 0.0005)
            {
                cout << "actual DB " << highestDualBound << " more or equal to UB " << actualBestUB << " gap: " << gap_highLB_UB << endl;
                isCgOK = false;
                break;
            }
        }
    } while(lowest_gap >= epsilon && !inputData->isExcTimeLimitReached());
    
    
    if(isCgOK)
    {
        if(!isDummyUsed || (isDummyUsed && dummyVal < 1e-10))
        {
            finalLB = primalBound;
            
            // if dummy is used, the heuristic rounding is performed only at the end of the successfull CG
            // otherwise is execute at every CG iteration, hence at CG end is not needed
            if(isDummyUsed)
            {
                if(computeIntegerSolution(integerSolCost, finalLB, true))
                {
                    inputData->printBestIntegerSolution();
                    
                    cout << " new integer solution while branching : try lagrangian fixing  on root" << endl;
                    cerr << " new integer solution while branching : try lagrangian fixing  on root" << endl;
                    tryLagrangianProbingOnNode(branchingRoot, cplex);
                }
            }
        }
        else
        {
            // CG ended successfully but dummy variable is still used -> hence the master problem is infeasible
            isCgOK = false;
        }
    }
    
    return(isCgOK);
}

// return the final lower bound of the CG, or -1 if CG was not successfull
double ShortPathHeur::getFinalLB()
{
    if(isCgOK)
    {
        return finalLB;
    }
    return -1;
}

// return final variables values of CG, or variables at value -1 if CG was not successfull
VarValues ShortPathHeur::getFinalVarValues()
{
    if(isCgOK)
    {
        return originalVarValues;
    }
    
    int nA = inputData->apCard;
    int nT = inputData->timeSlotCard;
    int nK = inputData->cloudletCard;
    
    return VarValues(nA,nK,nT);
}

// given current z variable values compute original x and y variables values, or all variables at -1 if CG was not successfull
VarValues ShortPathHeur::getOriginalVarValuesFromZ()
{

    int nA = inputData->apCard;
    int nT = inputData->timeSlotCard;
    int nK = inputData->cloudletCard;
    int i,j,k,t;
    
    VarValues varValues(nA,nK,nT);
    
    if(this->isCgOK)
    {
        // give right values to original x and y variables to print final solution
        // for each AP
        for(i = 0; i < zval.size(); i++)
        {
            // for each path related to the AP
            for(j = 0; j < zval[i].size(); j++)
            {
                // check if the path is used
                if(zval[i][j] > 1e-6)
                {
                    // get corresponding assignment path
                    AssociationPath assocPath = this->addedPaths[i][j];
                    // for each hop in the path
                    for(t = 0; t < assocPath.path.size();t++)
                    {
                        int time = t;
                        int cloudlet = assocPath.path[t];
                        // x variable
                        varValues.xval[i][cloudlet][time] += zval[i][j];
                        // y variable of the hop (if exist)
                        if(t < assocPath.path.size()-1)
                        {
                            int nextTime = t + 1;
                            int nextCloud = assocPath.path[t+1];
                            varValues.yval[i][cloudlet][nextCloud][nextTime] += zval[i][j];
                        }
                    }
                }
            }
        }
        
        // check feasibility of solution
        for(i = 0; i < nA; i++)
        {
            for(t = 0; t < nT; t++)
            {
                double assing = this->dummyVal;
                for(k = 0; k < nK; k++)
                {
                    assing += varValues.xval[i][k][t];
                }
                if(assing <= 1-1e-6)
                {
                    cout << " AP " << i << " not assigned in time " << t << " = " << assing << endl;
                    exit(-1);
                }
            }
        }
    }
    
    return varValues;
}

// given the current dual variables, check if some variable
// check if it's possible to fix further assignments that lead to a LB higher than current best UB
void ShortPathHeur::_lagrangianProbing(DualVariables duals,
                                     vector<int> & xToFix,
                                      vector<bool> & xAreToZero,
                                     vector<Transition> & forbiddenTransition,
                                     double currPB, double bestUB)
{
    if(bestUB <= 0)
    {
        bestUB = inputData->getBestIntegerCost();
        if(bestUB <= 0)
            return;
    }
    
    int i,t,k1,k2;
    
    // for each AP a shortest path problem has to be solved on a graph whose arc weigths are given by dual variables
    for(i = 0; i < inputData->apCard; i++)
    {
        // get for each node TIME-CLOUDLET the corresponding shortest path cost label
        vector<vector<double>> labelForward =
            vector<vector<double>>(inputData->timeSlotCard, vector<double>(inputData->cloudletCard, 0.0));
        
        vector<vector<double>> labelBackward =
            vector<vector<double>>(inputData->timeSlotCard, vector<double>(inputData->cloudletCard, 0.0));
        
        // to remove transition, for each node TIME-CLOUDLET keep the mininum label cost before the cost of the node itself
        vector<vector<double>> labelForwardBeforeNode =
            vector<vector<double>>(inputData->timeSlotCard, vector<double>(inputData->cloudletCard, 0.0));
        
        vector<vector<int>> pred = vector<vector<int>>(inputData->timeSlotCard, vector<int>(inputData->cloudletCard, -1));;
        
        // for every time-slot greater than 0, check the precedent time-slot minimum cost
        for(t = 0; t < inputData->timeSlotCard; t++)
        {
            int tback = inputData->timeSlotCard - 1 - t;
            for(k1 = 0; k1 < inputData->cloudletCard; k1++)
            {
                // cost of the node
                labelForward[t][k1] = inputData->nodeDemand[i][t]*(inputData->beta*inputData->distanceApCloudlet[i][k1] + duals.lambda[t][k1]);
                labelBackward[tback][k1] = inputData->nodeDemand[i][tback]*(inputData->beta*inputData->distanceApCloudlet[i][k1] + duals.lambda[tback][k1]);
                
                if(t > 0)
                {
                    // get as previous node the node minimizing the current cost
                    // find min{j in K}(costPrev[j] + alpha*d[i][t]*l[k][j])
                    double minCostForward = -1;
                    double minCostBackward = -1;
                    int predTemp = -1;

                    bool forwardFound = false;
                    bool backwardFound = false;
                    for(k2 = 0; k2 < inputData->cloudletCard; k2++)
                    {
                        if(!forwardFound)
                        {
                            // check if in the previous time-slot the triplet AP-(TIME-1)-CLOUDLET was forbidden
                            // if FORBIDDEN -> skip the check
                            if(this->fixedVar_small[i][t-1][k2] != 0)
                            {
                                double currCost = labelForward[t-1][k2] + inputData->alpha*inputData->nodeDemand[i][t]*inputData->distanceCloudletCloudlet[k2][k1];
                                
                                // check if in the previous time-slot a triplet was FORCED
                                // if FORCED -> skip the check of all other time-slot-cloudlet and assign directly
                                if(this->fixedVar_small[i][t-1][k2] == 1)
                                {
                                    if(isTransitionPossible[i][k2][k1][t]){
                                        minCostForward = currCost;
                                        predTemp = k2;
                                    }
                                    else{
                                        minCostForward = -1;
                                        predTemp = -1;
                                    }
                                    forwardFound = true;
                                }
                                else if(isTransitionPossible[i][k2][k1][t] &&
                                        (currCost < minCostForward || minCostForward == -1))
                                {
                                    minCostForward = currCost;
                                    predTemp = k2;
                                }
                            }
                        }
                        
                        if(!backwardFound)
                        {
                            // check if in the previous time-slot the triplet AP-(TIME-1)-CLOUDLET was forbidden
                            // if FORBIDDEN -> skip the check
                            if(this->fixedVar_small[i][tback+1][k2] != 0)
                            {
                                double currCost = labelBackward[tback+1][k2] + inputData->alpha*inputData->nodeDemand[i][tback+1]*inputData->distanceCloudletCloudlet[k1][k2];
                                
                                // check if in the previous time-slot a triplet was FORCED
                                // if FORCED -> skip the check of all other time-slot-cloudlet and assign directly
                                if(this->fixedVar_small[i][tback + 1][k2] == 1)
                                {
                                    if(isTransitionPossible[i][k1][k2][tback+1])
                                        minCostBackward = currCost;
                                    else
                                        minCostBackward = -1;
                                    
                                    backwardFound = true;
                                }
                                else if(isTransitionPossible[i][k1][k2][tback+1] &&
                                        (currCost < minCostBackward || minCostBackward == -1))
                                {
                                    minCostBackward = currCost;
                                }
                            }
                        }
                    }
                    // update label
                    if(minCostForward >= 0){
                        labelForward[t][k1] += minCostForward;
                        labelForwardBeforeNode[t][k1] = minCostForward;
                        pred[t][k1] = predTemp;
                    }
                    else{
                        labelForward[t][k1] = DBL_MAX;
                        labelForwardBeforeNode[t][k1] = DBL_MAX;
                        pred[t][k1] = -1;
                    }
                    
                    if(minCostBackward>=0)
                        labelBackward[tback][k1] += minCostBackward;
                    else
                        labelBackward[tback][k1] = DBL_MAX;
                }
            }
        }
        // get the minimum cost path among all, found in the last time-slot label forward (or first time-slot label backward)
        double minRC = DBL_MAX;
        // store the nodes in the path
        list<int> minPath;
        int lastCloud = -1;
        for(k1 = 0; k1 < inputData->cloudletCard; k1++)
        {
            if(fixedVar_small[i][inputData->timeSlotCard-1][k1] == 1){
                minRC = labelForward[inputData->timeSlotCard-1][k1];
                lastCloud = k1;
                break;
            }
            else if(fixedVar_small[i][inputData->timeSlotCard-1][k1] == 2 &&
                    (labelForward[inputData->timeSlotCard-1][k1] < minRC || lastCloud == -1))
            {
                    minRC = labelForward[inputData->timeSlotCard-1][k1];
                    lastCloud = k1;
            }
        }
        // save minimum reduced cost path and the path itself
        minRC -= duals.eta[i];
        // to store the path follows the predecessor backward
        minPath.push_back(lastCloud);
        for(t = inputData->timeSlotCard - 1; t > 0; t--)
        {
            int currNode = minPath.front();
            minPath.insert(minPath.begin(), pred[t][currNode]);
        }
        
        // FIND FORBIDDEN ASSIGNMENT TIME-CLOUDLET
        // subtract constant eta dual to the reduced cost
        // for each pair TIME-CLOUDLET that is not yet fixed, check if, forcing the assignment
        // the sum of LB + reduced cost of corresponding path is higher than current UB
        // if true --> fix the assignment to 0
        list<int>::iterator it = minPath.begin();
        for(t = 0; t < inputData->timeSlotCard; t++)
        {
            for(k1 = 0; k1 < inputData->cloudletCard; k1++)
            {
                if(fixedVar_small[i][t][k1] == 2)
                {
                    // FIND FORBIDDEN TRANSITIONS IN SUCCESSIVE TIMES
                    if(t > 0)
                    {
                        for(k2 = 0; k2 < inputData->cloudletCard; k2++)
                            if(fixedVar_small[i][t-1][k2] != 0 && isTransitionPossible[i][k2][k1][t])
                            {
                                double reduced_cost_transition = labelForwardBeforeNode[t][k1] - labelForward[t-1][k2];
                                // check if including the current transition in the current primal bound give a solution worst than current UB
                                if(currPB + reduced_cost_transition > bestUB)
                                {
                                    // forbid transition
                                    forbiddenTransition.push_back(Transition(i, k2, k1, t));
                                }
                            }
                    }
                    
                    double pathCost = -1;
                    
                    // check if current node belong to the minimum cost path
                    bool isMinNode = (*it) != k1;
                    
                    if(isMinNode)
                    {
                        // FIND FORBIDDEN ASSIGNMENT TIME-CLOUDLET
                        if(t == 0)
                        {
                            pathCost = labelBackward[0][k1];
                        }
                        else if(t == inputData->timeSlotCard -1)
                        {
                            pathCost = labelForward[t][k1];
                        }
                        else
                        {
                            pathCost = labelForwardBeforeNode[t][k1] + labelBackward[t][k1];
                        }
                    }
                    else
                    {
                        // FIND SECOND BEST ASSIGNMENT TIME-CLOUDLET BY REDUCED COST
                        for(k2 = 0; k2 < inputData->cloudletCard; k2++) if(k2!=k1 && fixedVar_small[i][t][k2] == 2)
                        {
                            if(t == 0)
                            {
                                if(labelBackward[0][k2] < pathCost || pathCost == -1)
                                    pathCost = labelBackward[0][k2];
                            }
                            else if(t == inputData->timeSlotCard - 1)
                            {
                                if(labelForward[t][k2] < pathCost || pathCost == -1)
                                    pathCost = labelForward[t][k2];
                            }
                            else
                            {
                                double minTemp = labelBackward[t][k2] + labelForwardBeforeNode[t][k2];
                                if(minTemp < pathCost || pathCost == -1){
                                    pathCost = minTemp;
                                }
                            }
                        }
                    }
                    // get actual path cost subtracting the constant dual eta
                    pathCost -= duals.eta[i];
                    // get gap between lowest cost path among all and minimum cost path using current node
                    double deltaRC = pathCost - minRC;
                    // if using the current node leads to a solution worst than current Upper Bound, the current node does not have to be used
                    if(deltaRC/2.0 + currPB >= bestUB)
                    {
                        XvarIndexes xtemp = XvarIndexes(i, k1, t);
                        // cout << "new assignment fixed ";
                        // xtemp.print();
                        // cout << endl;
                        xToFix.push_back(_getXvarIndexInVector(xtemp));
                        xAreToZero.push_back(isMinNode);
                    }
                }
            }
            // update minimum cost path iterator to the next time-slot
            ++it;
        }
    }
}

int ShortPathHeur::_getXvarIndexInVector(XvarIndexes xvarToCheck)
{
    return xIndexUnlistMapping[xvarToCheck.ap][xvarToCheck.cloudlet][xvarToCheck.time];
}

// pricing
vector<AssociationPath> ShortPathHeur::getNewPathsPricing(DualVariables duals)
{
    vector<AssociationPath> assocPaths;
    
    int i,t,k1,k2;
    
    // I have two structures: for the current and the successive time-slot
    // storing for every cloudlet the path reaching the cloudlet and its cost
    
    vector<double> costCurrent(inputData->cloudletCard,0.0);
    vector<double> costSucc(inputData->cloudletCard,0.0);
    
    // for each AP a shortest path problem has to be solved on a graph whose arc weigths are given by dual variables
    for(i = 0; i < inputData->apCard; i++)
    {
        vector<vector<int>> pathCurrent(inputData->cloudletCard,vector<int>());
        vector<vector<int>> pathSucc(inputData->cloudletCard,vector<int>());
        
        // initialize cloudlet cost for the first time slice
        for(k1=0; k1 < inputData->cloudletCard; k1++) if(inputData->isApCloudletPossible[i][k1])
        {
            // source arcs weights = d[i][0]*(beta*m[i][k] + lambda[0][k])
            costCurrent[k1] = inputData->nodeDemand[i][0]*(inputData->beta*inputData->distanceApCloudlet[i][k1] + duals.lambda[0][k1]);
            // first hop of the path is the cloudlet itself
            pathCurrent[k1].push_back(k1);
        }
        
        // for every time-slot greater than 0, check the precedent time-slot minimum cost
        for(t = 1; t < inputData->timeSlotCard;t++)
        {
            for(k1 = 0; k1 < inputData->cloudletCard; k1++) if(inputData->isApCloudletPossible[i][k1])
            {
                // cost of the node
                costSucc[k1] = inputData->nodeDemand[i][t]*(inputData->beta*inputData->distanceApCloudlet[i][k1] + duals.lambda[t][k1]);
                
                // as previous node get the one minimizing the current cost
                // find min{j in K}(costPrev[j] + alpha*d[i][t]*l[k][j])
                
                double minCostCurrent = -1;
                int argminCostCurrent = -1;
                
                // find min{j in K}(costPrev[j] + alpha*d[i][t]*l[k][j])
                for(k2 = 0; k2 < inputData->cloudletCard; k2++) if(inputData->isApCloudletPossible[i][k2])
                {
                    // check if in the previous time-slot the triplet AP-(TIME-1)-CLOUDLET was forbidden
                    // if FORBIDDEN -> skip the check
                    if(this->fixedVar_small[i][t-1][k2] != 0)
                    {
                        double currCost = costCurrent[k2] + inputData->alpha*inputData->nodeDemand[i][t]*inputData->distanceCloudletCloudlet[k2][k1];
                        
                        // check if in the previous time-slot a triplet was FORCED
                        // if FORCED -> skip the check of all other time-slot-cloudlet and assign directly
                        if(this->fixedVar_small[i][t-1][k2] == 1)
                        {
                            // check if transition is possible
                            if(isTransitionPossible[i][k2][k1][t])
                            {
                                minCostCurrent = currCost;
                                argminCostCurrent = k2;
                            }
                            else{
                                // if not possible, no path passing on the current node is present
                                minCostCurrent = -1;
                                argminCostCurrent = -1;
                            }
                            break;
                        }
                        else if(isTransitionPossible[i][k2][k1][t] && (currCost < minCostCurrent || minCostCurrent == -1))
                        {
                            minCostCurrent = currCost;
                            argminCostCurrent = k2;
                        }
                    }
                }
                
                if(argminCostCurrent >= 0)
                {
                    costSucc[k1] += minCostCurrent;
                    pathSucc[k1] = pathCurrent[argminCostCurrent];
                    pathSucc[k1].push_back(k1);
                }
                else{
                    // no successor found, no path exist
                    costSucc[k1] = DBL_MAX;
                    pathSucc[k1].clear();
                }
            }
            // copy costSucc in costPrev
            for(k1=0; k1 < inputData->cloudletCard; k1++)
            {
                costCurrent[k1] = costSucc[k1];
                pathCurrent[k1] = pathSucc[k1];
            }
        }
        
        double reducedCost = -1;
        vector<int> minPath;
        
        for(k1=0; k1 < inputData->cloudletCard; k1++) if(inputData->isApCloudletPossible[i][k1])
        {
            // check if the assignment AP-|T|-CLOUDLET is forbidden --> if true skip the check
            if(this->fixedVar_small[i][inputData->timeSlotCard-1][k1]!=0)
            {
                // check if the assignment AP-|T|-CLOUDLET is forced --> if true assign and do not check all other assignment
                if(this->fixedVar_small[i][inputData->timeSlotCard-1][k1]==1)
                {
                    reducedCost = costCurrent[k1];
                    minPath = pathCurrent[k1];
                    break;
                }
                else if(costCurrent[k1] < reducedCost || reducedCost == -1)
                {
                    reducedCost = costCurrent[k1];
                    minPath = pathCurrent[k1];
                }
            }
        }
        
        // subtract constant eta dual to the reduced cost
        reducedCost -= duals.eta[i];
        
        // check if reduced cost is negative and if the path has not been already inserted
        if(reducedCost < -1e-6 && minPath.size() > 0) // && !checkIfPathAlreadyExist(i, minPath))
        {
            // add assignment path to the path to return
            AssociationPath newPath(i,reducedCost, minPath);
            assocPaths.push_back(newPath);
        }
    }
    return assocPaths;
}

bool ShortPathHeur::fixVariablesThisBranch(BranchingNode* node)
{
    bool result = true;
    // I have to fix the given variable list, and to free all other variables
    // I free all z variables
    // then for each original variable to fix in the list I fix again the corresponding z variables bounds
    
    for(int i = 0; i < z.getSize(); i++)
    {
        for(int p = 0 ; p < z[i].getSize(); p++)
        {
            // I only set z variables UB to 0 --> to free variable I change all UB to 1
            z[i][p].setUB(1.0);
        }
    }
    // re-initialize fixed variables structure to use during pricing
    fixedVar_small = vector<vector<vector<int8_t>>>(inputData->apCard, vector<vector<int8_t>>(inputData->timeSlotCard,
                                                                                vector<int8_t>(inputData->cloudletCard,2)));
    
    isTransitionPossible = vector<vector<vector<vector<bool>>>>(inputData->apCard, vector<vector<vector<bool>>>(inputData->cloudletCard, vector<vector<bool>>(inputData->cloudletCard, vector<bool>(inputData->timeSlotCard,true))));
    
    // printBranchTree(node);
    
    while(node)
    {
        // fix node variables
        list<bool>::iterator itB = node->atZero.begin();
        for(list<int>::iterator it = node->xVarIdxInList.begin(); it != node->xVarIdxInList.end(); ++it)
        {
            int xtemp = (*it);
            XvarIndexes xvarToFix = this->xIndexUnlist[xtemp];
            bool fixedAtZero = (*itB);
            
            // check consistency of branch fixing
            if((fixedVar_small[xvarToFix.ap][xvarToFix.time][xvarToFix.cloudlet] == 0 && !fixedAtZero) ||
               (fixedVar_small[xvarToFix.ap][xvarToFix.time][xvarToFix.cloudlet] == 1 && fixedAtZero) )
            {
                return false;
            }
            
            if(!fixedAtZero){
                // a variable has to be fixed to 1, check if for the same pair AP-TIME another assignment is already fixed
                for(int k = 0; k < inputData->cloudletCard; k++) if(k != xvarToFix.cloudlet)
                {
                    if(fixedVar_small[xvarToFix.ap][xvarToFix.time][k] == 1)
                        return false;
                }
            }
            
            fixVariable(xvarToFix, fixedAtZero?AtZero:AtOne);
            ++itB;
        }
        
        // forbid transitions
        for(list<Transition>::iterator it = node->forbiddenTransition.begin(); it != node->forbiddenTransition.end(); ++it)
        {
            Transition tempTrans = (*it);
            isTransitionPossible[tempTrans.ap][tempTrans.cloudletPrev][tempTrans.cloudletSucc][tempTrans.timeSucc] = false;
        }

        node = node->father;
    }
    return result;
}

// populate IloModel with variable and constraints
void ShortPathHeur::populateModel(IloEnv env, InputDataApCloudletAssoc* inputDataIn, bool bpSingleVarFix)
{
    int i,t,k;
    
    this->inputData = inputDataIn;
    
    xIndexUnlistMapping = vector<vector<vector<int>>>(inputData->apCard, vector<vector<int>>(inputData->cloudletCard,vector<int>(inputData->timeSlotCard,-1)));
    // save unlisted x indexes
    for(i=0; i < inputData->apCard;i++)
    {
        for(k=0; k< inputData->cloudletCard;k++)
        {
            for(t= 0; t <inputData->timeSlotCard;t++)
            {
                xIndexUnlist.push_back(XvarIndexes(i, k, t));
                xIndexUnlistMapping[i][k][t] = ((int) xIndexUnlist.size()) - 1;
            }
        }
    }
    
    this-> branchingRoot = new BranchingNode(); //(bpSingleVarFix);
    
    // LP master model
    this->model = IloModel(env);
    
    // boolean stating if last CG execution was successfull
    this->isCgOK = false;
    // boolean stating if last CG execution solution is integer
    this->isIntegerSolution = false;
    
    // dual variables structure
    this->dualVariables = DualVariables(inputData->apCard, inputData->cloudletCard, inputData->timeSlotCard);
    
    // empty path variables for each AP
    z = IloArray<IloNumVarArray>(env, inputData->apCard);
    for(i = 0; i < inputData->apCard; i++)
    {
        z[i] = IloNumVarArray(env, 0);
    }
    // corresponding variable values
    zval = vector< vector<double> >(inputData->apCard, vector<double>(0));
    // corresponding association paths
    addedPaths = vector< vector<AssociationPath> >(inputData->apCard,vector<AssociationPath>());
    
    // dummy variable
    dummy = IloNumVar(env, 0, 1, "dummy");
    dummyVal = 0.0;
    isDummyUsed = true;
    // get objective function contribution of the dummy variable
    double dummyCost = inputData->getDummyVarOFcost();
    
    // OF (initial with the dummy variable alone)
    objFun = IloMinimize(env, dummyCost*dummy, "objFun");
    model.add(objFun);
    
    // CONSTRAINTS
    
    // cloudlet capacity constraints
    // - sum{i in A} sum{p in Omega[i]} d[i,t]*x[i,k,t]*z[p] >= - U*C  forall t in T, forall k in K
    cloudletCapacity = IloArray<IloRangeArray>(env, inputData->timeSlotCard);
    for(t = 0; t < inputData->timeSlotCard; t++)
    {
        cloudletCapacity[t] = IloRangeArray(env, inputData->cloudletCard);
        for(k = 0; k < inputData->cloudletCard; k++)
        {
            char name[100] = "";
            snprintf(name, 100, "cloudletCapacity_%d_%d", t, k );
            cloudletCapacity[t][k] = IloRange(env,-inputData->U*inputData->cloudletCap,+IloInfinity,name);
        }
        model.add(cloudletCapacity[t]);
    }
    
    // for each AP, one and only one assignment path has to be chosen
    // dummy variable is used if no z variables are present
    // dummy + sum{p in Omega[i]} z[i] = 1 forall i in A
    valorizeZ = IloRangeArray(env, inputData->apCard);
    for(i = 0; i < inputData->apCard; i++)
    {
        char name[100] = "";
        snprintf(name, 100, "valorizeZ_%d", i);
        valorizeZ[i] = IloRange(env,1,dummy,1,name);
    }
    model.add(valorizeZ);
    
}

// constructor given initial structures for the assignment paths to add and variable to fix
ShortPathHeur::ShortPathHeur(IloEnv env, InputDataApCloudletAssoc* inputData, BranchFixedVariables fixedVariableIn,
                             vector<vector<AssociationPath>> pathToAdd, bool bpSingleVarFix)
{
    // populate model
    this->populateModel(env, inputData, bpSingleVarFix);
    
    // branchingRoot = new BranchingNode(bpSingleVarFix);
    
    // copy added paths and add corresponding z variables
    for(int i=0; i< pathToAdd.size();i++)
    {
        for(int p = 0; p < pathToAdd[i].size();p++)
        {
            // add column corresponding to the greedy solution
            this->addColumn(pathToAdd[i][p]);
            this->addedPaths[pathToAdd[i][p].nodeId].push_back(pathToAdd[i][p]);
        }
    }
    
    // copy fixed variable structure
    this->copyFixedVariables(&fixedVar_small, fixedVariableIn);
    
    // fix new z variables
    for(int i = 0; i < inputData->apCard;i++)
    {
        for(int t =0; t< inputData->timeSlotCard;t++)
        {
            for(int j = 0; j < inputData->cloudletCard;j++)
            {
                if(fixedVar_small[i][t][j] == 0)
                {
                    this->fixVariable(XvarIndexes(i, j, t), AtZero);
                }
                else if(fixedVar_small[i][t][j] == 1)
                {
                    this->fixVariable(XvarIndexes(i, j, t), AtOne);
                }
            }
        }
    }
}

// copy fixed variable data structure
void ShortPathHeur::copyFixedVariables(BranchFixedVariables * copyTo, BranchFixedVariables copyFrom)
{
    (*copyTo).resize(copyFrom.size());
    for(int i=0; i < copyFrom.size();i++)
    {
        (*copyTo)[i].resize(copyFrom[i].size());
        for(int j=0;j < copyFrom[i].size();j++)
        {
            (*copyTo)[i][j].resize(copyFrom[i][j].size());
            for(int k=0; k < copyFrom[i][j].size(); k++)
            {
                (*copyTo)[i][j][k] = copyFrom[i][j][k];
            }
        }
    }
}

// costructor without any initial assignment paths and fixed variable --> to initialize the problem an heuristic integer solution is computed
ShortPathHeur::ShortPathHeur(IloEnv env,InputDataApCloudletAssoc* inputData, bool bpSingleVarFix)
{
    // branchingRoot = new BranchingNode(bpSingleVarFix);
    
    this->populateModel(env, inputData,bpSingleVarFix);
    
    int i,t,k;
    
    // execute greedy algorithm to initialize master model
    Solution bestIntegerSolution = inputData->getBestIntegerSolution();
    if(bestIntegerSolution.finalCost!=-1)
    {
        // fix dummy variable to get value 0 if initial feasible solution is found
        dummy.setUB(0.0);
        isDummyUsed = false;
        
        VarValues greedyVarValues = bestIntegerSolution.variableValues;
        
        // for each AP build the unique association path created by the greedy heuristic
        for(i = 0; i < inputData->apCard;i++)
        {
            AssociationPath greedyPath;
            greedyPath.nodeId = i;
            
            // build assignment paths to create CG master columns
            for(t = 0; t < inputData->timeSlotCard;t++)
            {
                for(k = 0; k< inputData->cloudletCard;k++)
                {
                    if(greedyVarValues.xval[i][k][t] > 1e-6)
                    {
                        greedyPath.path.push_back(k);
                        break;
                    }
                }
            }
            // add column corresponding to the greedy solution
            this->addColumn(greedyPath);
            this->addedPaths[i].push_back(greedyPath);
        }
    }
    
    // forcedAssignCloudlets contains for each AP and each time-slot, the cloudlet that is forced to be assigned
    // or -1 if no cloudlet is forced
    fixedVar_small = vector<vector<vector<int8_t>>>(inputData->apCard, vector<vector<int8_t>>(inputData->timeSlotCard,vector<int8_t>(inputData->cloudletCard,2)));
    isTransitionPossible = vector<vector<vector<vector<bool>>>>(inputData->apCard, vector<vector<vector<bool>>>(inputData->cloudletCard, vector<vector<bool>>(inputData->cloudletCard, vector<bool>(inputData->timeSlotCard,true))));
}

// getter for the fixed variable data structure
BranchFixedVariables ShortPathHeur::getFixedVariables()
{
    return this->fixedVar_small;
}

// getter for the IloEnv
IloEnv ShortPathHeur::getEnv()
{
    return this->model.getEnv();
}

// check if fixed variables costraints are observed
bool ShortPathHeur::checkFixedVariableObserved(VarValues varValues,BranchFixedVariables fixedVar)
{
    bool result = true;
    
    for(int i=0; i < fixedVar.size();i++)
    {
        for(int t=0; t< fixedVar[i].size();t++)
        {
            for(int k=0; k< fixedVar[i][t].size();k++)
            {
                if(fixedVar[i][t][k]==0)
                {
                    if(varValues.xval[i][k][t] > 1e-6)
                    {
                        cout << " ERROR: var " << i << " " << k << " " << t << " NOT FIXED TO 0 : " << varValues.xval[i][k][t] << endl;
                        cerr << " ERROR: var " << i << " " << k << " " << t << " NOT FIXED TO 0 : " << varValues.xval[i][k][t] << endl;
                        
                        return false;
                    }
                }
                else if(fixedVar[i][t][k]==1)
                {
                    // double temp = varValues.xval[i][k][t];
                    if(varValues.xval[i][k][t] < 1 - 1e-6)
                    {
                        cout << " ERROR: var " << i << " " << k << " " << t << " NOT FIXED TO 1 : " << varValues.xval[i][k][t] << endl;
                        cerr << " ERROR: var " << i << " " << k << " " << t << " NOT FIXED TO 1 : " << varValues.xval[i][k][t] << endl;
                        
                        return false;
                    }
                }
            }
        }
    }
    return result;
}

// getter for the assignment paths added to the model as columns
vector< vector<AssociationPath> > ShortPathHeur::getAddedPaths()
{
    return this->addedPaths;
}

// check if current CG solution has integrality property
bool ShortPathHeur::hasIntegrality()
{
    return isIntegerSolution;
}

// given the current CG solution, compute an integer solution
bool ShortPathHeur::computeIntegerSolution(double & currentIntegerSolCos, double currentPB, bool isBP)
{
    bool isNewSol = false;
    
    if(this->isIntegerSolution)
    {
        // current CG solution is also integer, check if it is better than current best integer
        if(currentPB < currentIntegerSolCos || currentIntegerSolCos == -1)
        {
            Solution newIntegerSol;
            newIntegerSol.variableValues = originalVarValues;
            newIntegerSol.finalCost = currentPB;
            inputData->setBestIntegerSolution(newIntegerSol);
            
            currentIntegerSolCos = currentPB;
            
            isNewSol = true;
        }
    }
    else
    {
        // otherwise, from the current fraction PB solution compute a new integer solution via rounding
        Solution currentIntegerSol;
        
        currentIntegerSol = inputData->getIntegerSolutionViaRounding(originalVarValues);
        isNewSol = this->_checkNewIntegerSolution(currentIntegerSolCos, currentIntegerSol);
        
        if(isBP && (currentIntegerSolCos - currentPB)/currentPB > 0.001 )
        {
            currentIntegerSol = inputData->getIntegerSolutionViaRounding(originalVarValues, fixedVar_small);
            isNewSol = isNewSol || this->_checkNewIntegerSolution(currentIntegerSolCos, currentIntegerSol);
        }
    }
    return isNewSol;
}

DualVariables ShortPathHeur::getDuals(){
    return this->dualVariables;
}

void ShortPathHeur::tryLagrangianProbingOnNode(BranchingNode * branchToCheck, DualVariables nodeDuals, double LB)
{
    if(this->fixVariablesThisBranch(branchToCheck))
    {
        // double LB = cplex->getObjValue();
        double UB = inputData->getBestIntegerCost();
        
        // try forbidding new transitions and assignments
        vector<int> newxVarToForbid;
        vector<Transition> newTransitionToForbid;
        vector<bool> newXvarAreToZero;
        this->_lagrangianProbing(nodeDuals, newxVarToForbid, newXvarAreToZero,newTransitionToForbid, LB, UB);
        
        if(branchToCheck->level == 0)
            cout << "lagrangian fixing at root" << endl;
        cout << " lagrangian fixing new x var to fix " << newxVarToForbid.size() << endl;
        cout << " lagrangian fixing new transitions to fix " << newTransitionToForbid.size() << endl;
        branchToCheck->addForbiddenVarAndTransitions(newTransitionToForbid, newxVarToForbid, newXvarAreToZero);
    }
}

void ShortPathHeur::tryLagrangianProbingOnNode(BranchingNode* branchToCheck, IloCplex *cplex)
{
    if(this->fixVariablesThisBranch(branchToCheck))
    {
       if(cplex->solve())
       {
           DualVariables nodeDuals = DualVariables(inputData->apCard, inputData->cloudletCard, inputData->timeSlotCard);
           updateDualValues(cplex, nodeDuals);
           
           double LB = cplex->getObjValue();
           
           tryLagrangianProbingOnNode(branchToCheck, nodeDuals, LB);
       }
    }
}


bool ShortPathHeur::_checkNewIntegerSolution(double & currentIntegerSolCos, Solution newSol)
{
    if(newSol.finalCost != -1)
    {
        if(newSol.finalCost < currentIntegerSolCos || currentIntegerSolCos == -1)
        {
            currentIntegerSolCos = newSol.finalCost;
            inputData->setBestIntegerSolution(newSol);
            _checkPathIntegerSolution(newSol);
            return true;
        }
    }
    return false;
}

void ShortPathHeur::_checkPathIntegerSolution(Solution sol)
{
    // get integer solution path
    for(int i = 0; i < inputData->apCard; i++)
    {
        vector<int> assingPath;
        for(int t = 0; t < inputData->timeSlotCard; t++)
        {
            for(int k = 0; k < inputData->cloudletCard; k++)
            {
                // found the cloudlet assigned to the
                if(sol.variableValues.xval[i][k][t] > 1 - 1e-6)
                {
                    assingPath.push_back(k);
                    break;
                }
            }
        }
        if(assingPath.size() < inputData->timeSlotCard)
        {
            break;
        }
        
        // check if a column with the same path already exist
        // if true -> set it to the 1 in the solution
        // otherise -> add the new column to the model and set it to 1
        vector<int> tempPath;
        bool found = true;
        for(int p = 0; p < addedPaths[i].size(); p++)
        {
            tempPath = addedPaths[i][p].path;
            found = true;
            for(int t = 0; t < assingPath.size(); t++)
            {
                if(assingPath[t] != tempPath[t])
                {
                    found = false;
                    break;
                }
            }
            if(found)
            {
                break;
            }
        }
        if(!found)
        {
            // assignment path is not a column of the model -> add the column
            AssociationPath newPath = AssociationPath(i, -1, assingPath);
            addColumn(newPath);
            this->addedPaths[i].push_back(newPath);
        }
    }
}

void ShortPathHeur::printBranchTree(BranchingNode * nodeLabel)
{
    if(nodeLabel)
    {
        if(nodeLabel->father)
        {
            printBranchTree(nodeLabel->father);
        }
        if(nodeLabel->level > 0)
        {
            cout << "lev. " << nodeLabel->level << endl;
        }
        else{
            cout << "root" << endl;
        }
        
        if(nodeLabel->xVarIdxInList.size() > 0){
            list<bool>::iterator itB = nodeLabel->atZero.begin();
            for(list<int>::iterator it = nodeLabel->xVarIdxInList.begin(); it != nodeLabel->xVarIdxInList.end(); ++it)
            {
                xIndexUnlist[(*it)].print();
                cout << " at " << ((*itB)?"0":"1");
                cout << endl;
                ++itB;
            }
        }
    }
}

// recursive function to destroy tree of branching labels
void ShortPathHeur::destroyBranchingTree(BranchingNode* nodeToDestroy)
{
    if(nodeToDestroy)
    {
        if(nodeToDestroy->leftson)
        {
            destroyBranchingTree(nodeToDestroy->leftson);
        }
        if(nodeToDestroy->rightson)
        {
            destroyBranchingTree(nodeToDestroy->rightson);
        }
        delete nodeToDestroy;
        nodeToDestroy = nullptr;
    }
}
