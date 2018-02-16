//
//  ShortPathHeurBP.cpp
//  ap_cloudlet_association
//
//  Created by Marco Premoli on 14/02/2017.
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

#include "ShortPathHeurBP.hpp"

ShortPathHeurBP::~ShortPathHeurBP(){}

// constructor
ShortPathHeurBP::ShortPathHeurBP(IloEnv env, InputDataApCloudletAssoc* inputData, BranchingRule branchRule): rootCG(ShortPathHeur(env,inputData, branchRule != SplitAssignment ))
{
    this->inputData = inputData;
    // save environment to create new models in the branches
    this->env = env;
    
    this->branchingRule = branchRule;
    
    // this->singleFixedXinBranch = branchingRule != SplitAssignment;
}

bool ShortPathHeurBP::getSplitAssignment(VarValues fatherVarValues,BranchFixedVariables fixedVariables,list<int> & firstSet, list<int> & secondSet)
{
    bool splitDone = false;
    int ap = -1;
    int time = -1;
    
    // get most fractional x variable
    // int xmostFract = getMostFractionalVar(fatherVarValues);
    
    if(this->getMostDividedAssignment(fatherVarValues, ap, time) > 0) //(xmostFract!=-1)
    {
        // get corresponding triplet AP-cloudlet-time
        // XvarIndexes xvarMostFract = rootCG.xIndexUnlist[xmostFract];
        
        // ap = xvarMostFract.ap;
        // time = xvarMostFract.time;
        
        // sort all NON-ALREADY-FIXED x variables of the same pair AP-TIME by descending value
        vector<double> xval;
        vector<int> xidx;
        for(int i=0; i < rootCG.xIndexUnlist.size(); i++)
        {
            XvarIndexes tempx = rootCG.xIndexUnlist[i];
            if(tempx.ap == ap && tempx.time == time
               && fixedVariables[ap][time][tempx.cloudlet] != 0
               && fixedVariables[ap][time][tempx.cloudlet] != 1)
            {
                xidx.push_back(i);
                xval.push_back(fatherVarValues.xval[ap][tempx.cloudlet][time]);
            }
        }
        
        if(xval.size() <=0 )
        {
            return splitDone;
        }
        
        splitDone = true;
        if(xval.size() > 1){
            // sort x variables
            vector<size_t> xvalSort = InputDataApCloudletAssoc::sort_indexes(xval, true);
            
            // assign sorted x variable to returning sets alternating from one set to the other
            // ex. : the first-third-fitht x to the first set, the second-fourth-sixtht x to the second set
            for(int i = 0; i < xvalSort.size(); i++)
            {
                if(i % 2 == 0){
                    firstSet.push_back(xidx[(int) xvalSort[i]]);
                }
                else{
                    secondSet.push_back(xidx[(int) xvalSort[i]]);
                }
            }
        }
        else{
            firstSet.push_back(xidx[0]);
        }
    }
    
    return splitDone;
}

bool ShortPathHeurBP::addBranchesToQueueBranchNode(stack<BranchingNode*, vector<BranchingNode*>>* model_stack ,
                                  priority_queue<BranchingNode*, vector<BranchingNode*>, CompareBranchingNodes>* pq,
                                  bool usePriorityQueue, VarValues fatherVarValues, double fatherLB, int * nodes_to_do, BranchingNode* fatherNode, BranchFixedVariables fatherFixedVar, BranchingRule currentBranchingRule, vector<double> & all_nodes_lb)
{
    bool atLeastOneBranchAdded = false;
    
    if(currentBranchingRule != SplitAssignment)
    {
    
        int xToFixIdx = -1;
        if(branchingRule==HighestFractional)
            xToFixIdx = getHighestFractionalVar(fatherVarValues);
        else
            xToFixIdx = getMostFractionalVar(fatherVarValues);
        
        // if chosen variable is valid (i.e. if it exist) add two new branches in the queue, to fix the variable in both directions
        if(xToFixIdx > -1)
        {
            for(int d = 0; d < 2; d++)
            {
                bool isLeft = d==0;
                list<int> xToFixIdxList;
                list<bool> atZero;
                xToFixIdxList.push_back(xToFixIdx);
                atZero.push_back(isLeft);
                
                all_nodes_lb.push_back(fatherLB);
                BranchingNode* newNode = fatherNode->addNewBranch(xToFixIdxList, atZero, isLeft,
                                                                  fatherLB, (int) all_nodes_lb.size() - 1);
                
                if(!usePriorityQueue)
                    model_stack->push(newNode);
                else
                    pq->push(newNode);
                
                // keep track of the number of branching nodes
                (*nodes_to_do)++;
                atLeastOneBranchAdded = true;
            }
        }
    }
    else{
        list<int> firstSet;
        list<int> secondSet;
        if(getSplitAssignment(fatherVarValues, fatherFixedVar, firstSet, secondSet))
        {
            if(firstSet.size() > 0 || secondSet.size() > 0)
                atLeastOneBranchAdded = true;
            
            if(firstSet.size() > 0)
            {
                list<bool> atZero(firstSet.size(), true);
                // in left branch the first set of x variables are set to 0
                all_nodes_lb.push_back(fatherLB);
                BranchingNode* leftNode = fatherNode->addNewBranch(firstSet, atZero ,true, fatherLB, (int) all_nodes_lb.size() - 1);
                
                if(!usePriorityQueue)
                    model_stack->push(leftNode);
                else
                    pq->push(leftNode);
                (*nodes_to_do)++;
            }
            
            if(secondSet.size() > 0)
            {
                list<bool> atZero(secondSet.size(), true);
                // in right branch the second set of x variables are set to 0
                all_nodes_lb.push_back(fatherLB);
                
                BranchingNode* rightNode = fatherNode->addNewBranch(secondSet, atZero, false, fatherLB,  (int) all_nodes_lb.size() - 1);
                
                if(!usePriorityQueue)
                    model_stack->push(rightNode);
                else
                    pq->push(rightNode);
                (*nodes_to_do)++;
            }
        }
    }
    return atLeastOneBranchAdded;
}

bool ShortPathHeurBP::solve(ExplorationStrategy explorationStrategy)
{
    bool branchResult = false;
    
    // get unique IloCplex instance
    cplex = IloCplex(rootCG.model);
    cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Algorithm::Primal);
    
    // solve root of CG
    branchResult = rootCG.solve(XvarIndexes(), None,-1, &cplex, false);
    
    int evaluatedNodes = 0;
    
    // the branching start if and only if the CG terminated sucecssfully and a feasible integer solution is available
    if(branchResult && inputData->getBestIntegerCost() != -1)
    {
    
        // save CG LB
        double cgLB = rootCG.getFinalLB();
        inputData->setCGRootLB(cgLB);
        
        // print CG solution and first integer solution
        cout << " CG SOLUTION " << endl;
        cout << " final CG cost " << cgLB << endl;
        
        // print original variable values
        InputDataApCloudletAssoc::printVarValues(rootCG.getFinalVarValues());
        
        cout << "INTEGER SOLUTION AFTER CG" << endl;
        inputData->printBestIntegerSolution();
        
        cout << "end of CG in " << inputData->getCurrentExcTime() << " seconds " << endl;
        cerr << "end of CG in " << inputData->getCurrentExcTime() << " seconds " << endl;
        // start B&P
        
        cout << "BRANCH-AND-PRICE START" << endl;
        cerr << "BRANCH-AND-PRICE START" << endl;
        
        bool usePriorityQueue = (explorationStrategy == PriorityQueue);
        
        stack<BranchingNode*, vector<BranchingNode*>> model_stack;
        priority_queue<BranchingNode*, vector<BranchingNode*>, CompareBranchingNodes> pq;
        
        // root of B&P has no fixed variables and the root LB
        rootCG.branchingRoot->fatherLB = cgLB;
        
        // nodes of the branching to visit
        int nodes_to_do = 0;
        // percentage integrality gap
        double perc_int_gap_max = 0.0005;
        
        vector<double> all_nodes_lb;
        all_nodes_lb.push_back(-1);
        
        //rootCG.tryLagrangianFixingOnNode(rootCG.branchingRoot, &cplex);
        
        // try to get a variable to fix from the current node and create two more branching nodes, added to the stack (or priority queue)
        if(!this->addBranchesToQueueBranchNode(&model_stack, &pq, usePriorityQueue, rootCG.getFinalVarValues(), cgLB, &nodes_to_do, rootCG.branchingRoot, rootCG.getFixedVariables(), this->branchingRule, all_nodes_lb))
        {
            // no more nodes added at the root -> exit B&P
            return branchResult;
        }
        
        bool empty_stack = true;
        
        // check if queue is empty
        if(!usePriorityQueue)
            empty_stack = model_stack.empty();
        else
            empty_stack = pq.empty();
        
        // continue until branches queue is not empty
        while(!empty_stack && !inputData->isExcTimeLimitReached())
        {
            // current model to solve, taken as top of the queue
            BranchingNode* current_node;
            
            if(!usePriorityQueue)
            {
                current_node = model_stack.top();
                model_stack.pop(); // remove top
            }
            else
            {
                current_node = pq.top();
                pq.pop(); // remove top
            }
            nodes_to_do--;
            evaluatedNodes++;
            
            // after the pop check again if node father LB is higher than current UB (that may be changed after the current node was inserted)
            // for sure the son could not have a higher LB
            double father_gap = (inputData->getBestIntegerCost() - current_node->fatherLB)/current_node->fatherLB;
            
            if(father_gap >= perc_int_gap_max)
            {
                // CHECK CONSISTENCY OF FIXED VARIABLES OF CURRENT NODE BRANCH
                // WHILE CHANGING FIXED VARIABLES OF CURRENT MODEL EXECUTION
                if(rootCG.fixVariablesThisBranch(current_node))
                {
                    // SOLVE CG
                    bool modelResult = rootCG.solve( inputData->getBestIntegerCost() , &cplex);
                    
                    if(modelResult)
                    {
                        double finalLB = rootCG.getFinalLB();
                        double gap = (inputData->getBestIntegerCost() - finalLB)/finalLB;
                        if(gap >= perc_int_gap_max)
                        {
                            // lower bound of the CG is lower than the current best upper bound -> NOT PRUNE
                            // add the current model in the list of model, ordered by descending lower bound
                            
                            VarValues currentVarValues = rootCG.getFinalVarValues();
                            
                            // check if fixed variables are actually fixed
                            if(!ShortPathHeur::checkFixedVariableObserved(currentVarValues, rootCG.getFixedVariables()))
                            {
                                // something wrong... abort execution
                                exit(-1);
                            }
                            
                            // CHECK IF A STRENGHTING ON VARIABLE FIXING CAN BE DONE WITH LAGRANGIAN FIXING
                            if(gap < perc_int_gap_max*1.5)
                            {
                                cout << " branching node gap " << gap << " : try lagrangian fixing " << endl;
                                cerr << " branching node gap " << gap << " : try lagrangian fixing " << endl;
                                rootCG.tryLagrangianProbingOnNode(current_node, rootCG.getDuals(), finalLB);
                            }

                            // check what is the highest number of split assignment
                            // if it is less than a certain percentage of the total number of cloudlet branch on single variable (the most fractional)
                            int t_a, t_t;
                            int maxSplit = this->getMostDividedAssignment(currentVarValues, t_a, t_t);
                            if(maxSplit  == 1) // || gap <= perc_int_gap_max*2)
                            {
                                cout << " most fractional branching rule " << endl;
                                this->branchingRule = MostFractional;
                            }
                            else{
                                cout << " split assignment branching rule " << endl;
                                this->branchingRule = SplitAssignment;
                            }
                            
                            // ADD NEW BRANCH NODES WITH NEW FIXED VARIABLES
                            this->addBranchesToQueueBranchNode(&model_stack, &pq, usePriorityQueue, currentVarValues, finalLB, &nodes_to_do, current_node, rootCG.getFixedVariables(), this->branchingRule, all_nodes_lb);
                        }
                        else
                        {
                            // lower bound of the CG is greater than the current best upper bound -> PRUNE THE BRANCH
                            cout << "prune branch LB = " << finalLB << " best UB = " << inputData->getBestIntegerCost() << " GAP = " << gap << endl;
                            cerr << "prune branch LB = " << finalLB << " best UB = " << inputData->getBestIntegerCost() << " GAP = " << gap << endl;
                            // inputData->printVarValues( rootCG.getFinalVarValues());
                            // break;
                        }
                    }
                    else{
                        // CG was not successfull -> PRUNE THE BRANCH
                        cout << "prune branch infeasible" << endl;
                        cerr << "prune branch infeasible" << endl;
                    }
                }
                else
                {
                    cout << "prune branch variable fixing inconsistency " << endl;
                    cerr << "prune branch variable fixing inconsistency " << endl;
                    // exit(-1);
                }
                
                if(!usePriorityQueue)
                {
                    empty_stack = model_stack.empty();
                }
                else
                {
                    empty_stack = pq.empty();
                }
                cout << "nodes to evaluate in stack " << nodes_to_do << endl;
                cerr << "nodes to evaluate in stack " << nodes_to_do << endl;
            }
            
            all_nodes_lb[current_node->nodeNumber] = -1;
            if(current_node->fatherLB == inputData->getMinActiveBranchNodeLB())
            {
                // update minimum active branch node LB
                double current_min = -1;
                if(usePriorityQueue)
                {
                    current_min = (pq.top())->fatherLB;
                }
                else{

                    for(int n = 0; n < all_nodes_lb.size(); n++)
                    {
                        if( (all_nodes_lb[n] >= 0 && all_nodes_lb[n] < current_min)  ||
                           (all_nodes_lb[n] >= 0 && current_min == -1))
                        {
                            current_min = all_nodes_lb[n];
                        }
                    }
                }
                if(current_min != -1)
                {
                    inputData->setMinActiveBranchNodeLB(current_min);
                }
            }
            
            if(!usePriorityQueue)
            {
                double excTime = inputData->getCurrentExcTime();
                if(excTime > -1 && excTime > 1800)
                {
                    switchToPriorityQueue(model_stack, pq, usePriorityQueue);
                }
            }
            
        }
        
        cout << "last min DB of active branching node " << inputData->getMinActiveBranchNodeLB() << endl;
        cout << "number of branching nodes considered " << all_nodes_lb.size() - 1 << " evaluated " << evaluatedNodes << " todo " << nodes_to_do <<  endl;
    }
    cplex.end();
        
    return branchResult;
}

int ShortPathHeurBP::getMostDividedAssignment(VarValues cgVarValues,int & ap, int & time)
{
    bool found = false;
    
    int maxCardFractX = -1;
    ap = time = -1;
    
    for(int i = 0; i < this->inputData->apCard; i++)
    {
        for(int t=0; t < this->inputData->timeSlotCard; t++)
        {
            int tempFractXcard = 0;
            for(int k = 0; k < this->inputData->cloudletCard; k++)
            {
                double xval = cgVarValues.xval[i][k][t]; // get temporary variable with x value for debug purpose
                if(IS_FRACTIONAL(xval))
                {
                    tempFractXcard++;
                }
            }
            if(tempFractXcard > maxCardFractX)
            {
                ap = i;
                time = t;
                maxCardFractX = tempFractXcard;
                found = true;
            }
        }
    }
    
    return maxCardFractX;
}

int ShortPathHeurBP::getMostFractionalVar(VarValues cgVarValues)
{
    double mingap = 1;
    int maxXidx = -1;
    for(int i = 0; i < this->rootCG.xIndexUnlist.size(); i++)
    {
        XvarIndexes currentX = rootCG.xIndexUnlist[i];
        int ap = currentX.ap;
        int cloud = currentX.cloudlet;
        int time =  currentX.time;
        double xval = cgVarValues.xval[ap][cloud][time];
        if(IS_FRACTIONAL(xval) &&
           fabs(xval - 0.5) < mingap)
        {
            mingap = fabs(xval - 0.5);
            maxXidx = i;
            if(mingap == 0)
                return maxXidx;
        }
    }
    
    return maxXidx;
}

int ShortPathHeurBP::getHighestFractionalVar(VarValues cgVarValues)
{
    double maxValue = 1e-12;
    int maxXidx = -1;
    for(int i = 0; i < this->rootCG.xIndexUnlist.size(); i++)
    {
        XvarIndexes currentX = rootCG.xIndexUnlist[i];
        int ap = currentX.ap;
        int cloud = currentX.cloudlet;
        int time =  currentX.time;
        double xval = cgVarValues.xval[ap][cloud][time];
        if(IS_FRACTIONAL(xval) &&
           xval > maxValue)
        {
            maxValue = xval;
            maxXidx = i;
        }
    }
    
    return maxXidx;
}

void ShortPathHeurBP::switchToPriorityQueue(stack<BranchingNode*, vector<BranchingNode*>> & stackMod,
                           priority_queue<BranchingNode*, vector<BranchingNode*>, CompareBranchingNodes> & pq,
                           bool & usePriorityQueue)
{
    usePriorityQueue = true;
    
    while(!stackMod.empty()){
        BranchingNode* current_node = stackMod.top();
        pq.push(current_node);
        stackMod.pop();
    }
    
}
