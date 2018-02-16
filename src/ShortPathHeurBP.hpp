//
//  ShortPathHeurBP.hpp
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

#ifndef ShortPathHeurBP_hpp
#define ShortPathHeurBP_hpp

#include "ShortPathHeur.hpp"
#include <string>
#include <functional>
#include <queue>
#include <stack>
#include <stdlib.h> 

enum ExplorationStrategy{ DepthFirst, PriorityQueue };

enum BranchingRule{ HighestFractional, MostFractional, SplitAssignment };

class CompareBranchingNodes{
public:
    // to sort branching node in ascending LB order
    bool operator() ( BranchingNode* & node1,  BranchingNode* & node2)
    {
        return node1->fatherLB >= node2->fatherLB;
    }
};

class ShortPathHeurBP{
public:
    ~ShortPathHeurBP();
    ShortPathHeurBP(IloEnv, InputDataApCloudletAssoc*, BranchingRule);
    
    bool solve(ExplorationStrategy);
    
private:
    
    IloCplex cplex;
    
    InputDataApCloudletAssoc * inputData;
    ShortPathHeur rootCG;
    
    IloEnv env;
    
    BranchingRule branchingRule;
    // bool singleFixedXinBranch;
    
    // get x variable to fix in the next branch as variable with maximum fractional value
    // XvarIndexes getVariableToFix(VarValues cgVarValues);
    int getHighestFractionalVar(VarValues cgVarValues);
    int getMostFractionalVar(VarValues cgVarValues);
    int getMostDividedAssignment(VarValues,int & ap, int & time);
    bool getSplitAssignment(VarValues, BranchFixedVariables, list<int> & firstSet, list<int> & secondSet );
    
    bool addBranchesToQueueBranchNode(stack<BranchingNode*, vector<BranchingNode*>>* ,
                                      priority_queue<BranchingNode*, vector<BranchingNode*>, CompareBranchingNodes>*,
                                      bool, VarValues, double, int *, BranchingNode*, BranchFixedVariables,
                                      BranchingRule, vector<double> &);
    
     // bool checkFixedVariableConsistencyFromNode(BranchingNode*);
    
    void switchToPriorityQueue(stack<BranchingNode*, vector<BranchingNode*>> &,
                               priority_queue<BranchingNode*, vector<BranchingNode*>, CompareBranchingNodes> &,
                               bool &);
    
};


#endif /* ShortPathHeurBP_hpp */
