//
//  ShortPathHeur.hpp
//  ap_cloudlet_association
//
//  Created by Marco Premoli on 31/01/2017.
//  Copyright © 2017 Marco Premoli. All rights reserved.
//
#pragma once
#ifndef ShortPathHeur_hpp
#define ShortPathHeur_hpp

#include "InputDataApCloudletAssoc.hpp"
// #include "supportFunctions.hpp"
#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h>
#include <list>

enum FixingDirection{ AtZero, AtOne, None };

// branching node for branching binary tree
struct BranchingNode
{
    // bool singleVar;
    list<bool> atZero;
    list<int> xVarIdxInList;
    list<Transition> forbiddenTransition;
    
    double fatherLB;
    int level;
    
    int nodeNumber;
    
    BranchingNode * father;
    BranchingNode * leftson;
    BranchingNode * rightson;
    
    BranchingNode(): father(nullptr), leftson(nullptr), rightson(nullptr), fatherLB(-1), level(0), nodeNumber(0) {}
    
    BranchingNode* addNewBranch(list<int> xvarIdx, list<bool> varFixAtZero,
                                bool isleftSon, double lb, int nodeNumber = -1)
    {
        BranchingNode* newNode = new BranchingNode();
        newNode->father = this;
        newNode->level = this->level + 1;
        newNode->fatherLB = lb;
        newNode->nodeNumber = nodeNumber;
        
        for(list<int>::iterator it = xvarIdx.begin(); it != xvarIdx.end(); ++it)
        {
            newNode->xVarIdxInList.push_back(*it);
        }
        
        for(list<bool>::iterator it = varFixAtZero.begin(); it != varFixAtZero.end(); ++it)
        {
            newNode->atZero.push_back(*it);
        }
       
        if(isleftSon)
        {
            // check if it is left or right son
            this->leftson = newNode;
        }
        else
        {
            this->rightson = newNode;
        }
        
        return newNode;
    }
    
    void addForbiddenVarAndTransitions(vector<Transition> newForbiddenTransition, vector<int> newXvarToFix, vector<bool> newXvarAreAtZero)
    {
        // cout << "new forbidden assignment transitions: " << endl;
        for(int i =0; i < newForbiddenTransition.size(); i++)
        {
            // newForbiddenTransition[i].print();
            this->forbiddenTransition.push_back(newForbiddenTransition[i]);
        }
        
        addXvarToFix(newXvarToFix,newXvarAreAtZero);
    }
    
    void addXvarToFix(vector<int> newXvarToFix, vector<bool> newXvarAreFixedAtZero)
    {
        if(newXvarToFix.size() == newXvarAreFixedAtZero.size())
        {
            for(int i = 0; i < newXvarToFix.size(); i++)
            {
                this->xVarIdxInList.push_back(newXvarToFix[i]);
                this->atZero.push_back(newXvarAreFixedAtZero[i]);
            }
        }
    }
    
};

class ShortPathHeur
{
public:
    ShortPathHeur(IloEnv, InputDataApCloudletAssoc*, bool bpSingleVarFix = true);
    ShortPathHeur(IloEnv, InputDataApCloudletAssoc*, BranchFixedVariables, vector<vector<AssociationPath>>, bool bpSingleVarFix = true);
    ~ShortPathHeur();
    
    bool solve();
    bool solve(double actualBestUB,IloCplex*cplex);
    bool solve(XvarIndexes, FixingDirection, double actualBestUB);
    bool solve(XvarIndexes, FixingDirection, double actualBestUB, IloCplex* cplex, bool useDummy);
    
    VarValues getFinalVarValues();
    
    double getFinalLB();
    
    BranchFixedVariables getFixedVariables();
    // BranchFixedVariables getFixedVariables(BranchingNode*);
    
    IloEnv getEnv();
    
    static bool checkFixedVariableObserved(VarValues, BranchFixedVariables);
    
    vector< vector<AssociationPath> > getAddedPaths();
    
    bool hasIntegrality();
    
    bool fixVariablesThisBranch(BranchingNode*);
    
    IloModel model;
    
    // unlist the x variable in a single vector in order to sort by descending value
    vector<XvarIndexes> xIndexUnlist;
    
    BranchingNode * branchingRoot;
    
    void tryLagrangianProbingOnNode(BranchingNode*, DualVariables, double);
    void tryLagrangianProbingOnNode(BranchingNode*, IloCplex*);
    
    DualVariables getDuals();
    
private:
    bool isCgOK;
    double finalLB;
    
    InputDataApCloudletAssoc* inputData;
    
    VarValues originalVarValues;
    
    // VARIABLES
    
    // z variables, for each AP and for each path related to the AP
    IloArray<IloNumVarArray> z;
    vector< vector<double> > zval;
    
    IloNumVar dummy;
    double dummyVal;
    bool isDummyUsed;
    
    // OBJ FUN
    IloObjective objFun;
    
    // MASTER CONSTRAINTS
    IloArray<IloRangeArray> cloudletCapacity;
    
    IloRangeArray valorizeZ;
    
    DualVariables dualVariables;
    
    // add new column to the model related to the path given in input
    void addColumn(AssociationPath);
    
    void updateDualValues(IloCplex*, DualVariables&);
    
    // associations paths corresponding to a z variable
    vector< vector<AssociationPath> > addedPaths;
    
    vector<AssociationPath> getNewPathsPricing(DualVariables);
    
    void _lagrangianProbing(DualVariables, vector<int> &,vector<bool> &, vector<Transition> &, double currPB, double bestUB = -1);
    
    VarValues getOriginalVarValuesFromZ();
    
    bool isIntegerSolution;
    
    bool computeIntegerSolution(double & currentIntegerSolCos, double currentPB, bool isBP);
    
    /**** FOR BRANCHING PURPOSE  ****/
    
    BranchFixedVariables fixedVar_small;
    vector<vector<vector<vector<bool>>>> isTransitionPossible;
    
    void populateModel(IloEnv env, InputDataApCloudletAssoc* inputData, bool bpSingleVarFix);
    
    void fixVariable(XvarIndexes, FixingDirection);
    
    bool checkIfPathAlreadyExist(int ap, vector<int> path);
    
    void copyFixedVariables(BranchFixedVariables*, BranchFixedVariables);
    
    void destroyBranchingTree(BranchingNode* nodeToDestroy);
    
    void printBranchTree(BranchingNode* nodeLabel);
    
    bool _checkNewIntegerSolution(double &,Solution);
    
    void _checkPathIntegerSolution(Solution);
    
    vector<vector<vector<int>>> xIndexUnlistMapping;
    int _getXvarIndexInVector(XvarIndexes);
    
};

#endif /* ShortPathHeur_hpp */
