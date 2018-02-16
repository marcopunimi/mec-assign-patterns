//
//  InputDataApCloudletAssoc.hpp
//  ap_cloudlet_association
//
//  Created by Marco Premoli on 17/01/2017.
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
#pragma once
#ifndef InputDataApCloudletAssoc_hpp
#define InputDataApCloudletAssoc_hpp

#include "supportFunctions.hpp"

#include <stdio.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <float.h>
#include <algorithm>
#include <math.h>

using namespace std;

//struct PathNode{
//    int timeSlot;
//    int cloudletId;
//    
//    // constructor
//    PathNode() : timeSlot(-1), cloudletId(-1) {}
//    PathNode(int cloudlet, int time) : timeSlot(time), cloudletId(cloudlet) { }
//    
//};

struct AssociationPath{
    vector<int> path;
    int nodeId;
    double cost;
    
    void printPath()
    {
        cout << nodeId << " -> ";
        for(int t = 0; t < path.size();t++)
        {
            cout << "(" << path[t] << "," << t << ")";
            if (t + 1 < path.size())
            {
                cout << " - ";
            }
        }
        cout << endl;
    }
    
    AssociationPath(): nodeId(-1), cost(0.0) {}
    AssociationPath(int AP, double costIn, vector<int> newPath){
        this->nodeId = AP;
        this->cost = costIn;
        for(int t = 0; t < newPath.size();t++){
            this->path.push_back(newPath[t]);
        }
    }
};

enum ModelVariant { Binary, Continuous, BinaryHeuristic,
    ShortestPathHeuristic, GeneralizedAssignHeuristic,
    ShortestPathBranchAndPrice, ShortestPathHeuristicCyclic,
    BinaryCyclic, ContinuousCyclic
};

struct VarValues{
    vector< vector< vector< double> > > xval;
    
    vector< vector< vector< vector<double> > > > yval;
    
    VarValues(){}
    VarValues(int apCard, int cloudCard, int timeCard){
        int nA=apCard;
        int nK = cloudCard;
        int nT = timeCard;
        this->xval = vector< vector< vector< double> > >(nA, vector< vector< double> >(nK, vector< double>(nT, 0.0)));
        this->yval = vector< vector< vector< vector<double> > > >(nA, vector< vector< vector< double> > >(nK, vector< vector< double> >(nK, vector< double>(nT, 0.0))));
    }
};

// structure containing a feasible solution for the problem, in terms of variable values and final cost
struct Solution
{
    VarValues variableValues;
    double finalCost;
    double assignmentCost;
    double switchingCost;
    bool isBinary;
    bool is_cycle;
    
    Solution(): finalCost(-1), isBinary(false), assignmentCost(-1), switchingCost(-1), is_cycle(false) {}
    Solution(int apCard, int cloudCard, int timeCard, bool is_cycle = false): variableValues(VarValues(apCard,cloudCard,timeCard)) , finalCost(-1), isBinary(false), assignmentCost(-1), switchingCost(-1) , is_cycle(is_cycle) {}
    
};

// data structure to  store if a variable of the compact formulation x[i,k,t] is free or fixed to 0 or fixed to 1
typedef vector<vector<vector<int8_t>>> BranchFixedVariables;

struct XvarIndexes
{
    int ap;
    int cloudlet;
    int time;
    
    XvarIndexes(): ap(-1),cloudlet(-1),time(-1){}
    XvarIndexes(int apIn,int cloudIn,int timeIn): ap(apIn),cloudlet(cloudIn),time(timeIn){}
    
    bool isValid()
    {
        return (ap != -1 && cloudlet != -1 && time != -1);
    }
    
    void print()
    {
        cout << "var x " << ap << " " << cloudlet << " " << time;
    }
};

struct Transition
{
    int ap;
    int cloudletPrev;
    int cloudletSucc;
    int timeSucc;
    
    Transition(): ap(-1), cloudletPrev(-1),cloudletSucc(-1), timeSucc(-1) {}
    Transition(int apIn, int cloudPrevIn, int cloudSuccIn, int timeIn): ap(apIn), cloudletPrev(cloudPrevIn),cloudletSucc(cloudSuccIn), timeSucc(timeIn) {}
    
    bool isValid(){
        return ap!=-1 && cloudletPrev != -1 && cloudletSucc != -1 && timeSucc != -1;
    }
    void print()
    {
        cout << " assignment transition " << ap << " " << cloudletPrev << " " << cloudletSucc << " t " << timeSucc << endl;
    }
};

struct DualVariables
{
    vector< vector<double> > lambda;
    vector<double> eta;
    
    DualVariables(){}
    DualVariables(vector< vector<double> > lambda, vector<double> eta) : lambda(lambda), eta(eta) {};
    DualVariables(int apCard, int cloudletCard, int timeslotCard)
    {
        this->lambda = vector< vector< double> >(timeslotCard,vector<double>(cloudletCard,0.0));
        this->eta = vector<double>(apCard, 0.0);
    }
};

class InputDataApCloudletAssoc{
public:
    InputDataApCloudletAssoc(string filenameParameters,
                             string filenameDistances,
                             string filenameDemands,
                             int timelimit = -1,
                             bool betazero = false,
                             double limitAssignmentDistance = 0);
    ~InputDataApCloudletAssoc();
    
    int apCard;         // cardinality APs
    int cloudletCard;   // cardinality cloudlets
    int timeSlotCard;   // cardinality time-slots
    
    double alpha;       // alpha parameter
    double beta;        // beta parameter
    int U;              // U parameter
    double cloudletCap;   // maximum cloudlet capacity
    double cloudletCapPerc; // percentage extra cloudlet capacity w.r.t. peak demand in time (>1)
    
    vector< vector< double > > nodeDemand;  //demand of node in a time-slot, for each AP a vector of demands in time
    vector< vector< double > > timeDemand;  //demand of node in a time-slot, for each time slot a vector of AP demands
    
    // for each time-slot, a vector with AP indexes sorted in descending demand in the time-slot, to use in the greedy heuristic 
    vector< vector< size_t > > sortedApIdxDescendingDemand;
    
    vector< vector< double > > distanceCloudletCloudlet; // distance btw pair of cloudlets
    vector< vector< double > > distanceApCloudlet; // distance btw APs and cloudlets
    vector< vector<size_t> > cloudletIdxSortedByApDistance; //
    
    vector< vector< bool > > isApCloudletPossible; // only Ap-cloudlet within distance threshold

    // for each AP keep the number of facilities it can reach given the threshold on assignment distance
    vector<int> ap_number_cloudlet_reached;
    // sort AP per the number of reachable facilities, in ascending order
    vector<size_t> sorted_ap_number_cloudlet_reached;
    
    bool isBetaZero;
    
    void printAllData();
    
    static void printVarValues(VarValues);
    
    bool readOk;
    
    bool single_assignment;
    
    bool greedy_heuristic(bool is_cycle = false);
    Solution getIntegerSolutionViaRounding(VarValues varValues, bool is_cycle = false);
    Solution getIntegerSolutionViaRounding(VarValues varValues, BranchFixedVariables, bool is_cycle = false);
    
    //VarValues getGreedyVarValues();
    
    template <typename T>
    static vector<size_t> sort_indexes(const vector<T> &v, bool descending);
    
    double computeCostGivenVariableValues(VarValues, bool is_cycle = false);
    
    double getBestIntegerCost();
    Solution getBestIntegerSolution(bool is_cycle = false);
    void setBestIntegerSolution(Solution, bool is_cycle = false);
    void printBestIntegerSolution();
    
    double getDummyVarOFcost();
    
    void setCGRootLB(double newRootLB);
    void setMinActiveBranchNodeLB(double newLB);
    double getMinActiveBranchNodeLB();
    double getIntegralityGAP();
    
    bool checkFeasibilitySolution(VarValues varValues);
    
    void setInitialExecutionTime(time_t);
    
    bool isExcTimeLimitReached();
    bool isTimeLimited();
    double getCurrentExcTime();
    
    void readIntegerSolutionFromFile(string datafile);
    
    void printAssignmentMatrixBestIntegerSolution();
    
    void getAssignmentAndSwitchingCostBestIntegerSolution(double & assignCost, double & switchingCost);
    void getAssignmentAndSwitchingCostSolution(double & assignCost, double & switchingCost, Solution);
    void getAssignmentAndSwitchingCostVarValues(double & assignCost, double & switchingCost,
                                                VarValues, bool is_cycle);
    
private:
    Solution _getIntegerSolutionViaRounding(VarValues varValues, bool is_cycle = false);
    
    void setLimitAssignmentDistance(double);
    
    double cg_root_LB;
    double min_active_branch_node_LB;
    double integralityGAP;
    
    Solution bestIntegerSolution;
    
    int nearestCloudletNotExceedingCapacity(vector<double> current_usage,
                                            vector<double> distances,
                                            double ap_demand,
                                            double maxCapacity);
    
    int nearestCloudletNotExceedingCapacity(vector<double> current_usage,
                                            vector<double> distances,
                                            double ap_demand,
                                            double maxCapacity,
                                            vector<bool> isAssignmentPossible);
    
    void computeYvarValuesGivenXvarValues(VarValues &);
    
    // dummy variable OF cost
    double dummyCost;
    //VarValues greedyVarValues;
    
    time_t initExecutionTime;
    
    int excTimeLimit;
    bool _isTimeLimited;
    
    static vector<vector<int>> computeAssignmentMatrix(VarValues);
    static void printAssignmentMatrix(vector<vector<int>>);

    bool greedy_heuristic_betapositive(bool is_cycle = false);
    bool greedy_heuristic_betazero(bool is_cycle = false);
    bool _greedy_heuristic_betazero_single_assingment(int t,
                                                      vector<int> & previousAssociation,
                                                      vector<int> & first_time_association,
                                                      double & objFun,
                                                      VarValues & greedyVarValues,
                                                      bool timeforward,
                                                      bool is_cycle=false);
    
    bool fileExists(ifstream &, string);
    void readParameterFiles(ifstream &);
    void readDistanceFile(ifstream &);
    void readDemandsFile(ifstream &);
    

};

#endif /* InputDataApCloudletAssoc_hpp */
