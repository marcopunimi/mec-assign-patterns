//
//  ap_cloudlet_assoc_model.hpp
//  ap_cloudlet_association
//
//  Created by Marco Premoli on 17/01/2017.
//  Copyright © 2017 Marco Premoli. All rights reserved.
//
#pragma once
#ifndef ap_cloudlet_assoc_model_hpp
#define ap_cloudlet_assoc_model_hpp

#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h>
#include "InputDataApCloudletAssoc.hpp"
// #include "supportFunctions.hpp"

struct ApCloudletAssoc{
    int apId;
    int cloudletId;
    int timeslotId;
    
    ApCloudletAssoc():
    apId(-1), cloudletId(-1), timeslotId(-1) {};
    
    ApCloudletAssoc(int apId, int cloudletId, int timeslotId):
    apId(apId), cloudletId(cloudletId), timeslotId(timeslotId) {};
    
};

class Ap_cloudlet_assoc_model{
public:
    Ap_cloudlet_assoc_model(IloEnv, InputDataApCloudletAssoc*, ModelVariant variant, bool isCyclic = false);
    ~Ap_cloudlet_assoc_model();
    
    void solve(IloCplex);
    
    IloModel model;
private:
    
    void create_Model(IloEnv,ModelVariant variant, bool isCyclic = false);
    void create_model_single_assign(IloEnv,ModelVariant variant, bool isCyclic = false);
    
    InputDataApCloudletAssoc* inputData;
    VarValues varValues;
    bool cyclicModel;
    bool isBinary;
    
    // VARIABLES
    // x variables for each triplet AP, cloudlet, time-slot
    IloArray< IloArray< IloNumVarArray > > x;
    
    // y variables for each quadruple AP, cloudlet 1, cloudlet 2, time-slot
    IloArray< IloArray< IloArray< IloNumVarArray > > > y;
    
    // CONSTRAINTS
    IloObjective minimizeCost;
    
    IloArray< IloRangeArray > associateOne;
    IloArray< IloRangeArray > cloudletMaxCap;
    IloArray< IloArray< IloRangeArray > > binaryLink1;
    IloArray< IloArray< IloRangeArray > > binaryLink2;
    
    IloArray< IloArray< IloArray< IloRangeArray > > > binaryLink_former;
    
    IloArray< IloArray< IloRangeArray > > countLink1;
    IloArray< IloArray< IloRangeArray > > countLink2;
    
    ModelVariant variant;
    
    void addLinkingCostraintsBinaryVsCountinous();
    
    ApCloudletAssoc pickFractionalXvarHeuristic(IloCplex* cplex, vector< vector< bool > >, vector< vector< double > >);
    
    void printVarValues();
    void addIntegralityConstraints();
    void addIntegralityConstraintsX();
    void addIntegralityConstraintsY();
    
    
    void solveBinaryHeuristic(IloCplex cplex);
    void solveCplexAsIs(IloCplex cplex);
    
    /* single assignment */
    IloArray< IloIntVarArray >  x_fixed;
    IloRangeArray associateOneAllTime;
    
    
};

#endif /* ap_cloudlet_assoc_model_hpp */
