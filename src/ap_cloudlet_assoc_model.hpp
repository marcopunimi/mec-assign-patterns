//
//  ap_cloudlet_assoc_model.hpp
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
