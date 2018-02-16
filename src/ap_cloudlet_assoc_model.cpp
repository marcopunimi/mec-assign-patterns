//
//  ap_cloudlet_assoc_model.cpp
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

#include "ap_cloudlet_assoc_model.hpp"

Ap_cloudlet_assoc_model::~Ap_cloudlet_assoc_model(){};

void Ap_cloudlet_assoc_model::solve(IloCplex cplex)
{
    solveCplexAsIs(cplex);
}

void Ap_cloudlet_assoc_model::solveCplexAsIs(IloCplex cplex)
{
    double optimalSol;
    int i,t,k1,k2;
    time_t excTimeInit;
    std::time(&excTimeInit);
    bool isCplexSolved = cplex.solve();
    double excTime = computeTime(excTimeInit);
    
    if(isCplexSolved)
    {
        if(cplex.getStatus() != IloAlgorithm::Infeasible && cplex.getCplexStatus() != IloCplex::Infeasible)
        {
            optimalSol =  cplex.getObjValue();
            if(!this->isBinary)
                cout << "CONTINUOUS ";
            cout << "cplex ok with o.f. value " << optimalSol << endl;
            if(this->isBinary){
                cout << "best bound open nodes " << cplex.getBestObjValue() << endl;
                cout << "last MIP GAP " << cplex.getMIPRelativeGap() << endl;
                cout << "number of nodes " << cplex.getNnodes() << endl;
            }
            cout << "cplex time " << excTime  << endl;
            
            if(!inputData->single_assignment){
                // save variables value
                for(i = 0; i < inputData->apCard; i++){
                    for(k1 = 0; k1< inputData->cloudletCard;k1++){
                        for(t = 0; t< inputData->timeSlotCard; t++){
                            varValues.xval[i][k1][t] = cplex.getValue(x[i][k1][t]);
                        }
                    }
                }
                
                // save variables value
                for(i = 0; i < inputData->apCard; i++){
                    for(k1 = 0; k1< inputData->cloudletCard;k1++){
                        for(k2 = 0; k2< inputData->cloudletCard;k2++) if(k1!=k2){
                            for(t = 0; t< inputData->timeSlotCard; t++){
                                varValues.yval[i][k1][k2][t] = cplex.getValue(y[i][k1][k2][t]);
                            }
                        }
                    }
                }
            }
            else{
                for(i = 0; i < inputData->apCard; i++){
                    for(k1 = 0; k1< inputData->cloudletCard;k1++){
                        for(t = 0; t< inputData->timeSlotCard; t++){
                            varValues.xval[i][k1][t] = cplex.getValue(x_fixed[i][k1]);
                        }
                    }
                }
                for(i = 0; i < inputData->apCard; i++){
                    for(k1 = 0; k1< inputData->cloudletCard;k1++){
                        for(k2 = 0; k2< inputData->cloudletCard;k2++) if(k1!=k2){
                            for(t = 0; t< inputData->timeSlotCard; t++){
                                varValues.yval[i][k1][k2][t] = 0.0;
                            }
                        }
                    }
                }
            }
                
            Solution sol;
            if(!this->isBinary){
                // otherwise, from the current fraction PB solution compute a new integer solution via rounding
                //Solution currentIntegerSol;
                
                sol = inputData->getIntegerSolutionViaRounding(varValues);
                varValues = sol.variableValues;
            }
            else{
                sol = Solution(inputData->apCard,
                             inputData->cloudletCard,
                             inputData->timeSlotCard,
                             cyclicModel);
            }
            double assignCost = 0;
            double switchingCost = 0;
            inputData->getAssignmentAndSwitchingCostVarValues(assignCost, switchingCost, varValues, cyclicModel);
            cout << "assignCost cost " << assignCost << " switching cost " << switchingCost << endl;
            
            printVarValues();
            
            if(this->isBinary){
                sol.finalCost = optimalSol;
                sol.isBinary = true;
                sol.variableValues = varValues;
            }
            sol.assignmentCost = assignCost;
            sol.switchingCost = switchingCost;
            
            inputData->setBestIntegerSolution(sol);
            inputData->printBestIntegerSolution();
            inputData->printAssignmentMatrixBestIntegerSolution();
        }
        else{
            cout << "cplex status " << cplex.getStatus() << " " << cplex.getCplexStatus() << endl;
            cerr << "cplex status " << cplex.getStatus() << " " << cplex.getCplexStatus() << endl;
        }
    }
    else{
        cout << "cplex status " << cplex.getStatus() << " " << cplex.getCplexStatus() << endl;
        cerr << "cplex status " << cplex.getStatus() << " " << cplex.getCplexStatus() << endl;
    }
}

void Ap_cloudlet_assoc_model::solveBinaryHeuristic(IloCplex cplex)
{
    double optimalSol;
    int i,t,k1;
    int round = 0;
    double xval = 0.0;
    
    // keep track of already fixed AP association in each time-slot (each AP is connected to 1 and 1 only cloudlet in each time-slot)
    vector< vector<bool> > is_fixed_ap(inputData->apCard, vector<bool>(inputData->timeSlotCard, false));
    // keep track of cloudlet usage in every time-slot given by fixed AP associations
    vector< vector<double> > current_usage(inputData->cloudletCard, vector<double>(inputData->timeSlotCard,0.0));
    
    while(true){
        vector<ApCloudletAssoc> varToFix;
        
        bool isCplexSolved = cplex.solve();
        if(isCplexSolved){
            if(cplex.getStatus() != IloAlgorithm::Infeasible && cplex.getCplexStatus() != IloCplex::Infeasible)
            {
                optimalSol =  cplex.getObjValue();
                cout << "round " << round <<  " cplex ok with o.f. value " << optimalSol << endl;
                
                // fix x variable with solution value 1, if not already fixed
                for(i = 0; i < inputData->apCard; i++)  {
                    for(t = 0; t< inputData->timeSlotCard; t++) if(!is_fixed_ap[i][t]) {
                        for(k1 = 0; k1< inputData->cloudletCard;k1++){
                            xval = cplex.getValue(x[i][k1][t]);
                            
                            // if value approach 1, fix x variable value
                            if(xval > 1-1e-6){
                                
                                varToFix.push_back(ApCloudletAssoc(i,k1,t));
                                
                                is_fixed_ap[i][t] = true;
                                
                                current_usage[k1][t] += inputData->nodeDemand[i][t];
                                cout << "fixed var " << i << "-" << k1 << "-" << t << endl;
                                break;
                            }
                        }
                    }
                }
                
                ApCloudletAssoc fractional_to_round = pickFractionalXvarHeuristic(&cplex, is_fixed_ap, current_usage);
                if(fractional_to_round.apId == -1){
                    break;
                }
                else{
                    x[fractional_to_round.apId][fractional_to_round.cloudletId][fractional_to_round.timeslotId].setLB(1);
                    is_fixed_ap[fractional_to_round.apId][fractional_to_round.timeslotId] = true;
                    current_usage[fractional_to_round.cloudletId][fractional_to_round.timeslotId] +=
                    inputData->nodeDemand[fractional_to_round.apId][fractional_to_round.timeslotId];
                    cout << "fixed var " << fractional_to_round.apId << "-" << fractional_to_round.cloudletId << "-" << fractional_to_round.timeslotId << endl;
                }
                
                for(i = 0; i < varToFix.size(); i++){
                    x[varToFix[i].apId][varToFix[i].cloudletId][varToFix[i].timeslotId].setLB(1);
                }
                
            }
            else{
                cout << "cplex status " << cplex.getStatus() << " " << cplex.getCplexStatus() << endl;
                cerr << "cplex status " << cplex.getStatus() << " " << cplex.getCplexStatus() << endl;
                break;
            }
        }
        else{
            cout << "cplex status " << cplex.getStatus() << " " << cplex.getCplexStatus() << endl;
            cerr << "cplex status " << cplex.getStatus() << " " << cplex.getCplexStatus() << endl;
            break;
        }
        round++;
    }
    addIntegralityConstraints();
    solveCplexAsIs(cplex);
}

ApCloudletAssoc Ap_cloudlet_assoc_model::pickFractionalXvarHeuristic(IloCplex* cplex,vector< vector< bool > > is_fixed_ap, vector< vector< double > > current_cloudlet_usage)
{
    
    int i,k1,t;
    
    ApCloudletAssoc highest_fract = ApCloudletAssoc();
    // highest_fract.apId = -1;
    // highest_fract.cloudletID = -1;
    // highest_fract.timeslotID = -1;
    double highest_fract_value = 0.0;
    for(i = 0; i < inputData->apCard; i++) {
        for(t = 0; t< inputData->timeSlotCard; t++) if(!is_fixed_ap[i][t]) {
            for(k1 = 0; k1< inputData->cloudletCard;k1++) if( current_cloudlet_usage[k1][t] + inputData->nodeDemand[i][t] <= inputData->U*inputData->cloudletCap){
                double current_value = cplex->getValue(x[i][k1][t]);
                if(IS_FRACTIONAL(current_value) && current_value > highest_fract_value)
                {
                    highest_fract_value = current_value;
                    
                    highest_fract.apId = i;
                    highest_fract.cloudletId = k1;
                    highest_fract.timeslotId = t;
                }
            }
        }
    }
    
    return highest_fract;
}

void Ap_cloudlet_assoc_model::addIntegralityConstraintsX()
{
    IloEnv env = this->model.getEnv();
    int nA = this->inputData->apCard;
    int nK = this->inputData->cloudletCard;
    int i,k1;
    
    // convert x and y variables to real domain
    IloArray< IloArray< IloConversion > > conversionX(env, nA);
    for(i = 0; i < nA; i++){
        conversionX[i] = IloArray< IloConversion >(env, nK);
        for(k1 = 0; k1 <nK; k1++){
            conversionX[i][k1] = IloConversion(env, x[i][k1], ILOINT);
            model.add(conversionX[i][k1]);
        }
    }

}

void Ap_cloudlet_assoc_model::addIntegralityConstraintsY()
{
    IloEnv env = this->model.getEnv();
    int nA = this->inputData->apCard;
    int nK = this->inputData->cloudletCard;
    int i,k1,k2;
    
    IloArray< IloArray< IloArray< IloConversion > > > conversionY(env, nA);
    for(i =0;i<nA;i++){
        conversionY[i] = IloArray< IloArray< IloConversion > >(env,nK);
        for(k1=0;k1<nK;k1++){
            conversionY[i][k1] = IloArray< IloConversion >(env, nK);
            for(k2=0;k2<nK;k2++){
                conversionY[i][k1][k2] = IloConversion(env, y[i][k1][k2], ILOINT);
                model.add(conversionY[i][k1][k2]);
            }
        }
    }
}

void Ap_cloudlet_assoc_model::addIntegralityConstraints()
{
    addIntegralityConstraintsX();
    addIntegralityConstraintsY();
}

void Ap_cloudlet_assoc_model::printVarValues()
{
    InputDataApCloudletAssoc::printVarValues(this->varValues);
}

void Ap_cloudlet_assoc_model::create_Model(IloEnv env,ModelVariant variant, bool isCyclic)
{
    int i,k1,k2,t;
    
    int nA = this->inputData->apCard;
    int nT = this-> inputData->timeSlotCard;
    int nK = this->inputData->cloudletCard;
    
    // x variables for each triplet AP, cloudlet, time-slot
    x = IloArray< IloArray< IloNumVarArray > >(env, nA);
    for(i =0;i<nA;i++)
    {
        x[i] = IloArray< IloNumVarArray >(env,nK);
        for(k1 = 0; k1 < nK; k1++)
        {
            x[i][k1] = IloNumVarArray(env, nT);
            for(t = 0; t < nT; t++)
            {
                char name[100] = "";
                snprintf(name, 100, "x_%d_%d_%d", i, k1, t );
                x[i][k1][t] = IloNumVar(env,0,1,name);
            }
            model.add(x[i][k1]);
        }
    }
    
    // y variables for each quadruple AP, cloudlet 1, cloudlet 2, time-slot
    y = IloArray< IloArray< IloArray< IloNumVarArray > > >(env, nA);
    for(i =0;i<nA;i++)
    {
        y[i] = IloArray< IloArray< IloNumVarArray > >(env,nK);
        for(k1=0; k1 < nK; k1++)
        {
            y[i][k1] = IloArray< IloNumVarArray >(env, nK);
            for(k2=0; k2 < nK; k2++)
            {
                y[i][k1][k2] = IloNumVarArray(env, nT);
                for(t=0 ; t < nT; t++)
                {
                    char name[100] = "";
                    snprintf(name, 100, "y_%d_%d_%d_%d", i ,k1, k2, t );
                    y[i][k1][k2][t] = IloNumVar(env,0,1,name);
                }
                model.add(y[i][k1][k2]);
            }
        }
    }
    
    // CONSTRAINTS
    //minimize objMinCost: sum{i in A, t in T, k1 in K} ( beta*d[i,t]*m[i,k1]*x[i,k1,t] + alpha*sum{k2 in K} (d[i,t]*l[k1,k2]*y[i,k1,k2,t]) );
    IloExpr expr(env);
    for(i = 0; i < nA; i++)
    {
        for(t = 0; t < nT; t++)
        {
            for(k1=0; k1 < nK; k1++)
            {
                expr += inputData->beta* inputData->nodeDemand[i][t] * inputData->distanceApCloudlet[i][k1]* x[i][k1][t];
                for(k2=0; k2 < nK; k2++)
                {
                    expr += inputData->alpha * inputData->nodeDemand[i][t] * inputData->distanceCloudletCloudlet[k1][k2]* y[i][k1][k2][t];
                }
            }
        }
    }
    
    minimizeCost = IloMinimize(env, expr, "minimizeCost");
    model.add(minimizeCost);
    
    //s.t. associateOne{i in A, t in T}: sum{ k in K } x[i,k,t] = 1;
    
    associateOne = IloArray< IloRangeArray >(env, nA);
    // for all APs
    for(i = 0; i < nA; i++)
    {
        associateOne[i] = IloRangeArray(env, nT);
        // for all time-slots
        for(t = 0; t < nT; t++)
        {
            IloExpr expr(env);
            // sum all x variables of all cloudlet
            for(k1=0; k1 < nK; k1++)
            {
                expr += x[i][k1][t];
            }
            char name[100] = "";
            snprintf(name, 100, "associateOne_%d_%d",i, t );
            associateOne[i][t] = IloRange(env, 1, expr, 1, name );
        }
        model.add(associateOne[i]);
    }
    
    //s.t. cloudletCap{t in T, k in K}: sum{i in A} d[i,t]*x[i,k,t] <= U*C;
    cloudletMaxCap = IloArray< IloRangeArray >(env, nT);
    // for all time slots
    for(t = 0; t < nT; t++){
        cloudletMaxCap[t] = IloRangeArray(env, nK);
        // for all cloudlets
        for(k1=0; k1 < nK; k1++){
            IloExpr expr(env);
            // sum all demands of APs associated to cloudlet
            for(i=0; i < nA; i++){
                expr += inputData->nodeDemand[i][t]*x[i][k1][t];
            }
            char name[100] = "";
            snprintf(name, 100, "cloudletMaxCap_%d_%d",t,k1 );
            cloudletMaxCap[t][k1] = IloRange(env, -IloInfinity, expr, inputData->U*inputData->cloudletCap, name );
        }
        model.add(cloudletMaxCap[t]);
    }
    
    
    addLinkingCostraintsBinaryVsCountinous();
    
    for(i =0;i<nA;i++)
    {
        for(k1=0; k1 < nK; k1++)
        {
            if(!inputData->isApCloudletPossible[i][k1]){
                for(t = 0; t < nT; t++){
                    x[i][k1][t].setUB(0.0);
                }
            }
        }
    }
}

void Ap_cloudlet_assoc_model::create_model_single_assign(IloEnv env,ModelVariant variant, bool isCyclic)
{
    int nA = this->inputData->apCard;
    int nT = this-> inputData->timeSlotCard;
    int nK = this->inputData->cloudletCard;
    
    x_fixed = IloArray< IloIntVarArray >(env,nA);
    for(int i =0;i<nA;i++)
    {
        x_fixed[i] = IloIntVarArray(env,nK);
        for(int k1 = 0; k1 < nK; k1++)
        {
            char name[100] = "";
            snprintf(name, 100, "xf_%d_%d", i, k1 );
            x_fixed[i][k1] = IloIntVar(env,0,1,name);
            
        }
        model.add(x_fixed[i]);
    }
    
    IloExpr expr(env);
    for(int i = 0; i < nA; i++)
    {
        for(int t = 0; t < nT; t++)
        {
            for(int k1=0; k1 < nK; k1++)
            {
                expr += inputData->beta* inputData->nodeDemand[i][t] * inputData->distanceApCloudlet[i][k1]* x_fixed[i][k1];
            }
        }
    }
    minimizeCost = IloMinimize(env, expr, "minimizeCost");
    model.add(minimizeCost);
    
    associateOneAllTime = IloRangeArray(env,nA);
    for(int i =0;i<nA;i++)
    {
        IloExpr tempExpr(env);
        for(int k1 = 0; k1 < nK; k1++)
        {
            tempExpr += x_fixed[i][k1];
        }
        char name[100] = "";
        snprintf(name, 100, "associateOneAllTime_%d",i );
        associateOneAllTime[i] = IloRange(env, 1, tempExpr, 1, name );
    }
    model.add(associateOneAllTime);
    
    cloudletMaxCap = IloArray< IloRangeArray >(env, nT);
    // for all time slots
    for(int t = 0; t < nT; t++)
    {
        cloudletMaxCap[t] = IloRangeArray(env, nK);
        // for all cloudlets
        for(int k1=0; k1 < nK; k1++)
        {
            IloExpr expr(env);
            // sum all demands of APs associated to cloudlet
            for(int i=0; i < nA; i++)
            {
                expr += inputData->nodeDemand[i][t]*x_fixed[i][k1];
            }
            char name[100] = "";
            snprintf(name, 100, "cloudletMaxCap_%d_%d",t,k1 );
            cloudletMaxCap[t][k1] = IloRange(env, -IloInfinity, expr, inputData->U*inputData->cloudletCap, name );
        }
        model.add(cloudletMaxCap[t]);
    }
    
}

Ap_cloudlet_assoc_model::Ap_cloudlet_assoc_model(IloEnv env, InputDataApCloudletAssoc* inputData, ModelVariant variant, bool isCyclic)
{
    this->cyclicModel = isCyclic;
    this->isBinary = (variant==BinaryCyclic || variant == Binary);
    this->variant = variant;
    this->inputData = inputData;
    this->model = IloModel(env);
    
    int nA = this->inputData->apCard;
    int nT = this-> inputData->timeSlotCard;
    int nK = this->inputData->cloudletCard;
    
    
    varValues.xval = vector< vector< vector< double> > >(nA, vector< vector< double> >(nK, vector< double>(nT, 0.0)));
    varValues.yval = vector< vector< vector< vector<double> > > >(nA, vector< vector< vector< double> > >(nK, vector< vector< double> >(nK, vector< double>(nT, 0.0))));
    
    if(!inputData->single_assignment)
        create_Model(env, variant, isCyclic);
    else
        create_model_single_assign(env, variant, isCyclic);
};

void Ap_cloudlet_assoc_model::addLinkingCostraintsBinaryVsCountinous()
{    
    IloEnv env = model.getEnv();
    int nA = inputData->apCard;
    int nT = inputData->timeSlotCard;
    int nK = inputData->cloudletCard;
    int i,k1,k2,t;
    
    
    //s.t. countLink1{i in A, t in T, k1 in K }: x[i,k1,t] = sum{k2 in K} y[i,k2,k1,t];
    countLink1 = IloArray< IloArray< IloRangeArray > >(env, nA);
    for(i = 0; i < nA; i++){
        int ntimes = nT - 1;
        if(this->cyclicModel)
        {
            ntimes = nT;
        }
        countLink1[i] =IloArray< IloRangeArray >(env, ntimes);
        
        for(t = 0; t < ntimes; t++)
        {
            int tt = t + 1;
            if(this->cyclicModel)
            {
                tt = t;
            }
            
            countLink1[i][t] = IloRangeArray(env, nK);
            for(k1 = 0; k1 < nK; k1++ ){
                char name[100] = "";
                snprintf(name, 100, "countLink1_%d_%d_%d",i, tt, k1 );
                IloExpr expr(env);
                for(k2 = 0; k2 < nK; k2++ )
                {
                    expr += y[i][k2][k1][tt];
                }
                countLink1[i][t][k1] = IloRange(env,0, x[i][k1][tt] - expr ,0,name);
            }
            model.add(countLink1[i][t]);
        }
    }
    
    //s.t. countLink2{i in A, t in T, k1 in K : t < Nt }: x[i,k1,t] = sum{k2 in K} y[i,k1,k2,t+1];
    countLink2 = IloArray< IloArray< IloRangeArray > >(env, nA);
    for(i = 0; i < nA; i++)
    {
        int ntimes = nT - 1;
        if(this->cyclicModel)
        {
            ntimes = nT;
        }
        countLink2[i] = IloArray< IloRangeArray >(env, ntimes);
        
        for(t = 0; t < ntimes; t++)
        {
            int nextTime = t + 1;
            if(this->cyclicModel && nextTime == nT)
            {
                nextTime = 0;
            }
            countLink2[i][t] = IloRangeArray(env, nK);
            for(k1 = 0; k1 < nK; k1++ )
            {
                
                char name[100] = "";
                snprintf(name, 100, "countLink2_%d_%d_%d",i, t, k1 );
                IloExpr expr(env);
                for(k2 = 0; k2 < nK; k2++ )
                {
                    expr += y[i][k1][k2][nextTime];
                }
                countLink2[i][t][k1] = IloRange(env,0, x[i][k1][t] - expr ,0,name);
            }
            model.add(countLink2[i][t]);
        }
        
    }
    
    if(this->isBinary)
    {
        addIntegralityConstraintsX();
    }

}
