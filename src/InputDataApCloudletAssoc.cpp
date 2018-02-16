//
//  InputDataApCloudletAssoc.cpp
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

#include "InputDataApCloudletAssoc.hpp"

bool InputDataApCloudletAssoc::isTimeLimited(){
    return _isTimeLimited;
}

double InputDataApCloudletAssoc::getCurrentExcTime()
{
    double exc_time = -1;
    //if(_isTimeLimited){
       exc_time = computeTime(this->initExecutionTime);
    //}
    return exc_time;
}

// return TRUE if the time limit of execution is reached
bool InputDataApCloudletAssoc::isExcTimeLimitReached()
{
    // if execution time is limited
    if(_isTimeLimited){
        // check if current execution time is higher than time limit
        double exc_time = computeTime(this->initExecutionTime);
        return exc_time > excTimeLimit;
    }
    return false;
}

// read and store integer solution from file
// solution file format is a matrix of integers, as assignment matrix:
//  A X T = one row per AP, one column per time-slot, containing the index of the cloudlet to which
// assign the AP in the time-slot (from 1 to |K|)
void InputDataApCloudletAssoc::readIntegerSolutionFromFile(string datafile)
{
    if(this->readOk)
    {
        ifstream file(datafile.c_str(), ios::in);
        string line = "";
        stringstream iss;
        int lineNb = 0;
        string entry = "";
        
        if (!file.is_open())
        {
            cout << "not existing file " << datafile.c_str() << endl;
            cerr << "not existing file " << datafile.c_str() << endl;
            return;
        }
        
        VarValues values(apCard, cloudletCard, timeSlotCard);
        
        for(int i = 0 ; i < apCard ; i++){
            getline(file,line);
            lineNb++;
            iss.clear();
            iss.str(line);
            int t = 0;
            while (getline(iss, entry, ',') && t < timeSlotCard)
            {
                int cloud_temp = atoi(entry.c_str());
                values.xval[i][cloud_temp-1][t] = 1.0;
                t++;
            }
        }
        this->computeYvarValuesGivenXvarValues(values);
        
        Solution newsol(apCard, cloudletCard, timeSlotCard);
        newsol.variableValues = values;
        newsol.finalCost = this->computeCostGivenVariableValues(values);
        newsol.isBinary = true;
        
        this->setBestIntegerSolution(newsol);
        
        file.close();
    }
}

bool InputDataApCloudletAssoc::fileExists(ifstream & filmToCheck,
                                          string filename){
    if (!filmToCheck.is_open())
    {
        cout << "not existing file " << filename << endl;
        cerr << "not existing file " << filename << endl;
        return false;
    }
    return true;
}

void InputDataApCloudletAssoc::readParameterFiles(ifstream & fileParam){
    string line = "";
    
    // alpha parameter
    getline(fileParam,line);
    alpha = atof(line.c_str());
    
    // beta parameter
    getline(fileParam,line);
    beta = atof(line.c_str());
    
    // U parameter
    getline(fileParam,line);
    U = atoi(line.c_str());
    
    getline(fileParam,line);
    cloudletCapPerc = atof(line.c_str());
    
}

void InputDataApCloudletAssoc::readDistanceFile(ifstream & fileDistances)
{
    string line = "";
    int lineNb = 0;
    string entry = "";
    stringstream iss;
    int i,j;
    
    // cardinality APs
    getline(fileDistances, line);
    apCard = atoi(line.c_str());
    
    // cardinality cloudlets
    getline(fileDistances, line);
    cloudletCard = atoi(line.c_str());
    
    // resize distance matrix
    distanceCloudletCloudlet.resize(cloudletCard,vector< double >(cloudletCard,0.0));
    
    // resize distance matrix
    distanceApCloudlet.resize(apCard,vector< double >(cloudletCard,0.0));
    isApCloudletPossible.resize(apCard,vector< bool >(cloudletCard, true));
    
    // for each AP keep the number of facilities it can reach given the threshold on assignment distance
    ap_number_cloudlet_reached.resize(apCard, cloudletCard);
    
    // read N x K distance matrix (distance APs - facilities)
    for(i = 0 ; i < apCard ; i++){
        getline(fileDistances, line);
        lineNb++;
        iss.clear();
        iss.str(line);
        j = 0;
        while (getline(iss, entry, ',') && j < cloudletCard){
            double demandRead = atof(entry.c_str());
            distanceApCloudlet[i][j] = demandRead;
            j++;
        }
    }
    
    // read K x K distance matrix (distance facility-facility)
    for(i = 0 ; i < cloudletCard ; i++)
    {
        getline(fileDistances, line);
        lineNb++;
        iss.clear();
        iss.str(line);
        j = 0;
        while (getline(iss, entry, ',') && j < cloudletCard)
        {
            double distanceRead = atof(entry.c_str());
            distanceCloudletCloudlet[i][j] = distanceRead;
            j++;
        }
    }
    
    // sort AP per the number of reachable facilities, in ascending order
    sorted_ap_number_cloudlet_reached.resize(apCard);
    sorted_ap_number_cloudlet_reached = sort_indexes(ap_number_cloudlet_reached, false);
    
    // for each AP, sort cloudlet indexes by ascending distance
    cloudletIdxSortedByApDistance.resize(apCard);
    for(i = 0 ; i < apCard ; i++){
        cloudletIdxSortedByApDistance[i] = sort_indexes(distanceApCloudlet[i], false);
    }
}

void InputDataApCloudletAssoc::readDemandsFile(ifstream & fileDemands){
    string line = "";
    int i,j,t;
    int lineNb = 0;
    string entry = "";
    stringstream iss;
    
    // cardinality time-slots
    getline(fileDemands, line);
    timeSlotCard = atoi(line.c_str());
    
    // maximum cloudlet capacity
    // getline(fileDemands, line);
    // cloudletCap = atof(line.c_str());
    
    // resize demand matrix
    nodeDemand.resize(apCard,vector< double >(timeSlotCard,0.0));
    timeDemand.resize(timeSlotCard,vector< double >(apCard,0.0));
    
    // read N x T demand matrix
    for(i = 0 ; i < apCard ; i++){
        getline(fileDemands, line);
        lineNb++;
        iss.clear();
        iss.str(line);
        j = 0;
        while (getline(iss, entry, ',') && j < timeSlotCard){
            double demandRead = atof(entry.c_str());
            nodeDemand[i][j] = demandRead;
            j++;
        }
    }
    // timeDemand is nodeDemand transposed
    double maxTimeDemand = 0;
    double maxTimeDemandTemp;
    for( t = 0; t < timeSlotCard;t++){
        maxTimeDemandTemp = 0;
        for(i = 0; i < apCard; i++){
            timeDemand[t][i] = nodeDemand[i][t];
            maxTimeDemandTemp += nodeDemand[i][t];
        }
        if(maxTimeDemandTemp > maxTimeDemand){
            maxTimeDemand = maxTimeDemandTemp;
        }
    }
    cloudletCap = ceil((maxTimeDemand/cloudletCard)*cloudletCapPerc);
    cout << "capacity computed " << cloudletCap << endl;
    
    // for each time slot, get the ordered AP indexes for descending demand in the time slot
    sortedApIdxDescendingDemand.resize(timeSlotCard);
    for( t = 0; t<timeSlotCard;t++){
        vector<double> demand_in_time = timeDemand[t];
        vector<size_t> sorted_ap_idx = sort_indexes(demand_in_time, true);
        sortedApIdxDescendingDemand[t] = sorted_ap_idx;
    }
}

InputDataApCloudletAssoc::~InputDataApCloudletAssoc(){};
InputDataApCloudletAssoc::InputDataApCloudletAssoc(string filenameParameters,
                                               string filenameDistances,
                                               string filenameDemands,
                                               int timelimit,
                                               bool betaZero,
                                               double limitAssignmentDistance)
{
    std::time(&this->initExecutionTime);
    
    this->excTimeLimit = -1;
    this->_isTimeLimited = false;
    if(timelimit > 0){
        this->excTimeLimit = timelimit;
        this->_isTimeLimited = true;
    }
    
    ifstream fileParam(filenameParameters.c_str(), ios::in);
    ifstream fileDistances(filenameDistances.c_str(), ios::in);
    ifstream fileDemands(filenameDemands.c_str(), ios::in);
    
    dummyCost = -1;
    apCard = cloudletCard = timeSlotCard = 0;
    
    if(!fileExists(fileParam, filenameParameters.c_str())
       || !fileExists(fileDistances, filenameDistances.c_str())
       || !fileExists(fileDemands, filenameDemands.c_str()))
    {
        readOk = false;
        return;
    }
    
    /*
     input file format:
         * PARAM FILE
     - 1 double = alpha parameter
     - 1 double = beta parameter
     - 1 double = U parameter
     - 1 double > 1 = cloudlet percentage extra capacity w.r.t. peak demand in time
         * DISTANCE FILE
     - 1 int = AP cardinality
     - 1 int = cloudlet cardinality
     - matrix A x K double = distance AP to cloudlet
     - matrix K x K double = distance cloudlet to cloudlet
         * DEMAND FILE
     - 1 int = time slot cardinality
     - matrix N x T double = demand nodes in time
     */
    
    readParameterFiles(fileParam);
    readDistanceFile(fileDistances);
    readDemandsFile(fileDemands);
 
    // if beta has to be set to zero
    // also set the threshold on the assignment distance
    isBetaZero = false;
    if(betaZero){
        this->isBetaZero = true;
        this->beta = 0;
        this->setLimitAssignmentDistance(limitAssignmentDistance);
    }
    
    fileParam.close();
    fileDemands.close();
    fileDistances.close();
    
    cg_root_LB = -1;
    min_active_branch_node_LB = -1;
    integralityGAP = -1;
    
    readOk = true;
    
    bestIntegerSolution = Solution(apCard, cloudletCard, timeSlotCard);
    
    single_assignment = false;
};

void InputDataApCloudletAssoc::setLimitAssignmentDistance(double limit)
{
    for(int i = 0; i < this->distanceApCloudlet.size(); i++)
    {
        for(int j = 0; j < this->distanceApCloudlet[i].size(); j++)
        {
            if(distanceApCloudlet[i][j] >= limit){
                isApCloudletPossible[i][j] = false;
                ap_number_cloudlet_reached[i]--;
            }
        }
    }
}

void InputDataApCloudletAssoc::printAllData(){
    if(!readOk){
        cout << "data not read" << endl;
        return;
    }
    int i,j;
    cout << " AP size = " << apCard << endl <<         // cardinality APs
            " cloudlet size " << cloudletCard << endl <<   // cardinality cloudlets
            " time-slot size " << timeSlotCard << endl <<   // cardinality time-slots
            " alpha " << alpha << endl <<       // alpha parameter
            " beta " << beta << endl <<       // beta parameter
            " U " << U << endl <<       // U parameter
            " cloudlet capacity " << cloudletCap << endl;   // maximum cloudlet capacity
    
    cout << "node demand " << endl;
    for(i = 0; i < apCard; i++){
        for(j = 0; j < timeSlotCard; j++){
            cout << nodeDemand[i][j] << " ";
        }
        cout << endl;
    }
    cout << "distance cloudlet - cloudlet " << endl;
    for(i = 0; i < cloudletCard; i++){
        for(j = 0; j < cloudletCard; j++){
            cout << distanceCloudletCloudlet[i][j] << " ";
        }
        cout << endl;
    }
    cout << "distance AP - cloudlet " << endl;
    // distance btw APs and cloudlets
    for(i = 0; i < apCard; i++){
        for(j = 0; j < cloudletCard; j++){
            cout << distanceApCloudlet[i][j] << " ";
        }
        cout << endl;
    }
}

// given a vector of data, return the vector of sorted indexes of the given vector
// by descending or ascending order
template <typename T>
vector<size_t> InputDataApCloudletAssoc::sort_indexes(const vector<T> &v, bool isDescending) {
    
    // initialize original index locations
    vector<size_t> idx(v.size());
    for (size_t i = 0; i != idx.size(); ++i)
        idx[i] = i;
    
    // sort indexes based on comparing values in v
    if(isDescending)
        sort(idx.begin(), idx.end(),
             [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});
    else
        sort(idx.begin(), idx.end(),
             [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
    return idx;
}

int InputDataApCloudletAssoc::nearestCloudletNotExceedingCapacity(vector<double> current_usage,
                                        vector<double> distances,
                                        double ap_demand,
                                        double maxCapacity,
                                        vector<bool> isAssignmentPossible)
{    
    int nearestCloudlet = -1;
    
    double nearestDistance = 1e10;
    
    for(int k = 0; k < distances.size(); k++) if(isAssignmentPossible[k])
    {
        if(distances[k] < nearestDistance && current_usage[k] + ap_demand <= maxCapacity){
            nearestCloudlet = k;
            nearestDistance = distances[k];
        }
    }
    
    return(nearestCloudlet);
}

int InputDataApCloudletAssoc::nearestCloudletNotExceedingCapacity(vector<double> current_usage,
                                                                  vector<double> distances,
                                                                  double ap_demand,
                                                                  double maxCapacity)
{
    int nearestCloudlet = -1;
    
    double nearestDistance = 1e10;
    
    for(int k = 0; k < distances.size(); k++)
    {
        if(distances[k] < nearestDistance && current_usage[k] + ap_demand <= maxCapacity){
            nearestCloudlet = k;
            nearestDistance = distances[k];
        }
    }
    
    return(nearestCloudlet);
}

// to use only for integer solutions!!
double InputDataApCloudletAssoc::computeCostGivenVariableValues(VarValues varValues, bool is_cycle)
{
    double cost = 0.0;
    
    int cloud_first_time = -1;
    int prev_cloud = -1;
    for(int i=0; i< apCard; i++){
        for(int t = 0; t < timeSlotCard; t++){
            for(int k = 0; k < cloudletCard; k++){
                if(varValues.xval[i][k][t] > 1e-6)
                {
                    cost += beta*nodeDemand[i][t]*distanceApCloudlet[i][k]*varValues.xval[i][k][t];
                    if(t == 0){
                        cloud_first_time = k;
                    }
                    else if(prev_cloud != k && prev_cloud!=-1){
                        cost += alpha*nodeDemand[i][t]*distanceCloudletCloudlet[prev_cloud][k];
                    }
                    if (t == timeSlotCard - 1 && is_cycle && cloud_first_time!=k){
                        cost += alpha*nodeDemand[i][0]*distanceCloudletCloudlet[k][cloud_first_time];
                    }
                    prev_cloud = k;
                }
            }
        }
    }
    
    return cost;
}

Solution InputDataApCloudletAssoc::getIntegerSolutionViaRounding(VarValues varValues, bool is_cycle)
{
    return _getIntegerSolutionViaRounding(varValues, is_cycle);
}

Solution InputDataApCloudletAssoc::getIntegerSolutionViaRounding(VarValues varValues, BranchFixedVariables fixedVar, bool is_cycle)
{
    // get integer solution via rounding, with already fixed variables that cannot be unfix

    // if fixed variable structure is empty, return the standard rounding heuristic
    if(fixedVar.size() == 0)
    {
        return _getIntegerSolutionViaRounding(varValues);
    }
    
    Solution roundSolution;
    roundSolution.isBinary = true;
    roundSolution.finalCost = -1;
    
    cout << "ROUNDING HEURISTIC WITH FIXED VARIABLES" << endl;
    
    // otherwise run the following heuristic
    // for each time slot
    //      if AP-TIME assignment was forced, keep it fix and update the cloudlet residual capacity
    //      for each AP with fractional assignment, compute the cost gap between the 1st and 2nd choice assignments (i.e. two topmost assignment per x variable values)
    //      if AP assignment was at integrality but not forced, the 2nd choice is the nearest cloudlet except the forbidden ones and the assigned one
    //      cost gap choice: (i) pure distance AP-cloudlet vs (ii) distance times AP demand
    //      sort AP by descending cost gap
    //      assign AP to cloudlet following cost gap order, updating the cloudlet residual capacity
    //      if 1st and 2nd choice are no more available, get the nearest cloudlet with enough capacity not forbidden by previous fixing
    for(int t = 0; t < timeSlotCard; t++)
    {
        // cloudlet residual capacity
        vector<double> residualCap(cloudletCard, U*cloudletCap);
        // keep track of the APs that have forced assignment
        vector<bool> apFixed(apCard, false);
        
        // check if AP-CLOUDLET-TIME assignment is forced by fixed variables
        for(int i = 0; i < apCard; i++)
        {
            for(int k = 0; k < cloudletCard; k++)
            {
                // assignment AP-TIME-CLOUDLET is forced -> force assignment in rounding
                if(fixedVar[i][t][k] == 1)
                {
                    // fix assignment
                    varValues.xval[i][k][t] = 1;
                    
                    // update cloudlet residual capacity
                    residualCap[k] -= nodeDemand[i][t];
                    apFixed[i] = true;
                    
                    break;
                }
            }
            if(apFixed[i])
            {
                // JUST TO BE SURE, FIX ALL OTHER VARIABLES TO 0
                for(int k=0; k < cloudletCard; k++)
                {
                    if(fixedVar[i][t][k] != 1)
                    {
                        varValues.xval[i][k][t] = 0;
                    }
                }
            }
        }
        
        // FOR EVERY FRACTIONALLY ASSIGN AP, COMPUTE COST GAP BETWEEN 1ST AND 2ND CHOICE
        vector<double> costGap;  // gap of the considered AP
        vector<int> apIdx;      // ID of the AP considered
        vector<int> firstCloud; // ID of the first cloudlet choose for the considered AP
        vector<int> secondCloud; // ID of the second cloudlet choose for the considered AP
        
        for(int i = 0; i < apCard; i++)
        {
            // check if current AP has already been fixed to a cloudlet by a forced assignment
            if(!apFixed[i])
            {
                double firstChoiceAssign = 0;
                double secondChoiceAssign = 0;
                int firstChoiceCloud = -1;
                int secondChoiceCloud = -1;
                // get first and second choice assignment, i.e. the cloudlet associated to two topmost x variable values
                for(int k = 0; k < cloudletCard; k++)
                {
                    // check if assignment is not forced nor forbidden
                    if(fixedVar[i][t][k] == 2)
                    {
                        if(varValues.xval[i][k][t] > firstChoiceAssign)
                        {
                            // current first choice become second choice
                            secondChoiceAssign = firstChoiceAssign;
                            secondChoiceCloud = firstChoiceCloud;
                            // current value becomes first choice
                            firstChoiceAssign = varValues.xval[i][k][t];
                            firstChoiceCloud = k;
                        }
                        else if(varValues.xval[i][k][t] > secondChoiceAssign)
                        {
                            // update second choice
                            secondChoiceCloud = k;
                            secondChoiceAssign = varValues.xval[i][k][t];
                        }
                    }
                }
                double currentCostGap = 0;
                // if first choice was not found an error has occured
                if(firstChoiceAssign == -1)
                {
                    cout << " ERROR IN ROUNDING HEURISTIC WITH FIXED VARIABLES " << endl;
                    return roundSolution;
                }
                
                // if a second choice cloud has not been found
                // get the nearest cloudlet not forbidden
                if(secondChoiceCloud == -1 )
                {
                    for(int k = 0; k < cloudletCard; k++)
                    {
                        int currentCloud = (int) cloudletIdxSortedByApDistance[i][k];
                        if(currentCloud != firstChoiceCloud && fixedVar[i][t][currentCloud] == 2)
                        {
                            secondChoiceCloud = currentCloud;
                            break;
                        }
                    }
                }
                
                // check again if a second choice was found -> i.e. only a cloudlet assignment is possible but not forced
                if(secondChoiceCloud == -1 )
                {
                    secondChoiceCloud = firstChoiceCloud;
                    currentCostGap = 0;
                }
                else
                    currentCostGap = nodeDemand[i][t]*distanceApCloudlet[i][firstChoiceCloud] - nodeDemand[i][t]*distanceApCloudlet[i][secondChoiceCloud];
                
                // save cost gap and corresponding AP index
                firstCloud.push_back(firstChoiceCloud);
                secondCloud.push_back(secondChoiceCloud);
                costGap.push_back(currentCostGap);
                apIdx.push_back(i);
            }
        }
        
        if(costGap.size() > 0)
        {
            // SORT AP BY DESCENDING COST GAP
            vector<size_t> sortedAP = sort_indexes(costGap, true);
            for(int i = 0; i < sortedAP.size(); i++)
            {
                // get current AP
                int currentAP = apIdx[(int) sortedAP[i]];
                int currentFirst = firstCloud[(int) sortedAP[i]];
                int currentSecond = secondCloud[(int) sortedAP[i]];
                
                // get final assigned cloudlet
                int assignedCloud = -1;
                
                // try to assign AP to its first choice
                if(nodeDemand[currentAP][t] <= residualCap[currentFirst])
                {
                    assignedCloud = currentFirst;
                }
                else if(nodeDemand[currentAP][t] <= residualCap[currentSecond])
                {
                    assignedCloud = currentSecond;
                }
                else
                {
                    // first and second choice are unavailable, as third choice get the nearest uvailable and not-fixed cloudlet
                    for(int k = 0; k < cloudletCard; k++)
                    {
                        double currentCloud = (int) cloudletIdxSortedByApDistance[currentAP][k];
                        if(fixedVar[currentAP][t][currentCloud] == 2 &&
                           nodeDemand[currentAP][t] <= residualCap[currentCloud])
                        {
                            assignedCloud = currentCloud;
                            break;
                        }
                    }
                }
                
                // assign and update cloudlet residual capacity
                if(assignedCloud != -1)
                {
                    residualCap[assignedCloud] -= nodeDemand[currentAP][t];
                    varValues.xval[currentAP][assignedCloud][t] = 1.0;
                    // fix all other cloudlet assignment to 0
                    for(int k = 0; k < cloudletCard;k++) if(k != assignedCloud)
                    {
                            varValues.xval[currentAP][k][t] = 0.0;
                    }
                }
                else
                {
                    cout << "impossible to reach integer solution" << endl;
                    cerr << "impossible to reach integer solution" << endl;
                    return roundSolution;
                }
            }
        }
    }
    // COMPUTE Y VARIABLES VALUES A POSTERIORI
    this->computeYvarValuesGivenXvarValues(varValues);
    
    // compute total cost a posteriori
    roundSolution.finalCost = computeCostGivenVariableValues(varValues, is_cycle);
    
    // print integer cost and solution
    cout << "rounded integer cost " << roundSolution.finalCost << endl;
    cerr << "rounded integer cost " << roundSolution.finalCost << endl;
    
    // InputDataApCloudletAssoc::printVarValues(varValues);
    
    roundSolution.variableValues = varValues;
    
    checkFeasibilitySolution(roundSolution.variableValues);
    
    return roundSolution;
}

bool InputDataApCloudletAssoc::checkFeasibilitySolution(VarValues varValues)
{
    bool results = true;
    
    vector<vector<double>> cloudletCapacity(cloudletCard, vector<double>( timeSlotCard ,U*cloudletCap));
    
    // for the feasibility of the solution, check if cloudlet capacity are respected
    for(int i = 0; i < apCard; i++)
    {
        for(int t = 0; t < timeSlotCard; t++)
        {
            bool isApAssigned = false;
            for(int k = 0 ; k < cloudletCard; k++)
            {
                if(varValues.xval[i][k][t] > 1e-6)
                {
                    cloudletCapacity[k][t] -= nodeDemand[i][t];
                    
                    if(cloudletCapacity[k][t] < 0)
                    {
                        cout << "ERROR: CLOUDLET CAPACITY EXCEEDED " << k << " TIME " << t << " = " << cloudletCapacity[k][t] << endl;
                        cerr << "ERROR: CLOUDLET CAPACITY EXCEEDED " << k << " TIME " << t << " = " << cloudletCapacity[k][t] << endl;
                        return false;
                    }
                    isApAssigned = true;
                }
            }
            if(!isApAssigned)
            {
                cout << "ERROR: AP " << i << " NOT ASSIGNED AT TIME " << t << endl;
                cerr << "ERROR: AP " << i << " NOT ASSIGNED AT TIME " << t << endl;
                return false;
            }
        }
    }
    
    return results;
}

void InputDataApCloudletAssoc::computeYvarValuesGivenXvarValues(VarValues & varValues)
{
    // re-initialize y var. values
    varValues.yval = vector< vector< vector< vector<double> > > >(apCard, vector< vector< vector< double> > >(cloudletCard, vector< vector< double> >(cloudletCard, vector< double>(timeSlotCard, 0.0))));
    
    for(int i = 0; i < apCard; i++)
    {
        int cloudPrevTime = -1;
        // get cloudlet associated first time slot
        for(int k=0; k< cloudletCard;k++)
        {
            if(varValues.xval[i][k][0] > 1.0-1e-6)
            {
                cloudPrevTime = k;
                break;
            }
        }
        
        for(int t = 1; t < timeSlotCard;t++)
        {
            int cloudTime = -1;
            
            // get cloudlet associated current time slot
            for(int k=0; k< cloudletCard;k++)
            {
                if(varValues.xval[i][k][t] > 1.0-1e-6)
                {
                    cloudTime = k;
                    break;
                }
            }
            // if cloudlets differs, give value 1 to corresponding variable y
            if(cloudPrevTime!=cloudTime)
            {
                varValues.yval[i][cloudPrevTime][cloudTime][t] = 1.0;
            }
            // update cloudlet of previous time slot to continue looping
            cloudPrevTime = cloudTime;
        }
    }
}

Solution InputDataApCloudletAssoc::_getIntegerSolutionViaRounding(VarValues varValues, bool is_cycle)
{
    // given current fractional solution round one variable at a time to get integer solution
    // for each time slot
    //      for each AP
    //          compute highest x_{i,k}^t value
    //      sort AP by descending highest x value
    //      for each sorted AP
    //          associate it to cloudlet corresponding with  highest x value with enough residual capacity
    //          update cloudlet residual capacity
    // compute total cost a posteriori
    
    cout << "ROUNDING HEURISTIC" << endl;
    
    Solution integerSolution;
    integerSolution.isBinary = true;
    integerSolution.finalCost = -1;
    integerSolution.is_cycle = is_cycle;
    
    // for each time slot
    for(int t=0; t < timeSlotCard; t++)
    {
        // for each AP compute highest x_{i,k}^t value
        vector<double> maxX(apCard,0.0);
        vector<int> cloudMaxX(apCard,-1);
        
        // cloudlet residual capacity
        vector<double> cloudCap(cloudletCard, U*cloudletCap);
        
        for(int i =0; i < apCard; i++)
        {
            for(int k = 0; k < cloudletCard; k++) if(isApCloudletPossible[i][k])
            {
                if(varValues.xval[i][k][t] > maxX[i])
                {
                    maxX[i] = varValues.xval[i][k][t];
                    cloudMaxX[i] = k;
                    // if maximum value is greater than 0.5 the loop can end (the sum of all x variables in a time is 1)
                    if(maxX[i] > 0.5)
                        break;
                }
            }
        }
        
        int splitCoef = 3;
        int halfApIndex = floor(apCard/splitCoef);
        
        for(int j = 0; j < splitCoef; j++)
        {
            vector<double> maxXhalf;
            vector<int> cloudMaxXhalf;
           
            int startIter = halfApIndex*j;
            int endIter = halfApIndex*(j+1);
            if(j+1 == splitCoef){
                endIter = apCard;
            }
            for(int i = startIter; i < endIter; i++)
            {
                int currentAP = -1;
                if(this->isBetaZero)
                    currentAP = (int) sorted_ap_number_cloudlet_reached[i];
                else
                    currentAP = (int) sortedApIdxDescendingDemand[t][i];
                maxXhalf.push_back(maxX[ currentAP ]);
                cloudMaxXhalf.push_back(cloudMaxX[ currentAP ]);
            }
            
            if(maxXhalf.empty())
                break;
            
            // sort AP by descending highest x value
            vector<size_t> sorted_ap_idx = InputDataApCloudletAssoc::sort_indexes(maxXhalf, true);
            
            // sort AP by descending highest x value
            
            // for each sorted AP
            for(int i = 0; i < sorted_ap_idx.size(); i++)
            {
                int current_ap = -1 ;
                if(isBetaZero)
                    current_ap = (int) sorted_ap_number_cloudlet_reached[(int) sorted_ap_idx[i] + j*halfApIndex];
                else
                    current_ap = (int) sortedApIdxDescendingDemand[t][(int) sorted_ap_idx[i] + j*halfApIndex];
                
                // associate it to cloudlet corresponding with highest x value with enough residual capacity
                double currentMaxX = -1.0;
                int currentCloudMaxX = -1;
                
                double currentMaxXhalf = maxXhalf[(int) sorted_ap_idx[i]];
                int currentCloudMaxXhalf = cloudMaxXhalf[(int) sorted_ap_idx[i]];
                if( nodeDemand[current_ap][t] <= cloudCap[currentCloudMaxXhalf])
                {
                    currentMaxX = currentMaxXhalf;
                    currentCloudMaxX = currentCloudMaxXhalf;
                }
                
                if(currentCloudMaxX==-1)
                {
                    // cloudlet corresponding to maximum x variable value is not available
                    // explore all other cloudlets in descending distance to the AP order
                    for(int k=0; k< cloudletCard;k++)
                    {
                        double currentCloud = (int) cloudletIdxSortedByApDistance[current_ap][k];
                        if(nodeDemand[current_ap][t] <= cloudCap[currentCloud]
                           && varValues.xval[current_ap][currentCloud][t] > currentMaxX
                           && this->isApCloudletPossible[current_ap][currentCloud])
                        {
                            currentMaxX = varValues.xval[current_ap][currentCloud][t];
                            currentCloudMaxX = currentCloud;
                            if(currentMaxX > 0.5)
                                break;
                        }
                    }
                }
                
                if(currentCloudMaxX == -1)
                {
                    cout << "impossible to reach integer solution" << endl;
                    cerr << "impossible to reach integer solution" << endl;
                    return integerSolution;
                }
                else
                {
                    // update residual capacity
                    cloudCap[currentCloudMaxX] -= nodeDemand[current_ap][t];
                    
                    // update corresponding xvalues for the AP
                    for(int k=0; k< cloudletCard;k++)
                    {
                        varValues.xval[current_ap][k][t] = 0.0;
                    }
                    varValues.xval[current_ap][currentCloudMaxX][t] = 1.0;
                }
            }
        }
    }
    
    // COMPUTE Y VARIABLES VALUES A POSTERIORI
    this->computeYvarValuesGivenXvarValues(varValues);
    
    // compute total cost a posteriori
    integerSolution.finalCost = computeCostGivenVariableValues(varValues, is_cycle);
    
    // print integer cost and solution
    cout << "rounded integer cost " << integerSolution.finalCost << endl;
    cerr << "rounded integer cost " << integerSolution.finalCost << endl;
    
    // InputDataApCloudletAssoc::printVarValues(varValues);
    
    integerSolution.variableValues = varValues;
    
    checkFeasibilitySolution(integerSolution.variableValues);
    
    return integerSolution;
}

double InputDataApCloudletAssoc::getBestIntegerCost(){
    return this->bestIntegerSolution.finalCost;
}

bool InputDataApCloudletAssoc::greedy_heuristic_betapositive(bool is_cycle){
    // vector< vector< vector< double> > > xval;
    // vector< vector< vector< vector<double> > > > yval;
    
    bool isGreedyOk = true;
    
    int nA = apCard;
    int nT = timeSlotCard;
    int nK = cloudletCard;
    int i,t;
    
    double maxCapacity = U*cloudletCap;
    double objValue = 0.0;
    
    VarValues greedyVarValues(nA,nK,nT);
    
    vector<int> previousAssociation(nA, -1);
    vector<int> first_time_association(nA, -1);
    
    // for every time-slot
    // for each AP in descending demand order
    // associate to nearest cloudlet without exceeding capacity
    // possibly without changing association with cloudlet of previous time-slot
    for(t = 0; t < nT; t++){
        // keep track of cloudlet usage in every time-slot given by fixed AP associations
        vector<double> current_usage(cloudletCard, 0.0);
        
        // sort AP demand in descending order
        vector<double> demand_in_time = timeDemand[t];
        //for(i = 0; i < nA; i++){
        //    demand_in_time.push_back(nodeDemand[i][t]);
        //}
        vector<size_t> sorted_ap_idx = sortedApIdxDescendingDemand[t]; // sort_indexes(demand_in_time);
        for(i = 0; i < nA; i++){
            // associate to nearest cloudlet without exceeding capacity
            // possibly without changing association with cloudlet of previous time-slot
            int current_ap = -1;
            if(isBetaZero)
                current_ap = (int) sorted_ap_number_cloudlet_reached[i];
            else
                current_ap = (int) sorted_ap_idx[i];
            
            double current_demand = demand_in_time[current_ap];
            
            int chosen_cloudlet = -1;
            
            if(previousAssociation[current_ap]!=-1){
                if(current_usage[previousAssociation[current_ap]] + current_demand <= maxCapacity ){
                    chosen_cloudlet = previousAssociation[current_ap];
                }
            }
            
            if(chosen_cloudlet==-1)
                chosen_cloudlet = nearestCloudletNotExceedingCapacity(current_usage,
                                                                      distanceApCloudlet[current_ap],
                                                                      current_demand,
                                                                      maxCapacity,
                                                                      isApCloudletPossible[current_ap]);
            
            if(chosen_cloudlet!=-1){
                greedyVarValues.xval[current_ap][chosen_cloudlet][t] = 1.0;
                current_usage[chosen_cloudlet] += current_demand;
                cout << "associate AP " << current_ap << " - " << chosen_cloudlet << " - " << t << endl;
                
                // beta*d[i,t]*m[i,k1]*x[i,k1,t]
                objValue += beta*distanceApCloudlet[current_ap][chosen_cloudlet]*current_demand;
                
                if(chosen_cloudlet!=previousAssociation[current_ap] && t > 0)
                {
                    greedyVarValues.yval[current_ap][previousAssociation[current_ap]][chosen_cloudlet][t] = 1.0;
                    // alpha*sum{k2 in K} (d[i,t]*l[k1,k2]*y[i,k1,k2,t])
                    objValue += alpha*current_demand*distanceCloudletCloudlet[previousAssociation[current_ap]][chosen_cloudlet];
                }
                
                previousAssociation[current_ap] = chosen_cloudlet;
                
                if(t == 0){
                    first_time_association[current_ap] = chosen_cloudlet;
                }
                else if(t == timeSlotCard - 1 && is_cycle){
                    objValue += alpha*nodeDemand[current_ap][0]*distanceCloudletCloudlet[chosen_cloudlet][first_time_association[current_ap]];
                }
            }
            else{
                cout << "ERROR! No cloudlet found for AP " << current_ap << " time " << t << endl;
                cout << "GREEDY HEURISTIC TERMINATED " << endl;
                // exit(-1);
                isGreedyOk = false;
                return isGreedyOk;
            }
        }
    }
    cout << "GREEDY SOLUTION " << endl;
    cout << "final cost = " << objValue << endl;
    
    InputDataApCloudletAssoc::printVarValues(greedyVarValues);
    
    if(bestIntegerSolution.finalCost == -1 || bestIntegerSolution.finalCost > objValue){
        // save initial integer solution as best
        this->bestIntegerSolution.variableValues = greedyVarValues;
        this->bestIntegerSolution.isBinary = true;
        this->bestIntegerSolution.finalCost = objValue;
        this->bestIntegerSolution.is_cycle = is_cycle;
    }
    
    return isGreedyOk;
}

bool InputDataApCloudletAssoc::_greedy_heuristic_betazero_single_assingment(int t,
                                                                    vector<int> & previousAssociation,
                                                                    vector<int> & first_time_association,
                                                                    double & objValue,
                                                                    VarValues & greedyVarValues,
                                                                    bool timeforward,
                                                                    bool is_cycle)
{
    bool isGreedyOk = true;
    
    int nA = apCard;
    int i;
    
    double maxCapacity = U*cloudletCap;
    
    // keep track of cloudlet usage in every time-slot given by fixed AP associations
    vector<double> current_usage(cloudletCard, 0.0);
    
    // sort AP demand in descending order
    vector<double> demand_in_time = timeDemand[t];
    for(i = 0; i < nA; i++)
    {
        // associate to nearest cloudlet without exceeding capacity
        // possibly without changing association with cloudlet of previous time-slot
        int current_ap = (int) sorted_ap_number_cloudlet_reached[i];
        
        double current_demand = demand_in_time[current_ap];
        
        int chosen_cloudlet = -1;
        
        if(previousAssociation[current_ap]!=-1){
            if(current_usage[previousAssociation[current_ap]] + current_demand <= maxCapacity ){
                chosen_cloudlet = previousAssociation[current_ap];
            }
        }
        
        if(chosen_cloudlet==-1)
            chosen_cloudlet = nearestCloudletNotExceedingCapacity(current_usage,
                                                                  distanceApCloudlet[current_ap],
                                                                  current_demand,
                                                                  maxCapacity,
                                                                  isApCloudletPossible[current_ap]);
        
        if(chosen_cloudlet!=-1){
            greedyVarValues.xval[current_ap][chosen_cloudlet][t] = 1.0;
            current_usage[chosen_cloudlet] += current_demand;
            cout << "associate AP " << current_ap << " - " << chosen_cloudlet << " - " << t << endl;
            
            // beta*d[i,t]*m[i,k1]*x[i,k1,t]
            objValue += beta*distanceApCloudlet[current_ap][chosen_cloudlet]*current_demand;
            
            if(chosen_cloudlet!=previousAssociation[current_ap]
               && t > 0
               && previousAssociation[current_ap] != -1)
            {
                if(timeforward){
                    greedyVarValues.yval[current_ap][previousAssociation[current_ap]][chosen_cloudlet][t] = 1.0;
                    // alpha*sum{k2 in K} (d[i,t]*l[k1,k2]*y[i,k1,k2,t])
                    objValue += alpha*current_demand*distanceCloudletCloudlet[previousAssociation[current_ap]][chosen_cloudlet];
                }
                else{
                    greedyVarValues.yval[current_ap][chosen_cloudlet][previousAssociation[current_ap]][t] = 1.0;
                    // alpha*sum{k2 in K} (d[i,t]*l[k1,k2]*y[i,k1,k2,t])
                    objValue += alpha*current_demand*distanceCloudletCloudlet[chosen_cloudlet][previousAssociation[current_ap]];
                }
            }
            
            previousAssociation[current_ap] = chosen_cloudlet;
            
            if(t == 0){
                first_time_association[current_ap] = chosen_cloudlet;
            }
            else if(t == timeSlotCard - 1 && is_cycle && first_time_association[current_ap] != -1)
            {
                objValue += alpha*nodeDemand[current_ap][0]*distanceCloudletCloudlet[chosen_cloudlet][first_time_association[current_ap]];
            }
        }
        else{
            cout << "ERROR! No cloudlet found for AP " << current_ap << " time " << t << endl;
            cout << "GREEDY HEURISTIC TERMINATED " << endl;
            isGreedyOk = false;
        }
    }
    return isGreedyOk;
}

bool InputDataApCloudletAssoc::greedy_heuristic_betazero(bool is_cycle)
{
    bool isGreedyOk = true;
    
    int nA = apCard;
    int nT = timeSlotCard;
    int nK = cloudletCard;
    int i,t, k1;
    
    // double maxCapacity = U*coudletCap;
    double objValue = 0.0;
    
    VarValues greedyVarValues(nA,nK,nT);
    
    vector<int> previousAssociation(nA, -1);
    vector<int> first_time_association(nA, -1);
    
    // find time with maximum sum of traffic
    int time_max = -1;
    double max_traffic = DBL_MIN;
    for(t = 0; t < nT; t++){
        double current_traffic_time = 0;
        for(i = 0; i < nA; i++){
            current_traffic_time += timeDemand[t][i];
        }
        if(current_traffic_time > max_traffic){
            max_traffic = current_traffic_time;
            time_max = t;
        }
    }
    
    cout << "time with maximum traffic " << time_max << " " << max_traffic << endl;
    
    // find first a solution for the time-slot with maximum total traffic
    // then find solution for all preceding and succeding time-slots, in order
    // possibly without changing association with cloudlet of previous (or next) time-slot
    // for each AP in descending demand order
    // associate to nearest cloudlet without exceeding capacity
    
    isGreedyOk = _greedy_heuristic_betazero_single_assingment(time_max,
                                                              previousAssociation,
                                                              first_time_association,
                                                              objValue,
                                                              greedyVarValues,
                                                              true,
                                                              is_cycle);
    // find solution for all time-slot preceding time_max
    t = time_max-1;
    while(t >= 0 && isGreedyOk){
        isGreedyOk = _greedy_heuristic_betazero_single_assingment(t,
                                                          previousAssociation,
                                                          first_time_association,
                                                          objValue,
                                                          greedyVarValues,
                                                          false,
                                                          is_cycle);
        t--;
    }
    // find solution for all time-slots succeding time-max
    // restore previousAssociation vector to the assignments of time-max time-slot
    t = time_max + 1;
    for(i = 0; i < nA; i++){
        for(k1 = 0; k1 < nK; k1++){
            if(greedyVarValues.xval[i][k1][time_max] > 1-1e-3){
                previousAssociation[i] = k1;
            }
        }
    }
    while(t < nT && isGreedyOk){
        isGreedyOk = _greedy_heuristic_betazero_single_assingment(t,
                                                          previousAssociation,
                                                          first_time_association,
                                                          objValue,
                                                          greedyVarValues,
                                                          true,
                                                          is_cycle);
        t++;
    }
    
    if(!isGreedyOk){
        return false;
    }
    
    cout << "GREEDY SOLUTION " << endl;
    cout << "final cost = " << objValue << endl;
    
    InputDataApCloudletAssoc::printVarValues(greedyVarValues);
    
    if(bestIntegerSolution.finalCost == -1 || bestIntegerSolution.finalCost > objValue){
        // save initial integer solution as best
        this->bestIntegerSolution.variableValues = greedyVarValues;
        this->bestIntegerSolution.isBinary = true;
        this->bestIntegerSolution.finalCost = objValue;
        this->bestIntegerSolution.is_cycle = is_cycle;
    }
    
    return isGreedyOk;
}

// simplest greedy heuristic. no guarantee on optimality. guarantee on feasibility.
bool InputDataApCloudletAssoc::greedy_heuristic(bool is_cycle){
    if(!isBetaZero){
        return greedy_heuristic_betapositive(is_cycle);
    }
    else{
        return greedy_heuristic_betazero(is_cycle);
        // return greedy_heuristic_betapositive(is_cycle);
    }
}

Solution InputDataApCloudletAssoc::getBestIntegerSolution(bool is_cycle ){
    if(bestIntegerSolution.finalCost == -1){
        // try to get an integer solution via greedy heuristic
        if(!this->greedy_heuristic(is_cycle)){
            // greedy heuristic not successfull
            cout << "NO CURRENT BEST INTEGER SOLUTION " << endl;
        }
    }
    return this->bestIntegerSolution;

}

double InputDataApCloudletAssoc::getIntegralityGAP()
{
    return integralityGAP;
}

void InputDataApCloudletAssoc::setCGRootLB(double newRootLB)
{
    this->cg_root_LB = newRootLB;
    this->min_active_branch_node_LB = newRootLB;
    if(this->bestIntegerSolution.finalCost > 0 && newRootLB > 0)
    {
        this->integralityGAP = (bestIntegerSolution.finalCost - cg_root_LB)/cg_root_LB;
    }
}

void InputDataApCloudletAssoc::setMinActiveBranchNodeLB(double newLB)
{
    this->min_active_branch_node_LB = newLB;
    if(this->bestIntegerSolution.finalCost > 0 && min_active_branch_node_LB > 0)
    {
        this->integralityGAP = (bestIntegerSolution.finalCost - min_active_branch_node_LB)/min_active_branch_node_LB;
    }
}

double InputDataApCloudletAssoc::getMinActiveBranchNodeLB()
{
    return min_active_branch_node_LB;
}

void InputDataApCloudletAssoc::setBestIntegerSolution(Solution newBestIntegerSolution, bool is_cycle)
{
    if(bestIntegerSolution.finalCost == -1 ||
       (bestIntegerSolution.finalCost > -1 && newBestIntegerSolution.finalCost < bestIntegerSolution.finalCost))
    {
        std::stringstream to_print;
        to_print.str("");
        to_print << "FOUND NEW BEST INTEGER SOLUTION WITH COST = " << newBestIntegerSolution.finalCost;
        
        if(cg_root_LB > 0)
        {
            double gap_to_use = cg_root_LB;
            if(min_active_branch_node_LB > cg_root_LB)
            {
                gap_to_use = min_active_branch_node_LB;
            }
            integralityGAP = (newBestIntegerSolution.finalCost - gap_to_use)/gap_to_use;
            to_print << " LB = " << gap_to_use << " GAP =  " << integralityGAP;
        }
        
        to_print << " after exc. time = " << computeTime(this->initExecutionTime);
        
        cout << to_print.str() << endl;
        cerr << to_print.str() << endl;
        this->bestIntegerSolution = newBestIntegerSolution;
        bestIntegerSolution.is_cycle = is_cycle;
    }
}

void InputDataApCloudletAssoc::printBestIntegerSolution(){
    // check if a current best integer solution exists
    if(bestIntegerSolution.finalCost == -1){
        // try to get an integer solution via greedy heuristic
        if(!this->greedy_heuristic()){
            // greedy heuristic not successfull
            cout << "NO CURRENT BEST INTEGER SOLUTION " << endl;
        }
    }
    
    // check if a current best integer solution exists
    if(bestIntegerSolution.finalCost != -1){
        // print best integer solution cost and variables
        cout << " BEST INTEGER SOLUTION " << endl;
        cout << " BEST INTEGER COST  " << this->bestIntegerSolution.finalCost << endl;
        double assignCost = 0;
        double switchCost = 0;
        this->getAssignmentAndSwitchingCostSolution(assignCost, switchCost, bestIntegerSolution);
        cout << " BEST INTEGER ASSIGN COST " << assignCost << " SWITCHING COST " << switchCost << endl;
        InputDataApCloudletAssoc::printVarValues(this->bestIntegerSolution.variableValues);
    }
    else{
        cout << "NO CURRENT BEST INTEGER SOLUTION " << endl;
    }
    
    // printAssignmentMatrixBestIntegerSolution();
}

void InputDataApCloudletAssoc::printVarValues(VarValues varValues){
    int i,k1,k2,t;
    cout << "x val positive values" << endl;
    // print variables value
    for(i = 0; i < varValues.xval.size(); i++){
        for(k1 = 0; k1< varValues.xval[i].size();k1++){
            for(t = 0; t< varValues.xval[i][k1].size(); t++){
                double xval = varValues.xval[i][k1][t];
                if(xval > 1e-6){
                    cout << i << " " << k1 << " " << t << " = " << xval << endl;
                }
            }
        }
    }
    
    cout << "y val positive values" << endl;
    // print variables value
    for(i = 0; i < varValues.yval.size(); i++){
        for(k1 = 0; k1< varValues.yval[i].size();k1++){
            for(k2 = 0; k2< varValues.yval[i][k1].size();k2++) if(k1!=k2) {
                for(t = 0; t< varValues.yval[i][k1][k2].size(); t++){
                    double yval = varValues.yval[i][k1][k2][t];
                    if(yval > 1e-6){
                        cout << i << " " << k1 << " " << k2 << " " << t << " = " << yval << endl;
                    }
                }
            }
        }
    }
    cout << "end of y val" << endl;
}

double InputDataApCloudletAssoc::getDummyVarOFcost()
{
    if(dummyCost < 0)
    {
        dummyCost = 0;
        // cost of the use of dummy variable is the sum of all possible costs
        for(int i = 0; i < this->apCard;i++){
            for(int t = 0; t < this->timeSlotCard;t++){
                for(int k=0; k < this->cloudletCard;k++){
                    dummyCost += this->nodeDemand[i][t]*this->distanceApCloudlet[i][k];
                    if(k < this->cloudletCard-2){
                        for(int j = k+1; j < this->cloudletCard;j++) if(j!=k) {
                            dummyCost += this->nodeDemand[i][t]*this->distanceCloudletCloudlet[k][j];
                            dummyCost += this->nodeDemand[i][t]*this->distanceCloudletCloudlet[j][k];
                        }
                    }
                }
            }
        }
    }
    return dummyCost;
}

void InputDataApCloudletAssoc::setInitialExecutionTime(time_t initTime)
{
    this->initExecutionTime = initTime;
}

void InputDataApCloudletAssoc::printAssignmentMatrixBestIntegerSolution()
{
    if(bestIntegerSolution.finalCost > -1)
    {
        cout << "BEST INTEGER SOLUTION ASSIGNMENT MATRIX" << endl;
        this->printAssignmentMatrix(computeAssignmentMatrix(bestIntegerSolution.variableValues));
        cout << "END ASSIGNMENT MATRIX" << endl;
    }
}

vector<vector<int>> InputDataApCloudletAssoc::computeAssignmentMatrix(VarValues varValues){
    int apCard = (int) varValues.xval.size();
    int timeCard = (int) varValues.xval[0][0].size();
    int cloudCard = (int) varValues.xval[0].size();
    
    vector<vector<int>> assign_matrix(apCard, vector<int>(timeCard, -1));
    
    for(int i = 0; i < apCard; i++){
        for(int t = 0; t < timeCard; t++){
            for(int k = 0; k < cloudCard; k++){
                if(IS_GE_ONE(varValues.xval[i][k][t])){
                    assign_matrix[i][t] = k;
                    break;
                }
            }
        }
    }
    
    return assign_matrix;
}
void InputDataApCloudletAssoc::printAssignmentMatrix(vector<vector<int>> assign_matrix)
{
    for(int i = 0; i < assign_matrix.size(); i++){
        for(int t = 0; t < assign_matrix[i].size(); t++){
            cout << assign_matrix[i][t] + 1;
            if(t != assign_matrix[i].size() - 1)
            {
                cout << ",";
            }
        }
        cout << endl;
    }
}

void InputDataApCloudletAssoc::getAssignmentAndSwitchingCostVarValues(
                                            double & assignCost, double & switchingCost,
                                            VarValues variableValues, bool is_cycle)
{
    assignCost = 0;
    switchingCost = 0;
    
    int cloud_first_time = -1;
    int prev_cloud = -1;
    for(int i=0; i< apCard; i++)
    {
        for(int t = 0; t < timeSlotCard; t++)
        {
            for(int k = 0; k < cloudletCard; k++)
            {
                if(variableValues.xval[i][k][t] > 1e-6)
                {
                    assignCost += beta*nodeDemand[i][t]*distanceApCloudlet[i][k]*variableValues.xval[i][k][t];
                    if(t == 0)
                    {
                        cloud_first_time = k;
                    }
                    else if(prev_cloud != k && prev_cloud != -1)
                    {
                        switchingCost += alpha*nodeDemand[i][t]*distanceCloudletCloudlet[prev_cloud][k];
                    }
                    if (t == timeSlotCard - 1 && is_cycle && cloud_first_time!=k)
                    {
                        switchingCost += alpha*nodeDemand[i][0]*distanceCloudletCloudlet[k][cloud_first_time];
                    }
                    prev_cloud = k;
                }
            }
        }
    }
}

void InputDataApCloudletAssoc::getAssignmentAndSwitchingCostSolution(
                                            double & assignCost,
                                           double & switchingCost,
                                           Solution sol)
{
    getAssignmentAndSwitchingCostVarValues(assignCost, switchingCost, sol.variableValues, sol.is_cycle);
}

void InputDataApCloudletAssoc::getAssignmentAndSwitchingCostBestIntegerSolution(double & assignCost
                                                                                , double & switchingCost)
{
    // check if a current best integer solution exists
    if(bestIntegerSolution.finalCost == -1)
    {
        cout << "NO CURRENT BEST INTEGER SOLUTION " << endl;
    }
    else{
        if(bestIntegerSolution.assignmentCost != -1 &&
           bestIntegerSolution.switchingCost != -1)
        {
            assignCost = bestIntegerSolution.assignmentCost;
            switchingCost = bestIntegerSolution.switchingCost;
        }
        else{
            getAssignmentAndSwitchingCostSolution(assignCost, switchingCost, bestIntegerSolution);
        }
    }
}
