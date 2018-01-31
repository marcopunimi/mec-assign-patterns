//
//  main.hpp
//  ap_cloudlet_association
//
//  Created by Marco Premoli on 31/01/2017.
//  Copyright © 2017 Marco Premoli. All rights reserved.
//

#pragma once

#ifndef main_h
#define main_h

#include "ap_cloudlet_assoc_model.hpp"
#include <iterator>
#include "ShortPathHeurBP.hpp"
//#include "GenerAssignHeur.hpp"
//#include "ShortPathHeurCycle.hpp"

void solveCompleteModel(IloEnv, InputDataApCloudletAssoc*,ModelVariant, double epgap = 0.001, int timelimit = -1);

#endif /* main_h */
