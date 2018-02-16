//
//  main.hpp
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
