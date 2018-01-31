//
//  supportFunctions.hpp
//  ap_cloudlet_association
//
//  Created by Marco Premoli on 31/01/2017.
//  Copyright © 2017 Marco Premoli. All rights reserved.
//
#pragma once
#ifndef supportFunctions_hpp
#define supportFunctions_hpp

#include <ctime>
// #include <vector>
// #include <iostream>

using namespace std;

#define IS_FRACTIONAL(x) (x > 1e-6 && x < 1-1e-6)
#define IS_GE_ONE(x) (x > 1-1e-6)
#define IS_LE_ZERO(x) (x < 1e-6)

double computeTime(time_t initialTime);

#endif /* supportFunctions_hpp */
