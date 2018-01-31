//
//  supportFunctions.cpp
//  ap_cloudlet_association
//
//  Created by Marco Premoli on 01/02/2017.
//  Copyright © 2017 Marco Premoli. All rights reserved.
//

#include "supportFunctions.hpp"

double computeTime(time_t initialTime){
    time_t endTime;
    std::time(&endTime);
    return difftime(endTime, initialTime);
}
