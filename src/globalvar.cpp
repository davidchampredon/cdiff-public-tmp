//
//  globalvar.cpp
//  abm-cdiff
//
//  Created by David CHAMPREDON on 2017-06-19.
//  Copyright © 2017 David CHAMPREDON. All rights reserved.
//

#include "globalvar.hpp"


// ==== Random seeds ====

uint RANDOM_SEED	= 12345;

std::mt19937_64	RANDOM_GENERATOR(RANDOM_SEED);

uint ISO_WARD = 0;

uint HUGENUMBER = 999999;
uint NA_UINT	= 999999;
float NA_FLOAT	= 999999;

uint MAX_N_RELAPSES	= 5;
