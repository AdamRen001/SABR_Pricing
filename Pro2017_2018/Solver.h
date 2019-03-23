//
//  Bisection.h
//  Pro2017_2018
//
//  Created by 任智桂 on 2019/3/12.
//  Copyright © 2019 triple S. All rights reserved.
//

#ifndef Solver_h
#define Solver_h
#include <cmath>
#include "BSCallClass.h"
double Newton(BSCall *FCptr, double target,double x0,
              double epsilon, int MaxIter);
double Bisection(BSCall *FCptr,
                 double target,
                 double a,
                 double b,
                 double epsilon,
                 int MaxIter);
#endif /* Solver_h */
