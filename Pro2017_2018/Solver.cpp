//
//  Solver.cpp
//  Pro2017_2018
//
//  Created by 任智桂 on 2019/3/12.
//  Copyright © 2019 triple S. All rights reserved.
//

#include "Solver.h"
double Newton(BSCall *FCptr, double target,double x0,
              double epsilon, int MaxIter)
{
    double xn, xnplus1, funxn, derivxn, testaccuracy;
    xn=x0;
    testaccuracy=epsilon+10.0;
    int i=0;
    while ((i < MaxIter)&&(testaccuracy > epsilon))
    {
        funxn=FCptr->Pricing(xn);
        derivxn=FCptr->Vega(xn);
        xnplus1=xn-(funxn-target)/derivxn;
        testaccuracy = fabs(xnplus1-xn);
        xn=xnplus1;
        i++;
    }
    return xn;
}
double Bisection(BSCall *FCptr,double target,double a,
                 double b,double epsilon,int MaxIter)
{
    double x = 0.5*(a+b);
    double y = FCptr->Pricing(x);
    int i=0;
    while ((i < MaxIter)&&(fabs(y-target) > epsilon))
    {
        if (y<target)
            a = x;
        if (y>target)
            b = x;
        x = 0.5*(a+b);
        y = FCptr->Pricing(x);
    }
    return x;
}
