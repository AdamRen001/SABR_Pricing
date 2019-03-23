//
//  BSFormula.cpp
//  Pro2017_2018
//
//  Created by 任智桂 on 2019/3/12.
//  Copyright © 2019 triple S. All rights reserved.
//

#include "BSFormula.h"
#include "Normals.h"
#include <cmath>
using namespace std;
double BlackScholesCall(double f,
                        double K,
                        double sigma,
                        double D,
                        double t_ex)
{
    double d1 = (log(f/K)+ 0.5*sigma*sigma*t_ex)/(sigma*sqrt(t_ex));
    double d2 = d1 - sigma * sqrt(t_ex);
    double CallPrice = D * f * CumulativeNormal(d1) - K*CumulativeNormal(d2);
    return CallPrice;
}
double BlackScholesVega(double f,
                        double K,
                        double sigma,
                        double D,
                        double t_ex)
{
    double d1 = (log(f/K)+ 0.5*sigma*sigma*t_ex)/(sigma*sqrt(t_ex));
    double d2 = d1 - sigma * sqrt(t_ex);
    double phi_d1 = NormalDensity(d1);
    double phi_d2 = NormalDensity(d2);
    double vega = D*(f*phi_d1*((0.5*sigma*sigma*t_ex-log(f/K))/(sigma*sigma*sqrt(t_ex)))+K*phi_d2*((0.5*sigma*sigma*t_ex+log(f/K))/(sigma*sigma*sqrt(t_ex))));
    return vega;
}
