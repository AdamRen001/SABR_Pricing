//
//  SABRFormula.cpp
//  Pro2017_2018
//
//  Created by 任智桂 on 2019/3/12.
//  Copyright © 2019 triple S. All rights reserved.
//

#include "SABRFormula.h"
#include <cmath>
#include "Normals.h"
#include "TheFunction.h"
#include <iostream>
using namespace std;
double SABR_IV(double f,double K,
               double t_ex,double alpha,
               double beta,double rho,
               double nu)
{
    double SigmaB;
    if (fabs(K-f)<1e-6)
    {
        cout<<"ATM CASE!!!!!!!"<<endl;
        double first = (1-beta)*(1-beta)/24*alpha*alpha/pow(f,2-2*beta);
        double second = 0.25*rho*beta*alpha*nu/pow(f,1-beta);
        double third = (2 - 3*rho*rho)*nu*nu/24;
        double multiplier = 1 + (first+second+third)*t_ex;
        SigmaB = (alpha/pow(f,1-beta))*multiplier;
    }
    else
    {
        SigmaB = sigma_B(alpha,beta,rho,nu,K,f,t_ex);
    }
    return SigmaB;
}
double SABRCall(double f,double K,
                double D,double t_ex,
                double alpha,double beta,
                double rho,double nu,double SigmaB)
{
    double d1 = (log(f/K)+0.5*SigmaB*SigmaB*t_ex)/(SigmaB*sqrt(t_ex));
    double d2 = d1 - SigmaB*sqrt(t_ex);
    double CallPrice=D*(f*CumulativeNormal(d1)-K*CumulativeNormal(d2));
    return CallPrice;
}
double PayOff(double Spot,double Strike)
{
    return (Spot-Strike)>0.0?(Spot-Strike):0.0;
    //return std::max(Spot-Strike,0.0);
}
