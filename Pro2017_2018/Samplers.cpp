//
//  Samplers.cpp
//  Pro2017_2018
//
//  Created by 任智桂 on 2019/3/14.
//  Copyright © 2019 triple S. All rights reserved.
//
#include<cstdlib>
#include<cmath>
#include "Samplers.h"
double Uniform::getnumber()
{
    double x;
    x = a + rand()*(b-a)/static_cast<double>(RAND_MAX);
    return x;
}
double InverseFuncSample::getnumber()
{
    double u,result;
    u = rand()/static_cast<double>(RAND_MAX);
    result = -log(1-u)/lambda;
    return result;
}
double BoxMuller::getnumber()
{
    double u1,u2,result,R,theta;
    u1 = rand()/static_cast<double>(RAND_MAX);
    u2 = rand()/static_cast<double>(RAND_MAX);
    R = -2*log(u1);
    theta = 2*3.141592653589793238462643*u2;
    result = sqrt(R)*cos(theta);
    return result;
}
double DoubleExp::getnumber()
{
    double u1=0.0;
    double result;
    while (u1==0.0||u1==1.0)
    {
        u1 = rand()/static_cast<double>(RAND_MAX);
    }
    if ((u1>0)&&(u1<=0.5))
    {
        result=log(2*u1);
    }
    else
    {
        result=-log(2*(1-u1));
    }
    return result;
}
double Rej_Samp::getnumber()
{
    double x=1.0;
    double u1=1000.0;
    double fx=1.0;
    double gx=1.0;
    while (u1>(fx/(gx*C)))
    {
        x = pt->getnumber();
        u1 = rand()/static_cast<double>(RAND_MAX);
        fx = sqrt(2*3.1415)*exp(-0.5*x*x);
        gx = 0.5*exp(-abs(x));
    }
    return x;
}

