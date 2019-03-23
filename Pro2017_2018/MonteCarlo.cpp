//
//  MonteCarlo.cpp
//  Pro2017_2018
//
//  Created by 任智桂 on 2019/3/19.
//  Copyright © 2019 triple S. All rights reserved.
//

#include "MonteCarlo.h"
#include "Samplers.h"
#include "SABRFormula.h"
#include <cmath>
#include <iostream>
using namespace std;
double FirstEulerSample(double f0,double t_ex,
                        double alpha0,double beta,
                        double rho,double nu)
{
    const int N = 10;
    double h=t_ex/N;
    for (int i=0;i<N;i++)
    {
        BoxMuller bx;
        double z1 = bx.getnumber();
        double z2 = bx.getnumber();
        z1 = sqrt(1-rho*rho)*z1+rho*z2;
        f0 = f0*exp(alpha0*pow(f0,beta-1)*sqrt(h)*
                    z1-0.5*alpha0*alpha0*pow(f0,2*beta-2)*h);
        alpha0 = alpha0*exp(-0.5*nu*nu*h+nu*sqrt(h)*z2);
    }
    return f0;
}

double MonteCarlo(double &sd1,int NumberOfPaths,double alpha0,
                  double beta0,double rho0,double nu0,double Strike,
                  double f0,double t_ex,double D)
{
    double FT;
    double result=0.0;
    double *myvector;
    myvector = new double[NumberOfPaths];
    for (int i=0;i<NumberOfPaths;i++)
    {
        FT = FirstEulerSample(f0,t_ex,alpha0,beta0,rho0,nu0);
        myvector[i] = D*PayOff(FT, Strike);
        result += myvector[i];
    }
    result = result/NumberOfPaths;
    double temp=0.0;
    for (int i=0;i<NumberOfPaths;i++)
    {
        temp += (myvector[i]-result)*(myvector[i]-result);
    }
    sd1 = sqrt(temp/(NumberOfPaths-1.0));
    delete [] myvector;
    return result;
}
double Antethetic(double &sd2,int NumberOfPaths,double alpha,
                  double beta,double rho,double nu,double Strike,
                  double f0,double t_ex,double D)
{
    double FT,FT_2;
    double result=0.0;
    double initial_alpha = alpha;
    double initial_f0 = f0;
    double *myvector;
    double *myvector2;
    myvector = new double[NumberOfPaths];
    myvector2 = new double[NumberOfPaths];
    for (int i=0;i<NumberOfPaths;i++)
    {
        BoxMuller bx;
        double f0_2 = initial_f0;
        f0 = initial_f0;
        double alpha_2 = initial_alpha;
        alpha = initial_alpha;
        const int N = 100;
        double h=t_ex/N;
        for (int j=0;j<N;j++)
        {
            double z1 = bx.getnumber();
            double z2 = bx.getnumber();
            z1 = sqrt(1-rho*rho)*z1+rho*z2;
            double z3 = -z1;
            double z4 = -z2;
            f0 = f0*exp(alpha*pow(f0,beta-1)*
                        sqrt(h)*z1-0.5*alpha*alpha*pow(f0,2*beta-2)*h);
            f0_2 = f0_2*exp(alpha_2*pow(f0_2,beta-1)
                          *sqrt(h)*z3-0.5*alpha_2*alpha_2*pow(f0_2,2*beta-2)*h);
            alpha = alpha*exp(-0.5*nu*nu*h+nu*sqrt(h)*z2);
            alpha_2 = alpha_2*exp(-0.5*nu*nu*h+nu*sqrt(h)*z4);
        }
        FT = f0;
        FT_2 = f0_2;
        myvector[i] = D*PayOff(FT, Strike);//1.1
        myvector2[i] = D*PayOff(FT_2,Strike);
        result += (myvector[i]+myvector2[i])*0.5;
    }
    result = result/NumberOfPaths;
    double temp=0.0;
    for (int i=0;i<NumberOfPaths;i++)
    {
        temp += (myvector[i]-result)*(myvector[i]-result);
    }
    sd2 = sqrt(temp/(2*NumberOfPaths-1));
    delete [] myvector;
    delete [] myvector2;
    return result;
}
double ControlVariates(double &sd3,int NumberOfPaths,double alpha0,
                       double beta0,double rho0,double nu0,double Strike,
                       double f0,double t_ex,double D)
{
    double result=0.0;
    double FT;
    double *myvector,*X,*Y;
    myvector = new double[NumberOfPaths];
    X = new double[NumberOfPaths];
    Y = new double[NumberOfPaths];
    double Yi;
    double avg_X = 0.0;
    double avg_Y = 0.0;
    for (int i=0;i<NumberOfPaths;i++)
    {
        FT = FirstEulerSample(f0,t_ex,alpha0,beta0,rho0,nu0);
        X[i] = FT;
        avg_X += X[i];
        Y[i] = PayOff(FT, Strike);
        avg_Y += Y[i];
    }
    avg_X /= NumberOfPaths;
    avg_Y /= NumberOfPaths;
    double numerator=0.0;
    double denominator = 0.0;
    double b;
    for (int i=0;i<NumberOfPaths;i++)
    {
        numerator += (X[i]-avg_X)*(Y[i]-avg_Y);
        denominator += (X[i]-avg_X)*(X[i]-avg_X);
    }
    b = numerator/denominator;
    for (int i=0;i<NumberOfPaths;i++)
    {
        FT = X[i];
        Yi = PayOff(FT, Strike);
        myvector[i] = D*(Yi - b*(FT-f0));//1.1
        result += myvector[i];
    }
    result = result/NumberOfPaths;
    double temp=0.0;
    for (int i=0;i<NumberOfPaths;i++)
    {
        temp += (myvector[i]-result)*(myvector[i]-result);
    }
    sd3 = sqrt(temp/(NumberOfPaths-1));
    delete [] myvector;
    return result;
}
