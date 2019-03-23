//
//  MonteCarlo.h
//  Pro2017_2018
//
//  Created by 任智桂 on 2019/3/19.
//  Copyright © 2019 triple S. All rights reserved.
//

#ifndef MonteCarlo_h
#define MonteCarlo_h
double FirstEulerSample(double f0, double t_ex,double alpha0,
                        double beta,double rho,double nu);
double MonteCarlo(double &sd1,int NumberOfPaths,double alpha0,
                  double beta0,double rho0,double nu0,double Strike,
                  double f0,double t_ex,double D);
double Antethetic(double &sd2,int NumberOfPaths,double alpha0,
                  double beta0,double rho0,double nu0,double Strike,
                  double f0,double t_ex,double D);
double ControlVariates(double &sd3,int NumberOfPaths,double alpha0,
                       double beta0,double rho0,double nu0,double Strike,
                       double f0,double t_ex,double D);

#endif /* MonteCarlo_h */
