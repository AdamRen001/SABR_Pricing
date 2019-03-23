//
//  SABRFormula.h
//  Pro2017_2018
//
//  Created by 任智桂 on 2019/3/12.
//  Copyright © 2019 triple S. All rights reserved.
//

#ifndef SABRFormula_h
#define SABRFormula_h
double SABR_IV(double f,double K,
               double t_ex,double alpha,
               double beta,double rho,
               double nu);
double SABRCall(double f,double K,
                double D,double t_ex,
                double alpha,double beta,
                double rho,double nu,double sigmaB);
double PayOff(double Spot,double Strike);
#endif /* SABRFormula_h */
