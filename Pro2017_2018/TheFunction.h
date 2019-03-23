//
//  Minimize.h
//  Pro2017_2018
//
//  Created by 任智桂 on 2019/3/13.
//  Copyright © 2019 triple S. All rights reserved.
//

#ifndef TheFunction_h
#define TheFunction_h
#include <gsl/gsl_multimin.h>
double chi_z(double z, double rho);
double z_num(double alpha, double beta, double nu, double K, double f = 0.025);
double trans_alpha(double x);
double trans_beta(double x);
double trans_rho (double x);
double trans_nu(double x);
double sigma_B(double alpha, double beta,
               double rho, double nu, double K, double f = 0.025, double t_ex = 1.0);
double my_f(const gsl_vector *v, void *params);
double my_f_2(const gsl_vector *v, void *params);
double my_f_3(const gsl_vector *v, void *params);
#endif /* TheFunction_h */
