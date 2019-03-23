//
//  Minimize.cpp
//  Pro2017_2018
//
//  Created by 任智桂 on 2019/3/14.
//  Copyright © 2019 triple S. All rights reserved.
//

#include "TheFunction.h"
#include <cstdlib>
#include <gsl/gsl_multimin.h>
#include <cmath>
double chi_z(double z, double rho)
{
    double temp = log((sqrt(1 - 2 * rho * z + z * z) + z - rho) / (1 - rho));
    return temp;
}


double z_num(double alpha, double beta, double nu,
             double K, double f)
{
    double temp = nu / alpha * (pow((f * K), (0.5 - 0.5 * beta))) * log(f / K);
    return temp;
}
double trans_alpha(double x)
{
    return exp(x);
}
double trans_beta(double x)
{
    return (1 / (1 + exp(x)));
}
double trans_rho (double x)
{
    return (2 / (1 + exp(x)) - 1);
}
double trans_nu(double x)
{
    return fabs(x);
}
double sigma_B(double alpha, double beta,
               double rho, double nu, double K,
               double f, double t_ex)
{
    double z = z_num(alpha, beta, nu, K, f);
    double x_z = chi_z(z, rho);
    
    double numerator = 1 + ((1 - beta) * (1- beta) / 24 * alpha * alpha / pow((f*K), 1 - beta)
                            + 0.25 * rho * beta * nu * alpha  / pow((f*K), (0.5 - 0.5 * beta))
                            + ((2 - 3 * rho * rho) / 24 * nu * nu)) * t_ex;
    
    double denominator = x_z * pow((f*K),(0.5 - 0.5 * beta)) * (1 + pow((1 - beta), 2) / 24 * pow((log(f / K)), 2) + pow((1 - beta), 4) / 1920 * pow((log(f / K)), 4));
    
    double result =  alpha * z * numerator / denominator;
    return(result);
}
double my_f(const gsl_vector *v, void *params)
{
    double *sigma = new double[7];
    double *p = (double *)params;
    double strike[7] = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07};
    //double strike[7] = {0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01};
    double sum = 0.0;
    
    double alpha, beta, rho, nu;
    alpha = gsl_vector_get(v, 0);
    beta = gsl_vector_get(v, 1);
    rho = gsl_vector_get(v, 2);
    nu = gsl_vector_get(v, 3);
    
    for (int i = 0; i < 7; i++)
    {
        sigma[i] = sigma_B(trans_alpha(alpha), trans_beta(beta), trans_rho(rho), trans_nu(nu), strike[i]);
        sum += ((sigma[i] - p[i]) * (sigma[i] - p[i]))/(p[i]*p[i]);
    }
    
    delete[] sigma;
    return sum;
}
double my_f_2(const gsl_vector *v, void *params)
{
    double *sigma = new double[7];
    double *p = (double *)params;
    double strike[7] = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07};
    //double strike[7] = {0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01};
    double sum = 0.0;
    double alpha,rho, nu;
    alpha = gsl_vector_get(v, 0);
    rho = gsl_vector_get(v, 1);
    nu = gsl_vector_get(v, 2);
    double beta = 1.0;
    for (int i = 0; i < 7; i++)
    {
        sigma[i] = sigma_B(trans_alpha(alpha),beta,trans_rho(rho),trans_nu(nu), strike[i]);
        sum += ((sigma[i] - p[i]) * (sigma[i] - p[i]))/(p[i]*p[i]);
    }
    delete[] sigma;
    return sum ;
}
double my_f_3(const gsl_vector *v, void *params)
{
    double *sigma = new double[7];
    double *p = (double *)params;
    double strike[7] = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07};
    //double strike[7] = {0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01};
    double sum = 0.0;
    double alpha,rho, nu;
    alpha = gsl_vector_get(v, 0);
    rho = gsl_vector_get(v, 1);
    nu = gsl_vector_get(v, 2);
    double beta = 0.5;
    for (int i = 0; i < 7; i++)
    {
        sigma[i] = sigma_B(trans_alpha(alpha),beta,trans_rho(rho),trans_nu(nu), strike[i]);
        sum += ((sigma[i] - p[i]) * (sigma[i] - p[i]))/(p[i]*p[i]);
    }
    delete[] sigma;
    return sum ;
}
