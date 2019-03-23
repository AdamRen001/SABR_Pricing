//
//  SABRClass.cpp
//  Pro2017_2018
//
//  Created by 任智桂 on 2019/3/12.
//  Copyright © 2019 triple S. All rights reserved.
//

#include "SABRClass.h"
#include "SABRFormula.h"
#include "Samplers.h"
#include "TheFunction.h"
#include <iostream>
using namespace std;

SABR::SABR(double f_,
           double D_,double t_ex_):f(f_),t_ex(t_ex_),D(D_)
{}
double SABR::Implied_Vol(double K_) const
{
    return SABR_IV(f, K_, t_ex, alpha, beta, rho, nu);
}

double SABR::get_alpha(){return alpha;}
double SABR::get_beta(){return beta;}
double SABR::get_rho(){return rho;}
double SABR::get_nu(){return nu;}
//
void SABR::Fit_with_all()
{
    double par[7] = {0.2295, 0.2139, 0.2054, 0.1996, 0.1954, 0.1920, 0.1892};
    gsl_vector *x,*ss;
    x = gsl_vector_alloc(4);
    gsl_vector_set(x, 0, 1); //alpha
    gsl_vector_set(x, 1, 0.0);//beta
    gsl_vector_set(x, 2, 0.0);//rho
    gsl_vector_set(x, 3, 1);//nu
    
    // initial point setting
    ss = gsl_vector_alloc(4);
    gsl_vector_set_all(ss, 1.0);
    // initial point setting
    //ss = gsl_vector_alloc(4);
    //gsl_vector_set(ss,0,0.01); //alpha
    //gsl_vector_set(ss,1,0.0099); //beta
    //gsl_vector_set(ss,2,0.019); //rho
    //gsl_vector_set(ss,3,0.05); //nu
    
    // initialize method and iterate
    gsl_multimin_function minex_func;
    minex_func.n = 4;
    minex_func.f = my_f;
    minex_func.params = par;
    
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
    gsl_multimin_fminimizer *s = nullptr;
    s = gsl_multimin_fminimizer_alloc(T, 4);
    gsl_multimin_fminimizer_set(s, &minex_func,x,ss);
    size_t iter = 0;
    int status;
    double size;
    
    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        
        if(status)
            break;
        
        size = gsl_multimin_fminimizer_size(s);
        status = gsl_multimin_test_size(size, 1e-4);
        if (status == GSL_SUCCESS)
        {
            std::cout << "converged to minimum at \n";
        }
        
        cout << "iter = " << iter << ", alpha = " << trans_alpha(gsl_vector_get(s->x, 0))
             << ", beta = " << trans_beta(gsl_vector_get(s->x, 1))
             << ",rho = " << trans_rho(gsl_vector_get(s->x, 2))<< ", nu = "
             << trans_nu(gsl_vector_get(s->x, 3)) << ", f(para) = " << s->fval
             << ", size = " << size << std::endl;
    }
    while (status == GSL_CONTINUE && iter < 10000);
    alpha = trans_alpha(gsl_vector_get(s->x, 0));
    beta = trans_beta(gsl_vector_get(s->x,1));
    rho = trans_rho(gsl_vector_get(s->x,2));
    nu = trans_nu(gsl_vector_get(s->x,3));
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);
}
void SABR::Fit_with_fixed_beta(double b)
{
    double par[7] = {0.2295, 0.2139, 0.2054, 0.1996, 0.1954, 0.1920, 0.1892};
    gsl_vector *x,*ss;
    x = gsl_vector_alloc(3);
    gsl_vector_set(x, 0, 0.7); //alpha
    gsl_vector_set(x, 1, -0.5);//rho
    gsl_vector_set(x, 2, 0.1);//nu
    
    // step size= 1.0
    ss = gsl_vector_alloc(3);
    gsl_vector_set_all(ss, 0.01);
    
    // initialize method and iterate
    gsl_multimin_function minex_func;
    minex_func.n = 3;
    if (b-1.0<0.0)
        minex_func.f = my_f_3;
    else
        minex_func.f = my_f_2;
    minex_func.params = par;
    
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
    gsl_multimin_fminimizer *s = nullptr;
    s = gsl_multimin_fminimizer_alloc(T, 3);
    gsl_multimin_fminimizer_set(s, &minex_func,x,ss);
    size_t iter = 0;
    int status;
    double size;
    
    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        
        if(status)
            break;
        
        size = gsl_multimin_fminimizer_size(s);
        status = gsl_multimin_test_size(size, 1e-4);
        if (status == GSL_SUCCESS)
        {
            std::cout << "converged to minimum at \n";
        }
        
        cout << "iter = " << iter << ", alpha = " << trans_alpha(gsl_vector_get(s->x, 0))
             << ", beta = " << b
             << ",rho = " << trans_rho(gsl_vector_get(s->x, 1))<< ", nu = "
             << trans_nu(gsl_vector_get(s->x, 2)) << ", f(para) = " << s->fval
             << ", size = " << size <<endl;
    }
    while (status == GSL_CONTINUE && iter < 10000);
    alpha = trans_alpha(gsl_vector_get(s->x, 0));
    beta = b;
    rho = trans_rho(gsl_vector_get(s->x,1));
    nu = trans_nu(gsl_vector_get(s->x,2));
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);
}
