//
//  BSCallClass.cpp
//  Pro2017_2018
//
//  Created by 任智桂 on 2019/3/12.
//  Copyright © 2019 triple S. All rights reserved.
//

#include "BSCallClass.h"
#include "BSFormula.h"
BSCall::BSCall(double f_,double K_,
               double sigma_,double D_,double t_ex_):f(f_),K(K_),sigma(sigma_),D(D_),t_ex(t_ex_)
{}
double BSCall::Pricing(double Vol) const
{
    return BlackScholesCall(f,K,Vol,D, t_ex);
}
double BSCall::Vega(double Vol) const
{
    return BlackScholesVega(f,K,Vol,D,t_ex);
}


