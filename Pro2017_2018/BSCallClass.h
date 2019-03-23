//
//  BSCallClass.h
//  Pro2017_2018
//
//  Created by 任智桂 on 2019/3/12.
//  Copyright © 2019 triple S. All rights reserved.
//

#ifndef BSCallClass_h
#define BSCallClass_h
class BSCall
{
public:
    BSCall(double f_,double K_,double sigma_,double D_,double t_ex_);
    double Pricing(double Vol) const;
    double Vega(double Vol) const;
private:
    double f;
    double K;
    double sigma;
    double D;
    double t_ex;
};

#endif /* BSCallClass_h */
