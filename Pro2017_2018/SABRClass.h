//
//  SABRClass.h
//  Pro2017_2018
//
//  Created by 任智桂 on 2019/3/12.
//  Copyright © 2019 triple S. All rights reserved.
//

#ifndef SABRClass_h
#define SABRClass_h
class SABR{
public:
    SABR(double f_,double D_,double t_ex_);
    double Implied_Vol(double K_) const;
    double get_alpha();
    double get_beta();
    double get_rho();
    double get_nu();
    void Fit_with_fixed_beta(double b);
    void Fit_with_all();
private:
    double f;
    double K;
    double D;
    double t_ex;
    double alpha;
    double beta;
    double rho;
    double nu;
    double sigma;
};

#endif /* SABRClass_h */
