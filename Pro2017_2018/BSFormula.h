//
//  BlackScholesFormulas.h
//  Pro2017_2018
//
//  Created by 任智桂 on 2019/3/12.
//  Copyright © 2019 triple S. All rights reserved.
//

#ifndef BSFormula_h
#define BSFormula_h
double BlackScholesCall(double f,
                        double K,
                        double sigma,
                        double D,
                        double t_ex);
double BlackScholesVega(double f,
                        double K,
                        double sigma,
                        double D,
                        double t_ex);
#endif /* BSFormula_h */
