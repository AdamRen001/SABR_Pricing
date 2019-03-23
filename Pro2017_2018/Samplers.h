//
//  Sampler.h
//  Pro2017_2018
//
//  Created by 任智桂 on 2019/3/14.
//  Copyright © 2019 triple S. All rights reserved.
//

#ifndef Samplers_h
#define Samplers_h
class Sampler{
public:
    virtual double getnumber()=0;
};
class Uniform:public Sampler{
public:
    Uniform(double aa,double bb){a = aa;b = bb;}
    virtual double getnumber();
private:
    double a,b;
};
class InverseFuncSample:public Sampler{
public:
    InverseFuncSample(double L=1.0){lambda=L;}
    virtual double getnumber();
private:
    double lambda;
};
class BoxMuller:public Sampler{
public:
    BoxMuller(){};
    virtual double getnumber();
};
class DoubleExp:public Sampler{
public:
    DoubleExp(){};
    virtual double getnumber();
};
class Rej_Samp:public Sampler{
public:
    Rej_Samp(double cc,DoubleExp p){C=cc;pt=&p;}
    virtual double getnumber();
private:
    double C=1.3155;
    DoubleExp *pt;
};
#endif /* Samplers_h */
