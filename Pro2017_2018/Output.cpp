//
//  Output.cpp
//  Pro2017_2018
//
//  Created by 任智桂 on 2019/3/17.
//  Copyright © 2019 triple S. All rights reserved.
//

#include "Output.h"
#include "Solver.h"
#include "BSCallClass.h"
#include "BSFormula.h"
#include "SABRClass.h"
#include "SABRFormula.h"
#include "TheFunction.h"
#include "MonteCarlo.h"
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <fstream>
using namespace std;
void output_BSCall_IV()
{
    ofstream outfile1("Newton_Implied_Vol.txt");
    ofstream outfile2("Bisection_Implied_Vol.txt");
    double t_ex = 0.5;
    double f = 80.0;
    double D= 1.0;
    int MaxIter = 1000;
    double epsilon=0.001;
    double Strike[7] = {50.0,60.0,70.0,80.0,90.0,100.0,110.0};
    double MarketPrice[7]={30.0195,20.2170,11.2588,4.7348,1.4439,0.3861,0.1089};
    double Price1[7];
    double Price2[7];
    for (int i=0;i<7;i++)
    {
        //initial guess Vol=5.0;
        double Vol = 1.0;
        double Vol1;
        double Vol2;
        BSCall theCall(f,Strike[i],Vol,D,t_ex);
        Vol1 = Newton(&theCall,MarketPrice[i],Vol,epsilon,MaxIter);
        Vol2 = Bisection(&theCall, MarketPrice[i], 0.0, 5.0, epsilon,MaxIter);
        outfile1<<Vol1<<" ";
        outfile2<<Vol2<<" ";
        cout<<"Vol calculated By Newton\n"<<Vol1<<endl;
        cout<<"Vol calculated By Bisection\n"<<Vol2<<endl;
        Price1[i] = BlackScholesCall(f,Strike[i],Vol1,D,t_ex);
        Price2[i] = BlackScholesCall(f,Strike[i],Vol2,D,t_ex);
        cout<<"MarketPrice vs Price(Newton IV):"<<MarketPrice[i]<<" vs "<<Price1[i]<<endl;
        cout<<"MarketPrice vs Price(Bisection IV):"<<MarketPrice[i]<<" vs "<<Price2[i]<<endl;
    }
    outfile1.close();
    outfile2.close();
}
void output_FitSABR()
{
    double Strike2[7] = {0.01,0.02,0.03,0.04,0.05,0.06,0.07};
    double IV[7] = {0.2295,0.2139,0.2054,0.1996,0.1954,0.1920,0.1892};
    ofstream outfile3("alpha_beta_rho_mu.txt");
    double f = 0.025;
    double t_ex = 1.0;
    double D = 1.0;
    SABR saber(f,D,t_ex);
    SABR saber2(f,D,t_ex);
    SABR saber3(f,D,t_ex);
    saber.Fit_with_all();
    double alpha = saber.get_alpha();
    double beta = saber.get_beta();
    double rho = saber.get_rho();
    double nu = saber.get_nu();
    outfile3<<"alpha = "<<alpha<<" "<<"beta = "<<beta<<" "<<"rho = "<<rho<<" "<<"nu = "<<nu;
    saber2.Fit_with_fixed_beta(1.0);
    alpha = saber2.get_alpha();
    beta = saber2.get_beta();
    rho = saber2.get_rho();
    nu = saber2.get_nu();
    outfile3<<"\n";
    outfile3<<"Optimize when Beta fixed at 1.0"<<"\n";
    outfile3<<"alpha = "<<alpha<<" "<<"beta = "<<beta<<" "<<"rho = "<<rho<<" "<<"nu = "<<nu;
    
    saber3.Fit_with_fixed_beta(0.5);
    alpha = saber3.get_alpha();
    beta = saber3.get_beta();
    rho = saber3.get_rho();
    nu = saber3.get_nu();
    outfile3<<"\n";
    outfile3<<"Optimize when Beta fixed at 0.5"<<"\n";
    outfile3<<"alpha = "<<alpha<<" "<<"beta = "<<beta<<" "<<"rho = "<<rho<<" "<<"nu = "<<nu;
    outfile3.close();
    for (int i=0;i<7;i++)
    {
        double vol=saber.Implied_Vol(Strike2[i]);
        cout<<"True Implied Volatility: "<<IV[i]<<"\n"<<"Calculated Implied Volatility: "<<vol<<endl;
        cout<<endl;
    }
    cout<<"*****************\n";
    for (int i=0;i<7;i++)
    {
        double vol=saber2.Implied_Vol(Strike2[i]);
        cout<<"True Implied Volatility: "<<IV[i]<<"\n"<<"Calculated Implied Volatility: "<<vol<<endl;
        cout<<endl;
    }
    cout<<"*****************\n";
    for (int i=0;i<7;i++)
    {
        double vol=saber3.Implied_Vol(Strike2[i]);
        cout<<"True Implied Volatility: "<<IV[i]<<"\n"<<"Calculated Implied Volatility: "<<vol<<endl;
        cout<<endl;
    }
}
void output_Montecarlo()
{
    ofstream outfile("MonteCarloPrice.txt");
    double f = 1.1;
    double t_ex = 0.5;
    double D = 1.0;
    double sd1 = 0.0;
    double sd2 = 0.0;
    double sd3 = 0.0;
    double K = 1.1;
    int num = 2000;
    SABR saber(f,D,t_ex);
    double alpha = 0.2;
    double beta = 0.7;
    double rho = 0.0;
    double nu = 0.1;
    double sigmaB = SABR_IV(f, K, t_ex, alpha, beta, rho, nu);
    double result1 = MonteCarlo(sd1,num,alpha,
                                beta,rho,nu,K,
                                f,t_ex,D);
    double result2 = Antethetic(sd2, num, alpha,
                                beta,rho, nu, K,
                                f, t_ex, D);
    double result3 = ControlVariates(sd3, num, alpha,
                                     beta, rho, nu, K,
                                     f, t_ex, D);
    cout<<"MonteCarlo's Price:"<<result1<<endl;
    outfile<<"MonteCarlo's Price: "<<result1<<endl;
    cout<<"95% confidence interval is : \n["
        <<(result1-sd1*1.96/num)<<","<<(result1+sd1*1.96/num)<<"]"<<"\n";
    outfile<<"95% confidence interval is : \n["
           <<(result1-sd1*1.96/num)<<","<<(result1+sd1*1.96/num)<<"]"<<"\n";
    cout<<"\n";
    outfile<<"\n";
    cout<<"AntitheticMC's Price: "<<result2<<endl;
    outfile<<"AntitheticMC's Price: "<<result2<<endl;
    cout<<"95% confidence interval is : \n["
        <<(result2-sd2*1.96/num)<<","<<(result2+sd2*1.96/num)<<"]"<<"\n";
    outfile<<"95% confidence interval is : \n["
           <<(result2-sd2*1.96/num)<<","<<(result2+sd2*1.96/num)<<"]"<<"\n";
    cout<<"\n";
    outfile<<"\n";
    cout<<"ControlVariatesMC's Price: "<<result3<<endl;
    outfile<<"ControlVariatesMC's Price: "<<result3<<endl;
    cout<<"95% confidence interval is : \n["
        <<(result3-sd3*1.96/num)<<","<<(result3+sd3*1.96/num)<<"]"<<"\n";
    outfile<<"95% confidence interval is : \n["
           <<(result3-sd3*1.96/num)<<","<<(result3+sd3*1.96/num)<<"]"<<"\n";
    cout<<"\n";
    outfile<<"\n";
    cout<<"Pricing Formula's Price:"<<SABRCall(f,K, D,
                                               t_ex,alpha,beta,rho,nu,sigmaB)<<endl;
    outfile<<"Pricing Formula's Price:"<<SABRCall(f,K, D,
                                                  t_ex,alpha,beta,rho,nu,sigmaB)<<endl;
}
