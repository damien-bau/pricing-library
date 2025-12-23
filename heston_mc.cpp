/*
Author: Damien Bau
Description: Heston Stochastic Volatility Model option pricer using Monte Carlo simulations 
and Full Truncation Euler Discretisation
Date: October 2025
*/

#include <iostream>
#include <random>
#include <cmath>

using std::cout;
using std::string;
using std::max;
using std::sqrt;
using std::exp;





// Heston model parameters
struct Parameters{
    double kappa, theta, xi, rho;
    double v0, s0, r, T, K;
    int steps, paths;
    string option_type;
};

void generateCorrelatedNormals(double rho, double& X1, double& X2){
// Generates two correlated N(0,1) draws with correlation rho
    static std::mt19937 rng(std::random_device{}());
    static std::normal_distribution<double> norm(0., 1.);


    double Z1 = norm(rng);
    double Z2 = norm(rng);

    X1 = Z1;
    X2 = rho * Z1 + sqrt(1. - pow(rho,2.))*Z2;

}


double priceOption(const Parameters& par){

    double dt = par.T / par.steps; // time increments
    double payoff_sum = 0.;

    for (int i=0; i<par.paths; i++){
        double S  = par.s0;
        double v = par.v0;
        double K = par.K;

        for (int j=0; j< par.steps; j++){
            double X1, X2;
            generateCorrelatedNormals(par.rho, X1, X2);
            double v_max = max(v,0.);
            v += par.kappa*(par.theta - v)*dt + par.xi * sqrt(v_max)*X2*sqrt(dt);
            S *= exp((par.r - 0.5*v_max)*dt + sqrt(v_max)*X1*sqrt(dt));
        }

        double payoff;
        if (par.option_type == "call"){
            payoff = max(S-K, 0.);
        }
        else if (par.option_type == "put")
        {
            payoff = max(K-S, 0.);
        }
        else {
            cout << "Option type isn't valid." << '\n';
        }

        payoff_sum += payoff;
        
    }

    return exp(-par.r*par.T)*payoff_sum/par.paths;
}


int main(){

    // Define parameters
    Parameters parameters;
    parameters.kappa = 1.5;
    parameters.theta = 0.04;
    parameters.xi = 0.;
    parameters.rho = -0.6;
    parameters.v0 = 0.04;
    parameters.s0 = 100.0;
    parameters.r = 0.01;
    parameters.T = 1.0;
    parameters.K = 100.0;
    parameters.steps = 252;
    parameters.paths = 100000;
    parameters.option_type = "call";

    double price = priceOption(parameters);

    cout << "OPTION PRICE: " << price << '\n';

    return 0;
}