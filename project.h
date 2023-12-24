#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

class boyer_lindquist_metric{
    public:
    boyer_lindquist_metric(double a0, double M0) {
    // Initialize the parameters a and M
    a = a0;
    M = M0;
    }
    void compute_metric(double r, double th) {
    // Compute all the required metric components in one go

    
    rho = r*r + a*a*std::pow(std::cos(th),2);
    drho_dr = 2*r;
    drho_dth = a*a * 2 * std::cos(th) * (-std::sin(th));

    delta = (r*r + a*a - 2*M*r);
    ddelta_dr = 2*r - 2*M;
    ddelta_dth = 0;

    sigma = std::pow((r*r + a*a),2) - a*a * delta * std::pow(std::sin(th),2);
    dsigma_dr = 2*(r*r + a*a)*2*r - a*a*(2*r - 2*M)*std::pow(std::sin(th),2);
    dsigma_dth = - a*a * delta * 2 * std::sin(th) * std::cos(th);


    alpha = std::sqrt((rho*delta)/sigma);
    d_alpha_dr = 0.5 * std::pow((rho*delta)/sigma, -1.0/2.0) * ((delta/sigma)*drho_dr + (rho/sigma)*ddelta_dr + rho*delta* (-1.0 * std::pow(sigma, -2) * dsigma_dr));
    d_alpha_dth = 0.5 * std::pow((rho*delta)/sigma, -1.0/2.0) * ((delta/sigma)*drho_dth + (rho/sigma)*ddelta_dth + rho*delta* (-1.0 * std::pow(sigma, -2) * dsigma_dth));

    beta3 = (-2*M*a*r)/sigma;
    d_beta3_dr =  (-2*M*a)/sigma + (-2*M*a*r)*(-1.0 * std::pow(sigma, -2) * dsigma_dr);
    d_beta3_dth = (-2*M*a*r)*(-1.0 * std::pow(sigma, -2) * dsigma_dth);

    gamma11 = delta/rho;
    d_gamma11_dr = ddelta_dr/rho + delta*((-1.0 * std::pow(rho, -2) * drho_dr));
    d_gamma11_dth = ddelta_dth/rho + delta*((-1.0 * std::pow(rho, -2) * drho_dth));

    gamma22 = 1.0/rho;
    d_gamma22_dr = (-1.0 * std::pow(rho, -2) * drho_dr);
    d_gamma22_dth = (-1.0 * std::pow(rho, -2) * drho_dth);

    gamma33 = rho/(sigma*std::pow(std::sin(th),2));
    d_gamma33_dr = drho_dr/(sigma*std::pow(std::sin(th),2)) + (rho/std::pow(std::sin(th),2)) * (-1 * std::pow(sigma, -2) * dsigma_dr);
    d_gamma33_dth =  drho_dth/(sigma*std::pow(std::sin(th),2)) + (rho/std::pow(std::sin(th),2)) * (-1 * std::pow(sigma, -2) * dsigma_dth) + (rho/sigma)* (-2 * std::pow(std::sin(th), -3) * std::cos(th));

    g_00 = (2*M*r)/rho - 1;
    g_03 = ((-2*M*a*r)/rho) * std::pow(std::sin(th),2);
    g_11 = rho/delta;
    g_22 = rho;
    g_33 = (sigma/rho) * std::pow(std::sin(th),2);

    };
    double rho, delta, sigma, ut;
    double drho_dr, drho_dth, ddelta_dr, ddelta_dth, dsigma_dr, dsigma_dth;
    double a, M;
    double alpha, beta3;
    double gamma11, gamma22, gamma33; // components of upper gamma^ij
    double g_00, g_03, g_11, g_22, g_33; // components of lower g_\mu\nu
    double d_alpha_dr, d_beta3_dr, d_gamma11_dr, d_gamma22_dr, d_gamma33_dr;
    double d_alpha_dth, d_beta3_dth, d_gamma11_dth, d_gamma22_dth, d_gamma33_dth;
};