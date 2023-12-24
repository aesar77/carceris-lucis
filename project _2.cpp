#include <cmath>
#include <iostream>
#include <fstream>
#include "project.h"
#include "rk45.h"

int main() {

double a = 1.0;
double M = 1.0;

boyer_lindquist_metric metric(a,M);

auto dydx = [&metric] (double x, const std::vector<double> &y) {
metric.compute_metric(y[0], y[1]); // y[0] is r, y[1] is theta
// compute the right hand sides of equations (17)-(22)
     double ut = std::sqrt(metric.gamma11*y[3]*y[3] + metric.gamma22*y[4]*y[4] + metric.gamma33*y[5]*y[5])/metric.alpha;

     return std::vector<double>{
        metric.gamma11 * y[3]/ut,
        metric.gamma22 * y[4]/ut,
        metric.gamma33 * y[5]/ut - metric.beta3,
        -metric.alpha * ut * metric.d_alpha_dr + y[5] * metric.d_beta3_dr - (1.0/(2*ut)) *(y[3]*y[3]*metric.d_gamma11_dr + y[4]*y[4]*metric.d_gamma22_dr + y[5]*y[5]*metric.d_gamma33_dr),
        -metric.alpha * ut * metric.d_alpha_dth + y[5] * metric.d_beta3_dth - (1.0/(2*ut)) *(y[3]*y[3]*metric.d_gamma11_dth + y[4]*y[4]*metric.d_gamma22_dth + y[5]*y[5]*metric.d_gamma33_dth),
         0};
};

  rk45_dormand_prince rk45_dormand_prince(6, 1e-10, 1e-10);

  double h = 0.01;
  double x0 = 0.0;
  double r0 = 2;
  double th0 = M_PI/2.0;
  double phi0 = 0;
  double rho0 = r0*r0 + a*a * std::pow(std::cos(th0),2);
  double delta0 = r0*r0 + a*a - 2*M*r0;
  double sigma0 = std::pow((r0*r0 + a*a),2) - a*a * delta0 * std::pow(std::sin(th0),2);
  double ur0 = 0;
  double uth0 = std::sqrt(16);
  double uph0 = 1;

  std::vector<double> y0  = {r0, th0, phi0, ur0, uth0, uph0};

  auto stop = [r0](double x, const std::vector<double> &y) { return fabs(y[0] - r0) >= 0.9 ; };

  rk45_dormand_prince.integrate(dydx, stop, h, x0, y0);

      std::ofstream output_file("output2E.csv");
  for (int i = 0; i < rk45_dormand_prince.xs.size(); i++) {
    output_file << rk45_dormand_prince.xs[i] << "," << rk45_dormand_prince.result[i][0] << "," << rk45_dormand_prince.result[i][1] << "," << rk45_dormand_prince.result[i][2]<<
             "," << rk45_dormand_prince.result[i][3]<< "," << rk45_dormand_prince.result[i][4]<< "," << rk45_dormand_prince.result[i][5] << std::endl;
  }
  output_file.close();


};



