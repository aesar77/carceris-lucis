#include <cmath>
#include <iostream>
#include <fstream>
#include "project.h"
#include "rk45.h"
#include <bits/stdc++.h> 

//Note: The script runs for the 1-a Case, with theta0 = 85 and dx,dy = 0.25. It may take a while, about 10 minutes to run. 

int main() {

double a = 0.99;
double M = 1.0;
double Rin = 5*M;
double Rout = 20*M;
double r_H = M + std::sqrt(M*M - a*a);
std::vector<std::vector<double>> final_pars;
std::vector<std::vector<double>> brightness;

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

  rk45_dormand_prince rk45_dormand_prince(6, 1e-14, 1e-14);

  
  double D = 500.0;
  double xdim = 25.0;
  double ydim = 25.0;
  double bet;
  double alp;
  double dx = 0.25;
  double dy = 0.25;
  //double brightness;
  double omega;
  double xsc = -xdim; 
  double ysc = -ydim;

  for (double i=0; (i*dx < 2*xdim); i++){
    xsc += dx;
    bet = xsc/D;
    std::cout << xsc << std::endl;
    std::vector<double> ind_brightness;
    for (double j=0; (j*dy < 2*ydim); j++){ 
      ysc += dy;
      alp = ysc/D;

    double h = 0.01;
    double x0 = 0.0;
    double r0 = std::sqrt(D*D + std::pow((xsc),2) + std::pow((ysc),2));
    double th0 = (M_PI/180.0) * 85 - alp;
    double phi0 = bet;
    double checkpoint = 0.0;
    metric.compute_metric(r0, th0);
    double ur0 = std::cos(bet) * std::cos(alp) * (-std::sqrt(metric.g_11));
    double uth0 = std::sin(alp) * std::sqrt(metric.g_22);
    double uph0 = std::sin(bet) * std::cos(alp) * std::sqrt(metric.g_33);
    // double rho0 = r0*r0 + a*a * std::pow(std::cos(th0),2);
    // double delta0 = r0*r0 + a*a - 2*M*r0;
    // double sigma0 = std::pow((r0*r0 + a*a),2) - a*a * delta0 * std::pow(std::sin(th0),2);
    // double ur0 = std::cos(bet)*std::cos(alp) * (-std::sqrt(rho0/delta0));
    // double uth0 = std::sin(alp) * std::sqrt(rho0);
    // double uph0 = std::sin(bet)* std::cos(alp) * std::sqrt(sigma0*std::pow(std::sin(th0),2) / rho0);

    std::vector<double> y0  = {r0, th0, phi0, ur0, uth0, uph0};


    auto stop = [Rin, Rout, r0, x0, r_H, phi0](double x, const std::vector<double> &y) { return ((((y[0] > Rin) && (y[0]< Rout)) && (fabs(M_PI/2.0 - y[1]) < 1e-2 )) || y[0] < r_H * 1.02  || y[0] > 6*r0); };


    rk45_dormand_prince.integrate(dydx, stop, h, x0, y0);

   
    for (int i = 0; i < rk45_dormand_prince.result.size(); i++){
      if (((rk45_dormand_prince.result[i][0] > Rin) && (rk45_dormand_prince.result[i][0] < Rout)) && (fabs(M_PI/2.0 - rk45_dormand_prince.result[i][1]) < 1e-2 )){
      
      
      checkpoint = 1.0;
        //std::cout << ind_brightness.size() << std::endl;
      };
    };

    if (checkpoint){
      std::vector<double>f_p = rk45_dormand_prince.result.back(); //getting final values for all parameetrs
      //std::cout << rk45_dormand_prince.result.size() << std::endl;
      metric.compute_metric(f_p[0], f_p[1]);
      double ut_f = std::sqrt(metric.gamma11*f_p[3]*f_p[3] + metric.gamma22*f_p[4]*f_p[4] + metric.gamma33*f_p[5]*f_p[5])/metric.alpha;
      double u_tf = -std::pow(metric.alpha, 2)*ut_f + (-f_p[5]*metric.beta3);

      // double rho_f = f_p[0]*f_p[0] + a*a * std::pow(std::cos(f_p[1]),2);
      // double delta_f = (f_p[0]*f_p[0] + a*a - 2*M*f_p[0]);
      // double sigma_f = std::pow((f_p[0]*f_p[0] + a*a),2) - a*a * delta_f* std::pow(std::sin(f_p[1]),2);
      
      // double g00 = (2*M*f_p[0])/rho_f - 1;
      // double g03 = ((-2*M*a*f_p[0])/rho_f) * std::pow(std::sin(f_p[1]),2);
      // double g33 = (sigma_f/rho_f) * std::pow(std::sin(f_p[1]),2);

      omega = 1.0/(a + (std::pow(f_p[0], 3.0/2.0)/std::sqrt(M)));
      checkpoint = 1.0/(std::pow((1.0 + omega*(-f_p[5]/u_tf))/std::sqrt(-metric.g_00 - omega*omega*metric.g_33 - 2*omega*metric.g_03),3));

    };

    ind_brightness.push_back(checkpoint);

    //std::vector<double>f_p = rk45_dormand_prince.result.back(); //getting final values for all parameetrs
    //std::cout << rk45_dormand_prince.result.size() << std::endl;
    // if (f_p[0]  >= Rout || f_p[0] <= Rin){
    //   brightness = 0;
    // };
    // if ((f_p[0] > Rin) && (f_p[0] < Rout) && (f_p[1] = (M_PI/2))){
    //   metric.compute_metric(f_p[0], f_p[1]);
    //   double ut_f = std::sqrt(metric.gamma11*f_p[3]*f_p[3] + metric.gamma22*f_p[4]*f_p[4] + metric.gamma33*f_p[5]*f_p[5])/metric.alpha;


      // double rho_f = f_p[0]*f_p[0] + a*a * std:pow(std::cos(f_p[1]),2);
      // double delta_f = (f_p[0]*f_p[0] + a*a - 2*M*f_p[0]);
      // double sigma_f = std::pow((f_p[0]*f_p[0] + a*a),2) - a*a * delta_f* std::pow(std::sin(f_p[1]),2);
      
      // double g00 = (2*M*f_p[0])/rho_f - 1;
      // double g03 = ((-2*M*a*f_p[0])/rho_f) * std::pow(std::sin(f_p[1]),2);
      // double g33 = (sigma_f/rho_f) * std::pow(std::sin(f_p[1]),2);

      // omega = 1.0/(a + (std::pow(f_p[0], 3.0/2.0)/std::sqrt(M)));

      // brightness = 1.0 ;
    ;
    // f_p.push_back(xsc);
    // f_p.push_back(ysc);
    //f_p.push_back(brightness);
    //final_pars.push_back(f_p);
    
   };
   brightness.push_back(ind_brightness);
   ysc = -ydim;
   ind_brightness.clear();

  };

// std::ofstream output_file("output4.csv");
// for (int i = 0; i < final_pars.size(); i++) {
//   output_file << i << "," << final_pars[i][0] << "," <<  final_pars[i][2]<<
//           "," <<  final_pars[i][3]<< "," <<  final_pars[i][4]<< "," <<  final_pars[i][5] << "," <<  final_pars[i][6]<< "," <<  final_pars[i][7] << "," <<  final_pars[i][8] << std::endl;
// }
// output_file.close();

std::ofstream output_file("output_raytracing4.csv");
for (int i = 0; i < brightness.size(); i++) {
  for (int j = 0; j < brightness[0].size(); j++){
    if(j < brightness[0].size()-1){
  output_file  << brightness[j][i] << ",";
  }else{
    output_file  << brightness[j][i];
  };
};
output_file << "\n";
};
output_file.close();

// std::ofstream output_file("output5.csv");
// for (int i = 0; i < brightness.size(); i++) {
//   output_file << i << "," << brightness[i][:] << std::endl;
// }
// output_file.close();
return 0;
};

