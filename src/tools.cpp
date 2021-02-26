#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

MatrixXd Tools::CalculateStateTransitionMatrix( const double & dt ) {
  MatrixXd f(4,4);

  f <<  1, 0, dt, 0, 
        0, 1, 0, dt,
        0, 0, 1, 0,
        0, 0, 0, 1;

  return f;
}


MatrixXd Tools::CalculateProccessNosieCovMatrix( const double & dt ) {
  MatrixXd q(4,4);

  const double a = pow(dt,4)/4.0;
  const double b = pow(dt,2)/2.0;
  const double c = pow(dt,3)/2.0;
  const double d = pow(dt,2);
  
  constexpr double ax2 = pow(9.0, 2.0);
  constexpr double ay2 = pow(9.0, 2.0);

  q <<  a*ax2, 0, c*ax2, 0,
        0, a*ay2, 0, c*ay2,
        c*ax2, 0, d*ax2, 0,
        0, c*ay2, 0, d*ay2;

  return q;

}


VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
    const vector<VectorXd> &ground_truth) {

  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() != ground_truth.size()
      || estimations.size() == 0) {
    std::cout << "Invalid estimation or ground_truth data" << std::endl;
    return rmse;
  }

  // accumulate squared residuals
  for (unsigned int i=0; i < estimations.size(); ++i) {

    VectorXd residual = estimations[i] - ground_truth[i];

    // coefficient-wise multiplication
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  // calculate the mean
  rmse = rmse/estimations.size();

  // calculate the squared root
  rmse = rmse.array().sqrt();

  // return the result
  return rmse;
}


MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // pre-compute a set of terms to avoid repeated calculation
  float c1 = px*px+py*py;
  float c2 = sqrt(c1);
  float c3 = (c1*c2);

  // check division by zero
  if (fabs(c1) < 0.0001) {
    std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
    return Hj;
  }

  // compute the Jacobian matrix
  Hj <<   (px/c2), (py/c2), 0, 0,
          -(py/c1), (px/c1), 0, 0,
          py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

  return Hj;
}
