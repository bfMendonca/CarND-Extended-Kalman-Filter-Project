#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in ) {
  x_ = x_in;
  P_ = P_in;

  n_x_ = 4;
  n_z_ = 2;

  F_ = MatrixXd( n_x_, n_x_ );
  Q_ = MatrixXd( n_x_, n_x_ );
}

void KalmanFilter::Predict( const double & dt ) {
    x_ = F_*x_;
    P_ = F_*P_*F_.transpose() + Q_;
}

void KalmanFilter::Update( const MatrixXd &h, const VectorXd &z, const MatrixXd &r ) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  MatrixXd S;
  S = h*P_*h.transpose() + r;

    // calculate Kalman gain K;
  MatrixXd K( n_x_, n_z_ );
  K = P_*h.transpose()*S.inverse();

  // update state mean and covariance matrix
  VectorXd z_pred = h*x_;

  x_ = x_ + K * ( z - z_pred );

  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * h) * P_;
}

void KalmanFilter::UpdateEKF( const MatrixXd &h, const VectorXd &z, const MatrixXd &r ) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */

  MatrixXd S;
  S = h*P_*h.transpose() + r;

    // calculate Kalman gain K;
  MatrixXd K( n_x_, n_z_ );
  K = P_*h.transpose()*S.inverse();

  // update state mean and covariance matrix
  VectorXd z_pred(3);

  double px = x_(0);
  double py = x_(1);
  double vx = x_(2);
  double vy = x_(3);

  z_pred(0) = sqrt( pow( px, 2.0 ) + pow( py, 2.0 ) );
  z_pred(1) = atan2( py, px );
  z_pred(2) = (px*vx + py*vy)/z_pred(0);


  VectorXd y( ( z - z_pred ) );

  // angle normalization
  while ( y(1)> M_PI) y(1)-=2.*M_PI;
  while ( y(1)<-M_PI) y(1)+=2.*M_PI;

  x_ = x_ + K * y;

  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * h) * P_;
  
}
