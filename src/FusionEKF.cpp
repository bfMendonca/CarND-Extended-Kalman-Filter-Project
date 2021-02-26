#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;


  VectorXd xInit(4);  
  MatrixXd pInit(4, 4);
  pInit.setIdentity();

  pInit(0,0) = 0.15*0.15;
  pInit(1,1) = 0.15*0.15;
  pInit(2,2) = 9.0;
  pInit(3,3) = 9.0;


  ekf_.Init( xInit, pInit );
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  
  /**
   * Initialization
   */
  if (!is_initialized_ ) {

      static bool first = true;
      static int counter = 0;
      static double initPx = 0;
      static double initPy = 0;
      static double initVx = 0;
      static double initVy = 0;

      static double totalTime = 0.0;

      static double displacementX = 0.0;
      static double displacementY = 0.0;

      static double lastPx;
      static double lastPy;

      double newx; 
      double newy;


      //Acquiring measurements for estimating the intiial position and velocity of the target

      if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)  {
          newx = measurement_pack.raw_measurements_(0);
          newy = measurement_pack.raw_measurements_(1);
      }else if( measurement_pack.sensor_type_ == MeasurementPackage::RADAR ) {
          newx = measurement_pack.raw_measurements_(0)*cos(measurement_pack.raw_measurements_(1));
          newy = measurement_pack.raw_measurements_(0)*sin(measurement_pack.raw_measurements_(2));

      }

      double dt = 0;
      if( !first ) {
          dt = (measurement_pack.timestamp_ - previous_timestamp_)/1e6;
          double dx = newx - lastPx;
          double dy = newy - lastPy;
          
          displacementX += dx;
          displacementY += dy;
          totalTime += dt;
      }
      first = false;

      
      previous_timestamp_ = measurement_pack.timestamp_;
      lastPx = newx;
      lastPy = newy;
      
      if( ++counter == 3 ) {
          double initVx = displacementX/totalTime;
          double initVy = displacementY/totalTime;
  
          ekf_.x_(0) = lastPx;
          ekf_.x_(1) = lastPy;
          ekf_.x_(2) = initVx;
          ekf_.x_(3) = initVy;

          std::cout << ekf_.x_ << std::endl;

          is_initialized_ = true;
      }
                
  }

  if( !is_initialized_ ) {
      return;
  }

  /**
   * Prediction
   */
  const double dt = (measurement_pack.timestamp_ - previous_timestamp_)/1e6;
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.F_ = tools.CalculateStateTransitionMatrix( dt );
  ekf_.Q_ = tools.CalculateProccessNosieCovMatrix( dt ); 

  ekf_.Predict( dt );
  
  /**
   * Update
   */
  if ( measurement_pack.sensor_type_ == MeasurementPackage::RADAR ) {
    Hj_ = tools.CalculateJacobian( ekf_.x_);
    ekf_.UpdateEKF( Hj_, measurement_pack.raw_measurements_, R_radar_ );

  } else if( measurement_pack.sensor_type_ == MeasurementPackage::LASER ) {
    ekf_.Update( H_laser_, measurement_pack.raw_measurements_, R_laser_ );
  }
}
