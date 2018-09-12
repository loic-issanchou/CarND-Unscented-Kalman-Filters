#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI/8;//0.3925
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;

  //number of state vector elements 
  n_x_ = 5;

  //number of augmented state vector elements 
  n_aug_ = 7;

  //Matrix of sigma points predictions
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  time_us_ = 0.0;

  weights_ = VectorXd(2*n_aug_+1);

  lambda_ = 3 - n_x_;

  i = 0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_){
    
    x_ << 0.5, 0.5, 1, 1, 0.4;
  
    P_ << 0.15,    0, 0, 0, 0,
	    0, 0.15, 0, 0, 0, 
	    0,    0, 1, 0, 0,
	    0,    0, 0, 1, 0,
	    0,    0, 0, 0, 0.5;

    use_laser_ = true;
    use_radar_ = true;

    if (meas_package.sensor_type_ = MeasurementPackage::LASER){

      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);

    }else if (meas_package.sensor_type_ = MeasurementPackage::RADAR){

      float rho = meas_package.raw_measurements_(0);
      float phi = meas_package.raw_measurements_(1);

      x_(0) = rho * cos(phi);
      x_(1) = rho * sin(phi);
    }
    
    // set weights
    double weight_0_ = lambda_/(lambda_+n_aug_);
    weights_(0) = weight_0_;
    for (int i=1; i<2*n_aug_+1; i++) {  
      double weight_ = 0.5/(n_aug_+lambda_);
      weights_(i) = weight_;
    }

    time_us_ = meas_package.timestamp_;

    is_initialized_ = true;
    return;
  }

  //Calcul of time between two measurements
  float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
  time_us_ = meas_package.timestamp_;

  i += 1;
  cout << "i = " << i << "\n";

  Prediction(delta_t);
  
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
    UpdateRadar(meas_package);
  }else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
    UpdateLidar(meas_package);
  }

}

////////////////////
//// PREDICTION ////
////////////////////

//-- Generate Sigma Points - FUNCTION --//
void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out_) {

  //create sigma point matrix
  MatrixXd Xsig_ = MatrixXd(n_x_, 2 * n_x_ + 1);

  //calculate square root of P
  MatrixXd A_ = P_.llt().matrixL();


  //set first column of sigma point matrix
  Xsig_.col(0)  = x_;

  //set remaining sigma points
  for (int i = 0; i < n_x_; i++)
  {
    Xsig_.col(i+1)      = x_ + sqrt(lambda_+n_x_) * A_.col(i);
    Xsig_.col(i+1+n_x_) = x_ - sqrt(lambda_+n_x_) * A_.col(i);
  }

  *Xsig_out_ = Xsig_;

}

//-- Augmented Sigma Points - FUNCTION --//
void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out_) {

  //create augmented mean vector
  VectorXd x_aug_ = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug_.head(n_x_) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;

  //create augmented covariance matrix
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(n_x_,n_x_) = P_;
  P_aug_(5,5) = std_a_*std_a_;
  P_aug_(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L_ = P_aug_.llt().matrixL();

  //create augmented sigma points
  Xsig_aug_.col(0)  = x_aug_;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug_.col(i+1)       = x_aug_ + sqrt(lambda_+n_aug_) * L_.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_+n_aug_) * L_.col(i);
  }

  //write result
  *Xsig_out_ = Xsig_aug_;
}

//-- Sigma Point Prediction - FUNCTION --//
void UKF::SigmaPointPrediction(MatrixXd Xsig_aug_, double delta_t) {

  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x_      = Xsig_aug_(0,i);
    double p_y_      = Xsig_aug_(1,i);
    double v_        = Xsig_aug_(2,i);
    double yaw_      = Xsig_aug_(3,i);
    double yawd_     = Xsig_aug_(4,i);
    double nu_a_     = Xsig_aug_(5,i);
    double nu_yawdd_ = Xsig_aug_(6,i);

    //predicted state values
    double px_p_, py_p_;


    //avoid division by zero
    if (fabs(yawd_) > 0.001) {
        px_p_ = p_x_ + v_/yawd_ * ( sin(yaw_ + yawd_*delta_t) - sin(yaw_));
        py_p_ = p_y_ + v_/yawd_ * ( cos(yaw_) - cos(yaw_+yawd_*delta_t) );
    }
    else {
        px_p_ = p_x_ + v_*delta_t*cos(yaw_);
        py_p_ = p_y_ + v_*delta_t*sin(yaw_);
    }

    double v_p_ = v_;
    double yaw_p_ = yaw_ + yawd_*delta_t;
    double yawd_p_ = yawd_;

    //add noise
    px_p_ = px_p_ + 0.5*nu_a_*delta_t*delta_t * cos(yaw_);
    py_p_ = py_p_ + 0.5*nu_a_*delta_t*delta_t * sin(yaw_);
    v_p_ = v_p_ + nu_a_*delta_t;

    yaw_p_ = yaw_p_ + 0.5*nu_yawdd_*delta_t*delta_t;
    yawd_p_ = yawd_p_ + nu_yawdd_*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p_;
    Xsig_pred_(1,i) = py_p_;
    Xsig_pred_(2,i) = v_p_;
    Xsig_pred_(3,i) = yaw_p_;
    Xsig_pred_(4,i) = yawd_p_;
  }
}

//-- Predict Mean And Covariance - FUNCTION --//
void UKF::PredictMeanAndCovariance(VectorXd weights_) {

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff_ = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff_(3)> M_PI) x_diff_(3)-=2.*M_PI;
    while (x_diff_(3)<-M_PI) x_diff_(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff_ * x_diff_.transpose() ;
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  MatrixXd Xsig_ = MatrixXd(n_x_, 2*n_x_+1);
  GenerateSigmaPoints(&Xsig_);

  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2*n_aug_+1);
  AugmentedSigmaPoints(&Xsig_aug_);

  SigmaPointPrediction(Xsig_aug_, delta_t);

  PredictMeanAndCovariance(weights_);
}

////////////////////
////// UPDATE //////
////////////////////

//-- Update state - Same function for both Lidar or Radar - FUNCTION --//
void UKF::UpdateState(MeasurementPackage meas_package, MatrixXd Zsig_, VectorXd z_pred_, MatrixXd S_, double* NIS_, VectorXd z_, int n_z_) {
 

  //create matrix for cross correlation Tc
  MatrixXd Tc_ = MatrixXd(n_x_, n_z_);

  //calculate cross correlation matrix
  Tc_.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff_ = Zsig_.col(i) - z_pred_;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
      //angle normalization
      while (z_diff_(1)> M_PI) z_diff_(1)-=2.*M_PI;
      while (z_diff_(1)<-M_PI) z_diff_(1)+=2.*M_PI;
    }

    // state difference
    VectorXd x_diff_ = Xsig_pred_.col(i) - x_;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
      //angle normalization
      while (x_diff_(3)> M_PI) x_diff_(3)-=2.*M_PI;
      while (x_diff_(3)<-M_PI) x_diff_(3)+=2.*M_PI;
    }

    Tc_ = Tc_ + weights_(i) * x_diff_ * z_diff_.transpose();
  }

  //Kalman gain K;
  MatrixXd K_ = Tc_ * S_.inverse();

  //residual
  VectorXd z_diff_ = z_ - z_pred_;

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
    //angle normalization
    while (z_diff_(1)> M_PI) z_diff_(1)-=2.*M_PI;
    while (z_diff_(1)<-M_PI) z_diff_(1)+=2.*M_PI;
  }

  //update state mean and covariance matrix
  x_ = x_ + K_ * z_diff_;
  P_ = P_ - K_*S_*K_.transpose();

  // Calcul of NIS 
  *NIS_ = z_diff_.transpose() * S_.inverse() * z_diff_;
}


//-- Predict Lidar Measurement - FUNCTION --//
void UKF::PredictLidarMeasurement(VectorXd* z_out_, MatrixXd* S_out_, MatrixXd* Zsig_){

  //set measurement dimension, lidar can measure px and py
  int n_z_ = 2;

  MatrixXd Zsig_temp_ = MatrixXd(2, 2 * n_aug_ + 1);
  MatrixXd S_ = MatrixXd(2, 2);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x_ = Xsig_pred_(0,i);
    double p_y_ = Xsig_pred_(1,i);

    // measurement model
    Zsig_temp_(0,i) = p_x_;    //p_x_
    Zsig_temp_(1,i) = p_y_;    //p_y_
  }

  //mean predicted measurement
  VectorXd z_pred_ = VectorXd(n_z_);
  z_pred_.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred_ = z_pred_ + weights_(i) * Zsig_temp_.col(i);
  }

  //innovation covariance matrix S
  S_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff_ = Zsig_temp_.col(i) - z_pred_;

    S_ = S_ + weights_(i) * z_diff_ * z_diff_.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R_ = MatrixXd(n_z_,n_z_);
  R_ <<   std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;

  S_ = S_ + R_;

  *z_out_ = z_pred_;
  *S_out_ = S_;
  *Zsig_ = Zsig_temp_;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  VectorXd z_pred_ = VectorXd(2);
  MatrixXd S_ = MatrixXd(2, 2);
  //create matrix for sigma points in measurement space
  MatrixXd Zsig_ = MatrixXd(2, 2 * n_aug_ + 1);

  PredictLidarMeasurement(&z_pred_, &S_, &Zsig_);

  double NIS_lidar_;
  VectorXd z_ = VectorXd(2);   
  z_ = meas_package.raw_measurements_;
  //set measurement dimension, radar can measure p_x_ and p_y_
  int n_z_ = 2;
  UpdateState(meas_package, Zsig_, z_pred_, S_, &NIS_lidar_, z_, n_z_);

  cout << "NIS lidar = " << NIS_lidar_ << "\n";

}

//-- Predict Radar Measurement - FUNCTION --//
void UKF::PredictRadarMeasurement(VectorXd* z_out_, MatrixXd* S_out_, MatrixXd* Zsig_) {

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z_ = 3;

  MatrixXd Zsig_temp_ = MatrixXd(3, 2 * n_aug_ + 1);
  MatrixXd S_ = MatrixXd(3, 3);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x_ = Xsig_pred_(0,i);
    double p_y_ = Xsig_pred_(1,i);
    double v_  = Xsig_pred_(2,i);
    double yaw_ = Xsig_pred_(3,i);

    double v1_ = cos(yaw_)*v_;
    double v2_ = sin(yaw_)*v_;
    
    // measurement model
    Zsig_temp_(0,i) = sqrt(p_x_*p_x_ + p_y_*p_y_);                          //r
    Zsig_temp_(1,i) = atan2(p_y_,p_x_);                                     //phi
    Zsig_temp_(2,i) = (p_x_*v1_ + p_y_*v2_ ) / sqrt(p_x_*p_x_ + p_y_*p_y_); //r_dot
  }
  
  //mean predicted measurement
  VectorXd z_pred_ = VectorXd(n_z_);
  z_pred_.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred_ = z_pred_ + weights_(i) * Zsig_temp_.col(i);
  }

  //innovation covariance matrix S
  S_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff_ = Zsig_temp_.col(i) - z_pred_;
    
    //angle normalization
    while (z_diff_(1)> M_PI) z_diff_(1)-=2.*M_PI;
    while (z_diff_(1)<-M_PI) z_diff_(1)+=2.*M_PI;
    
    S_ = S_ + weights_(i) * z_diff_ * z_diff_.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R_ = MatrixXd(n_z_,n_z_);
  R_ <<   std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S_ = S_ + R_;

  *z_out_ = z_pred_;
  *S_out_ = S_;
  *Zsig_ = Zsig_temp_;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  VectorXd z_pred_ = VectorXd(3);
  MatrixXd S_ = MatrixXd(3, 3);
  //create matrix for sigma points in measurement space
  MatrixXd Zsig_ = MatrixXd(3, 2 * n_aug_ + 1);

  PredictRadarMeasurement(&z_pred_, &S_, &Zsig_);

  double NIS_radar_;
  VectorXd z_ = VectorXd(3);
  z_ = meas_package.raw_measurements_;
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z_ = 3;
  UpdateState(meas_package, Zsig_, z_pred_, S_, &NIS_radar_, z_, n_z_);

  cout << "NIS radar = " << NIS_radar_ << "\n";
}
