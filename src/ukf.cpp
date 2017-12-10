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
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;
  
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

  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  //time when the state is true, in us
  time_us_ = 0;

  // State dimension
  n_x_ = 5;
  // Augmented state dimension
  n_aug_ = 7;

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  weights_.tail(2*n_aug_).fill(0.5 / (lambda_+n_aug_));
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

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/

   if (!is_initialized_) {
    // first measurement
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates and initialize state
      float rho0 = meas_package.raw_measurements_(0);
      float phi0 = meas_package.raw_measurements_(1);
      float px0 = rho0 * cos(phi0);
      float py0 = rho0 * sin(phi0);
      x_ << px0, py0, 0, 0, 0;
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // initialize state
      float px0 = meas_package.raw_measurements_(0);
      float py0 = meas_package.raw_measurements_(1);
      x_ << px0, py0, 0, 0, 0;
    }

    // initialize state covariance matrix
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }


  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  // calculate the elapsed time
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  // predict the state and the state covariance matrix
  Prediction(dt);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    // Radar update
    UpdateRadar(meas_package);
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    // Laser update
    UpdateLidar(meas_package);
  }

  // print the output
  // cout << "x_ = " << x_ << endl;
  // cout << "P_ = " << P_ << endl;
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

  /*****************************************************************************
   *  Generate Sigma Points
   ****************************************************************************/

  // augmented mean state
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(n_x_) = x_;
  x_aug.tail(2).fill(0.0);

  // augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_ * std_a_;
  P_aug(n_x_+1, n_x_+1) = std_yawdd_ * std_yawdd_;

  // square root of covariance matrix
  MatrixXd L = P_aug.llt().matrixL();

  // augmented sigma points
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.col(0) = x_aug;
  Xsig_aug.block(0, 1, n_aug_, n_aug_) = x_aug.replicate(1, n_aug_)
                                         + sqrt(lambda_ + n_aug_) * L;
  Xsig_aug.block(0, n_aug_+1, n_aug_, n_aug_) = x_aug.replicate(1, n_aug_)
                                                - sqrt(lambda_ + n_aug_) * L;

  /*****************************************************************************
   *  Predict Sigma Points
   ****************************************************************************/

  double delta_t2 = delta_t * delta_t;

  for(int i = 0; i < 2*n_aug_ + 1; i++) {
    VectorXd x_aug = Xsig_aug.col(i);
    VectorXd f_k(n_x_);
    VectorXd nu_k(n_x_);

    // calculate f_k
    if (x_aug(4) > 0.001) {
        f_k <<  x_aug(2)/x_aug(4) * ( sin(x_aug(3)+x_aug(4)*delta_t) - sin(x_aug(3)) ),
                x_aug(2)/x_aug(4) * (-cos(x_aug(3)+x_aug(4)*delta_t) + cos(x_aug(3)) ),
                0,
                x_aug(4) * delta_t,
                0;
    } else { //avoid division by zero
        f_k <<  x_aug(2) * cos(x_aug(3)) * delta_t,
                x_aug(2) * sin(x_aug(3)) * delta_t,
                0,
                x_aug(4) * delta_t,
                0;
    }

    // calculate nu_k
    nu_k << 0.5 * delta_t2 * cos(x_aug(3)) * x_aug(5),
            0.5 * delta_t2 * sin(x_aug(3)) * x_aug(5),
            delta_t * x_aug(5),
            0.5 * delta_t2 * x_aug(6),
            delta_t * x_aug(6);

    // calculate predicted sigma point
    Xsig_pred_.col(i) = x_aug.head(n_x_) + f_k + nu_k;
  }

  /*****************************************************************************
   *  Predict Mean and Covariance
   ****************************************************************************/

  // predict mean state
  x_ = Xsig_pred_ * weights_;

  // predict state covariance
  P_.fill(0.0);
  for(int i = 0; i < 2*n_aug_ + 1; i++) {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
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

  /*****************************************************************************
   *  Predict Lidar Measurement
   ****************************************************************************/

  // measurement dimension - laser
  int n_z = 2;

  // measurement matrix  - laser
  MatrixXd H_laser = MatrixXd(n_z, n_x_);
  H_laser << 1, 0, 0, 0, 0,
             0, 1, 0, 0, 0;

  // transform sigma points into measurement space
  MatrixXd Zsig = H_laser * Xsig_pred_;

  // calculate mean predicted measurement
  VectorXd z_pred = Zsig * weights_;

  // calculate measurement covariance matrix S and cross correlation matrix Tc
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_laspx_ * std_laspx_, 0,
       0, std_laspy_ * std_laspy_;
  MatrixXd S = R;

  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  for(int i = 0; i < 2*n_aug_ + 1; i++) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) <-M_PI) x_diff(3) += 2.*M_PI;

    S += weights_(i) * z_diff * z_diff.transpose();
    Tc+= weights_(i) * x_diff * z_diff.transpose();
  }

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  // new measurement
  VectorXd z = meas_package.raw_measurements_;

  // calculate Kalman gain K
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z - z_pred;

  // update state mean and covariance matrix
  x_ += K * z_diff;
  P_ += -K * S * K.transpose();

  // calculate lidar NIS
  double NIS_laser = z_diff.transpose() * S.inverse() * z_diff;
  cout << "NIS Lidar = " << NIS_laser / 5.991 << endl;
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


  /*****************************************************************************
   *  Predict Radar Measurement
   ****************************************************************************/

  // measurement dimension - radar
  int n_z = 3;

  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // transform sigma point into measurement space
  for(int i = 0; i < 2*n_aug_ + 1; i++) {
      double px = Xsig_pred_(0,i);
      double py = Xsig_pred_(1,i);
      double v = Xsig_pred_(2,i);
      double yaw = Xsig_pred_(3,i);

      Zsig(0,i) = sqrt(px*px + py*py);
      Zsig(1,i) = atan2(py, px);
      Zsig(2,i) = (px*cos(yaw)*v + py*sin(yaw)*v) / sqrt(px*px + py*py);
  }

  // calculate mean predicted measurement
  VectorXd z_pred = Zsig * weights_;

  // calculate measurement covariance matrix S and cross correlation matrix Tc
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_ * std_radr_, 0, 0,
       0, std_radphi_ * std_radphi_, 0,
       0, 0, std_radrd_ * std_radrd_;
  MatrixXd S = R;

  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  for(int i = 0; i < 2*n_aug_ + 1; i++) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1) <-M_PI) z_diff(1) += 2.*M_PI;
    while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3) <-M_PI) x_diff(3) += 2.*M_PI;

    S += weights_(i) * z_diff * z_diff.transpose();
    Tc+= weights_(i) * x_diff * z_diff.transpose();
  }

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  // new measurement
  VectorXd z = meas_package.raw_measurements_;

  // calculate Kalman gain
  MatrixXd K = Tc * S.inverse();

  // residual
  VectorXd z_diff = z - z_pred;

  // angle normalization
  while (z_diff(1) > M_PI) z_diff(1) -= 2.*M_PI;
  while (z_diff(1) <-M_PI) z_diff(1) += 2.*M_PI;

  // update state mean and covariance matrix
  x_ += K * z_diff;
  P_ += -K * S * K.transpose();

  // calculate radar NIS
  double NIS_radar = z_diff.transpose() * S.inverse() * z_diff;
  cout << "NIS Radar = " << NIS_radar / 7.815 << endl;
}
