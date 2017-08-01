#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        /*MatrixXd &H_in, MatrixXd &R_in,*/ MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  // H_ = H_in;
  // R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  // process noise is modeled as N(0, Q) 0 mean Q variance hence mean ignored 
  // while calculating predicted location x_. 
  // Q_ covariance added while calculating predicted covariance.
  x_ = F_*x_; 
  P_ = F_*P_*F_.transpose() + Q_; 
}

void KalmanFilter::common_update_(const VectorXd &y, const MatrixXd &H, const MatrixXd &R) {
  MatrixXd Ht = H.transpose();

  // measurement covariance matrix is added to this update 
  MatrixXd S = H*P_*Ht + R;
  
  //
  MatrixXd K = P_*Ht*S.inverse();

  // update the model as a function of predicted location, error between measured and predicted and model matrix K
  // this determines how much error should correct the model 
  x_ = x_ + K*y;

  // update covariance
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K*H)*P_;
}

void KalmanFilter::Update(const VectorXd &z, const MatrixXd &H, const MatrixXd &R) {
  // error between measured and calculated
  VectorXd y = z - H*x_;
  common_update_(y, H, R);
}


void KalmanFilter::UpdateEKF(const VectorXd &z, const MatrixXd &H, const MatrixXd &R) {
  // error between measured and calculated 3x1
  VectorXd y = z - tools_.cartesian_to_polar(x_);
  if(y[1] > M_PI) {
    y[1] -= 2*M_PI;
  }

  if (y[1] < -M_PI) {
    y[1] += 2*M_PI;
  }
  common_update_(y, H, R);
}
