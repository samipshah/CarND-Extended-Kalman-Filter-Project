#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
 	VectorXd rmse(4);
	rmse << 0,0,0,0;

	if(estimations.size() == 0) {
	    return rmse;
	}
	
	if(estimations.size() != ground_truth.size()) {
	    return rmse;
	}

	// accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){
	    VectorXd residual = estimations[i] - ground_truth[i];
	    residual = residual.array()*residual.array();
	    rmse += residual;
	}

	//calculate the mean
	rmse = rmse /  estimations.size();

	//calculate the squared root
  rmse = rmse.array().sqrt();
	return rmse; 
}

VectorXd Tools::cartesian_to_polar(const VectorXd &z) {
  VectorXd ret = VectorXd(3);
  float px = z[0];
  float py = z[1];
  float vx = z[2];
  float vy = z[3];

  float rho = sqrt(px*px + py*py);
  float phi = atan2(py, px);
  float rho_dot = (px*vx + py*vy)/rho;
  ret << rho, phi, rho_dot;
  return ret;
}

VectorXd Tools::polar_to_cartesian(const VectorXd &z) {
  VectorXd ret = VectorXd(4);

  float rho = z[0];
  float phi = z[1];

  float px = rho*sin(phi);
  float py = (-1)*(rho*cos(phi));
  ret << px, py, 0, 0;
  return ret;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//pre-compute a set of terms to avoid repeated calculation
	float c1 = px*px+py*py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);

	//check division by zero
	if(fabs(c1) < 0.0001){
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		Hj << 0, 0, 0, 0,
				0, 0, 0, 0,
				0, 0, 0, 0; 
		return Hj;
	}

	//compute the Jacobian matrix
	Hj << (px/c2), (py/c2), 0, 0,
		  -(py/c1), (px/c1), 0, 0,
		  py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

	return Hj;
}
