#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	// ... your code here
    if (!estimations.size() || (ground_truth.size() != estimations.size())) {
        cout << "you are using an empty matrix" << endl;
        return rmse;
    }
    
    vector<VectorXd> rmse_temp(4);
    
	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i) {
        // ... your code here
		rmse_temp[i] = estimations[i] - ground_truth[i];
		rmse_temp[i] = rmse_temp[i].array().pow(2);
	}

	//calculate the mean
	// ... your code here
	cout << rmse_temp[0] << endl;
	for (int i=0; i < rmse_temp.size(); ++i) {
	    rmse(i) = rmse_temp[i].sum()/rmse_temp.size();
	}
	
	//calculate the squared root
	// ... your code here
	rmse = rmse.array().sqrt();
	
	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//check division by zero
	for(int i; i < x_state.size(); ++i) {
        if(!x_state[i]) {
            cout << "Error - Division by zero" << endl;
            return Hj;
        }
    }
    
	//compute the Jacobian matrix
    float H00 = px/sqrt(pow(px, 2) + pow(py, 2));
    float H01 = py/sqrt(pow(px, 2) + pow(py, 2));
    float H02 = 0;
    float H03 = 0;
    float H10 = -py/(pow(px, 2) + pow(py, 2));
    float H11 = px/(pow(px, 2) + pow(py, 2));
    float H12 = 0;
    float H13 = 0;
    float H20 = py*(vx*py - vy*px)/pow(pow(px, 2) + pow(py, 2), 1.5);
    float H21 = px*(vy*px - vx*py)/pow(pow(px, 2) + pow(py, 2), 1.5);
    float H22 = H00;
    float H23 = H01;
    
    Hj << H00, H01, H02, H03,
          H10, H11, H12, H13,
          H20, H21, H22, H23;
    
	return Hj;
}
