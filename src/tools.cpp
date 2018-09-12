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
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  if (estimations.size()==0){
	cout<<"Error - The estimation vector size is zero\n";
	return rmse;
  }
  else if (estimations.size() != ground_truth.size()){
	cout<<"Error - The estimation vector size and ground_truth vector size are not equal\n";
	return rmse;
  }

  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
	// ... your code here

        VectorXd sum = estimations[i] - ground_truth[i];
        sum = sum.array() * sum.array();
        rmse += sum;
  }

  //calculate the mean
  rmse = rmse/estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  return rmse;
}
