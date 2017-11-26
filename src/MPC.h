#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include <string>

using namespace std;

class MPC {
 public:
  MPC();

  virtual ~MPC();
	
	double prev_a = 0;

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
	void set_mpc_param(std::string json_file);
};

#endif /* MPC_H */
