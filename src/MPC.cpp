#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
#include <string>
#include <fstream>
#include "json.hpp"

using json = nlohmann::json;
using CppAD::AD;

// TUNABLE PARAMETERS
// These global variables are set by calling MPC:set_mpc_param() with a json config file as parameter

// Reference values for the MPC
double ref_cte;
double ref_eps;
double ref_v;

// Set the timestep length and duration
size_t N;
double dt;

// Cost factors
double cost_cte;
double cost_eps;
double cost_v;
double cost_delta_cur;
double cost_delta_diff;
double cost_a_cur;
double cost_a_diff;

// COMPUTED PARAMETERS
// These are set by calling MPC_set_mpc_param() and depend on the json config file passed as parameter

// Calculate start index in the state vector
size_t x_idx;
size_t y_idx;
size_t psi_idx;
size_t v_idx;
size_t cte_idx;
size_t eps_idx;
size_t delta_idx;
size_t a_idx;

// FIXED PARAMETERS

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

class FG_eval {
	 public:
		// Fitted polynomial coefficients
		Eigen::VectorXd coeffs;
		FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

		typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
		
		void operator()(ADvector& fg, const ADvector& vars) {
				// `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
				// NOTE: You'll probably go back and forth between this function and
				// the Solver function below.
				
				// ----------------------------
				// Cost for the reference state
				// ----------------------------
				
				// Initialise cost to zero
				fg[0] = 0;
				
				// Cost for distance from reference state -> brings current state to reference state
				for (size_t i=0; i < N; i++) {
						fg[0] += cost_cte * pow(vars[cte_idx + i] - ref_cte, 2);
						fg[0] += cost_eps * pow(vars[eps_idx + i] - ref_eps, 2);
						fg[0] += cost_v * pow(vars[v_idx + i] - ref_v, 2);		
				}
				
				// Cost for usage of actuators -> penalizes erratic actuation
				for (size_t i=0; i < N-1; i++) {
						fg[0] += cost_delta_cur * pow(vars[delta_idx + i], 2);
						fg[0] += cost_a_cur * pow(vars[a_idx + i], 2);
				}
				
				// Cost for value difference in actuation between time steps -> smoothes actuation
				for (size_t i=0; i < N-2; i++) {
						fg[0] += cost_delta_diff * pow(vars[delta_idx + i + 1] - vars[delta_idx + i], 2);
						fg[0] += cost_a_diff * pow(vars[a_idx + i + 1] - vars[a_idx + i], 2);
				}

				// -----------
				// Constraints
				// -----------
				
				// Set the inital state, use idx + 1 since 0 holds the cost function
				fg[x_idx + 1] = vars[x_idx];
				fg[y_idx + 1] = vars[y_idx];
				fg[psi_idx + 1] = vars[psi_idx];
				fg[v_idx + 1] = vars[v_idx];
				fg[cte_idx + 1] = vars[cte_idx];
				fg[eps_idx + 1] = vars[eps_idx];
				
				// Calculate the predicted states for N-1 time steps
				for (size_t i=0; i < N-1; i++) {
						// Temporary variables for the state at time t
						AD<double> x0 = vars[x_idx + i];
						AD<double> y0 = vars[y_idx + i];
						AD<double> psi0 = vars[psi_idx + i];
						AD<double> v0 = vars[v_idx + i];
						AD<double> cte0 = vars[cte_idx + i];
						AD<double> eps0 = vars[eps_idx + i];
	
						AD<double> delta0 = vars[delta_idx + i];
						AD<double> a0 = vars[a_idx + i];
								 
						AD<double> f0 = coeffs[0] + coeffs[1]*x0 + coeffs[2]*pow(x0,2) + coeffs[3]*pow(x0,3);
						AD<double> psides0 = CppAD::atan( coeffs[1] + 2* coeffs[2]*x0 + 3*coeffs[3]*pow(x0,2) );
	
						// Temporary variables for the state at time t+1
						AD<double> x1 = vars[x_idx + i + 1];
						AD<double> y1 = vars[y_idx + i + 1];
						AD<double> psi1 = vars[psi_idx + i + 1];
						AD<double> v1 = vars[v_idx + i + 1];
						AD<double> cte1 = vars[cte_idx + i + 1];
						AD<double> eps1 = vars[eps_idx + i + 1];
						
						// Calculate difference between actual and predicted state and store in fg
						fg[2 + x_idx + i] = x1 - (x0 + v0*CppAD::cos(psi0)*dt);
						fg[2 + y_idx + i] = y1 - (y0 + v0*CppAD::sin(psi0)*dt);
						fg[2 + psi_idx + i] = psi1 - (psi0 - v0*(delta0/Lf)*dt);
						fg[2 + v_idx + i] = v1 - (v0 + a0*dt);
						fg[2 + cte_idx + i] =	cte1 - ((f0 - y0) + (v0*CppAD::sin(eps0)*dt));
						fg[2 + eps_idx + i] =	eps1 - ((psi0 - psides0) - v0*(delta0/Lf)*dt);
				}
		}
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

void MPC::set_mpc_param(std::string json_file) {
		// Load parameters from json file
		std::ostringstream param_buf; 
		std::ifstream param_file(json_file); 
		param_buf << param_file.rdbuf(); 
		auto param = json::parse(param_buf.str());
		
		// Set reference values
		ref_cte = param["reference"]["cte"].get<double>();
		ref_eps = param["reference"]["eps"].get<double>();
		ref_v = param["reference"]["v"].get<double>();
		
		// Set time horizon
		N = (size_t) param["horizon"]["N"].get<int>();
		dt = param["horizon"]["dt"].get<double>();
		
		// Set cost factors
		cost_cte = param["cost"]["cte"].get<double>();
		cost_eps = param["cost"]["eps"].get<double>();
		cost_v = param["cost"]["v"].get<double>();
		cost_delta_cur = param["cost"]["delta_cur"].get<double>();
		cost_delta_diff = param["cost"]["delta_diff"].get<double>();
		cost_a_cur = param["cost"]["a_cur"].get<double>();
		cost_a_diff = param["cost"]["a_diff"].get<double>();
		
		// Set the dependend size variables
		x_idx = 0;
		y_idx = x_idx + N;
		psi_idx = y_idx + N;
		v_idx = psi_idx + N;
		cte_idx = v_idx + N;
		eps_idx = cte_idx + N;
		delta_idx = eps_idx + N;
		a_idx = delta_idx + N - 1;
}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
		bool ok = true;
		typedef CPPAD_TESTVECTOR(double) Dvector;

		// Model Variables
		// State: [x, y, psi, v, cte, eps]
		// Actuators: [delta, a]
		size_t n_vars = 6*N + 2*(N-1);
		
		// Constraints
		// [x, y, psi, v, cte, eps]
		size_t n_constraints = 6*N;

		// -------------------------------------------
		// Initial value of the independent variables.
		// -------------------------------------------
		// SHOULD BE 0 besides initial state.
		Dvector vars(n_vars);
		for (size_t i=0; i < n_vars; i++) {
			vars[i] = 0;
		}

		vars[x_idx] = state[0];
		vars[y_idx] = state[1];
		vars[psi_idx] = state[2];
		vars[v_idx] = state[3];
		vars[cte_idx] = state[4];
		vars[eps_idx] = state[5];
		
		// -----------------------------------------------
		// Lower and upper bounds of independent variables
		// -----------------------------------------------
		Dvector vars_lowerbound(n_vars);
		Dvector vars_upperbound(n_vars);

		// Bounds for [x, y, psi, v, cte, eps] not required, set to very high values
		for (size_t i=0; i < delta_idx; i++) {
				vars_upperbound[i] = 1.0e10;
				vars_lowerbound[i] = -1.0e10;
		}
  
		// Steering angle to max values of simulator
		for (size_t i=delta_idx; i < a_idx; i++)
		{
				vars_upperbound[i] = M_PI/8;
				vars_lowerbound[i] = -M_PI/8;
		}
  
		// Acceleration to -1 and 1
		for (size_t i=a_idx; i < n_vars; i++)
		{
				vars_upperbound[i] = 1.0;
				vars_lowerbound[i] = -1.0;
		}
  
		// ------------------------------------------
		// Lower and upper limits for the constraints
		// ------------------------------------------
		// Should be 0 besides initial state.
		Dvector constraints_lowerbound(n_constraints);
		Dvector constraints_upperbound(n_constraints);
		
		for (size_t i=0; i < n_constraints; i++) {
				constraints_lowerbound[i] = 0;
				constraints_upperbound[i] = 0;
		}
		// Set init state lower and upper limits
		constraints_lowerbound[x_idx] = state[0];
		constraints_lowerbound[y_idx] = state[1];
		constraints_lowerbound[psi_idx] = state[2];
		constraints_lowerbound[v_idx] = state[3];
		constraints_lowerbound[cte_idx] = state[4];
		constraints_lowerbound[eps_idx] = state[5];
		
		constraints_upperbound[x_idx] = state[0];
		constraints_upperbound[y_idx] = state[1];
		constraints_upperbound[psi_idx] = state[2];
		constraints_upperbound[v_idx] = state[3];
		constraints_upperbound[cte_idx] = state[4];
		constraints_upperbound[eps_idx] = state[5];

		// ----------------------------------------------
		// object that computes objective and constraints
		// ----------------------------------------------
		FG_eval fg_eval(coeffs);

		//
		// NOTE: You don't have to worry about these options
		//
		// options for IPOPT solver
		std::string options;
		// Uncomment this if you'd like more print information
		options += "Integer print_level  0\n";
		// NOTE: Setting sparse to true allows the solver to take advantage
		// of sparse routines, this makes the computation MUCH FASTER. If you
		// can uncomment 1 of these and see if it makes a difference or not but
		// if you uncomment both the computation time should go up in orders of
		// magnitude.
		options += "Sparse  true        forward\n";
		options += "Sparse  true        reverse\n";
		// NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
		// Change this as you see fit.
		options += "Numeric max_cpu_time          0.5\n";

		// place to return solution
		CppAD::ipopt::solve_result<Dvector> solution;

		// solve the problem
		CppAD::ipopt::solve<Dvector, FG_eval>(
				options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
				constraints_upperbound, fg_eval, solution);

		// Check some of the solution values
		ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

		// Cost
		auto cost = solution.obj_value;
		std::cout << "Cost " << cost << std::endl;

		// --------------------------------
		// Return the first actuator values 
		// --------------------------------
		vector<double> result;
		result.push_back(solution.x[delta_idx]);
		result.push_back(solution.x[a_idx]);
		
		for (size_t i=0; i < N-1; i++) {
				result.push_back(solution.x[x_idx + i + 1]);
				result.push_back(solution.x[y_idx + i + 1]);
		}
		return result;
}
