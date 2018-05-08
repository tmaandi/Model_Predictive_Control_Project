#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

/* Set the timestep length and duration */
size_t N = 10;
double dt = 0.05;

/* This is the length from front to CoG of the car in simulator */
const double Lf = 2.67;

/* Reference Speed */
double ref_v = 40;

double cte_w;
double psie_w;
double delta_w;
double latency;

unsigned int latency_steps = (unsigned int)(latency/dt);

/*
 * The solver takes all the state variables and actuator
 * variables in a singular vector. Thus, we should to establish
 * when one variable starts and another ends.
*/
size_t x_start = 0;
size_t y_start = x_start + N;
size_t psi_start = y_start + N;
size_t v_start = psi_start + N;
size_t cte_start = v_start + N;
size_t epsi_start = cte_start + N;
size_t delta_start = epsi_start + N;
size_t a_start = delta_start + N - 1;

class FG_eval
{
 public:
  /* Fitted polynomial coefficients */
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars)
  {
    /*
     * `fg` a vector of the cost constraints, `vars` is a
     * vector of variable values (state & actuators)
     * The cost is stored is the first element of `fg`.
     * Any additions to the cost should be added to `fg[0]`.
     */
     fg[0] = 0;

    /*
     * Update Cost
     */
    for (unsigned int t = 0; t < N; ++t)
    {
      /*
       * Adding weighted error based on reference lane position and
       * vehicle orientation to the cost
       */
      fg[0] += cte_w * CppAD::pow(vars[cte_start + t], 2);
      fg[0] += psie_w * CppAD::pow(vars[epsi_start + t], 2);
      /*
       * Adding weighted error based on reference velocity to the cost
       */
      fg[0] += 0.2 * CppAD::pow(vars[v_start + t] - ref_v, 2);

      if (t < (N - 1))
      {
        /*
         * Adding to the cost, weighted magnitude of control input per step
         * to penalize large control inputs
         */
        fg[0] += delta_w * CppAD::pow(vars[delta_start + t], 2);
        fg[0] += 5 * CppAD::pow(vars[a_start + t], 2);

        if (t < (N - 2))
        {
          /*
           * Adding to the cost, weighted change in the magnitude of control
           * inputs between two consecutive steps, penalizing the
           * high rate of change of control inputs
           */
          fg[0] += 100 * CppAD::pow((vars[delta_start + t + 1] - vars[delta_start + t]), 2);
          fg[0] += 100 * CppAD::pow((vars[a_start + t + 1] - vars[a_start + t]), 2);
        }
      }
    }

    /*
     * Setup Constraints
     */
    /* Initial constraints
     * We add 1 to each of the starting indices due to cost being located at
     * index 0 of `fg`.
     * This bumps up the position of all the other values.
     */
    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];

    for (unsigned int t = 1; t < N; t++)
    {
      AD<double> x1 = vars[x_start + t];
      AD<double> x0 = vars[x_start + t - 1];

      AD<double> y1 = vars[y_start + t];
      AD<double> y0 = vars[y_start + t - 1];

      AD<double> psi1 = vars[psi_start + t];
      AD<double> psi0 = vars[psi_start + t - 1];

      AD<double> v1 = vars[v_start + t];
      AD<double> v0 = vars[v_start + t - 1];

      AD<double> cte1 = vars[cte_start + t];
      AD<double> cte0 = vars[cte_start + t - 1];

      AD<double> epsi1 = vars[epsi_start + t];
      AD<double> epsi0 = vars[epsi_start + t - 1];

      AD<double> delta = vars[delta_start + t - 1];
      AD<double> a = vars[a_start + t - 1];

      /* Considering the actuation delayed due to
       * latency */
      if (latency_steps > 0)
      {
        if (t > latency_steps)
        {
          delta = vars[delta_start + t - 1 - latency_steps];
          a = vars[a_start + t - 1 - latency_steps];
        }
        else
        {
          delta = 0;
          a = 0;
        }
      }

      /* The idea here is to constraint this value to be 0. */
      fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + psi_start + t] = psi1 - (psi0 + v0 * delta * dt /Lf);
      fg[1 + v_start + t] = v1 - (v0 + a * dt);
      fg[1 + cte_start + t] = cte1 - (cte0 + v0 * CppAD::sin(epsi0) * dt);
      fg[1 + epsi_start + t] = epsi1 - (epsi0 + v0 * delta * dt /Lf);
    }
  }
};

/*
 * MPC class definition implementation.
 */
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs)
{
  bool ok = true;

  typedef CPPAD_TESTVECTOR(double) Dvector;

  /*
   * Number of model variables (states and inputs)
   * Here, N timesteps == N - 1 actuations
   */
  size_t n_vars = N * 6 + (N - 1) * 2;
  /* Number of constraints */
  size_t n_constraints = N * 6;

  /*
   * Initial value of the independent variables.
   * SHOULD BE 0 besides initial state
   */
  Dvector vars(n_vars);
  for (unsigned int i = 0; i < n_vars; i++)
  {
    vars[i] = 0;
  }

  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double v = state[3];
  double cte = state[4];
  double epsi = state[5];

  /* Set the initial variable values */
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  /*
   * Set all non-actuators upper and lowerlimits
   * to the max negative and positive values.
   */
  for (unsigned int i = 0; i < delta_start; i++)
  {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  /*
   * The upper and lower limits of delta are set to -25 and 25
   * degrees (values in radians).
   */
  for (unsigned int i = delta_start; i < a_start; i++)
  {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }
  /*
   * Acceleration/decceleration upper and lower limits
   */
  for (unsigned int i = a_start; i < n_vars; i++)
  {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }
  /*
   * Lower and upper limits for the constraints
   * Should be 0 besides initial state.
   */
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (unsigned int i = 0; i < n_constraints; i++)
  {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;

  /* object that computes objective and constraints */
  FG_eval fg_eval(coeffs);

  /* options for IPOPT solver */
  std::string options;
  /* Uncomment this if you'd like more print information */
  options += "Integer print_level  0\n";
  /*
   * NOTE: Setting sparse to true allows the solver to take advantage
   * of sparse routines, this makes the computation MUCH FASTER. If you
   * can uncomment 1 of these and see if it makes a difference or not but
   * if you uncomment both the computation time should go up in orders of
   * magnitude.
   */
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  /*
   * NOTE: Currently the solver has a maximum time limit of 5 seconds.
   * Change this as you see fit.
   */
  options += "Numeric max_cpu_time          5.0\n";

  /* place to return solution */
  CppAD::ipopt::solve_result<Dvector> solution;

  /* solve the optimization problem */
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  /* Check some of the solution values */
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  /* Cost */
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  /* Returning actuator commands and predicted trajectory */
  return {solution.x[delta_start], solution.x[a_start],
          solution.x[x_start + 1],   solution.x[x_start + 2],
          solution.x[x_start + 3],   solution.x[x_start + 4],
          solution.x[x_start + 5],   solution.x[x_start + 6],
          solution.x[y_start + 1],   solution.x[y_start + 2],
          solution.x[y_start + 3],   solution.x[y_start + 4],
          solution.x[y_start + 5],   solution.x[y_start + 6]};
}
