#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC
{
 public:
  MPC();

  virtual ~MPC();

  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};

extern double cte_w;
extern double psie_w;
extern double delta_w;

extern const double Lf;

extern double latency;

#endif /* MPC_H */
