
#include "bipopcmastrategy.h"
#include "opti_err.h"
#include <random>
#include <glog/logging.h>
#include <time.h>

namespace libcmaes
{
  BIPOPCMAStrategy::BIPOPCMAStrategy(FitFunc &func,
				     CMAParameters &parameters)
    :IPOPCMAStrategy(func,parameters),_lambda_def(parameters._lambda),_lambda_l(parameters._lambda)
  {
    std::random_device rd;
    _gen = std::mt19937(rd());
    _gen.seed(static_cast<uint64_t>(time(NULL)));
    _unif = std::uniform_real_distribution<>(0,1);
    _lambda_def = 4.0+ceil(3.0+log(_parameters._dim));
    _parameters._lambda = _lambda_def;
    _solutions = CMASolutions(_parameters);
  }

  BIPOPCMAStrategy::~BIPOPCMAStrategy()
  {
  }

  void BIPOPCMAStrategy::tell()
  {
    CMAStrategy::tell();
  }

  int BIPOPCMAStrategy::optimize()
  {
    std::array<int,2> budgets = {0,0}; // 0: r1, 1: r2
    CMASolutions best_run;
    for (int r=0;r<_parameters._nrestarts;r++)
      {
	while(budgets[0]>budgets[1])
	  {
	    r2();
	    reset_search_state();
	    CMAStrategy::optimize();
	    budgets[1] += _solutions._niter * _parameters._lambda;
	    capture_best_solution(best_run);
	  }
	if (r > 0) // use lambda_def on first call.
	  {
	    r1();
	    reset_search_state();
	  }
	CMAStrategy::optimize();
	budgets[0] += _solutions._niter * _parameters._lambda;
	capture_best_solution(best_run);
      }
    _solutions = best_run;
    if (_solutions._run_status >= 0)
      return OPTI_SUCCESS;
    else return OPTI_ERR_TERMINATION; // exact termination code is in _solutions._run_status.
  }

  void BIPOPCMAStrategy::r1()
  {
    _parameters._lambda = _lambda_l;
    lambda_inc(); // from IPOP.
    _lambda_l = _parameters._lambda;
  }

  void BIPOPCMAStrategy::r2()
  {
    double u = _unif(_gen);
    double ltmp = pow(0.5*(_lambda_l/_lambda_def),u);
    double nlambda = ceil(_lambda_def * ltmp);
    LOG_IF(INFO,!_parameters._quiet) << "Restart => lambda_s=" << nlambda << " / lambda_old=" << _parameters._lambda << " / lambda_l=" << _lambda_l << " / lambda_def=" << _lambda_def << std::endl;
    _parameters._lambda = nlambda;
  }
  
}
