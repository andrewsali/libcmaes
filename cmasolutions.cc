
#include "cmasolutions.h"
#include "opti_err.h"
#include <limits>
#include <iostream>

namespace libcmaes
{
  
  CMASolutions::CMASolutions(Parameters &p)
    :_hsig(1),_max_eigenv(0.0),_min_eigenv(0.0),_niter(0),_kcand(1),_eigeniter(0),_updated_eigen(true),_run_status(0),_elapsed_time(0)
  {
    try
      {
	_cov = dMat::Identity(p._dim,p._dim);
      }
    catch (std::bad_alloc &e)
      {
	_run_status = OPTI_ERR_OUTOFMEMORY;
	return;
      }
    if (p._x0min == p._x0max)
      {
	if (p._x0min == dVec::Constant(p._dim,std::numeric_limits<double>::min()))
	  _xmean = dVec::Random(p._dim) * 4.0; // initial mean randomly sampled from -4,4 in all dimensions.
	else _xmean = p._x0min;
      }
    else
      {
	_xmean = 0.5*(dVec::Random(p._dim) + dVec::Constant(p._dim,1.0)); // scale to [0,1].
	_xmean = (p._x0max - p._x0min)*_xmean + p._x0min; // scale to bounds.
      }
    if (static_cast<CMAParameters&>(p)._sigma_init > 0.0)
      _sigma = static_cast<CMAParameters&>(p)._sigma_init;
    else static_cast<CMAParameters&>(p)._sigma_init = _sigma = 1.0/static_cast<double>(p._dim); // XXX: sqrt(trace(cov)/dim)
    
    _psigma = dVec::Zero(p._dim);
    _pc = dVec::Zero(p._dim);
    _candidates.resize(p._lambda);
    _kcand = static_cast<int>(1.0+floor(0.1+p._lambda/4.0));
  }

  CMASolutions::~CMASolutions()
  {
  }

  void CMASolutions::update_best_candidates()
  {
    _best_candidates_hist.push_back(_candidates.at(0)); // supposed candidates is sorted.
    _k_best_candidates_hist.push_back(_candidates.at(_kcand));

    _bfvalues.push_back(_candidates.at(0)._fvalue);
    if (_bfvalues.size() > 20)
      _bfvalues.erase(_bfvalues.begin());

    // get median of candidate's scores, used in termination criteria (stagnation).
    double median = 0.0;
    size_t csize = _candidates.size();
    if (csize % 2 == 0)
      median = (_candidates[csize/2-1]._fvalue + _candidates[csize/2]._fvalue)/2.0;
    else median = _candidates[csize/2]._fvalue;
    _median_fvalues.push_back(median);
    if (_median_fvalues.size() > static_cast<size_t>(ceil(0.2*_niter+120+30*_xmean.size()/static_cast<double>(_candidates.size()))))
      _median_fvalues.erase(_median_fvalues.begin());
    
    //debug
    /*std::cerr << "ordered candidates:\n";
    for (size_t i=0;i<_candidates.size();i++)
      {
	std::cerr << _candidates.at(i)._fvalue << " / " << _candidates.at(i)._x.transpose() << std::endl;
	}*/
    //debug
  }

  void CMASolutions::update_eigenv(const dVec &eigenvalues,
				   const dMat &eigenvectors)
  {
    _max_eigenv = eigenvalues.maxCoeff();
    _min_eigenv = eigenvalues.minCoeff();
    _leigenvalues = eigenvalues;
    _leigenvectors = eigenvectors;
  }

  std::ostream& CMASolutions::print(std::ostream &out,
				    const int &verb_level) const
  {
    out << "best solution => f-value=" << best_candidate()._fvalue << " / sigma=" << _sigma << " / iter=" << _niter << " / elaps=" << _elapsed_time << "ms";
    if (verb_level)
      {
	out << "\ncovdiag=" << _cov.diagonal().transpose() << std::endl;
	out << "psigma=" << _psigma.transpose() << std::endl;
	out << "pc=" << _pc.transpose() << std::endl;
      }
    return out;
  }

  std::ostream& operator<<(std::ostream &out, const CMASolutions &cmas)
  {
    cmas.print(out,0);
    return out;
  }
  
}
