/**
 * CMA-ES, Covariance Matrix Adaptation Evolution Strategy
 * Copyright (c) 2014 INRIA
 * Author: Emmanuel Benazera <emmanuel.benazera@lri.fr>
 *
 * This file is part of libcmaes.
 *
 * libcmaes is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * libcmaes is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with libcmaes.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ERRSTATS_H
#define ERRSTATS_H

#include "cmaes.h"

namespace libcmaes
{
  /**
   * \brief profile likelihood as a set of points and values. 
   */
  class pli
  {
  public:
  pli(const int &k, const int &samplesize, const int &dim,
      const dVec &xm, const double &fvalue)
    :_k(k),_samplesize(samplesize),_fvaluem(dVec::Zero(2*samplesize+1)),_xm(dMat::Zero(2*samplesize+1,dim)),_min(0.0),_max(0.0)
      {
	_fvaluem[samplesize] = fvalue;
	_xm.row(samplesize) = xm.transpose();
      }
    ~pli() {};

    void setMinMax()
    {
      _min = _xm(0,_k);
      _max = _xm(2*_samplesize,_k);
      if (_min > _max)
	std::swap(_min,_max);
    }

    std::pair<double,double> getMinMax(const double &fvalue)
    {
      dMat::Index mindex[2];
      (_fvaluem.head(_samplesize) - dVec::Constant(_samplesize,fvalue)).cwiseAbs().minCoeff(&mindex[0]);
      (_fvaluem.tail(_samplesize) - dVec::Constant(_samplesize,fvalue)).cwiseAbs().minCoeff(&mindex[1]);
      double min = _xm(mindex[0],_k);
      double max = _xm(_samplesize + 1 + mindex[1],_k);
      if (min > max)
	std::swap(min,max);
      return std::pair<double,double>(min,max);
    }

    int _k;
    int _samplesize;
    dVec _fvaluem;
    dMat _xm;
    double _min;
    double _max;
  };

  template <class TGenoPheno=GenoPheno<NoBoundStrategy>>
  class errstats
    {
    public:
    static pli profile_likelihood(FitFunc &func,
				  CMAParameters<TGenoPheno> &parameters,
				  const CMASolutions &cmasol,
				  const int &k,
				  const bool &curve=false,
				  const int &samplesize=1000,
				  const double &fup=0.1,
				  const double &delta=0.1);

    static void profile_likelihood_search(FitFunc &func,
					  CMAParameters<TGenoPheno> &parameters,
					  pli &le,
					  const CMASolutions &cmasol,
					  const int &k,
					  const bool &neg,
					  const int &samplesize,
					  const double &fup,
					  const double &delta,
					  const bool &curve);
    
    static void take_linear_step(FitFunc &func,
				 const int &k,
				 const double &minfvalue,
				 const double &fup,
				 const bool &curve,
				 dVec &x,
				 double &dxk);

    static CMASolutions optimize_pk(FitFunc &func,
				    CMAParameters<TGenoPheno> &parameters,
				    const CMASolutions &cmasol,
				    const int &k,
				    const double &vk);
    
    // DEPRECATED
    static CMASolutions optimize_reduced_pk(FitFunc &func,
					    CMAParameters<TGenoPheno> &parameters,
					    const CMASolutions &cmasol,
					    const int &k,
					    const double &vk);
    };    
  
}

#endif
