/**
 * CMA-ES, Covariance Matrix Adaptation Evolution Strategy
 * Copyright (c) 2014 Inria
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

#ifndef CANDIDATE_H
#define CANDIDATE_H

#include "eo_matrix.h"
#include "cmaparameters.h"

namespace libcmaes
{
  /**
   * \brief candidate solution point, in function parameter space.
   */
  class Candidate
  {
  public:
    /**
     * \brief empty constructor.
     */
  Candidate():
    _fvalue(0.0) {}
    
    /**
     * \brief constructor.
     * @param fvalue function value
     * @param x function parameter vector
     */
  Candidate(const double &fvalue,
	    const dVec &x)
    :_fvalue(fvalue),_x(x)
    {}

  ~Candidate() {}

  double get_fvalue() const { return _fvalue; }

  dVec get_x_dvec() const { return _x; }

  const double* get_x() const { return _x.data(); }

  template<class TGenoPheno>
    dVec get_x_pheno(const CMAParameters<TGenoPheno> &p) const
    {
      dVec gx = p._gp.pheno(_x);
      return gx;
    }
  
  double _fvalue; /**< function value. */
  dVec _x; /**< function parameter vector. */
  };

}

#endif
