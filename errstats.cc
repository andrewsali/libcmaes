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

#include "errstats.h"
#include <llogging.h>
#include <iostream>

namespace libcmaes
{
  
  template <class TGenoPheno>
  pli errstats<TGenoPheno>::profile_likelihood(FitFunc &func,
					       const CMAParameters<TGenoPheno> &parameters,
					       CMASolutions &cmasol,
					       const int &k,
					       const bool &curve,
					       const int &samplesize,
					       const double &fup,
					       const double &delta)
  {
    dVec x = cmasol.best_candidate()._x;
    double minfvalue = cmasol.best_candidate()._fvalue;
    dVec phenox = parameters._gp.pheno(x);
    
    //debug
    //std::cout << "xk=" << x[k] << " / minfvalue=" << minfvalue << std::endl;
    //std::cout << "phenox=" << phenox << std::endl;
    //debug

    pli le(k,samplesize,parameters._dim,parameters._gp.pheno(x),minfvalue,fup,delta);

    errstats<TGenoPheno>::profile_likelihood_search(func,parameters,le,cmasol,k,false,samplesize,fup,delta,curve); // positive direction
    errstats<TGenoPheno>::profile_likelihood_search(func,parameters,le,cmasol,k,true,samplesize,fup,delta,curve);  // negative direction
    
    le.setErrMinMax();
    cmasol._pls.insert(std::pair<int,pli>(k,le));
    return le;
  }

  template <class TGenoPheno>
  void errstats<TGenoPheno>::profile_likelihood_search(FitFunc &func,
						       const CMAParameters<TGenoPheno> &parameters,
						       pli &le,
						       const CMASolutions &cmasol,
						       const int &k,
						       const bool &neg,
						       const int &samplesize,
						       const double &fup,
						       const double &delta,
						       const bool &curve)
  {
    int sign = neg ? -1 : 1;
    dVec x = cmasol.best_candidate()._x;
    double xk = x[k];
    double minfvalue = cmasol.best_candidate()._fvalue;
    CMASolutions citsol = cmasol;
    double dxk = sign * xk * 0.1;
    for (int i=0;i<samplesize;i++)
      {
	// get a new xk point.
	bool iterend = errstats<TGenoPheno>::take_linear_step(func,parameters,k,minfvalue,fup,curve,x,dxk);
		
	//debug
	//std::cout << "new xk point: " << x.transpose() << std::endl;
	//debug
	
	// minimize.
	CMASolutions ncitsol = errstats<TGenoPheno>::optimize_pk(func,parameters,citsol,k,x[k]);
	if (ncitsol._run_status < 0)
	  {
	    LOG(WARNING) << "profile likelihood linesearch: optimization error " << ncitsol._run_status << std::endl;
	    // pad and return.
	    /*for (int j=i+1;j<samplesize;j++)
	      {
		le._fvaluem[samplesize+sign*(1+j)] = le._fvaluem[samplesize+sign*i];
		le._xm.row(samplesize+sign*(1+j)) = le._xm.row(samplesize+sign*i);
	      }
	      return;*/
	  }
	else // update current point and solution.
	  {
	    citsol = ncitsol;
	    x = citsol.best_candidate()._x;
	    minfvalue = citsol.best_candidate()._fvalue;
	  }
	
	// store points.
	dVec phenobx = parameters._gp.pheno(citsol.best_candidate()._x);
	le._fvaluem[samplesize+sign*(1+i)] = citsol.best_candidate()._fvalue;
	le._xm.row(samplesize+sign*(1+i)) = phenobx.transpose();
	le._err[samplesize+sign*(1+i)] = ncitsol._run_status;
	
	if (!curve && iterend)
	  {
	    // pad and return.
	    for (int j=i+1;j<samplesize;j++)
	      {
		le._fvaluem[samplesize+sign*(1+j)] = citsol.best_candidate()._fvalue;
		le._xm.row(samplesize+sign*(1+j)) = phenobx.transpose();
		le._err[samplesize+sign*(1+j)] = ncitsol._run_status;
	      }
	    return;
	  }
      }
  }
						
  template <class TGenoPheno>
  bool errstats<TGenoPheno>::take_linear_step(FitFunc &func,
					      const CMAParameters<TGenoPheno> &parameters,
					      const int &k,
					      const double &minfvalue,
					      const double &fup,
					      const bool &curve,
					      dVec &x,
					      double &dxk)
  {
    static double fdiff_relative_increase = 0.1;
    double fdelta = 0.1 * fup;
    double threshold = minfvalue + fup;
    dVec xtmp = x;
    xtmp[k] += dxk;
    dVec phenoxtmp = parameters._gp.pheno(xtmp);
    double fvalue = func(phenoxtmp.data(),xtmp.size());
    double fdiff = fvalue - minfvalue;

    
    //debug
    /*std::cout << "xtmp=" << xtmp.transpose() << std::endl;
    std::cout << "phenoxtmp=" << phenoxtmp.transpose() << std::endl;
    std::cout << "dxk=" << dxk << " / threshold=" << threshold << " / fvalue=" << fvalue << " / fdiff=" << fdiff << " / fabs=" << fabs(fvalue-fup) << " / fdelta=" << fdelta << std::endl;*/
    //debug

    if (fdiff > threshold * fdiff_relative_increase) // decrease dxk
      {
	while((curve || fabs(fvalue-fup)>fdelta)
	      && fdiff > threshold * fdiff_relative_increase
	      && phenoxtmp[k] >= parameters._gp._boundstrategy.getPhenoLBound(k))
	  {
	    //std::cerr << "fvalue=" << fvalue << " / fdelta=" << fdelta << " / fup= " << fup << " / xtmpk=" << phenoxtmp[k] << " / dxk=" << dxk << " / fdiff=" << fdiff << " / thresh=" << threshold * fdiff_relative_increase << std::endl;//" / lbound=" << parameters._gp._boundstrategy.getLBound(k) << std::endl;
	    dxk /= 2.0;
	    xtmp[k] = x[k] + dxk;
	    phenoxtmp = parameters._gp.pheno(xtmp);
	    fvalue = func(phenoxtmp.data(),xtmp.size());
	    fdiff = fvalue - minfvalue;
	  }
      }
    else // increase dxk
      {
	while ((curve || fabs(fvalue-fup)>fdelta)
	       && fdiff < threshold * fdiff_relative_increase
	       && phenoxtmp[k] <= parameters._gp._boundstrategy.getPhenoUBound(k))
	  {
	    //std::cerr << "fvalue=" << fvalue << " / xtmpk=" << phenoxtmp[k] << " / ubound=" << parameters._gp._boundstrategy.getUBound(k) << std::endl;
	    dxk *= 2.0;
	    xtmp[k] = x[k] + dxk;
	    phenoxtmp =parameters._gp.pheno(xtmp);
	    fvalue = func(phenoxtmp.data(),xtmp.size());
	    fdiff = fvalue - minfvalue;
	  }
	dxk /= 2.0;
      }
    x[k] = xtmp[k]; // set value.
    dVec phenox = parameters._gp.pheno(x);
    return (fabs(fvalue-fup) < fdelta || phenox[k] < parameters._gp._boundstrategy.getPhenoLBound(k) || phenox[k] > parameters._gp._boundstrategy.getPhenoUBound(k));
  }

  template <class TGenoPheno>
  CMASolutions errstats<TGenoPheno>::optimize_vpk(FitFunc &func,
						  const CMAParameters<TGenoPheno> &parameters,
						  const CMASolutions &cmasol,
						  const std::vector<int> &k,
						  const std::vector<double> &vk)
  {
    CMASolutions ncmasol = cmasol;
    CMAParameters<TGenoPheno> nparameters = parameters;
    nparameters._quiet = true; //TODO: option.
    ncmasol._sigma = 0.0;
    for (size_t i=0;i<k.size();i++)
      {
	nparameters.set_fixed_p(k[i],vk[i]);
	nparameters._sigma_init = ncmasol._sigma = std::max(ncmasol._sigma,fabs(cmasol.best_candidate()._x[k[i]]-vk[i])); // XXX: possibly a better sigma selection ?
      }
    return cmaes(func,nparameters,CMAStrategy<CovarianceUpdate,TGenoPheno>::_defaultPFunc,nullptr,ncmasol); //TODO: explicitely set the initial covariance.
  }

  template <class TGenoPheno>
  CMASolutions errstats<TGenoPheno>::optimize_pk(FitFunc &func,
						 const CMAParameters<TGenoPheno> &parameters,
						 const CMASolutions &cmasol,
						 const int &k,
						 const double &vk)
  {
    std::vector<int> tk = {k};
    std::vector<double> tvk = {vk};
    return optimize_vpk(func,parameters,cmasol,tk,tvk);
  }

  template <class TGenoPheno>
  CMASolutions errstats<TGenoPheno>::optimize_reduced_pk(FitFunc &func,
							 CMAParameters<TGenoPheno> &parameters,
							 const CMASolutions &cmasol,
							 const int &k,
							 const double &vk)
  {
    // set re-arranged solution as starting point of the new optimization problem.
    CMASolutions ncmasol = cmasol;
    ncmasol.reset_as_fixed(k);

    // re-arrange function
    FitFunc nfunc = [func,k,vk](const double *x, const int N)
      {
	std::vector<double> nx;
	nx.reserve(N+1);
	for (int i=0;i<N+1;i++)
	  {
	    if (i < k)
	      nx.push_back(x[i]);
	    else if (i == k)
	      nx.push_back(vk);
	    else nx.push_back(x[i-1]);
	  }
	return func(&nx.front(),N+1);
      };
    
    // optimize and return result.
    CMAParameters<TGenoPheno> nparameters = parameters;
    nparameters.reset_as_fixed(k);
    return cmaes(nfunc,nparameters);
  }

  template <class TGenoPheno>
  contour errstats<TGenoPheno>::contour_points(FitFunc & func, const int &px, const int &py, const int &npoints, const double &fup,
					       const CMAParameters<TGenoPheno> &parameters,
					       CMASolutions &cmasol)
  {
    // find first two points.
    int samplesize = 10;
    pli plx,ply;
    if (!cmasol.get_pli(px,plx))
      {
	errstats<TGenoPheno>::profile_likelihood(func,parameters,cmasol,px,false,samplesize,fup);
	cmasol.get_pli(px,plx);
      }
    
    // find second two points.
    if (!cmasol.get_pli(py,ply))
      {
	errstats<TGenoPheno>::profile_likelihood(func,parameters,cmasol,py,false,samplesize,fup);
	cmasol.get_pli(py,ply);
      }

    double valx = cmasol.best_candidate()._x[px];
    double valy = cmasol.best_candidate()._x[py];
    
    // find upper y value for x parameter.
    CMAParameters<TGenoPheno> nparameters = parameters;
    nparameters.set_x0(cmasol.best_candidate()._x);
    nparameters.set_fixed_p(px,valx+plx._errmax);
    CMASolutions exy_up = cmaes(func,nparameters);
    std::cout << "exy_up=" << exy_up.best_candidate()._x.transpose() << std::endl;
    
    // find lower y value for x parameter.
    nparameters = parameters;
    nparameters.set_x0(cmasol.best_candidate()._x);
    nparameters.set_fixed_p(px,valx+plx._errmin);
    CMASolutions exy_lo = cmaes(func,nparameters);
    std::cout << "exy_lo=" << exy_lo.best_candidate()._x.transpose() << std::endl;
    
    // find upper x value for y parameter.
    //nparameters = parameters;
    TGenoPheno gp;
    //nparameters = CMAParameters<TGenoPheno>(parameters._dim,cmasol.best_candidate()._x.data(),0.1,-1,0,gp);
    nparameters = parameters;
    nparameters.set_x0(cmasol.best_candidate()._x);
    nparameters.set_fixed_p(py,valy+ply._errmax);
    CMASolutions eyx_up = cmaes<TGenoPheno>(func,nparameters);
    std::cout << "eyx_up=" << eyx_up.best_candidate()._x.transpose() << std::endl;
    
    // find lower x value for y parameter.
    nparameters = parameters;
    nparameters.set_x0(cmasol.best_candidate()._x);
    nparameters.set_fixed_p(py,valy+ply._errmin);
    CMASolutions eyx_lo = cmaes(func,nparameters);
    std::cout << "eyx_lo=" << eyx_lo.best_candidate()._x.transpose() << std::endl;
    
    contour c;
    c.add_point(valx+plx._errmin,exy_lo.best_candidate()._x[py]);
    c.add_point(eyx_lo.best_candidate()._x[px],valy+ply._errmin); 
    c.add_point(valx+plx._errmax,exy_up.best_candidate()._x[py]);
    c.add_point(eyx_up.best_candidate()._x[px],valy+ply._errmax);

    double scalx = 1.0/(plx._errmax - plx._errmin);
    double scaly = 1.0/(ply._errmax - ply._errmin);

    std::cout << "contour:" << c << std::endl;
    
    //TODO: more than 4 points.
    for (int i=4;i<npoints;i++)
      {
	std::cout << "=> generating point #" << i << std::endl;
	
	//TODO: - check on max budget, and return if exceeded.
	
	//TODO:- get most distant points.
	std::vector<std::pair<double,double>>::iterator idist1 = c._points.end()-1;
	std::vector<std::pair<double,double>>::iterator idist2 = c._points.begin();
	double dx = idist1->first - idist2->first;
	double dy = idist1->second - idist2->second;
	double bigdis = scalx*scalx*dx*dx + scaly*scaly*dy*dy;

	for (auto ipair = c._points.begin();ipair!=c._points.end()-1;ipair++)
	  {
	    dx = ipair->first - (ipair+1)->first;
	    dy = ipair->second - (ipair+1)->second;
	    double dist = scalx*scalx*dx*dx + scaly*scaly*dy*dy;
	    if (dist > bigdis)
	      {
		bigdis = dist;
		idist1 = ipair;
		idist2 = ipair+1;
	      }
	    //std::cout << "idist10=" << idist1->first << " -- idist20=" << idist2->first << std::endl;
	  }
	
	//TODO:- select mid-range point x and direction dir along the two axis of interest.
	double a1 = 0.5;
	double a2 = 0.5;
	double sca = 1.0;
	double xmidcr = a1*idist1->first + a2*idist2->first;
	double ymidcr = a1*idist1->second + a2*idist2->second;
	double xdir = idist2->second - idist1->second;
	double ydir = idist1->first - idist2->first;
	double scalfac = sca*std::max(fabs(xdir*scalx), fabs(ydir*scaly));
	double xdircr = xdir/scalfac;
	double ydircr = ydir/scalfac;
	std::vector<double> pmid = {xmidcr, ymidcr};
	std::vector<double> pdir = {xdircr, ydircr};
	std::vector<int> par = {px,py};
	
	//TODO: printout of pmid / pdir.
	std::cout << "idist10=" << idist1->first << " / idist11=" << idist1->second << " / idist20=" << idist2->first << " / idist21=" << idist2->second << std::endl;
	std::cout << "pmid0=" << pmid[0] << " / pmid1=" << pmid[1] << " / pdir0=" << pdir[0] << " / pdir1=" << pdir[1] << std::endl;
	
	//TODO: find crossing point from x with direction dir where function is equal to min + fup.
	fcross fc = errstats<TGenoPheno>::cross(func,parameters,cmasol,fup,par,pmid,pdir,parameters._ftolerance);

	//debug
	std::cout << "fcross point=" << fc._x.transpose() << std::endl;
	//debug

	if (idist2 == c._points.begin())
	  {
	    c.add_point(fc._x(par[0]),fc._x(par[1]));
	  }
	else
	  {
	    c.add_point(idist2,fc._x(par[0]),fc._x(par[1]));
	  }
      }

    //debug
    std::cout << "number of contour points=" << c._points.size() << std::endl;
    //debug
    
    return c;
  }

  template <class TGenoPheno>
  fcross errstats<TGenoPheno>::cross(FitFunc &func,
				     const CMAParameters<TGenoPheno> &parameters,
				     CMASolutions &cmasol,
				     const double &fup,
				     const std::vector<int> &par, const std::vector<double> &pmid,
				     const std::vector<double> &pdir, const double &ftol)
  {
    double aopt = 0.0;
    std::vector<double> alsb(3,0.0), flsb(3,0.0);
    std::vector<CMASolutions> cmasols;
    double aminsv = cmasol.best_candidate()._fvalue;
    double aim = aminsv + fup;
    double tla = ftol;
    double tlf = ftol*fup;
    int nfcn = 0;
    unsigned int maxitr = 15, ipt = 0;

    //debug
    std::cout << "aminsv=" << aminsv << " / aim=" << aim << " / tla=" << tla << " / tlf=" << tlf << std::endl;
    //debug
    
    // get a first optimized point.
    CMASolutions cmasol1 = errstats<TGenoPheno>::optimize_vpk(func,parameters,cmasol,par,pmid);
    alsb[0] = 0.0;
    flsb[0] = cmasol1.best_candidate()._fvalue;
    flsb[0] = std::max(flsb[0],aminsv+0.1*fup);
    nfcn += cmasol1._nevals;
    cmasols.push_back(cmasol1);
    ipt++;

    //debug
    std::cout << "contour / fvalue=" << cmasol1.best_candidate()._fvalue << " / optimized point1=" << cmasol1.best_candidate()._x.transpose() << std::endl;
    //debug
    
    // update aopt and get a second optimized point.
    aopt = sqrt(fup/(flsb[0]-aminsv))-1.0;
    if (aopt > 1.0)
      aopt = 1.0;
    else if (aopt < -0.5)
      aopt = 0.5;
    std::vector<double> pmiddir2 = {pmid[0]+aopt*pdir[0],pmid[1]+aopt*pdir[1]};
    CMASolutions cmasol2 = errstats<TGenoPheno>::optimize_vpk(func,parameters,cmasol,par,pmiddir2);
    alsb[1] = aopt;
    flsb[1] = cmasol2.best_candidate()._fvalue;
    nfcn += cmasol2._nevals;
    double dfda = (flsb[1]-flsb[0])/(alsb[1]-alsb[0]);
    cmasols.push_back(cmasol2);
    ipt++;

    //debug
    std::cout << "contour / fvalue=" << cmasol2.best_candidate()._fvalue << " / optimized point2=" << cmasol2.best_candidate()._x.transpose() << std::endl;
    //debug
    
    // study slope between the two points.
    if (dfda < 0.0)
      {

	//debug
	std::cout << "negative slope\n";
	//debug
	
	// if negative slope, update until positive.
	unsigned int maxlk = maxitr - ipt;
	for (unsigned int it=0;it<maxlk;it++)
	  {
	    alsb[0] = alsb[1];
	    flsb[0] = flsb[1];
	    aopt = alsb[0] + 0.2*(it+1);
	    std::vector<double> pmidt = {pmid[0]+aopt*pdir[0],pmid[1]+aopt*pdir[1]};
	    CMASolutions cmasolt = errstats<TGenoPheno>::optimize_vpk(func,parameters,cmasol,par,pmidt);
	    alsb[1] = aopt;
	    flsb[1] = cmasolt.best_candidate()._fvalue;
	    dfda = (flsb[1]-flsb[0])/(alsb[1]-alsb[0]);
	    nfcn += cmasolt._nevals;
	    cmasols.at(1) = cmasolt;
	    if (dfda > 0.0)
	      break;
	  }
      }

    // once positive slope, find a third point.
    std::cout << "positive slope, dfda=" << dfda << std::endl;
    aopt = alsb[1] + (aim-flsb[1])/dfda;
    double fdist = std::min(fabs(aim  - flsb[0]), fabs(aim  - flsb[1]));
    double adist = std::min(fabs(aopt - alsb[0]), fabs(aopt - alsb[1]));
    if (fabs(aopt) > 1.0)
      tla = ftol*fabs(aopt);
    if (adist < tla && fdist < tlf)
      {
	// return second optimized point.
	std::cout << "below tolerance, returning second optimized point\n";
	return fcross(aopt,cmasol2.best_candidate()._fvalue,
		      nfcn,cmasol2.best_candidate()._x); // XXX: returning the modified aopt, same as original code...
      }
    double bmin = std::min(alsb[0], alsb[1]) - 1.;
    if (aopt < bmin) aopt = bmin;
    double bmax = std::max(alsb[0], alsb[1]) + 1.;
    if (aopt > bmax) aopt = bmax;

    // get third point.
    std::vector<double> pmiddir3 = {pmid[0]+aopt*pdir[0],pmid[1]+aopt*pdir[1]};
    CMASolutions cmasol3 = errstats<TGenoPheno>::optimize_vpk(func,parameters,cmasol,par,pmiddir3);
    alsb[2] = aopt;
    flsb[2] = cmasol3.best_candidate()._fvalue;
    nfcn += cmasol3._nevals;
    cmasols.push_back(cmasol3);

    //debug
    std::cout << "contour / fvalue=" << cmasol3.best_candidate()._fvalue << " / optimized point3=" << cmasol3.best_candidate()._x.transpose() << std::endl;
    //debug
    
    // from three points < or > objective, decide what to do.
    double ecarmn = fabs(flsb[2] - aim);
    double ecarmx = 0.;
    unsigned int ibest = 2;
    unsigned int iworst = 0;
    unsigned int noless = 0;

    for(unsigned int i = 0; i < 3; i++)
      {
	double ecart = fabs(flsb[i] - aim);
	if(ecart > ecarmx) {
	  ecarmx = ecart;
	  iworst = i;
	}
	if(ecart < ecarmn)
	  {
	    ecarmn = ecart;
	    ibest = i;
	  }
	if(flsb[i] < aim) noless++;
      }
    
    // at least one on each side of AIM (contour)
    if(noless == 1 || noless == 2)
      {
	// XXX: could do parabola instead.
	//std::cout << "ecarmn=" << ecarmn << " / aim=" << aim << std::endl;
	int srefp = -1;
	//int sside = (ecarmn > 0.0) - (ecarmn < 0.0);
	bool bestbelow = (flsb[ibest] < aim);
	double srefecar = ecarmx;
	for (int i=0;i<3;i++)
	  {
	    if (i != ibest)
	      {
		if ((bestbelow && flsb[i] > aim)
		    || (!bestbelow && flsb[i] < aim))
		  {
		    double ecar = fabs(flsb[i]-aim);
		    if (ecar <= srefecar)
		      {
			srefecar = ecar;
			srefp = i;
		      }
		  }
	      }
	  }
	/*std::cout << "bestbelow=" << bestbelow << " / srefp=" << srefp << " / ibest=" << ibest << " / srefecar=" << srefecar << " / val=" << flsb[srefp]-aim << " / ecarmn=" << ecarmn << " / ecarmx=" << ecarmx << std::endl;

	std::cout << "a point on each side, performing linesearch\n";
	std::cout << "srefecar x=" << cmasols.at(srefp).best_candidate()._x.transpose() << std::endl;
	std::cout << "ibest x=" << cmasols.at(ibest).best_candidate()._x.transpose() << std::endl;*/
	int sdir0 = (cmasols.at(srefp).best_candidate()._x(par[0])-cmasols.at(ibest).best_candidate()._x(par[0]) > 0.0) ? 1 : -1;
	int sdir1 = (cmasols.at(srefp).best_candidate()._x(par[1])-cmasols.at(ibest).best_candidate()._x(par[1]) > 0.0) ? 1 : -1;
	/*std::cout << "ref sdir0=" << cmasols.at(srefp).best_candidate()._x(par[0])-cmasols.at(ibest).best_candidate()._x(par[0]) << " / sdir0=" << sdir0 << std::endl;
	std::cout << "ref sdir1=" << cmasols.at(srefp).best_candidate()._x(par[1])-cmasols.at(ibest).best_candidate()._x(par[1]) << " / sdir1=" << sdir1 << std::endl;*/
	dVec xstart;
	double nminfvalue;
	if (cmasols.at(ibest).best_candidate()._fvalue < aim)
	  {
	    xstart = cmasols.at(ibest).best_candidate()._x;
	    nminfvalue = cmasols.at(ibest).best_candidate()._fvalue;
	  }
	else
	  {
	    xstart = cmasols.at(srefp).best_candidate()._x;
	    nminfvalue = cmasols.at(srefp).best_candidate()._fvalue;
	  }
	double dxk0 = sdir0 * fabs(xstart(par[0])) * 0.1;
	double dxk1 = sdir1 * fabs(xstart(par[1])) * 0.1;
	
	//std::cout << "sdir0=" << sdir0 << " / sdir1=" << sdir1 << " / dxk0=" << dxk0 << " / dxk1=" << dxk1 << std::endl;
	CMASolutions citsol = cmasols.at(ibest);
	while (true)
	  {
	    // advance incrementally instead of more complicated linear_step procedure.
	    xstart(par[0]) += dxk0;
	    xstart(par[1]) += dxk1;

	    //debug
	    //std::cout << "xstart=" << xstart.transpose() << std::endl;
	    //debug
	    
	    std::vector<double> vxk = {xstart(par[0]),xstart(par[1])};
	    CMASolutions ncitsol = errstats<TGenoPheno>::optimize_vpk(func,parameters,citsol,par,vxk);
	    //std::cout << "optimization status=" << ncitsol._run_status << std::endl;
	    if (ncitsol._run_status < 0)
	      {
		LOG(WARNING) << "contour linesearch: optimization error " << ncitsol._run_status << std::endl;
	      }
	    nminfvalue = ncitsol.best_candidate()._fvalue;
	    if (nminfvalue > aim)
	      break;
	    citsol = ncitsol;
	  }
	//std::cout << "contour linesearch best point=" << citsol.best_candidate()._x.transpose() << std::endl;
	return fcross(0.0,citsol.best_candidate()._fvalue,nfcn,citsol.best_candidate()._x);
      }
    // if all three are above AIM, third point must be the closest to AIM, return it
    if(noless == 0 && ibest != 2)
      {
	std::cout << "all points above, returning third point as best\n";
	return fcross(aopt,flsb[2],
		      nfcn,cmasols[2].best_candidate()._x);
      }
    // if all three below and third is not best then the slope has again gone negative,
    // re-iterate and look for positive slope
    if(noless == 3 && ibest != 2)
      {
	//TODO.
	std::cout << "slope is again negative\n";
      }

    // in other case new straight line thru first two points
    std::cout << "new straight line through first two points\n";
    flsb[iworst] = flsb[2];
    alsb[iworst] = alsb[2];
    dfda = (flsb[1] - flsb[0])/(alsb[1] - alsb[0]);
    //TODO: restart from point 3...
    
    return fcross();
  }
  
  template class errstats<GenoPheno<NoBoundStrategy>>;
  template class errstats<GenoPheno<pwqBoundStrategy>>;
  template class errstats<GenoPheno<NoBoundStrategy,linScalingStrategy>>;
  template class errstats<GenoPheno<pwqBoundStrategy,linScalingStrategy>>;
}
