/**
* @file lrtsac.hpp
* @brief Ransac variant based on likelihood
* @author Clement Riu, Pascal Monasse
*
* Copyright (c) 2020-2021 Clement Riu
* Copyright (c) 2021 Pascal Monasse
* All rights reserved.
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>

#include "lrtsac.hpp"
#include "sampling.hpp"

// The algorithm is described in
// Automatic RANSAC by Likelihood Maximization
// C. Riu, V. Nozick, P. Monasse
// IPOL
// Comments below refer to this article.

/// Chi2 table: {{dof,p},F^{-1}(p)} with F the cumulative Chi2 distribution.
/// Used dof: 2+ (Line2P,Essential5P,PNP6P, Fundamental7P,Fundemental8P,9:???)
static std::map<std::pair<int, double>, double> chi2Table = {
    {{2+2, 0.90},7.779},{{2+2, 0.91},8.043},{{2+2, 0.92},8.337},{{2+2, 0.93},8.666},{{2+2, 0.94},9.044},{{2+2, 0.95},9.488},{{2+2, 0.96},10.026},{{2+2, 0.97},10.712},{{2+2, 0.98},11.668},{{2+2, 0.99},13.277},
    {{5+2, 0.90},12.017},{{5+2, 0.91},12.337},{{5+2, 0.92},12.691},{{5+2, 0.93},13.088},{{5+2, 0.94},13.540},{{5+2, 0.95},14.067},{{5+2, 0.96},14.703},{{5+2, 0.97},15.509},{{5+2, 0.98},16.622},{{5+2, 0.99},18.475},
    {{6+2, 0.90},13.3616},{{6+2, 0.91},13.6975},{{6+2, 0.92},14.0684},{{6+2, 0.93},14.4836},{{6+2, 0.94},14.9563},{{6+2, 0.95},15.5073},{{6+2, 0.96},16.1708},{{6+2, 0.97},17.0105},{{6+2, 0.98},18.1682},{{6+2, 0.99},20.0902 },
    {{7+2, 0.90},14.684},{{7+2, 0.91},15.034},{{7+2, 0.92},15.421},{{7+2, 0.93},15.854},{{7+2, 0.94},16.346},{{7+2, 0.95},16.919},{{7+2, 0.96},17.608},{{7+2, 0.97},18.480},{{7+2, 0.98},19.679},{{7+2, 0.99},21.666},
    {{8+2, 0.90},15.987},{{8+2, 0.91},16.352},{{8+2, 0.92},16.753},{{8+2, 0.93},17.203},{{8+2, 0.94},17.713},{{8+2, 0.95},18.307},{{8+2, 0.96},19.021},{{8+2, 0.97},19.922},{{8+2, 0.98},21.161},{{8+2, 0.99},23.209},
    {{9+2, 0.90},17.275},{{9+2, 0.91},17.653},{{9+2, 0.92},18.069},{{9+2, 0.93},18.533},{{9+2, 0.94},19.061},{{9+2, 0.95},19.675},{{9+2, 0.96},20.412},{{9+2, 0.97},21.342},{{9+2, 0.98},22.618},{{9+2, 0.99},24.725}
};

namespace orsa {
/// Constructor of the LRTSAC algorithm
/// \param model Model to evaluate with the data
LRTSac::LRTSac(const ModelEstimator *model) : RansacAlgorithm(model) {
    _sigmaMin = 0.25;
    _B = 100; // Should be adjusted to balance bailout test and processing time
    setHyperParameters();
}

/// Set all LRTSac hyper-paramters.
/// \param cpI   Confidence proba wrt type I error
/// \param cpIIB Confidence proba wrt type II error (bailout test)
/// \param cpIIT Confidence proba wrt type II error (limited iterations)
/// \param sigmaMax Parameter controlling the max error considered
void LRTSac::setHyperParameters(double cpI, double cpIIB,
                                double cpIIT, double sigmaMax,
                                bool reduceSigma) {
    setCpI(cpI);
    setCpIIB(cpIIB);
    setCpIIT(cpIIT);
    setSigmaMax(sigmaMax);
    setReduceSigma(reduceSigma);
}

/// Sets the value of cpI.
/// A value of 0 accepts anything.
/// \param cpI Confidence proba wrt type I error
void LRTSac::setCpI(double cpI) {
    assert(cpI<1);
    if (cpI > 0) {
        // minL computed by inverting chi2 cumulative distribution
        std::pair<int, double> param = {_model->nDegreeOfFreedom() + 2, cpI};
        if (chi2Table.find(param) == chi2Table.end())
            throw std::invalid_argument("LRTSac's chi2 value not tabulated");
        _minL = chi2Table[param]/(2*_model->NbData()); // (13)
    } else
        _minL = 0;
}

/// Sets the value of cpIIB.
/// A value of 1 disables the bailout test.
/// \param cpIIB Probability of type II-error in bailout test
void LRTSac::setCpIIB(double cpIIB) { // Set to 0 if useBailout == false.
    assert(cpIIB > 0);
    _cpIIB = cpIIB;
    if (_cpIIB > 1.0) {
        std::cerr<<"LRTSac cpIIB adjusted to not exceed 1"<<std::endl;
        _cpIIB = 1.0;
    }
}

/// Sets the value of cpIIT.
/// A value of 1 disables the dynamic adjustment of the iteration number.
/// \param cpIIT Confidence proba wrt type II error (limited iterations)
void LRTSac::setCpIIT(double cpIIT) {
    assert(cpIIT > 0);
    _cpIIT = cpIIT;
    if (_cpIIT > 1.0) {
        std::cerr<<"LRTSac cpIIT adjusted to not exceed 1"<<std::endl;
        _cpIIT = 1.0;
    }
}

//// Sets the value of sigmaMax.
/// \param sigmaMax Max error considered.
void LRTSac::setSigmaMax(double sigmaMax) {
    _sigmaMax = sigmaMax;
    if (sigmaMax < _sigmaMin)
        _sigmaMin = sigmaMax;
}

void LRTSac::setReduceSigma(bool reduceSigma) {
    _reduceSigma = reduceSigma;
}

/// Print info during the run
void LRTSac::printRunInfo(int iter, int T, double lambda, double eps,
                          double sigma, const std::vector<double> &Sigma)const{
    std::cout << "  L=" << lambda <<  " inliers=" << int(eps*_model->NbData())
              << " precision=" << sigma
              << " iter=" << iter << "/" << T << " Sigma={"
              << Sigma.front() << "..." << Sigma.back() << '}' << std::endl;
}

/// Fill set of sigma values.
/// Geometric progression from _sigmaMin to _sigmaMax.
void LRTSac::initSigma(std::vector<double>& Sigma) const {
    const double sigmaMultiplier=sqrt(2.0);
    Sigma.push_back(_sigmaMin);
    while(true) {
        double sigma = Sigma.back()*sigmaMultiplier;
        if(sigma>_sigmaMax) break;
        Sigma.push_back(sigma);
    }
    if(Sigma.back() != _sigmaMax)
        Sigma.push_back(_sigmaMax);
}

/// Computation of the log-likelihood function. Equation (10)
double LRTSac::likelihood(double eps, double sigma) const {
    double p=_model->pSigma(sigma), q=1-p;
    if (p<1.0e-10 || q<1.0e-10)
        throw std::domain_error("Error likelihood:\n"
                                "pSigma too small or too close to 1.");
    if (eps < p)
        return 0;

    return (eps==1)? -std::log(p):
                     eps*std::log(eps/p) + (1-eps)*std::log((1-eps)/q);
}

/// Bisection based on the likelihood function (in algorithm 3).
/// Find inlier ratio at given \a sigma to reach log-likelihood value \a L.
/// Granularity is 1/NbData.
double LRTSac::bisectLikelihood(double sigma, double L) const {
    double iMin=0, iMax=1;
    double LMin = likelihood(iMin, sigma);
    double LMax = likelihood(iMax, sigma);

    if (L <= LMin) return iMin;
    if (L >= LMax) return iMax;

    while((iMax-iMin)*_model->NbData() > 1.0) {
        double iMid = (iMin+iMax)*0.5;
        double LMid = likelihood(iMid, sigma);
        assert(LMin<=LMid && LMid<=LMax);
        if (L < LMid) {
            iMax = iMid;
            LMax = LMid;
        } else {
            iMin = iMid;
            LMin = LMid;
        }
    }
    return iMin;
}

/// Computation of the inlier ratios (\a eps) for each sigma (algorithm 4).
/// Early bailout may occur if the model is unlikely a better one.
/// \param model Current model to test
/// \param Sigma Set of possible values for sigma
/// \param[out] eps Inlier ratio for each value of sigma
/// \param[out] vpm Number of applied verifications
/// \param epsMin Min eps value for better model (used only for bailout)
/// \return Indicate whether eps is exact (no early bailout)
bool LRTSac::computeEps(const ModelEstimator::Model &model,
                        const std::vector<double> &Sigma,
                        std::vector<double> &eps, int &vpm,
                        const std::vector<double> &epsMin) const {
    const int n = _model->NbData();
    const double increment = 1.0/n;
    for (int j=0, bailCount=0; j<n; j++) {
        double error = _model->Error(model, j);
        vpm++;

        for (size_t i=0; i<Sigma.size(); i++)
            if (error <= Sigma[i] * Sigma[i])
                eps[i] += increment;

        if (_cpIIB<1 && ++bailCount==_B) {
            bailCount=0; // Round counter, cheaper than Euclidean division
            double tau = std::sqrt(-(std::log(1-_cpIIB)
                                     - std::log(std::floor(n/_B)))
                                   / (2*(j+1))); // (19)
            bool bailout = true;
            for (size_t i=0; bailout && i<Sigma.size(); i++)
                if (eps[i]*n >= (j+1)*(epsMin[i]-tau))
                    bailout = false;
            if (bailout)
                return false;
        }
    }
    return true;
}

/// Find sigma leading to best log-likelihood based on inlier ratios.
/// Algorigthm 2, line 7.
/// \param Sigma Set of possible values for sigma
/// \param eps Inlier ratio for each value of sigma
/// \param[out] L The highest log-likelihood
/// \param[out] bestEps Inlier ratio for best sigma
/// \return The value of the best sigma
double LRTSac::bestSigma(const std::vector<double> &Sigma,
                         const std::vector<double> &eps,
                         double &L, double &epsBest) const {
    double sigma=0;
    L = -1.0;
    for (size_t i=0; i<Sigma.size(); i++) {
        double lambda = likelihood(eps[i], Sigma[i]);
        if (lambda > L) {
            L = lambda;
            sigma = Sigma[i];
            epsBest = eps[i];
        }
    }
    return sigma;
}

/// Compute min epsilon for each sigma to reach log-likelihood \a L.
/// This may also reduce the maximum sigma (algorithm 3).
void LRTSac::computeEpsMin(std::vector<double> &Sigma,
                           std::vector<double> &epsMin,
                           double L) const {
    std::vector<double>::iterator it=Sigma.begin();
    for (int i=0; it!=Sigma.end(); ++it, ++i) {
        if (likelihood(1,*it) <= L)
            break;
        epsMin[i] = (likelihood(0,*it)>=L)? 0: bisectLikelihood(*it,L);
    }
    if(_reduceSigma)
        Sigma.erase(it, Sigma.end());
}

/// Minimal number of iterations for a given inlier ratio \a eps (algorithm 5)
double LRTSac::computeIter(double eps) const {
    double num = std::log(1-_cpIIT);
    double den = std::pow(eps, _model->SizeSample());
    den = std::log(1 - _cpIIB*den);
    if (den==0) // Happens if _cpIIB*den<<1
        num=1; // So as to return +infty
    return num/den;
}

/// Computes the best model with the LRTSAC algorithm (algorithm 2).
/// Return the maximum log-likelihood.
double LRTSac::run(RunResult& res, int nIterMax, bool verbose) const {
    const int nData = _model->NbData();
    const int sizeSample = _model->SizeSample();

    double Lbest = 0;
    res.sigma = _sigmaMax; res.T=0; res.vpm=0;
    if (nData <= sizeSample)
        return Lbest;

    if(nIterMax<=0) // Authorize "infinite" number of iterations
        nIterMax = std::numeric_limits<int>::max();

    // Computation of array of values for sigma
    std::vector<double> Sigma;
    initSigma(Sigma);

    std::vector<double> epsMin(Sigma.size(), 0.0);
    int T = nIterMax;
    if (_minL > 0) {
        if (_cpIIT<1 || _cpIIB<1 || _reduceSigma)
            computeEpsMin(Sigma, epsMin, _minL);
        if (_cpIIT<1)
            T = (int)std::min((double)T, computeIter(epsMin.front()));
    }

    if(verbose) {
        std::cout << "(init)";
        printRunInfo(-1, T, _minL, 0, res.sigma, Sigma);
    }

    int modelCount = 0;
    for (; res.T<T && !Sigma.empty(); res.T++) {
        std::vector<int> vSample(sizeSample);
        UniformSample(sizeSample, nData, &vSample);
        std::vector<ModelEstimator::Model> vModels;
        _model->Fit(vSample, &vModels);

        std::vector<ModelEstimator::Model>::const_iterator it;
        for (it=vModels.begin(); it!=vModels.end(); ++it) {
            ++modelCount;
            std::vector<double> eps(Sigma.size(),0); // Inlier ratios
            int vpm = 0;
            bool noBailout = computeEps(*it, Sigma, eps, vpm, epsMin);
            res.vpm += vpm;
            if(! noBailout)
                continue;

            double L, epsBest=0;
            double sigma = bestSigma(Sigma, eps, L, epsBest);
            if (L > Lbest) {
                Lbest = L;
                res.model = *it;
                res.sigma = sigma;

                if (_cpIIT<1 || _cpIIB<1 || _reduceSigma)
                    computeEpsMin(Sigma, epsMin, Lbest);
                if (_cpIIT<1 && !Sigma.empty())
                    T = (int)std::min((double)T, computeIter(epsMin.front()));

                if(verbose)
                    printRunInfo(res.T, T, Lbest, epsBest, sigma, Sigma);
            }
        }
    }
    res.vpm /= (double) modelCount;

    // Computation of inliers based on the best sigma found
    _model->FindInliers(res.model, res.sigma, res.vInliers);

    return Lbest;
}

/// Verify the runOutput metric is good enough.
bool LRTSac::satisfyingRun(double runOutput) const {
    return runOutput>=_minL;
}

}
