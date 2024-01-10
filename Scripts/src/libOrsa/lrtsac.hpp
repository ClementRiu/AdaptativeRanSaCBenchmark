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

#ifndef MMM_ORSA_LRTSAC_HPP
#define MMM_ORSA_LRTSAC_HPP

#include "ransac_algorithm.hpp"

// Model estimation using the LRTSAC algorithm.
// Original method:
// @inproceedings{cohen2015likelihood,
// title={The likelihood-ratio test and efficient robust estimation},
//        author={Cohen, Andrea and Zach, Christopher},
//        booktitle={Proceedings of ICCV},
//        pages={2282--2290},
//        year={2015}
// }


namespace orsa {
class LRTSac : public RansacAlgorithm {
public:
    explicit LRTSac(const ModelEstimator *model);

    /// Setters for LRTSac hyper-parameters.
    void setHyperParameters(double cpI = 0, double cpIIB = 0.95,
                            double cpIIT = 0.99, double sigmaMax = 16,
                            bool reduceSigma = true);

    void setCpI(double cpI);     ///< Type I confidence
    void setCpIIB(double cpIIB); ///< Type II confidence wrt early bailout
    void setCpIIT(double cpIIT); ///< Type II confidence wrt iterations
    void setSigmaMax(double sigmaMax);
    void setReduceSigma(bool reduceSigma);

    // Computes the best model based on the saved parameters.
    double run(RunResult& res, int nIterMax=0, bool verbose=false) const;

private:
    // Verifies the runOutput metric is good enough.
    bool satisfyingRun(double runOutput) const;
    double _cpIIB; ///< Confidence proba wrt type II error due to bailout
    double _minL; ///< Min log-likelihood deduced from alpha
    bool _reduceSigma; ///< Discard inactive upper values of Sigma

    int _B; ///< Frequency to check for early bailout.

protected:
    // Values of sigma to try.
    double _sigmaMin, _sigmaMax;
    double _cpIIT; ///< Confidence proba wrt type II error due to iterations

    void initSigma(std::vector<double>& Sigma) const;

private:
    void printRunInfo(int iter, int T, double lambda, double eps,
                      double sigma, const std::vector<double> &Sigma) const;
    double likelihood(double eps, double sigma) const;
    double bisectLikelihood(double sigma, double L) const;
    void computeEpsMin(std::vector<double> &Sigma, std::vector<double> &epsMin,
                       double L) const;
    bool computeEps(const ModelEstimator::Model &model,
                             const std::vector<double> &Sigma,
                             std::vector<double> &eps, int &vpm,
                             const std::vector<double> &epsMin) const;
    double computeIter(double eps) const;
    double bestSigma(const std::vector<double> &Sigma,
                     const std::vector<double> &eps,
                     double &L, double &epsBest) const;
};
}

#endif //MMM_ORSA_LRTSAC_HPP
