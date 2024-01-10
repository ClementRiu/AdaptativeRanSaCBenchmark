// @INPROCEEDINGS{5206678,
//  author={Jongmoo Choi and Medioni, Gerard},
//  booktitle={2009 IEEE Conference on Computer Vision and Pattern Recognition},
//  title={StaRSaC: Stable random sample consensus for parameter estimation},
//  year={2009},
//  volume={},
//  number={},
//  pages={675-682},
//  doi={10.1109/CVPR.2009.5206678}}
//
// Created by Clément Riu.
//

#ifndef MMM_ORSA_STARSAC_H
#define MMM_ORSA_STARSAC_H

#include "lrtsac.hpp"

namespace orsa {
class StarSac : public LRTSac {
public:
    /// Constructor
    explicit StarSac(const ModelEstimator *estimator);

    /// Generic implementation of RANSAC
    double run(RunResult &res, int nIterMax = 1000, bool verbose = false) const;

    void setHyperParameters(double cpIIT = 0.99, double sigmaMax = 16);

private:
    // Verifies the runOutput metric is good enough.
    bool satisfyingRun(double runOutput) const;

    const int _numRansacRun = 5;

    double computeVar(const std::vector<RunResult> &estimatedModels) const;

    void FindAllInliers(const ModelEstimator::Model &model,
                        const std::vector<double> &precisions,
                        std::vector<std::vector<int>> &inliers) const;

    int ransac(const std::vector<double> &precisions,
               std::vector<RunResult> &res,
               int nIterMax) const;
};

} // namespace orsa

#endif //MMM_ORSA_STARSAC_H
