// Adapted from https://github.com/vxl/vxl/tree/c6c899aaf9cbf7ccabfca1ba35b2248fb611ffbc
// Algorithm from MUSE: Robust Surface Fitting using Unbiased Scale Estimates, James V. Miller and Charles V. Stewart
// By Clément Riu
// 2021

#ifndef MUSE_H
#define MUSE_H

#include "ransac_algorithm.hpp"
#include "model_estimator.hpp"
#include "muse_lookup.hpp"

namespace orsa {

/// Model estimation with ORSA algorithm.
class Muse : public RansacAlgorithm {
public:
    /// Constructor
    explicit Muse(const ModelEstimator *estimator);

    //// Setters for Muse hyperparameters.
    void setHyperParameters(double cpIIT = 0.99,
                            double min_frac = 0.25,
                            double max_frac = 0.95,
                            double frac_inc = 0.05,
                            bool use_sk_refine = true);
    void setCpIIT(double cpIIT);

    void setMinFrac(double min_frac);

    void setMaxFrac(double max_frac);

    void setFracInc(double frac_inc);

    void setSKRefine(bool use_sk_refine);

    /// Generic implementation of Muse
    double run(RunResult &res, int nIterMax = 1000, bool verbose = false) const;

    bool satisfyingRun(double runOutput) const;

protected:

private:
    double _cpIIT;
    double _min_frac, _max_frac, _frac_inc;
    MuseTable *_table;
    bool _use_sk_refine;

    double computeObjective(std::vector<double> &residuals) const;
};

}  // namespace orsa

#endif
