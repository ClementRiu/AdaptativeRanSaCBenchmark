// https://github.com/danini/magsac
//@inproceedings{barath2019magsac,
//author = {Barath, Daniel and Matas, Jiri and Noskova, Jana},
//        title = {MAGSAC: marginalizing sample consensus},
//        booktitle = {Conference on Computer Vision and Pattern Recognition},
//        year = {2019},
//}
//
//@inproceedings{barath2019magsacplusplus,
//author = {Barath, Daniel and Noskova, Jana and Ivashechkin, Maksym and Matas, Jiri},
//        title = {MAGSAC++, a fast, reliable and accurate robust estimator},
//        booktitle = {arXiv preprint:1912.05909},
//        year = {2019},
//}
//
// Adapted by Clément Riu.
//

#ifndef MMM_ORSA_MAGSAC_H
#define MMM_ORSA_MAGSAC_H

#include "ransac_algorithm.hpp"
#include "model_estimator.hpp"

namespace orsa {

class ModelScore {
public:
    /* number of inliers, rectangular gain function */
    size_t inlier_number;
    /* MSAC scoring, truncated quadratic gain function */
    double score;
    /* The log probability of the model considering that the inliers are normally and the outliers are uniformly distributed */
    double probability;
    /* Iteration number when it is found */
    size_t iteration;

    ModelScore() : inlier_number(0), score(0), probability(0), iteration(0) {}

    inline bool operator<(const ModelScore &score_) {
        return score < score_.score;
    }
};

using Score = ModelScore;

class MAGSAC : public RansacAlgorithm {
public:
    enum Version {
        // The original version of MAGSAC. It works well, however, can be quite slow in many cases.
        MAGSAC_ORIGINAL,
        // The recently proposed MAGSAC++ algorithm which keeps the accuracy of the original MAGSAC but is often orders of magnitude faster.
        MAGSAC_PLUS_PLUS
    };

    MAGSAC(const ModelEstimator *estimator, const Version magsac_version_ = Version::MAGSAC_ORIGINAL);

    ~MAGSAC() {}

    //// Setters for RANSAC hyperparameters.
    void setHyperParameters(double confidence, double maximumThreshold = 10.0, size_t partitionNumber = 5,
                            const double fps = -1, const double refThreshold = 2.0);

    void setConfidence(double confidence);

    void setMaximumThreshold(double maximumThreshold);

    void setRefThreshold(double refThreshold);

    void setPartitionNumber(size_t partitionNumber);

    void setFPS(double fps_);

    // A function to run MAGSAC.
    double run(RunResult &res, int maxiter, bool verbose) const; // The score of the estimated model

    bool sigmaConsensus(const ModelEstimator::Model &model_,
                        ModelEstimator::Model &refined_model_,
                        ModelScore &score_,
                        const ModelScore &best_score_,
                        int &last_iteration_number,
                        std::vector<int> &inliersIdxsSaved,
                        std::vector<double> &weightsSaved) const;

    bool sigmaConsensusPlusPlus(const ModelEstimator::Model &model_,
                                ModelEstimator::Model &refined_model_,
                                ModelScore &score_,
                                const ModelScore &best_score_,
                                int &last_iteration_number,
                                std::vector<int> &inliersIdxsSaved,
                                std::vector<double> &weightsSaved) const;

    bool satisfyingRun(double runOutput) const;

    size_t number_of_irwls_iters;

    inline std::vector<double> getWeights() const { return weightsSaved_; };
private:
    static constexpr double getSigmaQuantile() {
        return 3.64;
    }

    void getModelQuality(const ModelEstimator::Model &model_, // The model parameter
                         double &marginalized_iteration_number_, // The marginalized iteration number to be calculated
                         double &score_) const;

    void getModelQualityPlusPlus(const ModelEstimator::Model &model_, // The model parameter
                                 double &score_,
                                 const double &previous_best_score_) const;

    double _confidence;

    Version magsac_version; // The version of MAGSAC used
    size_t mininum_iteration_number; // Minimum number of iteration before terminating
    double maximum_threshold; // The maximum sigma value
    size_t core_number; // Number of core used in sigma-consensus
    double _time_limit; // A time limit after the algorithm is interrupted
    double _desired_fps; // The desired FPS (TODO: not tested with MAGSAC)
    bool apply_post_processing; // Decides if the post-processing step should be applied
    double _log_confidence; // The logarithm of the required confidence
    size_t partition_number; // Number of partitions used to speed up sigma-consensus
    double interrupting_threshold; // A threshold to speed up MAGSAC by interrupting the sigma-consensus procedure whenever there is no chance of being better than the previous so-far-the-best model
    mutable std::vector<double> weightsSaved_;

};

} // namespace orsa

#endif //MMM_ORSA_MAGSAC_H
