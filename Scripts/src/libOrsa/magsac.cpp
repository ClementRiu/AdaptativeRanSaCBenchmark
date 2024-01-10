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

#include <limits>
#include <chrono>
#include <cmath>
#include <algorithm>

#include "magsac.hpp"
#include "sampling.hpp"
#include "gamma_values.cpp"

namespace orsa {

MAGSAC::MAGSAC(const ModelEstimator *estimator, const Version magsac_version_) :
        RansacAlgorithm(estimator),
        number_of_irwls_iters(1),
        magsac_version(magsac_version_),
        mininum_iteration_number(50),
        maximum_threshold(10.0),
        core_number(1),
        _time_limit(std::numeric_limits<double>::max()),
        _desired_fps(-1),
        apply_post_processing(true),
        _log_confidence(0),
        partition_number(10),
        interrupting_threshold(2.0) {
}

/// Setters for MAGSAC hyperparameters.
void
MAGSAC::setHyperParameters(const double confidence, const double maximumThreshold, const size_t partitionNumber,
                           const double fps, const double refThreshold) {
    setConfidence(confidence);
    setMaximumThreshold(maximumThreshold);
    setPartitionNumber(partitionNumber);
    setFPS(fps);
    setRefThreshold(refThreshold);
}

void MAGSAC::setConfidence(const double confidence) {
    _confidence = confidence;
    _log_confidence = log(1.0 - confidence);
}

void MAGSAC::setMaximumThreshold(const double maximumThreshold) {
    maximum_threshold = maximumThreshold;
}

void MAGSAC::setRefThreshold(const double refThreshold) {
    interrupting_threshold = refThreshold;
}

void MAGSAC::setPartitionNumber(const size_t partitionNumber) {
    partition_number = partitionNumber;
}

// A function to set a desired minimum frames-per-second (FPS) value.
void MAGSAC::setFPS(const double fps_) {
    _desired_fps = fps_; // The required FPS.
    // The time limit which the FPS implies
    _time_limit = fps_ <= 0 ?
                  std::numeric_limits<double>::max() :
                  1.0 / fps_;
}

double MAGSAC::run(RunResult &res, const int maxiter, bool verbose) const {
    // Initialize variables
    std::chrono::time_point<std::chrono::system_clock> start, end; // Variables for time measuring: start and end times
    std::chrono::duration<double> elapsed_seconds; // Variables for time measuring: elapsed time
    const int point_number = _model->NbData(); // Number of points
    const int sample_size = _model->SizeSample(); // The sample size required for the estimation
    size_t max_iteration = maxiter; // The maximum number of iterations initialized to the iteration limit
    size_t iteration = 0; // Current number of iterations
    ModelScore so_far_the_best_score; // The score of the current best model
    std::vector<int> minimal_sample(sample_size); // The sample used for the estimation

    std::vector<size_t> pool(point_number);
    for (size_t point_idx = 0; point_idx < (size_t) point_number; ++point_idx)
        pool[point_idx] = point_idx;

    if (point_number < sample_size) {
        fprintf(stderr,
                "There are not enough points for applying robust estimation. Minimum is %d; while %d are given.\n",
                sample_size, point_number);
        return false;
    }

    // Set the start time variable if there is some time limit set
    if (_desired_fps > -1)
        start = std::chrono::system_clock::now();

    std::vector<int> inliersIdxsToSave;
    std::vector<double> weightsToSave;

    // Main MAGSAC iteration
    while (mininum_iteration_number > iteration ||
           iteration < max_iteration) {
        // Increase the current iteration number
        ++iteration;
        // Sample a minimal subset
        std::vector<ModelEstimator::Model> models; // The set of estimated models
        // Try to select a minimal sample and estimate the implied model parameters
        // Get a minimal sample randomly
        UniformSample(sample_size, point_number, &minimal_sample);
        _model->Fit(minimal_sample, &models); // The estimated models

        // Select the so-far-the-best from the estimated models
        std::vector<ModelEstimator::Model>::const_iterator it;
        for (it = models.begin(); it != models.end(); ++it) {
            ModelScore score; // The score of the current model
            ModelEstimator::Model refined_model; // The refined model parameters

            // Apply sigma-consensus to refine the model parameters by marginalizing over the noise level sigma
            bool success;
            int last_iteration_number = 0;
            if (magsac_version == Version::MAGSAC_ORIGINAL) {
                success = sigmaConsensus(*it,
                                         refined_model,
                                         score,
                                         so_far_the_best_score,
                                         last_iteration_number,
                                         inliersIdxsToSave,
                                         weightsToSave);
            } else {
                success = sigmaConsensusPlusPlus(*it,
                                                 refined_model,
                                                 score,
                                                 so_far_the_best_score,
                                                 last_iteration_number,
                                                 inliersIdxsToSave,
                                                 weightsToSave);
            }
            // Continue if the model was rejected
            if (!success || score.score == -1)
                continue;

            // Save the iteration number when the current model is found
            score.iteration = iteration;

            // Update the best model parameters if needed
            if (so_far_the_best_score < score) {
                res.model = refined_model; // Update the best model parameters
                so_far_the_best_score = score; // Update the best model's score
                max_iteration = std::min(max_iteration,
                                         (size_t) last_iteration_number); // Update the max iteration number, but do not allow to increase
                if (verbose) {
                    std::cout << "Score: " << score.score << " - Iterations: " << iteration << " / "
                              << max_iteration << std::endl;
                }
                res.vInliers = inliersIdxsToSave;
                weightsSaved_ = weightsToSave;
            }
        }

        // Update the time parameters if a time limit is set
        if (_desired_fps > -1) {
            end = std::chrono::system_clock::now();
            elapsed_seconds = end - start;

            // Interrupt if the time limit is exceeded
            if (elapsed_seconds.count() > _time_limit)
                break;
        }
    }

    // Apply sigma-consensus as a post processing step if needed and the estimated model is valid
    if (apply_post_processing) {
        // TODO
    }

    res.T = iteration;
    res.sigma = maximum_threshold;
    res.vpm = point_number;
    return so_far_the_best_score.score;
}

bool MAGSAC::sigmaConsensus(const ModelEstimator::Model &model_,
                            ModelEstimator::Model &refined_model_,
                            ModelScore &score_,
                            const ModelScore &best_score_,
                            int &last_iteration_number,
                            std::vector<int> &inliersIdxsSaved,
                            std::vector<double> &weightsSaved) const {
    // Set up the parameters
    const double L = 1.05;
    const double k = getSigmaQuantile();
    const double threshold_to_sigma_multiplier = 1.0 / k;
    const size_t sample_size = _model->SizeSample();
    static auto comparator = [](std::pair<double, int> left, std::pair<double, int> right) {
        return left.first < right.first;
    };
    const int point_number = _model->NbData();
    double current_maximum_sigma = maximum_threshold;

    // Calculating the residuals
    std::vector<std::pair<double, size_t> > all_residuals;
    all_residuals.reserve(point_number);

    // If it is not the first run, consider the previous best and interrupt the validation when there is no chance of being better
    if (best_score_.inlier_number > 0) {
        // Number of inliers which should be exceeded
        int points_remaining = best_score_.inlier_number;

        // Collect the points which are closer than the threshold which the maximum sigma implies
        for (int point_idx = 0; point_idx < point_number; ++point_idx) {
            // Calculate the residual of the current point
            const double residual = sqrt(_model->Error(model_, point_idx));
            if (current_maximum_sigma > residual) {
                // Store the residual of the current point and its index
                all_residuals.emplace_back(std::make_pair(residual, point_idx));

                // Count points which are closer than a reference threshold to speed up the procedure
                if (residual < interrupting_threshold)
                    --points_remaining;
            }

            // Interrupt if there is no chance of being better
            // TODO: replace this part by SPRT test
            if (point_number - point_idx < points_remaining)
                return false;
        }

        // Store the number of really close inliers just to speed up the procedure
        // by interrupting the next verifications.
        score_.inlier_number = best_score_.inlier_number - points_remaining;
    } else {
        // The number of really close points
        size_t points_close = 0;

        // Collect the points which are closer than the threshold which the maximum sigma implies
        for (size_t point_idx = 0; point_idx < (size_t) point_number; ++point_idx) {
            // Calculate the residual of the current point
            const double residual = sqrt(_model->Error(model_, point_idx));
            if (current_maximum_sigma > residual) {
                // Store the residual of the current point and its index
                all_residuals.emplace_back(std::make_pair(residual, point_idx));

                // Count points which are closer than a reference threshold to speed up the procedure
                if (residual < interrupting_threshold)
                    ++points_close;
            }
        }

        // Store the number of really close inliers just to speed up the procedure
        // by interrupting the next verifications.
        score_.inlier_number = points_close;
    }

    std::vector<ModelEstimator::Model> sigma_models;
    std::vector<int> sigma_inliers;
    std::vector<double> final_weights;

    // The number of possible inliers
    const size_t possible_inlier_number = all_residuals.size();

    // Sort the residuals in ascending order
    std::sort(all_residuals.begin(), all_residuals.end(), comparator);

    // The maximum threshold is set to be slightly bigger than the distance of the
    // farthest possible inlier.
    current_maximum_sigma =
            all_residuals.back().first + std::numeric_limits<double>::epsilon();

    const double sigma_step = current_maximum_sigma / partition_number;

    last_iteration_number = 10000;

    score_.score = 0;

    // The weights calculated by each parallel process
    std::vector<std::vector<double>> point_weights_par(partition_number,
                                                       std::vector<double>(possible_inlier_number, 0));

    // If OpenMP is used, calculate things in parallel
//#ifdef USE_OPENMP
#pragma omp parallel for num_threads(core_number)
    for (int partition_idx = 0; partition_idx < partition_number; ++partition_idx) {
        // The maximum sigma value in the current partition
        const double max_sigma = (partition_idx + 1) * sigma_step;

        // Find the last element which has smaller distance than 'max_threshold'
        // Since the vector is ordered binary search can be used to find that particular element.
        const auto &last_element = std::upper_bound(all_residuals.begin(), all_residuals.end(),
                                                    std::make_pair(max_sigma, 0), comparator);
        const size_t sigma_inlier_number = last_element - all_residuals.begin();

        // Put the indices into a vector
        std::vector<int> sigma_inliers;
        sigma_inliers.reserve(sigma_inlier_number);

        // Store the points which are closer than the current sigma limit
        for (size_t relative_point_idx = 0; relative_point_idx < sigma_inlier_number; ++relative_point_idx)
            sigma_inliers.emplace_back(all_residuals[relative_point_idx].second);

        // Check if there are enough inliers to fit a model
        if (sigma_inliers.size() > sample_size) {
            // Estimating the model which the current set of inliers imply
            std::vector<ModelEstimator::Model> sigma_models;
            _model->Fit(sigma_inliers, &sigma_models);
            // If the estimation was successful calculate the implied probabilities
            if (sigma_models.size() == 1) {
                const double max_sigma_squared_2 = 2 * max_sigma * max_sigma;
                double residual_i_2, // The residual of the i-th point
                probability_i; // The probability of the i-th point

                // Iterate through all points to estimate the related probabilities
                for (size_t relative_point_idx = 0;
                     relative_point_idx < sigma_inliers.size(); ++relative_point_idx) {
                    // TODO: Replace with Chi-square instead of normal distribution
                    const size_t &point_idx = sigma_inliers[relative_point_idx];

                    // Calculate the residual of the current point
                    residual_i_2 = _model->Error(sigma_models[0], point_idx);
                    // Calculate the probability of the i-th point assuming Gaussian distribution
                    // TODO: replace by Chi-square distribution
                    probability_i = exp(-residual_i_2 / max_sigma_squared_2);

                    // Store the probability of the i-th point coming from the current partition
                    point_weights_par[partition_idx][relative_point_idx] += probability_i;


                }
            }
        }
    }
//#else
//    fprintf(stderr, "Not implemented yet.\n");
//#endif
    // The weights used for the final weighted least-squares fitting
    final_weights.reserve(possible_inlier_number);

    // Collect all points which has higher probability of being inlier than zero
    sigma_inliers.reserve(possible_inlier_number);
    for (size_t point_idx = 0; point_idx < possible_inlier_number; ++point_idx) {
        // Calculate the weight of the current point
        double weight = 0.0;
        for (size_t partition_idx = 0; partition_idx < partition_number; ++partition_idx)
            weight += point_weights_par[partition_idx][point_idx];

        // If the weight is approx. zero, continue.
        if (weight < std::numeric_limits<double>::epsilon())
            continue;

        // Store the index and weight of the current point
        sigma_inliers.emplace_back(all_residuals[point_idx].second);
        final_weights.emplace_back(weight);
    }

    // If there are fewer inliers than the size of the minimal sample interupt the procedure
    if (sigma_inliers.size() < sample_size)
        return false;

    // Estimate the model parameters using weighted least-squares fitting
    _model->Fit(sigma_inliers, &sigma_models, &(final_weights)[0]);

    bool is_model_updated = false;

    if (sigma_models.size() == 1) // && // If only a single model is estimated
//        estimator_.isValidModel(sigma_models.back(), // TODO
//                                points_,
//                                sigma_inliers,
//                                &(sigma_inliers)[0],
//                                interrupting_threshold,
//                                is_model_updated)) // and it is valid
    {
        // Return the refined model
        refined_model_ = sigma_models.back();

        // Calculate the score of the model and the implied iteration number
        double marginalized_iteration_number;

        getModelQuality(refined_model_, // The estimated model
                        marginalized_iteration_number, // The marginalized inlier ratio
                        score_.score); // The marginalized score

        if (marginalized_iteration_number < 0 || std::isnan(marginalized_iteration_number))
            last_iteration_number = std::numeric_limits<int>::max();
        else
            last_iteration_number = static_cast<int>(round(marginalized_iteration_number));
        inliersIdxsSaved = sigma_inliers;
        weightsSaved = final_weights;
        return true;
    }
    return false;
}

bool MAGSAC::sigmaConsensusPlusPlus(const ModelEstimator::Model &model_,
                                    ModelEstimator::Model &refined_model_,
                                    ModelScore &score_,
                                    const ModelScore &best_score_,
                                    int &last_iteration_number,
                                    std::vector<int> &inliersIdxsSaved,
                                    std::vector<double> &weightsSaved) const {
    // The degrees of freedom of the data from which the model is estimated.
    // E.g., for models coming from point correspondences (x1,y1,x2,y2), it is 4.
    const size_t degrees_of_freedom = _model->dataDegreeOfFreedom();
    // A 0.99 quantile of the Chi^2-distribution to convert sigma values to residuals
    const double k = getSigmaQuantile();
    // A multiplier to convert residual values to sigmas
    const double threshold_to_sigma_multiplier = 1.0 / k;
    // Calculating k^2 / 2 which will be used for the estimation and,
    // due to being constant, it is better to calculate it a priori.
    const double squared_k_per_2 = k * k / 2.0;
    // Calculating (DoF - 1) / 2 which will be used for the estimation and,
    // due to being constant, it is better to calculate it a priori.
    const double dof_minus_one_per_two = (degrees_of_freedom - 1.0) / 2.0;
    const double C = _model->C();
    // The size of a minimal sample used for the estimation
    const size_t sample_size = _model->SizeSample();
    // Calculating 2^(DoF - 1) which will be used for the estimation and,
    // due to being constant, it is better to calculate it a priori.
    static const double two_ad_dof = std::pow(2.0, dof_minus_one_per_two);
    // Calculating C * 2^(DoF - 1) which will be used for the estimation and,
    // due to being constant, it is better to calculate it a priori.
    static const double C_times_two_ad_dof = C * two_ad_dof;
    // Calculating the gamma value of (DoF - 1) / 2 which will be used for the estimation and,
    // due to being constant, it is better to calculate it a priori.
    static const double gamma_value = tgamma(dof_minus_one_per_two);
    // Calculating the upper incomplete gamma value of (DoF - 1) / 2 with k^2 / 2.
    const double gamma_k = _model->upperGamma();
    // Calculating the lower incomplete gamma value of (DoF - 1) / 2 which will be used for the estimation and,
    // due to being constant, it is better to calculate it a priori.
    static const double gamma_difference = gamma_value - gamma_k;
    // The number of points provided
    const int point_number = _model->NbData();
    // The manually set maximum inlier-outlier threshold
    double current_maximum_sigma = this->maximum_threshold;
    // Calculating the pairs of (residual, point index).
    std::vector<std::pair<double, size_t> > residuals;
    // Occupy the maximum required memory to avoid doing it later.
    residuals.reserve(point_number);

    // If it is not the first run, consider the previous best and interrupt the validation when there is no chance of being better
    if (best_score_.inlier_number > 0) {
        // Number of points close to the previous so-far-the-best model.
        // This model should have more inliers.
        int points_remaining = best_score_.inlier_number;

        // Collect the points which are closer than the threshold which the maximum sigma implies
        for (int point_idx = 0; point_idx < point_number; ++point_idx) {
            // Calculate the residual of the current point
            const double residual = std::sqrt(_model->Error(model_, point_idx));
            if (current_maximum_sigma > residual) {
                // Store the residual of the current point and its index
                residuals.emplace_back(std::make_pair(residual, point_idx));
                // all_residuals.emplace_back(std::make_pair(residual * threshold_to_sigma_multiplier, point_idx));

                // Count points which are closer than a reference threshold to speed up the procedure
                if (residual < interrupting_threshold)
                    --points_remaining;
            }

            // Interrupt if there is no chance of being better
            // TODO: replace this part by SPRT test
            if (point_number - point_idx < points_remaining)
                return false;
        }

        // Store the number of really close inliers just to speed up the procedure
        // by interrupting the next verifications.
        score_.inlier_number = best_score_.inlier_number - points_remaining;
    } else {
        // The number of really close points
        size_t points_close = 0;

        // Collect the points which are closer than the threshold which the maximum sigma implies
        for (size_t point_idx = 0; point_idx < point_number; ++point_idx) {
            // Calculate the residual of the current point
            const double residual = sqrt(_model->Error(model_, point_idx));
            if (current_maximum_sigma > residual) {
                // Store the residual of the current point and its index
                residuals.emplace_back(std::make_pair(residual, point_idx));

                // Count points which are closer than a reference threshold to speed up the procedure
                if (residual < interrupting_threshold)
                    ++points_close;
            }
        }

        // Store the number of really close inliers just to speed up the procedure
        // by interrupting the next verifications.
        score_.inlier_number = points_close;
    }

    // Models fit by weighted least-squares fitting
    std::vector<ModelEstimator::Model> sigma_models;
    // Points used in the weighted least-squares fitting
    std::vector<int> sigma_inliers;
    // Weights used in the the weighted least-squares fitting
    std::vector<double> sigma_weights;
    // Number of points considered in the fitting
    const size_t possible_inlier_number = residuals.size();
    // Occupy the memory to avoid doing it inside the calculation possibly multiple times
    sigma_inliers.reserve(possible_inlier_number);
    // Occupy the memory to avoid doing it inside the calculation possibly multiple times
    sigma_weights.reserve(possible_inlier_number);

    // Calculate 2 * \sigma_{max}^2 a priori
    const double squared_sigma_max_2 = current_maximum_sigma * current_maximum_sigma * 2.0;
    // Divide C * 2^(DoF - 1) by \sigma_{max} a priori
    const double one_over_sigma = C_times_two_ad_dof / current_maximum_sigma;
    // Calculate the weight of a point with 0 residual (i.e., fitting perfectly) a priori
    const double weight_zero = one_over_sigma * gamma_difference;

    // Initialize the polished model with the initial one
    ModelEstimator::Model polished_model = model_;
    // A flag to determine if the initial model has been updated
    bool updated = false;

    // Do the iteratively re-weighted least squares fitting
    for (size_t iterations = 0; iterations < number_of_irwls_iters; ++iterations) {
        // If the current iteration is not the first, the set of possibly inliers
        // (i.e., points closer than the maximum threshold) have to be recalculated.
        if (iterations > 0) {
            // The number of points close to the model
            size_t points_close = 0;
            // Remove everything from the residual vector
            residuals.clear();

            // Collect the points which are closer than the maximum threshold
            for (size_t point_idx = 0; point_idx < point_number; ++point_idx) {
                // Calculate the residual of the current point
                const double residual = sqrt(_model->Error(polished_model, point_idx));
                if (current_maximum_sigma > residual) {
                    // Store the residual of the current point and its index
                    residuals.emplace_back(std::make_pair(residual, point_idx));

                    // Count points which are closer than a reference threshold to speed up the procedure
                    if (residual < interrupting_threshold)
                        ++points_close;
                }
            }

            // Store the number of really close inliers just to speed up the procedure
            // by interrupting the next verifications.
            score_.inlier_number = points_close;

            // Number of points closer than the threshold
            const size_t possible_inlier_number = residuals.size();

            // Clear the inliers and weights
            sigma_inliers.clear();
            sigma_weights.clear();

            // Occupy the memory for the inliers and weights
            sigma_inliers.reserve(possible_inlier_number);
            sigma_weights.reserve(possible_inlier_number);
        }

        // Calculate the weight of each point
        for (const auto &[residual, idx] : residuals) {
            // The weight
            double weight = 0.0;
            // If the residual is ~0, the point fits perfectly and it is handled differently
            if (residual < std::numeric_limits<double>::epsilon())
                weight = weight_zero;
            else {
                // Calculate the squared residual
                const double squared_residual = residual * residual;
                // Get the position of the gamma value in the lookup table
                size_t x = round(precision_of_stored_gammas * squared_residual / squared_sigma_max_2);
                // Put the index of the point into the vector of points used for the least squares fitting
                sigma_inliers.emplace_back(idx);

                // If the sought gamma value is not stored in the lookup, return the closest element
                if (stored_gamma_number < x)
                    x = stored_gamma_number;

                // Calculate the weight of the point
                weight = one_over_sigma * (stored_gamma_values[x] - gamma_k);
            }

            // Store the weight of the point
            sigma_weights.emplace_back(weight);
        }

        // If there are fewer than the minimum point close to the model,
        // terminate.
        if (sigma_inliers.size() < sample_size)
            return false;

        // Estimate the model parameters using weighted least-squares fitting
        _model->Fit(sigma_inliers, &sigma_models, &(sigma_weights)[0]);

        // Estimate the model parameters using weighted least-squares fitting
        if (sigma_models.empty()) {
            if (iterations == 0) {
                return false;
            }
            // Otherwise, if the iteration was successfull at least one,
            // simply break it.
            break;
        }

        // Update the model parameters
        polished_model = sigma_models[0];
        // Clear the vector of models and keep only the best
        sigma_models.clear();
        // The model has been updated
        updated = true;
    }

    bool is_model_updated = false;

    if (updated) {
        // Return the refined model
        refined_model_ = polished_model;

        // Calculate the score of the model and the implied iteration number
        getModelQualityPlusPlus(refined_model_, // The estimated model
                                score_.score, // The marginalized score
                                best_score_.score); // The score of the previous so-far-the-best model

        // Update the iteration number
        last_iteration_number = _log_confidence /
                log(1.0 - std::pow(static_cast<double>(score_.inlier_number) / point_number, sample_size));
        inliersIdxsSaved = sigma_inliers;
        weightsSaved = sigma_weights;
        return true;
    }
    return false;
}

void MAGSAC::getModelQuality(const ModelEstimator::Model &model_, // The model parameter
                             double &marginalized_iteration_number_, // The marginalized iteration number to be calculated
                             double &score_) // The score to be calculated
const {
    // Set up the parameters
    const size_t sample_size = _model->SizeSample();
    const size_t point_number = _model->NbData();

    // Getting the inliers
    std::vector<std::pair<double, size_t>> all_residuals;
    all_residuals.reserve(point_number);

    double max_distance = 0;
    for (size_t point_idx = 0; point_idx < point_number; ++point_idx) {
        // Calculate the residual of the current point
        const double residual = std::sqrt(_model->Error(model_, point_idx)); // TODO
        // If the residual is smaller than the maximum threshold, add it to the set of possible inliers
        if (maximum_threshold > residual) {
            max_distance = std::max(max_distance, residual);
            all_residuals.emplace_back(std::make_pair(residual, point_idx));
        }
    }

    // Set the maximum distance to be slightly bigger than that of the farthest possible inlier
    max_distance = max_distance +
                   std::numeric_limits<double>::epsilon();

    // Number of possible inliers
    const size_t possible_inlier_number = all_residuals.size();

    // The extent of a partition
    const double threshold_step = max_distance / partition_number;

    // The maximum threshold considered in each partition
    std::vector<double> thresholds(partition_number);
    std::vector<double> thresholds_squared(partition_number);
    std::vector<double> thresholds_2_squared(partition_number);

    // Calculating the thresholds for each partition
    for (size_t i = 0; i < partition_number; ++i) {
        thresholds[i] = (i + 1) * threshold_step;
        thresholds_squared[i] = thresholds[i] * thresholds[i];
        thresholds_2_squared[i] = 2 * thresholds_squared[i];
    }

    double residual_i, // Residual of the i-th point
    residual_i_squared, // Squared residual of the i-th poin
    probability_i; // Probability of the i-th point given the model

    std::vector<double> inliers(partition_number, 0), // RANSAC score for each partition
    probabilities(partition_number, 1); // Probabilities for each partition
    for (size_t point_idx = 0; point_idx < possible_inlier_number; ++point_idx) {
        residual_i = all_residuals[point_idx].first;
        residual_i_squared = residual_i * residual_i;

        for (size_t i = 0; i < partition_number; ++i) {
            if (residual_i < thresholds[i]) {
                probability_i = 1.0 - residual_i_squared / thresholds_squared[i];
                ++inliers[i];
                probabilities[i] += probability_i;
            }
        }
    }

    score_ = 0;
    marginalized_iteration_number_ = 0.0;
    for (size_t i = 0; i < partition_number; ++i) {
        score_ += probabilities[i];
        marginalized_iteration_number_ +=
                _log_confidence / log(1.0 - std::pow(inliers[i] / point_number, sample_size));
    }
    marginalized_iteration_number_ = marginalized_iteration_number_ / partition_number;
}

void MAGSAC::getModelQualityPlusPlus(const ModelEstimator::Model &model_, // The model parameter
                                     double &score_,
                                     const double &previous_best_score_)// The score of the previous so-far-the-best model
const {
    // The degrees of freedom of the data from which the model is estimated.
    // E.g., for models coming from point correspondences (x1,y1,x2,y2), it is 4.
    const size_t degrees_of_freedom = _model->dataDegreeOfFreedom();
    // A 0.99 quantile of the Chi^2-distribution to convert sigma values to residuals
    const double k = getSigmaQuantile();
    // A multiplier to convert residual values to sigmas
    const double threshold_to_sigma_multiplier = 1.0 / k;
    // Calculating k^2 / 2 which will be used for the estimation and,
    // due to being constant, it is better to calculate it a priori.
    const double squared_k_per_2 = k * k / 2.0;
    // Calculating (DoF - 1) / 2 which will be used for the estimation and,
    // due to being constant, it is better to calculate it a priori.
    const double dof_minus_one_per_two = (degrees_of_freedom - 1.0) / 2.0;
    // Calculating (DoF + 1) / 2 which will be used for the estimation and,
    // due to being constant, it is better to calculate it a priori.
    const double dof_plus_one_per_two = (degrees_of_freedom + 1.0) / 2.0;
    const double C = _model->C();
    // Calculating 2^(DoF - 1) which will be used for the estimation and,
    // due to being constant, it is better to calculate it a priori.
    static const double two_ad_dof_minus_one = std::pow(2.0, dof_minus_one_per_two);
    // Calculating 2^(DoF + 1) which will be used for the estimation and,
    // due to being constant, it is better to calculate it a priori.
    static const double two_ad_dof_plus_one = std::pow(2.0, dof_plus_one_per_two);
    // Calculate the gamma value of k
    const double gamma_value_of_k = _model->upperGamma();
    // Calculate the lower incomplete gamma value of k
    const double lower_gamma_value_of_k = _model->lowerGamma();
    // The number of points provided
    const int point_number = _model->NbData();
    // The previous best loss
    const double previous_best_loss = 1.0 / previous_best_score_;
    // Convert the maximum threshold to a sigma value
    const double maximum_sigma = threshold_to_sigma_multiplier * maximum_threshold;
    // Calculate the squared maximum sigma
    const double maximum_sigma_2 = maximum_sigma * maximum_sigma;
    // Calculate \sigma_{max}^2 / 2
    const double maximum_sigma_2_per_2 = maximum_sigma_2 / 2.0;
    // Calculate 2 * \sigma_{max}^2
    const double maximum_sigma_2_times_2 = maximum_sigma_2 * 2.0;
    // Calculate the loss implied by an outlier
    const double outlier_loss = maximum_sigma * two_ad_dof_minus_one * lower_gamma_value_of_k;
    // Calculating 2^(DoF + 1) / \sigma_{max} which will be used for the estimation and,
    // due to being constant, it is better to calculate it a priori.
    const double two_ad_dof_plus_one_per_maximum_sigma = two_ad_dof_plus_one / maximum_sigma;
    // The loss which a point implies
    double loss = 0.0,
    // The total loss regarding the current model
    total_loss = 0.0;

    // Iterate through all points to calculate the implied loss
    for (size_t point_idx = 0; point_idx < point_number; ++point_idx) {
        // Calculate the residual of the current point
        const double residual = std::sqrt(_model->Error(model_, point_idx));

        // If the residual is smaller than the maximum threshold, consider it outlier
        // and add the loss implied to the total loss.
        if (maximum_threshold < residual)
            loss = outlier_loss;
        else // Otherwise, consider the point inlier, and calculate the implied loss
        {
            // Calculate the squared residual
            const double squared_residual = residual * residual;
            // Divide the residual by the 2 * \sigma^2
            const double squared_residual_per_sigma = squared_residual / maximum_sigma_2_times_2;
            // Get the position of the gamma value in the lookup table
            size_t x = round(precision_of_stored_incomplete_gammas * squared_residual_per_sigma);
            // If the sought gamma value is not stored in the lookup, return the closest element
            if (stored_incomplete_gamma_number < x)
                x = stored_incomplete_gamma_number;

            // Calculate the loss implied by the current point
            loss = maximum_sigma_2_per_2 * stored_lower_incomplete_gamma_values[x] +
                   squared_residual / 4.0 * (stored_complete_gamma_values[x] -
                                             gamma_value_of_k);
            loss = loss * two_ad_dof_plus_one_per_maximum_sigma;
        }

        // Update the total loss
        total_loss += loss;

        // Break the validation if there is no chance of being better than the previous
        // so-far-the-best model.
        if (previous_best_loss < total_loss)
            break;
    }
    // Calculate the score of the model from the total loss
    score_ = 1.0 / total_loss;
}

/// Verify the runOutput metric is good enough.
bool MAGSAC::satisfyingRun(double runOutput) const {
    return runOutput > 0;
}

} // namespace orsa
