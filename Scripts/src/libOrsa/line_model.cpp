/**
* @file line_model.cpp
* @brief Compute 2D line from plane points
* @author Clement Riu
*
* Copyright (c) 2020 Clement Riu
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

#include <iostream>
#include <random>


#include "line_model.hpp"

namespace orsa {

/// Constructor.
LineModel::LineModel(const Mat &data, bool symError, int nPoints, double trueInlierRatio, double trueSigma,
                     const int width, const int height, const unsigned int randomSeed)
: ModelEstimator(data, symError), _nPoints(nPoints),
  _trueInlierRatio(trueInlierRatio),
  _trueSigma(trueSigma), _width(width), _height(height), _area(width * height),
  _diameter(std::sqrt(std::pow(width, 2) + std::pow(height, 2))) {
    _generate2DArtificialData(randomSeed);
}

// Generation of random data according to the model:
// - _outliers from uniform background distribution.
// - _inliers from a random linear model.
void LineModel::_generate2DArtificialData(const unsigned int randomSeed) {
    //// Generation of the data:
    int nInliers = (int) (_trueInlierRatio * _nPoints);
    int nOutliers = _nPoints - nInliers;

    Mat _outliers(nOutliers, 2);
    Mat _inliers(nInliers, 2);


    // Random generation of outliers: elements are uniformly sampled in the considered range:
    for (int indexOutlier = 0; indexOutlier < nOutliers; indexOutlier++) {
        _outliers(indexOutlier, 0) = std::rand() % _width;
        _outliers(indexOutlier, 1) = std::rand() % _height;
    }

    // Generation of the model: the points are uniformly drawn from the
    // considered range and the associated line is computed:
    int x1 = std::rand() % _width;
    int y1 = std::rand() % _height;
    int x2, y2;
    do {
        x2 = std::rand() % _width;
        y2 = std::rand() % _height;
    } while (abs(x1 - x2) < 1);

    std::vector<Model> trueModels;
    Fit(x1, x2, y1, y2, &trueModels);
    Model trueModel = trueModels[0];

    std::cout << "The true model is: ";
    std::cout << trueModel << std::endl;

    // Random inliers: x-coord is drawn uniformly in range of "printable" xs.
    double xMin = 0.0;
    double xMax = (double) _width;
    if (trueModel(0, 0) > 1.0e-10) { // TODO: use global parameter for definition of "too small".
        xMin = std::ceil(std::max(-trueModel(0,1)/trueModel(0,0), 0.0));
        xMax = std::floor(std::min((_height-trueModel(0,1))/trueModel(0,0),
                                   (double)_width));
    }
    if (trueModel(0, 0) < -1.0e-10) { // TODO: use global parameter for definition of "too small".
        xMin = std::floor(std::min(-trueModel(0,1) / trueModel(0,0),
                                   (double)_width));
        xMax = std::ceil(std::max((_height-trueModel(0,1))/trueModel(0,0),
                                  0.0));
    }

    // Initialisation of the random generators: uniform to draw the x-coordinate and gaussian for the noise:
    std::default_random_engine generator;
    generator.seed(randomSeed);
    std::uniform_real_distribution<double> distributionUniform(xMin, xMax);
    std::normal_distribution<double> distributionGaussian(0.0, _trueSigma);

    // Generation of the inliers:
    for (int indexInlier = 0; indexInlier < nInliers; indexInlier++) {
        double x = distributionUniform(generator);
        double y = trueModel(0, 0) * x + trueModel(0, 1);

        // Some noise is added.
        double xNoisy = x + distributionGaussian(generator);
        double yNoisy = y + distributionGaussian(generator);

        if (xNoisy < 0 || xNoisy > _width)
            xNoisy = x;
        if (yNoisy < 0 || yNoisy > _height)
            yNoisy = y;

        _inliers(indexInlier, 0) = (int) std::round(xNoisy);
        _inliers(indexInlier, 1) = (int) std::round(yNoisy);
    }

    //// Merging outliers and inliers randomly:
    int emptyedOutliers = 0;
    int emptyedInliers = 0;
    int iter = 0;
    while (iter < _nPoints) {
        int currentIndex = emptyedOutliers + emptyedInliers;
        if (emptyedOutliers < nOutliers && emptyedInliers < nInliers) {
            int chosenDataSet = std::rand() % 2;
            if (chosenDataSet == 0) {
                this->data_(0, currentIndex) = _outliers(emptyedOutliers, 0);
                this->data_(1, currentIndex) = _outliers(emptyedOutliers, 1);
                _dataIsInlier.push_back(false);
                emptyedOutliers++;
            } else {
                this->data_(0, currentIndex) = _inliers(emptyedInliers, 0);
                this->data_(1, currentIndex) = _inliers(emptyedInliers, 1);
                _dataIsInlier.push_back(true);
                emptyedInliers++;
            }
        } else {
            if (emptyedOutliers >= nOutliers && emptyedInliers < nInliers) {
                this->data_(0, currentIndex) = _inliers(emptyedInliers, 0);
                this->data_(1, currentIndex) = _inliers(emptyedInliers, 1);
                _dataIsInlier.push_back(true);
                emptyedInliers++;
            } else {
                if (emptyedInliers >= nInliers && emptyedOutliers < nOutliers) {
                    this->data_(0, currentIndex) = _outliers(emptyedOutliers, 0);
                    this->data_(1, currentIndex) = _outliers(emptyedOutliers, 1);
                    emptyedOutliers++;
                    _dataIsInlier.push_back(false);

                } else {
                    throw std::runtime_error(
                            "PROBLEM DURING GENERATION OF ARTIFICIAL 2D DATASET:\nNot enough points generated");
                }
            }
        }
        iter++;
    }
}

/// Computes the models associated to indexed sample.
/// \param indices Indices of points to consider for model estimation.
/// \param models  Estimated model(s) from sampled point.
void LineModel::Fit(const std::vector<int> &indices, std::vector<Model> *lines) const {
    if (indices.size() < 2) {
        return;
    }
    if (indices.size() == 2) { // Could be removed as the least-squared solution is the same.
        int x1 = (int) this->data_(0, indices[0]);
        int x2 = (int) this->data_(0, indices[1]);
        int y1 = (int) this->data_(1, indices[0]);
        int y2 = (int) this->data_(1, indices[1]);
        Fit(x1, x2, y1, y2, lines);
    } else {
        double sumXi = 0, sumXiSquared = 0, sumXiYi = 0, sumYi = 0;
        int nInliers = (int) indices.size();
        std::vector<int>::const_iterator it = indices.begin();
        for (; it != indices.end(); ++it) {
            //            for (int i = 0; i < nInliers; i++) {
            int index = *it;
            double xi = this->data_(0, index);
            double yi = this->data_(1, index);
            sumXi += xi;
            sumXiSquared += xi * xi;
            sumXiYi += xi * yi;
            sumYi += yi;
        }
        Model model(1, 2);
        double den = nInliers * sumXiSquared - sumXi * sumXi;
        model(0, 0) = (nInliers * sumXiYi - sumXi * sumYi) / den;
        model(0, 1) = (-sumXi * sumXiYi + sumYi * sumXiSquared) / den;
        lines->push_back(model);
    }
}

/// Computes the models associated to indexed sample.
/// \param x1, x2, y1, y2 coordinates of the 2D points to consider for model estimation.
/// \param models  Estimated model(s) from sampled point.
void LineModel::Fit(const int x1, const int x2, const int y1, const int y2, std::vector<Model> *lines) const {
    Model model(1, 2);
    if (abs(x1 - x2) < 1) {
        model(0, 0) = -1;
        model(0, 1) = -1;
    } else {
        model(0, 0) = (double) (y2 - y1) / (double) (x2 - x1);
        model(0, 1) = ((double) y1 - model(0, 0) * x1);
    }
    lines->push_back(model);
}

/// Computes the square error of a correspondence wrt \a model.
/// \param model The model to evaluate.
/// \param index The point index stored in the Kernel.
/// \param[out] Not used.
double LineModel::Error(const Model &line, int index, int *side) const {
    int x = (int) this->data_(0,index);
    int y = (int) this->data_(1,index);
    double num = (line(0,0)*x+line(0,1) - y) * (line(0,0)*x+line(0,1) - y);
    double den = 1 + line(0,0)*line(0,0);
    if (side) *side = 0;
    return num / den;
}

/// Computes the square error of a correspondence wrt \a model.
/// \param model The model to evaluate.
/// \param index The point index stored in the Kernel.
/// \param[out] Not used.
double LineModel::Error(const Model &line, Mat testData, int *side) const {
    int x = (int) testData(0,0);
    int y = (int) testData(1,0);
    double num = (line(0,0)*x+line(0,1) - y) * (line(0,0)*x+line(0,1) - y);
    double den = 1 + line(0,0)*line(0, 0);
    if (side) *side = 0;
    return num / den;
}

/// Computes the area of the inlier region according to the AC-RANSAC
/// paper for a given error margin. The inlier region is over-estimated to
/// ease computation. This is a point-to-line error.
/// \param sigma The error margin.
double LineModel::pSigma(const double sigma) const {
    return 2 * sigma * _diameter / (double) _area;
}

} // namespace orsa
