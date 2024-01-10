/**
* @file line_model.hpp
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

#ifndef MMM_ORSA_LINE_MODEL_HPP
#define MMM_ORSA_LINE_MODEL_HPP

#include "model_estimator.hpp"

namespace orsa {

    class LineModel : public ModelEstimator {
    public:
        LineModel(const Mat &data, bool symError, int nPoints, double trueInlierRatio, double trueSigma, int width,
                  int height, unsigned int randomSeed);

        /// 2 points are required to compute a line.
        int SizeSample() const { return 2; }

        /// Up to 1 line is computed from a sample of 2 points.
        int NbModels() const { return 1; }

        /// Degree of freedom of the model
        int nDegreeOfFreedom() const { return 2; };
        /// Degree of freedom of the data (for magsac++)
        int dataDegreeOfFreedom() const {return 2;};
        /// C(n) = 1/(n/2 * Gamma(n/2)) where n is DoF (for magsac++)
        double C() const {return 1;};
        /// UG((DoF - 1) /2, k**2/2) (k is sigma quantile) (for magsac++)
        double upperGamma() const {return 0.0004832383458165898;};
        /// LG((DoF + 1) /2, k**2/2) (k is sigma quantile) (for magsac++)
        double lowerGamma() const {return 0.8825696646186669357;};

        /// Distance used to distinguish inlier/outlier is to a line
        virtual bool DistToPoint() const { return false; }

        /// Computes the area of the inlier region according to the AC-RANSAC
        /// paper for a given error margin. The inlier region is over-estimated to
        /// ease computation. This is a point-to-line error.
        /// \param sigma The error margin.
        double pSigma(double sigma) const;

        /// Computes the models associated to indexed sample.
        /// \param indices Indices of points to consider for model estimation.
        /// \param models  Estimated model(s) from sampled point.
        void Fit(const std::vector<int> &indices, std::vector<Model> *lines) const;

        /// Computes the models associated to indexed sample.
        /// \param x1, x2, y1, y2 coordinates of the 2D points to consider for model estimation.
        /// \param models  Estimated model(s) from sampled point.
        void Fit(int x1, int x2, int y1, int y2, std::vector<Model> *lines) const;

        /// Distance to line.
        double Error(const Model &line, int index, int *side = 0) const;

        double Error(const Model &line, Mat testData, int *side = 0) const;

        /// Returns the fundamental matrix. For Lines the output is the input.
        /// Necessary for essential matrix.
        inline Model toPixelSpace(const Model &L) const { return L; }

    private:
        //// Private parameters:
        // Number of points to generate.
        int _nPoints;

        // Parameters of the generation.
        double _trueInlierRatio;
        double _trueSigma;

        // Size of the range.
        int _width, _height;
        // Area of the region.
        const int _area;
        // Diameter of the region.
        const double _diameter;


        // Wether a point is an inlier or outlier.
        std::vector<bool> _dataIsInlier;

        //// Private functions:
        // Generation of random data according to the model:
        // - _outliers from uniform background distribution.
        // - _inliers from a random linear model.
        void _generate2DArtificialData(unsigned int randomSeed);
    };


} // namespace orsa

#endif //MMM_ORSA_LINE_MODEL_HPP
