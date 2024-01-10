/**
 * @file siftMatch.hpp
 * @brief SIFT extraction and matching
 * @author Lionel Moisan, Pascal Monasse, Pierre Moulon
 * 
 * Copyright (c) 2011 Lionel Moisan, Pascal Monasse, Pierre Moulon
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

#ifndef SIFT_MATCH_HPP
#define SIFT_MATCH_HPP

#include "libImage/image.hpp"

#include "third_party/sift_anatomy/lib_sift_anatomy.h"
#include "third_party/sift_anatomy/lib_matching.h"
#include "libOrsa/match.hpp"

#include <algorithm>
#include <istream>
#include <ostream>
#include <vector>


namespace utility {

/// Rectangle in image
class Geometry {
public:
    int x0, x1, y0, y1;

    bool inside(double x, double y) const {
        int ix = static_cast<int>(x), iy = static_cast<int>(y);
        return (x0 <= ix && ix < x1 && y0 <= iy && iy < y1);
    }

    bool inside(keypoint k) const {
        return inside(k.y, k.x); // Weird SiftAnatomy's coordinate system...
    }
};

/// Output @geo.
std::ostream &operator<<(std::ostream &str, Geometry &geo);

/// Input @geo. Format: wxh+x0+y0, eg 100x100+0+0
std::istream &operator>>(std::istream &str, Geometry &geo);

/// SIFT matches
void SIFT(const Image<unsigned char> &im1,
          const Image<unsigned char> &im2,
          std::vector<Match> &vec_matchings,
          float fMatchRatio = 0.6f, Geometry *rect = 0);

/// Remove multiple "same position" matches
template <typename M>
void rm_duplicates(std::vector<M> &m) {
    std::sort(m.begin(), m.end());
    typename std::vector<M>::iterator end = std::unique(m.begin(), m.end());
    if (end != m.end()) {
        std::cout << "Remove " << std::distance(end, m.end())
                  << "/" << m.size() << " duplicate matches, "
                  << "keeping " << std::distance(m.begin(), end) << std::endl;
        m.erase(end, m.end());
    }
}

} // namespace utility

#endif
