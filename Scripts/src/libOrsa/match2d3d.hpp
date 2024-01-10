//
// Created by riuclement on 9/2/21.
//

#ifndef MMM_ORSA_MATCH2D3D_HPP
#define MMM_ORSA_MATCH2D3D_HPP


#include "libNumerics/matrix.h"
#include "eigen/Eigen/Eigen"
#include <vector>
#include <fstream>
#include <sstream>
#include <string>

/// Save matching position between two points.
struct Match2D3D {
    Match2D3D() {}

    Match2D3D(float x2D, float y2D, float x3D, float y3D, float z3D) {
        p2D(0) = x2D;
        p2D(1) = y2D;
        p3D(0) = x3D;
        p3D(1) = y3D;
        p3D(2) = z3D;
    }

    Match2D3D(Eigen::Vector2d p2D, Eigen::Vector3d p3D)
            : p2D(p2D), p3D(p3D) {}

    Eigen::Vector2d p2D;
    Eigen::Vector3d p3D;

    /**
    * Transform into matrix where each column is (x1 y1 x2 y2 z2)^T.
    * \param m The matches that we transform.
    * \return Matrix 5xn where n is the size of vector \a m.
    */
    static libNumerics::matrix<double> toMat(const std::vector<Match2D3D> &m) {
        libNumerics::matrix<double> M(5, static_cast<int>(m.size()));
        std::vector<Match2D3D>::const_iterator it = m.begin();
        for (int j = 0; it != m.end(); ++it, ++j) {
            M(0, j) = it->p2D(0);
            M(1, j) = it->p2D(1);
            M(2, j) = it->p3D(0);
            M(3, j) = it->p3D(1);
            M(4, j) = it->p3D(2);
        }
        return M;
    }

    /**
    * Load the corresponding matches from file.
    * \param nameFile   The file where matches were saved.
    * \param vec_match  The loaded corresponding points.
    * \return bool      True if everything was ok, otherwise false.
    */
    static bool loadMatch2D3D(const char *nameFile, std::vector<Match2D3D> &vec_match) {
        vec_match.clear();
        std::ifstream f(nameFile);
        while (f.good()) {
            std::string str;
            std::getline(f, str);
            if (f.good()) {
                std::istringstream s(str);
                Match2D3D m;
                s >> m;
                if (!s.fail())
                    vec_match.push_back(m);
            }
        }
        return !vec_match.empty();
    }

    /**
    * Save the corresponding matches to file.
    * \param nameFile   The file where matches will be saved.
    * \param vec_match  The matches that we want to save.
    * \return bool True if everything was ok, otherwise false.
    */
    static bool saveMatch2D3D(const char *nameFile, const std::vector<Match2D3D> &vec_match) {
        std::ofstream f(nameFile);
        if (f.is_open()) {
            std::vector<Match2D3D>::const_iterator it = vec_match.begin();
            for (; it != vec_match.end(); ++it)
                f << *it;
        }
        return f.is_open();
    }

    /// Lexicographical ordering of matches. Used to remove duplicates.
    friend bool operator<(const Match2D3D &m1, const Match2D3D &m2) {
        if (m1.p2D(0) < m2.p2D(0)) return true;
        if (m1.p2D(0) > m2.p2D(0)) return false;

        if (m1.p2D(1) < m2.p2D(1)) return true;
        if (m1.p2D(1) > m2.p2D(1)) return false;

        if (m1.p3D(0) < m2.p3D(0)) return true;
        if (m1.p3D(0) > m2.p3D(0)) return false;

        if (m1.p3D(1) < m2.p3D(1)) return true;
        if (m1.p3D(1) > m2.p3D(1)) return false;

        return (m1.p3D(2) < m2.p3D(2));
    }

    /// Comparison Operator
    friend bool operator==(const Match2D3D &m1, const Match2D3D &m2) {
        return (m1.p2D.isApprox(m2.p2D) && m1.p3D.isApprox(m2.p3D));
    }

    friend std::ostream &operator<<(std::ostream &os, const Match2D3D &m) {
        return os << m.p2D(0) << " "
                  << m.p2D(1) << " "
                  << m.p3D(0) << " "
                  << m.p3D(1) << " "
                  << m.p3D(2) << std::endl;
    }

    friend std::istream &operator>>(std::istream &in, Match2D3D &m) {
        return in >> m.p2D(0) >> m.p2D(1) >> m.p3D(0) >> m.p3D(1) >> m.p3D(2);
    }
};

#endif //MMM_ORSA_MATCH2D3D_HPP
