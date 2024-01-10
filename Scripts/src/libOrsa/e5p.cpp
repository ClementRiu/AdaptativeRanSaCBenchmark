#include "e5p.h"

/// Script extracted from MAGSAC implementation.
// Multiply two degree one polynomials of variables x, y, z.
// E.g. p1 = a[0]x + a[1]y + a[2]z + a[3]
// Output order: x^2 xy y^2 xz yz z^2 x y z 1 (GrevLex)
Eigen::Matrix<double, 1, 10> multiplyDegOnePoly(
        const Eigen::RowVector4d &a,
        const Eigen::RowVector4d &b) {
    Eigen::Matrix<double, 1, 10> output;
    // x^2
    output(0) = a(0) * b(0);
    // xy
    output(1) = a(0) * b(1) + a(1) * b(0);
    // y^2
    output(2) = a(1) * b(1);
    // xz
    output(3) = a(0) * b(2) + a(2) * b(0);
    // yz
    output(4) = a(1) * b(2) + a(2) * b(1);
    // z^2
    output(5) = a(2) * b(2);
    // x
    output(6) = a(0) * b(3) + a(3) * b(0);
    // y
    output(7) = a(1) * b(3) + a(3) * b(1);
    // z
    output(8) = a(2) * b(3) + a(3) * b(2);
    // 1
    output(9) = a(3) * b(3);
    return output;
}

/// Script extracted from MAGSAC implementation.
// Multiply a 2 deg poly (in x, y, z) and a one deg poly in GrevLex order.
// x^3 x^2y xy^2 y^3 x^2z xyz y^2z xz^2 yz^2 z^3 x^2 xy y^2 xz yz z^2 x y z 1
Eigen::Matrix<double, 1, 20> multiplyDegTwoDegOnePoly(
        const Eigen::Matrix<double, 1, 10> &a,
        const Eigen::RowVector4d &b) {
    Eigen::Matrix<double, 1, 20> output;
    // x^3
    output(0) = a(0) * b(0);
    // x^2y
    output(1) = a(0) * b(1) + a(1) * b(0);
    // xy^2
    output(2) = a(1) * b(1) + a(2) * b(0);
    // y^3
    output(3) = a(2) * b(1);
    // x^2z
    output(4) = a(0) * b(2) + a(3) * b(0);
    // xyz
    output(5) = a(1) * b(2) + a(3) * b(1) + a(4) * b(0);
    // y^2z
    output(6) = a(2) * b(2) + a(4) * b(1);
    // xz^2
    output(7) = a(3) * b(2) + a(5) * b(0);
    // yz^2
    output(8) = a(4) * b(2) + a(5) * b(1);
    // z^3
    output(9) = a(5) * b(2);
    // x^2
    output(10) = a(0) * b(3) + a(6) * b(0);
    // xy
    output(11) = a(1) * b(3) + a(6) * b(1) + a(7) * b(0);
    // y^2
    output(12) = a(2) * b(3) + a(7) * b(1);
    // xz
    output(13) = a(3) * b(3) + a(6) * b(2) + a(8) * b(0);
    // yz
    output(14) = a(4) * b(3) + a(7) * b(2) + a(8) * b(1);
    // z^2
    output(15) = a(5) * b(3) + a(8) * b(2);
    // x
    output(16) = a(6) * b(3) + a(9) * b(0);
    // y
    output(17) = a(7) * b(3) + a(9) * b(1);
    // z
    output(18) = a(8) * b(3) + a(9) * b(2);
    // 1
    output(19) = a(9) * b(3);
    return output;
}

/// Script extracted from MAGSAC implementation.
// Shorthand for multiplying the Essential matrix with its transpose.
Eigen::Matrix<double, 1, 10> computeEETranspose(
        const Eigen::Matrix<double, 1, 4> nullSpace[3][3],
        int i,
        int j) {
    return multiplyDegOnePoly(nullSpace[i][0], nullSpace[j][0]) +
           multiplyDegOnePoly(nullSpace[i][1], nullSpace[j][1]) +
           multiplyDegOnePoly(nullSpace[i][2], nullSpace[j][2]);
}

/// Script extracted from MAGSAC implementation.
// Builds the trace constraint: EEtE - 1/2 trace(EEt)E = 0
Eigen::Matrix<double, 9, 20> getTraceConstraint(
        const Eigen::Matrix<double, 1, 4> nullSpace[3][3]) {
    Eigen::Matrix<double, 9, 20> traceConstraint;

    // Compute EEt.
    Eigen::Matrix<double, 1, 10> eet[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            eet[i][j] = 2 * computeEETranspose(nullSpace, i, j);
        }
    }

    // Compute the trace.
    const Eigen::Matrix<double, 1, 10> trace = eet[0][0] + eet[1][1] + eet[2][2];

    // Multiply EEt with E.
    for (auto i = 0; i < 3; i++) {
        for (auto j = 0; j < 3; j++) {
            traceConstraint.row(3 * i + j) = multiplyDegTwoDegOnePoly(eet[i][0], nullSpace[0][j]) +
                                             multiplyDegTwoDegOnePoly(eet[i][1], nullSpace[1][j]) +
                                             multiplyDegTwoDegOnePoly(eet[i][2], nullSpace[2][j]) -
                                             0.5 * multiplyDegTwoDegOnePoly(trace, nullSpace[i][j]);
        }
    }

    return traceConstraint;
}

/// Script extracted from MAGSAC implementation.
Eigen::Matrix<double, 1, 20> getDeterminantConstraint(
        const Eigen::Matrix<double, 1, 4> nullSpace[3][3]) {
    // Singularity constraint.
    return multiplyDegTwoDegOnePoly(
            multiplyDegOnePoly(nullSpace[0][1], nullSpace[1][2]) -
            multiplyDegOnePoly(nullSpace[0][2], nullSpace[1][1]),
            nullSpace[2][0]) +
           multiplyDegTwoDegOnePoly(
                   multiplyDegOnePoly(nullSpace[0][2], nullSpace[1][0]) -
                   multiplyDegOnePoly(nullSpace[0][0], nullSpace[1][2]),
                   nullSpace[2][1]) +
           multiplyDegTwoDegOnePoly(
                   multiplyDegOnePoly(nullSpace[0][0], nullSpace[1][1]) -
                   multiplyDegOnePoly(nullSpace[0][1], nullSpace[1][0]),
                   nullSpace[2][2]);
}

/// Script extracted from MAGSAC implementation.
Eigen::Matrix<double, 10, 20> buildConstraintMatrix(
        const Eigen::Matrix<double, 1, 4> nullSpace[3][3]) {
    Eigen::Matrix<double, 10, 20> constraintMatrix;
    constraintMatrix.block<9, 20>(0, 0) = getTraceConstraint(nullSpace);
    constraintMatrix.row(9) = getDeterminantConstraint(nullSpace);
    return constraintMatrix;
}
