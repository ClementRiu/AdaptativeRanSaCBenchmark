PROJECT(libNumerics)

SET(SRC
        cubicRoots.h
        ccmath_svd.cpp
        homography.h
        matrix.cpp matrix.h
        numerics.cpp numerics.h
        vector.cpp)

ADD_LIBRARY(Numerics ${SRC})

IF(UNIX)
    SET_TARGET_PROPERTIES(Numerics PROPERTIES
                          COMPILE_FLAGS "-std=c++11 -Wall -Wextra")
ENDIF(UNIX)

#----------
# UNIT TEST
#----------
UNIT_TEST(cubicRoots Numerics)
UNIT_TEST(ccmath_svd Numerics)
