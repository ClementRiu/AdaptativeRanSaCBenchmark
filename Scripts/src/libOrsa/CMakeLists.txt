SET(SRC conditioning.hpp conditioning.cpp
        match.hpp
        match2d3d.hpp
        model_estimator.hpp model_estimator.cpp
        homography_model.hpp homography_model.cpp
        fundamental_model.hpp fundamental_model.cpp
        essential_model.hpp essential_model.cpp
        e5p.h e5p.cpp
        line_model.hpp line_model.cpp
        pnp_model.hpp pnp_model.cpp
        pnp_poly.hpp pnp_poly.cpp
        sampling.hpp sampling.cpp
        ransac_algorithm.hpp ransac_algorithm.cpp
        ransac.hpp ransac.cpp
        starsac.hpp starsac.cpp
        muse.hpp muse.cpp
        muse_lookup.hpp muse_lookup.cpp
        orsa.hpp orsa.cpp
        orsa_fast.hpp orsa_fast.cpp
        lrtsac.hpp lrtsac.cpp
        magsac.hpp magsac.cpp
        gamma_values.cpp
        magsac_metrics.hpp magsac_metrics.cpp)

ADD_LIBRARY(Orsa ${SRC})

IF(UNIX)
    SET_TARGET_PROPERTIES(Orsa PROPERTIES
                          COMPILE_FLAGS "-std=c++11 -Wall -Wextra")
ENDIF(UNIX)

IF(MSVC)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
ENDIF(MSVC)

ADD_SUBDIRECTORY(libNumerics)

TARGET_LINK_LIBRARIES(Orsa Numerics)

#----------
# UNIT TEST
#----------

UNIT_TEST(fundamental_model Orsa)
UNIT_TEST(fundamental_orsa Orsa)
UNIT_TEST(homography_model Orsa)
UNIT_TEST(homography_orsa Orsa)
