FIND_PACKAGE(OpenMP)
if (OPENMP_FOUND)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif (OPENMP_FOUND)

# Executable demo to estimate the fundamental matrix from SIFT matches
ADD_EXECUTABLE(residual_distrib
        residual_distrib.cpp
        apply_transform.hpp apply_transform.cpp
        read_write_experiment.hpp read_write_experiment.cpp
        )
TARGET_LINK_LIBRARIES(residual_distrib Orsa Utilities SiftAnatomy)

ADD_EXECUTABLE(outlier_gen_distrib
        outlier_gen_distrib.cpp
        apply_transform.hpp apply_transform.cpp
        read_write_experiment.hpp read_write_experiment.cpp
        )
TARGET_LINK_LIBRARIES(outlier_gen_distrib Orsa Utilities SiftAnatomy)

ADD_EXECUTABLE(residual_distrib_pnp
        residual_distrib_pnp.cpp
        apply_transform.hpp apply_transform.cpp
        read_write_experiment.hpp read_write_experiment.cpp
        )
TARGET_LINK_LIBRARIES(residual_distrib_pnp Orsa Utilities SiftAnatomy)

ADD_EXECUTABLE(outlier_gen_distrib_pnp
        outlier_gen_distrib_pnp.cpp
        apply_transform.hpp apply_transform.cpp
        read_write_experiment.hpp read_write_experiment.cpp
        )
TARGET_LINK_LIBRARIES(outlier_gen_distrib_pnp Orsa Utilities SiftAnatomy)

# Executable demo to estimate the fundamental matrix from SIFT matches
ADD_EXECUTABLE(generate_artificial_dataset
        generate_artificial_dataset.cpp
        apply_transform.hpp apply_transform.cpp
        read_write_experiment.hpp read_write_experiment.cpp
        )
TARGET_LINK_LIBRARIES(generate_artificial_dataset Orsa Utilities SiftAnatomy)

# Executable demo to estimate the fundamental matrix from SIFT matches
ADD_EXECUTABLE(pnp_generate_artificial_dataset
        pnp_generate_artificial_dataset.cpp
        apply_transform.hpp apply_transform.cpp
        read_write_experiment.hpp read_write_experiment.cpp
        )
TARGET_LINK_LIBRARIES(pnp_generate_artificial_dataset Orsa Utilities SiftAnatomy)

ADD_EXECUTABLE(general_experiment
        general_experiment.cpp
        apply_transform.hpp apply_transform.cpp
        read_write_experiment.hpp read_write_experiment.cpp
        )
TARGET_LINK_LIBRARIES(general_experiment Orsa Utilities SiftAnatomy)

ADD_EXECUTABLE(pnp_experiment
        pnp_experiment.cpp
        apply_transform.hpp apply_transform.cpp
        read_write_experiment.hpp read_write_experiment.cpp
        )
TARGET_LINK_LIBRARIES(pnp_experiment Orsa Utilities SiftAnatomy)

ADD_EXECUTABLE(new_algo_experiment
        new_algo_experiment.cpp
        read_write_experiment.hpp read_write_experiment.cpp
        )
TARGET_LINK_LIBRARIES(new_algo_experiment Orsa Utilities SiftAnatomy)

IF (UNIX)
    SET_TARGET_PROPERTIES(generate_artificial_dataset PROPERTIES
            COMPILE_FLAGS "-std=c++11 -Wall -Wextra")
    SET_TARGET_PROPERTIES(pnp_generate_artificial_dataset PROPERTIES
            COMPILE_FLAGS "-std=c++11 -Wall -Wextra")
    SET_TARGET_PROPERTIES(general_experiment PROPERTIES
            COMPILE_FLAGS "-std=c++11 -Wall -Wextra")
    SET_TARGET_PROPERTIES(pnp_experiment PROPERTIES
            COMPILE_FLAGS "-std=c++11 -Wall -Wextra")
    SET_TARGET_PROPERTIES(new_algo_experiment PROPERTIES
            COMPILE_FLAGS "-std=c++11 -Wall -Wextra")
ENDIF (UNIX)

ADD_EXECUTABLE(time_experiment
        time_experiment.cpp
        )
TARGET_LINK_LIBRARIES(time_experiment Orsa Utilities SiftAnatomy)

ADD_EXECUTABLE(lrt_experiment
        lrt_experiment.cpp
        read_write_experiment.hpp read_write_experiment.cpp
        apply_transform.hpp apply_transform.cpp
        )
TARGET_LINK_LIBRARIES(lrt_experiment Orsa Utilities SiftAnatomy)

ADD_EXECUTABLE(evaluate_artificial
        evaluate_artificial.cpp
        read_write_experiment.hpp read_write_experiment.cpp
        apply_transform.hpp apply_transform.cpp
        )
TARGET_LINK_LIBRARIES(evaluate_artificial Orsa Utilities SiftAnatomy)


IF (UNIX)
    SET_TARGET_PROPERTIES(time_experiment PROPERTIES
            COMPILE_FLAGS "-std=c++11 -Wall -Wextra")
    SET_TARGET_PROPERTIES(lrt_experiment PROPERTIES
            COMPILE_FLAGS "-std=c++11 -Wall -Wextra")
    SET_TARGET_PROPERTIES(evaluate_artificial PROPERTIES
            COMPILE_FLAGS "-std=c++11 -Wall -Wextra")
ENDIF (UNIX)
