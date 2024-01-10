FIND_PACKAGE(OpenMP)
if (OPENMP_FOUND)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif (OPENMP_FOUND)

ADD_EXECUTABLE(demo demo.cpp)
TARGET_LINK_LIBRARIES(demo Orsa Utilities SiftAnatomy)

ADD_EXECUTABLE(demo_PNP demo_PNP.cpp)
TARGET_LINK_LIBRARIES(demo_PNP Orsa Utilities SiftAnatomy)

# Utility to overlay points and lines on images (used by IPOL demo)
ADD_EXECUTABLE(put_epipolar put_epipolar.cpp ../utilities/Rect.cpp)
TARGET_LINK_LIBRARIES(put_epipolar image)

IF (UNIX)
    SET_TARGET_PROPERTIES(demo PROPERTIES
            COMPILE_FLAGS "-std=c++11 -Wall -Wextra")
    SET_TARGET_PROPERTIES(demo_PNP PROPERTIES
            COMPILE_FLAGS "-std=c++11 -Wall -Wextra")
    SET_TARGET_PROPERTIES(put_epipolar PROPERTIES
            COMPILE_FLAGS "-std=c++11 -Wall -Wextra")
ENDIF (UNIX)

# Show small epipolar lines
ADD_EXECUTABLE(show_epi show_epi.cpp
        )
TARGET_LINK_LIBRARIES(show_epi Orsa Utilities SiftAnatomy)

IF (UNIX)
    SET_TARGET_PROPERTIES(show_epi PROPERTIES
            COMPILE_FLAGS "-std=c++11 -Wall -Wextra")
ENDIF (UNIX)

# Read and print an homography
ADD_EXECUTABLE(read_and_save_homography read_and_save_homography.cpp)
TARGET_LINK_LIBRARIES(read_and_save_homography Orsa Utilities SiftAnatomy)

IF (UNIX)
    SET_TARGET_PROPERTIES(read_and_save_homography PROPERTIES
            COMPILE_FLAGS "-std=c++11 -Wall -Wextra")
ENDIF (UNIX)
