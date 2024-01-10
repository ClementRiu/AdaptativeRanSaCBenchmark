SET(SRC artificial_output.hpp artificial_output.cpp
        cmdLine.h
        data_handler.hpp
        fundamental_graphical_output.hpp fundamental_graphical_output.cpp
        graphical_output.hpp graphical_output.cpp
        homography_graphical_output.hpp homography_graphical_output.cpp
        histogram.hpp
        metrics.hpp metrics.cpp
        Rect.hpp Rect.cpp
        siftMatch.hpp siftMatch.cpp
        usac_file_handler.hpp usac_file_handler.cpp
        warping.hpp warping.cpp
        )

ADD_LIBRARY(Utilities ${SRC})

IF(UNIX)
    SET_TARGET_PROPERTIES(Utilities PROPERTIES
            COMPILE_FLAGS "-std=c++11 -Wall -Wextra")
ENDIF(UNIX)

ADD_SUBDIRECTORY(libImage)

TARGET_LINK_LIBRARIES(Utilities image)
