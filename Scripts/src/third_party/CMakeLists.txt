#---------------------
#Third party libraries
#---------------------

add_subdirectory(CppUnitLite)
add_subdirectory(sift_anatomy)
add_subdirectory(eigen)

install(FILES CMakeLists.txt DESTINATION src/third_party)
