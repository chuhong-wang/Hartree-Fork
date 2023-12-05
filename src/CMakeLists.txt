add_library (src molecule.cpp point.cpp scf.cpp)

target_include_directories(src PUBLIC "${PROJECT_SOURCE_DIR}/include" "${PROJECT_SOURCE_DIR}/third_party")