add_library(BiCGStab bicgstab.cpp)
target_include_directories(BiCGStab PUBLIC ${PROJECT_SOURCE_DIR}/include_bicgstab)