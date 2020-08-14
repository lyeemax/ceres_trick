cmake_minimum_required(VERSION 3.15)
project(ceres_debug_trick)

set(CMAKE_CXX_STANDARD 14)
find_package(OpenCV REQUIRED)
find_package(Eigen3 3.3.4 REQUIRED)
find_package(Ceres 2.0 REQUIRED)
find_package(PCL REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIRS})
include_directories(${CERES_INCLUDE_DIRS})
include_directories(${PCL_INCLUDE_DIRS})
include_directories(io)
include_directories(factor)
include_directories(utils)

set(LIBS ${OpenCV_LIBS} ${CERES_LIBRARIES} ${PCL_LIBRARIES})

add_executable(curve_fitting src/curve_fitting.cpp factor/curve_factor.h )
target_link_libraries(curve_fitting ${CERES_LIBRARIES})

add_executable(sphere  src/g2o_sphere.cpp factor/posegraph_factor.h utils/quaternion_utils.h factor/PoseLocalParameter.h)
target_link_libraries(sphere ${CERES_LIBRARIES})

add_executable(sphere_analytic  src/g2o_sphere_analytical.cpp factor/posegraph_factor.h utils/quaternion_utils.h factor/PoseLocalParameter.h)
target_link_libraries(sphere_analytic ${CERES_LIBRARIES})

add_executable(slam2d src/slam2d.cpp factor/posegraph2d_factor.h)
target_link_libraries(slam2d  ${CERES_LIBRARIES})

add_executable(tum_demo src/tum_rgbd.cpp factor/projection_factor.h factor/PoseLocalParameter.h)
target_link_libraries(tum_demo ${LIBS})