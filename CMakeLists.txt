cmake_minimum_required(VERSION 3.10)
project(RayTracing)

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(OpenCV REQUIRED)

set(CMAKE_CXX_STANDARD 17)

add_executable(RayTracing main.cpp dip.cpp Object.hpp Vector.cpp Vector.hpp global.hpp Triangle.hpp Scene.cpp
        Scene.hpp Light.hpp AreaLight.hpp BVH.cpp BVH.hpp Bounds3.hpp Ray.hpp Material.hpp Intersection.hpp Texture.hpp 
        Renderer.cpp Renderer.hpp)
target_link_libraries(RayTracing ${OpenCV_LIBRARIES})
