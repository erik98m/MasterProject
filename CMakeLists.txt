cmake_minimum_required(VERSION 3.12)
project(masteroppgave LANGUAGES CXX)


set(SOURCEFILES src/main.cpp src/fitness.cpp src/neighbor.cpp)

add_executable(hillclimber ${SOURCEFILES})
set_target_properties(hillclimber PROPERTIES CXX_STANDARD 17)

# === Set up Paradiseo library ===
set(PARADISEO_SRC_DIR ${CMAKE_SOURCE_DIR}/paradiseo)
set(PARADISEO_BUILD_DIR ${PARADISEO_SRC_DIR}/build)

add_library(paradiseo INTERFACE)
target_include_directories(paradiseo INTERFACE ${PARADISEO_SRC_DIR}/eo/src)
target_include_directories(paradiseo INTERFACE ${PARADISEO_SRC_DIR}/mo/src)
target_include_directories(paradiseo INTERFACE ${PARADISEO_SRC_DIR}/problems)
target_compile_definitions(paradiseo INTERFACE WITH_EIGEN=1)
target_include_directories(paradiseo INTERFACE /usr/include/eigen3)
target_link_directories(paradiseo INTERFACE ${PARADISEO_BUILD_DIR}/lib)
target_link_libraries(paradiseo INTERFACE -leo -leoutils -les)

target_link_libraries(hillclimber PRIVATE paradiseo)

add_definitions(-DCMAKE_EXPORT_COMPILE_COMMANDS=ON)