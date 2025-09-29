cmake_minimum_required(VERSION 3.6)
set(NAME GammaCube)
set(CMAKE_CXX_STANDARD 17)
project(${NAME})

set (CMAKE_PREFIX_PATH /GEANT4/)
set (CMAKE_PREFIX_PATH /opt/GEANT4/geant4-v11.3.0-install)

option(WITH_GEANT4_UIVIS "Build project with Geant4 UI and Vis drivers" ON)
if (WITH_GEANT4_UIVIS)
    find_package(Geant4 REQUIRED ui_all vis_all)
else ()
    find_package(Geant4 REQUIRED)
endif ()

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)


file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/src/Flux/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/src/Flux/*.hh)

add_executable(${NAME} GammaCube.cc ${sources} ${headers})
target_link_libraries(${NAME} ${Geant4_LIBRARIES})
