PROJECT( mt )

##############################################################################
# POLICIES ###################################################################
##############################################################################

cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)

if (POLICY CMP0057)
  cmake_policy(SET CMP0057 NEW)
endif()

##############################################################################
# OPTIONS ####################################################################
##############################################################################

OPTION(useFFTW "Turn FTTW on and off" ON)
OPTION(usePython2 "Turn Python2.7 on and off" ON)
MESSAGE("useFFTW = ${useFFTW}")
MESSAGE("usePython2 = ${usePython2}")

##############################################################################
# SEARCH FOR SOURCE FILES ####################################################
##############################################################################

file(GLOB_RECURSE SRC_PATH "src/*.cpp")
MESSAGE("Sources: \n ${SRC_PATH}")
file(GLOB_RECURSE QT_RESOURCES "src/*.qrc")
MESSAGE("Sources: \n ${QT_RESOURCES}")

# sources files
Set (PROJECT_SRCS
    ${SRC_PATH}
)

##############################################################################
# FIND LIBARIES ##############################################################
##############################################################################

##############here starts shark compilation dependent setup##########################

#map default build type to release
IF(NOT CMAKE_BUILD_TYPE)
       SET(CMAKE_BUILD_TYPE "Release")
ENDIF(NOT CMAKE_BUILD_TYPE)
IF((CMAKE_BUILD_TYPE MATCHES "debug") OR (CMAKE_BUILD_TYPE MATCHES "DEBUG"))
    SET(CMAKE_BUILD_TYPE "Debug")
ENDIF((CMAKE_BUILD_TYPE MATCHES "debug") OR (CMAKE_BUILD_TYPE MATCHES "DEBUG"))
IF((CMAKE_BUILD_TYPE MATCHES "release") OR (CMAKE_BUILD_TYPE MATCHES "RELEASE"))
    SET(CMAKE_BUILD_TYPE "Release")
ENDIF((CMAKE_BUILD_TYPE MATCHES "release") OR (CMAKE_BUILD_TYPE MATCHES "RELEASE"))
#enable or disable debugging
IF(CMAKE_BUILD_TYPE MATCHES "Debug")
    ADD_DEFINITIONS(-DDEBUG )
    SET( LINK_LIBRARIES
        /home/dith/builds/boost_1_55_0/stage/lib/libboost_system.a;/home/dith/builds/boost_1_55_0/stage/lib/libboost_date_time.a;/home/dith/builds/boost_1_55_0/stage/lib/libboost_filesystem.a;/home/dith/builds/boost_1_55_0/stage/lib/libboost_program_options.a;/home/dith/builds/boost_1_55_0/stage/lib/libboost_serialization.a;/home/dith/builds/boost_1_55_0/stage/lib/libboost_thread.a;/home/dith/builds/boost_1_55_0/stage/lib/libboost_unit_test_framework.a
        /usr/local/lib/libshark_debug.a
    )
ENDIF (CMAKE_BUILD_TYPE MATCHES "Debug")
IF(CMAKE_BUILD_TYPE MATCHES "Release")
    ADD_DEFINITIONS(-DNDEBUG )
    ADD_DEFINITIONS(-DBOOST_UBLAS_NDEBUG )
    SET( LINK_LIBRARIES
        /home/dith/builds/boost_1_55_0/stage/lib/libboost_system.a;/home/dith/builds/boost_1_55_0/stage/lib/libboost_date_time.a;/home/dith/builds/boost_1_55_0/stage/lib/libboost_filesystem.a;/home/dith/builds/boost_1_55_0/stage/lib/libboost_program_options.a;/home/dith/builds/boost_1_55_0/stage/lib/libboost_serialization.a;/home/dith/builds/boost_1_55_0/stage/lib/libboost_thread.a;/home/dith/builds/boost_1_55_0/stage/lib/libboost_unit_test_framework.a
        /usr/local/lib/libshark.a
    )
ENDIF (CMAKE_BUILD_TYPE MATCHES "Release")

#handle GotoBLAS
OPTION(OPT_ENABLE_GOTOBLAS "Enable GOTOBLAS as Linear Algebra library" )
mark_as_advanced(
    OPT_ENABLE_GOTOBLAS
)

IF ( OPT_ENABLE_GOTOBLAS )
    SET(LINK_LIBRARIES
        ${LINK_LIBRARIES}
    )
    ADD_DEFINITIONS( -DSHARK_USE_GOTOBLAS)
    include_directories ( /include/ )
ENDIF ( OPT_ENABLE_GOTOBLAS )

ADD_DEFINITIONS( -DSHARK_USE_OPENMP)

####################here Ends shark dependent setup########################

#find_package(Shark REQUIRED)
#include(${SHARK_USE_FILE})
#set_property(TARGET ExampleProject PROPERTY CXX_STANDARD 11)

SET(EIGEN3_INCLUDE_DIR SYSTEM /usr/include/eigen3)

find_package(VTK REQUIRED)
if (VTK_FOUND)
        include(${VTK_USE_FILE})
        #include_directories(${VTK_INCLUDE_DIRS})
        message("VTK found, VTK_USE_FILE: ${VTK_USE_FILE}")
else()
        message(ERROR "VTK not found")
endif()

IF (useFFTW)
    SET(FFTW_LINKER_FLAGS "-lfftw3_threads -lfftw3")
ENDIF (useFFTW)

IF (usePython2)
    SET(PYTHON27_LINKER_FLAGS "-lpython2.7")
    SET(PYTHON27_INCLUDE_DIR SYSTEM /usr/include/python2.7)
ENDIF (usePython2)

find_package(Qt5Widgets REQUIRED)

set(CMAKE_AUTOMOC ON)

qt5_add_resources(UI_RESOURCES src/application.qrc)

ADD_EXECUTABLE( mt ${PROJECT_SRCS} ${UI_RESOURCES} )

target_compile_features(mt PRIVATE cxx_range_for)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
##############################################################################
# SET FLAGS DIRECTORIES AND DEFINITIONS ######################################
##############################################################################

set(CMAKE_CXX_FLAGS " -O3 -Wall -Wextra -fopenmp -Wno-sign-compare")
SET(CMAKE_EXE_LINKER_FLAGS "-lm -ltiff ${FFTW_LINKER_FLAGS} ${PYTHON27_LINKER_FLAGS}")

MESSAGE("Linker flags: ${CMAKE_EXE_LINKER_FLAGS}")

INCLUDE_DIRECTORIES (
    /usr/include
    /usr/local/include
    /home/dith/builds/boost_1_55_0
    ${PYTHON27_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIR}
)

##############################################################################
# TARGET LINK LIBARIES #######################################################
##############################################################################

TARGET_LINK_LIBRARIES( mt ${LINK_LIBRARIES} ${VTK_LIBRARIES} ${SHARK_LIBRARIES} Qt5::Widgets )
