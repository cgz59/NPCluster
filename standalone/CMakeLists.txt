CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(npcluster)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")

if(APPLE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -DNDEBUG -msse4.2 -std=c++0x")
else(APPLE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -DNDEBUG -msse4.2 -std=c++0x")
endif(APPLE)


# OS Definitions
IF(UNIX)
  IF(APPLE)   #OSX
    ADD_DEFINITIONS(-DOS_MAC)

    SET(CMAKE_MACOSX_RPATH ON)
    SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
    SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

    LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
    IF("${isSystemDir}" STREQUAL "-1")
      SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    ENDIF("${isSystemDir}" STREQUAL "-1")
  ELSE(APPLE) #Linux
    ADD_DEFINITIONS(-DOS_LNX)
  ENDIF()
ELSE(${UNIX}) #Windows
  ADD_DEFINITIONS(-DOS_WIN -DNOMINMAX)
ENDIF()

# Version information
#INCLUDE(${CMAKE_MODULE_PATH}/Version.cmake)

include(cotire)

# Get R information
find_package(R REQUIRED)

#set(JNI_SOURCE_FILES
#	src/jni/dr_app_beagle_multidimensionalscaling_NativeMDSSingleton.cpp
#	src/factory.cpp
#	)

#add_library(mds_jni SHARED ${JNI_SOURCE_FILES})
#target_link_libraries(mds_jni mds_opencl)
#target_link_libraries(mds_jni ${TBB_LIBRARIES})

##
# Build benchmark
##

include_directories(${CMAKE_SOURCE_DIR}/../src)

#set(R_HOME /Library/Frameworks/R.framework/Resources)
#set(R_LIB_INC 
##	/Users/msuchard/Library/R/3.2
#	/Library/Frameworks/R.framework/Versions/3.2/Resources
#	)

include_directories(
	${R_HOME}/include
	${R_LIBDIR}/Rcpp/include
	${R_LIBDIR}/RInside/include
	${R_LIBDIR}/BH/include
	${R_LIBDIR}/RcppNT2/include
)

link_directories(
 	${R_HOME}/lib 
	${R_SHAREDLIBDIR}
	${R_LIBDIR}/RInside/lib
	${R_LIBDIR}/BoostCompiledLibraries/lib
	${R_LIBDIR}/RcppParallel/lib
)

add_executable(benchmark src/benchmark.cpp)
target_link_libraries(benchmark ${R_USED_LIBS})
target_link_libraries(benchmark tbb RInside Boost)
cotire(benchmark)



#target_link_libraries(benchmark R)
#add_library(tbb SHARED IMPORTED)
#add_library(R SHARED IMPORTED)
#target_link_libraries(benchmark mds_jni)
#target_link_libraries(benchmark ${TBB_LIBRARIES})
#target_link_libraries(benchmark boost_program_options)

