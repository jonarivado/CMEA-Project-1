# Install script for directory: /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/unsupported/Eigen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen_install")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump.exe")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE FILE FILES
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/unsupported/Eigen/AdolcForward"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/unsupported/Eigen/AlignedVector3"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/unsupported/Eigen/ArpackSupport"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/unsupported/Eigen/AutoDiff"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/unsupported/Eigen/BVH"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/unsupported/Eigen/EulerAngles"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/unsupported/Eigen/FFT"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/unsupported/Eigen/IterativeSolvers"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/unsupported/Eigen/KroneckerProduct"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/unsupported/Eigen/LevenbergMarquardt"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/unsupported/Eigen/MatrixFunctions"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/unsupported/Eigen/MoreVectorization"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/unsupported/Eigen/MPRealSupport"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/unsupported/Eigen/NonLinearOptimization"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/unsupported/Eigen/NumericalDiff"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/unsupported/Eigen/OpenGLSupport"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/unsupported/Eigen/Polynomials"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/unsupported/Eigen/Skyline"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/unsupported/Eigen/SparseExtra"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/unsupported/Eigen/SpecialFunctions"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/unsupported/Eigen/Splines"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE DIRECTORY FILES "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/unsupported/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen-prefix/src/Eigen-build/unsupported/Eigen/CXX11/cmake_install.cmake")

endif()

