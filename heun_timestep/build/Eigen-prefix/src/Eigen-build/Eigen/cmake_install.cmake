# Install script for directory: /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen" TYPE FILE FILES
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/Cholesky"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/CholmodSupport"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/Core"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/Dense"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/Eigen"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/Eigenvalues"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/Geometry"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/Householder"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/IterativeLinearSolvers"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/Jacobi"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/LU"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/MetisSupport"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/OrderingMethods"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/PaStiXSupport"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/PardisoSupport"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/QR"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/QtAlignedMalloc"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/SPQRSupport"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/SVD"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/Sparse"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/SparseCholesky"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/SparseCore"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/SparseLU"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/SparseQR"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/StdDeque"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/StdList"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/StdVector"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/SuperLUSupport"
    "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/UmfPackSupport"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen" TYPE DIRECTORY FILES "/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

