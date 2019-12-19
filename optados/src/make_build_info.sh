#!/bin/bash

BUILD=`git describe --tags`
COMPILER=`grep SYSTEM ../make.system | head -n 1 | awk '{print $3}'` 
BUILD_TYPE=`grep BUILD ../make.system | head -n 1 | awk '{print $3}'`
COMMS_ARCH=`grep COMMS_ARCH ../make.system | head -n 1 | awk '{print $3}'`
SOURCE_TIME=`git log | grep Date | head -n 1 |awk '{print $5}'`
SOURCE_DATE=`git log | grep Date | head -n 1 |awk '{print $2, $4, $3, $6}'`
COMPILE_TIME=`date +"%R %Z"`
COMPILE_DATE=`date +"%a %e %b %Y"`

echo "module od_build" > build.f90
echo " implicit none" >> build.f90
echo " ">> build.f90
echo " private ! unless otherwise stated" >> build.f90
echo " " >> build.f90
echo " type, public ::  build_info_type" >> build.f90 
echo "  character(len=20) :: build='"${BUILD}"'" >> build.f90
echo "  character(len=20) :: compiler='"${COMPILER}"'" >> build.f90
echo "  character(len=20) :: build_type='"${BUILD_TYPE}"'" >> build.f90
echo "  character(len=20) :: comms_arch='"${COMMS_ARCH}"'" >> build.f90
echo "  character(len=20) :: source_time='"${SOURCE_TIME}"'" >> build.f90
echo "  character(len=20) :: source_date='"${SOURCE_DATE}"'" >> build.f90
echo "  character(len=20) :: compile_date='"${COMPILE_DATE}"'" >> build.f90
echo "  character(len=20) :: compile_time='"${COMPILE_TIME}"'" >> build.f90
echo " end type build_info_type" >> build.f90

echo " type(build_info_type), public, save :: build_info" >> build.f90

echo "endmodule od_build" >> build.f90
