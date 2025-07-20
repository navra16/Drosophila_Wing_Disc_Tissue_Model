################################################################################
# Automatically-generated file. Do not edit!
################################################################################
CXX := gcc
NVCC := nvcc
CFLAGS := -pthread -std=c++11 -Wall -Wextra
MFLAGS := -f
NVCCFLAGS := -std=c++11

current_dir := $(shell pwd)
LIBS:=  -lpugixml -L/$(current_dir)/pugixml/lib64
#-lgsl -lgslcblas

ILIBS_cuda8 = -I/opt/linux/centos/7.x/x86_64/pkgs/cuda/8.0/include/
ILIBS_cuda9 := -I/opt/linux/centos/7.x/x86_64/pkgs/cuda/9.1/include/

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../AreaTriangles.cu \
../LinearSprings.cu \
../VolumeComp.cu \
../VolumeSprings.cu \
../LineTensionSprings.cu \
../NodeAdvance.cu \
../System.cu \
../Utilities.cpp \
../StrainTensor.cu \
../SystemBuilder.cpp \
../Storage.cpp \
../main.cpp


# this is a variable
OBJS += \
./AreaTriangles.o \
./LinearSprings.o \
./VolumeComp.o \
./VolumeSprings.o \
./LineTensionSprings.o \
./NodeAdvance.o \
./System.o \
./Utilities.o \
./StrainTensor.o \
./SystemBuilder.o \
./Storage.o \
./main.o


CPP_DEPS += \
./AreaTriangles.d \
./LinearSprings.d \
./VolumeComp.d \
./VolumeSprings.d \
./LineTensionSprings.d \
./NodeAdvance.d \
./System.d \
./Utilities.d \
./StrainTensor.d \
./SystemBuilder.d \
./Storage.d \
./main.d

#need o have ILIBS2
#cpp files
%.o : ./%.cpp 
	 $(CXX) $(CFLAGS) $(ILIBS_cuda8) $(LIBS) -o $@ -c $^ 

	
#cuda files
%.o : ./%.cu 
	$(NVCC) $(NVCCFLAGS) $(ILIBS_cuda9) -dc -o $@ $^
