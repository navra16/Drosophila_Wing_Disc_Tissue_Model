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

ILIBSCPP := -I/afs/crc.nd.edu/x86_64_linux/c/cuda/8.0/include/
ILIBS1 := -I/afs/crc.nd.edu/x86_64_linux/c/cuda/9.2/include/

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../AreaTriangles.cu \
../BendingTriangles.cu \
../MemRepulsionSprings.cu\
../LinearSprings.cu \
../LJSprings.cu \
../LJSprings_LJ.cu \
../VolumeComp.cu \
../VolumeSprings.cu \
../NodeAdvance.cu \
../AreaTrianglesEnergy.cu \
../BendingTrianglesEnergy.cu \
../LinearSpringsEnergy.cu \
../MemRepulsionEnergy.cu \
../System.cu \
../Edgeswap_test.cpp \
../SystemBuilder.cpp \
../Storage.cpp \
../main.cpp


# this is a variable
OBJS += \
./AreaTriangles.o \
./BendingTriangles.o \
./MemRepulsionSprings.o\
./LinearSprings.o \
./LJSprings.o \
./LJSprings_LJ.o \
./VolumeComp.o \
./VolumeSprings.o \
./NodeAdvance.o \
./AreaTrianglesEnergy.o \
./BendingTrianglesEnergy.o \
./LinearSpringsEnergy.o \
./MemRepulsionEnergy.o \
./System.o \
./Edgeswap_test.o \
./SystemBuilder.o \
./Storage.o \
./main.o


CPP_DEPS += \
./AreaTriangles.d \
./BendingTriangles.d \
./MemRepulsionSprings.d\
./LinearSprings.d \
./LJSprings.d \
./LJSPrings_LJ.d \
./VolumeComp.d \
./VolumeSprings.d \
./NodeAdvance.d \
./AreaTrianglesEnergy.d \
./BendingTrianglesEnergy.d \
./LinearSpringsEnergy.d \
./MemRepulsionEnergy.d \
./System.d \
./Edgeswap_test.d \
./SystemBuilder.d \
./Storage.d \
./main.d

#need o have ILIBS2
#cpp files
%.o : ./%.cpp 
	 $(CXX) $(CFLAGS) $(ILIBSCPP) $(LIBS) -o $@ -c $^ 

	
#cuda files
%.o : ./%.cu 
	$(NVCC) $(NVCCFLAGS) $(ILIBS1) -dc -o $@ $^
