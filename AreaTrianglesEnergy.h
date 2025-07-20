#ifndef AREATRIANGLESENERGY_H_
#define AREATRIANGLESENERGY_H_

#include "SystemStructures.h"

// Include the header file "SystemStructures.h" which contains the necessary structure definitions.

// The following returns a double value that is the energy of a triangle I believe
double ComputeAreaTrianglesEnergy(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    AreaTriangleInfoVecs& areaTriangleInfoVecs);

struct AreaEnergyFunctor {
	// Declare a struct named AreaEnergyFunctor.

    double area_0;
    double spring_constant;
    double* locXAddr;
    double* locYAddr;
    double* locZAddr;
    int* idKey;

    // Declare member variables of the struct.

    __host__ __device__ AreaEnergyFunctor(
        double& _area_0,
        double& _spring_constant,
        double* _locXAddr,
        double* _locYAddr,
        double* _locZAddr,
        int* _idKey):
        area_0(_area_0),
        spring_constant(_spring_constant),
        locXAddr(_locXAddr),
        locYAddr(_locYAddr),
        locZAddr(_locZAddr),
        idKey(_idKey) {}

    // Declare a constructor for the struct.
    // The constructor takes references to doubles and pointers as arguments and initializes the member variables with the provided values.

    //hand in counting iterator and id of triangle
    __device__ double operator()(const Tuuuu &u4) {
        // Declare a member function named operator() that takes an argument of type Tuuuu and returns a double.

        int counter = thrust::get<0>(u4);
        int place = 3 * counter;
        // Retrieve the first element of the tuple u4 and assign it to the variable counter.
        // Calculate the value of place as 3 times the value of counter.

        int id_i = thrust::get<1>(u4);
        int id_j = thrust::get<2>(u4);
        int id_k = thrust::get<3>(u4);
        // Retrieve the second, third, and fourth elements of the tuple u4 and assign them to the variables id_i, id_j, and id_k, respectively.

        CVec3 ri = thrust::make_tuple<double>(locXAddr[id_i], locYAddr[id_i], locZAddr[id_i]);
        CVec3 rj = thrust::make_tuple<double>(locXAddr[id_j], locYAddr[id_j], locZAddr[id_j]);
        CVec3 rk = thrust::make_tuple<double>(locXAddr[id_k], locYAddr[id_k], locZAddr[id_k]);
        // Create CVec3 tuples ri, rj, and rk using the values from locXAddr, locYAddr, and locZAddr arrays based on the indices id_i, id_j, and id_k.

        CVec3 rkj = CVec3_subtract(rk, rj);
        CVec3 rij = CVec3_subtract(ri, rj);
        // Calculate the differences rkj and rij by subtracting rj from rk and ri, respectively.

        double area_current = sqrt( CVec3_dot( CVec3_cross(rkj, rij), CVec3_cross(rkj, rij) ) )/2;
        // Calculate the current area of the triangle using cross products and dot products of rkj and rij vectors.

        double energy =  spring_constant/2 * (area_current - area_0) * (area_current - area_0) / area_0;
        // Calculate the energy of the triangle using the spring_constant, area_current, and area_0.

        return energy;
        // Return the calculated energy.
    };
};

#endif //AREATRIANGLESENERGY_H_
