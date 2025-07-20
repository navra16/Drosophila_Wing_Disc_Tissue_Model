#ifndef VOLUMECOMP_H_
#define VOLUMECOMP_H_ 

#include "SystemStructures.h"

void ComputeVolume(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    LJInfoVecs& ljInfoVecs);
    
// Functor for computing volume using the provided node locations and triangle information.
struct VolumeCompFunctor {
    
    double spring_constant;  // The spring constant used in volume calculations.
    double* locXAddr; // Raw pointer to the X-coordinates of node locations.
    double* locYAddr; // Raw pointer to the Y-coordinates of node locations.
    double* locZAddr; // Raw pointer to the Z-coordinates of node locations.

	__host__ __device__ VolumeCompFunctor(
        double& _spring_constant,
        double* _locXAddr,
        double* _locYAddr,
        double* _locZAddr):

        spring_constant(_spring_constant),
        locXAddr(_locXAddr),
        locYAddr(_locYAddr),
        locZAddr(_locZAddr){}

	//hand in counting iterator and id of two edges and preferred length
	// Function for computing the volume using three nodes (r1, r2, r3) of a triangle.
    // The input is a tuple u4 containing (counter, r1, r2, r3).
    __device__ double operator()(const Tuuuu &u4) {
        		
        //counter ranges from 0 to num_edges. 
        // Extract the elements from the input tuple.
        int counter = thrust::get<0>(u4); // Counter ranges from 0 to num_edges.
        int r1 = thrust::get<1>(u4); // Index of the first node of the triangle.
        int r2 = thrust::get<2>(u4); // Index of the second node of the triangle.
        int r3 = thrust::get<3>(u4); // Index of the third node of the triangle.

        // Check if all nodes of the triangle are valid (not equal to INT_MAX).
        if (r1 != INT_MAX && r2 != INT_MAX && r3 != INT_MAX){
        // Get the X, Y, Z coordinates of the nodes of the triangle.
            double r1x = locXAddr[r1];
        double r1y = locYAddr[r1];
        double r1z = locZAddr[r1];
        double r2x = locXAddr[r2];
        double r2y = locYAddr[r2];
        double r2z = locZAddr[r2];
        double r3x = locXAddr[r3];
        double r3y = locYAddr[r3];
        double r3z = locZAddr[r3];

        // Compute the components of the normal vector (N) for the triangle.
            double N1 = (r2y - r1y)*(r3z - r1z) - (r3y - r1y)*(r2z - r1z);
        double N2 = -(r2x - r1x)*(r3z - r1z) + (r3x - r1x)*(r2z - r1z);
        double N3 = (r2x - r1x)*(r3y - r1y) - (r3x - r1x)*(r2y - r1y);
        
        
        // Compute the norm of the normal vector.
            double normN = sqrt(N1*N1 + N2*N2 + N3*N3);
        // Normalize the components of the normal vector.
            N1 = N1/normN;
        N2 = N2/normN;
        N3 = N3/normN;
        
        // Compute the dot product between the position vector (r1) and the normalized normal vector (N).
            double r1_dot_N = r1x*N1 + r1y*N2 + r1z*N3;
        
        // Compute the cross products between the nodes of the triangle to determine the signed volume of the tetrahedron.
            double r1cr2x = r1y*r2z - r2y*r1z;
        double r1cr2y = -r1x*r2z + r2x*r1z;
        double r1cr2z = r1x*r2y - r2x*r1y;
        double r2cr3x = r2y*r3z - r3y*r2z;
        double r2cr3y = -r2x*r3z + r3x*r2z;
        double r2cr3z = r2x*r3y - r3x*r2y;
        double r3cr1x = r3y*r1z - r1y*r3z;
        double r3cr1y = -r3x*r1z + r1x*r3z;
        double r3cr1z = r3x*r1y - r1x*r3y;

        // Compute the dot product between the normalized normal vector (N) and the sum of cross products.
            double NN = N1*(r1cr2x + r2cr3x + r3cr1x) + N2*(r1cr2y + r2cr3y + r3cr1y) + N3*(r1cr2z + r2cr3z + r3cr1z);
        
        
        // Compute the signed volume of the tetrahedron using the dot products.
            double volume = r1_dot_N*sqrt(NN*NN);


        return volume;
        
        }
        else{
            
            // If any of the nodes of the triangle are invalid, return a volume of 0.0.
            double volume = 0.0;
            
            return volume;
        }


    }
};

#endif

//In summary, the VolumeComp.h header file defines the ComputeVolume function and a VolumeCompFunctor struct. The ComputeVolume function's purpose is to compute the volume of a structure based on its node locations and triangle information. The VolumeCompFunctor is used as a functor for the computation, where it takes a tuple of (counter, r1, r2, r3) representing three nodes forming a triangle and calculates the signed volume of the tetrahedron formed by these nodes. The functor utilizes the provided node locations to perform the necessary calculations, and the ComputeVolume function allows for the overall volume calculation for the entire structure.