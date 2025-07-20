
#include "System.h"
#include "SystemStructures.h"
#include "VolumeComp.h"

void ComputeVolume(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    LJInfoVecs& ljInfoVecs) {  
    
    // Create a counting iterator for triangle IDs from 0 to num_triangles.
    thrust::counting_iterator<int> triangleIdBegin(0);
    // Note: The triangleIdEnd is not used in the code.

    // Calculate the current total volume of the system using thrust::transform_reduce.
    // The VolumeCompFunctor is used to compute the volume for each triangle.
    // It takes the counting iterator and the indices of three nodes (r1, r2, r3) of each triangle as input.
    // The output of this transform_reduce is the sum of all the computed volumes, which represents the current total volume.
    generalParams.current_total_volume = thrust::transform_reduce(  
        thrust::make_zip_iterator(
            thrust::make_tuple(
                triangleIdBegin,
                coordInfoVecs.triangles2Nodes_1.begin(),
                coordInfoVecs.triangles2Nodes_2.begin(),
                coordInfoVecs.triangles2Nodes_3.begin()
                )
            ),
        thrust::make_zip_iterator( 
            thrust::make_tuple(
                triangleIdBegin,
                coordInfoVecs.triangles2Nodes_1.begin(),
                coordInfoVecs.triangles2Nodes_2.begin(), 
                coordInfoVecs.triangles2Nodes_3.begin()
                )
            ) + coordInfoVecs.num_triangles,
        // VolumeCompFunctor is used to calculate the volume of each triangle.
        VolumeCompFunctor(
            linearSpringInfoVecs.spring_constant, 
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()) 
            ),
        0.0, 
        thrust::plus<double>() 
        ); 
        // This sum is the part without the absolute value and factor of (1/6) in the formula.

    // Calculate the true current total volume by taking the absolute value of the computed current total volume.
    generalParams.true_current_total_volume = sqrt(generalParams.current_total_volume*generalParams.current_total_volume)/6.0;
    
    // Calculate the volume energy using the volume spring constant and the deviation of the current total volume from the equilibrium volume.
    // The result is stored in generalParams.volume_energy.
    generalParams.volume_energy = generalParams.volume_spring_constant*(generalParams.true_current_total_volume - generalParams.eq_total_volume)*
                                        (generalParams.true_current_total_volume - generalParams.eq_total_volume)/
                                        (2.0*generalParams.Rmin*generalParams.Rmin*generalParams.Rmin*generalParams.eq_total_volume);

};


//Note: The code computes the total volume of a system using the VolumeCompFunctor and performs additional calculations related to the system's volume energy.