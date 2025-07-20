#include "System.h"
#include "VolumeSprings.h"



// WARNING: function must not reset coordInfoVecs.nodeForceX etc.
void ComputeVolumeSprings(
    CoordInfoVecs& coordInfoVecs,
	  LinearSpringInfoVecs& linearSpringInfoVecs, 
	  CapsidInfoVecs& capsidInfoVecs,
    GeneralParams& generalParams,
    AuxVecs& auxVecs) {    
    
    // The purpose of this function is to compute volume spring forces for the nodes in the system.
    // The volume springs are used to maintain the volume of a structure, and the forces help regulate
    // its shape.

    // Create a counting iterator starting from 0 and ending at generalParams.maxNodeCount.
    thrust::counting_iterator<int> begin(0);

    // Use thrust::transform to compute the volume spring forces for each node.
    thrust::transform(  
        // Input range: a zip iterator containing tuples of (index, id_bucket, nodeForceX, nodeForceY, nodeForceZ).
        thrust::make_zip_iterator(
            thrust::make_tuple(
                begin,
                auxVecs.id_bucket.begin(),
                coordInfoVecs.nodeForceX.begin(),
                coordInfoVecs.nodeForceY.begin(),
                coordInfoVecs.nodeForceZ.begin()
                )
            ),
        
        // Input range: a zip iterator containing tuples of (index, id_bucket, nodeForceX, nodeForceY, nodeForceZ)
        // incremented by generalParams.maxNodeCount.
        thrust::make_zip_iterator(
            thrust::make_tuple(
                begin,
                auxVecs.id_bucket.begin(),
                coordInfoVecs.nodeForceX.begin(),
                coordInfoVecs.nodeForceY.begin(),
                coordInfoVecs.nodeForceZ.begin()
                )
            ) + generalParams.maxNodeCount,

        // Output range: a zip iterator containing tuples of (nodeForceX, nodeForceY, nodeForceZ) for each node.
        thrust::make_zip_iterator(
            thrust::make_tuple(
                coordInfoVecs.nodeForceX.begin(),
                coordInfoVecs.nodeForceY.begin(),
                coordInfoVecs.nodeForceZ.begin())),

        // Functor that computes the volume spring forces for each node.
        VolumeSpringFunctor(
            generalParams.current_total_volume, // Current total volume of the structure.
            generalParams.true_current_total_volume, // True current total volume of the structure.
            generalParams.eq_total_volume, // Equilibrium total volume of the structure.
            generalParams.volume_spring_constant, // Spring constant for volume springs.
            coordInfoVecs.num_triangles, // Number of triangles in the structure.
            generalParams.Rmin, // Minimum distance used in volume calculations.

                // Raw pointers to the triangle information and node locations.
            thrust::raw_pointer_cast(coordInfoVecs.triangles2Nodes_1.data()),
            thrust::raw_pointer_cast(coordInfoVecs.triangles2Nodes_2.data()),
            thrust::raw_pointer_cast(coordInfoVecs.triangles2Nodes_3.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()), 
            
                      
            // Raw pointers to auxiliary vectors for the volume spring calculations.
            thrust::raw_pointer_cast(auxVecs.id_value_expanded.data()),
            thrust::raw_pointer_cast(auxVecs.keyBegin.data()),
            thrust::raw_pointer_cast(auxVecs.keyEnd.data())
            )
    );
            
            
    // The following loop prints the force exerted by volume springs on node 36 (index 35) after
    // the computation is done.
    for (int i = 0; i < generalParams.maxNodeCount; i++){
        std::cout<<"Force from volume on node 36 = "<<coordInfoVecs.nodeForceX[35]<<" "
            <<coordInfoVecs.nodeForceY[35]<<" "<<coordInfoVecs.nodeForceZ[35]<<std::endl;
    }
};

