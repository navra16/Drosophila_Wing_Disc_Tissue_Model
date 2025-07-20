#include "System.h"
#include "SystemStructures.h"
#include "AreaTrianglesEnergy.h"

// Compute the total energy of area triangles.
double ComputeAreaTrianglesEnergy(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    AreaTriangleInfoVecs& areaTriangleInfoVecs) {

    // Create a counting iterator for element IDs.
    thrust::counting_iterator<int> elemId(0);
    
        
    // Compute the energy of each area triangle and sum them using transform_reduce.
    areaTriangleInfoVecs.area_triangle_energy = thrust::transform_reduce( 
			thrust::make_zip_iterator(
				thrust::make_tuple(
          elemId,
					coordInfoVecs.triangles2Nodes_1.begin(),
					coordInfoVecs.triangles2Nodes_2.begin(),
					coordInfoVecs.triangles2Nodes_3.begin())),
			thrust::make_zip_iterator(
				thrust::make_tuple(
                    elemId,
					coordInfoVecs.triangles2Nodes_1.begin(),
					coordInfoVecs.triangles2Nodes_2.begin(),
					coordInfoVecs.triangles2Nodes_3.begin())) + coordInfoVecs.num_triangles,
      AreaEnergyFunctor( 
          areaTriangleInfoVecs.initial_area,
          areaTriangleInfoVecs.spring_constant,
          thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
          thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
          thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()),
				  thrust::raw_pointer_cast(areaTriangleInfoVecs.tempNodeIdUnreduced.data())),
      0.0, thrust::plus<double>() );
      
      // Energy of triangles to be returned. Have to see why this is needed <<insert reason here>>.     
      return areaTriangleInfoVecs.area_triangle_energy;
};