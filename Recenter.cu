#include "SystemStructures.h"
#include "System.h" 
#include "Recenter.h" 
 
void Recentering(
    double old_center_x,
    double old_center_y,
    double old_center_z,
	CoordInfoVecs& coordInfoVecs,
	GeneralParams& generalParams) {

    
		//At this point, the previous node location is the same as the current node, 
		//we can therefore use previous node locations to update nodeLoc. 
		
		 
		  
		thrust::transform( 
			thrust::make_zip_iterator(
				thrust::make_tuple(
					coordInfoVecs.isNodeFixed.begin(),
					coordInfoVecs.nodeLocX.begin(),
					coordInfoVecs.nodeLocY.begin(),
					coordInfoVecs.nodeLocZ.begin())),
			thrust::make_zip_iterator(
				thrust::make_tuple(
					coordInfoVecs.isNodeFixed.begin(),
					coordInfoVecs.nodeLocX.begin(),
					coordInfoVecs.nodeLocY.begin(),
					coordInfoVecs.nodeLocZ.begin())) + generalParams.maxNodeCount,
			//second vector begin
			thrust::make_zip_iterator(
				thrust::make_tuple(
					coordInfoVecs.nodeForceX.begin(),
					coordInfoVecs.nodeForceY.begin(),
					coordInfoVecs.nodeForceZ.begin())),
			//save result in third vector to test values
			thrust::make_zip_iterator(
				thrust::make_tuple(
					coordInfoVecs.displacementx.begin(),
					coordInfoVecs.displacementy.begin(),
					coordInfoVecs.displacementz.begin())),
				SaxpyFunctorPrimary(
				generalParams.dt,
				generalParams.nodeMass,
				generalParams.maxNodeCount,
				domainParams.originMaxX,
				domainParams.originMaxY,
				domainParams.originMaxZ));
				
		//now that nodeLoc is different, we can calculate change and then set previous location
		//to the current location. 
	
}; 
