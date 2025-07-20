
#include "System.h"
#include "SystemStructures.h"
#include "LinearSpringsEnergy.h"

double ComputeLinearSpringsEnergy(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    LJInfoVecs& ljInfoVecs) {

    
        thrust::counting_iterator<int> edgeIdBegin(0);
        thrust::counting_iterator<int> edgeIdEnd(coordInfoVecs.num_edges);


    linearSpringInfoVecs.linear_spring_energy=
    thrust::transform_reduce(  
        thrust::make_zip_iterator(
            thrust::make_tuple(
                edgeIdBegin,
                coordInfoVecs.edges2Nodes_1.begin(),
                coordInfoVecs.edges2Nodes_2.begin(),
                linearSpringInfoVecs.edge_initial_length.begin())),
        thrust::make_zip_iterator( 
            thrust::make_tuple(
                edgeIdEnd,
                coordInfoVecs.edges2Nodes_1.end(),
                coordInfoVecs.edges2Nodes_2.end(), 
                linearSpringInfoVecs.edge_initial_length.end())),
        LinearSpringEnergyFunctor(
            linearSpringInfoVecs.spring_constant, 
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()),

            thrust::raw_pointer_cast(linearSpringInfoVecs.tempNodeIdUnreduced.data())),
            //thrust::raw_pointer_cast(linearSpringInfoVecs.tempNodeForceXUnreduced.data()),
            //thrust::raw_pointer_cast(linearSpringInfoVecs.tempNodeForceYUnreduced.data()),
            //thrust::raw_pointer_cast(linearSpringInfoVecs.tempNodeForceZUnreduced.data())),
        0.0, thrust::plus<double>() ); 
        //std::cout<<"linear energy from energy.cu: "<< linearSpringInfoVecs.linear_spring_energy<<std::endl;
        
        printf("current_len = %g, zero_len = %g\n", length_current, length_zero);
        return linearSpringInfoVecs.linear_spring_energy; 

};