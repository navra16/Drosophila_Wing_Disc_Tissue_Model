#pragma once
/******************************************************************************************************
 *  Spontaneous-strain engine for planar / weak-curvature spring meshes                               *
 *                                                                                                    *
 *  Implements the axi-symmetric growth tensor lambda(x,t) used in the Science Advances               *
 *  2024 wing-eversion paper (Eqs 22-29) by Furhmann et. al.                                                             *
 *                                                                                                    *
 *      1. buildVertexLambda()     lambda  at every vertex (radial, circumf., through-thick)          *
 *      2. updateEdgeRestLengths()     lambda-projected edge rest lengths (Eq 26)                     *
 *      3. updatePreferredAngles()     optional anisotropic spontaneous curvature  (not used here)    *
 *                                                                                                    *
 *  Author: Navaira Sherwani, 2025                                                                    *
 ******************************************************************************************************/

#include <thrust/device_vector.h>
#include "SystemStructures.h"
#include "System.h"


/* ============================================================================= */
/*  GPU container that stores three diagonal entries of ? at every vertex        */
/* ============================================================================= */
struct LambdaField {
    
    thrust::device_vector<CVec3> e_h, e_R, e_phi; // 3x1 basis at v_i     

    thrust::device_vector<double> lam_rr;   ///< radial component   lambda_rr
    thrust::device_vector<double> lam_pp;   ///< circumf. component lambda_ff
    thrust::device_vector<double> lam_ss;   ///< thickness          lambda_ss
    
    thrust::device_vector<Mat_3x3> lam_alpha; // full lambda tensor 3x3
    
    thrust::device_vector<double> rho;      ///< normalized radius at vertex i 
    void resize(std::size_t N)
    {    
        rho.resize(N);
        lam_rr.resize(N);
        lam_pp.resize(N);
        lam_ss.resize(N);
        e_h.resize(N);// e_h_y.resize(N); e_h_z.resize(N);
        e_R.resize(N);// e_R_y.resize(N); e_R_z.resize(N);
        e_phi.resize(N);// e_phi_y.resize(N); e_phi_z.resize(N);
        lam_alpha.resize(N);
    }
};

/* ============================================================================= */
/*  Public interface (all runs on GPU; wrapper functions are in this namespace)  */
/* ============================================================================= */
namespace StrainTensorGPU { 

    /** Build the basis vectors and lambda field at the current pseudo-time fraction tFrac (T/Tf) at each vertex. */
    void buildVertexLambda(GeneralParams& gp,
                           CoordInfoVecs& coord,
                           LambdaField&         field,
                           double               tFrac);

    /** Update linear-spring rest lengths with edge-wise ?-projection
        (plain diagonal average). */
    void updateEdgeRestLengths(CoordInfoVecs&      coord,
                               GeneralParams& gp,
                               const LambdaField&        field,
                               LinearSpringInfoVecs&     lsInfo, 
                               int targetLayer);


    //--------------------------------------------------------------
    //  Preferred dihedral angles  ?0  for every bending triangle
    //--------------------------------------------------------------
    void updatePreferredAngles(                     //  <<<  HOST  >>>
            BendingTriangleInfoVecs&  btiv,   // t2e*, theta0
            const CoordInfoVecs& coord,     // e2n*, edgeLen
            const LambdaField& field,
            const GeneralParams& gp,
            const LinearSpringInfoVecs& lsInfo);     // gp.thickness

} // namespace StrainTensorGPU


