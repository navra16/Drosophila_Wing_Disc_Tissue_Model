//#ifndef VOLUMESPRINGS_H_
//#define VOLUMESPRINGS_H_
//
//#include "SystemStructures.h"
//
///**
// * ComputeVolumeSprings - L&L-compliant LOCAL volume conservation
// * 
// * Instead of a single global constraint E = k_v*(V_total - V0_total)^2,
// * this applies per-prism volume conservation:
// * 
// *     E = sum_p  k_v * (V_p - V0_p)^2
// * 
// * The force on node n is:
// *     F_n = sum_{p containing n}  -2*k_v*(V_p - V0_p) * dV_p/dr_n
// * 
// * Each prism has its own prefactor, so forces on node n depend only 
// * on the volumes of prisms that contain n — purely local coupling,
// * consistent with L&L §2-4.
// * 
// * Requires: prismInfoVecs.eq_prism_volume[p] to be pre-computed
// * (call ComputeEquilibriumPrismVolumes once during initialization).
// */
//void ComputeVolumeSprings(
//    CoordInfoVecs& coordInfoVecs,
//    LinearSpringInfoVecs& linearSpringInfoVecs,
//    CapsidInfoVecs& capsidInfoVecs,
//    GeneralParams& generalParams,
//    PrismInfoVecs& prismInfoVecs);
//
///**
// * ComputeEquilibriumPrismVolumes - Store per-prism rest volumes
// * 
// * Call ONCE after initial relaxation (when the mesh is in its 
// * stress-free reference state) to record V0_p for each prism.
// * These are stored in prismInfoVecs.eq_prism_volume.
// */
//void ComputeEquilibriumPrismVolumes(
//    CoordInfoVecs& coordInfoVecs,
//    GeneralParams& generalParams,
//    PrismInfoVecs& prismInfoVecs);
//
//#endif // VOLUMESPRINGS_H_


// Nav commented out the bottom. The bottom code is for a global volume conservation. It works. It just isnt right for us. 03/08/26
#ifndef VOLUMESPRINGS_H_
#define VOLUMESPRINGS_H_

#include "SystemStructures.h"

/**
 * ComputeVolumeSprings - Apply volume-conserving forces to all nodes
 * 
 * Implements the force: F_n = -dE_vol/dr_n = -2*k_v*(V - V0)*dV/dr_n
 * 
 * The gradient dV/dr_n is computed by summing contributions from all
 * prisms that contain node n, decomposed into tetrahedra.
 *
 * Forces are accumulated via atomicAdd (matching LinearSprings pattern)
 * to avoid overwriting forces from other springs.
 */
void ComputeVolumeSprings(
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    CapsidInfoVecs& capsidInfoVecs,
    GeneralParams& generalParams,
    PrismInfoVecs& prismInfoVecs);

// Note: The VolumeSpringAtomicFunctor is defined internally in VolumeSprings.cu.
// The old VolumeSpringPrismFunctor (thrust::transform-based) has been removed 
// because it caused a force overwrite bug when combined with LinearSprings' 
// atomicAdd pattern.

#endif // VOLUMESPRINGS_H_
