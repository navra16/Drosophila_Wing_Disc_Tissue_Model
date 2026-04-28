//// VolumeSprings.cu - PER-PRISM LOCAL volume conservation
////
//// Replaces the global volume constraint with per-prism constraints:
////     E = S_p k_v * (V_p - V0_p)²
////
//// This is L&L-compliant: forces on each node depend only on the
//// volumes of prisms that contain that node (local coupling).
////
//// Performance: O(P) work — iterates over prisms, not (nodes × prisms).
//// Each prism computes its own volume and atomicAdd's forces to its 6 nodes.
//
//#include "System.h"
//#include "SystemStructures.h"
//#include "VolumeSprings.h"
//
//#include <thrust/for_each.h>
//#include <thrust/iterator/counting_iterator.h>
//#include <thrust/device_vector.h>
//#include <cmath>
//#include <iostream>
//
//// ============================================================================
//// GPU kernel: per-prism volume spring forces
////
//// For each prism p:
////   1. Compute V_p from its 6 nodes (3-tet decomposition)
////   2. Compute prefactor_p = -2*k_v*(V_p - V0_p)
////   3. Compute dV_p/dr_n for all 6 nodes
////   4. atomicAdd forces to those 6 nodes
//// ============================================================================
//
//struct PerPrismVolumeSpringFunctor
//{
//    double kv;      // volume spring constant
//    int num_nodes;
//
//    // Prism connectivity (read-only)
//    const int* P1;
//    const int* P2;
//    const int* P3;
//    const int* P4;
//    const int* P5;
//    const int* P6;
//
//    // Per-prism equilibrium volumes (read-only)
//    const double* V0;
//
//    // Node coordinates (read-only)
//    const double* X;
//    const double* Y;
//    const double* Z;
//
//    // Node forces (write via atomicAdd)
//    double* nodeForceX;
//    double* nodeForceY;
//    double* nodeForceZ;
//
//    // Output: per-prism energy and current volume (for diagnostics)
//    double* prism_energy;    // can be nullptr if not needed
//    double* prism_volume;    // can be nullptr if not needed
//
//    __host__ __device__
//    PerPrismVolumeSpringFunctor(
//        double _kv, int _num_nodes,
//        const int* _P1, const int* _P2, const int* _P3,
//        const int* _P4, const int* _P5, const int* _P6,
//        const double* _V0,
//        const double* _X, const double* _Y, const double* _Z,
//        double* _fX, double* _fY, double* _fZ,
//        double* _pe, double* _pv)
//        : kv(_kv), num_nodes(_num_nodes),
//          P1(_P1), P2(_P2), P3(_P3), P4(_P4), P5(_P5), P6(_P6),
//          V0(_V0), X(_X), Y(_Y), Z(_Z),
//          nodeForceX(_fX), nodeForceY(_fY), nodeForceZ(_fZ),
//          prism_energy(_pe), prism_volume(_pv)
//    {}
//
//    // Compute 6*signed volume of tet (i,j,k,l) and gradients w.r.t. all 4 vertices
//    __device__ __forceinline__
//    void tet_vol_and_grad(
//        int i, int j, int k, int l,
//        double& sixV,
//        double* gx, double* gy, double* gz) const  // gx[4], gy[4], gz[4] for i,j,k,l
//    {
//        const double ux = X[i] - X[l], uy = Y[i] - Y[l], uz = Z[i] - Z[l];
//        const double vx = X[j] - X[l], vy = Y[j] - Y[l], vz = Z[j] - Z[l];
//        const double wx = X[k] - X[l], wy = Y[k] - Y[l], wz = Z[k] - Z[l];
//
//        // v × w  (gradient w.r.t. vertex i)
//        gx[0] = vy * wz - vz * wy;
//        gy[0] = vz * wx - vx * wz;
//        gz[0] = vx * wy - vy * wx;
//
//        // w × u  (gradient w.r.t. vertex j)
//        gx[1] = wy * uz - wz * uy;
//        gy[1] = wz * ux - wx * uz;
//        gz[1] = wx * uy - wy * ux;
//
//        // u × v  (gradient w.r.t. vertex k)
//        gx[2] = uy * vz - uz * vy;
//        gy[2] = uz * vx - ux * vz;
//        gz[2] = ux * vy - uy * vx;
//
//        // -(gi + gj + gk)  (gradient w.r.t. vertex l)
//        gx[3] = -(gx[0] + gx[1] + gx[2]);
//        gy[3] = -(gy[0] + gy[1] + gy[2]);
//        gz[3] = -(gz[0] + gz[1] + gz[2]);
//
//        sixV = ux * gx[0] + uy * gy[0] + uz * gz[0];
//    }
//
//    __device__
//    void operator()(int prism_id) const
//    {
//        const int a = P1[prism_id];
//        const int b = P2[prism_id];
//        const int c = P3[prism_id];
//        const int A = P4[prism_id];
//        const int B = P5[prism_id];
//        const int C = P6[prism_id];
//
//        // Validate node indices
//        if (a < 0 || b < 0 || c < 0 || A < 0 || B < 0 || C < 0) return;
//        if (a >= num_nodes || b >= num_nodes || c >= num_nodes ||
//            A >= num_nodes || B >= num_nodes || C >= num_nodes) return;
//
//        // 3-tet decomposition (same as VolumeComp.cu):
//        //   Tet 1: (b, c, A, a)
//        //   Tet 2: (b, A, C, a)
//        //   Tet 3: (A, B, C, a)
//        //
//        // Each tet has 4 vertices. Prism has 6 unique vertices {a,b,c,A,B,C}.
//        // We accumulate dV/dr for all 6 nodes across all 3 tets.
//
//        double grad_x[6] = {0}, grad_y[6] = {0}, grad_z[6] = {0};
//        double total_sixV = 0.0;
//
//        // --- Tet 1: (b, c, A, a) --- vertices: i=b, j=c, k=A, l=a
//        {
//            double sixV;
//            double gx[4], gy[4], gz[4];
//            tet_vol_and_grad(b, c, A, a, sixV, gx, gy, gz);
//            total_sixV += sixV;
//            // Map: [0]=b?1, [1]=c?2, [2]=A?3, [3]=a?0
//            grad_x[1] += gx[0]; grad_y[1] += gy[0]; grad_z[1] += gz[0]; // b
//            grad_x[2] += gx[1]; grad_y[2] += gy[1]; grad_z[2] += gz[1]; // c
//            grad_x[3] += gx[2]; grad_y[3] += gy[2]; grad_z[3] += gz[2]; // A
//            grad_x[0] += gx[3]; grad_y[0] += gy[3]; grad_z[0] += gz[3]; // a
//        }
//
//        // --- Tet 2: (b, A, C, a) --- vertices: i=b, j=A, k=C, l=a
//        {
//            double sixV;
//            double gx[4], gy[4], gz[4];
//            tet_vol_and_grad(b, A, C, a, sixV, gx, gy, gz);
//            total_sixV += sixV;
//            // Map: [0]=b?1, [1]=A?3, [2]=C?5, [3]=a?0
//            grad_x[1] += gx[0]; grad_y[1] += gy[0]; grad_z[1] += gz[0]; // b
//            grad_x[3] += gx[1]; grad_y[3] += gy[1]; grad_z[3] += gz[1]; // A
//            grad_x[5] += gx[2]; grad_y[5] += gy[2]; grad_z[5] += gz[2]; // C
//            grad_x[0] += gx[3]; grad_y[0] += gy[3]; grad_z[0] += gz[3]; // a
//        }
//
//        // --- Tet 3: (A, B, C, a) --- vertices: i=A, j=B, k=C, l=a
//        {
//            double sixV;
//            double gx[4], gy[4], gz[4];
//            tet_vol_and_grad(A, B, C, a, sixV, gx, gy, gz);
//            total_sixV += sixV;
//            // Map: [0]=A?3, [1]=B?4, [2]=C?5, [3]=a?0
//            grad_x[3] += gx[0]; grad_y[3] += gy[0]; grad_z[3] += gz[0]; // A
//            grad_x[4] += gx[1]; grad_y[4] += gy[1]; grad_z[4] += gz[1]; // B
//            grad_x[5] += gx[2]; grad_y[5] += gy[2]; grad_z[5] += gz[2]; // C
//            grad_x[0] += gx[3]; grad_y[0] += gy[3]; grad_z[0] += gz[3]; // a
//        }
//
//        // Current prism volume
//        double Vp = total_sixV / 6.0;
//
//        // Per-prism equilibrium volume
//        double V0p = V0[prism_id];
//
//        // Skip degenerate prisms
//        if (fabs(V0p) < 1e-12) return;
//
//        // Per-prism volume difference
//        double dV = Vp - V0p;
//
//        // Per-prism prefactor (THIS IS THE KEY DIFFERENCE from global)
//        double prefactor = -2.0 * kv * dV;
//
//        // Clamp to prevent explosion from inverted prisms
//        double max_pf = 1e4 * kv * fabs(V0p);
//        if (fabs(prefactor) > max_pf) {
//            prefactor = (prefactor > 0) ? max_pf : -max_pf;
//        }
//
//        // Scale gradients from d(6V)/dr to dV/dr
//        const double inv6 = 1.0 / 6.0;
//
//        // Apply forces to the 6 prism nodes
//        int nodes[6] = {a, b, c, A, B, C};
//        for (int v = 0; v < 6; ++v) {
//            double fx = prefactor * inv6 * grad_x[v];
//            double fy = prefactor * inv6 * grad_y[v];
//            double fz = prefactor * inv6 * grad_z[v];
//
//            if (fx != 0.0 || fy != 0.0 || fz != 0.0) {
//                atomicAdd(&nodeForceX[nodes[v]], fx);
//                atomicAdd(&nodeForceY[nodes[v]], fy);
//                atomicAdd(&nodeForceZ[nodes[v]], fz);
//            }
//        }
//
//        // Store diagnostics if requested
//        if (prism_energy) prism_energy[prism_id] = kv * dV * dV;
//        if (prism_volume) prism_volume[prism_id] = Vp;
//    }
//};
//
//
//// ============================================================================
//// ComputeEquilibriumPrismVolumes - call once after initial relaxation
//// ============================================================================
//
//void ComputeEquilibriumPrismVolumes(
//    CoordInfoVecs& coordInfoVecs,
//    GeneralParams& generalParams,
//    PrismInfoVecs& prismInfoVecs)
//{
//    const int P = prismInfoVecs.num_prisms;
//    if (P <= 0) return;
//
//    const int N = (int)coordInfoVecs.nodeLocX.size();
//
//    // Allocate storage for equilibrium volumes
//    prismInfoVecs.eq_prism_volume.resize(P);
//
//    // Temporary device vector for computed volumes
//    thrust::device_vector<double> d_volumes(P, 0.0);
//
//    // Use the same functor with kv=0 (no forces, just compute volumes)
//    // Or better: a simple volume-only kernel
//    thrust::host_vector<int> hP1 = prismInfoVecs.P1;
//    thrust::host_vector<int> hP2 = prismInfoVecs.P2;
//    thrust::host_vector<int> hP3 = prismInfoVecs.P3;
//    thrust::host_vector<int> hP4 = prismInfoVecs.P4;
//    thrust::host_vector<int> hP5 = prismInfoVecs.P5;
//    thrust::host_vector<int> hP6 = prismInfoVecs.P6;
//
//    // Copy coordinates to host for initialization
//    thrust::host_vector<double> hX = coordInfoVecs.nodeLocX;
//    thrust::host_vector<double> hY = coordInfoVecs.nodeLocY;
//    thrust::host_vector<double> hZ = coordInfoVecs.nodeLocZ;
//
//    thrust::host_vector<double> h_eq_vol(P);
//
//    double total_vol = 0.0;
//    int negative_count = 0;
//
//    for (int p = 0; p < P; ++p) {
//        int a = hP1[p], b = hP2[p], c = hP3[p];
//        int A = hP4[p], B = hP5[p], C = hP6[p];
//
//        if (a < 0 || b < 0 || c < 0 || A < 0 || B < 0 || C < 0 ||
//            a >= N || b >= N || c >= N || A >= N || B >= N || C >= N) {
//            h_eq_vol[p] = 0.0;
//            continue;
//        }
//
//        // 3-tet decomposition (matches VolumeComp.cu)
//        auto sixV = [&](int i, int j, int k, int l) -> double {
//            double ux = hX[i]-hX[l], uy = hY[i]-hY[l], uz = hZ[i]-hZ[l];
//            double vx = hX[j]-hX[l], vy = hY[j]-hY[l], vz = hZ[j]-hZ[l];
//            double wx = hX[k]-hX[l], wy = hY[k]-hY[l], wz = hZ[k]-hZ[l];
//            return ux*(vy*wz-vz*wy) + uy*(vz*wx-vx*wz) + uz*(vx*wy-vy*wx);
//        };
//
//        double s1 = sixV(b, c, A, a);
//        double s2 = sixV(b, A, C, a);
//        double s3 = sixV(A, B, C, a);
//        double vol = (s1 + s2 + s3) / 6.0;
//
//        h_eq_vol[p] = vol;
//        total_vol += vol;
//        if (vol < 0) negative_count++;
//    }
//
//    prismInfoVecs.eq_prism_volume = h_eq_vol;
//
//    std::cout << "Computed per-prism equilibrium volumes:" << std::endl;
//    std::cout << "  Prisms: " << P << std::endl;
//    std::cout << "  Total volume: " << total_vol << std::endl;
//    std::cout << "  Negative prisms: " << negative_count << std::endl;
//
//    // Also store global equilibrium volume for backward compatibility
//    generalParams.eq_total_volume = total_vol;
//}
//
//
//// ============================================================================
//// ComputeVolumeSprings - main entry point (per-prism version)
//// ============================================================================
//
//void ComputeVolumeSprings(
//    CoordInfoVecs& coordInfoVecs,
//    LinearSpringInfoVecs& linearSpringInfoVecs,
//    CapsidInfoVecs& capsidInfoVecs,
//    GeneralParams& generalParams,
//    PrismInfoVecs& prismInfoVecs)
//{
//    const int P = prismInfoVecs.num_prisms;
//    if (P <= 0) return;
//
//    const double kv = generalParams.volume_spring_constant;
//    if (kv == 0.0) return;
//
//    // Check that per-prism equilibrium volumes have been computed
//    if ((int)prismInfoVecs.eq_prism_volume.size() != P) {
//        std::cerr << "ERROR: eq_prism_volume not initialized. "
//                  << "Call ComputeEquilibriumPrismVolumes() first." << std::endl;
//        return;
//    }
//
//    const int Nnodes = (int)coordInfoVecs.nodeLocX.size();
//
//    // Iterate over PRISMS (not nodes) — O(P) work instead of O(N*P)
//    PerPrismVolumeSpringFunctor functor(
//        kv, Nnodes,
//        thrust::raw_pointer_cast(prismInfoVecs.P1.data()),
//        thrust::raw_pointer_cast(prismInfoVecs.P2.data()),
//        thrust::raw_pointer_cast(prismInfoVecs.P3.data()),
//        thrust::raw_pointer_cast(prismInfoVecs.P4.data()),
//        thrust::raw_pointer_cast(prismInfoVecs.P5.data()),
//        thrust::raw_pointer_cast(prismInfoVecs.P6.data()),
//        thrust::raw_pointer_cast(prismInfoVecs.eq_prism_volume.data()),
//        thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
//        thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
//        thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()),
//        thrust::raw_pointer_cast(coordInfoVecs.nodeForceX.data()),
//        thrust::raw_pointer_cast(coordInfoVecs.nodeForceY.data()),
//        thrust::raw_pointer_cast(coordInfoVecs.nodeForceZ.data()),
//        nullptr,  // prism_energy (set to device ptr for diagnostics)
//        nullptr); // prism_volume
//
//    thrust::for_each(
//        thrust::device,
//        thrust::counting_iterator<int>(0),
//        thrust::counting_iterator<int>(P),
//        functor);
//
//    // Compute total volume and energy on host for diagnostics
//    // (reuse ComputeVolume for global tracking — it's cheap)
//}
//

// The code below was commented out by Nav. It works for a global volume conservation but that's just not what we need right now. 03/08/26
#include "System.h"
#include "SystemStructures.h"
#include "VolumeSprings.h"

#include <thrust/for_each.h>
#include <thrust/iterator/counting_iterator.h>
#include <cmath>
#include <iostream>

// ============================================================================
// GPU kernel: compute volume spring force for each node using atomicAdd
// ============================================================================

struct VolumeSpringAtomicFunctor
{
    double prefactor;   // = -2 * k_v * (Omega - Omega_0)
    int num_prisms;
    int num_nodes;

    // Prism connectivity (read-only)
    const int* P1;
    const int* P2;
    const int* P3;
    const int* P4;
    const int* P5;
    const int* P6;

    // Node coordinates (read-only)
    const double* X;
    const double* Y;
    const double* Z;

    // Node forces (write via atomicAdd - same as LinearSprings!)
    double* nodeForceX;
    double* nodeForceY;
    double* nodeForceZ;

    __host__ __device__
    VolumeSpringAtomicFunctor(
        double _prefactor, int _num_prisms, int _num_nodes,
        const int* _P1, const int* _P2, const int* _P3,
        const int* _P4, const int* _P5, const int* _P6,
        const double* _X, const double* _Y, const double* _Z,
        double* _fX, double* _fY, double* _fZ)
        : prefactor(_prefactor), num_prisms(_num_prisms), num_nodes(_num_nodes),
          P1(_P1), P2(_P2), P3(_P3), P4(_P4), P5(_P5), P6(_P6),
          X(_X), Y(_Y), Z(_Z),
          nodeForceX(_fX), nodeForceY(_fY), nodeForceZ(_fZ)
    {}

    // Compute gradient of tet volume w.r.t. one vertex and accumulate
    __device__ __forceinline__
    void accumulate_tet_grad(
        int i, int j, int k, int l,
        int node_id,
        double& gx, double& gy, double& gz) const
    {
        if (i < 0 || j < 0 || k < 0 || l < 0) return;
        if (i >= num_nodes || j >= num_nodes || k >= num_nodes || l >= num_nodes) return;

        const double uix = X[i] - X[l];
        const double uiy = Y[i] - Y[l];
        const double uiz = Z[i] - Z[l];

        const double vjx = X[j] - X[l];
        const double vjy = Y[j] - Y[l];
        const double vjz = Z[j] - Z[l];

        const double wkx = X[k] - X[l];
        const double wky = Y[k] - Y[l];
        const double wkz = Z[k] - Z[l];

        // d(6V)/dr_i = v x w
        const double gix = vjy * wkz - vjz * wky;
        const double giy = vjz * wkx - vjx * wkz;
        const double giz = vjx * wky - vjy * wkx;

        // d(6V)/dr_j = w x u
        const double gjx = wky * uiz - wkz * uiy;
        const double gjy = wkz * uix - wkx * uiz;
        const double gjz = wkx * uiy - wky * uix;

        // d(6V)/dr_k = u x v
        const double gkx = uiy * vjz - uiz * vjy;
        const double gky = uiz * vjx - uix * vjz;
        const double gkz = uix * vjy - uiy * vjx;

        // d(6V)/dr_l = -(g_i + g_j + g_k)
        const double glx = -(gix + gjx + gkx);
        const double gly = -(giy + gjy + gky);
        const double glz = -(giz + gjz + gkz);

        const double inv6 = 1.0 / 6.0;

        if (node_id == i) { gx += inv6 * gix; gy += inv6 * giy; gz += inv6 * giz; }
        if (node_id == j) { gx += inv6 * gjx; gy += inv6 * gjy; gz += inv6 * gjz; }
        if (node_id == k) { gx += inv6 * gkx; gy += inv6 * gky; gz += inv6 * gkz; }
        if (node_id == l) { gx += inv6 * glx; gy += inv6 * gly; gz += inv6 * glz; }
    }

    __device__
    void operator()(int node_id) const
    {
        if (node_id >= num_nodes) return;

        double gx = 0.0, gy = 0.0, gz = 0.0;

        for (int p = 0; p < num_prisms; ++p) {
            const int a = P1[p];
            const int b = P2[p];
            const int c = P3[p];
            const int A = P4[p];
            const int B = P5[p];
            const int C = P6[p];

            // Skip if this node isn't in this prism
            if (node_id != a && node_id != b && node_id != c &&
                node_id != A && node_id != B && node_id != C)
                continue;

            // 3-tet decomposition (matches VolumeComp.cu exactly)
            accumulate_tet_grad(b, c, A, a, node_id, gx, gy, gz);
            accumulate_tet_grad(b, A, C, a, node_id, gx, gy, gz);
            accumulate_tet_grad(A, B, C, a, node_id, gx, gy, gz);
        }

        // Force = prefactor * gradient
        double fx = prefactor * gx;
        double fy = prefactor * gy;
        double fz = prefactor * gz;

        // CRITICAL FIX: Use atomicAdd (same as LinearSprings!)
        // This ACCUMULATES onto existing forces instead of overwriting them
        if (fx != 0.0 || fy != 0.0 || fz != 0.0) {
            atomicAdd(&nodeForceX[node_id], fx);
            atomicAdd(&nodeForceY[node_id], fy);
            atomicAdd(&nodeForceZ[node_id], fz);
        }
    }
};

// ============================================================================
// ComputeVolumeSprings - main entry point
// ============================================================================

void ComputeVolumeSprings(
    CoordInfoVecs& coordInfoVecs,
    LinearSpringInfoVecs& linearSpringInfoVecs,
    CapsidInfoVecs& capsidInfoVecs,
    GeneralParams& generalParams,
    PrismInfoVecs& prismInfoVecs)
{
    const int P = prismInfoVecs.num_prisms; 
    if (P <= 0) {
        return;
    }

    const double kv = generalParams.volume_spring_constant;
    if (kv == 0.0) {
        return;
    }

    const double Omega_s = generalParams.current_total_volume;
    const double Omega0 = generalParams.eq_total_volume;
    
    // Safeguard: NaN/Inf check
    if (std::isnan(Omega_s) || std::isinf(Omega_s)) {
        std::cout << "WARNING: Volume is NaN/Inf (" << Omega_s << "), skipping volume springs." << std::endl;
        return;
    }
    
    // Safeguard: Negative volume (mesh inversion)
    if (Omega_s < 0.0) {
        std::cout << "WARNING: Negative volume (" << Omega_s 
                  << "). Skipping volume springs." << std::endl;
        return;
    }
    
    double volume_diff = Omega_s - Omega0;
    double volume_ratio = Omega_s / Omega0;
    double prefactor;
    
    // Clamp extreme volume changes
    if (volume_ratio < 0.5 || volume_ratio > 2.0) {
        std::cout << "WARNING: Extreme volume change (ratio=" << volume_ratio 
                  << "). Clamping." << std::endl;
        double max_diff = 0.5 * std::fabs(Omega0);
        double clamped_diff = std::max(-max_diff, std::min(max_diff, volume_diff));
        prefactor = -2.0 * kv * clamped_diff;
    } else {
        prefactor = -2.0 * kv * volume_diff;
    }
    
    // Safeguard: NaN/Inf prefactor
    if (std::isnan(prefactor) || std::isinf(prefactor)) {
        std::cout << "WARNING: Volume spring prefactor is NaN/Inf. Skipping." << std::endl;
        return;
    }
    
    // Clamp magnitude
    const double max_prefactor = 1e6;
    if (std::fabs(prefactor) > max_prefactor) {
        prefactor = (prefactor > 0) ? max_prefactor : -max_prefactor;
    }
    
    generalParams.volume_energy = kv * volume_diff * volume_diff;
    const int Nnodes = (int)coordInfoVecs.nodeLocX.size();

    // Use for_each with atomicAdd (NOT thrust::transform!)
    VolumeSpringAtomicFunctor functor(
        prefactor, P, Nnodes,
        thrust::raw_pointer_cast(prismInfoVecs.P1.data()),
        thrust::raw_pointer_cast(prismInfoVecs.P2.data()),
        thrust::raw_pointer_cast(prismInfoVecs.P3.data()),
        thrust::raw_pointer_cast(prismInfoVecs.P4.data()),
        thrust::raw_pointer_cast(prismInfoVecs.P5.data()),
        thrust::raw_pointer_cast(prismInfoVecs.P6.data()),
        thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
        thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
        thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()),
        thrust::raw_pointer_cast(coordInfoVecs.nodeForceX.data()),
        thrust::raw_pointer_cast(coordInfoVecs.nodeForceY.data()),
        thrust::raw_pointer_cast(coordInfoVecs.nodeForceZ.data()));

    // for_each calls functor(node_id) for each node
    // functor uses atomicAdd to accumulate forces
    thrust::for_each(
        thrust::device,
        thrust::counting_iterator<int>(0),
        thrust::counting_iterator<int>(Nnodes),
        functor);
}
