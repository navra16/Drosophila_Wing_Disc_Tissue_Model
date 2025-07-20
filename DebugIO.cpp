#include "DebugIO.h"
#include "StrainTensor.h"           // where LambdaField is defined
#include <thrust/host_vector.h>
#include <fstream>

void dumpLambdaCSV(const LambdaField& field,
                   const std::string& fileName,
                   std::size_t maxToWrite)
{
    thrust::host_vector<double> h_rr = field.lam_rr;
    thrust::host_vector<double> h_pp = field.lam_pp;
    thrust::host_vector<double> h_ss = field.lam_ss;

    std::size_t N = h_rr.size();
    if (maxToWrite && maxToWrite < N) N = maxToWrite;

    std::ofstream out(fileName);
    out << "idx,lam_rr,lam_pp,lam_ss\n";
    for (std::size_t i = 0; i < N; ++i)
        out << i << ',' << h_rr[i] << ',' << h_pp[i] << ',' << h_ss[i] << '\n';
}
