#ifndef GRADIENTRELAX_H_
#define GRADIENTRELAX_H_

#include "System.h"        // for class System and Solve_Forces()
#include "NodeAdvance.h"   // for AdvancePositions(...)
  
/// Run overdamped (gradient-descent) relaxation on the system
/// until forces fall below generalParams.tol. Returns number of iterations.
int relaxUntilConverged(System& system);

#endif // GRADIENTRELAX_H_
