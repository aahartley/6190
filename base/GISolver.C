
#include "GISolver.h"



using namespace pba;





pba::GISolver pba::CreateGISolverSubstep( pba::GISolver& s, int nbsteps )
{
   return GISolver( new GISolverSubstep(s, nbsteps) );
}


pba::GISolver pba::CreateLeapFrogSolver( pba::GISolver& A, pba::GISolver&  B )
{
   return pba::GISolver( new pba::LeapFrogSolver(A, B) );
}

pba::GISolver pba::CreateForwardEulerSolver( pba::GISolver& A, pba::GISolver& B )
{
   return GISolver( new ForwardEulerSolver(A, B) );
}

pba::GISolver pba::CreateBackwardEulerSolver( pba::GISolver& A, pba::GISolver& B )
{
   return GISolver( new BackwardEulerSolver(A, B) );
}

