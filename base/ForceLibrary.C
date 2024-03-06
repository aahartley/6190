#include "ForceLibrary.h"
#include <cstdlib>
#include <iostream>

using namespace pba;




void AccumulatingForce::compute( DynamicalState& s, const double dt )
{
#pragma omp parallel for
   for( size_t i=0;i<s->nb();i++ )
   {
      s->set_accel(i,Vector(0,0,0));
   }

   for( size_t i=0;i<forces.size();i++ )
   {
      forces[i]->compute(s,dt);
   }
}


// void AccumulatingForce::compute( RigidBodyState& s, const double dt )
// {
// #pragma omp parallel for
//    for( size_t i=0;i<s->nb();i++ )
//    {
//       s->set_accel(i,Vector(0,0,0));
//    }

//    for( size_t i=0;i<forces.size();i++ )
//    {
//       forces[i]->compute(s,dt);
//    }
// }


// void AccumulatingForce::compute( SoftBodyState& s, const double dt )
// {
// #pragma omp parallel for
//    for( size_t i=0;i<s->nb();i++ )
//    {
//       s->set_accel(i,Vector(0,0,0));
//    }

//    for( size_t i=0;i<forces.size();i++ )
//    {
//       forces[i]->compute(s,dt);
//    }
// }


// void AccumulatingForce::compute( SPHState& s, const double dt )
// {
// #pragma omp parallel for
//    for( size_t i=0;i<s->nb();i++ )
//    {
//       s->set_accel(i,Vector(0,0,0));
//    }

//    for( size_t i=0;i<forces.size();i++ )
//    {
//       forces[i]->compute(s,dt);
//    }
// }








void AccumulatingForce::add( Force& f )
{
   forces.push_back(f);
}


pba::Force pba::CreateAccumulatingForce()
{
   return pba::Force( new pba::AccumulatingForce() );
}



void AccumulatingGravityForce::compute( pba::DynamicalState& s, const double dt )
{
#pragma omp parallel for
   for( size_t i=0;i<s->nb();i++ )
   {
      s->set_accel(i, s->accel(i) + (G)); 
   }
}


pba::Force pba::CreateAccumulatingGravityForce( const Vector& G )
{
   return pba::Force( new pba::AccumulatingGravityForce(G) );
}