//-------------------------------------------------------
//
//  ExplicitDynamics.C
//
//  Various Dynamics Models Using Explicit Solvers
//
//  Copyright (c) 2017 Jerry Tessendorf
//
//
//--------------------------------------------------------




#include "ExplicitDynamics.h"



using namespace pba;

AdvancePosition::AdvancePosition( DynamicalState& pq ):
  PQ (pq)
  {}

void AdvancePosition::init(){}

void AdvancePosition::solve(const double dt)
{
#pragma omp parallel for
   for( size_t i=0;i<PQ->nb();i++ )
   {
      PQ->set_pos( i, PQ->pos(i) + PQ->vel(i)*dt );
   }
}




AdvanceVelocity::AdvanceVelocity( DynamicalState& pq, Force& f ):
  PQ (pq),
  force (f)
  {}

void AdvanceVelocity::init(){}

void AdvanceVelocity::solve(const double dt)
{
   force->compute(PQ, dt);
#pragma omp parallel for
   for( size_t i=0;i<PQ->nb();i++ )
   {
      PQ->set_vel( i, PQ->vel(i) + PQ->accel(i)*dt );
   }
}


pba::GISolver pba::CreateAdvancePosition( pba::DynamicalState& pq )
{
   return pba::GISolver( new pba::AdvancePosition(pq) );
}

pba::GISolver pba::CreateAdvanceVelocity( pba::DynamicalState& pq, pba::Force& f )
{
   return pba::GISolver( new pba::AdvanceVelocity(pq,f) );
}



AdvancePositionWithCollisions::AdvancePositionWithCollisions( DynamicalState& pq, CollisionHandler& cs ):
  PQ (pq),
  CS (cs)
  {}



void AdvancePositionWithCollisions::init(){}
void AdvancePositionWithCollisions::solve(const double dt)
{
#pragma omp parallel for
   for( size_t i=0;i<PQ->nb();i++ )
   {
      PQ->set_pos( i, PQ->pos(i) + PQ->vel(i)*dt );
   }
   CS.handle_collisions( dt, PQ );
}


pba::GISolver pba::CreateAdvancePosition( DynamicalState& pq, CollisionHandler& cs )
{
   return pba::GISolver( new pba::AdvancePositionWithCollisions(pq, cs) );

}











AdvanceRandomVelocity::AdvanceRandomVelocity( DynamicalState& pq, Force& f ):
  PQ (pq),
  force (f)
  {}

void AdvanceRandomVelocity::init(){}

void AdvanceRandomVelocity::solve(const double dt)
{
   force->compute(PQ, dt);
#pragma omp parallel for
   for( size_t i=0;i<PQ->nb();i++ )
   {
      PQ->set_vel( i, PQ->accel(i) );
   }
}


pba::GISolver pba::CreateAdvanceRandomVelocity( pba::DynamicalState& pq, pba::Force& f )
{
   return pba::GISolver( new pba::AdvanceRandomVelocity(pq,f) );
}


