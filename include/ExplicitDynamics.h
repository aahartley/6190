 //-------------------------------------------------------
 //
 //  ExplicitDynamics.h
 //
 //  Various Dynamics Models Using Explicit Solvers
 //
 //  Copyright (c) 2017 Jerry Tessendorf
 //
 //
 //--------------------------------------------------------
 #ifndef ____PBA_EXPLICITDYNAMICS_H____
 #define ____PBA_EXPLICITDYNAMICS_H____
  
  
 #include "DynamicalState.h"
 #include "Force.h"
 #include "GISolver.h"
 #include "CollisionHandler.h"
  
  
  
 namespace pba
 {
  
  
 //   Solvers with no collisions or constraints
  
 class AdvancePosition : public GISolverBase
 {
   public:
  
     AdvancePosition( DynamicalState& pq );
    ~AdvancePosition(){}
  
  
     void init();
     void solve(const double dt);
  
   private:
  
     DynamicalState PQ;
  
 };
  
  
 class AdvanceVelocity : public GISolverBase
 {
   public:
  
     AdvanceVelocity( DynamicalState& pq, Force& f );
    ~AdvanceVelocity(){}
  
  
     void init();
     void solve(const double dt);
  
   private:
  
     DynamicalState PQ;
     Force force;
  
 };
  
  
 GISolver CreateAdvancePosition( DynamicalState& pq );
 GISolver CreateAdvanceVelocity( DynamicalState& pq, Force& f );
  
  
 //   Solvers that involve collisions
  
 class AdvancePositionWithCollisions : public GISolverBase 
 {
   public:
  
     AdvancePositionWithCollisions( DynamicalState& pq, CollisionHandler& coll );
    ~AdvancePositionWithCollisions(){}
  
  
     void init();
     void solve(const double dt);
  
   private:
  
     DynamicalState PQ;
     CollisionHandler& CS;
  
 };
  
 GISolver CreateAdvancePosition( DynamicalState& pq, CollisionHandler& cs );
  

  
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
 class AdvanceRandomVelocity : public GISolverBase
 {
   public:
  
     AdvanceRandomVelocity( DynamicalState& pq, Force& f );
    ~AdvanceRandomVelocity(){}
  
  
     void init();
     void solve(const double dt);
  
   private:
  
     DynamicalState PQ;
     Force force;
  
 };
  
  
 GISolver CreateAdvanceRandomVelocity( DynamicalState& pq, Force& f );
  
  
 }
  
  
 #endif