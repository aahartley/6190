//-------------------------------------------------------
 //
 //  ForceLibrary.h
 //
 //  A collection of force implementations
 //
 //  Copyright (c) 2017 Jerry Tessendorf
 //
 //
 //--------------------------------------------------------
#ifndef ____PBA_FORCELIBRARY_H____
 #define ____PBA_FORCELIBRARY_H____
  
 #include "Force.h"
  
 namespace pba
 {
  
  
  
 class AccumulatingForce : public ForceBase
 {
   public:
  
     AccumulatingForce() {}
  
    ~AccumulatingForce(){}
  
     void compute( DynamicalState& s, const double dt );
    //  void compute( RigidBodyState& s, const double dt );
    //  void compute( SoftBodyState& s, const double dt );
    //  void compute( SPHState& s, const double dt );
  
     void add( Force& f );
  
  
   private:
  
     std::vector<Force> forces;
  
 };
  
 Force CreateAccumulatingForce();
  
  
  
  class AccumulatingGravityForce : public pba::ForceBase
{

  public:

    AccumulatingGravityForce( const Vector& g ) : 
    G(g)
    {}

   ~AccumulatingGravityForce(){};

    void compute( pba::DynamicalState& pq, const double dt );
    // void compute( RigidBodyState& s, const double dt )
    // {
    //    DynamicalState ss = std::dynamic_pointer_cast<DynamicalStateData, RigidBodyStateData>(s);
    //    compute(ss,dt);
    // }
    // void compute( SoftBodyState& s, const double dt )
    // {
    //    DynamicalState ss = std::dynamic_pointer_cast<DynamicalStateData, SoftBodyStateData>(s);
    //    compute(ss,dt);
    // }
    // void compute( SPHState& s, const double dt )
    // {
    //    DynamicalState ss = std::dynamic_pointer_cast<DynamicalStateData, SPHStateData>(s);
    //    compute(ss,dt);
    // }



    void set_strength(const double v ){ G *= v/G.magnitude(); };
    void set_strength(const Vector& v ){ G = v; };
    const double get_strength() const { return G.magnitude(); };
    const Vector& get_strengthv() const { return G; };

  private:

    Vector G;
};

pba::Force CreateAccumulatingGravityForce( const Vector& G );
  
  
  

  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
 }
 #endif