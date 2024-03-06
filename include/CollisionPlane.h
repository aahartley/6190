//*******************************************************************
 //
 //   CollisionPlane.h
 //
 // Class for a flat plane for collisions
 //
 //  Copyright (c) 2017 Jerry Tessendorf
 //
 //
 //*******************************************************************
#ifndef __PBA_COLLISIONPLANE_H__
 #define __PBA_COLLISIONPLANE_H__
  
 #include "Vector.h"
  #include <iostream>
  #include <algorithm>
  
 namespace pba{
  
 class CollisionPlane
 {
   public:
    CollisionPlane(){P0=Vector(0,0,0); normal = Vector(0,0,0);}
    CollisionPlane( const Vector& normal, const Vector& p0  );
    ~CollisionPlane(){}
    
    bool hit( const Vector& X0, const Vector& Xu, const Vector& V, const double dt, Vector& XH_cand, double& dtH_cand, float radius ) const;
  
    // // If returns true, then XH and dtH are filled with the hit point and hit time.
    // bool hit( const Vector& XS, const Vector VS, const double& dt, Vector& XH, double& dtH ) const;
    // Takes in hit data and returns reflected position and velocity
    void handle( const Vector& XS, const Vector& VS, const double& dt, const Vector& XH, const double& dtH, Vector& XR, Vector& VR, float cs, float cr ) const;  

    const Vector& getP0(){return P0;} 

   private:

     Vector P0;
     Vector normal;
  
 };
  
 }
  
 #endif