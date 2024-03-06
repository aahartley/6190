  
  //-------------------------------------------------------
 //
 //  CollisionGeometryLibrary.h
 //
 // Create collision geometry.
 //
 //  Copyright (c) 2017 Jerry Tessendorf
 //
 //
 //--------------------------------------------------------
 #ifndef __PBA_COLLISIONGEOMETRYLIBRARY_H__
 #define __PBA_COLLISIONGEOMETRYLIBRARY_H__
  
 #include "Vector.h"
 #include "Color.h"
 #include "CollisionSurface.h"
 #include <string>
 #include <vector>
  
  
  
  
 using namespace std;
  
 namespace pba{
  
  
  
 pba::CollisionSurface GenerateCollisionCube();
  
  
 pba::CollisionSurface ObjCollisionSurface( const string filename );
  
  
 };
  
  
 #endif
  