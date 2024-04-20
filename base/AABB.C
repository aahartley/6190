//*******************************************************************
//
//   AABB.C
//
//  Axis Aligned Bounding Box for intersection testing
//
//  Copyright (c) 2017 Jerry Tessendorf
//
//
//*******************************************************************

#include "AABB.h"
#include <iostream>

using namespace std;

using namespace pba;


AABB::AABB(const Vector& lc, const Vector& uc) :
  llc (lc),
  urc (uc)
  {}

AABB::~AABB(){}
    
const bool AABB::intersects( const AABB& aabb ) const
{
   if( llc[0] > aabb.urc[0] ){ return false; }
   if( llc[1] > aabb.urc[1] ){ return false; }
   if( llc[2] > aabb.urc[2] ){ return false; }
   if( aabb.llc[0] > urc[0] ){ return false; }
   if( aabb.llc[1] > urc[1] ){ return false; }
   if( aabb.llc[2] > urc[2] ){ return false; }
   return true;
}

AABB AABB::Intersection( const AABB& aabb ) const
{
   Vector nllc = llc;
   Vector nurc = urc;
   if( !intersects(aabb) ){ return AABB(nllc,nllc); }
   for( int i=0;i<2;i++ )
   {
      if( llc[i] < aabb.llc[i] ){ nllc[i] = aabb.llc[i]; }
      if( urc[i] > aabb.urc[i] ){ nurc[i] = aabb.urc[i]; }
   }
   return AABB(nllc,nurc);
}

AABB AABB::Union( const AABB& aabb ) const
{
   Vector nllc = llc;
   Vector nurc = urc;
   for( int i=0;i<3;i++ )
   {
      if( nllc[i] > aabb.llc[i] ){ nllc[i] = aabb.llc[i]; }
      if( nurc[i] < aabb.urc[i] ){ nurc[i] = aabb.urc[i]; }
   }
   return AABB(nllc,nurc);
}


void AABB::split( const int component,  AABB& aabb1, AABB& aabb2 ) const
{
   Vector center = (llc+urc)*0.5;
   Vector split_plane = llc;
   split_plane[component] = center[component];
   aabb1.llc = split_plane;
   aabb1.urc = urc;
   split_plane = urc;
   split_plane[component] = center[component];
   aabb2.llc = llc;
   aabb2.urc = split_plane;
}

const double AABB::volume() const
{
   Vector diff = urc-llc;
   return diff[0]*diff[1]*diff[2];
}

const bool AABB::isInside( const Vector& P ) const
{
   for( int i=0;i<3;i++ )
   {
      if( P[i] < llc[i] ){ return false; }
      if( P[i] > urc[i] ){ return false; }
   }
   return true;
}

const double AABB::intersect( const Vector& start, const Vector& direction ) const
{
   // Smits method ala Shirley document
   double tmin, tmax, tymin, tymax, tzmin, tzmax;
   double divx = 1.0/direction[0];
   double divy = 1.0/direction[1];
   double divz = 1.0/direction[2];
   if (divx >= 0) 
   {
      tmin = (llc[0] - start[0]) * divx;
      tmax = (urc[0] - start[0]) * divx;
   }
   else 
   {
      tmin = (urc[0] - start[0]) * divx;
      tmax = (llc[0] - start[0]) * divx;
   }
   if (divy >= 0) 
   {
      tymin = (llc[1] - start[1]) * divy;
      tymax = (urc[1] - start[1]) * divy;
   }
   else 
   {
      tymin = (urc[1] - start[1]) * divy;
      tymax = (llc[1] - start[1]) * divy;
   }
   if ( (tmin > tymax) || (tymin > tmax) ) { return -1.0; }
   if (tymin > tmin) { tmin = tymin; }
   if (tymax < tmax) { tmax = tymax; }
   if (divz >= 0) 
   {
      tzmin = (llc[2] - start[2]) * divz;
      tzmax = (urc[2] - start[2]) * divz;
   }
   else 
   {
      tzmin = (urc[2] - start[2]) * divz;
      tzmax = (llc[2] - start[2]) * divz;
   }
   if ( (tmin > tzmax) || (tzmin > tmax) ) { return -1.0; }
   if (tzmin > tmin) { tmin = tzmin; }
   if( tmin < 0 || std::isnan(tmin) ){ tmin = -1.0; }
   return tmin;
}



const double AABB::intersect( const Vector& P, const Vector& V, const Vector& W, const Matrix& rot, const Vector& com ) const
{
   return -1.0;

   /*
   double ommag = W.magnitude();
   Vector om = W.unitvector();



   // first plane
   Vector normal(1,0,0);
   Vector P0(llc);
   if( ommag == 0.0 )
   {
      t = (normal*(P-P0))/(normal*V);
   }
   else
   {
      // midpoint root finding
      double t0 = 0;
      double f0 = normal*( P-P0 );
      if( f0 == 0.0 ){ t = 0.0; std::cout << "f0 zero\n"; return true; }
      double t1 = tmax;
      Matrix U = pba::rotation( om, ommag*t1 ) * rot;
      double f1 = normal*(U*(P-com) + com - V*t1 - P0 );
      if( f0*f1 > 0.0 ){ return false; }
      if( f1 == 0.0 ){ t = tmax; return true; }
      bool converged = false;
      while(!converged)
      {
         double tt = (t0+t1)*0.5;
         U = pba::rotation(om, ommag*tt ) * rot;
         double ff = normal*(U*(P-com) + com - V*tt - P0 );
         if( ff == 0.0 ){ t = tt; return true; }
         if( ff*f0 < 0.0 )
         {
            f1 = ff;
            t1 = tt;
         }
         else
         {
            f0 = ff;
            t0 = tt;
         }
         if( std::fabs(t0-t1)/tmax < 0.0001 ){ t = (t0+t1)*0.5; converged = true; }
      }
   }

   if( t*tmax <= 0.0 ){  return false; }
   if( std::fabs(t) > std::fabs(tmax) ) { return false; }

   Vector X = P-com;
   if( ommag > 0.0 )
   {
      X = pba::rotation(om, ommag*t) * X;
   }
   X += com - V*t;
   X -= P0;


   double e1X = X*e1;
   double e2X = X*e2;

   double a = (e1X*(e2*e2) - e2X*(e1*e2))/det;
   if( a < 0.0 ){ return false; }
   if( a > 1.0 ){ return false; }

   double b = (e2X*(e1*e1) - e1X*(e1*e2))/det;
   if( b < 0.0 ){ return false; }
   if( b > 1.0 ){ return false; }
   if( a+b > 1.0 ){ return false; }

   */

}


