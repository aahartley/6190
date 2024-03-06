//*******************************************************************
//
//   AABB.h
//
//  Axis Aligned Bounding Box for intersection testing
//
//  Copyright (c) 2017 Jerry Tessendorf
//
//
//*******************************************************************


#ifndef ___PBA_AABB_H____
#define ___PBA_AABB_H____


#include "Vector.h"
#include "Matrix.h"

namespace pba {

/*!
  AABB is an axis-aligned bounding box, suitable to designating rectangular regions of space.
 */
class AABB
{
  public:
    AABB(const Vector& lc, const Vector& uc);
   ~AABB();
    
    //! Test whether the input AABB intersects this AABB
    const bool intersects( const AABB& aabb ) const;
    //! Return an AABB for the region of intersection of the input AABB and this AABB.
    AABB Intersection( const AABB& aabb ) const;
    //! Return the union AABB for the input AABB and this AABB
    AABB Union( const AABB& aabb ) const;
    //! Divide this AABB into two equal AABBs along the component direction
    void split( const int component, AABB& aabb1, AABB& aabb2 ) const;
    //! Return the volume of this AABB
    const double volume() const;
    //! Test whether the input point is inside this AABB
    const bool isInside( const Vector& P ) const;
    //! Return the distance at which the input ray intersects this AABB
    const double intersect( const Vector& start, const Vector& direction ) const;
    //! Not ready to use
    const double intersect( const Vector& P, const Vector& V, const Vector& W, const Matrix& rot, const Vector& com ) const;
    //! Return the lower left corner of this AABB
    const Vector& LLC() const { return llc; }
    //! Return the upper right corner of this AABB
    const Vector& URC() const { return urc; }


  private:

    Vector llc;
    Vector urc;
};



}

#endif
