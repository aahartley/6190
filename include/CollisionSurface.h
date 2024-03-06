 //*******************************************************************
 //
 //   CollisionSurface.h
 //
 // A collection of collision geometry
 //
 //  Copyright (c) 2017 Jerry Tessendorf
 //
 //
 //*******************************************************************
#ifndef __PBA_COLLISIONSURFACE_H__
 #define __PBA_COLLISIONSURFACE_H__
  
 #include "CollisionTriangle.h"
 #include "CollisionPlane.h"
//  #include "RigidBodyState.h"
//  #include "SoftBodyState.h"
  
 #include <vector>
 #include <memory>
  
 namespace pba{
  
  
 struct CollisionData
 {
   double hit_time;
   CollisionTriangle tri;
   CollisionPlane plane;
   bool hit_tri;
   bool hit_plane;
   size_t hit_index;
   Vector XH; 
 };
  
 // The nudge parameter is used because the intersection calcularion is usually just a little bit of for rigid bodies, so the nudge is a correction factor.
  
  
 class CollisionSurfaceRaw
 {
   public:
  
     CollisionSurfaceRaw();
    ~CollisionSurfaceRaw(){}
  
     void addTriangle( const CollisionTriangle& t );
     void addPlane (const CollisionPlane& p);
     void clear() { tri_elements.clear(); plane_elements.clear();}
     size_t triangle_size() const { return tri_elements.size(); }
     size_t plane_size() const { return plane_elements.size();}
  
     bool hit( const Vector& X0, const Vector& XU, const Vector& V, const double tmax, CollisionData& t, float radius ) const;
  
    //  bool hit( const RigidBodyState& s, const size_t i, const double tmax, CollisionData& t ) const ;
  
    //  bool hit( const SoftBodyState& s, const size_t i, const double tmax, CollisionData& t ) const ;
  
  
  
     CollisionTriangle& get_triangle(size_t i) { return tri_elements[i]; }
     const CollisionTriangle& get_triangle(size_t i) const { return tri_elements[i]; }

     CollisionPlane& get_plane(size_t i) { return plane_elements[i]; }
     const CollisionPlane& get_plane(size_t i) const { return plane_elements[i]; }
  
  
     void toggle_visible(){ visible = !visible; }
     bool is_visible() const { return visible; }
     void toggle_wireframe(){ wireframe = !wireframe; }
     bool use_wireframe() const { return wireframe; }
     void toggle_points(){ points = !points; }
     bool use_points() const { return points; }
  
     void set_coeff_restitution( const double v ){ coeff_of_restitution = v; }
     const double& coeff_restitution() const { return coeff_of_restitution; }
  
     void set_coeff_sticky( const double v ){ coeff_of_sticky = v; }
     const double& coeff_sticky() const { return coeff_of_sticky; }
  
  
     void decay();
  
   private:
  
     bool visible;
     bool wireframe;
     bool points;
  
     double coeff_of_restitution;
     double coeff_of_sticky;
  
     std::vector<CollisionTriangle> tri_elements;
     std::vector<CollisionPlane> plane_elements;
 
 };
  
  
 typedef std::shared_ptr<CollisionSurfaceRaw> CollisionSurface;
  
 CollisionSurface makeCollisionSurface();
  
  
  
 }
  
 #endif