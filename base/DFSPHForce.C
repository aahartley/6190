//-------------------------------------------------------
//
//  DFSPHForce.C
//
//  A collection of force implementations for DFSPH
//
//  Copyright (c) 2019 Jerry Tessendorf
//
//
//--------------------------------------------------------

#include "DFSPHForce.h"
#include "OccupancyVolume.h"
#include <cmath>
#include <iostream>

namespace pba
{

DFSPHForce::DFSPHForce() :
  dynamic_viscosity(0.01f)
{

}

void DFSPHForce::compute(SPHState& s, const double dt )
{
#pragma omp parallel for
   for( size_t p=0;p<s->nb();p++)
   {
      float density = s->get_float_attr( "density", p );
      float mass = s->mass(p);
      float kinematic_viscosity = dynamic_viscosity / density;
      float h = s->get_radius();
      const int dimensions = 3; //make part of solver or state?
      const Vector& P = s->pos(p);
      const Vector& V = s->vel(p);
      size_t pindex = s->index(P);
      size_t i,j,k;
      s->anti_index(pindex, i,j,k);
      const std::vector<size_t>& cells = s->cell_contents(i,j,k);

      Vector laplacian{};
      Vector viscosity_force{}; 

// PARTICLES IN SAME CELL 
      for(size_t a=0;a<cells.size();a++)
      {
         //pid != p???
        size_t pid = cells[a];
        Vector v_ij = V - s->vel(pid);
        Vector x_ij = P - s->pos(pid);
        float pid_density = s->get_float_attr("density", pid);
        float pid_mass = s->mass(pid);

        laplacian += ( (pid_mass / pid_density) * 
                        ( (v_ij * x_ij) / ((x_ij.magnitude() * x_ij.magnitude()) + 0.01 * (h*h)) ) ) * s->grad_weight(pid, P); 
      

      }

// PARTICLES IN NEIGHBORING CELLS 
      std::vector<size_t> neighbors;
      s->neighbor_cells(i,j,k,neighbors);
      for(size_t a=0;a<neighbors.size();a++)
      {
        size_t pid = neighbors[a];
	    Vector v_ij = V - s->vel(pid);
        Vector x_ij = P - s->pos(pid);
        float pid_density = s->get_float_attr("density", pid);
        float pid_mass = s->mass(pid);

        laplacian += ( (pid_mass / pid_density) *
                         ( (v_ij * x_ij) / ((x_ij.magnitude() * x_ij.magnitude()) + 0.01 * (h*h)) ) ) * s->grad_weight(pid, P);
	
	 
      }
      laplacian = (2*(dimensions+2) * laplacian);
      viscosity_force = (mass * kinematic_viscosity) * laplacian;
      //#pragma omp critical
      s->set_accel(p, s->accel(p) + (viscosity_force)/s->mass(p)); // / by mass?
   } 
}



Force CreateDFSPHForce()
{
   return pba::Force( new pba::DFSPHForce() );
}


}
