//*******************************************************************
//
//   OccupancyVolume.C
//
//
//  Occupancy Volume to track binning of particles
//
//  Copyright (c) 2019 Jerry Tessendorf
//
//
//*******************************************************************


#include "OccupancyVolume.h"
#include <iostream>

namespace pba {

OccupancyVolume::OccupancyVolume( const AABB& aabb, const double h ) :
   bounds(aabb),
   cellsize(h),
   nx(0),
   ny(0),
   nz(0)
{
   compute_size();
}
  

OccupancyVolume::~OccupancyVolume(){}

OccupancyVolume::OccupancyVolume(const OccupancyVolume& o ) :
   bounds(o.bounds),
   cellsize(o.cellsize),
   nx(o.nx),
   ny(o.ny),
   nz(o.nz),
   contents(o.contents),
   center_of_mass(o.center_of_mass),
   total_mass(o.total_mass),
   LX(o.LX)
{}




OccupancyVolume& OccupancyVolume::operator=(const OccupancyVolume& o )
{
   bounds = o.bounds;
   cellsize = o.cellsize;
   nx = o.nx;
   ny = o.ny;
   nz = o.nz;
   contents = o.contents;
   center_of_mass = o.center_of_mass;
   total_mass = o.total_mass;
   LX = o.LX;
   return *this;
}



   
size_t OccupancyVolume::nbx() const { return nx; }
size_t OccupancyVolume::nby() const { return ny; }
size_t OccupancyVolume::nbz() const { return nz; }
size_t OccupancyVolume::cell_nb(size_t i, size_t j, size_t k) const 
{
   size_t ii = index(i,j,k);
   if (ii < contents.size())
   {
      return contents[ii].size();
   }
   return 0;
}

const std::vector<size_t>& OccupancyVolume::cell_contents(size_t i, size_t j, size_t k) const
{
   size_t ii = index(i,j,k);
   if (ii < contents.size())
   {
      return contents[ii];
   }
   return contents[0];
}

Vector OccupancyVolume::cell_center_of_mass(size_t i, size_t j, size_t k) const
{
   size_t ii = index(i,j,k);
   if (ii < center_of_mass.size())
   {
      return center_of_mass[ii];
   }
   return Vector(0,0,0);	
}

double OccupancyVolume::cell_total_mass(size_t i, size_t j, size_t k) const
{
   size_t ii = index(i,j,k);
   if (ii < total_mass.size())
   {
      return total_mass[ii];
   }
   return 0.0;	
}


size_t OccupancyVolume::index( size_t i, size_t j, size_t k ) const
{
   if(i>=nx){ return nx*ny*nz; }
   if(j>=ny){ return nx*ny*nz; }
   if(k>=nz){ return nx*ny*nz; }
   return i + nx*(j+ny*k);
}

size_t OccupancyVolume::index( const Vector& P ) const
{
   Vector X = P - bounds.LLC();
   for(size_t q=0;q<3;q++)
   {
      if(X[q] < 0.0 ){ return nx*ny*nz;}
   }
   for(size_t q=0;q<3;q++)
   {
      X[q] /= LX[q];
   }
   for(size_t q=0;q<3;q++)
   {
      if( X[q] < 0.0 ){ return nx*ny*nz; }
   }
   size_t i = (size_t)(X[0]*nx);
   size_t j = (size_t)(X[1]*ny);
   size_t k = (size_t)(X[2]*nz);
   return index(i,j,k);
}



void OccupancyVolume::populate(const DynamicalStateData& state)
{
   contents.clear();
   center_of_mass.clear();
   total_mass.clear();
   contents.resize(nx*ny*nz);
   center_of_mass.resize(nx*ny*nz);
   total_mass.resize(nx*ny*nz);


   Vector zero(0,0,0);
#pragma omp parallel for
   for(size_t i=0;i<nx*ny*nz;i++)
   {
      center_of_mass[i] = zero;
      total_mass[i] = 0.0;
   }

   // Step 2: populate contents
#pragma omp parallel for
   for(size_t i=0;i<state.nb();i++)
   {
      size_t ii = index(state.pos(i));
      //const Vector& P = state.pos(i);
//std::cout << "Particle " << i << " at cell " << ii << "    " << P.X() << " " << P.Y() << " " << P.Z() << std::endl;
      if(ii < nx*ny*nz)
      {
#pragma omp critical
{
         contents[ii].push_back(i);
         center_of_mass[ii] += state.pos(i) * state.mass(i);
         total_mass[ii] += state.mass(i);
//std::cout << "Particle " << i << " pushed to cell " << ii << "\n";
}
      }
//      else
//      {
//#pragma omp critical
//{
//std::cout << "Particle " << i << " NOT PUSHED\n";
//}
//      }
   }

   // Step 3 normalize
#pragma omp parallel for
   for(size_t i=0;i<nx*ny*nz;i++)
   {
      if( total_mass[i] > 0 )
      {
         center_of_mass[i] /= total_mass[i];
      }
   }
}


void OccupancyVolume::anti_index( const size_t ind, size_t& i, size_t& j, size_t& k ) const
{
   k = ind/(nx*ny);
   j = (ind-k*nx*ny)/nx;
   i = (ind-k*nx*ny-j*nx);
}

void OccupancyVolume::neighbor_cells( size_t i, size_t j, size_t k, std::vector<size_t>& neighbors ) const
{
   neighbors.clear();
   for(size_t kk=k-1;kk<=k+1;kk++)
   {
      for(size_t jj=j-1;jj<=j+1;jj++)
      {
         for(size_t ii=i-1;ii<=i+1;ii++)
	 {
	    if( (ii!=i) || (jj!=j) || (kk!=k) )
	    {   
               size_t ind = index(ii,jj,kk);
	       if( ind < nx*ny*nz )
	       {
	          for(size_t pp=0;pp<contents[ind].size();pp++)
		  {
	             neighbors.push_back(contents[ind][pp]);
		  }
	       }
	    }
	 }
      }
   }

}

void OccupancyVolume::distant_cells( size_t i, size_t j, size_t k, std::vector<size_t>& neighbors ) const
{
   neighbors.clear();
   for(size_t p=0;p<nx*ny*nz;p++)
   {
      size_t ii,jj,kk;
      anti_index(p, ii,jj,kk);
      if( fabs(i-ii) > 1 || fabs(j-jj) > 1 || fabs(k-kk) > 1 )
      {
	 if( cell_nb(ii,jj,kk) > 0 )
	 { 
            neighbors.push_back(p);
	 }
      }
   }
}



void OccupancyVolume::compute_size()
{
   LX = bounds.URC()-bounds.LLC();
   nx = (size_t)LX.X()/cellsize + 1;
   ny = (size_t)LX.Y()/cellsize + 1;
   nz = (size_t)LX.Z()/cellsize + 1;
   std::cout << "Occupancy Volume Construction:\n";
   std::cout << "\tLLC        " << bounds.LLC().X() << " " << bounds.LLC().Y() << " " << bounds.LLC().Z() << "\n";
   std::cout << "\tURC        " << bounds.URC().X() << " " << bounds.URC().Y() << " " << bounds.URC().Z() << "\n";
   std::cout << "\tlengths    " << LX.X() << " " << LX.Y() << " " << LX.Z() << "\n";
   std::cout << "\tDimensions " << nx << " " << ny << " " << nz << "\n\n";
   contents.clear();
   center_of_mass.clear();
   total_mass.clear();
   contents.resize(nx*ny*nz);
   center_of_mass.resize(nx*ny*nz);
   total_mass.resize(nx*ny*nz);
}


void OccupancyVolume::set_cellsize(const double& v )
{
   cellsize = v;
   compute_size();
}


void OccupancyVolume::set_bounds( const Vector& llc, const Vector& urc )
{
   bounds = AABB(llc, urc);
   compute_size();
}

const Vector& OccupancyVolume::get_bounds_llc() const { return bounds.LLC(); }
const Vector& OccupancyVolume::get_bounds_urc() const { return bounds.URC(); }


void OccupancyVolume::recompute_bounds(const DynamicalStateData& state)
{
   if( state.nb() == 0 )
   {
      set_bounds( -Vector(cellsize,cellsize,cellsize), Vector(cellsize,cellsize,cellsize) );
      return;
   }

   Vector llc = state.pos(0);
   Vector urc = state.pos(0);
   for(size_t i=1;i<state.nb();i++)
   {
      const Vector& P = state.pos(i);
      if( P.X() < llc.X()){ llc[0] = P.X(); }
      if( P.X() > urc.X()){ urc[0] = P.X(); }
      if( P.Y() < llc.Y()){ llc[0] = P.Y(); }
      if( P.Y() > urc.Y()){ urc[0] = P.Y(); }
      if( P.Z() < llc.Z()){ llc[0] = P.Z(); }
      if( P.Z() > urc.Z()){ urc[0] = P.Z(); }
   }
   llc -= Vector(cellsize,cellsize,cellsize);
   urc += Vector(cellsize, cellsize, cellsize);
   set_bounds(llc,urc);
}

}

