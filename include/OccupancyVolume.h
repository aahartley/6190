//*******************************************************************
//
//   OccupancyVolume.h
//
//
//  Occupancy Volume to track binning of particles
//
//  Copyright (c) 2019 Jerry Tessendorf
//
//
//*******************************************************************


#ifndef ___PBA_OCCUPANCY_VOLUME_H____
#define ___PBA_OCCUPANCY_VOLUME_H____


#include "AABB.h"
#include "DynamicalState.h"

namespace pba {



/*!
  Occupancy volume is an axis-aligned 3D grid of cells with each cell containing:
  (1) A list of particles in that cell
  (2) The total mass of the particles in the cell
  (3) The center of mass of the particles in the cell
 */
class OccupancyVolume 
{
  public:
    OccupancyVolume( const AABB& aabb, const double h );
    OccupancyVolume( const OccupancyVolume& o );
   ~OccupancyVolume();
  
    OccupancyVolume& operator=(const OccupancyVolume& o );

    void populate(const DynamicalStateData& pq);

  
    size_t nbx() const;
    size_t nby() const;
    size_t nbz() const;
    size_t cell_nb(size_t i, size_t j, size_t k) const;
    const std::vector<size_t>& cell_contents(size_t i, size_t j, size_t k) const;
    Vector cell_center_of_mass(size_t i, size_t j, size_t k) const;
    double cell_total_mass(size_t i, size_t j, size_t k) const;

    void neighbor_cells( size_t i, size_t j, size_t k, std::vector<size_t>& neighbors ) const;
    void distant_cells( size_t i, size_t j, size_t k, std::vector<size_t>& neighbors ) const;

    size_t index( const Vector& P ) const;
    size_t index( size_t i, size_t j, size_t k ) const;
    void anti_index( const size_t ind, size_t& i, size_t& j, size_t& k ) const;

    const double& get_cellsize() const { return cellsize; }
    void set_cellsize(const double& v );

    void set_bounds( const Vector& llc, const Vector& urc );
    const Vector& get_bounds_llc() const;
    const Vector& get_bounds_urc() const;

    void recompute_bounds(const DynamicalStateData& pq);
  private:

    AABB bounds;
    double cellsize;
    size_t nx, ny, nz;
    std::vector< std::vector<size_t> > contents;
    std::vector< Vector > center_of_mass;
    std::vector< double > total_mass;
    pba::Vector LX;

    void compute_size();



};



}

#endif
