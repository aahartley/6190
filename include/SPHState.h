//-------------------------------------------------------
//
//  SPHState.h
//
//  Container for data associated with the dynamics
//  degrees of freedom of sph system.
//
//  Copyright (c) 2019 Jerry Tessendorf
//
//
//--------------------------------------------------------

#ifndef ____PBA_SPHSTATE_H____
#define ____PBA_SPHSTATE_H____

#include "DynamicalState.h"
#include "OccupancyVolume.h"
#include <iostream>
#include <fstream>
#include <algorithm>

namespace pba
{


class SPHStateData : public DynamicalStateData, public OccupancyVolume
{
  public:

    SPHStateData( const AABB& bounds, const double h, const std::string& nam = "SPHDataNoName" );
    SPHStateData( const SPHStateData& d );
   ~SPHStateData();

    SPHStateData& operator= ( const SPHStateData& d );

    const float get_radius() const { return radius; }
    void set_radius( const float& v );

    const float get_particle_radius() const { return particle_radius; }
    const float get_density0() const { return density0; }
    void set_density0(const float d0);
    const int get_maxIter() const { return maxIter; }
    void set_maxIter(int maxi);
    const float get_meps() const { return m_eps;}
    const float get_mMaxError() const { return m_maxError;}
    void set_mMaxError(int mmx);
    const float get_maxError() const { return maxError;}
    void set_maxError(int mx);


    const float weight( size_t p, const Vector& P ) const;
    const Vector grad_weight( size_t p, const Vector& P ) const;

    void compute_density();
    void compute_predicted_density(const size_t p, const double dt);
    void compute_divergence();
    void compute_density_derivative(size_t p);
    void compute_factor();
    void adjust_density_for_divergence(const double dt);
    void populate();

    float average_density(); //const;
    float average_density_derivative(); //const;
    float average_predicted_density(); //const;
    float average_divergence() const;
    float max_velocity() const;

    std::ofstream file;
    int densiter =0;
    int densavgiter =0;
    int pdensiter =0;
    int pdensavgiter =0;
    int ddensiter =0;
    int ddensavgiter =0;
  private:

    float radius;
    float particle_radius;
    float density0;
    float m_eps;
    float maxError; //divergence
    float m_maxError; //density
    int maxIter;



};



typedef std::shared_ptr<SPHStateData> SPHState;

SPHState CreateSPH( const AABB& bounds, const double h, const std::string& nam = "SPHDataNoName" );

SPHState copy( const SPHState d );




}
#endif
