//-------------------------------------------------------
//
//  DFSPHForce.h
//
//  A collection of force implementations for DFSPH
//
//  Copyright (c) 2019 Jerry Tessendorf
//
//
//--------------------------------------------------------

#ifndef ____PBA_DFSPHFORCE_H____
#define ____PBA_DFSPHFORCE_H____

#include "Force.h"
#include "SPHState.h"

namespace pba
{


//! All of the sph forces in one object
class DFSPHForce : public ForceBase
{
  public:

    DFSPHForce();
   ~DFSPHForce(){};

    void compute( SPHState& s, const double dt );

    const float& get_dynamic_viscosity() const { return dynamic_viscosity; }

    void set_dynamic_viscosity(const float v) { dynamic_viscosity = v; }


  private:

    float dynamic_viscosity;

};

Force CreateDFSPHForce();






}
#endif
