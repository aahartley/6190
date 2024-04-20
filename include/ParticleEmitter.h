//*******************************************************************
//
//   ParticleEmitter.h
//
//  Emit particles
//
//  Copyright (c) 2017 Jerry Tessendorf
//
//
//*******************************************************************

#ifndef __PBA_PARTICLEEMITTER_H__
#define __PBA_PARTICLEEMITTER_H__

#include "Vector.h"
#include "Color.h"
#include "SPHState.h"

namespace pba{

class ParticleEmitter 
{
  public:

    ParticleEmitter( const Vector& loc, const Vector& velocity, const double R, const double rms_vel );
   ~ParticleEmitter(){}

    virtual void emit( Vector& p, Vector& v, Color& c );
    void emitCube(SPHState& state, int numParticlesPerAxis, int distance, const Vector& center);
    int emission_rate() const { return rate; }
    void set_emission_rate(const int r){ rate = r; }

  private:

    Vector location;
    Vector mean_velocity;
    double radius;
    double velocity_rms;
    int rate;
};





}

#endif
