//*******************************************************************
//
//   ParticleEmitter.C
//
//  Emit particles
//
//  Copyright (c) 2017 Jerry Tessendorf
//
//
//*******************************************************************

#include "ParticleEmitter.h"
#include <cmath>
#include <stdlib.h>


using namespace pba;

ParticleEmitter::ParticleEmitter( const Vector& loc, const Vector& vel, const double R, const double rms_vel ) :
   location      (loc),
   mean_velocity (vel),
   radius        (R),
   velocity_rms  (rms_vel),
   rate(0)
{
   //prn.setSeed(57575);
}

void ParticleEmitter::emit( Vector& p, Vector& v, Color& c )
{
   float s = 1.0 - 2.0*drand48();
   float ss = std::sqrt( 1.0-s*s );
   float phi = drand48()*2.0*3.14159265;
   p = Vector( ss*std::cos(phi), s, ss*std::sin(phi) );
   v = p*velocity_rms + mean_velocity;
   //v.set(0,0,0);
   p = location + p*radius*std::pow(drand48(),1.0/3.0);
   c = Color( drand48(), drand48(), drand48(), 1.0 );
}


void ParticleEmitter::emitCube(SPHState& state, int numParticlesPerAxis, int distance, const Vector& center)
{
   float diameter = state->get_particle_radius()*distance;

   int numPoints = numParticlesPerAxis * numParticlesPerAxis * numParticlesPerAxis;
   state->add(numPoints);
   int i = 0;

   for (int x = 0; x < numParticlesPerAxis; x++)
   {
      for (int y = 0; y < numParticlesPerAxis; y++)
      {
            for (int z = 0; z < numParticlesPerAxis; z++)
            {
               float tx = x / (numParticlesPerAxis - 1.f);
               float ty = y / (numParticlesPerAxis - 1.f);
               float tz = z / (numParticlesPerAxis - 1.f);

               float px = (tx - 0.5f) * diameter + center.X();
               float py = (ty - 0.5f) * diameter + center.Y();
               float pz = (tz - 0.5f) * diameter + center.Z();
               state->set_pos(i, Vector(px,py,pz));
               state->set_vel(i, Vector(5.0,0,5.0));
               state->set_mass(i, state->get_float_attr("volume", i) * state->get_density0());
               state->set_ci(i, Color(0,0,1,1));
               i++;
            }
      }
   }


}



