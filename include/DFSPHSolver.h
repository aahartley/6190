//-------------------------------------------------------
//
//  DFSPHSPHSolver.h
//
//  Solvers for DFSPH Dynamics
//
//  DFSPH functions are based off the 2015 paper: https://animation.rwth-aachen.de/media/papers/2015-SCA-DFSPH.pdf
//  2017 paper, where the equations are written differently (still equivalent): https://animation.rwth-aachen.de/media/papers/2017-TVCG-ViscousDFSPH.pdf
//
//  Copyright (c) 2024 Jerry Tessendorf
//
//
//--------------------------------------------------------

#ifndef ____PBA_DFSPHSOLVER_H____
#define ____PBA_DFSPHSOLVER_H____

#include "SPHState.h"
#include "Force.h"
#include "DFSPHForce.h"
#include "GISolver.h"
#include "CollisionHandler.h"
#include <chrono>


namespace pba
{


class DFSPHSolver : public GISolverBase
{
  public:
    DFSPHSolver(SPHState& pq, Force& f, double vclamp, double aclamp, DFSPHForce& sphf);
    ~DFSPHSolver(){}
    
    void init();
    void solve(const double dt);

    void advance_velocity();
    void advance_position();
    void correct_density_error();
    void correct_divergence_error();
    void get_timestep();

  private:
    SPHState PQ;
    Force force;
    DFSPHForce sphforce;
    float velocity_clamp;
    float acceleration_clamp;
    float dt;


};

GISolver CreateDFSPHSolver( SPHState& pq, Force& f, float vel_clamp, float accel_clamp, DFSPHForce& sphforce );


class DFSPHSolverWithCollisions : public GISolverBase
{
  public:
    DFSPHSolverWithCollisions(SPHState& pq, Force& f, double vclamp, double aclamp, DFSPHForce& sphf, ElasticCollisionHandler& coll);
    ~DFSPHSolverWithCollisions(){}
    
    void init();
    void solve(const double dt);

    void advance_velocity();
    void advance_position();
    void correct_density_error();
    void density_solve_iteration(float& avg_density_error);
    void correct_divergence_error();
    void divergence_solve_iteration(float& avg_density_error);
    void compute_pressure_acc(size_t p, std::string type);
    float compute_aij_pj(size_t p);
    void get_timestep();
    void fakecs();



  private:
    SPHState PQ;
    Force force;
    DFSPHForce sphforce;
    ElasticCollisionHandler& CS;
    float velocity_clamp;
    float acceleration_clamp;
    float user_dt;
    float dt;
    float dt_accum;


};

GISolver CreateDFSPHSolver( SPHState& pq, Force& f, float vel_clamp, float accel_clamp, DFSPHForce& sphforce, ElasticCollisionHandler& cs );



}

#endif