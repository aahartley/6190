//-------------------------------------------------------
//
//  DFSPHSolver.C
//
//  Solvers for DFSPH Dynamics
//
//  Copyright (c) 2019 Jerry Tessendorf
//
//
//--------------------------------------------------------


#include "DFSPHSolver.h"
#include "ForceLibrary.h"

using namespace pba;
using namespace std;

DFSPHSolver::DFSPHSolver(SPHState& pq, Force& f, double vclamp, double aclamp, DFSPHForce& sphf) :
  PQ(pq),
  force(f),
  sphforce(sphf),
  velocity_clamp(vclamp),
  acceleration_clamp(aclamp),
  dt(1.0/48.0)
{}

void DFSPHSolver::init() 
{
    PQ->populate();
    PQ->compute_density();
    PQ->compute_factor();
}

void DFSPHSolver::solve(const double deltaTime)
{
    force->compute(PQ, dt);
    get_timestep();
    advance_velocity();
    correct_density_error();
    advance_position();
    PQ->populate();
    PQ->compute_density();
    PQ->compute_factor();
    correct_divergence_error();

}

void DFSPHSolver::correct_density_error()
{
//    float densAvg = PQ->average_density();
//    float densRest = 1000;
//    int iter = 0;
//    float threshold = 100.f;
//    while((((densAvg - densRest) > threshold) || iter < 2)&& iter < 100)
//    {
//       //predict dens
//       PQ->compute_predicted_density(dt);
//       #pragma omp parallel for
//       for(size_t p = 0; p < PQ->nb(); p++)
//       {
//          const Vector P = PQ->pos(p);
//          const Vector V = PQ->vel(p);
//          float density = PQ->get_float_attr("density", p);
//          size_t pindex = PQ->index(P);
//          size_t i,j,k;
//          PQ->anti_index(pindex, i,j,k);
//          float ki = ((PQ->get_float_attr("predicted_density", p) - densRest) / (dt * dt)) * PQ->get_float_attr("factor", p);

//          Vector sum{};
// //cell
//          const std::vector<size_t>& cells = PQ->cell_contents(i,j,k); 
//          for(size_t a=0;a<cells.size();a++)
//          {
//             size_t pid = cells[a]; 
//             float kj = ((PQ->get_float_attr("predicted_density", pid) - densRest) / (dt * dt)) * PQ->get_float_attr("factor", pid);
//             float density_j = PQ->get_float_attr("density", pid);
//             sum += PQ->mass(pid) * ( (ki / density) + (kj / density_j) ) * PQ->grad_weight(pid, P);
//          }
// //neighbors
//          std::vector<size_t> neighbors;
//          PQ->neighbor_cells(i,j,k,neighbors);
//          for(size_t a=0;a<neighbors.size();a++)
//          {
//             size_t pid = cells[a]; 
//             float kj = ((PQ->get_float_attr("predicted_density", pid) - densRest) / (dt * dt)) * PQ->get_float_attr("factor", pid);
//             float density_j = PQ->get_float_attr("density", pid);
//             sum += PQ->mass(pid) * ( (ki / density) + (kj / density_j) ) * PQ->grad_weight(pid, P);
//          }

//          PQ->set_vel(p, V -  (dt * sum));

//       }
//       //what to do about new avg??
//       densAvg = PQ->average_predicted_density();
//       iter++;
      
//    }
}

void DFSPHSolver::correct_divergence_error()
{
//    PQ->compute_density_derivative(); //use from prev?
//    float densDAvg = PQ->average_density_derivative();
//    int iter = 0;
//    float threshold = 10.f;
//    while((densDAvg > threshold || iter < 1) && iter < 100)   {
//       //predict dens
//       PQ->compute_density_derivative();
//       #pragma omp parallel for
//       for(size_t p = 0; p < PQ->nb(); p++)
//       {
//          const Vector P = PQ->pos(p);
//          const Vector V = PQ->vel(p);
//          float density = PQ->get_float_attr("density", p);
//          size_t pindex = PQ->index(P);
//          size_t i,j,k;
//          PQ->anti_index(pindex, i,j,k);
//          float ki = (1.0/dt) * PQ->get_float_attr("density_derivative",p) * PQ->get_float_attr("factor", p);

//          Vector sum{};
// //cell
//          const std::vector<size_t>& cells = PQ->cell_contents(i,j,k); 
//          for(size_t a=0;a<cells.size();a++)
//          {
//             size_t pid = cells[a]; 
//             float kj = (1.0/dt) * PQ->get_float_attr("density_derivative",pid) * PQ->get_float_attr("factor", pid);
//             float density_j = PQ->get_float_attr("density", pid);
//             sum += PQ->mass(pid) * ( (ki / density) + (kj / density_j) ) * PQ->grad_weight(pid, P);
//          }
// //neighbors
//          std::vector<size_t> neighbors;
//          PQ->neighbor_cells(i,j,k,neighbors);
//          for(size_t a=0;a<neighbors.size();a++)
//          {
//             size_t pid = cells[a]; 
//             float kj = (1.0/dt) * PQ->get_float_attr("density_derivative",pid) * PQ->get_float_attr("factor", pid);
//             float density_j = PQ->get_float_attr("density", pid);
//             sum += PQ->mass(pid) * ( (ki / density) + (kj / density_j) ) * PQ->grad_weight(pid, P);
//          }

//          PQ->set_vel(p, V -  (dt * sum));

//       }
//       //what to do about new avg??
//       //PQ->compute_density_derivative();
//       densDAvg = PQ->average_density_derivative();
//       iter++;
      
//    }
}

void DFSPHSolver::get_timestep()
{
    double pdiameter = PQ->get_radius() / 2.0;
    float max_vel = PQ->max_velocity();
    float lambda = 0.4f;
    double max_dt;
    if(max_vel != 0)
       max_dt = lambda * (pdiameter/max_vel);
    else max_dt = 1.0/48.0;   
    if(dt > max_dt) dt = max_dt;
    if(dt > 1.0/48.0) dt = 1.0/48.0;
    //dt=1.0/48.0;

}

void DFSPHSolver::advance_velocity()
{
#pragma omp parallel for
    for( size_t i=0;i<PQ->nb();i++ )
    {    
      Vector A = PQ->accel(i);
      float Amag = A.magnitude();
      if(Amag > acceleration_clamp)
      {
         A *= acceleration_clamp/Amag;
      }
      Vector V = PQ->vel(i) + A*dt;
      float Vmag = V.magnitude();
      if(Vmag > velocity_clamp)
      {
         V *= velocity_clamp/Vmag;
      }
      PQ->set_vel( i, V );
   }
}

void DFSPHSolver::advance_position()
{
#pragma omp parallel for
   for( size_t i=0;i<PQ->nb();i++ )
   {
      PQ->set_pos( i, PQ->pos(i) + PQ->vel(i)*dt );
   }
}

pba::GISolver pba::CreateDFSPHSolver( SPHState& pq, Force& f, float vclamp, float aclamp, DFSPHForce& sphforce )
{
   return GISolver( new DFSPHSolver( pq, f, vclamp, aclamp, sphforce ) );
}



DFSPHSolverWithCollisions::DFSPHSolverWithCollisions(SPHState& pq, Force& f, double vclamp, double aclamp, DFSPHForce& sphf, ElasticCollisionHandler& coll) :
  PQ(pq),
  force(f),
  sphforce(sphf),
  CS(coll),
  velocity_clamp(vclamp),
  acceleration_clamp(aclamp),
  dt(1.0/48.0)
{  
}

void DFSPHSolverWithCollisions::init()
{
    PQ->populate();
    PQ->compute_density();
    PQ->compute_factor();
}

void DFSPHSolverWithCollisions::solve(const double userdt)
{
    user_dt = userdt;
    //dt = user_dt;
    correct_divergence_error();
    #pragma omp parallel for
    for(size_t p = 0; p < PQ->nb(); p++)
    {
      PQ->set_accel(p, Vector(0,0,0));
    }

    force->compute(PQ, dt);
    get_timestep();
    advance_velocity();
    correct_density_error();
    advance_position();
    DynamicalState ss =  std::dynamic_pointer_cast<DynamicalStateData, SPHStateData>(PQ);
    CS.handle_collisions( dt, ss );
    PQ->populate();
    PQ->compute_density();
    PQ->compute_factor();


}

void DFSPHSolverWithCollisions::correct_density_error()
{

   #pragma omp parallel for
   for(size_t p = 0; p < PQ->nb(); p++)
   {
      PQ->compute_predicted_density(p, dt); //density solver, above equation 12
      const float s_i = 1.0f - PQ->get_float_attr("predicted_density", p); //find deviation of density
      const float residuum = std::min(s_i, 0.0f); //ensures predicted density its not negative?
      PQ->set_attr("density_pressure", p, -residuum * (PQ->get_float_attr("factor",p)* (1.0f/(dt*dt))) ); //pressure to find acc
   }

   int iter = 0;
   float average_density_error = 0.0f;
   bool check = false;

   while((!check || iter < 2) && iter < PQ->get_maxIter())   
   {
      check = true;
      average_density_error = 0.0f;
      density_solve_iteration(average_density_error);

      const float eta = 0.01 * PQ->get_density0() * PQ->get_mMaxError(); //0.01 == m_maxError,  this equals 0.1??
      check = check && (average_density_error <= eta);
      iter++;
   }  

   std::cout << "DE dens error: " <<average_density_error << '\n';

   #pragma omp parallel for
   for(size_t p = 0; p < PQ->nb(); p++)
   {
      compute_pressure_acc(p, "density"); // ∑_j m_j( κ_vi / ρ_i + κ_vj / ρj)∇Wij, in the 2 solvers
      Vector V = PQ->vel(p);
      PQ->set_vel(p, V - dt * PQ->get_vector_attr("pressure_acc", p) );
   }


//    PQ->compute_predicted_density(dt);
//    float densAvg = PQ->average_predicted_density();
//    //std::cout<< "avg: " << densAvg << ' ' << dt << '\n';
//    float densRest = 1000;
//    int iter = 0;
//    float threshold = 0.1;
//    while(((std::abs(densAvg - densRest) > threshold) || iter < 2)&& iter < 100)
//    {
//       //predict dens
//       //PQ->compute_predicted_density(dt);
//       #pragma omp parallel for
//       for(size_t p = 0; p < PQ->nb(); p++)
//       {
//          const Vector P = PQ->pos(p);
//          const Vector V = PQ->vel(p);
//          float density = PQ->get_float_attr("density", p);
//          if(density == 0) std::cout << "zero danger\n";
//          size_t pindex = PQ->index(P);
//          size_t i,j,k;
//          PQ->anti_index(pindex, i,j,k);
//          float ki = ((PQ->get_float_attr("predicted_density", p) - densRest) / (dt * dt)) * PQ->get_float_attr("factor", p);
//          if(ki < 0.25) ki =0.25;
//          //#pragma omp critical
//          //std::cout << "iter:"<<iter<< " particle:"<<p<<" ki: \t" << ki <<'\n';

//          Vector sum{0,0,0};
// //cell
//          const std::vector<size_t>& cells = PQ->cell_contents(i,j,k); 
//          for(size_t a=0;a<cells.size();a++)
//          {
//             size_t pid = cells[a]; 
//             float kj = ((PQ->get_float_attr("predicted_density", pid) - densRest) / (dt * dt)) * PQ->get_float_attr("factor", pid);
//             if(kj < 0.25) kj = 0.25;
//             //#pragma omp critical
//             //std::cout << "iter:"<<iter<< " particle:"<<p<<" kj: \t" << kj << '\n';

//             float density_j = PQ->get_float_attr("density", pid);
//             if(density_j == 0) std::cout << "zero danger\n";
//             sum += PQ->mass(pid) * ( (ki / density) + (kj / density_j) ) * PQ->grad_weight(pid, P);
         
//          }
// //neighbors
//          std::vector<size_t> neighbors;
//          PQ->neighbor_cells(i,j,k,neighbors);
//          for(size_t a=0;a<neighbors.size();a++)
//          {
//             size_t pid = neighbors[a]; 
//             float kj = ((PQ->get_float_attr("predicted_density", pid) - densRest) / (dt * dt)) * PQ->get_float_attr("factor", pid);
//             if(kj < 0.25) kj=0.25;
//             float density_j = PQ->get_float_attr("density", pid);
//             if(density_j == 0) std::cout << "zero danger\n";
//             sum += PQ->mass(pid) * ( (ki / density) + (kj / density_j) ) * PQ->grad_weight(pid, P);
//          }
//          //#pragma omp critical
//          //std::cout << "iter:"<<iter<< " particle:"<<p<<" sum: \t" << sum.X() << ' ' << sum.Y() << ' ' << sum.Z() << '\n';
//          //#pragma omp critical
//          PQ->set_vel(p, V -  (dt * sum));
//          //#pragma omp critical
//          //std::cout << "iter:"<<iter<<" particle:"<< p <<" vel: \t" << PQ->vel(p).X() << ' ' << PQ->vel(p).Y() << ' ' << PQ->vel(p).Z() << '\n';

//       }
//       PQ->compute_predicted_density(dt);
//       densAvg = PQ->average_predicted_density();
//       //std::cout << "predicted density in loop: " << densAvg << '\n';
//       iter++;
      
//    }
}


void DFSPHSolverWithCollisions::density_solve_iteration(float& average_density_error)
{
   float density_error = 0.0f;


   #pragma omp parallel for
   for(size_t p = 0; p < PQ->nb(); p++)
   {
      compute_pressure_acc(p, "density"); // ∑_j m_j( κ_vi / ρ_i + κ_vj / ρ_j)∇Wij
   }   
   #pragma omp  parallel for reduction(+:density_error)
   for(size_t p = 0; p < PQ->nb(); p++)
   {
      //forces to correct density error
      float aij_pj = compute_aij_pj(p);//∑_j m_j( f_Pi / ρ_i - f_Pj / ρ_j)∇Wij
      aij_pj *= dt * dt;

      const float predicted_density = PQ->get_float_attr("predicted_density", p);
      const float s_i = 1.0f - predicted_density;

      float density_pressure = PQ->get_float_attr("density_pressure", p);
      const float residuum = std::min(s_i - aij_pj, 0.0f); //pdens - forces to correct

      //density pressure - half of the error instead of pressure - rest density
      density_pressure = std::max(density_pressure - 0.5f * (s_i - aij_pj) * (PQ->get_float_attr("factor",p) * (1.0f/(dt*dt))), 0.0f );
      PQ->set_attr("density_pressure", p, density_pressure);

      density_error -= PQ->get_density0() * residuum;

   }

   average_density_error = density_error / PQ->nb();
   // std::cout << "dens error: " <<average_density_error << '\n';
}


void DFSPHSolverWithCollisions::correct_divergence_error()
{

   #pragma omp parallel for
   for(size_t p = 0; p < PQ->nb(); p++)
   {
      PQ->compute_density_derivative(p); // equation 9: Dρi/Dt = ∑j mj (vi − vj )∇Wij 
      float density_derivative = PQ->get_float_attr("density_derivative", p);
      density_derivative = std::max(density_derivative, 0.0f); //divergence free so we want density to never be negative

      int num_neighbors = 0;
      const Vector P = PQ->pos(p);
      size_t pindex = PQ->index(P);
      size_t i,j,k;
      PQ->anti_index(pindex, i,j,k);
      const std::vector<size_t>& cells = PQ->cell_contents(i,j,k);
      std::vector<size_t> neighbors;
      PQ->neighbor_cells(i,j,k,neighbors);
      num_neighbors += (int)cells.size()-1; //i dont think they are counting the particle itself
      num_neighbors += (int)neighbors.size();
      //do not solve for particles without enough influence
      if(num_neighbors < 20)
      {
         //PQ->set_attr("density_derivative", p, 0.0f);
         density_derivative = 0.0f;
      }
      //float k_v = PQ->get_float_attr("density_derivative", p) * (PQ->get_float_attr("factor", p) * (1.0f/dt));
      float k_v = density_derivative * (PQ->get_float_attr("factor", p) * (1.0f/dt)); //pressure value used to compute acc
      PQ->set_attr("divergence_pressure", p, k_v);

   }

   int iter = 0;

   float average_density_error = 0.0f;
   bool check = false;

   while((!check || iter < 1) && iter < PQ->get_maxIter())   
   {
      check = true;
      
      average_density_error = 0.0f;
      divergence_solve_iteration(average_density_error);
      //this jsut equals 1/dt?? 
      float eta = (1.0f / dt) * PQ->get_maxError() * 0.01f * PQ->get_density0();  // maxError is given in percent, 0.1 == max error
      check = check && (average_density_error <= eta);
      
      iter++;
   }
   std::cout << "DivE dens error: " <<average_density_error << '\n';

   #pragma omp parallel for 
   for(size_t p = 0; p < PQ->nb(); p++)
   {
      compute_pressure_acc(p, "divergence"); // ∑_j m_j( κ_vi / ρ_i + κ_vj / ρj)∇Wij, in the 2 solvers
      Vector V = PQ->vel(p);
      PQ->set_vel(p, V - dt * PQ->get_vector_attr("pressure_acc", p) );
   }

//    float densDAvg = PQ->average_density_derivative();
//    int iter = 0;
//    float threshold = 1.f * (1.0/dt);
//    threshold = 10;
//    while(std::abs(densDAvg > threshold || iter < 1) && iter < 100)   {
//       //predict dens
//       //PQ->compute_density_derivative();
//       #pragma omp parallel for
//       for(size_t p = 0; p < PQ->nb(); p++)
//       {
//          const Vector P = PQ->pos(p);
//          const Vector V = PQ->vel(p);
//          float density = PQ->get_float_attr("density", p);
//          if(density == 0) std::cout << "zero danger\n";
//          size_t pindex = PQ->index(P);
//          size_t i,j,k;
//          PQ->anti_index(pindex, i,j,k);
//          float ki = (1.0/dt) * PQ->get_float_attr("density_derivative",p) * PQ->get_float_attr("factor", p);
//          if(ki < 0.5) ki =0.5; ki *=0.5;
//          Vector sum{0,0,0};
// //cell
//          const std::vector<size_t>& cells = PQ->cell_contents(i,j,k);
//          for(size_t a=0;a<cells.size();a++)
//          {
//             size_t pid = cells[a]; 
//             float kj = (1.0/dt) * PQ->get_float_attr("density_derivative",pid) * PQ->get_float_attr("factor", pid);
//             if(kj < 0.5) kj=0.5; kj*=0.5;
//             float density_j = PQ->get_float_attr("density", pid);
//             if(density_j == 0) std::cout << "zero danger\n";

//             sum += PQ->mass(pid) * ( (ki / density) + (kj / density_j) ) * PQ->grad_weight(pid, P);
//          }
// //neighbors
//          std::vector<size_t> neighbors;
//          PQ->neighbor_cells(i,j,k,neighbors);
//          for(size_t a=0;a<neighbors.size();a++)
//          {
//             size_t pid = neighbors[a]; 
//             float kj = (1.0/dt) * PQ->get_float_attr("density_derivative",pid) * PQ->get_float_attr("factor", pid);
//             if(kj < 0.5) kj=0.5; kj*=0.5;
//             float density_j = PQ->get_float_attr("density", pid);
//             if(density_j == 0) std::cout << "zero danger\n";
//             sum += PQ->mass(pid) * ( (ki / density) + (kj / density_j) ) * PQ->grad_weight(pid, P);
//          }
//          //#pragma omp critical
//          //#pragma omp critical
//          //std::cout << "diver iter:"<<iter<<" particle:"<< p <<" vel: \t" << PQ->vel(p).X() << ' ' << PQ->vel(p).Y() << ' ' << PQ->vel(p).Z() << '\n';
//          PQ->set_vel(p, V -  (dt * sum));

//       }
//       PQ->compute_density_derivative();
//       densDAvg = PQ->average_density_derivative();
//       //std::cout << "density derivative in loop: " << densDAvg << '\n';
//       iter++;
      
//    }
}

void DFSPHSolverWithCollisions::divergence_solve_iteration(float& average_density_error)
{
   float density_error = 0.0f;

   #pragma omp parallel for
   for(size_t p = 0; p < PQ->nb(); p++)
   {
      compute_pressure_acc(p, "divergence"); // ∑_j m_j( κ_vi / ρ_i + κ_vj / ρ_j)∇Wij
   }
   #pragma omp parallel for reduction(+:density_error)
   for(size_t p = 0; p < PQ->nb(); p++)
   {
      //forces to correct density error for pedicted dens??
      float aij_pj = compute_aij_pj(p); //∑_j m_j( f_Pi / ρ_i - f_Pj / ρ_j)∇Wij
      aij_pj *= dt; // equation says dt^2???

      const float density_derivative = PQ->get_float_attr("density_derivative", p);
      const float s_i = -density_derivative;

      float divergence_pressure = PQ->get_float_attr("divergence_pressure", p);
      float residuum = std::min(s_i - aij_pj, 0.0f); //-dd - forces to correct dens, should be -?

      int num_neighbors = 0;
      const Vector P = PQ->pos(p);
      size_t pindex = PQ->index(P);
      size_t i,j,k;
      PQ->anti_index(pindex, i,j,k);
      const std::vector<size_t>& cells = PQ->cell_contents(i,j,k);
      std::vector<size_t> neighbors;
      PQ->neighbor_cells(i,j,k,neighbors);
      num_neighbors += (int)cells.size()-1; //i dont think they are counting the particle itself
      num_neighbors += (int)neighbors.size();

      if(num_neighbors < 20)
         residuum = 0.0f;
      divergence_pressure = std::max(divergence_pressure - 0.5f*(s_i - aij_pj) * (PQ->get_float_attr("factor",p)*(1.0f/dt)), 0.0f); 
      PQ->set_attr("divergence_pressure", p, divergence_pressure);

      density_error -= PQ->get_density0() * residuum; // - * 1000, adds % of error?

   }

   average_density_error = density_error / PQ->nb();


}

void DFSPHSolverWithCollisions::compute_pressure_acc(size_t p, std::string type)
{
   Vector pressure_acci_i(0.f,0.f,0.f);
   float ki;
   if(type == "divergence")
      ki = PQ->get_float_attr("divergence_pressure", p);
   else if(type == "density")
      ki = PQ->get_float_attr("density_pressure", p);
   else std::cout <<"ERROR_compute_pressure_acc\n";

   const Vector P = PQ->pos(p);
   size_t pindex = PQ->index(P);
   size_t i,j,k;
   PQ->anti_index(pindex, i,j,k);
   const std::vector<size_t>& cells = PQ->cell_contents(i,j,k);
   std::vector<size_t> neighbors;
   PQ->neighbor_cells(i,j,k,neighbors);

//cell
   for(size_t a=0;a<cells.size();a++)
   {
      size_t pid = cells[a]; 
      if(pid != p)
      {
         float kj;
         if(type == "divergence")
            kj = PQ->get_float_attr("divergence_pressure", pid);
         else if(type == "density")
            kj = PQ->get_float_attr("density_pressure", pid);
         float psum = ki + PQ->get_density0()/PQ->get_density0()* kj;
         if(fabs(psum) > PQ->get_meps())//m_eps = (1.0e-5); b/c we set pressure to 0 for certain scenarios
         {
            Vector grad_pj = PQ->get_float_attr("volume", pid) * PQ->grad_weight(pid, P);
            pressure_acci_i += psum * grad_pj;
         }
      }
   }
//neighbors
   for(size_t a=0;a<neighbors.size();a++)
   {
      size_t pid = neighbors[a]; 
   
      float kj;
      if(type == "divergence")
         kj = PQ->get_float_attr("divergence_pressure", pid);
      else if(type == "density")
         kj = PQ->get_float_attr("density_pressure", pid);
      float psum = ki + PQ->get_density0()/PQ->get_density0()* kj;
      if(fabs(psum) > PQ->get_meps()) //m_eps = (1.0e-5);
      {
         Vector grad_pj = PQ->get_float_attr("volume", pid) * PQ->grad_weight(pid, P);
         pressure_acci_i += psum * grad_pj;
      }
      
   }
   PQ->set_attr("pressure_acc", p, pressure_acci_i);

}
//forces to correct density error p* - p0, rhs equation 12
float DFSPHSolverWithCollisions::compute_aij_pj(size_t p)
{
   float aij_pj = 0.f;

   const Vector P = PQ->pos(p);
   Vector pressure_acc_i = PQ->get_vector_attr("pressure_acc", p);

   size_t pindex = PQ->index(P);
   size_t i,j,k;
   PQ->anti_index(pindex, i,j,k);
   const std::vector<size_t>& cells = PQ->cell_contents(i,j,k);
   std::vector<size_t> neighbors;
   PQ->neighbor_cells(i,j,k,neighbors);

//cell
   for(size_t a=0;a<cells.size();a++)
   {
      size_t pid = cells[a]; 
      if(pid != p)
      {
         Vector pressure_acc_j = PQ->get_vector_attr("pressure_acc", pid);
         aij_pj += (pressure_acc_i - pressure_acc_j) * PQ->grad_weight(pid, P); // dot with gradient
      }
   }
//neighbors
   for(size_t a=0;a<neighbors.size();a++)
   {
      size_t pid = neighbors[a]; 
      Vector pressure_acc_j = PQ->get_vector_attr("pressure_acc", pid);
      aij_pj += (pressure_acc_i - pressure_acc_j) * PQ->grad_weight(pid, P);
      
   }
   aij_pj *= PQ->get_float_attr("volume", p);
   return aij_pj;
}

void DFSPHSolverWithCollisions::get_timestep()
{
    double pdiameter = PQ->get_radius() / 2.0;
    float max_vel = PQ->max_velocity();
    float lambda = 0.4f;
    double max_dt;
    if(max_vel != 0)
       max_dt = lambda * (pdiameter/max_vel);
    else max_dt = user_dt;   
    //if(dt > max_dt) 

    dt = max_dt;
    if(dt > user_dt || dt==0) dt = user_dt;
   //  if(user_dt != dt)
   //    std::cout << "dt: " << dt << '\n';

   //  dt_accum += dt;
   //  if(dt_accum >= user_dt)
   //  {
   //    dt_accum = 0;
   //    std::cout << "dt: " << dt << '\n';
   //  }
   //  else
   //  {
   //    std::cout << "user_dt too big, dt used: " << dt <<  ", time left: " << user_dt - dt_accum << '\n';
   //  }

}

void DFSPHSolverWithCollisions::advance_velocity()
{
#pragma omp parallel for
    for( size_t i=0;i<PQ->nb();i++ )
    {    
      Vector A = PQ->accel(i);
      float Amag = A.magnitude();
      if(Amag > acceleration_clamp)
      {
         A *= acceleration_clamp/Amag;
      }
      Vector V = PQ->vel(i) + A*dt;
      float Vmag = V.magnitude();
      if(Vmag > velocity_clamp)
      {
         V *= velocity_clamp/Vmag;
      }
      PQ->set_vel( i, V );
      //#pragma omp critical
      //std::cout << "particle:"<< i <<" vel: \t" << PQ->vel(i).X() << ' ' << PQ->vel(i).Y() << ' ' << PQ->vel(i).Z() << '\n';

   }
}

void DFSPHSolverWithCollisions::advance_position()
{
#pragma omp parallel for
   for( size_t i=0;i<PQ->nb();i++ )
   {
      PQ->set_pos( i, PQ->pos(i) + PQ->vel(i)*dt );
   }
}

pba::GISolver pba::CreateDFSPHSolver( SPHState& pq, Force& f, float vclamp, float aclamp, DFSPHForce& sphforce, ElasticCollisionHandler& coll )
{
   return GISolver( new DFSPHSolverWithCollisions( pq, f, vclamp, aclamp, sphforce, coll ) );
}
