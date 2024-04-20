//-------------------------------------------------------
//
//  DFSPHSPHSolver.C
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
//sim loop does this already --- no need call
//requirements for sim to start
void DFSPHSolverWithCollisions::init()
{
    PQ->populate();
    PQ->compute_density();
    PQ->compute_factor();
}

//naive collisions to keep particles in box
void DFSPHSolverWithCollisions::fakecs()
{
   float length = 10 * (0.025*2);
   float botx = -(length/2)*2;
   float WIDTHx = std::abs(botx) ;

   float boty = -length;
   float WIDTHy = std::abs(boty) ;

   float botz = -(length/2);
   float WIDTHz = std::abs(botz) ;

   float radius = PQ->get_particle_radius();
   float coef = 0.9;
   #pragma omp parallel for
    for(int i = 0; i < PQ->nb(); i++)
    {
        if(PQ->pos(i).X() >= WIDTHx-radius)
        { 
            PQ->set_pos(i, Vector((WIDTHx-radius), PQ->pos(i).Y(), PQ->pos(i).Z()));
            PQ->set_vel(i, Vector((PQ->vel(i).X() * -coef),PQ->vel(i).Y(),PQ->vel(i).Z()));
            //PQ->set_vel(i, Vector(0,PQ->vel(i).Y(),PQ->vel(i).Z()));

        }
        if(PQ->pos(i).X() <= botx+radius)
        {
            PQ->set_pos(i, Vector((botx+radius), PQ->pos(i).Y(), PQ->pos(i).Z()));
            PQ->set_vel(i, Vector((PQ->vel(i).X() * -coef),PQ->vel(i).Y(),PQ->vel(i).Z()));
            //PQ->set_vel(i, Vector(0,PQ->vel(i).Y(),PQ->vel(i).Z()));

        }
        if(PQ->pos(i).Y() >= WIDTHy-radius)
        { 
            PQ->set_pos(i, Vector(PQ->pos(i).X(), (WIDTHy-radius), PQ->pos(i).Z()));
            PQ->set_vel(i, Vector(PQ->vel(i).X(),(PQ->vel(i).Y() * -coef),PQ->vel(i).Z()));
            //PQ->set_vel(i, Vector(PQ->vel(i).X(),0,PQ->vel(i).Z()));

        }
        if(PQ->pos(i).Y() <= boty+radius)
        {
            //std::cout << "botplane hit\n";
            PQ->set_pos(i, Vector(PQ->pos(i).X(), (boty+radius), PQ->pos(i).Z()));
            PQ->set_vel(i, Vector(PQ->vel(i).X(),(PQ->vel(i).Y() * -coef),PQ->vel(i).Z()));
            //PQ->set_vel(i, Vector(PQ->vel(i).X(),0,PQ->vel(i).Z()));

        }
         if(PQ->pos(i).Z() >= WIDTHz-radius)
        { 
            PQ->set_pos(i, Vector(PQ->pos(i).X(), PQ->pos(i).Y(),(WIDTHz-radius)));
            PQ->set_vel(i, Vector(PQ->vel(i).X(),PQ->vel(i).Y(),(PQ->vel(i).Z() * -coef)));
            //PQ->set_vel(i, Vector(PQ->vel(i).X(),PQ->vel(i).Y(),0));

        }
        if(PQ->pos(i).Z() <= botz+radius)
        {
            PQ->set_pos(i, Vector(PQ->pos(i).X(), PQ->pos(i).Y(),(botz+radius)));
            PQ->set_vel(i, Vector(PQ->vel(i).X(),PQ->vel(i).Y(),(PQ->vel(i).Z() * -coef)));
            //PQ->set_vel(i, Vector(PQ->vel(i).X(),PQ->vel(i).Y(),0));

        }


    }
}

void DFSPHSolverWithCollisions::solve(const double userdt)
{
    user_dt = userdt;
    //dt = user_dt;

    //Occupancy Grid
    PQ->populate(); 
    PQ->compute_density();
    //Equ (8)
    PQ->compute_factor();

    //PPE 2 (Alg. 2): Dρi/Dt = 0
    correct_divergence_error();

    //Non-pressure forces (viscosity, gravity)
    //TODO: surface tension?
    force->compute(PQ, dt);

    //CFL Condition
    get_timestep();

    // ___ Euler integration
    advance_velocity();

    //correct_divergence_error();
    // PPE 1 (Alg. 3): ρ∗i - ρ0 = 0
    correct_density_error();

    // ___ Euler integration
    advance_position();

    //Boundary Handling
    fakecs();

    //auto start2 = std::chrono::high_resolution_clock::now();

   //  auto end2 = std::chrono::high_resolution_clock::now();
   //  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2);
   //  std::cout << "Time taken by function: " << duration.count() << " ms" << std::endl;

    DynamicalState ss =  std::dynamic_pointer_cast<DynamicalStateData, SPHStateData>(PQ);
    //CS.handle_collisions( dt, ss );



}

// PPE 1 (Alg. 3): ρ∗i - ρ0 = 0
void DFSPHSolverWithCollisions::correct_density_error()
{

   #pragma omp parallel for
   for(size_t p = 0; p < PQ->nb(); p++)
   {
      //Equation above (12)
      PQ->compute_predicted_density(p, dt);

      ////const float s_i = 1.0f - PQ->get_float_attr("predicted_density", p); //find deviation of density (1 is used b/c pdens had p/p0)
      //ρ∗i - ρ0
      const float s_i =  PQ->get_float_attr("predicted_density", p) - PQ->get_density0(); 

      ////si > 0 = pdens is < dens0 || pdens <=0, if si < 0 = pdens is > dens0
      //r = b - A * p
      //r = s_i - A * density_pressure
      ////const float residuum = std::min(s_i, 0.0f); // if s_i is neg, that means compression. clamp to 0 dont worry about expansion ,i.e pdens is negative
      //clamp source -- if source is neg == expansion, only solve for compression
      const float residuum = std::max(s_i, 0.0f); 

      ////PQ->set_attr("density_pressure", p, -residuum * (PQ->get_float_attr("factor",p)* (1.0f/(dt*dt))) ); //pressure stiffness to find acc(ki), -residuum to counteract negation prev
      //Under Equation (9)
      //pressure stifness parameter ki = 1/∆t^2 (ρ∗i - ρ0) αi 
      PQ->set_attr("density_pressure", p, residuum * (PQ->get_float_attr("factor",p)* (1.0f/(dt*dt))) ); 

   }

   int iter = 0;
   float average_density_error = 0.0f;
   bool check = false;

   while((!check || iter < 2) && iter < PQ->get_maxIter())   
   {
      check = true;
      average_density_error = 0.0f;
      //Jacobi Iteration to solve for ki
      density_solve_iteration(average_density_error);

      //0.01 == m_maxError,  m_maxError is given as a percent. p0*maxerror = 0.01% * density0
      const float eta = PQ->get_density0() * PQ->get_mMaxError() * 0.01; 
      check = check && (average_density_error <= eta);
      iter++;
   }  

   std::cout << "DE dens error: " <<average_density_error << '\n';

   #pragma omp parallel for
   for(size_t p = 0; p < PQ->nb(); p++)
   {
      //End of section 3.2 is Fpi,total = -mi ∑_j m_j( κ_vi / ρ_i + κ_vj / ρj)∇Wij
      //this is acc =  -∑_j m_j( κ_vi / ρ_i + κ_vj / ρj)∇Wij
      compute_pressure_acc(p, "density"); 
      Vector V = PQ->vel(p);
      //pressure_acc has negative sign accounted for, unlike psuedo-code in paper
      PQ->set_vel(p, V + dt * PQ->get_vector_attr("pressure_acc", p));
      //PQ->set_vel(p, V - dt * PQ->get_vector_attr("pressure_acc", p) );

      //#pragma omp critical
      //{
         // if(std::isnan(PQ->vel(p).X()) || std::isinf(PQ->vel(p).X())) std::cout << "dense PQ->vel(p).X() x bad" << "P: "<< p <<'\n';
         // if(std::isnan(PQ->vel(p).Y()) || std::isinf(PQ->vel(p).Y())) std::cout << "dense PQ->vel(p).Y() y bad" << "P: "<< p <<'\n';
         // if(std::isnan(PQ->vel(p).Z()) || std::isinf(PQ->vel(p).Z())) std::cout << "dense PQ->vel(p).Y() z bad" << "P: "<< p <<'\n';
      //}
   }

}

//Jacobi Iteration to solve for ki
void DFSPHSolverWithCollisions::density_solve_iteration(float& average_density_error)
{
   float density_error = 0.0f;


   #pragma omp parallel for
   for(size_t p = 0; p < PQ->nb(); p++)
   {
      //End of section 3.2 is Fpi,total = -mi ∑_j m_j( κ_vi / ρ_i + κ_vj / ρj)∇Wij
      //this is acc =  -∑_j m_j( κ_vi / ρ_i + κ_vj / ρj)∇Wij
      compute_pressure_acc(p, "density"); 
   }   
   #pragma omp  parallel for reduction(+:density_error)
   for(size_t p = 0; p < PQ->nb(); p++)
   {
      //RHS Equation (9) -- forces to correct density error (ρ∗i - ρ0)
      //∑_j m_j( f_Pi / m_i - f_Pj / m_i)∇Wij
      //which is ∑_j m_j( Acc_Pi - Acc_Pj)∇Wij
      float aij_pj = compute_aij_pj(p);
      //add the timestep
      aij_pj *= dt * dt;
      // #pragma omp critical
      // std::cout << "aijpj : " << aij_pj<< '\n';

      const float predicted_density = PQ->get_float_attr("predicted_density", p);
      //const float s_i = 1.0f - predicted_density;  //si > 0 = pdens is < dens0 || pdens <=0, if si < 0 = pdens is > dens0
      const float s_i =  PQ->get_float_attr("predicted_density", p) - PQ->get_density0(); //find deviation of density (1 is used b/c pdens had p/p0)


      float density_pressure = PQ->get_float_attr("density_pressure", p);
      ////if si > 0 no compression - forces to crrerct <= 0 ==clamp to 0
      ////else si<0 compression - foreces to crrect >=0
      ////const float residuum = std::min(s_i - aij_pj, 0.0f); //pdens - forces to correct dens error (r=negative)

      //Equation (9)
      //(ρ∗i - ρ0) = 1/∆t^2  ∑_j m_j( f_Pi / m_i - f_Pj / m_i)∇Wij (in my derivations it should be -(1/∆t^2), typo in paper?)
      //(ρ∗i - ρ0) = -1/∆t^2  ∑_j m_j( f_Pi / m_i - f_Pj / m_i)∇Wij
      //(ρ∗i - ρ0) + 1/∆t^2  ∑_j m_j( f_Pi / m_i - f_Pj / m_i)∇Wij = error
      const float residuum = std::max(s_i + aij_pj, 0.0f); 


      ////density pressure - half of the error instead of pressure - rest density
      ////if ki ==0, 0 - 0.5(+ - (- or 0)) == clamp to 0 no compression
      ////if compression , + - 0.5(- - (+)) == > 0, ki grows
      ////density_pressure = std::max(density_pressure - 0.5f * (s_i - aij_pj) * (PQ->get_float_attr("factor",p) * (1.0f/(dt*dt))), 0.0f );

      //adjust ki by half of the error during jacobi, clamp to 0 for no negative pressure
      //ki never needs to go down? check
      density_pressure = std::max(density_pressure + 0.5f * (s_i + aij_pj) * (PQ->get_float_attr("factor",p) * (1.0f/(dt*dt))), 0.0f );

      PQ->set_attr("density_pressure", p, density_pressure);

      //density_error -= PQ->get_density0() * residuum; //prev method used volume, scale to mass
      density_error +=  residuum;

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
      density_derivative = std::max(density_derivative, 0.0f); //divergence free so we want density to never be negative? positive derivative means chance of compression?
      PQ->set_attr("density_derivative", p, density_derivative);

      int num_neighbors = 0;
      const Vector P = PQ->pos(p);
      size_t pindex = PQ->index(P);
      size_t i,j,k;
      PQ->anti_index(pindex, i,j,k);
      const std::vector<size_t>& cells = PQ->cell_contents(i,j,k);
      std::vector<size_t> neighbors;
      PQ->neighbor_cells(i,j,k,neighbors);
      num_neighbors += (int)cells.size(); //i dont think they are counting the particle itself -1 or not to -1
      num_neighbors += (int)neighbors.size();
      //do not solve for particles without enough influence
      if(num_neighbors < 20)
      {
         //PQ->set_attr("density_derivative", p, 0.0f); //which one??????????????
         density_derivative = 0.0f;
         PQ->set_attr("density_derivative", p, density_derivative);

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
      //this jsut equals 1/dt?? 1/0.005 is 200
      //we want density error dp/dt - forces*dt = 0 to be less than the rate of change of time * a fraction of the rest denstiy
      //0.1% * 1000=1
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
      PQ->set_vel(p, V + dt * PQ->get_vector_attr("pressure_acc", p) );
      #pragma omp critical
      {
         // if(std::isnan(PQ->vel(p).X()) || std::isinf(PQ->vel(p).X())) std::cout << "dive PQ->vel(p).X() x bad\n";
         // if(std::isnan(PQ->vel(p).Y()) || std::isinf(PQ->vel(p).Y())) std::cout << "dive PQ->vel(p).Y() x bad\n";
         // if(std::isnan(PQ->vel(p).Z()) || std::isinf(PQ->vel(p).Z())) std::cout << "dive PQ->vel(p).Z() x bad\n";
      }
   }

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
      aij_pj *= dt; // equation says -dt??
      //#pragma omp critical
      //std::cout << "aijpj: " << aij_pj << '\n';
      const float density_derivative = PQ->get_float_attr("density_derivative", p);
      //const float s_i = -density_derivative;
      const float s_i = density_derivative;

      float divergence_pressure = PQ->get_float_attr("divergence_pressure", p);
      //if -dp/dt then + - (- or 0) == clamp to 0
      //else - - +
      //float residuum = std::min(s_i - aij_pj, 0.0f); //-dd - forces to correct dens, should be -?
      float residuum = std::max(s_i + aij_pj, 0.0f); //-dd - forces to correct dens, should be -?

      int num_neighbors = 0;
      const Vector P = PQ->pos(p);
      size_t pindex = PQ->index(P);
      size_t i,j,k;
      PQ->anti_index(pindex, i,j,k);
      const std::vector<size_t>& cells = PQ->cell_contents(i,j,k);
      std::vector<size_t> neighbors;
      PQ->neighbor_cells(i,j,k,neighbors);
      num_neighbors += (int)cells.size(); //i dont think they are counting the particle itself
      num_neighbors += (int)neighbors.size();
      // #pragma omp critical
      // std::cout << "neighbors " << num_neighbors <<'\n';
      if(num_neighbors < 20)
         residuum = 0.0f;
      //divergence_pressure = std::max(divergence_pressure - 0.5f*(s_i - aij_pj) * (PQ->get_float_attr("factor",p)*(1.0f/dt)), 0.0f); 
      divergence_pressure = std::max(divergence_pressure + 0.5f*(s_i + aij_pj) * (PQ->get_float_attr("factor",p)*(1.0f/dt)), 0.0f); 

      PQ->set_attr("divergence_pressure", p, divergence_pressure);
      //if(std::isnan(residuum)) std::cout << "res is nan\n";
      //density_error -= PQ->get_density0() * residuum; // - * 1000, adds % of error?
      density_error += residuum; // - * 1000, adds % of error?


   }

   average_density_error = density_error / PQ->nb();


}

void DFSPHSolverWithCollisions::compute_pressure_acc(size_t p, std::string type)
{
   Vector pressure_acci_i(0.f,0.f,0.f);
   float di = PQ->get_float_attr("density", p);
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
         float dj = PQ->get_float_attr("density", pid);

     // if(pid != p)
      //{
         float kj;
         if(type == "divergence")
            kj = PQ->get_float_attr("divergence_pressure", pid);
         else if(type == "density")
            kj = PQ->get_float_attr("density_pressure", pid);
         //float psum = (ki/di) + (kj/dj);
            float psum = (ki) + PQ->get_density0()/PQ->get_density0()* (kj);

         if(fabs(psum) > PQ->get_meps())//m_eps = (1.0e-5); b/c we set pressure to 0 for certain scenarios
         {
            //Vector grad_pj = -PQ->get_float_attr("volume", pid) * PQ->grad_weight(pid, P); //-volume to cancel out - volume in force for acc
            Vector grad_pj = PQ->get_float_attr("volume", pid) * PQ->grad_weight(pid, P)  * 1000; //-volume to cancel out - volume in force for acc

            pressure_acci_i += psum * grad_pj;
         }
    //  }
   }
//neighbors
   for(size_t a=0;a<neighbors.size();a++)
   {
      size_t pid = neighbors[a]; 
      float dj = PQ->get_float_attr("density", pid);

      float kj;
      if(type == "divergence")
         kj = PQ->get_float_attr("divergence_pressure", pid);
      else if(type == "density")
         kj = PQ->get_float_attr("density_pressure", pid);
      //float psum = (ki/di) + (kj/dj);
       //float psum = (ki) + (kj);
            float psum = (ki) + PQ->get_density0()/PQ->get_density0()* (kj);

      if(fabs(psum) > PQ->get_meps()) //m_eps = (1.0e-5);
      {
         //Vector grad_pj = -PQ->get_float_attr("volume", pid) * PQ->grad_weight(pid, P);
         Vector grad_pj = PQ->get_float_attr("volume", pid) * PQ->grad_weight(pid, P) *1000  ; //-volume to cancel out - volume in force for acc

         pressure_acci_i += psum * grad_pj;
      }
      
   }
   PQ->set_attr("pressure_acc", p, -pressure_acci_i );
   #pragma omp critical
   {
   if(std::isnan(pressure_acci_i.X()) || std::isinf(pressure_acci_i.X())) std::cout << "pressure_acci_i x bad\n";
   if(std::isnan(pressure_acci_i.Y()) || std::isinf(pressure_acci_i.Y())) std::cout << "pressure_acci_i y bad\n";
   if(std::isnan(pressure_acci_i.Z()) || std::isinf(pressure_acci_i.Z())) std::cout << "pressure_acci_i z bad\n";
   // if(pressure_acci_i.X() < 0)  std::cout << "pressure_acci_i x <=0 " << pressure_acci_i.X() << '\n';
   // if(pressure_acci_i.Y() < 0) std::cout << "pressure_acci_i y <=0 " << pressure_acci_i.Y() << '\n';
   // if(pressure_acci_i.Z() < 0) std::cout << "pressure_acci_i z <=0 " << pressure_acci_i.Z() << '\n';
   }
}
//forces to correct density error p* - p0, rhs equation 12
float DFSPHSolverWithCollisions::compute_aij_pj(size_t p)
{
   float aij_pj = 0.f;
   float di = PQ->get_float_attr("density",p);
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
      //if(pid != p)
      //{
         Vector pressure_acc_j = PQ->get_vector_attr("pressure_acc", pid);
         //aij_pj += ((pressure_acc_i/di) - (pressure_acc_j/di)) * PQ->grad_weight(pid, P); // dot with gradient
         aij_pj += PQ->mass(pid) * ((pressure_acc_i) - (pressure_acc_j)) * PQ->grad_weight(pid, P); // dot with gradient

      //}
   }
//neighbors
   for(size_t a=0;a<neighbors.size();a++)
   {
      size_t pid = neighbors[a]; 
      Vector pressure_acc_j = PQ->get_vector_attr("pressure_acc", pid) ;
      //aij_pj += ((pressure_acc_i/di) - (pressure_acc_j/di)) * PQ->grad_weight(pid, P);
      aij_pj += PQ->mass(pid) * ((pressure_acc_i) - (pressure_acc_j)) * PQ->grad_weight(pid, P); // dot with gradient
      
   }
   if(std::isnan(aij_pj) || std::isinf(aij_pj)) std::cout << "aij_pj bad\n";
   // if(aij_pj > 0 || aij_pj < 0)
   // {
   //    #pragma omp critical
   //    std::cout << "aij_pj: " << aij_pj << '\n';
   // }
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
    if(dt > 0.005) dt = 0.005;
    else if (dt < 0.0001) dt = 0.0001;
    std::cout << "Dt: " << dt << '\n';
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
