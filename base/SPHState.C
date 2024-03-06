//-------------------------------------------------------
//
//  SPHState.C
//
//  Container for data associated with the dynamics
//  degrees of freedom of a single rigid body.
//
//  Copyright (c) 2017 Jerry Tessendorf
//
//
//--------------------------------------------------------


#include "SPHState.h"
#include <memory>



using namespace pba;



SPHStateData::SPHStateData( const AABB& bounds, const double h, const std::string& nam ) :
  DynamicalStateData(nam+"SPHStateData"),
  OccupancyVolume( bounds, 2.0*h ),
  radius(h),
  density0(1000.f),
  particle_radius(radius/4.0f),
  m_eps(1.0e-5f),
  maxError(0.1f),
  m_maxError(0.01f),
  maxIter(2)
{
   float particle_diamter = particle_radius * 2;
   float mV = 0.8f * particle_diamter * particle_diamter * particle_diamter;

   float default_value = 0.0f;
   Vector default_vector(0.f,0.f,0.f);

   create_attr("density", default_value);
   create_attr("predicted_density", default_value);
   create_attr("density_derivative", default_value);
   create_attr("pressure", default_value );
   create_attr("divergence", default_value );
   create_attr("factor", default_value);
   create_attr("volume", mV);
   create_attr("density_pressure", default_value);
   create_attr("divergence_pressure", default_value);
   create_attr("pressure_acc", default_vector);
   
}
 



SPHStateData::SPHStateData( const SPHStateData& d ) :
  DynamicalStateData(d),
  OccupancyVolume(d),
  radius  (d.radius)
{}
  


SPHStateData::~SPHStateData(){}

SPHStateData& SPHStateData::operator= ( const SPHStateData& d )
{
  DynamicalStateData::operator=(d);
  OccupancyVolume::operator=(d);
  radius = d.radius;
  return *this;
}




SPHState pba::CreateSPH( const AABB& bounds, const double h, const std::string& nam )
{
   return SPHState(new SPHStateData(bounds, h, nam)); 
}

SPHState pba::copy( const SPHState d )
{
   SPHState rbd = SPHState(new SPHStateData(*d));
   return rbd;
}

// const float SPHStateData::weight(size_t p, const Vector& P)const 
// {
//    Vector X = (P- pos(p))/radius;
// 	double q = X.magnitude();

// 	return (1.0/(radius*radius*radius*M_PI)*(1.0 - 1.5*q*q + 0.75*q*q*q));
// }
// // Gradient for cubic spline kernel function
// const Vector SPHStateData::grad_weight(size_t p, const Vector& P) const
//  {
//    Vector X = (P- pos(p))/radius;
// 	double q = X.magnitude();
//    if ((P - pos(p)).magnitude() == 0)
// 		return Vector(1.0, 1.0, 1.0);
// 	else 
// 		return (P-pos(p))*(1.0/(radius*radius*radius*radius*M_PI*(P-pos(p)).magnitude()))*(- 3.0*q + 2.25*q*q);
// }
// const float SPHStateData::weight( size_t p, const Vector& P ) const
// {
//    //Vector X = (pos(p)-P)/radius;
//    Vector X = (P- pos(p))/radius;
//    float Xmag = X.magnitude();
//    if( Xmag >= 2.0 ){ return 0.0; }
//    float scale = 1.0/(3.14159265 * radius*radius*radius);
//    if( Xmag >= 1.0 )
//    {
//       float d = 2.0 - Xmag;
//       return scale*d*d*d/4.0;
//    }
//    return scale*( 1.0 - 1.5*Xmag*Xmag + 0.75*Xmag*Xmag*Xmag );
// }

// const float SPHStateData::weight( size_t p, const Vector& P ) const
// {
//    //Vector X = (pos(p)-P)/radius;
//    Vector X = (P- pos(p))/radius;
//    float Xmag = X.magnitude();
//    if( Xmag >2.0 ){ return 0.0; }
//    float scale = 1.0/(3.14159265 * radius*radius*radius);
//    if( Xmag >1.0 && Xmag <= 2.0 )
//    {
//       float d = 2.0 - Xmag;
//       return scale*d*d*d/4.0;
//    }
//    return scale*( 1.0 - 1.5*Xmag*Xmag + 0.75*Xmag*Xmag*Xmag );
// }

// //smoothing kernel
// const float SPHStateData::weight( size_t p, const Vector& P) const
// {

//     //cubic spline from Smoothed Particle HydrodynamicsTechniques for the Physics Based Simulation of Fluids and Solids
//     //Dan Koschier1, Jan Bender2, Barbara Solenthaler3, and Matthias Teschner4
//     //normalization factors: sigma dimensions=2  40/7 *pi *h^2            d=3 8/pi*h^3
//     //smmothing length == kernel support radius == h, h == 4 * particle radius, h = 0.3 m
//     Vector X = (P - pos(p));
//     float q = (1.0f/radius) * (X.magnitude());
//     float sigma_3 = 8.0f / (3.14159265 * (radius * radius * radius)); 
//     if(q >= 0 && q <= 1.0f/2.0f) 
//     {
//         //std::cout << q << ' ' << sigma_2 * ((6 * (std::pow(q, 3) - std::pow(q, 2))) + 1) << '\n';
//         return sigma_3 * ((6 * ((q * q * q) - (q * q))) + 1);
//     }
//     else if (q > 1.0f/2.0f && q <= 1) 
//     {
//         //std::cout << q << ' ' << sigma_2 * (2 * std::pow((1 - q), 3)) << '\n';
//         return sigma_3 * ( 2 * ((1 - q) * (1 - q) * (1 - q)) );
//     }
//     else 
//     {
//         //std::cout << 0 << '\n';
//         return 0;
//     }

// }

const float SPHStateData::weight(size_t p, const Vector& P) const
{
			float res = 0.0;
         Vector X = (P - pos(p));
         float x = X.magnitude();
         float h3 = radius * radius * radius;
         float m_k = static_cast<float>(8.0) / (M_PI*h3);
			const float q = x / radius;
			if (q <= 1.0)
			{
				if (q <= 0.5)
				{
					const float q2 = q*q;
					const float q3 = q2*q;
					res = m_k * (static_cast<float>(6.0)*q3 - static_cast<float>(6.0)*q2 + static_cast<float>(1.0));
				}
				else
				{
					res = m_k * (static_cast<float>(2.0)*pow(static_cast<float>(1.0) - q, static_cast<float>(3.0)));
				}
			}
			return res;
}

const Vector SPHStateData::grad_weight(size_t p, const Vector& P) const
{
			Vector res; 
         Vector X = (P- pos(p));
         const float x = X.magnitude();
			const float q = x/ radius;
         float h3 = radius * radius * radius;
         float m_l = static_cast<float>(48.0) / (M_PI*h3);
			if ((x > 1.0e-9) && (q <= 1.0))
			{
				Vector gradq = X / x;
				gradq /= radius;
				if (q <= 0.5)
				{
					res = m_l*q*((float) 3.0*q - static_cast<float>(2.0))*gradq;
				}
				else
				{
					const float factor = static_cast<float>(1.0) - q;
					res = m_l*(-factor*factor)*gradq;
				}
			}
			else
				res.set(0,0,0);

			return res;
}

// const Vector SPHStateData::grad_weight( size_t p, const Vector& P ) const
// {
//    Vector value(0,0,0);
//    //Vector X = (pos(p)-P)/radius;
//    Vector X = (P- pos(p))/radius;   
//    float Xmag = X.magnitude();
//    if( Xmag >= 2.0 ){ return value; }
//    if( Xmag == 0.0 ){ return value; }
//    X = X/Xmag;
//    float scale = 1.0/(3.14159265 * radius*radius*radius*radius);
//    if( Xmag >= 1.0 )
//    {
//       float d = 2.0 - Xmag;
//       value = X*( -3.0*scale*d*d/4.0 );
//    }
//    else
//    {
//       value = X*( 3.0*scale*Xmag*Xmag*( 0.75 - Xmag )  );
//    }
//    return value;
// }

// const Vector SPHStateData::grad_weight( size_t p, const Vector& P ) const
// {
//    Vector value(0,0,0);
//    //Vector X = (pos(p)-P)/radius;
//    Vector X = (P- pos(p))/radius;   
//    float Xmag = X.magnitude();
//    if( Xmag > 2.0 ){ return value; }
//    if( Xmag == 0.0 ){ return value; }
//    X = X/Xmag;
//    float scale = 1.0/(3.14159265 * radius*radius*radius*radius);
//    if( Xmag > 1.0 && Xmag <= 2.0)
//    {
//       float d = 2.0 - Xmag;
//       value = X*( -3.0*scale*d*d/4.0 );
//    }
//    else
//    {
//       value = X*( 3.0*scale*Xmag*Xmag*( 0.75 - Xmag )  );
//    }
//    return value;
// }

// const Vector SPHStateData::grad_weight(size_t p, const Vector& P) const
// {
//     //∇W = (∂W/∂x, ∂W/∂y)
//     Vector X = P - pos(p);
//     Vector gradient{};
//     float q = (1.0f/radius) * X.magnitude();
//     float sigma_3 = 8.0f / (3.14159265 * (radius * radius * radius * radius)); 
//     float dx, dy, dz;
//     if(q >= 0 && q <= 1.0f/2.0f) 
//     {
//         dx = (sigma_3 * ( (18 * (q*q)) - (12 * q) )  *  X.X() ) / (radius * X.magnitude()) ;
//         dy = (sigma_3 * ( (18 * (q*q)) - (12 * q) )  *  X.Y() ) / (radius * X.magnitude()) ;
//         dz = (sigma_3 * ( (18 * (q*q)) - (12 * q) )  *  X.Z() ) / (radius * X.magnitude()) ;

//     }
//     else if (q > 1.0f/2.0f && q <= 1)
//     { 
//         dx = (sigma_3 * (-6 * ((1 - q) * (1-q)) )  *  X.X() ) / (radius * X.magnitude()) ;
//         dy = (sigma_3 * (-6 * ((1 - q) * (1-q)) )  * X.Y() ) / (radius * X.magnitude()) ;
//         dz = (sigma_3 * (-6 * ((1 - q) * (1-q)) )  * X.Z() )  / (radius * X.magnitude()) ;

//     }
//     else 
//     {
//         dx = 0;
//         dy = 0;
//         dz = 0;
//     }
//     gradient.set(dx, dy, dz);
//     //std::cout << "gradient: " << gradient.X() << ' ' << gradient.Y() << '\n';
//     //std::cout << gradient.magnitude() << '\n';
//     return gradient;
// }

void SPHStateData::compute_density()
{


   file.open("./densf/dens"+std::to_string(densiter)+".txt");
   file << densiter << '\n';
#pragma omp parallel for
   for( size_t p=0;p<nb();p++)
   {
      float density = 0.0;
      const Vector P = pos(p);
      size_t pindex = index(P);
      size_t i,j,k;
      anti_index(pindex, i,j,k);

      const std::vector<size_t>& cells = cell_contents(i,j,k); 
      for(size_t a=0;a<cells.size();a++)
      {
         size_t pid = cells[a]; 

         density += get_float_attr("volume", pid) * weight(pid,P);
         //density += mass(pid) * weight(pid,P);

         
      }

      std::vector<size_t> neighbors;
      neighbor_cells(i,j,k,neighbors);
      for(size_t a=0;a<neighbors.size();a++)
      {
         size_t pid = neighbors[a];
         density += get_float_attr("volume", pid) * weight(pid,P) ;
         //density += mass(pid) * weight(pid,P);
      }
      density *= density0;
      #pragma omp critical
      file << density << "\n";
      set_attr("density", p, density);
//#pragma omp critical
//      {
//	 std::cout << "Particle " << p << std::endl;
//	 std::cout << "\tposition " << P.X() << " " << P.Y() << " " << P.Z() << std::endl;
//         std::cout << "\tpindex " << pindex << std::endl;
//         std::cout << "\tijk " << i << " " << j << " " << k << std::endl;
//         std::cout << "\tNumber of cells " << cells.size() << std::endl;
//         std::cout << "\tNumber of neighbors " << neighbors.size() << std::endl;
//         std::cout << "\tMass " << mass(p) << std::endl;
//         std::cout << "\tDensity " << density << std::endl;
//      }
   } 
   densiter++;
   file.close();


}

void SPHStateData::compute_predicted_density(size_t p, const double dt)
{

   float density = get_float_attr("density", p);
   float pdensity = 0.0f;
   const Vector P = pos(p);
   const Vector V = vel(p);
   size_t pindex = index(P);
   size_t i,j,k;
   anti_index(pindex, i,j,k);

   const std::vector<size_t>& cells = cell_contents(i,j,k); 
   for(size_t a=0;a<cells.size();a++)
   {
      size_t pid = cells[a]; 
      if(pid != p)
      {
         const Vector V_j = vel(pid);
         pdensity += ((V-V_j) * grad_weight(pid,P));
      }
   }

   std::vector<size_t> neighbors;
   neighbor_cells(i,j,k,neighbors);
   for(size_t a=0;a<neighbors.size();a++)
   {
      size_t pid = neighbors[a];
      const Vector V_j = vel(pid);
      pdensity += ((V-V_j) * grad_weight(pid,P));     
   }

   pdensity *= get_float_attr("volume", p);

   pdensity = (density / density0) + dt *pdensity;
   // #pragma omp critical
   // std::cout << "pdens: " << pdensity << '\n';
   set_attr("predicted_density", p, pdensity);
   


}

void SPHStateData::compute_density_derivative(size_t p)
{


      float density_change = 0.0;
      const Vector P = pos(p);
      const Vector V = vel(p);
      size_t pindex = index(P);
      size_t i,j,k;
      anti_index(pindex, i,j,k);

      const std::vector<size_t>& cells = cell_contents(i,j,k); 
      for(size_t a=0;a<cells.size();a++)
      {
         size_t pid = cells[a]; 
         if(pid != p)
         {
            const Vector V_j = vel(pid);
            density_change += ((V-V_j) * grad_weight(pid,P));
         }
      }

      std::vector<size_t> neighbors;
      neighbor_cells(i,j,k,neighbors);
      for(size_t a=0;a<neighbors.size();a++)
      {
         size_t pid = neighbors[a];
         const Vector V_j = vel(pid);
         density_change += ((V-V_j) * grad_weight(pid,P));     
      }
      density_change *= get_float_attr("volume",p); //all fluid have constant volume

      // density_change = std::max(density_change, 0.0f);
      // #pragma omp critical
      // std::cout << "Dens deriv: " << density_change << '\n';
      set_attr("density_derivative", p, density_change);
   
 
}
//equ. 7 and 8
void SPHStateData::compute_factor()
{
#pragma omp parallel for   
   for( size_t p=0;p<nb();p++)
   {
      float factor = 0;
      float sum_grad_p = 0.0f; // sum of grad pi and pj

      Vector grad_p_i(0,0,0); //pressure gradient of ith particle
      const Vector P = pos(p);
      size_t pindex = index(P);
      size_t i,j,k;
      anti_index(pindex, i,j,k);

      const std::vector<size_t>& cells = cell_contents(i,j,k);
      for(size_t a=0;a<cells.size();a++)
      {         
         size_t pid = cells[a]; 
         if(pid != p)
         {
            //const Vector grad_p_j = mass(pid) * grad_weight(pid, P);
            const Vector grad_p_j = get_float_attr("volume",pid) * grad_weight(pid, P);
            //const Vector grad_p_j = get_float_attr("volume",pid) * grad_weight(pid, P) * density0;

            sum_grad_p += grad_p_j * grad_p_j; // this equals magnitude squared
            grad_p_i += grad_p_j;
         }
      }
      std::vector<size_t> neighbors;
      neighbor_cells(i,j,k,neighbors);
      for(size_t a=0;a<neighbors.size();a++)
      {
         size_t pid = neighbors[a];
         //const Vector grad_p_j = mass(pid) * grad_weight(pid, P);
         const Vector grad_p_j = get_float_attr("volume",pid) * grad_weight(pid, P) ;
         //const Vector grad_p_j = get_float_attr("volume",pid) * grad_weight(pid, P) * density0;

         sum_grad_p += grad_p_j * grad_p_j;
         grad_p_i += grad_p_j;
      }

      sum_grad_p += grad_p_i * grad_p_i;

      sum_grad_p *= density0;

      if(sum_grad_p > m_eps)
      {
         //factor  =  get_float_attr("density",p) / sum_grad_p;
         factor  =  1.0 / sum_grad_p;
      }
      else
         factor = 0.f;
      set_attr("factor", p, factor);
   }
}
// void SPHStateData::compute_factor()
// {
// #pragma omp parallel for   
//    for( size_t p=0;p<nb();p++)
//    {
//       float factor = 0;
//       Vector lsum{};
//       const Vector P = pos(p);
//       size_t pindex = index(P);
//       size_t i,j,k;
//       anti_index(pindex, i,j,k);

//       const std::vector<size_t>& cells = cell_contents(i,j,k);
//       for(size_t a=0;a<cells.size();a++)
//       {         
//          size_t pid = cells[a]; 

//          lsum += (mass(pid) * grad_weight(pid,P));
//       }
//       std::vector<size_t> neighbors;
//       neighbor_cells(i,j,k,neighbors);
//       for(size_t a=0;a<neighbors.size();a++)
//       {
//          size_t pid = neighbors[a];
//          lsum += (mass(pid) * grad_weight(pid,P));

//       }

//       float rsum = 0;
//       const std::vector<size_t>& cells2 = cell_contents(i,j,k);
//       for(size_t a=0;a<cells2.size();a++)
//       {         
//          size_t pid = cells2[a];  
//          rsum += (mass(pid) * grad_weight(pid,P)).magnitude() * 
//                  (mass(pid) * grad_weight(pid,P)).magnitude();
      
//       }
//       std::vector<size_t> neighbors2;
//       neighbor_cells(i,j,k,neighbors2);
//       for(size_t a=0;a<neighbors2.size();a++)
//       {
//          size_t pid = neighbors2[a];
//          rsum += (mass(pid) * grad_weight(pid,P)).magnitude() * 
//                  (mass(pid) * grad_weight(pid,P)).magnitude();
//       }
//       float density = get_float_attr("density", p);
//       float denom = ( (lsum.magnitude() * lsum.magnitude()) + rsum );
//       if(denom > m_eps)
//          factor  =  1.0 / denom;
//       else
//          factor = 0;
//       set_attr("factor", p, factor);
//    }
// }

void SPHStateData::compute_divergence()
{
   #pragma omp parallel for
   for( size_t p=0;p<nb();p++)
   {
      const Vector P = pos(p);
      const Vector V = vel(p);
      size_t pindex = index(P);
      size_t i,j,k;
      anti_index(pindex, i,j,k);

      float divergence = 0.0;
      const std::vector<size_t>& cells = cell_contents(i,j,k);
      for(size_t a=0;a<cells.size();a++)
      {
         size_t pid = cells[a];
	 const Vector VV = vel(pid);
         divergence += mass(pid) * ( (VV-V) * grad_weight(pid,P) )/get_float_attr("density",pid);
      }

      std::vector<size_t> neighbors;
      neighbor_cells(i,j,k,neighbors);
      for(size_t a=0;a<neighbors.size();a++)
      {
         size_t pid = neighbors[a];
	 const Vector VV = vel(pid);
         divergence += mass(pid) * ( (VV-V) * grad_weight(pid,P) )/get_float_attr("density",pid);
      }

      set_attr("divergence", p, divergence);
   } 
}


void SPHStateData::adjust_density_for_divergence(const double dt)
{
   compute_divergence();

#pragma omp parallel for
   for( size_t p=0;p<nb();p++)
   {
      float density = 0.0;
      const Vector P = pos(p);
      size_t pindex = index(P);
      size_t i,j,k;
      anti_index(pindex, i,j,k);

      const std::vector<size_t>& cells = cell_contents(i,j,k);
      for(size_t a=0;a<cells.size();a++)
      {
         size_t pid = cells[a];
         density += mass(pid) * weight(pid,P) * std::exp( get_float_attr("divergence",pid)*dt );
      }

      std::vector<size_t> neighbors;
      neighbor_cells(i,j,k,neighbors);
      for(size_t a=0;a<neighbors.size();a++)
      {
         size_t pid = neighbors[a];
         density += mass(pid) * weight(pid,P) * std::exp( get_float_attr("divergence",pid)*dt );
      }

      set_attr("density", p, density);
   } 




//#pragma omp parallel for
//   for( size_t p=0;p<nb();p++)
//   {
//      float density = get_float_attr("density",p) * std::exp( get_float_attr("divergence",p) * dt );
//      set_attr("density", p, density);
//   }
}

void SPHStateData::populate()
{
   OccupancyVolume::populate(*this);
}


void SPHStateData::set_radius( const float& v )
{
   radius = v;
   particle_radius = radius/4.0f;
   float particle_diamter = particle_radius * 2;
   float mV = 0.8f * particle_diamter * particle_diamter * particle_diamter;
   #pragma omp parall for
   for(int p = 0; p < nb(); p++)
   {
      set_attr("volume", p, mV);

   }
   set_cellsize(2.0*radius); 
}

void SPHStateData::set_maxIter(int maxi)
{
   if(maxi < 2) maxi = 2;
   maxIter = maxi;
}

void SPHStateData::set_mMaxError(int mmx)
{
   if(mmx < 0.01f ) mmx = 0.01f;
   m_maxError = mmx;
}

void SPHStateData::set_maxError(int mx)
{
   if(mx < 0.1f) mx = 0.1f;
   maxError = mx;
}

void SPHStateData::set_density0(float d0)
{
   if(d0 == 0) d0 =1;
   density0 = d0;
}


float SPHStateData::average_density() //const
{
   //file.open("./densavgf/densavg"+std::to_string(densavgiter)+".txt");  
   //file << densavgiter << '\n';
   float density = 0.0;
   for( size_t p=0;p<nb();p++)
   {
      density += get_float_attr("density",p);
   }
   density /= nb();
   //file << density << '\n';
   //file.close();
   densavgiter++;
   return density;
}

float SPHStateData::average_predicted_density() //const
{
   //file.open("./pdensavgf/pdensavg"+std::to_string(pdensavgiter)+".txt");  
   //file << pdensavgiter << '\n';
   float density = 0.0;
   for( size_t p=0;p<nb();p++)
   {
      density += get_float_attr("predicted_density",p);
   }
   density /= nb();
  // file << density << '\n';
   //file.close();
   pdensavgiter++;
   return density;
}

float SPHStateData::average_density_derivative() //const
{
   //file.open("./ddensavgf/ddensavg"+std::to_string(ddensavgiter)+".txt");  
   //file << ddensavgiter << '\n';
   float density = 0.0;
   for( size_t p=0;p<nb();p++)
   {
      density += get_float_attr("density_derivative",p);
   }
   density /= nb();
   //file << density << '\n';
   //file.close();
   ddensavgiter++;
   return density;
}

float SPHStateData::average_divergence() const
{
   float div = 0.0;
   for( size_t p=0;p<nb();p++)
   {
      div += get_float_attr("divergence",p);
   }
   div /= nb();
   return div;
}

float SPHStateData::max_velocity() const
{
   float vel_m = vel(0).magnitude();
   for(size_t i = 1; i < nb(); i++)
   {
      float curr_vel = vel(i).magnitude();
      if (curr_vel > vel_m) vel_m = curr_vel;
   }
   return vel_m;
}
