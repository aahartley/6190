//*******************************************************************
//
//   CollisionTriangle.C
//
// Class for a triangle for collisions
//
//  Copyright (c) 2017 Jerry Tessendorf
//
//
//*******************************************************************

#include "CollisionTriangle.h"
#include "LinearAlgebra.h"
#include <cmath>
#include <iostream>

using namespace pba;


CollisionTriangleRaw::CollisionTriangleRaw( const Vector& p0, const Vector& p1, const Vector& p2  ) :
  P0(p0),
  P1(p1),
  P2(p2),
  color(Color(0.5,0.5,0.5,0)),
  hitcolor(Color(0,0,0,0)),
  litcolor(Color(0.75,0.75,0.75,0)),
  decayrate(0.01),
  visible(true)
    {
    e1 = p1-p0;
    e2 = p2-p0;
    normal = (e1^e2).unitvector();
    det = (e1*e1)*(e2*e2) - (e1*e2)*(e1*e2);

    double xllc = p0.X();
    double yllc = p0.Y();
    double zllc = p0.Z();
    if( p1.X() < xllc ){ xllc = p1.X(); }
    if( p2.X() < xllc ){ xllc = p2.X(); }
    if( p1.Y() < yllc ){ yllc = p1.Y(); }
    if( p2.Y() < yllc ){ yllc = p2.Y(); }
    if( p1.Z() < zllc ){ zllc = p1.Z(); }
    if( p2.Z() < zllc ){ zllc = p2.Z(); }

    double xurc = p0.X();
    double yurc = p0.Y();
    double zurc = p0.Z();
    if( p1.X() > xurc ){ xurc = p1.X(); }
    if( p2.X() > xurc ){ xurc = p2.X(); }
    if( p1.Y() > yurc ){ yurc = p1.Y(); }
    if( p2.Y() > yurc ){ yurc = p2.Y(); }
    if( p1.Z() > zurc ){ zurc = p1.Z(); }
    if( p2.Z() > zurc ){ zurc = p2.Z(); }


  }



   
bool CollisionTriangleRaw::hit( const Vector& P, const Vector& V, const double tmax, double& t ) 
{
   double f0 = normal*(P-P0);
   if( f0 == 0.0 ){ t = 0.0; return true; }
   double f1 = normal*(P-P0-V*tmax);
   if( f0*f1 > 0.0 ){ return false; }
   t = f0 / (normal*V);


   if( t*tmax < 0.0 ){ return false; }
   if( std::fabs(t) >= std::fabs(tmax) ) { return false; }
   if( (tmax-t)/tmax < 1.0e-6 ) { return false; }

   Vector X = P - V*t - P0;
   return is_in_triangle(X);
}


bool CollisionTriangleRaw::is_in_triangle( const Vector& X )
{
   double e1X = X*e1;
   double e2X = X*e2;

   double a = (e1X*(e2*e2) - e2X*(e1*e2))/det;
   if( a < 0.0 ){ return false; }
   if( a > 1.0 ){ return false; }

   double b = (e2X*(e1*e1) - e1X*(e1*e2))/det;
   if( b < 0.0 ){ return false; }
   if( b > 1.0 ){ return false; }
   if( a+b > 1.0 ){ return false; }
   hitcolor = litcolor;
   return true;
}






bool CollisionTriangleRaw::hit( const Vector& P, const Vector& V, const double R, const double tmax, double& t )
{
    double sp = normal*(P-P0);
    if( std::fabs(sp) >= R ){ return false; }
    double nn = normal*V;
    double t0 = sp / nn;
    double t1 = t0 + R/nn;
    if( std::isnan(t1) ){ return false; }
    double t2 = t0 - R/nn;
    if( std::isnan(t2) ){ return false; }

    t = 0.0;

// First, check what it takes to move out of the plane

    bool case1 = false;
    if( (tmax*t1 > 0.0) && (std::fabs(t1) <= std::fabs(tmax)) )
    {
       Vector X = P - V*t1 - P0;
       if( (X*normal)*sp >= 0.0 )
       {
       double e1X = X*e1;
       double e2X = X*e2;
       double a = (e1X*(e2*e2) - e2X*(e1*e2))/det;
       if( a >= 0.0 && a <= 1.0 )
       {
             double b = (e2X*(e1*e1) - e1X*(e1*e2))/det;
             if( b >= 0.0 && b <= 1.0 && a+b <= 1.0 )
             {
                case1 = true;
             }
       }
       }
    }

    bool case2 = false;
    if( (tmax*t2 > 0.0) && (std::fabs(t2) <= std::fabs(tmax)) )
    {
       Vector X = P - V*t2 - P0;
       if( (X*normal)*sp >= 0.0 )
       {
       double e1X = X*e1;
       double e2X = X*e2;
       double a = (e1X*(e2*e2) - e2X*(e1*e2))/det;
       if( a >= 0.0 && a <= 1.0 )
       {
             double b = (e2X*(e1*e1) - e1X*(e1*e2))/det;
             if( b >= 0.0 && b <= 1.0 && a+b <= 1.0 )
             {
                case2 = true;
             }
       }
       }
    }


    if( case1 && case2 )
    {
       t = (std::fabs(t1) > std::fabs(t2) ) ? t1 : t2;
    }
    else if(case1)
    {
       t = t1;
    }
    else if(case2)
    {
       t = t2;
    }

// Second, check what it takes to move enough to have one or more edges be tangent

       // edge e1
       Vector Sperp = P-P0;
       double q = Sperp*e1/(e1*e1);
       if( q >= 0.0 && q <= 1.0 )
       {
          Vector Vperp = V - ((e1*V)/(e1*e1))*e1;
          Sperp = Sperp - q*e1;
          double radical = (Vperp*Sperp)*(Vperp*Sperp) + (Vperp*Vperp)*(Sperp*Sperp-R*R);
          if( radical >= 0.0 )
          {
             radical = std::sqrt(radical)/(Vperp*Vperp);
             double te0 = -(Vperp*Sperp)/(Vperp*Vperp);
             double tep = te0 + radical;
             double tem = te0 - radical;
             if( tep*tmax > 0.0 && std::fabs(tep) < std::fabs(tmax) )
             { 
                t = ( std::fabs(t) > std::fabs(tep) ) ? t : tep;
             } 
             if ( tem*tmax > 0.0 && std::fabs(tem) < std::fabs(tmax) )
             {
                t = ( std::fabs(t) > std::fabs(tem) ) ? t : tem;
             }
          }
       }

       // edge e2
       Sperp = P-P0;
       q = Sperp*e2/(e2*e2);
       if( q >= 0.0 && q <= 1.0 )
       {
          Vector Vperp = V - ((e2*V)/(e2*e2))*e2;
          Sperp = Sperp - q*e2;
          double radical = (Vperp*Sperp)*(Vperp*Sperp) + (Vperp*Vperp)*(Sperp*Sperp-R*R);
          if( radical >= 0.0 )
          {
             radical = std::sqrt(radical)/(Vperp*Vperp);
             double te0 = -(Vperp*Sperp)/(Vperp*Vperp);
             double tep = te0 + radical;
             double tem = te0 - radical;
             if( tep*tmax > 0.0 && std::fabs(tep) < std::fabs(tmax) )
             { 
                t = ( std::fabs(t) > std::fabs(tep) ) ? t : tep;
             } 
             if ( tem*tmax > 0.0 && std::fabs(tem) < std::fabs(tmax) )
             {
                t = ( std::fabs(t) > std::fabs(tem) ) ? t : tem;
             }
          }
       }

       // edge e3
       Vector e3 = P2-P1;
       Sperp = P-P1;
       q = Sperp*e3/(e3*e3);
       if( q >= 0.0 && q <= 1.0 )
       {
          Vector Vperp = V - ((e3*V)/(e3*e3))*e3;
          Sperp = Sperp - q*e3;
          double radical = (Vperp*Sperp)*(Vperp*Sperp) + (Vperp*Vperp)*(Sperp*Sperp-R*R);
          if( radical >= 0.0 )
          {
             radical = std::sqrt(radical)/(Vperp*Vperp);
             double te0 = -(Vperp*Sperp)/(Vperp*Vperp);
             double tep = te0 + radical;
             double tem = te0 - radical;
             if( tep*tmax > 0.0 && std::fabs(tep) < std::fabs(tmax) )
             { 
                t = ( std::fabs(t) > std::fabs(tep) ) ? t : tep;
             } 
             if ( tem*tmax > 0.0 && std::fabs(tem) < std::fabs(tmax) )
             {
                t = ( std::fabs(t) > std::fabs(tem) ) ? t : tem;
             }
          }
       }


// Check again each point:


       Sperp = P-P0;
       double radical = (V*Sperp)*(V*Sperp) + (V*V)*(R*R - Sperp*Sperp);
       if( radical >= 0.0 )
       {
             radical = std::sqrt(radical)/(V*V);
             double te0 = -(V*Sperp)/(V*V);
             double tep = te0 + radical;
             double tem = te0 - radical;
             if( tep*tmax > 0.0 && std::fabs(tep) < std::fabs(tmax) )
             { 
                t = ( std::fabs(t) > std::fabs(tep) ) ? t : tep;
             } 
             if ( tem*tmax > 0.0 && std::fabs(tem) < std::fabs(tmax) )
             {
                t = ( std::fabs(t) > std::fabs(tem) ) ? t : tem;
             }
       }

       Sperp = P-P1;
       radical = (V*Sperp)*(V*Sperp) + (V*V)*(R*R - Sperp*Sperp);
       if( radical >= 0.0 )
       {
             radical = std::sqrt(radical)/(V*V);
             double te0 = -(V*Sperp)/(V*V);
             double tep = te0 + radical;
             double tem = te0 - radical;
             if( tep*tmax > 0.0 && std::fabs(tep) < std::fabs(tmax) )
             { 
                t = ( std::fabs(t) > std::fabs(tep) ) ? t : tep;
             } 
             if ( tem*tmax > 0.0 && std::fabs(tem) < std::fabs(tmax) )
             {
                t = ( std::fabs(t) > std::fabs(tem) ) ? t : tem;
             }
       }


       Sperp = P-P2;
       radical = (V*Sperp)*(V*Sperp) + (V*V)*(R*R - Sperp*Sperp);
       if( radical >= 0.0 )
       {
             radical = std::sqrt(radical)/(V*V);
             double te0 = -(V*Sperp)/(V*V);
             double tep = te0 + radical;
             double tem = te0 - radical;
             if( tep*tmax > 0.0 && std::fabs(tep) < std::fabs(tmax) )
             { 
                t = ( std::fabs(t) > std::fabs(tep) ) ? t : tep;
             } 
             if ( tem*tmax > 0.0 && std::fabs(tem) < std::fabs(tmax) )
             {
                t = ( std::fabs(t) > std::fabs(tem) ) ? t : tem;
             }
       }



       if( t == 0.0 ){ return false; }


/*
   if( tmax*t1 > 0.0 && tmax*t2 > 0.0 )
   {
      if(  std::fabs(t1) < std::fabs(tmax) && std::fabs((t1-tmax)/tmax) >= 0.000001 )
      {
         t = t1;
      }
      if(  std::fabs(t2) < std::fabs(tmax) && std::fabs((t2-tmax)/tmax) >= 0.000001 )
      {
         t = ( std::fabs(t) > std::fabs(t2) ) ? t : t2;
      }
   }
   else if( (tmax*t1 > 0.0) && (std::fabs(t1) < std::fabs(tmax)) &&  std::fabs((t1-tmax)/tmax) >= 0.000001 )
   {
      t = t1;
   }
   else if( tmax*t2 > 0.0 && std::fabs(t2) < std::fabs(tmax) && std::fabs((t2-tmax)/tmax) >= 0.000001 )
   {
      t = t2;
   }
   else { return false; }

   Vector X = P - V*t - P0;

   double e1X = X*e1;
   double e2X = X*e2;

   double a = (e1X*(e2*e2) - e2X*(e1*e2))/det;
   if( a < 0.0 ){ return false; }
   if( a > 1.0 ){ return false; }

   double b = (e2X*(e1*e1) - e1X*(e1*e2))/det;
   if( b < 0.0 ){ return false; }
   if( b > 1.0 ){ return false; }
   if( a+b > 1.0 ){ return false; }
*/

   hitcolor = litcolor;
   return true;
}








void CollisionTriangleRaw::set_orientation( const Vector& P )
{
   if( (P-P0)*normal < 0.0 )
   {
      Vector temp = e2;
      e2 = e1;
      e1 = temp;
      temp = P2;
      P2 = P1;
      P1 = temp;
      normal = -normal; 
   }
}



void CollisionTriangleRaw::decay(){ hitcolor = hitcolor*std::exp(-decayrate); }

CollisionTriangle pba::makeCollisionTriangle(  const Vector& p0, const Vector& p1, const Vector& p2 )
{
   return CollisionTriangle( new CollisionTriangleRaw(p0,p1,p2) );
}

