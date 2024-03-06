

#include "CollisionPlane.h"



using namespace pba;


CollisionPlane::CollisionPlane( const Vector& norm, const Vector& p0 ) :
  normal(norm),
  P0(p0)
  {

  }



   
bool CollisionPlane::hit( const Vector& X0, const Vector& Xu, const Vector& V, const double dt, Vector& XH_cand, double& dtH_cand, float rad) const
{
    bool hit = false;
    float fxu = normal * ((Xu) - P0);
    float fx0 = normal * ((X0) - P0);

    // if(fxu >=0)
    // {
    //     fxu = (Xu - P0) * normal - rad;
    // }
    // else fxu = (Xu - P0) * normal + rad;

    // if(fx0 >=0)
    // {
    //     fx0 = (X0 - P0) * normal - rad;
    // }
    // else fx0 = (X0 - P0) * normal + rad;

    if((fxu == 0 || fx0 * fxu < 0 || (fxu <=0 && fx0 <=0))) hit = true;
    if(hit)
    {
        //std::cout << V.X() << ' ' << V.Y() << ' ' << V.Z() << '\n';
        if(V.isZero() == true)
        {
            std::cout << "Zero? " << V.isZero() << '\n';
            //float epsilon = 0.06041;
            XH_cand = X0 + (normal * rad); //- (normal*epsilon);
            dtH_cand = dt;
            return hit;
        }       
        // if(V.isZero())
        // {
        //    std::cout << V.X() << ' ' << V.Y() << ' ' << V.Z() << '\n';
        //    std::cout << "Zero v: "<< normal*V << '\n';
        // } 
        XH_cand = X0 + V * ( (normal * (P0 - X0)) / (normal * V) );
        dtH_cand = (normal * (P0-X0)) / (normal * V);
      
        //std::cout << "Dth: " << dtH_cand << '\n';
        if(dtH_cand > dt) std::cout << "dth bigger than dt!!?!!!!!??\n";
    }

    return hit;
}


void CollisionPlane::handle( const Vector& XS, const Vector& VS, const double& dt, const Vector& XH, const double& dtH, Vector& XR, Vector& VR, float cs, float cr ) const
{

    VR = (cs * VS) - ((cs + cr) * normal * (normal * VS));
    XR = XH + VR * (dt - dtH);
}



