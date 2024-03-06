#include "CollisionHandler.h"
#include "LinearAlgebra.h"
#include <cmath>
#include <iostream>


using namespace pba;


void CollisionHandler::set_collision_surface( CollisionSurface& c )
{
   surf = c;
//    if(tree){ delete tree; }
//    tree = new TraceTree( c->aabb().LLC(), c->aabb().URC(), 0, 20, 1 );
//    tree->addObject(surf);
//    tree->Divide();
}






// Check if positions imply a collision has already taken place within the allotted time.
// If so, backs the position up along the velocity direction to the point of impact, then
// does an elastic bounce
void ElasticCollisionHandler::handle_collisions( const double dt, DynamicalState& PQ ) 
{
    for( size_t i=0;i<PQ->nb();i++ )
    {
        Vector V0 = PQ->vel(i);
        Vector X0 = PQ->pos(i) - V0*dt;
        Vector XR = PQ->pos(i);
        Vector VR = V0;
        float radius = PQ->get_float_attr("radius", i);
        CollisionData data;
        double running_dt = dt;
        bool keep_checking_for_hits = true;
        int count = 0;
        while(keep_checking_for_hits && !(running_dt <= 0))
        {
            keep_checking_for_hits = false;
            if( surf->hit( X0, XR, V0, running_dt, data , radius))
            {
                keep_checking_for_hits = true;
                // Handle collision on plane with index pH, with the smalled dtH at hit point XH
                if(data.hit_plane)
                    surf->get_plane(data.hit_index).handle( X0, V0, dt, data.XH, data.hit_time, XR, VR, surf->coeff_sticky(), surf->coeff_restitution());
                X0 = data.XH;
                V0 = VR;
                running_dt = dt - data.hit_time;
                //std::cout << "hit " << count << " : "<< running_dt << '\n';
                count++;
                if( running_dt == 0.0 ){ keep_checking_for_hits = false; }
                //break;
            }
        }
        PQ->set_pos( i, XR );
        PQ->set_vel( i, VR );
    }
}
