#include "CollisionGeometryLibrary.h"

namespace pba
{

pba::CollisionSurface GenerateCollisionCube()
{
    CollisionSurface s = makeCollisionSurface();

    CollisionPlane bot(Vector(0,1,0), Vector(0,-2,0));
    CollisionPlane top(Vector(0,-1,0), Vector(0,2,0));
    CollisionPlane left(Vector(1,0,0), Vector(-2,0,0));
    CollisionPlane right(Vector(-1,0,0), Vector(2,0,0));
    CollisionPlane front(Vector(0,0,-1), Vector(0,0,2));
    CollisionPlane back(Vector(0,0,1), Vector(0,0,-2));



    s->addPlane(bot);
    s->addPlane(top);
    s->addPlane(left);
    s->addPlane(right);
    s->addPlane(front);
    s->addPlane(back);



    return s;
}


}