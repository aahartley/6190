
#include "MyThing1.h"
#include <cstdlib>
#include <GL/gl.h>   // OpenGL itself.
#include <GL/glu.h>  // GLU support library.
#include <GL/glut.h> // GLUT support library.
#include <iostream>



using namespace std;

using namespace pba;





MyThing1::MyThing1(const std::string nam) :
 PbaThingyDingy (nam),
 emit       (false)
{
    Reset();
    std::cout << name << " constructed\n";

    box = makeCollisionSurface();

    particles = CreateDynamicalState("ParticleState");
    particles->add(1);
    std::cout << "Emit: Total Points " << particles->nb() << std::endl;
    particles->set_pos(0,Vector(0,0,0));
    particles->set_vel(0,Vector(0,0,0));
    particles->set_ci(0,Color(1,1,0,1));
    particles->create_attr("radius",0.075f);
    force = CreateAccumulatingForce();
    gravityforce = CreateAccumulatingGravityForce(Vector(0,-9.81f,0));

    std::shared_ptr<AccumulatingForce> f = dynamic_pointer_cast<AccumulatingForce>(force); 
	f->add(gravityforce);
    GISolver a = CreateAdvancePosition(particles, collisions);
    GISolver b = CreateAdvanceVelocity(particles, force);
    solver = CreateForwardEulerSolver(a, b);

}

MyThing1::~MyThing1(){}

void MyThing1::Init( const std::vector<std::string>& args ) {}
    
void MyThing1::Display() 
{
    //    glBegin(GL_LINES);

    // for(int i = 0; i < box->plane_size(); i ++)
    // {
    //     CollisionPlane plan = box->get_plane(i);
    //     Vector v = plan.getP0();
    //     glVertex3f(v.X(), v.Y(), v.Z());
    //     glVertex3f(-v.X(), -v.Y(), -v.Z());

    // }
    // glEnd();
    // glBegin(GL_LINES);
    // glVertex3f(-2,-2,-2);
    // glVertex3f(2,-2,-2);
    // glVertex3f(-2,-2,2);
    // glVertex3f(2,-2,2);
    // glEnd();
    glColor3f(0,0,1);
    glBegin(GL_QUADS);
    glVertex3f(-2,-2,-2);
    glVertex3f(2,-2,-2);
    glVertex3f(2,-2,2);
    glVertex3f(-2,-2,2);
    glEnd();
   //glPointSize(5.0);
   //glBegin(GL_POINTS);
   for(size_t i = 0; i < particles->nb(); i++)
   {
    const Color& ci = particles->ci(i);
    const pba::Vector& v = particles->pos(i);
    glPushMatrix();
    glColor3f( ci.red(), ci.green(), ci.blue() );
    glTranslatef(v.X(), v.Y(),v.Z());
    glutSolidSphere(0.075, 30,30);
    glPopMatrix();
    //glVertex3f( v.X(), v.Y(), v.Z() );
   }
   
   //glEnd();
}

void MyThing1::Keyboard( unsigned char key, int x, int y )
{
       PbaThingyDingy::Keyboard(key,x,y);
       if( key == 'e' ){ emit = !emit; }
}


void MyThing1::solve()
{
    solver->solve(dt);
}

void MyThing1::Reset()
{

}

void MyThing1::Usage()
{
   PbaThingyDingy::Usage();
   cout << "=== " << name << " ===\n";
   cout << "e            toggle particle emission on/off\n";
}

void MyThing1::AddCollisionSurface(pba::CollisionSurface& s)
{
    std::cout << "Add CollisionSurface\n";
    box = s;
    s->set_coeff_restitution(0.5);
    //s->set_coeff_sticky(0.1);
    collisions.set_collision_surface(box);
}

pba::PbaThing pba::CreateMyThing1(){ return PbaThing( new MyThing1() ); }


