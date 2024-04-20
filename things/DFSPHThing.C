
#include "DFSPHThing.h"
#include <cstdlib>
#include <GL/gl.h>   // OpenGL itself.
#include <GL/glu.h>  // GLU support library.
#include <GL/glut.h> // GLUT support library.
#include <iostream>



using namespace std;

using namespace pba;





DFSPHThing::DFSPHThing(const std::string nam) :
 PbaThingyDingy (nam),
 emit       (false)
{
    Reset();
    std::cout << name << " constructed\n";

    box = makeCollisionSurface();
    AABB bounds (Vector(-2,-2,-2), Vector(2,2,2));
    float h = 0.3;
    particles = CreateSPH(bounds, h , "DFsph");
    // particles->add(1);
    // particles->set_pos(0,Vector(0,0,0));
    // particles->set_vel(0,Vector(0,0,0));
    // particles->set_ci(0,Color(1,1,0,1));
    // particles->create_attr("radius",0.075f);
    pba::ParticleEmitter emitter(Vector(0,0,0),Vector(0,0,0), 1,1);
    emitter.emitCube(particles, 10, 6, Vector(0,1,0));
    std::cout << "Emit: Total Points " << particles->nb() << std::endl;
    force = CreateAccumulatingForce();
    gravityforce = CreateAccumulatingGravityForce(Vector(0,-9.81f,0));
    sphforce = CreateDFSPHForce();

    std::shared_ptr<AccumulatingForce> f = dynamic_pointer_cast<AccumulatingForce>(force); 
	f->add(gravityforce);
    f->add(sphforce);
    // GISolver a = CreateAdvancePosition(particles, collisions);
    // GISolver b = CreateAdvanceVelocity(particles, force);
    //solver = CreateForwardEulerSolver(a, b);
    solver = CreateDFSPHSolver(particles, force, 1.0e+06, 1.0e+06, *dynamic_pointer_cast<DFSPHForce>(sphforce), collisions);
    solver->init();
}

DFSPHThing::~DFSPHThing(){}

void DFSPHThing::Init( const std::vector<std::string>& args ) {SetSimulationTimestep(1.0/48.0);}
    
void DFSPHThing::Display() 
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
    glColor3f(0,1,1);
    glBegin(GL_QUADS);
    glVertex3f(-2,-2,-2);
    glVertex3f(2,-2,-2);
    glVertex3f(2,-2,2);
    glVertex3f(-2,-2,2);
    glEnd();
   glPointSize(5.0);
   glBegin(GL_POINTS);
   for(size_t i = 0; i < particles->nb(); i++)
   {
    const Color& ci = particles->ci(i);
    const pba::Vector& v = particles->pos(i);
    //glPushMatrix();
    glColor3f( ci.red(), ci.green(), ci.blue() );
    //glTranslatef(v.X(), v.Y(),v.Z());
    //glutSolidSphere(0.075, 30,30);
    //glPopMatrix();
    glVertex3f( v.X(), v.Y(), v.Z() );
   }
   
   glEnd();
}

void DFSPHThing::Keyboard( unsigned char key, int x, int y )
{
       PbaThingyDingy::Keyboard(key,x,y);
       if( key == 'e' ){ emit = !emit; }
}


void DFSPHThing::solve()
{
    solver->solve(dt);
}

void DFSPHThing::Reset()
{

}

void DFSPHThing::Usage()
{
   PbaThingyDingy::Usage();
   cout << "=== " << name << " ===\n";
   cout << "e            toggle particle emission on/off\n";
}

void DFSPHThing::AddCollisionSurface(pba::CollisionSurface& s)
{
    std::cout << "Add CollisionSurface\n";
    box = s;
    s->set_coeff_restitution(0.5);
    //s->set_coeff_sticky(0.1);
    collisions.set_collision_surface(box);
}

pba::PbaThing pba::CreateDFSPHThing(){ return PbaThing( new DFSPHThing() ); }


