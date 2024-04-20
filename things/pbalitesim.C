

#include "PbaViewer.h"
#include "ScreenCapturePPM.h"
#include "MyThing.h"
#include "MyThing1.h"
#include "DFSPHThing.h"
#include "CollisionSurface.h"
#include "CollisionGeometryLibrary.h"


int main(int argc, char** argv)
{
   // Set up command line arguments, if any
   std::vector<std::string> args;
   for(int i=0;i<argc;i++)
   {
      args.push_back( argv[i] );
   }


   // Instantiate a viewer
   pba::PbaViewer* viewer = pba::CreateViewer();

  
   // Set up a simulation thing
   pba::PbaThing mything = pba::CreateDFSPHThing();
   //pba::PbaThing mything = pba::CreateMyThing1();

   pba::CollisionSurface cube = pba::GenerateCollisionCube();

   mything->AddCollisionSurface(cube);
   viewer->AddThing(mything);

   //Optional thing to capture the screen each window and write it to a file
   pba::PbaThing imager = pba::CreateScreenCapturePPM("./" + mything->Name());
   viewer->AddThing(imager);



   // Initialize viewer
   viewer->Init(args);
   // Run the (GLUT) main loop
   viewer->MainLoop();

}
