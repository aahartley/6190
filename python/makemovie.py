#!/usr/bin/python3


#
#
#
#  This is a script to convert ppm files to
#  jpg files, then a movie, then clean up the
#  jpg files.  Assuming you are created a 
#  collection of jpg files with names
#
#      some_base_name.XXXX.ppm
#
#  where XXXX is a zero-padded four digit frame
#  number, you can use this tool like this:
#
#       python ./makemovie.py some_base_name.*.ppm
#
#  when done, a file "test.mov" will be generated.
#
#  The ImageMagick routine "convert" is used in this script.
#
#





import os
import sys

if len(sys.argv) > 1:
    if "-h" in sys.argv:
        print( "makemovie.py usage:" )
        print( "\tmakemovie.py <list of files to convert to a movie>" )
        print( "Wild cards make be used in the list." )
        print( "The output movie will be named \"test.mp4\"." )
        print( "Example:" )
        print( "\t./makemovie.py *.ppm" )
        print( " " )
        print( "Before generating the movie, a collection of jpeg conversion" )
        print( "files are created using the ImageMagick tool \"convert\"." )
        print( "After the movie is generated, these files are removed." )
        print( " " )
        sys.exit()
        
frame = 1
for img in sys.argv[1:]:
    print( str(frame) + " " + img )
    padframe = str(frame)
    if frame < 1000:
        padframe = "0" + padframe
    if frame < 100:
        padframe = "0" + padframe
    if frame < 10:
        padframe = "0" + padframe
    outimage = "conversion." + padframe + ".jpg"
    cmd = "convert " + img + " " + outimage
    os.system(cmd)
    frame = frame + 1

if len(sys.argv[1:]) > 1:
    cmd = "ffmpeg -r 24 -i conversion.%04d.jpg -r 24 test.mp4"
    os.system(cmd)
    cmd = "rm -rf conversion.????.jpg"
    os.system(cmd)
    #cmd = "ffplay -loop 0 test.mov"
    #os.system(cmd)
else:
    print( "Not enough frames found. No movie generated." )

