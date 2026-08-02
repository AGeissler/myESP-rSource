#!/usr/bin/python3

# Version 1.1
# Initial version with only basic python3 facilities used.

# Date 12 September 2024
# What it does:
# This python3 script tests prj facilities by creating a simple model with two zones
# via sequences of commands. It proceeds in stages with user confirmation for each step.
# It tests the following facilities:
#   initial model registration,
#   creation of a rectangular zone (4mx3mx2.7m),
#   naming and attribution of surfaces,
#   addition of a window frame glazing and a door,
#   creation of a zone construction and operations file
#   and a model contents report.
#   duplicate the zone to the right, change partition construction and make connection.

# It generates shell scripts to do the work and offers to run them.
# Each function involves a fresh invocation of prj and each finishes by working
# back up the menu structure of prj.  It borrows ideas from John Allison's
# espy.

# USE:
# python3 ./create_test_model_2zone_scripts.py -r <root_name> -s <synopsis> -l <latitude> -g <longitude-diff> -a <altitude> -y <year>
# A -v in the command line echos the command sequence generated for each task.
# A -i in the command line switches to incremental processing asking the user for confirmation at each step.

# TODO: There is scope for using dictionaries to capture common keystroke sequences.

import sys, getopt, os
from subprocess import PIPE, run

verbose = False

def ask_YorN(s_desc):
# Part of APUF (Andy's Python Useful Functions).
# Ask yes or no question and return a logical.
# Argument is a string containing the question, the function automatically adds syntax guidance onto the end of this.
    while True:
        s_YorN=input(s_desc+' Enter "y" or "n", or "e" to exit: ')
        if s_YorN=='y' or s_YorN=='Y':
            return True
        elif s_YorN=='n' or s_YorN=='N':
            return False
        elif s_YorN=='e' or s_YorN=='E':
            sys.exit('Exiting.')
        elif s_YorN=='42':
            print('Nice reference you hoopy frood, but try again.')
        else:
            print('Unrecognised input, please try again.')

# Function to create a rectangular zone, name the surfaces and attribute
# MLC and boundary conditions.
def create_box_room(root,cfg_file,verbose):
   doc = "a simple room with frame window and door"
   origin = "0.0 0.0 0.0"
   plan = "4.0 3.0 2.7"
   orientation = "0.0"
   elevation = "0.0"

# Command to run within prj to build the initial shape of the zone.
   initial_cmd = [ "m", "c", "a", "a", "box", doc, "a", origin, plan, orientation, elevation ]

# Command to attribute all of the surface names.
   surf_name_cmd = [ "f", "*", "a", "*", "-", "front", "right", "back", "left", "ceiling", "floor" ]

# Command to attribute composition of wall surfaces where 'k' is cavity brick..
   wall_mlc_cmd = [ "*", "b", "a", "k", "a", "b", "c", "d", "-" ]

# Command to attribute composition of ceiling.
   ceil_mlc_cmd = [ "e", "e", "f", "f", "-" ]

# Command to attribute composition of floor where 'd' is entry_floor.
   floor_mlc_cmd = [ "f", "e", "h", "d", "-" ]

# Command to set all but the floor to ambient and then save the zone.
   set_ambient_cmd = [ "*", "c", "a", "-", "*", "f", "-", "f", "f", "f", "-", "1", "-", "-", "!", "", "y" ]
   set_exit_cmd = [ "-", "-", "-", "-", "-", "-" ]

# Join up the commands.
   cmd = initial_cmd + surf_name_cmd + wall_mlc_cmd + ceil_mlc_cmd + floor_mlc_cmd + set_ambient_cmd + set_exit_cmd
   cmd = "\n".join(cmd)
   if (verbose == True): print (cmd)
   fn = root+"_initial_steps.sh"
   f = open(fn, "a")
   f.write(cmd)
   f.write("\nXXX")
   f.close()
   msgos = "chmod a+x " + fn
   os.system(msgos)

# Change to the cfg folder.
#    os.chdir(root)
#    os.chdir("cfg")
   print (os.getcwd())
   print ('the cfg is '+cfg_file)
   print()

   print ('Initial rectangular zone script added')
#    return prj2

# Function to add a window to front wall of rectangular zone as well
# as a door in the back wall.
def add_window_door(root,cfg_file,verbose):
   frame_name = "s_frame"
   geometry = "1.0 0.9 2.0 1.5"
   glass_name = "s_glaze"
   door_name = "back_door"
   door_geometry = "0.2 0.8 2.1"

# Command to run within prj to get to the zone and the add surface option for.
   initial_cmd = [ "m", "c", "a", "a", "e", "+", "c", "a", "a", geometry, frame_name, "e", "i", "b", "-", "c", "y" ]

# Create the glazing.
   add_glazing = [ "+", "c", "g", "e", "0.1", glass_name, "d", "a", "e", "-", "b", "y" ]

# Create the door into the back (3rd) surface with door MLC and type undercut..
   add_door = [ "+", "c", "c", "b", door_geometry, door_name, "c", "a", "a", "-", "b", "y", "-" ]

# Save and work back up the menu structures.
   exit_cmd = [ "!", "", "y", "-", "-", "-", "!", "", "-", "-" ]
   cmd = initial_cmd + add_glazing + add_door + exit_cmd
   cmd = "\n".join(cmd)
   if verbose == True: print (cmd)
   fn = root+"_add_window_door.sh"
   f = open(fn, "w")
   f.write("#!/bin/sh")
   f.write("\n# adds window and door and attributes.")
   msgcd = "\ncd " + root + '/cfg'
   f.write(msgcd)
   msg = "\nprj -mode text -file " + cfg_file + "<<XXX\n"
   f.write(msg)
   f.write(cmd)
   f.write("\nXXX")
   f.close()
   msgos = "chmod a+x " + fn
   os.system(msgos)

   fdn = root+"_do_all.sh"
   f = open(fdn, "a")
   f.write("\necho "  "")
   msg = "\n./"+fn
   f.write(msg)
   f.write("\necho "  "\n")
   f.close()
   print ('Window frame and door script added')

# Function to add an overhang at top edge of sourth facade.
def add_overhang(root,cfg_file,verbose):
   ovh_name = "overhang"
   origin = "0.0 -0.5 2.6"
   blk_geo = "4.0 0.45 0.1"

# Commands to run within prj to get to the zone and the add the overhang.
   initial_cmd = [ "m", "c", "a", "a", "h", "a", "*", "a", "a", "a", origin, "b", blk_geo, "g", "e", "k", "f", ovh_name, ]
   second_cmd = [ "-", ">", "b", "a", "c", "a", "-", "-", "!", "", "y", "-", "-", "f", "c", "-", "!", "", "-", "-", ]

   cmd = initial_cmd + second_cmd
   cmd = "\n".join(cmd)
   if verbose == True: print (cmd)
   fn = root+"_add_overhang.sh"
   f = open(fn, "w")
   f.write("#!/bin/sh")
   f.write("\n# adds overhang at south facade.")
   msgcd = "\ncd " + root + '/cfg'
   f.write(msgcd)
   msg = "\nprj -mode text -file " + cfg_file + "<<XXX\n"
   f.write(msg)
   f.write(cmd)
   f.write("\nXXX")
   f.close()
   msgos = "chmod a+x " + fn
   os.system(msgos)

   fdn = root+"_do_all.sh"
   f = open(fdn, "a")
   f.write("\necho "  "")
   msg = "\n./"+fn
   f.write(msg)
   f.write("\necho "  "\n")
   f.close()
   print ('Overhang added')

# Function to create initial zone construction file.
def create_zone_construction_file(root,cfg_file,verbose):

# Command to run within prj to get to menu for zone constructions and
# then, since all attributes are known, save the constructions file.
   initial_cmd = [ "m", "c", "b", "a", "b", "b", ">", "", "-", "-", "-" ]

# Save cfg and work back up the menu structures.
   exit_cmd = [ "!", "", "-", "-" ]
   cmd = initial_cmd + exit_cmd
   cmd = "\n".join(cmd)
   if verbose == True: print (cmd)
   fn = root+"_create_constructions.sh"
   f = open(fn, "w")
   f.write("#!/bin/sh")
   f.write("\n# creates zone constructions from attributes.")
   msgcd = "\ncd " + root + '/cfg'
   f.write(msgcd)
   msg = "\nprj -mode text -file " + cfg_file + "<<XXX\n"
   f.write(msg)
   f.write(cmd)
   f.write("\nXXX")
   f.close()
   msgos = "chmod a+x " + fn
   os.system(msgos)

   fdn = root+"_do_all.sh"
   f = open(fdn, "a")
   f.write("\necho "  "")
   msg = "\n./"+fn
   f.write(msg)
   f.write("\necho "  "\n")
   f.close()

   print ('Created zone construction file script.')

# Function to create initial zone operation file (for a small office).
def create_zone_operation_file(root,cfg_file,verbose):

# Command to run within prj to get to menu for zone operations and
# then, set to a small office regime and set 0.5 ach.
   initial_cmd = [ "m", "c", "c", "a", "", "i", "" ]
   vent_doc = "Assumes 0.4 ach all days."
   infil_rate = "0.4"
   set_ach = [ "a", "", "a", vent_doc, "c", "b", "a", infil_rate, ">", "b", "a", infil_rate, ">", "b", "a", infil_rate, "-" ]

# Save opr and then cfg and work back up the menu structures.
   exit_cmd = [ ">", "", "y", "-", "-", "-", "!", "", "-", "-" ]
   cmd = initial_cmd + set_ach + exit_cmd
   cmd = "\n".join(cmd)
   if verbose == True: print (cmd)
   fn = root+"_create_operations.sh"
   f = open(fn, "w")
   f.write("#!/bin/sh")
   f.write("\n# creates zone operations.")
   msgcd = "\ncd " + root + '/cfg'
   f.write(msgcd)
   msg = "\nprj -mode text -file " + cfg_file + "<<XXX\n"
   f.write(msg)
   f.write(cmd)
   f.write("\nXXX")
   f.close()
   msgos = "chmod a+x " + fn
   os.system(msgos)

   fdn = root+"_do_all.sh"
   f = open(fdn, "a")
   f.write("\necho "  "")
   msg = "\n./"+fn
   f.write(msg)
   f.write("\necho "  "\n")
   f.close()
   
   print ('Created zone operations file script.')

# Function to duplicates zone, change right wall to partition, find near surfaces,
# rebuild construction files.
def dup_zone(root,cfg_file,verbose):
   room_name = "right"
   room_doc = "room on the right"
   offset = "4.1 0.0 0.0"

# Commands to run within prj to duplicate the zone.
   initial_cmd = [ "m", "c", "a", "*", "c", "a", offset, "n", room_name, room_doc, "a", "a", "", "a", "a", "f", ]

# Commands to adjust the right wall to be a partition, update the boundary etc.
   second_cmd = [ "b", "e", "n", "y", "b", "c", "f", "m", "-", "y", "-", "-", "!", "", "y", "-", "b", "#", "y", "b", "-", "-", "!", "", "-", "-", ]

   cmd = initial_cmd + second_cmd
   cmd = "\n".join(cmd)
   if verbose == True: print (cmd)
   fn = root+"_duplicate_zone.sh"
   f = open(fn, "w")
   f.write("#!/bin/sh")
   f.write("\n# duplicates zone and connects partition.")
   msgcd = "\ncd " + root + '/cfg'
   f.write(msgcd)
   msg = "\nprj -mode text -file " + cfg_file + "<<XXX\n"
   f.write(msg)
   f.write(cmd)
   f.write("\nXXX")
   f.close()
   msgos = "chmod a+x " + fn
   os.system(msgos)

   fdn = root+"_do_all.sh"
   f = open(fdn, "a")
   f.write("\necho "  "")
   msg = "\n./"+fn
   f.write(msg)
   f.write("\necho "  "\n")
   f.close()
   print ('Room duplicated')


# Function to set and initial february simulation parameter set.
def create_SPS(root,cfg_file,verbose):
   sps_name = "february"
   period_start = "1 2"
   period_end   = "28 2"
   initial_cmd = [ "m", "s", "a", "y", sps_name, "g", period_start, period_end, "-" ]
   exit_cmd = [ "!", "", "-", "-" ]
   cmd = initial_cmd + exit_cmd
   cmd = "\n".join(cmd)
   if verbose == True: print (cmd)
   fn = root+"_create_SPS.sh"
   f = open(fn, "w")
   f.write("#!/bin/sh")
   f.write("\n# creates a simulation parameter set.")
   msgcd = "\ncd " + root + '/cfg'
   f.write(msgcd)
   msg = "\nprj -mode text -file " + cfg_file + "<<XXX\n"
   f.write(msg)
   f.write(cmd)
   f.write("\nXXX")
   f.close()
   msgos = "chmod a+x " + fn
   os.system(msgos)

   fdn = root+"_do_all.sh"
   f = open(fdn, "a")
   f.write("\necho "  "")
   msg = "\n./"+fn
   f.write(msg)
   f.write("\necho "  "\n")
   f.close()
   
   print ('Set initial simulation parameter set script created.')


# Function to generate a markdown model contents report.
def create_QA(root,cfg_file,verbose):
   initial_cmd = [ "m", "u", "y", "a", "c", "h", "h", ">", "", "!", "-" ]
   exit_cmd = [ "!", "", "-", "-" ]
   cmd = initial_cmd + exit_cmd
   cmd = "\n".join(cmd)
   if verbose == True: print (cmd)
   fn = root+"_create_QA.sh"
   f = open(fn, "w")
   f.write("#!/bin/sh")
   f.write("\n# creates QA report.")
   msgcd = "\ncd " + root + '/cfg'
   f.write(msgcd)
   msg = "\nprj -mode text -file " + cfg_file + "<<XXX\n"
   f.write(msg)
   f.write(cmd)
   f.write("\nXXX")
   f.close()
   msgos = "chmod a+x " + fn
   os.system(msgos)

   fdn = root+"_do_all.sh"
   f = open(fdn, "a")
   f.write("\necho "  "")
   msg = "\n./"+fn
   f.write(msg)
   f.write("\necho "  "\n")
   f.close()
   
   print ('Created QA report script.')

# Place any other model manipulation routines here....

# The MAIN which takes in command line arguments. It does the initial
# registration of the ESP-r model and if the user agrees it then
# creates and attributes an initial box shaped zone.
def main(argv):

# These are default values in case the user misses out some command line options.
   root = "box"
   synopsis = "Initial model for testing"
   latitude = "55.86"
   longitude = "-4.25"
   altitude = "50.0"
   year = "2024"
   incremental = False
   verbose = False
   try:
      opts, args = getopt.getopt(argv,"hir:s:l:g:a:y:v",["root=","synopsis=","latitude=","longitude=","alt=","year="])
   except getopt.GetoptError:
      print ('create_test_model_2zone_scripts.py -r --root <root> -s --synopsis <synopsis> -l --latitude <latitude>')
      print (' -g --longitude <longitude-diff> -a --altitude <altitude> -y --year <year>')
      print (' -i forces incremental processing -v verbose mode echoing command sequences')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print ('create_test_model_2zone_scripts.py -r <root> -s <synopsis> -l <latitude> -g <longitude-diff> -a <altitude> -y <year>')
         sys.exit()
      elif opt in ("-i"):
         incremental = True
         print ('interactive mode')
      elif opt in ("-r", "--root"):
         root = arg
      elif opt in ("-s", "--synopsis"):
         synopsis = arg
      elif opt in ("-l", "--latitude"):
         latitude = arg
      elif opt in ("-g", "--longitude"):
         longitude = arg
      elif opt in ("-a", "--altitude"):
         altitude = arg
      elif opt in ("-y", "--year"):
         year = arg
      elif opt in ("-v"):
         verbose = True
         print ('setting verbose mode')

# Assign the cfg file name and echo the attributes either defaulted or supplied.
   cfg_file = root+'.cfg'
   print ('Model root is ', root)
   print ('cfg is ', cfg_file)
   print ('Synopsis is ', synopsis)
   print ('Lat Long are ', latitude,longitude)
   print ('Altitude is ', altitude)
   print ('Year is ', year)
   print ('incremental processing ', incremental)
   print ('verbose ', verbose)

# Create commands to run within prj to register new model and create the folders.
# For Linux we seem to need encoding='utf-8' rather than encoding='ascii'
   cmd = [ "e", root, synopsis, "", "n", "n", latitude, longitude, altitude, year, "b", "e","d", "", "-", "-", "-", ]
   cmd = "\n".join(cmd)   # Separate each of the command tokens with newline characters
   if (verbose == True): print (cmd)
   fn = root+"_initial_steps.sh"
   f = open(fn, "w")
   f.write("#!/bin/sh")
   f.write("\n# creates basic model.")
   msg = "\nprj -mode script " + "<<XXX\n"
   f.write(msg)
   f.write(cmd)
   f.write("\n")
   f.close()
   msgos = "chmod a+x " + fn
   os.system(msgos)
   print ('Initial registration of model '+root+'/cfg/'+cfg_file+' created.')

# Start high level script to drive others.
   fdn = root+"_do_all.sh"
   f = open(fdn, "w")
   f.write("#!/bin/sh")
   f.write("\n# creates basic model via individual scripts.")
   f.write("\necho ")
   msg = "\n./"+fn
   f.write(msg)
   f.write("\necho ")
   f.close()
   msgos = "chmod a+x " + fdn
   os.system(msgos)

# Ask user how to proceed.
   if incremental == True:
      l_reply = ask_YorN('Proceed to add a rectangular zone?')
      if ( l_reply == False ): sys.exit('Exiting.')
   print ()
   print(os.getcwd())

# Do initial shell and attribution of a rectangular room.
   create_box_room(root,cfg_file,verbose)

   if incremental == True:
      l_reply = ask_YorN('Add a window in front wall and a back door?')
      if ( l_reply == False ): sys.exit('Exiting.')
   print ('Processing window and door...')
   add_window_door(root,cfg_file,verbose)
   
   if incremental == True:
      l_reply = ask_YorN('Add an overhang on South facade?')
      if ( l_reply == False ): sys.exit('Exiting.')
   print ('Processing overhang...')
   add_overhang(root,cfg_file,verbose)

   if incremental == True:
      l_reply = ask_YorN('Create zone construction file?')
      if ( l_reply == False ): sys.exit('Exiting.')
   print ('Processing zone constructions...')
   create_zone_construction_file(root,cfg_file,verbose)

   if incremental == True:
     l_reply = ask_YorN('Create zone operations file?')
     if ( l_reply == False ): sys.exit('Exiting.')
   print ('Processing zone operations...')
   create_zone_operation_file(root,cfg_file,verbose)

   if incremental == True:
      l_reply = ask_YorN('Setup an initial simulation parameter set?')
      if ( l_reply == False ): sys.exit('Exiting.')
   print ('Processing simulaton parameter set...')
   create_SPS(root,cfg_file,verbose)

   if incremental == True:
      l_reply = ask_YorN('Duplicate zone and connect?')
      if ( l_reply == False ): sys.exit('Exiting.')
   print ('Processing zone duplication...')
   dup_zone(root,cfg_file,verbose)

   if incremental == True:
      l_reply = ask_YorN('Create model contents report?')
      if ( l_reply == False ): sys.exit('Exiting.')
   print ('Processing model contents report...')
   create_QA(root,cfg_file,verbose)

   if incremental == True:
      l_reply = ask_YorN('Create model from scripts?')
      if ( l_reply == False ): sys.exit('Exiting.')
   print ('Generating model...')
   fdn = "./"+root+"_do_all.sh"
   os.system(fdn)
#    prj5 = run(["./do_all.sh"],text=True,shell=True)
   print ('Generating model... done.')
#    return prj5

   if incremental == True:
      l_reply = ask_YorN('Open new model in prj?')
      if ( l_reply == False ): sys.exit('Exiting.')
   print ('Opening model...')
   os.chdir(root)
   os.chdir("cfg")
   fdn = "prj -file " + cfg_file
   os.system(fdn)

if __name__ == "__main__":
   main(sys.argv[1:])
