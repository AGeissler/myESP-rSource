*Geometry 1.1,GEN,mixTop # tag version, format, zone name
*date Wed Jan 15 08:34:35 2014  # latest file modification 
mixTop describes the "mixing zone" at top (outlet)
# tag, X co-ord, Y co-ord, Z co-ord
*vertex,0.00000,-0.75000,4.30000  #   1
*vertex,4.50000,-0.75000,4.30000  #   2
*vertex,4.50000,-0.05000,4.30000  #   3
*vertex,0.00000,-0.05000,4.30000  #   4
*vertex,0.00000,-0.75000,4.45000  #   5
*vertex,4.50000,-0.75000,4.45000  #   6
*vertex,4.50000,-0.05000,4.45000  #   7
*vertex,0.00000,-0.05000,4.45000  #   8
# 
# tag, number of vertices followed by list of associated vert
*edges,4,1,2,6,5  #  1
*edges,4,2,3,7,6  #  2
*edges,4,3,4,8,7  #  3
*edges,4,4,1,5,8  #  4
*edges,4,5,6,7,8  #  5
*edges,4,1,4,3,2  #  6
# 
# surf attributes:
#  surf name, surf position VERT/CEIL/FLOR/SLOP/UNKN
#  child of (surface name), useage (pair of tags) 
#  construction name, optical name
#  boundary condition tag followed by two data items
*surf,Wall-1,VERT,-,-,-,door,OPAQUE,ADIABATIC,0,0  #   1 ||< adiabatic
*surf,Wall-2,VERT,-,-,-,door,OPAQUE,ADIABATIC,0,0  #   2 ||< adiabatic
*surf,Wall-3,VERT,-,-,-,door,OPAQUE,ADIABATIC,0,0  #   3 ||< adiabatic
*surf,Wall-4,VERT,-,-,-,door,OPAQUE,ADIABATIC,0,0  #   4 ||< adiabatic
*surf,Top-5,CEIL,-,-,-,door,OPAQUE,ADIABATIC,0,0  #   5 ||< adiabatic
*surf,Base-6,FLOR,-,-,-,door,OPAQUE,ANOTHER,03,05  #   6 ||< Top-5:TheChannel
# 
*insol,3,0,0,0  # default insolation distribution
# 
# shading directives
*shad_calc,none  # no temporal shading requested
# 
*insol_calc,none  # no insolation requested
# 
*base_list,1,6,     3.15 0  # zone base list
