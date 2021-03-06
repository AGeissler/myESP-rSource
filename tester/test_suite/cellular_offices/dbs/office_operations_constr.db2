# composite construction db defined in ../dbs/office_vent_constr.db2
# based on materials db /usr/esru/esp-r/databases/constr_db3.materialdb
   13     # no of composites 
# layers  description   optics 
    5    extern_wall   OPAQ  OPAQUE              
# db ref  thick   db name & air gap R 
    6    0.1000  Lt brown brick
  211    0.1500  Glasswool
    0    0.0500  air  0.170 0.170 0.170
   22    0.1000  Aerated conc block
  108    0.0120  White ptd Gypboard
# layers  description   optics 
    3    insul_frame   OPAQ  OPAQUE              
# db ref  thick   db name & air gap R 
   46    0.0040  Grey cotd aluminium
  281    0.0800  Glass Fibre Quilt
   46    0.0040  Grey cotd aluminium
# layers  description   optics 
    1    door          OPAQ  OPAQUE              
# db ref  thick   db name & air gap R 
   69    0.0250  Oak (radial)
# layers  description   optics 
    1    mass_part     OPAQ  OPAQUE              
# db ref  thick   db name & air gap R 
   35    0.2400  Block white ptd inner (3% mc)
# layers  description   optics 
    1    fict          TRAN  SC_fictit           
# db ref  thick   db name & air gap R 
  245    0.0060  fict
# layers  description   optics 
    3    dbl_glz       TRAN  DCF7671_06nb        
# db ref  thick   db name & air gap R 
  242    0.0060  Plate glass
    0    0.0120  air  0.170 0.170 0.170
  242    0.0060  Plate glass
# layers  description   optics 
    4    roof          OPAQ  OPAQUE              
# db ref  thick   db name & air gap R 
   43    0.0030  Aluminium
    0    0.0250  air  0.170 0.170 0.170
  281    0.0800  Glass Fibre Quilt
   43    0.0030  Aluminium
# layers  description   optics 
    5    susp_floor    OPAQ  OPAQUE              
# db ref  thick   db name & air gap R 
  221    0.0060  Wilton
   67    0.0190  Chipboard
    0    0.0500  air  0.170 0.170 0.170
   32    0.1400  Heavy mix concrete
   42    0.0040  Steel
# layers  description   optics 
    5    susp_flr_re   OPAQ  OPAQUE              
# db ref  thick   db name & air gap R 
   42    0.0040  Steel
   32    0.1400  Heavy mix concrete
    0    0.0500  air  0.170 0.170 0.170
   67    0.0190  Chipboard
  221    0.0060  Wilton
# layers  description   optics 
    2    ceiling       OPAQ  OPAQUE              
# db ref  thick   db name & air gap R 
  211    0.1000  Glasswool
  150    0.0100  Ceiling (mineral)
# layers  description   optics 
    2    ceiling_rev   OPAQ  OPAQUE              
# db ref  thick   db name & air gap R 
  150    0.0100  Ceiling (mineral)
  211    0.1000  Glasswool
# layers  description   optics 
    5    gyp_blk_ptn   OPAQ  OPAQUE              
# db ref  thick   db name & air gap R 
  108    0.0130  White ptd Gypboard
    0    0.0500  air  0.170 0.170 0.170
   28    0.1000  Block inner (3% mc)
    0    0.0500  air  0.170 0.170 0.170
  108    0.0130  White ptd Gypboard
# layers  description   optics 
    3    gyp_gyp_ptn   OPAQ  OPAQUE              
# db ref  thick   db name & air gap R 
  108    0.0120  White ptd Gypboard
    0    0.0500  air  0.170 0.170 0.170
  108    0.0120  White ptd Gypboard
