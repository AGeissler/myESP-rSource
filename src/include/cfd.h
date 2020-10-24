C This include file should follow building.h in the source code.
C  NTCELX, NTCELY, NTCELZ, and MCEL1D define the resolution of the CFD
C  domain. These variables significantly influence the size of dfs and
C  bps and are set in building.h. The other variables have little
C  impact of the size of dfs and bps.

C  Maximum number of new mfs connections created for CFD domain (usual setting MCFND=20)
      INTEGER, PARAMETER :: MCFND=20
C  Maximum number of heat source (usual setting MNHS=36)
      INTEGER, PARAMETER :: MNHS=36
C  Maximum number of solid boundaries per zone
C (set equal to 2* or greater than MS in building.h, because of current
C method of specification of solid boundary conditions from building
C surfaces)
      INTEGER, PARAMETER :: MNSBZ=MS*2
C  Maximum number of zones with CFD (usual setting MNZ=4)
      INTEGER, PARAMETER :: MNZ=4
C  Maximum number of key volumes (usual setting MNVLS=120)
      INTEGER, PARAMETER :: MNVLS=120
C  Maximum number of contaminants that can be modelled (Set equal to
C  INTEGER, PARAMETER :: MCONTM in net_flow.h, if this is changed also change MCONTM)
      INTEGER, PARAMETER :: MCTM=4
C Maximum frequency of residuals plotting during CFD solutions
      INTEGER, PARAMETER :: MFRP=1000
C Maximum number of cells to potentially refine gridding
      INTEGER, PARAMETER :: MRFN=NTCELX*NTCELY*NTCELZ/1000
C Maximum number of times we might refine in total
      INTEGER, PARAMETER :: MRFNT=NTCELX*NTCELY*NTCELZ


