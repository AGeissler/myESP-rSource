C This file is part of the ESP-r system.
C Copyright Energy Systems Research Unit, University of
C Strathclyde, Glasgow Scotland, 2001.

C ESP-r is free software.  You can redistribute it and/or
C modify it under the terms of the GNU General Public
C License as published by the Free Software Foundation 
C (version 2 or later).

C ESP-r is distributed in the hope that it will be useful
C but WITHOUT ANY WARRANTY; without even the implied
C warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
C PURPOSE. See the GNU General Public License for more
C details.

C You should have received a copy of the GNU General Public
C License along with ESP-r. If not, write to the Free
C Software Foundation, Inc., 59 Temple Place, Suite 330,
C Boston, MA 02111-1307 USA.

C This source file contains:
C CFDCOMP - controls the editing of CFD input files.
C CFDGGRD - define domain gridding.
C CFDSVAR - controls the editing of CFD solution variables.
C CFDEDAT - controls the editing of CFD equations.
C CFDEDATC- controls the editing of contaminant relaxation factors.
C CFDBVAR - controls the editing of CFD boundary variables.
C EDAIRO  - controls the editing of CFD air flow opening boundary variables.
C MFCONF  - defines conflation of CFD and mfs 
C EDSLDB  - controls the editing of CFD solid boundary variables.
C EDSRC   - controls the editing of CFD sources - humidity, CO2, etc.
C EDBLK   - controls the editing of CFD blockages to flow
C HANDSHK - controls the editing of the handshaking mechanism.
C PIKCELS - controls the editing of the cells defining a boundary region.
C PIKFACE - controls the editing of the faceof the boundary region.
C TOPBOT  - define upper and lower curvilinear x-section for use in
C           gridding generation.
C GRIDDEF - define gridding via regions and specify number of cells
C           and cell distribution for each region.
C ASKMFNOD- asks for one or more mass flow nodes (copy of mfoutp.F).
C ASKMFCON- Asks for one or more mass flow connections (copy of mfoutp.F).  

C << TODO
C << The introduction of cancel options for easki calls should be
C << extended to ensure that data in a sequence of edits is only
C << instanciated if no cancels are called in the sequence. Currently
C << the code does not re-establish some of the earlier values in
C << a mixed integer and real sequence of edits.

C ************* CFDCOMP 
C CFDCOMP controls the editing of CFD input files.
C ITRU unit number for user output, IUF unit number for CFD input file.
C IER=0 indicates no error.

      SUBROUTINE CFDCOMP(izone,iuf,IER)
#include "building.h"
#include "model.h"
#include "geometry.h"
#include "cfd.h"
#include "espriou.h"
C #include "esprdbfile.h"
#include "prj3dv.h"
      
      integer lnblnk  ! function definition

      COMMON/SPAD/MMOD,LIMIT,LIMTTY
      common/OUTIN/IUOUT,IUIN
      common/gzonpik/izgfoc,nzg,nznog(mcom)
      COMMON/DEFLT3/DFCFD,DECMPDBFL,DICONDBFL
      common/cfdfil/LCFD(MCOM),IFCFD(MCOM)
      COMMON/ICFNOD/ICFD,ICP
      common/param2/TITLE(MNZ),CFTRFL(MNZ),LPHI(MNZ)
      common/GRIDEFN/origin,xgrid,ygrid,zgrid,zegrid,ortho
      common/ndcfd/ncfdnd,icfdnd(MNZ),NCONF
      integer ncfdnd,icfdnd,NCONF
      COMMON/CFDVIS/HAS_GEOM,ISHSB,ISHAO,IFACES,ISHBLK,ISHGEO,INITD
      logical HAS_GEOM
      character INITD*6

C 3D visualisation mode.
      COMMON/MODVIS/IVISMOD

      logical origin,xgrid,ygrid,zgrid,zegrid,ortho
      logical OK,dok,XST,CFDsave

      character OUTSTR*124,WORD*124,ETEXT*82
      CHARACTER ITEM(12)*34
      character ltmp*72,dtmp*72,LCFD*72
      CHARACTER*72 TITLE,CFTRFL,LPHI
      character*72 DFCFD,DECMPDBFL,DICONDBFL
      integer found,zngeo
      integer IWO  ! for radio button
      integer NITMS,INO ! max items and current menu item

#ifdef OSI
      integer iside,isize,ifont
#else
      integer*8 iside,isize,ifont
#endif

C For help messages
      character helpinsub*24 ! subroutine name
      character helptopic*24 ! string (unique) for topic
      integer nbhelp     ! number of help lines found

      helpinsub='CFDCOMP'  ! set for subroutine

      IVISMOD=-1

C ICFD=1 for the first CFD domain in the building, 2 for the second CFD
C domain, etc. ICFD is set independantly in VERMAN in prj.F, any updates
C here should be mirrored there as well.
      found=0
      do 234 iconf=1,nconf
        if ( icfdnd(iconf).eq.izone )then
          ICFD=iconf; found=1
        endif
  234 continue
      if(found.eq.0)then

C No existing domain for this zone, clear domain indepenent commons,
C set ICFD, and initialise simulation parameters.
C Do not increment NCONF yet, as user might not include the domain.
        CALL RESETCFD(IER)
        ICFD=nconf+1
        CALL CFDDEFLT
      endif
      ICP = izone
      TITLE(ICFD) = '...'

C Set the flag to indicate that the data has not been saved to the CFD input file.
      CFDsave=.false.

C Set flags to indiate the origin of CFD domain and gridding of axes
C have not yet been defined.
      origin=.false.; xgrid=.false.; ygrid=.false.; zgrid=.false.

C= Assume orthogonal domain.
      ortho=.true.

C Draw the zone and display the vertex numbers.
C << changed to use CFD visualisation subroutines instead >>
      call tstamp('>','PRJ: edit CFD domain') 
      nzg=1; nznog(1)=izone; izgfoc=izone
      itvno=0; itsnm=0

C Setup help text for main menu.
      helptopic='cfd_opening_menu'
      call gethelptext(helpinsub,helptopic,nbhelp)

C Get the name of the CFD input file.
      if(LCFD(izone)(1:2).eq.'  ')then
        if(netpth(1:2).eq.'  '.or.netpth(1:2).eq.'./')then
          WRITE(LCFD(izone),'(A,A4)')
     &      zname(izone)(1:lnzname(izone)),'.dfd'
        else
          WRITE(LCFD(izone),'(A,A,A,A4)') netpth(1:lnblnk(netpth)),'/',
     &      zname(izone)(1:lnzname(izone)),'.dfd'
        endif
      endif
      ltmp=LCFD(izone)
      CALL EASKS(ltmp,' CFD input file?','  ',72,DFCFD,'dfd file',
     &  IER,nbhelp)
      if(ltmp(1:2).ne.'  ')LCFD(izone)=ltmp

C If the CFD input file exists, read in its contents.
C This is a two-stage process: first open the file and get the conflation
C flag.  Then read in the data based on the conflation type.
      CALL ERPFREE(IUF,ISTAT)
      call FINDFIL(LCFD(izone),XST)
      if(XST)then
C set currentfile
        write(currentfile,'(a)')LCFD(izone)(1:LNBLNK(LCFD(izone)))
        CALL EFOPSEQ(IUF,LCFD(izone),1,IER)
        IF(IER.NE.0)goto 1
        CALL STRIPC(IUF,OUTSTR,0,ND,1,'dfd line 1',IER)
        if(OUTSTR(1:7).eq.'*DFS V2')then
          CALL ERPFREE(IUF,ISTAT)
          call DFDREAD(IZONE,0,iuout,IER)
          if (IER.ne.0) return
        elseif(OUTSTR(1:15).eq.'DFS DESCRIPTION')then
          CALL STRIPC(IUF,OUTSTR,0,ND,1,'dfd line 2',IER)
          K=0
          CALL EGETW(OUTSTR,K,WORD,'W','tag',IFLAG)
          if(WORD(1:11).eq.'*conflation')then

C Set the conflation flag.
            CALL EGETWI(OUTSTR,K,iv,0,3,'F','confla type',IER)
            IFCFD(izone)=iv
            CALL ERPFREE(IUF,ISTAT)

C Set up for reading in CFD input file.
            IBLD=0; IMFN=0
            if(IFCFD(izone).eq.1) IBLD=1
            if(IFCFD(izone).eq.2) IBLD=1
            if(IFCFD(izone).eq.3) IMFN=1
            CALL EFOPSEQ(IUF,LCFD(izone),1,IER)

C << Force reporting. >>
            CALL CFDDTA(2,iuout,izone,IBLD,IMFN,IER)

C Set flags to indicate that CFD origin and gridding has been defined.
            origin=.true.; xgrid=.true.; ygrid=.true.; zgrid=.true.
          else
            call usrmsg('DFD confl. tag missing in:',LCFD(izone),'W')
            CALL ERPFREE(IUF,ISTAT)
            return
          endif
        else
          call usrmsg('Not a recognised file in:',LCFD(izone),'W')
          CALL ERPFREE(IUF,ISTAT)
          return
        endif

C Default values for visualisation.
        if (origin.and.xgrid.and.ygrid.and.zgrid) then
          IVISMOD=2
          HAS_GEOM=.TRUE.
          ISHSB=-1; ISHAO=-1; IFACES=1; ISHBLK=-1; ISHGEO=1
        endif
      endif
      CALL ERPFREE(IUF,ISTAT)

      MODIFYVIEW=.TRUE.
      MODLEN=.TRUE.
      MODBND=.TRUE.

C Present menu of options for editing CFD input file.
   10 INO=-4
      ITEM(1) ='a Title: '
      write (ITEM(2),'(a,a)') '  ',TITLE(ICFD)(1:30)
      ITEM(3) ='  -------------------------- '
      if (IFCFD(izone).eq.0) then
        ITEM(4) ='a CFD coupling >> Off '
      else
        ITEM(4) ='a CFD coupling >> On  '
      endif
      ITEM(5) ='b Geometry and gridding '
      ITEM(6) ='c Solution variables '
      ITEM(7) ='d Boundary conditions '
      ITEM(8) ='  -------------------------- '
      ITEM(9) ='f Report domain details '
      ITEM(10) ='> Save CFD input file '
      ITEM(11)='? Help'
      ITEM(12)='- exit this menu'
      NITMS=12

      IF(MMOD.EQ.8)THEN

C Only display if domain geometry has been defined.
        if (origin.and.xgrid.and.ygrid.and.zgrid) then


          call redraw(IER)

        else

C Just display title.
          WRITE(ETEXT,'(2A)')'Domain: ',
     &      TITLE(izone)(1:lnblnk(TITLE(izone)))
          iside=1; isize=1; ifont=1
          call viewtext(ETEXT,iside,isize,ifont)
        endif
      
      ENDIF

C Help text for this menu.
  12  helptopic='cfd_opening_menu'
      call gethelptext(helpinsub,helptopic,nbhelp)
 
      CALL EMENU('Zone air movement description',ITEM,NITMS,INO)

      IF(INO.EQ.0)THEN

C Wrong pick.
        INO=-1
        GOTO 12
      ELSEIF(INO.EQ.1)THEN
        dtmp='Model of room air flow'
        ltmp=TITLE(ICFD)
        call EASKS(ltmp,'Description of CFD model:',' ',72,dtmp,
     &             'CFD description',IER,42)
        TITLE(ICFD)=ltmp
      ELSEIF(INO.EQ.4)THEN

C Conflation options.
        helptopic='cfd_conflation_options'
        call gethelptext(helpinsub,helptopic,nbhelp)
        dok=.true.
        call ASKOK(
     &      'Do you want to couple CFD with other simulation domains?',
     &                                               ' ',OK,dok,nbhelp)
        if(OK)then
          IFCFD(izone)=4
          helptopic='cfd_conflation_active'
          call gethelptext(helpinsub,helptopic,nbhelp)
          call PHELPD('CFD coupling activated',nbhelp,'-',0,0,IER)
        else 
          IFCFD(izone)=0
        endif
      ELSEIF(INO.EQ.5)THEN 

C Geometry and Gridding
        call CFDGGRD(izone,IER)

C Recalculate gridding in case it has changed.
        CALL GRID(ier)
        CALL NEW2OLD
      ELSEIF(INO.EQ.6)THEN 

C Solution Variables and Initial Conditions
C User to specify equations to solve and initial variables.
        call CFDSVAR(IER)
      ELSEIF(INO.EQ.7)THEN 

C Boundary conditions.
        call CFDBVAR(IZONE,IER)
      ELSEIF(INO.EQ.9)THEN

C Report model details. Do this by re-reading input file with trace output.
        CALL ERPFREE(IUF,ISTAT)
        call DFDREAD(IZONE,2,iuout,IER)
      ELSEIF(INO.EQ.10)THEN 

C Save Data to dfs input file.
        call DFDSV(IUF,IZONE,IER)
        if (IER.eq.0) CFDsave=.true.

C Help.
      ELSEIF(INO.EQ.(NITMS-1))THEN
        helptopic='cfd_opening_menu'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL PHELPD('CFD inputs',nbhelp,'-',0,0,IER)
      ELSEIF(INO.EQ.NITMS)THEN

C Exit CFD editing menu.
        if(.not.CFDsave)then
          dok=.false.
          CALL ASKOK('Data not yet saved to CFD input file.',
     &               'Are you sure you want to exit? ',OK,dok,nbhelp)
          IF(.NOT.OK)GOTO 10
        endif
        return
      ENDIF
      GOTO 10

C Error handling
    1 CALL USRMSG('Problem with CFD input file line:',OUTSTR,'W')
      call edisp(iuout,' Plese check your CFD input file.')
      return

      end

C ******************* CFDGGRD *******************
C CFDGGRD - define domain gridding.
C IUF unit number for CFD input file.
C IER=0 indicates no error.

      SUBROUTINE CFDGGRD(izone,IER)
#include "building.h"
#include "model.h"
#include "geometry.h"
#include "cfd.h"
#include "prj3dv.h"

      COMMON/SPAD/MMOD,LIMIT,LIMTTY
      COMMON/FILEP/IFIL
      common/OUTIN/IUOUT,IUIN
      integer menuchw,igl,igr,igt,igb,igw,igwh
      common/VIEWPX/menuchw,igl,igr,igt,igb,igw,igwh
      common/gzonpik/izgfoc,nzg,nznog(mcom)

C NREG(1,*) is X, NREG(2,*) is Y, NREG(3,*) is west, NREG(4,*) is east
C (the east (4,*) is for curvilinear gridding).
      common/GRIDFN/NCELX(MNREG,MNZ),NCELY(MNREG,MNZ),NCELZ(MNREG,MNZ),
     &  NCELZE(MNREG,MNZ),XREG(MNREG,MNZ),YREG(MNREG,MNZ),
     &  ZREG(MNREG,MNZ),ZREGE(MNREG,MNZ),Xplaw(MNREG,MNZ),
     &  Yplaw(MNREG,MNZ),Zplaw(MNREG,MNZ),Zplawe(MNREG,MNZ),NREG(4,MNZ)
      common/GRDVRTS/iorg(MNZ),ixend(MNZ),iyend(MNZ),izend(MNZ),
     &  izende(MNZ)

      common/GRIDEFN/origin,xgrid,ygrid,zgrid,zegrid,ortho
      COMMON/ICFNOD/ICFD,ICP
      common/grdmax/NTCX,NTCY,NTCZ

      common/EQTION/CALLU(MNZ),CALLV(MNZ),CALLW(MNZ),CALLT(MNZ),
     &             CALLC(MNZ),KEMDL(MNZ),BUOY(MNZ),BOUSSI(MNZ),
     &             ZEROT(MNZ),ZandKE(MNZ),MITzero(MNZ)

      COMMON/CFDVIS/HAS_GEOM,ISHSB,ISHAO,IFACES,ISHBLK,ISHGEO,INITD
      logical HAS_GEOM
      character INITD*6

C 3D visualisation mode.
      COMMON/MODVIS/IVISMOD
      common/param2/TITLE(MNZ),CFTRFL(MNZ),LPHI(MNZ)

C Arrays for sorting domains from (see code for comments).
      dimension XS(MTV),XSS(MTV),XSI(MTV),YS(MTV),YSS(MTV),YSI(MTV)
      dimension ZS(MTV),ZSI(MTV),ZSS(MTV)

      LOGICAL CALLU,CALLV,CALLW,CALLT,CALLC,KEMDL,BUOY,BOUSSI
      logical close,griderr
      logical origin,xgrid,ygrid,zgrid,zegrid,ortho
      LOGICAL ZEROT,ZandKE,MITzero,vclose

      character ITEM(14)*36
      character outs*124,outs248*248,ETEXT*82
      CHARACTER*72 TITLE,CFTRFL,LPHI
      integer zngeo

      integer iglib   ! if 1 then X11, if 2 then GTK, if 3 then text only.
      integer iorgt,ixendt,iyendt,izendt,izendet  ! for local editing
      integer NITMS,INO ! max items and current menu item

#ifdef OSI
      integer iigl,iigr,iigt,iigb,iigw,iigwh
      integer iiw1,iiw2,iiw3,iiw4,iimenu
      integer iside,isize,ifont
      integer ilf,igfw,igfh,ild,igdw,igdh
#else
      integer*8 iigl,iigr,iigt,iigb,iigw,iigwh
      integer*8 iiw1,iiw2,iiw3,iiw4,iimenu
      integer*8 iside,isize,ifont
      integer*8 ilf,igfw,igfh,ild,igdw,igdh
#endif

C For help messages
      character helpinsub*24 ! subroutine name
      character helptopic*24 ! string (unique) for topic
      integer nbhelp     ! number of help lines found

      helpinsub='CFDGGRD'  ! set for subroutine

C Assume axes defined and set initial values for general resolution
C of the X Y Z axis and visualisation options.
C << Todo: set these based on a higher level concept of
C << low/medium/high resolution and/or the size of the zone.
C -- Visualisation variables already defined in calling routine
      origin=.TRUE.
      ieorg=0; icellx=20; icelly=20; icellz=20
      itrc=0   ! silent feedback initially
c      ISHSB=1
c      ISHAO=1
c      ISHBLK=1
c      IFACES=1

C Help text for this menu.
      helptopic='cfd_gridding_setup'
      call gethelptext(helpinsub,helptopic,nbhelp)

C Ensure that G1 common blocks are correctly assigned.
      call georead(IFIL+1,LGEOM(izone),izone,1,IUOUT,IER)

C Check if we have vertex at ends of axis data.
      if (IORG(ICFD).eq.0) then

        call usrmsg('You will have to define the CFD axis relationship',
     &              'to the zone geometry via vertex numbers.','W')

C Display zone geometry with vertex numbers to aid in selection
C Restore standard graphics window variables.
        if(MMOD.EQ.8)then

C Setup and pass in parameters to win3d.
          iiw1=4; iiw2=1; iiw3=1; iiw4=3; iimenu=menuchw
          iigl=igl; iigr=igr; iigt=igt; iigb=igb; iigw=igw; iigwh=igwh
          ilf=2; ild=LIMTTY
          iglib = igraphiclib()  ! find out if X11 or GTK or text support only.
          if(iglib.eq.1)then
            CALL feedbox(iimenu,ilf,igfw,igfh)
            CALL opengdisp(iimenu,ild,ilf,igdw,igdh)
          endif
          CALL win3d(iimenu,iiw1,iiw2,iiw3,iiw4,
     &      iigl,iigr,iigt,iigb,iigw,iigwh)
          igl=int(iigl); igr=int(iigr); igt=int(iigt); igb=int(iigb)
          igw=int(iigw); igwh=int(iigwh)
          call startbuffer()

C Redraw zone and vertex points. Note, if user has returned
C from 3D grid visualisation then G1 common blocks will need
C to be refreshed.
          nzg=1; nznog(1)=izone; izgfoc=izone
          CALL ESCZONE(izone)
          CALL BNDOBJ(0,IER)
          CALL ERCZONE(izone)
          MODIFYVIEW=.TRUE.; MODLEN=.TRUE.; MODBND=.TRUE.
          itvno=0; itznm=1; itsnm=1
          CALL INLNST(1)
          CALL ADJVIEW(IER)
        endif

        write(outs,*) 'Coordinates in the zone...'
        call edisp(iuout,outs)
        write(outs,*) '         1      2      3      4      5      6'
        call edisp(iuout,outs)
        indisp=MIN0(32,NZTV(izone))
        write(outs248,'(a,32F7.3)') 'x axis',(X(j),j=1,indisp)
        call edisp248(iuout,outs248,100)
        write(outs248,'(a,32F7.3)') 'y axis',(Y(j),j=1,indisp)
        call edisp248(iuout,outs248,100)
        write(outs248,'(a,32F7.3)') 'z axis',(Y(j),j=1,indisp)
        call edisp248(iuout,outs248,100)

C Ask user for zone vertex indices of the origin of cfd domain
C followed by vertex at end of X axis, vertex at end of Y axis
C and the vertex for the upper Z axis point (usually above the
C Y axis end point).
        iorgt=IORG(ICFD)
        CALL EASKI(iorgt,' Specify vertex at origin of CFD domain. ',
     &    ' Which `v`? ',1,'F',NTV,'F',1,'orig vertex',IERI,nbhelp)
        if(ieri.eq.-3) then
          origin=.false.
          return
        endif
        ixendt=ixend(ICFD)
        CALL EASKI(ixendt,' Specify vertex at end of X axis. ',
     &    ' Which `v`? ',1,'F',NTV,'F',1,'x-axis vertex',IERI,nbhelp)
        if(ieri.eq.-3) then
          origin=.false.
          return
        endif

        iyendt=iyend(ICFD)
        CALL EASKI(iyendt,' Specify vertex at end of Y axis. ',
     &    ' Which `v`? ',1,'F',NTV,'F',1,'y-axis vertex',IERI,nbhelp)
        if(ieri.eq.-3) then
          origin=.false.
          return
        endif

        izendt=izend(ICFD)
        CALL EASKI(izendt,' Specify vertex at end of Z axis. ',
     &    ' Which `v`? ',1,'F',NTV,'F',1,'z-axis vertex',IERI,nbhelp)
        if(ieri.eq.-3) then
          origin=.false.
          return
        endif

C No cancels so instanciate variables.
        IORG(ICFD)=iorgt
        ixend(ICFD)=ixendt; iyend(ICFD)=iyendt; izend(ICFD)=izendt
        IVISMOD=2
        HAS_GEOM=.TRUE.
        ISHSB=-1; ISHAO=-1; IFACES=1; ISHBLK=-1
      endif

C Make sure axes are orthogonal, as dfs only supports Cartesian coordinates.
C First ensure that G1 common blocks are correctly assigned.
      call georead(IFIL+1,LGEOM(izone),izone,1,IUOUT,IER)

C Check X and Y.
      call ang3vtx(X(ixend(ICFD)),Y(ixend(ICFD)),Z(ixend(ICFD)),
     &  X(IORG(ICFD)),Y(IORG(ICFD)),Z(IORG(ICFD)),X(iyend(ICFD)),
     &  Y(iyend(ICFD)),Z(iyend(ICFD)),ang3)
      call eclose(ang3,90.00,0.1,close)
      if(.NOT.close) then
        call usrmsg(
     &    'X & Y axes are not orthogonal. Please check and respecify.',
     &    '  ','W')
        origin=.FALSE.; xgrid=.FALSE.
      endif

C Check X and Z.
      call ang3vtx(X(ixend(ICFD)),Y(ixend(ICFD)),Z(ixend(ICFD)),
     &  X(IORG(ICFD)),Y(IORG(ICFD)),Z(IORG(ICFD)),X(izend(ICFD)),
     &  Y(izend(ICFD)),Z(izend(ICFD)),ang3)
      call eclose(ang3,90.00,0.1,close)
      if(.NOT.close) then
        call usrmsg(
     &    'X & Z axes are not orthogonal. Please check and respecify.',
     &    '  ','W')
        origin=.FALSE.; ygrid=.FALSE.
      endif

C Check Y and Z.
      call ang3vtx(X(iyend(ICFD)),Y(iyend(ICFD)),Z(iyend(ICFD)),
     &  X(IORG(ICFD)),Y(IORG(ICFD)),Z(IORG(ICFD)),X(izend(ICFD)),
     &  Y(izend(ICFD)),Z(izend(ICFD)),ang3)
      call eclose(ang3,90.00,0.1,close)
      if(.NOT.close) then
        call usrmsg(
     &    'Y & Z axes are not orthogonal. Please check and respecify.',
     &    '  ','W')
        origin=.FALSE.; zgrid=.FALSE.
      endif

C Determine length of each axis.
      xdis=crowxyz(X(IORG(ICFD)),Y(IORG(ICFD)),Z(IORG(ICFD)),
     &  X(ixend(ICFD)),Y(ixend(ICFD)),Z(ixend(ICFD)))
      ydis=crowxyz(X(IORG(ICFD)),Y(IORG(ICFD)),Z(IORG(ICFD)),
     &  X(iyend(ICFD)),Y(iyend(ICFD)),Z(iyend(ICFD)))
      zdis=crowxyz(X(IORG(ICFD)),Y(IORG(ICFD)),Z(IORG(ICFD)),
     &  X(izend(ICFD)),Y(izend(ICFD)),Z(izend(ICFD)))

C Non orthogonal.
      if (.NOT.ortho) then

C Check Y and Ze.
        call ang3vtx(X(iyend(ICFD)),Y(iyend(ICFD)),Z(iyend(ICFD)),
     &    X(IORG(ICFD)),Y(IORG(ICFD)),Z(IORG(ICFD)),X(izende(ICFD)),
     &    Y(izende(ICFD)),Z(izende(ICFD)),ang3)
        call eclose(ang3,90.00,0.1,close)
        if(.NOT.close) then
          call usrmsg('Y & Ze axes are not orthogonal. ',
     &                'Please check and respecify.','W')
          INO=1
          goto 9
        endif
      endif

C Set-up data for menu.
   14 INO=-4

C Check for gridding errors.
      if (origin.AND.xgrid.AND.ygrid.AND.zgrid) then
        griderr=.FALSE.
      else
        griderr=.TRUE.
      endif

C Need to define counter for variable length menu.
      NCP=0

C Create menu text.
      ITEM(1) ='a Define origin and axes.'
      ITEM(2) ='  Vertex ids:'
      if (ortho) then
        write (ITEM(3),'(4(a,i2))') ' O=',IORG(ICFD),': Vx=',
     &    ixend(ICFD),', Vy=',iyend(ICFD),', Vz=',izend(ICFD)
      else
        write (ITEM(3),'(5(a,i2))') ' O=',IORG(ICFD),': Vx=',
     &    ixend(ICFD),', Vy=',iyend(ICFD),', Vz=',izend(ICFD),' Vze=',
     &    izende(ICFD)
      endif
      ITEM(4) =' ------------------------------ '
      ITEM(5) ='b Estimate regions from geometry'
      ITEM(6) =' Axis Regions Total cells'
      write (ITEM(7),'(a,2i6)') 'c X  ',NREG(1,ICFD),NTCX
      write (ITEM(8),'(a,2i6)') 'd Y  ',NREG(2,ICFD),NTCY
      write (ITEM(9),'(a,2i6)') 'e Z  ',NREG(3,ICFD),NTCZ
      if (.NOT.ortho) then
        write (ITEM(10),'(a,2i6)') 'f Ze ',NREG(4,ICFD),NTCZe
        NCP=1
      endif
      if (griderr) then
        ITEM(10+NCP) ='  Gridding incomplete or errors! '
      else
        ITEM(10+NCP) ='g Visualize gridding             '
      endif
      ITEM(11+NCP) =' ------------------------------ '
      if(itrc.eq.0)then
        ITEM(12+NCP) ='> feedback silent             '
      elseif(itrc.eq.1)then
        ITEM(12+NCP) ='> feedback brief              '
      elseif(itrc.eq.2)then
        ITEM(12+NCP) ='> feedback verbose            '
      endif
      ITEM(13+NCP)   ='? help                        '
      ITEM(14+NCP)   ='- exit                        '
      NITMS=14+NCP

      IF(MMOD.EQ.8)THEN

C Recalculate gridding in case this has changed.
        CALL INICNT
        CALL GRID(ier)
        call NEW2OLD

        call redraw(IER)
      ENDIF

C Help text for this menu.
      helptopic='cfd_gridding_setup'
      call gethelptext(helpinsub,helptopic,nbhelp)

C Display menu.
      CALL EMENU('Geometry and gridding',ITEM,NITMS,INO)
 9    if(INO.EQ.1)then

C Display zone geometry with vertex numbers to aid in selection
C Restore standard graphics window variables.
        if(MMOD.EQ.8)then

C Setup and pass in parameters to win3d.
          iiw1=4; iiw2=1; iiw3=1; iiw4=3; iimenu=menuchw
          iigl=igl; iigr=igr; iigt=igt; iigb=igb; iigw=igw; iigwh=igwh
          ilf=2; ild=LIMTTY
          iglib = igraphiclib()  ! find out if X11 or GTK or text support only.
          if(iglib.eq.1)then
            CALL feedbox(iimenu,ilf,igfw,igfh)
            CALL opengdisp(iimenu,ild,ilf,igdw,igdh)
          endif
          CALL win3d(iimenu,iiw1,iiw2,iiw3,iiw4,
     &      iigl,iigr,iigt,iigb,iigw,iigwh)
          igl=int(iigl); igr=int(iigr); igt=int(iigt); igb=int(iigb)
          igw=int(iigw); igwh=int(iigwh)
          call startbuffer()

C Redraw zone and vertex points. Note, if user has returned
C from 3D grid visualisation then G1 common blocks will need
C to be refreshed.
          nzg=1; nznog(1)=izone; izgfoc=izone
          call georead(IFIL+1,LGEOM(izone),izone,1,IUOUT,IER)
          CALL ESCZONE(izone)
          CALL BNDOBJ(0,IER)
          CALL ERCZONE(izone)
          MODIFYVIEW=.TRUE.; MODLEN=.TRUE.; MODBND=.TRUE.
          itvno=0; itsnm=0
          CALL INLNST(1)
          CALL ADJVIEW(IER)
        endif

C Specify origin and axes of CFD domain.
        write(outs,*) 'Coordinates in the zone...'
        call edisp(iuout,outs)
        write(outs,*) '         1      2      3      4      5      6'
        call edisp(iuout,outs)
        indisp=MIN0(32,NZTV(izone))
        write(outs248,'(a,32F7.3)') 'x axis',(X(j),j=1,indisp)
        call edisp248(iuout,outs248,100)
        write(outs248,'(a,32F7.3)') 'y axis',(Y(j),j=1,indisp)
        call edisp248(iuout,outs248,100)
        write(outs248,'(a,32F7.3)') 'z axis',(Y(j),j=1,indisp)
        call edisp248(iuout,outs248,100)
        iorgt=IORG(ICFD)
        CALL EASKI(iorgt,' Specify vertex at origin of CFD domain. ',
     &    ' Which `v`? ',1,'F',NTV,'F',1,'orig vertex',IERI,nbhelp)
        if(ieri.eq.-3) then
          origin=.false.
          return
        endif

        ixendt=ixend(ICFD)
        CALL EASKI(ixendt,' Specify vertex at end of X axis. ',
     &    ' Which `v`? ',1,'F',NTV,'F',1,'x-axis vertex',IERI,nbhelp)
        if(ieri.eq.-3) then
          origin=.false.
          return
        endif

        iyendt=iyend(ICFD)
        CALL EASKI(iyendt,' Specify vertex at end of Y axis. ',
     &    ' Which `v`? ',1,'F',NTV,'F',1,'y-axis vertex',IERI,nbhelp)
        if(ieri.eq.-3) then
          origin=.false.
          return
        endif

        izendt=izend(ICFD)
        CALL EASKI(izendt,' Specify vertex at end of Z axis. ',
     &    ' Which `v`? ',1,'F',NTV,'F',1,'z-axis vertex',IERI,nbhelp)
        if(ieri.eq.-3) then
          origin=.false.
          return
        endif

C No cancels so instanciate variables.
        IORG(ICFD)=iorgt
        ixend(ICFD)=ixendt; iyend(ICFD)=iyendt; izend(ICFD)=izendt

C Make sure axes are orthogonal, as dfs only supports Cartesian coordinates.
C Check X and Y.
        call ang3vtx(X(ixend(ICFD)),Y(ixend(ICFD)),Z(ixend(ICFD)),
     &    X(IORG(ICFD)),Y(IORG(ICFD)),Z(IORG(ICFD)),X(iyend(ICFD)),
     &    Y(iyend(ICFD)),Z(iyend(ICFD)),ang3)
        call eclose(ang3,90.00,0.1,close)
        if(.NOT.close) then
          call edisp(iuout,' X & Y axes are not orthogonal. Respecify.')
          goto 9
        endif

C Check X and Z.
        call ang3vtx(X(ixend(ICFD)),Y(ixend(ICFD)),Z(ixend(ICFD)),
     &    X(IORG(ICFD)),Y(IORG(ICFD)),Z(IORG(ICFD)),X(izend(ICFD)),
     &    Y(izend(ICFD)),Z(izend(ICFD)),ang3)
        call eclose(ang3,90.00,0.1,close)
        if(.NOT.close) then
          call edisp(iuout,' X & Z axes are not orthogonal. Respecify.')
          goto 9
        endif

C Check Y and Z.
        call ang3vtx(X(iyend(ICFD)),Y(iyend(ICFD)),Z(iyend(ICFD)),
     &    X(IORG(ICFD)),Y(IORG(ICFD)),Z(IORG(ICFD)),X(izend(ICFD)),
     &    Y(izend(ICFD)),Z(izend(ICFD)),ang3)
        call eclose(ang3,90.00,0.1,close)
        if(.NOT.close) then
          call edisp(iuout,' Y & Z axes are not orthogonal. Respecify.')
          goto 9
        endif

C Determine length of each axis.
        xdis=crowxyz(X(IORG(ICFD)),Y(IORG(ICFD)),Z(IORG(ICFD)),
     &    X(ixend(ICFD)),Y(ixend(ICFD)),Z(ixend(ICFD)))
        ydis=crowxyz(X(IORG(ICFD)),Y(IORG(ICFD)),Z(IORG(ICFD)),
     &    X(iyend(ICFD)),Y(iyend(ICFD)),Z(iyend(ICFD)))
        zdis=crowxyz(X(IORG(ICFD)),Y(IORG(ICFD)),Z(IORG(ICFD)),
     &    X(izend(ICFD)),Y(izend(ICFD)),Z(izend(ICFD)))

C Set flag to indicate that origin has been defined.
        origin=.true.

      elseif(INO.EQ.5)then

C Setup regions via scan of geometry.
C Scan X,Y,Z coordinates and then sort in ascending order. Make use of
C XS,YS,ZS for array to sort and XSI,YSI,ZSI for their link back to
C initial position. Use XSS,YSS,ZSS for the compacted sorted list (with
C duplicates removed).

C Take into account if the zone is offset from the site origin.
C Use the zone bounds (ZXMN ZYMN ZZMN). Double check the max on each
C of the axis.
C Rescan to ensure G1 common blocks are correctly assigned.
        call georead(IFIL+1,LGEOM(izone),izone,1,IUOUT,IER)
        xxmax=ZXMN(izone); yymax=ZYMN(izone); zzmax=ZZMN(izone)
        do 42 il=1,NZTV(izone)
          XS(il)=X(il); YS(il)=Y(il); ZS(il)=Z(il)
          if(X(il).gt.xxmax) xxmax=X(il)
          if(Y(il).gt.yymax) yymax=Y(il)
          if(Z(il).gt.zzmax) zzmax=Z(il)
          XSS(il)=ZXMN(izone); YSS(il)=ZYMN(izone); ZSS(il)=ZZMN(izone)
          XSI(il)=il; YSI(il)=il; ZSI(il)=il
  42    continue

C Sort each axis in ascending order.
C Debug.
C        write(6,*) 'maximums found ',xxmax,yymax,zzmax
C        write(6,*) 'max geom file ',ZXMX(izone),ZYMX(izone),ZZMX(izone)
C        write(6,*) 'array x',(X(j),j=1,NZTV(izone))
C        write(6,*) 'array xs',(XS(j),j=1,NZTV(izone))
        KFLAG = 2
        call SORTR(XS,XSI,NZTV(izone),KFLAG)
C        write(6,*) 'sorted array xs',(XS(j),j=1,NZTV(izone))
C        write(6,*) 'array xsi ',(xsi(j),j=1,NZTV(izone))
        call SORTR(YS,YSI,NZTV(izone),KFLAG)
        call SORTR(ZS,ZSI,NZTV(izone),KFLAG)

C Start at smallest value for each axis. tolx,toly,tolz is how small
C a cell can be at max resolution.
        vclose=.false.
        ilx=1
        tolx=(xdis/real(NTCELX)) * 0.2
        XSS(1)=XS(1)
        ily=1
        toly=(ydis/real(NTCELY)) * 0.2
        YSS(1)=YS(1)
        ilz=1
        ZSS(1)=ZS(1)
        tolz=(zdis/real(NTCELZ)) * 0.2

C Debug.
        if(itrc.gt.0)then
          write(outs,'(a,i3,F6.3,i3,F6.3,i3,F6.3)') 
     &     'Max cells & smallest cell dimension in each axis: ',
     &      NTCELX,tolx,NTCELY,toly,NTCELZ,tolz
          call edisp(iuout,outs)
        endif
        do 43 il=2,NZTV(izone)

C For X/Y/Z axis in turn see if next sorted point is sufficiently bigger
C than the last one to qualify as a new region boundary. Note that in
C less complex zones only the latter part of XS,YS,ZS array will have
C anything other than zero.
          vclose=.false.
          CALL ECLOSE(XS(il),XSS(ilx),tolx,vclose)
          if(vclose)then
            continue
          else
            if(XS(il).gt.XSS(ilx))then
              ilx=ilx+1
              XSS(ilx)=XS(il)
            endif
          endif
          vclose=.false.
          CALL ECLOSE(YS(il),YSS(ily),toly,vclose)
          if(vclose)then
            continue
          else
            if(YS(il).gt.YSS(ily))then
              ily=ily+1
              YSS(ily)=YS(il)
            endif
          endif
          vclose=.false.
          CALL ECLOSE(ZS(il),ZSS(ilz),tolz,vclose)
          if(vclose)then
            continue
          else
            if(ZS(il).gt.ZSS(ilz))then
              ilz=ilz+1
              ZSS(ilz)=ZS(il)
            endif
          endif
  43    continue

C Notify user of the number of regions in each axis and the width of each.
        if(itrc.gt.0)then
          call edisp(iuout,'  ')
          call edisp(iuout,'  regions   points at region extents...')
          write(outs,'(a,i2,10F7.3)') 'x axis ',ilx-1,(XSS(j),j=1,ilx)
          call edisp(iuout,outs)
          write(outs,'(a,10F7.3)')    '  reg-width ',
     &      ((XSS(j+1)-XSS(j)),j=1,ilx-1)
          call edisp(iuout,outs)
          write(outs,'(a,i2,10F7.3)') 'y axis ',ily-1,(YSS(j),j=1,ily)
          call edisp(iuout,outs)
          write(outs,'(a,10F7.3)')    '  reg-width ',
     &      ((YSS(j+1)-YSS(j)),j=1,ily-1)
          call edisp(iuout,outs)
          write(outs,'(a,i2,10F7.3)') 'z axis ',ilz-1,(ZSS(j),j=1,ilz)
          call edisp(iuout,outs)
          write(outs,'(a,10F7.3)')    '  reg-width ',
     &      ((ZSS(j+1)-ZSS(j)),j=1,ilz-1)
          call edisp(iuout,outs)
        endif

C Record number of regions in each axis (one less than scan above).
        if(NREG(1,ICFD).eq.0)NREG(1,ICFD)=ilx-1
        if(NREG(2,ICFD).eq.0)NREG(2,ICFD)=ily-1
        if(NREG(3,ICFD).eq.0)NREG(3,ICFD)=ilz-1

C Initial request for total cells in X axis.
C << could do better by asking region by region >>
        write(outs,'(a,i3,a,f6.2,a)') 
     &    'How many total cells along X axis (max is ',NTCELX,
     &    ' over ',xdis,'m) would you like to have?'
        CALL EASKI(icellx,outs,' ',1,'F',NTCELX,'F',10,
     &    'requested total X cells',IERI,nbhelp)
        if(ieri.eq.-3) then
          xgrid=.FALSE.
          return
        endif
        tolx=(xdis/real(icellx))
        ntcx=0
        if(ilx.ge.2)then
          do 44 ij=1,ilx-1
            xd = XSS(ij+1)-XSS(ij)
            XREG(ij,ICFD) = xd
            NCELX(ij,ICFD) = int(xd/tolx)
            if(NCELX(ij,ICFD).eq.0)NCELX(ij,ICFD) = 1
            Xplaw(ij,ICFD) = 1.0
            ntcx = ntcx + NCELX(ij,ICFD)
  44      continue
        else
          xd = XSS(2)-XSS(1)  ! in case of one region
          XREG(1,ICFD) = xd
          NCELX(1,ICFD) = int(xd/tolx)
          if(NCELX(1,ICFD).eq.0)NCELX(1,ICFD) = 1
          Xplaw(1,ICFD) = 1.0
          ntcx = ntcx + NCELX(1,ICFD)
        endif

C Debug.
        if(itrc.gt.1)then
          write(outs,*) 'xd tolx ',xd,tolx,int(xd/tolx)
          call edisp(iuout,outs)
          write(outs,*) 'ncell in x ',ntcx,(NCELX(j,icfd),j=1,ilx-1)
          call edisp(iuout,outs)
        endif

C Initial request for total cells in Y axis.
C << alternative asking region by region >>
        write(outs,'(a,i3,a,f6.2,a)') 
     &    'How many total cells along Y axis (max is ',NTCELY,
     &    ' over ',ydis,'m) would you like to have?'
        CALL EASKI(icelly,outs,' ',1,'F',NTCELY,'F',10,
     &    'requested total Y cells',IERI,nbhelp)
        if(ieri.eq.-3) then
          xgrid=.FALSE.
          return
        endif
        toly=(ydis/real(icelly))
        ntcy=0
        if(ily.ge.2)then
          do 45 ij=1,ily-1
            yd = YSS(ij+1)-YSS(ij)
            YREG(ij,ICFD) = yd
            NCELY(ij,ICFD) = int(yd/toly)
            if(NCELY(ij,ICFD).eq.0)NCELY(ij,ICFD) = 1
            Yplaw(ij,ICFD) = 1.0
            ntcy = ntcy + NCELY(ij,ICFD)
  45      continue
        else
          yd = YSS(2)-YSS(1)  ! in case of one region
          YREG(1,ICFD) = yd
          NCELY(1,ICFD) = int(yd/toly)
          if(NCELY(1,ICFD).eq.0)NCELY(1,ICFD) = 1
          Yplaw(1,ICFD) = 1.0
          ntcy = ntcy + NCELY(1,ICFD)
        endif

C Debug.
        if(itrc.gt.1)then
          write(outs,*) 'yd toly ',yd,toly,int(yd/toly)
          call edisp(iuout,outs)
          write(outs,*) 'ncell in y ',(NCELY(j,icfd),j=1,ily-1)
          call edisp(iuout,outs)
        endif

C Initial request for total cells in Z axis.
C << alternative asking region by region ? >>
        write(outs,'(a,i3,a,f6.2,a)') 
     &    'How many total cells along Z axis (max is ',NTCELZ,
     &    ' over ',zdis,'m) would you like to have?'
        CALL EASKI(icellz,outs,' ',1,'F',NTCELZ,'F',10,
     &    'requested total Z cells',IERI,nbhelp)
        if(ieri.eq.-3) then
          xgrid=.FALSE.
          return
        endif
        tolz=(zdis/real(icellz))
        ntcz=0
        if(ilz.ge.2)then
          do 46 ij=1,ilz-1
            zd = ZSS(ij+1)-ZSS(ij)
            ZREG(ij,ICFD) = zd
            NCELZ(ij,ICFD) = int(zd/tolz)
            if(NCELZ(ij,ICFD).eq.0)NCELZ(ij,ICFD) = 1
            Zplaw(ij,ICFD) = 1.0
            ntcz = ntcz + NCELZ(ij,ICFD)
  46      continue
        else
          zd = ZSS(2)-ZSS(1)  ! in case of one region
          ZREG(1,ICFD) = zd
          NCELZ(1,ICFD) = int(zd/tolz)
          if(NCELZ(1,ICFD).eq.0)NCELZ(1,ICFD) = 1
          Zplaw(1,ICFD) = 1.0
          ntcz = ntcz + NCELZ(1,ICFD)
        endif

C Debug.
        if(itrc.gt.1)then
          write(outs,*)'Z distance ',zd,' tollerance ',tolz,int(zd/tolz)
          call edisp(iuout,outs)
          write(outs,*) 'ncell in z ',(NCELZ(j,icfd),j=1,ilz-1)
          call edisp(iuout,outs)
        endif

C Set flags to indicate that CFD origin and gridding have been defined.
        origin=.true.; xgrid=.true.; ygrid=.true.; zgrid=.true.

C Advise the user what to do next.
        helptopic='cfd_gridding_auto'
        call gethelptext(helpinsub,helptopic,nbhelp)
        call PHELPD('CFD what to do next',nbhelp,'-',0,0,IER)

      elseif(INO.EQ.7)then

C Grid X-axis.
        call GRIDDEF(1,NREG,NTCX,xdis,NCELX,XREG,Xplaw,ier)
        if (ier.ne.0) then
          griderr=.TRUE.; xgrid=.FALSE.
        else
          xgrid=.TRUE.
        endif
      elseif(INO.EQ.8)then

C Grid Y-axis.
        call GRIDDEF(2,NREG,NTCY,ydis,NCELY,YREG,Yplaw,ier)
        if (ier.ne.0) then
          griderr=.TRUE.; ygrid=.FALSE.
        else
          ygrid=.TRUE.
        endif
      elseif(INO.EQ.9)then

C Have user subdivide Z-axis into regions and specify number of cells
C and cell distribution for each region. Additional questions asked
C if curvilinear.

C Ask whether Z is Orthogonal or curvilinear.
        CALL EASKABC('Options for Z gridding:',' ','orthogonal',
     &               'curvilinear','cancel',IWO,nbhelp)
        if(iwo.ne.3)then
          ortho=.TRUE.

C Grid Zwest-axis for both cases.
          call GRIDDEF(3,NREG,NTCZ,zdis,NCELZ,ZREG,Zplaw,ier)
          if (ier.ne.0) then
            griderr=.TRUE.; zgrid=.FALSE.
          else
            zgrid=.TRUE.
          endif
          if(iwo.eq.2)then

C Additional Z axis must be gridded.
            ortho=.FALSE.
            call usrmsg(
     &        'A non-orthogonal grid requires a Z distribution',
     &        'of points on the East face of the domain. ','W')
            if(ieorg.eq.0) ieorg=1
            CALL EASKI(ieorg,'Specify vertex at start of east Z axis.',
     &        'Which `v`? ',1,'F',NTV,'F',1,'ez axis vertex',
     &        IERI,nbhelp)
            if(ieri.eq.-3) then
              zegrid=.FALSE.
              return
            endif
            izendet=izende(ICFD)
            CALL EASKI(izendet,'Specify vertex at end of east Z axis.',
     &        'Which `v`? ',1,'F',NTV,'F',1,'ez axis vertex',
     &        IERI,nbhelp)
            if(ieri.eq.-3) then
              zegrid=.FALSE.
              return
            endif

C User did not cancel so instanciate values.
            izende(ICFD)=izendet

            zedis=crowxyz(X(ieorg),Y(ieorg),Z(ieorg),X(izende(ICFD)),
     &              Y(izende(ICFD)),Z(izende(ICFD)))

C Need to ask for vertices which represent the
C top and bottom boundaries of the curvilinear x-section.
            call TOPBOT(izone,'-',IER)

C Grid Zeast-axis.
            call GRIDDEF(4,NREG,NTCZe,zedis,NCELZe,ZREGe,Zplawe,ier)
            if (ier.ne.0) then
              griderr=.TRUE.; zegrid=.FALSE.
            else
              zegrid=.TRUE.
            endif
          endif
        endif

      elseif(INO.EQ.(10+NCP).and.(.NOT.griderr))then

C Visualize gridding.
        if(origin.and.xgrid.and.ygrid.and.zgrid)then
          call cgd(izone,ierr)

C << Functionality changed so that grid visualisation chosen in previous
C    menu is persistent to this and the root CFD menu. >>
        else
          CALL USRMSG(' Define origin of CFD domain ',
     &         ' and grid axes first ','W')
        endif

      elseif(INO.EQ.NITMS-1)then

        helptopic='cfd_gridding_setup'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL PHELPD('gridding',nbhelp,'-',0,0,IER)

      elseif(INO.EQ.NITMS-2)then

        itrc=itrc+1           ! toggle feedback level
        if(itrc.gt.2) itrc=0

      elseif(INO.EQ.NITMS)then

C Return to higher level menu.
        if (griderr) CALL USRMSG(
     &    'The current gridding scheme is incomplete!','  ','W')

C Check grid definition and current velocity solution requirements are 
C consistant.
        if (CALLU(ICFD).OR.CALLV(ICFD).OR.CALLW(ICFD)) then
          if (NTCX.eq.1) then
            CALLU(ICFD)=.FALSE.
          else
            CALLU(ICFD)=.TRUE.
          endif
          if (NTCY.eq.1) then
            CALLV(ICFD)=.FALSE.
          else
            CALLV(ICFD)=.TRUE.
          endif
          if (NTCZ.eq.1) then
            CALLW(ICFD)=.FALSE.
          else
            CALLW(ICFD)=.TRUE.
          endif
        endif
        return
      endif
      goto 14

      end
    

C ************************ CFDSVAR ************************
C CFDSVAR: controls the editing of CFD solution variables.
C IER=0 indicates no error.

      SUBROUTINE CFDSVAR(IER)
#include "cfd.h"

      COMMON/ICFNOD/ICFD,ICP
      common/grdmax/NTCX,NTCY,NTCZ

C Equations solved.
      common/EQTION/CALLU(MNZ),CALLV(MNZ),CALLW(MNZ),CALLT(MNZ),
     &             CALLC(MNZ),KEMDL(MNZ),BUOY(MNZ),BOUSSI(MNZ),
     &             ZEROT(MNZ),ZandKE(MNZ),MITzero(MNZ)
      common/EQTION3/CALLMA(MNZ),CALPOL(MCTM,MNZ),POLNAM(MCTM,MNZ),
     &               NCTM(MNZ),JHUMINDX(MNZ),URFC(MCTM)

C Solution methods.
      common/METHDS/ITURB(MNZ),IBUOY(MNZ)

      character ITEM(11+MCTM)*40, POLNAM*12
      LOGICAL CALLU,CALLV,CALLW,CALLT,CALLC,KEMDL,BUOY,BOUSSI
      LOGICAL ZEROT,ZandKE,MITzero
      LOGICAL CALPOL,CALLMA,OK,dok
      integer NITMS,INO ! max items and current menu item

C For help messages
      character helpinsub*24 ! subroutine name
      character helptopic*24 ! string (unique) for topic
      integer nbhelp     ! number of help lines found

      helpinsub='CFDSVAR'  ! set for subroutine

C Check gridding to set velocity solution.
      if (NTCX.gt.1) then
        CALLU(ICFD)=.TRUE.
      else
        CALLU(ICFD)=.FALSE.
      endif
      if (NTCY.gt.1) then
        CALLV(ICFD)=.TRUE.
      else
        CALLV(ICFD)=.FALSE.
      endif
      if (NTCZ.gt.1) then
        CALLW(ICFD)=.TRUE.
      else
        CALLW(ICFD)=.FALSE.
      endif

C Create menu showing 4 main solution categories and allow togelling 
C to activate various solution cats/ methods.  Once a solution category 
C has been activated then display associated state variables and allow 
C definition of initial values and relaxation factors.
C First clear menu.
 5    do 10 I=1,14
        ITEM(I)='  '
 10   continue

C Clear display area.
      call usrmsg('  ','  ','-')

C Now create menu.
      ITEM(1)= '  Velocity: solved'
      if (.NOT.CALLT(ICFD)) then
        ITEM(2)= 'a Temperature >> not solved'
      else
        ITEM(2)= 'a Temperature >> solved'
      endif
      if (ITURB(ICFD).eq.0) then
        ITEM(3)= 'b Turbulence >> none (laminar)'
      elseif (ITURB(ICFD).eq.1) then
        ITEM(3)= 'b Turbulence >> k-e turbulence'
      elseif (ITURB(ICFD).eq.2) then
        ITEM(3)= 'b Turbulence >> fixed eddy viscosity'
      elseif (ITURB(ICFD).eq.3) then
        ITEM(3)= 'b Turbulence >> MIT zero-equation'
      elseif (ITURB(ICFD).eq.4) then
        ITEM(3)= 'b Turbulence >> fixed visc -> k-e'
      endif
      if (IBUOY(ICFD).eq.0) then
        ITEM(4)= 'c Buoyancy >> not solved'
      elseif (IBUOY(ICFD).eq.1) then
        ITEM(4)= 'c Buoyancy >> ideal gas'
      elseif (IBUOY(ICFD).eq.2) then
        ITEM(4)= 'c Buoyancy >> Boussinesq approx'
      endif
      DO 121 ICTM=1,NCTM(ICFD)
        IF(CALPOL(ICTM,ICFD))THEN
          WRITE(ITEM(4+ICTM),'(A,1X,A,1X,A,1X)')CHAR(99+ICTM),
     &         POLNAM(ICTM,ICFD),'>> solved'
        ELSE
          WRITE(ITEM(4+ICTM),'(A,1X,A,1X,A,1X)')CHAR(99+ICTM),
     &         POLNAM(ICTM,ICFD),'>> not solved'
        ENDIF
 121  CONTINUE
      if (.NOT.CALLMA(ICFD)) then
        WRITE(ITEM(5+NCTM(ICFD)),'(2A)')CHAR(100+NCTM(ICFD)),
     &       ' Mean age of air >> not solved'
      else
        WRITE(ITEM(5+NCTM(ICFD)),'(2A)')CHAR(100+NCTM(ICFD)),
     &       ' Mean age of air >> solved'
      endif
      ITEM(6+NCTM(ICFD))= ' ------------------------------------- '
      ITEM(7+NCTM(ICFD))= '# Edit solution parameters'
      ITEM(8+NCTM(ICFD))= '* Reset to defaults'
      ITEM(9+NCTM(ICFD))= ' ------------------------------------- '
      ITEM(10+NCTM(ICFD))='? Help'
      ITEM(11+NCTM(ICFD))='- exit this menu'
      NITMS=11+NCTM(ICFD)

C Help text for this menu.
      helptopic='cfd_equations_to_solve'
      call gethelptext(helpinsub,helptopic,nbhelp)

C Display menu.
      CALL EMENU('Required solutions',ITEM,NITMS,INO)

C Check for pick of empty box.
      if(INO.EQ.2)then

C Toggle temperature solution.
        if (.NOT.CALLT(ICFD)) then
          CALLT(ICFD)=.TRUE.
        else
          CALLT(ICFD)=.FALSE.
        endif
      elseif(INO.EQ.3)then

C Toggle through turbulence methods.
        ITURB(ICFD)=ITURB(ICFD)+1
        if (ITURB(ICFD).gt.4) ITURB(ICFD)=0
      elseif(INO.EQ.4)then

C Toggle buoyancy.
        IBUOY(ICFD)=IBUOY(ICFD)+1
        if (IBUOY(ICFD).gt.2) IBUOY(ICFD)=0
      elseif(INO.GT.4.AND.INO.LE.(4+NCTM(ICFD)))then

C Contaminants
        IF(CALPOL(INO-4,ICFD))THEN
          CALPOL(INO-4,ICFD)=.FALSE.
        ELSE
          CALPOL(INO-4,ICFD)=.TRUE.
        ENDIF
      elseif(INO.EQ.5+NCTM(ICFD))then

C Mean age of air.
        if (CALLMA(ICFD)) then
          CALLMA(ICFD)=.false.
        else
          CALLMA(ICFD)=.true.
        endif
      elseif(INO.EQ.8+NCTM(ICFD))then

C Reset defaults.
        dok=.false.
        call askok('Resetting solution parameters to default values.',
     &             'Are you sure? (see help)',OK,dok,nbhelp)
        if (OK) call CFDDEFLT
      elseif(INO.eq.7+NCTM(ICFD))then

C Edit solution parameters.
        call CFDEDAT        
      elseif(INO.eq.(NITMS-1)) then

C Help.
        helptopic='cfd_equations_to_solve'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL PHELPD('solution variables',nbhelp,'-',0,0,IER)
      elseif(INO.eq.NITMS) then
        return
      endif
      goto 5

      end

C ************************ CFDEDAT ************************
C CFDEDAT: controls the editing of CFD equations.

      SUBROUTINE CFDEDAT
C #include "building.h"
#include "cfd.h"

C      common/OUTIN/IUOUT,IUIN
      COMMON/ICFNOD/ICFD,ICP
      common/grdmax/NTCX,NTCY,NTCZ

C Equations solved.
      common/EQTION/CALLU(MNZ),CALLV(MNZ),CALLW(MNZ),CALLT(MNZ),
     &             CALLC(MNZ),KEMDL(MNZ),BUOY(MNZ),BOUSSI(MNZ),
     &             ZEROT(MNZ),ZandKE(MNZ),MITzero(MNZ)
      common/EQTION3/CALLMA(MNZ),CALPOL(MCTM,MNZ),POLNAM(MCTM,MNZ),
     &               NCTM(MNZ),JHUMINDX(MNZ),URFC(MCTM)

C Solution methods.
      common/METHDS/ITURB(MNZ),IBUOY(MNZ)

C Initial values.
      COMMON/INITIA/UINIT(MNZ),VINIT(MNZ),WINIT(MNZ),PINIT(MNZ),
     &              TINIT(MNZ),TEINIT(MNZ),EDINIT(MNZ),POLINIT(MNZ,MCTM)

C Under-relaxation factors.
      COMMON/LINRFC/URFCU(MNZ),URFCV(MNZ),URFCW(MNZ),URFCP(MNZ),
     &              URFCT(MNZ),URFCK(MNZ),URFCE(MNZ),URFCVS(MNZ),
     &              URFCC(MNZ,MCTM)
      COMMON/LINRFC2/URFCU2(MNZ),URFCV2(MNZ),URFCW2(MNZ),URFCP2(MNZ),
     &              URFCT2(MNZ),URFCK2(MNZ),URFCE2(MNZ),URFCVS2(MNZ),
     &              URFCC2(MNZ,MCTM)

C Convergence criteria.
      common/param1/MAXITR(MNZ),IMONT(MNZ),JMONT(MNZ),KMONT(MNZ),
     &             IPPHI(MNZ),SRMAX(MNZ)

C Additional data needed for some models.
      COMMON/ZTURB/rMOOT(MNZ),nZtoKE(MNZ)
      COMMON/BUOYAN/BUOYA,BOUSSA,TBAR(MNZ)
 
      character ITEM(33)*40,ltmp*16,dtmp*16, POLNAM*12

      LOGICAL CALLU,CALLV,CALLW,CALLT,CALLC,KEMDL,BUOY,BOUSSI
      LOGICAL ZEROT,ZandKE,MITzero
      LOGICAL BUOYA,BOUSSA
      LOGICAL CALPOL,CALLMA
      integer NITMS,INO ! max items and current menu item

C For help messages
      character helpinsub*24 ! subroutine name
      character helptopic*24 ! string (unique) for topic
      integer nbhelp     ! number of help lines found

      helpinsub='CFDEDAT'  ! set for subroutine

C Check gridding to set velocity solution.
      if (NTCX.gt.1) then
        CALLU(ICFD)=.TRUE.
      else
        CALLU(ICFD)=.FALSE.
      endif
      if (NTCY.gt.1) then
        CALLV(ICFD)=.TRUE.
      else
        CALLV(ICFD)=.FALSE.
      endif
      if (NTCZ.gt.1) then
        CALLW(ICFD)=.TRUE.
      else
        CALLW(ICFD)=.FALSE.
      endif

C Create menu showing 4 main solution categories and allow togelling 
C to activate various solution cats/ methods.  Once a solution category 
C has been activated then display associated state variables and allow 
C definition of initial values and relaxation factors.
C First clear menu.
 5    INO=-4
      do 10 I=1,31
        ITEM(I)='  '
 10   continue

C Clear display area.
      call usrmsg('  ','  ','-')

C Now create menu.
      ITEM(1)= ' ----- Velocity in X direction ----- '
      if (CALLU(ICFD)) then
        write(ITEM(2),'(a,f6.3)')'a Initial value: ',UINIT(ICFD)
        write(ITEM(3),'(a,f5.2)')'b Relaxation factor 1: ',URFCU(ICFD)
        write(ITEM(4),'(a,f5.2)')'c Relaxation factor 2: ',URFCU2(ICFD)
      else
        ITEM(2)='  N/A (not solved)'
        ITEM(3)='  N/A (not solved)'
        ITEM(4)='  N/A (not solved)'
      endif
      ITEM(5)= ' ----- Velocity in Y direction ----- '
      if (CALLV(ICFD)) then
        write(ITEM(6),'(a,f6.3)')'d Initial value: ',VINIT(ICFD)
        write(ITEM(7),'(a,f5.2)')'e Relaxation factor 1: ',URFCV(ICFD)
        write(ITEM(8),'(a,f5.2)')'f Relaxation factor 2: ',URFCV2(ICFD)
      else
        ITEM(6)='  N/A (not solved)'
        ITEM(7)='  N/A (not solved)'
        ITEM(8)='  N/A (not solved)'
      endif
      ITEM(9)= ' ----- Velocity in Z direction ----- '
      if (CALLW(ICFD)) then
        write(ITEM(10),'(a,f6.3)')'g Initial value: ',WINIT(ICFD)
        write(ITEM(11),'(a,f5.2)')'h Relaxation factor 1: ',URFCW(ICFD)
        write(ITEM(12),'(a,f5.2)')'i Relaxation factor 2: ',URFCW2(ICFD)
      else
        ITEM(10)='  N/A (not solved)'
        ITEM(11)='  N/A (not solved)'
        ITEM(12)='  N/A (not solved)'
      endif
      ITEM(13)=' ----------- Temperature ----------- '
      if (CALLT(ICFD)) then
        ITEM(14)='j Temperature >> solved'
      else
        ITEM(14)='j Temperature >> not solved'
      endif
      write(ITEM(15),'(a,f6.3)')'k Initial value: ',TINIT(ICFD)
      write(ITEM(16),'(a,f5.2)')'l Relaxation factor 1: ',URFCT(ICFD)
      write(ITEM(17),'(a,f5.2)')'m Relaxation factor 2: ',URFCT2(ICFD)
      ITEM(18)=' ----------- Turbulence ------------ '
      if (ITURB(ICFD).eq.0) then
        ITEM(19)= 'n Turbulence >> none (laminar)'
      elseif (ITURB(ICFD).eq.1) then
        ITEM(19)= 'n Turbulence >> k-e turbulence'
        write(ITEM(20),'(a,E8.2E2,a,E8.2E2)')'o Initial; k: ',
     &                                 TEINIT(ICFD),' e: ',EDINIT(ICFD)
        write(ITEM(21),'(a,f5.2,a,f5.2)')'p Relax. factor 1; k:',
     &                                   URFCK(ICFD),' e: ',URFCE(ICFD)
        write(ITEM(22),'(a,f5.2,a,f5.2)')'q Relax. factor 2; k:',
     &                                 URFCK2(ICFD),' e: ',URFCE2(ICFD)
      elseif (ITURB(ICFD).eq.2) then
        ITEM(19)= 'n Turbulence >> fixed eddy viscosity'
        write (ITEM(20),'(a,f6.2)') 'o Eddy viscosity: ',rMOOT(ICFD)
      elseif (ITURB(ICFD).eq.3) then
        ITEM(19)= 'n Turbulence >> MIT zero-equation'
      elseif (ITURB(ICFD).eq.4) then
        ITEM(19)= 'n Turbulence >> fixed visc -> k-e'
        write(ITEM(20),'(a,f6.2)')'o Eddy viscosity: ',rMOOT(ICFD)
        write(ITEM(21),'(a,i5)')'p Transition on iteration: ',
     &                                                     nZtoKE(ICFD)
        write(ITEM(22),'(a)')'q Toggle to k-e to set parameters'
      endif
      ITEM(23)=' ------------- Buoyancy ------------ '
      if (IBUOY(ICFD).eq.0) then
        ITEM(24)= 'r Buoyancy >> not solved'
      elseif (IBUOY(ICFD).eq.1) then
        ITEM(24)= 'r Buoyancy >> ideal gas'
      elseif (IBUOY(ICFD).eq.2) then
        ITEM(24)= 'r Buoyancy >> Boussinesq approx.'
        write(ITEM(25),'(a,f6.2)')'s Reference temperature: ',
     &                                                       TBAR(ICFD)
      endif
      ITEM(26)=' ------- Convergence criteria ------ '
      write(ITEM(27),'(a,i5)')  't Maximum iterations: ',MAXITR(ICFD)
      write(ITEM(28),'(a,f7.5)')'u Max sum of residuals: ',SRMAX(ICFD)
      write(ITEM(29),'(a,3i4)') 'v Monitoring cell (i,j,k): ',
     &                        IMONT(ICFD)-1,JMONT(ICFD)-1,KMONT(ICFD)-1
      ITEM(30)=' ----------------------------------- '
      write(ITEM(31),'(a)')'+ Contaminants'
      ITEM(32)='? Help'
      ITEM(33)='- Exit this menu'

      NITMS=33

C Help text for this menu.
      helptopic='cfd_equations_to_edit'
      call gethelptext(helpinsub,helptopic,nbhelp)

C Display menu.
      CALL EMENU('Solution parameters',ITEM,NITMS,INO)

C Check for pick of empty box.
      if(INO.EQ.2)then

C Ask for initial value of velocity in X direction.
        VAL=UINIT(ICFD)
        CALL EASKR(VAL,'Initial value for X velocity?','(m/s)',-1.0,
     &             'W',1.0,'W',0.001,'U init',IER,nbhelp)
        if (IER.eq.0) UINIT(ICFD)=VAL
      elseif(INO.EQ.3)then

C Ask for relaxation factor 1 for velocity in X direction.
        VAL=URFCU(ICFD)
        CALL EASKR(VAL,'Relaxation factor 1 (see help) for',
     &   'X velocity?',0.0,'F',1.5,'W',0.5,'U relax 1',IER,nbhelp)
        if (IER.eq.0) URFCU(ICFD)=VAL
      elseif(INO.EQ.4)then

C Ask for relaxation factor 2 for velocity in X direction.
        VAL=URFCU2(ICFD)
        CALL EASKR(VAL,'Relaxation factor 2 (see help) for',
     &   'X velocity?',0.0,'F',1.5,'W',0.05,'U relax 2',IER,nbhelp)
        if (IER.eq.0) URFCU2(ICFD)=VAL
      elseif(INO.EQ.6)then

C Ask for initial value of velocity in Y direction.
        VAL=VINIT(ICFD)
        CALL EASKR(VAL,'Initial value for Y velocity?','(m/s)',-1.0,
     &             'W',1.0,'W',0.001,'V init',IER,nbhelp)
        if (IER.eq.0) VINIT(ICFD)=VAL
      elseif(INO.EQ.7)then

C Ask for relaxation factor 1 for velocity in Y direction.
        VAL=URFCV(ICFD)
        CALL EASKR(VAL,'Relaxation factor 1 (see help) for',
     &   'Y velocity?',0.0,'F',1.5,'W',0.5,'V relax 1',IER,nbhelp)
        if (IER.eq.0) URFCV(ICFD)=VAL
      elseif(INO.EQ.8)then

C Ask for relaxation factor 2 for velocity in Y direction.
        VAL=URFCV2(ICFD)
        CALL EASKR(VAL,'Relaxation factor 2 (see help) for',
     &   'Y velocity?',0.0,'F',1.5,'W',0.05,'V relax 2',IER,nbhelp)
        if (IER.eq.0) URFCV2(ICFD)=VAL
      elseif(INO.eq.10)then

C Ask for initial value of velocity in Z direction.
        VAL=WINIT(ICFD)
        CALL EASKR(VAL,'Initial value for Z velocity?','(m/s)',-1.0,
     &            'W',1.0,'W',0.001,'W init',IER,nbhelp)
        if (IER.eq.0) WINIT(ICFD)=VAL
      elseif(INO.EQ.11)then

C Ask for relaxation factor 1 for velocity in Z direction.
        VAL=URFCW(ICFD)
        CALL EASKR(VAL,'Relaxation factor 1 (see help) for',
     &   'Z velocity?',0.0,'F',1.5,'W',0.5,'W relax 1',IER,nbhelp)
        if (IER.eq.0) URFCW(ICFD)=VAL
      elseif(INO.EQ.12)then

C Ask for relaxation factor 2 for velocity in Z direction.
        VAL=URFCW2(ICFD)
        CALL EASKR(VAL,'Relaxation factor 2 (see help) for',
     &   'Z velocity?',0.0,'F',1.5,'W',0.05,'W relax 2',IER,nbhelp)
        if (IER.eq.0) URFCW2(ICFD)=VAL
      elseif(INO.eq.14)then      

C (De)activate temperature solution.
        if (.NOT.CALLT(ICFD)) then
          CALLT(ICFD)=.TRUE.
        else
          CALLT(ICFD)=.FALSE.
        endif
      elseif(INO.eq.15)then

C Ask for initial value of temperature.
        VAL=TINIT(ICFD)
        CALL EASKR(VAL,'Initial value for temperature?','(degrees C)',
     &                         10.,'W',30.,'W',20.,'T init',IER,nbhelp)
        if (IER.eq.0) TINIT(ICFD)=VAL
      elseif(INO.EQ.16)then

C Ask for relaxation factor 1 for temperature.
        VAL=URFCT(ICFD)
        CALL EASKR(VAL,'Relaxation factor 1 (see help) for',
     &   'temperature?',0.0,'F',1.5,'W',1.0,'T relax 1',IER,nbhelp)
        if (IER.eq.0) URFCT(ICFD)=VAL
      elseif(INO.EQ.17)then

C Ask for relaxation factor 2 for temperature.
        VAL=URFCT2(ICFD)
        CALL EASKR(VAL,'Relaxation factor 2 (see help) for',
     &   'temperature?',0.0,'F',1.5,'W',0.25,'T relax 2',IER,nbhelp)
        if (IER.eq.0) URFCT2(ICFD)=VAL
      elseif(INO.eq.19)then   

C Toggle through turbulence methods.
        ITURB(ICFD)=ITURB(ICFD)+1
        if (ITURB(ICFD).gt.4) ITURB(ICFD)=0
      elseif(INO.EQ.20)then

C Edit turbulence data according to current method.
        if (ITURB(ICFD).eq.1) then

C Ask for initial values for k and epsilon.
          VAL=TEINIT(ICFD)
          CALL EASKE(VAL,'Initial value for turbulent energy?',' ',
     &                      0.0,'W',0.1,'W',0.005,'TE init',IER,nbhelp)
          if (IER.eq.0) TEINIT(ICFD)=VAL
          VAL=EDINIT(ICFD)
          CALL EASKE(VAL,'Initial value for energy dissipation?',' ',
     &                      0.0,'W',0.1,'W',0.005,'ED init',IER,nbhelp)
          if (IER.eq.0) EDINIT(ICFD)=VAL
        elseif ((ITURB(ICFD).eq.2).or.(ITURB(ICFD).eq.4)) then

C Ask for fixed eddy viscosity.
          VAL=rMOOT(ICFD)
          CALL EASKR(VAL,'Fixed value for eddy viscosity?',' ',0.0,'F',
     &               500.0,'F',90.0,'eddy viscos',IER,nbhelp)
          if (IER.eq.0) rMOOT(ICFD)=VAL
        endif
      elseif(INO.EQ.21)then

C Edit turbulence data according to current method.
        if (ITURB(ICFD).eq.1) then

C Ask for relaxation factor 1 for k and epsilon.
          VAL=URFCK(ICFD)
          CALL EASKR(VAL,'Relaxation factor 1 (see help) for',
     & 'turbulent energy?',0.0,'F',1.5,'W',1.0,'TE relax 1',IER,nbhelp)
          if (IER.eq.0) URFCK(ICFD)=VAL
          VAL=URFCE(ICFD)
          CALL EASKR(VAL,'Relaxation factor 1 (see help) for',
     &      'energy dissipation?',0.0,'F',1.5,'W',1.0,'ED relax 1',IER,
     &                                                          nbhelp)
          if (IER.eq.0) URFCE(ICFD)=VAL
        elseif (ITURB(ICFD).eq.4) then

C Ask for iteration to transition from fixed eddy viscosity to k-e.
          IVAL=nZtoKE(ICFD)
          CALL EASKI(IVAL,'Iteration number to transition to',
     &      'k-e model?',1,'F',1000,'W',10,'fev->ke iter',IER,nbhelp)
          if (IER.eq.0) nZtoKE(ICFD)=IVAL
        endif
      elseif(INO.EQ.22)then

C Edit turbulence data according to current method.
        if (ITURB(ICFD).eq.1) then

C Ask for relaxation factor 2 for k and epsilon.
          VAL=URFCK2(ICFD)
          CALL EASKR(VAL,'Relaxation factor 2 (see help) for',
     &       'turbulent energy?',0.0,'F',1.5,'W',0.05,'TE relax 2',IER,
     &                                                          nbhelp)
          if (IER.eq.0) URFCK2(ICFD)=VAL
          VAL=URFCE2(ICFD)
          CALL EASKR(VAL,'Relaxation factor 2 (see help) for',
     &     'energy dissipation?',0.0,'F',1.5,'W',0.05,'ED relax 2',IER,
     &                                                          nbhelp)
          if (IER.eq.0) URFCE2(ICFD)=VAL
        elseif (ITURB(ICFD).eq.4) then

C Toggle to k-e.
          ITURB(ICFD)=1
        endif
      elseif(INO.EQ.24)then

C Toggle buoyancy.
        IBUOY(ICFD)=IBUOY(ICFD)+1
        if (IBUOY(ICFD).gt.2) IBUOY(ICFD)=0
      elseif(INO.EQ.25)then
        if (IBUOY(ICFD).eq.2) then

C Ask for reference temperature.
          VAL=TBAR(ICFD)
          CALL EASKR(VAL,'Reference temperature?','(degrees C)?',10.,
     &                               'W',30.,'W',20.,'Tbar',IER,nbhelp)
          if (IER.eq.0) TBAR(ICFD)=VAL
        endif
      elseif(INO.EQ.27)then

C Ask for maximum number of iterations.
        helptopic='cfd_convergence_crit'
        call gethelptext(helpinsub,helptopic,nbhelp)
        IVAL=MAXITR(ICFD)
        CALL EASKI(IVAL,'Maximum number of iterations to perform?',' ',
     &             1,'F',10000,'W',3750,'CFD iters',IER,nbhelp)
        if (IER.eq.0) MAXITR(ICFD)=IVAL
      elseif (INO.eq.28) then

C Ask for maximum residual.
        helptopic='cfd_convergence_crit'
        call gethelptext(helpinsub,helptopic,nbhelp)
        VAL=SRMAX(ICFD)
        CALL EASKR(VAL,' ','Maximum equation residual?',0.,'F',1.,'F',
     &                                     0.02,'max resid',IER,nbhelp)
        if (IER.eq.0) SRMAX(ICFD)=VAL
      elseif (INO.eq.29) then

C Monitoring cell.
        write (ltmp,'(3i4)') IMONT(ICFD)-1,JMONT(ICFD)-1,KMONT(ICFD)-1
        write (dtmp,'(a)') ' 1  1  1'
        call EASKS(ltmp,'Monitoring cell (i,j,k)? ','  ',16,
     &                                    dtmp,'monit cell',IER,nbhelp)
        K=0
        call EGETWI(ltmp,K,IVAL,0,NTCELX,'W','I mon cell',IER)
        if (IER.eq.0) IMONT(ICFD)=IVAL+1
        call EGETWI(ltmp,K,IVAL,0,NTCELY,'W','J mon cell',IER)
        if (IER.eq.0) JMONT(ICFD)=IVAL+1
        call EGETWI(ltmp,K,IVAL,0,NTCELZ,'W','K mon cell',IER)
        if (IER.eq.0) KMONT(ICFD)=IVAL+1
      elseif (INO.eq.NITMS-2) then

C Check that contaminants have been defined.
        if (NCTM(ICFD).gt.0) then

C Go to contaminant equation setup menu.
          CALL CFDEDATC
        else

C Display popup.
          helptopic='cfd_no_contam_defined'
          call gethelptext(helpinsub,helptopic,nbhelp)
          CALL PHELPD('No contaminants defined',nbhelp,'-',0,0,IER)
        endif        
      elseif(INO.eq.(NITMS-1)) then

C Help.
        helptopic='cfd_equations_to_edit'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL PHELPD('Solution parameters',nbhelp,'-',0,0,IER)
      elseif(INO.eq.NITMS) then
        return
      endif
      goto 5

      end

C ************************* CFDEDATC **************************
C CFDEDATC is a continuation of CFDEDAT, controlling editing
C of CFD contaminant equations.

      SUBROUTINE CFDEDATC
C #include "building.h"
#include "cfd.h"

C      common/OUTIN/IUOUT,IUIN
      COMMON/ICFNOD/ICFD,ICP

C Equations solved.
      common/EQTION3/CALLMA(MNZ),CALPOL(MCTM,MNZ),POLNAM(MCTM,MNZ),
     &               NCTM(MNZ),JHUMINDX(MNZ),URFC(MCTM)

C Under-relaxation factors.
      COMMON/LINRFC/URFCU(MNZ),URFCV(MNZ),URFCW(MNZ),URFCP(MNZ),
     &              URFCT(MNZ),URFCK(MNZ),URFCE(MNZ),URFCVS(MNZ),
     &              URFCC(MNZ,MCTM)
      COMMON/LINRFC2/URFCU2(MNZ),URFCV2(MNZ),URFCW2(MNZ),URFCP2(MNZ),
     &              URFCT2(MNZ),URFCK2(MNZ),URFCE2(MNZ),URFCVS2(MNZ),
     &              URFCC2(MNZ,MCTM)

      character ITEM(3+3*NCTM(ICFD))*40,POLNAM*12,tstr*13,tstr1*30,
     &  tstr2*30

      LOGICAL CALPOL,CALLMA
      integer NITMS,INO ! max items and current menu item

C For help messages
      character helpinsub*24 ! subroutine name
      character helptopic*24 ! string (unique) for topic
      integer nbhelp     ! number of help lines found

      helpinsub='CFDEDATC'  ! set for subroutine

      NITMS=3+3*NCTM(ICFD)
C Debug
C      write(6,*)'NITMS: ',NITMS

C First clear menu.
 5    INO=-4
      do I=1,NITMS
        ITEM(I)='  '
      enddo

C Now create menu.
      do ICTM=1,NCTM(ICFD)
        IND=(ICTM-1)*3+1
C Debug
C        write(6,*)'IND: ',IND
        IF(CALPOL(ICTM,ICFD))THEN
          WRITE(ITEM(IND),'(A,1X,A,1X,A,1X)')CHAR(96+IND),
     &       POLNAM(ICTM,ICFD)(1:lnblnk(POLNAM(ICTM,ICFD))),'>> solved'
        ELSE
          WRITE(ITEM(IND),'(A,1X,A,1X,A,1X)')CHAR(96+IND),
     &   POLNAM(ICTM,ICFD)(1:lnblnk(POLNAM(ICTM,ICFD))),'>> not solved'
        ENDIF
        IND=(ICTM-1)*3+2
C Debug
C        write(6,*)'IND: ',IND
        write(ITEM(IND),'(a,1x,a,f5.2)')CHAR(96+IND),
     &                         'Relaxation factor 1: ',URFCC(ICFD,ICTM)
        IND=(ICTM-1)*3+3
C Debug
C        write(6,*)'IND: ',IND
        write(ITEM(IND),'(a,1x,a,f5.2)')CHAR(96+IND),
     &                        'Relaxation factor 2: ',URFCC2(ICFD,ICTM)
      enddo
      ITEM(NITMS-2)=' ------------------------------ '
      ITEM(NITMS-1)='? Help'
      ITEM(NITMS)=  '- Exit this menu'

C Help text for this menu.
      helptopic='cfd_contam_equ_to_edit'
      call gethelptext(helpinsub,helptopic,nbhelp)

C Display menu.
      CALL EMENU('Contaminant relaxation factors',ITEM,NITMS,INO)

C Process INO to get ICTM and option 1, 2 or 3.
      ICTM=(INO+2)/3
      IND=INO-(ICTM-1)*3
C Debug
C      write(6,*)'ICTM: ',ICTM
C      write(6,*)'IND: ',IND
C Create temporary strings.
      write(tstr,'(2a)')POLNAM(ICTM,ICFD)(1:lnblnk(POLNAM(ICTM,ICFD))),
     &                                                              '?'
      write(tstr1,'(a,i1,a)')'contam ',ICTM,' relax 1'
      write(tstr2,'(a,i1,a)')'contam ',ICTM,' relax 2'
      
C First check for Help or Exit.
      IF(INO.eq.(NITMS-1))THEN

C Help text for this menu.
        helptopic='cfd_contam_equ_to_edit'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL PHELPD('Solution parameters (contam.)',nbhelp,'-',0,0,IER)
      ELSEIF(INO.eq.NITMS)THEN

C Exit.
        return
      ELSEIF(IND.eq.1)THEN

C Toggle contaminant equation.
        if (CALPOL(ICTM,ICFD)) then
          CALPOL(ICTM,ICFD)=.false.
        else
          CALPOL(ICTM,ICFD)=.true.
        endif
      ELSEIF(IND.eq.2)THEN

C Ask for relaxation factor 1.
        VAL=URFCC(ICFD,ICTM)
        CALL EASKR(VAL,'Relaxation factor 1 (see help) for',
     &                       tstr,0.0,'F',1.5,'W',1.0,tstr1,IER,nbhelp)
        if (IER.eq.0) URFCC(ICFD,ICTM)=VAL
      ELSEIF(IND.eq.3)THEN

C Ask for relaxation factor 2.
        VAL=URFCC2(ICFD,ICTM)
        CALL EASKR(VAL,'Relaxation factor 2 (see help) for',
     &                      tstr,0.0,'F',1.5,'W',0.25,tstr2,IER,nbhelp)
        if (IER.eq.0) URFCC2(ICFD,ICTM)=VAL
      ENDIF

      goto 5
      END      

C ************************** CFDBVAR **************************
C CFDBVAR controls the editing of CFD boundary variables.
C ITRU unit number for user output.
C IER=0 indicates no error.

      SUBROUTINE CFDBVAR(IZONE,IER)
#include "building.h"
#include "cfd.h"
#include "epara.h"
#include "prj3dv.h"

      COMMON/SPAD/MMOD,LIMIT,LIMTTY

      COMMON/ICFNOD/ICFD,ICP
      COMMON/OUTIN/IUOUT,IUIN   
      common/cfdconf/ICFBLD(MNZ),ICFMFS(MNZ)
      common/KEYVOLS/NVOL(MNZ),IVOLF(MNVLS,MNZ),IVCELLS(MNVLS,MNZ,2),
     &               JVCELLS(MNVLS,MNZ,2),KVCELLS(MNVLS,MNZ,2)
      common/KEYVOLN/VOLNAME(MNVLS,MNZ),VCsurf(MNVLS,MNZ)
      common/KEYVDAT/IVTYPE(MNVLS,MNZ),VOLTemp(MNVLS,MNZ),
     &          VOLHeat(MNVLS,MNZ),IVConfl(MNVLS,MNZ),VOLHum(MNVLS,MNZ),
     &          VOLCO2(MNVLS,MNZ),VOLVel(MNVLS,MNZ),VOLDir(MNVLS,MNZ,2),
     &          VOLArea(MNVLS,MNZ),VOLPres(MNVLS,MNZ),
     &          VOLPol(MCTM,MNVLS,MNZ)
      common/GRIDEFN/origin,xgrid,ygrid,zgrid,zegrid,ortho
      COMMON/CFDVIS/HAS_GEOM,ISHSB,ISHAO,IFACES,ISHBLK,ISHGEO,INITD
      logical HAS_GEOM
      character INITD*6
      common/param2/TITLE(MNZ),CFTRFL(MNZ),LPHI(MNZ)
      COMMON/NDMAP/NOPEN(MNZ),MFNODE(MCFND,MNZ),IOPENi(MCFND,MNZ),
     &             IOPENf(MCFND,MNZ),JOPENi(MCFND,MNZ),
     &             JOPENf(MCFND,MNZ),KOPENi(MCFND,MNZ),
     &             KOPENf(MCFND,MNZ),FIXM(MCFND,MNZ),
     &             FIXT(MCFND,MNZ),FIXC(MCFND,MNZ),
     &             FIXK(MCFND,MNZ),FIXE(MCFND,MNZ),
     &             IWOPEN(MCFND,MNZ),ICFDCN(MCFND,MNZ),
     &             ICNACT(MCFND,MNZ),IVOLNOP(MCFND,MNZ)
      COMMON/Sbdary/NSB(MNZ),ISBi(MNSBZ,MNZ),ISBf(MNSBZ,MNZ),
     &              JSBi(MNSBZ,MNZ),JSBf(MNSBZ,MNZ),
     &              KSBi(MNSBZ,MNZ),KSBf(MNSBZ,MNZ),
     &              ISUFLC(MNSBZ,MNZ),IWSB(MNSBZ,MNZ),SSB(MNSBZ,MNZ),
     &              SSBHC(MNSBZ,MNZ),IVOLNSB(MNSBZ,MNZ),
     &              ITCtype(MNSBZ,MNZ),icTREF(MNSBZ,MNZ)
      common/blksso/NBLK(MNZ),INBLK(MNVLS,MNZ),NSSO(MNZ),
     &          INSSO(MNVLS,MNZ)

      dimension IDEL(MNVLS)

      character ITEM(MNVLS+8)*40
      character KEY*1,t24*24
      character VOLNAME*12,VCsurf*12,t12*12,d12*12
      character NAMES(MNVLS)*12
      character zngeo*3,etext*82
      CHARACTER*72 TITLE,CFTRFL,LPHI
      character outs*124
      logical origin,xgrid,ygrid,zgrid,zegrid,ortho,ok,dok,show

C      integer iglib   ! if 1 then X11, if 2 then GTK, if 3 then text only.
      integer INB     ! for radio button
      integer NITMS,INO,IAIN,IIND ! max items and current menu item
      real dum(1)

#ifdef OSI
      integer iside,isize,ifont
#else
      integer*8 iside,isize,ifont
#endif

C For help messages
      character helpinsub*24 ! subroutine name
      character helptopic*24 ! string (unique) for topic
      integer nbhelp     ! number of help lines found

      helpinsub='CFDBVAR'  ! set for subroutine

C << Note: for more complex domains grow this menu needs
C << to be one that can support multiple pages.

C << Changed to allow paging, example taken from vertex selection menu >>
C << EDVERT in edgeo.F                                                 >>

C C Initialise zone vertex menu size variables based on window size. 
C INO is the menu position, NITMS the current number of menu lines.
      MHEAD=2
      MCTL=5
      ILEN=NVOL(ICFD)
      IPACT=CREATE
      CALL EKPAGE(IPACT)

C Initial menu entry setup.
 5    IER=0
      ILEN=NVOL(ICFD)

C Create menu text.
      M=MHEAD
      ITEM(1)='1 Infer from geometry'
      ITEM(2)=' --------------------'
      do 10 I=1,ILEN
        IF(I.GE.IST.AND.(I.LE.(IST+MIFULL)))THEN
          M=M+1
          call EMKEY(I,KEY,IER)
          t24='  '
          if (IVTYPE(I,ICFD).eq.1) then
            t24='Solid  | Temp'
          elseif (IVTYPE(I,ICFD).eq.2) then
            t24='Solid  | Heat'
          elseif (IVTYPE(I,ICFD).eq.3) then
            t24='Solid  | Symmetrical'
          elseif (IVTYPE(I,ICFD).eq.4) then
            t24='Solid  | Conflated'
          elseif (IVTYPE(I,ICFD).eq.5) then
            t24='Solid  | Conflated'
          elseif (IVTYPE(I,ICFD).eq.6) then
            t24='Solid  | Conflated'
          elseif (IVTYPE(I,ICFD).eq.10) then
            t24='Opening| Pressure'
            IF (ICFMFS(ICFD).EQ.1)T24='Opening| Mass flow nwk'
          elseif (IVTYPE(I,ICFD).eq.11) then
            t24='Opening| Velocity'
            IF (ICFMFS(ICFD).EQ.1)T24='Opening| Mass flow nwk'
          elseif (IVTYPE(I,ICFD).eq.12) then
            t24='Opening| Zero gradient'
            IF (ICFMFS(ICFD).EQ.1)T24='Opening| Mass flow nwk'
          elseif (IVTYPE(I,ICFD).eq.13) then
            t24='Opening| Mass flow nwk'
            IF (ICFMFS(ICFD).EQ.1)T24='Opening| Mass flow nwk'
          elseif (IVTYPE(I,ICFD).eq.20) then
            t24='Source | Contaminant'
          elseif (IVTYPE(I,ICFD).eq.30) then
            t24='Blockage'
          else
            t24='UNKNOWN'
          endif
          write (ITEM(M),'(a,1x,3a)') KEY,VOLNAME(I,ICFD),':',t24
        endif
 10   continue

C Number of actual items displayed.
      NITMS=M+MCTL

C If a long list include page facility text.      
      IF(IPFLG.EQ.0)THEN  
        ITEM(M+1)=' --------------------'
      ELSE
        WRITE(ITEM(M+1),'(a,i2,a,i2,a)')'0 Page: ',IPM,' of ',MPM,
     &                                                      ' --------'
      ENDIF
      ITEM(M+2)='+ Add/ delete a boundary definition'
      ITEM(M+3)='! View current gridding'
      ITEM(M+4)='? Help'
      ITEM(M+5)='- exit'

      IF(MMOD.EQ.8)THEN

C Recalculate gridding in case this has changed.
        CALL INICNT
        CALL GRID(ier)
        call NEW2OLD

        call redraw(ier)
      ENDIF

C Help text for this menu.
      helptopic='cfd_boundary_variables'
      call gethelptext(helpinsub,helptopic,nbhelp)

C Display menu.
      INO=-1
      CALL EMENU('Boundary conditions',ITEM,NITMS,INO)

C Do things in reverse order here.
      if (INO.eq.NITMS) then

C Return
        return
      elseif (INO.eq.(NITMS-1)) then

C Help.
        helptopic='cfd_boundary_variables'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL PHELPD('boundary conditions',nbhelp,'-',0,0,IER)
      elseif (INO.eq.(NITMS-2)) then

C View the current gridding.
        if (origin.and.xgrid.AND.ygrid.AND.zgrid) then
          call cgd(izone,ierr)
        else
          CALL USRMSG(' Define origin of CFD domain ',
     &         ' and grid axes first ','W')
        endif

      elseif (INO.eq.(NITMS-3)) then

C Add a new/ delete a boundary definition.
        call easkabc('Add/ delete a boundary condition?','  ','Add',
     &               'Delete','Cancel',IAD,nbhelp)
        if (IAD.eq.1) then
          INB=1
          call EASKATOG('Which type of boundary?','  ',
     &       'Air flow opening','Solid surface','Source',
     &       'Blockage','Cancel',' ',' ',INB,nbhelp)

C Create boundary condition and call relevant editing routine.
          dum=(/0/)
          show=.false.
          t12='            '
          if (INB.eq.1) then

C Air flow opening.
            write (t12,'(a,i3.3)') 'Open',NVOL(ICFD)+1
            d12='  '
            call EASKS(t12,'Boundary condition name?',' ',12,d12,
     &        'BC name',IER,nbhelp)
            if (IER.eq.0) then
              call st2name(t12,d12)
            endif
            call MKCFDBC(show,11,0,0,0,0,0,0,0,d12,t12,dum,1,
     &        20.0,0.0,0.0,0.0,0.0,IER)
            call EDAIRO(NVOL(ICFD),IER)
            call NEW2OLD            
          elseif (INB.eq.2) then

C A solid boundary has been added.
            write (t12,'(a,i3.3)') 'Wall',NVOL(ICFD)+1
            d12='  '
            call EASKS(t12,'Boundary condition name?',' ',12,d12,
     &        'BC name',IER,nbhelp)
            if (IER.eq.0) then
              call st2name(t12,d12)
            endif
            call MKCFDBC(show,1,0,0,0,0,0,0,0,d12,t12,dum,1,
     &        20.0,0.0,0.0,0.0,0.0,IER)
            call EDSLDB(IZONE,NVOL(ICFD),IER)
            call NEW2OLD
            
          elseif (INB.eq.3) then

C A contaminate source has been added.
            write (t12,'(a,i3.3)') 'Source',NVOL(ICFD)+1
            d12='  '
            call EASKS(t12,'Boundary condition name?',' ',12,d12,
     &        'BC name',IER,nbhelp)
            if (IER.eq.0) then
              call st2name(t12,d12)
            endif
            call MKCFDBC(show,20,9,0,0,0,0,0,0,d12,t12,dum,1,
     &        0.0,0.0,0.0,0.0,0.0,IER)
            call EDSRC(NVOL(ICFD),IER)
            call NEW2OLD
            
          elseif (INB.eq.4) then

C A blockage has been added.
            write (t12,'(a,i3.3)') 'Block',NVOL(ICFD)+1
            d12='  '
            call EASKS(t12,'Boundary condition name?',' ',12,d12,
     &        'BC name',IER,nbhelp)
            if (IER.eq.0) then
              call st2name(t12,d12)
            endif
            call MKCFDBC(show,30,9,0,0,0,0,0,0,d12,t12,dum,1,
     &        0.0,0.0,0.0,0.0,0.0,IER)
            call EDBLK(NVOL(ICFD),IER)
            call NEW2OLD
          endif
        elseif (IAD.eq.2) then

C Delete a boundary condition. First, ask which one.
          if (NVOL(ICFD).gt.0) then
            INDEL=1
            do 233 I=1,NVOL(ICFD)
              NAMES(I)=VOLNAME(I,ICFD)
 233        continue
            INVOL=NVOL(ICFD)
            call EPICKS(INDEL,IDEL,' ',' ',12,INVOL,NAMES,
     &        'Key volume name',IER,nbhelp)
            CALL DELCFDBC(IDEL(1))
          endif
        endif

      elseif (INO.eq.(NITMS-4)) then
 
C If there are enough items allow paging control via EKPAGE.
        IF(IPFLG.EQ.1)THEN
          IPACT=EDIT
          CALL EKPAGE(IPACT)
        ENDIF
          
      elseif (INO.le.NVOL(ICFD).and.INO.gt.MHEAD) then

C Edit volume identified by KEYIND.
        CALL KEYIND(NITMS,INO,IAIN,IO)

C Existing volume, edit according to type.
        if (IVTYPE(IAIN,ICFD).ge.1.AND.IVTYPE(IAIN,ICFD).le.6) then
          call EDSLDB(IZONE,IAIN,IER)
        elseif (IVTYPE(IAIN,ICFD).ge.10.AND.IVTYPE(IAIN,ICFD).le.13)
     &                                                             then
          call EDAIRO(IAIN,IER)
        elseif (IVTYPE(IAIN,ICFD).eq.20) then
          call EDSRC(IAIN,IER)
        elseif (IVTYPE(IAIN,ICFD).eq.30) then
          call EDBLK(IAIN,IER)
        endif
      elseif (INO.eq.1) then

C Infer boundary conditions from geometry.
        CALL GENCFDBC(IER)
        if (IER.ne.0) then
          CALL EDISP(IUOUT,'Boundary conditions not generated.')
        endif

      endif
      goto 5

      end

C ************************** DELCFDBC **************************
C DELCFDBC deletes a CFD boundary condition, index ibc, from data
C structures.

      SUBROUTINE DELCFDBC(ibc)

#include "building.h"
#include "cfd.h"

      COMMON/ICFNOD/ICFD,ICP
      common/KEYVOLS/NVOL(MNZ),IVOLF(MNVLS,MNZ),IVCELLS(MNVLS,MNZ,2),
     &               JVCELLS(MNVLS,MNZ,2),KVCELLS(MNVLS,MNZ,2)
      common/KEYVOLN/VOLNAME(MNVLS,MNZ),VCsurf(MNVLS,MNZ)
      character VOLNAME*12,VCsurf*12
      common/KEYVDAT/IVTYPE(MNVLS,MNZ),VOLTemp(MNVLS,MNZ),
     &          VOLHeat(MNVLS,MNZ),IVConfl(MNVLS,MNZ),VOLHum(MNVLS,MNZ),
     &          VOLCO2(MNVLS,MNZ),VOLVel(MNVLS,MNZ),VOLDir(MNVLS,MNZ,2),
     &          VOLArea(MNVLS,MNZ),VOLPres(MNVLS,MNZ),
     &          VOLPol(MCTM,MNVLS,MNZ)

C Move all definitions greater than selected volume up one slot.
      do I=ibc,NVOL(ICFD)
        if (I.lt.NVOL(ICFD)) then
          IVOLF(I,ICFD)=IVOLF(I+1,ICFD)
          IVCELLS(I,ICFD,1)=IVCELLS(I+1,ICFD,1)
          IVCELLS(I,ICFD,2)=IVCELLS(I+1,ICFD,2)
          JVCELLS(I,ICFD,1)=JVCELLS(I+1,ICFD,1)
          JVCELLS(I,ICFD,2)=JVCELLS(I+1,ICFD,2)
          KVCELLS(I,ICFD,1)=KVCELLS(I+1,ICFD,1)
          KVCELLS(I,ICFD,2)=KVCELLS(I+1,ICFD,2)
          VOLNAME(I,ICFD)=VOLNAME(I+1,ICFD)
          VCsurf(I,ICFD)=VCsurf(I+1,ICFD)
          IVTYPE(I,ICFD)=IVTYPE(I+1,ICFD)
          VOLTemp(I,ICFD)=VOLTemp(I+1,ICFD)
          VOLHeat(I,ICFD)=VOLHeat(I+1,ICFD)
          IVConfl(I,ICFD)=IVConfl(I+1,ICFD)
          VOLHum(I,ICFD)=VOLHum(I+1,ICFD)
          VOLCO2(I,ICFD)=VOLCO2(I+1,ICFD)
          VOLVel(I,ICFD)=VOLVel(I+1,ICFD)
          VOLDir(I,ICFD,1)=VOLDir(I+1,ICFD,1)
          VOLDir(I,ICFD,2)=VOLDir(I+1,ICFD,2)
          VOLArea(I,ICFD)=VOLArea(I+1,ICFD)
          VOLPres(I,ICFD)=VOLPres(I+1,ICFD)
          do J=1,MCTM
            VOLPol(J,I,ICFD)=VOLPol(J,I,ICFD)
          enddo
        else
          IVOLF(I,ICFD)=0
          IVCELLS(I,ICFD,1)=0; IVCELLS(I,ICFD,2)=0
          JVCELLS(I,ICFD,1)=0; JVCELLS(I,ICFD,2)=0
          KVCELLS(I,ICFD,1)=0; KVCELLS(I,ICFD,2)=0
          VOLNAME(I,ICFD)=' '; VCsurf(I,ICFD)=' '
          IVTYPE(I,ICFD)=0
          VOLTemp(I,ICFD)=0.; VOLHeat(I,ICFD)=0.
          IVConfl(I,ICFD)=0
          VOLHum(I,ICFD)=0.; VOLCO2(I,ICFD)=0.
          VOLVel(I,ICFD)=0.
          VOLDir(I,ICFD,1)=0.; VOLDir(I,ICFD,2)=0.
          VOLArea(I,ICFD)=0.; VOLPres(I,ICFD)=0.
          do J=1,MCTM
            VOLPol(J,I,ICFD)=0.
          enddo
        endif
      enddo
      NVOL(ICFD)=NVOL(ICFD)-1
      CALL NEW2OLD

      RETURN
      END

C ************************** MKCFDBC **************************
C GENSGLBC creates a single CFD boundary condition. Input arguments are:

C itype (IVTYPE): BC type
C solid 
C 1 = temperature
C 2 = heat flux
C 3 = symmetrical
C 4 = temperature, conflated
C 5 = heat flux, conflated
C 6 = symmetrical, conflated
C air flow opening
C 10 = pressure
C 11 = velocity
C 12 = zero gradient
C 13 = mass
C source
C 20 = source
C blockage
C 30 = blockage

C ifac (IVOLF): face
C 1 = west
C 2 = east
C 3 = south
C 4 = north
C 5 = low
C 6 = high
C 7 = whole
C 8 = block
C 9 = source

C i1,i2,j1,j2,k1,k2 (IVCELLS,JVCELLS,KVCELLS):
C start and end cell indices for i,j, and k directions

C nam (VOLNAME): volume name
C confsurf (VCsurf): conflated surface name (solid only)
C pol(npol) (VOLPol): contaminant concentrations (source only)

C dat1,dat2,dat3,dat4,dat5: depending on itype, data for the BC
C itype=1:  VOLTemp,-      ,-      ,-      ,-
C itype=2:  VOLHeat,-      ,-      ,-      ,-
C itype=3:  -      ,-      ,-      ,-      ,-
C itype=4:  VOLTemp,IVConfl,-      ,-      ,-
C itype=5:  VOLHeat,IVConfl,-      ,-      ,-
C itype=6:  -      ,IVConfl,-      ,-      ,-
C itype=10: VOLTemp,VOLArea,VOLPres,-      ,-
C itype=11: VOLTemp,VOLArea,VOLVel ,VOLDir1,VOLDir2
C itype=12: VOLTemp,VOLArea,-      ,-      ,-
C itype=13: VOLTemp,VOLArea,VOLVel ,-      ,-
C itype=20: VOLHeat,VOLHum ,VOLCO2 ,-      ,-
C itype=30: VOLHeat,-      ,-      ,-      ,-

C If a blockage or source is created, this will first ensure the cells
C do not already contain a source or a blockage. If so, returns silently
C with IER=2.

C If show is true, will highlight the boundary condition in the 3D model
C display and print text feedback, then pause for half a second.

C Error codes:
C IER=1: unrecognised itype
C IER=2: cells already contain a blockage or source (silent)

      SUBROUTINE MKCFDBC(show,itype,ifac,i1,i2,j1,j2,k1,k2,nam,confsurf,
     &  pol,npol,dat1,dat2,dat3,dat4,dat5,IER)

#include "building.h"
#include "cfd.h"

      logical show
      integer itype,ifac,i1,i2,j1,j2,k1,k2,npol
      character nam*12,confsurf*12
      real pol(npol),dat1,dat2,dat3,dat4,dat5

      COMMON/OUTIN/IUOUT,IUIN      
      COMMON/ICFNOD/ICFD,ICP
      common/KEYVOLS/NVOL(MNZ),IVOLF(MNVLS,MNZ),IVCELLS(MNVLS,MNZ,2),
     &               JVCELLS(MNVLS,MNZ,2),KVCELLS(MNVLS,MNZ,2)
      common/KEYVOLN/VOLNAME(MNVLS,MNZ),VCsurf(MNVLS,MNZ)
      character VOLNAME*12,VCsurf*12
      common/KEYVDAT/IVTYPE(MNVLS,MNZ),VOLTemp(MNVLS,MNZ),
     &          VOLHeat(MNVLS,MNZ),IVConfl(MNVLS,MNZ),VOLHum(MNVLS,MNZ),
     &          VOLCO2(MNVLS,MNZ),VOLVel(MNVLS,MNZ),VOLDir(MNVLS,MNZ,2),
     &          VOLArea(MNVLS,MNZ),VOLPres(MNVLS,MNZ),
     &          VOLPol(MCTM,MNVLS,MNZ)
      COMMON/NDMAP/NOPEN(MNZ),MFNODE(MCFND,MNZ),IOPENi(MCFND,MNZ),
     &             IOPENf(MCFND,MNZ),JOPENi(MCFND,MNZ),
     &             JOPENf(MCFND,MNZ),KOPENi(MCFND,MNZ),
     &             KOPENf(MCFND,MNZ),FIXM(MCFND,MNZ),
     &             FIXT(MCFND,MNZ),FIXC(MCFND,MNZ),
     &             FIXK(MCFND,MNZ),FIXE(MCFND,MNZ),
     &             IWOPEN(MCFND,MNZ),ICFDCN(MCFND,MNZ),
     &             ICNACT(MCFND,MNZ),IVOLNOP(MCFND,MNZ)
      COMMON/Sbdary/NSB(MNZ),ISBi(MNSBZ,MNZ),ISBf(MNSBZ,MNZ),
     &              JSBi(MNSBZ,MNZ),JSBf(MNSBZ,MNZ),
     &              KSBi(MNSBZ,MNZ),KSBf(MNSBZ,MNZ),
     &              ISUFLC(MNSBZ,MNZ),IWSB(MNSBZ,MNZ),SSB(MNSBZ,MNZ),
     &              SSBHC(MNSBZ,MNZ),IVOLNSB(MNSBZ,MNZ),
     &              ITCtype(MNSBZ,MNZ),icTREF(MNSBZ,MNZ)
      common/blksso/NBLK(MNZ),INBLK(MNVLS,MNZ),NSSO(MNZ),
     &          INSSO(MNVLS,MNZ)
      common/KEYVblkd/IVblkd(NTCELX,NTCELY,NTCELZ)
      COMMON/CFDVIS/HAS_GEOM,ISHSB,ISHAO,IFACES,ISHBLK,ISHGEO,INITD
      logical HAS_GEOM
      character INITD*6

      character outs*124,stype*8,sfac*15
      logical xit

      if (itype.eq.1) then
        NVOL(ICFD)=NVOL(ICFD)+1
        n=NVOL(ICFD)
        IVTYPE(n,ICFD)=itype
        IVOLF(n,ICFD)=ifac
        IVCELLS(n,ICFD,1)=i1; IVCELLS(n,ICFD,2)=i2
        JVCELLS(n,ICFD,1)=j1; JVCELLS(n,ICFD,2)=j2
        KVCELLS(n,ICFD,1)=k1; KVCELLS(n,ICFD,2)=k2
        VOLNAME(n,ICFD)=nam
        VCsurf(n,ICFD)=confsurf
        VOLTemp(n,ICFD)=dat1
      elseif (itype.eq.2) then
        NVOL(ICFD)=NVOL(ICFD)+1
        n=NVOL(ICFD)
        IVTYPE(n,ICFD)=itype
        IVOLF(n,ICFD)=ifac
        IVCELLS(n,ICFD,1)=i1; IVCELLS(n,ICFD,2)=i2
        JVCELLS(n,ICFD,1)=j1; JVCELLS(n,ICFD,2)=j2
        KVCELLS(n,ICFD,1)=k1; KVCELLS(n,ICFD,2)=k2
        VOLNAME(n,ICFD)=nam
        VCsurf(n,ICFD)=confsurf
        VOLHeat(n,ICFD)=dat1
      elseif (itype.eq.3) then
        NVOL(ICFD)=NVOL(ICFD)+1
        n=NVOL(ICFD)
        IVTYPE(n,ICFD)=itype
        IVOLF(n,ICFD)=ifac
        IVCELLS(n,ICFD,1)=i1; IVCELLS(n,ICFD,2)=i2
        JVCELLS(n,ICFD,1)=j1; JVCELLS(n,ICFD,2)=j2
        KVCELLS(n,ICFD,1)=k1; KVCELLS(n,ICFD,2)=k2
        VOLNAME(n,ICFD)=nam
        VCsurf(n,ICFD)=confsurf
      elseif (itype.eq.4) then
        NVOL(ICFD)=NVOL(ICFD)+1
        n=NVOL(ICFD)
        IVTYPE(n,ICFD)=itype
        IVOLF(n,ICFD)=ifac
        IVCELLS(n,ICFD,1)=i1; IVCELLS(n,ICFD,2)=i2
        JVCELLS(n,ICFD,1)=j1; JVCELLS(n,ICFD,2)=j2
        KVCELLS(n,ICFD,1)=k1; KVCELLS(n,ICFD,2)=k2
        VOLNAME(n,ICFD)=nam
        VCsurf(n,ICFD)=confsurf
        VOLTemp(n,ICFD)=dat1
        IVconfl(n,ICFD)=int(dat2)
      elseif (itype.eq.5) then
        NVOL(ICFD)=NVOL(ICFD)+1
        n=NVOL(ICFD)
        IVTYPE(n,ICFD)=itype
        IVOLF(n,ICFD)=ifac
        IVCELLS(n,ICFD,1)=i1; IVCELLS(n,ICFD,2)=i2
        JVCELLS(n,ICFD,1)=j1; JVCELLS(n,ICFD,2)=j2
        KVCELLS(n,ICFD,1)=k1; KVCELLS(n,ICFD,2)=k2
        VOLNAME(n,ICFD)=nam
        VCsurf(n,ICFD)=confsurf
        VOLHeat(n,ICFD)=dat1
        IVconfl(n,ICFD)=int(dat2)
      elseif (itype.eq.6) then
        NVOL(ICFD)=NVOL(ICFD)+1
        n=NVOL(ICFD)
        IVTYPE(n,ICFD)=itype
        IVOLF(n,ICFD)=ifac
        IVCELLS(n,ICFD,1)=i1; IVCELLS(n,ICFD,2)=i2
        JVCELLS(n,ICFD,1)=j1; JVCELLS(n,ICFD,2)=j2
        KVCELLS(n,ICFD,1)=k1; KVCELLS(n,ICFD,2)=k2
        VOLNAME(n,ICFD)=nam
        VCsurf(n,ICFD)=confsurf
        IVconfl(n,ICFD)=int(dat2)
      elseif (itype.eq.10) then
        NVOL(ICFD)=NVOL(ICFD)+1
        n=NVOL(ICFD)
        IVTYPE(n,ICFD)=itype
        IVOLF(n,ICFD)=ifac
        IVCELLS(n,ICFD,1)=i1; IVCELLS(n,ICFD,2)=i2
        JVCELLS(n,ICFD,1)=j1; JVCELLS(n,ICFD,2)=j2
        KVCELLS(n,ICFD,1)=k1; KVCELLS(n,ICFD,2)=k2
        VOLNAME(n,ICFD)=nam
        VOLTemp(n,ICFD)=dat1
        VOLArea(n,ICFD)=dat2
        VOLPres(n,ICFD)=dat3
      elseif (itype.eq.11) then
        NVOL(ICFD)=NVOL(ICFD)+1
        n=NVOL(ICFD)
        IVTYPE(n,ICFD)=itype
        IVOLF(n,ICFD)=ifac
        IVCELLS(n,ICFD,1)=i1; IVCELLS(n,ICFD,2)=i2
        JVCELLS(n,ICFD,1)=j1; JVCELLS(n,ICFD,2)=j2
        KVCELLS(n,ICFD,1)=k1; KVCELLS(n,ICFD,2)=k2
        VOLNAME(n,ICFD)=nam
        VOLTemp(n,ICFD)=dat1
        VOLArea(n,ICFD)=dat2
        VOLVel(n,ICFD)=dat3
        VOLDir(n,ICFD,1)=dat4
        VOLDir(n,ICFD,2)=dat5
      elseif (itype.eq.12) then
        NVOL(ICFD)=NVOL(ICFD)+1
        n=NVOL(ICFD)
        IVTYPE(n,ICFD)=itype
        IVOLF(n,ICFD)=ifac
        IVCELLS(n,ICFD,1)=i1; IVCELLS(n,ICFD,2)=i2
        JVCELLS(n,ICFD,1)=j1; JVCELLS(n,ICFD,2)=j2
        KVCELLS(n,ICFD,1)=k1; KVCELLS(n,ICFD,2)=k2
        VOLNAME(n,ICFD)=nam
        VOLTemp(n,ICFD)=dat1
        VOLArea(n,ICFD)=dat2
      elseif (itype.eq.13) then
        NVOL(ICFD)=NVOL(ICFD)+1
        n=NVOL(ICFD)
        IVTYPE(n,ICFD)=itype
        IVOLF(n,ICFD)=ifac
        IVCELLS(n,ICFD,1)=i1; IVCELLS(n,ICFD,2)=i2
        JVCELLS(n,ICFD,1)=j1; JVCELLS(n,ICFD,2)=j2
        KVCELLS(n,ICFD,1)=k1; KVCELLS(n,ICFD,2)=k2
        VOLNAME(n,ICFD)=nam
        VOLTemp(n,ICFD)=dat1
        VOLArea(n,ICFD)=dat2
        VOLVel(n,ICFD)=dat3
      elseif (itype.eq.20) then

        NVOL(ICFD)=NVOL(ICFD)+1
        n=NVOL(ICFD)
        IVTYPE(n,ICFD)=itype
        IVOLF(n,ICFD)=ifac
        IVCELLS(n,ICFD,1)=i1; IVCELLS(n,ICFD,2)=i2
        JVCELLS(n,ICFD,1)=j1; JVCELLS(n,ICFD,2)=j2
        KVCELLS(n,ICFD,1)=k1; KVCELLS(n,ICFD,2)=k2
        VOLNAME(n,ICFD)=nam
        VOLPol(1:size(pol),n,ICFD)=pol
        VOLHeat(n,ICFD)=dat1
        VOLHum(n,ICFD)=dat2
        VOLCO2(n,ICFD)=dat3
      elseif (itype.eq.30) then

        NVOL(ICFD)=NVOL(ICFD)+1
        n=NVOL(ICFD)
        IVTYPE(n,ICFD)=itype
        IVOLF(n,ICFD)=ifac
        IVCELLS(n,ICFD,1)=i1; IVCELLS(n,ICFD,2)=i2
        JVCELLS(n,ICFD,1)=j1; JVCELLS(n,ICFD,2)=j2
        KVCELLS(n,ICFD,1)=k1; KVCELLS(n,ICFD,2)=k2
        VOLNAME(n,ICFD)=nam
        VOLHeat(n,ICFD)=dat1
        do i=i1,i2
          do j=j1,j2
            do k=k1,k2              
              IVblkd(i,j,k)=1
            enddo
          enddo
        enddo
      else
        IER=1
        write(outs,'(a,i2)')'MKCFDBC: unrecognised type ',itype
        call EDISP(IUOUT,outs)
        goto 999
      endif
      call NEW2OLD

      if (show) then
        if (itype.ge.1 .and. itype.le.6) then
          stype='solid   '
          lstype=5          
          ISHSB=NSB(ICFD); ISHAO=-1; ISHBLK=-1
        elseif (itype.ge.10 .and. itype.le.13) then
          stype='opening '
          lstype=7
          ISHSB=-1; ISHAO=NOPEN(ICFD); ISHBLK=-1
        elseif (itype.eq.20) then
          stype='source  '
          lstype=6
          ISHSB=-1; ISHAO=-1; ISHBLK=-1 ! << Need to implement sources
        elseif (itype.eq.30) then
          stype='blockage'
          lstype=8
          ISHSB=-1; ISHAO=-1; ISHBLK=NBLK(ICFD)
        endif

        if (ifac.eq.1) then 
          sfac=' on west face '
          lsfac=14
        elseif (ifac.eq.2) then
          sfac=' on east face '
          lsfac=14
        elseif (ifac.eq.3) then
          sfac=' on south face '
          lsfac=15
        elseif (ifac.eq.4) then
          sfac=' on north face '
          lsfac=15
        elseif (ifac.eq.5) then
          sfac=' on low face '
          lsfac=13
        elseif (ifac.eq.6) then
          sfac=' on high face '
          lsfac=14
        else
          sfac=' '
          lsfac=1
        endif

        write(outs,'(6a)')' Created ',stype(1:lstype),' BC',
     &    sfac(1:lsfac),' called ',nam(1:lnblnk(nam))
        call EDISP(IUOUT,outs)
        call redraw(IER)
        call pausems(500)
      endif

  999 RETURN
      END

C ************************** GENCFDBC **************************
C GENCFDBC attempts to infer CFD boundary conditions from zone geometry.
C IER=0 indicates no error.

      SUBROUTINE GENCFDBC(IER)

#include "building.h"
#include "model.h"
#include "geometry.h"
#include "cfd.h"

#ifdef OSI
      integer iicol
#else
      integer*8 iicol
#endif

      COMMON/FILEP/IFIL
      COMMON/OUTIN/IUOUT,IUIN      
      common/GRDVRTS/iorg(MNZ),ixend(MNZ),iyend(MNZ),izend(MNZ),
     &  izende(MNZ)
      COMMON/GEOM/XP(ntcelx),YP(ntcely),ZP(ntcelz),
     1            DXEP(ntcelx),DXPW(ntcelx),DYNP(ntcely),DYPS(ntcely),
     2            DZHP(ntcelz),DZPL(ntcelz),
     3            SEW(ntcelx),SNS(ntcely),SHL(ntcelz),
     4            XU(ntcelx),YV(ntcely),ZW(ntcelz)
      COMMON/ICFNOD/ICFD,ICP
      common/cfdconf/ICFBLD(MNZ),ICFMFS(MNZ)
      common/KEYVOLS/NVOL(MNZ),IVOLF(MNVLS,MNZ),IVCELLS(MNVLS,MNZ,2),
     &               JVCELLS(MNVLS,MNZ,2),KVCELLS(MNVLS,MNZ,2)
      common/KEYVOLN/VOLNAME(MNVLS,MNZ),VCsurf(MNVLS,MNZ)
      character VOLNAME*12,VCsurf*12
      common/KEYVDAT/IVTYPE(MNVLS,MNZ),VOLTemp(MNVLS,MNZ),
     &          VOLHeat(MNVLS,MNZ),IVConfl(MNVLS,MNZ),VOLHum(MNVLS,MNZ),
     &          VOLCO2(MNVLS,MNZ),VOLVel(MNVLS,MNZ),VOLDir(MNVLS,MNZ,2),
     &          VOLArea(MNVLS,MNZ),VOLPres(MNVLS,MNZ),
     &          VOLPol(MCTM,MNVLS,MNZ)
      COMMON/NDMAP/NOPEN(MNZ),MFNODE(MCFND,MNZ),IOPENi(MCFND,MNZ),
     &             IOPENf(MCFND,MNZ),JOPENi(MCFND,MNZ),
     &             JOPENf(MCFND,MNZ),KOPENi(MCFND,MNZ),
     &             KOPENf(MCFND,MNZ),FIXM(MCFND,MNZ),
     &             FIXT(MCFND,MNZ),FIXC(MCFND,MNZ),
     &             FIXK(MCFND,MNZ),FIXE(MCFND,MNZ),
     &             IWOPEN(MCFND,MNZ),ICFDCN(MCFND,MNZ),
     &             ICNACT(MCFND,MNZ),IVOLNOP(MCFND,MNZ)
      COMMON/Sbdary/NSB(MNZ),ISBi(MNSBZ,MNZ),ISBf(MNSBZ,MNZ),
     &              JSBi(MNSBZ,MNZ),JSBf(MNSBZ,MNZ),
     &              KSBi(MNSBZ,MNZ),KSBf(MNSBZ,MNZ),
     &              ISUFLC(MNSBZ,MNZ),IWSB(MNSBZ,MNZ),SSB(MNSBZ,MNZ),
     &              SSBHC(MNSBZ,MNZ),IVOLNSB(MNSBZ,MNZ),
     &              ITCtype(MNSBZ,MNZ),icTREF(MNSBZ,MNZ)
      COMMON/CFDVIS/HAS_GEOM,ISHSB,ISHAO,IFACES,ISHBLK,ISHGEO,INITD
      logical HAS_GEOM
      character INITD*6
      COMMON/C24/IZSTOCN(MCOM,MS)

      character outs*124,t6*6,t12*12,facstr*1
      logical ok,dok
      logical bumin,bumax,bvmin,bvmax,bwmin,bwmax,doafnasc
      real umin,umax,vmin,vmax,wmin,wmax,r1,r2,u,v,w
      real sumin,sumax,svmin,svmax,swmin,swmax
      integer snum

C For help messages
      character helpinsub*24 ! subroutine name
      character helptopic*24 ! string (unique) for topic
      integer nbhelp     ! number of help lines found

      helpinsub='GENCFDBC'  ! set for subroutine

C Remember visualisation settings.
      IISHSB=ISHSB; IISHAO=ISHAO; IISHBLK=ISHBLK
      IISHGEO=ISHGEO; IIFACES=IFACES

C Ensure parent-child relationships are defined.
      call GEOREAD(IFIL+1,LGEOM(ICP),ICP,1,IUOUT,IER)
      call suredgeadj(0,'c',ICP,IER)

C If there are already some defined, ask user if they want to clear
C existing boundary conditions.
      if (NVOL(ICFD).gt.0) then
        outs='Do you want to clear all existing boundary conditions?'
        dok=.true.
        ok=dok
        helptopic='cfdbc_infer_delete_exist'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL ASKOK(outs,' ',ok,dok,nbhelp)
        if (ok) then
          do I=1,NVOL(ICFD)
            IVOLF(I,ICFD)=0
            IVCELLS(I,ICFD,1)=0; IVCELLS(I,ICFD,2)=0
            JVCELLS(I,ICFD,1)=0; JVCELLS(I,ICFD,2)=0
            KVCELLS(I,ICFD,1)=0; KVCELLS(I,ICFD,2)=0
            VOLNAME(I,ICFD)=' '; VCsurf(I,ICFD)=' '
            IVTYPE(I,ICFD)=0
            VOLTemp(I,ICFD)=0.; VOLHeat(I,ICFD)=0.
            IVConfl(I,ICFD)=0
            VOLHum(I,ICFD)=0.; VOLCO2(I,ICFD)=0.
            VOLVel(I,ICFD)=0.
            VOLDir(I,ICFD,1)=0.; VOLDir(I,ICFD,2)=0.
            VOLArea(I,ICFD)=0.; VOLPres(I,ICFD)=0.
          enddo
          NVOL(ICFD)=0
        endif
      endif

C Make sure CFD - geometry transforms have been initialised.
      write(t6,'(2I3.3)')ICFD,ICP
      if (INITD.ne.t6) then
        if (HAS_GEOM) then
          CALL INIT_CFDTRANS(IER)
        else
          CALL INIT_CFDTRANS_NOGEOM(IER)
        endif
      endif
      if (IER.ne.0) then
        write(outs,'(2a)')
     &    'Error initialising coordinate transforms. ',
     &    'Check gridding.'
        CALL EDISP(IUOUT,outs)
        goto 999
      endif

C Grab the X, Y, and Z coordinates of the domain extremes (in CFD axes)
      umin=0.0; vmin=0.0; wmin=0.0
      CALL CFDTRANS(2,VCOORD(ICP,ixend(ICFD),1),
     &                VCOORD(ICP,ixend(ICFD),2),
     &                VCOORD(ICP,ixend(ICFD),3),umax,r1,r2,IER)
      CALL CFDTRANS(2,VCOORD(ICP,iyend(ICFD),1),
     &                VCOORD(ICP,iyend(ICFD),2),
     &                VCOORD(ICP,iyend(ICFD),3),r1,vmax,r2,IER)
      CALL CFDTRANS(2,VCOORD(ICP,izend(ICFD),1),
     &                VCOORD(ICP,izend(ICFD),2),
     &                VCOORD(ICP,izend(ICFD),3),r1,r2,wmax,IER)
      if (IER.ne.0) then
        write(outs,'(2a)')
     &    'Error transforming coordinates for domain axes. ',
     &    'Check gridding.'
        CALL EDISP(IUOUT,outs)
        goto 999
      endif

C Now, loop through surfaces in the zone and create solid and opening
C BCs.
      if (NZSUR(ICP).gt.0) then
        write(outs,'(a)')
     &    '*** Scanning surfaces for solid and opening BCs ***'
        CALL EDISP(IUOUT,outs)

        ok=.true.
        do ISUR=1,NZSUR(ICP)
          ICNN=IZSTOCN(ICP,ISUR)

C If this surface is a child, skip it; children are processed along with
C their parent.
          if (iparent(ICNN).ne.0) CYCLE

          ifac=0
          CALL GENEDGBC(isur,umax,vmax,wmax,ifac,ibctyp,isumin,isumax,
     &      isvmin,isvmax,iswmin,iswmax,doafnasc,ISTAT)

          if (ISTAT.eq.0) then
            CALL EDISP(IUOUT,'* No problems *')
          elseif (ISTAT.eq.4 .or. ISTAT.eq.6) then
            CALL EDISP(IUOUT,'* Surface skipped, but no problems *')            
          else
            CALL EDISP(IUOUT,'! Problem encountered !')
            ok=.false.
          endif
        enddo
        if (.not.ok) then
          CALL EDISP(IUOUT,
     &      'Warning: boundaries may not be fully defined.')
          CALL EDISP(IUOUT,
     &      'Please check text feedback above.')
        endif
      else        
        CALL EDISP(IUOUT,'*** No surfaces detected ***')
      endif

C Restore visualisation settings.
      ISHSB=IISHSB; ISHAO=IISHAO; ISHBLK=IISHBLK
      ISHGEO=IISHGEO; IFACES=IIFACES
      call redraw(IER)

  999 RETURN
      END

C ************ GENEDGBC
C GENEDGBC is a recursive sub used in the auto-generation of CFD BCs
C inferred from geometry. It is recursive because child surfaces can
C also be parents.

C Given a surface index isur, the extent of each axis ?max, and
C optionally the face index ifac, this subroutine:
C    finds the surface extents and, if ifac=0, determines what boundary
C    face the surface lies on,
C    checks the usage attributes to determine what boundary conditions
C    are needed (returned in ibctyp),
C    checks if the surface is rectangular and oriented with the CFD
C    domain axes,
C    finds the CFD grid cells that best match the extents of the surface
C    (returned in is[u/v/w][min/max]),
C    checks for child surfaces, runs this subroutine on them, and
C    confirms that we can handle the parent-child topology,
C    creates BC(s) for child surface(s) if applicable,
C    creates BC(s) for parent surface if applicable,

C ifac = 0: don't know yet, check
C ifac = 1: West
C ifac = 2: East
C ifac = 3: South
C ifac = 4: North
C ifac = 5: Low
C ifac = 6: High
C ifac = 7: not on a boundary face, surface normal in the u direction
C ifac = 8: not on a boundary face, surface normal in the v direction
C ifac = 9: not on a boundary face, surface normal in the w direction

C ibctyp = 1:  solid bc(s) for the whole surface
C ibctyp = 2:  single opening bc for the whole surface
C ibctyp = 3:  opening bc for bottom row of cells, solid bc for the rest
C ibctyp = 4:  opening bc for middle row of cells, solid bcs for the rest
C ibctyp = 5:  opening bc for 25% of height, at bottom edge, solid bc for
C              the rest
C ibctyp = 6:  opening bc for 25% of width, at vertical centre line,
C              solid bcs for rest
C ibctyp = 7:  opening bcs for 12.5% of height, at top and bottom edges,
C              solid bc for rest
C ibctyp = 8:  opening bc for 75% of height and width, at bottom edge,
C              solid bcs for rest
C ibctyp = 9:  for frames; opening bc as close to middle row of parent
C              and child as we can get
C ibctyp = 10: for frames; opening bc for 20% of child width, as close
C              to top edge as we can get

C ISTAT = 1: Fatal error, already reported.
C ISTAT = 2: Surface is not rectangular or aligned.
C ISTAT = 3: More than 1 child.
C ISTAT = 4: Surface is too narrow.
C ISTAT = 5: Problem with usage tags.
C ISTAT = 6: Surface is not on a boundary.

C doafnasc is currently unused, pending further development.

      RECURSIVE SUBROUTINE GENEDGBC(isur,umax,vmax,wmax,ifac,ibctyp,
     &  isumin,isumax,isvmin,isvmax,iswmin,iswmax,doafnasc,ISTAT)
#include "building.h"
#include "model.h"
#include "geometry.h"
#include "cfd.h"

      integer isur,ifac,ibctyp,isumin,isumax,isvmin,isvmax,iswmin,iswmax
      real umax,vmax,wmax,umin,vmin,wmin
      logical doafnasc

      COMMON/OUTIN/IUOUT,IUIN
      COMMON/ICFNOD/ICFD,ICP
      COMMON/CFDVIS/HAS_GEOM,ISHSB,ISHAO,IFACES,ISHBLK,ISHGEO,INITD
      logical HAS_GEOM
      character INITD*6
      COMMON/C3/IC1(MCON),IE1(MCON),ICT(MCON),IC2(MCON),IE2(MCON)
      COMMON/C24/IZSTOCN(MCOM,MS)
      common/KEYVOLS/NVOL(MNZ),IVOLF(MNVLS,MNZ),IVCELLS(MNVLS,MNZ,2),
     &               JVCELLS(MNVLS,MNZ,2),KVCELLS(MNVLS,MNZ,2)
      common/KEYVOLN/VOLNAME(MNVLS,MNZ),VCsurf(MNVLS,MNZ)
      character VOLNAME*12,VCsurf*12
      common/KEYVDAT/IVTYPE(MNVLS,MNZ),VOLTemp(MNVLS,MNZ),
     &          VOLHeat(MNVLS,MNZ),IVConfl(MNVLS,MNZ),VOLHum(MNVLS,MNZ),
     &          VOLCO2(MNVLS,MNZ),VOLVel(MNVLS,MNZ),VOLDir(MNVLS,MNZ,2),
     &          VOLArea(MNVLS,MNZ),VOLPres(MNVLS,MNZ),
     &          VOLPol(MCTM,MNVLS,MNZ)
      COMMON/NDMAP/NOPEN(MNZ),MFNODE(MCFND,MNZ),IOPENi(MCFND,MNZ),
     &             IOPENf(MCFND,MNZ),JOPENi(MCFND,MNZ),
     &             JOPENf(MCFND,MNZ),KOPENi(MCFND,MNZ),
     &             KOPENf(MCFND,MNZ),FIXM(MCFND,MNZ),
     &             FIXT(MCFND,MNZ),FIXC(MCFND,MNZ),
     &             FIXK(MCFND,MNZ),FIXE(MCFND,MNZ),
     &             IWOPEN(MCFND,MNZ),ICFDCN(MCFND,MNZ),
     &             ICNACT(MCFND,MNZ),IVOLNOP(MCFND,MNZ)
      COMMON/Sbdary/NSB(MNZ),ISBi(MNSBZ,MNZ),ISBf(MNSBZ,MNZ),
     &              JSBi(MNSBZ,MNZ),JSBf(MNSBZ,MNZ),
     &              KSBi(MNSBZ,MNZ),KSBf(MNSBZ,MNZ),
     &              ISUFLC(MNSBZ,MNZ),IWSB(MNSBZ,MNZ),SSB(MNSBZ,MNZ),
     &              SSBHC(MNSBZ,MNZ),IVOLNSB(MNSBZ,MNZ),
     &              ITCtype(MNSBZ,MNZ),icTREF(MNSBZ,MNZ)
      common/blksso/NBLK(MNZ),INBLK(MNVLS,MNZ),NSSO(MNZ),
     &          INSSO(MNVLS,MNZ)

      logical ok,dok,dochild(MCHILD),dothisur,lumped,show,isrect      
      logical bumin,bumax,bvmin,bvmax,bwmin,bwmax
      character outs*124,t12*12,facstr*1,snam*12
      real u,v,w,dum(1),sumin,sumax,svmin,svmax,swmin,swmax
      integer iparext(4),iparca,icldext(MCHILD,4),icldtotca,isclose(4),
     &        icldbctyp(MCHILD),icldverts(MV)

C For help messages
      character helpinsub*24 ! subroutine name
      character helptopic*24 ! string (unique) for topic
      integer nbhelp     ! number of help lines found

      helpinsub='GENCFDBC'  ! set for calling subroutine

      ISTAT=0
      dum=(/0/)
      show=.true.
      icnn=IZSTOCN(ICP,isur)

C Check if this surface is a child. Unless it is also a parent, it will
C be created by the invocation of this subroutine for its parent.
      if (IPARENT(icnn).gt.0) then
        if (NBCHILD(icnn).gt.0) then
          dothisur=.true.
        else          
          dothisur=.false.
        endif
      else
        dothisur=.true.
      endif

      snam=SSNAME(ICNN) ! surface name

C Establish surface extents, and if not already known determine what
C face the surface lies on if any (south, west, low, etc).
C b?min/max  = is on the face matching the min/max extent of the domain
C              in the ? direction (logical).
C s?min/max  = min/max extent of the surface in the ? direction in CFD
C              domain coords (real).
      bumin=.TRUE.; bumax=.TRUE.
      bvmin=.TRUE.; bvmax=.TRUE.
      bwmin=.TRUE.; bwmax=.TRUE.
      iminmax=6
      umin=0.0; vmin=0.0; wmin=0.0
      sumin=umax; sumax=umin
      svmin=vmax; svmax=vmin
      swmin=wmax; swmax=wmin
      IER=0
      do I=1,NZNVER(ICNN)
        IVRT=NZJVN(ICNN,I)
        CALL CFDTRANS(2,VCOORD(ICP,IVRT,1),VCOORD(ICP,IVRT,2),
     &    VCOORD(ICP,IVRT,3),u,v,w,IER)
        if (IER.ne.0) then
          write(outs,'(2a,i3.3,1x,2a)')'Error transforming ',
     &      'coordinates for surface ',ISUR,snam(1:lnblnk(snam)),
     &      '. Check gridding.'
          CALL EDISP(IUOUT,outs)
          iminmax=-1
          EXIT
        endif
        if (u.lt.sumin) sumin=u
        if (u.gt.sumax) sumax=u
        if (v.lt.svmin) svmin=v
        if (v.gt.svmax) svmax=v
        if (w.lt.swmin) swmin=w
        if (w.gt.swmax) swmax=w

        if (ifac.eq.0) then
          if (bumin) then
            call eclose(umin,u,0.01,ok)
            if (.not.ok) then
              bumin=.FALSE.
              iminmax=iminmax-1
            endif
          endif
          if (bumax) then
            call eclose(umax,u,0.01,ok)
            if (.not.ok) then
              bumax=.FALSE.
              iminmax=iminmax-1
            endif
          endif
          if (bvmin) then
            call eclose(vmin,v,0.01,ok)
            if (.not.ok) then
              bvmin=.FALSE.
              iminmax=iminmax-1
            endif
          endif
          if (bvmax) then
            call eclose(vmax,v,0.01,ok)
            if (.not.ok) then
              bvmax=.FALSE.
              iminmax=iminmax-1
            endif
          endif
          if (bwmin) then
            call eclose(wmin,w,0.01,ok)
            if (.not.ok) then
              bwmin=.FALSE.
              iminmax=iminmax-1
            endif
          endif
          if (bwmax) then
            call eclose(wmax,w,0.01,ok)
            if (.not.ok) then
              bwmax=.FALSE.
              iminmax=iminmax-1
            endif
          endif
        endif
      enddo

      if (ifac.eq.0) then
        if (iminmax.eq.0) then

C Surface is not on a boundary face. Next surface.
          write(outs,'(a,i3.3,1x,2a)')'Surface ',ISUR,
     &      snam(1:lnblnk(snam)),
     &      ' is not on a boundary face, skipping.'
          CALL EDISP(IUOUT,outs)
          ISTAT=6
          goto 999
        elseif (iminmax.eq.1) then

C Surface is on a boundary. Continue.           
          write(outs,'(a,i3.3,1x,2a)')'Surface ',ISUR,
     &     snam(1:lnblnk(snam)),' is on a boundary face, processing ...'
          CALL EDISP(IUOUT,outs)
          if (bumin) then ! west
            ifac=1
          elseif (bumax) then ! east
            ifac=2
          elseif (bvmin) then ! south
            ifac=3
          elseif (bvmax) then ! north
            ifac=4
          elseif (bwmin) then ! low
            ifac=5
          elseif (bwmax) then ! high
            ifac=6
          endif
        elseif (iminmax.gt.1) then

C More than one boundary face. Notify the user, then skip to the next
C surface.
          write(outs,'(a,i3.3,1x,3a)')'Surface ',ISUR,
     &      snam(1:lnblnk(snam)),
     &    ' seems to be on more than one boundary face. Check gridding.'
          CALL EDISP(IUOUT,outs)
          ISTAT=1
          goto 999
        else

C Something has gone wrong, presumably already reported. Skip this
C surface.
          ISTAT=1
          goto 999
        endif
      endif

C Now make sure the surface is rectangular and oriented parallel to the
C CFD domain axes.
C If the surface has more than 4 vertices, we ensure it is rectangular
C by checking all vertices to see if they share at least one coordinate
C in the plane of the boundary face with the surface extents.
C We ensure it is oriented correctly by checking that at least 4 of the
C vertices share both coordinates in the plane of the boundary face with
C the surface extents.
      if (NZNVER(icnn).gt.4) then

C First, check if this surface is a parent. If it is, we first need to
C establish that the child is rectangular - wait until later.
        if (NBCHILD(icnn).gt.0) then
          isrect=.false.

        else
          icorners=0
          do I=1,NZNVER(ICNN)
            IVRT=NZJVN(ICNN,I)
            CALL CFDTRANS(2,VCOORD(ICP,IVRT,1),VCOORD(ICP,IVRT,2),
     &        VCOORD(ICP,IVRT,3),u,v,w,IER)
            if (ifac.eq.1 .or. ifac.eq.2) then
              call ECLOSE(v,svmin,0.0001,ok)
              if (.not.ok) call ECLOSE(v,svmax,0.0001,ok)
              call ECLOSE(w,swmin,0.0001,dok)
              if (.not.dok) call ECLOSE(w,swmax,0.0001,dok)
              if (ok.and.dok) icorners=icorners+1
              if (ok.or.dok) CYCLE
            elseif (ifac.eq.3 .or. ifac.eq.4) then
              call ECLOSE(u,sumin,0.0001,ok)
              if (.not.ok) call ECLOSE(u,sumax,0.0001,ok)
              call ECLOSE(w,swmin,0.0001,dok)
              if (.not.dok) call ECLOSE(w,swmax,0.0001,dok)
              if (ok.and.dok) icorners=icorners+1
              if (ok.or.dok) CYCLE
            elseif (ifac.eq.5 .or. ifac.eq.6) then
              call ECLOSE(u,sumin,0.0001,ok)
              if (.not.ok) call ECLOSE(u,sumax,0.0001,ok)
              call ECLOSE(v,svmin,0.0001,dok)
              if (.not.dok) call ECLOSE(v,svmax,0.0001,dok)
              if (ok.and.dok) icorners=icorners+1
              if (ok.or.dok) CYCLE
            endif
            ISTAT=2
            EXIT
          enddo
          if (ISTAT.eq.2) then
            CALL EDISP(IUOUT,
     &        ' This surface is not rectangular, skipping.')
            goto 999
          endif
          if (icorners.lt.4) then
            ISTAT=2
            CALL EDISP(IUOUT,
     &        ' This surface is not aligned with CFD axes, skipping.')
            goto 999
          endif
          isrect=.true.
        endif

C If the surface has less than 4 vertices, it definitely isn't
C rectanglar.
      elseif (NZNVER(icnn).lt.4) then
        ISTAT=2
        CALL EDISP(IUOUT,
     &    ' This surface is not rectangular, skipping.')
        goto 999

C If it has only 4 vertices, ensure that each one shares both
C coordinates in the plane of the boundary face with the surface
C extents.
      else
        do I=1,NZNVER(ICNN)
          IVRT=NZJVN(ICNN,I)
          CALL CFDTRANS(2,VCOORD(ICP,IVRT,1),VCOORD(ICP,IVRT,2),
     &      VCOORD(ICP,IVRT,3),u,v,w,IER)
          if (ifac.eq.1 .or. ifac.eq.2) then
            call ECLOSE(v,svmin,0.0001,ok)
            if (.not.ok) call ECLOSE(v,svmax,0.0001,ok)
            call ECLOSE(w,swmin,0.0001,dok)
            if (.not.dok) call ECLOSE(w,swmax,0.0001,dok)
            if (ok.and.dok) CYCLE
          elseif (ifac.eq.3 .or. ifac.eq.4) then
            call ECLOSE(u,sumin,0.0001,ok)
            if (.not.ok) call ECLOSE(u,sumax,0.0001,ok)
            call ECLOSE(w,swmin,0.0001,dok)
            if (.not.dok) call ECLOSE(w,swmax,0.0001,dok)
            if (ok.and.dok) CYCLE
          elseif (ifac.eq.5 .or. ifac.eq.6) then
            call ECLOSE(u,sumin,0.0001,ok)
            if (.not.ok) call ECLOSE(u,sumax,0.0001,ok)
            call ECLOSE(v,svmin,0.0001,dok)
            if (.not.dok) call ECLOSE(v,svmax,0.0001,dok)
            if (ok.and.dok) CYCLE
          endif
          ISTAT=2
          EXIT
        enddo
        if (ISTAT.eq.2) then
          CALL EDISP(IUOUT,
     &     ' This surface is not rectangular and/or aligned, skipping.')
          goto 999
        endif
        isrect=.true.
      endif

C Initialise face variables.
      if (ifac.eq.1) then
        IFACES=3
        facstr='W'
      elseif (ifac.eq.2) then
        IFACES=3
        facstr='E'
      elseif (ifac.eq.3) then
        IFACES=2
        facstr='S'
      elseif (ifac.eq.4) then
        IFACES=2
        facstr='N'
      elseif (ifac.eq.5) then
        IFACES=4
        facstr='L'
      elseif (ifac.eq.6) then
        IFACES=4
        facstr='H'
      endif

C Sniff around usage flags to determine what this surface is. This code
C is designed to be as consistent as possible with flow network
C autogeneration (scan_flow_use in mfprb3.F). If the usage flags in
C geometry files are updated, this code shouild be updated accordingly.
      ibctyp=0
      t12=SSUSE(icnn,1)
      if (t12(1:4).eq.'WALL') then
        ibctyp=1

      elseif (t12(1:4).eq.'DOOR' .or. t12(3:6).eq.'DOOR') then
        t12=SSUSE(ICNN,2)
        if (t12(1:6).eq.'CLOSED') then
          ibctyp=1
        elseif (t12(1:8).eq.'UNDERCUT') then
          ibctyp=3
        elseif (t12(1:4).eq.'OPEN') then
          ibctyp=2
        elseif (t12(1:5).eq.'BIDIR') then
          ibctyp=2
        else
          goto 51
        endif

      elseif (t12(1:5).eq.'FRAME' .or. t12(3:7).eq.'FRAME') then
        t12=SSUSE(ICNN,2)
        if (t12(1:6).eq.'CLOSED') then
          ibctyp=1
        elseif (t12(1:5).eq.'CRACK') then
          ibctyp=9
        elseif (t12(1:5).eq.'VENT') then
          ibctyp=10
        else
          goto 51
        endif        

      elseif (t12(1:6).eq.'WINDOW' .or. t12(3:8).eq.'WINDOW') then
        t12=SSUSE(ICNN,2)
        if (t12(1:6).eq.'CLOSED') then
          ibctyp=1
        elseif (t12(1:5).eq.'CRACK') then
          ibctyp=4
        elseif (t12(1:4).eq.'OPEN') then
          CALL ZSURLEHI(ICP,ISUR,XYDIS,ZDIS,i1,i2,i3,i4)
          if (XYDIS/ZDIS.gt.1.0) then
            ibctyp=6
          else
            ibctyp=5
          endif
        elseif (t12(1:4).eq.'SASH') then
          ibctyp=7
        elseif (t12(1:5).eq.'BIDIR') then
          ibctyp=2
        else
          goto 51
        endif

      elseif (t12(1:5).eq.'GRILL') then
        t12=SSUSE(ICNN,2)
        if (t12(1:6).eq.'CLOSED') then
          ibctyp=1
        elseif (t12(1:5).eq.'CRACK') then
          ibctyp=4
        elseif (t12(1:5).eq.'INLET') then
          ibctyp=2
        elseif (t12(1:7).eq.'EXTRACT') then
          ibctyp=2
        elseif (t12(1:4).eq.'OPEN') then
          ibctyp=2
        else
          goto 51
        endif

      elseif (t12(1:4).eq.'ROOF') then
        ibctyp=1

      elseif (t12(1:4).eq.'FICT') then
        t12=SSUSE(ICNN,2)
        if (t12(1:6).eq.'CLOSED') then
          ibctyp=1
        elseif (t12(1:5).eq.'CRACK') then
          ibctyp=3
        elseif (t12(1:4).eq.'OPEN') then
          CALL ZSURLEHI(ICP,ISUR,XYDIS,ZDIS,i1,i2,i3,i4)
          if (XYDIS/ZDIS.gt.1.0) then
            ibctyp=6
          else
            ibctyp=5
          endif
        elseif (t12(1:5).eq.'BIDIR') then
          ibctyp=8
        else
          goto 51
        endif

      elseif (t12(1:5).eq.'FLOOR') then
        ibctyp=1

      elseif (t12(1:5).eq.'FURNI') then
        ibctyp=1

      elseif (t12(1:5).eq.'PARTN') then
        ibctyp=1

      elseif (t12(1:5).eq.'STRUC') then
C << Extend to scan for these not on boundaries, and create
C blockages/sources >>
        ibctyp=1

      elseif (t12(1:7).eq.'ITEQUIP') then
C << Extend to scan for these not on boundaries, and create
C blockages/sources >>
        ibctyp=1

      elseif (t12(1:7).eq.'FIXTURE') then
C << Extend to scan for these not on boundaries, and create
C blockages/sources >>
        ibctyp=1

      elseif (t12(1:6).eq.'PLANTS') then
C << Extend to scan for these not on boundaries, and create
C blockages/sources >>
        ibctyp=1

      elseif (t12(1:1).eq.'-') then
C Assume that a surface with no usage tag is a solid surface.
        ibctyp=1
        write(outs,'(2a)')' This surface has no usage tag, assuming a',
     &    ' solid surface.'
        CALL EDISP(IUOUT,outs)
      else
        goto 51
      endif

C If usage tag is not recognised, notify user and skip to the next
C surface.
      goto 52
   51 write(outs,'(3a)')' The usage tag ',t12(1:lnblnk(t12)),
     &  ' is not recognised, skipping this surface.'
      CALL EDISP(IUOUT,outs)
      ISTAT=5
      goto 999

C Find grid cells to best match the extents of the surface. Note that
C this will find indices on the velocity grid. This is negatively offset
C from the scalar grid, so for example XU(1) is the negative U edge of
C the scalar cell XP(1), and the positive U edge is XU(2).
C So the actual (scalar) grid cells to use are is?min and is?max-1.
   52 CALL FDCFDPT(2,sumin,svmin,swmin,isumin,isvmin,iswmin,IER)
      CALL FDCFDPT(2,sumax,svmax,swmax,isumax,isvmax,iswmax,IER)
      if (IER.ne.0) then
        write(outs,'(a)')
     &    ' Error finding cells to match geometry. Check gridding.'
        CALL EDISP(IUOUT,outs)
        ISTAT=1
        goto 999
      endif
      if (ifac.eq.1 .or. ifac.eq.2) then
        if (ifac.eq.1) then
          isumax=isumin+1
        else
          isumin=isumax-1
        endif
        if (isvmin.eq.isvmax) then
          write(outs,'(2a)')' This surface is too narrow in the Y',
     &      ' direction, no BC will be created.'
          CALL EDISP(IUOUT,outs)
          ISTAT=4
          goto 999
        elseif (iswmin.eq.iswmax) then
          write(outs,'(2a)')' This surface is too narrow in the Z',
     &      ' direction, no BC will be created.'
          CALL EDISP(IUOUT,outs)
          ISTAT=4
          goto 999
        endif
        iparext(1:4)=(/isvmin,isvmax,iswmin,iswmax/)
        iparca=(isvmax-isvmin)*(iswmax-iswmin)
      elseif (ifac.eq.3 .or. ifac.eq.4) then
        if (ifac.eq.3) then
          isvmax=isvmin+1
        else
          isvmin=isvmax-1
        endif
        if (isumin.eq.isumax) then
          write(outs,'(2a)')' This surface is too narrow in the X',
     &      ' direction, no BC will be created.'
          CALL EDISP(IUOUT,outs)
          ISTAT=4
          goto 999
        elseif (iswmin.eq.iswmax) then
          write(outs,'(2a)')' This surface is too narrow in the Z',
     &      ' direction, no BC will be created.'
          CALL EDISP(IUOUT,outs)
          ISTAT=4
          goto 999
        endif
        iparext(1:4)=(/isumin,isumax,iswmin,iswmax/)
        iparca=(isumax-isumin)*(iswmax-iswmin)
      elseif (ifac.eq.5 .or. ifac.eq.6) then
        if (ifac.eq.5) then
          iswmax=iswmin+1
        else
          iswmin=iswmax-1
        endif
        if (isumin.eq.isumax) then
          write(outs,'(2a)')' This surface is too narrow in the X',
     &      ' direction, no BC will be created.'
          CALL EDISP(IUOUT,outs)
          ISTAT=4
          goto 999
        elseif (isvmin.eq.isvmax) then
          write(outs,'(2a)')' This surface is too narrow in the Y',
     &      ' direction, no BC will be created.'
          CALL EDISP(IUOUT,outs)
          ISTAT=4
          goto 999
        endif
        iparext(1:4)=(/isumin,isumax,isvmin,isvmax/)
        iparca=(isumax-isumin)*(isvmax-isvmin)
      endif

C Now, check for children.
      if (NBCHILD(icnn).gt.0) then
        write(outs,'(a)')' Child surfaces detected.'
        CALL EDISP(IUOUT,outs)
        icldtotca=0
        icldvert=0
        ncldverts=0
        do icld=1,NBCHILD(icnn)
          icldcnn=ICHILD(icnn,icld)
          icldsur=IE1(icldcnn)
          t12=SSNAME(icldcnn)
          write(outs,'(a,i3.3,1x,2a)')' Processing child surface ',
     &      icldsur,t12(1:lnblnk(t12)),' ...'
          CALL EDISP(IUOUT,outs)

C Call this subroutine for the child.
          CALL GENEDGBC(icldsur,umax,vmax,wmax,ifac,iibctyp,iisumin,
     &      iisumax,iisvmin,iisvmax,iiswmin,iiswmax,doafnasc,ISTAT)

C Check for problems with the child.
          if (ISTAT.eq.1) then
            write(outs,'(a)')
     &        ' Error with child, parent will also be skipped.'
            CALL EDISP(IUOUT,outs)
            goto 999
          endif
          if (ISTAT.eq.2) then
            write(outs,'(a)')
     &        ' Child is not rectangular, parent will also be skipped.'
            CALL EDISP(IUOUT,outs)
            goto 999
          endif
          if (ISTAT.ge.3 .and. ISTAT.le.5) then ! can't do child, but can still do parent
            dochild(icld)=.false.
            CYCLE
          else
            dochild(icld)=.true.
          endif
          icldbctyp(icld)=iibctyp

C Store the extents of the child (in cells), and the total cells covered
C by all children.
          if (ifac.eq.1) then
            iisumax=iisumin+1
            icldext(icld,1:4)=(/iisvmin,iisvmax,iiswmin,iiswmax/)
            icldtotca=icldtotca+(iisvmax-iisvmin)*(iiswmax-iiswmin)
          elseif (ifac.eq.2) then
            iisumax=iisumin+1
            icldext(icld,1:4)=(/iisvmin,iisvmax,iiswmin,iiswmax/)
            icldtotca=icldtotca+(iisvmax-iisvmin)*(iiswmax-iiswmin)
          elseif (ifac.eq.3) then
            iisvmax=iisvmin+1
            icldext(icld,1:4)=(/iisumin,iisumax,iiswmin,iiswmax/)
            icldtotca=icldtotca+(iisumax-iisumin)*(iiswmax-iiswmin)
          elseif (ifac.eq.4) then
            iisvmax=iisvmin+1
            icldext(icld,1:4)=(/iisumin,iisumax,iiswmin,iiswmax/)
            icldtotca=icldtotca+(iisumax-iisumin)*(iiswmax-iiswmin)
          elseif (ifac.eq.5) then
            iiswmax=iiswmin+1
            icldext(icld,1:4)=(/iisumin,iisumax,iisvmin,iisvmax/)
            icldtotca=icldtotca+(iisumax-iisumin)*(iisvmax-iisvmin)
          elseif (ifac.eq.6) then
            iiswmax=iiswmin+1
            icldext(icld,1:4)=(/iisumin,iisumax,iisvmin,iisvmax/)
            icldtotca=icldtotca+(iisumax-iisumin)*(iisvmax-iisvmin)
          endif

C Assemble an array of vertices in child surfaces, for checking
C rectangularity.
          do i=1,NZNVER(icldcnn)
            ivrt=NZJVN(icldcnn,i)            
            ncldverts=ncldverts+1
            icldverts(ncldverts)=ivrt
          enddo
        enddo

C We now have the extents of all children, and have confirmed that they
C are all rectangular. We now need to ensure that the parent is also
C rectangular. To do this, use a procedure similar to childless surfaces
C with more than 4 vertices, with the exception that any vertex that is
C in one of the children must be ommited from the check.
        if (.not.isrect) then
          icorners=0
          do i=1,NZNVER(icnn)
            ivrt=NZJVN(icnn,i)
            ok=.false.
            do j=1,ncldverts
              if (ivrt.eq.icldverts(j)) then
                ok=.true.
                EXIT
              endif
            enddo
            if (ok) CYCLE
            CALL CFDTRANS(2,VCOORD(ICP,IVRT,1),VCOORD(ICP,IVRT,2),
     &        VCOORD(ICP,IVRT,3),u,v,w,IER)
            if (ifac.eq.1 .or. ifac.eq.2) then
              call ECLOSE(v,svmin,0.0001,ok)
              if (.not.ok) call ECLOSE(v,svmax,0.0001,ok)
              call ECLOSE(w,swmin,0.0001,dok)
              if (.not.dok) call ECLOSE(w,swmax,0.0001,dok)
              if (ok.and.dok) icorners=icorners+1
              if (ok.or.dok) CYCLE
            elseif (ifac.eq.3 .or. ifac.eq.4) then
              call ECLOSE(u,sumin,0.0001,ok)
              if (.not.ok) call ECLOSE(u,sumax,0.0001,ok)
              call ECLOSE(w,swmin,0.0001,dok)
              if (.not.dok) call ECLOSE(w,swmax,0.0001,dok)
              if (ok.and.dok) icorners=icorners+1
              if (ok.or.dok) CYCLE
            elseif (ifac.eq.5 .or. ifac.eq.6) then
              call ECLOSE(u,sumin,0.0001,ok)
              if (.not.ok) call ECLOSE(u,sumax,0.0001,ok)
              call ECLOSE(v,svmin,0.0001,dok)
              if (.not.dok) call ECLOSE(v,svmax,0.0001,dok)
              if (ok.and.dok) icorners=icorners+1
              if (ok.or.dok) CYCLE
            endif
            ISTAT=2
            EXIT
          enddo
          if (ISTAT.eq.2) then
            CALL EDISP(IUOUT,
     &        ' This surface is not rectangular, skipping.')
            goto 999
          endif
          if (icorners.lt.4) then
            ISTAT=2
            CALL EDISP(IUOUT,
     &        ' This surface is not aligned with CFD axes, skipping.')
            goto 999
          endif
          isrect=.true.
        endif

C Now, we try to establish if the parent is worth representing on its
C own, or should be lumped in with the child(ren).

C Rules are (within the discretized CFD grid):
C 1: If the whole extent of the parent is covered by its child(ren),
C    then the parent will not be represented with a solid BC.
C 2: Otherwise, provided there is a single child, if any edge of the
C    child is more than 1 cell away from its neighbouring parent edge,
C    the parent will be represented with a solid BC.
C 3: Otherwise, highlight the surface and ask the user.
        if (icldtotca.ge.iparca) then ! all cells of parent taken up by children
          write(outs,'(2a)')
     &  ' Parent framing is narrow; no solid BC will be created for it.'
          CALL EDISP(IUOUT,outs)
          lumped=.true.
        else          
          if (NBCHILD(icnn).eq.1) then

C Check proximity of child edges to parent edges.
            if (icldext(1,1).eq.iparext(1)) then
              isclose(1)=1
            elseif (icldext(1,1)-1.eq.iparext(1)) then
              isclose(1)=2
            else
              isclose(1)=0
            endif
            if (icldext(1,2).eq.iparext(2)) then
              isclose(2)=1
            elseif (icldext(1,2)+1.eq.iparext(2)) then
              isclose(2)=2
            else
              isclose(2)=0
            endif
            if (icldext(1,3).eq.iparext(3)) then
              isclose(3)=1
            elseif (icldext(1,3)-1.eq.iparext(3)) then
              isclose(3)=2
            else
              isclose(3)=0
            endif
            if (icldext(1,4).eq.iparext(4)) then
              isclose(4)=1
            elseif (icldext(1,4)+1.eq.iparext(4)) then
              isclose(4)=2
            else
              isclose(4)=0
            endif

C If all edges of the child are close to the edges of the present
C surface, ask user if they want to lump the present surface in with its
C child. Do not ask if ...
            if (isclose(1).gt.0 .and. isclose(2).gt.0 .and.
     &          isclose(3).gt.0 .and. isclose(4).gt.0) then

C ... the present surface is a frame and we need space for an opening
C BC.
              if ((ibctyp.eq.9 .or. ibctyp.eq.10) .or.

C ... the child is also a parent (in which case the BCs would have been
C created in the invocation of GENEDGBC for that surface).
     &            NBCHILD(icldcnn).gt.0) then
                lumped=.false.
                write(outs,'(3a)')' Parent framing is narrow; BCs ',
     &            'created for it may not well represent surface ',
     &            'position and area.'
                call EDISP(IUOUT,outs)
                write(outs,'(a)')
     &            ' Please check after autogeneration is complete.'
                call EDISP(IUOUT,outs)
                goto 70
              endif

              ISHSB=0; ISHAO=0; ISHBLK=0; ISHGEO=-ISUR
              call redraw(IER)
              IR=2
              helptopic='cfdbc_parentframing'
              call gethelptext(helpinsub,helptopic,nbhelp)
   61         CALL EASKABC(
     &         'Do you want to represent parent surface with solid BCs?'
     &         ,'(see help for details)','yes','no','preview',IR,nbhelp)
              ok=.false.
              if (IR.eq.1) then
                lumped=.false.
              elseif (IR.eq.2) then
                lumped=.true.

C Extend the extents of the child to those of the parent.
                icldext(1,1)=iparext(1)
                icldext(1,2)=iparext(2)
                icldext(1,3)=iparext(3)
                icldext(1,4)=iparext(4)
              elseif (IR.eq.3) then ! preview

C First, show the parent surface highlighted.
                ISHSB=0; ISHAO=0; ISHBLK=0; ISHGEO=-isur
                CALL redraw(IER)
                CALL continuebox('Parent surface.',' ','ok')

C Create the boundary conditions that would be created for the parent.
                icreated=0
                if (isclose(1).eq.2) then
                  if (ifac.eq.1 .or. ifac.eq.2) then
                    i1=isumin; i2=isumax-1
                    j1=isvmin; j2=isvmin
                    k1=iswmin; k2=iswmax-1
                  else
                    i1=isumin; i2=isumin
                    j1=isvmin; j2=isvmax-1
                    k1=iswmin; k2=iswmax-1
                  endif
                  write (t12,'(2a,i3.3,a,i3.3)')
     &              facstr,'sld',NSB(ICFD)+1,'_s',ISUR
                  call MKCFDBC(.false.,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &              SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
                  icreated=icreated+1

                endif
                if (isclose(2).eq.2) then
                  if (ifac.eq.1 .or. ifac.eq.2) then
                    i1=isumin;   i2=isumax-1
                    j1=isvmax-1; j2=isvmax-1
                    k1=iswmin;   k2=iswmax-1
                  else
                    i1=isumax-1; i2=isumax-1
                    j1=isvmin;   j2=isvmax-1
                    k1=iswmin;   k2=iswmax-1
                  endif
                  write (t12,'(2a,i3.3,a,i3.3)')
     &              facstr,'sld',NSB(ICFD)+1,'_s',ISUR
                  call MKCFDBC(.false.,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &              SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
                  icreated=icreated+1
                endif
                if (isclose(3).eq.2) then
                  if (ifac.eq.5 .or. ifac.eq.6) then
                    if (isclose(1).eq.2) then                      
                      i1=isumin+1
                    else
                      i1=isumin
                    endif
                    if (isclose(2).eq.2) then
                      i2=isumax-2
                    else
                      i2=isumax-1
                    endif
                    j1=isvmin; j2=isvmin
                    k1=iswmin; k2=iswmax-1
                  else
                    if (isclose(1).eq.2) then                      
                      i1=isumin+1
                    else
                      i1=isumin
                    endif
                    if (isclose(2).eq.2) then
                      i2=isumax-2
                    else
                      i2=isumax-1
                    endif
                    j1=isvmin; j2=isvmax-1
                    k1=iswmin; k2=iswmin
                  endif
                  write (t12,'(2a,i3.3,a,i3.3)')
     &              facstr,'sld',NSB(ICFD)+1,'_s',ISUR
                  call MKCFDBC(.false.,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &              SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
                  icreated=icreated+1
                endif
                if (isclose(4).eq.2) then
                  if (ifac.eq.5 .or. ifac.eq.6) then
                    if (isclose(1).eq.2) then                      
                      i1=isumin+1
                    else
                      i1=isumin
                    endif
                    if (isclose(2).eq.2) then
                      i2=isumax-2
                    else
                      i2=isumax-1
                    endif
                    j1=isvmax-1; j2=isvmax-1
                    k1=iswmin;   k2=iswmax-1
                  else
                    if (isclose(1).eq.2) then                      
                      i1=isumin+1
                    else
                      i1=isumin
                    endif
                    if (isclose(2).eq.2) then
                      i2=isumax-2
                    else
                      i2=isumax-1
                    endif
                    j1=isvmin;   j2=isvmax-1
                    k1=iswmax-1; k2=iswmax-1
                  endif
                  write (t12,'(2a,i3.3,a,i3.3)')
     &              facstr,'sld',NSB(ICFD)+1,'_s',ISUR
                  call MKCFDBC(.false.,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &              SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
                  icreated=icreated+1
                endif

                if (icreated.eq.1) then
                  ISHSB=NSB(ICFD); ISHAO=0; ISHBLK=0; ISHGEO=-9999
                  call redraw(IER)
                  call continuebox(
     &              'Preview of parent surface boundary condition.',
     &              ' ','ok')
                else

C Show the new boundary conditions highlighted, one by one, with a half
C second delay.
   62             do i=1,icreated
                    ISHSB=NSB(ICFD)-i+1; ISHAO=0; ISHBLK=0; ISHGEO=-9999
                    call redraw(IER)
                    if (i.lt.icreated) call pausems(500)
                  enddo
                  IR=1
                  call EASKAB(
     &              'Preview of parent surface boundary conditions.',
     &              ' ','ok','show again',IR,nbhelp)
                  if (IR.eq.2) goto 62
                endif

C Delete BC(s) created for preview.
                do i=1,icreated
                  call DELCFDBC(NVOL(ICFD))
                enddo

C Show child surface.
                ISHSB=0; ISHAO=0; ISHBLK=0; ISHGEO=-icldsur
                CALL redraw(IER)
                CALL continuebox('Child surface.',' ','ok')

C Preview child BC with parent BC(s).
   63           if (ifac.eq.1 .or. ifac.eq.2) then
                  i1=isumin; i2=isumax-1
                  if (isclose(1).eq.2) then
                    j1=iparext(1)+1
                  else
                    j1=iparext(1)
                  endif
                  if (isclose(2).eq.2) then
                    j2=iparext(2)-2
                  else
                    j2=iparext(2)-1
                  endif
                  if (isclose(3).eq.2) then
                    k1=iparext(3)+1
                  else
                    k1=iparext(3)
                  endif
                  if (isclose(4).eq.2) then
                    k2=iparext(4)-2
                  else
                    k2=iparext(4)-1
                  endif
                elseif (ifac.eq.3 .or. ifac.eq.4) then
                  if (isclose(1).eq.2) then
                    i1=iparext(1)+1
                  else
                    i1=iparext(1)
                  endif
                  if (isclose(2).eq.2) then
                    i2=iparext(2)-2
                  else
                    i2=iparext(2)-1
                  endif
                  j1=isvmin; j2=isvmax-1
                  if (isclose(3).eq.2) then
                    k1=iparext(3)+1
                  else
                    k1=iparext(3)
                  endif
                  if (isclose(4).eq.2) then
                    k2=iparext(4)-2
                  else
                    k2=iparext(4)-1
                  endif
                elseif (ifac.eq.5 .or. ifac.eq.6) then
                  if (isclose(1).eq.2) then
                    i1=iparext(1)+1
                  else
                    i1=iparext(1)
                  endif
                  if (isclose(2).eq.2) then
                    i2=iparext(2)-2
                  else
                    i2=iparext(2)-1
                  endif
                  if (isclose(3).eq.2) then
                    j1=iparext(3)+1
                  else
                    j1=iparext(3)
                  endif
                  if (isclose(4).eq.2) then
                    j2=iparext(4)-2
                  else
                    j2=iparext(4)-1
                  endif
                  k1=iswmin; k2=iswmax-1
                endif
                write (t12,'(2a,i3.3,a,i3.3)')
     &            facstr,'sld',NSB(ICFD)+1,'_s',icldsur
                call MKCFDBC(.false.,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &            SSNAME(icldcnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
                ISHSB=NSB(ICFD); ISHAO=0; ISHBLK=0; ISHGEO=-9999
                call redraw(IER)
                call continuebox(
     &            'Preview of child surface BC with parent BC(s).',
     &            ' ','ok')

C Delete preview BC.
                call DELCFDBC(NVOL(ICFD))

C Preview child BC without parent BC(s).
                if (ifac.eq.1 .or. ifac.eq.2) then
                  i1=isumin;     i2=isumax-1
                  j1=iparext(1); j2=iparext(2)-1
                  k1=iparext(3); k2=iparext(4)-1
                elseif (ifac.eq.3 .or. ifac.eq.4) then
                  i1=iparext(1); i2=iparext(2)-1
                  j1=isvmin;     j2=isvmax-1
                  k1=iparext(3); k2=iparext(4)-1
                elseif (ifac.eq.5 .or. ifac.eq.6) then
                  i1=iparext(1); i2=iparext(2)-1
                  j1=iparext(3); j2=iparext(4)-1
                  k1=iswmin;     k2=iswmax-1
                endif
                write (t12,'(2a,i3.3,a,i3.3)')
     &            facstr,'sld',NSB(ICFD)+1,'_s',icldsur
                call MKCFDBC(.false.,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &            SSNAME(icldcnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
                ISHSB=NSB(ICFD); ISHAO=0; ISHBLK=0; ISHGEO=-9999
                call redraw(IER)
                call EASKAB(
     &            'Preview of child surface BC without parent BC(s).',
     &              ' ','ok','show again',IR,nbhelp)

C Delete preview BC.
                call DELCFDBC(NVOL(ICFD))

                if (IR.eq.2) goto 63

                ok=.true.
              else
                ok=.true.
              endif
              if (ok) goto 61
   70         CONTINUE
            else
C Otherwise, create solid BC(s) for parent.
              lumped=.false.
            endif

C If there is more than one child, we don't handle it.
          else
            write(outs,'(2a)')
     &        ' Parent has more than one child and is not narrow,',
     &        ' skipping parent and children.'
            call EDISP(IUOUT,outs)
            ISTAT=3
            goto 999
          endif
        endif

C ibctyp = 1:  solid BC(s) for the whole surface
C ibctyp = 2:  single opening BC for the whole surface
C ibctyp = 3:  opening BC for bottom row of cells, solid BC for the rest
C ibctyp = 4:  opening BC for middle row of cells, solid BCs for the rest
C ibctyp = 5:  opening BC for 25% of height, at bottom edge, solid BC for
C              the rest
C ibctyp = 6:  opening BC for 25% of width, at vertical centre line,
C              solid BCs for rest
C ibctyp = 7:  opening BCs for 12.5% of height, at top and bottom edges,
C              solid BC for rest
C ibctyp = 8:  opening BC for 75% of height and width, at bottom edge,
C              solid BCs for rest
C ibctyp = 9:  for frames; opening BC as close to middle row of parent
C              and child as we can get
C ibctyp = 10: for frames; opening BC for 20% of child width, as close
C              to top edge as we can get

Create child(ren) BCs.
        do icld=1,NBCHILD(icnn)
          if (.not.dochild(icld)) CYCLE
          icldcnn=ICHILD(icnn,icld)
          icldsur=IE1(icldcnn)
          iibctyp=icldbctyp(icld)

C If this child is also a parent, it will have already been created in
C its own invocation of this function; do not create BCs for it.
          if (NBCHILD(icldcnn).gt.0) CYCLE

          write(outs,'(a,i3.3,2a)')' Creating BC(s) for surface ',
     &      icldsur,SSNAME(icldcnn)(1:lnblnk(SSNAME(icldcnn))),' ...'
          ISHSB=0; ISHAO=0; ISHBLK=0; IFACES=1; ISHGEO=-icldsur
          call redraw(IER)
          call pausems(500)
          ISHGEO=-9999

C Single solid BC for the whole surface.
          if (iibctyp.eq.1) then
            if (ifac.eq.1 .or. ifac.eq.2) then
              i1=isumin;          i2=isumax-1
              j1=icldext(icld,1); j2=icldext(icld,2)-1
              k1=icldext(icld,3); k2=icldext(icld,4)-1
            elseif (ifac.eq.3 .or. ifac.eq.4) then
              i1=icldext(icld,1); i2=icldext(icld,2)-1
              j1=isvmin;          j2=isvmax-1
              k1=icldext(icld,3); k2=icldext(icld,4)-1
            elseif (ifac.eq.5 .or. ifac.eq.6) then
              i1=icldext(icld,1); i2=icldext(icld,2)-1
              j1=icldext(icld,3); j2=icldext(icld,4)-1
              k1=iswmin;          k2=iswmax-1
            endif
            write (t12,'(2a,i3.3,a,i3.3)')
     &        facstr,'sld',NSB(ICFD)+1,'_s',icldsur
            call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &        SSNAME(icldcnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)

C Single opening BC for the whole surface.
          elseif (iibctyp.eq.2) then
            if (ifac.eq.1 .or. ifac.eq.2) then
              i1=isumin;          i2=isumax-1
              j1=icldext(icld,1); j2=icldext(icld,2)-1
              k1=icldext(icld,3); k2=icldext(icld,4)-1
            elseif (ifac.eq.3 .or. ifac.eq.4) then
              i1=icldext(icld,1); i2=icldext(icld,2)-1
              j1=isvmin;          j2=isvmax-1
              k1=icldext(icld,3); k2=icldext(icld,4)-1
            elseif (ifac.eq.5 .or. ifac.eq.6) then
              i1=icldext(icld,1); i2=icldext(icld,2)-1
              j1=icldext(icld,3); j2=icldext(icld,4)-1
              k1=iswmin;          k2=iswmax-1
            endif
            write (t12,'(2a,i3.3,a,i3.3)')
     &        facstr,'opn',NOPEN(ICFD)+1,'_s',icldsur
            call MKCFDBC(show,11,ifac,i1,i2,j1,j2,k1,k2,t12,
     &        SSNAME(icldcnn),dum,1,20.0,0.0,0.01,0.0,0.0,IER)

C Opening BC for bottom row of cells, solid BC for the rest.
          elseif (iibctyp.eq.3) then
            if (ifac.eq.1 .or. ifac.eq.2) then
              i1=isumin;          i2=isumax-1
              j1=icldext(icld,1); j2=icldext(icld,2)-1
              k1=icldext(icld,3); k2=icldext(icld,3)
            elseif (ifac.eq.3 .or. ifac.eq.4) then
              i1=icldext(icld,1); i2=icldext(icld,2)-1
              j1=isvmin;          j2=isvmax-1
              k1=icldext(icld,3); k2=icldext(icld,3)
            elseif (ifac.eq.5 .or. ifac.eq.6) then
              i1=icldext(icld,1); i2=icldext(icld,2)-1
              j1=icldext(icld,3); j2=icldext(icld,3)
              k1=iswmin;          k2=iswmax-1
            endif
            write (t12,'(2a,i3.3,a,i3.3)')
     &        facstr,'opn',NOPEN(ICFD)+1,'_s',icldsur
            call MKCFDBC(show,11,ifac,i1,i2,j1,j2,k1,k2,t12,
     &        SSNAME(icldcnn),dum,1,20.0,0.0,0.01,0.0,0.0,IER)

C Is there any more of the surface for a solid BC?
            if (icldext(icld,4)-1.gt.icldext(icld,3)) then
              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;            i2=isumax-1
                j1=icldext(icld,1);   j2=icldext(icld,2)-1
                k1=icldext(icld,3)+1; k2=icldext(icld,4)-1
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=icldext(icld,1);   i2=icldext(icld,2)-1
                j1=isvmin;            j2=isvmax-1
                k1=icldext(icld,3)+1; k2=icldext(icld,4)-1
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=icldext(icld,1);   i2=icldext(icld,2)-1
                j1=icldext(icld,3)+1; j2=icldext(icld,4)-1
                k1=iswmin;            k2=iswmax-1
              endif
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'sld',NSB(ICFD)+1,'_s',icldsur
              call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icldcnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)

            endif

C Opening BC for middle row of cells, solid BCs for the rest.
          elseif (iibctyp.eq.4) then

C Find middle row.
            irows=icldext(icld,4)-icldext(icld,3)
            imidrow=icldext(icld,3)+irows/2

            if (ifac.eq.1 .or. ifac.eq.2) then
              i1=isumin;          i2=isumax-1
              j1=icldext(icld,1); j2=icldext(icld,2)-1
              k1=imidrow;         k2=imidrow
            elseif (ifac.eq.3 .or. ifac.eq.4) then
              i1=icldext(icld,1); i2=icldext(icld,2)-1
              j1=isvmin;          j2=isvmax-1
              k1=imidrow;         k2=imidrow
            elseif (ifac.eq.5 .or. ifac.eq.6) then
              i1=icldext(icld,1); i2=icldext(icld,2)-1
              j1=imidrow;         j2=imidrow
              k1=iswmin;          k2=iswmax-1
            endif
            write (t12,'(2a,i3.3,a,i3.3)')
     &        facstr,'opn',NOPEN(ICFD)+1,'_s',icldsur
            call MKCFDBC(show,11,ifac,i1,i2,j1,j2,k1,k2,t12,
     &        SSNAME(icldcnn),dum,1,20.0,0.0,0.01,0.0,0.0,IER)

C Is there any more of the surface below for a solid BC?
            if (icldext(icld,3).lt.imidrow) then
              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;          i2=isumax-1
                j1=icldext(icld,1); j2=icldext(icld,2)-1
                k1=icldext(icld,3); k2=imidrow-1
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=icldext(icld,1); i2=icldext(icld,2)-1
                j1=isvmin;          j2=isvmax-1
                k1=icldext(icld,3); k2=imidrow-1
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=icldext(icld,1); i2=icldext(icld,2)-1
                j1=icldext(icld,3); j2=imidrow-1
                k1=iswmin;          k2=iswmax-1
              endif
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'sld',NSB(ICFD)+1,'_s',icldsur
              call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icldcnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
            endif

C Is there any more of the surface above for a solid BC?
            if ((icldext(icld,4)-1).gt.imidrow) then
              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;          i2=isumax-1
                j1=icldext(icld,1); j2=icldext(icld,2)-1
                k1=imidrow+1;       k2=icldext(icld,4)-1
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=icldext(icld,1); i2=icldext(icld,2)-1
                j1=isvmin;          j2=isvmax-1
                k1=imidrow+1;       k2=icldext(icld,4)-1
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=icldext(icld,1); i2=icldext(icld,2)-1
                j1=imidrow+1;       j2=icldext(icld,4)-1
                k1=iswmin;          k2=iswmax-1
              endif
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'sld',NSB(ICFD)+1,'_s',icldsur
              call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icldcnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
            endif

C Opening BC for 25% of height, at bottom edge, solid BC for the rest.
          elseif (iibctyp.eq.5) then

C Find 25% of height, rounded, at least 1 row.
            rrows=real(icldext(icld,4)-icldext(icld,3))
            iqtr=nint(rrows/4.0)
            if (iqtr.eq.0) iqtr=1

            if (ifac.eq.1 .or. ifac.eq.2) then
              i1=isumin;          i2=isumax-1
              j1=icldext(icld,1); j2=icldext(icld,2)-1
              k1=icldext(icld,3); k2=icldext(icld,3)+iqtr-1
            elseif (ifac.eq.3 .or. ifac.eq.4) then
              i1=icldext(icld,1); i2=icldext(icld,2)-1
              j1=isvmin;          j2=isvmax-1
              k1=icldext(icld,3); k2=icldext(icld,3)+iqtr-1
            elseif (ifac.eq.5 .or. ifac.eq.6) then
              i1=icldext(icld,1); i2=icldext(icld,2)-1
              j1=icldext(icld,3); j2=icldext(icld,3)+iqtr-1
              k1=iswmin;          k2=iswmax-1
            endif
            write (t12,'(2a,i3.3,a,i3.3)')
     &        facstr,'opn',NOPEN(ICFD)+1,'_s',icldsur
            call MKCFDBC(show,11,ifac,i1,i2,j1,j2,k1,k2,t12,
     &        SSNAME(icldcnn),dum,1,20.0,0.0,0.01,0.0,0.0,IER)

C Is there any more of the surface for a solid BC?
            if ((icldext(icld,4)-1).gt.(icldext(icld,3)+iqtr)) then
              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;               i2=isumax-1
                j1=icldext(icld,1);      j2=icldext(icld,2)-1
                k1=icldext(icld,3)+iqtr; k2=icldext(icld,4)-1
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=icldext(icld,1);      i2=icldext(icld,2)-1
                j1=isvmin;               j2=isvmax-1
                k1=icldext(icld,3)+iqtr; k2=icldext(icld,4)-1
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=icldext(icld,1);      i2=icldext(icld,2)-1
                j1=icldext(icld,3)+iqtr; j2=icldext(icld,4)-1
                k1=iswmin;               k2=iswmax-1
              endif
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'sld',NSB(ICFD)+1,'_s',icldsur
              call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icldcnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
            endif

C Opening BC for 25% of width, at vertical centre line, solid BCs for
C rest.
          elseif (iibctyp.eq.6) then

C Find 25% of width, rounded, at least 1 column, on vertical centre
C line.
            icols=icldext(icld,2)-icldext(icld,1)
            rcols=real(icols)
            iqtr=nint(rcols/4.0)
            if (iqtr.eq.0) iqtr=1
            imidcol=icldext(icld,1)+icols/2
            if (mod(iqtr,2).eq.0) then
              imidqtrmin=imidcol-iqtr/2
              imidqtrmax=imidcol+iqtr/2-1
            else
              imidqtrmin=imidcol-iqtr/2
              imidqtrmax=imidcol+iqtr/2
            endif

            if (ifac.eq.1 .or. ifac.eq.2) then
              i1=isumin;          i2=isumax-1
              j1=imidqtrmin;      j2=imidqtrmax
              k1=icldext(icld,3); k2=icldext(icld,4)-1
            elseif (ifac.eq.3 .or. ifac.eq.4) then
              i1=imidqtrmin;      i2=imidqtrmax
              j1=isvmin;          j2=isvmax-1
              k1=icldext(icld,3); k2=icldext(icld,4)-1
            elseif (ifac.eq.5 .or. ifac.eq.6) then
              i1=imidqtrmin;      i2=imidqtrmax
              j1=icldext(icld,3); j2=icldext(icld,4)-1
              k1=iswmin;          k2=iswmax-1
            endif
            write (t12,'(2a,i3.3,a,i3.3)')
     &        facstr,'opn',NOPEN(ICFD)+1,'_s',icldsur
            call MKCFDBC(show,11,ifac,i1,i2,j1,j2,k1,k2,t12,
     &        SSNAME(icldcnn),dum,1,20.0,0.0,0.01,0.0,0.0,IER)

C Is there any more of the surface to the left for a solid BC?
            if ((icldext(icld,1)).lt.imidqtrmin) then
              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;          i2=isumax-1
                j1=icldext(icld,1); j2=imidqtrmin-1
                k1=icldext(icld,3); k2=icldext(icld,4)-1
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=icldext(icld,1); i2=imidqtrmin-1
                j1=isvmin;          j2=isvmax-1
                k1=icldext(icld,3); k2=icldext(icld,4)-1
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=icldext(icld,1); i2=imidqtrmin-1
                j1=icldext(icld,3); j2=icldext(icld,4)-1
                k1=iswmin;          k2=iswmax-1
              endif
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'sld',NSB(ICFD)+1,'_s',icldsur
              call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icldcnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
            endif

C Is there any more of the surface to the right for a solid BC?
            if ((icldext(icld,2)-1).gt.imidqtrmax) then
              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;          i2=isumax-1
                j1=imidqtrmax+1;    j2=icldext(icld,2)-1
                k1=icldext(icld,3); k2=icldext(icld,4)-1
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=imidqtrmax+1;    i2=icldext(icld,2)-1
                j1=isvmin;          j2=isvmax-1
                k1=icldext(icld,4); k2=icldext(icld,4)-1
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=imidqtrmax+1;    i2=icldext(icld,2)-1
                j1=icldext(icld,3); j2=icldext(icld,4)-1
                k1=iswmin;          k2=iswmax-1
              endif
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'sld',NSB(ICFD)+1,'_s',icldsur
              call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icldcnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
            endif

C Opening BCs for 12.5% of height, at top and bottom edges, solid BC for
C rest.
          elseif (iibctyp.eq.7) then

C Find 12.5% of height, rounded, at least 1 row.
            irows=icldext(icld,4)-icldext(icld,3)
            rrows=real(irows)
            i8th=nint(rrows/8.0)
            if (i8th.eq.0) i8th=1

C Is there enough space for both openings?
            ok=.true.
   80       if (irows.lt.(i8th*2)) then
              if (i8th.gt.1) then
                i8th=i8th-1
              else
                write(outs,'(3a)')' Not enough space for openings at',
     &            ' top and bottom of this sash window; only one',
     &            ' opening will be created.'
                call EDISP(IUOUT,outs)
                ok=.false.
              endif
              goto 80
            endif

C Opening at bottom edge.
            if (ifac.eq.1 .or. ifac.eq.2) then
              i1=isumin;          i2=isumax-1
              j1=icldext(icld,1); j2=icldext(icld,2)-1
              k1=icldext(icld,3); k2=icldext(icld,3)+i8th-1
            elseif (ifac.eq.3 .or. ifac.eq.4) then
              i1=icldext(icld,1); i2=icldext(icld,2)-1
              j1=isvmin;          j2=isvmax-1
              k1=icldext(icld,3); k2=icldext(icld,3)+i8th-1
            elseif (ifac.eq.5 .or. ifac.eq.6) then
              i1=icldext(icld,1); i2=icldext(icld,2)-1
              j1=icldext(icld,3); j2=icldext(icld,3)+i8th-1
              k1=iswmin;          k2=iswmax-1
            endif
            write (t12,'(2a,i3.3,a,i3.3)')
     &        facstr,'opn',NOPEN(ICFD)+1,'_s',icldsur
            call MKCFDBC(show,11,ifac,i1,i2,j1,j2,k1,k2,t12,
     &        SSNAME(icldcnn),dum,1,20.0,0.0,0.01,0.0,0.0,IER)

C Is there room for the top opening?
            if (ok) then
              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;               i2=isumax-1
                j1=icldext(icld,1);      j2=icldext(icld,2)-1
                k1=icldext(icld,4)-i8th; k2=icldext(icld,4)-1
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=icldext(icld,1);      i2=icldext(icld,2)-1
                j1=isvmin;               j2=isvmax-1
                k1=icldext(icld,4)-i8th; k2=icldext(icld,4)-1
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=icldext(icld,1);      i2=icldext(icld,2)-1
                j1=icldext(icld,4)-i8th; j2=icldext(icld,4)-1
                k1=iswmin;               k2=iswmax-1
              endif
              itop=icldext(icld,4)-i8th-1
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'opn',NOPEN(ICFD)+1,'_s',icldsur
              call MKCFDBC(show,11,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icldcnn),dum,1,20.0,0.0,0.01,0.0,0.0,IER)
            else
              itop=icldext(icld,4)-1
            endif

C Is there any more of the surface for a solid BC?
            if ((icldext(icld,3)+i8th).le.itop) then
              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;               i2=isumax-1
                j1=icldext(icld,1);      j2=icldext(icld,2)-1
                k1=icldext(icld,3)+i8th; k2=itop
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=icldext(icld,1);      i2=icldext(icld,2)-1
                j1=isvmin;               j2=isvmax-1
                k1=icldext(icld,3)+i8th; k2=itop
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=icldext(icld,1);      i2=icldext(icld,2)-1
                j1=icldext(icld,3)+i8th; j2=itop
                k1=iswmin;               k2=iswmax-1
              endif
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'sld',NSB(ICFD)+1,'_s',icldsur
              call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icldcnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
            endif

C Opening BC for 75% of height and width, at bottom edge, solid BCs for
C rest.            
          elseif (iibctyp.eq.8) then

C Find 75% of height, rounded, at least 1 row, at bottom edge.
            irows=icldext(icld,4)-icldext(icld,3)
            rrows=real(irows)
            i3qtrhgt=nint(3.0*rrows/4.0)
            if (i3qtrhgt.eq.0) i3qtrhgt=1
            ibtm3qtrhgtmin=icldext(icld,3)
            ibtm3qtrhgtmax=icldext(icld,3)+i3qtrhgt-1

C Find 75% of width, rounded, at least 1 column, on vertical centre
C line.
            icols=icldext(icld,2)-icldext(icld,1)
            rcols=real(icols)
            i3qtrwdt=nint(3.0*rcols/4.0)
            if (i3qtrwdt.eq.0) i3qtrwdt=1
            imidcol=icldext(icld,1)+icols/2
            if (mod(i3qtrwdt,2).eq.0) then
              imid3qtrwdtmin=imidcol-i3qtrwdt/2
              imid3qtrwdtmax=imidcol+i3qtrwdt/2-1
            else
              imid3qtrwdtmin=imidcol-i3qtrwdt/2
              imid3qtrwdtmax=imidcol+i3qtrwdt/2
            endif

C Opening BC.
            if (ifac.eq.1 .or. ifac.eq.2) then
              i1=isumin;         i2=isumax-1
              j1=imid3qtrwdtmin; j2=imid3qtrwdtmax
              k1=ibtm3qtrhgtmin; k2=ibtm3qtrhgtmax
            elseif (ifac.eq.3 .or. ifac.eq.4) then
              i1=imid3qtrwdtmin; i2=imid3qtrwdtmax
              j1=isvmin;         j2=isvmax-1
              k1=ibtm3qtrhgtmin; k2=ibtm3qtrhgtmax
            elseif (ifac.eq.5 .or. ifac.eq.6) then
              i1=imid3qtrwdtmin; i2=imid3qtrwdtmax
              j1=ibtm3qtrhgtmin; j2=ibtm3qtrhgtmax
              k1=iswmin;         k2=iswmax-1
            endif
            write (t12,'(2a,i3.3,a,i3.3)')
     &        facstr,'opn',NOPEN(ICFD)+1,'_s',icldsur
            call MKCFDBC(show,11,ifac,i1,i2,j1,j2,k1,k2,t12,
     &        SSNAME(icldcnn),dum,1,20.0,0.0,0.01,0.0,0.0,IER)

C Is there room left at the top?
            if (ibtm3qtrhgtmax.lt.(icldext(icld,4)-1)) then
              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;           i2=isumax-1
                j1=icldext(icld,1);  j2=icldext(icld,2)-1
                k1=ibtm3qtrhgtmax+1; k2=icldext(icld,4)-1
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=icldext(icld,1);  i2=icldext(icld,2)-1
                j1=isvmin;           j2=isvmax-1
                k1=ibtm3qtrhgtmax+1; k2=icldext(icld,4)-1
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=icldext(icld,1);  i2=icldext(icld,2)-1
                j1=ibtm3qtrhgtmax+1; j2=icldext(icld,4)-1
                k1=iswmin;           k2=iswmax-1
              endif
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'sld',NSB(ICFD)+1,'_s',icldsur
              call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icldcnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
              ok=.true.
            else
              ok=.false.
            endif

C Is there room at the left?
            if (imid3qtrwdtmin.gt.icldext(icld,1)) then
              if (ok) then
                itop=icldext(icld,4)-2
              else
                itop=icldext(icld,4)-1
              endif
              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;          i2=isumax-1
                j1=icldext(icld,1); j2=imid3qtrwdtmin-1
                k1=icldext(icld,3); k2=itop
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=icldext(icld,1); i2=imid3qtrwdtmin-1
                j1=isvmin;          j2=isvmax-1
                k1=icldext(icld,3); k2=itop
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=icldext(icld,1); i2=imid3qtrwdtmin-1
                j1=icldext(icld,3); j2=itop
                k1=iswmin;          k2=iswmax-1
              endif
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'sld',NSB(ICFD)+1,'_s',icldsur
              call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icldcnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
            endif

C Is there room at the right?
            if (imid3qtrwdtmax.lt.(icldext(icld,2)-1)) then
              if (ok) then
                itop=icldext(icld,4)-2
              else
                itop=icldext(icld,4)-1
              endif
              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;           i2=isumax-1
                j1=imid3qtrwdtmax+1; j2=icldext(icld,2)-1
                k1=icldext(icld,3);  k2=itop
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=imid3qtrwdtmax+1; i2=icldext(icld,2)-1
                j1=isvmin;           j2=isvmax-1
                k1=icldext(icld,3);  k2=itop
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=imid3qtrwdtmax+1; i2=icldext(icld,2)-1
                j1=icldext(icld,3);  j2=itop
                k1=iswmin;           k2=iswmax-1
              endif
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'sld',NSB(ICFD)+1,'_s',icldsur
              call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icldcnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
            endif

C 9 and 10 are for frames.
          elseif (iibctyp.eq.9 .or. iibctyp.eq.10) then
            write(6,'(a)')
     &' This surface is a frame, but has no children. Check usage tags.'
            call EDISP(IUOUT,outs)
          endif
        enddo

C Create parent surface BC.
        if (dothisur .and. .not.lumped) then

          write(outs,'(a,i3.3,2a)')' Creating BC(s) for surface ',
     &      isur,SSNAME(icnn)(1:lnblnk(SSNAME(icnn))),' ...'
          ISHSB=0; ISHAO=0; ISHBLK=0; IFACES=1; ISHGEO=-isur
          call redraw(IER)
          call pausems(500)
          ISHGEO=-9999

C Solid BCs drawn around child.
          if (ibctyp.eq.1) then

C Is there space on the top?
            irows=iparext(4)-icldext(1,4)
            if (irows.gt.0) then
              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;       i2=isumax-1
                j1=iparext(1);   j2=iparext(2)-1
                k1=icldext(1,4); k2=iparext(4)-1
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=iparext(1);   i2=iparext(2)-1
                j1=isvmin;       j2=isvmax-1
                k1=icldext(1,4); k2=iparext(4)-1
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=iparext(1);   i2=iparext(2)-1
                j1=icldext(1,4); j2=iparext(4)-1
                k1=iswmin;       k2=iswmax-1
              endif
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'sld',NSB(ICFD)+1,'_s',isur
              call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
              itop=iparext(4)-1-irows
            else
              itop=iparext(4)-1
            endif

C Is there space on the bottom?
            irows=icldext(1,3)-iparext(3)
            if (irows.gt.0) then
              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;     i2=isumax-1
                j1=iparext(1); j2=iparext(2)-1
                k1=iparext(3); k2=icldext(1,3)-1
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=iparext(1); i2=iparext(2)-1
                j1=isvmin;     j2=isvmax-1
                k1=iparext(3); k2=icldext(1,3)-1
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=iparext(1); i2=iparext(2)-1
                j1=iparext(3); j2=icldext(1,3)-1
                k1=iswmin;     k2=iswmax-1
              endif
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'sld',NSB(ICFD)+1,'_s',isur
              call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
              ibot=iparext(3)+irows
            else
              ibot=iparext(3)
            endif

C Is there space on the left?
            icols=icldext(1,1)-iparext(1)
            if (icols.gt.0) then
              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;     i2=isumax-1
                j1=iparext(1); j2=icldext(1,1)-1
                k1=ibot;       k2=itop
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=iparext(1); i2=icldext(1,1)-1
                j1=isvmin;     j2=isvmax-1
                k1=ibot;       k2=itop
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=iparext(1); i2=icldext(1,1)-1
                j1=ibot;       j2=itop
                k1=iswmin;     k2=iswmax-1
              endif
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'sld',NSB(ICFD)+1,'_s',isur
              call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
            endif

C Is there space on the right?
            icols=iparext(2)-icldext(1,2)
            if (icols.gt.0) then
              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;       i2=isumax-1
                j1=icldext(1,2); j2=iparext(2)-1
                k1=ibot;         k2=itop
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=icldext(1,2); i2=iparext(2)-1
                j1=isvmin;       j2=isvmax-1
                k1=ibot;         k2=itop
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=icldext(1,2); i2=iparext(2)-1
                j1=ibot;         j2=itop
                k1=iswmin;       k2=iswmax-1
              endif
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'sld',NSB(ICFD)+1,'_s',isur
              call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
            endif

C Parent is not compatible with ibctyp=2 to 8.
          elseif (ibctyp.ge.2 .and. ibctyp.le.8) then
            write(outs,'(2a)')' Parent surface BCs could not be',
     &        ' created based on usage flags.'
            call EDISP(IUOUT,outs)
            ISTAT=5
            goto 999

C Frame with crack. Use first edge where there is space (in order of
C top, bottom, left, right) for crack, other edges where there is space
C use a solid BC.
          elseif (ibctyp.eq.9) then

            ok=.false.
C Is there space on the top?
            irows=iparext(4)-icldext(1,4)
            if (irows.gt.0) then
              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;       i2=isumax-1
                j1=iparext(1);   j2=iparext(2)-1
                k1=icldext(1,4); k2=icldext(1,4)
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=iparext(1);   i2=iparext(2)-1
                j1=isvmin;       j2=isvmax-1
                k1=icldext(1,4); k2=icldext(1,4)
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=iparext(1);   i2=iparext(2)-1
                j1=icldext(1,4); j2=icldext(1,4)
                k1=iswmin;       k2=iswmax-1
              endif
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'opn',NOPEN(ICFD)+1,'_s',isur
              call MKCFDBC(show,11,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icnn),dum,1,20.0,0.0,0.01,0.0,0.0,IER)
              ok=.true.

C Any space left over that needs filling with a solid bc?
              if (irows.gt.1) then
                if (ifac.eq.1 .or. ifac.eq.2) then
                  i1=isumin;         i2=isumax-1
                  j1=iparext(1);     j2=iparext(2)-1
                  k1=icldext(1,4)+1; k2=iparext(4)-1
                elseif (ifac.eq.3 .or. ifac.eq.4) then
                  i1=iparext(1);     i2=iparext(2)-1
                  j1=isvmin;         j2=isvmax-1
                  k1=icldext(1,4)+1; k2=iparext(4)-1
                elseif (ifac.eq.5 .or. ifac.eq.6) then
                  i1=iparext(1);     i2=iparext(2)-1
                  j1=icldext(1,4)+1; j2=iparext(4)-1
                  k1=iswmin;         k2=iswmax-1
                endif
                write (t12,'(2a,i3.3,a,i3.3)')
     &            facstr,'sld',NSB(ICFD)+1,'_s',isur
                call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &            SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
              endif

              itop=iparext(4)-1-irows
            else
              itop=iparext(4)-1
            endif

C Is there space on the bottom?
            irows=icldext(1,3)-iparext(3)
            if (irows.gt.0) then
              if (ok) then ! opening already created?
                if (ifac.eq.1 .or. ifac.eq.2) then
                  i1=isumin;     i2=isumax-1
                  j1=iparext(1); j2=iparext(2)-1
                  k1=iparext(3); k2=icldext(1,3)-1
                elseif (ifac.eq.3 .or. ifac.eq.4) then
                  i1=iparext(1); i2=iparext(2)-1
                  j1=isvmin;     j2=isvmax-1
                  k1=iparext(3); k2=icldext(1,3)-1
                elseif (ifac.eq.5 .or. ifac.eq.6) then
                  i1=iparext(1); i2=iparext(2)-1
                  j1=iparext(3); j2=icldext(1,3)-1
                  k1=iswmin;     k2=iswmax-1
                endif
                write (t12,'(2a,i3.3,a,i3.3)')
     &            facstr,'sld',NSB(ICFD)+1,'_s',isur
                call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &            SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
              else
                if (ifac.eq.1 .or. ifac.eq.2) then
                  i1=isumin;         i2=isumax-1
                  j1=iparext(1);     j2=iparext(2)-1
                  k1=icldext(1,3)-1; k2=icldext(1,3)-1
                elseif (ifac.eq.3 .or. ifac.eq.4) then
                  i1=iparext(1);     i2=iparext(2)-1
                  j1=isvmin;         j2=isvmax-1
                  k1=icldext(1,3)-1; k2=icldext(1,3)-1
                elseif (ifac.eq.5 .or. ifac.eq.6) then
                  i1=iparext(1);     i2=iparext(2)-1
                  j1=icldext(1,3)-1; j2=icldext(1,3)-1
                  k1=iswmin;         k2=iswmax-1
                endif
                write (t12,'(2a,i3.3,a,i3.3)')
     &            facstr,'opn',NOPEN(ICFD)+1,'_s',isur
                call MKCFDBC(show,11,ifac,i1,i2,j1,j2,k1,k2,t12,
     &            SSNAME(icnn),dum,1,20.0,0.0,0.01,0.0,0.0,IER)
                ok=.true.

C Any space left over that needs filling with a solid bc?
                if (irows.gt.1) then
                  if (ifac.eq.1 .or. ifac.eq.2) then
                    i1=isumin;     i2=isumax-1
                    j1=iparext(1); j2=iparext(2)-1
                    k1=iparext(3); k2=icldext(1,3)-2
                  elseif (ifac.eq.3 .or. ifac.eq.4) then
                    i1=iparext(1); i2=iparext(2)-1
                    j1=isvmin;     j2=isvmax-1
                    k1=iparext(3); k2=icldext(1,3)-2
                  elseif (ifac.eq.5 .or. ifac.eq.6) then
                    i1=iparext(1); i2=iparext(2)-1
                    j1=iparext(3); j2=icldext(1,3)-2
                    k1=iswmin;     k2=iswmax-1
                  endif
                  write (t12,'(2a,i3.3,a,i3.3)')
     &              facstr,'sld',NSB(ICFD)+1,'_s',isur
                  call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &              SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
                endif
              endif
              ibot=iparext(3)+irows
            else
              ibot=iparext(3)
            endif

C Is there space on the left?
            icols=icldext(1,1)-iparext(1)
            if (icols.gt.0) then
              if (ok) then
                if (ifac.eq.1 .or. ifac.eq.2) then
                  i1=isumin;     i2=isumax-1
                  j1=iparext(1); j2=icldext(1,1)-1
                  k1=ibot;       k2=itop
                elseif (ifac.eq.3 .or. ifac.eq.4) then
                  i1=iparext(1); i2=icldext(1,1)-1
                  j1=isvmin;     j2=isvmax-1
                  k1=ibot;       k2=itop
                elseif (ifac.eq.5 .or. ifac.eq.6) then
                  i1=iparext(1); i2=icldext(1,1)-1
                  j1=ibot;       j2=itop
                  k1=iswmin;     k2=iswmax-1
                endif
                write (t12,'(2a,i3.3,a,i3.3)')
     &            facstr,'sld',NSB(ICFD)+1,'_s',isur
                call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &            SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
              else
                if (ifac.eq.1 .or. ifac.eq.2) then
                  i1=isumin;         i2=isumax-1
                  j1=icldext(1,1)-1; j2=icldext(1,1)-1
                  k1=ibot;           k2=itop
                elseif (ifac.eq.3 .or. ifac.eq.4) then
                  i1=icldext(1,1)-1; i2=icldext(1,1)-1
                  j1=isvmin;         j2=isvmax-1
                  k1=ibot;           k2=itop
                elseif (ifac.eq.5 .or. ifac.eq.6) then
                  i1=icldext(1,1)-1; i2=icldext(1,1)-1
                  j1=ibot;           j2=itop
                  k1=iswmin;         k2=iswmax-1
                endif
                write (t12,'(2a,i3.3,a,i3.3)')
     &            facstr,'opn',NOPEN(ICFD)+1,'_s',isur
                call MKCFDBC(show,11,ifac,i1,i2,j1,j2,k1,k2,t12,
     &            SSNAME(icnn),dum,1,20.0,0.0,0.01,0.0,0.0,IER)

C Any space left over that needs filling with a solid bc?
                if (icols.gt.1) then
                  if (ifac.eq.1 .or. ifac.eq.2) then
                    i1=isumin;     i2=isumax-1
                    j1=iparext(1); j2=icldext(1,1)-2
                    k1=ibot;       k2=itop
                  elseif (ifac.eq.3 .or. ifac.eq.4) then
                    i1=iparext(1); i2=icldext(1,1)-2
                    j1=isvmin;     j2=isvmax-1
                    k1=ibot;       k2=itop
                  elseif (ifac.eq.5 .or. ifac.eq.6) then
                    i1=iparext(1); i2=icldext(1,1)-2
                    j1=ibot;       j2=itop
                    k1=iswmin;     k2=iswmax-1
                  endif
                  write (t12,'(2a,i3.3,a,i3.3)')
     &              facstr,'sld',NSB(ICFD)+1,'_s',isur
                  call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &              SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
                endif
              endif
            endif

C Is there space on the right?
            icols=iparext(2)-icldext(1,2)
            if (icols.gt.0) then
              if (ok) then
                if (ifac.eq.1 .or. ifac.eq.2) then
                  i1=isumin;       i2=isumax-1
                  j1=icldext(1,2); j2=iparext(2)-1
                  k1=ibot;         k2=itop
                elseif (ifac.eq.3 .or. ifac.eq.4) then
                  i1=icldext(1,2); i2=iparext(2)-1
                  j1=isvmin;       j2=isvmax-1
                  k1=ibot;         k2=itop
                elseif (ifac.eq.5 .or. ifac.eq.6) then
                  i1=icldext(1,2); i2=iparext(2)-1
                  j1=ibot;         j2=itop
                  k1=iswmin;       k2=iswmax-1
                endif
                write (t12,'(2a,i3.3,a,i3.3)')
     &            facstr,'sld',NSB(ICFD)+1,'_s',isur
                call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &            SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
              else
                if (ifac.eq.1 .or. ifac.eq.2) then
                  i1=isumin;       i2=isumax-1
                  j1=icldext(1,2); j2=icldext(1,2)
                  k1=ibot;         k2=itop
                elseif (ifac.eq.3 .or. ifac.eq.4) then
                  i1=icldext(1,2); i2=icldext(1,2)
                  j1=isvmin;       j2=isvmax-1
                  k1=ibot;         k2=itop
                elseif (ifac.eq.5 .or. ifac.eq.6) then
                  i1=icldext(1,2); i2=icldext(1,2)
                  j1=ibot;         j2=itop
                  k1=iswmin;       k2=iswmax-1
                endif
                write (t12,'(2a,i3.3,a,i3.3)')
     &            facstr,'opn',NOPEN(ICFD)+1,'_s',isur
                call MKCFDBC(show,11,ifac,i1,i2,j1,j2,k1,k2,t12,
     &            SSNAME(icnn),dum,1,20.0,0.0,0.01,0.0,0.0,IER)

C Any space left over that needs filling with a solid bc?
                if (icols.gt.1) then
                  if (ifac.eq.1 .or. ifac.eq.2) then
                    i1=isumin;         i2=isumax-1
                    j1=icldext(1,2)+1; j2=iparext(2)-1
                    k1=ibot;           k2=itop
                  elseif (ifac.eq.3 .or. ifac.eq.4) then
                    j1=icldext(1,2)+1; j2=iparext(2)-1
                    j1=isvmin;         j2=isvmax-1
                    k1=ibot;           k2=itop
                  elseif (ifac.eq.5 .or. ifac.eq.6) then
                    j1=icldext(1,2)+1; j2=iparext(2)-1
                    j1=ibot;           j2=itop
                    k1=iswmin;         k2=iswmax-1
                  endif
                  write (t12,'(2a,i3.3,a,i3.3)')
     &              facstr,'sld',NSB(ICFD)+1,'_s',isur
                  call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &              SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
                endif
              endif
            endif

C Frame with a trickle vent. Use first edge where there is space (in
C order of top, bottom, left, right) for vent, other edges where there
C is space use a solid BC.
          elseif (ibctyp.eq.10) then

c Find 20% of frame inside width, rounded, at least 1 cell.            
            icols=icldext(icld,2)-icldext(icld,1)
            rcols=real(icols)
            iwdt=nint(rcols/5.0)
            if (iwdt.eq.0) iwdt=1
            icols=iparext(2)-iparext(1)

            ok=.false.

C Is there space on the top?
            irows=iparext(4)-icldext(1,4)
            if (irows.gt.0) then

C Find cells to put the trickle vent in the middle.
              imidcol=icldext(icld,1)+icols/2
              if (mod(iwdt,2).eq.0) then
                iwdtmin=imidcol-iwdt/2
                iwdtmax=imidcol+iwdt/2-1
              else
                iwdtmin=imidcol-iwdt/2
                iwdtmax=imidcol+iwdt/2
              endif

              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;       i2=isumax-1
                j1=iwdtmin;      j2=iwdtmax
                k1=icldext(1,4); k2=icldext(1,4)
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=iwdtmin;      i2=iwdtmax
                j1=isvmin;       j2=isvmax-1
                k1=icldext(1,4); k2=icldext(1,4)
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=iwdtmin;      i2=iwdtmax
                j1=icldext(1,4); j2=icldext(1,4)
                k1=iswmin;       k2=iswmax-1
              endif
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'opn',NOPEN(ICFD)+1,'_s',isur
              call MKCFDBC(show,11,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icnn),dum,1,20.0,0.0,0.01,0.0,0.0,IER)
              ok=.true.

              if (iparext(1).lt.iwdtmin) then
                if (ifac.eq.1 .or. ifac.eq.2) then
                  i1=isumin;       i2=isumax-1
                  j1=iparext(1);   j2=iwdtmin-1
                  k1=icldext(1,4); k2=icldext(1,4)
                elseif (ifac.eq.3 .or. ifac.eq.4) then
                  i1=iparext(1);   i2=iwdtmin-1
                  j1=isvmin;       j2=isvmax-1
                  k1=icldext(1,4); k2=icldext(1,4)
                elseif (ifac.eq.5 .or. ifac.eq.6) then
                  i1=iparext(1);   i2=iwdtmin-1
                  j1=icldext(1,4); j2=icldext(1,4)
                  k1=iswmin;       k2=iswmax-1
                endif
                write (t12,'(2a,i3.3,a,i3.3)')
     &            facstr,'sld',NSB(ICFD)+1,'_s',isur
                call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &            SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
              endif

              if ((iparext(2)-1).gt.iwdtmax) then
                if (ifac.eq.1 .or. ifac.eq.2) then
                  i1=isumin;       i2=isumax-1
                  j1=iwdtmax+1;    j2=iparext(2)-1
                  k1=icldext(1,4); k2=icldext(1,4)
                elseif (ifac.eq.3 .or. ifac.eq.4) then
                  i1=iwdtmax+1;    i2=iparext(2)-1
                  j1=isvmin;       j2=isvmax-1
                  k1=icldext(1,4); k2=icldext(1,4)
                elseif (ifac.eq.5 .or. ifac.eq.6) then
                  i1=iwdtmax+1;    i2=iparext(2)-1
                  j1=icldext(1,4); j2=icldext(1,4)
                  k1=iswmin;       k2=iswmax-1
                endif
                write (t12,'(2a,i3.3,a,i3.3)')
     &            facstr,'sld',NSB(ICFD)+1,'_s',isur
                call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &            SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
              endif

C Any space left over that needs filling with a solid bc?
              if (irows.gt.1) then
                if (ifac.eq.1 .or. ifac.eq.2) then
                  i1=isumin;         i2=isumax-1
                  j1=iparext(1);     j2=iparext(2)-1
                  k1=icldext(1,4)+1; k2=iparext(4)-1
                elseif (ifac.eq.3 .or. ifac.eq.4) then
                  i1=iparext(1);     i2=iparext(2)-1
                  j1=isvmin;         j2=isvmax-1
                  k1=icldext(1,4)+1; k2=iparext(4)-1
                elseif (ifac.eq.5 .or. ifac.eq.6) then
                  i1=iparext(1);     i2=iparext(2)-1
                  j1=icldext(1,4)+1; j2=iparext(4)-1
                  k1=iswmin;         k2=iswmax-1
                endif
                write (t12,'(2a,i3.3,a,i3.3)')
     &            facstr,'sld',NSB(ICFD)+1,'_s',isur
                call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &            SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
              endif

              itop=iparext(4)-1-irows
            else
              itop=iparext(4)-1
            endif

C Is there space on the bottom?
            irows=icldext(1,3)-iparext(3)
            if (irows.gt.0) then
              if (ok) then ! opening already created?
                if (ifac.eq.1 .or. ifac.eq.2) then
                  i1=isumin;     i2=isumax-1
                  j1=iparext(1); j2=iparext(2)-1
                  k1=iparext(3); k2=icldext(1,3)-1
                elseif (ifac.eq.3 .or. ifac.eq.4) then
                  i1=iparext(1); i2=iparext(2)-1
                  j1=isvmin;     j2=isvmax-1
                  k1=iparext(3); k2=icldext(1,3)-1
                elseif (ifac.eq.5 .or. ifac.eq.6) then
                  i1=iparext(1); i2=iparext(2)-1
                  j1=iparext(3); j2=icldext(1,3)-1
                  k1=iswmin;     k2=iswmax-1
                endif
                write (t12,'(2a,i3.3,a,i3.3)')
     &            facstr,'sld',NSB(ICFD)+1,'_s',isur
                call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &            SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
              else

C Find cells to put the trickle vent in the middle.
                imidcol=icldext(icld,1)+icols/2
                if (mod(iwdt,2).eq.0) then
                  iwdtmin=imidcol-iwdt/2
                  iwdtmax=imidcol+iwdt/2-1
                else
                  iwdtmin=imidcol-iwdt/2
                  iwdtmax=imidcol+iwdt/2
                endif

                if (ifac.eq.1 .or. ifac.eq.2) then
                  i1=isumin;         i2=isumax-1
                  j1=iwdtmin;        j2=iwdtmax
                  k1=icldext(1,3)-1; k2=icldext(1,3)-1
                elseif (ifac.eq.3 .or. ifac.eq.4) then
                  i1=iwdtmin;        i2=iwdtmax
                  j1=isvmin;         j2=isvmax-1
                  k1=icldext(1,3)-1; k2=icldext(1,3)-1
                elseif (ifac.eq.5 .or. ifac.eq.6) then
                  i1=iwdtmin;        i2=iwdtmax
                  j1=icldext(1,3)-1; j2=icldext(1,3)-1
                  k1=iswmin;         k2=iswmax-1
                endif
                write (t12,'(2a,i3.3,a,i3.3)')
     &            facstr,'opn',NOPEN(ICFD)+1,'_s',isur
                call MKCFDBC(show,11,ifac,i1,i2,j1,j2,k1,k2,t12,
     &            SSNAME(icnn),dum,1,20.0,0.0,0.01,0.0,0.0,IER)
                ok=.true.

                if (iparext(1).lt.iwdtmin) then
                  if (ifac.eq.1 .or. ifac.eq.2) then
                    i1=isumin;         i2=isumax-1
                    j1=iparext(1);     j2=iwdtmin-1
                    k1=icldext(1,3)-1; k2=icldext(1,3)-1
                  elseif (ifac.eq.3 .or. ifac.eq.4) then
                    i1=iparext(1);     i2=iwdtmin-1
                    j1=isvmin;         j2=isvmax-1
                    k1=icldext(1,3)-1; k2=icldext(1,3)-1
                  elseif (ifac.eq.5 .or. ifac.eq.6) then
                    i1=iparext(1);     i2=iwdtmin-1
                    j1=icldext(1,3)-1; j2=icldext(1,3)-1
                    k1=iswmin;         k2=iswmax-1
                  endif
                  write (t12,'(2a,i3.3,a,i3.3)')
     &              facstr,'sld',NSB(ICFD)+1,'_s',isur
                  call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &              SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
                endif

                if ((iparext(2)-1).gt.iwdtmax) then
                  if (ifac.eq.1 .or. ifac.eq.2) then
                    i1=isumin;         i2=isumax-1
                    j1=iwdtmax+1;      j2=iparext(2)-1
                    k1=icldext(1,3)-1; k2=icldext(1,3)-1
                  elseif (ifac.eq.3 .or. ifac.eq.4) then
                    i1=iwdtmax+1;      i2=iparext(2)-1
                    j1=isvmin;         j2=isvmax-1
                    k1=icldext(1,3)-1; k2=icldext(1,3)-1
                  elseif (ifac.eq.5 .or. ifac.eq.6) then
                    i1=iwdtmax+1;      i2=iparext(2)-1
                    j1=icldext(1,3)-1; j2=icldext(1,3)-1
                    k1=iswmin;         k2=iswmax-1
                  endif
                  write (t12,'(2a,i3.3,a,i3.3)')
     &              facstr,'sld',NSB(ICFD)+1,'_s',isur
                  call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &              SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
                endif

C Any space left over that needs filling with a solid bc?
                if (irows.gt.1) then
                  if (ifac.eq.1 .or. ifac.eq.2) then
                    i1=isumin;     i2=isumax-1
                    j1=iparext(1); j2=iparext(2)-1
                    k1=iparext(3); k2=icldext(1,3)-2
                  elseif (ifac.eq.3 .or. ifac.eq.4) then
                    i1=iparext(1); i2=iparext(2)-1
                    j1=isvmin;     j2=isvmax-1
                    k1=iparext(3); k2=icldext(1,3)-2
                  elseif (ifac.eq.5 .or. ifac.eq.6) then
                    i1=iparext(1); i2=iparext(2)-1
                    j1=iparext(3); j2=icldext(1,3)-2
                    k1=iswmin;     k2=iswmax-1
                  endif
                  write (t12,'(2a,i3.3,a,i3.3)')
     &              facstr,'sld',NSB(ICFD)+1,'_s',isur
                  call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &              SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
                endif
              endif
              ibot=iparext(3)+irows
            else
              ibot=iparext(3)
            endif

            irows=iparext(4)-iparext(3)

C Is there space on the left?
            icols=icldext(1,1)-iparext(1)
            if (icols.gt.0) then
              if (ok) then
                if (ifac.eq.1 .or. ifac.eq.2) then
                  i1=isumin;     i2=isumax-1
                  j1=iparext(1); j2=icldext(1,1)-1
                  k1=ibot;       k2=itop
                elseif (ifac.eq.3 .or. ifac.eq.4) then
                  i1=iparext(1); i2=icldext(1,1)-1
                  j1=isvmin;     j2=isvmax-1
                  k1=ibot;       k2=itop
                elseif (ifac.eq.5 .or. ifac.eq.6) then
                  i1=iparext(1); i2=icldext(1,1)-1
                  j1=ibot;       j2=itop
                  k1=iswmin;     k2=iswmax-1
                endif
                write (t12,'(2a,i3.3,a,i3.3)')
     &            facstr,'sld',NSB(ICFD)+1,'_s',isur
                call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &            SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
              else

C Find cells to put the trickle vent in the middle.
                imidrow=icldext(icld,3)+irows/2
                if (mod(iwdt,2).eq.0) then
                  iwdtmin=imidrow-iwdt/2
                  iwdtmax=imidrow+iwdt/2-1
                else
                  iwdtmin=imidrow-iwdt/2
                  iwdtmax=imidrow+iwdt/2
                endif

                if (ifac.eq.1 .or. ifac.eq.2) then
                  i1=isumin;         i2=isumax-1
                  j1=icldext(1,1)-1; j2=icldext(1,1)-1
                  k1=iwdtmin;        k2=iwdtmax
                elseif (ifac.eq.3 .or. ifac.eq.4) then
                  i1=icldext(1,1)-1; i2=icldext(1,1)-1
                  j1=isvmin;         j2=isvmax-1
                  k1=iwdtmin;        k2=iwdtmax
                elseif (ifac.eq.5 .or. ifac.eq.6) then
                  i1=icldext(1,1)-1; i2=icldext(1,1)-1
                  j1=iwdtmin;        j2=iwdtmax
                  k1=iswmin;         k2=iswmax-1
                endif
                write (t12,'(2a,i3.3,a,i3.3)')
     &            facstr,'opn',NOPEN(ICFD)+1,'_s',isur
                call MKCFDBC(show,11,ifac,i1,i2,j1,j2,k1,k2,t12,
     &            SSNAME(icnn),dum,1,20.0,0.0,0.01,0.0,0.0,IER)
                ok=.true.

                if (ibot.lt.iwdtmin) then
                  if (ifac.eq.1 .or. ifac.eq.2) then
                    i1=isumin;         i2=isumax-1
                    j1=icldext(1,1)-1; j2=icldext(1,1)-1
                    k1=ibot;           k2=iwdtmin-1
                  elseif (ifac.eq.3 .or. ifac.eq.4) then
                    i1=icldext(1,1)-1; i2=icldext(1,1)-1
                    j1=isvmin;         j2=isvmax-1
                    k1=ibot;           k2=iwdtmin-1
                  elseif (ifac.eq.5 .or. ifac.eq.6) then
                    i1=icldext(1,1)-1; i2=icldext(1,1)-1
                    j1=ibot;           j2=iwdtmin-1
                    k1=iswmin;         k2=iswmax-1
                  endif
                  write (t12,'(2a,i3.3,a,i3.3)')
     &              facstr,'sld',NSB(ICFD)+1,'_s',isur
                  call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &              SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
                endif

                if (itop.gt.iwdtmax) then
                  if (ifac.eq.1 .or. ifac.eq.2) then
                    i1=isumin;         i2=isumax-1
                    j1=icldext(1,1)-1; j2=icldext(1,1)-1
                    k1=iwdtmax+1;      k2=itop
                  elseif (ifac.eq.3 .or. ifac.eq.4) then
                    i1=icldext(1,1)-1; i2=icldext(1,1)-1
                    j1=isvmin;         j2=isvmax-1
                    k1=iwdtmax+1;      k2=itop
                  elseif (ifac.eq.5 .or. ifac.eq.6) then
                    i1=icldext(1,1)-1; i2=icldext(1,1)-1
                    j1=iwdtmax+1;      j2=itop
                    k1=iswmin;         k2=iswmax-1
                  endif
                  write (t12,'(2a,i3.3,a,i3.3)')
     &              facstr,'sld',NSB(ICFD)+1,'_s',isur
                  call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &              SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
                endif

C Any space left over that needs filling with a solid bc?
                if (icols.gt.1) then
                  if (ifac.eq.1 .or. ifac.eq.2) then
                    i1=isumin;     i2=isumax-1
                    j1=iparext(1); j2=icldext(1,1)-2
                    k1=ibot;       k2=itop
                  elseif (ifac.eq.3 .or. ifac.eq.4) then
                    i1=iparext(1); i2=icldext(1,1)-2
                    j1=isvmin;     j2=isvmax-1
                    k1=ibot;       k2=itop
                  elseif (ifac.eq.5 .or. ifac.eq.6) then
                    i1=iparext(1); i2=icldext(1,1)-2
                    j1=ibot;       j2=itop
                    k1=iswmin;     k2=iswmax-1
                  endif
                  write (t12,'(2a,i3.3,a,i3.3)')
     &              facstr,'sld',NSB(ICFD)+1,'_s',isur
                  call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &              SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
                endif
              endif
            endif

C Is there space on the right?
            icols=iparext(2)-icldext(1,2)
            if (icols.gt.0) then
              if (ok) then
                if (ifac.eq.1 .or. ifac.eq.2) then
                  i1=isumin;       i2=isumax-1
                  j1=icldext(1,2); j2=iparext(2)-1
                  k1=ibot;         k2=itop
                elseif (ifac.eq.3 .or. ifac.eq.4) then
                  i1=icldext(1,2); i2=iparext(2)-1
                  j1=isvmin;       j2=isvmax-1
                  k1=ibot;         k2=itop
                elseif (ifac.eq.5 .or. ifac.eq.6) then
                  i1=icldext(1,2); i2=iparext(2)-1
                  j1=ibot;         j2=itop
                  k1=iswmin;       k2=iswmax-1
                endif
                write (t12,'(2a,i3.3,a,i3.3)')
     &            facstr,'sld',NSB(ICFD)+1,'_s',isur
                call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &            SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
              else

C Find cells to put the trickle vent in the middle.
                imidrow=icldext(icld,3)+irows/2
                if (mod(iwdt,2).eq.0) then
                  iwdtmin=imidrow-iwdt/2
                  iwdtmax=imidrow+iwdt/2-1
                else
                  iwdtmin=imidrow-iwdt/2
                  iwdtmax=imidrow+iwdt/2
                endif

                if (ifac.eq.1 .or. ifac.eq.2) then
                  i1=isumin;       i2=isumax-1
                  j1=icldext(1,2); j2=icldext(1,2)
                  k1=iwdtmin;      k2=iwdtmax
                elseif (ifac.eq.3 .or. ifac.eq.4) then
                  i1=icldext(1,2); i2=icldext(1,2)
                  j1=isvmin;       j2=isvmax-1
                  k1=iwdtmin;      k2=iwdtmax
                elseif (ifac.eq.5 .or. ifac.eq.6) then
                  i1=icldext(1,2); i2=icldext(1,2)
                  j1=iwdtmin;      j2=iwdtmax
                  k1=iswmin;       k2=iswmax-1
                endif
                write (t12,'(2a,i3.3,a,i3.3)')
     &            facstr,'opn',NOPEN(ICFD)+1,'_s',isur
                call MKCFDBC(show,11,ifac,i1,i2,j1,j2,k1,k2,t12,
     &            SSNAME(icnn),dum,1,20.0,0.0,0.01,0.0,0.0,IER)
                ok=.true.

                if (ibot.lt.iwdtmin) then
                  if (ifac.eq.1 .or. ifac.eq.2) then
                    i1=isumin;       i2=isumax-1
                    j1=icldext(1,2); j2=icldext(1,2)
                    k1=ibot;         k2=iwdtmin-1
                  elseif (ifac.eq.3 .or. ifac.eq.4) then
                    i1=icldext(1,2); i2=icldext(1,2)
                    j1=isvmin;       j2=isvmax-1
                    k1=ibot;         k2=iwdtmin-1
                  elseif (ifac.eq.5 .or. ifac.eq.6) then
                    i1=icldext(1,2); i2=icldext(1,2)
                    j1=ibot;         j2=iwdtmin-1
                    k1=iswmin;       k2=iswmax-1
                  endif
                  write (t12,'(2a,i3.3,a,i3.3)')
     &              facstr,'sld',NSB(ICFD)+1,'_s',isur
                  call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &              SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
                endif

                if (itop.gt.iwdtmax) then
                  if (ifac.eq.1 .or. ifac.eq.2) then
                    i1=isumin;       i2=isumax-1
                    j1=icldext(1,2); j2=icldext(1,2)
                    k1=iwdtmax+1;    k2=itop
                  elseif (ifac.eq.3 .or. ifac.eq.4) then
                    i1=icldext(1,2); i2=icldext(1,2)
                    j1=isvmin;       j2=isvmax-1
                    k1=iwdtmax+1;    k2=itop
                  elseif (ifac.eq.5 .or. ifac.eq.6) then
                    i1=icldext(1,2); i2=icldext(1,2)
                    j1=iwdtmax+1;    j2=itop
                    k1=iswmin;       k2=iswmax-1
                  endif
                  write (t12,'(2a,i3.3,a,i3.3)')
     &              facstr,'sld',NSB(ICFD)+1,'_s',isur
                  call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &              SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
                endif

C Any space left over that needs filling with a solid bc?
                if (icols.gt.1) then
                  if (ifac.eq.1 .or. ifac.eq.2) then
                    i1=isumin;         i2=isumax-1
                    j1=icldext(1,2)+1; j2=iparext(2)-1
                    k1=ibot;           k2=itop
                  elseif (ifac.eq.3 .or. ifac.eq.4) then
                    j1=icldext(1,2)+1; j2=iparext(2)-1
                    j1=isvmin;         j2=isvmax-1
                    k1=ibot;           k2=itop
                  elseif (ifac.eq.5 .or. ifac.eq.6) then
                    j1=icldext(1,2)+1; j2=iparext(2)-1
                    j1=ibot;           j2=itop
                    k1=iswmin;         k2=iswmax-1
                  endif
                  write (t12,'(2a,i3.3,a,i3.3)')
     &              facstr,'sld',NSB(ICFD)+1,'_s',isur
                  call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &              SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
                endif
              endif
            endif
          endif
        endif
      else

C This surface is not a parent; create BC(s).
        if (dothisur) then

          write(outs,'(a,i3.3,2a)')' Creating BC(s) for surface ',
     &      isur,SSNAME(icnn)(1:lnblnk(SSNAME(icnn))),' ...'
          ISHSB=0; ISHAO=0; ISHBLK=0; IFACES=1; ISHGEO=-isur
          call redraw(IER)
          call pausems(500)
          ISHGEO=-9999

C "frame with a crack" usage tags (ibctyp 9) might be used to put a
C crack in a solid surface - in this case, switch to a default crack
C treatment (ibctyp 4).
          if (ibctyp.eq.9) ibctyp=4

C Single solid BC for the whole surface.
          if (ibctyp.eq.1) then
            if (ifac.eq.1 .or. ifac.eq.2) then
              i1=isumin;     i2=isumax-1
              j1=iparext(1); j2=iparext(2)-1
              k1=iparext(3); k2=iparext(4)-1
            elseif (ifac.eq.3 .or. ifac.eq.4) then
              i1=iparext(1); i2=iparext(2)-1
              j1=isvmin;     j2=isvmax-1
              k1=iparext(3); k2=iparext(4)-1
            elseif (ifac.eq.5 .or. ifac.eq.6) then
              i1=iparext(1); i2=iparext(2)-1
              j1=iparext(3); j2=iparext(4)-1
              k1=iswmin;     k2=iswmax-1
            endif
            write (t12,'(2a,i3.3,a,i3.3)')
     &        facstr,'sld',NSB(ICFD)+1,'_s',isur
            call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &        SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)

C Single opening BC for the whole surface.
          elseif (ibctyp.eq.2) then
            if (ifac.eq.1 .or. ifac.eq.2) then
              i1=isumin;     i2=isumax-1
              j1=iparext(1); j2=iparext(2)-1
              k1=iparext(3); k2=iparext(4)-1
            elseif (ifac.eq.3 .or. ifac.eq.4) then
              i1=iparext(1); i2=iparext(2)-1
              j1=isvmin;     j2=isvmax-1
              k1=iparext(3); k2=iparext(4)-1
            elseif (ifac.eq.5 .or. ifac.eq.6) then
              i1=iparext(1); i2=iparext(2)-1
              j1=iparext(3); j2=iparext(4)-1
              k1=iswmin;     k2=iswmax-1
            endif
            write (t12,'(2a,i3.3,a,i3.3)')
     &        facstr,'opn',NOPEN(ICFD)+1,'_s',isur
            call MKCFDBC(show,11,ifac,i1,i2,j1,j2,k1,k2,t12,
     &        SSNAME(icnn),dum,1,20.0,0.0,0.01,0.0,0.0,IER)

C Opening BC for bottom row of cells, solid BC for the rest.
          elseif (ibctyp.eq.3) then
            if (ifac.eq.1 .or. ifac.eq.2) then
              i1=isumin;     i2=isumax-1
              j1=iparext(1); j2=iparext(2)-1
              k1=iparext(3); k2=iparext(3)
            elseif (ifac.eq.3 .or. ifac.eq.4) then
              i1=iparext(1); i2=iparext(2)-1
              j1=isvmin;     j2=isvmax-1
              k1=iparext(3); k2=iparext(3)
            elseif (ifac.eq.5 .or. ifac.eq.6) then
              i1=iparext(1); i2=iparext(2)-1
              j1=iparext(3); j2=iparext(3)
              k1=iswmin;     k2=iswmax-1
            endif
            write (t12,'(2a,i3.3,a,i3.3)')
     &        facstr,'opn',NOPEN(ICFD)+1,'_s',isur
            call MKCFDBC(show,11,ifac,i1,i2,j1,j2,k1,k2,t12,
     &        SSNAME(icnn),dum,1,20.0,0.0,0.01,0.0,0.0,IER)

C Is there any more of the surface for a solid BC?
            if (iparext(4)-1.gt.iparext(3)) then
              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;       i2=isumax-1
                j1=iparext(1);   j2=iparext(2)-1
                k1=iparext(3)+1; k2=iparext(4)-1
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=iparext(1);   i2=iparext(2)-1
                j1=isvmin;       j2=isvmax-1
                k1=iparext(3)+1; k2=iparext(4)-1
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=iparext(1);   i2=iparext(2)-1
                j1=iparext(3)+1; j2=iparext(4)-1
                k1=iswmin;       k2=iswmax-1
              endif
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'sld',NSB(ICFD)+1,'_s',isur
              call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)

            endif

C Opening BC for middle row of cells, solid BCs for the rest.
          elseif (ibctyp.eq.4) then

C Find middle row.
            irows=iparext(4)-iparext(3)
            imidrow=iparext(3)+irows/2

            if (ifac.eq.1 .or. ifac.eq.2) then
              i1=isumin;          i2=isumax-1
              j1=icldext(icld,1); j2=icldext(icld,2)-1
              k1=imidrow;         k2=imidrow
            elseif (ifac.eq.3 .or. ifac.eq.4) then
              i1=icldext(icld,1); i2=icldext(icld,2)-1
              j1=isvmin;          j2=isvmax-1
              k1=imidrow;         k2=imidrow
            elseif (ifac.eq.5 .or. ifac.eq.6) then
              i1=icldext(icld,1); i2=icldext(icld,2)-1
              j1=imidrow;         j2=imidrow
              k1=iswmin;          k2=iswmax-1
            endif
            write (t12,'(2a,i3.3,a,i3.3)')
     &        facstr,'opn',NOPEN(ICFD)+1,'_s',isur
            call MKCFDBC(show,11,ifac,i1,i2,j1,j2,k1,k2,t12,
     &        SSNAME(icnn),dum,1,20.0,0.0,0.01,0.0,0.0,IER)

C Is there any more of the surface below for a solid BC?
            if (iparext(3).lt.imidrow) then
              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;     i2=isumax-1
                j1=iparext(1); j2=iparext(2)-1
                k1=iparext(3); k2=imidrow-1
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=iparext(1); i2=iparext(2)-1
                j1=isvmin;     j2=isvmax-1
                k1=iparext(3); k2=imidrow-1
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=iparext(1); i2=iparext(2)-1
                j1=iparext(3); j2=imidrow-1
                k1=iswmin;     k2=iswmax-1
              endif
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'sld',NSB(ICFD)+1,'_s',isur
              call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
            endif

C Is there any more of the surface above for a solid BC?
            if ((iparext(4)-1).gt.imidrow) then
              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;     i2=isumax-1
                j1=iparext(1); j2=iparext(2)-1
                k1=imidrow+1;  k2=iparext(4)-1
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=iparext(1); i2=iparext(2)-1
                j1=isvmin;     j2=isvmax-1
                k1=imidrow+1;  k2=iparext(4)-1
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=iparext(1); i2=iparext(2)-1
                j1=imidrow+1;  j2=iparext(4)-1
                k1=iswmin;     k2=iswmax-1
              endif
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'sld',NSB(ICFD)+1,'_s',isur
              call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
            endif

C Opening BC for 25% of height, at bottom edge, solid BC for the rest.
          elseif (ibctyp.eq.5) then

C Find 25% of height, rounded, at least 1 row.
            rrows=real(iparext(4)-iparext(3))
            iqtr=nint(rrows/4.0)
            if (iqtr.eq.0) iqtr=1

            if (ifac.eq.1 .or. ifac.eq.2) then
              i1=isumin;     i2=isumax-1
              j1=iparext(1); j2=iparext(2)-1
              k1=iparext(3); k2=iparext(3)+iqtr-1
            elseif (ifac.eq.3 .or. ifac.eq.4) then
              i1=iparext(1); i2=iparext(2)-1
              j1=isvmin;     j2=isvmax-1
              k1=iparext(3); k2=iparext(3)+iqtr-1
            elseif (ifac.eq.5 .or. ifac.eq.6) then
              i1=iparext(1); i2=iparext(2)-1
              j1=iparext(3); j2=iparext(3)+iqtr-1
              k1=iswmin;     k2=iswmax-1
            endif
            write (t12,'(2a,i3.3,a,i3.3)')
     &        facstr,'opn',NOPEN(ICFD)+1,'_s',isur
            call MKCFDBC(show,11,ifac,i1,i2,j1,j2,k1,k2,t12,
     &        SSNAME(icnn),dum,1,20.0,0.0,0.01,0.0,0.0,IER)

C Is there any more of the surface for a solid BC?
            if ((iparext(4)-1).gt.(iparext(3)+iqtr)) then
              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;          i2=isumax-1
                j1=iparext(1);      j2=iparext(2)-1
                k1=iparext(3)+iqtr; k2=iparext(4)-1
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=iparext(1);      i2=iparext(2)-1
                j1=isvmin;          j2=isvmax-1
                k1=iparext(3)+iqtr; k2=iparext(4)-1
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=iparext(1);      i2=iparext(2)-1
                j1=iparext(3)+iqtr; j2=iparext(4)-1
                k1=iswmin;          k2=iswmax-1
              endif
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'sld',NSB(ICFD)+1,'_s',isur
              call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
            endif

C Opening BC for 25% of width, at vertical centre line, solid BCs for
C rest.
          elseif (ibctyp.eq.6) then

C Find 25% of width, rounded, at least 1 column, on vertical centre
C line.
            icols=iparext(2)-iparext(1)
            rcols=real(icols)
            iqtr=nint(rcols/4.0)
            if (iqtr.eq.0) iqtr=1
            imidcol=iparext(1)+icols/2
            if (mod(iqtr,2).eq.0) then
              imidqtrmin=imidcol-iqtr/2
              imidqtrmax=imidcol+iqtr/2-1
            else
              imidqtrmin=imidcol-iqtr/2
              imidqtrmax=imidcol+iqtr/2
            endif

            if (ifac.eq.1 .or. ifac.eq.2) then
              i1=isumin;     i2=isumax-1
              j1=imidqtrmin; j2=imidqtrmax
              k1=iparext(3); k2=iparext(4)-1
            elseif (ifac.eq.3 .or. ifac.eq.4) then
              i1=imidqtrmin; i2=imidqtrmax
              j1=isvmin;     j2=isvmax-1
              k1=iparext(3); k2=iparext(4)-1
            elseif (ifac.eq.5 .or. ifac.eq.6) then
              i1=imidqtrmin; i2=imidqtrmax
              j1=iparext(3); j2=iparext(4)-1
              k1=iswmin;     k2=iswmax-1
            endif
            write (t12,'(2a,i3.3,a,i3.3)')
     &        facstr,'opn',NOPEN(ICFD)+1,'_s',isur
            call MKCFDBC(show,11,ifac,i1,i2,j1,j2,k1,k2,t12,
     &        SSNAME(icnn),dum,1,20.0,0.0,0.01,0.0,0.0,IER)

C Is there any more of the surface to the left for a solid BC?
            if ((iparext(1)).lt.imidqtrmin) then
              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;     i2=isumax-1
                j1=iparext(1); j2=imidqtrmin-1
                k1=iparext(3); k2=iparext(4)-1
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=iparext(1); i2=imidqtrmin-1
                j1=isvmin;     j2=isvmax-1
                k1=iparext(3); k2=iparext(4)-1
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=iparext(1); i2=imidqtrmin-1
                j1=iparext(3); j2=iparext(4)-1
                k1=iswmin;     k2=iswmax-1
              endif
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'sld',NSB(ICFD)+1,'_s',isur
              call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
            endif

C Is there any more of the surface to the right for a solid BC?
            if ((iparext(2)-1).gt.imidqtrmax) then
              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;       i2=isumax-1
                j1=imidqtrmax+1; j2=iparext(2)-1
                k1=iparext(3);   k2=iparext(4)-1
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=imidqtrmax+1; i2=iparext(2)-1
                j1=isvmin;       j2=isvmax-1
                k1=iparext(4);   k2=iparext(4)-1
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=imidqtrmax+1; i2=iparext(2)-1
                j1=iparext(3);   j2=iparext(4)-1
                k1=iswmin;       k2=iswmax-1
              endif
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'sld',NSB(ICFD)+1,'_s',isur
              call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
            endif

C Opening BCs for 12.5% of height, at top and bottom edges, solid BC for
C rest.
          elseif (ibctyp.eq.7) then

C Find 12.5% of height, rounded, at least 1 row.
            irows=iparext(4)-iparext(3)
            rrows=real(irows)
            i8th=nint(rrows/8.0)
            if (i8th.eq.0) i8th=1

C Is there enough space for both openings?
            ok=.true.
   90       if (irows.lt.(i8th*2)) then
              if (i8th.gt.1) then
                i8th=i8th-1
              else
                write(outs,'(3a)')' Not enough space for openings at',
     &            ' top and bottom of this sash window; only one',
     &            ' opening will be created.'
                call EDISP(IUOUT,outs)
                ok=.false.
              endif
              goto 90
            endif

C Opening at bottom edge.
            if (ifac.eq.1 .or. ifac.eq.2) then
              i1=isumin;     i2=isumax-1
              j1=iparext(1); j2=iparext(2)-1
              k1=iparext(3); k2=iparext(3)+i8th-1
            elseif (ifac.eq.3 .or. ifac.eq.4) then
              i1=iparext(1); i2=iparext(2)-1
              j1=isvmin;     j2=isvmax-1
              k1=iparext(3); k2=iparext(3)+i8th-1
            elseif (ifac.eq.5 .or. ifac.eq.6) then
              i1=iparext(1); i2=iparext(2)-1
              j1=iparext(3); j2=iparext(3)+i8th-1
              k1=iswmin;     k2=iswmax-1
            endif
            write (t12,'(2a,i3.3,a,i3.3)')
     &        facstr,'opn',NOPEN(ICFD)+1,'_s',isur
            call MKCFDBC(show,11,ifac,i1,i2,j1,j2,k1,k2,t12,
     &        SSNAME(icnn),dum,1,20.0,0.0,0.01,0.0,0.0,IER)

C Is there room for the top opening?
            if (ok) then
              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;          i2=isumax-1
                j1=iparext(1);      j2=iparext(2)-1
                k1=iparext(4)-i8th; k2=iparext(4)-1
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=iparext(1);      i2=iparext(2)-1
                j1=isvmin;          j2=isvmax-1
                k1=iparext(4)-i8th; k2=iparext(4)-1
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=iparext(1);      i2=iparext(2)-1
                j1=iparext(4)-i8th; j2=iparext(4)-1
                k1=iswmin;          k2=iswmax-1
              endif
              itop=iparext(4)-i8th-1
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'opn',NOPEN(ICFD)+1,'_s',isur
              call MKCFDBC(show,11,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icnn),dum,1,20.0,0.0,0.01,0.0,0.0,IER)
            else
              itop=iparext(4)-1
            endif

C Is there any more of the surface for a solid BC?
            if ((iparext(3)+i8th).le.itop) then
              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;          i2=isumax-1
                j1=iparext(1);      j2=iparext(2)-1
                k1=iparext(3)+i8th; k2=itop
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=iparext(1);      i2=iparext(2)-1
                j1=isvmin;          j2=isvmax-1
                k1=iparext(3)+i8th; k2=itop
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=iparext(1);      i2=iparext(2)-1
                j1=iparext(3)+i8th; j2=itop
                k1=iswmin;          k2=iswmax-1
              endif
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'sld',NSB(ICFD)+1,'_s',isur
              call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
            endif

C Opening BC for 75% of height and width, at bottom edge, solid BCs for
C rest.            
          elseif (ibctyp.eq.8) then

C Find 75% of height, rounded, at least 1 row, at bottom edge.
            irows=iparext(4)-iparext(3)
            rrows=real(irows)
            i3qtrhgt=nint(3.0*rrows/4.0)
            if (i3qtrhgt.eq.0) i3qtrhgt=1
            ibtm3qtrhgtmin=iparext(3)
            ibtm3qtrhgtmin=iparext(3)+i3qtrhgt-1

C Find 75% of width, rounded, at least 1 column, on vertical centre
C line.
            icols=iparext(2)-iparext(1)
            rcols=real(icols)
            i3qtrwdt=nint(3.0*rcols/4.0)
            if (i3qtrwdt.eq.0) i3qtrwdt=1
            imidcol=iparext(1)+icols/2
            if (mod(i3qtrwdt,2).eq.0) then
              imid3qtrwdtmin=imidcol-i3qtrwdt/2
              imid3qtrwdtmax=imidcol+i3qtrwdt/2-1
            else
              imid3qtrwdtmin=imidcol-i3qtrwdt/2
              imid3qtrwdtmax=imidcol+i3qtrwdt/2
            endif

C Opening BC.
            if (ifac.eq.1 .or. ifac.eq.2) then
              i1=isumin;         i2=isumax-1
              j1=imid3qtrwdtmin; j2=imid3qtrwdtmax
              k1=ibtm3qtrhgtmin; k2=ibtm3qtrhgtmax
            elseif (ifac.eq.3 .or. ifac.eq.4) then
              i1=imid3qtrwdtmin; i2=imid3qtrwdtmax
              j1=isvmin;         j2=isvmax-1
              k1=ibtm3qtrhgtmin; k2=ibtm3qtrhgtmax
            elseif (ifac.eq.5 .or. ifac.eq.6) then
              i1=imid3qtrwdtmin; i2=imid3qtrwdtmax
              j1=ibtm3qtrhgtmin; j2=ibtm3qtrhgtmax
              k1=iswmin;         k2=iswmax-1
            endif
            write (t12,'(2a,i3.3,a,i3.3)')
     &        facstr,'opn',NOPEN(ICFD)+1,'_s',isur
            call MKCFDBC(show,11,ifac,i1,i2,j1,j2,k1,k2,t12,
     &        SSNAME(icnn),dum,1,20.0,0.0,0.01,0.0,0.0,IER)

C Is there room left at the top?
            if (ibtm3qtrhgtmax.lt.(iparext(4)-1)) then
              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;           i2=isumax-1
                j1=iparext(1);       j2=iparext(2)-1
                k1=ibtm3qtrhgtmax+1; k2=iparext(4)-1
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=iparext(1);       i2=iparext(2)-1
                j1=isvmin;           j2=isvmax-1
                k1=ibtm3qtrhgtmax+1; k2=iparext(4)-1
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=iparext(1);       i2=iparext(2)-1
                j1=ibtm3qtrhgtmax+1; j2=iparext(4)-1
                k1=iswmin;           k2=iswmax-1
              endif
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'sld',NSB(ICFD)+1,'_s',isur
              call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
              ok=.true.
            else
              ok=.false.
            endif

C Is there room at the left?
            if (imid3qtrwdtmin.gt.iparext(1)) then
              if (ok) then
                itop=iparext(4)-2
              else
                itop=iparext(4)-1
              endif
              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;     i2=isumax-1
                j1=iparext(1); j2=imid3qtrwdtmin-1
                k1=iparext(3); k2=itop
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=iparext(1); i2=imid3qtrwdtmin-1
                j1=isvmin;     j2=isvmax-1
                k1=iparext(3); k2=itop
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=iparext(1); i2=imid3qtrwdtmin-1
                j1=iparext(3); j2=itop
                k1=iswmin;     k2=iswmax-1
              endif
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'sld',NSB(ICFD)+1,'_s',isur
              call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
            endif

C Is there room at the right?
            if (imid3qtrwdtmax.lt.(iparext(2)-1)) then
              if (ok) then
                itop=iparext(4)-2
              else
                itop=iparext(4)-1
              endif
              if (ifac.eq.1 .or. ifac.eq.2) then
                i1=isumin;           i2=isumax-1
                j1=imid3qtrwdtmax+1; j2=iparext(2)-1
                k1=iparext(3);       k2=itop
              elseif (ifac.eq.3 .or. ifac.eq.4) then
                i1=imid3qtrwdtmax+1; i2=iparext(2)-1
                j1=isvmin;           j2=isvmax-1
                k1=iparext(3);       k2=itop
              elseif (ifac.eq.5 .or. ifac.eq.6) then
                i1=imid3qtrwdtmax+1; i2=iparext(2)-1
                j1=iparext(3);       j2=itop
                k1=iswmin;           k2=iswmax-1
              endif
              write (t12,'(2a,i3.3,a,i3.3)')
     &          facstr,'sld',NSB(ICFD)+1,'_s',isur
              call MKCFDBC(show,4,ifac,i1,i2,j1,j2,k1,k2,t12,
     &          SSNAME(icnn),dum,1,20.0,0.0,0.0,0.0,0.0,IER)
            endif

C 9 and 10 are for frames.
          elseif (ibctyp.eq.9 .or. ibctyp.eq.10) then
            write(6,'(a)')
     &' This surface is a frame, but has no children. Check usage tags.'
            call EDISP(IUOUT,outs)
          endif
        endif
      endif

  999 if (ISTAT.gt.0 .and. ibctyp.gt.1) doafnasc=.false.

      RETURN
      END












C ************************** EDAIRO **************************
C EDAIRO controls the editing of CFD air flow opening boundary variables.
C IER=0 indicates no error.

      SUBROUTINE EDAIRO(IAO,IER)
C#include "building.h"
#include "cfd.h"

      COMMON/ICFNOD/ICFD,ICP
      COMMON/SPAD/MMOD,LIMIT,LIMTTY

      common/cfdconf/ICFBLD(MNZ),ICFMFS(MNZ)
      common/KEYVOLS/NVOL(MNZ),IVOLF(MNVLS,MNZ),IVCELLS(MNVLS,MNZ,2),
     &               JVCELLS(MNVLS,MNZ,2),KVCELLS(MNVLS,MNZ,2)
      common/KEYVOLN/VOLNAME(MNVLS,MNZ),VCsurf(MNVLS,MNZ)
      common/KEYVDAT/IVTYPE(MNVLS,MNZ),VOLTemp(MNVLS,MNZ),
     &          VOLHeat(MNVLS,MNZ),IVConfl(MNVLS,MNZ),VOLHum(MNVLS,MNZ),
     &          VOLCO2(MNVLS,MNZ),VOLVel(MNVLS,MNZ),VOLDir(MNVLS,MNZ,2),
     &          VOLArea(MNVLS,MNZ),VOLPres(MNVLS,MNZ),
     &          VOLPol(MCTM,MNVLS,MNZ)
      common/param2/TITLE(MNZ),CFTRFL(MNZ),LPHI(MNZ)
      CHARACTER*72 TITLE,CFTRFL,LPHI

      COMMON/NDMAP/NOPEN(MNZ),MFNODE(MCFND,MNZ),IOPENi(MCFND,MNZ),
     &             IOPENf(MCFND,MNZ),JOPENi(MCFND,MNZ),
     &             JOPENf(MCFND,MNZ),KOPENi(MCFND,MNZ),
     &             KOPENf(MCFND,MNZ),FIXM(MCFND,MNZ),
     &             FIXT(MCFND,MNZ),FIXC(MCFND,MNZ),
     &             FIXK(MCFND,MNZ),FIXE(MCFND,MNZ),
     &             IWOPEN(MCFND,MNZ),ICFDCN(MCFND,MNZ),
     &             ICNACT(MCFND,MNZ),IVOLNOP(MCFND,MNZ)
      COMMON/CFDVIS/HAS_GEOM,ISHSB,ISHAO,IFACES,ISHBLK,ISHGEO,INITD
      logical HAS_GEOM
      character INITD*6

      character ITMSS(13)*50
      character face*6,ltmp*17,t12*12,d12*12
      character VOLNAME*12, VCsurf*12,t20*20
      integer IWO  ! for radio button
      integer NITMS,INO ! max items and current menu item
      character zngeo*3,etext*82

#ifdef OSI
      integer iside,isize,ifont
#else
      integer*8 iside,isize,ifont
#endif

C For help messages
      character helpinsub*24 ! subroutine name
      character helptopic*24 ! string (unique) for topic
      integer nbhelp     ! number of help lines found

      helpinsub='EDAIRO'  ! set for subroutine

C Adjust each cell by -1 as dfs adds fictitious cells around the domain.
      Ii=IVCELLS(IAO,ICFD,1)-1; If=IVCELLS(IAO,ICFD,2)-1
      Ji=JVCELLS(IAO,ICFD,1)-1; Jf=JVCELLS(IAO,ICFD,2)-1
      Ki=KVCELLS(IAO,ICFD,1)-1; Kf=KVCELLS(IAO,ICFD,2)-1

C Check if new opening.
      if (Ii.lt.0.or.If.lt.0.or.
     &    Ji.lt.0.or.Jf.lt.0.or.
     &    Ki.lt.0.or.Kf.lt.0) then 
        t20='Air opening'
        if(IVOLF(IAO,ICFD).ne.0)IFACE=IVOLF(IAO,ICFD)
        call PIKCELS(t20,VOLNAME(IAO,ICFD),Ii,If,Ji,Jf,Ki,Kf,
     &    IFACE,IER)
        IVOLF(IAO,ICFD)=IFACE
C        write(6,*) 'face a',IAO,ICFD,IFACE
        IVCELLS(IAO,ICFD,1)=Ii+1; IVCELLS(IAO,ICFD,2)=If+1
        JVCELLS(IAO,ICFD,1)=Ji+1; JVCELLS(IAO,ICFD,2)=Jf+1
        KVCELLS(IAO,ICFD,1)=Ki+1; KVCELLS(IAO,ICFD,2)=Kf+1
      endif

C Generate text for location index
 5    FACE='  '
      if (IVOLF(IAO,ICFD).eq.1) then
        FACE='West'
      elseif (IVOLF(IAO,ICFD).eq.2) then
        FACE='East'
      elseif (IVOLF(IAO,ICFD).eq.3) then
        FACE='South'
      elseif (IVOLF(IAO,ICFD).eq.4) then
        FACE='North'
      elseif (IVOLF(IAO,ICFD).eq.5) then
        FACE='Low'
      elseif (IVOLF(IAO,ICFD).eq.6) then
        FACE='High'
      elseif (IVOLF(IAO,ICFD).eq.7) then
        FACE='Whole'
      elseif (IVOLF(IAO,ICFD).eq.8) then
        FACE='Block'
      elseif (IVOLF(IAO,ICFD).eq.9) then
        FACE='Source'
      endif

C Set boundary type string.
      ltmp='  '
      if (IVTYPE(IAO,ICFD).eq.10) then
        ltmp='Pressure'
      elseif (IVTYPE(IAO,ICFD).eq.11) then
        ltmp='Velocity'
      elseif (IVTYPE(IAO,ICFD).eq.12) then
        ltmp='Zero gradient'
      elseif (IVTYPE(IAO,ICFD).eq.13) then
        ltmp='Mass flow'
      endif
      IF (ICFMFS(ICFD).EQ.1)THEN
        LTMP='from flow network'
      ENDIF

C Create a menu showing definitions for different conflation schemes.
      write (ITMSS(1),'(a,a)')   'a Name: ',VOLNAME(IAO,ICFD)
      write (ITMSS(2),'(a,a)')   'b Type: ',ltmp
      write (ITMSS(3),'(a)')     ' -----------------------------'
      write (ITMSS(4),'(a,a)')   '         Is  If  Js  Jf  Ks  Kf '
      write (ITMSS(5),'(a,6(i3,1x))') 'c Cells:',Ii,If,Ji,Jf,Ki,Kf
      write (ITMSS(7),'(a)')     ' -----------------------------'
      if (IVTYPE(IAO,ICFD).eq.10) then
        write(ITMSS(6),'(a)')   '  '
        write(ITMSS(8),'(a,f6.2)')'e Pressure (Pa):',VOLPres(IAO,ICFD)
      elseif (IVTYPE(IAO,ICFD).eq.11) then
        write(ITMSS(6),'(a,a)')   '  Face: ',FACE
        write(ITMSS(8),'(a,f6.2,a,2f6.2)')'e Flow rate (kg/s):',
     &        VOLVel(IAO,ICFD),', Dir (deg):',(VOLDir(IAO,ICFD,I),I=1,2)
      elseif (IVTYPE(IAO,ICFD).eq.12) then
        write(ITMSS(6),'(a,a)')   '  Face: ',FACE
        write(ITMSS(8),'(a)')   '  '
      elseif (IVTYPE(IAO,ICFD).eq.13) then
        write(ITMSS(6),'(a)')   '  '
        write(ITMSS(8),'(a,f6.2)')'e Mass (kg/s):',VOLVel(IAO,ICFD)
      endif
      write(ITMSS(9),'(a,f6.2)')'f Temperature (degC):',
     &                                             VOLTemp(IAO,ICFD)
      if (VOLArea(IAO,ICFD).gt.0.) then
        write(ITMSS(10),'(a,f6.2)')'g Real area (m^2):',
     &                                             VOLArea(IAO,ICFD)
      else
        write(ITMSS(10),'(a)')'g Real area (m^2): cell area'
      endif
      ITMSS(11)=' -----------------------------'
      ITMSS(12)='? Help'
      ITMSS(13)='- exit'
      NITMS=13

      IF(MMOD.EQ.8)THEN
C Find index in old commons of opening specified by index IAO in new commons
        do 638 I=1,NOPEN(ICFD)
          if (IAO.eq.IVOLNOP(I,ICFD)) then
            ISHAO_TMP=ISHAO
            ISHGEO_TMP=ISHGEO
            ISHAO=I
            ISHGEO=0
            CALL redraw(IER)
            ISHAO=ISHAO_TMP
            ISHGEO=ISHGEO_TMP
          endif
638     continue
      ENDIF

C Help text for this menu.
      helptopic='cfd_air_boundary'
      call gethelptext(helpinsub,helptopic,nbhelp)

C Display menu.
      INO=-1
      CALL EMENU('Air flow opening edit',ITMSS,NITMS,INO)

      if (INO.eq.1) then
        t12=VOLNAME(IAO,ICFD)
        d12='  '
        call EASKS(t12,'Opening name? ','  ',12,d12,'open name',
     &    IER,nbhelp)
        if (IER.eq.0) then
          call st2name(t12,d12)
          VOLNAME(IAO,ICFD)=d12
        endif
      elseif (INO.eq.2) then

C Choose an opening type.
        IWO=1
        CALL EASKATOG('Opening types: (see help)',' ','pressure',
     &    'velocity','zero gradient','from flow network',
     &    'cancel','  ','  ',IWO,nbhelp)
        IF(IWO.EQ.4)THEN
          ICFMFS(ICFD)=1
          CALL MFCONF(ICFD)
        ENDIF

C Set BC type.
        if (IWO.ge.1.AND.IWO.le.3) THEN
          IVTYPE(IAO,ICFD)=IWO+9
          ICFMFS(ICFD)=0
        ENDIF

C Reset location to whole of volume for non face choices.
        if (IWO.eq.1.OR.IWO.eq.3) then
          IVOLF(IAO,ICFD)=7
        else
          t20='Air opening'
          call PIKFACE(t20,VOLNAME(IAO,ICFD),Ii,If,Ji,Jf,Ki,Kf,
     &      IVOLF(IAO,ICFD),IER)
        endif
      elseif (INO.eq.5) then
        t20='Air opening'
        if(IVOLF(IAO,ICFD).ne.0)IFACE=IVOLF(IAO,ICFD)
        call PIKCELS(t20,VOLNAME(IAO,ICFD),Ii,If,Ji,Jf,Ki,Kf,
     &    IFACE,IER)
        IVOLF(IAO,ICFD)=IFACE
C        write(6,*) 'face b',IAO,ICFD,IFACE
        IVCELLS(IAO,ICFD,1)=Ii+1; IVCELLS(IAO,ICFD,2)=If+1
        JVCELLS(IAO,ICFD,1)=Ji+1; JVCELLS(IAO,ICFD,2)=Jf+1
        KVCELLS(IAO,ICFD,1)=Ki+1; KVCELLS(IAO,ICFD,2)=Kf+1
      elseif (INO.eq.8) then
        if (IVTYPE(IAO,ICFD).eq.10) then
          VAL=VOLPres(IAO,ICFD)
          call EASKR(VAL,' ',' Pressure (Pa) of cells?',
     &      0.,'W',10000.,'W',50.0,'pressure',IER,nbhelp)
          if (IER.eq.0) VOLPres(IAO,ICFD)=VAL
        elseif (IVTYPE(IAO,ICFD).eq.11) then
          VAL=VOLVel(IAO,ICFD)
          call EASKR(VAL,' ',
     &      'Mass flow rate (kg/s) at chosen face of cells?',
     &       -50.,'W',50.,'W',0.5,'flow rate',IER,nbhelp)
          if (IER.eq.0) VOLVel(IAO,ICFD)=VAL
        elseif (IVTYPE(IAO,ICFD).eq.13) then
          VAL=VOLVel(IAO,ICFD)
          call EASKR(VAL,' ','Mass flow rate (kg/s) entering cells?',
     &      -10.,'W',10.,'W',0.02,'mass flow',IER,nbhelp)
          if (IER.eq.0) VOLVel(IAO,ICFD)=VAL
        endif
      elseif (INO.eq.9) then
        VAL=VOLTemp(IAO,ICFD)
        call EASKR(VAL,' ','Temperature (degC) of air entering domian?',
     &     -50.,'W',50.,'W',20.,'temperature',IER,nbhelp)
        if (IER.eq.0) VOLTemp(IAO,ICFD)=VAL
      elseif (INO.eq.9) then
        VAL=VOLArea(IAO,ICFD)
        call EASKR(VAL,' ','Real area (m^2) of opening?',
     &    -0.1,'W',50.,'W',0.,'area',IER,nbhelp)
        if (IER.eq.0) VOLArea(IAO,ICFD)=VAL
      elseif (INO.eq.(NITMS-1)) then

C Help.
        helptopic='cfd_air_boundary'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL PHELPD('boundary conditions',nbhelp,'-',0,0,IER)
      elseif (INO.eq.NITMS) then

C Return.
        return
      endif

      goto 5
      end

C ************************** MFCONF **************************
C MFCONF defines conflation of CFD and mfs 

      SUBROUTINE MFCONF(ICFD)
#include "cfd.h"
#include "net_flow.h"

      common/KEYVOLS/NVOL(MNZ),IVOLF(MNVLS,MNZ),IVCELLS(MNVLS,MNZ,2),
     &               JVCELLS(MNVLS,MNZ,2),KVCELLS(MNVLS,MNZ,2)
      common/KEYVOLN/VOLNAME(MNVLS,MNZ),VCsurf(MNVLS,MNZ)
      common/KEYVDAT/IVTYPE(MNVLS,MNZ),VOLTemp(MNVLS,MNZ),
     &          VOLHeat(MNVLS,MNZ),IVConfl(MNVLS,MNZ),VOLHum(MNVLS,MNZ),
     &          VOLCO2(MNVLS,MNZ),VOLVel(MNVLS,MNZ),VOLDir(MNVLS,MNZ,2),
     &          VOLArea(MNVLS,MNZ),VOLPres(MNVLS,MNZ),
     &          VOLPol(MCTM,MNVLS,MNZ)
      COMMON/CFDMFS1/ICFDNOD(MNZ),ICFDCNN(MNVLS,MNZ)
      DIMENSION ICPK(MCNN),INPK(MNOD)

      character VOLNAME*12, VCsurf*12
      CHARACTER PROMPT1*72, PROMPT2*72

C For help messages
      character helpinsub*24 ! subroutine name
      character helptopic*24 ! string (unique) for topic
      integer nbhelp     ! number of help lines found

      helpinsub='MFCONF'  ! set for subroutine

C Note: use verions of askmfnod and askmfcon defined in edcfd.F
C rather than the ones defined in esrures/mfoutp.F

C Select mfs node representative of CFD domain
      helptopic='cfd_mass_flow_node'
      call gethelptext(helpinsub,helptopic,nbhelp)
      INOPT=0
      NNPK=1
      PROMPT1='Which node is represetative of CFD domain'
      PROMPT2='  '
      call ASKMFNOD(INOPT,NNPK,INPK,PROMPT1,PROMPT2,nbhelp)
      if (NNPK.gt.0.and.INPK(1).gt.0) then
        ICFDNOD(ICFD)=INPK(1)
      endif

C Select mfs connections representative of all CFD openings
      INODE=INPK(1)
      DO 100 I=1,NVOL(ICFD)
        IF(IVTYPE(I,ICFD).GE.10.AND.IVTYPE(I,ICFD).LE.13)THEN
          WRITE(PROMPT1,'(2A)')
     &    'Which connection is represented by opening ',VOLNAME(I,ICFD)
          CALL ASKMFCON(INODE,NNPK,ICPK,PROMPT1,PROMPT2)
          if (NNPK.gt.0.and.ICPK(1).gt.0) then
            ICFDCNN(I,ICFD)=ICPK(1)
          endif
        ENDIF
 100  CONTINUE

      RETURN
      END


C ************************** EDSLDB **************************
C EDSLDB controls the editing of CFD solid boundary variables.
C IER=0 indicates no error.

      SUBROUTINE EDSLDB(IZONE,ISB,IER)
#include "building.h"
C #include "model.h"
#include "geometry.h"
#include "cfd.h"
C #include "prj3dv.h"

      COMMON/OUTIN/IUOUT,IUIN
      COMMON/SPAD/MMOD,LIMIT,LIMTTY
      COMMON/ICFNOD/ICFD,ICP
      COMMON/GEOM/XP(ntcelx),YP(ntcely),ZP(ntcelz),
     1            DXEP(ntcelx),DXPW(ntcelx),DYNP(ntcely),DYPS(ntcely),
     2            DZHP(ntcelz),DZPL(ntcelz),
     3            SEW(ntcelx),SNS(ntcely),SHL(ntcelz),
     4            XU(ntcelx),YV(ntcely),ZW(ntcelz)

      common/METHDS/ITURB(MNZ),IBUOY(MNZ)
      common/KEYVOLS/NVOL(MNZ),IVOLF(MNVLS,MNZ),IVCELLS(MNVLS,MNZ,2),
     &               JVCELLS(MNVLS,MNZ,2),KVCELLS(MNVLS,MNZ,2)
      common/KEYVOLN/VOLNAME(MNVLS,MNZ),VCsurf(MNVLS,MNZ)
      common/KEYVDAT/IVTYPE(MNVLS,MNZ),VOLTemp(MNVLS,MNZ),
     &          VOLHeat(MNVLS,MNZ),IVConfl(MNVLS,MNZ),VOLHum(MNVLS,MNZ),
     &          VOLCO2(MNVLS,MNZ),VOLVel(MNVLS,MNZ),VOLDir(MNVLS,MNZ,2),
     &          VOLArea(MNVLS,MNZ),VOLPres(MNVLS,MNZ),
     &          VOLPol(MCTM,MNVLS,MNZ)
      COMMON/CFDVIS/HAS_GEOM,ISHSB,ISHAO,IFACES,ISHBLK,ISHGEO,INITD
      logical HAS_GEOM
      character INITD*6
      common/param2/TITLE(MNZ),CFTRFL(MNZ),LPHI(MNZ)

      COMMON/Sbdary/NSB(MNZ),ISBi(MNSBZ,MNZ),ISBf(MNSBZ,MNZ),
     &              JSBi(MNSBZ,MNZ),JSBf(MNSBZ,MNZ),
     &              KSBi(MNSBZ,MNZ),KSBf(MNSBZ,MNZ),
     &              ISUFLC(MNSBZ,MNZ),IWSB(MNSBZ,MNZ),SSB(MNSBZ,MNZ),
     &              SSBHC(MNSBZ,MNZ),IVOLNSB(MNSBZ,MNZ),
     &              ITCtype(MNSBZ,MNZ),icTREF(MNSBZ,MNZ)

      character ITMSS(13)*40
      character face*6,ltmp*16,t12*12,d12*12,outs*124
      character VOLNAME*12, VCsurf*12
      CHARACTER t20*20,ETEXT*82,zngeo*3
      character TITLE*72,CFTRFL*72,LPHI*72
      integer IWO  ! for radio button
      integer NITMS,INO ! max items and current menu item

C For help messages
      character helpinsub*24 ! subroutine name
      character helptopic*24 ! string (unique) for topic
      integer nbhelp     ! number of help lines found

#ifdef OSI
      integer iside,isize,ifont
#else
      integer*8 iside,isize,ifont
#endif

      helpinsub='EDSLDB'  ! set for subroutine

C Adjust each cell by -1 as dfs adds fictitious cells around the domain.
      Ii=IVCELLS(ISB,ICFD,1)-1; If=IVCELLS(ISB,ICFD,2)-1
      Ji=JVCELLS(ISB,ICFD,1)-1; Jf=JVCELLS(ISB,ICFD,2)-1
      Ki=KVCELLS(ISB,ICFD,1)-1; Kf=KVCELLS(ISB,ICFD,2)-1
      itrc=0   ! silent feedback initially
      XZ=0.; YZ=0.; ZZ=0.

C Check if new solid boundary.
      if (Ii.lt.0.or.If.lt.0.or.
     &    Ji.lt.0.or.Jf.lt.0.or.
     &    Ki.lt.0.or.Kf.lt.0) then 
        t20='Solid boundary'
        if(IVOLF(ISB,ICFD).ne.0)IFACE=IVOLF(ISB,ICFD)
        call PIKCELS(t20,VOLNAME(ISB,ICFD),Ii,If,Ji,Jf,Ki,Kf,
     &    IFACE,IER)
        IVOLF(ISB,ICFD)=IFACE
C        write(6,*) 'face c',ISB,ICFD,IFACE
        IVCELLS(ISB,ICFD,1)=Ii+1; IVCELLS(ISB,ICFD,2)=If+1
        JVCELLS(ISB,ICFD,1)=Ji+1; JVCELLS(ISB,ICFD,2)=Jf+1
        KVCELLS(ISB,ICFD,1)=Ki+1; KVCELLS(ISB,ICFD,2)=Kf+1
      endif

C Generate text for location index
 5    FACE='  '
      if (IVOLF(ISB,ICFD).eq.1) then
        FACE='West'
      elseif (IVOLF(ISB,ICFD).eq.2) then
        FACE='East'
      elseif (IVOLF(ISB,ICFD).eq.3) then
        FACE='South'
      elseif (IVOLF(ISB,ICFD).eq.4) then
        FACE='North'
      elseif (IVOLF(ISB,ICFD).eq.5) then
        FACE='Low'
      elseif (IVOLF(ISB,ICFD).eq.6) then
        FACE='High'
      elseif (IVOLF(ISB,ICFD).eq.7) then
        FACE='Whole'
      elseif (IVOLF(ISB,ICFD).eq.8) then
        FACE='Block'
      elseif (IVOLF(ISB,ICFD).eq.9) then
        FACE='Source'
      endif

C Set boundary type string.
      ltmp='  '
      if (IVTYPE(ISB,ICFD).eq.1) then
        ltmp='Temperature'
      elseif (IVTYPE(ISB,ICFD).eq.2) then
        ltmp='Heat flux'
      elseif (IVTYPE(ISB,ICFD).eq.3) then
        ltmp='Symmetrical'
      elseif (IVTYPE(ISB,ICFD).eq.4) then
        ltmp='Conflated  '
      elseif (IVTYPE(ISB,ICFD).eq.5) then
        ltmp='Conflated  '
      elseif (IVTYPE(ISB,ICFD).eq.6) then
        ltmp='Conflated  '
      endif

C Create a menu showing definitions for different conflation schemes.
      write (ITMSS(1),'(a,a)')   'a Name: ',VOLNAME(ISB,ICFD)
      write (ITMSS(2),'(a,a)')   'b Type: ',ltmp
      write (ITMSS(3),'(a)')     ' -----------------------------'
      write (ITMSS(4),'(a,a)')   '         Is  If  Js  Jf  Ks  Kf '
      write (ITMSS(5),'(a,6(i3,1x))') 'c Cells:',Ii,If,Ji,Jf,Ki,Kf

C If not silent echo some XYZ coordinates in space. Need to offset
C indicies to get correct position.
C      write(6,*) 'XU Ii Ii+1 If If+1 If+2 If+3 ',Ii,If,
C     &  XU(Ii),XU(Ii+1),XU(If),XU(If+1),XU(If+2),XU(If+3)
C      write(6,*) 'XP Ii Ii+1 If If+1 If+2 ',
C     &  XP(Ii),XP(Ii+1),XP(If),XP(If+1),XP(If+2)
C      write(6,*) 'DXEP Ii Ii+1 If If+1 If+2 If+3',Ii,If,
C     &  DXEP(Ii),DXEP(Ii+1),DXEP(If),DXEP(If+1),DXEP(If+2)
      call edisp(iuout,
     &  'Cell X start & end     Y start & end     Z start & end')
      write(outs,'(6f10.4)') XU(Ii+1),XU(If+2),YV(Ji+1),YV(Jf+2),
     &  ZW(Ki+1),ZW(Kf+2)
      call edisp(iuout,outs)
      write (ITMSS(6),'(a,a)')   '  Face: ',FACE
      write (ITMSS(7),'(a)')     ' -----------------------------'
      if (IVTYPE(ISB,ICFD).eq.1) then
        write(ITMSS(8),'(a,f6.2)')'e Temperature (degC): ',
     &    VOLTemp(ISB,ICFD)
        write(ITMSS(9),'(a)')     '  '
      elseif (IVTYPE(ISB,ICFD).eq.2) then
        write(ITMSS(8),'(a,f7.0)')'e Heat flux (W) ',VOLHeat(ISB,ICFD)
        write(ITMSS(9),'(a)')     '  '
      elseif (IVTYPE(ISB,ICFD).eq.3) then
        write(ITMSS(8),'(a)')   '  '
        write(ITMSS(9),'(a)')   '  '
      elseif (IVTYPE(ISB,ICFD).eq.4) then
        write(ITMSS(8),'(a,a)')'e Conflated to surface: ',
     &                                              VCsurf(ISB,ICFD)
        ivcfl=IVConfl(ISB,ICFD)
        if(ivcfl.eq.0)then
          write(ITMSS(9),'(a)')'f Handshaking: not yet defined'
        elseif(ivcfl.eq.1)then
          write(ITMSS(9),'(a)')'f Handshaking: One-way log-law CFD'
        elseif(ivcfl.eq.3)then
          write(ITMSS(9),'(a)')'f Handshaking: One-way Yuan    CFD'
        elseif(ivcfl.eq.4)then
          write(ITMSS(9),'(a)')'f Handshaking: One-way Yuan   BSim'
        elseif(ivcfl.eq.5)then
          write(ITMSS(9),'(a)')'f Handshaking: One-way Yuan CFD avg'
        elseif(ivcfl.eq.6)then
          write(ITMSS(9),'(a)')'f Handshaking: One-way log-law BSim'
        elseif(ivcfl.eq.7)then
          write(ITMSS(9),'(a)')'f Handshaking: One-way log-law CFD avg'
        elseif(ivcfl.eq.8)then
          write(ITMSS(9),'(a)')'f Handshaking: One-way log-law CFD loc'
        elseif(ivcfl.eq.9)then
          write(ITMSS(9),'(a)')'f Handshaking: Two-way log-law CFD'
        elseif(ivcfl.eq.11)then
          write(ITMSS(9),'(a)')'f Handshaking: Two-way Yuan CFD'
        elseif(ivcfl.eq.12)then
          write(ITMSS(9),'(a)')'f Handshaking: Two-way Yuan CFD loc'
        elseif(ivcfl.eq.13)then
          write(ITMSS(9),'(a)')'f Handshaking: Two-way log-law CFD'
        elseif(ivcfl.eq.14)then
          write(ITMSS(9),'(a)')'f Handshaking: Two-way log-law local'
        else
          write(ITMSS(9),'(a,i4)')'f Handshaking: UNKNOWN type',ivcfl
        endif
C        write(ITMSS(9),'(a,i4)')'f Handshaking: ',IVConfl(ISB,ICFD)
      endif
      ITMSS(10) =' -----------------------------'
      if(itrc.eq.0)then
        ITMSS(11) ='> feedback silent             '
      elseif(itrc.eq.1)then
        ITMSS(11) ='> feedback brief              '
      elseif(itrc.eq.2)then
        ITMSS(11) ='> feedback verbose            '
      endif
      
      ITMSS(12)='? Help'
      ITMSS(13)='- exit'
      NITMS=13

      IF(MMOD.EQ.8)THEN
C Find index in old commons of boundary specified by index ISB in new commons
        do 195 I=1,NSB(ICFD)
          if (ISB.eq.IVOLNSB(I,ICFD)) then
            ISHSB_TMP=ISHSB
            ISHGEO_TMP=ISHGEO
            ISHSB=I
            ISHGEO=0
            CALL redraw(IER)
            ISHSB=ISHSB_TMP
            ISHGEO=ISHGEO_TMP
          endif
195     continue
      ENDIF

C Help text for this menu.
      helptopic='cfd_solid_boundary_edit'
      call gethelptext(helpinsub,helptopic,nbhelp)

C Display menu.
      INO=-1
      CALL EMENU('Solid boundary edit',ITMSS,NITMS,INO)
      if (INO.eq.1) then

C Get boundary name.
        t12=VOLNAME(ISB,ICFD)
        d12='  '
        call EASKS(t12,'Boundary name? ','  ',12,d12,'bnd name',
     &    IER,nbhelp)
        if (IER.eq.0) then
          call st2name(t12,d12)
          VOLNAME(ISB,ICFD)=d12
        endif
      elseif (INO.eq.2) then

C Choose an opening type.
        IWO=1
        CALL EASKATOG('Boundary types:',' ','temperature','heat flux',
     &    'symmetrical','building surface','cancel','  ','  ',
     &    IWO,nbhelp)

C Try to remember non conflated option if building surface chosen.
        if (IWO.eq.4) then
          if (ITURB(ICFD).eq.0) then
            call usrmsg('Building surface boundary conditions are only',
     &                'applicable in turbulent flows. ','W')
            call usrmsg('Please select a method for modelling  ',
     &         'turbulence via the solution variables menu pick. ','W')
          else
            IVTYPE(ISB,ICFD)=IVTYPE(ISB,ICFD)+3
          endif
        else
          IVTYPE(ISB,ICFD)=IWO
        endif
      elseif (INO.eq.5) then
        t20='Solid boundary'
        if(IVOLF(ISB,ICFD).ne.0)IFACE=IVOLF(ISB,ICFD)
        call PIKCELS(t20,VOLNAME(ISB,ICFD),Ii,If,Ji,Jf,Ki,Kf,
     &    IFACE,IER)
        IVOLF(ISB,ICFD)=IFACE
C        write(6,*) 'face d',ISB,ICFD,IFACE
        IVCELLS(ISB,ICFD,1)=Ii+1; IVCELLS(ISB,ICFD,2)=If+1
        JVCELLS(ISB,ICFD,1)=Ji+1; JVCELLS(ISB,ICFD,2)=Jf+1
        KVCELLS(ISB,ICFD,1)=Ki+1; KVCELLS(ISB,ICFD,2)=Kf+1
      elseif (INO.eq.8) then
        if (IVTYPE(ISB,ICFD).eq.1) then
          VAL=VOLTemp(ISB,ICFD)
          call EASKR(VAL,'temperature (degC) at chosen face of cells?',
     &      ' ',-100.,'W',100.,'W',10.,'temperature',IER,nbhelp)
          if (IER.eq.0) VOLTemp(ISB,ICFD)=VAL
        elseif (IVTYPE(ISB,ICFD).eq.2) then
          VAL=VOLHeat(ISB,ICFD)
          call EASKR(VAL,'Heat flux (W) at chosen face of cells?',
     &      ' ',-1000.,'W',1000.,'W',10.,'heat flux',IER,nbhelp)
          if (IER.eq.0) VOLHeat(ISB,ICFD)=VAL
        elseif (IVTYPE(ISB,ICFD).eq.4) then
          CALL EASKSUR(IZONE,IS,'-','Select a surface.',' ',IER)
          VCsurf(ISB,ICFD)=SNAME(IZONE,IS)
        endif
      elseif (INO.eq.9) then
        if (IVTYPE(ISB,ICFD).eq.4) then
          call HANDSHK(IHST,IER)
          if (IER.eq.0) IVConfl(ISB,ICFD)=IHST
        endif
      elseif (INO.eq.(NITMS-1)) then

C Help.
        helptopic='cfd_solid_boundary_edit'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL PHELPD('boundary conditions',nbhelp,'-',0,0,IER)
      elseif (INO.eq.NITMS) then

C Return.
        return
      endif

      goto 5
      end

C ************************** EDSRC **************************
C EDSRC controls the editing of CFD sources - humidity, CO2, etc.
C IER=0 indicates no error.

      SUBROUTINE EDSRC(ISB,IER)
#include "cfd.h"
#include "net_flow.h"
#include "net_flow_data.h"

      COMMON/SPAD/MMOD,LIMIT,LIMTTY
      COMMON/ICFNOD/ICFD,ICP
      common/KEYVOLS/NVOL(MNZ),IVOLF(MNVLS,MNZ),IVCELLS(MNVLS,MNZ,2),
     &               JVCELLS(MNVLS,MNZ,2),KVCELLS(MNVLS,MNZ,2)
      common/KEYVOLN/VOLNAME(MNVLS,MNZ),VCsurf(MNVLS,MNZ)
      common/KEYVDAT/IVTYPE(MNVLS,MNZ),VOLTemp(MNVLS,MNZ),
     &          VOLHeat(MNVLS,MNZ),IVConfl(MNVLS,MNZ),VOLHum(MNVLS,MNZ),
     &          VOLCO2(MNVLS,MNZ),VOLVel(MNVLS,MNZ),VOLDir(MNVLS,MNZ,2),
     &          VOLArea(MNVLS,MNZ),VOLPres(MNVLS,MNZ),
     &          VOLPol(MCTM,MNVLS,MNZ)
      common/KEYCASGN/IDcasgn(MNVLS,MNZ),Fcasgn(MNVLS,MNZ)
      common/EQTION3/CALLMA(MNZ),CALPOL(MCTM,MNZ),POLNAM(MCTM,MNZ),
     &               NCTM(MNZ),JHUMINDX(MNZ),URFC(MCTM)
      COMMON/SCHMTT/SCHMT(MCTM),GFM(MCTM),VCRIT(MCTM),TBOIL(MCTM),
     &       TCRIT(MCTM),ISCHMT(MCTM)
      COMMON/CTDFAF/ICTDFAF,SRCE(MNVLS,MCTM,MNZ),ICCSRC
      COMMON/CONTM0/NCONTM,NOCNTM,CONTMNAM(MCONTM)
      COMMON/CONTM5/SPMSUP(MSPMNO,MCSD),SSLINK2(MSPMNO,MNOD),
     &SPMTYP(MSPMNO),SSNAME(MSPMNO),NSPMNO,SSLINK1(MSPMNO,MCONTM)
      common/cfdconf/ICFBLD(MNZ),ICFMFS(MNZ)
      COMMON/CFDMFS1/ICFDNOD(MNZ),ICFDCNN(MNVLS,MNZ)
      COMMON/CFDVIS/HAS_GEOM,ISHSB,ISHAO,IFACES,ISHBLK,ISHGEO,INITD
      logical HAS_GEOM
      character INITD*6
      common/param2/TITLE(MNZ),CFTRFL(MNZ),LPHI(MNZ)
      common/pophelp/h(60)   ! needed for dynamic popup
      character H*72

      character ITMSS(14+MCTM)*60,ITUM(MCTM+1)*50,WORD*50
      character ltmp*26,t12*12,d12*12,t20*20,FACE*7,HOLD*32,WORD1*12
      character*12 VOLNAME,VCsurf,SRCNAM,POLNAM,CONTMNAM,SSNAME,SRCE
      CHARACTER CNAM*12,CNAME*3,NNAME*3,ITMSS1(14+MCTM)*60
      CHARACTER*60 ITMSS2(14+MCTM)*60,ITMSS3(14+MCTM)*60
      character outs*124
      LOGICAL CALPOL,CALLMA,OK,DOK
      INTEGER SPMTYP,SSLINK1,SSLINK2
      integer NITMS,INO,NITMSS1,IIO,NITMSS2,IXO,NITUM,INUO,NITMSS3 ! max items and current menu item
      character TITLE*72,CFTRFL*72,LPHI*72
      CHARACTER ETEXT*82,zngeo*3

#ifdef OSI
      integer iside,isize,ifont
#else
      integer*8 iside,isize,ifont
#endif

C For help messages
      character helpinsub*24 ! subroutine name
      character helptopic*24 ! string (unique) for topic
      integer nbhelp     ! number of help lines found

      helpinsub='EDSRC'  ! set for subroutine

C Initialise SRCE 
      IF(ICTDFAF.NE.1)THEN
        DO 420 ICTM=1,NCTM(ICFD)
          DO 520 II=1,NVOL(ICFD)
            SRCE(II,ICTM,ICFD)=' '
 520      CONTINUE
 420    CONTINUE
      ENDIF

C Adjust each cell by -1 as dfs adds fictitious cells around the domain.
      Ii=IVCELLS(ISB,ICFD,1)-1; If=IVCELLS(ISB,ICFD,2)-1
      Ji=JVCELLS(ISB,ICFD,1)-1; Jf=JVCELLS(ISB,ICFD,2)-1
      Ki=KVCELLS(ISB,ICFD,1)-1; Kf=KVCELLS(ISB,ICFD,2)-1

C Check if new source.
      if (Ii.lt.0.or.If.lt.0.or.
     &    Ji.lt.0.or.Jf.lt.0.or.
     &    Ki.lt.0.or.Kf.lt.0) then 
        t20='Source definition'
        if(IVOLF(ISB,ICFD).ne.0)IFACE=IVOLF(ISB,ICFD)
        call PIKCELS(t20,VOLNAME(ISB,ICFD),Ii,If,Ji,Jf,Ki,Kf,
     &    IFACE,IER)
        IVOLF(ISB,ICFD)=IFACE
C        write(6,*) 'face e',ISB,ICFD,IFACE
        IVCELLS(ISB,ICFD,1)=Ii+1; IVCELLS(ISB,ICFD,2)=If+1
        JVCELLS(ISB,ICFD,1)=Ji+1; JVCELLS(ISB,ICFD,2)=Jf+1
        KVCELLS(ISB,ICFD,1)=Ki+1; KVCELLS(ISB,ICFD,2)=Kf+1
      endif

C Generate text for location index
 5    FACE='  '
      if (IVOLF(ISB,ICFD).eq.1) then
        FACE='West   '
      elseif (IVOLF(ISB,ICFD).eq.2) then
        FACE='East   '
      elseif (IVOLF(ISB,ICFD).eq.3) then
        FACE='South  '
      elseif (IVOLF(ISB,ICFD).eq.4) then
        FACE='North  '
      elseif (IVOLF(ISB,ICFD).eq.5) then
        FACE='Low    '
      elseif (IVOLF(ISB,ICFD).eq.6) then
        FACE='High   '
      elseif (IVOLF(ISB,ICFD).eq.7) then
        FACE='Whole  '
      elseif (IVOLF(ISB,ICFD).eq.8) then
        FACE='Block  '
      elseif (IVOLF(ISB,ICFD).eq.9) then
        FACE='Source '
      endif

C Set boundary type string.
      ltmp=FACE
 144  if (IVTYPE(ISB,ICFD).eq.20)then
        ltmp='Static source'
        IF(ICTDFAF.EQ.1)LTMP='From ctm network'
      endif

C Create a menu showing definitions for different sources.
      write (ITMSS(1),'(a,a)')   'a Name: ',VOLNAME(ISB,ICFD)
      write (ITMSS(2),'(a,a)')   'b Type:',ltmp
      write (ITMSS(3),'(a)')     ' -----------------------------'
      write (ITMSS(4),'(a,a)')   '         Is  If  Js  Jf  Ks  Kf '
      write (ITMSS(5),'(a,6(i3,1x))') 'c Cells:',Ii,If,Ji,Jf,Ki,Kf
      write (ITMSS(6),'(a)')     ' -----------------------------'
      write(ITMSS(7),'(a,f8.0)')'e Heat flux(W)',VOLHeat(ISB,ICFD)
      write(ITMSS(8),'(a,i3)') 'f Cas gain index  ',IDCasgn(ISB,ICFD)
      write(ITMSS(9),'(a,f4.2)')'g Cas gain fraction ',FCasgn(ISB,ICFD)
      write (ITMSS(10),'(a)')    '  Contaminant | source strength '
      WRITE (ITMSS(11),'(A)')    '     name     | (kg/s) or name  '
      IF(ICTDFAF.EQ.1)THEN
        DO 141 ICTM=1,NCTM(ICFD)
          CALPOL(ICTM,ICFD)=.TRUE.
          WORD1=SRCE(ISB,ICTM,ICFD)
          IF(WORD1(:1).EQ.' ')WORD1='none'
          IF(JHUMINDX(ICFD).NE.ICTM)THEN
            WRITE(ITMSS(11+ICTM),'(A,2X,A,1X,A)')CHAR(103+ICTM),
     &      POLNAM(ICTM,ICFD),WORD1
          ELSE
            WRITE(ITMSS(11+ICTM),'(A,2X,A,1X,2A)')CHAR(103+ICTM),
     &      POLNAM(ICTM,ICFD),WORD1,' (humidity) '
          ENDIF
 141    CONTINUE   
      ELSE
        DO 131 ICTM=1,NCTM(ICFD)
          CALPOL(ICTM,ICFD)=.TRUE.
          IF(JHUMINDX(ICFD).NE.ICTM)THEN
            WRITE(ITMSS(11+ICTM),'(A,1X,A,1X,F12.10)')CHAR(103+ICTM),
     &      POLNAM(ICTM,ICFD),VOLPOL(ICTM,ISB,ICFD)
          ELSE
            WRITE(ITMSS(11+ICTM),'(A,1X,A,1X,F12.10,A)')CHAR(103+ICTM),
     &      POLNAM(ICTM,ICFD),VOLPOL(ICTM,ISB,ICFD),' (humidity) '
          ENDIF
 131    CONTINUE   
      ENDIF
      ITMSS(12+NCTM(ICFD)) =' -----------------------------'
      ITMSS(13+NCTM(ICFD)) = '+ Add/Delete Contaminant'
      ITMSS(14+NCTM(ICFD))='? Help'
      ITMSS(15+NCTM(ICFD))='- exit'
      NITMS=15+NCTM(ICFD)

C Help text for this menu.
  7   helptopic='cfd_source_boundary_edit'
      call gethelptext(helpinsub,helptopic,nbhelp)

C If in graphic model show the current grid.
      IF(MMOD.EQ.8)THEN

C Recalculate gridding in case this has changed.
        CALL INICNT
        CALL GRID(ier)
        call NEW2OLD

        call redraw(ier)
      ENDIF

C Display menu.
      INO=-1
      IF(VOLNAME(ISB,ICFD)(1:4).EQ.'SRC_')THEN
        CALL PHELPD('Warning',nbhelp,'-',0,0,IER)
      ENDIF
      CALL EMENU('Sources edit',ITMSS,NITMS,INO)
      write(6,*)INO
      if (INO.lt.1) then
        goto 7
      elseif (ITMSS(INO)(1:1).eq.' ') then
        goto 7
      elseif (INO.eq.1) then

C Get boundary name.
        t12=VOLNAME(ISB,ICFD)
        d12='  '
        call EASKS(t12,'Source name? ','  ',12,d12,'src name',
     &    IER,nbhelp)
        if (IER.eq.0) then
          call st2name(t12,d12)
          VOLNAME(ISB,ICFD)=d12
        endif

C Set coupling with contaminants network on/off
      elseif (INO.eq.2) then
        CALL EASKABC('How to specify coundary condition',' ','static',
     &    'from contaminants network','cancel',IHOWL,nbhelp)
        IF(IHOWL.EQ.1)THEN
          ICTDFAF=0
        ELSEIF(IHOWL.EQ.2)THEN

C Contaminants currently only work for CFD domain 1.
C If the user tries to input contaminants in any other domain,
C inform them and do not allow it.
          IF(ICFD.NE.1)GOTO 756

C Take contaminant definitions from contaminant flow network
          DOK=.TRUE.
          CALL ASKOK('Link CFD with contaminant network (see help) ?',
     &      ' ',OK,dok,nbhelp)
          IF(.NOT.OK)GOTO 144
          IF(ICFMFS(ICFD).NE.1)THEN
            CALL PHELPD(' ',nbhelp,'-',0,0,IER)
          ENDIF
          ICTDFAF=1
          NCTM(ICFD)=NCONTM
          IF(ICCSRC.NE.1)THEN
            ICCSRC=1

C Add one more source type boundary condition for each opening type
C boundary conditions and define the source strength to be from 
            NOBC=NVOL(ICFD)
            DO 3056 IV=1,NOBC
              IF(IVTYPE(IV,ICFD).GE.10.AND.IVTYPE(IV,ICFD).LE.13)THEN
                NVOL(ICFD)=NVOL(ICFD)+1
                NUMVOL=NVOL(ICFD)
                WRITE(VOLNAME(NUMVOL,ICFD),'(2A)')'SRC_',
     &          VOLNAME(IV,ICFD)(1:8)
                IVTYPE(NUMVOL,ICFD)=20
                IVOLF(NUMVOL,ICFD)=9
                DO 3057 L=1,2
                  IVCELLS(NUMVOL,ICFD,L)=IVCELLS(IV,ICFD,L)
                  JVCELLS(NUMVOL,ICFD,L)=JVCELLS(IV,ICFD,L)
                  KVCELLS(NUMVOL,ICFD,L)=KVCELLS(IV,ICFD,L)
 3057           CONTINUE
                VOLHEAT(NUMVOL,ICFD)=0.0
                IDCASGN(NUMVOL,ICFD)=0
                FCASGN(NUMVOL,ICFD)=0.0
                DO 3058 ICTM=1,NCTM(ICFD)
                  VOLPOL(ICTM,NUMVOL,ICFD)=0.0
                  CALL INTSTR(ICTM,CNAM,ILEN,IER)
                  IF(ILEN.EQ.1)THEN
                    WRITE(CNAME,'(2A)')'00',CNAM(:1)
                  ELSEIF(ILEN.EQ.2)THEN                
                    WRITE(CNAME,'(2A)')'0',CNAM(:2)
                  ELSEIF(ILEN.EQ.3)THEN
                    WRITE(CNAME,'(A)')CNAM(:3)
                  ENDIF
                  ICNN=ICFDCNN(IV,ICFD)
                  IF(NODPS(ICNN).EQ.ICFDNOD(ICFD))THEN
                    IEXTNOD=NODNE(ICNN)
                  ELSEIF(NODNE(ICNN).EQ.ICFDNOD(ICFD))THEN
                    IEXTNOD=NODPS(ICNN)
                  ELSE
                    GOTO 3059
                  ENDIF
                  CALL INTSTR(IEXTNOD,CNAM,ILEN,IER)
                  IF(ILEN.EQ.1)THEN
                    WRITE(NNAME,'(2A)')'00',CNAM(:1)
                  ELSEIF(ILEN.EQ.2)THEN                
                    WRITE(NNAME,'(2A)')'0',CNAM(:2)
                  ELSEIF(ILEN.EQ.3)THEN
                    WRITE(NNAME,'(A)')CNAM(:3)
                  ENDIF
                  WRITE(SRCE(NUMVOL,ICTM,ICFD),'(3A)')'SRC',CNAME,NNAME
 3059             CONTINUE
 3058           CONTINUE
              ENDIF
 3056       CONTINUE
          ENDIF


C Ask which contaminant (or none) is moisture 
          ITMSS1(1) = ' name:  '
          ITMSS1(2) ='         '
          DO 2122 ICTM=1,NCTM(ICFD)
            WRITE(ITMSS1(2+ICTM),'(2A,A12)')CHAR(96+ICTM),' ',
     &      CONTMNAM(ICTM)
 2122     CONTINUE
          ITMSS1(NCTM(ICFD)+3) ='# none of the above'
          ITMSS1(NCTM(ICFD)+4) =' ------------------------'
          ITMSS1(NCTM(ICFD)+5) ='? help                   '
          ITMSS1(NCTM(ICFD)+6) ='- exit                   '
          NITMSS1=NCTM(ICFD)+6
 2121     IIO=-1

C Help text for this menu.
          CALL EMENU('Moisture is ... ',ITMSS1,NITMSS1,IIO)
          IF(IIO.EQ.NCTM(ICFD)+6)THEN
            IF(JHUMINDX(ICFD).EQ.0)THEN
              call edisp(iuout,'None of the contaminants is moisture')
              GOTO 2121
            ENDIF
          ELSEIF(IIO.EQ.NCTM(ICFD)+5)THEN
            call edisp(iuout,'Choose a moisture contaminant.')
            GOTO 2121
          ELSEIF(IIO.EQ.NCTM(ICFD)+3)THEN
            JHUMINDX(ICFD)=0
          ELSEIF(IIO.GT.2.AND.IIO.LE.NCTM(ICFD)+2)THEN
            JHUMINDX(ICFD)=IIO-2
          ELSE
            GOTO 2121
          ENDIF

C Take contaminant name from contaminant network
          DO 789 ICTM=1,NCTM(ICFD)
            POLNAM(ICTM,ICFD)=CONTMNAM(ICTM)

C Ask for information to calculate Schmidt number
            helptopic='cfd_contam_distribution'
            call gethelptext(helpinsub,helptopic,nbhelp)
            CALL EASKABC(
     &        'How to specify turbulent Schmidt #(see help)?',
     &        POLNAM(ICTM,ICFD),'parameters','manual','default',
     &        ISCHM,nbhelp)
            IF(ISCHM.EQ.1)THEN
              ISCHMT(ICTM)=1
              CALL EASKS(HOLD,
     &          'enter mol mass(g/mol),crit vol(cm^3/mol)',
     &          'boiling pt (K) and crit pt (K)',32,
     &          ' 44.0 94.0 217.0 304.0',' ',IER,nbhelp)
              K=0
              CALL EGETWR(HOLD,K,VL1,1.,1000.,'W','molar mass',IER)
              CALL EGETWR(HOLD,K,VL2,1.,1000.,'W','crit vol',IER)
              CALL EGETWR(HOLD,K,VL3,100.,500.,'W','T_boil',IER)
              CALL EGETWR(HOLD,K,VL4,100.,1000.,'W','T_crit',IER)
              GFM(ICTM)=VL1
              VCRIT(ICTM)=VL2; TBOIL(ICTM)=VL3; TCRIT(ICTM)=VL4
            ELSEIF(ISCHM.EQ.2)THEN
              ISCHMT(ICTM)=0
              call EASKR(VAL,'What is Turbulent Schmidt Number for',
     &          POLNAM(ICTM,ICFD),0.1,'W',2.5,'W',1.,'schmidt #',
     &          IER,nbhelp)
              SCHMT(ICTM)=VAL
            ELSEIF(ISCHM.EQ.3)THEN
              ISCHMT(ICTM)=0; SCHMT(ICTM)=1.0
              IF(JHUMINDX(ICFD).EQ.ICTM)SCHMT(ICTM)=0.59
            ENDIF
 789      CONTINUE
          helptopic='cfd_contam_sinks'
          call gethelptext(helpinsub,helptopic,nbhelp)
          DOK=.TRUE.
          CALL ASKOK('Link contaminant network source/sink with this',
     &      'CFD source type boundary condition now ?',OK,dok,nbhelp)
          IF(OK)THEN
            DO 3040 ICTM=1,NCTM(ICFD)
              IX=2
              ITMSS2(1)= ' Source / sink name:'
              DO 3050 ISPMNO=1,NSPMNO
                IF(SSLINK1(ISPMNO,ICTM).NE.0)THEN
                  WRITE(ITMSS2(IX),'(A,1X,A)')CHAR(96+IX-1),
     &            SSNAME(SSLINK1(ISPMNO,ICTM))
                  IX=IX+1
                ENDIF
 3050         CONTINUE
              ITMSS2(IX)  = ' '
              ITMSS2(IX+1)= '0 none'
              ITMSS2(IX+2)= ' ------------------'
              ITMSS2(IX+3)= '? help'
              ITMSS2(IX+4)= '- exit'
              NITMSS2=IX+4
              WRITE(outs,'(2A)')  'Sources and sinks for contaminant ',
     &                       CONTMNAM(ICTM)
              call edisp(iuout,outs)
              call edisp(iuout,
     &          'are displayed. Choose which one represents this')
              call edisp(iuout,
     &          'source type CFD boundary condition (if at all)')
 3041         IXO=-1
              WRITE(WORD,'(2A)')'Choose source type for contaminant ',
     &                       CONTMNAM(ICTM)
              IF(IX.EQ.2)GOTO 3042
              CALL EMENU(WORD,ITMSS2,NITMSS2,IXO)
              IF(IXO.EQ.NITMSS2)THEN
                CONTINUE
              ELSEIF(IXO.EQ.NITMSS2-1)THEN
                helptopic='cfd_contam_sinks'
                call gethelptext(helpinsub,helptopic,nbhelp)
                CALL PHELPD('Choose source/sink',nbhelp,'-',0,0,IER)
                GOTO 3041
              ELSEIF(IXO.EQ.NITMSS2-3)THEN
                SRCE(ISB,ICTM,ICFD)=' '
              ELSEIF(IXO.GT.1.AND.IXO.LT.NITMSS2-3)THEN
                SRCE(ISB,ICTM,ICFD)=SSNAME(SSLINK1(IXO-1,ICTM))
              ELSE
                GOTO 3041
              ENDIF
 3040       CONTINUE
            GOTO 3044
 3042       call edisp(iuout,
     &        'No sources/sinks defined for this contaminant')
 3044       continue
          ENDIF
        ENDIF
        GOTO 144

C Define CFD grid cells for source type boundary condition
      elseif (INO.eq.5) then
        helptopic='cfd_source_boundary_edit'
        call gethelptext(helpinsub,helptopic,nbhelp)
        t20='Source definition'
        if(IVOLF(ISB,ICFD).ne.0)IFACE=IVOLF(ISB,ICFD)
        call PIKCELS(t20,VOLNAME(ISB,ICFD),Ii,If,Ji,Jf,Ki,Kf,
     &    IFACE,IER)
        IVOLF(ISB,ICFD)=IFACE
C        write(6,*) 'face f',ISB,ICFD,IFACE
        IVCELLS(ISB,ICFD,1)=Ii+1; IVCELLS(ISB,ICFD,2)=If+1
        JVCELLS(ISB,ICFD,1)=Ji+1; JVCELLS(ISB,ICFD,2)=Jf+1
        KVCELLS(ISB,ICFD,1)=Ki+1; KVCELLS(ISB,ICFD,2)=Kf+1
      elseif (INO.eq.7) then
        helptopic='cfd_source_boundary_edit'
        call gethelptext(helpinsub,helptopic,nbhelp)
        VAL=VOLHeat(ISB,ICFD)
        call EASKR(VAL,'Heat flux (W) in source volume?',
     &    ' ',-1000.,'W',1000.,'W',10.,'heat flux',IER,nbhelp)
        if (IER.eq.0) VOLHeat(ISB,ICFD)=VAL
      elseif (INO.eq.8) then

C<< Should ideally popup a list of casual gain types from zone
C operations file here with some regard for standalone running of dfs
        helptopic='cfd_source_boundary_edit'
        call gethelptext(helpinsub,helptopic,nbhelp)
        IVAL=IDCasgn(ISB,ICFD)
        call EASKI(IVAL,'Which casual gain index represents heat ',
     &  'generation?',0,'W',3,'W',1,'casgn index',IERI,nbhelp)
        if(ieri.eq.-3) goto 7
        if (IER.eq.0) IDCasgn(ISB,ICFD)=IVAL
      elseif (INO.eq.9) then
        helptopic='cfd_source_boundary_edit'
        call gethelptext(helpinsub,helptopic,nbhelp)
        VAL=FCasgn(ISB,ICFD)
        call EASKR(VAL,'What fraction of casual gain is generated? ',
     &  ' ',0.,'W',1.,'W',1.,'casgn fraction',IER,nbhelp)
        if (IER.eq.0) FCasgn(ISB,ICFD)=VAL
      elseif (INO.eq.13+NCTM(ICFD)) then
        helptopic='cfd_contam_sinks'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL EASKABC('Add/Delete contaminant source',' ','add',
     &      'delete','cancel',IOWL,nbhelp)
        IF(IOWL.EQ.1)THEN

C Contaminants currently only work for CFD domain 1.
C If the user tries to input contaminants in any other domain,
C inform them and do not allow it.
          IF(ICFD.NE.1)GOTO 756

          IF(ICTDFAF.EQ.1)THEN
            helptopic='cfd_contam_delete_warn'
            call gethelptext(helpinsub,helptopic,nbhelp)
            CALL PHELPD('Not allowed',nbhelp,'-',0,0,IER)
          ELSE
            helptopic='cfd_contam_sinks'
            call gethelptext(helpinsub,helptopic,nbhelp)
            SRCNAM=' '
            d12='Contaminant'
            call EASKS(SRCNAM,'Contaminant name? ','  ',12,d12,
     &      'cntmnt name',IER,nbhelp)
            NCTM(ICFD)=NCTM(ICFD)+1; ICTM=NCTM(ICFD)
            POLNAM(NCTM(ICFD),ICFD)=SRCNAM
            call EASKR(VAL,'What is source strength (kg/s) ',
     &      ' ',0.,'W',10.,'W',1.,'src strenght',IER,nbhelp)
            VOLPOL(NCTM(ICFD),ISB,ICFD)=VAL
            CALPOL(NCTM(ICFD),ICFD)=.TRUE.
            CALL ASKOK(' ','Is this contaminant water',OK,dok,0)
            IF(OK)JHUMINDX(ICFD)=ICTM

C Ask for information to calculate Schmidt number
            helptopic='cfd_contam_schmidt_no'
            call gethelptext(helpinsub,helptopic,nbhelp)
            CALL EASKABC('How to specify turbulent Schmidt #(see help)?'
     &       ,' ','parameters','manual','default',ISCHM,nbhelp)
            IF(ISCHM.EQ.1)THEN
              ISCHMT(ICTM)=1
              CALL EASKS(HOLD,'enter mol mass(g/mol),crit vol(cm^3/mol)'
     &        ,'boiling pt (K) and crit pt (K)',32
     &        ,' 44.0 94.0 217.0 304.0',' ',IER,nbhelp)
              K=0
              CALL EGETWR(HOLD,K,VL1,1.,1000.,'W','molar mass',IER)
              CALL EGETWR(HOLD,K,VL2,1.,1000.,'W','crit vol',IER)
              CALL EGETWR(HOLD,K,VL3,100.,500.,'W','T_boil',IER)
              CALL EGETWR(HOLD,K,VL4,100.,1000.,'W','T_crit',IER)
              GFM(ICTM)=VL1
              VCRIT(ICTM)=VL2; TBOIL(ICTM)=VL3; TCRIT(ICTM)=VL4
            ELSEIF(ISCHM.EQ.2)THEN
              ISCHMT(ICTM)=0
              call EASKR(VAL,'What is Turbulent Schmidt Number ',
     &           ' ',0.1,'W',2.5,'W',1.,'schmidt #',IER,nbhelp)
              SCHMT(ICTM)=VAL
            ELSEIF(ISCHM.EQ.3)THEN
              ISCHMT(ICTM)=0; SCHMT(ICTM)=1.0
              IF(JHUMINDX(ICFD).EQ.ICTM)SCHMT(ICTM)=0.59
            ENDIF
          ENDIF

C Delete Contaminant
        ELSEIF(IOWL.EQ.2)THEN
          IF(ICTDFAF.EQ.1)THEN
            helptopic='cfd_contam_delete_warn'
            call gethelptext(helpinsub,helptopic,nbhelp)
            CALL PHELPD('Not allowed',nbhelp,'-',0,0,IER)
          ELSE
            WRITE(ITUM(1),'(a)')' Contaminant name: '
            DO 132 ICTM=1,NCTM(ICFD)
              WRITE(ITUM(1+ICTM),'(3A)')
     &        CHAR(96+ICTM),' ',POLNAM(ICTM,ICFD)
 132        CONTINUE
            NITUM=1+NCTM(ICFD)
            INUO=-1
            CALL EMENU('Choose contaminant to delete',ITUM,NITUM,INUO)
            DO 133 ICTM=INUO-1,NCTM(ICFD)-1
              POLNAM(ICTM,ICFD)=POLNAM(ICTM+1,ICFD)
              CALPOL(ICTM,ICFD)=CALPOL(ICTM+1,ICFD)
              SCHMT(ICTM)=SCHMT(ICTM+1)
              ISCHMT(ICTM)=ISCHMT(ICTM+1)
              GFM(ICTM)=GFM(ICTM+1)
              VCRIT(ICTM)=VCRIT(ICTM+1)
              TBOIL(ICTM)=TBOIL(ICTM+1)
              TCRIT(ICTM)=TCRIT(ICTM+1)

C Reset humidity / water index if applicable
              IF(JHUMINDX(ICFD).EQ.INUO)THEN
                JHUMINDX(ICFD)=0
              ELSEIF(JHUMINDX(ICFD).GT.INUO)THEN
                JHUMINDX(ICFD)=JHUMINDX(ICFD)-1
              ENDIF

C Delete this contaminant from all key volumes
              DO 135 I=1,NVOL(ICFD)
                VOLPOL(ICTM,I,ICFD)=VOLPOL(ICTM+1,I,ICFD)
 135          CONTINUE
 133        CONTINUE
            POLNAM(NCTM(ICFD),ICFD)=' '
            CALPOL(NCTM(ICFD),ICFD)=.FALSE.
            DO 134 I=1,NVOL(ICFD)
              IF(IVTYPE(I,ICFD).EQ.20)VOLPOL(NCTM(ICFD),I,ICFD)=0.
 134        CONTINUE
            NCTM(ICFD)=NCTM(ICFD)-1
          ENDIF
        ENDIF
        IF(ICTDFAF.NE.1)GOTO 144

C Edit contaminants 
      elseif (INO.gt.11.and.INO.lt.12+NCTM(ICFD)) then
        ICTM=INO-11
        IF(ICTDFAF.EQ.1)THEN
          helptopic='cfd_link_contam_sink'
          call gethelptext(helpinsub,helptopic,nbhelp)
          CALL PHELPD('Contaminants Edit',10,'-',0,0,IER)
        ELSE
          helptopic='cfd_link_contam_sink'
          call gethelptext(helpinsub,helptopic,nbhelp)
          SRCNAM=POLNAM(ICTM,ICFD)
          d12='Contaminant'
          call EASKS(SRCNAM,'Contaminant name? ','  ',12,d12,
     &    'cntmnt name',IER,nbhelp)
          POLNAM(ICTM,ICFD)=SRCNAM
        ENDIF
        
        VAL=VOLPOL(ICTM,ISB,ICFD)
        call EASKR(VAL,'What is source strength (kg/s) ',
     &    ' ',0.,'W',10.,'W',1.,'src strenght',IER,nbhelp)
        VOLPOL(ICTM,ISB,ICFD)=VAL
        CALL ASKOK(' ','Is this contaminant water(humidity)',
     &    OK,dok,nbhelp)
        IF(OK)THEN
          JHUMINDX(ICFD)=ICTM
        ELSE
          IF(JHUMINDX(ICFD).EQ.ICTM)JHUMINDX(ICFD)=0
        ENDIF

C Ask for information to calculate Schmidt number
        helptopic='cfd_contam_schmidt_no'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL EASKABC('How to specify turbulent Schmidt # (see help) ?'
     &     ,' ','parameters','manual','default',ISCHM,nbhelp)
        IF(ISCHM.EQ.1)THEN
          ISCHMT(ICTM)=1
          CALL EASKS(HOLD,'enter mol mass (g/mol) crit vol (cm^3/mol)'
     &      ,'boiling pt (K) and crit pt (K)',32
     &      ,' 44.0 94.0 217.0 304.0',' ',IER,nbhelp)
          K=0
          CALL EGETWR(HOLD,K,VL1,1.,1000.,'W','molar mass',IER)
          CALL EGETWR(HOLD,K,VL2,1.,1000.,'W','crit vol',IER)
          CALL EGETWR(HOLD,K,VL3,100.,500.,'W','T_boil',IER)
          CALL EGETWR(HOLD,K,VL4,100.,1000.,'W','T_crit',IER)
          GFM(ICTM)=VL1
          VCRIT(ICTM)=VL2; TBOIL(ICTM)=VL3; TCRIT(ICTM)=VL4
        ELSEIF(ISCHM.EQ.2)THEN
          ISCHMT(ICTM)=0
          call EASKR(VAL,'What is Turbulent Schmidt Number ',
     &         ' ',0.5,'W',1.5,'W',1.,'schmidt #',IER,nbhelp)
          SCHMT(ICTM)=VAL
        ELSEIF(ISCHM.EQ.3)THEN
          ISCHMT(ICTM)=0; SCHMT(ICTM)=1.0
          IF(JHUMINDX(ICFD).EQ.ICTM)SCHMT(ICTM)=0.59
        ENDIF

C Now get source sink information from contaminants network
        IF(ICTDFAF.EQ.1)THEN
          helptopic='cfd_link_contam_sink'
          call gethelptext(helpinsub,helptopic,nbhelp)
          DOK=.FALSE.
          CALL ASKOK('Link contaminant network source/sink with this',
     &      'CFD source type boundary condition now ?',OK,dok,nbhelp)
          IF(OK)THEN
            IX=2
            ITMSS3(1)= ' Source / sink name:'
            DO 2050 ISPMNO=1,NSPMNO
              IF(SSLINK1(ISPMNO,ICTM).NE.0)THEN
                WRITE(ITMSS3(IX),'(A,1X,A)')CHAR(96+IX-1),
     &          SSNAME(SSLINK1(ISPMNO,ICTM))
                IX=IX+1
              ENDIF
 2050       CONTINUE
            ITMSS3(IX)  = ' '
            ITMSS3(IX+1)= '0 none'
            ITMSS3(IX+2)= ' ------------------'
            ITMSS3(IX+3)= '? help'
            ITMSS3(IX+4)= '- exit'
            NITMSS3=IX+4
            WRITE(outs,'(2A)')  'Sources and sinks for contaminant ',
     &                       CONTMNAM(ICTM)
            call edisp(iuout,outs)
            call edisp(iuout,
     &        'are displayed. Choose which one represents this')
            call edisp(iuout,
     &        'source type CFD boundary condition (if at all)')
 2041       IXO=-1
            WRITE(WORD,'(2A)')'Choose source type for contaminant ',
     &                     CONTMNAM(ICTM)
            IF(IX.EQ.2)GOTO 3043
            CALL EMENU(WORD,ITMSS3,NITMSS3,IXO)
            IF(IXO.EQ.NITMSS3)THEN
              CONTINUE
            ELSEIF(IXO.EQ.NITMSS3-1)THEN
              helptopic='cfd_contam_sinks'
              call gethelptext(helpinsub,helptopic,nbhelp)
              CALL PHELPD('Choose source/sink',nbhelp,'-',0,0,IER)
              GOTO 2041
            ELSEIF(IXO.EQ.NITMSS3-3)THEN
                SRCE(ISB,ICTM,ICFD)=' '
            ELSEIF(IXO.GT.1.AND.IXO.LT.NITMSS3-3)THEN
              SRCE(ISB,ICTM,ICFD)=SSNAME(SSLINK1(IXO-1,ICTM))
            ELSE
              GOTO 2041
            ENDIF
          ENDIF
          GOTO 3045
 3043     call edisp(iuout,
     &      'No sources/sinks defined for this contaminant')
 3045     continue
        ENDIF
      elseif (INO.eq.(NITMS-1)) then

C Create help for 'sources edit' menu
        helptopic='cfd_source_boundary_edit'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL PHELPD('sources edit help',nbhelp,'-',0,0,IER)
      elseif (INO.eq.NITMS) then

C Return.
        return
      else
        goto 7
      endif
        goto 5

C Prevent user from entering contaminants.
  756 write(H(1),'(a)')'Currently contaminants only work correctly in'
      write(H(2),'(a)')'CFD domain 1. This limitation will probably be'
      write(H(3),'(a)')'addressed in a future release of ESP-r.'
      CALL PHELPD('Contaminants',3,'-',0,0,IER)
      goto 144
      end


C ************************** EDBLK **************************
C EDBLK controls the editing of CFD blockages to flow
C IER=0 indicates no error.

      SUBROUTINE EDBLK(ISB,IER)
#include "cfd.h"

      COMMON/ICFNOD/ICFD,ICP
      COMMON/SPAD/MMOD,LIMIT,LIMTTY

      common/KEYVOLS/NVOL(MNZ),IVOLF(MNVLS,MNZ),IVCELLS(MNVLS,MNZ,2),
     &               JVCELLS(MNVLS,MNZ,2),KVCELLS(MNVLS,MNZ,2)
      common/KEYVOLN/VOLNAME(MNVLS,MNZ),VCsurf(MNVLS,MNZ)
      common/KEYVDAT/IVTYPE(MNVLS,MNZ),VOLTemp(MNVLS,MNZ),
     &          VOLHeat(MNVLS,MNZ),IVConfl(MNVLS,MNZ),VOLHum(MNVLS,MNZ),
     &          VOLCO2(MNVLS,MNZ),VOLVel(MNVLS,MNZ),VOLDir(MNVLS,MNZ,2),
     &          VOLArea(MNVLS,MNZ),VOLPres(MNVLS,MNZ),
     &          VOLPol(MCTM,MNVLS,MNZ)
      COMMON/GEOM/XP(ntcelx),YP(ntcely),ZP(ntcelz),
     1            DXEP(ntcelx),DXPW(ntcelx),DYNP(ntcely),DYPS(ntcely),
     2            DZHP(ntcelz),DZPL(ntcelz),
     3            SEW(ntcelx),SNS(ntcely),SHL(ntcelz),
     4            XU(ntcelx),YV(ntcely),ZW(ntcelz)
      common/param2/TITLE(MNZ),CFTRFL(MNZ),LPHI(MNZ)
      CHARACTER*72 TITLE,CFTRFL,LPHI

      common/blksso/NBLK(MNZ),INBLK(MNVLS,MNZ),NSSO(MNZ),
     &          INSSO(MNVLS,MNZ)
      COMMON/CFDVIS/HAS_GEOM,ISHSB,ISHAO,IFACES,ISHBLK,ISHGEO,INITD
      logical HAS_GEOM
      character INITD*6
 
      character ITMSS(14)*50,FACE*6
      character ltmp*16,t12*12,d12*12,t20*20
      character VOLNAME*12, VCsurf*12
      character zngeo*3,etext*82

      logical ERROR
      integer NITMS,INO ! max items and current menu item

#ifdef OSI
      integer iside,isize,ifont
#else
      integer*8 iside,isize,ifont
#endif

C For help messages
      character helpinsub*24 ! subroutine name
      character helptopic*24 ! string (unique) for topic
      integer nbhelp     ! number of help lines found

      helpinsub='EDBLK'  ! set for subroutine
      helptopic='cfd_blockage_menu'
      call gethelptext(helpinsub,helptopic,nbhelp)

C Adjust each cell by -1 as dfs adds fictitious cells around the domain.
      Ii=IVCELLS(ISB,ICFD,1)-1; If=IVCELLS(ISB,ICFD,2)-1
      Ji=JVCELLS(ISB,ICFD,1)-1; Jf=JVCELLS(ISB,ICFD,2)-1
      Ki=KVCELLS(ISB,ICFD,1)-1; Kf=KVCELLS(ISB,ICFD,2)-1

C Check if new blockage.
      if (Ii.lt.0.or.If.lt.0.or.
     &    Ji.lt.0.or.Jf.lt.0.or.
     &    Ki.lt.0.or.Kf.lt.0) then 
        t20='Blockage'
        if(IVOLF(ISB,ICFD).ne.0)IFACE=IVOLF(ISB,ICFD)
        call PIKCELS(t20,VOLNAME(ISB,ICFD),Ii,If,Ji,Jf,Ki,Kf,
     &    IFACE,IER)
        IVOLF(ISB,ICFD)=IFACE
C        write(6,*) 'face g',ISB,ICFD,IFACE
        IVCELLS(ISB,ICFD,1)=Ii+1; IVCELLS(ISB,ICFD,2)=If+1
        JVCELLS(ISB,ICFD,1)=Ji+1; JVCELLS(ISB,ICFD,2)=Jf+1
        KVCELLS(ISB,ICFD,1)=Ki+1; KVCELLS(ISB,ICFD,2)=Kf+1
      endif


C Generate text for location index. (not yet used in the interface)
 5    FACE='  '
      if (IVOLF(ISB,ICFD).eq.1) then
        FACE='West'
      elseif (IVOLF(ISB,ICFD).eq.2) then
        FACE='East'
      elseif (IVOLF(ISB,ICFD).eq.3) then
        FACE='South'
      elseif (IVOLF(ISB,ICFD).eq.4) then
        FACE='North'
      elseif (IVOLF(ISB,ICFD).eq.5) then
        FACE='Low'
      elseif (IVOLF(ISB,ICFD).eq.6) then
        FACE='High'
      elseif (IVOLF(ISB,ICFD).eq.7) then
        FACE='Whole'
      elseif (IVOLF(ISB,ICFD).eq.8) then
        FACE='Block'
      elseif (IVOLF(ISB,ICFD).eq.9) then
        FACE='Source'
      endif

C Set boundary type string.
      ltmp='  '
      if (IVTYPE(ISB,ICFD).eq.30) then
        ltmp='Blockage'
      endif

C Create a menu showing definitions for different conflation schemes.
      xd=XU(If+2)-XU(Ii+1)  ! need to add 1 to indicies
      yd=YV(Jf+2)-YV(Ji+1)
      zd=ZW(Kf+2)-ZW(Ki+1)

C Debug echo the coordinates of the blockage.
      write(6,*) 'cords X are:',XU(Ii+1),XU(If+2),xd
      write(6,*) 'cords Y are:',YV(Ji+1),YV(Jf+2),yd
      write(6,*) 'cords Z are:',ZW(Ki+1),ZW(Kf+2),zd

      write (ITMSS(1),'(a,a)')   'a Name: ',VOLNAME(ISB,ICFD)
      write (ITMSS(2),'(3a)')   '  Type: ',ltmp
      write (ITMSS(3),'(a)')     ' -----------------------------'
      write (ITMSS(4),'(a,a)')   '         Is  If  Js  Jf  Ks  Kf '
      write (ITMSS(5),'(a,6(i4))') 'c Cells:',Ii,If,Ji,Jf,Ki,Kf
      write (ITMSS(6),'(a)')     ' -----------------------------'
      write (ITMSS(7),'(a,f6.3,a,f6.3,a,f5.3)') '  Z cords: ',XU(Ii+1),
     &                  ' to ',XU(If+2),' delta ',xd
      write (ITMSS(8),'(a,f6.3,a,f6.3,a,f5.3)') '  Y cords: ',YV(Ji+1),
     &                  ' to ',YV(Jf+2),' delta ',yd
      write (ITMSS(9),'(a,f6.3,a,f6.3,a,f5.3)') '  Z cords: ',ZW(Ki+1),
     &                  ' to ',ZW(Kf+2),' delta ',zd
      write(ITMSS(10),'(a)')     ' -----------------------------'
      write(ITMSS(11),'(a,f6.2)')'e Heat flux (W/m^2):',
     &                  VOLHeat(ISB,ICFD)
      ITMSS(12) =' -----------------------------'
      ITMSS(13)='? Help'
      ITMSS(14)='- exit'
      NITMS=14

      IF(MMOD.EQ.8)THEN
        do 639 I=1,NBLK(ICFD)
          if (ISB.eq.INBLK(I,ICFD)) then
            ISHBLK_TMP=ISHBLK
            ISHGEO_TMP=ISHGEO
            ISHBLK=I
            ISHGEO=0
            CALL redraw(IER)
            ISHBLK=ISHBLK_TMP
            ISHGEO=ISHGEO_TMP
          endif
639     continue   
      ENDIF

C Help text for this menu.
      helptopic='cfd_blockage_menu'
      call gethelptext(helpinsub,helptopic,nbhelp)

C Display menu.
      INO=-1
      CALL EMENU('Blockage edit',ITMSS,NITMS,INO)
      if (INO.eq.1) then

C Get boundary name.
        t12=VOLNAME(ISB,ICFD)
        d12=' '
        call EASKS(t12,'Blockage name? ','  ',12,d12,'blk name',
     &    IER,nbhelp)
        if (IER.eq.0) then
          call st2name(t12,d12)
          VOLNAME(ISB,ICFD)=d12
        endif
      elseif (INO.eq.5) then
 79     t20='Blockage'
        if(IVOLF(ISB,ICFD).ne.0)IFACE=IVOLF(ISB,ICFD)
        call PIKCELS(t20,VOLNAME(ISB,ICFD),Ii,If,Ji,Jf,Ki,Kf,
     &    IFACE,IER)
        IVOLF(ISB,ICFD)=IFACE
C        write(6,*) 'face h',ISB,ICFD,IFACE
        IVCELLS(ISB,ICFD,1)=Ii+1; IVCELLS(ISB,ICFD,2)=If+1
        JVCELLS(ISB,ICFD,1)=Ji+1; JVCELLS(ISB,ICFD,2)=Jf+1
        KVCELLS(ISB,ICFD,1)=Ki+1; KVCELLS(ISB,ICFD,2)=Kf+1

C Check for illegal volume (vol must be 1 cells thick in all directions).
        ERROR=.false.
        if (IVCELLS(ISB,ICFD,1).gt.IVCELLS(ISB,ICFD,2)) then
          call usrmsg('Blockage cell order reversed',
     &                'in the x direction.','W')
          ERROR=.true.
        endif
        if (JVCELLS(ISB,ICFD,1).gt.JVCELLS(ISB,ICFD,2)) then
          call usrmsg('Blockage cell order reversed',
     &                'in the y direction.','W')
          ERROR=.true.
        endif
        if (KVCELLS(ISB,ICFD,1).gt.KVCELLS(ISB,ICFD,2)) then
          call usrmsg('Blockage cell order reversed',
     &                'in the z direction.','W')
          ERROR=.true.
        endif
        if (ERROR) goto 79
      elseif (INO.eq.11) then
        VAL=VOLHeat(ISB,ICFD)
        call EASKR(VAL,'Heat flux (W/m^2) on blockage surfaces?',
     &    ' ',-1000.,'W',1000.,'W',10.,'heat flux',IER,nbhelp)
        if (IER.eq.0) VOLHeat(ISB,ICFD)=VAL
      elseif (INO.eq.(NITMS-1)) then

C Help.
        helptopic='cfd_blockage_menu'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL PHELPD('blocakage edit',nbhelp,'-',0,0,IER)
      elseif (INO.eq.NITMS) then

C Return.
        return
      endif

      goto 5
      end


C ************************** HANDSHK **************************
C HANDSHK controls the editing of the handshaking mechanism.

      SUBROUTINE HANDSHK(IHST,IER)
#include "cfd.h"

      COMMON/ICFNOD/ICFD,ICP
      common/METHDS/ITURB(MNZ),IBUOY(MNZ)

      dimension IHSID(13)
      character ITMSS(13)*60
      integer NITMS,INO ! max items and current menu item

C Set up menu - display only applicable handshaking mechanisms based on 
C current turbulence modelling choice.
C Return if not k-e or MIT.
      IER=0
      if (ITURB(ICFD).ne.1.AND.ITURB(ICFD).ne.3) then
        call usrmsg('Building surface boundary conditions are only',
     &             'applicable for k-e or MIT turbulence models. ','W')
        IER=1
        return
      endif
      ITMSS(1)=   'Direction|Wall function|HTCs / Qsurf     |Tref'
      if (ITURB(ICFD).eq.1) then
        ITMSS(2)= 'a One-way|log-law n/a  |CFD              |n/a'
        IHSID(2)=1
        ITMSS(3)= 'b One-way|Yuan    n/a  |CFD              |n/a'
        IHSID(3)=3
        ITMSS(4)= 'c One-way|Yuan    n/a  |A*HC*(Tsurf-Tref)|BSim '
        IHSID(4)=4
        ITMSS(5)= 'd One-way|Yuan    n/a  |A*HC*(Tsurf-Tref)|CFD ave.'
        IHSID(5)=5
        ITMSS(6)= 'e One-way|log-law n/a  |A*HC*(Tsurf-Tref)|BSim'
        IHSID(6)=6
        ITMSS(7)= 'f One-way|log-law n/a  |A*HC*(Tsurf-Tref)|CFD ave.'
        IHSID(7)=7
        ITMSS(8)= 'g One-way|log-law local|A*HC*(Tsurf-Tref)|CFD local'
        IHSID(8)=8
        ITMSS(9)= 'h Two-way|log-law n/a  |CFD              |n/a'
        IHSID(9)=9
        ITMSS(10)='i Two-way|Yuan    n/a  |CFD              |n/a'
        IHSID(10)=11
        ITMSS(11)='j Two-way|Yuan    local|A*HC*(Tsurf-Tref)|CFD local'
        IHSID(11)=12
        ITMSS(12)='k Two-way|log-law n/a  |CFD              |n/a'
        IHSID(12)=13
        ITMSS(12)='l Two-way|log-law local|A*HC*(Tsurf-Tp)  |local'
        IHSID(12)=14
        NITMS=12
      else
        ITMSS(2)= 'a One-way|n/a n/a      |CFD              |n/a'
        IHSID(2)=2
        ITMSS(3)= 'b Two-way|n/a n/a      |CFD              |n/a'
        IHSID(3)=10
        NITMS=3
      endif

C Display menu and then decode returned value.
      CALL EMENU('Handshaking mechanism',ITMSS,NITMS,INO)
      if (INO.eq.0) then
        ier=1
        return
      else
        IHST=IHSID(INO)
      endif

      return
      end


C ************************** PIKCELS **************************
C PIKCELS controls the editing of the cells defining a boundary region.

      SUBROUTINE PIKCELS(BTYPE,BNAME,Ii,If,Ji,Jf,Ki,Kf,IFACE,IER)

#include "cfd.h"
      integer iCountWords

      COMMON/OUTIN/IUOUT,IUIN
      common/grdmax/NTCX,NTCY,NTCZ

C XP YP ZP are the centre of the cells.
C XU YV ZW are the start of the cells.
      COMMON/GEOM/XP(ntcelx),YP(ntcely),ZP(ntcelz),
     1            DXEP(ntcelx),DXPW(ntcelx),DYNP(ntcely),DYPS(ntcely),
     2            DZHP(ntcelz),DZPL(ntcelz),
     3            SEW(ntcelx),SNS(ntcely),SHL(ntcelz),
     4            XU(ntcelx),YV(ntcely),ZW(ntcelz)
      COMMON/ALL/NI,NJ,NK,NIM1,NJM1,NKM1,NIM2,NJM2,NKM2

      character BTYPE*20,BNAME*12,s30*30,dflt*30,lface*7,t30*30
      character outs*124
      integer IFACE   ! for radio button

C For help messages
      character helpinsub*24 ! subroutine name
      character helptopic*24 ! string (unique) for topic
      integer nbhelp     ! number of help lines found

      helpinsub='PIKCELS'  ! set for subroutine

 5    IER=0
      t30='  '
      
C Ask for face first if necessary.
      helptopic='cfd_cell_selection'
      call gethelptext(helpinsub,helptopic,nbhelp)
      if (IFACE.eq.1) then
        lface='West'
      elseif (IFACE.eq.2) then
        lface='East'
      elseif (IFACE.eq.3) then
        lface='South'
      elseif (IFACE.eq.4) then
        lface='North'
      elseif (IFACE.eq.5) then
        lface='Low'
      elseif (IFACE.eq.6) then
        lface='High'
      else
        lface='UNKNOWN'
      endif
      if (IFACE.lt.7) then
        write (t30,'(2a)') 'Current face: ',lface
        IFACE=1
        CALL EASKATOG('Boundary face is on which side?',t30,'West',
     &    'East','South','North','Low','High',' ',IFACE,nbhelp)
      endif
      IONE=1

C Remind the user of the cell coordinates for the particular
C fact they they selected.
      XZ=0.; YZ=0.; ZZ=0.
      if (IFACE.eq.1.or.IFACE.eq.2) then

C List out the Y-Z coordinates for west face and east face.
        call edisp(iuout,'Grid in the Y axis... ')
        call edisp(iuout,'Cell index, Y start & end,   distance(m)')
        DO 102 I=2,(NJ-1)
          yd=YV(I+1)-YV(I)
          write(outs,'(i7,3f10.4)') I-1,YV(I)+YZ,YV(I+1)+YZ,yd 
          call edisp(iuout,outs)
  102   CONTINUE
        call edisp(iuout,' ')
        call edisp(iuout,'Grid in the Z axis... ')
        call edisp(iuout,'Cell index, Z start & end,   distance(m)')
        DO 103 I=2,(NK-1)
          zd=ZW(I+1)-ZW(I)
          write(outs,'(i7,3f10.4)') I-1,ZW(I)+ZZ,ZW(I+1)+ZZ,zd
          call edisp(iuout,outs)
  103   CONTINUE
      elseif (IFACE.eq.3.or.IFACE.eq.4) then

C List out the X-Z coordinates for south and north face.
        call edisp(iuout,'Grid in the X axis... ')
        call edisp(iuout,'Cell index, X start & end,   distance(m)')
        DO 106 I=2,(NI-1)
          xd=XU(I+1)-XU(I)
          write(outs,'(i7,3f10.4)') I-1,XU(I)+XZ,XU(I+1)+XZ,xd
          call edisp(iuout,outs)
  106   CONTINUE
        call edisp(iuout,' ')
        call edisp(iuout,'Grid in the Z axis... ')
        call edisp(iuout,'Cell index, Z start & end,   distance(m)')
        DO 107 I=2,(NK-1)
          zd=ZW(I+1)-ZW(I)
          write(outs,'(i7,4f10.4)') I-1,ZW(I)+ZZ,ZW(I+1)+ZZ,zd
          call edisp(iuout,outs)
  107   CONTINUE
      elseif (IFACE.eq.5.or.IFACE.eq.6) then

C List out X and Y cell coordinates for low and high face.
        call edisp(iuout,'Grid in the X axis... ')
        call edisp(iuout,'Cell index, X start & end,   distance(m)')
        DO 110 I=2,(NI-1)
          xd=XU(I+1)-XU(I)
          write(outs,'(i7,3f10.4)') I-1,XU(I)+XZ,XU(I+1)+XZ,xd
          call edisp(iuout,outs)
  110   CONTINUE
        call edisp(iuout,'  ')
        call edisp(iuout,'Grid in the Y axis... ')
        call edisp(iuout,'Cell index, Y start & end,   distance(m)')
        DO 111 I=2,(NJ-1)
          yd=YV(I+1)-YV(I)
          write(outs,'(i7,3f10.4)') I-1,YV(I)+YZ,YV(I+1)+YZ,yd 
          call edisp(iuout,outs)
  111   CONTINUE
      endif

C Setup Ii and If and dflt and t30 for each of the faces.
      if (IFACE.eq.1) then

C West face.
        Ii=1; If=1
        write(dflt,'(6i4)') Ii,If,IONE,NTCY,IONE,NTCZ
        t30='(current face: west)'
      elseif (IFACE.eq.2) then

C East face.
        Ii=NTCX; If=NTCX
        write(dflt,'(6i4)') Ii,If,IONE,NTCY,IONE,NTCZ
        t30='(current face: east)'
      elseif (IFACE.eq.3) then

C South face.
        Ji=1; Jf=1
        write(dflt,'(6i4)') IONE,NTCX,Ji,Jf,IONE,NTCZ
        t30='(current face: south)'
      elseif (IFACE.eq.4) then

C North face.
        Ji=NTCY; Jf=NTCY
        write(dflt,'(6i4)') IONE,NTCX,Ji,Jf,IONE,NTCZ
        t30='(current face: north)'
      elseif (IFACE.eq.5) then

C Low face.
        Ki=1; Kf=1
        write(dflt,'(6i4)') IONE,NTCX,IONE,NTCY,Ki,Kf
        t30='(current face: low)'
      elseif (IFACE.eq.6) then

C High face.
        Ki=NTCZ; Kf=NTCZ
        write(dflt,'(6i4)') IONE,NTCX,IONE,NTCY,Ki,Kf
        t30='(current face: high)'
      else

C Unknown face.
        write(dflt,'(6i4)') IONE,NTCX,IONE,NTCY,IONE,NTCZ
        t30='(current face: UNKNOWN)'
      endif
      if (Ii.lt.0.or.If.lt.0.or.
     &    Ji.lt.0.or.Jf.lt.0.or.
     &    Ki.lt.0.or.Kf.lt.0) then 
        s30=dflt
      else

C Ask for cells.
        write(s30,'(6i4)') Ii,If,Ji,Jf,Ki,Kf
      endif
 10   call EASKS(s30,'Boundary cells (Ii If Ji Jf Ki Kf)?',
     &  t30,30,dflt,'bndry cells',IER,nbhelp)
      NITMS = iCountWords(s30)
      if (NITMS.ne.6) goto 10

C Read six cell locations.
      K=0
      call EGETWI(s30,K,IVAL,1,NTCX,'F','x cells start',IER)
      if (IER.eq.0) Ii=IVAL
      call EGETWI(s30,K,IVAL,1,NTCX,'F','x cells end',IER)
      if (IER.eq.0) If=IVAL
      call EGETWI(s30,K,IVAL,1,NTCY,'F','y cells start',IER)
      if (IER.eq.0) Ji=IVAL
      call EGETWI(s30,K,IVAL,1,NTCY,'F','y cells end',IER)
      if (IER.eq.0) Jf=IVAL
      call EGETWI(s30,K,IVAL,1,NTCZ,'F','z cells start',IER)
      if (IER.eq.0) Ki=IVAL
      call EGETWI(s30,K,IVAL,1,NTCZ,'F','z cells end',IER)
      if (IER.eq.0) Kf=IVAL

      if (IER.ne.0) goto 10
      
C Check for incompatable face/ cell definition
      if (IFACE.eq.1) then
        if (Ii.ne.1.or.If.ne.1) IER=2
      elseif (IFACE.eq.2) then
        if (Ii.ne.NTCX.or.If.ne.NTCX) IER=2
      elseif (IFACE.eq.3) then
        if (Ji.ne.1.or.Jf.ne.1) IER=2
      elseif (IFACE.eq.4) then
        if (Ji.ne.NTCY.or.Jf.ne.NTCY) IER=2
      elseif (IFACE.eq.5) then
        if (Ki.ne.1.or.Kf.ne.1) IER=2
      elseif (IFACE.eq.6) then
        if (Ki.ne.NTCZ.or.Kf.ne.NTCZ) IER=2
      endif
      if (IER.ne.0) then
        call usrmsg('Your cells are not compatible with the',
     &              'chosen face.  Please redefine.','W') 
        goto 5
      endif

      return
      end

C ************************** PIKFACE **************************
C PIKFACE controls the editing of the face of the boundary region.

      SUBROUTINE PIKFACE(BTYPE,BNAME,Ii,If,Ji,Jf,Ki,Kf,IFACE,IER)
C Check which two start/ end locations are the same.  Keep count of 
C number of matches so as confusion can be removed.

      common/grdmax/NTCX,NTCY,NTCZ

      character BTYPE*20,BNAME*12

C For help messages
      character helpinsub*24 ! subroutine name
      character helptopic*24 ! string (unique) for topic
      integer nbhelp     ! number of help lines found

      helpinsub='PIKFACE'  ! set for subroutine

 5    IER=0

C Check that face is needed.
      if (IFACE.eq.8.or.IFACE.eq.9) return

C Try to guess which face.
      IMTCH=0
      if (Ii.eq.If) then
        if (Ii.eq.1) then

C West face.
          IMTCH=IMTCH+1; IFACE=1
        elseif (Ii.eq.NTCX) then

C East face.
          IMTCH=IMTCH+1; IFACE=2
        endif
      endif
      if (Ji.eq.Jf) then
        if (Ji.eq.1) then

C South face.
          IMTCH=IMTCH+1; IFACE=3
        elseif (Ji.eq.NTCY) then

C North face.
          IMTCH=IMTCH+1; IFACE=4
        endif
      endif
      if (Ki.eq.Kf) then
        if (Ki.eq.1) then

C Base face.
          IMTCH=IMTCH+1; IFACE=5
        elseif (Ki.eq.NTCZ) then

C Top face.
          IMTCH=IMTCH+1; IFACE=6
        endif
      endif
      helptopic='cfd_face_selection'
      call gethelptext(helpinsub,helptopic,nbhelp)

C Check for conflicts.
      if (IMTCH.eq.0) then
        call usrmsg ('Cells must be at a boundary face.','  ','W')
        call PIKCELS(BTYPE,BNAME,Ii,If,Ji,Jf,Ki,Kf,IFACE,IER)
      elseif (IMTCH.gt.1) then
        call usrmsg ('Cannot resolve boundary face.',
     &               'Please choose from the following list.','W')
        IFACE=1
        CALL EASKATOG('Boundary face is on which side?',' ','West',
     &    'East','South','North','Low','High',' ',IFACE,nbhelp)
      endif

C Check for incompatable face/ cell definition
      if (IFACE.eq.1) then
        if (Ii.ne.1.or.If.ne.1) IER=2
      elseif (IFACE.eq.2) then
        if (Ii.ne.NTCX.or.If.ne.NTCX) IER=2
      elseif (IFACE.eq.3) then
        if (Ji.ne.1.or.Jf.ne.1) IER=2
      elseif (IFACE.eq.4) then
        if (Ji.ne.NTCY.or.Jf.ne.NTCY) IER=2
      elseif (IFACE.eq.5) then
        if (Ki.ne.1.or.Kf.ne.1) IER=2
      elseif (IFACE.eq.6) then
        if (Ki.ne.NTCZ.or.Kf.ne.NTCZ) IER=2
      endif
      if (IER.ne.0) then
        call usrmsg('Your cells are not compatible with the',
     &              'chosen face.  Please redefine.','W') 
        goto 5
      endif

      return
      end

C ************* TOPBOT 
C TOPBOT: Define upper and lower curvilinear x-section for use in
C gridding generation.
      SUBROUTINE TOPBOT(ICOMP,act,IER)
#include "building.h"
#include "geometry.h"
#include "prj3dv.h"
      integer iCountWords

      COMMON/SPAD/MMOD,LIMIT,LIMTTY
      COMMON/OUTIN/IUOUT,IUIN
      COMMON/gzonpik/izgfoc,nzg,nznog(mcom)

C Curvilinear boundary vertex numbers.
      common/CLVTXNO/ivxhigh(MV),ivxlow(MV),ihigh,ilow

      logical clkok,found

      DIMENSION COG1(3),COG2(3),COG3(3),JJVN(MV)
      integer ISTRW

      CHARACTER act*1,HOLD*96,outs*124
#ifdef OSI
      integer iix,iiy,iixx,iiyy,ixd,iyd,iik
#else
      integer*8 iix,iiy,iixx,iiyy,ixd,iyd,iik
#endif

C For help messages
      character helpinsub*24 ! subroutine name
      character helptopic*24 ! string (unique) for topic
      integer nbhelp     ! number of help lines found

      helpinsub='TOPBOT'  ! set for subroutine

C Tell the user the rules of selection.
      ier=0
      helptopic='cfd_curvilinear_grid'
      call gethelptext(helpinsub,helptopic,nbhelp)
      CALL PHELPD('gen start popup',nbhelp,'-',0,0,IER)

C Redraw the zone with vertices.
      MODIFYVIEW=.TRUE.
      MODBND=.TRUE.
      CALL INLNST(1)
      ITVNO=0
      itsnm=0
      nzg=1
      nznog(1)=ICOMP
      izgfoc=ICOMP
      CALL ADJVIEW(IER)

C Ask for and display the top defining vertices (west to east).
C Allow user to type in the list of vertices as a string.  The
C number of items is the number of vertices and then parse out
C the individual vertex numbers from the list.
   41 HOLD=' '
      if(MMOD.ne.8)then
        CALL EASKS(HOLD,' Give high vertices (west to east)','  ',
     &    72,' ','associated vertices',IER,nbhelp)
      else
        ISTRW=72
        CALL EASKSCMD(HOLD,' Give high vertices (west to east)',
     &   ' ','via mouse', clkok,ISTRW,' 1 2 3','assoc vert+mouse',
     &   IER,nbhelp)
        if(clkok)then
          ijvn=0
          call edisp(iuout,
     &        'Select points via cursor...type `e` to finish.')
  46      CALL trackview(iik,iixx,iiyy)
          if(iik.eq.69.or.iik.eq.101)goto 47
          found=.false.
          do 45 i=1,NTV
            COG1(1)=X(I); COG1(2)=Y(I); COG1(3)=Z(I)
            CALL VECTRN(COG1,TEMAT,COG2,IER)
            CALL VECTRN(COG2,ETSMAT,COG3,IER)
            call u2pixel(COG3(1),COG3(2),iix,iiy)
            ixd=iix-iixx; iyd=iiy-iiyy
            if(abs(ixd).lt.5.and.abs(iyd).lt.5)then
              if(found)then
                call edisp(iuout,'Close points...try again.')
                goto 46
              endif
              WRITE(outs,'(a,i5,a,i5,a,i3)')' The point @ x=',iixx,
     &          ' & y=',iiyy,' matches vertex ',i
              call edisp(iuout,outs)
              found=.true.

C << Should Check that points are in the same plane and that X is increasing. >>
              ijvn=ijvn+1; jjvn(ijvn)=i
              goto 46
            endif
  45      continue
          if(.NOT.found)goto 46
  47      if(ijvn.lt.2)then
            call edisp(iuout,'Not enough points.')
            goto 41
          endif

C This write assumes not-very-complex zones and surfaces.
C << needs to support lots more edges >>
          HOLD=' '
          WRITE(HOLD,'(32I3)')(JJVN(J),J=1,ijvn)
          CALL EASKS(HOLD,' Associated high vertices (confirm):',
     &      '  ',96,' ','associated vertices',IIER,nbhelp)
        endif
      endif
      NV = iCountWords(HOLD)
      K=0
      DO 94 J=1,NV
        CALL EGETWI(HOLD,K,JV,1,NTV,'W','vertex list',IER)
        ivxhigh(J)=JV
   94 CONTINUE
      ihigh=NV

C Ask for and display the bottom defining vertices (west to east).
   51 HOLD=' '
      if(MMOD.ne.8)then
        CALL EASKS(HOLD,' Give low vertices (west to east)','  ',
     &    72,' ','associated vertices',IER,nbhelp)
      else
        ISTRW=72
        CALL EASKSCMD(HOLD,' Give low vertices (west to east)',
     &   ' ','via mouse', clkok,ISTRW,' 1 2 3','assoc vert+mouse',
     &   IER,nbhelp)
        if(clkok)then
          ijvn=0
          call edisp(iuout,
     &        'Select points via cursor...type `e` to finish.')
  56      CALL trackview(iik,iixx,iiyy)
          if(iik.eq.69.or.iik.eq.101)goto 57
          found=.false.
          do 55 i=1,NTV
            COG1(1)=X(I); COG1(2)=Y(I); COG1(3)=Z(I)
            CALL VECTRN(COG1,TEMAT,COG2,IER)
            CALL VECTRN(COG2,ETSMAT,COG3,IER)
            call u2pixel(COG3(1),COG3(2),iix,iiy)
            ixd=iix-iixx; iyd=iiy-iiyy
            if(abs(ixd).lt.5.and.abs(iyd).lt.5)then
              if(found)then
                call edisp(iuout,'Close points...try again.')
                goto 56
              endif
              WRITE(outs,'(a,i5,a,i5,a,i3)')' The point @ x=',iixx,
     &          ' & y=',iiyy,' matches vertex ',i
              call edisp(iuout,outs)
              found=.true.
              ijvn=ijvn+1; jjvn(ijvn)=i
              goto 56
            endif
  55      continue
          if(.NOT.found)goto 56
  57      if(ijvn.lt.2)then
            call edisp(iuout,'Not enough points.')
            goto 51
          endif

C This write statement assumes not-very-complex zones and polygons jjvn. 
          HOLD=' '
          WRITE(HOLD,'(32I3)')(JJVN(J),J=1,ijvn)
          CALL EASKS(HOLD,' Associated low vertices (confirm):',
     &      '  ',96,' ','associated vertices',IIER,nbhelp)
        endif
      endif
      NV = iCountWords(HOLD)
      K=0
      DO 95 J=1,NV
        CALL EGETWI(HOLD,K,JV,1,NTV,'W','vertex list',IER)
        ivxlow(J)=JV
   95 CONTINUE
      ilow=NV
      return
      end

C ******************* GRIDDEF *******************
C GRIDDEF - define gridding via regions and specify number of cells
C and cell distribution for each region.

      SUBROUTINE GRIDDEF(IAX,NRG,NTC,TOTDIS,NCEL,REG,plaw,ierr)
#include "building.h"
#include "geometry.h"
#include "cfd.h"
#include "epara.h"

      COMMON/SPAD/MMOD,LIMIT,LIMTTY
      common/OUTIN/IUOUT,IUIN
      COMMON/MFTRA/IMFTU
      COMMON/SHOUT/ICOUT
      integer menuchw,igl,igr,igt,igb,igw,igwh
      COMMON/VIEWPX/menuchw,igl,igr,igt,igb,igw,igwh
      COMMON/ICFNOD/ICFD,ICP

      dimension NCEL(MNREG,MNZ),plaw(MNREG,MNZ)
      dimension NRG(4,MNZ),REG(MNREG,MNZ)

      logical OK,DOK,GRMODE

      character ITEM(MNREG+5)*44,KEY*1,HOLD*24,outs*124,TITLE*30

      integer icellt    ! for local editing.
      integer MVERT,INO ! max items and current menu item

#ifdef OSI
      integer iigl,iigr,iigt,iigb,iigw,iigwh
      integer iiw1,iiw2,iiw3,iiw4,iimenu
#else
      integer*8 iigl,iigr,iigt,iigb,iigw,iigwh
      integer*8 iiw1,iiw2,iiw3,iiw4,iimenu
#endif

C For help messages
      character helpinsub*24 ! subroutine name
      character helptopic*24 ! string (unique) for topic
      integer nbhelp     ! number of help lines found

      helpinsub='GRIDDEF'  ! set for subroutine

C Set additional output units to stdout. Then redirect warning 
C messages to stderr in case of rule script program control.
      IMFTU=IUOUT
      IF(MMOD.EQ.-6)then
        ICOUT=0
      else
        ICOUT=IUOUT
      endif

C Reset the display bounds for grid display.
      GRMODE=.false.
      IF(MMOD.EQ.8)THEN
        call startbuffer()

C Setup and pass in parameters to win3d.
        iiw1=6; iiw2=6; iiw3=2; iiw4=4; iimenu=menuchw
        iigl=igl; iigr=igr; iigt=igt; iigb=igb; iigw=igw; iigwh=igwh
        CALL win3d(iimenu,iiw1,iiw2,iiw3,iiw4,
     &    iigl,iigr,iigt,iigb,iigw,iigwh)
        igl=int(iigl); igr=int(iigr); igt=int(iigt); igb=int(iigb)
        igw=int(iigw); igwh=int(iigwh)
        GRMODE=.true.
      ENDIF

C Reset error flag.
      ierr=0

C Initialise menu size variables based on window size. 
C IVERT is the menu position, MVERT the current number of menu lines.
 5    INO=-4
      MHEAD=5
      MCTL=4
      ILEN=NRG(IAX,ICFD)
      IPACT=CREATE
      CALL EKPAGE(IPACT)

C Check current defined distance and total no of cells.
      CTDIS=0.0
      NTC=0
      do 8 I=1,NRG(IAX,ICFD)
        CTDIS=CTDIS+REG(I,ICFD)
        NTC=NTC+ABS(NCEL(I,ICFD))
 8    continue
      ierr=0
      if (NTC.gt.1) then
        if (abs(CTDIS-TOTDIS).gt.0.001)then
          if (IAX.eq.1) then
            write(outs,'(a,f6.3,a,f6.3,a)') 'The current length ',
     &       CTDIS,' & TOTDIS',TOTDIS,' are not close for X axis.'
          elseif(IAX.eq.2)then
            write(outs,'(a,f6.3,a,f6.3,a)') 'The current length ',
     &        CTDIS,' & TOTDIS',TOTDIS,' are not close for Y axis.'
          elseif(IAX.eq.3)then
            write(outs,'(a,f6.3,a,f6.3,a)') 'The current length ',
     &        CTDIS,' & TOTDIS',TOTDIS,' are not close for Z axis.'
          elseif(IAX.eq.4)then
            write(outs,'(a,f6.3,a,f6.3,a)') 'The current length',
     &        CTDIS,' & TOTDIS',TOTDIS,' are not close for Z eaxis.'
          endif
          call edisp(iuout,outs)
          ierr=1
        endif
      endif

C Draw grid in 1D.
      if (GRMODE) then
        xoff=0.; yoff=0.
        if (IAX.eq.1) then
          TITLE='X-axis gridding'
        elseif (IAX.eq.2) then
          TITLE='Y-axis gridding'
        elseif (IAX.eq.3) then
          TITLE='Z-axis gridding'
        elseif (IAX.eq.4) then
          TITLE='Ze-axis gridding'
        endif
        call VGRID1D(xoff,yoff,TOTDIS,CTDIS,NRG,IAX,ICFD,NCEL,REG,plaw,
     &               TITLE,'Length (m)')
      endif

C Create menu displaying data for all regions.
      write (ITEM(1),'(a,i4)')' Number of regions: ',NRG(IAX,ICFD)
      write (ITEM(2),'(a,f8.3)') ' Total domain length:    ',TOTDIS
      write (ITEM(3),'(a,f8.3)') ' Current defined length: ',CTDIS
      ITEM(4) = ' ------------------------------------'
      ITEM(5) = '   Region  Cells  Length  P-law coeff  sym'
      M=MHEAD
      do 10 L=1,ILEN
        IF(L.GE.IST.AND.(L.LE.(IST+MIFULL)))THEN
          M=M+1
          CALL EMKEY(L,KEY,IER)
          if (NCEL(L,ICFD).lt.0) then
            write(ITEM(M),'(a,2i7,2f9.3,6x,a)')KEY,L,abs(NCEL(L,ICFD)),
     &                                   REG(L,ICFD),plaw(L,ICFD),'yes'
          else
            write(ITEM(M),'(a,2i7,2f9.3,6x,a)')KEY,L,abs(NCEL(L,ICFD)),
     &                                   REG(L,ICFD),plaw(L,ICFD),' no'
          endif
        endif
 10   CONTINUE

C If a long list include page facility text.      
      IF(IPFLG.EQ.0)THEN  
        ITEM(M+1) = ' ------------------------------------'
      ELSE
        WRITE(ITEM(M+1),15)IPM,MPM 
 15     FORMAT   ('0 page: ',I2,' of ',I2,' ---------')
      ENDIF
      ITEM(M+2) = '+ add/delete region'
      ITEM(M+3)  ='? help                              '
      ITEM(M+4)  ='- exit                              '

C Number of actual items displayed.
      MVERT=M+MCTL

C Help text for this menu.
      helptopic='cfd_cells_in_region'
      call gethelptext(helpinsub,helptopic,nbhelp)

      if (IAX.eq.1) then
        CALL EMENU('x-axis gridding',ITEM,MVERT,INO)
      elseif (IAX.eq.2) then
        CALL EMENU('y-axis gridding',ITEM,MVERT,INO)
      elseif (IAX.eq.3) then
        CALL EMENU('z-axis gridding',ITEM,MVERT,INO)
      elseif (IAX.eq.4) then
        CALL EMENU('ze-axis gridding',ITEM,MVERT,INO)
      endif

      if (INO.GT.MHEAD.AND.INO.LT.(MVERT-MCTL+1)) then

C Edit a regions data.
        CALL KEYIND(MVERT,INO,IFOC,IO)
        write (HOLD,'(i4,2f8.3)') abs(NCEL(IFOC,ICFD)),REG(IFOC,ICFD),
     &                            plaw(IFOC,ICFD)
        write(outs,'(2a)')'Please specify: number of cells; length of',
     &                    ' region; and power law coefficient.'

        call EASKS(HOLD,outs,' ',32,' 1 1.0 1.0 ','cell data',
     &    IER,nbhelp)
        K=0
        if (IAX.eq.1) then
          call EGETWI(HOLD,K,ic,1,NTCELX,'F','no of x cells',IIER)
C Here check the total number of cells along x-axis 
          numcel=0 
          do 20 ihelp=1,NRG(1,ICFD)
C Exclude current region from cell counting 
            if(ihelp.ne.IFOC) then 
              numcel=numcel+NCEL(ihelp,ICFD)
            endif
 20       continue
          if((ic+numcel).gt.(NTCELX-2)) then 
            iier=1 
            write(outs,'(a26,i4,a16,i4)')
     &        'The total No of cells now:'
     &        ,ic+numcel,' is bigger than ',NTCELX-2  
            CALL USRMSG(outs,
     &    ' which is the maximum allowed number for x-axis ','W')       
          endif 

        elseif (IAX.eq.2) then
          call EGETWI(HOLD,K,ic,1,NTCELY,'F','no of y cells',IIER)
C Here check the total number of cells along y-axis 
          numcel=0 
          do 30 ihelp=1,NRG(2,ICFD)
C Exclude current region from cell counting 
            if(ihelp.ne.IFOC) then 
              numcel=numcel+NCEL(ihelp,ICFD)
            endif
 30       continue
          if((ic+numcel).gt.(NTCELY-2)) then 
            iier=1 
            write(outs,'(a26,i4,a16,i4)')
     &        'The total No of cells now:'
     &        ,ic+numcel,' is bigger than ',NTCELY-2  
            CALL USRMSG(outs,
     &    ' which is the maximum allowed number for y-axis ','W')       
          endif 

        else
          call EGETWI(HOLD,K,ic,1,NTCELZ,'F','no of z cells',IIER)
C Here check the total number of cells along z-axis 
          numcel=0 
          do 40 ihelp=1,NRG(3,ICFD)
C Exclude current region from cell counting 
            if(ihelp.ne.IFOC) then 
              numcel=numcel+NCEL(ihelp,ICFD)
            endif
 40       continue
          if((ic+numcel).gt.(NTCELZ-2)) then 
            iier=1 
            write(outs,'(a26,i4,a16,i4)')
     &        'The total No of cells now:'
     &        ,ic+numcel,' is bigger than ',NTCELZ-2  
            CALL USRMSG(outs,
     &    ' which is the maximum allowed number for z-axis ','W')       
          endif 
        endif

C Only copy data if no errors.
        if (iier.eq.0) NCEL(IFOC,ICFD)=ic

C If more that one cell then ask if distribution is symmetrical.
        if (NCEL(IFOC,ICFD).gt.1) then
          call EGETWR(HOLD,K,ci,0.0,TOTDIS,'W','length of reg',IIER)
          if (iier.eq.0) REG(IFOC,ICFD)=ci
          call EGETWR(HOLD,K,ci,-99.0,99.0,'F','pwr law coef',IIER)
          if (iier.eq.0) plaw(IFOC,ICFD)=ci
          dok=.true.
          CALL ASKOK(' ',' Symmetrical gridding?',OK,dok,nbhelp)
          if (OK) NCEL(IFOC,ICFD) = -1*NCEL(IFOC,ICFD)
        else
          plaw(IFOC,ICFD)=1.0
        endif

      elseif (INO.eq.(MVERT-3)) then

C Change page.
      elseif (INO.eq.(MVERT-2)) then

C Add delete regions.
        call EASKABC('Region operations:','  ','add region',
     &    'delete region','cancel',IRO,nbhelp)
        if (IRO.eq.1) then
          if (NRG(IAX,ICFD)+1.gt.MNREG) then
            call usrmsg ('Maximum number of regions exceeded',
     &                   'returning to menu','W')
          else
            NRG(IAX,ICFD)=NRG(IAX,ICFD)+1
            IFOC=NRG(IAX,ICFD)

C Ask for length of current region (supply TOTDIS-CTDIS as default)
            CALL EASKABC('Specify region length:',' ','by vertices',
     &                   'by size','cancel',IW,nbhelp)
            if(iw.eq.1)then

C Locate regions by vertices.
              HOLD='  1   2'
              call EASKS(HOLD,'Specify vertex at start & end of region',
     &                   ' ',32,' 1 2 ','vertices',IER,nbhelp)
              K=0
              call EGETWI(HOLD,K,is,1,NTV,'F','vertex',IIER)
              call EGETWI(HOLD,K,ie,1,NTV,'F','vertex',IIER)
              if (IAX.eq.1) then
                REG(IFOC,ICFD) = abs(X(ie)-X(is))
              elseif (IAX.eq.2) then
                REG(IFOC,ICFD) = abs(Y(ie)-Y(is))
              else
                REG(IFOC,ICFD) = abs(Z(ie)-Z(is))
              endif
            elseif(iw.eq.2)then
              VAL=TOTDIS-CTDIS
              CALL EASKR(VAL,' ','Length of current region (m) ? ',
     &          0.0,'F',TOTDIS,'W',1.,'reg length',IER,nbhelp)
              REG(IFOC,ICFD) = VAL
            endif

C Got length of region now ask for no of cells and power law and symmetry.
C Offer the user a cancel option.
            if (IW.ne.3) then
              icell=NCEL(IFOC,ICFD)
              icellt=icell
              CALL EASKI(icellt,' Specify number of cells in region',
     &          ' ',1,'F',NTCELX,'F',10,'single region',IERI,nbhelp)
              if(ieri.eq.-3) then
                goto 5
              else
                icell=icellt
                NCEL(IFOC,ICFD) = icell
              endif

C If the number of cells is greater than 1 then ask for details, otherwise 
C set default size and power law.
              if(NCEL(IFOC,ICFD).gt.1) then
                VAL=1.0
                CALL EASKR(VAL,'  ','Power-law coefficient?',
     &              -4.0,'W',4.0,'W',1.0,'power coeff',IER,nbhelp)
                plaw(IFOC,ICFD)=VAL
                dok=.true.
                CALL ASKOK('  ','Symmetrical gridding?',OK,dok,nbhelp)
                if(OK) NCEL(IFOC,ICFD) = -1*NCEL(IFOC,ICFD)
              else
                plaw(IFOC,ICFD)=1.0
              endif
            endif
          endif
        elseif (IRO.eq.2) then
          IDEL=NRG(IAX,ICFD)
          CALL EASKI(IDEL,' ','Which region to delete ? ',
     &      0,'F',NRG(IAX,ICFD),'F',1,'del reg',IERI,nbhelp)
          if(ieri.eq.-3) goto 5

          if (IDEL.gt.0) then
            do 100 ID=IDEL,NRG(IAX,ICFD)
              if ((ID+1).le.MNREG) then
                REG(ID,ICFD)=REG(ID+1,ICFD)
                plaw(ID,ICFD)=plaw(ID+1,ICFD)
                NCEL(ID,ICFD)=NCEL(ID+1,ICFD)
              endif
 100        continue
            NRG(IAX,ICFD)=NRG(IAX,ICFD)-1
          endif
        endif
      elseif (INO.eq.(MVERT-1)) then

C Help.
        helptopic='cfd_cells_in_region'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL PHELPD('CFD inputs',nbhelp,'-',0,0,IER)
      elseif (INO.eq.(MVERT)) then
        call startbuffer()
        return
      endif
      goto 5
      end

C ************************ ASKMFNOD ************************
C This subroutine is copied here from esrures/mfoutp.F
C ASKMFNOD: Asks for one or more mass flow nodes.
C IOPT - restrict node types on display:
C  IOPT=0 - display all nodes
C  IOPT=1 - display internal nodes only.
C PROMPT1 and PROMPT2 are text prompt strings
C NH is number of help strings.

      SUBROUTINE ASKMFNOD(IOPT,NNPK,INPK,PROMPT1,PROMPT2,NH)
#include "net_flow.h"
#include "net_flow_data.h"

      DIMENSION INPK(MNOD)

      CHARACTER VERT(MNOD)*12
      CHARACTER PROMPT1*72, PROMPT2*72

C Set up array of menu strings.
      do 10 I=1,NNOD
        VERT(I)=NDNAM(I)
 10   continue
      call EPICKS(NNPK,INPK,PROMPT1,PROMPT2,
     &  12,NNOD,VERT,'Node select',IER,NH)

      return
      end

C ************************ ASKMFCON ************************
C This subroutine is copied here from esrures/mfoutp.F
C ASKMFCON: Asks for one or more mass flow connections.  If INODE 
C is greater than zero then user specified list of connections
C associated with this node. If INODE = zero then allow display
C of all connections in the network. If INODE is less than zero
C then silently return all connections associated with the node
C (absolute value of INODE).

      SUBROUTINE ASKMFCON(INODE,NNPK,ICPK,PROMPT1,PROMPT2)
#include "net_flow.h"
#include "net_flow_data.h"

      DIMENSION IPK(MCNN),ICPK(MCNN),ICID(MCNN)

      CHARACTER VRT(MCNN)*43
      CHARACTER PROMPT1*72, PROMPT2*72

C Set up array of menu strings.
      NITMS=0
      do 10 I=1,NCNN
        IP=NODPS(I); IN=NODNE(I); IC=ITPCON(I)
        if (abs(INODE).gt.0) then
          if (abs(INODE).eq.IP.or.abs(INODE).eq.IN) then

C A node has been specified -> only list its connections.
            NITMS=NITMS+1; ICID(NITMS)=I
            write(VRT(NITMS),'(5a)')
     &        NDNAM(IP),'->',NDNAM(IN),' via ',CMNAM(IC)
          endif
        else

C List all connections.
          NITMS=NITMS+1; ICID(NITMS)=I
          write(VRT(NITMS),'(5a)')
     &      NDNAM(IP),'->',NDNAM(IN),' via ',CMNAM(IC)
        endif
 10   continue

      if (INODE.ge.0) then
        if (NITMS.lt.NNPK) NNPK=NITMS
        call EPICKS(NNPK,IPK,PROMPT1,PROMPT2,
     &    43,NITMS,VRT,'Connection selection',IER,0)

C Decode menu picks to connections.
        do 20 I=1,NNPK
          ICPK(I)=ICID(IPK(I))
 20     continue
      else

C Act silently.
        NNPK=NITMS
        do 40 I=1,NNPK
          ICPK(I)=ICID(I)
 40     continue
      endif

C Debug.
C      do 30 I=1,NNPK
C        IP=NODPS(ICPK(I))
C        IN=NODNE(ICPK(I))
C        IC=ITPCON(ICPK(I))
C        write(6,'(5a)') NDNAM(IP),'->',NDNAM(IN),' via ',CMNAM(IC)
C 30   continue

      return
      end
