C This file is part of the ESP-r system.
C Copyright Energy Systems Research Unit, University of
C Strathclyde, Glasgow Scotland, 2001-.

C ESP-r is free software.  You can redistribute it and/or
C modify it under the terms of the GNU General Public
C License as published by the Free Software Foundation 
C (version 2 orlater).

C ESP-r is distributed in the hope that it will be useful
C but WITHOUT ANY WARRANTY; without even the implied
C warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
C PURPOSE. See the GNU General Public License for more
C details.


C cmmon3dv.F contains the following Fortran subroutines:
C  CADJVIEW: Adjusts the viewing transforms & boundaries if necessary
C            using global data structures.
C  CDRAWESP: Loads the geometric description(s), converts into viewing 
c            coords, clips the image and draws it.
C  CHGAZI:  Called from C code with an increment or decrement of azimuth.
C  CHGELEV: Called from C code with an increment or decrement of elevation.
C  CHGEYE:  Called from C code with current contents of the ray2_ and
C           image_ c structures.
C  DRAWOBS: Checks for obstructions, loads, converts and draws them.
C  DRAWVIS: draws any visual entities within zone as slightly grey lines.

C  DRWSEN:  Loads, converts and draws MRT sensors (rec bodies).
C  EGRNDR:  Draws ground topology.
C  DSGRID:  Draws a dotted site grid with current scaling factors.
C  EMKVIEW: Constructs a 'viewer' format file.
C  wiresymbol: draws a user defined symbol in wireframe
C  REDRAW : Entry point for redrawing model display.
C  TMPMENU: A temporary menu that allows users to reposition the view.

C ************* CADJVIEW 
C Takes the current status of view and boundary mods and adjusts
C the viewing transforms & object boundaries if necessary before
C displaying the zone.
C IER=0 OK, IER=1 problem. 
C izgfoc > zero represents the zone which is currently being edited.
C If focussname is true than only labels highlighted surfaces.

      SUBROUTINE CADJVIEW(focussname,IER)
#include "building.h"
#include "model.h"
#include "prj3dv.h"
      
      integer lnblnk  ! function definition

C Passed parameter.
      logical focussname ! if true only label highlighted surfaces.
      integer ier

      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      COMMON/SPAD/MMOD,LIMIT,LIMTTY
      COMMON/gzonpik/izgfoc,nzg,nznog(mcom)
      common/initv/initvt,EYEMI(3),VIEWMI(3),ANGI
      INTEGER :: initvt
      integer menuchw,igl,igr,igt,igb,igw,igwh
      COMMON/VIEWPX/menuchw,igl,igr,igt,igb,igw,igwh

C Common HOLDVW remembers linescale parameters, so that we can display
C something else (say a graph) and then return to static view bounds in
C the same position as we left it.
      COMMON/HOLDVW/iglhld,igbhld,Xaddhld,Yaddhld,scahld

      character ETEXT*82
#ifdef OSI
      integer iside,isize,ifont     ! passed to viewtext
      integer numberofzones
      integer JITZNM,JITSNM,JITVNO,JITOBS,JITVIS,JITVOBJ,JITSNR,JITGRD
      integer JITORG,JITBND,JITDSP,JITHLS,JITHLZ,JITPPSW
      integer jizgfoc,jnzg,jnznog,jnznogv
      integer igwid,igheight          ! for use with axiscale
      integer iigl,iigr,iigt,iigb,iigw,iigwh
      integer iiw1,iiw2,iiw3,iiw4,iimenu
#else
      integer*8 iside,isize,ifont     ! passed to viewtext
      integer*8 numberofzones
      integer*8 JITZNM,JITSNM,JITVNO,JITOBS,JITVIS,JITVOBJ,JITSNR
      integer*8 JITGRD,JITORG,JITBND,JITDSP,JITHLS,JITHLZ,JITPPSW
      integer*8 jizgfoc,jnzg,jnznog,jnznogv
      integer*8 igwid,igheight          ! for use with axiscale
      integer*8 iigl,iigr,iigt,iigb,iigw,iigwh
      integer*8 iiw1,iiw2,iiw3,iiw4,iimenu
#endif

      real EVX,EVY,EVZ,VX,VY,VZ,EAN,GDIS

      SMALL=0.01

C Tell wireframe control the current number of selected zones.
C If only drawing a ground topology set updwire manually.
      if(itdsp.eq.6)then
        numberofzones=1
        if(mmod.eq.8)call updwire(numberofzones)  ! pass local integer
      else
        numberofzones=nzg
        if(mmod.eq.8)call updwire(numberofzones)  ! pass local integer
      endif

C Also tell wireframe control about the current views. If the
C user has previously set initvt to 1 by editing the viewing
C parameters then use EYEMI & VIEWMI rather than EYEM VIEWM.
      EVX=EYEM(1); EVY=EYEM(2); EVZ=EYEM(3)
      if(initvt.eq.0)then
        VX=VIEWM(1); VY=VIEWM(2); VZ=VIEWM(3)
      else
        VX=VIEWMI(1); VY=VIEWMI(2); VZ=VIEWMI(3)
      endif
      EAN=ANG; GDIS=GRDIS
      JITZNM=ITZNM; JITSNM=ITSNM; JITVNO=ITVNO
      JITOBS=ITOBS; JITVIS=ITVIS; JITVOBJ=ITVOBJ
      JITSNR=ITSNR; JITGRD=ITGRD
      JITORG=ITORG; JITBND=ITBND; JITDSP=ITDSP
      JITHLS=ITHLS; JITHLZ=ITHLZ; JITPPSW=ITPPSW

      call curviews(EVX,EVY,EVZ,VX,VY,VZ,EAN,JITZNM,JITSNM,JITVNO,
     &  JITOBS,JITVIS,JITVOBJ,JITSNR,JITGRD,JITORG,GDIS,JITBND,
     &  JITDSP,JITHLS,JITHLZ,JITPPSW)

C Tell wireframe how many zones to include in the image.
      jizgfoc=izgfoc; jnzg=nzg
      call pushgzonpik(jizgfoc,jnzg)
      if(nzg.gt.0)then
        do 42 ij=1,nzg
          jnznog=ij; jnznogv=nznog(ij)
          call pushnznog(jnznog,jnznogv)
  42    continue
      endif

C If altered then do a bound check on the current zone and update display.
C If static views ITBND=0 then do not update VIEWM when calling BNDOBJ.
      IF(MODIFYVIEW)THEN
        if(izgfoc.NE.0) CALL ESCZONE(izgfoc)
        IF(MODBND)THEN
          CALL BNDOBJ(0,IER)
          MODLEN=.TRUE.
        ENDIF

C If viewpoint or bounds changed then recalculate viewing parameters.
        IF(MODLEN)THEN
          DIS = (VIEWM(1)-EYEM(1))**2 + (VIEWM(2)-EYEM(2))**2 +
     &          (VIEWM(3)-EYEM(3))**2
          IF(DIS.GE.SMALL)THEN
            HITH=1.0; YON=1300.0
          ELSE
            ier=1
            CALL USRMSG(' ',
     &       ' The eye and viewed position are too close!','W')
            RETURN
          ENDIF
          CALL LENS(IER)
        ENDIF

C Clear current viewing box.
        if(mmod.eq.8)CALL startbuffer()

C Setup and pass in parameters to win3d.
C Add small margins on all sides of the graphic window.
        iiw1=4; iiw2=4; iiw3=2; iiw4=1; iimenu=menuchw
        iigl=igl; iigr=igr; iigt=igt; iigb=igb; iigw=igw; iigwh=igwh
        if(mmod.eq.8)then
          CALL win3d(iimenu,iiw1,iiw2,iiw3,iiw4,
     &      iigl,iigr,iigt,iigb,iigw,iigwh)
        else
          CALL win3dwwc(iimenu,iiw1,iiw2,iiw3,iiw4,
     &      iigl,iigr,iigt,iigb,iigw,iigwh)
        endif
        igl=int(iigl); igr=int(iigr); igt=int(iigt); igb=int(iigb)
        igw=int(iigw); igwh=int(iigwh)
        igwid=igw; igheight=igwh

C If optimal view bounds, reset scaling ratios.
        IF(ITBND.EQ.1)THEN
          CALL SITE2D(SXMX,SXMN,SYMX,SYMN,ier)
          CALL axiscale(igwid,igheight,SXMN,SXMX,SYMN,SYMX,xsc,ysc,sca,
     &      Xadd,Yadd)
          call linescale(iigl,Xadd,sca,iigb,Yadd,sca)
          iglhld=igl; Xaddhld=Xadd; igbhld=igb; Yaddhld=Yadd; scahld=sca ! remember values

C If static bounds, reset linescale in case we have returned from
C something else that changed it.
        ELSE
          iigl=igl; iigb=igb
          call linescale(iigl,Xaddhld,scahld,iigb,Yaddhld,scahld)
        ENDIF

C If view or scaling has changed, reset clipping planes.
        if (MODLEN .or. ITBND.EQ.1) call INICLP(ier)

C Draw the selected bodies based on the passed focussname.
        CALL CDRAWESP(focussname,IER)
        if (ier.ne.0) return

C Show feature buttons.
        if(mmod.eq.8) call redrawbuttons()

C Remind user of which model even when looking at one zone.
        WRITE(ETEXT,'(2A)')'Model: ',modeltitle(1:lnblnk(modeltitle))
        iside=1; isize=1; ifont=1
        if(mmod.eq.8)then
          call viewtext(etext,iside,isize,ifont)
        else
C          call viewtextwwc(etext,iside,isize,ifont)
        endif

C Recover information begin edited and reset modifyview flag.
        if(izgfoc.NE.0) CALL ERCZONE(izgfoc)
        MODIFYVIEW=.FALSE.
        MODLEN=.FALSE.
        MODBND=.FALSE.
      ENDIF

      RETURN
      END


C ******************** CDRAWESP 
C Loads the geometric description(s) from common blocks,
C converts into viewing coords, clips the image and draws it.
C ier is returned as non-zero if there has been a problem.

      SUBROUTINE CDRAWESP(focussname,ier)
#include "building.h"
#include "model.h"
#include "geometry.h"
#include "prj3dv.h"
#include "esprdbfile.h"
#include "material.h"
      
      integer lnblnk  ! function definition

C Passed parameters.
      logical focussname ! if true only label highlighted surfaces.
      integer ier

      COMMON/gzonpik/izgfoc,nzg,nznog(mcom)
      COMMON/SPAD/MMOD,LIMIT,LIMTTY
      COMMON/OUTIN/IUOUT,IUIN,IEOUT
C      COMMON/FILEP/IFIL
      integer ifs,itfs,imfs
      COMMON/GFONT/IFS,ITFS,IMFS
      integer menuchw,igl,igr,igt,igb,igw,igwh
      COMMON/VIEWPX/menuchw,igl,igr,igt,igb,igw,igwh

      integer IZSTOCN
      COMMON/C24/IZSTOCN(MCOM,MS)
      common/appcols/mdispl,nifgrey,ncset,ngset,nzonec

      COMMON/FOPENED/CFGOK,MLDBOK,MATDBOK,CTLOK,OPTKOK,CFCDBOK
      LOGICAL  ISTOK,CFGOK,MLDBOK,MATDBOK,CTLOK,OPTKOK,CFCDBOK
      logical greyok,found
      integer ier2,ier3,ier4,ier6  ! traps

      DIMENSION  AX(MPL),AY(MPL),AZ(MPL),IANXT(MPL)
      DIMENSION  BX(MPL),BY(MPL),BZ(MPL),IBNXT(MPL)
      DIMENSION  SBBOX(3,2),COG1(3),COG2(3)
      DIMENSION  CG(3),VN(3),ipoints(6,2)

      CHARACTER ETEXT*82
      CHARACTER temp*16,DESCRC*25,stemp*13,temp20*21
      CHARACTER MLCN*32,t13*13
      integer izsize,issize

C Local variables to pass to edline and axiscale.
#ifdef OSI
      integer iside,isize,ifont       ! passed to viewtext
      integer iix,iiy,iixlast,iiylast,iixc,iiyc,iid1,iid2,iid3,iid4
      integer iicol,ipixw,ipixh       ! passed to textsizeatxy
#else
      integer*8 iside,isize,ifont     ! passed to viewtext
      integer*8 iix,iiy,iixlast,iiylast,iixc,iiyc,iid1,iid2,iid3,iid4
      integer*8 iicol,ipixw,ipixh     ! passed to textsizeatxy
#endif

C If libraries and monitor support greyscale remember this.
      greyok=.false.
      if(nifgrey.gt.4)then
        greyok=.true.
      endif

C Set font for displaying zone & surface names etc.
      issize=4; izsize=5  ! for surfaces and zone name

      IER=0; ier2=0; ier3=0; ier4=0; ier6=0

C Refresh the dialogue box if terminal type 8.
      IF(MMOD.EQ.8)CALL USRMSG(' ',' ','-')

C If site grid or ground is asked for draw now.
      if(ITGRD.EQ.0)then
        CALL DSGRID(0.0,GRDIS,5,IER2)
        if(ier2.ne.0)then
          call to_session('CDRAWESP: Error in drawing grid.')
          ier=-2
        endif
      endif
      if(ITDSP.EQ.5.or.ITDSP.EQ.6)then
        call EGRNDR(IER3)
        if(ier3.ne.0)then
          call to_session('CDRAWESP: Error in drawing ground.')
          ier=-3
        endif
      endif
      IF(ITDSP.EQ.6)return

C Remind user of which project.
      WRITE(etext,'(2A)')'Model: ',modeltitle(1:lnblnk(modeltitle))
      iside=1; isize=1; ifont=1
      if(mmod.eq.8)then
        call viewtext(etext,iside,isize,ifont)
      else
C        call viewtextwwc(etext,iside,isize,ifont)
      endif

C The first task is to set iz to zero and istok to true to force the
C code to draw the site origin symbol.
      IZ=0
      ISTOK=.TRUE.

  100 CONTINUE

C If the origin of the site is to be shown include this in the image
C and then go on to the normal objects. After this point ISTOK is .FALSE.
      if(ITORG.EQ.0.and.ISTOK)then
        if(mmod.eq.8)CALL ORIGESP
        ISTOK=.FALSE.
      elseif(ITORG.EQ.1.and.ISTOK)then
        ISTOK=.FALSE.
      endif

      IZ=IZ+1
      IF(IZ.GT.nzg)RETURN

C If obstructions are to be drawn transform attributes into a
C gen body plotted before going on to the surfaces in the zone itself. 
C Also draw any visual entities included in the model (e.g. to preview 
C what will be passed to Radiance).
      if(ITDSP.EQ.0.or.ITDSP.eq.5.or.ITDSP.eq.7)then
        CALL DRAWOBS(nznog(IZ),ier4)
        if(ier4.ne.0)then
          call to_session('CDRAWESP: Error drawing obstructions.')
          ier=-4
        endif
      endif
      if(ITDSP.EQ.0.or.ITDSP.eq.8)then
        CALL DRAWVIS(nznog(IZ),ier4)
        if(ier4.ne.0)then
          call to_session('CDRAWESP: Error drawing visual entities.')
          ier=-4
        endif
      endif


C If zone names are to be displayed place near the centre of zone.
C For clarity use larger font and then restore font to normal size.
C Transform the site COG coords into eyepoint and then screen coords.
C Find the pixels required for the name and centre text.
C Note: temp character t13 is 1 char wider than zname().
      IF(ITZNM.EQ.0)THEN
        COG1(1)=ZCOG(nznog(IZ),1)
        COG1(2)=ZCOG(nznog(IZ),2)
        COG1(3)=ZCOG(nznog(IZ),3)
        CALL VECTRN(COG1,TSMAT,COG2,IER6)
        if(ier6.ne.0)then
          call to_session('CDRAWESP: Error in name coord transforms.')
          ier=-6
        endif
        call u2pixel(COG2(1),COG2(2),iix,iiy)

C Note call to textpixwidth corrupts t13 so need to
C re-create t13 after textpixwidth call.
        t13=' '; iicol=0
        lnzn=lnblnk(zname(nznog(IZ)))
        write(t13,'(A)') zname(nznog(IZ))(1:lnzn)
        if(mmod.eq.8)then
          call winfnt(izsize)
          call textpixwidth(t13,ipixw,ipixh)
        else
          ipixw=lnblnk(t13)*7  ! assume 7 pixels wide
          ipixh=7              ! assume 7 pixels high
        endif

        write(t13,'(A)')zname(nznog(IZ))(1:lnzn)
        iid4 = iix - (ipixw/2)

C Check for clipping.
        call pixel2u(iid4,iiy,xx,xy)
        call CLIPPT(xx,xy,COG2(3),iclp1)
        call pixel2u(iix+(ipixw/2),iiy-ipixh,xx,xy)
        call CLIPPT(xx,xy,COG2(3),iclp2)
        if (iclp1.eq.0 .and. iclp2.eq.0) then
          iicol=0
          if(mmod.eq.8)then
            call textatxy(iid4,iiy,t13,'-',iicol)
          else
            call textatxywwc(iid4,iiy,t13,'-',iicol)
          endif
        endif
        if (mmod.eq.8) CALL winfnt(IMFS)
      ENDIF

      if(mmod.eq.8) call forceflush()

C If there are no surfaces yet defined then draw the vertices with 
C label and loop back to the next body.
      IF(NZSUR(nznog(IZ)).EQ.0.AND.NZTV(nznog(IZ)).GT.1)THEN
        DO 202 IV=1,NZTV(nznog(IZ))
          COG1(1)=szcoords(nznog(IZ),IV,1)
          COG1(2)=szcoords(nznog(IZ),IV,2)
          COG1(3)=szcoords(nznog(IZ),IV,3)
          CALL VECTRN(COG1,TSMAT,COG2,IER)
          call CLIPPT(COG2(1),COG2(2),COG2(3),iclp) ! check clipping
          if (iclp.eq.0) then
            call u2pixel(COG2(1),COG2(2),iix,iiy)
            if(mmod.eq.8) call ecirc(iix,iiy,3,1)
            if(mmod.eq.8) CALL VERTLBL(iix,iiy,COG2(3),IV,ier)
          endif
  202   CONTINUE
        if(mmod.eq.8) call forceflush()
        GOTO 100
      ELSEIF(NZSUR(nznog(IZ)).GE.1.AND.NZTV(nznog(IZ)).GT.1)THEN

C Check for unlinked vertices, if so then draw them before proceeding
C to draw the surfaces. For each vertex, loop through each surface
C and its vertex list and see how many hits.
        found=.false.
        DO 204 IZV=1,NZTV(nznog(IZ))
          IHIT=0
          DO 8792 IS=1,NZSUR(nznog(IZ))
            DO 8794 IV=1,isznver(nznog(IZ),IS)
              IF(IZV.EQ.iszjvn(nznog(IZ),IS,IV))IHIT=IHIT+1
 8794       CONTINUE
 8792     CONTINUE

C If less than 2 hits draw the vertex as an open circle, if no hits
C draw it as a solid circle.
          IF(IHIT.LT.2)THEN
            COG1(1)=szcoords(nznog(IZ),IZV,1)
            COG1(2)=szcoords(nznog(IZ),IZV,2)
            COG1(3)=szcoords(nznog(IZ),IZV,3)
            CALL VECTRN(COG1,TSMAT,COG2,IER)
            call CLIPPT(COG2(1),COG2(2),COG2(3),iclp) ! check clipping
            if (iclp.eq.0) then
              found=.true.
              call u2pixel(COG2(1),COG2(2),iix,iiy)
              IF(IHIT.EQ.0)THEN
                if(mmod.eq.8) CALL ecirc(iix,iiy,3,1)
              ELSEIF(IHIT.GT.0.AND.IHIT.LT.2)THEN
                if(mmod.eq.8) CALL ecirc(iix,iiy,3,0)
              ENDIF
              if(mmod.eq.8) CALL VERTLBL(iix,iiy,COG2(3),IZV,ier)
            endif
          ENDIF
  204   CONTINUE

C If any were found then draw a legend for the meaning of the vertex
C symbols if in graphic mode.
        if(found .and. mmod.eq.8)then
          CALL winfnt(issize)
          iid1=igl+10; iid2=igb+5;
          CALL ecirc(iid1,iid2,3,1)
          write(temp,'(A)')'Unlinked vertex'
          iid1=igl+20; iid2=igb+10;
          iicol=0
          call textatxy(iid1,iid2,temp,'-',iicol)
          iid1=igl+130; iid2=igb+5;
          if(mmod.eq.8) CALL ecirc(iid1,iid2,3,0)
          write(temp20,'(A)')'Single-linked vertex'
          iid1=igl+140; iid2=igb+10;
          call textatxy(iid1,iid2,temp20,'-',iicol)
          CALL winfnt(IMFS)
        endif
        if(mmod.eq.8) call forceflush()
      ENDIF

C Continue on with any surfaces in the zone. If zone is zero then
C we are drawing grid.
      DO 1000 I=1,NZSUR(nznog(IZ))
        isur=I
        if(nznog(IZ).ne.0)then
          icc=IZSTOCN(nznog(IZ),isur)
        else
          icc=0
        endif

C If the configuration file has been read and there is a filter for
C which surfaces should be displayed test this here. Locate the
C connection (IC) which relates to this surface if a real zone.
        IF(CFGOK)THEN
          if(iz.ne.0)then
            CALL SURADJ(nznog(IZ),isur,IE,TMP,IZC,ISC,IC,DESCRC)
            IF(ITDSP.EQ.2.AND.IE.NE.0)GOTO 1000
            IF(ITDSP.EQ.3.AND.IE.NE.3)GOTO 1000
            IF(ITDSP.EQ.4.AND.(IE.EQ.0.OR.IE.EQ.3))GOTO 1000
          endif
        ENDIF

C Copy polygon for Surface I into structure A for compatibility with
C viewer format 'holes' in surfaces.
C NAP   = Number of vertex points on surface
C IANXT = Index of 'next' point
C IAPNT = Pointer to first vertex of polygon
        IAPNT = 1
        NAP = isznver(nznog(IZ),I)
        DO 1100 J = 1,isznver(nznog(IZ),I)
          K = iszjvn(nznog(IZ),I,J)
          if(k.gt.0)then
            AX(J) = szcoords(nznog(IZ),K,1)
            AY(J) = szcoords(nznog(IZ),K,2)
            AZ(J) = szcoords(nznog(IZ),K,3)
            IANXT(J) = J + 1
          else
            AX(J) = 0.0; AY(J) = 0.0; AZ(J) = 0.0 ! if JVN malformed
          endif
 1100   CONTINUE

        IANXT(isznver(nznog(IZ),I)) = IAPNT

C Transform surface polygon to screen co-ordinates
C Take structure A multiply by TSMAT return structure B
C TSMAT = Model to Screen Matrix
        CALL MATPOL(NAP,IAPNT,AX,AY,AZ,IANXT,TSMAT,
     &              SBBOX,NBP,IBPNT,BX,BY,BZ,IBNXT,IERR)

C Generate clipping flags and clip geometry in eye coords.
C If ISTAT =  0 : totally inside frustrum
C If ISTAT =  1 : totally outside frustrum
C If ISTAT = -1 : straddles frustrum
        CALL CLIPFL(NBP,BX,BY,BZ,ISTAT)
        IF (ISTAT .EQ. 1) THEN
          GOTO  1000
        ELSEIF (ISTAT .EQ. -1) THEN
          CALL CUTPOL(NB,NBP,IBPNT,BX,BY,BZ,IBNXT,ISTAT)
        else
          NB=1
        ENDIF

        iixc=0
        iiyc=0

C Set line style for this surface. linsty=2 is double width line. 
        if(ICC.ne.0)then
          if(LINSTY(ICC).eq.2)then
            IPEN= -305
          else
            IPEN= 1
          endif
        else
          IPEN= 1
        endif

C Check for hilights. ITHLS=1 is composition, ITHLS=2 trn:opq,
C ITHLS = 3 partial attributes.
        if(ITHLS.eq.1.and.ITHLZ.gt.0)then
          WRITE(MLCN,'(A)') mlcname(ITHLZ)
          if(ICC.gt.0)then
            lnssmlc=lnblnk(SMLCN(nznog(IZ),isur))
            if(MLCN(1:lnssmlc).eq.SMLCN(nznog(IZ),isur)(1:lnssmlc))then
              IPEN= -305
            else
              IPEN= 1
            endif
          else
            IPEN= 1
          endif
        elseif(ITHLS.eq.2)then
          if(ICC.gt.0)then
            if(ITHLZ.eq.1.and.SOTF(nznog(IZ),isur)(1:4).eq.'OPAQ')then
              IPEN= -305
            else
              IPEN= 1
            endif
            if(ITHLZ.eq.2.and.SOTF(nznog(IZ),isur)(1:4).eq.'TRAN')then
              IPEN= -305
            else
              IPEN= 1
            endif
          else
            IPEN= 1
          endif
          if(ITHLZ.eq.2.and.SOTF(nznog(IZ),I)(1:4).ne.'OPAQ')then
            IPEN= -305
          else
            IPEN= 1
          endif
        elseif(ITHLS.eq.3)then
          if(ICC.gt.0)then
            if(SMLCN(nznog(IZ),isur)(1:3).eq.'UNK'.or.
     &         zboundarytype(nznog(IZ),isur,1).eq.-1)then
              IPEN= -305
            else
              IPEN= 1
            endif
          else
            IPEN=1
          endif
        endif

C Debug.
C        write(6,*)'icc linsty ipen ithls ithlz',icc,linsty(icc),ipen,
C     &    ithls,ithlz,ISTOK

C Count the edges in case enhanced edge drawing required.
        iedge=0
        DO 300 J = 1,NB

C Move to the beginning and transform screen to pixel coords.
          IP = IABS(IBPNT)
          call u2pixel(BX(IP),BY(IP),iix,iiy)

C If vertex are to be labeled do it in the smallest font, also place
C a 2pixel square dot at the vertex so if it is embedded in a line
C it can be seen. Then move cursor to vertex and remember position.
          IF(ITVNO.EQ.0.AND.NAP.EQ.NBP)THEN
            IV=iszjvn(nznog(IZ),I,IP)
            if(mmod.eq.8)then
              call esymbol(iix,iiy,8,1)
            else
              call esymbolwwc(iix,iiy,8,1)
            endif
            CALL VERTLBL(iix,iiy,BZ(IP),IV,ier)
          ENDIF
          iixlast=iix   ! cast for edline call
          iiylast=iiy

          IP1 = IBNXT(IP)

 351      CONTINUE
          iedge=iedge+1

C Ensure edge is less than MV.
          if(iedge.gt.MV) goto 300

          call u2pixel(BX(IP1),BY(IP1),iix,iiy)
          if(IPEN.eq.1)then

C Normal (not highlighted) edge being drawn, if NAP = NBP then check
C common/G8 for hints as to alternate greys or dots.

C Debug.
C            write(6,*) 'nap nbp nbedgdup iedge j',nap,nbp,nbedgdup(icc),
C     &        iedge,j

            if(NAP.ne.NBP)then
              if(mmod.eq.8)then
                call eswline(iixlast,iiylast,iix,iiy)
              else
                call eswlinewwc(iixlast,iiylast,iix,iiy)
              endif
            elseif(icc.eq.0)then
              if(mmod.eq.8)then
                call eswline(iixlast,iiylast,iix,iiy)
              else
                call eswlinewwc(iixlast,iiylast,iix,iiy)
              endif
            else

C If there will be diagonal lines drawn (duplicate edges) check if the
C current edge is marked. If dotted line is drawn then don't bother
C with the co-planer.
              found=.false.
              if(nbedgdup(icc).ne.0)then
                if(iedgdup(icc,iedge).eq.0)then
                  found=.false.
                else
                  LD=4
                  if(mmod.eq.8)then
                    call edline(iixlast,iiylast,iix,iiy,LD)
                  else
                    call edlinewwc(iixlast,iiylast,iix,iiy,LD)
                  endif
                  found=.true.
                endif
              endif

C If this edge is shared with a co-planer surface with the same construction
C and greys can be done reset forground to 50% grey otherwise draw black,
C unless duplicate line has already been drawn.
              if(imatshr(icc,iedge).eq.0)then
                if(.NOT.found)then
                  if(mmod.eq.8)then
                    call eswline(iixlast,iiylast,iix,iiy)
                  else
                    call eswlinewwc(iixlast,iiylast,iix,iiy)
                  endif
                endif
              else
                if(greyok)then
                  iicol=nifgrey-3
                  if(mmod.eq.8) call winscl('i',iicol)
                endif
                if(mmod.eq.8)then
                  call eswline(iixlast,iiylast,iix,iiy)
                else
                  call eswlinewwc(iixlast,iiylast,iix,iiy)
                endif
                if(greyok)then
                  iicol=0
                  if(mmod.eq.8) call winscl('-',iicol)
                endif
              endif
            endif
          elseif(IPEN.eq.-305)then
            if(mmod.eq.8)then
              call edwline(iixlast,iiylast,iix,iiy)
            else
              call edwlinewwc(iixlast,iiylast,iix,iiy)
            endif
          else
            if(mmod.eq.8)then
              call eswline(iixlast,iiylast,iix,iiy)
            else
              call eswlinewwc(iixlast,iiylast,iix,iiy)
            endif
          endif
          iixc=iixc+iix
          iiyc=iiyc+iiy
          if(mmod.eq.8) call forceflush()

C If vertex are to be labeled do it in the smallest font and
C move cursor back to vertex and remember position.
          IF(ITVNO.EQ.0.AND.NAP.EQ.NBP)THEN
            IV=iszjvn(nznog(IZ),I,IP1)
            if(mmod.eq.8)then
              call esymbol(iix,iiy,8,1)
            else
              call esymbolwwc(iix,iiy,8,1)
            endif
            CALL VERTLBL(iix,iiy,BZ(IP1),IV,ier)
            iixlast=iix   ! cast for edline call
            iiylast=iiy
          ENDIF

C Repeat until next point is the first, remember position.
          IP1 = IBNXT(IP1)
          iixlast=iix   ! cast for edline call
          iiylast=iiy
          IF(IP1 .NE. IP) GOTO 351

C Complete the polygon.
          iedge=iedge+1
          call u2pixel(BX(IP),BY(IP),iix,iiy)
          if(IPEN.eq.1)then

C Normal (not highlighted) edge being drawn, if NAP = NBP then check
C common/G8 for hints as to alternate greys or dots.
            if(NAP.ne.NBP)then
              if(mmod.eq.8)then
                call eswline(iixlast,iiylast,iix,iiy)
              else
                call eswlinewwc(iixlast,iiylast,iix,iiy)
              endif
            elseif(icc.eq.0)then
              if(mmod.eq.8)then
                call eswline(iixlast,iiylast,iix,iiy)
              else
                call eswlinewwc(iixlast,iiylast,iix,iiy)
              endif
            else

C If there will be diagonal lines drawn (duplicate edges) check if the
C current edge is marked. If dotted line is drawn then don't bother
C with the co-planer.
              found=.false.
              if(nbedgdup(icc).ne.0)then
                if(iedgdup(icc,iedge).eq.0)then
                  found=.false.
                else
                  LD=4
                  if(mmod.eq.8)then
                    call edline(iixlast,iiylast,iix,iiy,LD)
                  else
                    call edlinewwc(iixlast,iiylast,iix,iiy,LD)
                  endif
                  found=.true.
                endif
              endif

C If this edge is shared with a co-planer surface with the same construction
C and greys can be done reset forground to 50% grey otherwise draw black,
C unless duplicate line has already been drawn.
              if(imatshr(icc,iedge).eq.0)then
                if(.NOT.found)then
                  if(mmod.eq.8)then
                    call eswline(iixlast,iiylast,iix,iiy)
                  else
                    call eswlinewwc(iixlast,iiylast,iix,iiy)
                  endif
                endif
              else
                if(greyok)then
                  iicol=nifgrey-3
                  if(mmod.eq.8) call winscl('i',iicol)
                endif
                if(mmod.eq.8)then
                  call eswline(iixlast,iiylast,iix,iiy)
                else
                  call eswlinewwc(iixlast,iiylast,iix,iiy)
                endif
                if(greyok)then
                  iicol=0
                  if(mmod.eq.8) call winscl('-',iicol)
                endif
              endif
            endif
          elseif(IPEN.eq.-305)then
            if(mmod.eq.8)then
              call edwline(iixlast,iiylast,iix,iiy)
            else
              call edwlinewwc(iixlast,iiylast,iix,iiy)
            endif
          else
            if(mmod.eq.8)then
              call eswline(iixlast,iiylast,iix,iiy)
            else
              call eswlinewwc(iixlast,iiylast,iix,iiy)
            endif
          endif
          iixc=iixc+iix
          iiyc=iiyc+iiy
          if(mmod.eq.8) call forceflush()
  300   CONTINUE

C If surface names are to be displayed place near the centre of polygon.
C If the surface is not drawn (iixc=0 & iiyc=0) do not label it. If the
C surface has been clipped then do not label it. For clarity drop to
C the next smaller font for the surface name and then restore font
C to normal size.
C Patch for cfg is to only label highlighted surfaces.
        IF(ITSNM.EQ.0)THEN
C          write(6,*) 'cdrawesp_iixc iiyc nap nbp ',iixc,iiyc,nap,nbp
          IF(iixc.GT.0.AND.iiyc.GT.0.AND.(NAP.EQ.NBP))THEN
            if(icc.ne.0)then
              if(mmod.eq.8)CALL winfnt(issize)
              write(stemp,'(A)')SNAME(nznog(IZ),isur)
              call textpixwidth(stemp,ipixw,ipixh)    ! corrupts stemp
              write(stemp,'(A)')SNAME(nznog(IZ),isur)

C Logic that follows is to only label surfaces that are highlighted. 
              if(focussname)then
                if(LINSTY(icc).eq.2)then
                  CG(1)=SURCOG(nznog(IZ),I,1)
                  CG(2)=SURCOG(nznog(IZ),I,2)
                  CG(3)=SURCOG(nznog(IZ),I,3)
                  VN(1)=SURVN(nznog(IZ),I,1)
                  VN(2)=SURVN(nznog(IZ),I,2)
                  VN(3)=SURVN(nznog(IZ),I,3)
                  CALL VECTRN(VN,TSMAT,COG2,IER)
                  call arrow(CG,VN,0.3,0.1,ipoints,'a',2)
                  call u2pixel(COG2(1),COG2(2),iixc,iiyc)
                  write(stemp,'(A)')SNAME(nznog(IZ),I)
                  CALL winfnt(issize)
                  iixc=iixc+5

C Check clipping.
                  call CLIPST(stemp,iixc,iiyc,COG2(3),iclp)
                  if (iclp.eq.0) then
                    iicol=0
                    call textatxy(iixc,iiyc,stemp,'-',iicol)
                  endif
                  CALL winfnt(IMFS)
                endif
              else

C Normal display of surface names as in cdrawesp.
C If surface normals to be displayed draw an arrow and
C offset the surface name depending on orientation.
                if(itsnr.eq.0)then
                  CG(1)=SURCOG(nznog(IZ),isur,1)
                  CG(2)=SURCOG(nznog(IZ),isur,2)
                  CG(3)=SURCOG(nznog(IZ),isur,3)
                  VN(1)=SURVN(nznog(IZ),isur,1)
                  VN(2)=SURVN(nznog(IZ),isur,2)
                  VN(3)=SURVN(nznog(IZ),isur,3)
                  CALL VECTRN(VN,TSMAT,COG2,IER)
                  call u2pixel(COG2(1),COG2(2),iixc,iiyc)

C Check for clipping.
                  call CLIPPT(COG2(1),COG2(2),COG2(3),iclp1)
                  call pixel2u(iixc+ipixw+5,iiyc-ipixh,xx,xy)
                  call CLIPPT(xx,xy,COG2(3),iclp2)

                  if (iclp1.eq.0 .and. iclp2.eq.0) then
                    call arrow(CG,VN,0.3,0.1,ipoints,'a',2)
                    IF(SVFC(nznog(IZ),isur).EQ.'VERT'.or.
     &                SVFC(nznog(IZ),isur).EQ.'CEIL'.or.
     &                SVFC(nznog(IZ),isur).EQ.'FLOR')THEN
                      iixc=iixc+5
                    ELSE
                      iixc=iixc+5
                    ENDIF
                    iicol=0
                    if(mmod.eq.8)then
                      call textatxy(iixc,iiyc,stemp,'-',iicol)
                    else
                      call textatxywwc(iixc,iiyc,stemp,'-',iicol)
                    endif
                  endif
                else
                  CG(1)=SURCOG(nznog(IZ),isur,1)
                  CG(2)=SURCOG(nznog(IZ),isur,2)
                  CG(3)=SURCOG(nznog(IZ),isur,3)
                  CALL VECTRN(CG,TSMAT,COG2,IER)
                  call u2pixel(COG2(1),COG2(2),iixc,iiyc)                
                  IF(SVFC(nznog(IZ),isur).EQ.'VERT')THEN

C Check for clipping.
                    call pixel2u(iixc,iiyc+3,xx,xy)
                    call CLIPPT(xx,xy,COG2(3),iclp1)
                    call pixel2u(iixc+ipixw+7,iiyc+3-ipixh,xx,xy)
                    call CLIPPT(xx,xy,COG2(3),iclp2)

C Draw arrow and horizontal line.
                    if (iclp1.eq.0 .and. iclp2.eq.0) then
                      iid1=iixc+3; iid2=iiyc-3
                      if(mmod.eq.8)then
                        call eswline(iixc,iiyc,iid1,iid2)
                      else
                        call eswlinewwc(iixc,iiyc,iid1,iid2)
                      endif
                      iid1=iixc+3; iid2=iiyc+3
                      if(mmod.eq.8)then
                        call eswline(iixc,iiyc,iid1,iid2)
                      else
                        call eswlinewwc(iixc,iiyc,iid1,iid2)
                      endif
                      iid1=iixc+7
                      if(mmod.eq.8)then
                        call eswline(iixc,iiyc,iid1,iiyc)
                      else
                        call eswlinewwc(iixc,iiyc,iid1,iiyc)
                      endif
                      iixc=iixc+7; iiyc=iiyc+3
                      found=.true.
                    else
                      found=.false.
                    endif
                  ELSEIF(SVFC(nznog(IZ),isur).EQ.'CEIL')THEN
                  
C Check for clipping.
                    call pixel2u(iixc-3,iiyc,xx,xy)
                    call CLIPPT(xx,xy,COG2(3),iclp1)
                    call pixel2u(iixc+ipixw+7,iiyc-ipixh,xx,xy)
                    call CLIPPT(xx,xy,COG2(3),iclp2)

C Draw arrow to surface then up and horizontal to the text.
                    if (iclp1.eq.0 .and. iclp2.eq.0) then
                      iid1=iixc+3; iid2=iiyc-3
                      if(mmod.eq.8)then
                        call eswline(iixc,iiyc,iid1,iid2)
                      else
                        call eswlinewwc(iixc,iiyc,iid1,iid2)
                      endif
                      iid1=iixc-3; iid2=iiyc-3
                      if(mmod.eq.8)then
                        call eswline(iixc,iiyc,iid1,iid2)
                      else
                        call eswlinewwc(iixc,iiyc,iid1,iid2)
                      endif
                      iid1=iiyc-5
                      if(mmod.eq.8)then
                        call eswline(iixc,iiyc,iixc,iid1)
                      else
                        call eswlinewwc(iixc,iiyc,iixc,iid1)
                      endif
                      iid1=iiyc-5; iid2=iixc+7; iid3=iiyc-5
                      if(mmod.eq.8)then
                        call eswline(iixc,iid1,iid2,iid3)
                      else
                        call eswlinewwc(iixc,iid1,iid2,iid3)
                      endif
                      iixc=iixc+7
                      found=.true.
                    else
                      found=.false.
                    endif
                  ELSEIF(SVFC(nznog(IZ),isur).EQ.'FLOR')THEN

C Check for clipping.
                    call pixel2u(iixc-3,iiyc+8,xx,xy)
                    call CLIPPT(xx,xy,COG2(3),iclp1)
                    call pixel2u(iixc+ipixw+7,iiyc+8-ipixh,xx,xy)
                    call CLIPPT(xx,xy,COG2(3),iclp2)

C Draw arrow to surface then down and horizontal to the text.
                    if (iclp1.eq.0 .and. iclp2.eq.0) then
                      iid1=iixc+3; iid2=iiyc+3
                      if(mmod.eq.8)then
                        call eswline(iixc,iiyc,iid1,iid2)
                      else
                        call eswlinewwc(iixc,iiyc,iid1,iid2)
                      endif
                      iid1=iixc-3; iid2=iiyc+3
                      if(mmod.eq.8)then
                        call eswline(iixc,iiyc,iid1,iid2)
                      else
                        call eswlinewwc(iixc,iiyc,iid1,iid2)
                      endif
                      iid1=iiyc+5
                      if(mmod.eq.8)then
                        call eswline(iixc,iiyc,iixc,iid1)
                      else
                        call eswlinewwc(iixc,iiyc,iixc,iid1)
                      endif
                      iid1=iiyc+5; iid2=iixc+7; iid3=iiyc+5
                      if(mmod.eq.8)then
                        call eswline(iixc,iid1,iid2,iid3)
                      else
                        call eswlinewwc(iixc,iid1,iid2,iid3)
                      endif
                      iixc=iixc+7; iiyc=iiyc+8
                      found=.true.
                    else
                      found=.false.
                    endif
                  ELSE

C Check for clipping.
                    call CLIPPT(COG2(1),COG2(2),COG2(3),iclp1)
                    call pixel2u(iixc+ipixw,iiyc-ipixh,xx,xy)
                    call CLIPPT(xx,xy,COG2(3),iclp2)
                    if (iclp1.eq.0 .and. iclp2.eq.0) then
                      found=.true.
                    else
                      found=.false.
                    endif
                  ENDIF
                  if (found) then
                    iicol=0
                    if(mmod.eq.8)then
                      call textatxy(iixc,iiyc,stemp,'-',iicol)
                    else
                      call textatxywwc(iixc,iiyc,stemp,'-',iicol)
                    endif
                  endif
                  if (mmod.eq.8) CALL winfnt(IMFS)
                endif
              endif
            endif
          ENDIF
        ENDIF
        if(mmod.eq.8) call forceflush()
 1000 CONTINUE

C Loop back again to see if there is more.
      GOTO 100

      END


C ********* chgazi
C Called from C code in esp-r.c or esru_x.c with an increment or decrement
C to the viewing azimuth (icazi). If not needed then provide a dummy. 

      subroutine chgazi(icazi,ifrlk)
#include "building.h"
#include "prj3dv.h"

C Passed paremter.
      integer icazi
      
      PI = 4.0 * ATAN(1.0)
      tdis= crow(eyem,viewm)
      call LN2AZ(VIEWM(1),VIEWM(2),VIEWM(3),eyem(1),eyem(2),eyem(3),
     &  az,el)
      aznew=az + float(icazi)

      RAD = PI/180.
      RYAZI = aznew*RAD
      RSALT = el*RAD
      Z = tdis*SIN(RSALT)
      XYDIS = tdis*COS(RSALT)
      IF (XYDIS .LT. 1E-6)THEN
        X = 0.
        Y = 0.
      ELSE
        X = XYDIS*SIN(RYAZI)
        Y = XYDIS*COS(RYAZI)
      ENDIF
      EYEM(1)=X+VIEWM(1)
      EYEM(2)=Y+VIEWM(2)
      EYEM(3)=Z+VIEWM(3)
      MODIFYVIEW=.TRUE.
      MODLEN=.TRUE.

C If we are in a freelook, do minimal rendering to save lag.
      if (ifrlk.eq.1) then
        ITZNMP=ITZNM
        ITZNM=1  ! turn off zone names
        ITSNMP=ITSNM
        ITSNM=1  ! turn off surface names
        ITORGP=ITORG
        ITORG=1  ! turn off site origin
        ITGRDP=ITGRD
        ITGRD=1  ! turn off grid
        ITDSPP=ITDSP
        ITDSP=1  ! surfaces only
      endif

      call redraw(IER)

      if (ifrlk.eq.1) then
        ITZNM=ITZNMP
        ITSNM=ITSNMP
        ITORG=ITORGP
        ITGRD=ITGRDP
        ITDSP=ITDSPP
      endif

      return
      end

C ********* chgelev
C Called from C code in esp-r.c or esru_x.c with an increment or decrement
C the viewing elevation (icelev). If not needed then provide a dummy. 

      subroutine chgelev(icelev,ifrlk)
#include "building.h"
#include "prj3dv.h"

C Passed parameter.
      integer icelev
      
      PI = 4.0 * ATAN(1.0)
      tdis= crow(eyem,viewm)
      call ln2az(VIEWM(1),VIEWM(2),VIEWM(3),eyem(1),eyem(2),eyem(3),
     &  az,el)
      elevnew=el + float(icelev)

      RAD = PI/180.
      RYAZI = az*RAD
      RSALT = elevnew*RAD
      Z = tdis*SIN(RSALT)
      XYDIS = tdis*COS(RSALT)
      IF (XYDIS .LT. 1E-6)THEN
        X = 0.
        Y = 0.
      ELSE
        X = XYDIS*SIN(RYAZI)
        Y = XYDIS*COS(RYAZI)
      ENDIF
      EYEM(1)=X+VIEWM(1)
      EYEM(2)=Y+VIEWM(2)
      EYEM(3)=Z+VIEWM(3)
      MODIFYVIEW=.TRUE.
      MODLEN=.TRUE.

C If we are in a freelook, do minimal rendering to save lag.
      if (ifrlk.eq.1) then
        ITZNMP=ITZNM
        ITZNM=1  ! turn off zone names
        ITSNMP=ITSNM
        ITSNM=1  ! turn off surface names
        ITORGP=ITORG
        ITORG=1  ! turn off site origin
        ITGRDP=ITGRD
        ITGRD=1  ! turn off grid
        ITDSPP=ITDSP
        ITDSP=1  ! surfaces only
      endif

      call redraw(IER)

      if (ifrlk.eq.1) then
        ITZNM=ITZNMP
        ITSNM=ITSNMP
        ITORG=ITORGP
        ITGRD=ITGRDP
        ITDSP=ITDSPP
      endif

      return
      end

C ********* chgpan
C Called from C code in esp-r.c or esru_x.c with x and y
C increments to the viewed point. This effectivley implements mouse
C panning. If not needed then provide a dummy. 

      subroutine chgpan(ix,iy,ifrlk)
#include "building.h"
#include "prj3dv.h"

C Passed paremter.
      integer ix,iy
      real vec(3)

C Convert viewed point to screen coordinates.
      CALL VECTRN(VIEWM,TSMAT,vec,IER)

C Scale movement increments by view angle.
      idx=-ix*10
      idy=iy*10

C Subtract pan increments.
      vec(1)=vec(1)+idx
      vec(2)=vec(2)+idy

C Convert back to model coordinates and save.
      CALL VECTRN(vec,STMMAT,VIEWM,IER)

C Make sure view bounds are static, and grid and origin are off.
      if (ITBND.ne.0) ITBND=0

      MODIFYVIEW=.TRUE.
      MODLEN=.TRUE.
      
C If we are in a freelook, do minimal rendering to save lag.
      if (ifrlk.eq.1) then
        ITZNMP=ITZNM
        ITZNM=1  ! turn off zone names
        ITSNMP=ITSNM
        ITSNM=1  ! turn off surface names
        ITORGP=ITORG
        ITORG=1  ! turn off site origin
        ITGRDP=ITGRD
        ITGRD=1  ! turn off grid
        ITDSPP=ITDSP
        ITDSP=1  ! surfaces only
      endif

      call redraw(IER)

      if (ifrlk.eq.1) then
        ITZNM=ITZNMP
        ITSNM=ITSNMP
        ITORG=ITORGP
        ITGRD=ITGRDP
        ITDSP=ITDSPP
      endif

      return
      end

C ********* chgzoom
C Called from C code in esp-r.c or esru_x.c with 1 or 2.
C 1 means zoom out, 2 means zoom in. Do this by adjusting view angle.

      subroutine chgzoom(imode)
#include "building.h"
#include "prj3dv.h"

C Passed paremter.
      integer imode

C Increment view angle by 5 percent.
      if (imode.eq.1) then
        ANG=ANG-ANG*0.05
      elseif (imode.eq.2) then
        ANG=ANG+ANG*0.05
      endif
      
C Make sure view bounds are static.
      if (ITBND.ne.0) ITBND=0

      MODIFYVIEW=.TRUE.
      MODLEN=.TRUE.
      call redraw(IER)

      return
      end

C ********* optview
C Called from C code in esp-r.c or esru_x.c. It sets optimum
C view bounds.

      subroutine optview
#include "building.h"
#include "prj3dv.h"

C Activate optimum view bounds, as well as origin and grid.
      if (ITBND.ne.1) ITBND=1    
      if (ITORG.ne.0) ITORG=0 
      if (ITGRD.ne.0) ITGRD=0

      MODIFYVIEW=.TRUE.
      MODLEN=.TRUE.
      MODBND=.TRUE.
      call redraw(IER)

      return
      end

C ********* chgeye
C Called from C code in esru_x.c with current contents of the ray2_ and
C image_ c structures.  Integers are specially treated based on ifdef.

      subroutine chgeye(EVX,EVY,EVZ,VX,VY,VZ,EAN,JITZNM,JITSNM,JITVNO,
     &   JITOBS,JITVIS,JITVOBJ,JITSNR,JITGRD,JITORG,GDIS,JITBND,JITDSP,
     &   JITHLS,JITHLZ,JITPPSW)
#include "building.h"
#include "prj3dv.h"

C Passed parameters.
      real EVX,EVY,EVZ,VX,VY,VZ,EAN,GDIS

C Depending on computer type set integer size of passed parameters.
#ifdef OSI
      integer JITZNM,JITSNM,JITVNO,JITOBS,JITVIS,JITVOBJ,JITSNR,JITGRD
      integer JITORG,JITBND,JITDSP,JITHLS,JITHLZ,JITPPSW
#else
      integer*8 JITZNM,JITSNM,JITVNO,JITOBS,JITVIS,JITVOBJ,JITSNR
      integer*8 JITGRD,JITORG,JITBND,JITDSP,JITHLS,JITHLZ,JITPPSW
#endif     
      
C Debug
C      write(6,*) 'chgeye passed vals ',EVX,EVY,EVZ,VX,VY,VZ,EAN,
C     &  JITZNM,JITSNM,JITVNO,JITOBS,'snr',JITSNR,JITGRD,JITORG,GDIS,
C     &  JITBND,JITDSP,JITHLS,JITHLZ,JITPPSW
C      write(6,*) 'chgeye com before  ',eyem(1),eyem(2),eyem(3),
C     &  VIEWM(1),VIEWM(2),VIEWM(3),ANG,HANG
C      write(6,*) 'ray2 com ',ITDSP,ITBND,ITEPT,ITZNM,ITSNM,
C     &  ITVNO,ITORG,'snr',ITSNR,ITOBS,ITHLS,ITHLZ,ITGRD,GRDIS,ITPPSW
      eyem(1)=EVX
      eyem(2)=EVY
      eyem(3)=EVZ
      VIEWM(1)=VX
      VIEWM(2)=VY
      VIEWM(3)=VZ
      ANG=EAN
      HANG=ANG/2.0

C To avoid implicit casts also check machine type.
#ifdef OSI
      ITDSP=JITDSP; ITBND=JITBND; ITPPSW=JITPPSW
      ITZNM=JITZNM; ITSNM=JITSNM; ITVNO=JITVNO
      ITORG=JITORG; ITSNR=JITSNR; ITOBS=JITOBS; ITVIS=JITVIS;
      ITVOBJ=JITVOBJ; ITHLS=JITHLS; ITHLZ=JITHLZ; ITGRD=JITGRD
#else
      ITDSP=int(JITDSP); ITBND=int(JITBND); ITPPSW=int(JITPPSW)
      ITZNM=int(JITZNM); ITSNM=int(JITSNM); ITVNO=int(JITVNO)
      ITORG=int(JITORG); ITVIS=int(JITVIS); ITVOBJ=int(JITVOBJ)
      ITSNR=int(JITSNR); ITOBS=int(JITOBS)
      ITHLS=int(JITHLS); ITHLZ=int(JITHLZ); ITGRD=int(JITGRD)
#endif     

      GRDIS=GDIS

C Debug
C      write(6,*) 'chgeye com after  ',eyem(1),eyem(2),eyem(3),
C     &  VIEWM(1),VIEWM(2),VIEWM(3),ANG,HANG
C      write(6,*) 'ray2 aft ',ITDSP,ITBND,ITEPT,ITZNM,ITSNM,
C     &  ITVNO,ITORG,'snr',ITSNR,ITOBS,ITHLS,ITHLZ,ITGRD,GRDIS,ITPPSW

      MODIFYVIEW=.TRUE.
      MODLEN=.TRUE.

C Use same logic as in EVSET.
      if(ITPPSW.eq.0)then

C If perspective view do the normal wireframe calls.
        if(ITHLS.eq.1)then

C If ITHLS set to 1 then ask about which construction.
C But this is disabled because epkmlc tries to use the current
C main menu structure which causes a GTK fault. Need an alternative
C list management for this.
          call usrmsg('Hilight by construction is not yet working',
     &                'this will be added in a later version.','P')
C          CALL EPKMLC(ISEL,'Select a construction to hilight.',' ',IER)
C          ITHLZ=ISEL
C          MODIFYVIEW=.TRUE.
        endif
        CALL INLNST(1)
        call redraw(IER)
        MODLEN=.FALSE.
        MODBND=.FALSE.
        MODIFYVIEW=.FALSE.
      elseif(ITPPSW.eq.1)then
        call PLELEV('P')
      elseif(ITPPSW.eq.2)then
        call PLELEV('S')
      elseif(ITPPSW.eq.3)then
        call PLELEV('E')
      endif

      return
      end

C ********* chgsun
C Called from C code in esru_x.ca and esp-r.c requesting solar view. 
C If isunhour zero then do initial setup of date and time and ask user what
C else they want to do (previous next animation close). If isunhour one
C step to next hour. If isunhour two then present an animation. 

      subroutine chgsun(isunhour)
#include "building.h"
#include "prj3dv.h"

C Passed parameter.
      integer isunhour

      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      COMMON/PREC8/SLAT,SLON
      integer ifs,itfs,imfs
      COMMON/GFONT/IFS,ITFS,IMFS
      COMMON/SPAD/MMOD,LIMIT,LIMTTY

C Local variables
      character outs*124,ETEXT*60,DESCRH*5,DESCRD*5,DESCRJ*5
      integer NTS,IMO,IDO,IJDAY,IT,IER
      logical doanim
      real STIME,X1,Y1,Z1
      DIMENSION MTHNAM(12)
      CHARACTER MTHNAM*3

#ifdef OSI
      integer ivt2,ivt3,ivt4  ! for use with viewtext
#else
      integer*8 ivt2,ivt3,ivt4  ! for use with viewtext
#endif

C Remember some values.
      save STIME,NTS,IMO,IDO  ! remember for next call

      DATA MTHNAM/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug',
     &            'Sep','Oct','Nov','Dec'/

C Draw the zones and any obstructions based on sun
C position. Confirm the lat, long, request day and time,
C compute positon and update the view.
      doanim=.false.
      write(outs,'(A,F6.2,A,F6.2)')
     &  ' Current latitude = ',SLAT,
     &  ', longitude difference from time meridian=',SLON
      call edisp(iuout,' ')
      call edisp(iuout,outs)
      if(isunhour.eq.0)then

C Ask for initial date and time and present view.
        call edisp(iuout,' ')
        call edisp(iuout,' Date and time for the view... ')
        NTS=1
        IMO=6
        IDO=1
        STIME=1.
  67    CALL ASKTIM(2,NTS,IMO,IDO,IJDAY,STIME,IT,IER)
  68    CALL EAZALTS(STIME,IJDAY,SLAT,SLON,ISUNUP,SAZI,SALT)
        IF(SALT.GT.0.0)THEN
          CALL ANGXYZ(SAZI,SALT,X1,Y1,Z1)
          EYEM(1)=X1
          EYEM(2)=Y1
          EYEM(3)=Z1
          ITSNM=1
          MODIFYVIEW=.TRUE.
          MODLEN=.TRUE.
          call redraw(IER)

C Generate a heading for the view.
          CALL EDTIME(STIME,DESCRH,DESCRD,DESCRJ,TIMER)
          WRITE(ETEXT,'(A,I2,1X,A3,A,A5,A,F6.1,A,F5.1)')
     &      ' View: ',IDO,MTHNAM(IMO),' @',DESCRH,
     &      '  azim:',SAZI,' elev:',SALT
          ivt2=2; ivt3=1; ivt4=IFS
          if(mmod.eq.8)then
            call viewtext(ETEXT,ivt2,ivt3,ivt4)
          else
            call viewtextwwc(ETEXT,ivt2,ivt3,ivt4)
          endif

C Ask user what they want to do next.
          CALL EASKMBOX(' View choices: ',' ','previous hour',
     &      'next hour ','animation','cancel',
     &      ' ',' ',' ',' ',IW,0)
          if(IW.eq.1)then
            STIME=STIME-((60.0/FLOAT(NTS))/60.0)
            goto 68
          elseif(IW.eq.2)then
            STIME=STIME+((60.0/FLOAT(NTS))/60.0)
            goto 68
          elseif(IW.eq.3)then
            doanim=.true.
          elseif(IW.eq.4)then
            ITSNM=0
            return
          endif
        else
          CALL EASKMBOX(' View point below horizon:',' ',
     &        'specify another time','exit',
     &        ' ',' ',' ',' ',' ',' ',IW,0)
          if(IW.eq.1)goto 67
          ITSNM=0
          return
        endif
      endif

      if(isunhour.eq.1)then

C User asked for next hour, increment and display unless below horizon.
  69    STIME=STIME+((60.0/FLOAT(NTS))/60.0)
        CALL EAZALTS(STIME,IJDAY,SLAT,SLON,ISUNUP,SAZI,SALT)
        IF(SALT.GT.0.0)THEN
          CALL ANGXYZ(SAZI,SALT,X1,Y1,Z1)
          EYEM(1)=X1
          EYEM(2)=Y1
          EYEM(3)=Z1
          ITSNM=1
          MODIFYVIEW=.TRUE.
          MODLEN=.TRUE.
          call redraw(IER)

C Generate a heading for the view.
          CALL EDTIME(STIME,DESCRH,DESCRD,DESCRJ,TIMER)
          WRITE(ETEXT,'(A,I2,1X,A3,A,A5,A,F6.1,A,F5.1)')
     &      ' View: ',IDO,MTHNAM(IMO),' @',DESCRH,
     &      '  azim:',SAZI,' elev:',SALT
          ivt2=2; ivt3=1; ivt4=IFS
          if(mmod.eq.8)then
            call viewtext(ETEXT,ivt2,ivt3,ivt4)
            call pauses(1)   ! pause before returning
          else
            call viewtextwwc(ETEXT,ivt2,ivt3,ivt4)
          endif
          ITSNM=0
          return
        else
          if(STIME.LT.11.0) goto 69  ! try next hour
          ITSNM=0
          return
        endif
      elseif(isunhour.eq.2)then

C User asked for animation.
        doanim=.true.
      endif

      if(doanim)then
        if(STIME.lt.3.0) STIME = 3.0  ! if not set
 70     CALL EAZALTS(STIME,IJDAY,SLAT,SLON,ISUNUP,SAZI,SALT)
        NTS=4  ! do animation at 15 minute steps
        IF(SALT.GT.0.0)THEN
          CALL ANGXYZ(SAZI,SALT,X1,Y1,Z1)
          EYEM(1)=X1
          EYEM(2)=Y1
          EYEM(3)=Z1
          ITSNM=1
          MODIFYVIEW=.TRUE.
          MODLEN=.TRUE.
          call redraw(IER)

C Generate a heading for the view.
          CALL EDTIME(STIME,DESCRH,DESCRD,DESCRJ,TIMER)
          WRITE(ETEXT,'(A,I2,1X,A3,A,A5,A,F6.1,A,F5.1)')
     &      ' View: ',IDO,MTHNAM(IMO),' @',DESCRH,
     &      '  azim:',SAZI,' elev:',SALT
          ivt2=2
          ivt3=1
          ivt4=IFS
          if(mmod.eq.8)then
            CALL viewtext(ETEXT,ivt2,ivt3,ivt4)
            STIME=STIME+((60.0/FLOAT(NTS))/60.0)
            call forceflush()
            call pauses(1)
          else
            CALL viewtextwwc(ETEXT,ivt2,ivt3,ivt4)
            STIME=STIME+((60.0/FLOAT(NTS))/60.0)
          endif
          goto 70
        else
          if(STIME.LT.11.0)then
            goto 70  ! try next hour
          else
            ITSNM=0
            return
          endif
        endif
      endif

      return
      end

C ********* chgzonpik
C Called from C code in esp-r.c with parameters representing the
C current state of izgfoc and nzg in the gzonpik_ c structure. Sparate
C calls are made to chgzonpikarray is used to fill the nznog array

      subroutine chgzonpik(jizgfoc,jnzg)
#include "building.h"

      COMMON/gzonpik/izgfoc,nzg,nznog(mcom)

#ifdef OSI
      integer jizgfoc,jnzg  ! for use with viewtext
#else
      integer*8 jizgfoc,jnzg  ! for use with viewtext
#endif

C To avoid implicit casts also check machine type.
#ifdef OSI
        izgfoc=jizgfoc
        nzg=jnzg
#else
        izgfoc=int(jizgfoc)
        nzg=int(jnzg)
#endif     

C Debug.
C      write(6,*) 'common gzonpik ',izgfoc,nzg,nznog(1),nznog(2),
C     &    nznog(3),nznog(4),nznog(5)
      return
      end

C ********* chgzonpikarray
C Called from C code in esp-r.c with parameters representing
C one element of nznog array (index and value) in the gzonpik_ c structure.
C This is tedious, but does ensure that array elements are properly passed
C from C to Fortran.

      subroutine chgzonpikarray(jnznog,jnznogv)
#include "building.h"

      COMMON/gzonpik/izgfoc,nzg,nznog(mcom)

#ifdef OSI
      integer jnznog,jnznogv
#else
      integer*8 jnznog,jnznogv
#endif

C To avoid implicit casts also check machine type. Increment the
C index jnznog to account for fortran array indexing starts from
C one and C starts from zero.
#ifdef OSI
      nznog(jnznog+1)=jnznogv
#else
      nznog(int(jnznog)+1)=int(jnznogv)
#endif     
      return
      end

C ******************** SVDSOPT ********************
C Save graphics display options into common block RAY2SV.

      subroutine SVDSOPT
#include "building.h"
#include "prj3dv.h"

      COMMON/RAY2SV/ITDSPsv,ITBNDsv,ITEPTsv,ITZNMsv,ITSNMsv,ITVNOsv,
     &      ITORGsv,ITSNRsv,ITOBSsv,ITVISsv,ITVOBJsv,ITHLSsv,ITHLZsv,
     &      ITGRDsv,GRDISsv,ITPPSWsv

      ITDSPsv=ITDSP
      ITBNDsv=ITBND
      ITEPTsv=ITEPT
      ITZNMsv=ITZNM
      ITSNMsv=ITSNM
      ITVNOsv=ITVNO
      ITORGsv=ITORG
      ITSNRsv=ITSNR
      ITOBSsv=ITOBS
      ITVISsv=ITVIS
      ITVOBJsv=ITVOBJ
      ITHLSsv=ITHLS
      ITHLZsv=ITHLZ
      ITGRDsv=ITGRD
      GRDISsv=GRDIS
      ITPPSWsv=ITPPSW

      return
      end

C ******************** RCDSOPT ********************
C Recover graphics display options from common block RAY2SV.

      subroutine RCDSOPT
#include "building.h"
#include "prj3dv.h"

      COMMON/RAY2SV/ITDSPsv,ITBNDsv,ITEPTsv,ITZNMsv,ITSNMsv,ITVNOsv,
     &      ITORGsv,ITSNRsv,ITOBSsv,ITVISsv,ITVOBJsv,ITHLSsv,ITHLZsv,
     &      ITGRDsv,GRDISsv,ITPPSWsv

      ITDSP=ITDSPsv
      ITBND=ITBNDsv
      ITEPT=ITEPTsv
      ITZNM=ITZNMsv
      ITSNM=ITSNMsv
      ITVNO=ITVNOsv
      ITORG=ITORGsv
      ITSNR=ITSNRsv
      ITOBS=ITOBSsv
      ITVIS=ITVISsv
      ITVOBJ=ITVOBJsv
      ITHLS=ITHLSsv
      ITHLZ=ITHLZsv
      ITGRD=ITGRDsv
      GRDIS=GRDISsv
      ITPPSW=ITPPSWsv

      return
      end

C **************DRAWOBS 
C Checks to see if obstructions are associated with zone (IFOC).
C If so takes common block BG1 and then converts into viewing coords, 
C clips the image and draws it.
C IER is returned as non-zero if there is a problem.
C Obstructions are in global commons so re-scan of files probably not needed.

      SUBROUTINE DRAWOBS(IFOC,ier)
#include "building.h"
#include "model.h"

C geometry.h provides commons G0, G2, GS5, GS6.
#include "geometry.h"
#include "prj3dv.h"

C Passed parameters.
      integer ifoc,ier
      
      COMMON/OUTIN/IUOUT,IUIN,IEOUT
C      COMMON/FILEP/IFIL
      COMMON/SPAD/MMOD,LIMIT,LIMTTY
      integer ifs,itfs,imfs
      COMMON/GFONT/IFS,ITFS,IMFS
      common/appcols/mdispl,nifgrey,ncset,ngset,nzonec

      COMMON/GB1/XB(12),YB(12),ZB(12),JVNB(6,4)

C Local arrays.
      real AX,AY,AZ,BX,BY,BZ,SBBOX,COG1,COG2
      integer IANXT,IBNXT
      DIMENSION  AX(MV),AY(MV),AZ(MV),IANXT(MV)
      DIMENSION  BX(MV),BY(MV),BZ(MV),IBNXT(MV)
      DIMENSION  SBBOX(3,2),COG1(3),COG2(3)

      CHARACTER btemp*12

#ifdef OSI
      integer iix,iiy,iixlast,iiylast
      integer iicol,ibsize        ! passed to textsizeatxy
#else
      integer*8 iix,iiy,iixlast,iiylast
      integer*8 iicol,ibsize       ! passed to textsizeatxy
#endif     

C If libraries and monitor support greyscale remember this.
      if(nifgrey.gt.4)then
        iicol=nifgrey-3
      endif

C Set font for displaying block names.
      if(IFS.eq.0.or.IFS.eq.4)then
        ibsize=0; isize=0
      endif
      if(IFS.eq.1.or.IFS.eq.5)then
        ibsize=1; isize=1
      endif
      if(IFS.eq.2.or.IFS.eq.6)then
        ibsize=2; isize=2
      endif
      if(IFS.eq.3.or.IFS.eq.7)then
        ibsize=2; isize=2
      endif

C If no obstruction for this zone don't bother reading file.
C If iobs() is 1 then scan zone obstructions file. If iob()
C is 2 then re-read the version 1.1 zone geometry file.
      IF(IOBS(ifoc).EQ.0)THEN
        RETURN
      ELSEIF(IOBS(ifoc).EQ.1)THEN
        continue
      ELSEIF(IOBS(ifoc).EQ.2)THEN
        continue
      ENDIF

      DO 301 IB=1,nbobs(ifoc)
        if(BLOCKTYP(ifoc,IB)(1:4).eq.'obs ')then
          CALL CNVBLK(XOB(ifoc,IB),YOB(ifoc,IB),ZOB(ifoc,IB),
     &      DXOB(ifoc,IB),DYOB(ifoc,IB),DZOB(ifoc,IB),
     &      BANGOB(ifoc,IB,1))
        elseif(BLOCKTYP(ifoc,IB)(1:4).eq.'obs3')then
          CALL CNVBLK3A(XOB(ifoc,IB),YOB(ifoc,IB),ZOB(ifoc,IB),
     &      DXOB(ifoc,IB),DYOB(ifoc,IB),DZOB(ifoc,IB),
     &      BANGOB(ifoc,IB,1),BANGOB(ifoc,IB,2),BANGOB(ifoc,IB,3))
        elseif(BLOCKTYP(ifoc,IB)(1:4).eq.'obsp')then
          call CNVBLKP(ifoc,IB) ! convert obsp type.
        endif

C If block names are to be displayed place near first vertex.
C Transform the location into eyepoint and then screen coords.
        IF(ITSNM.EQ.0)THEN
          COG1(1)=XB(1)
          COG1(2)=YB(1)
          COG1(3)=ZB(1)
          CALL VECTRN(COG1,TSMAT,COG2,IER)
          call u2pixel(COG2(1),COG2(2),iix,iiy)
          btemp='           '
          write(btemp,'(a)') 
     &      BLOCKNAME(IFOC,IB)(1:lnblnk(BLOCKNAME(IFOC,IB)))
          if(mmod.eq.8)then
            call winfnt(isize)
            call CLIPST(btemp,iix,iiy,COG2(3),iclp)
            if (iclp.eq.0) then
              call textatxy(iix,iiy,btemp,'-',iicol)
              iicol=0
              call winscl('-',iicol)
            endif
            call winfnt(IFS)
          else
            call textsizeatxywwc(iix,iiy,btemp,ibsize,'-',iicol)
          endif
        ENDIF

C If editing an obstruction, highlight it.
        if(ITOBS.eq.IB)then
          IPEN= -305
        else
          IPEN= 1
        endif

C Continue on with surfaces in the block.
C Copy polygon for obstruction surface I into structure AX AY AZ
C where IANXT = Index of 'next' point. Clear IANXT & IBNXT arrays.
        do 1101 J = 1,MV
          IANXT(J) = 0
          IBNXT(J) = 0
 1101   continue

C There are 4 initial vertices per block side and 6 sides.
        NAP=4
        DO 1000 I=1,6
          DO 1100 J = 1,NAP
            K = JVNB(I,J)
            AX(J) = XB(K)
            AY(J) = YB(K)
            AZ(J) = ZB(K)
            IANXT(J) = J + 1
 1100     CONTINUE
          IANXT(4) = 1

C Transform block polygon to EYE co-ordinates where
C TEMAT = Model to Eye Matrix
          CALL MATPOLS(NAP,AX,AY,AZ,IANXT,TSMAT,SBBOX,
     &                 NBP,BX,BY,BZ,IBNXT,IERR)

C Generate clipping flags and clip geometry in eye coords.
C If ISTAT =  0 : totally inside frustrum
C If ISTAT =  1 : totally outside frustrum
C If ISTAT = -1 : straddles frustrum
          IBPNT=1
          call CLIPSUR(NBP,BX,BY,BZ,ISTAT)
          IF (ISTAT .EQ. 1) THEN
            GOTO  1000
          ELSEIF (ISTAT .EQ. -1) THEN
            CALL CUTSUR(NB1,NBP,IBPNT,BX,BY,BZ,IBNXT,ISTAT)
          else
            NB1=1
          ENDIF
            
C If there are no edges to draw then jump.
          if(NB1.eq.0)goto 1000

C Multiply each vertex by ETSMAT (Eye to Screen Matrix) to transform
C to screen coords. First move to the beginning point (as passed from CUTSUR).
          IP = IABS(IBPNT)
          call u2pixel(BX(IP),BY(IP),iix,iiy)
          iixlast=iix   ! cast for edline call
          iiylast=iiy

C If IP1 is ever negative then jump to a later point.
          IP1 = IBNXT(IP)
          if(IP1.lt.1)IP1 = IBNXT(IP+1)
          if(IP1.lt.1)IP1 = IBNXT(IP+2)

 351      CONTINUE
          call u2pixel(BX(IP1),BY(IP1),iix,iiy)
          if(mmod.eq.8)then
            if(IPEN.eq.1)then
              call eswline(iixlast,iiylast,iix,iiy)
            elseif(IPEN.eq.-305)then
              call edwline(iixlast,iiylast,iix,iiy)
            else
              call eswline(iixlast,iiylast,iix,iiy)
            endif
          else
            if(IPEN.eq.1)then
              call eswlinewwc(iixlast,iiylast,iix,iiy)
            elseif(IPEN.eq.-305)then
              call edwlinewwc(iixlast,iiylast,iix,iiy)
            else
              call eswlinewwc(iixlast,iiylast,iix,iiy)
            endif
          endif
          if(mmod.eq.8) call forceflush()

C Repeat until next point is the first, remember position and skip
C past -1 values of IBNXT, save IP1 as IPT for testing.
          IPT = IP1
          IP1 = IBNXT(IP1)
          if(IP1.lt.1)IP1 = IBNXT(IPT+1)
          if(IP1.lt.1)IP1 = IBNXT(IPT+2)
          iixlast=iix   ! cast for edline call
          iiylast=iiy
          IF(IP1 .NE. IP) GOTO 351

C Complete the polygon.
          call u2pixel(BX(IP),BY(IP),iix,iiy)
          if(mmod.eq.8)then
            if(IPEN.eq.1)then
              call eswline(iixlast,iiylast,iix,iiy)
            elseif(IPEN.eq.-305)then
              call edwline(iixlast,iiylast,iix,iiy)
            else
              call eswline(iixlast,iiylast,iix,iiy)
            endif
          else
            if(IPEN.eq.1)then
              call eswlinewwc(iixlast,iiylast,iix,iiy)
            elseif(IPEN.eq.-305)then
              call edwlinewwc(iixlast,iiylast,iix,iiy)
            else
              call eswlinewwc(iixlast,iiylast,iix,iiy)
            endif
          endif
          if(mmod.eq.8) call forceflush()
 1000   CONTINUE
  301 CONTINUE

      RETURN
      END

C ************** DRAWVIS
C Draws any visual entities within zone (IFOC) as slightly grey lines.
C It converts visual entities from V1.1 Geometry file into common block BG1 
C and then converts into viewing coords, clips the image and draws it.
C If ITVIS (from prj3dv.h) is non-zero that entity is highlighted.
C If ITVIS is zero but ITVOBJ is non-zero then we need to highlight
C each of the entities associated with the compound object.
C IER is returned as non-zero if there is a problem.

      SUBROUTINE DRAWVIS(IFOC,ier)
#include "building.h"
#include "model.h"

C geometry.h provides commons G0, G2, GS5, GS6.
#include "geometry.h"
#include "prj3dv.h"

C Passed parameters.
      integer ifoc,ier
      
      COMMON/OUTIN/IUOUT,IUIN,IEOUT
C      COMMON/FILEP/IFIL
      COMMON/SPAD/MMOD,LIMIT,LIMTTY
      integer ifs,itfs,imfs
      COMMON/GFONT/IFS,ITFS,IMFS
      common/appcols/mdispl,nifgrey,ncset,ngset,nzonec

      COMMON/GB1/XB(12),YB(12),ZB(12),JVNB(6,4)

C Local arrays.
      real AX,AY,AZ,BX,BY,BZ,SBBOX,COG1,COG2
      integer IANXT,IBNXT
      DIMENSION  AX(MV),AY(MV),AZ(MV),IANXT(MV)
      DIMENSION  BX(MV),BY(MV),BZ(MV),IBNXT(MV)
      DIMENSION  SBBOX(3,2),COG1(3),COG2(3)
      integer ihighlight(MB)

      CHARACTER btemp*12
      logical greyok

#ifdef OSI
      integer iix,iiy,iixlast,iiylast,iicol,ibsize
#else
      integer*8 iix,iiy,iixlast,iiylast,iicol,ibsize
#endif     

C If libraries and monitor support greyscale then reset forground
C to 50% grey when drawing visual entities.
      greyok=.false.
      if(nifgrey.gt.4)then
        iicol=nifgrey-3
        if(mmod.eq.8) call winscl('i',iicol)
        greyok=.true.
      endif

C Set font for displaying block names.
      if(IFS.eq.0.or.IFS.eq.4)then
        ibsize=0; isize=0
      endif
      if(IFS.eq.1.or.IFS.eq.5)then
        ibsize=1; isize=1
      endif
      if(IFS.eq.2.or.IFS.eq.6)then
        ibsize=2; isize=2
      endif
      if(IFS.eq.3.or.IFS.eq.7)then
        ibsize=2; isize=2
      endif

C      write(6,*) 'drawvis ',ifoc,nbvis(ifoc),itvis,NBVOBJ(ifoc),itvobj

C Clear local list of possibly associated visual entities for compound.
      do loop=1,MB
        ihighlight(loop)=0
      enddo

C Visual entities attributes will be in memory (done during intial
C scan of the model.
C If no visaul entities for this zone don't bother.
      if(nbvis(ifoc).eq.0.and.NBVOBJ(ifoc).eq.0) RETURN

C If there are any compound objects and ITVOBJ is non-zero.
      if(NBVOBJ(ifoc).gt.0.and.ITVOBJ.gt.0)then
        DO I=1,NBVOBJ(ifoc)
          if(I.eq.ITVOBJ)then  ! matching focus compound
C Debug      WRITE(6,'(I3,4a)') I,' ',
C     &        VOBJNAME(ifoc,I)(1:lnblnk(VOBJNAME(ifoc,I))),'  ',
C     &        VOBJDESC(ifoc,I)(1:lnblnk(VOBJDESC(ifoc,I)))
            do J=1,NBVOBJLIST(ifoc,I)
              do K=1,nbvis(ifoc)
                if(VOBJLIST(ifoc,I,J)(1:12).eq.
     &             VISNAME(ifoc,K)(1:12))then
                  nbhighlight=nbhighlight+1
                  ihighlight(k)=1  ! mark it to highlight
                endif
              enddo  ! of K
            enddo  ! of J
          endif
        ENDDO    ! of I
      endif

C Draw all of the visual entities and highlight if ITVIS has
C been set.
      DO 301 IB=1,nbvis(ifoc)
        if(VISTYP(ifoc,IB)(1:4).eq.'vis ')then
          CALL CNVBLK(XOV(ifoc,IB),YOV(ifoc,IB),ZOV(ifoc,IB),
     &      DXOV(ifoc,IB),DYOV(ifoc,IB),DZOV(ifoc,IB),
     &      BANGOV(ifoc,IB,1))
        elseif(VISTYP(ifoc,IB)(1:4).eq.'vis3')then
          CALL CNVBLK3A(XOV(ifoc,IB),YOV(ifoc,IB),ZOV(ifoc,IB),
     &      DXOV(ifoc,IB),DYOV(ifoc,IB),DZOV(ifoc,IB),
     &      BANGOV(ifoc,IB,1),BANGOV(ifoc,IB,2),BANGOV(ifoc,IB,3))
        elseif(VISTYP(ifoc,IB)(1:4).eq.'visp')then
          call CNVVISP(ifoc,IB) ! convert visp type.
        endif

C If visual names are to be displayed place near first vertex.
C Transform the location into eyepoint and then screen coords.
        IF(ITSNM.EQ.0)THEN
          COG1(1)=XB(1)
          COG1(2)=YB(1)
          COG1(3)=ZB(1)
          CALL VECTRN(COG1,TSMAT,COG2,IER)
          call u2pixel(COG2(1),COG2(2),iix,iiy)
          btemp='           '
          write(btemp,'(a)') 
     &      VISNAME(IFOC,IB)(1:lnblnk(VISNAME(IFOC,IB)))
          if(mmod.eq.8)then
            call winfnt(isize)
            call CLIPST(btemp,iix,iiy,COG2(3),iclp)
            if (iclp.eq.0) then
              call textatxy(iix,iiy,btemp,'-',iicol)
            endif
            call winfnt(IFS)
          else
            call textsizeatxywwc(iix,iiy,btemp,ibsize,'-',iicol)
          endif
        ENDIF

C If editing an visual entity, highlight it. ITVIS would have been
C set to the index of the visual entity before this subroutine was
C called
        if(ITVIS.eq.IB.or.ihighlight(ib).eq.1)then
          IPEN= -305
        else
          IPEN= 1
        endif

C Continue on with surfaces in the block.
C Copy polygon for visual surface I into structure AX AY AZ
C where IANXT = Index of 'next' point. Clear IANXT & IBNXT arrays.
        do 1101 J = 1,MV
          IANXT(J) = 0
          IBNXT(J) = 0
 1101   continue

C There are 4 initial vertices per block side and 6 sides.
        NAP=4
        DO 1000 I=1,6
          DO 1100 J = 1,NAP
            K = JVNB(I,J)
            AX(J) = XB(K)
            AY(J) = YB(K)
            AZ(J) = ZB(K)
            IANXT(J) = J + 1
 1100     CONTINUE
          IANXT(4) = 1

C Transform block polygon to EYE co-ordinates where
C TEMAT = Model to Eye Matrix
          CALL MATPOLS(NAP,AX,AY,AZ,IANXT,TSMAT,SBBOX,
     &                 NBP,BX,BY,BZ,IBNXT,IERR)

C Generate clipping flags and clip geometry in eye coords.
C If ISTAT =  0 : totally inside frustrum
C If ISTAT =  1 : totally outside frustrum
C If ISTAT = -1 : straddles frustrum
          call CLIPSUR(NBP,BX,BY,BZ,ISTAT)
          IBPNT=1
          IF (ISTAT .EQ. 1) THEN
            GOTO  1000
          ELSEIF (ISTAT .EQ. -1) THEN
            CALL CUTSUR(NB1,NBP,IBPNT,BX,BY,BZ,IBNXT,ISTAT)
          else
            NB1=1
          ENDIF
            
C If there are no edges to draw then jump.
          if(NB1.eq.0)goto 1000

C Multiply each vertex by ETSMAT (Eye to Screen Matrix) to transform
C to screen coords. First move to the beginning point (as passed from CUTSUR).
          IP = IABS(IBPNT)
          call u2pixel(BX(IP),BY(IP),iix,iiy)
          iixlast=iix   ! cast for edline call
          iiylast=iiy

C If IP1 is ever negative then jump to a later point.
          IP1 = IBNXT(IP)
          if(IP1.lt.1)IP1 = IBNXT(IP+1)
          if(IP1.lt.1)IP1 = IBNXT(IP+2)

 351      CONTINUE
          call u2pixel(BX(IP1),BY(IP1),iix,iiy)

          if(mmod.eq.8)then
            if(IPEN.eq.1)then
              call eswline(iixlast,iiylast,iix,iiy)
            elseif(IPEN.eq.-305)then
              call edwline(iixlast,iiylast,iix,iiy)
            else
              call eswline(iixlast,iiylast,iix,iiy)
            endif
          else
            if(IPEN.eq.1)then
              call eswlinewwc(iixlast,iiylast,iix,iiy)
            elseif(IPEN.eq.-305)then
              call edwlinewwc(iixlast,iiylast,iix,iiy)
            else
              call eswlinewwc(iixlast,iiylast,iix,iiy)
            endif
          endif
          if(mmod.eq.8) call forceflush()

C Repeat until next point is the first, remember position and skip
C past -1 values of IBNXT, save IP1 as IPT for testing.
          IPT = IP1
          IP1 = IBNXT(IP1)
          if(IP1.lt.1)IP1 = IBNXT(IPT+1)
          if(IP1.lt.1)IP1 = IBNXT(IPT+2)
          iixlast=iix   ! cast for edline call
          iiylast=iiy
          IF(IP1 .NE. IP) GOTO 351

C Complete the polygon.
          call u2pixel(BX(IP),BY(IP),iix,iiy)

          if(mmod.eq.8)then
            if(IPEN.eq.1)then
              call eswline(iixlast,iiylast,iix,iiy)
            elseif(IPEN.eq.-305)then
              call edwline(iixlast,iiylast,iix,iiy)
            else
              call eswline(iixlast,iiylast,iix,iiy)
            endif
          else
            if(IPEN.eq.1)then
              call eswlinewwc(iixlast,iiylast,iix,iiy)
            elseif(IPEN.eq.-305)then
              call edwlinewwc(iixlast,iiylast,iix,iiy)
            else
              call eswlinewwc(iixlast,iiylast,iix,iiy)
            endif
          endif
          if(mmod.eq.8) call forceflush()
 1000   CONTINUE
  301 CONTINUE

C Reset color before returning.
      if(greyok)then
        iicol=0
        if(mmod.eq.8) call winscl('-',iicol)
      endif
      if(mmod.eq.8) call forceflush()

      RETURN
      END

C **************** DRWSEN
C Loads the geometric description within an MRT sensor file
C into common block BG1 and then converts into viewing coords,
C clips the image and draws it.
C IER is returned as non-zero if there is a problem.

      SUBROUTINE DRWSEN(iz,ier)
#include "building.h"
#include "prj3dv.h"
#include "geometry.h"

C Passed parameter.
      integer ier
      
      COMMON/SPAD/MMOD,LIMIT,LIMTTY
      integer ifs,itfs,imfs
      COMMON/GFONT/IFS,ITFS,IMFS
      common/appcols/mdispl,nifgrey,ncset,ngset,nzonec

      COMMON/GB1/XB(12),YB(12),ZB(12),JVNB(6,4)

C Local arrays.
      real AX,AY,AZ,BX,BY,BZ,SBBOX,COG1,COG2
      integer IANXT,IBNXT
      DIMENSION  AX(MV),AY(MV),AZ(MV),IANXT(MV)
      DIMENSION  BX(MV),BY(MV),BZ(MV),IBNXT(MV)
      DIMENSION  SBBOX(3,2),COG1(3),COG2(3)
      logical greyok
#ifdef OSI
      integer iix,iiy,iixlast,iiylast
      integer iicol,ibsize        ! passed to textsizeatxy
#else
      integer*8 iix,iiy,iixlast,iiylast
      integer*8 iicol,ibsize      ! passed to textsizeatxy
#endif     

      CHARACTER temp*7
      logical isfront

C If not in graphic mode return.
      if(MMOD.lt.8)return

C If libraries and monitor support greyscale remember this.
      if(nifgrey.gt.4)then
        iicol=nifgrey-3
      endif

C Set font for displaying block names.
      if(IFS.eq.0.or.IFS.eq.4)then
        ibsize=0; isize=0
      endif
      if(IFS.eq.1.or.IFS.eq.5)then
        ibsize=1; isize=1
      endif
      if(IFS.eq.2.or.IFS.eq.6)then
        ibsize=2; isize=2
      endif
      if(IFS.eq.3.or.IFS.eq.7)then
        ibsize=2; isize=2
      endif

      if(NCUB(iz).eq.0)return
      DO 301 IB=1,NCUB(iz)
        CALL CNVBLK(XOC(IB),YOC(IB),ZOC(IB),DXC(IB),DYC(IB),DZC(IB),
     &    CANG(IB))

C Place sensor names near first vertex.
C Transform the location into eyepoint and then screen coords.
        COG1(1)=XB(1)
        COG1(2)=YB(1)
        COG1(3)=ZB(1)
        CALL VECTRN(COG1,TSMAT,COG2,IER)
        call u2pixel(COG2(1),COG2(2),iix,iiy)
        write(temp,'(A)')CUBN(IB)
        if(mmod.eq.8)then
          call winfnt(isize)
          call CLIPST(temp,iix,iiy,COG2(3),iclp)
          if (iclp.eq.0) then
            call textatxy(iix,iiy,temp,'-',iicol)
          endif
          call winfnt(IFS)
        else
          call textsizeatxywwc(iix,iiy,temp,ibsize,'-',iicol)
        endif
        iicol=0
        if(mmod.eq.8) call winscl('-',iicol)

C Continue on with surfaces in the sensor.
C Copy polygon for sensor surface I into structure AX AY AZ
C where IANXT = Index of 'next' point. Clear IANXT & IBNXT arrays.
        do 1101 J = 1,MV
          IANXT(J) = 0
          IBNXT(J) = 0
 1101   continue

C There are 4 initial vertices per block side and 6 sides.
        NAP=4
        DO 1000 I=1,6

c I=3 is the north side on an unrotated cube, which is the way the
C sensor points (for the purposes of radiance viewpoints).
          if (I.eq.3) then
            isfront=.true.
          else
            isfront=.false.
          endif

          DO 1100 J = 1,NAP
            K = JVNB(I,J)
            AX(J) = XB(K)
            AY(J) = YB(K)
            AZ(J) = ZB(K)
            IANXT(J) = J + 1
 1100     CONTINUE
          IANXT(4) = 1

C Transform block polygon to EYE co-ordinates.
          CALL MATPOLS(NAP,AX,AY,AZ,IANXT,TSMAT,SBBOX,
     &                 NBP,BX,BY,BZ,IBNXT,IERR)

C Generate clipping flags and clip geometry in eye coords.
C If ISTAT =  0 : totally inside frustrum
C If ISTAT =  1 : totally outside frustrum
C If ISTAT = -1 : straddles frustrum
          call CLIPSUR(NBP,BX,BY,BZ,ISTAT)
          IBPNT=1
          IF (ISTAT .EQ. 1) THEN
            GOTO  1000
          ELSEIF (ISTAT .EQ. -1) THEN
            CALL CUTSUR(NB1,NBP,IBPNT,BX,BY,BZ,IBNXT,ISTAT)
          else
            NB1=1
          ENDIF

C If there are no edges to draw then jump.
          if(NB1.eq.0)goto 1000

C Multiply each vertex by ETSMAT (Eye to Screen Matrix).
C First move to the beginning (as passed back from CUTSUR).
          IP = IABS(IBPNT)
          call u2pixel(BX(IP),BY(IP),iix,iiy)
          iixlast=iix   ! cast for edline call
          iiylast=iiy

C If IP1 is ever negative then jump to a later point.
          IP1 = IBNXT(IP)
          if(IP1.lt.1)IP1 = IBNXT(IP+1)
          if(IP1.lt.1)IP1 = IBNXT(IP+2)

 351      CONTINUE
          call u2pixel(BX(IP1),BY(IP1),iix,iiy)
          if (isfront) then
            if(mmod.eq.8)then
              call edwline(iixlast,iiylast,iix,iiy)
            else
              call edwlinewwc(iixlast,iiylast,iix,iiy)
            endif
          else
            if(mmod.eq.8)then
              call eswline(iixlast,iiylast,iix,iiy)
            else
              call eswlinewwc(iixlast,iiylast,iix,iiy)
            endif
          endif
          if(mmod.eq.8) call forceflush()

C Repeat until next point is the first, remember position and skip
C past any -1 values of IBNXT, save IP1 as IPT for testing.
          IPT = IP1
          IP1 = IBNXT(IP1)
          if(IP1.lt.1)IP1 = IBNXT(IPT+1)
          if(IP1.lt.1)IP1 = IBNXT(IPT+2)
          iixlast=iix   ! cast for edline call
          iiylast=iiy
          IF(IP1 .NE. IP) GOTO 351

C Complete the polygon.
          call u2pixel(BX(IP),BY(IP),iix,iiy)
          if (isfront) then
            if(mmod.eq.8)then
              call edwline(iixlast,iiylast,iix,iiy)
            else
              call edwlinewwc(iixlast,iiylast,iix,iiy)
            endif
          else
            if(mmod.eq.8)then
              call eswline(iixlast,iiylast,iix,iiy)
            else
              call eswlinewwc(iixlast,iiylast,iix,iiy)
            endif
          endif
          if(mmod.eq.8) call forceflush()
 1000   CONTINUE
  301 CONTINUE

      RETURN
      END


C **************** DRWBB
C Previews a bounding box as described by current contents of
C common block GB1. Also shows labels passed. Similar to code in the
C subroutine DRWSEN.
C IER is returned as non-zero if there is a problem.

      SUBROUTINE DRWBB(bbname,ier)
#include "building.h"
#include "prj3dv.h"

C Passed parameter.
      character bbname*12  ! preview name
      integer ier
      
      COMMON/SPAD/MMOD,LIMIT,LIMTTY
      integer ifs,itfs,imfs
      COMMON/GFONT/IFS,ITFS,IMFS
      common/appcols/mdispl,nifgrey,ncset,ngset,nzonec

      COMMON/GB1/XB(12),YB(12),ZB(12),JVNB(6,4)

C Local arrays.
      real AX,AY,AZ,BX,BY,BZ,SBBOX,COG1,COG2
      integer IANXT,IBNXT
      DIMENSION  AX(MV),AY(MV),AZ(MV),IANXT(MV)
      DIMENSION  BX(MV),BY(MV),BZ(MV),IBNXT(MV)
      DIMENSION  SBBOX(3,2),COG1(3),COG2(3)
#ifdef OSI
      integer iix,iiy,iixlast,iiylast
      integer iicol,ibsize         ! passed to textsizeatxy
#else
      integer*8 iix,iiy,iixlast,iiylast
      integer*8 iicol,ibsize       ! passed to textsizeatxy
#endif     

      CHARACTER temp*13

C If not in graphic mode return.
      if(MMOD.lt.8)return

C If libraries and monitor support greyscale remember this.
      if(nifgrey.gt.4)then
        iicol=nifgrey-3
      endif

C Set font for displaying bb names.
      if(IFS.eq.0.or.IFS.eq.4)then
        ibsize=0; isize=0
      endif
      if(IFS.eq.1.or.IFS.eq.5)then
        ibsize=1; isize=1
      endif
      if(IFS.eq.2.or.IFS.eq.6)then
        ibsize=2; isize=2
      endif
      if(IFS.eq.3.or.IFS.eq.7)then
        ibsize=2; isize=2
      endif

C Place sensor names near first vertex.
C Transform the location into eyepoint and then screen coords.
      COG1(1)=XB(1)
      COG1(2)=YB(1)
      COG1(3)=ZB(1)
      CALL VECTRN(COG1,TSMAT,COG2,IER)
      call u2pixel(COG2(1),COG2(2),iix,iiy)
      write(temp,'(A)')bbname
      if(mmod.eq.8)then
        call winfnt(isize)
        call CLIPST(temp,iix,iiy,COG2(3),iclp)
        if (iclp.eq.0) then
          call textatxy(iix,iiy,temp,'-',iicol)
        endif
        call winfnt(IFS)
      else
        call textsizeatxywwc(iix,iiy,temp,ibsize,'-',iicol)
      endif
      iicol=0
      if(mmod.eq.8) call winscl('-',iicol)

C Continue on with surfaces in the sensor.
C Copy polygon for sensor surface I into structure AX AY AZ
C where IANXT = Index of 'next' point. Clear IANXT & IBNXT arrays.
      do 1101 J = 1,MV
        IANXT(J) = 0
        IBNXT(J) = 0
 1101 continue

C There are 4 initial vertices per block side and 6 sides.
      NAP=4
      DO 1000 I=1,6
        DO 1100 J = 1,NAP
          K = JVNB(I,J)
          AX(J) = XB(K)
          AY(J) = YB(K)
          AZ(J) = ZB(K)
          IANXT(J) = J + 1
 1100   CONTINUE
        IANXT(4) = 1

C Transform block polygon to EYE co-ordinates.
        CALL MATPOLS(NAP,AX,AY,AZ,IANXT,TSMAT,SBBOX,
     &               NBP,BX,BY,BZ,IBNXT,IERR)

C Generate clipping flags and clip geometry in eye coords.
C If ISTAT =  0 : totally inside frustrum
C If ISTAT =  1 : totally outside frustrum
C If ISTAT = -1 : straddles frustrum
        IBPNT=1
        call CLIPSUR(NBP,BX,BY,BZ,ISTAT)
        IF (ISTAT .EQ. 1) THEN
          GOTO  1000
        ELSEIF (ISTAT .EQ. -1) THEN
          CALL CUTSUR(NB1,NBP,IBPNT,BX,BY,BZ,IBNXT,ISTAT)
        else
          NB1=1
        ENDIF

C If there are no edges to draw then jump.
        if(NB1.eq.0)goto 1000

C Multiply each vertex by ETSMAT (Eye to Screen Matrix).
C First move to the beginning (as passed back from CUTSUR).
        IP = IABS(IBPNT)
        call u2pixel(BX(IP),BY(IP),iix,iiy)
        iixlast=iix   ! cast for edline call
        iiylast=iiy

C If IP1 is ever negative then jump to a later point.
        IP1 = IBNXT(IP)
        if(IP1.lt.1)IP1 = IBNXT(IP+1)
        if(IP1.lt.1)IP1 = IBNXT(IP+2)

 351    CONTINUE
        call u2pixel(BX(IP1),BY(IP1),iix,iiy)
        if(mmod.eq.8)then
          call eswline(iixlast,iiylast,iix,iiy)
          call forceflush()
        else
          call eswlinewwc(iixlast,iiylast,iix,iiy)
        endif

C Repeat until next point is the first, remember position and skip
C past any -1 values of IBNXT, save IP1 as IPT for testing.
        IPT = IP1
        IP1 = IBNXT(IP1)
        if(IP1.lt.1)IP1 = IBNXT(IPT+1)
        if(IP1.lt.1)IP1 = IBNXT(IPT+2)
        iixlast=iix   ! cast for edline call
        iiylast=iiy
        IF(IP1 .NE. IP) GOTO 351

C Complete the polygon.
        call u2pixel(BX(IP),BY(IP),iix,iiy)
        if(mmod.eq.8)then
          call eswline(iixlast,iiylast,iix,iiy)
          call forceflush()
        else
          call eswlinewwc(iixlast,iiylast,iix,iiy)
        endif
 1000 CONTINUE

      RETURN
      END

C **************** withinview
C Checks if a dummy surface at xyz is within viewing cone.
C Derived from  DRWSEN.
C IER is returned as non-zero if there is a problem.

      SUBROUTINE withinview(tx,ty,tz,cansee,ier)
#include "building.h"
#include "prj3dv.h"

C Passed parameter.
      real tx,ty,tz
      logical cansee
      integer ier
      
      COMMON/SPAD/MMOD,LIMIT,LIMTTY

C Local arrays.
      real AX,AY,AZ,BX,BY,BZ,SBBOX
      integer IANXT,IBNXT
      DIMENSION  AX(MV),AY(MV),AZ(MV),IANXT(MV)
      DIMENSION  BX(MV),BY(MV),BZ(MV),IBNXT(MV)
      DIMENSION  SBBOX(3,2)

C If not in graphic mode return.
      if(MMOD.lt.8)return
      cansee=.false.  ! assume we cannot see.

      tx2=tx+0.3; tx3=tx-0.3
      ty2=ty; ty3=ty
      tz2=tz+0.3; tz3=tz+0.3

C Continue with small triangle offset from tx ty tz.
C Copy polygon for sensor surface I into structure AX AY AZ
C where IANXT = Index of 'next' point. Clear IANXT & IBNXT arrays.
      do 1101 J = 1,3
        IANXT(J) = 0
        IBNXT(J) = 0
 1101 continue

C There are 3 initial vertices in polygon.
      NAP=3
      AX(1) = tx; AY(1) = ty; AZ(1) = tz
      AX(2) = tx2; AY(1) = ty2; AZ(1) = tz2
      AX(3) = tx3; AY(1) = ty3; AZ(1) = tz3
      IANXT(1) = 2; IANXT(2) = 3;IANXT(3) = 1

C Transform polygon to EYE co-ordinates.
      CALL MATPOLS(NAP,AX,AY,AZ,IANXT,TSMAT,SBBOX,
     &             NBP,BX,BY,BZ,IBNXT,IERR)

C Generate clipping flags and clip geometry in eye coords.
C If ISTAT =  0 : totally inside frustrum
C If ISTAT =  1 : totally outside frustrum
C If ISTAT = -1 : straddles frustrum
      call CLIPSUR(NBP,BX,BY,BZ,ISTAT)
      IF (ISTAT .EQ. 1) THEN
        cansee=.false.
      ELSEIF (ISTAT .EQ. -1) THEN
        cansee=.true.
      ELSEIF (ISTAT .EQ. 0) THEN
        cansee=.true.
      ENDIF
      ier=IERR

      RETURN
      END  ! of withinview

C ****** EGRNDR
C Draws ground topology held in GT5 and grndpl commons.
C IER is returned as non-zero if there has been a problem.

      SUBROUTINE EGRNDR(IER)
#include "building.h"
#include "prj3dv.h"

C Passed parameter.
      integer ier
      
      COMMON/SPAD/MMOD,LIMIT,LIMTTY
      integer ifs,itfs,imfs
      COMMON/GFONT/IFS,ITFS,IMFS
      common/grndpl/NGT,NGTV,XGT(MGTV),YGT(MGTV),ZGT(MGTV),JGVN(MGRT,8),
     &  NGVER(MGRT),IVEDGE(MGRT)
      COMMON/GT5/GSNAME(MGRT),GMLCN(MGRT)
      COMMON/RAY6G/LINSTYG(MGRT)

C Local arrays.
      real AX,AY,AZ,BX,BY,BZ,SBBOX,COG1,COG2
      integer IANXT,IBNXT
      DIMENSION  AX(MPL),AY(MPL),AZ(MPL),IANXT(MPL)
      DIMENSION  BX(MPL),BY(MPL),BZ(MPL),IBNXT(MPL)
      DIMENSION  SBBOX(3,2),COG1(3),COG2(3)

      CHARACTER GMLCN*32, GSNAME*6,stemp*7

C Local variables to pass to edline.
#ifdef OSI
      integer iix,iiy,iixlast,iiylast,iixc,iiyc,iid2,iicol
#else
      integer*8 iix,iiy,iixlast,iiylast,iixc,iiyc,iid2,iicol
#endif     

C Save the current font and switch to a smaller font for displaying
C surface names etc. Ixc,iyc keep track of surface/pixel COG.
      if(MMOD.lt.8)return
      issize=4  ! for surfaces
      LD=3
      if(ngtv.eq.0)return
      IF(ITVNO.EQ.0)then
        do 42 iv = 1,NGTV
          COG1(1)=XGT(IV)
          COG1(2)=YGT(IV)
          COG1(3)=ZGT(IV)
          CALL VECTRN(COG1,TSMAT,COG2,IER)
          call CLIPPT(COG2(1),COG2(2),COG2(3),iclp)
          if (iclp.eq.0) then
            call u2pixel(COG2(1),COG2(2),iix,iiy)
            CALL ecirc(iix,iiy,2,1)
            CALL VERTLBL(iix,iiy,COG2(3),IV,ier)
          endif
  42    continue
      endif
      if(ngt.eq.0)return
      do 43 igs = 1, NGT
        iixc=0
        iiyc=0

C Copy polygon for Surface I into structure A for compatibility with
C viewer format 'holes' in surfaces.
C NAP   = Number of vertex points on surface
C IANXT = Index of 'next' point
C IAPNT = Pointer to first vertex of polygon
        IAPNT = 1
        NAP = NGVER(igs)
        DO 1100 J = 1,NGVER(igs)
          K = JGVN(igs,J)
          AX(J) = XGT(K)
          AY(J) = YGT(K)
          AZ(J) = ZGT(K)
          IANXT(J) = J + 1
 1100   CONTINUE
        IANXT(NGVER(igs)) = IAPNT

C-----------------------------------------------------------------------
C Transform surface polygon as above.
        CALL MATPOL(NAP,IAPNT,AX,AY,AZ,IANXT,TSMAT,
     &             SBBOX,NBP,IBPNT,BX,BY,BZ,IBNXT, IERR)
        CALL CLIPFL(NBP,BX,BY,BZ,ISTAT)
        IF (ISTAT .EQ. 1) THEN
          GOTO  43
        ELSEIF (ISTAT .EQ.-1) THEN
          CALL  CUTPOL(NB,NBP,IBPNT,BX,BY,BZ,IBNXT,ISTAT)
        else
          NB=1
        ENDIF

        DO 302 JB = 1,NB

C Move to first point.
          IP = IABS(IBPNT)
          call u2pixel(BX(IP),BY(IP),iix,iiy)
          iixlast=iix   ! cast for edline call
          iiylast=iiy
          iixc=iixc+iix
          iiyc=iiyc+iiy

C Draw succeeding points until last vertex.
          IP1 = IBNXT(IP)
 451      CONTINUE
          call u2pixel(BX(IP1),BY(IP1),iix,iiy)
          if(linstyg(igs).eq.3)then
            if(mmod.eq.8)then
              call edline(iixlast,iiylast,iix,iiy,2)
            else
              call edlinewwc(iixlast,iiylast,iix,iiy,2)
            endif
          elseif(linstyg(igs).eq.2)then
            if(mmod.eq.8)then
              call edwline(iixlast,iiylast,iix,iiy)
            else
              call edwlinewwc(iixlast,iiylast,iix,iiy)
            endif
          elseif(linstyg(igs).eq.1)then
            if(mmod.eq.8)then
              call eswline(iixlast,iiylast,iix,iiy)
            else
              call eswlinewwc(iixlast,iiylast,iix,iiy)
            endif
          else
            if(mmod.eq.8)then
              call edline(iixlast,iiylast,iix,iiy,LD)
            else
              call edlinewwc(iixlast,iiylast,iix,iiy,LD)
            endif
          endif
          iixc=iixc+iix
          iiyc=iiyc+iiy

          IP1 = IBNXT(IP1)
          iixlast=iix   ! cast for call to edline
          iiylast=iiy
          IF(IP1 .NE. IP) GOTO 451

C Complete ground surface.
          call u2pixel(BX(IP),BY(IP),iix,iiy)
          if(linstyg(igs).eq.3)then
            if(mmod.eq.8)then
              call edline(iixlast,iiylast,iix,iiy,2)
            else
              call edlinewwc(iixlast,iiylast,iix,iiy,2)
            endif
          elseif(linstyg(igs).eq.2)then
            if(mmod.eq.8)then
              call edwline(iixlast,iiylast,iix,iiy)
            else
              call edwlinewwc(iixlast,iiylast,iix,iiy)
            endif
          elseif(linstyg(igs).eq.1)then
            if(mmod.eq.8)then
              call eswline(iixlast,iiylast,iix,iiy)
            else
              call eswlinewwc(iixlast,iiylast,iix,iiy)
            endif
          else
            if(mmod.eq.8)then
              call edline(iixlast,iiylast,iix,iiy,LD)
              call forceflush()
            else
              call edlinewwc(iixlast,iiylast,iix,iiy,LD)
            endif
          endif

C If ground surface names are to be displayed place near the centre of polygon.
C If the surface is not drawn (iixc=0 & iyc=0) do not label it.  For clarity drop to
C the next smaller font for the surface name and then restore font
C to normal size.
          IF(ITSNM.EQ.0)THEN
            IF(iixc.GT.0.AND.iiyc.GT.0.AND.(NAP.EQ.NBP))THEN
              write(stemp,'(A)')GSNAME(igs)
              iixc=iixc/NBP
              iiyc=iiyc/NBP
              iid2=iiyc+3
              CALL winfnt(issize)
              call CLIPST(stemp,iixc,iid2,BZ(IP1),iclp)
              if (iclp.eq.0) then
                iicol=0
                if(mmod.eq.8)then
                  call textatxy(iixc,iid2,stemp,'-',iicol)
                else
                  call textatxywwc(iixc,iid2,stemp,'-',iicol)
                endif
              endif
              CALL winfnt(IFS)
            endif
          endif
  302   CONTINUE
  43  continue
      if(mmod.eq.8) call forceflush()

      return
      end

C ******** DSGRID
C Draws a dotted grid at a given reference height (RH), grid
C density (GD), and dot density (LD) within the viewport
C with the current scaling factors etc.
C IER is returned as non-zero if there is a problem.

      SUBROUTINE DSGRID(RH,GD,LD,IER)
#include "building.h"
#include "prj3dv.h"

C Passed parameters.
      real RH,GD
      integer LD,IER
      
      COMMON/SPAD/MMOD,LIMIT,LIMTTY
      common/appcols/mdispl,nifgrey,ncset,ngset,nzonec

C Explicit local variables
      real COG1,COG2
      dimension COG1(3),COG2(3)
      DIMENSION  X(2),Y(2),Z(2)
      real GRIDD,XD,DD
      real GXMN,GYMN,XINC,YINC
      logical close,greyok

C Local variables to pass to edline.
#ifdef OSI
      integer iix,iiy,iix1,iiy1,iicol
#else
      integer*8 iix,iiy,iix1,iiy1,iicol
#endif       

      IER=0

C If libraries and monitor support greyscale then reset forground
C to 50% grey when drawing dots.
      greyok=.false.
      if(nifgrey.gt.4)then
        iicol=nifgrey-3
        if(mmod.eq.8) call winscl('i',iicol)
        greyok=.true.
      endif

      GRIDD=GD

C If GRIDD is zero compute default grid spacing.
C Use diagonal distance across the site.
      CALL ECLOSE(GRIDD,0.0,0.01,close)
      IF(close)THEN
        XD=(XMX-XMN)**2 + (YMX-YMN)**2  ! get diagonal
        DD=SQRT(XD)
        IF((DD.LE.5.))THEN
          GRIDD=0.5
        ELSEIF((DD.GT.5.).AND.(DD.LE.50.))THEN
          GRIDD=2.0
        ELSEIF((DD.GT.50.).AND.(DD.LE.100.))THEN
          GRIDD=5.0
        ELSEIF((DD.GT.100.).AND.(DD.LE.200.))THEN
          GRIDD=10.0
        ELSEIF((DD.GT.200.).AND.(DD.LE.500.))THEN
          GRIDD=25.0
        ELSEIF(DD.GT.500.)THEN
          GRIDD=50.0
        ENDIF
      ENDIF

C Make sure that the gridding will pass over the 0.0 point.
      GXMN=XMN
      IF(XMN-AINT(XMN).LT.GRIDD)GXMN=AINT(XMN)-1.

C Loop up in one direction.
      XINC=GXMN
  100 CONTINUE
      IF(XINC.LT.XMX+GRIDD)THEN

C Transform one end of line.
        COG1(1)=XINC
        COG1(2)=YMN-GRIDD
        COG1(3)=RH
        CALL VECTRN(COG1,TSMAT,COG2,IER)
        X(1)=COG2(1); Y(1)=COG2(2); Z(1)=COG2(3)

C Transform other end of line.
        COG1(1)=XINC
        COG1(2)=YMX+GRIDD
        COG1(3)=RH
        CALL VECTRN(COG1,TSMAT,COG2,IER)
        X(2)=COG2(1); Y(2)=COG2(2); Z(2)=COG2(3)

C Clip the line.
        CALL CLIPLIN(X,Y,Z,ISTAT)
        if (ISTAT.eq.1) GOTO 101
        if (ISTAT.eq.-1) then
          CALL CUTLIN(X,Y,Z,ISTAT)
          if (ISTAT.eq.-1) goto 101
        endif

C Draw the dotted line, first take calculated pixels and
C cast to local variable to pass to edline. Only do this if
C the length of the line is more than one pixel.
        call u2pixel(X(1),Y(1),iix,iiy)
        call u2pixel(X(2),Y(2),iix1,iiy1)

        if(iix.eq.iix1.and.iiy.eq.iiy1)then
          continue
        else
          if(mmod.eq.8)then
            call edline(iix,iiy,iix1,iiy1,LD)
          else
            call edlinewwc(iix,iiy,iix1,iiy1,LD)
          endif
        endif

C Move to next grid position.
  101   XINC=XINC+GRIDD
        GOTO 100
      ENDIF
      if(mmod.eq.8) call forceflush()

C Loop up in other direction.
C Make sure that the gridding will pass over the 0.0 point.
      GYMN=YMN
      IF(YMN-AINT(YMN).LT.GRIDD)GYMN=AINT(YMN)-1.
      YINC=GYMN
  110 CONTINUE
      IF(YINC.LT.YMX+GRIDD)THEN

C Transform one end of line.
        COG1(1)=XMN-GRIDD
        COG1(2)=YINC
        COG1(3)=RH
        CALL VECTRN(COG1,TSMAT,COG2,IER)
        X(1)=COG2(1); Y(1)=COG2(2); Z(1)=COG2(3)

C Transform other end of line.
        COG1(1)=XMX+GRIDD
        COG1(2)=YINC
        COG1(3)=RH
        CALL VECTRN(COG1,TSMAT,COG2,IER)
        X(2)=COG2(1); Y(2)=COG2(2); Z(2)=COG2(3)

C Clip the line.
        CALL CLIPLIN(X,Y,Z,ISTAT)
        if (ISTAT.eq.1) GOTO 111
        if (ISTAT.eq.-1) then
          CALL CUTLIN(X,Y,Z,ISTAT)
          if (ISTAT.eq.-1) GOTO 111
        endif

        call u2pixel(X(1),Y(1),iix,iiy)
        call u2pixel(X(2),Y(2),iix1,iiy1)

C Draw the dotted line after casting to local variable.
        if(iix.eq.iix1.and.iiy.eq.iiy1)then
          continue
        else
          if(mmod.eq.8)then
            call edline(iix,iiy,iix1,iiy1,LD)
          else
            call edlinewwc(iix,iiy,iix1,iiy1,LD)
          endif
        endif

C Move to next grid position.
  111   YINC=YINC+GRIDD
        GOTO 110
      ENDIF
      if(greyok)then
        iicol=0
        if(mmod.eq.8) call winscl('-',iicol)
      endif
      if(mmod.eq.8) call forceflush()

      RETURN
      END

C ******************* EMKVIEW 
C Constructs a 'viewer' format file. IUO is the file unit number.
C CFGOK is passed in as .TRUE. if the user can select zones.
C IER is returned as non-zero if there is a problem.

      SUBROUTINE EMKVIEW(IUO,CFGOK,IER)
#include "building.h"
#include "model.h"
#include "geometry.h"
#include "help.h"

C Passed parameters.
      integer IUO,IER
      logical CFGOK
      
      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      COMMON/FILEP/IFIL
      
      integer ncomp,ncon
      COMMON/C1/NCOMP,NCON
      COMMON/GB1/XB(12),YB(12),ZB(12),JVNB(6,4)
      common/grndpl/NGT,NGTV,XGT(MGTV),YGT(MGTV),ZGT(MGTV),JGVN(MGRT,8),
     &  NGVER(MGRT),IVEDGE(MGRT)
      COMMON/GT5/GSNAME(MGRT),GMLCN(MGRT)
      LOGICAL OK

C Use JGVNT for ground surfaces and JVNZ to hold inverted order for zone.
      DIMENSION IVALS(MCOM),JGVNT(MGRT,8),JVNZ(MS,MV)
      CHARACTER GMLCN*32,GSNAME*6,outs*124

      IER=0
      IUF=IFIL+1
      helpinsub='common3dv'  ! set for subroutine

C Write to previously opened viewer file.

C Present a list of zone names to pick from if the configuration file
C has been read in, otherwise ask for zones one by one.
  106 IF(CFGOK)THEN
        helptopic='zones_for_inclusion'
        call gethelptext(helpinsub,helptopic,nbhelp)
        INPIC=NCOMP
        CALL EPICKS(INPIC,IVALS,' ','Which zones to include: ',
     &    12,NCOMP,zname,' zone list',IER,nbhelp)
      ELSE
        helptopic='zones_for_inclusion'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL EASKI(INPIC,' ','Number of zone to include?',
     &       1,'W',MCOM,'W',1,'no of zones',IERI,nbhelp)
        if(ieri.eq.-3) goto 105
        IF(INPIC.EQ.0)GOTO 105

        DO 79 IZZ=1,INPIC

C Pass ICOMP in as a positive so that 'ALL' is included in list.
          ICOMP=1
          CALL EASKGEOF('Pick zones to include in view',CFGOK,ICOMP,
     &      '-',34,IER)
          IVALS(IZZ)=ICOMP
   79   CONTINUE
      ENDIF

C << shift to global variables >>
      IF(INPIC.EQ.0)GOTO 105
      DO 83 IZ=1,INPIC
        IF(IVALS(IZ).GT.0)THEN
          call georead(IUF,LGEOM(IVALS(IZ)),IVALS(IZ),1,iuout,IER)
          IF(IER.NE.0)THEN
            helptopic='zones_scan_fault'
            call gethelptext(helpinsub,helptopic,nbhelp)
            CALL EASKOK('Problem while scanning geometry.',
     &                 'Retry?',OK,nbhelp)
            if(ok)goto 106
            goto 83
          ENDIF

C Backup JVN to JVNZ
          do 121 J=1,nsur
            K=NVER(J)
            do 122 L=1,K
              JVNZ(J,L)=JVN(J,L)
  122       continue
  121     continue

C Write surface information to the viewer file. First invert edges.
          DO 120 J=1,NSUR
            K=NVER(J)
            DO 130 L=1,K/2
              ITEMP=JVNZ(J,L)
              JVNZ(J,L)=JVNZ(J,K+1-L)
              JVNZ(J,K+1-L)=ITEMP
  130       CONTINUE
  120     CONTINUE

C Write to output file.
          WRITE(IUO,18)zname(IVALS(IZ))
   18     FORMAT('GEN',1X,A12)
          WRITE(IUO,'(2I7)')NTV,NSUR
          DO 40 J=1,NTV
            WRITE(IUO,'(3F8.3)')X(J),Y(J),Z(J)
   40     CONTINUE
          DO 50 J=1,NSUR
            WRITE(IUO,'(I5,1X,99I4)')NVER(J),(JVNZ(J,K),K=1,NVER(J))
   50     CONTINUE
        ENDIF

C Check if obstructions should be included in image.
        if(IOBS(IZ).EQ.1)then
          write(outs,'(a,a)') ' Using: ',ZOBS(IZ)
          CALL USRMSG(' ',outs,'-')
        elseif(IOBS(IZ).EQ.2)THEN
          continue  ! any obstructions are in the geometry file.
        endif
        if(nbobs(iz).gt.0)then
          DO 301 IB=1,nbobs(iz)

C Conversion using cnvblk or cvnblk3a and common GB1.
            if(BLOCKTYP(iz,IB)(1:4).eq.'obs ')then
              CALL CNVBLK(XOB(iz,IB),YOB(iz,IB),ZOB(iz,IB),
     &          DXOB(iz,IB),DYOB(iz,IB),DZOB(iz,IB),BANGOB(iz,IB,1))
            elseif(BLOCKTYP(iz,IB)(1:4).eq.'obs3')then
              CALL CNVBLK3A(XOB(iz,IB),YOB(iz,IB),ZOB(iz,IB),
     &          DXOB(iz,IB),DYOB(iz,IB),DZOB(iz,IB),BANGOB(iz,IB,1),
     &          BANGOB(iz,IB,2),BANGOB(iz,IB,3))
            elseif(BLOCKTYP(iz,IB)(1:4).eq.'obsp')then
              call CNVBLKP(iz,IB) ! convert obsp type.
            endif

C Write obstruction information to the viewer file. First invert edges.
            DO 320 J=1,6
              K=4
              DO 330 L=1,K/2
                ITEMP=JVNB(J,L)
                JVNB(J,L)=JVNB(J,K+1-L)
                JVNB(J,K+1-L)=ITEMP
  330         CONTINUE
  320       CONTINUE

C Write obstructions to output file.
            WRITE(IUO,218)zname(IVALS(IZ))
  218       FORMAT('GEN',1X,A12,' obstructions')
            WRITE(IUO,'(a)')'  8  6'
            DO 240 J=1,8
              WRITE(IUO,'(3F8.3)')XB(J),YB(J),ZB(J)
  240       CONTINUE
            DO 250 J=1,6
              WRITE(IUO,'(a,99I4)') '  4 ',(JVNB(J,K),K=1,4)
  250       CONTINUE
  301     CONTINUE
        ENDIF
   83 CONTINUE

C If there is a ground description ask if to be included and then
C add at the bottom of the file.
      if(NGTV.gt.0)then
        helptopic='include_ground_in_view'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL EASKOK(' ','Include ground topography?',OK,nbhelp)
        if(OK)then

C Write surface information to the viewer file. First invert edges.
C Use temporary array for this.
          DO 142 J=1,NGT
            DO 143 L=1,NGVER(J)
              JGVNT(J,L)=JGVN(J,L)
 143        continue
 142      continue
          DO 144 J=1,NGT
            K=NGVER(J)
            DO 145 L=1,K/2
              ITEMP=JGVNT(J,L)
              JGVNT(J,L)=JGVNT(J,K+1-L)
              JGVNT(J,K+1-L)=ITEMP
  145       CONTINUE
  144     CONTINUE

C Write to output file.
          do 242,JG=1,NGT
            WRITE(IUO,'(a,1x,a)')'PLA',GSNAME(JG)
            WRITE(IUO,'(I7)')NGVER(JG)
            DO 243 J=1,NGVER(JG)
              WRITE(IUO,'(3F8.3)')XGT(JGVNT(JG,J)),YGT(JGVNT(JG,J)),
     &          ZGT(JGVNT(JG,J))
  243       continue
  242     continue
        ENDIF

        endif

C Free the viewer file.
  105 CONTINUE
      CALL ERPFREE(IUO,ISTAT)

      RETURN
      END ! of EMKVIEW 

C ******* wiresymbol ******
C Passed a 3D coordinate and draws a symbol at that
C location in a wireframe view with choice of colour.
C Notice passing of isymbol=32 or 33 and instead
C use symbol 8 with size 0 or 2.

      subroutine wiresymbol(xcord,ycord,zcord,icolour,isymbol)
#include "building.h"
#include "prj3dv.h"

C Parameters
      real xcord,ycord,zcord  ! location
      integer icolour ! colour index in the winscl zone scale
      integer isymbol ! symbol index to pass to esymbol
                      ! 8 is 2x3 block 32 is 2x2 block
                      ! 33 is 4x4 block
      common/appcols/mdispl,nifgrey,ncset,ngset,nzonec
      COMMON/SPAD/MMOD,LIMIT,LIMTTY

      real COG1,COG2
      dimension COG1(3),COG2(3)
      logical greyok

#ifdef OSI
      integer iix,iiy,iicol    ! passed to c code
#else
      integer*8 iix,iiy,iicol  ! passed to c code
#endif

      greyok=.false.      ! see if colour available
      if(nifgrey.gt.4)then
        greyok=.true.
      endif

      COG1(1)=xcord   ! setup for translation into viewing coords
      COG1(2)=ycord
      COG1(3)=zcord
      CALL VECTRN(COG1,TSMAT,COG2,IER)  ! convert to screen coords
      call CLIPPT(COG2(1),COG2(2),COG2(3),iclp)  ! check clipping
      if (iclp.ne.0) return
      call u2pixel(COG2(1),COG2(2),iix,iiy)  
      iicol=icolour                    ! set user color
      if(mmod.eq.8)then
        if(greyok)call winscl('z',iicol)
        if(isymbol.eq.32)then
          call esymbol(iix,iiy,8,0)  ! draw 2x2 dot
        elseif(isymbol.eq.33)then
          call esymbol(iix,iiy,8,2)  ! draw 4x4 dot
        else
          call esymbol(iix,iiy,isymbol,1)  ! draw the symbol
        endif
      else
        if(isymbol.eq.32)then
          call esymbolwwc(iix,iiy,8,0)  ! draw 2x2 dot
        elseif(isymbol.eq.33)then
          call esymbolwwc(iix,iiy,8,2)  ! draw 4x4 dot
        else
          call esymbolwwc(iix,iiy,isymbol,1)  ! draw the symbol
        endif
      endif
      iicol=0                          ! return to black
      if(mmod.eq.8)then
        if(greyok)call winscl('-',iicol)
        call forceflush()  ! ensure it is drawn
      endif
      return
      end

C ************* REDRAW
C Refreshes the 3D model display. Any view mode that takes notice
C of the change azimuth and elevation buttons should be implemented
C here.
C IVISMOD is the visualisation mode:
C -1 = blank
C 1 = building geometry (standard)
C 2 = CFD domain
C 3 = 3-D gridding (grd module)

      subroutine redraw(IER)

      common/SPAD/MMOD,LIMIT,LIMTTY
      COMMON/C6/INDCFG
      COMMON/MODVIS/IVISMOD
      logical focussname ! if true only label highlighted surfaces.

C Check we're in graphic mode, just in case.
      if (MMOD.ne.8) return

C Check we have a model to display.
      if ((IVISMOD.eq.1.or.IVISMOD.eq.3).and.
     &  (indcfg.ne.1.and.indcfg.ne.3)) goto 999

      if (IVISMOD.eq.-1) then
        call startbuffer()
      elseif (IVISMOD.eq.1) then
        focussname=.false.
        CALL cadjview(focussname,IER)
      elseif (IVISMOD.eq.2) then
        CALL CFDVIEW(IER)
      elseif (IVISMOD.eq.3) then
        i1=-1
        i2=-1
        CALL GRAAPH(i1,i2)
      endif

 999  return
      end


C ************* TMPMENU
C This will display a small menu with an option to continue. Otherwise
C it behaves like any other menu subroutine. The intention is to call
C this during mouse-related operations, to allow users to reposition the
C view without interfering with the operation in progress. Once this
C subroutine returns, the operation can resume.

      SUBROUTINE TMPMENU

      character*10 cont(1)

      cont(1)='- continue'
      i=0
      do while (i.eq.0)
        call emenu('paused ...',cont,1,i)
      enddo

      return
      end
