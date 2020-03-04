C mfprb3.F includes subroutines which either support the drawing of
C 3D flow entities or user dialogs establishing postion.

C  mfwfdraw - draws 3D flow network overlay on wireframe from commons passed colour.
C  symbol_flow_nd - draws a flow node symbol passed name type colour & position
C  symbol_flow_cmp - draws a flow component symbol passed name type colour & position
C  draw_flow_link - draws a flow connection between two coordinates
C  draw_connect_assoc - draws the nodes & components associated
C                       with flow connection IFCN.
C  draw_surf_flow_symbol - places flow symbol after user chooses surface USE 

C mfcompposition - user interactions eastablishing flow component position.
C mfnodeposition - user interactions for flow node position.
C scan_flow_use - uses create_zone_nodes and create_surface_node_cmp to
C   fills flow common blocks from zone & surface attributes and creates
C   both legacy and a 3D flow files while drawing wireframe overlay.
C create_surface_node_cmp - scans surface USE attributes related to mass flows,
C   extracts position information, fills flow common blocks.
C isuseflowrelated - tests whether surface attributes are flow related.

C ************************* MFWFDRAW
C MFWFDRAW Fluid flow file: draw common block contents.
C colour is either 'r' red, 'b' blue or 'g' grey50 for components.
      SUBROUTINE MFWFDRAW(col)

#include "building.h"
#include "net_flow.h"
#include "net_flow_data.h"

C      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      COMMON/SPAD/MMOD,LIMIT,LIMTTY
      character col*1
      CHARACTER name*12,twochar*2
      integer POS,NEG,CMPIS
      real xc,yc,zc,xw,yw,zw
      logical inview,inv2,inv3,iswater
      character altcol*1,ccol*1

      if(MMOD.lt.8) return  ! in text mode do not draw
      altcol='b'  ! blue alternative colour

C Node data.
      IF(NNOD.gt.0)then
        DO INOD=1,NNOD
          if(NNOD.gt.50)then
            write(name,'(a,i3.3)') 'n:',inod
          else
            name=NDNAM(INOD)
          endif
          xw=HNOD(INOD,1);yw=HNOD(INOD,2); zw=HNOD(INOD,3)

C Check within viewing angle.
          call withinview(xw,yw,zw,inview,ier)
          if(inview)then
            if(NDTYP(INOD).eq.0)then
              if(NDFLD(INOD).eq.1)then
                call symbol_flow_nd(name,'iu',col,xw,yw,zw)
              else
                call symbol_flow_nd(name,'iu',altcol,xw,yw,zw)
              endif
              call pausems(30)
            elseif(NDTYP(INOD).eq.1)then
              if(NDFLD(INOD).eq.1)then
                call symbol_flow_nd(name,'ik',col,xw,yw,zw)
              else
                call symbol_flow_nd(name,'ik',altcol,xw,yw,zw)
              endif
              call pausems(30)
            elseif(NDTYP(INOD).eq.2)then
              call symbol_flow_nd(name,'bk',col,xw,yw,zw)
              call pausems(30)
            elseif(NDTYP(INOD).eq.3)then
              call symbol_flow_nd(name,'bw',col,xw,yw,zw)
              call pausems(30)
            endif
          endif
        ENDDO
      endif

C Component data.
      IF(NCMP.gt.0)then
        DO ICMP=1,NCMP
          if(NCMP.gt.50)then
            write(name,'(a,i3.3)') 'c:',icmp
          else
            name=CMNAM(ICMP)
          endif
          call get_component_code(ICMP,twochar)
          
          xc=HCMP(ICMP,1,1);yc=HCMP(ICMP,1,2); zc=HCMP(ICMP,1,3)
          call withinview(xc,yc,zc,inview,ier)
          if(inview)then
            call ECLOSE(SUPCMP(ICMP,1),2.0,0.0001,iswater)
            IF(ITPCMP(ICMP).EQ. 10) THEN
              if(.NOT.iswater)then
                call symbol_flow_cmp(name,'bx',col,xc,yc,zc) ! power law flow (box)
              else
                call symbol_flow_cmp(name,'bx',altcol,xc,yc,zc) ! power law flow (box)
              endif
              call pausems(30)
            ELSE IF(ITPCMP(ICMP).EQ. 11) THEN
              continue  ! Self regulating vent
            ELSE IF(ITPCMP(ICMP).EQ. 12) THEN
              if(.NOT.iswater)then
                call symbol_flow_cmp(name,'bx',col,xc,yc,zc) ! power law vol flow (box)
              else
                call symbol_flow_cmp(name,'bx',altcol,xc,yc,zc) ! power law vol flow (box)
              endif
              call pausems(30)
            ELSE IF(ITPCMP(ICMP).EQ. 15) THEN
              if(.NOT.iswater)then
                call symbol_flow_cmp(name,'bx',col,xc,yc,zc) ! power law mass flow (box)
              else
                call symbol_flow_cmp(name,'bx',altcol,xc,yc,zc) ! power law mass flow (box)
              endif
              call pausems(30)
            ELSE IF(ITPCMP(ICMP).EQ. 17) THEN
              if(.NOT.iswater)then
                call symbol_flow_cmp(name,'bx',col,xc,yc,zc) ! power law mass flow (box)
              else
                call symbol_flow_cmp(name,'bx',altcol,xc,yc,zc) ! power law mass flow (box)
              endif
              call pausems(30)
            ELSE IF(ITPCMP(ICMP).EQ. 20) THEN
              if(.NOT.iswater)then
                call symbol_flow_cmp(name,'bx',col,xc,yc,zc) ! quadratic law vol flow (box)
              else
                call symbol_flow_cmp(name,'bx',altcol,xc,yc,zc) ! quadratic law vol flow (box)
              endif
              call pausems(30)
            ELSE IF(ITPCMP(ICMP).EQ. 25) THEN
              if(.NOT.iswater)then
                call symbol_flow_cmp(name,'bx',col,xc,yc,zc) ! quadratic law mass flow (box)
              else
                call symbol_flow_cmp(name,'bx',altcol,xc,yc,zc) ! quadratic law mass flow (box)
              endif
              call pausems(30)
            ELSE IF(ITPCMP(ICMP).EQ. 30) THEN
              if(.NOT.iswater)then
                call symbol_flow_cmp(name,'bx',col,xc,yc,zc) ! constant vol flow (box)
              else
                call symbol_flow_cmp(name,'bx',altcol,xc,yc,zc) ! constant vol flow (box)
              endif
              call pausems(30)
            ELSE IF(ITPCMP(ICMP).EQ. 35) THEN
              if(.NOT.iswater)then
                call symbol_flow_cmp(name,'bx',col,xc,yc,zc) ! constant mass flow (box)
              else
                call symbol_flow_cmp(name,'bx',altcol,xc,yc,zc) ! constant mass flow (box)
              endif
              call pausems(30)
            ELSE IF(ITPCMP(ICMP).EQ. 40) THEN
              if(.NOT.iswater)then
                call symbol_flow_cmp(name,'or',col,xc,yc,zc) ! Common orifice flow
              else
                call symbol_flow_cmp(name,'or',altcol,xc,yc,zc) ! Common orifice flow
              endif
              call pausems(30)
            ELSE IF(ITPCMP(ICMP).EQ. 50) THEN
              if(.NOT.iswater)then
                call symbol_flow_cmp(name,'cd',col,xc,yc,zc) ! Laminar pipe vol
              else
                call symbol_flow_cmp(name,'cd',altcol,xc,yc,zc) ! Laminar pipe vol
              endif
              call pausems(30)
            ELSE IF(ITPCMP(ICMP).EQ.110) THEN
              call symbol_flow_cmp(name,'or',col,xc,yc,zc) ! Specific air flow
              call pausems(30)
            ELSE IF(ITPCMP(ICMP).EQ.120) THEN
              call symbol_flow_cmp(name,'cr',col,xc,yc,zc) ! Specific air flow crack
              call pausems(30)
            ELSE IF(ITPCMP(ICMP).EQ.130) THEN
              call symbol_flow_cmp(name,'bi',col,xc,yc,zc) ! Specific air flow door
              call pausems(30)
            ELSE IF(ITPCMP(ICMP).EQ.210) THEN
              if(.NOT.iswater)then
                call symbol_flow_cmp(name,'cd',col,xc,yc,zc) ! General flow conduit
              else
                call symbol_flow_cmp(name,'cd',altcol,xc,yc,zc) ! General flow conduit
              endif
              call pausems(30)
            ELSE IF(ITPCMP(ICMP).EQ.211) THEN
              call symbol_flow_cmp(name,'bx',col,xc,yc,zc) ! Cowls and roof outlets (box)
              call pausems(30)
            ELSE IF(ITPCMP(ICMP).EQ.220) THEN
              call symbol_flow_cmp(name,'cd',col,xc,yc,zc) ! Conduit ending in conv 3 leg
              call pausems(30)
            ELSE IF(ITPCMP(ICMP).EQ.230) THEN
              call symbol_flow_cmp(name,'cd',col,xc,yc,zc) ! Conduit starts in diverg 3 leg
              call pausems(30)
            ELSE IF(ITPCMP(ICMP).EQ.240) THEN
              call symbol_flow_cmp(name,'cd',col,xc,yc,zc) ! Conduit ending in conv 4 leg
              call pausems(30)
            ELSE IF(ITPCMP(ICMP).EQ.250) THEN
              call symbol_flow_cmp(name,'cd',col,xc,yc,zc) ! Conduit starts in diverg 4 leg
              call pausems(30)
            ELSE IF(ITPCMP(ICMP).EQ.310) THEN
              if(.NOT.iswater)then
                call symbol_flow_cmp(name,'bx',col,xc,yc,zc) ! General flow inducer (box)
              else
                call symbol_flow_cmp(name,'bx',altcol,xc,yc,zc) ! General flow inducer (box)
              endif
              call pausems(30)
            ELSE IF(ITPCMP(ICMP).EQ.410) THEN
              if(.NOT.iswater)then
                call symbol_flow_cmp(name,'bx',col,xc,yc,zc) ! General flow corrector (box)
              else
                call symbol_flow_cmp(name,'bx',altcol,xc,yc,zc) ! General flow corrector (box)
              endif
              call pausems(30)
            ELSE IF(ITPCMP(ICMP).EQ.420) THEN
              if(.NOT.iswater)then
                call symbol_flow_cmp(name,'bx',col,xc,yc,zc) ! Corrector polynomial (box)
              else
                call symbol_flow_cmp(name,'bx',altcol,xc,yc,zc) ! Corrector polynomial (box)
              endif
              call pausems(30)
            ELSE IF(ITPCMP(ICMP).EQ.460) THEN
              call symbol_flow_cmp(name,'bx',col,xc,yc,zc) ! Fixed flow rates controller (box)
              call pausems(30)
            ELSE IF(ITPCMP(ICMP).EQ.500) THEN
              call symbol_flow_cmp(name,'bx',col,xc,yc,zc) ! multi-configuration (box)
              call pausems(30)
            endif
          endif
        ENDDO
      endif

C Connections drawn as dotted lines from node -> component -> node.
      IF(NCNN.gt.0)then
        DO ICNN=1,NCNN
          POS=NODPS(ICNN); NEG=NODNE(ICNN)
          if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
            ccol='r'
          else
            ccol='b'
          endif
          CMPIS=ITPCON(ICNN)
          xc=HCMP(CMPIS,1,1); yc=HCMP(CMPIS,1,2); zc=HCMP(CMPIS,1,3)

C Check if all 3 points are within view.
          call withinview(xc,yc,zc,inview,ier)
          call withinview(HNOD(POS,1),HNOD(POS,2),HNOD(POS,3),inv2,ier)
          if(inview.and.inv2)then
            call draw_flow_link(HNOD(POS,1),HNOD(POS,2),HNOD(POS,3),
     &        xc,yc,zc,ccol)
          endif
          call withinview(HNOD(NEG,1),HNOD(NEG,2),HNOD(NEG,3),inv3,ier)
          if(inview.and.inv3)then
            call draw_flow_link(xc,yc,zc,
     &        HNOD(NEG,1),HNOD(NEG,2),HNOD(NEG,3),ccol)
          endif
        ENDDO
      endif

  100 RETURN

      END

C ********* symbol_flow_nd
C Draw a symbol for a flow node with optional name and type.
C type='iu' internal unknown, 'ik' internal known
C type='bw' boundary wind, 'bk' boundary known

      subroutine symbol_flow_nd(name,type,col,xw,yw,zw)
#include "building.h"
#include "prj3dv.h"
      COMMON/SPAD/MMOD,LIMIT,LIMTTY
      CHARACTER name*12,type*2,col*1
      real xw,yw,zw

      integer ifs,itfs,imfs
      COMMON/GFONT/IFS,ITFS,IMFS
      DIMENSION  COG1(3),COG2(3),COG3(3)
      character t13*13
      integer izsize,issize
#ifdef OSI
      integer iix,iiy,iid4,iicol
#else
      integer*8 iix,iiy,iid4,iicol
#endif

      if(MMOD.lt.8) return  ! in text mode do not draw
      issize=4; izsize=5  ! for surfaces and zone name

C Get pixel position to draw at.
      COG1(1)=xw; COG1(2)=yw; COG1(3)=zw
      CALL VECTRN(COG1,TEMAT,COG2,IER)
      CALL VECTRN(COG2,ETSMAT,COG3,IER)
      call u2pixel(COG3(1),COG3(2),iix,iiy)

C Draw circles depending on type. (powderblue is 27) (navy blue is 34)
      if(col.eq.'r') iicol=0   ! red
      if(col.eq.'b') iicol=13  ! royal blue
      if(col.eq.'g') iicol=32  ! grey50
      call winscl('z',iicol)
      if(type.eq.'bw')then
        CALL ecirc(iix,iiy,5,0)
        CALL ecirc(iix,iiy,3,0)
      elseif(type.eq.'bk')then
        CALL ecirc(iix,iiy,5,0)
        CALL ecirc(iix,iiy,4,0)
        CALL ecirc(iix,iiy,3,0)
      elseif(type.eq.'iu')then
        CALL ecirc(iix,iiy+8,4,1)
      elseif(type.eq.'ik')then
        CALL ecirc(iix,iiy,5,0)
        CALL ecirc(iix,iiy,3,1)
      endif

C Write name for other than 'iu'.
      if(type.eq.'iu')then
        write(t13,'(A)')name(1:lnblnk(name))
        CALL winfnt(issize)
        iid4 = iix +8  ! offset to right of circle
        if(col.eq.'r') iicol=0   ! red
        if(col.eq.'g') iicol=32  ! grey50
        call textatxy(iid4,iiy+8,t13,'z',iicol)
        CALL winfnt(IFS)
        call forceflush()
      else
        write(t13,'(A)')name(1:lnblnk(name))
        CALL winfnt(issize)
        iid4 = iix +8  ! offset to right of circle
        if(col.eq.'r') iicol=0   ! red
        if(col.eq.'g') iicol=32  ! grey50
        call textatxy(iid4,iiy+5,t13,'z',iicol)
        CALL winfnt(IFS)
        call forceflush()
      endif
      iicol=0
      call winscl('-',iicol)

      return
      end

C ********* symbol_flow_cmp
C Draw a symbol for a flow component with name and symbol type.
C type='bx' box, type='cd' conduit type='cr' crack, 'or' orifice,
C type='bi' bi-directional
C col is 'r' red  'b' blue  'g' grey
      subroutine symbol_flow_cmp(name,type,col,xw,yw,zw)
#include "building.h"
#include "prj3dv.h"
      COMMON/SPAD/MMOD,LIMIT,LIMTTY
      CHARACTER name*12,type*2,col*1
      real xw,yw,zw

      integer ifs,itfs,imfs
      COMMON/GFONT/IFS,ITFS,IMFS
      DIMENSION  COG1(3),COG2(3),COG3(3)
      character t13*13
      integer izsize,issize
#ifdef OSI
      integer iix,iiy,iid4,iicol
#else
      integer*8 iix,iiy,iid4,iicol
#endif

      if(MMOD.lt.8) return  ! in text mode do not draw
      issize=4; izsize=5  ! for surfaces and zone name

C Get pixel position to draw at.
      COG1(1)=xw; COG1(2)=yw; COG1(3)=zw
      CALL VECTRN(COG1,TEMAT,COG2,IER)
      CALL VECTRN(COG2,ETSMAT,COG3,IER)
      call u2pixel(COG3(1),COG3(2),iix,iiy)

C Draw symbol depending on type.
      if(col.eq.'r') iicol=0   ! red
      if(col.eq.'b') iicol=13  ! royal blue
      if(col.eq.'g') iicol=32  ! grey50
      call winscl('z',iicol)
      if(type.eq.'bx')then
        call esymbol(iix,iiy,38,1) ! box symbol
      elseif(type.eq.'cr')then
        call esymbol(iix,iiy,34,1) ! crack symbol
      elseif(type.eq.'cd')then
        call esymbol(iix,iiy,37,1) ! conduit symbol
      elseif(type.eq.'or')then
        call esymbol(iix,iiy,36,1) ! orifice symbol
      elseif(type.eq.'bi')then
        call esymbol(iix,iiy+5,35,1) ! bidir symbol
      endif

C Write name to the right of the symbol.
      t13=' '
      if(type.eq.'bx')then
        write(t13,'(A)')name(1:lnblnk(name))
        CALL winfnt(issize)
        iid4 = iix +8
        if(col.eq.'r') iicol=0   ! red
        if(col.eq.'b') iicol=13  ! royal blue
        if(col.eq.'g') iicol=32  ! grey50
        call textatxy(iid4,iiy+5,t13,'z',iicol)
        CALL winfnt(IFS)
        call forceflush()
      elseif(type.eq.'cr')then
        write(t13,'(A)')name(1:lnblnk(name))
        CALL winfnt(issize)
        iid4 = iix +10
        if(col.eq.'r') iicol=0   ! red
        if(col.eq.'b') iicol=13  ! royal blue
        if(col.eq.'g') iicol=32  ! grey50
        call textatxy(iid4,iiy+5,t13,'z',iicol)
        CALL winfnt(IFS)
        call forceflush()
      elseif(type.eq.'cd')then
        write(t13,'(A)')name(1:lnblnk(name))
        CALL winfnt(issize)
        iid4 = iix +10
        if(col.eq.'r') iicol=0   ! red
        if(col.eq.'b') iicol=13  ! royal blue
        if(col.eq.'g') iicol=32  ! grey50
        call textatxy(iid4,iiy+5,t13,'z',iicol)
        CALL winfnt(IFS)
        call forceflush()
      elseif(type.eq.'or')then
        write(t13,'(A)')name(1:lnblnk(name))
        CALL winfnt(issize)
        iid4 = iix +8
        if(col.eq.'r') iicol=0   ! red
        if(col.eq.'b') iicol=13  ! royal blue
        if(col.eq.'g') iicol=32  ! grey50
        call textatxy(iid4,iiy+5,t13,'z',iicol)
        CALL winfnt(IFS)
        call forceflush()
      elseif(type.eq.'bi')then
        write(t13,'(A)')name(1:lnblnk(name))
        CALL winfnt(issize)
        iid4 = iix +8
        if(col.eq.'r') iicol=0   ! red
        if(col.eq.'b') iicol=13  ! royal blue
        if(col.eq.'g') iicol=32  ! grey50
        call textatxy(iid4,iiy+5,t13,'z',iicol)
        CALL winfnt(IFS)
        call forceflush()
      endif
      iicol=0
      call winscl('-',iicol)

      return
      end

C ********* draw_flow_link
C Draw_flow_link places dotted line between a flow node and
C component in such as way as to not obscure node and/or component.

      subroutine draw_flow_link(xs,ys,zs,xe,ye,ze,col)
#include "building.h"
#include "prj3dv.h"
      COMMON/SPAD/MMOD,LIMIT,LIMTTY
C      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      real xs,ys,zs,xe,ye,ze
      character col*1
      real tdis
C      character outs*124
      DIMENSION COG1(3),COG2(3),COG3(3)
#ifdef OSI
      integer iix1,iiy1,iix2,iiy2,iicol
#else
      integer*8 iix1,iiy1,iix2,iiy2,iicol
#endif

      if(MMOD.lt.8) return  ! in text mode do not draw

C Find distance and report.
      tdis= crowxyz(xs,ys,zs,xe,ye,ze)
C      write(outs,'(a,3f7.3,a,3f7.3,a,f9.4)') ' Distance between ',
C     &   xs,ys,zs,' & ',xe,ye,ze,' =',tdis
C      call edisp(iuout,outs)
      if(tdis.gt.10)then
        vdis=0.3
      elseif(tdis.gt.2.and.tdis.le.10)then
        vdis=0.2
      else
        vdis=0.1
      endif

C Use ratio calculation.
      r2 = tdis - vdis
      r1 = vdis
      x3 = ((r2 * xs) + (r1 * xe))/tdis
      y3 = ((r2 * ys) + (r1 * ye))/tdis
      z3 = ((r2 * zs) + (r1 * ze))/tdis
      COG1(1)=x3; COG1(2)=y3; COG1(3)=z3
      CALL VECTRN(COG1,TEMAT,COG2,IER)
      CALL VECTRN(COG2,ETSMAT,COG3,IER)
      call u2pixel(COG3(1),COG3(2),iix1,iiy1)

C Now do the other end.
      r2 = tdis - vdis
      r1 = vdis
      x3 = ((r2 * xe) + (r1 * xs))/tdis
      y3 = ((r2 * ye) + (r1 * ys))/tdis
      z3 = ((r2 * ze) + (r1 * zs))/tdis
      COG1(1)=x3; COG1(2)=y3; COG1(3)=z3
      CALL VECTRN(COG1,TEMAT,COG2,IER)
      CALL VECTRN(COG2,ETSMAT,COG3,IER)
      call u2pixel(COG3(1),COG3(2),iix2,iiy2)

C Draw a dotted line (in red=0) (powderblue 27) (navy blue 34) (midnight blue 1).
      if(col.eq.'r') iicol=0   ! red
      if(col.eq.'b') iicol=13  ! royal blue
      if(col.eq.'g') iicol=32  ! grey50
C      iicol=0  ! red
      call winscl('z',iicol)
      call edline(iix1,iiy1,iix2,iiy2,2)
      call forceflush()
      iicol=0
      call winscl('-',iicol)
      return
      end

C Based on the surface USE attributes draw relevant symbol
C on the wire-frame. If act 'v' show names.
      subroutine draw_surf_flow_symbol(loop,act)
#include "building.h"
#include "geometry.h"
#include "net_flow.h"
#include "net_flow_data.h"
#include "esprdbfile.h"
C esprdbfile.h supplies the following:
C LPRFDB,IPRODB (for event profile database)
#include "prj3dv.h"

C Passed parameters.
      integer loop     ! the connection number for the surface
      character act*1  ! '-' try not to show component name

      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      COMMON/C3/IC1(MCON),IE1(MCON),ICT(MCON),IC2(MCON),IE2(MCON)
      COMMON/AFN/IAIRN,LAPROB,ICAAS(MCOM)
      INTEGER :: iairn,icaas
      character LAPROB*72

      CHARACTER outs*144,name*12
      CHARACTER SN*12,ZN*12,SNO*12,ZNO*12
      CHARACTER zsn*28
      logical dupedges
      logical inthenetwork
      logical isrelated ! true if surface attributes are flow related
      integer loop2
      integer llpos,lrpos,ulpos,urpos     ! closest to BB corners for parent
      real edgelower,edgeupper  ! Z of current surface lower & upper edges.
      real xc,yc,zc

C Cases to ignore.
      iz=IC1(loop); is=IE1(loop)     ! get current zone & surface.
      write(SN,'(a)') SNAME(iz,is)
      write(ZN,'(a)') ZNAME(iz)
      lnsn=lnblnk(SN)
      lnzn=lnblnk(ZN)
      write(zsn,'(3a)') ZN(1:lnzn),':',SN(1:lnsn)
      lnzsn=lnblnk(zsn)  ! remember length of combined zone:surface name

      if(ICT(loop).eq.3)then
        izo=IC2(loop); iso=IE2(loop) ! other side zone & surface.
        write(SNO,'(a)') SNAME(izo,iso)
        write(ZNO,'(a)') ZNAME(izo)
        lnsno=lnblnk(SNO)
        lnzno=lnblnk(ZNO)
      else
        izo=0; iso=0  ! no other side to check.
      endif
      call isuseflowrelated(iz,is,isrelated)
      if(.NOT.isrelated)then
        return
      endif

C Under what circumstances do we not want to draw symbols?  << NOT YET DECIDED >>
C      if(IAIRN.eq.3) return

      PI = 4.0 * ATAN(1.0)

      name='  '
      inthenetwork=.false.
      xc=0.0; yc=0.0; zc=0.0
      if(IAIRN.eq.3)then

C If IAIRN = 3 loop through the components looking for one
C which references the same zone and surface (or the other
C side if a partition) and use the name and position
C in the flow network. Otherwise use the approximate location
C and truncated name. Jump out when inthenetwork becomes true.
        do loop2=1,NCMP
          lnczn=lnblnk(CMPASSOC(loop2,1))
          lncsn=lnblnk(CMPASSOC(loop2,2))
          if(ZN(1:lnzn).eq.CMPASSOC(loop2,1)(1:lnczn).and.
     &       SN(1:lnsn).eq.CMPASSOC(loop2,2)(1:lncsn))then
            write(name,'(a)') CMNAM(loop2)(1:lnblnk(CMNAM(loop2)))
            xc=HCMP(loop2,1,1); yc=HCMP(loop2,1,2) ! recover
            zc=HCMP(loop2,1,3)
            inthenetwork=.true.
            EXIT
          endif
          if(ICT(loop).eq.3)then  ! Check match with other surface.
            if(ZNO(1:lnzno).eq.CMPASSOC(loop2,1)(1:lnczn).and.
     &         SNO(1:lnsno).eq.CMPASSOC(loop2,2)(1:lncsn))then
              write(name,'(a)') CMNAM(loop2)(1:lnblnk(CMNAM(loop2)))
              xc=HCMP(loop2,1,1); yc=HCMP(loop2,1,2) ! recover
              zc=HCMP(loop2,1,3)
              inthenetwork=.true.
              EXIT
            endif
          endif
        enddo  ! of loop2
      endif

C Note: llpos, lrpos etc are the position in the iszjvn list at each
C of the corners of the surface.
      CALL ZSURLEHI(IZ,IS,XYMAX,ZMAX,llpos,lrpos,ulpos,urpos)
      if(llpos.gt.0.and.ulpos.gt.0)then
        edgelower=szcoords(iz,iszjvn(iz,is,llpos),3)
        edgeupper=szcoords(iz,iszjvn(iz,is,ulpos),3)
      endif

C For DOOR variants.
      if(SUSE(iz,is,1)(1:4).eq.'DOOR'.or.
     &   SUSE(iz,is,1)(1:6).eq.'P-DOOR'.or.
     &   SUSE(iz,is,1)(1:6).eq.'H-DOOR'.or.
     &   SUSE(iz,is,1)(1:6).eq.'H-DOOR')then
        if(SUSE(iz,is,2)(1:6).eq.'CLOSED')then

C Door crack - crack component, assume width = 2mm,
C length = perimeter, on centre line.
          call zsurfprm(iz,is,dupedges,perim)
          if(inthenetwork)then
            continue
          else
            write(name,'(a)') 'DoCr'
            xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
            zc=SURCOG(iz,is,3)
          endif
          write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),' ',name,
     &      ' door crack width 2mm len',perim,'m COG @XYZ',
     &      xc,yc,zc
          if(act.eq.'v') call edisp(iuout,outs)
          if(inthenetwork)then
            call symbol_flow_cmp(name,'cr','r',xc,yc,zc)
          else
            call symbol_flow_cmp(name,'cr','g',xc,yc,zc)
          endif
          call pausems(30)
        elseif(SUSE(iz,is,2)(1:8).eq.'UNDERCUT')then

C Undercut at door - common orifice at lowest horizontal edge.
          if(inthenetwork)then
            continue
          else
            write(name,'(a)') 'DoUc'
            xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
            zc=edgelower
          endif
          write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),' ',name,
     &    ' door undercut width ',XYMAX,'m @XYZ',xc,yc,zc
          if(act.eq.'v') call edisp(iuout,outs)
          if(inthenetwork)then
            call symbol_flow_cmp(name,'or','r',xc,yc,zc)
          else
            call symbol_flow_cmp(name,'or','g',xc,yc,zc)
          endif
          call pausems(30)
        elseif(SUSE(iz,is,2)(1:4).eq.'OPEN')then

C Open door - common orifice component, (type 40) IVALCM(10) MF040I
          if(inthenetwork)then
            continue
          else
            write(name,'(a)') 'DoOP'
            xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
            zc=SURCOG(iz,is,3)
          endif
          write(outs,'(4a,f6.1,a,f6.1,a,f7.3,a,3f6.2)') 
     &      zsn(1:lnzsn),' ',name,
     &      ' door open within surface of width ',XYMAX,' height',ZMAX,
     &      ' area',SNA(iz,is),'m^2 COG @XYZ',xc,yc,zc
          if(act.eq.'v') call edisp(iuout,outs)
          if(inthenetwork)then
            call symbol_flow_cmp(name,'or','r',xc,yc,zc)
          else
            call symbol_flow_cmp(name,'or','g',xc,yc,zc)
          endif
          call pausems(30)
        elseif(SUSE(iz,is,2)(1:5).eq.'BIDIR'.or.
     &         SUSE(iz,is,2)(1:10).eq.'AJAR-BIDIR')then

C Bidirectional door - door component,  IVALCM(14)=130 MF130I
          if(inthenetwork)then
            continue
          else
            if(SUSE(iz,is,2)(1:5).eq.'BIDIR')then
              write(name,'(a)') 'DoBi'
            elseif(SUSE(iz,is,2)(1:10).eq.'AJAR-BIDIR')then
              write(name,'(a)') 'DABi'
            endif
            xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
            zc=SURCOG(iz,is,3)
          endif
          VAL3=ABS(ZCOG(iz,3)-edgelower)  ! base to zone COG
          write(outs,'(4a,f6.1,a,f6.1,a,3f7.3,a,f6.3)') 
     &      zsn(1:lnzsn),' ',name,
     &      ' door bidir in surface w= ',XYMAX,
     &      ' h=',ZMAX,' coef 0.5 base @XYZ',xc,yc,zc,
     &      ' delta to zone COG',VAL3
          if(act.eq.'v') call edisp(iuout,outs)
          if(inthenetwork)then
            call symbol_flow_cmp(name,'bi','r',xc,yc,zc)
          else
            call symbol_flow_cmp(name,'bi','g',xc,yc,zc)
          endif
          call pausems(30)
        endif
      elseif(SUSE(iz,is,1)(1:5).eq.'FRAME'.or.
     &       SUSE(iz,is,1)(1:7).eq.'F-FRAME')then

C For FRAME vairants.
        if(SUSE(iz,is,2)(1:6).eq.'CLOSED'.or.
     &     SUSE(iz,is,2)(1:5).eq.'CRACK')then

C Frame crack - crack component, assume width = 1mm,
C length = frame inside perimeter, on centre line.
          call zsurfprm(iz,is,dupedges,perim)
          if(inthenetwork)then
            continue
          else
            write(name,'(a)') 'FrCr'
            xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
            zc=SURCOG(iz,is,3)
          endif
          write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),' ',name,
     &      ' frame crack width 1mm len',perim,'m COG @XYZ',
     &      xc,yc,zc
          if(act.eq.'v') call edisp(iuout,outs)
          if(inthenetwork)then
            call symbol_flow_cmp(name,'cr','r',xc,yc,zc)
          else
            call symbol_flow_cmp(name,'cr','g',xc,yc,zc)
          endif
          call pausems(30)
        elseif(SUSE(iz,is,2)(1:4).eq.'VENT')then

C Frame trickle vent - air flow opening, assume width = 0.25* frame inside width,
C height = 15mm, at highest horizontal edge. OK
          if(inthenetwork)then
            continue
          else
            write(name,'(a)') 'FrVe'
            xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
            zc=edgeupper
          endif
          trwidth=XYMAX*0.6
          if(trwidth.gt.0.7) trwidth=0.7
          if(trwidth.lt.0.4.and.XYMAX.lt.0.4) trwidth=XYMAX
          write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),'  ',name,
     &      ' frame trickle vent width ',trwidth,' height 15mm @XYZ',
     &      xc,yc,zc
          if(act.eq.'v') call edisp(iuout,outs)
          if(inthenetwork)then
            call symbol_flow_cmp(name,'or','r',xc,yc,zc)
          else
            call symbol_flow_cmp(name,'or','g',xc,yc,zc)
          endif
          call pausems(30)
        endif
      elseif(SUSE(iz,is,1)(1:5).eq.'GRILL')then

C For GRILLs.
        if(SUSE(iz,is,2)(1:6).eq.'CLOSED'.or.
     &     SUSE(iz,is,2)(1:5).eq.'CRACK')then

C Grill crack - crack component, assume width = 1mm,
C length = frame perimeter, on centre line.
          call zsurfprm(iz,is,dupedges,perim)
          if(inthenetwork)then
            continue
          else
            write(name,'(a)') 'GrCR'
            xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
            zc=SURCOG(iz,is,3)
          endif
          write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),' ',name,
     &      ' grill crack width 1mm length',perim,'m COG @XYZ',
     &      xc,yc,zc
          if(act.eq.'v') call edisp(iuout,outs)
          if(inthenetwork)then
            call symbol_flow_cmp(name,'cr','r',xc,yc,zc)
          else
            call symbol_flow_cmp(name,'cr','g',xc,yc,zc)
          endif
          call pausems(30)
        elseif(SUSE(iz,is,2)(1:5).eq.'INLET')then

C Grill inlet - a conduit (type 210) IVALCM(15), as in MF210I
          call zsurfprm(iz,is,dupedges,perim)
          if(inthenetwork)then
            continue
          else
            write(name,'(a)') 'GrIN'
            xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
            zc=SURCOG(iz,is,3)
          endif
          write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),' ',name,
     &      ' grill inlet hydr diam',
     &      4.0*SNA(iz,is)/perim,'m COG @XYZ',xc,yc,zc
          if(act.eq.'v') call edisp(iuout,outs)
          if(inthenetwork)then
            call symbol_flow_cmp(name,'cd','r',xc,yc,zc)
          else
            call symbol_flow_cmp(name,'cd','g',xc,yc,zc)
          endif
          call pausems(30)
        elseif(SUSE(iz,is,2)(1:7).eq.'EXTRACT')then

C Grill extract - constant vol flow comp. (type 30) IVALCM(8), as in MF030I
          call zsurfprm(iz,is,dupedges,perim)
          if(inthenetwork)then
            continue
          else
            write(name,'(a)') 'GrEX'
            xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
            zc=SURCOG(iz,is,3)
          endif
          write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),' ',name,
     &      ' grill extract hydr diam',
     &      4.0*SNA(iz,is)/perim,'m COG @XYZ',xc,yc,zc
          if(act.eq.'v') call edisp(iuout,outs)
          if(inthenetwork)then
            call symbol_flow_cmp(name,'bx','r',xc,yc,zc)
          else
            call symbol_flow_cmp(name,'bx','g',xc,yc,zc)
          endif
          call pausems(30)
        elseif(SUSE(iz,is,2)(1:4).eq.'OPEN')then

C Grill opening - common orifice component, (type 40) IVALCM(10) MF040I
          call zsurfprm(iz,is,dupedges,perim)
          if(inthenetwork)then
            continue
          else
            write(name,'(a)') 'GrOP'
            xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
            zc=SURCOG(iz,is,3)
          endif
          write(outs,'(4a,f6.1,a,3f6.2)') zsn(1:lnzsn),' ',name,
     &      ' grill orifice area',SNA(iz,is),'m^2 COG @XYZ',
     &      xc,yc,zc
          if(act.eq.'v') call edisp(iuout,outs)
          if(inthenetwork)then
            call symbol_flow_cmp(name,'or','r',xc,yc,zc)
          else
            call symbol_flow_cmp(name,'or','g',xc,yc,zc)
          endif
          call pausems(30)
        endif

      elseif(SUSE(iz,is,1)(1:6).eq.'WINDOW'.or.
     &    SUSE(iz,is,1)(1:8).eq.'D-WINDOW'.or.
     &    SUSE(iz,is,1)(1:8).eq.'S-WINDOW'.or.
     &    SUSE(iz,is,1)(1:8).eq.'C-WINDOW')then
        if(SUSE(iz,is,2)(1:6).eq.'CLOSED'.or.
     &     SUSE(iz,is,2)(1:5).eq.'CRACK')then

C Window crack - crack component, assume width = 1mm,
C length = window perimeter, at surface COG.
          call zsurfprm(iz,is,dupedges,perim)
          if(inthenetwork)then
            continue
          else
            write(name,'(a)') 'WiCr'
            xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
            zc=SURCOG(iz,is,3)
          endif
          write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),' ',name,
     &      ' window crack width 1mm length',perim,' @XYZ',xc,yc,zc
          if(act.eq.'v') call edisp(iuout,outs)
          if(inthenetwork)then
            call symbol_flow_cmp(name,'cr','r',xc,yc,zc)
          else
            call symbol_flow_cmp(name,'cr','g',xc,yc,zc)
          endif
          call pausems(30)
        elseif(SUSE(iz,is,2)(1:4).eq.'OPEN')then

C Window open - common orifice component, (type 40) IVALCM(10) MF040I
          if(inthenetwork)then
            continue
          else
            write(name,'(a)') 'WiOp'
          endif
          if(inthenetwork)then
            continue
          else
            xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
            zc=SURCOG(iz,is,3)
          endif
          write(outs,'(4a,f6.1,a,f6.1,a,f6.2,a,3f7.3)') 
     &      zsn(1:lnzsn),' ',name,
     &      ' window open within surface width ',XYMAX,' height',ZMAX,
     &      ' area',SNA(iz,is),'m^2 COG @Z',xc,yc,zc
          if(act.eq.'v') call edisp(iuout,outs)
          if(inthenetwork)then
            call symbol_flow_cmp(name,'or','r',xc,yc,zc)
          else
            call symbol_flow_cmp(name,'or','g',xc,yc,zc)
          endif
          call pausems(30)
        elseif(SUSE(iz,is,2)(1:4).eq.'SASH')then

C Window sash - 2 common orifice component, (type 40) IVALCM(10) MF040I
          if(inthenetwork)then
            continue
          else
            write(name,'(a)') 'WiSL'
            xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
            zc=edgelower
          endif
          write(outs,'(4a,f6.1,a,f6.1,a,f6.2,a,3f6.2)') 
     &      zsn(1:lnzsn),' ',name,
     &      ' within surface of width ',XYMAX,' height',ZMAX,
     &      ' area',SNA(iz,is),'m2 lower @XYZ',xc,yc,zc
          if(act.eq.'v') call edisp(iuout,outs)
          if(inthenetwork)then
            call symbol_flow_cmp(name,'or','r',xc,yc,zc)
          else
            call symbol_flow_cmp(name,'or','g',xc,yc,zc)
          endif
          call pausems(30)

C Now create the upper orifice.
          if(inthenetwork)then
            continue
          else
            write(name,'(a)') 'WiSU'
            xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
            zc=edgeupper
          endif
          write(outs,'(4a,f6.1,a,f6.1,a,f6.2,a,3f6.2)') 
     &      zsn(1:lnzsn),' ',name,
     &      ' width ',XYMAX,' height',ZMAX*0.125,
     &      ' area',(XYMAX*(ZMAX*0.125)),'m2 upper @XYZ',xc,yc,zc
          if(act.eq.'v') call edisp(iuout,outs)
          if(inthenetwork)then
            call symbol_flow_cmp(name,'or','r',xc,yc,zc)
          else
            call symbol_flow_cmp(name,'or','g',xc,yc,zc)
          endif
          call pausems(30)
        elseif(SUSE(iz,is,2)(1:5).eq.'BIDIR'.or.
     &         SUSE(iz,is,2)(1:10).eq.'AJAR-BIDIR')then

C Window bidirectional - specific door component width & height
          if(inthenetwork)then
            continue
          else
            if(SUSE(iz,is,2)(1:5).eq.'BIDIR')then
              write(name,'(a)') 'WiBi'
            elseif(SUSE(iz,is,2)(1:10).eq.'AJAR-BIDIR')then
              write(name,'(a)') 'WABi'
            endif
            xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
            zc=SURCOG(iz,is,3)
          endif
          VAL3=ABS(ZCOG(iz,3)-edgelower)  ! base to zone COG
          if(SUSE(iz,is,2)(1:5).eq.'BIDIR')then
            write(outs,'(4a,f6.1,a,f6.1,a,3f7.3,a,f6.3)') 
     &      zsn(1:lnzsn),' ',name,
     &      ' bidir width ',XYMAX,' height',ZMAX,
     &      ' coef 0.6 @XYZ',xc,yc,zc,
     &      ' delta base to zone COG',VAL3
          elseif(SUSE(iz,is,2)(1:10).eq.'AJAR-BIDIR')then
            write(outs,'(4a,f6.1,a,3f7.3,a,f6.3)') 
     &      zsn(1:lnzsn),' ',name,
     &      ' bidir width 100mm height',ZMAX,
     &      ' coef 0.6 @XYZ',xc,yc,zc,
     &      ' delta base to zone COG',VAL3
          endif
          if(act.eq.'v') call edisp(iuout,outs)
          if(inthenetwork)then
            call symbol_flow_cmp(name,'bi','r',xc,yc,zc)
          else
            call symbol_flow_cmp(name,'bi','g',xc,yc,zc)
          endif
          call pausems(30)
        endif

      elseif(SUSE(iz,is,1)(1:4).eq.'FICT')then
        if(SUSE(iz,is,2)(1:6).eq.'CLOSED'.or.
     &     SUSE(iz,is,2)(1:5).eq.'CRACK')then

C Fict crack - crack component, assume width = 1mm,
C length = opening perimeter, at lowest horizontal edge.
          call zsurfprm(iz,is,dupedges,perim)
          if(inthenetwork)then
            continue
          else
            write(name,'(a)') 'FiCR'
            xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
            zc=edgelower
          endif
          write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),' ',name,
     &      ' fict 1mm crack length ',perim,' @XYZ',xc,yc,zc
          if(act.eq.'v') call edisp(iuout,outs)
          if(inthenetwork)then
            call symbol_flow_cmp(name,'cr','r',xc,yc,zc)
          else
            call symbol_flow_cmp(name,'cr','g',xc,yc,zc)
          endif
          call pausems(30)
        elseif(SUSE(iz,is,2)(1:4).eq.'OPEN')then

C Fict opening - common orifice component, (type 40) IVALCM(10) MF040I
          if(inthenetwork)then
            continue
          else
            write(name,'(a)') 'FiOP'
          endif
          if((XYMAX/ZMAX).gt.1.0)then
            if(inthenetwork)then
              continue
            else
              xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
              zc=SURCOG(iz,is,3)
            endif
            write(outs,'(4a,f6.1,a,f6.1,a,f6.2,a,3f7.3)') 
     &      zsn(1:lnzsn),' ',name,
     &      ' fict open width within ',XYMAX,' height',ZMAX,
     &      ' area',SNA(iz,is),'m^2 COG @Z',xc,yc,zc
          else
            if(inthenetwork)then
              continue
            else
              xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
              zc=SURCOG(iz,is,3)
            endif
            write(outs,'(4a,f6.1,a,f6.1,a,f6.2,a,3f7.3)') 
     &      zsn(1:lnzsn),' ',name,
     &      ' fict open width within ',XYMAX,' height',ZMAX,
     &      ' area',SNA(iz,is),'m^2 COG @Z',xc,yc,zc
          endif
          if(act.eq.'v') call edisp(iuout,outs)
          if(inthenetwork)then
            call symbol_flow_cmp(name,'or','r',xc,yc,zc)
          else
            call symbol_flow_cmp(name,'or','g',xc,yc,zc)
          endif
          call pausems(30)
        elseif(SUSE(iz,is,2)(1:5).eq.'BIDIR'.or.
     &         SUSE(iz,is,2)(1:10).eq.'AJAR-BIDIR')then

C Fict bi-directional - door component IVALCM(14)=130 MF130I
          if(inthenetwork)then
            continue
          else
            if(SUSE(iz,is,2)(1:5).eq.'BIDIR')then
              write(name,'(a)') 'FiBI'
            elseif(SUSE(iz,is,2)(1:10).eq.'AJAR-BIDIR')then
              write(name,'(a)') 'FABI'
            endif
            xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
            zc=edgelower
          endif
          write(outs,'(4a,f6.1,a,f6.1,a,3f7.3)') 
     &      zsn(1:lnzsn),' ',name,
     &      ' fict bidir within width ',XYMAX,' height',ZMAX,
     &      ' coef 0.6 base @XYZ',xc,yc,zc
          if(act.eq.'v') call edisp(iuout,outs)
          if(inthenetwork)then
            call symbol_flow_cmp(name,'bi','r',xc,yc,zc)
          else
            call symbol_flow_cmp(name,'bi','g',xc,yc,zc)
          endif
          call pausems(30)
        endif
      endif
      return
      end

C ***************** draw_node_assoc
C draw_connect_assoc draws the connections & components associated
C with flow node X.
      subroutine draw_node_assoc(IFNOD)
#include "building.h"
#include "net_flow.h"
#include "net_flow_data.h"
#include "model.h"
#include "geometry.h"
#include "prj3dv.h"

      integer lnblnk  ! function definition

      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      COMMON/FILEP/IFIL
      COMMON/SPAD/MMOD,LIMIT,LIMTTY
      COMMON/AFN/IAIRN,LAPROB,ICAAS(MCOM)
      INTEGER :: iairn,icaas
      CHARACTER LAPROB*72
      COMMON/gzonpik/izgfoc,nzg,nznog(mcom)
      integer ncomp,ncon
      COMMON/C1/NCOMP,NCON

      CHARACTER outs*124
      logical inview
      logical found
      CHARACTER name*12,zn*12,sn*12,twochar*2
      character altcol*1
      character col*1
      real xc,yc,zc,xw,yw,zw

C If 3D and graphic mode draw the zones associated with the
C connection if they are internal known/unknown as well as 
c adjacent connections (in grey) and highlight the current connection.
      nbzones=0; izn1=0; izn2=0; iposn=0; inegn=0; nzg=0
C If 3D and graphic mode draw the node and adjacent connections.
C If associated with a specific zone also read and draw it.
      if(MMOD.eq.8)then  ! in text mode do not draw
        found=.false.
        DO loop=1,NCOMP
          if(ICAAS(loop).gt.0)then
            i=ICAAS(loop)
            if(NDNAM(i).eq.NDNAM(IFNOD))then
              write(outs,'(3A)')zname(loop)(1:lnzname(loop)),
     &          ' is currently linked to: ',NDNAM(I)
              call edisp(iuout,' ')
              call edisp(iuout,outs)
              CALL INLNST(1)
              nzg=1
              nznog(1)=loop
              izgfoc=loop
              call georead(IFIL+1,LGEOM(loop),loop,0,iuout,IER)
              MODIFYVIEW=.TRUE.; MODBND=.TRUE.
              ITORG=1; ITGRD=1; ITBND=1  ! do not show origin or grid
              CALL redraw(IER)
              found=.true.
              exit   ! no need to loop any more
            endif
          endif
        enddo  ! of loop

C If not found then it is a boundary node. Add logic to draw
C it and the adjacent zone.
        if(.NOT.found)then

C Check associated components for links to a zone.
          nbzones=0
          DO ICON=NCNN,1,-1
            IF(NODPS(ICON).EQ.IFNOD.OR.NODNE(ICON).EQ.IFNOD)THEN
              zn=cmpassoc(itpcon(icon),1)  ! get associated zone
              lnzn=lnblnk(zn)
              sn=cmpassoc(itpcon(icon),2)  ! get associated surface
              lnsn=lnblnk(sn)
              WRITE(outs,'(9A)')
     &          NDNAM(NODPS(ICON)),' --> ',NDNAM(NODNE(ICON)),
     &          ' via ',CMNAM(ITPCON(ICON)),' associated with ',
     &          zn(1:lnzn),':',sn(1:lnsn)
              CALL EDISP(IUOUT,' ')
              CALL EDISP(IUOUT,OUTS)
              do loop=1,NCOMP
                 lntzn=lnblnk(zname(loop))
                 if(zn(1:lnzn).eq.zname(loop)(1:lntzn))then
                   nbzones=nbzones+1
                   nzg=nbzones
                   nznog(nbzones)=loop
                   izgfoc=loop
                 endif
               enddo
               if(nbzones.eq.1)then
                  call georead(IFIL+1,LGEOM(izgfoc),izgfoc,0,iuout,IER)
                  CALL INLNST(1)
                  MODIFYVIEW=.TRUE.; MODBND=.TRUE.
                  ITORG=1; ITGRD=1; ITBND=1 ! do not show origin or grid
                  CALL redraw(IER)
               elseif(nbzones.gt.1)then
                  CALL INLNST(1)
                  MODIFYVIEW=.TRUE.; MODBND=.TRUE.
                  ITORG=1; ITGRD=1; ITBND=1  ! do not show origin or grid
                  CALL redraw(IER)
               endif
             endif
          ENDDO
        endif
        col='r'     ! red for current node.
        altcol='b'  ! blue alternative colour
        name=NDNAM(IFNOD)
        xw=HNOD(IFNOD,1);yw=HNOD(IFNOD,2); zw=HNOD(IFNOD,3)

C Check within viewing angle and draw node symbol.
        call withinview(xw,yw,zw,inview,ier)
        if(inview)then
          if(NDTYP(IFNOD).eq.0)then
            if(NDFLD(IFNOD).eq.1)then
              call symbol_flow_nd(name,'iu',col,xw,yw,zw)
            else
              call symbol_flow_nd(name,'iu',altcol,xw,yw,zw)
            endif
            call pausems(30)
          elseif(NDTYP(IFNOD).eq.1)then
            if(NDFLD(IFNOD).eq.1)then
              call symbol_flow_nd(name,'ik',col,xw,yw,zw)
            else
              call symbol_flow_nd(name,'ik',altcol,xw,yw,zw)
            endif
            call pausems(30)
          elseif(NDTYP(IFNOD).eq.2)then
            call symbol_flow_nd(name,'bk',col,xw,yw,zw)
            call pausems(30)
          elseif(NDTYP(IFNOD).eq.3)then
            call symbol_flow_nd(name,'bw',col,xw,yw,zw)
            call pausems(30)
          endif
                
C Display associated connections.
          DO ICON=NCNN,1,-1
            IF(NODPS(ICON).EQ.IFNOD.OR.NODNE(ICON).EQ.IFNOD)THEN
              zn=cmpassoc(itpcon(icon),1)  ! get associated zone
              lnzn=lnblnk(zn)
              sn=cmpassoc(itpcon(icon),2)  ! get associated surface
              lnsn=lnblnk(sn)
              WRITE(OUTS,'(A,I3,10A)') 
     &          'Used by connection... ',ICON,' ',
     &          NDNAM(NODPS(ICON)),' --> ',NDNAM(NODNE(ICON)),
     &          ' via ',CMNAM(ITPCON(ICON)),' associated with ',
     &          zn(1:lnzn),':',sn(1:lnsn)
              CALL EDISP(IUOUT,OUTS)

C Draw the associated components.
              name=CMNAM(ITPCON(ICON))

              IFCN=ITPCON(ICON)
              call get_component_code(IFCN,twochar)

              xc=HCMP(ITPCON(ICON),1,1); yc=HCMP(ITPCON(ICON),1,2)
              zc=HCMP(ITPCON(ICON),1,3)
              IF(ITPCMP(ITPCON(ICON)).EQ. 10) THEN
                call symbol_flow_cmp(name,twochar,'b',xc,yc,zc)
              ELSE
                call symbol_flow_cmp(name,twochar,col,xc,yc,zc)
              ENDIF

C Draw the connection betwen node and component (but not to the other node).
              call draw_flow_link(xw,yw,zw,xc,yc,zc,col)
            endif
          ENDDO
        endif    ! node within viewing area
      endif      ! of graphic mode
      
      end

C ***************** draw_connect_assoc
C draw_connect_assoc draws the nodes & components associated
C with flow connection IFCN.
      subroutine draw_connect_assoc(IFCN)
#include "building.h"
#include "net_flow.h"
#include "net_flow_data.h"
#include "model.h"
#include "geometry.h"
#include "prj3dv.h"

      integer lnblnk  ! function definition

      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      COMMON/FILEP/IFIL
      COMMON/AFN/IAIRN,LAPROB,ICAAS(MCOM)
      INTEGER :: iairn,icaas
      CHARACTER LAPROB*72
      COMMON/gzonpik/izgfoc,nzg,nznog(mcom)
      integer ncomp,ncon
      COMMON/C1/NCOMP,NCON

      CHARACTER outs*124
      logical inview
      CHARACTER name*12,zn*12,sn*12,twochar*2
      character altcol*1
      character col*1
      real xc,yc,zc,xw,yw,zw

C If 3D and graphic mode draw the zones associated with the
C connection if they are internal known/unknown as well as 
c adjacent connections (in grey) and highlight the current connection.
      nbzones=0; izn1=0; izn2=0; iposn=0; inegn=0; nzg=0
      if(NDTYP(NODPS(IFCN)).eq.0.or.NDTYP(NODPS(IFCN)).eq.1)then
        do loop=1,NCOMP
          if(ICAAS(loop).gt.0)then
            i=ICAAS(loop)
            if(NDNAM(i).eq.NDNAM(NODPS(IFCN)))then
              write(outs,'(3A)')zname(loop)(1:lnzname(loop)),
     &          ' is currently linked to +node: ',NDNAM(I)
              call edisp(iuout,' ')
              call edisp(iuout,outs)
              nbzones=nbzones+1
              nzg=nzg+1
              nznog(nbzones)=loop
              izgfoc=loop
              izn1=loop; iposn=NODPS(IFCN) ! remember zone and pos node
              CYCLE  ! no need to check further
            endif
          endif
        enddo
      else
        iposn=NODPS(IFCN)
      endif
      if(NDTYP(NODNE(IFCN)).eq.0.or.NDTYP(NODNE(IFCN)).eq.1)then
        do loop=1,NCOMP
          if(ICAAS(loop).gt.0)then
            i=ICAAS(loop)
            if(NDNAM(i).eq.NDNAM(NODNE(IFCN)))then
              write(outs,'(3A)')zname(loop)(1:lnzname(loop)),
     &          ' is currently linked to -node: ',NDNAM(I)
              call edisp(iuout,' ')
              call edisp(iuout,outs)
              nbzones=nbzones+1  ! node positive
              nzg=nzg+1
              nznog(nbzones)=loop
              izgfoc=loop
              if(nbzones.eq.1)izn1=loop   ! remember zone
              if(nbzones.eq.2)then
                izn2=loop
                izgfoc=0    ! set to draw multiple zones
              endif
              inegn=NODNE(IFCN)           ! and negative node
              CYCLE  ! no need to check further
            endif
          endif
        enddo
      else
        inegn=NODNE(IFCN)
      endif
      if(nbzones.eq.1)then
        call georead(IFIL+1,LGEOM(izn1),izn1,0,iuout,IER)
        CALL INLNST(1)
        MODIFYVIEW=.TRUE.; MODBND=.TRUE.
        ITORG=1; ITGRD=1; ITBND=1  ! do not show origin or grid
        CALL redraw(IER)
      elseif(nbzones.gt.1)then
        CALL INLNST(1)
        MODIFYVIEW=.TRUE.; MODBND=.TRUE.
        ITORG=1; ITGRD=1; ITBND=1  ! do not show origin or grid
        CALL redraw(IER)
      endif

C Draw the positive node.
      col='r'     ! red for current node.
      altcol='b'  ! blue alternative colour
      name=NDNAM(iposn)
      xw=HNOD(iposn,1);yw=HNOD(iposn,2); zw=HNOD(iposn,3)

C Check within viewing angle and draw node symbol.
      call withinview(xw,yw,zw,inview,ier)
      if(inview)then
        if(NDTYP(iposn).eq.0)then
          if(NDFLD(iposn).eq.1)then
            call symbol_flow_nd(name,'iu',col,xw,yw,zw)
          else
            call symbol_flow_nd(name,'iu',altcol,xw,yw,zw)
          endif
          call pausems(30)
        elseif(NDTYP(iposn).eq.1)then
          if(NDFLD(iposn).eq.1)then
            call symbol_flow_nd(name,'ik',col,xw,yw,zw)
          else
            call symbol_flow_nd(name,'ik',altcol,xw,yw,zw)
          endif
          call pausems(30)
        elseif(NDTYP(iposn).eq.2)then
          call sym bol_flow_nd(name,'bk',col,xw,yw,zw)
          call pausems(30)
        elseif(NDTYP(iposn).eq.3)then
          call symbol_flow_nd(name,'bw',col,xw,yw,zw)
          call pausems(30)
        endif

C Draw connections associated with positive node.
        DO ICON=NCNN,1,-1
          IF(NODPS(ICON).EQ.iposn.OR.NODNE(ICON).EQ.iposn)THEN
            if(ICON.eq.IFCN)then
              col='r'     ! red for current connection.
            else
              col='g'     ! grey for current connection.
            endif
            zn=cmpassoc(itpcon(icon),1)  ! get associated zone
            lnzn=lnblnk(zn)
            sn=cmpassoc(itpcon(icon),2)  ! get associated surface
            lnsn=lnblnk(sn)
            WRITE(OUTS,'(A,I3,10A)') 
     &        'Connection... ',ICON,' ',
     &        NDNAM(NODPS(ICON)),' --> ',NDNAM(NODNE(ICON)),
     &        ' via ',CMNAM(ITPCON(ICON)),' associated with ',
     &        zn(1:lnzn),':',sn(1:lnsn)
            CALL EDISP(IUOUT,OUTS)

C Draw the associated components.
            name=CMNAM(ITPCON(ICON))
            xc=HCMP(ITPCON(ICON),1,1); yc=HCMP(ITPCON(ICON),1,2)
            zc=HCMP(ITPCON(ICON),1,3)

            IFCNT=ITPCON(ICON)
            call get_component_code(IFCNT,twochar)
            IF(ITPCMP(ITPCON(ICON)).EQ. 10) THEN
              call symbol_flow_cmp(name,twochar,'b',xc,yc,zc)
            ELSE
              call symbol_flow_cmp(name,twochar,col,xc,yc,zc)
            ENDIF

C Draw the connection betwen node and component (but not to the other node).
            call draw_flow_link(xw,yw,zw,xc,yc,zc,col) ! comp to positive
          endif
        enddo  ! of ICON
      endif

C Draw negative node.
      col='r'     ! red for current node.
      altcol='b'  ! blue alternative colour
      name=NDNAM(inegn)
      xw=HNOD(inegn,1);yw=HNOD(inegn,2); zw=HNOD(inegn,3)

C Check within viewing angle and draw node symbol.
      call withinview(xw,yw,zw,inview,ier)
      if(inview)then
        if(NDTYP(inegn).eq.0)then
          if(NDFLD(inegn).eq.1)then
            call symbol_flow_nd(name,'iu',col,xw,yw,zw)
          else
            call symbol_flow_nd(name,'iu',altcol,xw,yw,zw)
          endif
          call pausems(30)
        elseif(NDTYP(inegn).eq.1)then
          if(NDFLD(inegn).eq.1)then
            call symbol_flow_nd(name,'ik',col,xw,yw,zw)
          else
            call symbol_flow_nd(name,'ik',altcol,xw,yw,zw)
          endif
          call pausems(30)
        elseif(NDTYP(inegn).eq.2)then
          call symbol_flow_nd(name,'bk',col,xw,yw,zw)
          call pausems(30)
        elseif(NDTYP(inegn).eq.3)then
          call symbol_flow_nd(name,'bw',col,xw,yw,zw)
          call pausems(30)
        endif

C Draw connections (grey) associated with negative node.
        DO ICON=NCNN,1,-1
          IF(NODPS(ICON).EQ.inegn.OR.NODNE(ICON).EQ.inegn)THEN
            if(ICON.eq.IFCN)then
              col='r'     ! red for current connection.
            else
              col='g'     ! red for current connection.
            endif
            zn=cmpassoc(itpcon(icon),1)  ! get associated zone
            lnzn=lnblnk(zn)
            sn=cmpassoc(itpcon(icon),2)  ! get associated surface
            lnsn=lnblnk(sn)
            WRITE(OUTS,'(A,I3,10A)') 
     &        'Connection... ',ICON,' ',
     &        NDNAM(NODPS(ICON)),' --> ',NDNAM(NODNE(ICON)),
     &        ' via ',CMNAM(ITPCON(ICON)),' associated with ',
     &        zn(1:lnzn),':',sn(1:lnsn)
            CALL EDISP(IUOUT,OUTS)

C Draw the associated components.
            name=CMNAM(ITPCON(ICON))
            xc=HCMP(ITPCON(ICON),1,1); yc=HCMP(ITPCON(ICON),1,2)
            zc=HCMP(ITPCON(ICON),1,3)

            IFCNT=ITPCON(ICON)
            call get_component_code(IFCNT,twochar)
            IF(ITPCMP(ITPCON(ICON)).EQ. 10) THEN
              call symbol_flow_cmp(name,twochar,'b',xc,yc,zc)
            ELSE
              call symbol_flow_cmp(name,twochar,col,xc,yc,zc)
            ENDIF

C Draw the connection betwen node and component (but not to the other node).
            call draw_flow_link(xw,yw,zw,xc,yc,zc,col) ! comp to negative
          endif
        enddo  ! of ICON
      endif

C Force highlight both legs of connection plus component.
      col='r'     ! red for current connection.
      xw=HNOD(iposn,1);yw=HNOD(iposn,2); zw=HNOD(iposn,3)
      xc=HCMP(ITPCON(IFCN),1,1); yc=HCMP(ITPCON(IFCN),1,2)
      zc=HCMP(ITPCON(IFCN),1,3)
      call get_component_code(IFCN,twochar)
      call symbol_flow_cmp(name,twochar,col,xc,yc,zc)
      call draw_flow_link(xw,yw,zw,xc,yc,zc,col) ! comp to positive
      xw=HNOD(inegn,1);yw=HNOD(inegn,2); zw=HNOD(inegn,3)
      call draw_flow_link(xw,yw,zw,xc,yc,zc,col) ! comp to negative
      return
      end    ! of draw_connect_assoc

C ***************** get_component_code
C get_component_code returns a 2 character code for the component
C type to use when drawing components.
C with flow connection IFCN.
      subroutine get_component_code(IFCN,twochar)
#include "building.h"
#include "net_flow.h"
#include "net_flow_data.h"

      integer IFCN    ! flow component
      character twochar*2

      twochar='--'
      IF(ITPCMP(IFCN).EQ. 10) THEN
        twochar='bx'
      ELSEIF(ITPCMP(IFCN).EQ. 12) THEN
        twochar='bx' ! power law vol flow (box)
      ELSEIF(ITPCMP(IFCN).EQ. 15) THEN
        twochar='bx' ! power law mass flow (box)
      ELSEIF(ITPCMP(IFCN).EQ. 17) THEN
        twochar='bx' ! power law mass flow (box)
      ELSEIF(ITPCMP(IFCN).EQ. 20) THEN
        twochar='bx' ! quadratic law vol flow (box)
      ELSEIF(ITPCMP(IFCN).EQ. 25) THEN
        twochar='bx' ! quadratic law mass flow (box)
      ELSEIF(ITPCMP(IFCN).EQ. 30) THEN
        twochar='bx' ! constant vol flow (box)
      ELSEIF(ITPCMP(IFCN).EQ. 35) THEN
        twochar='bx' ! constant mass flow (box)
      ELSEIF(ITPCMP(IFCN).EQ. 40) THEN
        twochar='or' ! Common orifice flow
      ELSEIF(ITPCMP(IFCN).EQ. 50) THEN
        twochar='cd' ! Laminar pipe vol
      ELSEIF(ITPCMP(IFCN).EQ. 110) THEN
        twochar='or' ! Specific air flow
      ELSEIF(ITPCMP(IFCN).EQ. 120) THEN
        twochar='cr' ! Specific air flow crack
      ELSEIF(ITPCMP(IFCN).EQ. 130) THEN
        twochar='bi' ! Specific air flow door
      ELSEIF(ITPCMP(IFCN).EQ. 210) THEN
        twochar='cd' ! General flow conduit
      ELSEIF(ITPCMP(IFCN).EQ. 211) THEN
        twochar='bx' ! Cowls and roof outlets (box)
      ELSEIF(ITPCMP(IFCN).EQ. 220) THEN
        twochar='cd' ! Conduit ending in conv 3 leg
      ELSEIF(ITPCMP(IFCN).EQ. 230) THEN
        twochar='cd' ! Conduit starts in diverg 3 leg
      ELSEIF(ITPCMP(IFCN).EQ. 240) THEN
        twochar='cd' ! Conduit ending in conv 4 leg
      ELSEIF(ITPCMP(IFCN).EQ. 250) THEN
        twochar='cd' ! Conduit starts in diverg 4 leg
      ELSEIF(ITPCMP(IFCN).EQ. 310) THEN
        twochar='bx' ! General flow inducer (box)
      ELSEIF(ITPCMP(IFCN).EQ. 410) THEN
        twochar='bx' ! General flow corrector (box)
      ELSEIF(ITPCMP(IFCN).EQ. 420) THEN
        twochar='bx' ! Corrector polynomial (box)
      ELSEIF(ITPCMP(IFCN).EQ. 460) THEN
        twochar='bx' ! Fixed flow rates controller (box)
      ELSEIF(ITPCMP(IFCN).EQ. 500) THEN
        twochar='bx' ! multi-configuration (box)
      ENDIF
      return
      end  ! of get_component_code

C ***************** mfcompposition
C mfcompposition does user interactions for flow component position.
C passed component index and returns associated zone & surface, the
C position and several useful surface derived attributes.

C A future extension to make inferences and set surface flow attributes
C would save lots of duplication of effort.

      SUBROUTINE mfcompposition(ICMP,iz,iso,cx,cy,cz,carea,cperim,
     &  cwidth,cheight,IER)
#include "building.h"
#include "model.h"
#include "geometry.h"
#include "net_flow.h"
#include "net_flow_data.h"
#include "prj3dv.h"

C Parameters passed
      integer ICMP        ! index of component
      real cx,cy,cz       ! returned position
      real carea,cperim   ! returned surface area and perimeter
      real cwidth,cheight ! returned surface width & height
      integer IER  ! 0=use all returns, 1 use only position,  3 not 3D file

      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      COMMON/FILEP/IFIL
      COMMON/AFN/IAIRN,LAPROB,ICAAS(MCOM)
      INTEGER :: iairn,icaas
      CHARACTER LAPROB*72
      COMMON/FOPENED/CFGOK,MLDBOK,MATDBOK,CTLOK,OPTKOK,CFCDBOK
      LOGICAL        CFGOK,MLDBOK,MATDBOK,CTLOK,OPTKOK,CFCDBOK
      common/gzonpik/izgfoc,nzg,nznog(mcom)

      CHARACTER outs*124,HOLD*36,MSG*48
      real edgelower,edgeupper  ! Z of current surface lower & upper edges.
      integer llpos,lrpos,ulpos,urpos     ! closest to BB corners for parent
      logical dupedges

C For help messages
      character helpinsub*24 ! subroutine name
      character helptopic*24 ! string (unique) for topic
      integer nbhelp     ! number of help lines found

      helpinsub='mfcompposition'  ! set for subroutine

      IER=0

C Synopsis of flow component.
      WRITE(outs,'(A)') LTPCMP(ICMP)
      call edisp(iuout,outs)
      helptopic='comp_position'
      call gethelptext(helpinsub,helptopic,nbhelp)

C If not 3D network exit.
      if(IAIRN.lt.3)then
        ier=3
        return
      endif

C Position in space if IAIRN=3.
      write(hold,'(f9.3,f9.3,f8.3,a)') HCMP(ICMP,1,1),HCMP(ICMP,1,2),
     &  HCMP(ICMP,1,3),'   '
      write(MSG,'(3a)') 'Component: ',CMNAM(ICMP),' position in space:'
 342  CALL EASKS2CMD(HOLD,MSG,' ',
     &  'accept current','via a specific surface',iclkok,36,
     &  ' 1. 0. 0.','position',iser,nbhelp)
      if(iclkok.eq.1)then    ! notice accept current so just return
        cx=HCMP(ICMP,1,1)
        cy=HCMP(ICMP,1,2)
        cz=HCMP(ICMP,1,3)
        carea=0.; cperim=0.; cwidth=0.; cheight=0.
        ier=1
        call usrmsg(' ',' ','-')
        return
      elseif(iclkok.eq.2)then

C Ask for zone and surface.
        CALL EASKGEOF('Select an existing zone, else exit.',
     &    CFGOK,IZ,'-',34,IER)
        if(IZ.GT.0)then

C General image option flags. Rescan the zone geometry file.
          call georead(IFIL+1,LGEOM(IZ),IZ,1,iuout,IER)
          CALL ESCZONE(IZ)
          CALL BNDOBJ(0,IER)
          CALL ERCZONE(IZ)
          ITDSP=1; ITBND=1; ITEPT=0; ITZNM=0; ITSNM=0
          ITVNO=1; ITORG=1; ITBND=1; ITSNR=1; ITGRD=1
          GRDIS=0.0; ITPPSW=0

C Reset all surface lines to standard width and display the zone.
          MODIFYVIEW=.TRUE.; MODBND=.TRUE.
          CALL INLNST(1)
          nzg=1
          nznog(1)=IZ
          izgfoc=IZ
          CALL redraw(IER)

          ISO=1
          CALL EPMENSV
          CALL EASKSUR(IZ,ISO,'A','Select surface ',
     &      'to associate with this flow entity.',IER)
          CALL EPMENRC

C Note: llpos, lrpos etc are the position in the iszjvn list at each
C of the corners of the surface.
          call zsurfprm(iz,iso,dupedges,perim)
          CALL ZSURLEHI(IZ,ISO,XYMAX,ZMAX,llpos,lrpos,ulpos,urpos)
          if(llpos.gt.0.and.ulpos.gt.0)then
            edgelower=szcoords(iz,iszjvn(iz,iso,llpos),3)
            edgeupper=szcoords(iz,iszjvn(iz,iso,ulpos),3)
          endif
          write(outs,'(a,3f8.3)') 'Surface COG is:',
     &      SURCOG(iz,iso,1),SURCOG(iz,iso,2),SURCOG(iz,iso,3)
          call edisp(iuout,outs)
          write(outs,'(a,3f8.3)') 'Surface upper edge is:',
     &      SURCOG(iz,iso,1),SURCOG(iz,iso,2),edgeupper
          call edisp(iuout,outs)
          write(outs,'(a,3f8.3)') 'Surface lower edge is:',
     &      SURCOG(iz,iso,1),SURCOG(iz,iso,2),edgelower
          call edisp(iuout,outs)
          write(outs,'(a,f7.1,f6.2)') 'Surface area & perimeter:',
     &      SNA(iz,iso),perim
          call edisp(iuout,outs)
          write(outs,'(a,f7.1,f6.2)') 'Surface width & height:',
     &      XYMAX,ZMAX
          call edisp(iuout,outs)
          carea=SNA(iz,iso)  ! instanciate suggested values
          cperim=perim
          cwidth=XYMAX
          cheight=ZMAX

C Present options COG, upper edge, lower edge or manual edit.
C In the future remember the COG/upper/lower/manual decision.
          CALL EASKABCD('Position options:',' ',
     &      'surface COG','surface upper edge',
     &      'surface lower edge','manual edit',IW,nbhelp)
          if(iw.eq.1)then
            cx=SURCOG(iz,iso,1)
            cy=SURCOG(iz,iso,2)
            cz=SURCOG(iz,iso,3)
            ier=0
            call usrmsg(' ',' ','-')
            return
          elseif(iw.eq.2)then
            cx=SURCOG(iz,iso,1)
            cy=SURCOG(iz,iso,2)
            cz=edgeupper
            ier=0
            call usrmsg(' ',' ','-')
            return
          elseif(iw.eq.3)then
            cx=SURCOG(iz,iso,1)
            cy=SURCOG(iz,iso,2)
            cz=edgelower
            ier=0
            call usrmsg(' ',' ','-')
            return
          else
            CALL EASKS(HOLD,'Position in space?',
     &        ' ',36,' 1.0  0.0  0.0 ','position',IER,nbhelp)
          endif
        else

C No selection so just return current position.
          cx=HCMP(ICMP,1,1)
          cy=HCMP(ICMP,1,2)
          cz=HCMP(ICMP,1,3)
          carea=0.; cperim=0.; cwidth=0.; cheight=0.
          ier=1
          call usrmsg(' ',' ','-')
          return
        endif

C Extract XYZ values from the string.
        K=0
        CALL EGETWR(HOLD,K,cx,-999.9,999.9,'W','X',IER)
        CALL EGETWR(HOLD,K,cy,-99.9,999.9,'W','Y',IER)
        CALL EGETWR(HOLD,K,cz,-99.9,999.0,'W','Z',IER)
        if(ier.ne.0)goto 342
        ier=0
        call usrmsg(' ',' ','-')
        return
      else

C Extract XYZ values from the string.
        K=0
        CALL EGETWR(HOLD,K,cx,-999.9,999.9,'W','X',IER)
        CALL EGETWR(HOLD,K,cy,-99.9,999.9,'W','Y',IER)
        CALL EGETWR(HOLD,K,cz,-99.9,999.0,'W','Z',IER)
        if(ier.ne.0)goto 342
        ier=1  ! signal only position changed
      endif
      call usrmsg(' ',' ','-')
      RETURN
      END

C ***************** mfnodeposition
C mfnodeposition does user interactions for flow node position.
C passed node index and returns associated zone & surface, the
C position and several useful surface derived attributes.
C If boundary wind or fixed pressure offset from the surface.

      SUBROUTINE mfnodeposition(INOD,iz,iso,cx,cy,cz,carea,cperim,
     &  cwidth,cheight,IER)
#include "building.h"
#include "model.h"
#include "geometry.h"
#include "net_flow.h"
#include "net_flow_data.h"
#include "prj3dv.h"

C Parameters passed
      integer INOD        ! index of node
      real cx,cy,cz       ! returned position
      real carea,cperim   ! returned surface area and perimeter
      real cwidth,cheight ! returned surface width & height
      integer IER  ! 0=use all returns, 1 use only position,  3 not 3D file

      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      COMMON/FILEP/IFIL
      COMMON/AFN/IAIRN,LAPROB,ICAAS(MCOM)
      INTEGER :: iairn,icaas
      CHARACTER LAPROB*72
      COMMON/FOPENED/CFGOK,MLDBOK,MATDBOK,CTLOK,OPTKOK,CFCDBOK
      LOGICAL        CFGOK,MLDBOK,MATDBOK,CTLOK,OPTKOK,CFCDBOK
      common/gzonpik/izgfoc,nzg,nznog(mcom)

      CHARACTER outs*124,HOLD*36
      real edgelower,edgeupper  ! Z of current surface lower & upper edges.
      integer llpos,lrpos,ulpos,urpos     ! closest to BB corners for parent
      logical dupedges

C For help messages
      character helpinsub*24 ! subroutine name
      character helptopic*24 ! string (unique) for topic
      integer nbhelp     ! number of help lines found

      helpinsub='mfnodeposition'  ! set for subroutine

      IER=0

C Synopsis of flow node.
      if(NDTYP(INOD).eq.0)then
        WRITE(outs,'(2A)') NDNAM(INOD),' internal unknown pressure'
      elseif(NDTYP(INOD).eq.1)then
        WRITE(outs,'(2A)') NDNAM(INOD),' internal known pressure'
      elseif(NDTYP(INOD).eq.2)then
        WRITE(outs,'(2A)') NDNAM(INOD),' external known pressure'
      elseif(NDTYP(INOD).eq.3)then
        WRITE(outs,'(2A)') NDNAM(INOD),' external wind pressure'
      endif
      call edisp(iuout,outs)
      helptopic='node_position'
      call gethelptext(helpinsub,helptopic,nbhelp)

C If not 3D network exit.
      if(IAIRN.lt.3)then
        ier=3
        return
      endif

C For calculating offset.
      PI = 4.0 * ATAN(1.0); RAD = PI/180.0;

C Position in space if IAIRN=3.
      write(hold,'(f9.3,f9.3,f8.3,a)') HNOD(INOD,1),HNOD(INOD,2),
     &  HNOD(INOD,3),'   '
 342  CALL EASKS2CMD(HOLD,'Node: position in space:',' ',
     &  'accept current','via a specific zone or surface',iclkok,36,
     &  ' 1. 0. 0.','position',iser,nbhelp)
      if(iclkok.eq.1)then    ! notice accept current so just return
        cx=HNOD(INOD,1)
        cy=HNOD(INOD,2)
        cz=HNOD(INOD,3)
        carea=0.; cperim=0.; cwidth=0.; cheight=0.
        ier=1
        call usrmsg(' ',' ','-')
        return
      elseif(iclkok.eq.2)then

C Ask for zone and optionally for surface.
        CALL EASKGEOF('Select an existing zone, else exit.',
     &    CFGOK,IZ,'-',34,IER)
        if(IZ.GT.0)then

C General image option flags. Rescan the zone geometry file.
          call georead(IFIL+1,LGEOM(IZ),IZ,1,iuout,IER)
          CALL ESCZONE(IZ)
          CALL BNDOBJ(0,IER)
          CALL ERCZONE(IZ)
          ITDSP=1; ITBND=1; ITEPT=0; ITZNM=0; ITSNM=0
          ITVNO=1; ITORG=1; ITBND=1; ITSNR=1; ITGRD=1
          GRDIS=0.0; ITPPSW=0

C Reset all surface lines to standard width and display the zone.
          MODIFYVIEW=.TRUE.; MODBND=.TRUE.
          CALL INLNST(1)
          nzg=1
          nznog(1)=IZ
          izgfoc=IZ
          CALL redraw(IER)
          MODIFYVIEW=.TRUE.; MODBND=.TRUE.

          if(NDTYP(INOD).eq.0)then

C If node is internal unknown then set the zone COG and return.
            cx=ZCOG(iz,1); cy=ZCOG(iz,2); cz=ZCOG(iz,3)
            carea=0.; cperim=0.; cwidth=0.; cheight=0.
            ier=1
            call usrmsg(' ',' ','-')
            return
          elseif(NDTYP(INOD).eq.1)then

C If node is internal known then set the zone COG and return.
            cx=ZCOG(iz,1); cy=ZCOG(iz,2); cz=ZCOG(iz,3)
            carea=0.; cperim=0.; cwidth=0.; cheight=0.
            ier=1
            call usrmsg(' ',' ','-')
            return
          endif

C For boundary known and wind pressure request a surface.
          ISO=1
          CALL EPMENSV
          CALL EASKSUR(IZ,ISO,'A','Select surface ',
     &      'to associate with this node.',IER)
          CALL EPMENRC

C Calculate offset.
          V1=0.9; AZ=spazi(iz,iso); EL=spelv(iz,iso)
          RYAZI = AZ*RAD; RSALT = EL*RAD
          z3 = V1*SIN(RSALT)
          XYDIS = V1*COS(RSALT)
          IF (XYDIS .LT. 1E-6)THEN
            x3 = 0.0; y3 = 0.0
          ELSE
            x3 = XYDIS*SIN(RYAZI); y3 = XYDIS*COS(RYAZI)
          ENDIF

C Note: llpos, lrpos etc are the position in the iszjvn list at each
C of the corners of the surface.
          call zsurfprm(iz,iso,dupedges,perim)
          CALL ZSURLEHI(IZ,ISO,XYMAX,ZMAX,llpos,lrpos,ulpos,urpos)
          if(llpos.gt.0.and.ulpos.gt.0)then
            edgelower=szcoords(iz,iszjvn(iz,iso,llpos),3)
            edgeupper=szcoords(iz,iszjvn(iz,iso,ulpos),3)
          endif
          write(outs,'(a,3f8.3)') 'Surface COG is:',
     &      SURCOG(iz,iso,1),SURCOG(iz,iso,2),SURCOG(iz,iso,3)
          call edisp(iuout,outs)
          write(outs,'(a,3f8.3)') 'Surface upper edge is:',
     &      SURCOG(iz,iso,1),SURCOG(iz,iso,2),edgeupper
          call edisp(iuout,outs)
          write(outs,'(a,3f8.3)') 'Surface lower edge is:',
     &      SURCOG(iz,iso,1),SURCOG(iz,iso,2),edgelower
          call edisp(iuout,outs)
          write(outs,'(a,f7.1,f6.2)') 'Surface area & perimeter:',
     &      SNA(iz,iso),perim
          call edisp(iuout,outs)
          write(outs,'(a,f7.1,f6.2)') 'Surface width & height:',
     &      XYMAX,ZMAX
          call edisp(iuout,outs)
          carea=SNA(iz,iso) ! instanciate suggested values
          cperim=perim
          cwidth=XYMAX
          cheight=ZMAX

C Present options COG, upper edge, lower edge or manual edit.
          CALL EASKABCD(
     &      'Position options (node will be 1m away from surface:',
     &      ' ','surface COG','surface upper edge',
     &      'surface lower edge','manual edit',IW,nbhelp)
          if(iw.eq.1)then
            cx=x3+SURCOG(iz,iso,1)
            cy=y3+SURCOG(iz,iso,2)
            cz=z3+SURCOG(iz,iso,3)
            ier=0
            call usrmsg(' ',' ','-')
            return
          elseif(iw.eq.2)then
            cx=x3+SURCOG(iz,iso,1)
            cy=y3+SURCOG(iz,iso,2)
            cz=z3+edgeupper
            ier=0
            call usrmsg(' ',' ','-')
            return
          elseif(iw.eq.3)then
            cx=x3+SURCOG(iz,iso,1)
            cy=y3+SURCOG(iz,iso,2)
            cz=z3+edgelower
            ier=0
            call usrmsg(' ',' ','-')
            return
          else
            CALL EASKS(HOLD,'Node: position in space (no offset)?',
     &        ' ',36,' 1.0  0.0  0.0 ','position',IER,nbhelp)
          endif
        else

C No selection so just return current position.
          cx=HNOD(INOD,1)
          cy=HNOD(INOD,2)
          cz=HNOD(INOD,3)
          carea=0.; cperim=0.; cwidth=0.; cheight=0.
          ier=1
          call usrmsg(' ',' ','-')
          return
        endif

C Extract XYZ values from the string.
        K=0
        CALL EGETWR(HOLD,K,cx,-999.9,999.9,'W','X',IER)
        CALL EGETWR(HOLD,K,cy,-99.9,999.9,'W','Y',IER)
        CALL EGETWR(HOLD,K,cz,-99.9,999.0,'W','Z',IER)
        if(ier.ne.0)goto 342
        ier=0
        call usrmsg(' ',' ','-')
        return
      else

C Extract XYZ values from the string.
        K=0
        CALL EGETWR(HOLD,K,cx,-999.9,999.9,'W','X',IER)
        CALL EGETWR(HOLD,K,cy,-99.9,999.9,'W','Y',IER)
        CALL EGETWR(HOLD,K,cz,-99.9,999.0,'W','Z',IER)
        if(ier.ne.0)goto 342
        ier=1  ! signal only position changed
      endif
      call usrmsg(' ',' ','-')
      RETURN
      END

C ***********************************
C scan_flow_use scanns surface USE attributes related to mass flows,
C extracts position information, fills flow common blocks and creates
C both legacy and a 3D flow file while drawing over wireframe.

      subroutine scan_flow_use()
#include "building.h"
#include "model.h"
#include "geometry.h"
#include "net_flow.h"
#include "net_flow_data.h"
#include "esprdbfile.h"
C esprdbfile.h supplies the following:
C LPRFDB,IPRODB (for event profile database)
#include "prj3dv.h"
      INTEGER, PARAMETER :: MNWKTYP=6  ! from gnetwk.h
      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      COMMON/FILEP/IFIL
      integer ncomp,ncon
      COMMON/C1/NCOMP,NCON
      COMMON/C3/IC1(MCON),IE1(MCON),ICT(MCON),IC2(MCON),IE2(MCON)
      COMMON/AFN/IAIRN,LAPROB,ICAAS(MCOM)
      INTEGER :: iairn,icaas
      character LAPROB*72
      COMMON/MFLOWIT/fndegc,imix

C To signal to MFWRIT that globals are available from graphic network.
C Currently set .FALSE.
      logical haveglobal
      COMMON/MFLOW11/haveglobal

      integer linkatpartition  ! have we done component on side of partition
      COMMON/linkpart/linkatpartition(MCON)

C Documentation shared with graphic flow network.
      COMMON/NWKSTR/LEGNWKNAM,NWKNAM,NWKDSC,NWKTYPSTR(MNWKTYP)
      CHARACTER LEGNWKNAM*72,NWKNAM*72,NWKDSC*72,NWKTYPSTR*12
      COMMON/NWKTYP/INWKTYP,vergnf
      INTEGER inwktyp
      REAL vergnf   ! 1.0 is 3D variant of ASCII network

      DIMENSION IVALS(MCOM)
      CHARACTER outs*144
      CHARACTER prompt*72,fs*1
      CHARACTER LTMP*72
      logical unixok
      integer loop
      integer ulpos,urpos     ! closest to BB corners for parent

C For help messages
      character helpinsub*24 ! subroutine name
      character helptopic*24 ! string (unique) for topic
      integer nbhelp     ! number of help lines found

      helpinsub='scan_flow_use'  ! set for subroutine

C Set folder separator (fs) to \ or / as required.
      call isunix(unixok)
      if(unixok)then
        fs = char(47)
      else
        fs = char(92)
      endif

      WRITE(OUTS,'(A)')
     &  'Which zones to include in the fluid flow network?'
      IZN=NCOMP

C Clear zone->node mapping and then generate their names & Z values.
      NNOD=0; NCMP=0; NCNN=0  ! Reset
      DO ICZ=1,IZN
        ICAAS(ICZ)=0
      ENDDO  ! of ICZ
      DO loop=1,NCON
        linkatpartition(loop)=0  ! connection not yet done
      ENDDO ! of loop
      CALL EPICKS(IZN,IVALS,' ',OUTS,12,NCOMP,zname,'zone list',
     &  IER,nbhelp)

      IF(IZN.EQ.0)THEN
        CALL USRMSG('No zones selected, returning.',' ','W')
        ier=1
        return
      ENDIF

C Temporarily use same file unit as profiles db.
      NDNAM(0)='    '
      fndegc=20.0
      imix=1
      vergnf=1.0  ! signal 3D flow file variant
      inwktyp=2   ! to match graphic flow network type
      IF(INWKTYP.EQ.0)THEN

C Set the network type (default is plant).
        INWKTYP=1
        NWKTYPSTR(1)='HVAC'
      ENDIF
      NWKTYPSTR(2)='Flow'
      NWKTYPSTR(3)='Electrical'
      NWKTYPSTR(4)='Control'
      NWKTYPSTR(5)='Hygroscopic'
      NWKTYPSTR(6)='PrimitivePt'

      if(IPRODB.eq.IFIL+6)then
        IUM=IPRODB
      else
        IUM=IFIL+6
      endif
      CALL ERPFREE(IUM,ISTAT)
      CALL MFCDAT

C Attempt to read in a flow file. If file has yet to be named then base
C it on cfgroot and place it in the netpth folder (differienciate between
C unix and non-unix machine types).
      if(LAPROB(1:2).eq.'  '.or.LAPROB(1:4).eq.'UNKN')then
        if(unixok)then
          if(netpth(1:2).eq.'  '.or.netpth(1:2).eq.'./')then
            WRITE(LAPROB,'(2A)')cfgroot(1:lnblnk(cfgroot)),'.afn'
            WRITE(NWKNAM,'(2A)')cfgroot(1:lnblnk(cfgroot)),'.3dafn'
            WRITE(LEGNWKNAM,'(2A)')cfgroot(1:lnblnk(cfgroot)),'.afn'
          else
            WRITE(LAPROB,'(4a)') netpth(1:lnblnk(netpth)),fs,
     &      cfgroot(1:lnblnk(cfgroot)),'.afn'
            WRITE(NWKNAM,'(4a)') netpth(1:lnblnk(netpth)),fs,
     &      cfgroot(1:lnblnk(cfgroot)),'.3dafn'
            WRITE(LEGNWKNAM,'(4a)') netpth(1:lnblnk(netpth)),fs,
     &      cfgroot(1:lnblnk(cfgroot)),'.afn'
          endif
        else
          if(netpth(1:2).eq.'  '.or.(ichar(netpth(1:1)).eq.46.and.
     &       ichar(netpth(2:2)).eq.92))then
            WRITE(LAPROB,'(2A)')cfgroot(1:lnblnk(cfgroot)),'.afn'
            WRITE(NWKNAM,'(2A)')cfgroot(1:lnblnk(cfgroot)),'.3dafn'
            WRITE(LEGNWKNAM,'(2A)')cfgroot(1:lnblnk(cfgroot)),'.afn'
          else
            WRITE(LAPROB,'(4a)') netpth(1:lnblnk(netpth)),fs,
     &      cfgroot(1:lnblnk(cfgroot)),'.afn'
            WRITE(NWKNAM,'(4a)') netpth(1:lnblnk(netpth)),fs,
     &      cfgroot(1:lnblnk(cfgroot)),'.3dafn'
            WRITE(LEGNWKNAM,'(4a)') netpth(1:lnblnk(netpth)),fs,
     &      cfgroot(1:lnblnk(cfgroot)),'.afn'
          endif
        endif
      else
        if(unixok)then
          if(netpth(1:2).eq.'  '.or.netpth(1:2).eq.'./')then
            WRITE(NWKNAM,'(2A)')cfgroot(1:lnblnk(cfgroot)),'.3dafn'
            WRITE(LEGNWKNAM,'(2A)')cfgroot(1:lnblnk(cfgroot)),'.afn'
          else
            WRITE(NWKNAM,'(4a)') netpth(1:lnblnk(netpth)),fs,
     &      cfgroot(1:lnblnk(cfgroot)),'.3dafn'
            WRITE(LEGNWKNAM,'(4a)') netpth(1:lnblnk(netpth)),fs,
     &      cfgroot(1:lnblnk(cfgroot)),'.afn'
          endif
        else
          if(netpth(1:2).eq.'  '.or.(ichar(netpth(1:1)).eq.46.and.
     &       ichar(netpth(2:2)).eq.92))then
            WRITE(NWKNAM,'(2A)')cfgroot(1:lnblnk(cfgroot)),'.3dafn'
            WRITE(LEGNWKNAM,'(2A)')cfgroot(1:lnblnk(cfgroot)),'.afn'
          else
            WRITE(NWKNAM,'(4a)') netpth(1:lnblnk(netpth)),fs,
     &      cfgroot(1:lnblnk(cfgroot)),'.3dafn'
            WRITE(LEGNWKNAM,'(4a)') netpth(1:lnblnk(netpth)),fs,
     &      cfgroot(1:lnblnk(cfgroot)),'.afn'
          endif
        endif
      endif
   78 LTMP=LAPROB
      CALL EASKS(LTMP,'Fluid flow network file?',
     &  '(non-graphic version)',72,' ','mass flow file',IER,nbhelp)
      if(LTMP(1:2).ne.'  '.and.LTMP(1:4).ne.'UNKN')then
        LAPROB=LTMP
      else
        goto 78
      endif
   79 LTMP=NWKNAM
      CALL EASKS(LTMP,'3D Fluid flow network file?',
     &  '(non-graphic 3D version)',72,' ','mass flow file',IER,nbhelp)
      if(LTMP(1:2).ne.'  '.and.LTMP(1:4).ne.'UNKN')then
        NWKNAM=LTMP
      else
        goto 79
      endif
   80 LTMP=NWKDSC
      CALL EASKS(LTMP,'3D Fluid flow description?',
     &  '(<72 char)',72,' ','mass flow descr',IER,nbhelp)
      if(LTMP(1:2).ne.'  '.and.LTMP(1:4).ne.'UNKN')then
        NWKDSC=LTMP
      else
        goto 80
      endif
      CALL EFOPSEQ(IUM,LAPROB,4,IER)
      CALL USRMSG(
     &  'Generating flow nodes for selected zones ...',' ','-')
      DO IZT=1,IZN

C Setup default internal flow node for zone IZT.
        call create_zone_node(IZT)
        call pausems(100)
      ENDDO

C Loop through all the surfaces in the model. Based on the use
C attributes build up a list of likely components.
      call edisp(iuout,' ')
      call edisp(iuout,
     & 'Scanning surface attributes for boundary nodes & connections.')
      PI = 4.0 * ATAN(1.0); RAD = PI/180.0;
      do loop=1,NCON
        iz=IC1(loop); is=IE1(loop)  ! get current zone & surface.

C Check if partition and linkatpartition already set.
C If so we do not need to create this component.
        if(ICT(loop).eq.3.and.linkatpartition(loop).gt.0)then
          CYCLE  ! no need to create anything for this surface
        endif

C Cases to ignore.
        if(SUSE(iz,is,1)(1:1).eq.'-'.or.
     &     SUSE(iz,is,1)(1:4).eq.'WALL'.or.
     &     SUSE(iz,is,1)(1:5).eq.'FLOOR'.or.
     &     SUSE(iz,is,1)(1:5).eq.'FURNI'.or.
     &     SUSE(iz,is,1)(1:7).eq.'ITEQUIP'.or.
     &     SUSE(iz,is,1)(1:5).eq.'PARTN'.or.
     &     SUSE(iz,is,1)(1:4).eq.'ROOF'.or.
     &     SUSE(iz,is,1)(1:5).eq.'STRUC'.or.
     &     SUSE(iz,is,1)(1:7).eq.'FIXTURE'.or.
     &     SUSE(iz,is,1)(1:6).eq.'PLANTS')then
          CYCLE
        endif

C Create a default boundary flow node or component based on surface attributes.
        call create_surface_node_cmp(iz,is)

      enddo
      haveglobal=.FALSE.  ! global preferences not yet available
      CALL MFWRIT(IUM)
      CALL ERPFREE(IUM,ISTAT)

C List what we have thus far.
      CALL MFLIST(iuout)

C Write 3D attributes and legacy version of flow networks.
      CALL ERPFREE(IUM,ISTAT)
      CALL EFOPSEQ(IUM,NWKNAM,4,IER)
      call MF3DWRIT(IUM)
      CALL ERPFREE(IUM,ISTAT)
      lnnwknam=lnblnk(NWKNAM)
      write(LEGNWKNAM,'(2a)') NWKNAM(1:lnnwknam-5),'.afn'
      CALL EFOPSEQ(IUM,LEGNWKNAM,4,IER)
      IAIRN=1  ! set so legacy format is written
      CALL MFWRIT(IUM)
      CALL ERPFREE(IUM,ISTAT)
      IAIRN=3  ! reset to 3D attributes format

      call edisp(iuout,' ')
      call edisp(iuout,
     &  'Next task is to confirm links between the nodes & zones.')

C Link nodes and zones.
      DO 33 IZ=1,NCOMP
        if(ICAAS(IZ).gt.0)then
          write(outs,'(3A)')zname(IZ)(1:lnzname(IZ)),
     &      ' is currently linked to: ',NDNAM(ICAAS(IZ))
          helptopic='link_node_to_zone'
          call gethelptext(helpinsub,helptopic,nbhelp)
          CALL EASKABC(outs,' ','accept',
     &      'select another','free link',iacc,nbhelp)
          if(iacc.eq.1)then
            goto 33
          elseif(iacc.eq.2)then
            IC=0
            call ASKRNOD(' available nodes','-',IC,IER)
            ICAAS(IZ)=IC
          else
            ICAAS(IZ)=0
          endif
        else
          helptopic='link_node_to_zone'
          call gethelptext(helpinsub,helptopic,nbhelp)
          write(outs,'(A,A)') zname(IZ)(1:lnzname(IZ)),
     &                        ' has no mass flow node.'
          CALL EASKAB(outs,' ','accept','select a node',iacc,nbhelp)
          if(iacc.eq.1)then
            ICAAS(IZ)=0
          elseif(iacc.eq.2)then
            IC=0
            call ASKRNOD(' available nodes','-',IC,IER)
            ICAAS(IZ)=IC
          endif
        endif
  33  CONTINUE
      IAIRN=3  ! so it shows up in cfg file.

C Update the configuration file so that connections and mass flow
C network are know.
      call edisp(iuout,' ')
      call edisp(iuout,' Updating configuration mass flow links...')
      CALL EMKCFG('-',IER)

      return
      end

C Subroutine to setup common blocks for a default internal flow
C node for a zone. Calling code assumed to call MFWRIT & MF3DWRIT
C after this subroutine returns.
      subroutine create_zone_node(izone)
#include "building.h"
#include "model.h"
#include "geometry.h"
#include "net_flow.h"
#include "net_flow_data.h"
C LPRFDB,IPRODB (for event profile database)
#include "prj3dv.h"
      INTEGER, PARAMETER :: MNWKTYP=6  ! from gnetwk.h
      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      COMMON/AFN/IAIRN,LAPROB,ICAAS(MCOM)
      INTEGER :: iairn,icaas
      character LAPROB*72
      COMMON/MFLOWIT/fndegc,imix

C Zone fluid.
C  zfldK - conductivity;
C  zfldD - density;
C  zfldC - specific heat capacity;
C  zflsA - total shortwave absorptivity;
C  zSWA - shortwave absorption (W).
      COMMON/zfluid/znotair(mcom),zfldK,zfldD,zfldC,zfldA,
     &              zSWAp(mcom),zSWAf(mcom)
      real zfldK,zfldD,zfldC,zfldA,zSWAp,zSWAf
      LOGICAL znotair

      CHARACTER outs*144,name*12
      logical unixok

C Set folder separator (fs) to \ or / as required.
      call isunix(unixok)

      IZ=izone
      NNOD=NNOD+1

C Get zone information. For the zone associated with the current
C surface instanciate pointer from zone to the node.
      ICAAS(IZ)=NNOD

C Set node type to internal and fill in the data on volume, position
C and temperature, setting the node name to the zone name.
      NDTYP(NNOD)=0; NDFLD(NNOD)=1 ! internal unkn air
      HNOD(NNOD,1)=ZCOG(IZ,1); HNOD(NNOD,2)=ZCOG(IZ,2)
      HNOD(NNOD,3)=ZCOG(IZ,3)
      SUPNOD(NNOD,2)=VOL(IZ); TNOD(NNOD)=fndegc
      SUPNOD(NNOD,1)=1  ! initial assumption flow node as air
      if(znotair(iz))then
        SUPNOD(NNOD,1)=2  ! mark flow node as water
        NDFLD(NNOD)=2
      else
        SUPNOD(NNOD,1)=1  ! mark flow node as air
        NDFLD(NNOD)=1
      endif
      NDNAM(NNOD)=ZNAME(IZ)
      NODASSOC(NNOD,1)=ZNAME(IZ); NODASSOC(NNOD,2)='-'
      WRITE(OUTS,'(A,I3,3A)')'Node :',NNOD,' name: ',
     &  NDNAM(NNOD), ' type: internal, unknown pressure '
      CALL edisp(iuout,OUTS)
      WRITE(OUTS,'(A,3F9.3,A,F10.3)') ' XYZ(m): ',
     &  ZCOG(IZ,1),ZCOG(IZ,2),ZCOG(IZ,3), ' volume (m^3): ',
     &  SUPNOD(NNOD,2)
      CALL edisp(iuout,OUTS)
      write(name,'(a)') NDNAM(NNOD)
      if(znotair(iz))then
        call symbol_flow_nd(name,'iu','b',
     &    ZCOG(IZ,1),ZCOG(IZ,2),ZCOG(IZ,3))
      else
        call symbol_flow_nd(name,'iu','r',
     &    ZCOG(IZ,1),ZCOG(IZ,2),ZCOG(IZ,3))
      endif

      return
      end  ! of create_zone_node


C ***********************************
C create_surface_node_cmp scans surface USE attributes related to mass flows,
C extracts position information, fills flow common blocks.

      subroutine create_surface_node_cmp(izone,isurf)
#include "building.h"
#include "model.h"
#include "geometry.h"
#include "net_flow.h"
#include "net_flow_data.h"
#include "esprdbfile.h"
C esprdbfile.h supplies the following:
C LPRFDB,IPRODB (for event profile database)
#include "prj3dv.h"
      INTEGER, PARAMETER :: MNWKTYP=6  ! from gnetwk.h
      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      integer ncomp,ncon
      COMMON/C1/NCOMP,NCON
      COMMON/C3/IC1(MCON),IE1(MCON),ICT(MCON),IC2(MCON),IE2(MCON)
      integer IZSTOCN
      COMMON/C24/IZSTOCN(MCOM,MS)
      COMMON/AFN/IAIRN,LAPROB,ICAAS(MCOM)
      INTEGER :: iairn,icaas
      character LAPROB*72
      COMMON/MFLOWIT/fndegc,imix

      integer linkatpartition  ! have we done component on side of partition
      COMMON/linkpart/linkatpartition(MCON)

C Zone fluid.
C  zfldK - conductivity;
C  zfldD - density;
C  zfldC - specific heat capacity;
C  zflsA - total shortwave absorptivity;
C  zSWA - shortwave absorption (W).
      COMMON/zfluid/znotair(mcom),zfldK,zfldD,zfldC,zfldA,
     &              zSWAp(mcom),zSWAf(mcom)
      real zfldK,zfldD,zfldC,zfldA,zSWAp,zSWAf
      LOGICAL znotair

      DIMENSION IVALS2(MCOM)
      CHARACTER outs*144,name*12
      character CXSTR*78
      CHARACTER prompt*72,hold*32,fs*1,SN*12,zsn*28
      character ccol*1
      logical ok,unixok
      logical dupedges
      integer loop
      integer llpos,lrpos,ulpos,urpos     ! closest to BB corners for parent
      integer CMPIS  ! index of component we just created
      integer POS,NEG
      logical iswater
      real edgelower,edgeupper  ! Z of current surface lower & upper edges.
      real trwidth   ! width of trickle vent

C For help messages
      character helpinsub*24 ! subroutine name
      character helptopic*24 ! string (unique) for topic
      integer nbhelp     ! number of help lines found

      helpinsub='scan_flow_use'  ! set for subroutine

C Set folder separator (fs) to \ or / as required.
      call isunix(unixok)
      if(unixok)then
        fs = char(47)
      else
        fs = char(92)
      endif
      IF(izone.eq.0.or.isurf.EQ.0)THEN
        CALL USRMSG('Zone or surface index zero, returning.',' ','W')
        ier=1
        return
      ENDIF

C Temporarily use same file unit as profiles db.
      NDNAM(0)='    '
      fndegc=20.0
      imix=1

C For izone and isurf process based on the use
C attributes build up a list of likely components.
      PI = 4.0 * ATAN(1.0); RAD = PI/180.0;
      iz=izone; is=isurf
      loop=IZSTOCN(izone,isurf)

      write(SN,'(a)') SNAME(izone,isurf)  ! setup surface name context
      lsn=lnblnk(SNAME(izone,isurf))
      write(zsn,'(3a)') zname(iz)(1:lnzname(iz)),
     &  ':',SNAME(izone,isurf)(1:lsn)
      lnzsn=lnblnk(zsn)  ! remember length of combined zone:surface name

C Check if partition and linkatpartition already set.
C If so we do not need to create this component.
      if(ICT(loop).eq.3.and.linkatpartition(loop).gt.0)then
        return  ! no need to create anything for this surface
      endif
      name='  '
      call conxinfo(1,loop,CXSTR)

C Note: llpos, lrpos etc are the position in the iszjvn list at each
C of the corners of the surface.
      CALL ZSURLEHI(IZ,IS,XYMAX,ZMAX,llpos,lrpos,ulpos,urpos)
      if(llpos.gt.0.and.ulpos.gt.0)then
        edgelower=szcoords(iz,iszjvn(izone,isurf,llpos),3)
        edgeupper=szcoords(iz,iszjvn(izone,isurf,ulpos),3)
      endif

C Figure out the delta for the boundary node would be if
C 0.9m along normal vector.
      V1=0.9; AZ=spazi(izone,isurf); EL=spelv(izone,isurf)
      RYAZI = AZ*RAD; RSALT = EL*RAD
      z3 = V1*SIN(RSALT)
      XYDIS = V1*COS(RSALT)
      IF (XYDIS .LT. 1E-6)THEN
        x3 = 0.0; y3 = 0.0
      ELSE
        x3 = XYDIS*SIN(RYAZI); y3 = XYDIS*COS(RYAZI)
      ENDIF

C Cases to ignore.
      if(SUSE(iz,is,1)(1:1).eq.'-'.or.
     &   SUSE(iz,is,1)(1:4).eq.'WALL'.or.
     &   SUSE(iz,is,1)(1:5).eq.'FLOOR'.or.
     &   SUSE(iz,is,1)(1:5).eq.'FURNI'.or.
     &   SUSE(iz,is,1)(1:7).eq.'ITEQUIP'.or.
     &   SUSE(iz,is,1)(1:5).eq.'PARTN'.or.
     &   SUSE(iz,is,1)(1:4).eq.'ROOF'.or.
     &   SUSE(iz,is,1)(1:5).eq.'STRUC'.or.
     &   SUSE(iz,is,1)(1:7).eq.'FIXTURE'.or.
     &   SUSE(iz,is,1)(1:6).eq.'PLANTS')then
        return
      endif

C For DOOR variants.
      if(SUSE(iz,is,1)(1:4).eq.'DOOR'.or.
     &   SUSE(iz,is,1)(1:6).eq.'P-DOOR'.or.
     &   SUSE(iz,is,1)(1:6).eq.'H-DOOR'.or.
     &   SUSE(iz,is,1)(1:6).eq.'H-DOOR')then
        if(SUSE(iz,is,2)(1:6).eq.'CLOSED')then

C Door crack - crack component, assume width = 2mm,
C length = perimeter, on centre line. OK
          call zsurfprm(iz,is,dupedges,perim)
          write(name,'(a,i2.2,a,i3.3)') 'DoCrz',iz,':',is
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=SURCOG(iz,is,3)
          write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),'  ',name,
     &      ' door crack width 2mm len',perim,'m COG @XYZ',
     &      xc,yc,zc
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          NCMP=NCMP+1 ! use logic in MF120I
          CMNAM(NCMP)=name
          CMPASSOC(NCMP,1)=ZNAME(IZ); CMPASSOC(NCMP,2)=SN
          ITPCMP(NCMP)=IVALCM(13); LTPCMP(NCMP)=LVALCM(13)  !crack (120)
          ISDCMP(NCMP)=3; ISDCNN(NCMP)=0
          SUPCMP(NCMP,1)=1.
          SUPCMP(NCMP,2)=2.0/1000.0
          SUPCMP(NCMP,3)=perim
          NWPCMP(NCMP)=0  ! start and end the same
          HCMP(NCMP,1,1)=xc; HCMP(NCMP,1,2)=yc; HCMP(NCMP,1,3)=zc
          call symbol_flow_cmp(name,'cr','r',xc,yc,zc)
          call pausems(100)

C Report on Delta Height with associated .
          write(outs,'(a,f7.3)') 'Delta Zone with surface COG is',
     &      SURCOG(iz,is,3)-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call pausems(100)
          if(ICT(loop).eq.0)then
            write(name,'(a,i2.2,a,i3.3)') 'BD-Cr',iz,':',is
            xw=x3+SURCOG(iz,is,1);yw=y3+SURCOG(iz,is,2)
            zw=z3+SURCOG(iz,is,3)
            write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),'  ',name,
     &        ' Wind boundary crack azim',spazi(iz,is),
     &        ' @XYZ',xw,yw,zw
            call edisp(iuout,outs)
            NNOD=NNOD+1; NDTYP(NNOD)=3; NDFLD(NNOD)=1
            HNOD(NNOD,1)=xw; HNOD(NNOD,2)=yw; HNOD(NNOD,3)=zw
            SUPNOD(NNOD,2)=spazi(iz,is)
            SUPNOD(NNOD,1)=1.0
            TNOD(NNOD)=0.0
            write(NDNAM(NNOD),'(a)') name
            NODASSOC(NNOD,1)=ZNAME(IZ); NODASSOC(NNOD,2)=SN
            call symbol_flow_nd(name,'bw','r',xw,yw,zw)
            call pausems(100)

C Establish connection from this boundary node to the component
C and then to the adjacent zone.
            POS=NNOD; NEG=ICAAS(IZ); CMPIS=NCMP; NCNN=NCNN+1
            if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
              ccol='r'
            else
              ccol='b'
            endif

C Reality checks here...
            NODPS(NCNN)=POS; NODNE(NCNN)=NEG
            HGTPS(NCNN)=0.0   ! boundary always in line with crack
            HGTNE(NCNN)=SURCOG(iz,is,3)-ZCOG(IZ,3)
            NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
            ITPCON(NCNN)=CMPIS ! use the just-defined crack
            write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &        NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
            call edisp(iuout,outs)
            call draw_flow_link(xw,yw,zw,xc,yc,zc,ccol)
            call draw_flow_link(xc,yc,zc,
     &        HNOD(NEG,1),HNOD(NEG,2),HNOD(NEG,3),ccol)
          elseif(ICT(loop).eq.3)then

C Define connection with adjacent zone. Check if linkatpartition is still
C zero. If so mark this connection and update the other surface as well.
            if(linkatpartition(loop).eq.0)then
              POS=ICAAS(IC1(loop))  ! the current zone node is positive
              NEG=ICAAS(IC2(loop)) ! the other zone
              if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
                ccol='r'
              else
                ccol='b'
              endif
              CMPIS=NCMP; NCNN=NCNN+1
              linkatpartition(loop)=1
              icoth=IZSTOCN(IC2(loop),IE2(loop))
              linkatpartition(icoth)=1

C Reality checks here...
              NODPS(NCNN)=POS; NODNE(NCNN)=NEG
              HGTPS(NCNN)= HCMP(NCMP,1,3)-ZCOG(IZ,3)  ! diff comp & zone
              HGTNE(NCNN)= HCMP(NCMP,1,3)-ZCOG(IC2(loop),3)  ! diff comp & other zone
              NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
              ITPCON(NCNN)=CMPIS ! use the just-defined crack
              write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &          NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
              call edisp(iuout,outs)
              call draw_flow_link(ZCOG(IZ,1),ZCOG(IZ,2),ZCOG(IZ,3),
     &          HCMP(NCMP,1,1),HCMP(NCMP,1,2),HCMP(NCMP,1,3),ccol)
              call draw_flow_link(ZCOG(IC2(loop),1),ZCOG(IC2(loop),2),
     &          ZCOG(IC2(loop),3),HCMP(NCMP,1,1),HCMP(NCMP,1,2),
     &          HCMP(NCMP,1,3),ccol)
            endif
          else
            continue  ! << for similar creat another node ....
          endif
        elseif(SUSE(iz,is,2)(1:8).eq.'UNDERCUT')then

C Undercut at door - common orifice,  (type 40) IVALCM(10) MF040I
C Ask width = ~10mm, length = door width at lowest horizontal edge.
C And Cd of ~0.7 (checked against measurements in Gross & Habermann paper).
          write(name,'(a,i2.2,a,i3.3)') 'DoUcz',iz,':',is
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=edgelower
          write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),'  ',name,
     &      ' door undercut width ',XYMAX,'m @XYZ',xc,yc,zc
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          NCMP=NCMP+1
          CMNAM(NCMP)=name
          CMPASSOC(NCMP,1)=ZNAME(IZ); CMPASSOC(NCMP,2)=SN
          ITPCMP(NCMP)=IVALCM(10); LTPCMP(NCMP)=LVALCM(10) ! common orifice (40)
          ISDCMP(NCMP)=3; ISDCNN(NCMP)=0
          SUPCMP(NCMP,1)=1.
          write(prompt,'(3a,f7.5)') 'Undercut area (m^2) for ',
     &      CMNAM(NCMP),' is ',SUPCMP(NCMP,2)
          if(SUPCMP(NCMP,2).lt.0.001)then
            VAL2=0.015                ! initial 15mm high undercut
          else
            VAL2=SUPCMP(NCMP,2)/XYMAX  ! recover hight of undercut
          endif
          CALL EASKR(VAL2,prompt,'height of undercut (m)?',0.,'W',0.,
     &      '-',1.,'undercut (m) in comp',IER,nbhelp)
          SUPCMP(NCMP,2)=VAL2*XYMAX
          if(SUPCMP(NCMP,3).lt.0.001) SUPCMP(NCMP,3)=0.70
          write(prompt,'(3a)')
     &      'The discharge coef of ',CMNAM(NCMP),': '
          VAL2=SUPCMP(NCMP,3)
          CALL EASKR(VAL2,' ',prompt,
     &      0.,'W',0.,'-',1.,'discharge coef',IER,nbhelp)
          SUPCMP(NCMP,3)=VAL2
          NWPCMP(NCMP)=0  ! start and end the same
          HCMP(NCMP,1,1)=xc; HCMP(NCMP,1,2)=yc; HCMP(NCMP,1,3)=zc
          call symbol_flow_cmp(name,'or','r',xc,yc,zc)
          call pausems(100)
          if(ICT(loop).eq.0)then

C Define boundary point ~1m from facade.
            write(name,'(a,i2.2,a,i3.3)') 'BW-Cr',iz,':',is
            xw=x3+SURCOG(iz,is,1);yw=y3+SURCOG(iz,is,2)
            zw=z3+edgelower
            write(outs,'(2a,f6.1,a,3f7.3)') name,
     &        ' Wind boundary undercut azim',spazi(iz,is),
     &        ' @XYZ',xw,yw,zw
            call edisp(iuout,outs)
            NNOD=NNOD+1; NDTYP(NNOD)=3; NDFLD(NNOD)=1
            HNOD(NNOD,1)=xw; HNOD(NNOD,2)=yw; HNOD(NNOD,3)=zw
            SUPNOD(NNOD,2)=spazi(iz,is)
            SUPNOD(NNOD,1)=1.0
            TNOD(NNOD)=0.0
            write(NDNAM(NNOD),'(a)') name
            NODASSOC(NNOD,1)=ZNAME(IZ); NODASSOC(NNOD,2)=SN
            call symbol_flow_nd(name,'bw','r',xw,yw,zw)
            call pausems(100)

C Establish connection from this boundary node to the component
C and then to the adjacent zone.
            POS=NNOD; NEG=ICAAS(IZ); CMPIS=NCMP; NCNN=NCNN+1
            if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
              ccol='r'
            else
              ccol='b'
            endif

C Reality checks here...
            NODPS(NCNN)=POS; NODNE(NCNN)=NEG
            HGTPS(NCNN)=0.0   ! boundary always in line with component
            HGTNE(NCNN)=edgelower-ZCOG(IZ,3)
            NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
            ITPCON(NCNN)=CMPIS ! use the just-defined crack
            write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &        NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
            call edisp(iuout,outs)

C Draw from boundary to component and then from component to the
C associated zone.
            call draw_flow_link(xw,yw,zw,xc,yc,zc,ccol)
            call draw_flow_link(xc,yc,zc,
     &        HNOD(NEG,1),HNOD(NEG,2),HNOD(NEG,3),ccol)
          elseif(ICT(loop).eq.3)then

C Define connection with adjacent zone. Check if linkatpartition is still
C zero. If so mark this connection and update the other surface as well.
            if(linkatpartition(loop).eq.0)then
              POS=ICAAS(IC1(loop))  ! the current zone node is positive
              NEG=ICAAS(IC2(loop)) ! the other zone
              if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
                ccol='r'
              else
                ccol='b'
              endif
              CMPIS=NCMP; NCNN=NCNN+1
              linkatpartition(loop)=1
              icoth=IZSTOCN(IC2(loop),IE2(loop))
              linkatpartition(icoth)=1

C Reality checks here...
              NODPS(NCNN)=POS; NODNE(NCNN)=NEG
              HGTPS(NCNN)= HCMP(NCMP,1,3)-ZCOG(IZ,3)  ! diff comp & zone
              HGTNE(NCNN)= HCMP(NCMP,1,3)-ZCOG(IC2(loop),3)  ! diff comp & other zone
              NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
              ITPCON(NCNN)=CMPIS ! use the just-defined crack
              write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &          NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
              call edisp(iuout,outs)
              call draw_flow_link(ZCOG(IZ,1),ZCOG(IZ,2),ZCOG(IZ,3),
     &          HCMP(NCMP,1,1),HCMP(NCMP,1,2),HCMP(NCMP,1,3),ccol)
              call draw_flow_link(ZCOG(IC2(loop),1),ZCOG(IC2(loop),2),
     &          ZCOG(IC2(loop),3),HCMP(NCMP,1,1),HCMP(NCMP,1,2),
     &          HCMP(NCMP,1,3),ccol)
            endif
          else
            continue  ! << for similar creat another node ....
          endif

        elseif(SUSE(iz,is,2)(1:4).eq.'OPEN')then

C Open door - common orifice component, (type 40) IVALCM(10) MF040I
          write(name,'(a,i2.2,a,i3.3)') 'DoOPz',iz,':',is
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=SURCOG(iz,is,3)
          write(prompt,'(4a,f6.1,a,f6.1,a,f7.3,a,3f6.2)') zsn(1:lnzsn),
     &      '  ',name,' door open width ',XYMAX,' height',ZMAX,
     &      ' area',SNA(iz,is),'m^2 COG @XYZ',xc,yc,zc
          NCMP=NCMP+1 ! use logic in MF110I
          CMNAM(NCMP)=name
          CMPASSOC(NCMP,1)=ZNAME(IZ); CMPASSOC(NCMP,2)=SN
          ITPCMP(NCMP)=IVALCM(10); LTPCMP(NCMP)=LVALCM(10) ! common orifice (40)
          ISDCMP(NCMP)=3; ISDCNN(NCMP)=0
          SUPCMP(NCMP,1)=1.
          VAL2=100.0              ! 100% initial guess of area
          CALL EASKR(VAL2,prompt,'Surface area % for opening?',
     &      0.,'W',1.,'-',100.,'% area',IER,nbhelp)
          SUPCMP(NCMP,2)=SNA(iz,is)*(VAL2*0.001)
          if(SUPCMP(NCMP,3).lt.0.001) SUPCMP(NCMP,3)=0.6  ! initial Cd
          write(prompt,'(3a)')
     &      'The discharge coef of ',CMNAM(NCMP),':'
          VAL2=SUPCMP(NCMP,3)
          CALL EASKR(VAL2,' ',prompt,
     &      0.,'W',0.,'-',1.,'discharge coef',IER,nbhelp)
          SUPCMP(NCMP,3)=VAL2
          NWPCMP(NCMP)=0  ! start and end the same
          HCMP(NCMP,1,1)=xc; HCMP(NCMP,1,2)=yc; HCMP(NCMP,1,3)=zc

C Report on Delta Height with associated surf COG.
          write(outs,'(a,f7.3)') 'Delta Zone with surface COG is',
     &      SURCOG(iz,is,3)-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call symbol_flow_cmp(name,'or','r',xc,yc,zc)
          call pausems(100)
          if(ICT(loop).eq.0)then

C Define boundary point 1m from facade.
            write(name,'(a,i2.2,a,i3.3)') 'BW-Op',iz,':',is
            xw=x3+SURCOG(iz,is,1);yw=y3+SURCOG(iz,is,2)
            zw=z3+SURCOG(iz,is,3)
            write(outs,'(2a,f6.1,a,3f7.3)') name,
     &        ' Wind boundary opening azim',spazi(iz,is),
     &        ' @XYZ',xw,yw,zw
            call edisp(iuout,outs)
            NNOD=NNOD+1; NDTYP(NNOD)=3; NDFLD(NNOD)=1
            HNOD(NNOD,1)=xw; HNOD(NNOD,2)=yw; HNOD(NNOD,3)=zw
            SUPNOD(NNOD,2)=spazi(iz,is)
            SUPNOD(NNOD,1)=1.0
            TNOD(NNOD)=0.0
            write(NDNAM(NNOD),'(a)') name
            NODASSOC(NNOD,1)=ZNAME(IZ); NODASSOC(NNOD,2)=SN
            call symbol_flow_nd(name,'bw','r',xw,yw,zw)
            call pausems(100)

C Establish connection from this boundary node to the orifice
C and then to the adjacent zone.
            POS=NNOD; NEG=ICAAS(IZ); CMPIS=NCMP; NCNN=NCNN+1
            if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
              ccol='r'
            else
              ccol='b'
            endif

C Reality checks here...
            NODPS(NCNN)=POS; NODNE(NCNN)=NEG
            HGTPS(NCNN)=0.0   ! boundary always in line with component
            HGTNE(NCNN)=SURCOG(iz,is,3)-ZCOG(IZ,3)
            NDSCNN(NCNN,1)=0
            NDSCNN(NCNN,2)=0
            ITPCON(NCNN)=CMPIS ! use the just-defined orifice
            write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &        NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
            call edisp(iuout,outs)

C Draw from boundary to component and then from component to the
C associated zone.
            call draw_flow_link(xw,yw,zw,xc,yc,zc,ccol)
            call draw_flow_link(xc,yc,zc,
     &        HNOD(NEG,1),HNOD(NEG,2),HNOD(NEG,3),ccol)
          elseif(ICT(loop).eq.3)then

C Define connection with adjacent zone. Check if linkatpartition is still
C zero. If so mark this connection and update the other surface as well. 
            if(linkatpartition(loop).eq.0)then
              POS=ICAAS(IC1(loop))  ! the current zone node is positive
              NEG=ICAAS(IC2(loop)) ! the other zone
              if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
                ccol='r'
              else
                ccol='b'
              endif
              CMPIS=NCMP; NCNN=NCNN+1
              linkatpartition(loop)=1
              icoth=IZSTOCN(IC2(loop),IE2(loop))
              linkatpartition(icoth)=1

C Reality checks here...
              NODPS(NCNN)=POS; NODNE(NCNN)=NEG
              HGTPS(NCNN)= HCMP(NCMP,1,3)-ZCOG(IZ,3)  ! diff comp & zone
              HGTNE(NCNN)= HCMP(NCMP,1,3)-ZCOG(IC2(loop),3)  ! diff comp & other zone
              NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
              ITPCON(NCNN)=CMPIS ! use the just-defined crack
              write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &          NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
              call edisp(iuout,outs)
              call draw_flow_link(ZCOG(IZ,1),ZCOG(IZ,2),ZCOG(IZ,3),
     &          HCMP(NCMP,1,1),HCMP(NCMP,1,2),HCMP(NCMP,1,3),ccol)
              call draw_flow_link(ZCOG(IC2(loop),1),ZCOG(IC2(loop),2),
     &          ZCOG(IC2(loop),3),HCMP(NCMP,1,1),HCMP(NCMP,1,2),
     &          HCMP(NCMP,1,3),ccol)
            endif
          else
            continue  ! << for similar creat another node ....
          endif

        elseif(SUSE(iz,is,2)(1:5).eq.'BIDIR'.or.
     &         SUSE(iz,is,2)(1:10).eq.'AJAR-BIDIR')then

C Bidirectional door - door component,  IVALCM(14)=130 MF130I
C assume width = door width, height = door height, discharge factor = 0.5.
C If AJAR-BIDIR limit width to 200mm.
C Need to draw component at surf COG but the real base is at bottom
C edge.
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=SURCOG(iz,is,3)
          if(SUSE(iz,is,2)(1:5).eq.'BIDIR')then
            write(name,'(a,i2.2,a,i3.3)') 'DoBiz',iz,':',is
            VAL1=XYMAX; VAL2=ZMAX
            VAL3=ABS(ZCOG(iz,3)-edgelower)  ! base to zone COG
            VAL4=SURCOG(iz,is,3)-ZCOG(IZ,3) ! surf COG to zone COG
            write(outs,'(4a,f6.1,a,f6.1,a,3f7.3,a,f6.3)') zsn(1:lnzsn),
     &      '  ',name,' door bidir w= ',XYMAX,
     &      ' h=',ZMAX,' coef 0.5 base @XYZ',xc,yc,zc,
     &      ' delta to zone COG',VAL3
          elseif(SUSE(iz,is,2)(1:10).eq.'AJAR-BIDIR')then

C Bidirectional door - door component,  IVALCM(14)=130 MF130I
C assume width = 100mm, height = door height, discharge factor = 0.5.
            write(name,'(a,i2.2,a,i3.3)') 'DABiz',iz,':',is
            VAL1=0.1; VAL2=ZMAX   ! 100mm wide
            VAL3=ABS(ZCOG(iz,3)-edgelower)  ! base to zone COG
            VAL4=SURCOG(iz,is,3)-ZCOG(IZ,3) ! surf COG to zone COG
            write(outs,'(5a,f6.1,a,3f7.3,a,f6.3)') zsn(1:lnzsn),
     &      '  ',name,' ajar door bidir w= 100mm',
     &      ' h=',ZMAX,' coef 0.5 base @XYZ',xc,yc,zc,
     &      ' delta to zone COG',VAL3
          endif
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          NCMP=NCMP+1 ! use logic in MF040I
          CMNAM(NCMP)=name
          CMPASSOC(NCMP,1)=ZNAME(IZ); CMPASSOC(NCMP,2)=SN
          ITPCMP(NCMP)=IVALCM(14); LTPCMP(NCMP)=LVALCM(14) ! bi-directional (130)
          ISDCMP(NCMP)=5; ISDCNN(NCMP)=0
          SUPCMP(NCMP,1)=1.
          SUPCMP(NCMP,2)=VAL1
          SUPCMP(NCMP,3)=VAL2
          SUPCMP(NCMP,4)=VAL3
          if(SUPCMP(NCMP,5).lt.0.001) SUPCMP(NCMP,5)=0.5
          write(prompt,'(3a)') 'The discharge coefficient of ',
     &      CMNAM(NCMP),':'
          VAL2=SUPCMP(NCMP,5)
          CALL EASKR(VAL2,' ',prompt,0.,'W',0.,'-',10.,
     &      'discharge coef in comp',IER,nbhelp)
          SUPCMP(NCMP,5)=VAL2
          NWPCMP(NCMP)=0  ! start and end the same
          HCMP(NCMP,1,1)=xc; HCMP(NCMP,1,2)=yc; HCMP(NCMP,1,3)=zc
          call symbol_flow_cmp(name,'bi','r',xc,yc,zc)
          call pausems(100)

C Report on Delta Height between surface COG and zone COG.
          write(outs,'(a,f7.3)') 'Delta Zone with base is',VAL4
          call edisp(iuout,outs)
          call pausems(100)
          if(ICT(loop).eq.0)then
            write(name,'(a,i2.2,a,i3.3)') 'BW-Bi',iz,':',is
            xw=x3+SURCOG(iz,is,1);yw=y3+SURCOG(iz,is,2)
            zw=z3+SURCOG(iz,is,3)
            write(outs,'(2a,f6.1,a,3f7.3)') name,
     &        ' Wind boundary bidir azim',spazi(iz,is),
     &        ' @XYZ',xw,yw,zw
            call edisp(iuout,outs)
            NNOD=NNOD+1; NDTYP(NNOD)=3; NDFLD(NNOD)=1
            HNOD(NNOD,1)=xw; HNOD(NNOD,2)=yw; HNOD(NNOD,3)=zw
            SUPNOD(NNOD,2)=spazi(iz,is)
            SUPNOD(NNOD,1)=1.0
            TNOD(NNOD)=0.0
            write(NDNAM(NNOD),'(a)') name
            NODASSOC(NNOD,1)=ZNAME(IZ); NODASSOC(NNOD,2)=SN
            call symbol_flow_nd(name,'bw','r',xw,yw,zw)
            call pausems(100)

C Establish connection from this boundary node to the bi-directional
C and then to the adjacent zone.
            POS=NNOD; NEG=ICAAS(IZ); CMPIS=NCMP; NCNN=NCNN+1
            if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
              ccol='r'
            else
              ccol='b'
            endif

C Reality checks here...
            NODPS(NCNN)=POS; NODNE(NCNN)=NEG
            HGTPS(NCNN)=0.0   ! boundary always in line with component
            HGTNE(NCNN)=VAL4   ! bi-dir assume surf COG to zone COG
            NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
            ITPCON(NCNN)=CMPIS ! use the just-defined orifice
            write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &        NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
            call edisp(iuout,outs)
            call draw_flow_link(xw,yw,zw,xc,yc,zc,ccol)
            call draw_flow_link(xc,yc,zc,
     &        HNOD(NEG,1),HNOD(NEG,2),HNOD(NEG,3),ccol)
          elseif(ICT(loop).eq.3)then

C Define connection with adjacent zone. Check if linkatpartition is still
C zero. If so mark this connection and update the other surface as well. 
            if(linkatpartition(loop).eq.0)then
              POS=ICAAS(IC1(loop))  ! the current zone node is positive
              NEG=ICAAS(IC2(loop)) ! the other zone
              if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
                ccol='r'
              else
                ccol='b'
              endif
              CMPIS=NCMP; NCNN=NCNN+1
              linkatpartition(loop)=1
              icoth=IZSTOCN(IC2(loop),IE2(loop))
              linkatpartition(icoth)=1

C Reality checks here...
              NODPS(NCNN)=POS; NODNE(NCNN)=NEG
              HGTPS(NCNN)= HCMP(NCMP,1,3)-ZCOG(IZ,3)  ! diff comp & zone
              HGTNE(NCNN)= HCMP(NCMP,1,3)-ZCOG(IC2(loop),3)  ! diff comp & other zone
              NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
              ITPCON(NCNN)=CMPIS ! use the just-defined crack
              write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &          NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
              call edisp(iuout,outs)
              call draw_flow_link(ZCOG(IZ,1),ZCOG(IZ,2),ZCOG(IZ,3),
     &          HCMP(NCMP,1,1),HCMP(NCMP,1,2),HCMP(NCMP,1,3),ccol)
              call draw_flow_link(ZCOG(IC2(loop),1),ZCOG(IC2(loop),2),
     &          ZCOG(IC2(loop),3),HCMP(NCMP,1,1),HCMP(NCMP,1,2),
     &          HCMP(NCMP,1,3),ccol)
            endif
          else
            continue  ! << for similar creat another node ....
          endif
        endif
      elseif(SUSE(iz,is,1)(1:5).eq.'FRAME'.or.
     &       SUSE(iz,is,1)(1:7).eq.'F-FRAME')then

C For FRAME variants.
        if(SUSE(iz,is,2)(1:6).eq.'CLOSED'.or.
     &     SUSE(iz,is,2)(1:5).eq.'CRACK')then

C Frame crack - crack component, assume width = 2mm,
C length = frame inside perimeter, on centre line. OK
          call zsurfprm(iz,is,dupedges,perim)
          write(name,'(a,i2.2,a,i3.3)') 'FrCrz',iz,':',is
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=SURCOG(iz,is,3)
          write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),'  ',name,
     &      ' frame crack width 2mm len',perim,'m COG @XYZ',
     &      xc,yc,zc
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          NCMP=NCMP+1 ! use logic in MF120I
          CMNAM(NCMP)=name
          CMPASSOC(NCMP,1)=ZNAME(IZ); CMPASSOC(NCMP,2)=SN
          ITPCMP(NCMP)=IVALCM(13); LTPCMP(NCMP)=LVALCM(13)  !crack
          ISDCMP(NCMP)=3; ISDCNN(NCMP)=0
          SUPCMP(NCMP,1)=1.
          SUPCMP(NCMP,2)=2.0/1000.0
          SUPCMP(NCMP,3)=perim
          NWPCMP(NCMP)=0  ! start and end the same
          HCMP(NCMP,1,1)=xc; HCMP(NCMP,1,2)=yc; HCMP(NCMP,1,3)=zc
          call symbol_flow_cmp(name,'cr','r',xc,yc,zc)
          call pausems(100)

C Report on Delta Height with associated surf COG.
          write(outs,'(a,f7.3)') 'Delta Zone with surface COG is',
     &      SURCOG(iz,is,3)-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call pausems(100)
          if(ICT(loop).eq.0)then
            write(name,'(a,i2.2,a,i3.3)') 'BW-Cr',iz,':',is
            xw=x3+SURCOG(iz,is,1);yw=y3+SURCOG(iz,is,2)
            zw=z3+SURCOG(iz,is,3)
            write(outs,'(2a,f6.1,a,3f7.3)') name,
     &        ' Wind boundary crack azim',spazi(iz,is),
     &        ' @XYZ',xw,yw,zw
            call edisp(iuout,outs)
            NNOD=NNOD+1; NDTYP(NNOD)=3; NDFLD(NNOD)=1
            HNOD(NNOD,1)=xw; HNOD(NNOD,2)=yw; HNOD(NNOD,3)=zw
            SUPNOD(NNOD,2)=spazi(iz,is)
            SUPNOD(NNOD,1)=1.0
            TNOD(NNOD)=0.0
            write(NDNAM(NNOD),'(a)') name
            NODASSOC(NNOD,1)=ZNAME(IZ); NODASSOC(NNOD,2)=SN
            call symbol_flow_nd(name,'bw','r',xw,yw,zw)
            call pausems(100)

C Establish connection from this boundary node to the component
C and then to the adjacent zone.
            POS=NNOD; NEG=ICAAS(IZ); CMPIS=NCMP; NCNN=NCNN+1
            if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
              ccol='r'
            else
              ccol='b'
            endif

C Reality checks here...
            NODPS(NCNN)=POS; NODNE(NCNN)=NEG
            HGTPS(NCNN)=0.0   ! boundary always in line with crack
            HGTNE(NCNN)=SURCOG(iz,is,3)-ZCOG(IZ,3)
            NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
            ITPCON(NCNN)=CMPIS ! use the just-defined crack
            write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &        NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
            call edisp(iuout,outs)
            call draw_flow_link(xw,yw,zw,xc,yc,zc,ccol)
            call draw_flow_link(xc,yc,zc,
     &        HNOD(NEG,1),HNOD(NEG,2),HNOD(NEG,3),ccol)
          elseif(ICT(loop).eq.3)then

C Define connection with adjacent zone. Check if linkatpartition is still
C zero. If so mark this connection and update the other surface as well.
            if(linkatpartition(loop).eq.0)then
              POS=ICAAS(IC1(loop))  ! the current zone node is positive
              NEG=ICAAS(IC2(loop)) ! the other zone
              if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
                ccol='r'
              else
                ccol='b'
              endif
              CMPIS=NCMP; NCNN=NCNN+1
              linkatpartition(loop)=1
              icoth=IZSTOCN(IC2(loop),IE2(loop))
              linkatpartition(icoth)=1

C Reality checks here...
              NODPS(NCNN)=POS; NODNE(NCNN)=NEG
              HGTPS(NCNN)= HCMP(NCMP,1,3)-ZCOG(IZ,3)  ! diff comp & zone
              HGTNE(NCNN)= HCMP(NCMP,1,3)-ZCOG(IC2(loop),3)  ! diff comp & other zone
              NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
              ITPCON(NCNN)=CMPIS ! use the just-defined crack
              write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &          NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
              call edisp(iuout,outs)
              call draw_flow_link(ZCOG(IZ,1),ZCOG(IZ,2),ZCOG(IZ,3),
     &          HCMP(NCMP,1,1),HCMP(NCMP,1,2),HCMP(NCMP,1,3),ccol)
              call draw_flow_link(ZCOG(IC2(loop),1),ZCOG(IC2(loop),2),
     &          ZCOG(IC2(loop),3),HCMP(NCMP,1,1),HCMP(NCMP,1,2),
     &          HCMP(NCMP,1,3),ccol)
            endif
          else
            continue  ! << for similar creat another node ....
          endif

        elseif(SUSE(iz,is,2)(1:4).eq.'VENT')then

C Frame trickle vent - common orifice, assume width = 0.6* frame inside width,
C but not more than 0.7m and not less than 0.4m unless the frame is less wide,
C height = 15mm, at highest horizontal edge.
          write(name,'(a,i2.2,a,i3.3)') 'FrVez',iz,':',is
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=edgeupper
          trwidth=XYMAX*0.6
          if(trwidth.gt.0.7) trwidth=0.7
          if(trwidth.lt.0.4.and.XYMAX.lt.0.4) trwidth=XYMAX
          write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),'  ',name,
     &      ' frame trickle vent width ',trwidth,' height 15mm @XYZ',
     &      xc,yc,zc
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          NCMP=NCMP+1
          CMNAM(NCMP)=name
          CMPASSOC(NCMP,1)=ZNAME(IZ); CMPASSOC(NCMP,2)=SN
          ITPCMP(NCMP)=IVALCM(10); LTPCMP(NCMP)=LVALCM(10) ! common orifice (40)
          ISDCMP(NCMP)=3; ISDCNN(NCMP)=0
          SUPCMP(NCMP,1)=1.
          SUPCMP(NCMP,2)=trwidth*0.015  ! area assumes 15mm high
          if(SUPCMP(NCMP,3).lt.0.001) SUPCMP(NCMP,3)=0.34
          write(prompt,'(3a)')
     &      'The discharge coef of ',CMNAM(NCMP),':'
          VAL2=SUPCMP(NCMP,3)
          CALL EASKR(VAL2,' ',prompt,
     &      0.,'W',0.,'-',1.,'discharge coef in comp',IER,nbhelp)
          SUPCMP(NCMP,3)=VAL2
          NWPCMP(NCMP)=0  ! start and end the same
          HCMP(NCMP,1,1)=xc; HCMP(NCMP,1,2)=yc; HCMP(NCMP,1,3)=zc
          call symbol_flow_cmp(name,'or','r',xc,yc,zc)
          call pausems(100)

C Report on Delta Height with associated surf upper.
          write(outs,'(a,f7.3)') 'Delta Zone with surf top is',
     &      edgeupper-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call pausems(100)
          if(ICT(loop).eq.0)then
            write(name,'(a,i2.2,a,i3.3)') 'BW-Ve',iz,':',is
            xw=x3+SURCOG(iz,is,1);yw=y3+SURCOG(iz,is,2)
            zw=z3+edgeupper
            write(outs,'(2a,f6.1,a,3f7.3)') name,
     &        ' Wind boundary azim',spazi(iz,is),
     &        ' @XYZ',xw,yw,zw
            call edisp(iuout,outs)
            NNOD=NNOD+1; NDTYP(NNOD)=3; NDFLD(NNOD)=1
            HNOD(NNOD,1)=xw; HNOD(NNOD,2)=yw; HNOD(NNOD,3)=zw
            SUPNOD(NNOD,2)=spazi(iz,is)
            SUPNOD(NNOD,1)=1.0
            TNOD(NNOD)=0.0
            write(NDNAM(NNOD),'(a)') name
            NODASSOC(NNOD,1)=ZNAME(IZ); NODASSOC(NNOD,2)=SN
            call symbol_flow_nd(name,'bw','r',xw,yw,zw)
            call pausems(100)

C Establish connection from this boundary node to the vent
C and then to the adjacent zone.
            POS=NNOD; NEG=ICAAS(IZ); CMPIS=NCMP; NCNN=NCNN+1
            if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
              ccol='r'
            else
              ccol='b'
            endif

C Reality checks here...
            NODPS(NCNN)=POS; NODNE(NCNN)=NEG
            HGTPS(NCNN)=0.0   ! boundary always in line with crack
            HGTNE(NCNN)=edgeupper-ZCOG(IZ,3)
            NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
            ITPCON(NCNN)=CMPIS ! use the just-defined crack
            write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &        NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
            call edisp(iuout,outs)
            call draw_flow_link(xw,yw,zw,xc,yc,zc,ccol)
            call draw_flow_link(xc,yc,zc,
     &        HNOD(NEG,1),HNOD(NEG,2),HNOD(NEG,3),ccol)
          elseif(ICT(loop).eq.3)then

C Define connection with adjacent zone. Check if linkatpartition is still
C zero. If so mark this connection and update the other surface as well.
            if(linkatpartition(loop).eq.0)then
              POS=ICAAS(IC1(loop))  ! the current zone node is positive
              NEG=ICAAS(IC2(loop)) ! the other zone
              if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
                ccol='r'
              else
                ccol='b'
              endif
              CMPIS=NCMP; NCNN=NCNN+1
              linkatpartition(loop)=1
              icoth=IZSTOCN(IC2(loop),IE2(loop))
              linkatpartition(icoth)=1

C Reality checks here...
              NODPS(NCNN)=POS; NODNE(NCNN)=NEG
              HGTPS(NCNN)= HCMP(NCMP,1,3)-ZCOG(IZ,3)  ! diff comp & zone
              HGTNE(NCNN)= HCMP(NCMP,1,3)-ZCOG(IC2(loop),3)  ! diff comp & other zone
              NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
              ITPCON(NCNN)=CMPIS ! use the just-defined crack
              write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &          NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
              call edisp(iuout,outs)
              call draw_flow_link(ZCOG(IZ,1),ZCOG(IZ,2),ZCOG(IZ,3),
     &          HCMP(NCMP,1,1),HCMP(NCMP,1,2),HCMP(NCMP,1,3),ccol)
              call draw_flow_link(ZCOG(IC2(loop),1),ZCOG(IC2(loop),2),
     &          ZCOG(IC2(loop),3),HCMP(NCMP,1,1),HCMP(NCMP,1,2),
     &          HCMP(NCMP,1,3),ccol)
            endif
          else
            continue  ! << for similar creat another node ....
          endif

        endif
      elseif(SUSE(iz,is,1)(1:5).eq.'GRILL')then

C For GRILLs.
        if(SUSE(iz,is,2)(1:6).eq.'CLOSED'.or.
     &     SUSE(iz,is,2)(1:5).eq.'CRACK')then

C Grill crack - crack component, assume width = 2mm,
C length = frame inside perimeter, on centre line.
          call zsurfprm(iz,is,dupedges,perim)
          write(name,'(a,i2.2,a,i3.3)') 'GrCRz',iz,':',is
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=SURCOG(iz,is,3)
          write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),'  ',name,
     &      ' grill crack width 2mm length',perim,'m COG @XYZ',
     &      xc,yc,zc
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          NCMP=NCMP+1 ! use logic in MF120I
          CMNAM(NCMP)=name
          CMPASSOC(NCMP,1)=ZNAME(IZ); CMPASSOC(NCMP,2)=SN
          ITPCMP(NCMP)=IVALCM(13); LTPCMP(NCMP)=LVALCM(13)  !crack (120)
          ISDCMP(NCMP)=3; ISDCNN(NCMP)=0
          SUPCMP(NCMP,1)=1.
          SUPCMP(NCMP,2)=2.0/1000.0
          SUPCMP(NCMP,3)=perim
          NWPCMP(NCMP)=0  ! start and end the same
          HCMP(NCMP,1,1)=xc; HCMP(NCMP,1,2)=yc; HCMP(NCMP,1,3)=zc
          call symbol_flow_cmp(name,'cr','r',xc,yc,zc)
          call pausems(100)

C Report on Delta Height with associated surf COG.
          write(outs,'(a,f7.3)') 'Delta Zone with surface COG is',
     &      SURCOG(iz,is,3)-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call pausems(100)
          if(ICT(loop).eq.0)then
            write(name,'(a,i2.2,a,i3.3)') 'BW-Cr',iz,':',is
            xw=x3+SURCOG(iz,is,1);yw=y3+SURCOG(iz,is,2)
            zw=z3+SURCOG(iz,is,3)
            write(outs,'(2a,f6.1,a,3f7.3)') name,
     &        ' Wind boundary crack azim',spazi(iz,is),
     &        ' @XYZ',xw,yw,zw
            call edisp(iuout,outs)
            NNOD=NNOD+1; NDTYP(NNOD)=3; NDFLD(NNOD)=1
            HNOD(NNOD,1)=xw; HNOD(NNOD,2)=yw; HNOD(NNOD,3)=zw
            SUPNOD(NNOD,2)=spazi(iz,is)
            SUPNOD(NNOD,1)=1.0
            TNOD(NNOD)=0.0
            write(NDNAM(NNOD),'(a)') name
            NODASSOC(NNOD,1)=ZNAME(IZ); NODASSOC(NNOD,2)=SN
            call symbol_flow_nd(name,'bw','r',xw,yw,zw)
            call pausems(100)

C Establish connection from this boundary node to the crack
C and then to the adjacent zone.
            POS=NNOD; NEG=ICAAS(IZ); CMPIS=NCMP; NCNN=NCNN+1
            if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
              ccol='r'
            else
              ccol='b'
            endif

C Reality checks here...
            NODPS(NCNN)=POS; NODNE(NCNN)=NEG
            HGTPS(NCNN)=0.0   ! boundary always in line with crack
            HGTNE(NCNN)=SURCOG(iz,is,3)-ZCOG(IZ,3)
            NDSCNN(NCNN,1)=0
            NDSCNN(NCNN,2)=0
            ITPCON(NCNN)=CMPIS ! use the just-defined crack
            write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &        NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
            call edisp(iuout,outs)
            call draw_flow_link(xw,yw,zw,xc,yc,zc,ccol)
            call draw_flow_link(xc,yc,zc,
     &        HNOD(NEG,1),HNOD(NEG,2),HNOD(NEG,3),ccol)
          elseif(ICT(loop).eq.3)then

C Define connection with adjacent zone. Check if linkatpartition is still
C zero. If so mark this connection and update the other surface as well.
            if(linkatpartition(loop).eq.0)then
              POS=ICAAS(IC1(loop))  ! the current zone node is positive
              NEG=ICAAS(IC2(loop)) ! the other zone
              if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
                ccol='r'
              else
                ccol='b'
              endif
              CMPIS=NCMP; NCNN=NCNN+1
              linkatpartition(loop)=1
              icoth=IZSTOCN(IC2(loop),IE2(loop))
              linkatpartition(icoth)=1

C Reality checks here...
              NODPS(NCNN)=POS; NODNE(NCNN)=NEG
              HGTPS(NCNN)= HCMP(NCMP,1,3)-ZCOG(IZ,3)  ! diff comp & zone
              HGTNE(NCNN)= HCMP(NCMP,1,3)-ZCOG(IC2(loop),3)  ! diff comp & other zone
              NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
              ITPCON(NCNN)=CMPIS ! use the just-defined crack
              write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &          NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
              call edisp(iuout,outs)
              call draw_flow_link(ZCOG(IZ,1),ZCOG(IZ,2),ZCOG(IZ,3),
     &          HCMP(NCMP,1,1),HCMP(NCMP,1,2),HCMP(NCMP,1,3),ccol)
              call draw_flow_link(ZCOG(IC2(loop),1),ZCOG(IC2(loop),2),
     &          ZCOG(IC2(loop),3),HCMP(NCMP,1,1),HCMP(NCMP,1,2),
     &          HCMP(NCMP,1,3),ccol)
            endif
          else
            continue  ! << for similar creat another node ....
          endif

        elseif(SUSE(iz,is,2)(1:5).eq.'INLET')then

C A grill inlet could be a volume flow rate component or a conduit.
          write(outs,'(3a)')'Grill inlet for ',zsn(1:lnzsn),' can be'
 242      CALL EASKAB(outs,'a','constant volume flow (30)',
     &      'or conduit (210)',iacc,nbhelp)
          if(iacc.eq.1)then

C     SUPCMP(ICMP,1) - fluid type (1=air, 2=water)
C     SUPCMP(ICMP,2) - volume flow rate m^3/s
            write(name,'(a,i2.2,a,i3.3)') 'GrMEI',iz,':',is
            xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
            zc=SURCOG(iz,is,3)
            write(outs,'(4a,3f7.3)') zsn(1:lnzsn),'  ',name,
     &        ' grill mechanical extract COG @XYZ',xc,yc,zc
            call edisp(iuout,outs)
            NCMP=NCMP+1 ! use logic in MF030I
            CMNAM(NCMP)=name
            CMPASSOC(NCMP,1)=ZNAME(IZ); CMPASSOC(NCMP,2)=SN
            ITPCMP(NCMP)=IVALCM(8); LTPCMP(NCMP)=LVALCM(8) ! constant vol (30)
            ISDCMP(NCMP)=2; ISDCNN(NCMP)=0
            NWPCMP(NCMP)=0  ! start and end the same
            HCMP(NCMP,1,1)=xc; HCMP(NCMP,1,2)=yc; HCMP(NCMP,1,3)=zc
            SUPCMP(NCMP,1)=1.

C Offer shortcut definition (same as in mfprb2.F MF030I.
            write(prompt,'(3a,F9.6,a)') 'The flow of ',
     &        CMNAM(NCMP),' (m^3/s) is ',SUPCMP(NCMP,2),'.'
            write(prompt,'(3a)') 'Derive flow through ',CMNAM(NCMP),
     &        ' (m^3/s) via:'
            CALL EASKABCD(prompt,' ','direct edit','l/s/person',
     &        'm^3/hr','ach (use volume of node)',IW,nbhelp)
            if(iw.eq.1)then

C Flow m^3/s (zero flow allowed).
              VAL=SUPCMP(NCMP,2)
              write(prompt,'(3a)') 'Flow through ',CMNAM(NCMP),
     &           ' (m^3/s)?'
  220         CALL EASKR(VAL,' ',prompt,
     &          0.,'W',0.,'-',1.,'flow in comp',IER,nbhelp)
              IF(IER.NE.0) GOTO 220
              SUPCMP(NCMP,2)=VAL
            elseif(iw.eq.2)then

C Derive from litres/second/person and nb of people.
              write(prompt,'(3a,F8.5,a)') 'The flow of ',
     &          CMNAM(NCMP),' (m^3/s) is ',SUPCMP(NCMP,2),'.'
              hold = '  6.0    1.0 '
  243         CALL EASKS(hold,prompt,
     &          'litres/sec/person & number people?',
     &          32,' 6.0 1 ','l/s/p',IER,nbhelp)
              K=0
              CALL EGETWR(hold,K,VAL2,0.1,99.,'W','l/s/p',IER)
              CALL EGETWR(hold,K,VAL3,0.1,999.,'W','num people',IER)
              IF(IER.NE.0) GOTO 243

              VAL= (val2/1000.0) * val3
              write(prompt,'(3a,F8.5,a)') 'The flow through ',
     &          CMNAM(NCMP),' (m^3/s) is ',VAL,'.'
              call easkok(prompt,'Is this OK?',OK,nbhelp)
              if (.NOT.OK) goto 242
              SUPCMP(NCMP,2)=VAL
            elseif(iw.eq.3)then

C Convert any existing value to m3/hr
              write(prompt,'(3a,F8.5,a)') 'The flow of ',
     &          CMNAM(NCMP),' (m^3/s) is ',SUPCMP(NCMP,2),'.'
              VAL2=SUPCMP(NCMP,2)*3600.0
              write(prompt,'(3a)') 'Flow through ',CMNAM(NCMP),
     &          ' (m^3/hr)?'
  221         CALL EASKR(VAL2,' ',prompt,
     &          0.,'W',0.,'-',10.,'flow m3/hr in comp',IER,nbhelp)
              IF(IER.NE.0) GOTO 221
              if(val2.gt.0.0) VAL = (val2/3600.0)
              write(prompt,'(3a,F8.5,a)') 'The flow through ',
     &          CMNAM(NCMP),' (m^3/s) is ',VAL,'.'
              call easkok(prompt,'Is this OK?',OK,nbhelp)
              if (.NOT.OK) goto 242
              SUPCMP(NCMP,2)=VAL
            elseif(iw.eq.4)then

C Based on ach against the volume of a particular node.
              write(prompt,'(3a,F8.5,a)') 'The flow via ',
     &          CMNAM(NCMP),' (m^3/s) is ',SUPCMP(NCMP,2),'.'
              call edisp(iuout,prompt)
              call ASKRNOD(' volume from which node','-',IV,IER)
              write(prompt,*)'Associated node volume is ',SUPNOD(IV,2)
              call edisp(iuout,prompt)
              val2=1.0
              write(prompt,'(3a)') 'Air changes for ',CMNAM(NCMP),
     &          ' (ach)?'
  222         CALL EASKR(VAL2,' ',prompt,
     &          0.,'W',0.,'-',1.,'flow ach in comp',IER,nbhelp)
              IF(IER.NE.0) GOTO 222
              VAL = (SUPNOD(IV,2) * VAL2)/3600.0
              write(prompt,'(3a,F8.5,a)') 'The flow through ',
     &          CMNAM(NCMP),' (m^3/s) is ',VAL,'.'
              call easkok(prompt,'Is this OK?',OK,nbhelp)
              if (.NOT.OK) goto 242
              SUPCMP(NCMP,2)=VAL
            endif

            write(outs,'(2a,f6.1,a,3f7.3)') name,
     &        ' grill inlet volume flow',SUPCMP(NCMP,2),
     &        ' COG @XYZ',xc,yc,zc
            call edisp(iuout,outs)
            call ECLOSE(SUPCMP(NCMP,1),2.0,0.0001,iswater)
            if(.NOT.iswater)then
              call symbol_flow_cmp(name,'bx','r',xc,yc,zc)
            else
              call symbol_flow_cmp(name,'bx','b',xc,yc,zc)
            endif
            call pausems(30)

C << deal with start/end points likely being different >>
          else

C Grill inlet - a conduit (type 210) IVALCM(15), as in MF210I
C     SUPCMP(ICMP,1) - fluid type (1=air, 2=water)
C     SUPCMP(ICMP,2) - conduit hydraulic diameter (m)
C     SUPCMP(ICMP,3) - cross-sectional area (m^2)
C     SUPCMP(ICMP,4) - conduit length (m)  ?? get from distance ??
C     SUPCMP(ICMP,5) - absolute wall roughness (m)
C     SUPCMP(ICMP,6) - sum of local dynamic loss factors (-)
            call zsurfprm(iz,is,dupedges,perim)
            write(name,'(a,i2.2,a,i3.3)') 'GrINz',iz,':',is
            xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
            zc=SURCOG(iz,is,3)
            write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),'  ',name,
     &        ' grill inlet hydr diam',
     &        4.0*SNA(iz,is)/perim,'m COG @XYZ',xc,yc,zc
            call edisp(iuout,' ')
            call edisp(iuout,outs)
            NCMP=NCMP+1 ! use logic in MF110I
            CMNAM(NCMP)=name
            CMPASSOC(NCMP,1)=ZNAME(IZ); CMPASSOC(NCMP,2)=SN
            ITPCMP(NCMP)=IVALCM(15); LTPCMP(NCMP)=LVALCM(15) ! conduit (210)
            ISDCMP(NCMP)=6; ISDCNN(NCMP)=0
            SUPCMP(NCMP,1)=1.
            SUPCMP(NCMP,2)=4.0*SNA(iz,is)/perim
            SUPCMP(NCMP,3)=SNA(iz,is)
            SUPCMP(NCMP,4)=1.0
            SUPCMP(NCMP,5)=1.00000E-04
            SUPCMP(NCMP,6)=6.00000E-03
            NWPCMP(NCMP)=0  ! start and end the same
            HCMP(NCMP,1,1)=xc; HCMP(NCMP,1,2)=yc; HCMP(NCMP,1,3)=zc
            call ECLOSE(SUPCMP(NCMP,1),2.0,0.0001,iswater)
            if(.NOT.iswater)then
              call symbol_flow_cmp(name,'cd','r',xc,yc,zc)
            else
              call symbol_flow_cmp(name,'cd','b',xc,yc,zc)
            endif
            call pausems(50)
          endif

C Report on Delta Height with associated surf COG.
          write(outs,'(a,f7.3)') 'Delta Zone with surface COG is',
     &      SURCOG(iz,is,3)-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call pausems(50)
          if(ICT(loop).eq.0)then
            write(name,'(a,i2.2,a,i3.3)') 'BW-IN',iz,':',is
            xw=x3+SURCOG(iz,is,1);yw=y3+SURCOG(iz,is,2)
            zw=z3+SURCOG(iz,is,3)
            write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),'  ',name,
     &        ' Wind boundary inlet azim',spazi(iz,is),
     &        ' @XYZ',xw,yw,zw
            call edisp(iuout,outs)
            NNOD=NNOD+1; NDTYP(NNOD)=3; NDFLD(NNOD)=1
            HNOD(NNOD,1)=xw; HNOD(NNOD,2)=yw; HNOD(NNOD,3)=zw
            SUPNOD(NNOD,2)=spazi(iz,is)
            SUPNOD(NNOD,1)=1.0
            TNOD(NNOD)=0.0
            write(NDNAM(NNOD),'(a)') name
            NODASSOC(NNOD,1)=ZNAME(IZ); NODASSOC(NNOD,2)=SN
            call symbol_flow_nd(name,'bw','r',xw,yw,zw)
            call pausems(50)

C Establish connection from this boundary node to the orifice
C and then to the adjacent zone.
            POS=NNOD; NEG=ICAAS(IZ); CMPIS=NCMP; NCNN=NCNN+1
            if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
              ccol='r'
            else
              ccol='b'
            endif

C Reality checks here...
            NODPS(NCNN)=POS; NODNE(NCNN)=NEG
            HGTPS(NCNN)=0.0   ! boundary always in line with crack
            HGTNE(NCNN)=SURCOG(iz,is,3)-ZCOG(IZ,3)
            NDSCNN(NCNN,1)=0
            NDSCNN(NCNN,2)=0
            ITPCON(NCNN)=CMPIS ! use the just-defined crack
            write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &        NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
            call edisp(iuout,outs)
            call draw_flow_link(xw,yw,zw,xc,yc,zc,ccol)
            call draw_flow_link(xc,yc,zc,
     &        HNOD(NEG,1),HNOD(NEG,2),HNOD(NEG,3),ccol)
          elseif(ICT(loop).eq.3)then

C Define connection with adjacent zone. As in INLET the other zone is
C positive and this zone is negative. Check if linkatpartition is still
C zero. If so mark this connection and update the other surface as well.
            if(linkatpartition(loop).eq.0)then
              NEG=ICAAS(IC1(loop))  ! the current zone node
              POS=ICAAS(IC2(loop))  ! the other zone is positive
              if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
                ccol='r'
              else
                ccol='b'
              endif
              CMPIS=NCMP; NCNN=NCNN+1
              linkatpartition(loop)=1
              icoth=IZSTOCN(IC2(loop),IE2(loop))
              linkatpartition(icoth)=1

C Reality checks here...
              NODPS(NCNN)=POS; NODNE(NCNN)=NEG
              HGTPS(NCNN)= HCMP(NCMP,1,3)-ZCOG(IZ,3)  ! diff comp & zone
              HGTNE(NCNN)= HCMP(NCMP,1,3)-ZCOG(IC2(loop),3)  ! diff comp & other zone
              NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
              ITPCON(NCNN)=CMPIS ! use the just-defined crack
              write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &          NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
              call edisp(iuout,outs)
              call draw_flow_link(ZCOG(IZ,1),ZCOG(IZ,2),ZCOG(IZ,3),
     &          HCMP(NCMP,1,1),HCMP(NCMP,1,2),HCMP(NCMP,1,3),ccol)
              call draw_flow_link(ZCOG(IC2(loop),1),ZCOG(IC2(loop),2),
     &          ZCOG(IC2(loop),3),HCMP(NCMP,1,1),HCMP(NCMP,1,2),
     &          HCMP(NCMP,1,3),ccol)
            endif
          else
            continue  ! << for similar creat another node ....
          endif

        elseif(SUSE(iz,is,2)(1:7).eq.'EXTRACT')then

C Grill extract - constant vol flow comp. (type 30) IVALCM(8), as in MF030I
          write(name,'(a,i2.2,a,i3.3)') 'GrEXz',iz,':',is
          call zsurfprm(iz,is,dupedges,perim)
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=SURCOG(iz,is,3)
          write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),'  ',name,
     &      ' grill extract hydr diam',
     &      4.0*SNA(iz,is)/perim,'m COG @XYZ',xc,yc,zc
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          NCMP=NCMP+1 ! use logic in MF110I
          CMNAM(NCMP)=name
          CMPASSOC(NCMP,1)=ZNAME(IZ); CMPASSOC(NCMP,2)=SN
          ITPCMP(NCMP)=IVALCM(8); LTPCMP(NCMP)=LVALCM(8) ! constant vol (30)
          ISDCMP(NCMP)=2; ISDCNN(NCMP)=0
          NWPCMP(NCMP)=0  ! start and end the same
          HCMP(NCMP,1,1)=xc; HCMP(NCMP,1,2)=yc; HCMP(NCMP,1,3)=zc
          SUPCMP(NCMP,1)=1.

C Offer shortcut definition.
          write(prompt,'(3a,F9.6,a)') 'The flow via ',
     &      CMNAM(NCMP),' (m^3/s) is ',SUPCMP(NCMP,2),'.'
          call edisp(iuout,prompt)
  42      write(prompt,'(3a)') 'Derive flow through ',CMNAM(NCMP),
     &      ' (m^3/s) via:'
          CALL EASKABCD(prompt,' ','direct edit','l/s/person',
     &      'm^3/hr','ach (zone volume)',IW,nbhelp)
          if(iw.eq.1)then

C Flow m^3/s (zero flow allowed).
            VAL=SUPCMP(NCMP,2)
            write(prompt,'(3a)') 'Flow through ',CMNAM(NCMP),
     &        ' (m^3/s)?'
   20       CALL EASKR(VAL,' ',prompt,0.,'W',0.,'-',1.,
     &        'flow in comp',IER,nbhelp)
            IF(IER.NE.0) GOTO 20
            SUPCMP(NCMP,2)=VAL
          elseif(iw.eq.2)then

C Derive from litres/second/person and nb of people.
            write(prompt,'(3a,F8.5,a)') ' The current flow of ',
     &        CMNAM(NCMP),' (m^3/s) is ',SUPCMP(NCMP,2),'.'
            call edisp(iuout,prompt)
            hold = '  6.0    1.0 '
  43        CALL EASKS(hold,' litres/sec/person & number of people?',
     &          '(help for synopsis)',32,' 6.0 1 ','l/s/p',IER,nbhelp)
            K=0
            CALL EGETWR(hold,K,VAL2,0.1,99.,'W','l/s/p',IER)
            CALL EGETWR(hold,K,VAL3,0.1,999.,'W','number people',IER)
            IF(IER.NE.0) GOTO 43

            VAL= (val2/1000.0) * val3
            write(prompt,'(3a,F8.5,a)') ' The flow via ',
     &        CMNAM(NCMP),' (m^3/s) is ',VAL,'.'
            call easkok(prompt,'Is this OK?',OK,nbhelp)
            if (.NOT.OK) goto 42
            SUPCMP(NCMP,2)=VAL
          elseif(iw.eq.3)then

C Convert any existing value to m3/hr
            write(prompt,'(3a,F8.5,a)') ' The flow via ',
     &        CMNAM(NCMP),' (m^3/s) is ',SUPCMP(NCMP,2),'.'
            call edisp(iuout,prompt)
            VAL2=SUPCMP(NCMP,2)*3600.0
            write(prompt,'(3a)') 'Flow through ',CMNAM(NCMP),
     &        ' (m^3/hr)?'
   21       CALL EASKR(VAL2,' ',prompt,0.,'W',0.,'-',100.,
     &        'flow m3/hr in comp',IER,nbhelp)
            IF(IER.NE.0) GOTO 21
            if(val2.gt.0.0) VAL = (val2/3600.0)
            write(prompt,'(3a,F8.5,a)') ' The flow through ',
     &        CMNAM(NCMP),' (m^3/s) is ',VAL,'.'
            call easkok(prompt,'Is this OK?',OK,nbhelp)
            if (.NOT.OK) goto 42
            SUPCMP(NCMP,2)=VAL
          elseif(iw.eq.4)then

C Based on ach against the volume of a particular zone.
            write(prompt,'(3a,F8.5,a)') ' The flow via ',
     &        CMNAM(NCMP),' (m^3/s) is ',SUPCMP(NCMP,2),'.'
            call edisp(iuout,prompt)
            CALL EPICKS(IZN,IVALS2,' ',OUTS,12,NCOMP,zname,
     &        'zone list',IER,nbhelp)
            write(prompt,*)'Associated volume is ',VOL(IVALS2(1))
            call edisp(iuout,prompt)
            if(SUPCMP(NCMP,2).lt.0.001) val2=1.0
            write(prompt,'(3a)') 'Air changes for ',CMNAM(NCMP),
     &        ' (ach)?'
   22       CALL EASKR(VAL2,' ',prompt,
     &        0.,'W',0.,'-',1.,'flow ach in comp',IER,nbhelp)
            IF(IER.NE.0) GOTO 22
            VAL = (SUPNOD(NCMP,2) * VAL2)/3600.0
            write(prompt,'(3a,F8.5,a)') ' The flow through ',
     &        CMNAM(NCMP),' (m^3/s) is ',VAL,'.'
            call easkok(prompt,'Is this OK?',OK,nbhelp)
            if (.NOT.OK) goto 42
            SUPCMP(NCMP,2)=VAL
          endif
          call ECLOSE(SUPCMP(NCMP,1),2.0,0.0001,iswater)
          if(.NOT.iswater)then
            call symbol_flow_cmp(name,'bx','r',xc,yc,zc)
          else
            call symbol_flow_cmp(name,'bx','b',xc,yc,zc)
          endif
          call pausems(30)

C Report on Delta Height with associated surf COG.
          write(outs,'(a,f7.3)') 'Delta Zone with surface COG is',
     &      SURCOG(iz,is,3)-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call pausems(30)
          if(ICT(loop).eq.0)then
            write(name,'(a,i2.2,a,i3.3)') 'BW-EX',iz,':',is
            xw=x3+SURCOG(iz,is,1);yw=y3+SURCOG(iz,is,2)
            zw=z3+SURCOG(iz,is,3)
            write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),'  ',name,
     &        ' Wind boundary extract azim',spazi(iz,is),
     &        ' @XYZ',xw,yw,zw
            call edisp(iuout,outs)
            NNOD=NNOD+1; NDTYP(NNOD)=3; NDFLD(NNOD)=1
            HNOD(NNOD,1)=xw; HNOD(NNOD,2)=yw; HNOD(NNOD,3)=zw
            SUPNOD(NNOD,2)=spazi(iz,is)
            SUPNOD(NNOD,1)=1.0
            TNOD(NNOD)=0.0
            if(znotair(iz))then
              SUPNOD(NNOD,1)=2  ! mark flow node as water
              NDFLD(NNOD)=2
            else
              SUPNOD(NNOD,1)=1  ! mark flow node as air
              NDFLD(NNOD)=1
            endif
            write(NDNAM(NNOD),'(a)') name
            NODASSOC(NNOD,1)=ZNAME(IZ); NODASSOC(NNOD,2)=SN
            call symbol_flow_nd(name,'bw','r',xw,yw,zw)
            call pausems(30)

C Establish connection from zone to this boundary node via extract
C and then to the adjacent zone.
            POS=ICAAS(IZ); NEG=NNOD; CMPIS=NCMP; NCNN=NCNN+1
            if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
              ccol='r'
            else
              ccol='b'
            endif

C Reality checks here...
            NODPS(NCNN)=POS; NODNE(NCNN)=NEG
            HGTPS(NCNN)=0.0   ! an extract has zero height differences
            HGTNE(NCNN)=0.0   ! an extract has zero height differences
            NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
            ITPCON(NCNN)=CMPIS ! use the just-defined crack
            write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &        NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
            call edisp(iuout,outs)
            call draw_flow_link(xw,yw,zw,xc,yc,zc,ccol)
            call draw_flow_link(xc,yc,zc,
     &        HNOD(NEG,1),HNOD(NEG,2),HNOD(NEG,3),ccol)
          elseif(ICT(loop).eq.3)then

C Define connection with adjacent zone. As this is an extract
C the current zone is positive. Check if linkatpartition is still
C zero. If so mark this connection and update the other surface as well.
            if(linkatpartition(loop).eq.0)then
              POS=ICAAS(IC1(loop))  ! the current zone node is positive
              NEG=ICAAS(IC2(loop)) ! the other zone
              if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
                ccol='r'
              else
                ccol='b'
              endif
              CMPIS=NCMP; NCNN=NCNN+1
              linkatpartition(loop)=1
              icoth=IZSTOCN(IC2(loop),IE2(loop))
              linkatpartition(icoth)=1

C Reality checks here...
              NODPS(NCNN)=POS; NODNE(NCNN)=NEG
              HGTPS(NCNN)= HCMP(NCMP,1,3)-ZCOG(IZ,3)  ! diff comp & zone
              HGTNE(NCNN)= HCMP(NCMP,1,3)-ZCOG(IC2(loop),3)  ! diff comp & other zone
              NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
              ITPCON(NCNN)=CMPIS ! use the just-defined crack
              write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &          NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
              call edisp(iuout,outs)
              call draw_flow_link(ZCOG(IZ,1),ZCOG(IZ,2),ZCOG(IZ,3),
     &          HCMP(NCMP,1,1),HCMP(NCMP,1,2),HCMP(NCMP,1,3),ccol)
              call draw_flow_link(ZCOG(IC2(loop),1),ZCOG(IC2(loop),2),
     &          ZCOG(IC2(loop),3),HCMP(NCMP,1,1),HCMP(NCMP,1,2),
     &          HCMP(NCMP,1,3),ccol)
            endif
          else
            continue  ! << for similar creat another node ....
          endif

        elseif(SUSE(iz,is,2)(1:4).eq.'OPEN')then

C Grill opening - common orifice component, (type 40) IVALCM(10) MF040I
          call zsurfprm(iz,is,dupedges,perim)
          write(name,'(a,i2.2,a,i3.3)') 'GrOPz',iz,':',is
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=SURCOG(iz,is,3)
          write(outs,'(4a,f6.1,a,3f6.2)') zsn(1:lnzsn),'  ',name,
     &      ' grill orifice area',SNA(iz,is),'m^2 COG @XYZ',
     &      xc,yc,zc
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          NCMP=NCMP+1 ! use logic in MF040I
          CMNAM(NCMP)=name
          CMPASSOC(NCMP,1)=ZNAME(IZ); CMPASSOC(NCMP,2)=SN
          ITPCMP(NCMP)=IVALCM(10); LTPCMP(NCMP)=LVALCM(10) ! common orifice
          ISDCMP(NCMP)=3; ISDCNN(NCMP)=0
          SUPCMP(NCMP,1)=1.
          SUPCMP(NCMP,2)=SNA(iz,is)
          if(SUPCMP(NCMP,3).lt.0.001) SUPCMP(NCMP,3)=0.65
          write(prompt,'(3a)')
     &      'The discharge coef of ',CMNAM(NCMP),':'
          VAL2=SUPCMP(NCMP,3)
          CALL EASKR(VAL2,' ',prompt,
     &      0.,'W',0.,'-',10.,'discharge coef',IER,nbhelp)
          SUPCMP(NCMP,3)=VAL2
          NWPCMP(NCMP)=0  ! start and end the same
          HCMP(NCMP,1,1)=xc; HCMP(NCMP,1,2)=yc; HCMP(NCMP,1,3)=zc
          call symbol_flow_cmp(name,'or','r',xc,yc,zc)
          call pausems(30)

C Report on Delta Height with associated surf COG.
          write(outs,'(a,f7.3)') 'Delta Zone with surface COG is',
     &      SURCOG(iz,is,3)-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call pausems(30)
          if(ICT(loop).eq.0)then
            write(name,'(a,i2.2,a,i3.3)') 'BW-Op',iz,':',is
            xw=x3+SURCOG(iz,is,1);yw=y3+SURCOG(iz,is,2)
            zw=z3+SURCOG(iz,is,3)
            write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),'  ',name,
     &      ' Wind boundary opening azim',spazi(iz,is),
     &      ' @XYZ',xw,yw,zw
            call edisp(iuout,outs)
            NNOD=NNOD+1; NDTYP(NNOD)=3; NDFLD(NNOD)=1
            HNOD(NNOD,1)=xw; HNOD(NNOD,2)=yw; HNOD(NNOD,3)=zw
            SUPNOD(NNOD,2)=spazi(iz,is)
            SUPNOD(NNOD,1)=1.0
            TNOD(NNOD)=0.0
            if(znotair(iz))then
              SUPNOD(NNOD,1)=2  ! mark flow node as water
              NDFLD(NNOD)=2
            else
              SUPNOD(NNOD,1)=1  ! mark flow node as air
              NDFLD(NNOD)=1
            endif
            write(NDNAM(NNOD),'(a)') name
            NODASSOC(NNOD,1)=ZNAME(IZ); NODASSOC(NNOD,2)=SN
            call symbol_flow_nd(name,'bw','r',xw,yw,zw)
            call pausems(30)

C Establish connection from this boundary node to the orifice
C and then to the adjacent zone.
            POS=NNOD; NEG=ICAAS(IZ); CMPIS=NCMP; NCNN=NCNN+1
            if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
              ccol='r'
            else
              ccol='b'
            endif

C Reality checks here...
            NODPS(NCNN)=POS; NODNE(NCNN)=NEG
            HGTPS(NCNN)=0.0   ! boundary always in line with crack
            HGTNE(NCNN)=SURCOG(iz,is,3)-ZCOG(IZ,3)
            NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
            ITPCON(NCNN)=CMPIS ! use the just-defined crack
            write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &        NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
            call edisp(iuout,outs)
            call draw_flow_link(xw,yw,zw,xc,yc,zc,ccol)
            call draw_flow_link(xc,yc,zc,
     &        HNOD(NEG,1),HNOD(NEG,2),HNOD(NEG,3),ccol)
          elseif(ICT(loop).eq.3)then

C Define connection with adjacent zone. Check if linkatpartition is still
C zero. If so mark this connection and update the other surface as well.
            if(linkatpartition(loop).eq.0)then
              POS=ICAAS(IC1(loop))  ! the current zone node is positive
              NEG=ICAAS(IC2(loop)) ! the other zone
              if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
                ccol='r'
              else
                ccol='b'
              endif
              CMPIS=NCMP; NCNN=NCNN+1
              linkatpartition(loop)=1
              icoth=IZSTOCN(IC2(loop),IE2(loop))
              linkatpartition(icoth)=1

C Reality checks here...
              NODPS(NCNN)=POS; NODNE(NCNN)=NEG
              HGTPS(NCNN)= HCMP(NCMP,1,3)-ZCOG(IZ,3)  ! diff comp & zone
              HGTNE(NCNN)= HCMP(NCMP,1,3)-ZCOG(IC2(loop),3)  ! diff comp & other zone
              NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
              ITPCON(NCNN)=CMPIS ! use the just-defined crack
              write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &          NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
              call edisp(iuout,outs)
              call draw_flow_link(ZCOG(IZ,1),ZCOG(IZ,2),ZCOG(IZ,3),
     &          HCMP(NCMP,1,1),HCMP(NCMP,1,2),HCMP(NCMP,1,3),ccol)
              call draw_flow_link(ZCOG(IC2(loop),1),ZCOG(IC2(loop),2),
     &          ZCOG(IC2(loop),3),HCMP(NCMP,1,1),HCMP(NCMP,1,2),
     &          HCMP(NCMP,1,3),ccol)
            endif
          else
            continue  ! << for similar creat another node ....
          endif

        endif
      elseif(SUSE(iz,is,1)(1:6).eq.'WINDOW'.or.
     &    SUSE(iz,is,1)(1:8).eq.'D-WINDOW'.or.
     &    SUSE(iz,is,1)(1:8).eq.'S-WINDOW'.or.
     &    SUSE(iz,is,1)(1:8).eq.'C-WINDOW')then
        if(SUSE(iz,is,2)(1:6).eq.'CLOSED'.or.
     &     SUSE(iz,is,2)(1:5).eq.'CRACK')then

C Window crack - crack component, assume width = 1mm,
C length = window perimeter, at surface COG.
          call zsurfprm(iz,is,dupedges,perim)
          write(name,'(a,i2.2,a,i3.3)') 'WiCrz',iz,':',is
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=SURCOG(iz,is,3)
          write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),'  ',name,
     &      ' window crack width 2mm length',perim,' @XYZ',xc,yc,zc
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          NCMP=NCMP+1 ! use logic in MF120I
          CMNAM(NCMP)=name
          CMPASSOC(NCMP,1)=ZNAME(IZ); CMPASSOC(NCMP,2)=SN
          ITPCMP(NCMP)=IVALCM(13); LTPCMP(NCMP)=LVALCM(13) !crack
          ISDCMP(NCMP)=3; ISDCNN(NCMP)=0
          SUPCMP(NCMP,1)=1.
          SUPCMP(NCMP,2)=2.0/1000.0
          SUPCMP(NCMP,3)=perim
          NWPCMP(NCMP)=0  ! start and end the same
          HCMP(NCMP,1,1)=xc; HCMP(NCMP,1,2)=yc; HCMP(NCMP,1,3)=zc
          call symbol_flow_cmp(name,'cr','r',xc,yc,zc)
          call pausems(30)

C Report on Delta Height with associated surf lower.
          write(outs,'(a,f7.3)') 
     &      'Delta Zone with surface lower edge is',
     &      SURCOG(iz,is,3)-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call pausems(30)
          if(ICT(loop).eq.0)then
            write(name,'(a,i2.2,a,i3.3)') 'BW-Cr',iz,':',is
            xw=x3+SURCOG(iz,is,1);yw=y3+SURCOG(iz,is,2)
            zw=z3+SURCOG(iz,is,3)
            write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),'  ',name,
     &        ' Wind boundary crack azim',spazi(iz,is),
     &        ' @XYZ',xw,yw,zw
            call edisp(iuout,outs)
            NNOD=NNOD+1; NDTYP(NNOD)=3; NDFLD(NNOD)=1
            HNOD(NNOD,1)=xw; HNOD(NNOD,2)=yw; HNOD(NNOD,3)=zw
            SUPNOD(NNOD,2)=spazi(iz,is)
            SUPNOD(NNOD,1)=1.0
            TNOD(NNOD)=0.0
            write(NDNAM(NNOD),'(a)') name
            NODASSOC(NNOD,1)=ZNAME(IZ); NODASSOC(NNOD,2)=SN
            call symbol_flow_nd(name,'bw','r',xw,yw,zw)
            call pausems(30)

C Establish connection from this boundary node to the crack
C and then to the adjacent zone.
            POS=NNOD; NEG=ICAAS(IZ); CMPIS=NCMP; NCNN=NCNN+1
            if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
              ccol='r'
            else
              ccol='b'
            endif

C Reality checks here...
            NODPS(NCNN)=POS; NODNE(NCNN)=NEG
            HGTPS(NCNN)=0.0   ! boundary always in line with crack
            HGTNE(NCNN)=SURCOG(iz,is,3)-ZCOG(IZ,3)
            NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
            ITPCON(NCNN)=CMPIS ! use the just-defined crack
            write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &        NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
            call edisp(iuout,outs)
            call draw_flow_link(xw,yw,zw,xc,yc,zc,ccol)
            call draw_flow_link(xc,yc,zc,
     &        HNOD(NEG,1),HNOD(NEG,2),HNOD(NEG,3),ccol)
          elseif(ICT(loop).eq.3)then

C Define connection with adjacent zone. Check if linkatpartition is still
C zero. If so mark this connection and update the other surface as well.
            if(linkatpartition(loop).eq.0)then
              POS=ICAAS(IC1(loop))  ! the current zone node is positive
              NEG=ICAAS(IC2(loop)) ! the other zone
              if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
                ccol='r'
              else
                ccol='b'
              endif
              CMPIS=NCMP; NCNN=NCNN+1
              linkatpartition(loop)=1
              icoth=IZSTOCN(IC2(loop),IE2(loop))
              linkatpartition(icoth)=1

C Reality checks here...
              NODPS(NCNN)=POS; NODNE(NCNN)=NEG
              HGTPS(NCNN)= HCMP(NCMP,1,3)-ZCOG(IZ,3)  ! diff comp & zone
              HGTNE(NCNN)= HCMP(NCMP,1,3)-ZCOG(IC2(loop),3)  ! diff comp & other zone
              NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
              ITPCON(NCNN)=CMPIS ! use the just-defined crack
              write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &          NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
              call edisp(iuout,outs)
              call draw_flow_link(ZCOG(IZ,1),ZCOG(IZ,2),ZCOG(IZ,3),
     &          HCMP(NCMP,1,1),HCMP(NCMP,1,2),HCMP(NCMP,1,3),ccol)
              call draw_flow_link(ZCOG(IC2(loop),1),ZCOG(IC2(loop),2),
     &          ZCOG(IC2(loop),3),HCMP(NCMP,1,1),HCMP(NCMP,1,2),
     &          HCMP(NCMP,1,3),ccol)
            endif
          else
            continue  ! << for similar creat another node ....
          endif

        elseif(SUSE(iz,is,2)(1:4).eq.'OPEN')then

C Window open - common orifice component, (type 40) IVALCM(10) MF040I
C Most windows do not fully open so ask user the % of surfaces area.
          write(name,'(a,i2.2,a,i3.3)') 'WiOpz',iz,':',is
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=SURCOG(iz,is,3)
          write(prompt,'(4a,f6.1,a,f6.1,a,f6.2,a,3f7.3)') 
     &      zsn(1:lnzsn),'  ',name,
     &      ' width ',XYMAX,' height',ZMAX,
     &      ' area',SNA(iz,is),'m^2 COG @Z',xc,yc,zc
          VAL2=25.0              ! 25% initial guess of area
          CALL EASKR(VAL2,prompt,'Area % for flow opening?',
     &      0.,'W',1.,'-',100.,'% area',IER,nbhelp)
          VAL=SNA(iz,is)*(VAL2*0.001)
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          NCMP=NCMP+1 ! use logic in MF040I
          CMNAM(NCMP)=name
          CMPASSOC(NCMP,1)=ZNAME(IZ); CMPASSOC(NCMP,2)=SN
          ITPCMP(NCMP)=IVALCM(10); LTPCMP(NCMP)=LVALCM(10) ! common orifice
          ISDCMP(NCMP)=3; ISDCNN(NCMP)=0
          SUPCMP(NCMP,1)=1.
          SUPCMP(NCMP,2)=VAL
          if(SUPCMP(NCMP,3).lt.0.001) SUPCMP(NCMP,3)=0.5
          write(prompt,'(3a)')
     &      'The discharge coef. for ',CMNAM(NCMP),':'
          VAL2=SUPCMP(NCMP,3)
          CALL EASKR(VAL2,' ',prompt,
     &      0.,'W',0.,'-',10.,'discharge coef',IER,nbhelp)
          SUPCMP(NCMP,3)=VAL2
          call edisp(iuout,outs)
          NWPCMP(NCMP)=0  ! start and end the same
          HCMP(NCMP,1,1)=xc; HCMP(NCMP,1,2)=yc; HCMP(NCMP,1,3)=zc
          call symbol_flow_cmp(name,'or','r',xc,yc,zc)
          call pausems(30)

C Report on Delta Height with associated surf COG.
          write(outs,'(a,f7.3)') 'Delta Zone with surface COG is',
     &      SURCOG(iz,is,3)-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call pausems(30)
          if(ICT(loop).eq.0)then
            write(name,'(a,i2.2,a,i3.3)') 'BW-Op',iz,':',is
            xw=x3+SURCOG(iz,is,1);yw=y3+SURCOG(iz,is,2)
            zw=z3+SURCOG(iz,is,3)
            write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),'  ',name,
     &        ' Wind boundary opening azim',spazi(iz,is),
     &        ' @XYZ',xw,yw,zw
            call edisp(iuout,outs)
            NNOD=NNOD+1; NDTYP(NNOD)=3; NDFLD(NNOD)=1
            HNOD(NNOD,1)=xw; HNOD(NNOD,2)=yw; HNOD(NNOD,3)=zw
            SUPNOD(NNOD,2)=spazi(iz,is)
            SUPNOD(NNOD,1)=1.0
            TNOD(NNOD)=0.0
            write(NDNAM(NNOD),'(a)') name
            NODASSOC(NNOD,1)=ZNAME(IZ); NODASSOC(NNOD,2)=SN
            call symbol_flow_nd(name,'bw','r',xw,yw,zw)
            call pausems(30)

C Establish connection from this boundary node to the orifice
C and then to the adjacent zone.
            POS=NNOD; NEG=ICAAS(IZ); CMPIS=NCMP; NCNN=NCNN+1
            if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
              ccol='r'
            else
              ccol='b'
            endif

C Reality checks here...
            NODPS(NCNN)=POS; NODNE(NCNN)=NEG
            HGTPS(NCNN)=0.0   ! boundary always in line with crack
            HGTNE(NCNN)=SURCOG(iz,is,3)-ZCOG(IZ,3)
            NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
            ITPCON(NCNN)=CMPIS ! use the just-defined crack
            write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &        NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
            call edisp(iuout,outs)
            call draw_flow_link(xw,yw,zw,xc,yc,zc,ccol)
            call draw_flow_link(xc,yc,zc,
     &        HNOD(NEG,1),HNOD(NEG,2),HNOD(NEG,3),ccol)
          elseif(ICT(loop).eq.3)then

C Define connection with adjacent zone. Check if linkatpartition is still
C zero. If so mark this connection and update the other surface as well.
            if(linkatpartition(loop).eq.0)then
              POS=ICAAS(IC1(loop))  ! the current zone node is positive
              NEG=ICAAS(IC2(loop)) ! the other zone
              if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
                ccol='r'
              else
                ccol='b'
              endif
              CMPIS=NCMP; NCNN=NCNN+1
              linkatpartition(loop)=1
              icoth=IZSTOCN(IC2(loop),IE2(loop))
              linkatpartition(icoth)=1

C Reality checks here...
              NODPS(NCNN)=POS; NODNE(NCNN)=NEG
              HGTPS(NCNN)= HCMP(NCMP,1,3)-ZCOG(IZ,3)  ! diff comp & zone
              HGTNE(NCNN)= HCMP(NCMP,1,3)-ZCOG(IC2(loop),3)  ! diff comp & other zone
              NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
              ITPCON(NCNN)=CMPIS ! use the just-defined crack
              write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &          NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
              call edisp(iuout,outs)
              call draw_flow_link(ZCOG(IZ,1),ZCOG(IZ,2),ZCOG(IZ,3),
     &          HCMP(NCMP,1,1),HCMP(NCMP,1,2),HCMP(NCMP,1,3),ccol)
              call draw_flow_link(ZCOG(IC2(loop),1),ZCOG(IC2(loop),2),
     &          ZCOG(IC2(loop),3),HCMP(NCMP,1,1),HCMP(NCMP,1,2),
     &          HCMP(NCMP,1,3),ccol)
            endif
          else
            continue  ! << for similar creat another node ....
          endif

        elseif(SUSE(iz,is,2)(1:4).eq.'SASH')then

C Window sash - 2 common orifice component, (type 40) IVALCM(10) MF040I
C For lower opening, use surface width and ask user for height.
C For upper opening, use surface width and ask user for height.
          write(name,'(a,i2.2,a,i3.3)') 'WiSLz',iz,':',is
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=edgelower
          write(prompt,'(2a,f6.1,a,f6.1,a,f7.3,a,3f6.2)') 
     &      zsn(1:lnzsn),' width ',XYMAX,' height',ZMAX,
     &      ' area',SNA(iz,is),'m2 sill @XYZ',xc,yc,zc
          VAL2=ZMAX*0.125
          CALL EASKR(VAL2,prompt,'Height (m) for lower opening?',
     &      0.,'W',1.,'-',100.,'height for lower',IER,nbhelp)
          VAL=XYMAX*VAL2
          write(outs,'(4a,f6.1,a,f6.1,a,f6.2,a,3f6.2)') 
     &      zsn(1:lnzsn),'  ',name,
     &      ' width ',XYMAX,' height',VAL2,
     &      ' area',VAL,'m2 lower @XYZ',xc,yc,zc
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          NCMP=NCMP+1 ! use logic in MF040I
          CMNAM(NCMP)=name
          CMPASSOC(NCMP,1)=ZNAME(IZ); CMPASSOC(NCMP,2)=SN
          ITPCMP(NCMP)=IVALCM(10); LTPCMP(NCMP)=LVALCM(10)  ! common orifice opening
          ISDCMP(NCMP)=3; ISDCNN(NCMP)=0
          SUPCMP(NCMP,1)=1.
          SUPCMP(NCMP,2)=VAL
          if(SUPCMP(NCMP,3).lt.0.001) SUPCMP(NCMP,3)=0.5
          write(prompt,'(3a)')
     &      'The discharge coef of ',CMNAM(NCMP),':'
          VAL2=SUPCMP(NCMP,3)
          CALL EASKR(VAL2,' ',prompt,0.,'W',0.,'-',10.,
     &      'discharge coef',IER,nbhelp)
          SUPCMP(NCMP,3)=VAL2
          NWPCMP(NCMP)=0  ! start and end the same
          HCMP(NCMP,1,1)=xc; HCMP(NCMP,1,2)=yc; HCMP(NCMP,1,3)=zc
          call symbol_flow_cmp(name,'or','r',xc,yc,zc)
          call pausems(30)

C Report on Delta Height with associated surf lower.
          write(outs,'(a,f7.3)') 
     &      'Delta Zone with surface lower edge is',
     &      edgelower-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call pausems(30)

C If at boundary setup the lower boundary node.
          if(ICT(loop).eq.0)then
            write(name,'(a,i2.2,a,i3.3)') 'BW-SL',iz,':',is
            xw=x3+SURCOG(iz,is,1);yw=y3+SURCOG(iz,is,2)
            zw=z3+edgelower
            write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),'  ',name,
     &      ' Wind boundary sash azim',spazi(iz,is),
     &      ' @XYZ',xw,yw,zw
            call edisp(iuout,outs)
            NNOD=NNOD+1; NDTYP(NNOD)=3; NDFLD(NNOD)=1
            HNOD(NNOD,1)=xw; HNOD(NNOD,2)=yw; HNOD(NNOD,3)=zw
            SUPNOD(NNOD,2)=spazi(iz,is)
            SUPNOD(NNOD,1)=1.0
            TNOD(NNOD)=0.0
            write(NDNAM(NNOD),'(a)') name
            NODASSOC(NNOD,1)=ZNAME(IZ); NODASSOC(NNOD,2)=SN
            call symbol_flow_nd(name,'bw','r',xw,yw,zw)
            call pausems(30)

C Establish connection from this boundary node to the lower orifice
C and then to the adjacent zone.
            POS=NNOD; NEG=ICAAS(IZ); CMPIS=NCMP; NCNN=NCNN+1
            if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
              ccol='r'
            else
              ccol='b'
            endif

C Reality checks here...
            NODPS(NCNN)=POS; NODNE(NCNN)=NEG
            HGTPS(NCNN)=0.0   ! boundary always in line with orifice
            HGTNE(NCNN)=edgelower-ZCOG(IZ,3)
            NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
            ITPCON(NCNN)=CMPIS ! use the just-defined lower orifice
            write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &        NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
            call edisp(iuout,outs)
            call draw_flow_link(xw,yw,zw,xc,yc,zc,ccol)
            call draw_flow_link(xc,yc,zc,
     &        HNOD(NEG,1),HNOD(NEG,2),HNOD(NEG,3),ccol)
          elseif(ICT(loop).eq.3)then

C Define connection with adjacent zone. Check if linkatpartition is still
C zero. If so mark this connection and update the other surface as well.
            if(linkatpartition(loop).eq.0)then
              POS=ICAAS(IC1(loop))  ! the current zone node is positive
              NEG=ICAAS(IC2(loop)) ! the other zone
              if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
                ccol='r'
              else
                ccol='b'
              endif
              CMPIS=NCMP; NCNN=NCNN+1
              linkatpartition(loop)=1
              icoth=IZSTOCN(IC2(loop),IE2(loop))
              linkatpartition(icoth)=1

C Reality checks here...
              NODPS(NCNN)=POS; NODNE(NCNN)=NEG
              HGTPS(NCNN)= HCMP(NCMP,1,3)-ZCOG(IZ,3)  ! diff comp & zone
              HGTNE(NCNN)= HCMP(NCMP,1,3)-ZCOG(IC2(loop),3)  ! diff comp & other zone
              NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
              ITPCON(NCNN)=CMPIS ! use the just-defined crack
              write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &          NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
              call edisp(iuout,outs)
              call draw_flow_link(ZCOG(IZ,1),ZCOG(IZ,2),ZCOG(IZ,3),
     &          HCMP(NCMP,1,1),HCMP(NCMP,1,2),HCMP(NCMP,1,3),ccol)
              call draw_flow_link(ZCOG(IC2(loop),1),ZCOG(IC2(loop),2),
     &          ZCOG(IC2(loop),3),HCMP(NCMP,1,1),HCMP(NCMP,1,2),
     &          HCMP(NCMP,1,3),ccol)
            endif
          else
            continue  ! << for similar creat another node ....
          endif

C Now create the upper orifice.
          write(name,'(a,i2.2,a,i3.3)') 'WiSUz',iz,':',is
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=edgeupper
          write(prompt,'(2a,f6.1,a,f6.1,a,f7.3,a,3f6.2)') 
     &      zsn(1:lnzsn),' width ',XYMAX,' height',ZMAX,
     &      ' area',SNA(iz,is),'m2 head @XYZ',xc,yc,zc
          VAL2=ZMAX*0.125
          CALL EASKR(VAL2,prompt,'Height (m) for upper opening?',
     &      0.,'W',1.,'-',100.,'height for upper',IER,nbhelp)
          VAL=XYMAX*VAL2
          write(outs,'(4a,f6.1,a,f6.1,a,f6.2,a,3f6.2)') 
     &      zsn(1:lnzsn),'  ',name,
     &      ' width ',XYMAX,' height',VAL2,
     &      ' area',VAL,'m2 upper @XYZ',xc,yc,zc
          NCMP=NCMP+1 ! use logic in MF040I
          CMNAM(NCMP)=name
          CMPASSOC(NCMP,1)=ZNAME(IZ); CMPASSOC(NCMP,2)=SN
          ITPCMP(NCMP)=IVALCM(10); LTPCMP(NCMP)=LVALCM(10)  ! common orifice opening
          ISDCMP(NCMP)=3; ISDCNN(NCMP)=0
          SUPCMP(NCMP,1)=1.
          SUPCMP(NCMP,2)=VAL
          if(SUPCMP(NCMP,3).lt.0.001) SUPCMP(NCMP,3)=0.5
          write(prompt,'(3a)')
     &      'The discharge coef of ',CMNAM(NCMP),':'
          VAL2=SUPCMP(NCMP,3)
          CALL EASKR(VAL2,' ',prompt,
     &      0.,'W',0.,'-',10.,'discharge coef',IER,nbhelp)
          SUPCMP(NCMP,3)=VAL2
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          NWPCMP(NCMP)=0  ! start and end the same
          HCMP(NCMP,1,1)=xc; HCMP(NCMP,1,2)=yc; HCMP(NCMP,1,3)=zc
          call symbol_flow_cmp(name,'or','r',xc,yc,zc)
          call pausems(30)

C Report on Delta Height with associated surf upper.
          write(outs,'(a,f7.3)') 'Delta Zone with surf upper is',
     &      edgeupper-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call pausems(30)

          if(ICT(loop).eq.0)then

C Now the boundary at the upper orifice.
            write(name,'(a,i2.2,a,i3.3)') 'BW-SU',iz,':',is
            xw=x3+SURCOG(iz,is,1);yw=y3+SURCOG(iz,is,2)
            zw=z3+edgeupper
            write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),'  ',name,
     &      ' Wind boundary sash azim',spazi(iz,is),
     &      ' @XYZ',xw,yw,zw
            call edisp(iuout,outs)
            NNOD=NNOD+1; NDTYP(NNOD)=3; NDFLD(NNOD)=1
            HNOD(NNOD,1)=xw; HNOD(NNOD,2)=yw; HNOD(NNOD,3)=zw
            SUPNOD(NNOD,2)=spazi(iz,is)
            SUPNOD(NNOD,1)=1.0
            TNOD(NNOD)=0.0
            write(NDNAM(NNOD),'(a)') name
            NODASSOC(NNOD,1)=ZNAME(IZ); NODASSOC(NNOD,2)=SN
            call symbol_flow_nd(name,'bw','r',xw,yw,zw)
            call pausems(30)

C Establish connection from this boundary node to the upper orifice
C and then to the adjacent zone.
            POS=NNOD; NEG=ICAAS(IZ); CMPIS=NCMP; NCNN=NCNN+1
            if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
              ccol='r'
            else
              ccol='b'
            endif

C Reality checks here...
            NODPS(NCNN)=POS; NODNE(NCNN)=NEG
            HGTPS(NCNN)=0.0   ! boundary always in line with orifice
            HGTNE(NCNN)=edgeupper-ZCOG(IZ,3)
            NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
            ITPCON(NCNN)=CMPIS ! use the just-defined lower orifice
            write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &        NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
            call edisp(iuout,outs)
            call draw_flow_link(xw,yw,zw,xc,yc,zc,ccol)
            call draw_flow_link(xc,yc,zc,
     &        HNOD(NEG,1),HNOD(NEG,2),HNOD(NEG,3),ccol)
          elseif(ICT(loop).eq.3)then

C Define connection with adjacent zone. Check if linkatpartition is still
C zero. If so mark this connection and update the other surface as well.
            if(linkatpartition(loop).eq.0)then
              POS=ICAAS(IC1(loop))  ! the current zone node is positive
              NEG=ICAAS(IC2(loop)) ! the other zone
              if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
                ccol='r'
              else
                ccol='b'
              endif
              CMPIS=NCMP; NCNN=NCNN+1
              linkatpartition(loop)=1
              icoth=IZSTOCN(IC2(loop),IE2(loop))
              linkatpartition(icoth)=1

C Reality checks here...
              NODPS(NCNN)=POS; NODNE(NCNN)=NEG
              HGTPS(NCNN)= HCMP(NCMP,1,3)-ZCOG(IZ,3)  ! diff comp & zone
              HGTNE(NCNN)= HCMP(NCMP,1,3)-ZCOG(IC2(loop),3)  ! diff comp & other zone
              NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
              ITPCON(NCNN)=CMPIS ! use the just-defined crack
              write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &          NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
              call edisp(iuout,outs)
              call draw_flow_link(ZCOG(IZ,1),ZCOG(IZ,2),ZCOG(IZ,3),
     &          HCMP(NCMP,1,1),HCMP(NCMP,1,2),HCMP(NCMP,1,3),ccol)
              call draw_flow_link(ZCOG(IC2(loop),1),ZCOG(IC2(loop),2),
     &          ZCOG(IC2(loop),3),HCMP(NCMP,1,1),HCMP(NCMP,1,2),
     &          HCMP(NCMP,1,3),ccol)
            endif
          else
            continue  ! << for similar creat another node ....
          endif

        elseif(SUSE(iz,is,2)(1:5).eq.'BIDIR'.or.
     &         SUSE(iz,is,2)(1:10).eq.'AJAR-BIDIR')then

C Window bidirectional - specific door component width & height
C based on the surface width & height. discharge factor = 0.6,
C Draw on the centre line. The VAL4 is Surf COG to Zone COG.
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=SURCOG(iz,is,3)
          if(SUSE(iz,is,2)(1:5).eq.'BIDIR')then
            write(name,'(a,i2.2,a,i3.3)') 'WiBiz',iz,':',is
            VAL1=XYMAX; VAL2=ZMAX
            VAL3=ABS(ZCOG(iz,3)-edgelower)  ! base to zone COG
            VAL4=SURCOG(iz,is,3)-ZCOG(IZ,3)        ! surf COG to zone COG
            write(outs,'(4a,f6.1,a,f6.1,a,3f7.3,a,f6.3)')
     &        zsn(1:lnzsn),'  ',name,
     &        ' bidir width ',XYMAX,' height',ZMAX,
     &        ' coef 0.6 @XYZ',xc,yc,zc,
     &        ' delta base to zone COG',VAL3
          elseif(SUSE(iz,is,2)(1:10).eq.'AJAR-BIDIR')then

C Window bidirectional ajar - specific door component width 100mm &
C height based on the surface height. discharge factor = 0.6,
C Draw on the centre line. The VAL4 is Surf COG to Zone COG.
            write(name,'(a,i2.2,a,i3.3)') 'WABiz',iz,':',is
            VAL1=0.1; VAL2=ZMAX
            VAL3=ABS(ZCOG(iz,3)-edgelower)  ! base to zone COG
            VAL4=SURCOG(iz,is,3)-ZCOG(IZ,3)        ! surf COG to zone COG
            write(outs,'(4a,f6.1,a,3f7.3,a,f6.3)')
     &        zsn(1:lnzsn),'  ',name,' bidir ajar w=100mm h=',ZMAX,
     &        ' coef 0.6 @XYZ',xc,yc,zc,
     &        ' delta base to zone COG',VAL3
          endif
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          NCMP=NCMP+1 ! use logic in MF040I
          CMNAM(NCMP)=name
          CMPASSOC(NCMP,1)=ZNAME(IZ); CMPASSOC(NCMP,2)=SN
          ITPCMP(NCMP)=IVALCM(14); LTPCMP(NCMP)=LVALCM(14) ! bi-directional
          ISDCMP(NCMP)=5; ISDCNN(NCMP)=0
          SUPCMP(NCMP,1)=1.
          SUPCMP(NCMP,2)=VAL1
          SUPCMP(NCMP,3)=VAL2
          SUPCMP(NCMP,4)=VAL3
          if(SUPCMP(NCMP,5).lt.0.001) SUPCMP(NCMP,5)=0.6
          write(prompt,'(3a)') 'The discharge coefficient of ',
     &      CMNAM(NCMP),':'
          VAL2=SUPCMP(NCMP,5)
          CALL EASKR(VAL2,' ',prompt,0.,'W',0.,'-',10.,
     &      'discharge coef in comp',IER,nbhelp)
          SUPCMP(NCMP,5)=VAL2
          NWPCMP(NCMP)=0  ! start and end the same
          HCMP(NCMP,1,1)=xc; HCMP(NCMP,1,2)=yc; HCMP(NCMP,1,3)=zc
          call symbol_flow_cmp(name,'bi','r',xc,yc,zc)
          call pausems(30)

C Report on Delta Height with associated door base.
          write(outs,'(a,f7.3)') 
     &      'Delta Zone with base of door is',VAL4
          call edisp(iuout,outs)
          call pausems(30)
          if(ICT(loop).eq.0)then
            write(name,'(a,i2.2,a,i3.3)') 'BW-Bi',iz,':',is
            xw=x3+SURCOG(iz,is,1);yw=y3+SURCOG(iz,is,2)
            zw=z3+SURCOG(iz,is,3)
            write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),'  ',name,
     &      ' Wind boundary opening azim',spazi(iz,is),
     &      ' @XYZ',xw,yw,zw
            call edisp(iuout,outs)
            NNOD=NNOD+1; NDTYP(NNOD)=3; NDFLD(NNOD)=1
            HNOD(NNOD,1)=xw; HNOD(NNOD,2)=yw; HNOD(NNOD,3)=zw
            SUPNOD(NNOD,2)=spazi(iz,is)
            SUPNOD(NNOD,1)=1.0
            TNOD(NNOD)=0.0
            write(NDNAM(NNOD),'(a)') name
            NODASSOC(NNOD,1)=ZNAME(IZ); NODASSOC(NNOD,2)=SN
            call symbol_flow_nd(name,'bw','r',xw,yw,zw)
            call pausems(30)

C Establish connection from this boundary node to the bi-directional
C and then to the adjacent zone.
            POS=NNOD; NEG=ICAAS(IZ); CMPIS=NCMP; NCNN=NCNN+1
            if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
              ccol='r'
            else
              ccol='b'
            endif

C Reality checks here...
            NODPS(NCNN)=POS; NODNE(NCNN)=NEG
            HGTPS(NCNN)=0.0   ! boundary always in line with component
            HGTNE(NCNN)=VAL4   ! bi-dir assume surf COG to zone COG
            NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
            ITPCON(NCNN)=CMPIS ! use the just-defined orifice
            write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &        NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
            call edisp(iuout,outs)
            call draw_flow_link(xw,yw,zw,xc,yc,zc,ccol)
            call draw_flow_link(xc,yc,zc,
     &        HNOD(NEG,1),HNOD(NEG,2),HNOD(NEG,3),ccol)
          elseif(ICT(loop).eq.3)then

C Define connection with adjacent zone. Check if linkatpartition is still
C zero. If so mark this connection and update the other surface as well.
            if(linkatpartition(loop).eq.0)then
              POS=ICAAS(IC1(loop))  ! the current zone node is positive
              NEG=ICAAS(IC2(loop)) ! the other zone
              if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
                ccol='r'
              else
                ccol='b'
              endif
              CMPIS=NCMP; NCNN=NCNN+1
              linkatpartition(loop)=1
              icoth=IZSTOCN(IC2(loop),IE2(loop))
              linkatpartition(icoth)=1

C Reality checks here...
              NODPS(NCNN)=POS; NODNE(NCNN)=NEG
              HGTPS(NCNN)= HCMP(NCMP,1,3)-ZCOG(IZ,3)  ! diff comp & zone
              HGTNE(NCNN)= HCMP(NCMP,1,3)-ZCOG(IC2(loop),3)  ! diff comp & other zone
              NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
              ITPCON(NCNN)=CMPIS ! use the just-defined crack
              write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &          NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
              call edisp(iuout,outs)
              call draw_flow_link(ZCOG(IZ,1),ZCOG(IZ,2),ZCOG(IZ,3),
     &          HCMP(NCMP,1,1),HCMP(NCMP,1,2),HCMP(NCMP,1,3),ccol)
              call draw_flow_link(ZCOG(IC2(loop),1),ZCOG(IC2(loop),2),
     &          ZCOG(IC2(loop),3),HCMP(NCMP,1,1),HCMP(NCMP,1,2),
     &          HCMP(NCMP,1,3),ccol)
            endif
          else
            continue  ! <<  for similar creat another node ....
          endif

        endif
      elseif(SUSE(iz,is,1)(1:4).eq.'FICT')then
        if(SUSE(iz,is,2)(1:6).eq.'CLOSED'.or.
     &     SUSE(iz,is,2)(1:5).eq.'CRACK')then

C Fict crack - crack component, assume width = 1mm,
C length = window width, at lowest horizontal edge.
          write(name,'(a,i2.2,a,i3.3)') 'FiCRz',iz,':',is
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=edgelower
          write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),'  ',name,
     &      ' fict crack width ',XYMAX,' height 2mm @XYZ',xc,yc,zc
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          NCMP=NCMP+1 ! use logic in MF120I
          CMNAM(NCMP)=name
          CMPASSOC(NCMP,1)=ZNAME(IZ); CMPASSOC(NCMP,2)=SN
          ITPCMP(NCMP)=IVALCM(13); LTPCMP(NCMP)=LVALCM(13)  !crack
          ISDCMP(NCMP)=3; ISDCNN(NCMP)=0
          SUPCMP(NCMP,1)=1.
          SUPCMP(NCMP,2)=2.0/1000.0
          SUPCMP(NCMP,3)=XYMAX
          NWPCMP(NCMP)=0  ! start and end the same
          HCMP(NCMP,1,1)=xc; HCMP(NCMP,1,2)=yc; HCMP(NCMP,1,3)=zc
          call symbol_flow_cmp(name,'cr','r',xc,yc,zc)
          call pausems(30)

C Report on Delta Height with associated surf lower.
          write(outs,'(a,f7.3)') 
     &      'Delta Zone with surface lower edge is',
     &      edgelower-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call pausems(30)
          if(ICT(loop).eq.0)then
            write(name,'(a,i2.2,a,i3.3)') 'BW-Cr',iz,':',is
            xw=x3+SURCOG(iz,is,1);yw=y3+SURCOG(iz,is,2)
            zw=z3+edgelower
            write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),'  ',name,
     &      ' Wind boundary fict azim',spazi(iz,is),
     &      ' @XYZ',xw,yw,zw
            call edisp(iuout,outs)
            NNOD=NNOD+1; NDTYP(NNOD)=3; NDFLD(NNOD)=1
            HNOD(NNOD,1)=xw; HNOD(NNOD,2)=yw; HNOD(NNOD,3)=zw
            SUPNOD(NNOD,2)=spazi(iz,is)
            SUPNOD(NNOD,1)=1.0
            TNOD(NNOD)=0.0
            write(NDNAM(NNOD),'(a)') name
            NODASSOC(NNOD,1)=ZNAME(IZ); NODASSOC(NNOD,2)=SN
            call symbol_flow_nd(name,'bw','r',xw,yw,zw)
            call pausems(30)

C Establish connection from this boundary node to the crack
C and then to the adjacent zone.
            POS=NNOD; NEG=ICAAS(IZ); CMPIS=NCMP; NCNN=NCNN+1
            if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
              ccol='r'
            else
              ccol='b'
            endif

C Reality checks here...
            NODPS(NCNN)=POS; NODNE(NCNN)=NEG
            HGTPS(NCNN)=0.0   ! boundary always in line with crack
            HGTNE(NCNN)=edgelower-ZCOG(IZ,3)
            NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
            ITPCON(NCNN)=CMPIS ! use the just-defined crack
            write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd',
     &        NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
            call edisp(iuout,outs)
            call draw_flow_link(xw,yw,zw,xc,yc,zc,ccol)
            call draw_flow_link(xc,yc,zc,
     &        HNOD(NEG,1),HNOD(NEG,2),HNOD(NEG,3),ccol)
          elseif(ICT(loop).eq.3)then

C Define connection with adjacent zone. Check if linkatpartition is still
C zero. If so mark this connection and update the other surface as well.
            if(linkatpartition(loop).eq.0)then
              POS=ICAAS(IC1(loop))  ! the current zone node is positive
              NEG=ICAAS(IC2(loop))  ! the other zone
              if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
                ccol='r'
              else
                ccol='b'
              endif
              CMPIS=NCMP; NCNN=NCNN+1
              linkatpartition(loop)=1
              icoth=IZSTOCN(IC2(loop),IE2(loop))
              linkatpartition(icoth)=1

C Reality checks here...
              NODPS(NCNN)=POS; NODNE(NCNN)=NEG
              HGTPS(NCNN)= HCMP(NCMP,1,3)-ZCOG(IZ,3)  ! diff comp & zone
              HGTNE(NCNN)= HCMP(NCMP,1,3)-ZCOG(IC2(loop),3)  ! diff comp & other zone
              NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
              ITPCON(NCNN)=CMPIS ! use the just-defined crack
              write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &          NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
              call edisp(iuout,outs)
              call draw_flow_link(ZCOG(IZ,1),ZCOG(IZ,2),ZCOG(IZ,3),
     &          HCMP(NCMP,1,1),HCMP(NCMP,1,2),HCMP(NCMP,1,3),ccol)
              call draw_flow_link(ZCOG(IC2(loop),1),ZCOG(IC2(loop),2),
     &          ZCOG(IC2(loop),3),HCMP(NCMP,1,1),HCMP(NCMP,1,2),
     &          HCMP(NCMP,1,3),ccol)
            endif
          else
            continue  ! << for similar creat another node ....
          endif

        elseif(SUSE(iz,is,2)(1:4).eq.'OPEN')then

C Fict opening - common orifice component, (type 40) IVALCM(10) MF040I
C Most flow openings probably need to be a fraction of the physical
C opening to moderate flows. As user the % of surface area.
          write(name,'(a,i2.2,a,i3.3)') 'FiOPz',iz,':',is
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=SURCOG(iz,is,3)
          write(prompt,'(4a,f6.1,a,f6.1,a,f6.2,a,3f7.3)') 
     &      zsn(1:lnzsn),'  ',name,
     &      ' width ',XYMAX,' height',ZMAX,
     &      ' area',SNA(iz,is),'m^2 COG @Z',xc,yc,zc
          VAL2=25.0              ! 25% initial guess of area
          CALL EASKR(VAL2,prompt,'Area % for flow opening?',
     &      0.,'W',1.,'-',100.,'% area',IER,nbhelp)
          VAL=SNA(iz,is)*(VAL2*0.001)
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          NCMP=NCMP+1 ! use logic in MF040I
          CMNAM(NCMP)=name
          CMPASSOC(NCMP,1)=ZNAME(IZ); CMPASSOC(NCMP,2)=SN
          ITPCMP(NCMP)=IVALCM(10); LTPCMP(NCMP)=LVALCM(10) ! common orifice
          ISDCMP(NCMP)=3; ISDCNN(NCMP)=0
          SUPCMP(NCMP,1)=1.
          SUPCMP(NCMP,2)=VAL
          if(SUPCMP(NCMP,3).lt.0.001) SUPCMP(NCMP,3)=0.5
          write(prompt,'(3a)')
     &      'The discharge coef of ',CMNAM(NCMP),':'
          VAL2=SUPCMP(NCMP,3)
          CALL EASKR(VAL2,' ',prompt,
     &      0.,'W',0.,'-',10.,'discharge coef',IER,nbhelp)
          SUPCMP(NCMP,3)=VAL2
          NWPCMP(NCMP)=0  ! start and end the same
          HCMP(NCMP,1,1)=xc; HCMP(NCMP,1,2)=yc; HCMP(NCMP,1,3)=zc
          call symbol_flow_cmp(name,'or','r',xc,yc,zc)
          call pausems(30)

C Report on Delta Height with associated surf COG.
          write(outs,'(a,f7.3)') 'Delta Zone with surface COG is',
     &      SURCOG(iz,is,3)-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call pausems(30)
          if(ICT(loop).eq.0)then
            write(name,'(a,i2.2,a,i3.3)') 'BW-Op',iz,':',is
            xw=x3+SURCOG(iz,is,1);yw=y3+SURCOG(iz,is,2)
            zw=z3+SURCOG(iz,is,3)
            write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),'  ',name,
     &      ' Wind boundary opening azim',spazi(iz,is),
     &      ' @XYZ',xw,yw,zw
            call edisp(iuout,outs)
            NNOD=NNOD+1; NDTYP(NNOD)=3; NDFLD(NNOD)=1
            HNOD(NNOD,1)=xw; HNOD(NNOD,2)=yw; HNOD(NNOD,3)=zw
            SUPNOD(NNOD,2)=spazi(iz,is)
            SUPNOD(NNOD,1)=1.0; TNOD(NNOD)=0.0
            if(znotair(iz))then
              SUPNOD(NNOD,1)=2  ! mark flow node as water
              NDFLD(NNOD)=2
            else
              SUPNOD(NNOD,1)=1  ! mark flow node as air
              NDFLD(NNOD)=1
            endif
            write(NDNAM(NNOD),'(a)') name
            NODASSOC(NNOD,1)=ZNAME(IZ); NODASSOC(NNOD,2)=SN
            call symbol_flow_nd(name,'bw','r',xw,yw,zw)
            call pausems(30)

C Establish connection from this boundary node to the orifice
C and then to the adjacent zone.
            POS=NNOD; NEG=ICAAS(IZ); CMPIS=NCMP; NCNN=NCNN+1
            if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
              ccol='r'
            else
              ccol='b'
            endif

C Reality checks here...
            NODPS(NCNN)=POS; NODNE(NCNN)=NEG
            HGTPS(NCNN)=0.0   ! boundary always in line with orifice
            HGTNE(NCNN)=SURCOG(iz,is,3)-ZCOG(IZ,3)
            NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
            ITPCON(NCNN)=CMPIS ! use the just-defined crack
            write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &        NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
            call edisp(iuout,outs)
            call draw_flow_link(xw,yw,zw,xc,yc,zc,ccol)
            call draw_flow_link(xc,yc,zc,
     &        HNOD(NEG,1),HNOD(NEG,2),HNOD(NEG,3),ccol)
          elseif(ICT(loop).eq.3)then

C Define connection with adjacent zone. Check if linkatpartition is still
C zero. If so mark this connection and update the other surface as well.
            if(linkatpartition(loop).eq.0)then
              POS=ICAAS(IC1(loop))  ! the current zone node is positive
              NEG=ICAAS(IC2(loop)) ! the other zone
              if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
                ccol='r'
              else
                ccol='b'
              endif
              CMPIS=NCMP; NCNN=NCNN+1
              linkatpartition(loop)=1
              icoth=IZSTOCN(IC2(loop),IE2(loop))
              linkatpartition(icoth)=1

C Reality checks here...
              NODPS(NCNN)=POS; NODNE(NCNN)=NEG
              HGTPS(NCNN)= HCMP(NCMP,1,3)-ZCOG(IZ,3)  ! diff comp & zone
              HGTNE(NCNN)= HCMP(NCMP,1,3)-ZCOG(IC2(loop),3)  ! diff comp & other zone
              NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
              ITPCON(NCNN)=CMPIS ! use the just-defined crack
              write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &          NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
              call edisp(iuout,outs)
              call draw_flow_link(ZCOG(IZ,1),ZCOG(IZ,2),ZCOG(IZ,3),
     &          HCMP(NCMP,1,1),HCMP(NCMP,1,2),HCMP(NCMP,1,3),ccol)
              call draw_flow_link(ZCOG(IC2(loop),1),ZCOG(IC2(loop),2),
     &          ZCOG(IC2(loop),3),HCMP(NCMP,1,1),HCMP(NCMP,1,2),
     &          HCMP(NCMP,1,3),ccol)
            endif
          else
            continue  ! << for similar creat another node ....
          endif

        elseif(SUSE(iz,is,2)(1:5).eq.'BIDIR'.or.
     &         SUSE(iz,is,2)(1:10).eq.'AJAR-BIDIR')then

C Fict bi-directional - door component IVALCM(14)=130 MF130I
C Make is slightly smaller than the opening.
C width = 0.75* width, height = 0.75*height, discharge factor = 0.4
C Need to draw component at surf COG but the real base is at bottom
C edge.
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=edgelower
          if(SUSE(iz,is,2)(1:5).eq.'BIDIR')then
            VAL1=XYMAX*0.75; VAL2=ZMAX*0.75
            VAL3=ABS(ZCOG(iz,3)-edgelower)  ! base to zone COG
            VAL4=SURCOG(iz,is,3)-ZCOG(IZ,3)        ! surf COG to zone COG
            write(name,'(a,i2.2,a,i3.3)') 'FiBIz',iz,':',is
            write(outs,'(4a,f6.1,a,f6.1,a,3f7.3,a,f6.3)') 
     &      zsn(1:lnzsn),'  ',name,
     &      ' fict bidir in w= ',XYMAX,' h=',ZMAX,
     &      ' coef 0.6 base @XYZ',xc,yc,zc,
     &      ' delta base to zone COG',VAL3
          elseif(SUSE(iz,is,2)(1:10).eq.'AJAR-BIDIR')then

C Fict ajar bi-directional - door component IVALCM(14)=130 MF130I
C Make width = 0.1m height = height, discharge factor = 0.4
            VAL1=0.1; VAL2=ZMAX
            VAL3=ABS(ZCOG(iz,3)-edgelower)  ! base to zone COG
            VAL4=SURCOG(iz,is,3)-ZCOG(IZ,3)        ! surf COG to zone COG
            write(name,'(a,i2.2,a,i3.3)') 'FiBIz',iz,':',is
            write(outs,'(4a,f6.1,a,f6.1,a,3f7.3,a,f6.3)') 
     &      zsn(1:lnzsn),'  ',name,
     &      ' fict bidir in w= ',XYMAX,' h=',
     &      ZMAX,' coef 0.6 base @XYZ',xc,yc,zc,
     &      ' delta base to zone COG',VAL3
          endif
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          NCMP=NCMP+1 ! use logic in MF040I
          CMNAM(NCMP)=name
          CMPASSOC(NCMP,1)=ZNAME(IZ); CMPASSOC(NCMP,2)=SN
          ITPCMP(NCMP)=IVALCM(14); LTPCMP(NCMP)=LVALCM(14) ! bi-directional
          ISDCMP(NCMP)=5; ISDCNN(NCMP)=0
          SUPCMP(NCMP,1)=1.
          SUPCMP(NCMP,2)=VAL1
          SUPCMP(NCMP,3)=VAL2
          SUPCMP(NCMP,4)=VAL3
          if(SUPCMP(NCMP,5).lt.0.001) SUPCMP(NCMP,5)=0.5
          write(prompt,'(3a)') 'The discharge coefficient of ',
     &      CMNAM(NCMP),':'
          VAL2=SUPCMP(NCMP,5)
          CALL EASKR(VAL2,' ',prompt,0.,'W',0.,'-',10.,
     &      'discharge coef in comp',IER,nbhelp)
          SUPCMP(NCMP,5)=VAL2
          NWPCMP(NCMP)=0  ! start and end the same
          HCMP(NCMP,1,1)=xc; HCMP(NCMP,1,2)=yc; HCMP(NCMP,1,3)=zc
          call symbol_flow_cmp(name,'bi','r',xc,yc,VAL3)
          call pausems(30)

C Report on Delta Height with associated surf COG.
          write(outs,'(a,f7.3)') 'Delta Zone with surface COG is',VAL4
          call edisp(iuout,outs)
          call pausems(30)
          if(ICT(loop).eq.0)then
            write(name,'(a,i2.2,a,i3.3)') 'BW-BI',iz,':',is
            xw=x3+SURCOG(iz,is,1);yw=y3+SURCOG(iz,is,2)
            zw=z3+SURCOG(iz,is,3)
            write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),'  ',name,
     &      ' Wind boundary fict azim',spazi(iz,is),
     &      ' @XYZ',xw,yw,zw
            call edisp(iuout,outs)
            NNOD=NNOD+1; NDTYP(NNOD)=3; NDFLD(NNOD)=1
            HNOD(NNOD,1)=xw; HNOD(NNOD,2)=yw; HNOD(NNOD,3)=zw
            SUPNOD(NNOD,2)=spazi(iz,is)
            SUPNOD(NNOD,1)=1.0
            TNOD(NNOD)=0.0
            write(NDNAM(NNOD),'(a)') name
            NODASSOC(NNOD,1)=ZNAME(IZ); NODASSOC(NNOD,2)=SN
            call symbol_flow_nd(name,'bw','r',xw,yw,zw)
            call pausems(30)

C Establish connection from this boundary node to the bi-directional
C and then to the adjacent zone.
            POS=NNOD; NEG=ICAAS(IZ); CMPIS=NCMP; NCNN=NCNN+1
            if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
              ccol='r'
            else
              ccol='b'
            endif

C Reality checks here...
            NODPS(NCNN)=POS; NODNE(NCNN)=NEG
            HGTPS(NCNN)=0.0   ! boundary always in line with component
            HGTNE(NCNN)=VAL4   ! bi-dir assume surf COG to zone COG
            NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
            ITPCON(NCNN)=CMPIS ! use the just-defined orifice
            write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &        NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
            call edisp(iuout,outs)
            call draw_flow_link(xw,yw,zw,xc,yc,VAL3,ccol)
            call draw_flow_link(xc,yc,zc,
     &        HNOD(NEG,1),HNOD(NEG,2),HNOD(NEG,3),ccol)
          elseif(ICT(loop).eq.3)then

C Define connection with adjacent zone. Check if linkatpartition is still
C zero. If so mark this connection and update the other surface as well.
            if(linkatpartition(loop).eq.0)then
              POS=ICAAS(IC1(loop))  ! the current zone node is positive
              NEG=ICAAS(IC2(loop)) ! the other zone
              if(NDFLD(POS).eq.1.and.NDFLD(NEG).eq.1)then
                ccol='r'
              else
                ccol='b'
              endif
              CMPIS=NCMP; NCNN=NCNN+1
              linkatpartition(loop)=1
              icoth=IZSTOCN(IC2(loop),IE2(loop))
              linkatpartition(icoth)=1

C Reality checks here...
              NODPS(NCNN)=POS; NODNE(NCNN)=NEG
              HGTPS(NCNN)= HCMP(NCMP,1,3)-ZCOG(IZ,3)  ! diff comp & zone
              HGTNE(NCNN)= HCMP(NCMP,1,3)-ZCOG(IC2(loop),3)  ! diff comp & other zone
              NDSCNN(NCNN,1)=0; NDSCNN(NCNN,2)=0
              ITPCON(NCNN)=CMPIS ! use the just-defined crack
              write(outs,'(a,i3,6a)') 'add connection ',ncnn,' +nd ',
     &          NDNAM(POS),' to ',NDNAM(NEG),' via cmp',CMNAM(CMPIS)
              call edisp(iuout,outs)
              call draw_flow_link(ZCOG(IZ,1),ZCOG(IZ,2),ZCOG(IZ,3),
     &          HCMP(NCMP,1,1),HCMP(NCMP,1,2),HCMP(NCMP,1,3),ccol)
              call draw_flow_link(ZCOG(IC2(loop),1),ZCOG(IC2(loop),2),
     &          ZCOG(IC2(loop),3),HCMP(NCMP,1,1),HCMP(NCMP,1,2),
     &          HCMP(NCMP,1,3),ccol)
            endif
          else
            continue  ! << for similar creat another node ....
          endif
        endif
      endif

      return
      end

C Subroutine to update position of an flow node.
      subroutine update_node_pos(indexofnode,X,Y,Z,sup1,sup2)

      return
      end

C *************************
C Subroutine to update position of a flow component based on
C current surface USE attributes.
      subroutine update_cmp_pos(izone,isurf,indexofcmp)
#include "building.h"
#include "geometry.h"
#include "net_flow.h"
#include "net_flow_data.h"
#include "prj3dv.h"
      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      COMMON/C3/IC1(MCON),IE1(MCON),ICT(MCON),IC2(MCON),IE2(MCON)
      integer IZSTOCN
      COMMON/C24/IZSTOCN(MCOM,MS)
      COMMON/MFLOWIT/fndegc,imix

      integer linkatpartition  ! have we done component on side of partition
      COMMON/linkpart/linkatpartition(MCON)

      CHARACTER outs*144,name*12
      character CXSTR*78
      CHARACTER zsn*28
      logical dupedges
      integer loop
      integer llpos,lrpos,ulpos,urpos     ! closest to BB corners for parent
      logical iswater
      real edgelower,edgeupper  ! Z of current surface lower & upper edges.

C Set folder separator (fs) to \ or / as required.
      IF(izone.eq.0.or.isurf.EQ.0.or.indexofcmp.eq.0)THEN
        CALL USRMSG('Zone or surface index zero, returning.',' ','W')
        ier=1
        return
      ENDIF

C Temporarily use same file unit as profiles db.
      NDNAM(0)='    '
      fndegc=20.0
      imix=1

      PI = 4.0 * ATAN(1.0); RAD = PI/180.0;
      iz=izone; is=isurf; icmp=indexofcmp
      loop=IZSTOCN(izone,isurf)

      lsn=lnblnk(SNAME(izone,isurf))  ! setup surface name contex
      write(zsn,'(3a)') zname(iz)(1:lnzname(iz)),
     &  ':',SNAME(izone,isurf)(1:lsn)
      lnzsn=lnblnk(zsn)  ! remember length of combined zone:surface name

C Check if partition and linkatpartition already set.
C If so we do not need to create this component.
      if(ICT(loop).eq.3.and.linkatpartition(loop).gt.0)then
        return  ! no need to create anything for this surface
      endif
      name='  '
      call conxinfo(1,loop,CXSTR)

C Note: llpos, lrpos etc are the position in the iszjvn list at each
C of the corners of the surface.
      CALL ZSURLEHI(IZ,IS,XYMAX,ZMAX,llpos,lrpos,ulpos,urpos)
      if(llpos.gt.0.and.ulpos.gt.0)then
        edgelower=szcoords(iz,iszjvn(izone,isurf,llpos),3)
        edgeupper=szcoords(iz,iszjvn(izone,isurf,ulpos),3)
      endif

C Figure out the delta for the boundary node would be if
C 0.9m along normal vector.
      V1=0.9; AZ=spazi(iz,is); EL=spelv(iz,is)
      RYAZI = AZ*RAD; RSALT = EL*RAD

C Cases to ignore.
      if(SUSE(iz,is,1)(1:1).eq.'-'.or.
     &   SUSE(iz,is,1)(1:4).eq.'WALL'.or.
     &   SUSE(iz,is,1)(1:5).eq.'FLOOR'.or.
     &   SUSE(iz,is,1)(1:5).eq.'FURNI'.or.
     &   SUSE(iz,is,1)(1:7).eq.'ITEQUIP'.or.
     &   SUSE(iz,is,1)(1:5).eq.'PARTN'.or.
     &   SUSE(iz,is,1)(1:4).eq.'ROOF'.or.
     &   SUSE(iz,is,1)(1:5).eq.'STRUC'.or.
     &   SUSE(iz,is,1)(1:7).eq.'FIXTURE'.or.
     &   SUSE(iz,is,1)(1:6).eq.'PLANTS')then
        return
      endif

C For DOOR variants.
      if(SUSE(iz,is,1)(1:4).eq.'DOOR'.or.
     &   SUSE(iz,is,1)(1:6).eq.'P-DOOR'.or.
     &   SUSE(iz,is,1)(1:6).eq.'H-DOOR'.or.
     &   SUSE(iz,is,1)(1:6).eq.'H-DOOR')then
        if(SUSE(iz,is,2)(1:6).eq.'CLOSED')then

C Door crack - crack component, assume width = 1mm,
C length = perimeter, on centre line. OK
          call zsurfprm(iz,is,dupedges,perim)
          write(name,'(a,i2.2,a,i3.3)') 'DoCrz',iz,':',is
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=SURCOG(iz,is,3)
          write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),' ',name,
     &      ' door crack width 2mm len',perim,'m COG @XYZ',
     &      xc,yc,zc
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          HCMP(ICMP,1,1)=xc; HCMP(ICMP,1,2)=yc; HCMP(ICMP,1,3)=zc
          call symbol_flow_cmp(name,'cr','r',xc,yc,zc)
          call pausems(100)

C Report on Delta Height with associated surf COG.
          write(outs,'(a,f7.3)') 'Delta Zone with surface COG is',
     &      SURCOG(iz,is,3)-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call pausems(100)

        elseif(SUSE(iz,is,2)(1:8).eq.'UNDERCUT')then

C Undercut at door - common orifice, ask width = ~10mm, 
C length = door width at lowest horizontal edge. And Cd
C of 0.7 (checked against measurements in Gross & Habermann paper).
          write(name,'(a,i2.2,a,i3.3)') 'DoUcz',iz,':',is
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=edgelower
          write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),' ',name,
     &      ' door undercut ~10mm x ',XYMAX,'m @XYZ',xc,yc,zc
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          HCMP(ICMP,1,1)=xc; HCMP(ICMP,1,2)=yc; HCMP(ICMP,1,3)=zc
          call symbol_flow_cmp(name,'or','r',xc,yc,zc)
          call pausems(100)

C Report on Delta Height with associated surf lower.
          write(outs,'(a,f7.3)') 
     &      'Delta Zone with surface lower edge is',
     &      edgelower-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call pausems(100)

        elseif(SUSE(iz,is,2)(1:4).eq.'OPEN')then

C Open door - common orifice component, (type 40) IVALCM(10) MF040I
          write(name,'(a,i2.2,a,i3.3)') 'DoOPz',iz,':',is
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=SURCOG(iz,is,3)
          write(outs,'(4a,f6.1,a,f6.1,a,f7.3,a,3f6.2)') zsn(1:lnzsn),
     &      ' ',name,' door open width ',XYMAX,' height',ZMAX,
     &      ' area',SNA(iz,is),'m^2 COG @XYZ',xc,yc,zc
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          HCMP(ICMP,1,1)=xc; HCMP(ICMP,1,2)=yc; HCMP(ICMP,1,3)=zc

C Report on Delta Height with associated surf COG.
          write(outs,'(a,f7.3)') 'Delta Zone with surface COG is',
     &      SURCOG(iz,is,3)-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call symbol_flow_cmp(name,'or','r',xc,yc,zc)
          call pausems(100)

        elseif(SUSE(iz,is,2)(1:5).eq.'BIDIR'.or.
     &         SUSE(iz,is,2)(1:10).eq.'AJAR-BIDIR')then

C Bidirectional door - door component,  IVALCM(14)=130 MF130I
C assume width = door width, height = door height, discharge factor = 0.5.
C Need to draw component at surf COG but the real base is at bottom
C edge.
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=SURCOG(iz,is,3)
          if(SUSE(iz,is,2)(1:5).eq.'BIDIR')then
            write(name,'(a,i2.2,a,i3.3)') 'DoBiz',iz,':',is
            VAL3=ABS(ZCOG(iz,3)-edgelower)  ! base to zone COG
            VAL4=SURCOG(iz,is,3)-ZCOG(IZ,3)        ! surf COG to zone COG
            write(outs,'(4a,f6.1,a,f6.1,a,3f7.3,a,f6.3)') zsn(1:lnzsn),
     &      ' ',name,' door bidir width ',XYMAX,
     &      ' height',ZMAX,' coef 0.5 base @XYZ',xc,yc,zc,
     &      ' delta base to zone COG',VAL3
          elseif(SUSE(iz,is,2)(1:10).eq.'AJAR-BIDIR')then
            write(name,'(a,i2.2,a,i3.3)') 'DABiz',iz,':',is
            VAL3=ABS(ZCOG(iz,3)-edgelower)  ! base to zone COG
            VAL4=SURCOG(iz,is,3)-ZCOG(IZ,3)        ! surf COG to zone COG
            write(outs,'(4a,a,f6.1,a,3f7.3,a,f6.3)') zsn(1:lnzsn),
     &      ' ',name,' door ajar bidir width 100mm',
     &      ' height',ZMAX,' coef 0.5 base @XYZ',xc,yc,zc,
     &      ' delta base to zone COG',VAL3
          endif
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          HCMP(ICMP,1,1)=xc; HCMP(ICMP,1,2)=yc; HCMP(ICMP,1,3)=zc
          call symbol_flow_cmp(name,'bi','r',xc,yc,zc)
          call pausems(100)

C Report on Delta Height between surface COG and zone COG.
          write(outs,'(a,f7.3)') 'Delta Zone with base is',VAL4
          call edisp(iuout,outs)
          call pausems(100)
        endif
      elseif(SUSE(iz,is,1)(1:5).eq.'FRAME'.or.
     &       SUSE(iz,is,1)(1:7).eq.'F-FRAME')then

C For FRAME variants.
        if(SUSE(iz,is,2)(1:6).eq.'CLOSED'.or.
     &     SUSE(iz,is,2)(1:5).eq.'CRACK')then

C Frame crack - crack component, assume width = 1mm,
C length = frame inside perimeter, on centre line. OK
          call zsurfprm(iz,is,dupedges,perim)
          write(name,'(a,i2.2,a,i3.3)') 'FrCrz',iz,':',is
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=SURCOG(iz,is,3)
          write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),' ',name,
     &      ' frame crack width 2mm len',perim,'m COG @XYZ',
     &      xc,yc,zc
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          HCMP(ICMP,1,1)=xc; HCMP(ICMP,1,2)=yc; HCMP(ICMP,1,3)=zc
          call symbol_flow_cmp(name,'cr','r',xc,yc,zc)
          call pausems(100)

C Report on Delta Height with associated surf COG.
          write(outs,'(a,f7.3)') 'Delta Zone with surface COG is',
     &      SURCOG(iz,is,3)-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call pausems(100)

        elseif(SUSE(iz,is,2)(1:4).eq.'VENT')then

C Frame trickle vent - air flow opening, calculate width,
C height = 10mm, at highest horizontal edge.
          write(name,'(a,i2.2,a,i3.3)') 'FrVez',iz,':',is
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=edgeupper
          trwidth=XYMAX*0.6
          if(trwidth.gt.0.7) trwidth=0.7
          if(trwidth.lt.0.4.and.XYMAX.lt.0.4) trwidth=XYMAX
          write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),'  ',name,
     &      ' frame trickle vent width ',trwidth,' height 15mm @XYZ',
     &      xc,yc,zc
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          HCMP(ICMP,1,1)=xc; HCMP(ICMP,1,2)=yc; HCMP(ICMP,1,3)=zc
          call symbol_flow_cmp(name,'or','r',xc,yc,zc)
          call pausems(100)

C Report on Delta Height with associated surf upper.
          write(outs,'(a,f7.3)') 'Delta Zone with surf top is',
     &      edgeupper-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call pausems(100)
        endif
      elseif(SUSE(iz,is,1)(1:5).eq.'GRILL')then

C For GRILLs.
        if(SUSE(iz,is,2)(1:6).eq.'CLOSED'.or.
     &     SUSE(iz,is,2)(1:5).eq.'CRACK')then

C Grill crack - crack component, assume width = 1mm,
C length = frame inside perimeter, on centre line.
          call zsurfprm(iz,is,dupedges,perim)
          write(name,'(a,i2.2,a,i3.3)') 'GrCRz',iz,':',is
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=SURCOG(iz,is,3)
          write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),' ',name,
     &      ' grill crack width 2mm length',perim,'m COG @XYZ',
     &      xc,yc,zc
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          HCMP(ICMP,1,1)=xc; HCMP(ICMP,1,2)=yc; HCMP(ICMP,1,3)=zc
          call symbol_flow_cmp(name,'cr','r',xc,yc,zc)
          call pausems(100)

C Report on Delta Height with associated surf COG.
          write(outs,'(a,f7.3)') 'Delta Zone with surface COG is',
     &      SURCOG(iz,is,3)-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call pausems(100)

        elseif(SUSE(iz,is,2)(1:5).eq.'INLET')then

          write(name,'(a,i2.2,a,i3.3)') 'Grill ',iz,':',is
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=SURCOG(iz,is,3)
          write(outs,'(4a,3f7.3)') zsn(1:lnzsn),' ',name,
     &      ' grill mechanical extract COG @XYZ',xc,yc,zc
          call edisp(iuout,outs)
          HCMP(ICMP,1,1)=xc; HCMP(ICMP,1,2)=yc; HCMP(ICMP,1,3)=zc
          call ECLOSE(SUPCMP(ICMP,1),2.0,0.0001,iswater)
          if(.NOT.iswater)then
            call symbol_flow_cmp(name,'cd','r',xc,yc,zc)
          else
            call symbol_flow_cmp(name,'cd','b',xc,yc,zc)
          endif
          call pausems(50)

C Report on Delta Height with associated surf COG.
          write(outs,'(a,f7.3)') 'Delta Zone with surface COG is',
     &      SURCOG(iz,is,3)-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call pausems(50)

        elseif(SUSE(iz,is,2)(1:7).eq.'EXTRACT')then

C Grill extract - constant vol flow comp. (type 30) IVALCM(8), as in MF030I
          write(name,'(a,i2.2,a,i3.3)') 'GrEXz',iz,':',is
          call zsurfprm(iz,is,dupedges,perim)
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=SURCOG(iz,is,3)
          write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),' ',name,
     &      ' grill extract hydr diam',
     &      4.0*SNA(iz,is)/perim,'m COG @XYZ',xc,yc,zc
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          HCMP(ICMP,1,1)=xc; HCMP(ICMP,1,2)=yc; HCMP(ICMP,1,3)=zc

          call ECLOSE(SUPCMP(ICMP,1),2.0,0.0001,iswater)
          if(.NOT.iswater)then
            call symbol_flow_cmp(name,'bx','r',xc,yc,zc)
          else
            call symbol_flow_cmp(name,'bx','b',xc,yc,zc)
          endif
          call pausems(30)

C Report on Delta Height with associated surf COG.
          write(outs,'(a,f7.3)') 'Delta Zone with surface COG is',
     &      SURCOG(iz,is,3)-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call pausems(30)

        elseif(SUSE(iz,is,2)(1:4).eq.'OPEN')then

C Grill opening - common orifice component, (type 40) IVALCM(10) MF040I
          call zsurfprm(iz,is,dupedges,perim)
          write(name,'(a,i2.2,a,i3.3)') 'GrOPz',iz,':',is
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=SURCOG(iz,is,3)
          write(outs,'(4a,f6.1,a,3f6.2)') zsn(1:lnzsn),' ',name,
     &      ' grill orifice area',SNA(iz,is),'m^2 COG @XYZ',
     &      xc,yc,zc
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          HCMP(ICMP,1,1)=xc; HCMP(ICMP,1,2)=yc; HCMP(ICMP,1,3)=zc
          call symbol_flow_cmp(name,'or','r',xc,yc,zc)
          call pausems(30)

C Report on Delta Height with associated surf COG.
          write(outs,'(a,f7.3)') 'Delta Zone with surface COG is',
     &      SURCOG(iz,is,3)-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call pausems(30)

        endif
      elseif(SUSE(iz,is,1)(1:6).eq.'WINDOW'.or.
     &    SUSE(iz,is,1)(1:8).eq.'D-WINDOW'.or.
     &    SUSE(iz,is,1)(1:8).eq.'S-WINDOW'.or.
     &    SUSE(iz,is,1)(1:8).eq.'C-WINDOW')then
        if(SUSE(iz,is,2)(1:6).eq.'CLOSED'.or.
     &     SUSE(iz,is,2)(1:5).eq.'CRACK')then

C Window crack - crack component, assume width = 1mm,
C length = window perimeter, at surface COG.
          call zsurfprm(iz,is,dupedges,perim)
          write(name,'(a,i2.2,a,i3.3)') 'WiCrz',iz,':',is
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=SURCOG(iz,is,3)
          write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),' ',name,
     &      ' window crack width 2mm length',XYMAX,' @XYZ',xc,yc,zc
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          HCMP(ICMP,1,1)=xc; HCMP(ICMP,1,2)=yc; HCMP(ICMP,1,3)=zc
          call symbol_flow_cmp(name,'cr','r',xc,yc,zc)
          call pausems(30)

C Report on Delta Height with associated surf lower.
          write(outs,'(a,f7.3)') 
     &      'Delta Zone with surface lower edge is',
     &      SURCOG(iz,is,3)-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call pausems(30)

        elseif(SUSE(iz,is,2)(1:4).eq.'OPEN')then

C Window open - common orifice component, (type 40) IVALCM(10) MF040I
          write(name,'(a,i2.2,a,i3.3)') 'WiOpz',iz,':',is
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=SURCOG(iz,is,3)
          write(outs,'(4a,f6.1,a,f6.1,a,f6.2,a,3f7.3)') 
     &      zsn(1:lnzsn),' ',name,
     &      ' window within surface ',SNA(iz,is),'m^2 COG @Z',
     &      xc,yc,zc
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          HCMP(ICMP,1,1)=xc; HCMP(ICMP,1,2)=yc; HCMP(ICMP,1,3)=zc
          call symbol_flow_cmp(name,'or','r',xc,yc,zc)
          call pausems(30)

C Report on Delta Height with associated surf COG.
          write(outs,'(a,f7.3)') 'Delta Zone with surface COG is',
     &      SURCOG(iz,is,3)-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call pausems(30)

        elseif(SUSE(iz,is,2)(1:4).eq.'SASH')then

C Window sash - 2 common orifice component, (type 40) IVALCM(10) MF040I
          write(name,'(a,i2.2,a,i3.3)') 'WiSLz',iz,':',is
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=edgelower
          write(outs,'(4a,f6.1,a,f6.1,a,f6.2,a,3f6.2)') 
     &      zsn(1:lnzsn),' ',name,
     &      ' within surface with width ',XYMAX,' height',ZMAX,
     &      ' area',SNA(iz,is),'m2 sill @XYZ',xc,yc,zc
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          HCMP(ICMP,1,1)=xc; HCMP(ICMP,1,2)=yc; HCMP(ICMP,1,3)=zc
          call symbol_flow_cmp(name,'or','r',xc,yc,zc)
          call pausems(30)

C Report on Delta Height with associated surf lower.
          write(outs,'(a,f7.3)') 
     &      'Delta Zone with surface lower edge is',
     &      edgelower-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call pausems(30)

C Now the upper orifice.
          write(name,'(a,i2.2,a,i3.3)') 'WiSUz',iz,':',is
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=edgeupper
          write(outs,'(4a,f6.1,a,f6.1,a,f6.2,a,3f6.2)') 
     &      zsn(1:lnzsn),' ',name,
     &      ' within surface with width ',XYMAX,' height',ZMAX,
     &      ' area',SNA(iz,is),'m2 upper @XYZ',xc,yc,zc
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          HCMP(ICMP,1,1)=xc; HCMP(ICMP,1,2)=yc; HCMP(ICMP,1,3)=zc
          call symbol_flow_cmp(name,'or','r',xc,yc,zc)
          call pausems(30)

C Report on Delta Height with associated surf upper.
          write(outs,'(a,f7.3)') 
     &      'Delta Zone with head of upper sash is',
     &      edgeupper-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call pausems(30)

        elseif(SUSE(iz,is,2)(1:5).eq.'BIDIR'.or.
     &         SUSE(iz,is,2)(1:10).eq.'AJAR-BIDIR')then

C Window bidirectional - specific door component width & height
C based on the surface width & height. discharge factor = 0.6,
C Draw on the centre line. The VAL4 is Surf COG to Zone COG.
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=SURCOG(iz,is,3)
          if(SUSE(iz,is,2)(1:5).eq.'BIDIR')then
            write(name,'(a,i2.2,a,i3.3)') 'WiBiz',iz,':',is
            VAL3=ABS(ZCOG(iz,3)-edgelower)  ! base to zone COG
            VAL4=SURCOG(iz,is,3)-ZCOG(IZ,3)        ! surf COG to zone COG
            write(outs,'(4a,f6.1,a,f6.1,a,3f7.3,a,f6.3)')
     &      zsn(1:lnzsn),' ',name,' bidir width ',XYMAX,' height',ZMAX,
     &      ' coef 0.6 @XYZ',xc,yc,zc,
     &      ' delta base to zone COG',VAL3
          elseif(SUSE(iz,is,2)(1:10).eq.'AJAR-BIDIR')then
            write(name,'(a,i2.2,a,i3.3)') 'WABiz',iz,':',is
            VAL3=ABS(ZCOG(iz,3)-edgelower)  ! base to zone COG
            VAL4=SURCOG(iz,is,3)-ZCOG(IZ,3)        ! surf COG to zone COG
            write(outs,'(4a,a,f6.1,a,3f7.3,a,f6.3)')
     &      zsn(1:lnzsn),' ',name,' bidir width 100mm',' height',ZMAX,
     &      ' coef 0.6 @XYZ',xc,yc,zc,
     &      ' delta base to zone COG',VAL3
          endif
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          HCMP(ICMP,1,1)=xc; HCMP(ICMP,1,2)=yc; HCMP(ICMP,1,3)=zc
          call symbol_flow_cmp(name,'bi','r',xc,yc,zc)
          call pausems(30)

C Report on Delta Height with associated surf lower.
          write(outs,'(a,f7.3)') 
     &      'Delta Zone with sill of lower sash is',VAL4
          call edisp(iuout,outs)
          call pausems(30)
        endif
      elseif(SUSE(iz,is,1)(1:4).eq.'FICT')then
        if(SUSE(iz,is,2)(1:6).eq.'CLOSED'.or.
     &     SUSE(iz,is,2)(1:5).eq.'CRACK')then

C Fict crack - crack component, assume width = 1mm,
C length = window width, at lowest horizontal edge.
          write(name,'(a,i2.2,a,i3.3)') 'FiCRz',iz,':',is
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=edgelower
          write(outs,'(4a,f6.1,a,3f7.3)') zsn(1:lnzsn),' ',name,
     &      ' fict crack width ',XYMAX,' height 2mm @XYZ',xc,yc,zc
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          HCMP(ICMP,1,1)=xc; HCMP(ICMP,1,2)=yc; HCMP(ICMP,1,3)=zc
          call symbol_flow_cmp(name,'cr','r',xc,yc,zc)
          call pausems(30)

C Report on Delta Height with associated surf lower.
          write(outs,'(a,f7.3)') 
     &      'Delta Zone with surface lower edge is',
     &      edgelower-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call pausems(30)

        elseif(SUSE(iz,is,2)(1:4).eq.'OPEN')then

C Fict opening - common orifice component, (type 40) IVALCM(10) MF040I
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=SURCOG(iz,is,3)
          write(name,'(a,i2.2,a,i3.3)') 'FiOPz',iz,':',is
          write(outs,'(4a,f6.1,a,f6.1,a,f6.2,a,3f7.3)') 
     &      zsn(1:lnzsn),' ',name,
     &      ' within surface with width ',XYMAX,' height',ZMAX,
     &      ' area',SNA(iz,is),'m^2 COG @Z',xc,yc,zc
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          HCMP(ICMP,1,1)=xc; HCMP(ICMP,1,2)=yc; HCMP(ICMP,1,3)=zc
          call symbol_flow_cmp(name,'or','r',xc,yc,zc)
          call pausems(30)

C Report on Delta Height with associated surf COG.
          write(outs,'(a,f7.3)') 'Delta Zone with surface COG is',
     &      SURCOG(iz,is,3)-ZCOG(IZ,3)
          call edisp(iuout,outs)
          call pausems(30)

        elseif(SUSE(iz,is,2)(1:5).eq.'BIDIR'.or.
     &         SUSE(iz,is,2)(1:10).eq.'AJAR-BIDIR')then

C Fict bi-directional - door component IVALCM(14)=130 MF130I
C Draw component at surf COG but the real base is at bottom edge.
          xc=SURCOG(iz,is,1); yc=SURCOG(iz,is,2) ! remember
          zc=edgelower
          if(SUSE(iz,is,2)(1:5).eq.'BIDIR')then
            VAL3=ABS(ZCOG(iz,3)-edgelower)  ! base to zone COG
            VAL4=SURCOG(iz,is,3)-ZCOG(IZ,3)        ! surf COG to zone COG
            write(name,'(a,i2.2,a,i3.3)') 'FiBIz',iz,':',is
            write(outs,'(4a,f6.1,a,f6.1,a,3f7.3,a,f6.3)') 
     &      zsn(1:lnzsn),' ',name,
     &      ' fict bidir with surface with width ',XYMAX,' height',ZMAX,
     &      ' coef 0.6 base @XYZ',xc,yc,zc,
     &      ' delta base to zone COG',VAL3
          elseif(SUSE(iz,is,2)(1:10).eq.'AJAR-BIDIR')then
            VAL3=ABS(ZCOG(iz,3)-edgelower)  ! base to zone COG
            VAL4=SURCOG(iz,is,3)-ZCOG(IZ,3)        ! surf COG to zone COG
            write(name,'(a,i2.2,a,i3.3)') 'FABIz',iz,':',is
            write(outs,'(4a,f6.1,a,3f7.3,a,f6.3)') 
     &      zsn(1:lnzsn),' ',name,
     &      ' fict bidir width 100mm within surface with height',
     &      ZMAX,' coef 0.6 base @XYZ',xc,yc,zc,
     &      ' delta base to zone COG',VAL3
          endif
          call edisp(iuout,' ')
          call edisp(iuout,outs)
          HCMP(ICMP,1,1)=xc; HCMP(ICMP,1,2)=yc; HCMP(ICMP,1,3)=zc
          call symbol_flow_cmp(name,'bi','r',xc,yc,VAL3)
          call pausems(30)

C Report on Delta Height with associated surf COG.
          write(outs,'(a,f7.3)') 'Delta Zone with surface COG is',VAL4
          call edisp(iuout,outs)
          call pausems(30)
        endif
      endif

      return
      end ! of update_cmp_pos

C Subroutine to update position of a flow connection.
      subroutine update_cnn_pos(indexofcnn)

      return
      end

C *************************
C isuseflowrelated tests whether surface attributes are flow related.
      subroutine isuseflowrelated(izone,isurf,isrelated)
#include "building.h"
#include "geometry.h"
#include "net_flow.h"
#include "net_flow_data.h"

C Parameters passed.
      integer izone,isurf  ! the zone and surface index
      logical isrelated    ! true is USE associated with flow

      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      integer IZSTOCN
      COMMON/C24/IZSTOCN(MCOM,MS)
      COMMON/MFLOWIT/fndegc,imix

      character CXSTR*78
      integer loop

C Set folder separator (fs) to \ or / as required.
      IF(izone.eq.0.or.isurf.EQ.0)THEN
        CALL USRMSG('Zone or surface index zero, returning.',' ','W')
        isrelated=.false.
        return
      ENDIF

C Temporarily use same file unit as profiles db.
      isrelated=.false.
      NDNAM(0)='    '
      fndegc=20.0
      imix=1

      loop=IZSTOCN(izone,isurf)
      call conxinfo(1,loop,CXSTR)

C Cases to ignore.
      if(SUSE(izone,isurf,1)(1:1).eq.'-'.or.
     &   SUSE(izone,isurf,1)(1:4).eq.'WALL'.or.
     &   SUSE(izone,isurf,1)(1:5).eq.'FLOOR'.or.
     &   SUSE(izone,isurf,1)(1:5).eq.'FURNI'.or.
     &   SUSE(izone,isurf,1)(1:7).eq.'ITEQUIP'.or.
     &   SUSE(izone,isurf,1)(1:5).eq.'PARTN'.or.
     &   SUSE(izone,isurf,1)(1:4).eq.'ROOF'.or.
     &   SUSE(izone,isurf,1)(1:5).eq.'STRUC'.or.
     &   SUSE(izone,isurf,1)(1:7).eq.'FIXTURE'.or.
     &   SUSE(izone,isurf,1)(1:6).eq.'PLANTS')then
        isrelated=.false.
        return
      endif

C For DOOR variants.
      if(SUSE(izone,isurf,1)(1:4).eq.'DOOR'.or.
     &   SUSE(izone,isurf,1)(1:6).eq.'P-DOOR'.or.
     &   SUSE(izone,isurf,1)(1:6).eq.'H-DOOR'.or.
     &   SUSE(izone,isurf,1)(1:6).eq.'H-DOOR')then
        if(SUSE(izone,isurf,2)(1:6).eq.'CLOSED')then
          isrelated=.true.
        elseif(SUSE(izone,isurf,2)(1:8).eq.'UNDERCUT')then
          isrelated=.true.
        elseif(SUSE(izone,isurf,2)(1:4).eq.'OPEN')then
          isrelated=.true.
        elseif(SUSE(izone,isurf,2)(1:5).eq.'BIDIR'.or.
     &         SUSE(izone,isurf,2)(1:10).eq.'AJAR-BIDIR')then
          isrelated=.true.
        endif
      elseif(SUSE(izone,isurf,1)(1:5).eq.'FRAME'.or.
     &       SUSE(izone,isurf,1)(1:7).eq.'F-FRAME')then

C For FRAME vairants.
        if(SUSE(izone,isurf,2)(1:6).eq.'CLOSED')then
          isrelated=.true.
        elseif(SUSE(izone,isurf,2)(1:5).eq.'CRACK')then
          isrelated=.true.
        elseif(SUSE(izone,isurf,2)(1:4).eq.'VENT')then
          isrelated=.true.
        endif
      elseif(SUSE(izone,isurf,1)(1:5).eq.'GRILL')then

C For GRILLs.
        if(SUSE(izone,isurf,2)(1:6).eq.'CLOSED')then
          isrelated=.true.
        elseif(SUSE(izone,isurf,2)(1:5).eq.'CRACK')then
          isrelated=.true.
        elseif(SUSE(izone,isurf,2)(1:5).eq.'INLET')then
          isrelated=.true.
        elseif(SUSE(izone,isurf,2)(1:7).eq.'EXTRACT')then
          isrelated=.true.
        elseif(SUSE(izone,isurf,2)(1:4).eq.'OPEN')then
          isrelated=.true.
        endif
      elseif(SUSE(izone,isurf,1)(1:6).eq.'WINDOW'.or.
     &    SUSE(izone,isurf,1)(1:8).eq.'D-WINDOW'.or.
     &    SUSE(izone,isurf,1)(1:8).eq.'S-WINDOW'.or.
     &    SUSE(izone,isurf,1)(1:8).eq.'C-WINDOW')then
        if(SUSE(izone,isurf,2)(1:6).eq.'CLOSED')then
          isrelated=.true.
        elseif(SUSE(izone,isurf,2)(1:5).eq.'CRACK')then
          isrelated=.true.
        elseif(SUSE(izone,isurf,2)(1:4).eq.'OPEN')then
          isrelated=.true.
        elseif(SUSE(izone,isurf,2)(1:4).eq.'SASH')then
          isrelated=.true.
        elseif(SUSE(izone,isurf,2)(1:5).eq.'BIDIR'.or.
     &         SUSE(izone,isurf,2)(1:10).eq.'AJAR-BIDIR')then
          isrelated=.true.
        endif
      elseif(SUSE(izone,isurf,1)(1:4).eq.'FICT')then
        if(SUSE(izone,isurf,2)(1:6).eq.'CLOSED')then
          isrelated=.true.
        elseif(SUSE(izone,isurf,2)(1:5).eq.'CRACK')then
          isrelated=.true.
        elseif(SUSE(izone,isurf,2)(1:4).eq.'OPEN')then
          isrelated=.true.
        elseif(SUSE(izone,isurf,2)(1:5).eq.'BIDIR'.or.
     &         SUSE(izone,isurf,2)(1:10).eq.'AJAR-BIDIR')then
          isrelated=.true.
        endif
      endif

      return
      end ! of isuseflowrelated

C Subroutine to compare position of a flow connection and its associated
C surface. Called from flow facilities as a reality check.
      subroutine compare_cmp_pos(indexofcmp,iaz,ias,tdis)
#include "building.h"
#include "geometry.h"
#include "net_flow.h"
#include "net_flow_data.h"
#include "prj3dv.h"
      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      integer ncomp,ncon
      COMMON/C1/NCOMP,NCON

C Passed parameters:
      integer indexofcmp  ! passed in component index
      integer iaz,ias     ! index of matching zone and surface

      CHARACTER msg*124,zn*12,sn*12
      integer loop,loop2

      iaz=0; ias=0  ! Zero returned values
      if(indexofcmp.eq.0) return
      write(SN,'(a)') CMPASSOC(indexofcmp,2)  ! associated surface
      write(ZN,'(a)') CMPASSOC(indexofcmp,1)  ! associated zone
      if(ZN(1:1).eq.'-'.and.SN(1:1).eq.'-') return

C Use string matching to find the zone and surface index for the
C associated zone and surface names associated with this component.
      lnczn=lnblnk(ZN)
      lncsn=lnblnk(SN)
      do loop=1,NCOMP
        if(zname(loop)(1:lnzname(loop)).eq.ZN(1:lnczn))then
          do loop2=1,NZSUR(loop)
            lsn=lnblnk(SNAME(loop,loop2))
            if(SN(1:lncsn).eq.SNAME(loop,loop2)(1:lsn))then
              iaz=loop; ias=loop2    ! to pass back
              xc=SURCOG(loop,loop2,1); yc=SURCOG(loop,loop2,2) ! remember
              zc=SURCOG(loop,loop2,3)
              tdis= crowxyz(xc,yc,zc,HCMP(indexofcmp,1,1),
     &          HCMP(indexofcmp,1,2),HCMP(indexofcmp,1,3))
              if(tdis.gt.1.0)then
                write(msg,'(a,3f7.3,a,3f7.3,a,f7.3)')
     &           'The distance between component @ ',
     &           HCMP(indexofcmp,1,1),
     &           HCMP(indexofcmp,1,2),HCMP(indexofcmp,1,3),
     &           ' and surface COG @ ',xc,yc,zc,' is ',
     &           tdis
                call edisp(iuout,msg)
              endif
              return  ! found match return indices and distance.
            endif
          enddo
        endif
      enddo
      return
      end
