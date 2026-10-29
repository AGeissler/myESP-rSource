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


C The file emeta.f is a collection of support facilities for
C reading and creating esp-r model meta files:
C  silentmodel creates a model from passed parameters.
C  silentzone: creates a zone based on parameters
C              passed to it.
C  silentread: Reads esp-r meta file zones geometry data.
C  metawrite:  Creates an esp-r meta file. Data input is from
C              the current zones common blocks.

C ************* silentmodel
C Creates a model (registration level) from passed parameters.
C The parameters are:
C action (8 char) as follows:
C    'new' creates model folders based on information
C       in root and mpath parameters.
C    'within' assumes that the current folder will be used.
C root   (32 char) the root name of the model
C mpath  (72 char) the path to the folder with cfg file.
C weather (32 char) the weather file name (no path included)
C simact (6 char) assessments to be carried out '------' if
C   none requested, otherwise follow pattern in edipv.F.

C Current functionality is to update the weather & site data.

      subroutine silentmodel(action,root,mpath,weather,simact,ier)
#include "building.h"
#include "model.h"
#include "site.h"
#include "espriou.h"
#include "esprdbfile.h"
#include "seasons.h"
#include "ipvdata.h"

      integer lnblnk  ! function definition

      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      common/FILEP/IFIL

C Where ESP-r was installed (as recorded when it was compiled).
      common/deflt4/dinstpath
      character dinstpath*60

C Weather data.
      COMMON/CLMSET/ICYEAR,ICDNGH,CLAT,CLONG

      COMMON/SET1/IYEAR,IBDOY,IEDOY,IFDAY,IFTIME
      INTEGER :: IYEAR,IBDOY,IEDOY,IFDAY,IFTIME

      common/XMLDATA/XMLMETRICS(200),inumXMLmetrics
      character XMLMETRICS*248

C Simulation parameter sets.
      common/spfldat/nsset,isset,isstup,isbnstep,ispnstep,issave,isavgh
      INTEGER :: nsset,isset,isstup,isbnstep,ispnstep,issave,isavgh

C Passed in parameters.
      character action*8,root*32,mpath*72,menu*72
      character fs*1           ! file separator
      character weather*32     ! file name for the location (no path)
      character subpath*72
      character simact*6 ! action for creating assessments.

      character outs*124,OUTSTR*124
      logical unixok,XST,ok
      logical haveground ! set true if ground temp file has been found
      character inpxmlfl*144 !to write out the input.xml file
      character llclmdb*144
      integer llt,lndbp

      IER=0

C Determine operating system.
      call isunix(unixok)
      if(unixok)then
        fs = char(47)
      else
        fs = char(92)
      endif
      subpath=' '
      haveground=.false.

C If the action is `new` then setup model folder(s) and the registration
C level configuration file via call to pregist in `sn` mode.
C If the action is `within` do not create new folders for the model.
      call edisp(iuout,'  ')
      if(action(1:3).eq.'new')then
        call edisp(iuout,'You will be asked a [YES RESTART] question.')
        call edisp(iuout,'Click on the [CONTINUE] button instead!')
        write(menu,'(a)') modeltitle(1:lnblnk(modeltitle))
        call pregist('sn',root,mpath,menu,ier)
      elseif(action(1:6).eq.'within')then

C Assume the model folders exist when calling pregist.
        write(menu,'(a)') modeltitle(1:lnblnk(modeltitle))
        call pregist('sw',root,mpath,menu,ier)
      endif

C If there is a line with the root name of a climate file look
C for a match in the current climate location and if there is
C one set the name of that file, scan it and find its site. If
C there is no weather file then set the default seasons.
      if(weather(1:4).eq.'none')then

C Set default early winter, spring, summer, autumn, late winter periods.
        CALL EDAY(9,1,ia1wins);   CALL EDAY(15,1,ia1winf)
        CALL EDAY(6,3,ia1sprs);   CALL EDAY(12,3,ia1sprf)
        CALL EDAY(11,7,iasums);   CALL EDAY(17,7,iasumf)
        CALL EDAY(2,10,ia2sprs);  CALL EDAY(8,10,ia2sprf)
        CALL EDAY(20,11,ia2wins); CALL EDAY(26,11,ia2winf)

C Default season definitions.
        CALL EDAY(1,1,is1wins);   CALL EDAY(28,2,is1winf)
        CALL EDAY(1,11,is2wins);  CALL EDAY(31,12,is2winf)
        CALL EDAY(1,3,is1sprs);   CALL EDAY(30,4,is1sprf)
        CALL EDAY(1,9,is2sprs);   CALL EDAY(31,10,is2sprf)
        CALL EDAY(1,5,is1sums);   CALL EDAY(31,8,is1sumf)
      else

C Assume that weather is in the standard location, expand the file
C name and check if it exists via the usual pattern. If not
C check if the climate file exists via findwhichdbpath.
        if(ICLIM.eq.0)then
          ICLIM=IFIL
        endif
        llt=lnblnk(weather)
        write(LCLIM,'(a)') weather(1:llt) ! Set LCLIM
        call CLMOPB(0,ISTAT)
        if(ISTAT.ge.0)then
          call CLMRDBMD(IER)              ! Read site and year.
          CALL ERPFREE(ICLIM,ISTAT)
          ipathclim=2
        else

C Use findwhichdbpath and then expand weather file path for checking.
          call findwhichdbpath('clm',weather,ier)
          llt=lnblnk(LCLIM)
          lndbp=lnblnk(standardclmpath)
          if(ipathclim.eq.0.or.ipathclim.eq.1)then
            llclmdb=LCLIM
          elseif(ipathclim.eq.2)then
            write(llclmdb,'(3a)') standardclmpath(1:lndbp),fs,
     &        LCLIM(1:lnblnk(LCLIM))
          endif
          CALL ERPFREE(ICLIM,ISTAT)
          call CLMOPB(0,ISTAT)
          write(currentfile,'(a)') llclmdb(1:lnblnk(llclmdb))
          if(ISTAT.ge.0)then
            call CLMRDBMD(IER)            ! Read site and year.
            CALL ERPFREE(ICLIM,ISTAT)
          else

C Not found, use the default weather file as a backstop.
            write(currentfile,'(a)') DCLIM(1:lnblnk(DCLIM))
            LCLIM=DCLIM
            CALL ERPFREE(ICLIM,ISTAT)
            call CLMOPB(0,ISTAT)
            if(ISTAT.ge.0)then
              call CLMRDBMD(IER)
              CALL ERPFREE(ICLIM,ISTAT)
              LCLIM=DCLIM
              ipathclim=2  ! assume default weather is in standard folder
            endif
          endif
        endif

C Take year (IYEAR) from the META file.
        sitelat=CLAT
        sitelongdif=CLONG
        write(outs,'(A,F6.1,A,F5.1,a,i4,a,2F7.3,a)')
     &    'The weather lat is ',clat,' long diff is ',
     &    clong,' Year is ',IYEAR,
     &    '. Previously used ',sitelat,sitelongdif,'.'
        call edisp(iuout,outs)

C If there is not a model specific season definition rescan
C the `climatelist` file. Check if this climate is in the
C list. If not, instantiate season and typical start and end dates.
        if(ihaveseason.gt.0)then
          continue
        else
          INQUIRE (FILE=cdblfil,EXIST=XST)
          if(XST)then
            IUF=IFIL+2
            call scancdblist(IUF,llclmdb,'p',ok,ier)
            if(ok)then
              continue
            else

C Set default early winter, spring, summer, autumn, late winter periods.
              CALL EDAY(9,1,ia1wins);   CALL EDAY(15,1,ia1winf)
              CALL EDAY(6,3,ia1sprs);   CALL EDAY(12,3,ia1sprf)
              CALL EDAY(11,7,iasums);   CALL EDAY(17,7,iasumf)
              CALL EDAY(2,10,ia2sprs);  CALL EDAY(8,10,ia2sprf)
              CALL EDAY(20,11,ia2wins); CALL EDAY(26,11,ia2winf)

C Default season definitions.
              CALL EDAY(1,1,is1wins);   CALL EDAY(28,2,is1winf)
              CALL EDAY(1,11,is2wins);  CALL EDAY(31,12,is2winf)
              CALL EDAY(1,3,is1sprs);   CALL EDAY(30,4,is1sprf)
              CALL EDAY(1,9,is2sprs);   CALL EDAY(31,10,is2sprf)
              CALL EDAY(1,5,is1sums);   CALL EDAY(31,8,is1sumf)
            endif
          endif
        endif
      endif

C Fall through if no climate located.
  74  if(nsset.gt.0)then

C Instantiate the IPV data based on value of simact.
        call ipvdatinit(simact)

        if(nipvassmt.eq.1.or.nipvassmt.eq.3.or.nipvassmt.eq.5)then

C If one of the standard number of IPV assessments transfer the relevant
C data into the simulation parameter set data structure.
          call ipv2simpar(simact)
        endif
      endif
  
      CALL EMKCFG('-',IER)

C Also write the xml input file if there is a metric defined in the flat file
      if(inumXMLmetrics.gt.0)then
        if(unixok)then
          fs = char(47)
        else
          fs = char(92)
        endif

C Write the XML input file in the same location as the cfg file.
        write(inpxmlfl,'(a)') 'input.xml'
        iuf=ifil+2
        call efopseq(iuf,inpxmlfl,3,ier)
        write(iuf,'(a)',iostat=ios) 
     &    '<?xml version="1.0" encoding="UTF-8"?>'
        write(iuf,'(a)',iostat=ios)'<configuration>'
        write(iuf,'(a)',iostat=ios) 
     &    '  <apply_style_sheet>false</apply_style_sheet>'
        write(iuf,'(a)',iostat=ios) 
     &    '  <dump_all_data>false</dump_all_data>'  
        write(iuf,'(a)',iostat=ios) 
     &    '  <hierarchy>flat</hierarchy>'  
        write(iuf,'(a)',iostat=ios) 
     &    '  <link_style_sheet>false</link_style_sheet>'  
        write(iuf,'(a)',iostat=ios) 
     &    '  <output_dictionary>true</output_dictionary>'  
        write(iuf,'(a)',iostat=ios) 
     &'  <report_startup_period_data>false</report_startup_period_data>'
        write(iuf,'(a)',iostat=ios) 
     &    '  <save_to_disk>false</save_to_disk>'  
        write(iuf,'(a)',iostat=ios) 
     &    '  <time_step_averaging>true</time_step_averaging>'  
        write(iuf,'(a)',iostat=ios) 
     &    '  <enable_xml_wildcards>true</enable_xml_wildcards>'

C Loop around the possible defined xml outputs
        do 200 iwrxmlout=1,inumXMLmetrics
        write(iuf,'(3a)',iostat=ios)'  <step_variable>', 
     &    XMLMETRICS(iwrxmlout)(1:lnblnk(XMLMETRICS(iwrxmlout))),
     &   '</step_variable>'
 200    continue       
        write(iuf,'(a)',iostat=ios)'</configuration>'
        call erpfree(iuf,istat)
      endif

C Other logic here....

      return

C Error trap. Set ier as two.
 1000 WRITE(outs,774)ISTAT
  774 FORMAT(' Error ',I7,' reading Climate db information.')
      call edisp(iuout,outs)
      call edisp(iuout,' ')
      ier=2
      RETURN

      end

C ************* silentzone
C silentzone creates a zone based on meta parameters passed
C to it.
C hasconstr (string array) the name of the MLC to use for each surface
C hasoptic (string array) the name of the optics to use for each surface
C use integer*3 directives for zone use pattern.,
c usefile char*32 file name (w/o path) of pattern operation file.
C IER=0 OK.
      subroutine silentzone(ICOMP,metaver,hasconstr,hasoptic,hasparent,
     &  hasuse,use,usefile,literalctl,IER)
#include "building.h"
#include "model.h"
#include "geometry.h"
#include "esprdbfile.h"
#include "material.h"
#include "schedule.h"

      integer lnblnk  ! function definition

C Parameters
      integer icomp    ! focus zone index
      integer metaver  ! one for 1.1 two for 1.2 three for 1.3
      character hasconstr*32
      dimension hasconstr(MS)  ! construction to use for each surface
      character hasoptic*32
      dimension hasoptic(MS)   ! optics to use for each surface
      character hasparent*12
      dimension hasparent(MS)
      character hasuse*12
      dimension hasuse(MS,2)
      character use*8,usefile*32  ! usage pattern directives
      dimension use(3)

      common/caleni/nbdaytype,nbcaldays(MDTY),icalender(365)
      INTEGER NBDAYTYPE,NBCALDAYS,ICALENDER
      COMMON/FILEP/IFIL
      COMMON/OUTIN/IUOUT,IUIN,IEOUT

      integer ncomp,ncon
      COMMON/C1/NCOMP,NCON
      COMMON/C3/IC1(MCON),IE1(MCON),ICT(MCON),IC2(MCON),IE2(MCON)
      integer IZSTOCN
      COMMON/C24/IZSTOCN(MCOM,MS)
      common/user/browse
      common/cctlnm/ctldoc,lctlf
      integer icascf
      common/cctl/icascf(mcom)
      common/deflt4/dinstpath
      character dinstpath*60 ! install folder as recorded when ESP-r was compiled
      COMMON/OVRWT/AUTOVR
      logical AUTOVR
      common/SLNTFULL/ifullysilent
      integer ifullysilent !flag to drive the silent model creation without
                           !questions about the operation files data

C hasname (string array) the surface name for each surface (so filsur does not wipe)
      character hasname*12
      dimension hasname(MS)  ! name to use for each surface
      DIMENSION XX(MS),YY(MS)
      character fs*1
      character ctldoc*248,LCTLF*72,zd*64
      character DFILE*72,CFILE*72,OFILE*72,TFILE*72
      character act*2
      character subpath*84,sfile*144
      character icact*8,iwatt*8,afact*8,iflow*8,vflow*8
      character sn*12,constr*32,opt*12,outs*124
      LOGICAL browse,OKC,XST,unixok
      logical willneedtmc  !  if a tmc file will be required.
      logical closev       !  mass is close to vertical
      LOGICAL QUIET
      logical anunknownmlc ! if true then there was an unknown MLC
      logical literalctl

C Strings for surface attributes to pass to insrec.
      character rsuse1*12,rsuse2*12
      character rsparent*12
      integer newsurf   ! index of surface when adding door/glazing
      real angr         ! for rotation

C If browsing then user cannot create a new zone.
      if(browse)then
        call usrmsg('Cannot update model while in browse',
     &    'mode, you must `own` the model!','W')
        return
      endif

C Initialise flag
      ifullysilent=0
      anunknownmlc=.false.

C Set folder separator (fs) to \ or / as required.
      willneedtmc=.false.
      call isunix(unixok)
      if(unixok)then
        fs = char(47)
      else
        fs = char(92)
      endif

C Clear local coordinate array.
      DO 345 IW=1,MS
        XX(iw)=0.0
        YY(iw)=0.0
 345  continue

C Creation of a new zone, begin with length of zone name.
      lnzname(ICOMP)=lnblnk(zname(ICOMP))

C << Unix vs Windows needs to be updated >>
      if(zonepth(1:2).eq.'  '.or.zonepth(1:2).eq.'./')then
        WRITE(DFILE,'(2a)') zname(ICOMP)(1:lnzname(ICOMP)),'.geo'
        WRITE(CFILE,'(2a)') zname(ICOMP)(1:lnzname(ICOMP)),'.con'
        WRITE(TFILE,'(2a)') zname(ICOMP)(1:lnzname(ICOMP)),'.tmc'
        WRITE(OFILE,'(2a)') zname(ICOMP)(1:lnzname(ICOMP)),'.opr'
      else
        WRITE(DFILE,'(4a)') zonepth(1:lnblnk(zonepth)),fs,
     &    zname(ICOMP)(1:lnzname(ICOMP)),'.geo'
        WRITE(CFILE,'(4a)') zonepth(1:lnblnk(zonepth)),fs,
     &    zname(ICOMP)(1:lnzname(ICOMP)),'.con'
        WRITE(TFILE,'(4a)') zonepth(1:lnblnk(zonepth)),fs,
     &    zname(ICOMP)(1:lnzname(ICOMP)),'.tmc'
        WRITE(OFILE,'(4a)') zonepth(1:lnblnk(zonepth)),fs,
     &    zname(ICOMP)(1:lnzname(ICOMP)),'.opr'
      endif
      LGEOM(ICOMP)=DFILE

      write(zd,'(4a)')
     &  zname(ICOMP)(1:lnzname(ICOMP)),' auto-generated for shape ',
     &  shape(ICOMP)(1:lnblnk(shape(ICOMP))),' & attributes.'
      zdesc(ICOMP)=zd
      lnzdesc(ICOMP)=lnblnk(zd)  ! update the length of this string.

C Choices based on the value of shape.
      IER=0
      if(shape(ICOMP)(1:3).eq.'box')then
        ishape=1
      elseif(shape(ICOMP)(1:7).eq.'extrude')then
        ishape=2
      elseif(shape(ICOMP)(1:4).eq.'poly')then
        ishape=3
      else
        ier=1
        return
      endif
      IF(ishape.eq.1)THEN
        XOCB=zorigin(icomp,1)
        YOCB=zorigin(icomp,2)
        ZOCB=zorigin(icomp,3)
        DXCB=zsize(icomp,1)
        DYCB=zsize(icomp,2)
        DZCB=zsize(icomp,3)

C Pending rotation angle and set point to origin of the zone.
        ANGR=rotateit(ICOMP,1)
        rotateit(ICOMP,2)=XOCB
        rotateit(ICOMP,3)=YOCB

C Convert box into a gen description. Pass in critical dimensions
C and expect the data back via common G1.
        CALL ERECC(XOCB,YOCB,ZOCB,DXCB,DYCB,DZCB,ANGR)
        rotateit(ICOMP,1)=0.0
        rotateit(ICOMP,2)=0.0
        rotateit(ICOMP,3)=0.0
        NSUR=6
        NZSUR(ICOMP)=NSUR
        IUZBASEA(icomp)=0
        IBASES(ICOMP,1)=6
        IZBASELIST(ICOMP)=1
      ELSEIF(ishape.eq.2)THEN

C User begins with an extruded shape.
        Z1=zorigin(icomp,1)
        Z2=zorigin(icomp,2)
        NW=nbwalls(ICOMP)

        DO 45 IW=1,NW
          XX(IW)=szcoords(icomp,IW,1)
          YY(IW)=szcoords(icomp,IW,2)
   45   CONTINUE

C Assume that the origin of the zone is 1st coord.
        AR=rotateit(ICOMP,1)
        rotateit(ICOMP,2)=XX(1)
        rotateit(ICOMP,3)=YY(1)

C Convert into a gen body and rotate if required. Pass in 
C critical dimensions to eregc and expect the data back
C via common G1. ESCROT also works via G1.
        CALL EREGC(NW,Z1,Z2,XX,YY)
        IF(AR.LT.-.01.OR.AR.GT..01)then
          x1=XX(1)
          y1=YY(1)
          CALL ESCROT(AR,x1,y1)
          rotateit(ICOMP,1)=0.0
          rotateit(ICOMP,2)=0.0
          rotateit(ICOMP,3)=0.0
        endif
        NSUR=NW+2
        NZSUR(ICOMP)=NSUR
        IUZBASEA(icomp)=0
        IBASES(ICOMP,1)=NSUR
        IZBASELIST(ICOMP)=1
      ELSEIF(ishape.eq.3)THEN

C For a poly shaped zone use similar logic to reading geometry file.
C Instantiate nsur and ntv for use in escrot and filsur.
        NSUR=nbwalls(ICOMP)
        NZSUR(ICOMP)=nbwalls(ICOMP)
        NTV=nztv(icomp)
        AR=rotateit(ICOMP,1)
        DO 62 I=1,nztv(icomp)
          X(I)=szcoords(ICOMP,I,1)
          Y(I)=szcoords(ICOMP,I,2)
          Z(I)=szcoords(ICOMP,I,3)
   62   CONTINUE

C Remember first zone point as rotation point.
        rotateit(ICOMP,2)=X(1)
        rotateit(ICOMP,3)=Y(1)
        DO 10 I=1,NZSUR(ICOMP)
          NVER(I)=isznver(ICOMP,I)
          DO 12 KV=1,NVER(I)
            JVN(I,KV)=iszjvn(ICOMP,I,KV)
   12     CONTINUE
   10   CONTINUE
        IUZBASEA(icomp)=0
        IBASES(ICOMP,1)=NZSUR(ICOMP)  ! initial guess, need to check
        IZBASELIST(ICOMP)=1

C Rotate if pending. Pass in critical dimensions to escrot
C and expect the data back via common G1.
        if(AR.LT.-.01.OR.AR.GT..01)then
          x1=X(1)
          y1=Y(1)
          CALL ESCROT(AR,x1,y1)
        endif

      ENDIF

C Begin with default assumptions for each surface then overwrite. The
C subroutine filsur updates variables in common block G5.
C Because filsur clears the surface name rememeber in hasname first.
      do 322 ICC=1,NZSUR(ICOMP)
        write(hasname(icc),'(a)') sname(icomp,icc)  ! remember surface name
 322  continue
      CALL FILSUR(ICOMP,0)

C Update the connection list.
      ICCC=NCON
      DO 32 ICC=1,NZSUR(ICOMP)
        ICCC=ICCC+1
        IC1(ICCC)=ICOMP
        IE1(ICCC)=ICC
        ICT(ICCC)=zboundarytype(icomp,icc,1)  ! use boundarytype
        IC2(ICCC)=zboundarytype(icomp,icc,2)
        IE2(ICCC)=zboundarytype(icomp,icc,3)

C Instantiate the surface constructions.
        write(sname(icomp,icc),'(a)') hasname(icc)
        write(smlcn(icomp,icc),'(a)') hasconstr(icc)
        write(sparent(icomp,icc),'(a)') hasparent(icc)
        write(suse(icomp,icc,1),'(a)') hasuse(icc,1)
        write(suse(icomp,icc,2),'(a)') hasuse(icc,2)

C Loop thru all of the known constructions to find the matching construction.
C Also set smlcindex for this surface.
        imatch=0
        smlcindex(icomp,icc)=0  ! assume no matching MLC          
        do 5 ii=1,nmlc
          lnssmlc=lnblnk(SMLCN(icomp,icc))
          if(SMLCN(icomp,icc)(1:lnssmlc).eq.
     &       mlcname(ii)(1:lnmlcname(ii)))then
            imatch=ii
            smlcindex(icomp,icc)=ii   ! remember MLC index     
          endif
  5     continue

        if(metaver.eq.2.or.metaver.eq.3)then

C Take the optics from the META file. If there are unknown constructions
C mark anunknownmlc so construction files are not created.
          write(SOTF(icomp,icc),'(a)') hasoptic(icc)
          if(hasoptic(icc)(1:4).EQ.'OPAQ')then
            continue
          else
            willneedtmc=.true.
          endif
          if(imatch.eq.0) then
            anunknownmlc=.true.   ! mark that one MLC is unknown
          endif
        elseif(metaver.eq.1)then

C See if it is opaque or transparent and copy out the relevant
C section of the construction attributes.
          if(imatch.eq.0) then

C If we do not have a match it might be because we are referencing a local
C database (which does not yet exist because it has not been copied from
C the original model. In this case << ? >>
            call edisp(iuout,'Warning: no matching MLC defined!')
            SOTF(icomp,icc)='UNKNOWN'
            anunknownmlc=.true.   ! mark that one MLC is unknown
          else
            if(mlctype(imatch)(1:4).EQ.'OPAQ')then
              SOTF(icomp,icc)='OPAQUE'
            else
              WRITE(SOTF(icomp,icc),'(A)') mlcoptical(imatch)(1:12)
              willneedtmc=.true.
            endif
          endif
        endif

C Use zboundarytype to instantiate sother.
        IF(zboundarytype(icomp,icc,1).EQ.0)THEN
          SOTHER(icomp,icc,1)='EXTERIOR'
        ELSEIF(zboundarytype(icomp,icc,1).EQ.1)THEN
          SOTHER(icomp,icc,1)='SIMILAR'
        ELSEIF(zboundarytype(icomp,icc,1).EQ.2)THEN
          SOTHER(icomp,icc,1)='CONSTANT'
        ELSEIF(zboundarytype(icomp,icc,1).EQ.3)THEN
          SOTHER(icomp,icc,1)='ANOTHER'
        ELSEIF(zboundarytype(icomp,icc,1).EQ.4)THEN
          SOTHER(icomp,icc,1)='GROUND'
        ELSEIF(zboundarytype(icomp,icc,1).EQ.5)THEN
          SOTHER(icomp,icc,1)='ADIABATIC'
        ELSEIF(zboundarytype(icomp,icc,1).EQ.6)THEN
          SOTHER(icomp,icc,1)='BASESIMP'
        ELSEIF(zboundarytype(icomp,icc,1).EQ.7)THEN
          SOTHER(icomp,icc,1)='IDENT_CEN'
        ELSE

C Reset to represent UNKNOWN
          SOTHER(icomp,icc,1)='UNKNOWN'
        ENDIF
        write(SOTHER(icomp,icc,2),'(i2.2)') zboundarytype(icomp,icc,2)
        write(SOTHER(icomp,icc,3),'(i3.3)') zboundarytype(icomp,icc,3)

C Set connection based variables.
        IZSTOCN(icomp,icc)=iccc
   32 CONTINUE
      NCON=ICCC

      CTYPE(icomp)='GEN '
      NDP(ICOMP)=3
      IDPN(ICOMP,1)=0
      IDPN(ICOMP,2)=0
      IDPN(ICOMP,3)=0
      NZSUR(ICOMP)=NSUR  ! update nzsur() it is needed by zgupdate.
      NZTV(ICOMP)=NTV

C Update the G7 common blocks and then if assign ZBASEA.
      call zgupdate(0,ICOMP,ier)
      if(ishape.eq.1)then
        ZBASEA(icomp)= SNA(icomp,IBASES(ICOMP,1))
      elseif(ishape.eq.2)then
        ZBASEA(icomp)= SNA(icomp,IBASES(ICOMP,1))
      elseif(ishape.eq.3)then

C For poly shape loop to find the floor in order to set up the floor area
        lastlist=0
        do 222 iSurface=1,NSUR
          if(SVFC(icomp,iSurface).eq.'FLOR')then
            if(SUSE(icomp,iSurface,1)(1:5).eq.'FURNI')then
              continue  ! ignore horizontal surfaces marked as furniture.
            elseif(SUSE(Icomp,iSurface,1)(1:6).eq.'REVEAL')then
              continue  ! ignore horizontal surfaces marked as reveal.
            else

C If surface `flor` not included in the list add it.
              lastlist=lastlist+1
              if(lastlist.le.6)then
                IBASES(ICOMP,lastlist)=iSurface
                ZBASEA(ICOMP)=ZBASEA(ICOMP)+SNA(ICOMP,iSurface)
              endif	
            endif
          endif
 222    continue
  43  continue
      write(outs,'(a,f6.2,2a)') 'Base area estimated at ',
     &    ZBASEA(ICOMP),'m^2 for ',zname(icomp)
      call edisp(iuout,outs)
      if(lastlist.gt.0)then
        iuzbasea(icomp)=0           ! signal area from orientation scan
         izbaselist(icomp)=lastlist  ! remember how many items in list.
      endif
      endif

C Save this to file before passing into the geometry editing facility.
C Use gversion 1.1.
      NCOMP=NCOMP+1  ! Temporarily update NCOMP for geowrite use.
      if(igupgrade.lt.2)then
        igupgrade=2
        gversion(icomp) =1.1
      endif
      call geowrite(IFIL+2,LGEOM(ICOMP),ICOMP,iuout,3,IER)
      IF(IER.EQ.1)THEN
        call usrmsg('Problem creating geometry file.','giving up','W')
        return
      ENDIF

C Read in control file if defined. Note: the control file will not yet
C know about this new zone so decrement NCOMP before read and re-establish
C after reading control file. If meta file included literal blocks of
C control file then no need to adjust.
      OKC=.false.
      if(LCTLF(1:1).ne.' '.and.LCTLF(1:4).ne.'UNKN')OKC=.true.
      if(literalctl) OKC=.false.
      if(OKC)then
        NCOMP=NCOMP-1  ! Temporarily decrement NCOMP.
        ICTLF=IFIL+1
        CALL ERPFREE(ICTLF,ISTAT)
        call FINDFIL(LCTLF,XST)
        IF(XST) CALL EZCTLR(ICTLF,0,IUOUT,IER)
        NCOMP=NCOMP+1  ! Re-establish NCOMP.
      endif

      NCCODE(ICOMP)=NCOMP
      if(OKC)then
        icascf(NCOMP)=0
        call usrmsg(' updating control for additional zone...',' ','-')
        call CTLWRT(ICTLF,IER)
        call usrmsg(' ',' ','-')
      endif

C Update the G7 common blocks and then the configuration file.
C Use the 's' parameter to emkcfg to avoid asking for cnn file name.
      call zgupdate(0,ICOMP,ier)
      CALL EMKCFG('s',IER)

C Next, if there are mass pairs to add make standard calls to insert them.
      if(znbmass(icomp).gt.0)then
        do 79 imu=1,znbmass(icomp)
          call eclose(zdatamass(icomp,imu,1),90.0,0.01,closev)
          write(SN,'(a)')
     &      ztextmass(icomp,imu,1)(1:lnblnk(ztextmass(icomp,imu,1)))
          write(constr,'(a)')
     &      ztextmass(icomp,imu,2)(1:lnblnk(ztextmass(icomp,imu,2)))
          write(opt,'(a)')
     &      ztextmass(icomp,imu,3)(1:lnblnk(ztextmass(icomp,imu,3)))
          if(closev)then
            call addmass(ICOMP,'VM','s',zdatamass(icomp,imu,2),
     &        zdatamass(icomp,imu,3),zdatamass(icomp,imu,4),
     &        zdatamass(icomp,imu,7),zdatamass(icomp,imu,5),
     &        zdatamass(icomp,imu,6),SN,constr,OPT,INVT)
          else
            call addmass(ICOMP,'HM','s',zdatamass(icomp,imu,2),
     &        zdatamass(icomp,imu,3),zdatamass(icomp,imu,4),
     &        zdatamass(icomp,imu,7),zdatamass(icomp,imu,5),
     &        zdatamass(icomp,imu,6),SN,constr,OPT,INVT)
          endif

          call geowrite(IFIL+2,LGEOM(ICOMP),ICOMP,iuout,3,IER)
  79    continue
      endif

C Predefined objects - task to decide how these might be represented
C in meta file (if at all).
C << work to be done >>

      CALL EMKCFG('s',IER)

C If no errors then proceed to create constructions files. The
C file name is already known, assume that attribution is complete
C so write the header of construction file so something exists and then
C call edcon to complete the process. If the file already exists
C it should be overwritten.
      if(anunknownmlc)then   ! if a MLC is unknown do not bother with edcon
        continue
      else
        IUF=IFIL+2
        LTHRM(ICOMP)=CFILE
        CALL EFOPSEQ(IUF,LTHRM(ICOMP),3,IER)
        WRITE(IUF,31,IOSTAT=IOS,ERR=146)
     &    zname(ICOMP)(1:lnzname(ICOMP)),CFILE(1:lnblnk(CFILE))
  31    FORMAT('# thermophysical properties of ',a,' defined in ',a,/,
     &  '# no of |air |surface(from geo)| multilayer construction',/,
     &  '# layers|gaps|  no.  name      | database name ')
        CALL ERPFREE(IUF,ISTAT)
      endif

      if(willneedtmc)then
        ITW(icomp)=1
        LTWIN(ICOMP)=TFILE
      endif
      if(anunknownmlc)then   ! if a MLC is unknown do not bother with edcon
        call usrmsg('Problem autogenerating zone construction',
     &    'files. A surface construction was unknown.','W')
      else
        QUIET=.true.
        CALL EDCON(0,iuout,ICOMP,QUIET,IER)  ! should auto create files.
        if(ier.ne.0)then
          call usrmsg('Problem autogenerating zone construction',
     &      'files. Please check','W')
        endif
      endif

C Operational details for this zone need to be handled. One
C option is to rely on a file within the pattern folder for
C the initial description (usefile) and use(3) parameters.
      if(use(1)(1:7).eq.'pattern')then
        LPROJ(ICOMP)=OFILE
        call isunix(unixok)
        if(unixok)then
          fs = char(47)
        else
          fs = char(92)
        endif
        write(subpath,'(6a)',iostat=ios,err=1)
     &    dinstpath(1:lnblnk(dinstpath)),fs,'training',fs,'pattern',fs
        write(sfile,'(2a)')subpath(1:lnblnk(subpath)),
     &    usefile(1:lnblnk(usefile))

C Debug.
C        write(outs,*) 'pattern file is ',sfile(1:lnblnk(sfile))
C        call edisp(iuout,outs)

C Clear commons.
        DO 401 IDTY=1,NBDAYTYPE
          NAC(IDTY)=0
          NCAS(IDTY)=0
 401    CONTINUE
        ctlstr(icomp,1)='no control of air flow  '
        ctlstr(icomp,2)='no control of air flow  '
        ctlstr(icomp,3)='no control of air flow  '
        oprdesc(icomp)='no operations notes (yet)'
        lodlabel(icomp,1)='Occupt'
        lodlabel(icomp,2)='Lights'
        lodlabel(icomp,3)='Equipt'
        lodlabel(icomp,4)='Other '
        lodlabel(icomp,5)='Ann.El'
        lodlabel(icomp,6)='N/A   '
        lodlabel(icomp,7)='N/A   '

C Fill common blocks for new zone based on information in the
C pattern file.

C Set to current zone operation file format.
        ip3ver(icomp)=21
        afact = use(2)
        iflow='ach-ach '
        vflow='ach-ach '
        isource=0
        ifullysilent=1
        idest=icomp
        AUTOVR=.true.
        call COPYAIRFL(isource,idest,afact,iflow,vflow,sfile,
     &    itrc,'s',ier)

        isource=0
        idest=icomp
        icact = use(3)
        iwatt='W-W     '
        itrc=1
        call IMPORTZCASG(isource,idest,icact,iwatt,sfile,itrc,'s',ier)

C Just in case sort the data.
        call PROCESSOLDCAS(ICOMP,0,iuout)
        IUO=IFIL+2
        CALL EMKOPER(IUO,LPROJ(ICOMP),ICOMP,IER)
        AUTOVR=.false.
      endif

C Update configuration file to know about the construction files.
      CALL EMKCFG('s',IER)

      RETURN

  146 if(IOS.eq.2)then
        CALL USRMSG(' No permission to write constructions',' ','W')
      else
        CALL USRMSG(' File write error in constructions',' ','W')
      endif
      return

C Errors when creating string buffers.
   1  write(outs,*) 'error creating subpath: ',subpath
      call edisp(iuout,outs)
      ier=2
      return

      END

C ******************** SILENTREAD
C SILENTREAD reads an esp-r meta file zones data. Data input is
C based on 3rd party generators. If act is '?' then return a
C menu and descriptive text block for the file. Otherwise scan
C and act on the information within the file.

       SUBROUTINE SILENTREAD(IUNIT,LFILE,act,silentreturndirec,IER)
       integer MSZ   ! number of META zones array sizes, edit to
                     ! match MCOM in building.h
       PARAMETER (MSZ=82)

#include "building.h"
#include "model.h"
#include "site.h"
#include "geometry.h"
#include "esprdbfile.h"
#include "espriou.h"
#include "control.h"
#include "seasons.h"
#include "ipvdata.h"

      integer lnblnk  ! function definition

      COMMON/OUTIN/IUOUT,IUIN,IEOUT

      integer ncomp,ncon
      COMMON/C1/NCOMP,NCON
      COMMON/FILEP/IFIL

C Data from site.h:
C groundrefl - External ground reflectivity.
C groundreflmonth  - (12) Monthly values of ground reflectivity
C groundreflmodel  - Type of ground reflectivity model
C           1= constant albedo
C           2= simple model (monthly albedo, number of days with snow on ground
C           3= advanced model (monthly albedo, snow depth read from file)
C snowgroundrefl - Snow reflectivity
C dayswithsnow - (12) Monthly values of number of days with snow on ground
C SNFNAM - *72 Name of the file containing hourly snow depth information
      COMMON/PREC8/SLAT,SLON

      COMMON/SET1/IYEAR,IBDOY,IEDOY,IFDAY,IFTIME
      INTEGER :: IYEAR,IBDOY,IEDOY,IFDAY,IFTIME

C Obstruction blocks via geometry.h.

C Simulation parameter sets.
      common/spfldat/nsset,isset,isstup,isbnstep,ispnstep,issave,isavgh
      INTEGER :: nsset,isset,isstup,isbnstep,ispnstep,issave,isavgh

C IPV description via ipvdata.h.
      common/cctlnm/ctldoc,lctlf
      character LCTLF*72,CTLDOC*248

      integer icascf
      common/cctl/icascf(mcom)

C XMLMETRICS are the metrics that could be included in the input.xml
C in case the flag: <dump_all_data> is false
      common/XMLDATA/XMLMETRICS(200),inumXMLmetrics
      character XMLMETRICS*248

C The model root and mpath passed to silentmodel call.
      character metaroot*32,metampath*72
      common/cadmeta/metaroot,metampath

C iaplic(1) toggle for shading; iaplic(2) toggle for insolation.
C   toggle for assessment where:
C   1 is 'all_applicable', 0 is manual selection of surfaces.
C nsurfcalc nb of shaded surfaces, lstsfcalc() list of applicable surfaces.
C nsurfinso nb of insolation sources, isurfinso() list of insolation sources.
      common/ishdirec/iaplic(MCOM,2),nsurfcalc(MCOM),lstsfcalc(MCOM,MS),
     &       nsurfinso(MCOM),isurfinso(MCOM,MS)

      COMMON/FOPENED/CFGOK,MLDBOK,MATDBOK,CTLOK,OPTKOK,CFCDBOK
      LOGICAL CFGOK,MLDBOK,MATDBOK,CTLOK,OPTKOK,CFCDBOK

      common/CWEC_SOL/iSlr_half_hr_flg
      integer iSlr_half_hr_flg      ! 0 hour-centered; 1 half-hour centered

C Images.
      character imgfmt*4  ! GIF XBMP TIF JPG
      character imgfoc*4  ! FZON FNET FCTL FDFS FPLN
      character limgfil*72  ! file name (extend to 144 char)
      character imgdoc*248  ! text associated with image
      common/imagf/imgfmt(MIMG),imgfoc(MIMG),limgfil(MIMG),imgdoc(MIMG)

      integer noimg  ! number of images
      integer iton   ! 0 if images not yet shown, 1 if yes
      common/imagfi/noimg,iton
      
C Radiance processor cores.
      integer radcores
      common/radcor/radcores

C  zfldK - conductivity; zfldD - density; zfldC - specific heat capacity;
C  zflsA - total shortwave absorptivity; zSWA  - shortwave absorption (W).
      COMMON/zfluid/znotair(mcom),zfldK,zfldD,zfldC,zfldA,
     &              zSWAp(mcom),zSWAf(mcom)
      real zfldK,zfldD,zfldC,zfldA,zSWAp,zSWAf
      LOGICAL znotair

C Initial view.
      common/initv/initvt,EYEMI(3),VIEWMI(3),ANGI
      REAL :: angi,EYEMI,VIEWMI
      INTEGER :: initvt

      integer inumXMLmetrics !counts for the number of xml metrics defined

C << update outstr etc to deal with complex surface edge jvn lists >>
      CHARACTER OUTSTR*248,WORD*32,WORD2*20,outs*124,outs2*124
      CHARACTER tmpvfc*4,tother*15,tother1*15,tother2*15
      character dstmp*24
      character word248*248

      character act*1  ! action requested.
      character simact*6 ! action for creating assessments.
      character*84 odir    ! where are folders

C At the head of the file is a title and a description
C modeltitle (char 72) and modeldocblock (char 248) from model.h which can be
C used by scanning code to identify the contents of a silent file.
      character weather*32     ! file name for the location (no path)
      character*(*) lfile      ! name of the file
      character fs*1           ! file separator
      character doit*124

C For passing to silentmodel.
      character actions*8,root*32,mpath*72

C Same as local variable actions which is passed back to calling code.
      character silentreturndirec*8
      character text*72,literalfile*72,pth*148,literalpath*72

      logical unixok
      real valx,valy,valz      ! local variable for extrude co-ords.

C Paremeters passed to:silentzone
      character use*8,szuse*8,usefile*32  ! usage pattern directives
      dimension use(3),szuse(MSZ,3),usefile(MSZ)

      character hasconstr*32,szhasconstr*32
      dimension hasconstr(MS),szhasconstr(MSZ,MS)  ! construction to use for each surface
      character hasoptic*32,szhasoptic*32
      dimension hasoptic(MS),szhasoptic(MSZ,MS)    ! optics to use for each surface
      character hasparent*12,szhasparent*12
      dimension hasparent(MS),szhasparent(MSZ,MS)  ! construction to use for each surface
      character hasuse*12,szhasuse*12
      dimension hasuse(MS,2),szhasuse(MSZ,MS,2)  ! construction to use for each surface

C Site related local variables.
      logical havesite       ! true if tokens included
      character hourlysnowfile*72  ! same as SNFNAM
      logical havehourlysnowfile   ! true if set
      logical literalctl

      real hoursGTM   ! hours before or after GTM
      real groundtemps(12,MGRDP)
      real GVA(12)  ! temporary array to hold a dozen values.
      integer groundtempsets
      logical havegroundtempsets
      logical havesimparameters
      LOGICAL IGDCVS,IGDCNC,IGDNDC,IGDTAQ
      integer simstartup   ! use for isstup
      integer simzonetimestep ! use for isbnstep
      integer simplanttimestep ! use for ispnstep
      integer simsavelevel  ! use for issave
      integer i_ctl_link  !flag to indicate the zone has a basic control
      real ht_Setpoint,cl_Setpoint !heating and cooling set-points
      dimension i_ctl_link(MCOM),ht_Setpoint(MCOM),cl_Setpoint(MCOM)
      integer isilentncf !silent decide how many control functions needed
      integer lsn  ! length of currentfile
      integer itrc
      integer iva
      dimension iva(MS)  ! for use with writing arrays of numbers

      integer metaver  ! one for 1.1 two for 1.2 three for 1.3

      IER=0; itrc=0

C Determine operating system
      call isunix(unixok)
      if(unixok)then
        fs = char(47)
      else
        fs = char(92)
      endif

C Set icomp to zero if there are no zones in the model.
      if(ncomp.eq.0)then
        icomp=0
      endif

C Get the current folder and display options to the user.
C first list any files with .cfg in name.
      odir=' '
      call usrdir(odir)
      literalfile=' '; literalctl=.false.

C Clear the sz array structures. The logic is to scan all of the
C zone details and then process the zones after *end is detected
C in the file.
      if(act(1:1).ne.'?')then
        actions=' '; root=' '; mpath=' '
        modeltitle=' '; modeldocblock=' '; weather=' '
        hourlysnowfile=' '
        simact='------'   ! assume no simulation directives

        groundreflmodel=0
        skyview=0.0; groundview=0.0; buildingview=0.0
        havesite=.false.; havehourlysnowfile=.false.
        havegroundtempsets=.false.; groundtempsets=0
        havesimparameters=.false.
        simstartup=0; simzonetimestep=0; simplanttimestep=0
        simsavelevel=0; isilentncf=0

C For each of the possible ESP-r zones, clear values.
        do 41 iz=1,MCOM
          nztv(iz) = 0; nbwalls(iz) = 0
          zname(iz) = ' '; shape(iz) = ' '
          rotateit(iz,1) = 0.0; rotateit(iz,2) = 0.0
          rotateit(iz,3) = 0.0
          zorigin(iz,1) = 0.0; zorigin(iz,2) = 0.0
          zorigin(iz,3) = 0.0
          zsize(iz,1) = 0.0; zsize(iz,2) = 0.0; zsize(iz,3) = 0.0
          znbmass(iz) = 0
          do 39 isu=1,MS
            zboundarytype(iz,isu,1)=0; zboundarytype(iz,isu,2)=0
            zboundarytype(iz,isu,3)=0
            sname(iz,isu)=' '
            isznver(iz,isu)=0
            do 54 ivu=1,MV
              iszjvn(iz,isu,ivu)=0
  54        continue
  39      continue
          do 59 imu=1,4
            zdatamass(iz,imu,1)= 0.0; zdatamass(iz,imu,2)= 0.0
            zdatamass(iz,imu,3)= 0.0; zdatamass(iz,imu,4)= 0.0
            zdatamass(iz,imu,5)= 0.0; zdatamass(iz,imu,6)= 0.0
            zdatamass(iz,imu,7)= 0.0
            ztextmass(iz,imu,1)= ' '; ztextmass(iz,imu,2)= ' '
            ztextmass(iz,imu,3)= ' '
  59      continue
  41    continue

C For each of the META specific arrays, clear values.
        do 42 iz=1,MSZ
          szuse(iz,1) = ' '; szuse(iz,2) = ' '; szuse(iz,3) = ' '
          usefile(iz) = ' '
          nbobs(iz) = 0
          ht_Setpoint(iz)= 0.0; cl_Setpoint(iz)= 0.0
          i_ctl_link(iz)=0
          do 55 ibu=1,MB
            XOB(iz,ibu) = 0.0; YOB(iz,ibu) = 0.0; ZOB(iz,ibu) = 0.0
            DXOB(iz,ibu) = 0.0; DYOB(iz,ibu) = 0.0; DZOB(iz,ibu) = 0.0
            BANGOB(iz,ibu,1) = 0.0; BANGOB(iz,ibu,2) = 0.0
            BANGOB(iz,ibu,3) = 0.0
            OPOB(iz,ibu) = 1.0
            BLOCKNAME(iz,ibu) = ' '; LNBLOCKNAME(iz,ibu) = 0
            BLOCKMAT(iz,ibu) = ' '; LNBLOCKMAT(iz,ibu) = 0
            BLOCKTYP(iz,ibu) = 'obs '
            do 56 ibe=1,8
              XBP(iz,ibu,ibe)=0.0; YBP(iz,ibu,ibe)=0.0
              ZBP(iz,ibu,ibe)=0.0
  56        continue
  55      continue
          do 53 isu=1,MTV   ! clear the coordinates
            szcoords(iz,isu,1)=0.0; szcoords(iz,isu,2)=0.0
            szcoords(iz,isu,3)=0.0
  53      continue
          do 52 isu=1,MS
            szhasconstr(iz,isu)=' '  ! clear MLC name and optics name
            szhasoptic(iz,isu)=' '
            szhasparent(iz,isu)=' '
            szhasuse(iz,isu,1)=' '
            szhasuse(iz,isu,2)=' '
  52      continue
  42    continue
      endif

      NS=0     ! temporary array for counting surfaces.
      nsz=0
      weather='none'
      OUTSTR=' '
      metaver =2  ! initial assumption

C Initialise meta file. and set currentfile.
      CALL EFOPSEQ(IUNIT,LFILE,1,IER)
      IF(IER.LT.0)THEN
        write(outs,'(3a)') 'Meta file ',LFILE(1:lnblnk(LFILE)),
     &      ' could not be opened.'
        call edisp(iuout,outs)
        IER=1
        RETURN
      endif
      write(currentfile,'(a)') LFILE(1:lnblnk(LFILE))

C Read header lines from file, the 1.1 version looks like:
C *silent_input
C The version 1.2 (which includes surface optical attributes looks like:
C *meta_file and version 1.3 which has full set of surface attribures
C looks like *meta_file_v3
      CALL LSTRIPC(IUNIT,OUTSTR,99,ND,1,'meta line 1',IER)
      IF(IER.NE.0)goto 1001
      if(OUTSTR(1:13).eq.'*silent_input')then
        metaver =1
      elseif(OUTSTR(1:13).eq.'*meta_file_v3')then
        metaver =3
      elseif(OUTSTR(1:10).eq.'*meta_file')then
        metaver =2
      else

C If we have reached this position then it is an older geometry
C file so close the file and set error message.
        CALL ERPFREE(IUNIT,ios)
        ier=2
        call usrmsg('Not a silent instruction file.',LFILE,'W')
        return
      endif

C Read data line. If end of file found then jump
C to 1000 to see if there are zones to process.
  43  CALL LSTRIPC(IUNIT,OUTSTR,99,ND,1,'silent tags',IER)
      IF(IER.NE.0)goto 1000
      K=0
      CALL EGETW(OUTSTR,K,WORD,'W','*vertex or *surface tags',IER)
      if(WORD(1:5).eq.'*Date'.or.WORD(1:5).eq.'*date')then
        CALL EGETRM(OUTSTR,K,dstmp,'W','date stamp',IER)
        goto 43
      elseif(WORD(1:5).eq.'*menu')then  ! depreciated
        CALL EGETRM(OUTSTR,K,modeltitle,'W','meta title',IER)
        goto 43
      elseif(WORD(1:6).eq.'*title')then
        CALL EGETRM(OUTSTR,K,modeltitle,'W','meta title',IER)
        goto 43
      elseif(WORD(1:4).eq.'*doc')then
        CALL EGETRM(OUTSTR,K,modeldocblock,'W','meta documentation',
     &    IER)
        if(act(1:1).eq.'?')then
          return
        else
          goto 43
        endif
      elseif(WORD(1:7).eq.'*action')then

C An action line is like this:
C *action,new,boxer,single,/tmp/boxer
        CALL EGETW(OUTSTR,K,actions,'W','action string',IER)
        silentreturndirec=actions
        if(actions(1:3).eq.'new'.or.actions(1:6).eq.'within')then
          CALL EGETW(OUTSTR,K,root,'W','root name',IER)
          CALL EGETW(OUTSTR,K,WORD,'W','distributed',IER)
          CALL EGETRM(OUTSTR,K,mpath,'W','model folder',IER)

C Debug.
          call isunix(unixok)
          if(unixok)then
            write(6,*)  's action ',actions,' root ',root
            write(6,*)  's mpath ',mpath

C If a new model create the main model folder, zones doc and ctl
C folders in case there are literal file contents to extract 
C from the META file.
            if(actions(1:3).eq.'new')then
              write(doit,'(2a)') 'mkdir -p ',mpath(1:lnblnk(mpath))
              call runit(doit,'-')
              write(doit,'(4a)') 'mkdir -p ',mpath(1:lnblnk(mpath)),fs,
     &          'zones'
              call runit(doit,'-')
              write(doit,'(4a)') 'mkdir -p ',mpath(1:lnblnk(mpath)),fs,
     &          'cfg'
              call runit(doit,'-')
              write(doit,'(4a)') 'mkdir -p ',mpath(1:lnblnk(mpath)),fs,
     &          'ctl'
              call runit(doit,'-')
              write(doit,'(4a)') 'mkdir -p ',mpath(1:lnblnk(mpath)),fs,
     &          'doc'
              call runit(doit,'-')
              write(doit,'(4a)') 'mkdir -p ',mpath(1:lnblnk(mpath)),fs,
     &          'dbs'
              call runit(doit,'-')
	      write(literalpath,'(3a)') mpath(1:lnblnk(mpath)),fs,'cfg/'
              write(6,*) literalpath(1:lnblnk(literalpath))
            endif
          endif
        endif
        goto 43
      elseif(WORD(1:7).eq.'*zonpth')then
        CALL EGETRM(OUTSTR,K,zonepth,'W','path to zones',IER)
        goto 43
      elseif(WORD(1:7).eq.'*netpth')then
        CALL EGETRM(OUTSTR,K,netpth,'W','path to nets',IER)
        goto 43
      elseif(WORD(1:7).eq.'*ctlpth')then
        CALL EGETRM(OUTSTR,K,ctlpth,'W','path to ctl',IER)
        goto 43
      elseif(WORD(1:7).eq.'*mscpth')then
        CALL EGETRM(OUTSTR,K,mscpth,'W','path to miscel files',IER)
        goto 43
      elseif(WORD(1:7).eq.'*imgpth')then
        CALL EGETRM(OUTSTR,K,imgpth,'W','path to img',IER)
        goto 43
      elseif(WORD(1:7).eq.'*tmppth')then
        CALL EGETRM(OUTSTR,K,tmppth,'W','path to temporary files',IER)
        goto 43
      elseif(WORD(1:7).eq.'*docpth')then
        CALL EGETRM(OUTSTR,K,docpth,'W','path to documents',IER)
        goto 43
      elseif(WORD(1:7).eq.'*radpth')then
        CALL EGETRM(OUTSTR,K,radpth,'W','path to radiance',IER)
        goto 43
      elseif(WORD(1:7).eq.'*dbspth')then
        CALL EGETRM(OUTSTR,K,dbspth,'W','path to loc databases',IER)
        goto 43

C Databases. There may or may not be tags for specific databases. If
C none is supplied assume the default names apply. If local database
C then literal blocks may be included that will need to be copied into
C the dbs folder.
      elseif(WORD(1:4).eq.'*mat'.or.WORD(1:7).eq.'*stdmat')then
C        call scan_database_names(IUNIT,foundmould,foundsbem,IRNCMD,ier)
        iuf=ifil+2
        call scan_embedded_database(IUNIT,iuf,literalpath,ier)

        write(weather,'(a)') lclim(1:32)  ! Save weather file name.
        goto 43
      elseif(OUTSTR(1:8).eq.'* Ground'.or.
     &       OUTSTR(1:8).eq.'* GROUND'.or.
     &       WORD(1:8).eq.'*ground ')then
        call scan_ground(IUNIT,ITRC,IGDCVS,IGDCNC,IGDNDC,IGDTAQ,ier)
        if(ier.eq.0) goto 43
      elseif(WORD(1:9).eq.'*calename'.or.WORD(1:9).eq.'*calentag'.or.
     &       WORD(1:7).eq.'*daytag')then
        call scan_day_types(IUNIT,idty,ier)
        if(ier.eq.0) goto 43
      elseif(WORD(1:4).eq.'*sps' .or. WORD(1:6).eq. '*mysps')then
        call scan_sps(IUNIT,ier)     ! Process simulation paramater sets.
        goto 43
      elseif(WORD(1:5).eq.'*year')then
        CALL EGETWI(OUTSTR,K,IYEAR,1900,2051,'W','year',IER)
        goto 43
      elseif(WORD(1:8).eq.'*seasons')then
        call scan_season(IUNIT,ier)  ! Parse seasons attributes.
        goto 43
      elseif(WORD(1:8).eq.'*cooling')then
        CALL EGETWR(OUTSTR,K,cl_Setpoint(nsz),-102.,102.,'W',
     &       'cool Setp',IER)
        goto 43
      elseif(WORD(1:5).eq.'*cord')then
        nztv(nsz)=nztv(nsz)+1   ! increment number of coordinates in zone
        nbc=nztv(nsz)
        CALL EGETWR(OUTSTR,K,valx,0.,0.,'-','cord-x',IER)
        CALL EGETWR(OUTSTR,K,valy,0.,0.,'-','cord-y',IER)
        szcoords(nsz,nbc,1)=valx
        szcoords(nsz,nbc,2)=valy
        if(shape(nsz)(1:4).eq.'poly')then
          CALL EGETWR(OUTSTR,K,valz,0.,0.,'-','cord-z',IER)
          szcoords(nsz,nbc,3)=valz
        else
          szcoords(nsz,nbc,3)=0.0
        endif
        goto 43
      elseif(WORD(1:5).eq.'*ctl ')then
        CALL EGETRM(OUTSTR,K,text,'W','control file',IER)
        LCTLF=text(1:72)
        literalfile=text(1:72)
        literalctl=.true.
        goto 43
      elseif(WORD(1:5).eq.'*door')then

C znbdoor & zhasdoor are depreciated. Just jump to 43.
        goto 43
      elseif(WORD(1:4).eq.'*opr')then
        CALL EGETRM(OUTSTR,K,text,'W','operation file',IER)
        LPROJ(nsz)=text(1:72)
        literalfile=text(1:72)
        goto 43
      elseif(WORD(1:9).eq.'*end_zone')then
        goto 43
      elseif(WORD(1:5).eq.'*end ')then

C Jump to process the collected data.
        goto 1000
      elseif(WORD(1:6).eq.'*glaze')then

C Depreciated. Continue to 43.
        goto 43
      elseif(WORD(1:19).eq.'*ground_refl_annual')then
        CALL EGETWR(OUTSTR,K,groundrefl,0.,1.,'W',
     &    'Ground reflectance annual value',IER)
        groundreflmodel=1  ! if annual value provided set constant albedo model
        goto 43
      elseif(WORD(1:20).eq.'*ground_refl_monthly')then

C Ground reflectivity model, no-snow monthly albedos, snow albedo
        CALL EGETWI(OUTSTR,K,groundreflmodel,1,3,'F',
     &    'ground refl. model',IER)
        DO 441 I=1,12
          CALL EGETWR(OUTSTR,K,groundreflmonth(I),0.,1.,'W',
     &      'monthly gr. refl.',IER)
  441   CONTINUE
        CALL EGETWR(OUTSTR,K,snowgroundrefl,0.,1.,'W','snow refl.',IER)
        goto 43
      elseif(WORD(1:8).eq.'*heating')then
        CALL EGETWR(OUTSTR,K,ht_Setpoint(nsz),-102.,102.,'W',
     &       'heat setp',IER)
        goto 43
      elseif(WORD(1:18).eq.'*hourly_snow_depth')then

C If ground refl model is 3 then this file has hourly data.
        CALL EGETRM(OUTSTR,K,hourlysnowfile,'W','snow file',IER)
        havehourlysnowfile=.true.
        goto 43
      elseif(WORD(1:14).eq.'*ideal_control')then
        i_ctl_link(nsz)=1
        goto 43
      elseif(WORD(1:4).eq.'*ipv')then
        goto 43
      elseif(WORD(1:5).eq.'*img ')then

C An image file. If this is a version 4 file then there will be
C a subsequent (long) line with documentation for the images.
        noimg=noimg+1
        CALL EGETW(OUTSTR,K,WORD2,'W','img format',IFLAG)
        write(imgfmt(noimg),'(a4)') WORD2(1:4)
        CALL EGETW(OUTSTR,K,WORD2,'W','focus',IFLAG)
        write(imgfoc(noimg),'(a4)') WORD2(1:4)
        CALL EGETRM(OUTSTR,K,limgfil(noimg),'W','image file',IER)
        goto 43
      elseif(WORD(1:6).eq.'*imdoc')then
        if(noimg.gt.0)then
          CALL EGETRM(OUTSTR,K,imgdoc(noimg),'W','image docum',IER)
        endif
        goto 43
      elseif(WORD(1:4).eq.'*ihc')then
        CALL EGETRM(OUTSTR,K,text,'W','htc file',IER)
        ihc(nsz)=1
        LHCCO(nsz)=text(1:72)
        literalfile=text(1:72)
        goto 43
      elseif(WORD(1:4).eq.'*ivf')then
        CALL EGETRM(OUTSTR,K,text,'W','viewfactor file',IER)
        ivf(nsz)=1
        LVIEW(nsz)=text(1:72)
        literalfile=text(1:72)
        goto 43
      elseif(WORD(1:8).eq.'*latlong')then
        CALL EGETWR(OUTSTR,K,sitelat,-90.,90.,'W','Latitude',IER)
        SLAT=sitelat
        CALL EGETWR(OUTSTR,K,sitelongdif,-15.,15.,'W',
     &    'Longitude dif from time meridian',IER)
        SLON=sitelongdif
        goto 43   ! read another tag
      elseif(WORD(1:8).eq.'*sitealt')then
        call egetwr(outstr,k,sitealt,-100.,3000.,'W',
     &    'site altitude',IER)
        if(ier.eq.0) goto 43
      elseif(WORD(1:6).eq.'*site ')then
        CALL EGETWI(outstr,k,siteexposureindex,1,8,'W',
     &    'Site exposure index',IER)
        CALL EGETWR(OUTSTR,K,groundrefl,0.,1.,'W',
     &    'Ground reflect',IER)
        if(siteexposureindex.eq.8)then  ! Read 3 values on next line.
          CALL LSTRIPC(IUC,OUTSTR,3,ND,1,'header tags',IER)
          K=0
          CALL EGETWR(OUTSTR,K,skyview,0.,1.,'W','sky view',IER)
          CALL EGETWR(OUTSTR,K,groundview,0.,1.,'W','sky view',IER)
          CALL EGETWR(OUTSTR,K,buildingview,0.,1.,'W','sky view',IER)

C Range checking.
          IF(ABS(skyview+groundview+buildingview-1.).GT.0.01) 
     &      CALL USRMSG(
     &      ' Total exposure out of range in',OUTSTR,'W')
          goto 43   ! read another tag
        endif
        goto 43   ! read another tag
      elseif(WORD(1:12).eq.'*slr_half_hr')then ! Solar radiation (half-hour data) flag.
        call egetwi(outstr,k,iSlr_half_hr_flg,0,1,'W',
     &    'solar half-hour flag',IER)
        if(ier.eq.0) goto 43
      elseif(WORD(1:5).eq.'*list')then

C The current value of NS should still point to the surface data
C that was just scanned. The *list line holds the number of
C edges and then the index of each coordinate (as in the normal
C geometry file).
        CALL EGETWI(OUTSTR,K,ival1,3,MV,'W','nb of edges',IER)
        isznver(nsz,NS)=ival1

C Now proceed to read vertices on one or more lines.
        DO 12 KV=1,ival1
          CALL EGETWI(OUTSTR,K,ival2,0,MTV,'F','vertex',IERV)
          IF(IERV.NE.0) THEN
            CALL STRIPC(IUNIT,OUTSTR,0,ND,0,'vertex XYZ',IER)
            K=0
            CALL EGETWI(OUTSTR,K,ival2,0,MTV,'F','vertex',IERV)
          ENDIF
          iszjvn(nsz,NS,KV)=ival2
   12   CONTINUE
        goto 43
      elseif(WORD(1:5).eq.'*mass')then

C The *mass keyword is followed by 9 tokens as follows:
        znbmass(nsz)=znbmass(nsz)+1
        nma=znbmass(nsz)
        CALL EGETW(OUTSTR,K,WORD,'W','mass vertical or horiz',IER)
        if(WORD(1:2).eq.'VM'.or.WORD(1:2).eq.'vm')then
          zdatamass(nsz,nma,1)= 90.0  ! signal vertical
        elseif(WORD(1:2).eq.'HM'.or.WORD(1:2).eq.'hm')then
          zdatamass(nsz,nma,1)= 0.0   ! signal horizontal
        endif
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','mass X org',IER)
        zdatamass(nsz,nma,2)=val1
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','mass Y org',IER)
        zdatamass(nsz,nma,3)=val1
        CALL EGETWR(OUTSTR,K,val1, -99., 99.,'W','mass Z org',IER)
        zdatamass(nsz,nma,4)=val1
        CALL EGETWR(OUTSTR,K,val1,   0., 99.,'W','mass X dis',IER)
        zdatamass(nsz,nma,5)=val1
        CALL EGETWR(OUTSTR,K,val1,   0., 99.,'W','mass Z dis',IER)
        zdatamass(nsz,nma,6)=val1
        CALL EGETWR(OUTSTR,K,val1,-359.,359.,'W','mass rotate',IER)
        zdatamass(nsz,nma,7)=val1
        CALL EGETW(OUTSTR,K,WORD,'W','mass name',IFLAG)
        write(ztextmass(nsz,nma,1),'(a)') WORD(1:lnblnk(WORD))
        CALL EGETW(OUTSTR,K,WORD,'W','mass constr name',IFLAG)
        write(ztextmass(nsz,nma,2),'(a)') WORD(1:lnblnk(WORD))
        write(ztextmass(nsz,nma,3),'(a)') 'OPAQUE'   ! set as opaque
        goto 43
      elseif(WORD(1:16).eq.'*monthly_profile')then

C Ground monthly temperature profiles. 2nd item is the number of profiles.
C And this is followed by one line each for each profile.
        CALL EGETWI(OUTSTR,K,groundtempsets,1,MGRDP,'F','grn prfls',IER)
        IF(groundtempsets.eq.0) goto 43
        DO 28 IGRDP=1,groundtempsets
          CALL EGETWRA(IUNIT,GVA,12,-15.,99.,'W','Gr tmp',IER)
          DO 29 J=1,12
            groundtemps(J,IGRDP)=GVA(J)
   29     CONTINUE
  28    CONTINUE
        havegroundtempsets=.true.
        goto 43
      elseif(WORD(1:18).eq.'*monthly_snow_days')then

C Number of days with snow on the ground (if ground refl model 2)
        DO 442 I=1,12
          CALL EGETWI(OUTSTR,K,dayswithsnow(I),0,31,'W',
     &      'days with snow in each month',IER)
  442   CONTINUE
        goto 43
      elseif(WORD(1:8).eq.'*nbwalls')then

C Number of parent surfaces has a meaning for extruded shapes and set to 6 for
C box shapes.
        if(shape(nsz)(1:7).eq.'extrude')then
          CALL EGETWI(OUTSTR,K,nbw,1,MS,'W','nb parent surfs',IER)
          nbwalls(nsz)=nbw
        elseif(shape(nsz)(1:4).eq.'poly')then
          CALL EGETWI(OUTSTR,K,nbw,1,MS,'W','nb parent surfs',IER)
          nbwalls(nsz)=nbw
        else
          nbwalls(nsz)=4
        endif
        goto 43
      elseif(WORD(1:8).eq.'*nbsurfs')then

C Number of parent surfaces has a meaning for extruded shapes and set to 6 for
C box shapes (typically used with poly type zones).
        if(shape(nsz)(1:7).eq.'extrude')then
          CALL EGETWI(OUTSTR,K,nbw,1,MS,'W','nb parent surfs',IER)
          nbwalls(nsz)=nbw
        elseif(shape(nsz)(1:4).eq.'poly')then
          CALL EGETWI(OUTSTR,K,nbw,1,MS,'W','nb parent surfs',IER)
          nbwalls(nsz)=nbw
        else
          nbwalls(nsz)=4
        endif
        goto 43
      elseif(WORD(1:9).eq.'*contents')then
        CALL EGETRM(OUTSTR,K,lmodelqa,'W','model contents',IER)
        literalfile=lmodelqa(1:72)
        goto 43
      elseif(WORD(1:6).eq.'*notes')then

C Project notes file, read and then go on to fabric/plant description.
        CALL EGETRM(OUTSTR,K,lmodellog,'W','notes file',IER)
        literalfile=lmodellog(1:72)
        goto 43
      elseif(WORD(1:10).eq.'*shad_calc')then

C Shading calculation instructions.
        CALL EGETW(OUTSTR,K,WORD,'W','shad directive',IFLAG)
        if(WORD(1:4).eq.'none')then
          iaplic(nsz,1)=0
          nsurfcalc(nsz)=0
        else
          if(WORD(1:14).eq.'all_applicable')then
            iaplic(nsz,1)=1
          elseif(WORD(1:4).eq.'list')then
            iaplic(nsz,1)=0
          endif

C Read number of surfaces and then the list
          CALL EGETWI(OUTSTR,K,iv,0,MS,'-','nsurfcalc',IER)
          nsurfcalc(nsz)=iv
          IRVA=nsurfcalc(nsz)
          CALL EGETWIA(IUNIT,IVA,IRVA,0,MS,'F','lstsfcalc',IER)
          DO 342 ks=1,nsurfcalc(nsz)
            lstsfcalc(nsz,ks)=IVA(ks)
 342      CONTINUE
        endif
        goto 43
      elseif(WORD(1:11).eq.'*insol_calc')then

C Insolation calculation instructions.
        CALL EGETW(OUTSTR,K,WORD,'W','insolation directive',IFLAG)
        if(WORD(1:4).eq.'none')then
          iaplic(nsz,2)=0
          nsurfinso(nsz)=0
        else
          if(WORD(1:14).eq.'all_applicable')then
            iaplic(nsz,2)=1
          elseif(WORD(1:4).eq.'list')then
            iaplic(nsz,2)=0
          endif

C Applicable surfaces are the current exterior facing surfaces.
          CALL EGETWI(OUTSTR,K,iv,0,MS,'-','nsurfinso',IER)
          nsurfinso(nsz)=iv
          IRVA=nsurfinso(nsz)
          CALL EGETWIA(IUNIT,IVA,IRVA,0,MS,'F','isurfinso',IER)
          DO 344 ks=1,nsurfinso(nsz)
            isurfinso(nsz,ks)=IVA(ks)
 344      CONTINUE
        endif
        goto 43
      elseif(WORD(1:8).eq.'*radcore')then
        CALL EGETWI(OUTSTR,K,radcores,1,4,'W','radiance cores',IER)
        goto 43
      elseif(WORD(1:19).eq.'*surf_heat_transfer')then
        call egetwi(outstr,k,ihct,1,8,'W','inside hc method',IER)
        call egetwi(outstr,k,icorexhct,1,12,'W','ext hc method',IER)
        if(ier.eq.0) goto 43
      elseif(WORD(1:11).eq.'*solar_grid')then

C If there will be shading or insulation then the density of the
C grid is found after the token *solar_grid.
        CALL EGETWI(OUTSTR,K,NOX(nsz),4,MOX,'F','opq grid X',IER)
        CALL EGETWI(OUTSTR,K,NOZ(nsz),4,MOZ,'F','opq grid Z',IER)
        goto 43

C Start of a literal block from a model file. If literalfile is
C not blank then use metatoasci to create that file. The folders will
C not yet exist so create.
      elseif(WORD(1:14).eq.'*start_literal')then
        lnl=lnblnk(literalfile)
        if(lnl.gt.1)then
          iuf=ifil+2
          write(pth,'(2a)') literalpath(1:lnblnk(literalpath)),
     &      literalfile(1:lnblnk(literalfile))
          write(6,*) pth(1:lnblnk(pth))
          call metatoasci(IUNIT,iuf,pth,IER)
          literalfile=' '
        endif
        goto 43
      elseif(WORD(1:5).eq.'*obs3')then

C A 3 rotation obstruction block to be associated with the current zone.
C << how do we detect OPOB? >>
        nbobs(nsz)=nbobs(nsz)+1
        nob=nbobs(nsz)
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','X org',IER)
        XOB(nsz,nob)=val1
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','Y org',IER)
        YOB(nsz,nob)=val1
        CALL EGETWR(OUTSTR,K,val1, -99., 99.,'W','Z org',IER)
        ZOB(nsz,nob)=val1
        CALL EGETWR(OUTSTR,K,val1,   0., 99.,'W','X dis',IER)
        DXOB(nsz,nob)=val1
        CALL EGETWR(OUTSTR,K,val1,   0., 99.,'W','Y dis',IER)
        DYOB(nsz,nob)=val1
        CALL EGETWR(OUTSTR,K,val1,   0., 99.,'W','Z dis',IER)
        DZOB(nsz,nob)=val1
        CALL EGETWR(OUTSTR,K,val1,-359.,359.,'W','rot a',IER)
        BANGOB(nsz,nob,1)=val1
        CALL EGETWR(OUTSTR,K,val1,-180.,180.,'W','rot b',IER)
        BANGOB(nsz,nob,2)=val1
        CALL EGETWR(OUTSTR,K,val1,-180.,180.,'W','rot c',IER)
        BANGOB(nsz,nob,3)=val1
        CALL EGETW(OUTSTR,K,WORD,'W','obs blk name',IFLAG)
        write(BLOCKNAME(nsz,nob),'(a)') WORD(1:lnblnk(WORD))
        LNBLOCKNAME(nsz,nob)=lnblnk(BLOCKNAME(nsz,nob))
        CALL EGETW(OUTSTR,K,WORD,'W','obs mat name',IFLAG)
        write(BLOCKMAT(nsz,nob),'(a)') WORD(1:lnblnk(WORD))
        LNBLOCKMAT(nsz,nob)=lnblnk(BLOCKMAT(nsz,nob))
        BLOCKTYP(nsz,nob)='obs3'
        if(iobs(nsz).eq.0) iobs(nsz)=2  ! mark where obstructions held
        goto 43
      elseif(WORD(1:5).eq.'*obsp')then

C A general polygon obstruction to be associated with the current zone.
C The first line includes (current fixed) integer number of vertices
C followed by number of faces and the name and material.
C The 2nd line has the first 4 coordinates and the 3rd line has the
C next 4 coordinates. The surface and edge ordering is as in a
C standard obstruction block when converted into a GB1 common block.
C << how do we detect OPOB? >>
        nbobs(nsz)=nbobs(nsz)+1
        nob=nbobs(nsz)
        CALL EGETWI(OUTSTR,K,ival,8,8,'F','obs nb vertices',IER)
        CALL EGETWI(OUTSTR,K,ival,6,6,'F','pbs nb faces',IER)
        CALL EGETW(OUTSTR,K,WORD,'W','obs blk name',IFLAG)
        write(BLOCKNAME(nsz,nob),'(a)') WORD(1:lnblnk(WORD))
        LNBLOCKNAME(nsz,nob)=lnblnk(BLOCKNAME(nsz,nob))
        CALL EGETW(OUTSTR,K,WORD,'W','obs mat name',IFLAG)
        write(BLOCKMAT(nsz,nob),'(a)') WORD(1:lnblnk(WORD))
        LNBLOCKMAT(nsz,nob)=lnblnk(BLOCKMAT(nsz,nob))
        BLOCKTYP(nsz,nob)='obsp'
        if(iobs(nsz).eq.0) iobs(nsz)=2  ! mark where obstructions held

        CALL LSTRIPC(IUNIT,OUTSTR,99,ND,1,'first 4 coordinates',IER)
        IF(IER.NE.0)goto 1000
        K=0
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','XBP 1',IER)
        XBP(nsz,nob,1)=val1
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','YBP 1',IER)
        YBP(nsz,nob,1)=val1
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','ZBP 1',IER)
        ZBP(nsz,nob,1)=val1
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','XBP 2',IER)
        XBP(nsz,nob,2)=val1
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','YBP 2',IER)
        YBP(nsz,nob,2)=val1
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','ZBP 2',IER)
        ZBP(nsz,nob,2)=val1
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','XBP 3',IER)
        XBP(nsz,nob,3)=val1
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','YBP 3',IER)
        YBP(nsz,nob,3)=val1
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','ZBP 3',IER)
        ZBP(nsz,nob,3)=val1
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','XBP 4',IER)
        XBP(nsz,nob,4)=val1
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','YBP 4',IER)
        YBP(nsz,nob,4)=val1
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','ZBP 4',IER)
        ZBP(nsz,nob,4)=val1

        CALL LSTRIPC(IUNIT,OUTSTR,99,ND,1,'2nd 4 coordinates',IER)
        IF(IER.NE.0)goto 1000
        K=0
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','XBP 5',IER)
        XBP(nsz,nob,5)=val1
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','YBP 5',IER)
        YBP(nsz,nob,5)=val1
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','ZBP 5',IER)
        ZBP(nsz,nob,5)=val1
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','XBP 6',IER)
        XBP(nsz,nob,6)=val1
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','YBP 6',IER)
        YBP(nsz,nob,6)=val1
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','ZBP 6',IER)
        ZBP(nsz,nob,6)=val1
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','XBP 7',IER)
        XBP(nsz,nob,7)=val1
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','YBP 7',IER)
        YBP(nsz,nob,7)=val1
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','ZBP 7',IER)
        ZBP(nsz,nob,7)=val1
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','XBP 8',IER)
        XBP(nsz,nob,8)=val1
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','YBP 8',IER)
        YBP(nsz,nob,8)=val1
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','ZBP 8',IER)
        ZBP(nsz,nob,8)=val1
        goto 43

      elseif(WORD(1:5).eq.'*obs ')then

C An obstruction block to be associated with the current zone.
C << how do we detect OPOB? >>
        nbobs(nsz)=nbobs(nsz)+1
        nob=nbobs(nsz)
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','X org',IER)
        XOB(nsz,nob)=val1
        CALL EGETWR(OUTSTR,K,val1,-999.,998.,'W','Y org',IER)
        YOB(nsz,nob)=val1
        CALL EGETWR(OUTSTR,K,val1, -99., 99.,'W','Z org',IER)
        ZOB(nsz,nob)=val1
        CALL EGETWR(OUTSTR,K,val1,   0., 99.,'W','X dis',IER)
        DXOB(nsz,nob)=val1
        CALL EGETWR(OUTSTR,K,val1,   0., 99.,'W','Y dis',IER)
        DYOB(nsz,nob)=val1
        CALL EGETWR(OUTSTR,K,val1,   0., 99.,'W','Z dis',IER)
        DZOB(nsz,nob)=val1
        CALL EGETWR(OUTSTR,K,val1,-359.,359.,'W','rot a',IER)
        BANGOB(nsz,nob,1)=val1
        CALL EGETW(OUTSTR,K,WORD,'W','obs blk name',IFLAG)
        write(BLOCKNAME(nsz,nob),'(a)') WORD(1:lnblnk(WORD))
        LNBLOCKNAME(nsz,nob)=lnblnk(BLOCKNAME(nsz,nob))
        CALL EGETW(OUTSTR,K,WORD,'W','obs mat name',IFLAG)
        write(BLOCKMAT(nsz,nob),'(a)') WORD(1:lnblnk(WORD))
        LNBLOCKMAT(nsz,nob)=lnblnk(BLOCKMAT(nsz,nob))
        BLOCKTYP(nsz,nob)='obs '
        if(iobs(nsz).eq.0) iobs(nsz)=2  ! mark where obstructions held
        goto 43
      elseif(WORD(1:7).eq.'*origin')then
        if(shape(nsz)(1:3).eq.'box')then
          CALL EGETWR(OUTSTR,K,val1,0.,0.,'-','box X origin',IER)
          zorigin(nsz,1)=val1
          CALL EGETWR(OUTSTR,K,val1,0.,0.,'-','box Y origin',IER)
          zorigin(nsz,2)=val1
          CALL EGETWR(OUTSTR,K,val1,0.,0.,'-','box Z origin',IER)
          zorigin(nsz,3)=val1
        elseif(shape(nsz)(1:7).eq.'extrude')then
          CALL EGETWR(OUTSTR,K,val1,0.,0.,'-','extrude Z base',IER)
          zorigin(nsz,1)=val1
          CALL EGETWR(OUTSTR,K,val1,0.,0.,'-','extrude Z top',IER)
          zorigin(nsz,2)=val1
          CALL EGETWR(OUTSTR,K,val1,0.,0.,'-','extrude ignore',IER)
          zorigin(nsz,3)=val1
        elseif(shape(nsz)(1:4).eq.'poly')then
          zorigin(nsz,1)=0.0
          zorigin(nsz,2)=0.0
          zorigin(nsz,3)=0.0
        endif
        goto 43
      elseif(WORD(1:7).eq.'*rotate')then

C Scan for pending X and Y coordinates to rotate around if there are 3 items on the line.
        CALL EGETWR(OUTSTR,K,rot,0.,0.,'-','rotation angle',IER)
        rotateit(nsz,1)=rot
        if(ND.gt.2)then
          CALL EGETWR(OUTSTR,K,rot,0.,0.,'-','rotation point X',IER)
          rotateit(nsz,2)=rot
          CALL EGETWR(OUTSTR,K,rot,0.,0.,'-','rotation point Y',IER)
          rotateit(nsz,3)=rot
        else
          rotateit(nsz,2)=0.0  ! assume site origin to rotate around
          rotateit(nsz,3)=0.0
        endif
        goto 43
      elseif(WORD(1:18).eq.'*previous_rotation')then

C Deprecated feature.
        goto 43
      elseif(WORD(1:15).eq.'*sim_parameters')then

C Get the assessment independant simulation parameters.
        CALL EGETWI(OUTSTR,K,simstartup,0,300,'W','startup',IER)
        isstup=simstartup
        CALL EGETWI(OUTSTR,K,simzonetimestep,1,60,'W','zone_ts',IER)
        isbnstep=simzonetimestep
        CALL EGETWI(OUTSTR,K,simplanttimestep,0,100,'W','plant_ts',IER)
        ispnstep=simplanttimestep
        CALL EGETWI(OUTSTR,K,simsavelevel,0,6,'W','save_lv',IER)
        issave=simsavelevel
        havesimparameters=.true.
        goto 43
      elseif(WORD(1:9).eq.'*site_loc')then

C The *site_loc line includes latitute, longitude difference from time
C meridian, time zone (hours before or after GMT).
        CALL EGETWR(OUTSTR,K,sitelat,-90.,90.,'W','Latitude',IER)
        CALL EGETWR(OUTSTR,K,sitelongdif,-15.,15.,'W',
     &    'Longitude dif from time meridian',IER)
        CALL EGETWR(OUTSTR,K,hoursGTM,-12.,12.,'W','Hours +-GMT',IER)
        havesite=.true.
        goto 43
      elseif(WORD(1:9).eq.'*site_exp')then

C The *site_exp line includes site exposure and if
C siteexposureindex is 8 then read three more numbers
C for the general exposure to sky, ground and other buildings.
        CALL EGETWI(OUTSTR,K,siteexposureindex,1,8,'W',
     &    'Site exposure',IER)
        if(siteexposureindex.eq.8)then
          CALL EGETWR(OUTSTR,K,skyview,0.,1.,'W','Sky viewfactor',IER)
          CALL EGETWR(OUTSTR,K,groundview,0.,1.,'W',
     &      'Ground viewfactor',IER)
          CALL EGETWR(OUTSTR,K,buildingview,0.,1.,'W',
     &      'Building viewfactor',IER)
        endif
        havesite=.true.
        goto 43
      elseif(WORD(1:5).eq.'*site')then

C The *site line includes latitute, longitude difference from time
C meridian, time zone (hours before or after GMT), site exposure
C and if siteexposureindex is 8 then read three more numbers
C for the general exposure to sky, ground and other buildings.
C This format is included for older META files.
        CALL EGETWR(OUTSTR,K,sitelat,-90.,90.,'W','Latitude',IER)
        CALL EGETWR(OUTSTR,K,sitelongdif,-15.,15.,'W',
     &    'Longitude dif from time meridian',IER)
        CALL EGETWR(OUTSTR,K,hoursGTM,-12.,12.,'W','Hours +-GMT',IER)
        CALL EGETWI(OUTSTR,K,siteexposureindex,1,8,'W',
     &    'Site exposure',IER)
        if(siteexposureindex.eq.8)then
          CALL EGETWR(OUTSTR,K,skyview,0.,1.,'W','Sky viewfactor',IER)
          CALL EGETWR(OUTSTR,K,groundview,0.,1.,'W',
     &      'Ground viewfactor',IER)
          CALL EGETWR(OUTSTR,K,buildingview,0.,1.,'W',
     &      'Building viewfactor',IER)
        endif
        havesite=.true.
        goto 43
      elseif(WORD(1:5).eq.'*size')then
        if(shape(nsz)(1:3).eq.'box')then
          CALL EGETWR(OUTSTR,K,val1,0.,0.,'-','box length',IER)
          zsize(nsz,1)=val1
          CALL EGETWR(OUTSTR,K,val1,0.,0.,'-','box width',IER)
          zsize(nsz,2)=val1
          CALL EGETWR(OUTSTR,K,val1,0.,0.,'-','box height',IER)
          zsize(nsz,3)=val1
        elseif(shape(nsz)(1:7).eq.'extrude')then
          zsize(nsz,1)=0.0   ! this is ignored for extruded shape
          zsize(nsz,2)=0.0
          zsize(nsz,3)=0.0
        elseif(shape(nsz)(1:4).eq.'poly')then
          zsize(nsz,1)=0.0   ! this is ignored for extruded shape
          zsize(nsz,2)=0.0
          zsize(nsz,3)=0.0
        endif
        goto 43
      elseif(WORD(1:10).eq.'*start_xml')then
        inumXMLmetrics=0
        CALL LSTRIPC(IUNIT,OUTSTR,99,ND,1,'read metrics',IER)
        K=0
        CALL EGETW(OUTSTR,K,WORD248,'W','xml data',IER)
        do while(WORD248(1:8).ne.'*end_xml')
          inumXMLmetrics=1+inumXMLmetrics
          XMLMETRICS(inumXMLmetrics)=WORD248
          CALL LSTRIPC(IUNIT,OUTSTR,99,ND,1,'read metrics',IER)
          K=0
          CALL EGETW(OUTSTR,K,WORD248,'W','xml data',IER)
        end do
        if(WORD(1:8).eq.'*end_xml'.or.WORD248(1:8).eq.'*end_xml')then
          goto 43
        endif
        goto 43
      elseif(WORD(1:11).eq.'*start_zone')then

C If this is a new model increment icomp. If the model already has
C zones then set icomp to the next available zone.
        if(ncomp.eq.0)then
          icomp=icomp+1      ! new model, ncomp not yet set
        else
          if(icomp.eq.0)then
            icomp=ncomp+1    ! increment existing zone index
          elseif(icomp.gt.0)then
            icomp=icomp+1    ! increment index for subsequent zones
          endif
        endif
        nsz=nsz+1

C Scan for zone name, convert any unallowed characters via st2name.
        CALL EGETW(OUTSTR,K,WORD,'W','Zone name',IFLAG)
        call st2name(WORD,zname(nsz))
        goto 43
      elseif(WORD(1:8).eq.'*surface')then

C Expected order is the walls first then the top then the base and
C then any windows and then any doors.
        NS=NS+1
        if(metaver.eq.1.or.metaver.eq.2)then

C Surface name, allow for future spaces in name.
          CALL EGETP(OUTSTR,K,WORD,'W','surface name',IER)
          write(sname(nsz,NS),'(a)') WORD(1:lnblnk(WORD))

C Surface construction name, allow for spaces.
          CALL EGETP(OUTSTR,K,WORD,'W','surface construction',IER)
          write(szhasconstr(nsz,NS),'(a)') WORD(1:lnblnk(WORD))

C If metaver is two or three then also read surface optical property.
          if(metaver.eq.2)then
            CALL EGETP(OUTSTR,K,WORD,'W','surface optics',IER)
            write(szhasoptic(nsz,NS),'(a)') WORD(1:lnblnk(WORD))
          endif

C Read three indices representing boundary condition.
          CALL EGETWI(OUTSTR,K,ival1,-1,7,'W','connection type',IER)
          CALL EGETWI(OUTSTR,K,ival2,0,MCOM,'W','connection ic2',IER)
          CALL EGETWI(OUTSTR,K,ival3,0,MS,'W','connection ie2',IER)
          zboundarytype(nsz,NS,1)=ival1
          zboundarytype(nsz,NS,2)=ival2
          zboundarytype(nsz,NS,3)=ival3
          szhasparent(nsz,NS) = '-'
          szhasuse(nsz,NS,1) = '-'
          szhasuse(nsz,NS,2) = '-'
        elseif(metaver.eq.3)then

C Surface attributes in a typical line:
C *surf,door,VERT,-,DOOR,UNDERCUT,door,OPAQUE,ANOTHER,3,6 # >|< door in coridor

C (1) surf name (12 char, tbd - allow spaces)
C (2) surf posn (tags VERT, CEIL (faces down), FLOR (faces up), SLOP (other)
C (3) child of (name of surface, otherwise a -)
C (4-5) useage: a pair of tags for code complience and to assist with air
C  flow network creation. A full description is at the top of this source
C  file.
C (6) construction name (from MLC database)
C (7) optical name (OPAQUE or the name of the optical set name, (tbd allow spaces))
C (8-10) boundary condition (including indices from cnn file):
C   ANOTHER, i index of zone, j index of surface in that zone
C   EXTERIOR,i,j    where i is 0 and j is 0
C   SIMILAR,i,j     where i is offset degC and j is offset Watts
C   IDENT_CEN,i,j     where i is 0 and j is 0
C   CONSTANT,i,j     where i is degC and j is Watts
C   ADIABETIC,i,j     where i is 0 and j is 0
C   BASESIMP,i,j     where i is configuration index and j is % to this surface
C   GROUND_STD,i,j     where i is profile index and j is 0
C   GROUND_USR,i,j     where i is 0 and j is profile index
C   UNKNOWN,i,j     where i is 0 and j is 0
          CALL EGETP(OUTSTR,K,WORD,'W','surface name',IER)
          write(sname(nsz,NS),'(a)') WORD(1:lnblnk(WORD))
          CALL EGETW(OUTSTR,K,tmpvfc,'W','surface position',IER)
          IF(tmpvfc.EQ.'VERT'.OR.tmpvfc.EQ.'SLOP')then
            write(SVFC(nsz,NS),'(a)') tmpvfc
          ELSEIF(tmpvfc(1:4).EQ.'CEIL')then
            write(SVFC(nsz,NS),'(a)') 'CEIL'
            izsceil(nsz)=NS   ! identify as a ceiling
          ELSEIF(tmpvfc(1:4).EQ.'FLOR')then
            write(SVFC(nsz,NS),'(a)') 'FLOR'
            izsfloor(nsz)=NS  ! identify as a floor
          ENDIF
          CALL EGETP(OUTSTR,K,WORD,'W','surface parent',IER)
          write(szhasparent(nsz,NS),'(a)') WORD(1:lnblnk(WORD))
          write(sparent(nsz,NS),'(a)') WORD(1:lnblnk(WORD))
          CALL EGETP(OUTSTR,K,WORD,'W','surface use 1',IER)
          if(lnblnk(WORD).gt.12)then
            write(SUSE(nsz,NS,1),'(a)') '-'
            write(SUSE(nsz,NS,2),'(a)') '-'
            write(SMLCN(nsz,NS),'(a)') WORD(1:lnblnk(WORD))
          else
            write(SUSE(nsz,NS,1),'(a)') WORD(1:lnblnk(WORD))
            CALL EGETW(OUTSTR,K,WORD,'W','surface use 2',IER)
            write(SUSE(nsz,NS,2),'(a)') WORD(1:lnblnk(WORD))

C Surface construction name, allow for spaces.
            CALL EGETP(OUTSTR,K,WORD,'W','surface construction',IER)
          endif
          szhasuse(nsz,NS,1)=SUSE(nsz,NS,1)
          szhasuse(nsz,NS,2)=SUSE(nsz,NS,2)
          write(SMLCN(nsz,NS),'(a)') WORD(1:lnblnk(WORD))
          write(szhasconstr(nsz,NS),'(a)') WORD(1:lnblnk(WORD))
          CALL EGETP(OUTSTR,K,WORD,'W','surface optics',IER)
          write(SOTF(nsz,NS),'(a)') WORD(1:lnblnk(WORD))
          write(szhasoptic(nsz,NS),'(a)') WORD(1:lnblnk(WORD))

          CALL EGETW(OUTSTR,K,tother,'W','surface other 1',IER)
          CALL EGETW(OUTSTR,K,tother1,'W','surface other 2',IER)
          CALL EGETW(OUTSTR,K,tother2,'W','surface other 3',IER)
          if(tother(1:7).eq.'UNKNOWN')then
            write(sother(nsz,ns,1),'(a)') 'UNKNOWN'
            write(sother(nsz,ns,2),'(a)') '00'
            write(sother(nsz,ns,3),'(a)') '000'
            zboundarytype(nsz,NS,1)= -1
            zboundarytype(nsz,NS,2)=0
            zboundarytype(nsz,NS,3)=0
          elseif(tother(1:8).eq.'EXTERIOR')then
            write(sother(nsz,ns,1),'(a)') 'EXTERIOR'
            write(sother(nsz,ns,2),'(a)') '00'
            write(sother(nsz,ns,3),'(a)') '000'
            zboundarytype(nsz,NS,1)=0
            zboundarytype(nsz,NS,2)=0
            zboundarytype(nsz,NS,3)=0
          elseif(tother(1:9).eq.'ADIABATIC')then
            write(sother(nsz,ns,1),'(a)') 'ADIABATIC'
            write(sother(nsz,ns,2),'(a)') '00'
            write(sother(nsz,ns,3),'(a)') '000'
            zboundarytype(nsz,NS,1)=5
            zboundarytype(nsz,NS,2)=0
            zboundarytype(nsz,NS,3)=0
          elseif(tother(1:7).eq.'SIMILAR')then
            write(sother(nsz,ns,1),'(a)') 'SIMILAR'
            write(sother(nsz,ns,2),'(a)') tother1
            write(sother(nsz,ns,3),'(a)') tother2
            zboundarytype(nsz,NS,1)=1
            read(tother1,*) zboundarytype(nsz,NS,2)
            read(tother2,*) zboundarytype(nsz,NS,3)
          elseif(tother(1:8).eq.'CONSTANT')then
            write(sother(nsz,ns,1),'(a)') 'CONSTANT'
            write(sother(nsz,ns,2),'(a)') tother1
            write(sother(nsz,ns,3),'(a)') tother2
            zboundarytype(nsz,NS,1)=2
            read(tother1,*) zboundarytype(nsz,NS,2)
            read(tother2,*) zboundarytype(nsz,NS,3)
          elseif(tother(1:8).eq.'BASESIMP')then
            write(sother(nsz,ns,1),'(a)') 'BASESIMP'
            write(sother(nsz,ns,2),'(a)') tother1
            write(sother(nsz,ns,3),'(a)') tother2
            zboundarytype(nsz,NS,1)=6
            read(tother1,*) zboundarytype(nsz,NS,2)
            read(tother2,*) zboundarytype(nsz,NS,3)
          elseif(tother(1:6).eq.'GROUND')then
            write(sother(nsz,ns,1),'(a)') 'GROUND'
            write(sother(nsz,ns,2),'(a)') tother1
            write(sother(nsz,ns,3),'(a)') tother2
            zboundarytype(nsz,NS,1)=4
            read(tother1,*) zboundarytype(nsz,NS,2)
            read(tother2,*) zboundarytype(nsz,NS,3)
          elseif(tother(1:9).eq.'IDENT_CEN')then
            write(sother(nsz,ns,1),'(a)') 'IDENT_CEN'
            write(sother(nsz,ns,2),'(a)') tother1
            write(sother(nsz,ns,3),'(a)') tother2
            zboundarytype(nsz,NS,1)=7
            read(tother1,*) zboundarytype(nsz,NS,2)
            read(tother2,*) zboundarytype(nsz,NS,3)
          elseif(tother(1:7).eq.'ANOTHER')then
            write(sother(nsz,ns,1),'(a)') 'ANOTHER'
            write(sother(nsz,ns,2),'(a)') tother1
            write(sother(nsz,ns,3),'(a)') tother2
            zboundarytype(nsz,NS,1)=3
            read(tother1,*) zboundarytype(nsz,NS,2)
            read(tother2,*) zboundarytype(nsz,NS,3)
          else

C Assume partition so setup sother arrays with this in mind.
            write(sother(nsz,ns,1),'(a)') 'ANOTHER'
            write(sother(nsz,ns,2),'(a)') tother1
            write(sother(nsz,ns,3),'(a)') tother2
            zboundarytype(nsz,NS,1)=3
            read(tother1,*) zboundarytype(nsz,NS,2)
            read(tother2,*) zboundarytype(nsz,NS,3)
          endif
        endif

C Note: if shape=poly for the each of the declaired walls the
C following line should start with the token *list.
        goto 43
      elseif(WORD(1:6).eq.'*shape')then
        CALL EGETW(OUTSTR,K,WORD,'W','shape type',IFLAG)
        if(WORD(1:3).eq.'box')then
          write(shape(nsz),'(a)') 'box '
          nbwalls(nsz)=6
        elseif(WORD(1:7).eq.'extrude')then
          write(shape(nsz),'(a)') 'extrude'
        elseif(WORD(1:4).eq.'poly')then
          write(shape(nsz),'(a)') 'poly'
        endif
        NS=0     ! reset number of surfaces and clear arrays.
        goto 43

      elseif(WORD(1:6).eq.'*usage')then

C After *usage the next token is 'pattern' (tag indicating that
C an operations file in the training/pattern folder to use to
C get infiltration and casual gains). Note usefile is only
C 32 char buffer because other source code knowns the path
C to the training/pattern folder.
        CALL EGETW(OUTSTR,K,use(1),'W','pattern',IFLAG)
        if(use(1)(1:7).eq.'pattern')then
          szuse(nsz,1) = use(1)
          CALL EGETW(OUTSTR,K,use(2),'W','pattern infil',IFLAG)
          szuse(nsz,2) = use(2)
          CALL EGETW(OUTSTR,K,use(3),'W','pattern casual',IFLAG)
          szuse(nsz,3) = use(3)
          CALL EGETRM(OUTSTR,K,WORD,'W','pattern file',IER)
          write(usefile(nsz),'(a)') WORD(1:lnblnk(WORD))
        else
          call usrmsg('Usage 2nd parameter not recognised',' ','W')
        endif
        goto 43
      elseif(WORD(1:4).eq.'*vew')then
        CALL EGETWR(OUTSTR,K,EYEMI(1),0.,0.,'-','X cord',IER)
        CALL EGETWR(OUTSTR,K,EYEMI(2),0.,0.,'-','Y cord',IER)
        CALL EGETWR(OUTSTR,K,EYEMI(3),0.,0.,'-','Z cord',IER)
        CALL EGETWR(OUTSTR,K,VIEWMI(1),0.,0.,'-','X cord',IER)
        CALL EGETWR(OUTSTR,K,VIEWMI(2),0.,0.,'-','Y cord',IER)
        CALL EGETWR(OUTSTR,K,VIEWMI(3),0.,0.,'-','Z cord',IER)
        CALL EGETWR(OUTSTR,K,ANGI,0.,0.,'-','Z cord',IER)
        initvt=1
        goto 43
      elseif(WORD(1:15).eq.'*water_in_zones')then
        K=16
        CALL EGETWR(OUTSTR,K,zfldK,0.,0.,'-','wtr cond',IER)
        CALL EGETWR(OUTSTR,K,zfldD,0.,0.,'-','wtr dens',IER)
        CALL EGETWR(OUTSTR,K,zfldC,0.,0.,'-','wtr spec ht',IER)
        CALL EGETWR(OUTSTR,K,zfldA,0.,0.,'-','wtr sw abs',IER)
        goto 43
      elseif(WORD(1:16).eq.'*weather_station')then

C If there is a line with the root name of a climate file look
C for a match in the current climate location and if there is
C one set the name of that file, scan it and find its site.
        CALL EGETRM(OUTSTR,K,weather,'W','weather name',IER)
        goto 43
      elseif(WORD(1:9).eq.'*occupant')then
        goto 43   ! not yet implemented
      elseif(WORD(1:9).eq.'*lighting')then
        goto 43   ! not yet implemented
      elseif(WORD(1:10).eq.'*equipment')then
        goto 43   ! not yet implemented
      elseif(WORD(1:12).eq.'*environment')then
        goto 43   ! not yet implemented
      else
        goto 43
      endif

C Now close silent data file.
      CALL ERPFREE(IUNIT,ios)
      RETURN

C Process the data if the end of the file reached or if the
C *end mark was found.
 1000 CALL ERPFREE(IUNIT,ios)

C Create new model (registration). Note the onfiguration file will be in the
C folder /tmp/box/cfg if mpath is /tmp/box.
      write(metaroot,'(a)') root(1:lnblnk(root))
      write(metampath,'(a)') mpath(1:lnblnk(mpath))
      call silentmodel(actions,root,mpath,weather,simact,ier)

C Add in site information if tokens were included.
      if(havesite)then
        SLAT=sitelat
        SLON=sitelongdif
      endif
      if(havehourlysnowfile)then
        SNFNAM=hourlysnowfile
      endif
      if(havegroundtempsets)then
        NGRDP=groundtempsets
        do 128 IGRDP=1,NGRDP
          DO 129 J=1,12
            UGRDTP(J,IGRDP)=groundtemps(J,IGRDP)
  129     CONTINUE
  128   continue
      endif

C Just to be sure re-scan the MLC database so that the construction
C attributes array are filled.
      call opendb(ier)

C << location for further system level data >>

C If there are no zones and no *start_zone found return.
      if(ncomp.eq.0.and.icomp.eq.0)then
        return
      else
        if(ncomp.eq.0)then
          icomp=0      ! reset because it will be incrmented later
        elseif(ncomp.gt.0)then
          icomp=ncomp  ! for existing model set icomp to next avail zone
        endif
        do 63 isz=1,nsz

C Call silent zone and then loop back for more possible definiitons.
C Copy from sz* arrays into parameter arrays prior to silentzone call.
          use(1)= szuse(isz,1)
          use(2)= szuse(isz,2)
          use(3)= szuse(isz,3)
C        usefile(iz) = ' '
          do 49 isu=1,MS
            hasconstr(isu)=szhasconstr(isz,isu)
            hasoptic(isu)=szhasoptic(isz,isu)
            hasparent(isu)=szhasparent(isz,isu)
            hasuse(isu,1)=szhasuse(isz,isu,1)
            hasuse(isu,2)=szhasuse(isz,isu,2)
  49      continue
          icomp=icomp+1

C If there are obstructions fill in relevant common blocks. Set iobs()=2
C to signal that the common blocks have been instantiated and data should
C be included in the version 1.1 geometry file.
          if(nbobs(icomp).gt.0)then
            iobs(icomp)=2
          else
            iobs(icomp)=0
            nbobs(icomp)=0
          endif

C If a literal control file has been included do not bother with terse
C control description.
          if(literalctl)then
            continue
          else
            if(i_ctl_link(isz).gt.0)then
              isilentncf=1+isilentncf
              ncf=isilentncf
              IBSN(isilentncf,1)=0; IBSN(isilentncf,2)=0
              IBSN(isilentncf,3)=0; IBSN(isilentncf,4)=0
              IBAN(isilentncf,1)=0; IBAN(isilentncf,2)=0
              IBAN(isilentncf,3)=0
              NBCDT(isilentncf)=1
              IBCDV(isilentncf,1,1)=1; IBCDV(isilentncf,1,2)=365
              NBCDP(isilentncf,1)=1
              IBCTYP(isilentncf,1,1)=0
              IBCLAW(isilentncf,1,1)=1
              TBCPS(isilentncf,1,1)=0
              BMISCD(isilentncf,1,1,1)=7
              BMISCD(isilentncf,1,1,2)=999000
              BMISCD(isilentncf,1,1,3)=0.0
              BMISCD(isilentncf,1,1,4)=999000
              BMISCD(isilentncf,1,1,5)=0.0
              BMISCD(isilentncf,1,1,6)= ht_Setpoint(isz)
              BMISCD(isilentncf,1,1,7)=cl_Setpoint(isz)
              BMISCD(isilentncf,1,1,8)=0.0
              znctldoc='basic ideal controls'
              ICASCF(isz)=ncf
            endif
          endif

          write(outs,'(2a)') 'Processing ',zname(isz)
          call edisp(iuout,outs)

C Create zone based on a set of parameters (created above).
          call silentzone(ICOMP,metaver,hasconstr,hasoptic,hasparent,
     &      hasuse,use,usefile(isz),literalctl,IER)
  63    continue

        if(unixok)then
          fs = char(47)
        else
          fs = char(92)
        endif
        if(.NOT.literalctl)then  ! No control embedded in meta file.
          LN=max(1,LNBLNK(cfgroot))
          if(ctlpth(1:2).eq.'  '.or.ctlpth(1:2).eq.'./')then
            WRITE(LCTLF,'(2a)')cfgroot(1:ln),'.ctl'
          else
            WRITE(LCTLF,'(4a)') ctlpth(1:lnblnk(ctlpth)),fs,
     &         cfgroot(1:ln),'.ctl'
          endif
          ICTLF=IFIL+1
          call CTLWRT(ICTLF,IER)
        endif

C Update cfg file to know about the control file
        CALL EMKCFG('s',IER)
      endif
      RETURN

 1001 write(outs,'(3a)') 'SILENTREAD: conversion error in...',
     &  OUTSTR(1:50),'...'
      lsn=MIN0(lnblnk(currentfile),110)
      write(outs2,'(2a)') 'in: ',currentfile(1:lsn)
      call edisp(iuout,outs)
      call edisp(iuout,outs2)
      IER=1
      CALL ERPFREE(IUNIT,ios)
      RETURN

      END


C ******************** METAWRITE
C METAWRITE creates an esp-r meta file. Data input is from
C the current zones common blocks. If act is .... Otherwise create
C the file with minimum interaction.

      SUBROUTINE METAWRITE(ifilg,LFILE,act,IER)
      integer MSZ   ! number of META zones array sizes, edit to
                    ! match MCOM in building.h
      PARAMETER (MSZ=82)

#include "building.h"
#include "model.h"
#include "site.h"
#include "geometry.h"
#include "esprdbfile.h"
#include "MultiYear_simulations.h"
#include "espriou.h"
#include "seasons.h"
#include "ipvdata.h"
#include "FMI.h"
#include "agent.h"
#include "gremlin.h"

      integer lnblnk  ! function definition

      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      common/FILEP/IFIL

      integer ncomp,ncon
      COMMON/C1/NCOMP,NCON
      COMMON/C3/IC1(MCON),IE1(MCON),ICT(MCON),IC2(MCON),IE2(MCON)

C Obstruction blocks via geometry.h.
      common/rpath/path
      character path*72    ! where path folders
      character odir*84    ! where current working folder

C Simulation parameter sets.
      common/spfldat/nsset,isset,isstup,isbnstep,ispnstep,issave,isavgh
      INTEGER :: nsset,isset,isstup,isbnstep,ispnstep,issave,isavgh

C IPV description via ipvdata.h.

      common/cctlnm/ctldoc,lctlf
      character LCTLF*72,CTLDOC*248

      integer icascf
      common/cctl/icascf(mcom)

      integer IZSTOCN
      COMMON/C24/IZSTOCN(MCOM,MS)

C XMLMETRICS are the metrics that could be included in the input.xml
C in case the flag: <dump_all_data> is false
      common/XMLDATA/XMLMETRICS(200),inumXMLmetrics
      character XMLMETRICS*248

C iaplic(1) toggle for shading; iaplic(2) toggle for insolation.
C   toggle for assessment where:
C   1 is 'all_applicable', 0 is manual selection of surfaces.
C nsurfcalc nb of shaded surfaces, lstsfcalc() list of applicable surfaces.
C nsurfinso nb of insolation sources, isurfinso() list of insolation sources.
      common/ishdirec/iaplic(MCOM,2),nsurfcalc(MCOM),lstsfcalc(MCOM,MS),
     &       nsurfinso(MCOM),isurfinso(MCOM,MS)

C Images.
      character imgfmt*4  ! GIF XBMP TIF JPG
      character imgfoc*4  ! FZON FNET FCTL FDFS FPLN
      character limgfil*72  ! file name (extend to 144 char)
      character imgdoc*248  ! text associated with image
      common/imagf/imgfmt(MIMG),imgfoc(MIMG),limgfil(MIMG),imgdoc(MIMG)

      integer noimg  ! number of images
      integer iton   ! 0 if images not yet shown, 1 if yes
      common/imagfi/noimg,iton

      common/calena/calename,calentag(MDTY),calendayname(MDTY)
      character calename*32,calentag*12,calendayname*32

      integer radcores
      common/radcor/radcores

      COMMON/zfluid/znotair(mcom),zfldK,zfldD,zfldC,zfldA,
     &              zSWAp(mcom),zSWAf(mcom)
      real zfldK,zfldD,zfldC,zfldA,zSWAp,zSWAf
      LOGICAL znotair              ! zone can be air or water filled
      
      integer IHCT,ICOREXHCT
      COMMON/HCTHRY/IHCT,ICOREXHCT  ! Default for internal & external hc method.

      common/CWEC_SOL/iSlr_half_hr_flg
      integer iSlr_half_hr_flg      ! 0 hour-centered; 1 half-hour centered

C Ground topology.
      COMMON/GTFIL/GTGEOM
      character GTGEOM*72
      COMMON/GRND100/GRND3D
      LOGICAL GRND3D

C Initial view.
      common/initv/initvt,EYEMI(3),VIEWMI(3),ANGI
      REAL :: angi,EYEMI,VIEWMI
      INTEGER :: initvt

      logical dogrnd
 
      integer inumXMLmetrics !counts for the number of xml metrics defined

      CHARACTER outs*124,outsd*124,louts*496,loutsd*496
      character dstmp*24,OTHSTR*30

      character act*1      ! action requested.
      character simact*6   ! action for creating assessments.
      character subpath*72 ! not yet used

C At the head of the file is a menu string and a description
C modeltitle (char 72) and modeldocblock (char 248) which identifies
C the contents of the file via model.h.
      character weather*32     ! file name for the location (no path)
      character cpath*72,cfile*72 ! for finding the climate file name
      character*(*) lfile      ! name of the file
      character fs*1           ! file separator

C For passing to silentmodel.
      character actions*8,root*32,mpath*72

      logical unixok

C Paremeters passed to:silentzone
      character szuse*8,usefile*32  ! usage pattern directives
      dimension szuse(MSZ,3),usefile(MSZ)

C Site related local variables.
      logical havesnowdays     ! true if days with snow defined
      character hourlysnowfile*72  ! same as SNFNAM
      logical havehourlysnowfile   ! true if set
      real hoursGTM   ! hours before or after GTM
      real groundtemps(12,MGRDP)
      integer groundtempsets
      logical havegroundtempsets
      logical havesimparameters
      logical rotationknown  ! true if we need to remember a rotation
      logical :: XST
      integer simstartup   ! use for isstup
      integer simzonetimestep ! use for isbnstep
      integer simplanttimestep ! use for ispnstep
      integer simsavelevel  ! use for issave
      integer i_ctl_link  !flag to indicate the zone has a basic control
      real ht_Setpoint,cl_Setpoint !heating and cooling set-points
      dimension i_ctl_link(MCOM),ht_Setpoint(MCOM),cl_Setpoint(MCOM)
      integer isilentncf !silent decide how many control functions needed
      character tab*1
      character tokens*156,comment*76,aligned_str*156
      character left*62,leftd*62,right*42,rightd*42
      integer lnc  ! for length of climate file
      integer metaver  ! one for 1.1 two for 1.2 three for 1.3
      integer ifc

C  The following data statements carry the proportions of
C  surrounding buildings, sky vault and ground assumed
C  to be seen by a vertical surface for a particular site
C  defined by the Index of Exposure.
      DIMENSION PROPB(8),PROPS(8),PROPG(8)
      DATA PROPB/0.28,0.18,0.1,0.34,0.52,0.0,1.0,0.0/
      DATA PROPS/0.36,0.41,0.45,0.33,0.15,0.5,0.0,0.0/
      DATA PROPG/0.36,0.41,0.45,0.33,0.33,0.5,0.0,0.0/

      IER=0; XST=.false.

C Determine operating system
      call isunix(unixok)
      if(unixok)then
        fs = char(47)
      else
        fs = char(92)
      endif
      call dstamp(dstmp)

C Get the current folder and display options to the user.
C first list any files with .cfg in name.
      odir=' '
      call usrdir(odir)

C If subpath is less than 72 char with odir + path use that
C otherwise just use path.

C Debug...
C      write(6,*) 'usrdir is ',odir
C      write(6,*) 'path is ',path

      subpath=' '
      lnodir=lnblnk(odir)
      lnpath=lnblnk(path)
      if(lnodir+lnpath+1.le.72)then
        write(subpath,'(3a)') odir(1:lnodir),fs,
     &  path(1:lnpath)
      else
        write(subpath,'(a)') path(1:lnpath)
      endif

C Indicate progress.
      write(outs,'(2a)') 'configuration files in ',
     &  subpath(1:lnblnk(subpath))
      call edisp(iuout,outs)

C Instantiate the sz array structures. The logic is to scan all of the
C zone details and then write into the file.
      write(actions,'(a)') 'new'
      write(root,'(a)') cfgroot
      write(mpath,'(a)') cfgroot  ! model folder will take root name

C Currently set the block of text to the model title.
      write(modeldocblock,'(a)')  modeltitle(1:lnblnk(modeltitle))

C Site data.
C If standard or local set weather to lclim otherwise remove the
C path from the climate file and set the file name
C only (the first 32 char) into weather.
      if(ipathclim.eq.1.or.ipathclim.eq.2)then  ! standard or ../dbs
        lnc=lnblnk(lclim)
        if(lnc.lt.32)then
          write(weather,'(a)') lclim(1:lnc)
        else
          write(weather,'(a)') lclim(1:31)
        endif
      else
        call fdroot(LCLIM,cpath,cfile)
        lnc=lnblnk(cfile)
        if(lnc.lt.32)then
          write(weather,'(a)') cfile(1:lnc)
        else
          write(weather,'(a)') cfile(1:31)
        endif
      endif
      hoursGTM=0.0

C The hourlysnowfile is SNFNAM.
      write(hourlysnowfile,'(a)') SNFNAM

C Depending on the Index of Exposure the proportions of the
C surroundings are assigned new variable names.
      if(siteexposureindex.eq.8)then
        continue
      else
        buildingview=PROPB(siteexposureindex)
        skyview=PROPS(siteexposureindex)
        groundview=PROPG(siteexposureindex)
      endif

C Number of days with snow on the ground (if ground refl model 2)
      if(groundreflmodel.EQ.2)then
        havesnowdays=.true.
      endif

C Snow depth file (if ground reflectivity model 3)
      if(groundreflmodel.EQ.3)then
        write(hourlysnowfile,'(a)') SNFNAM(1:lnblnk(SNFNAM))
        havehourlysnowfile=.true.
      endif

      if(NGRDP.gt.0)then
        havegroundtempsets=.true.
        groundtempsets=NGRDP
        do 128 IGRDP=1,NGRDP
          DO 129 J=1,12
            groundtemps(J,IGRDP)=UGRDTP(J,IGRDP)
  129     CONTINUE
  128   continue
      endif

      dogrnd = .false.
      IF(NGRDP.gt.0.or.GRND3D)dogrnd=.true.
      IF(GTGEOM(1:7).ne.'UNKNOWN'.and.GTGEOM(1:2).ne.'  ')dogrnd=.true.

C Simulation parameter sets.
      if(nsset.gt.0)then
        havesimparameters=.true.
        simstartup=isstup
        simzonetimestep=isbnstep
        simplanttimestep=ispnstep
        simsavelevel=issave
        isilentncf=0   ! still to do
      endif

C For each of the zones instantiate values.
C << todo - add in zone description >>
      do 42 iz=1,ncomp
        shape(iz) = 'poly'    ! zone shape is poly
        szuse(iz,1) = ' '     ! szuse (skip for now)
        szuse(iz,2) = ' '
        szuse(iz,3) = ' '
        usefile(iz) = ' '       ! usefile (skip for now)

C Until we store all META component history the creation of a META
C file treats everything as polygons so that nbwalls is the same as
C nzsur. For some 3rd party sources nbwalls will be less than nzsur.
        nbwalls(iz) = NZSUR(IZ) ! surfaces in the zone

C Control data - ignore for now.
        ht_Setpoint(iz)= 0.0
        cl_Setpoint(iz)= 0.0
        i_ctl_link(iz)=0

C Obstructions - if the original model has obstructions then
C we will use them.

  42  continue

C Initialise geometry data file. and set currentfile.
      CALL EFOPSEQ(ifilg,LFILE,3,IER)
      IF(IER.LT.0)THEN
        write(outs,'(3a)') 'Meta file ',LFILE(1:lnblnk(LFILE)),
     &      ' could not be opened.'
        call edisp(iuout,outs)
        IER=1
        RETURN
      endif
      write(currentfile,'(a)') LFILE(1:lnblnk(LFILE))

      tab=','    ! make the separator a comma.
C Write header lines from file, the 1.1 version looks like:
C *silent_input
C      write(ifilg,'(a)',IOSTAT=ios,ERR=13) '*silent_input'

C The version 1.2 (which includes surface optical attributes looks like:
C *meta_file. version 1.3 which has full set of surface atrributes looks
C like *meta_file_v3
      metaver =3
      write(ifilg,'(a)',IOSTAT=ios,ERR=13) '*meta_file_v3'
      write(ifilg,'(3a)',IOSTAT=ios,ERR=13)'*title',tab,
     &  modeltitle(1:lnblnk(modeltitle))
      write(ifilg,'(3a)',IOSTAT=ios,ERR=13)'*doc',tab,
     &  modeldocblock(1:lnblnk(modeldocblock))
      write(ifilg,'(4a)',IOSTAT=IOS,ERR=13)'*date',tab,dstmp,
     &  '  # latest file modification '
      write(ifilg,'(9a)',IOSTAT=IOS,ERR=13)'*action',tab,
     &  actions(1:lnblnk(actions)),tab,
     &  root(1:lnblnk(root)),tab,'distributed',tab,
     &  mpath(1:lnblnk(mpath))

C Include the model folder tags as in the cfg file. This will be used when
C extracting to a set of standard ESP-r files.
      ifc=ifilg
      WRITE(IFC,'(A)',IOSTAT=IOS,ERR=13) '#'
      WRITE(IFC,'(A)',IOSTAT=IOS,ERR=13) '# ---- model folders ----'
      call write_folders(ifc)

C Do something about monthly snow data etc.
      if(havesnowdays)then
        WRITE(ifilg,804,IOSTAT=IOS,ERR=14)
     &    '*monthly_snow_days',tab,
     &    (dayswithsnow(I),I=1,12)
  804   FORMAT(2A,13(1X,I2))
      endif
      if(havehourlysnowfile)then
        write(ifilg,'(4a)',IOSTAT=IOS,ERR=13)'*hourly_snow_depth',tab,
     &    hourlysnowfile(1:lnblnk(hourlysnowfile)),' # snow depth file'
      endif

C Ground temperatures (monthly profiles). Header followed by one
C line of 12 values for each included profile.
      if(havegroundtempsets)then
        write(ifilg,'(2a,I3,a)',IOSTAT=IOS,ERR=13)
     &    '*monthly_profile',tab,groundtempsets,
     &    '  # nb of ground temperature profiles'
        DO 28 IGRDP=1,groundtempsets
          WRITE(ifilg,933,IOSTAT=IOS,ERR=14)
     &      (groundtemps(J,IGRDP),J=1,12)
  933     FORMAT(12F6.2)
  28    CONTINUE
      endif

C Write common data files related to this model. Use logic
C similar to subroutine write_databases but for local MLC
C database include literally. 
      ifc=ifilg
C      call write_databases(ifc,ier)

C Start of code block similar to emkcfg.F: write_databases
      ier=0               ! Reset so write errors detected on return.
      WRITE(IFC,'(A)',IOSTAT=IOS,ERR=2) '#'
      WRITE(IFC,'(A)',IOSTAT=IOS,ERR=2) '# ---- databases ----'
      if(ipathmat.eq.0.or.ipathmat.eq.1)then
        WRITE(IFC,'(2A)',IOSTAT=IOS,ERR=2) '*mat ',
     &      LFMAT(1:lnblnk(LFMAT))
        ltf=max(1,LNBLNK(LFMAT))
        INQUIRE (FILE=LFMAT(1:ltf),EXIST=XST)
        if (XST) THEN
          write(ifilg,'(a)') '*start_literal'
          iuf=ifil+2
          call ascitometa(iuf,LFMAT,ifilg,IER)
          write(ifilg,'(a)') '*end_literal'
        endif
      elseif(ipathmat.eq.2)then
        WRITE(IFC,'(2A)',IOSTAT=IOS,ERR=2) '*stdmat ',
     &      LFMAT(1:lnblnk(LFMAT))
      endif

      if(ipathcfc.eq.0.or.ipathcfc.eq.1)then
        WRITE(IFC,'(2A)',IOSTAT=IOS,ERR=2) '*cfcdb ',
     &      LCFCDB(1:lnblnk(LCFCDB))
      elseif(ipathcfc.eq.2)then
        WRITE(IFC,'(2A)',IOSTAT=IOS,ERR=2) '*stdcfcdb ',
     &      LCFCDB(1:lnblnk(LCFCDB))
      endif

      if(ipathmul.eq.0.or.ipathmul.eq.1)then
        WRITE(IFC,'(2A)',IOSTAT=IOS,ERR=2) '*mlc ',
     &      LFMUL(1:lnblnk(LFMUL))
        ltf=max(1,LNBLNK(LFMUL))
        INQUIRE (FILE=LFMUL(1:ltf),EXIST=XST)
        if (XST) THEN
          write(ifilg,'(a)') '*start_literal'
          iuf=ifil+2
          call ascitometa(iuf,LFMUL,ifilg,IER)
          write(ifilg,'(a)') '*end_literal'
        endif
      elseif(ipathmul.eq.2)then
        WRITE(IFC,'(2A)',IOSTAT=IOS,ERR=2) '*stdmlc ',
     &      LFMUL(1:lnblnk(LFMUL))
      endif

      if(ipathoptdb.eq.0.or.ipathoptdb.eq.1)then
        WRITE(IFC,'(2A)',IOSTAT=IOS,ERR=2) '*opt ',
     &      LOPTDB(1:lnblnk(LOPTDB))
        ltf=max(1,LNBLNK(LOPTDB))
        INQUIRE (FILE=LOPTDB(1:ltf),EXIST=XST)
        if (XST) THEN
          write(ifilg,'(a)') '*start_literal'
          iuf=ifil+2
          call ascitometa(iuf,LOPTDB,ifilg,IER)
          write(ifilg,'(a)') '*end_literal'
        endif
      elseif(ipathoptdb.eq.2)then
        WRITE(IFC,'(2A)',IOSTAT=IOS,ERR=2) '*stdopt ',
     &      LOPTDB(1:lnblnk(LOPTDB))
      endif

      if(ipathapres.eq.0.or.ipathapres.eq.1)then
        WRITE(IFC,'(2A)',IOSTAT=IOS,ERR=2) '*prs ',
     &      LAPRES(1:lnblnk(LAPRES))
      elseif(ipathapres.eq.2)then
        WRITE(IFC,'(2A)',IOSTAT=IOS,ERR=2) '*stdprs ',
     &      LAPRES(1:lnblnk(LAPRES))
      endif

      if(ipathprodb.eq.0.or.ipathprodb.eq.1)then
        WRITE(IFC,'(2A)',IOSTAT=IOS,ERR=2) '*evn ',
     &      LPRFDB(1:lnblnk(LPRFDB))
      elseif(ipathprodb.eq.2)then
        WRITE(IFC,'(2A)',IOSTAT=IOS,ERR=2) '*stdevn ',
     &      LPRFDB(1:lnblnk(LPRFDB))
      endif

      if(ipathclim.eq.0.or.ipathclim.eq.1)then
        WRITE(IFC,'(2A)',IOSTAT=IOS,ERR=2) '*clm ',
     &      LCLIM(1:lnblnk(LCLIM))
      elseif(ipathclim.eq.2)then
        WRITE(IFC,'(2A)',IOSTAT=IOS,ERR=2) '*stdclm ',
     &      LCLIM(1:lnblnk(LCLIM))
      endif

      if(ipathmsc.eq.0.or.ipathmsc.eq.1)then
        WRITE(IFC,'(2A)',IOSTAT=IOS,ERR=2) '*mscldb ',
     &      MCMPDBFL(1:lnblnk(MCMPDBFL))
      elseif(ipathmsc.eq.2)then
        WRITE(IFC,'(2A)',IOSTAT=IOS,ERR=2) '*stdmscldb ',
     &      MCMPDBFL(1:lnblnk(MCMPDBFL))
      endif
      if(ipathmould.eq.0.or.ipathmould.eq.1)then
        WRITE(IFC,'(2A)',IOSTAT=IOS,ERR=2) '*mould ',
     &    lfmould(1:lnblnk(lfmould))
      elseif(ipathmould.eq.2)then
        WRITE(IFC,'(2A)',IOSTAT=IOS,ERR=2) '*stdmould  ',
     &    lfmould(1:lnblnk(lfmould))
      endif

C Ground reflectivity model (no-snow monthly albedos & snow albedo
C if ground reflectivity model 2 or 3).
      IF(groundreflmodel.GT.1) WRITE(IFC,803) '*gref ',
     &   groundreflmodel,(groundreflmonth(I),I=1,12),snowgroundrefl
  803 FORMAT(A,I1,13(1X,F5.3))

C Number of days with snow on the ground (if ground reflectivity model 2).
      IF(groundreflmodel.EQ.2)  
     &  WRITE(IFC,805) '*snow',(dayswithsnow(I),I=1,12)
  805 FORMAT(A,12(1X,I2))

C Snow depth file (if ground reflectivity model 3).
      IF(groundreflmodel.EQ.3) WRITE(IFC,'(2A)',IOSTAT=IOS,ERR=2) 
     &  '*sndf ',SNFNAM(1:lnblnk(SNFNAM))

C Multi-year weather database.
      if (bMY_climates_defined) then 
        WRITE(IFC,'(2A)',IOSTAT=IOS,ERR=2) '*myclm ',
     &    cMY_climate_db_name(1:lnblnk(cMY_climate_db_name))
      endif

C Plant components database.
      if(ipathpcdb.eq.0.or.ipathpcdb.eq.1)then
        WRITE(IFC,'(2A)',IOSTAT=IOS,ERR=2) '*pdb ',
     &      LPCDB(1:lnblnk(LPCDB))
      elseif(ipathpcdb.eq.2)then
        WRITE(IFC,'(2A)',IOSTAT=IOS,ERR=2) '*stdpdb ',
     &      LPCDB(1:lnblnk(LPCDB))
      endif

C UK SBEM database.
      if(ipathsbem.eq.0.or.ipathsbem.eq.1)then
        WRITE(IFC,'(2a)',IOSTAT=IOS,ERR=2) '*sbem ',
     &    LSBEM(1:lnblnk(LSBEM))
      elseif(ipathsbem.eq.2)then
        WRITE(IFC,'(2a)',IOSTAT=IOS,ERR=2) '*stdsbem ',
     &    LSBEM(1:lnblnk(LSBEM))
      endif

C Predefined objects database.
      if(ipathpredef.eq.0.or.ipathpredef.eq.1)then
        WRITE(IFC,'(2a)',IOSTAT=IOS,ERR=2) '*predef ',
     &    LPREDEF(1:lnblnk(LPREDEF))
        ltf=max(1,LNBLNK(LPREDEF))
        INQUIRE (FILE=LPREDEF(1:ltf),EXIST=XST)
        if (XST) THEN
          write(ifilg,'(a)') '*start_literal'
          iuf=ifil+2
          call ascitometa(iuf,LPREDEF,ifilg,IER)
          write(ifilg,'(a)') '*end_literal'
        endif

      elseif(ipathpredef.eq.2)then
        WRITE(IFC,'(2a)',IOSTAT=IOS,ERR=2) '*stdpredef ',
     &    LPREDEF(1:lnblnk(LPREDEF))
      endif
C End of code block similar to emkcfg.F: write_databases

      
      WRITE(ifilg,'(A)',IOSTAT=IOS,ERR=13) '#'
      WRITE(ifilg,'(A)',IOSTAT=IOS,ERR=13) 
     &    '# ---- documentation and images ----'
      write(tokens,'(2a)') '*notes ',lmodellog(1:lnblnk(lmodellog))
      write(comment,'(a)') 'associated notes about the project'
      call align_comment(48,tokens,comment,aligned_str)
      write(ifilg,'(a)') aligned_str(1:lnblnk(aligned_str))

C If there is a notes file then literally copy it.
      ltf=max(1,LNBLNK(lmodellog))
      INQUIRE (FILE=lmodellog(1:ltf),EXIST=XST)
      if (XST) THEN
        write(ifilg,'(a)') '*start_literal'
        iuf=ifil+2
        call ascitometa(iuf,lmodellog,ifilg,IER)
        write(ifilg,'(a)') '*end_literal'
      endif
      
      if(lnblnk(lmodelqa).eq.0)then  ! Model contents.
        continue
      elseif(lmodelqa(1:7).eq.'UNKNOWN')then
        continue
      elseif(lmodelqa(1:2).eq.'  ')then
        continue
      else

C Contents file name set, check if the file exists.
        write(tokens,'(2a)') '*contents ',lmodelqa(1:lnblnk(lmodelqa))
        write(comment,'(a)') 'contents report for the project'
        call align_comment(48,tokens,comment,aligned_str)
        write(ifilg,'(a)') aligned_str(1:lnblnk(aligned_str))
        ltf=max(1,LNBLNK(lmodelqa))
        INQUIRE (FILE=lmodelqa(1:ltf),EXIST=XST)
        if (XST) THEN
          write(ifilg,'(a)') '*start_literal'
          iuf=ifil+2
          call ascitometa(iuf,lmodelqa,ifilg,IER)
          write(ifilg,'(a)') '*end_literal'
        endif
      endif

      if(noimg.gt.0)then
        do img=1,noimg
          WRITE(ifilg,'(a,a4,2x,a4,2x,a)',IOSTAT=IOS,ERR=13)'*img ',
     &      imgfmt(img),imgfoc(img),
     &      limgfil(img)(1:lnblnk(limgfil(img)))
          if(icfgv.gt.3)then
            WRITE(ifilg,'(2a)',IOSTAT=IOS,ERR=13)'*imdoc ',
     &        imgdoc(img)(1:lnblnk(imgdoc(img)))
          endif
        enddo
      endif

C Write the year and seasons.
      ifc=ifilg
      call write_seasons(ifc,ier)

C Also write the project calendar and day types.
      if(CALENAME(1:7).NE.'UNKNOWN')then
        ifc=ifilg
        call write_calendar(ifc,ier)
      endif

      WRITE(IFC,'(A)',IOSTAT=IOS,ERR=13) '#'
      WRITE(IFC,'(A)',IOSTAT=IOS,ERR=13) '# ---- site ----'
      write(tokens,'(a,f7.3,2x,f7.3)') '*latlong ',sitelat,sitelongdif
      write(comment,'(a)') 
     &  'latitude & longitude difference (from reference meridian)'
      call align_comment(48,tokens,comment,aligned_str)
      write(IFC,'(a)') aligned_str(1:lnblnk(aligned_str))

      write(tokens,'(a,i7,2x,f6.3)') '*site ',siteexposureindex,
     &  groundrefl
      if(siteexposureindex.EQ.1)then
        comment='exposure: typical city centre & ground reflectivity'
      elseif(siteexposureindex.EQ.2)then
        comment='exposure: typical urban & ground reflectivity'
      elseif(siteexposureindex.EQ.3)then
        comment='exposure: typical rural & ground reflectivity'
      elseif(siteexposureindex.EQ.4)then
        comment='city center, equal views to sky, ground & buildings'
      elseif(siteexposureindex.EQ.5)then
        comment='city center below mean height of other buildings'
      elseif(siteexposureindex.EQ.6)then
        comment='isolated rural & ground reflectivity'
      elseif(siteexposureindex.EQ.7)then
        comment='exposure: within a totally enclosed building'
      elseif(siteexposureindex.EQ.8)then
        comment='exposure: user defined & ground reflectivity '
      endif
      call align_comment(48,tokens,comment,aligned_str)
      write(IFC,'(a)') aligned_str(1:lnblnk(aligned_str))

      if(siteexposureindex.EQ.8)then
        write(tokens,'(3F8.3)') skyview,groundview,buildingview
        write(comment,'(a)') 
     &    'views to sky, ground & surrounding obstructions'
        call align_comment(48,tokens,comment,aligned_str)
        write(IFC,'(a)') aligned_str(1:lnblnk(aligned_str))
      endif
      write(tokens,'(a,f6.1)') '*sitealt ',sitealt
      write(comment,'(a)') 'site altitude (m)'
      call align_comment(48,tokens,comment,aligned_str)
      write(IFC,'(a)') aligned_str(1:lnblnk(aligned_str))

C Include ground information if associated attributes found.
      ifc=ifilg
      call write_ground(dogrnd,ifc)
      
      if(IHCT.eq.0) IHCT=1            ! if not defined assume A&H
      if(ICOREXHCT.eq.0) ICOREXHCT=2  ! if not defined assume MoWiTT
      write(tokens,'(a,2i5)') '*surf_heat_transfer ',IHCT,ICOREXHCT
      if(IHCT.eq.0.or.ICOREXHCT.eq.0)then
        write(comment,'(a)')'no methods for inside & outside htc'
      elseif(IHCT.eq.1.and.ICOREXHCT.eq.1)then
        write(comment,'(a)')'default methods for inside & outside htc'
      elseif(IHCT.eq.1.and.ICOREXHCT.eq.2)then
        write(comment,'(a)')'A&H inside MoWiTT outside htc'
      elseif(IHCT.eq.2.and.ICOREXHCT.eq.1)then
        write(comment,'(a)')'Khalifa rad inside McAdams outside htc'
      elseif(IHCT.eq.2.and.ICOREXHCT.eq.2)then
        write(comment,'(a)')'Khalifa rad inside MoWiTT outside htc'
      elseif(IHCT.eq.3.and.ICOREXHCT.eq.1)then
        write(comment,'(a)')'Khalifa no rad inside McAdams outsd htc'
      elseif(IHCT.eq.3.and.ICOREXHCT.eq.2)then
        write(comment,'(a)')'Khalifa no rad inside MoWiTT outsd htc'
      endif
      call align_comment(48,tokens,comment,aligned_str)
      write(IFC,'(a)') aligned_str(1:lnblnk(aligned_str))

      WRITE(IFC,'(A)',IOSTAT=IOS,ERR=13) '#'
      WRITE(IFC,'(A)',IOSTAT=IOS,ERR=13) 
     &  '# ---- model controls ----'
      write(tokens,'(2a)') '*ctl ',LCTLF(1:lnblnk(LCTLF))
      write(comment,'(a)') 'model control file'
      call align_comment(48,tokens,comment,aligned_str)
      write(IFC,'(a)') aligned_str(1:lnblnk(aligned_str))

C If there is a control file then literally copy it.
      ltf=max(1,LNBLNK(LCTLF))
      INQUIRE (FILE=LCTLF(1:ltf),EXIST=XST)
      if (XST) THEN
        write(IFC,'(a)') '*start_literal'
        iuf=ifil+2
        call ascitometa(iuf,LCTLF,ifilg,IER)
        write(ifilg,'(a)') '*end_literal'
      endif

C Simulation parameters.
      if(havesimparameters)then
        ifc=ifilg
        WRITE(IFC,'(A)',IOSTAT=IOS,ERR=13) '#'
        WRITE(IFC,'(A)',IOSTAT=IOS,ERR=13) 
     &    '# ---- simulation directives ----'
        write(tokens,'(a,i2)') '*radcore ',radcores
        write(comment,'(a)') 'cores available to Radiance'
        call align_comment(48,tokens,comment,aligned_str)
        write(IFC,'(a)') aligned_str(1:lnblnk(aligned_str))
        if(iSlr_half_hr_flg.eq.0)then
          write(tokens,'(a,i2)') '*slr_half_hr ',iSlr_half_hr_flg
          write(comment,'(a)') 'solar data hour centred'
          call align_comment(48,tokens,comment,aligned_str)
          write(IFC,'(a)') aligned_str(1:lnblnk(aligned_str))
        else
          write(tokens,'(a,i2)') '*slr_half_hr ',iSlr_half_hr_flg
          write(comment,'(a)') 'solar data half-hour centred'
          call align_comment(48,tokens,comment,aligned_str)
          write(IFC,'(a)') aligned_str(1:lnblnk(aligned_str))
        endif

        write(tokens,'(a,f4.1,3f7.1)') '*water_in_zones ',
     &    zfldK,zfldD,zfldC,zfldA
        write(comment,'(a)') 
     &    'water cond, density, spec heat, shortwave abs'
        call align_comment(48,tokens,comment,aligned_str)
        write(IFC,'(a)') aligned_str(1:lnblnk(aligned_str))
      
        if (is_FMU.and.FMISFL(1:7).ne.'UNKNOWN') then
          write(comment,'(a)') 'FMI specifications'
          write(tokens,'(2a)')'*FMI ',FMISFL(1:lnblnk(FMISFL))
          call align_comment(48,tokens,comment,aligned_str)
          write(IFC,'(a)') aligned_str(1:lnblnk(aligned_str))
        endif
        if (IS_AGT.and.AGTFIL(1:7).ne.'UNKNOWN') then  ! Agents.
          write(comment,'(a)') 'agent definitions'
          write(tokens,'(2a)') '*agent ',AGTFIL(1:lnblnk(AGTFIL))
          call align_comment(48,tokens,comment,aligned_str)
          write(IFC,'(a)') aligned_str(1:lnblnk(aligned_str))
        endif
        if (IS_GRM.and.GRMFIL(1:7).ne.'UNKNOWN') then  ! Gremlins.
          write(comment,'(a)') 'gremlin definitions'
          write(tokens,'(2a)') '*gremlin ',GRMFIL(1:lnblnk(GRMFIL))
          call align_comment(48,tokens,comment,aligned_str)
          write(IFC,'(a)') aligned_str(1:lnblnk(aligned_str))
        endif

        if(nsset.gt.0)then
          ifc=ifilg
          call write_parameter_sets(ifc,ier)
        endif
      endif

C Write out data for each zone in the model. The nsz loop
C is per zone, the nsc loop is for each coordinate.
      if(ncomp.gt.0)then
        WRITE(IFC,'(A)',IOSTAT=IOS,ERR=13) '#'
        WRITE(IFC,'(A)',IOSTAT=IOS,ERR=13) 
     &      '# ---- building attributes ----'
        if(initvt.eq.1)then  ! user view preference
          WRITE(IFC,'(a,6f8.1,f5.1)',IOSTAT=IOS,ERR=13)'*vew ',
     &      EYEMI(1),EYEMI(2),EYEMI(3),
     &      VIEWMI(1),VIEWMI(2),VIEWMI(3),ANGI
        endif
        do 142 nsz=1,ncomp
          write(ifilg,'(4a)',IOSTAT=IOS,ERR=13)
     &      '*start_zone',tab,zname(nsz),'  # zone name'
          write(ifilg,'(4a)',IOSTAT=IOS,ERR=13)
     &      '*shape',tab,shape(nsz),'  # polygon enclosure'

C Write out any pending rotate directives (usually only used by 3rd party).
          call eclose3(rotateit(nsz,1),rotateit(nsz,2),
     &      rotateit(nsz,3),0.0,0.0,0.0,0.01,rotationknown)
          if(rotationknown)then
            continue
          else
            write(ifilg,'(a,f6.2,a,f8.3,a,f8.3,a,f8.3,a)',IOSTAT=IOS,
     &        ERR=13) '*rotate ',
     &        rotateit(nsz,1),tab,rotateit(nsz,2),tab,
     &        rotateit(nsz,3),'  # pending rotation angle X Y'
          endif

          do 143 nsc=1,NZTV(nsz)
            write(outs,'(a,3F8.4)',IOSTAT=IOS,ERR=13)
     &        '*cord ',szcoords(nsz,nsc,1),szcoords(nsz,nsc,2),
     &        szcoords(nsz,nsc,3)
            call SDELIM(outs,outsd,'C',IW)
            write(ifilg,'(2a,i3)',IOSTAT=IOS,ERR=13)
     &        outsd(1:lnblnk(outsd)),'  # X Y Z for ',nsc
 143      continue

C The origin is ignored in poly type enclosure so do not bother with *origin.

C The size is ignored in poly type enclosure so do not bother with *size.

C Number of surfaces.
          if(shape(nsz)(1:3).eq.'box')then
            write(ifilg,'(2a,i3,a)',IOSTAT=IOS,ERR=13)
     &        '*nbwalls',tab,nbwalls(nsz),'  # number of surfaces'
          elseif(shape(nsz)(1:7).eq.'extrude')then
            write(ifilg,'(2a,i3,a)',IOSTAT=IOS,ERR=13)
     &        '*nbwalls',tab,nbwalls(nsz),'  # number of surfaces'
          else
            write(ifilg,'(2a,i3,a)',IOSTAT=IOS,ERR=13)
     &        '*nbsurfs',tab,nbwalls(nsz),'  # number of surfaces'
          endif

C << consider adding more attributes and putting each attribute
C << on a separate line with a tag.
          do 144 isu=1,NZSUR(nsz)
            call OTHERINFO(nsz,isu,OTHSTR)
            icc=IZSTOCN(nsz,isu)    ! get the connection
            lnsn=lnblnk(SNAME(nsz,isu))
            lnsc=lnblnk(SMLCN(nsz,isu))
            lnso=lnblnk(SOTF(nsz,isu))

            if(metaver.eq.1)then
              write(outs,'(5a,3I4)',IOSTAT=IOS,ERR=13)
     &          '*surface ',SNAME(nsz,isu)(1:lnsn),' ',
     &          SMLCN(nsz,isu)(1:lnsc),' ',ICT(icc),IC2(icc),IE2(icc)
              call SDELIM(outs,outsd,'C',IW)
              write(ifilg,'(2a)',IOSTAT=IOS,ERR=13)
     &          outsd(1:lnblnk(outsd)),'  # surface attributes '
            elseif(metaver.eq.2)then
              write(outs,'(7a,3I4)',IOSTAT=IOS,ERR=13)
     &          '*surface ',SNAME(nsz,isu)(1:lnsn),' ',
     &          SMLCN(nsz,isu)(1:lnsc),' ',SOTF(nsz,isu)(1:lnso),' ',
     &          ICT(icc),IC2(icc),IE2(icc)
              call SDELIM(outs,outsd,'C',IW)
              write(ifilg,'(2a)',IOSTAT=IOS,ERR=13)
     &          outsd(1:lnblnk(outsd)),'  # surface attributes '
            elseif(metaver.eq.3)then

C Follow pattern in egeometry.F ~line 2727.
              lso=lnblnk(SOTHER(nsz,isu,1))
              lso2=lnblnk(SOTHER(nsz,isu,2))
              if(lso2.eq.0)then
                SOTHER(nsz,isu,2)='000'; lso2=3
              endif
              lso3=lnblnk(SOTHER(nsz,isu,3))
              if(lso3.eq.0)then
                SOTHER(nsz,isu,3)='000'; lso3=3
              endif
              lsml=lnblnk(SMLCN(nsz,isu))
              if(lsml.eq.0)then
                SMLCN(nsz,isu)='UNKNOWN'; lsml=7
              endif
              lspa=lnblnk(SPARENT(nsz,isu))
              if(lspa.gt.12) lspa=12
              lsna=lnblnk(SNAME(nsz,isu))
              lnsname=lsna   ! update surface name length common
              lsot=lnblnk(SOTF(nsz,isu))
              loth=lnblnk(OTHSTR)
              luse1=lnblnk(SUSE(nsz,isu,1))
              luse2=lnblnk(SUSE(nsz,isu,2))
              if(SVFC(nsz,isu)(1:4).eq.'    ') SVFC(nsz,isu)='UNKN'
              if(SOTF(nsz,isu)(1:4).eq.'    ') SOTF(nsz,isu)='UNKN'
              if(SPARENT(nsz,isu)(1:2).eq.'  ') SPARENT(nsz,isu)='-'
              WRITE(left,'(10a)',IOSTAT=ios,ERR=13) '*surface ',
     &          SNAME(nsz,isu)(1:lsna),' ',SVFC(nsz,isu),' ',
     &          SPARENT(nsz,isu)(1:lspa),' ',SUSE(nsz,isu,1)(1:luse1),
     &          ' ',SUSE(nsz,isu,2)(1:luse2)
              call SDELIM(left,leftd,'C',IW)  ! pack tokens before SMLCN
              WRITE(right,'(7a)',IOSTAT=ios,ERR=13)
     &          SOTF(nsz,isu)(1:lsot),' ',SOTHER(nsz,isu,1)(1:lso),' ',
     &          SOTHER(nsz,isu,2)(1:lso2),' ',SOTHER(nsz,isu,3)(1:lso3)
              call SDELIM(right,rightd,'C',IW) ! pack tokens after SMLCN
              write(ifilg,'(6a,i3,2a)',IOSTAT=IOS,ERR=13)
     &          leftd(1:lnblnk(leftd)),',',SMLCN(nsz,isu)(1:lsml),',',
     &          rightd(1:lnblnk(rightd)),'  # ',isu,' ',OTHSTR(1:loth)
            endif

            write(louts,'(a,124I4)',IOSTAT=ios,ERR=13) '*list ',
     &        isznver(nsz,isu),(iszjvn(nsz,isu,J),J=1,isznver(nsz,isu))
            call SDELIM(louts,loutsd,'C',IW)
            write(ifilg,'(2a,i2)',IOSTAT=IOS,ERR=13)
     &        loutsd(1:lnblnk(loutsd)),'  # ',isu
 144      continue

C Shading calculation directives in the form of:
C *shad_calc,all_applic,20 # list of surfs
C << may need to adapt 82i4 format statement? >>
          if(nsurfcalc(nsz).eq.0)then
            write(ifilg,'(a)',IOSTAT=IOS,ERR=13) '# shading directives'
            write(ifilg,'(a)',IOSTAT=ios,ERR=13)
     &        '*shad_calc,none  # no temporal shading requested'
          elseif(nsurfcalc(nsz).gt.0)then
            write(ifilg,'(a)',IOSTAT=IOS,ERR=13) '# shading directives'
            if(iaplic(nsz,1).eq.1)then
              write(ifilg,'(a,i3,a)',IOSTAT=ios,ERR=13)
     &        '*shad_calc,all_applicable ',nsurfcalc(nsz),
     &        ' # list of surfs'
              write(ifilg,'(82i4)',IOSTAT=ios,ERR=14)
     &        (lstsfcalc(nsz,J),J=1,nsurfcalc(nsz))
            elseif(iaplic(nsz,1).eq.0)then
              write(ifilg,'(a,i4,a)',IOSTAT=ios,ERR=13)
     &        '*shad_calc,list',nsurfcalc(nsz),
     &        ' # list of surfs'
              write(ifilg,'(82i4)',IOSTAT=ios,ERR=14)
     &        (lstsfcalc(nsz,J),J=1,nsurfcalc(nsz))
            endif
          endif

C Insolation calculation directives.
C *insol_calc,all_applic,2 # insolation sources
C << may need to adapt 82i4 format statement >>
          write(ifilg,'(a)',IOSTAT=IOS,ERR=13)
     &      '# insolation directives'
          if(nsurfinso(nsz).eq.0)then
            write(ifilg,'(a)',IOSTAT=ios,ERR=13)
     &        '*insol_calc,none  # no insolation requested'
          elseif(nsurfinso(nsz).gt.0)then
            if(iaplic(nsz,2).eq.1)then
              write(ifilg,'(a,i3,a)',IOSTAT=ios,ERR=13)
     &        '*insol_calc,all_applicable ',nsurfinso(nsz),
     &        ' # insolation sources'
              write(ifilg,'(82i4)',IOSTAT=ios,ERR=14)
     &        (isurfinso(nsz,J),J=1,nsurfinso(nsz))
            elseif(iaplic(nsz,2).eq.0)then
              write(ifilg,'(a,i4,a)',IOSTAT=ios,ERR=13)
     &        '*insol_calc,list',nsurfinso(nsz),
     &        ' # insolation sources'
              write(ifilg,'(82i4)',IOSTAT=ios,ERR=14)
     &        (isurfinso(nsz,J),J=1,nsurfinso(nsz))
            endif
          endif

C Obstruction information lines.  Write out slightly different
C format if type 'obs' and 'obs3'
          if(nbobs(nsz).gt.0)then
            write(ifilg,'(a)',IOSTAT=IOS,ERR=13)
     &        '#  *obs = solar obstructions'
            write(ifilg,'(2a,i2,i3,a)',IOSTAT=IOS,ERR=13)
     &        '*solar_grid',tab,NOX(nsz),NOZ(nsz),
     &        ' # solar gridding density'
            do 55 ib=1,nbobs(nsz)
              if(BLOCKTYP(nsz,ib)(1:4).eq.'obs ')then
                write(outs,'(2a,7F9.3,1X,A,1X,A)',IOSTAT=ios,ERR=13)
     &            '*obs',tab,XOB(nsz,ib),YOB(nsz,ib),ZOB(nsz,ib),
     &            DXOB(nsz,ib),DYOB(nsz,ib),DZOB(nsz,ib),
     &            BANGOB(nsz,ib,1),BLOCKNAME(nsz,ib),
     &            BLOCKMAT(nsz,ib)
                call SDELIM(outs,outsd,'C',IW)
                write(ifilg,'(2a,i2)',IOSTAT=IOS,ERR=13)
     &            outsd(1:lnblnk(outsd)),'  # block ',ib
              elseif(BLOCKTYP(nsz,ib)(1:4).eq.'obs3')then
                write(outs,'(2a,9F9.3,1X,A,1X,A)',IOSTAT=ios,ERR=13)
     &            '*obs3',tab,XOB(nsz,ib),YOB(nsz,ib),ZOB(nsz,ib),
     &            DXOB(nsz,ib),DYOB(nsz,ib),DZOB(nsz,ib),
     &            BANGOB(nsz,ib,1),BANGOB(nsz,ib,2),BANGOB(nsz,ib,3),
     &            BLOCKNAME(nsz,ib),BLOCKMAT(nsz,ib)
                call SDELIM(outs,outsd,'C',IW)
                write(ifilg,'(2a,i2)',IOSTAT=IOS,ERR=13)
     &            outsd(1:lnblnk(outsd)),'  # block ',ib
              elseif(BLOCKTYP(nsz,ib)(1:4).eq.'obsp')then
                WRITE(outs,'(2a,A,1X,A)',IOSTAT=ios,ERR=13)
     &          '*obsp',' 8 6 ',BLOCKNAME(nsz,ib),BLOCKMAT(nsz,ib)
                call SDELIM(outs,outsd,'C',IW)
                write(ifilg,'(2a,i2,a)',IOSTAT=IOS,ERR=13)
     &            outsd(1:lnblnk(outsd)),'  # block ',ib,
     &            ' coords follow:'

                WRITE(outs,'(12F9.3)',IOSTAT=ios,ERR=13)
     &            XBP(nsz,ib,1),YBP(nsz,ib,1),ZBP(nsz,ib,1),
     &            XBP(nsz,ib,2),YBP(nsz,ib,2),ZBP(nsz,ib,2),
     &            XBP(nsz,ib,3),YBP(nsz,ib,3),ZBP(nsz,ib,3),
     &            XBP(nsz,ib,4),YBP(nsz,ib,4),ZBP(nsz,ib,4)
                call SDELIM(outs,outsd,'C',IW)
                write(ifilg,'(A)',IOSTAT=IOS,ERR=13)
     &            outsd(1:lnblnk(outsd))

                WRITE(outs,'(12F9.3)',IOSTAT=ios,ERR=13)
     &            XBP(nsz,ib,5),YBP(nsz,ib,5),ZBP(nsz,ib,5),
     &            XBP(nsz,ib,6),YBP(nsz,ib,6),ZBP(nsz,ib,6),
     &            XBP(nsz,ib,7),YBP(nsz,ib,7),ZBP(nsz,ib,7),
     &            XBP(nsz,ib,8),YBP(nsz,ib,8),ZBP(nsz,ib,8)
                call SDELIM(outs,outsd,'C',IW)
                write(ifilg,'(A)',IOSTAT=IOS,ERR=13)
     &            outsd(1:lnblnk(outsd))
              endif
  55        continue
          endif

C Control definitions. Scan the current control data to determine if the
C current zone has an ideal control - if so set i_ctl_link() to one.
C Next need to determine the number of periods - if one period then
C straightforware to get heating and cooling setpoints and periods. If
C multiple periods then there might be a setback. If free floating then
C this needs to be represented.

          if(icascf(nsz).gt.0) i_ctl_link(nsz)=1

C Idea control information is work in progress.
          write(ifilg,'(2a,F6.2,a)',IOSTAT=ios,ERR=13)
     &      '*heating',tab,ht_Setpoint(nsz),
     &      '  # ideal heating stepoint'
          write(ifilg,'(2a,F6.2,a)',IOSTAT=ios,ERR=13)
     &      '*cooling',tab,cl_Setpoint(nsz),
     &      '  # ideal heating stepoint'
          write(ifilg,'(2a,i3,a)',IOSTAT=ios,ERR=13)
     &      '*ideal_control',tab,i_ctl_link(nsz),
     &      '  # ideal control linkage'

C Zone operations. Include the name and then copy in line-by-line
C between 'start_literal' and 'end_literal'.
          write(tokens,'(2a)') '*opr ',LPROJ(nsz)(1:LNBLNK(LPROJ(nsz)))
          write(comment,'(a)') 'schedules'
          call align_comment(48,tokens,comment,aligned_str)
          write(ifilg,'(a)') aligned_str(1:lnblnk(aligned_str))
          write(ifilg,'(a)') '*start_literal'
          iuf=ifil+2
          call ascitometa(iuf,LPROJ(nsz),ifilg,IER)
          write(ifilg,'(a)') '*end_literal'

          if(IVF(nsz).EQ.1)then
            write(tokens,'(2a)') '*ivf ',
     &        LVIEW(nsz)(1:LNBLNK(LVIEW(nsz)))
            write(comment,'(a)') 'view factors'
            call align_comment(48,tokens,comment,aligned_str)
            write(ifilg,'(a)') aligned_str(1:lnblnk(aligned_str))
            ltf=max(1,LNBLNK(LVIEW(nsz)))
            INQUIRE (FILE=LVIEW(nsz)(1:ltf),EXIST=XST)
            if (XST) THEN
              write(ifilg,'(a)') '*start_literal'
              iuf=ifil+2
              call ascitometa(iuf,LVIEW(nsz),ifilg,IER)
              write(ifilg,'(a)') '*end_literal'
            endif
          endif
          if(IHC(nsz).EQ.1)then
            write(tokens,'(2a)') '*ihc ',
     &        LHCCO(nsz)(1:lnblnk(LHCCO(nsz)))
            write(comment,'(a)')
     &        'convective heat transfer coefficient file'
            call align_comment(48,tokens,comment,aligned_str)
            write(ifilg,'(a)') aligned_str(1:lnblnk(aligned_str))
            ltf=max(1,LNBLNK(LHCCO(nsz)))
            INQUIRE (FILE=LHCCO(nsz)(1:ltf),EXIST=XST)
            if (XST) THEN
              write(ifilg,'(a)') '*start_literal'
              iuf=ifil+2
              call ascitometa(iuf,LHCCO(nsz),ifilg,IER)
              write(ifilg,'(a)') '*end_literal'
            endif
          endif
 
C Mark the end of zone information.
          write(ifilg,'(a)',IOSTAT=IOS,ERR=13) '*end_zone'
 142    continue
      endif

C Mark the end of model information.
      write(ifilg,'(a)',IOSTAT=IOS,ERR=13) '*end'

C Now close meta file.
  99  CALL ERPFREE(ifilg,ios)
      RETURN

C Error messages.
   13 if(IOS.eq.2)then
        CALL USRMSG('No permission to write ',LFILE,'W')
      else
        CALL USRMSG('File write error in ',LFILE,'W')
      endif
      IER=1
      GOTO 99
   14 if(IOS.eq.2)then
        CALL USRMSG('No prmission to write array in ',LFILE,'W')
      else
        CALL USRMSG('Long arrary write error in ',LFILE,'W')
      endif
      IER=1
      GOTO 99
    2 CALL USRMSG(' ','Problem writing database names!','W')
      ier=2             ! Return error code,
      goto 99

      END

C ******************** ascitometa ********************
C Read of an ascii model file into current META file.

      SUBROUTINE ascitometa(iunit,LFIL,ifilg,IER)
      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      character*(*) LFIL
      character outstr*360
      logical :: XST
      logical :: unixok
      integer ltf   ! position of last character in the string.

C Transcribe an ASCII file, if file name begins with "/" open
C it directly.
      call isunix(unixok)
      if(LFIL(1:1).eq.'/')then
        ltf=max(1,LNBLNK(LFIL))
        INQUIRE (FILE=LFIL(1:ltf),EXIST=XST)
        IF (XST) THEN
          LN=max(1,lnblnk(LFIL))
          CALL EFOPSEQ(IUNIT,LFIL,1,IER)
        else
          call usrmsg('ascitometa: problem opening file!',LFIL,'W')
          return
        endif
      else
        if(.NOT.unixok)then   ! For W10 use FPOPEN
          CALL ERPFREE(IUNIT,ISTAT)
          call FPOPEN(IUNIT,ISTAT,1,0,LFIL)
          if(ISTAT.ge.0)XST=.true.
          if(ISTAT.eq.-301)XST=.false.
        else
          CALL EFOPSEQ(IUNIT,LFIL,1,IER)
        endif
      endif
      IF(IER.LT.0)THEN
        call usrmsg(' Problem detected while trying to open',LFIL,'W')
        return
      ENDIF

C Read ASCII file.
    7 READ(IUNIT,'(a)',IOSTAT=IOS,END=102) OUTSTR
      write(ifilg,'(a)') OUTSTR(1:lnblnk(outstr))
      goto 7

C End of file encountered, close the file and return to menu.
  102 if(IOS.eq.2)then
        CALL USRMSG(' ','ascitometa: end of file problem!','-')
      else
C        CALL USRMSG(' ','ascitometa: end of file reached!','-')
      endif
      CALL ERPFREE(IUNIT,ISTAT)
      return
      end


C ******************** metatoasci ********************
C Write an ascii model file from current META file literal.

      SUBROUTINE metatoasci(ifilg,iunit,LFIL,IER)
      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      character*(*) LFIL
      character outstr*360
      logical :: XST
      logical :: unixok
      integer ltf   ! position of last character in the string.

C Transcribe literal section of META file into ASCII file.
      call isunix(unixok)
      if(LFIL(1:1).eq.'/')then
        ltf=max(1,LNBLNK(LFIL))
        INQUIRE (FILE=LFIL(1:ltf),EXIST=XST)
        IF (XST) THEN
          LN=max(1,lnblnk(LFIL))
          CALL EFOPSEQ(IUNIT,LFIL,4,IER)
        else
          call usrmsg('metatoasci: problem opening file!',LFIL,'W')
          return
        endif
      else
        if(.NOT.unixok)then   ! For W10 use FPOPEN
          CALL ERPFREE(IUNIT,ISTAT)
          call FPOPEN(IUNIT,ISTAT,4,0,LFIL)
          if(ISTAT.ge.0)XST=.true.
          if(ISTAT.eq.-301)XST=.false.
        else
          CALL EFOPSEQ(IUNIT,LFIL,4,IER)
        endif
      endif
      IF(IER.LT.0)THEN
        call usrmsg(' Problem detected while trying to open',LFIL,'W')
        return
      ENDIF

C Read line from META file until *end_literal.
    7 READ(ifilg,'(a)') OUTSTR
      lno=lnblnk(outstr)
      if(OUTSTR(1:12).eq.'*end_literal')then
        CALL ERPFREE(IUNIT,ISTAT)
        return
      else
        write(iunit,'(a)') OUTSTR(1:lnblnk(outstr))
        goto 7
      endif

      return
      end

C *************** scan_embedded_database ***********
C Scan database file names and listerial blocks in META file.

      subroutine scan_embedded_database(IUC,iuf,literalpath,ier)
#include "building.h"
#include "model.h"
#include "esprdbfile.h"
#include "MultiYear_simulations.h"
#include "sbem.h"

      character loutstr*248,WORD*20
      character lworking*144
      character text*72,literalfile*72,pth*148,literalpath*72
      logical unixok

C Since some file names could be long strings, LSTRIPC and LOUTSTR are used.
C The trigger tag might be '* DATAB' or '*datab' or '*mat' or '*stdmat'
      literalfile=' '
      call isunix(unixok)
      backspace(IUC)
  44  CALL LSTRIPC(IUC,LOUTSTR,99,ND,1,'db names',IER)
      K=0
      CALL EGETW(LOUTSTR,K,WORD,'W','db names',IFLAG)
      if(LOUTSTR(1:7).eq.'* DATAB'.or.LOUTSTR(1:6).eq.'*datab')then
        CALL LSTRIPC(IUC,LOUTSTR,99,ND,1,'db names',IER)
        K=0
        CALL EGETW(LOUTSTR,K,WORD,'W','db names',IFLAG)
      elseif(LOUTSTR(1:4).eq.'*mat'.or.LOUTSTR(1:7).eq.'*stdmat')then
        continue
      endif

C Deal with each of the possible database names until a tag does not match.
      if(WORD(1:14).eq.'*start_literal')then
        lnl=lnblnk(literalfile)
        if(lnl.gt.1)then
          write(pth,'(2a)') literalpath(1:lnblnk(literalpath)),
     &      literalfile(1:lnblnk(literalfile))
          write(6,*) pth(1:lnblnk(pth))
          call metatoasci(iuc,iuf,pth,IER)
          literalfile=' '
        endif
      elseif(WORD(1:4).eq.'*prm')then
        CALL EGETRM(LOUTSTR,K,lworking,'W','material db',IER)

C Process the file name following *prm and depending on its path
C set the LFMAT common block string via call to findwhichdbpath.
        call findwhichdbpath('mat',lworking,ier)

      elseif(WORD(1:4).eq.'*mat')then
        CALL EGETRM(LOUTSTR,K,lworking,'W','material db',IER)

C Process the file name following *mat and depending on its path
C set the LFMAT common block string via call to findwhichdbpath.
        call findwhichdbpath('mat',lworking,ier)
        literalfile=lworking(1:72)

      elseif(WORD(1:7).eq.'*stdmat')then
        CALL EGETRM(LOUTSTR,K,LFMAT,'W','material db',IER)
        ipathmat=2  ! standard folder for material database

      elseif(WORD(1:4).eq.'*mlc')then
        CALL EGETRM(LOUTSTR,K,lworking,'W','multilayer db',IER)

C Process the file name following *mlc and depending on its path
C set the LFMUL common block string via call to findwhichdbpath.
        call findwhichdbpath('mul',lworking,ier)
        literalfile=lworking(1:72)
        goto 44

      elseif(WORD(1:7).eq.'*stdmlc')then
        CALL EGETRM(LOUTSTR,K,LFMUL,'W','multilayer db',IER)
        ipathmul=2  ! standard folder for MLC database

      elseif(WORD(1:4).eq.'*opt')then
        CALL EGETRM(LOUTSTR,K,lworking,'W','optical db',IER)

C Process the file name following *opt and depending on its path
C set the LOPTDB common block string via call to findwhichdbpath.
        call findwhichdbpath('opt',lworking,ier)
        literalfile=lworking(1:72)
        goto 44

      elseif(WORD(1:7).eq.'*stdopt')then
        ipathoptdb=2  ! standard folder for optical database
        CALL EGETRM(LOUTSTR,K,LOPTDB,'W','optical db',IER)

      elseif(WORD(1:4).eq.'*prs')then
        CALL EGETRM(LOUTSTR,K,lworking,'W','pressure db',IER)

C Process the file name following *prs and depending on its path
C set the lapres common block string via call to findwhichdbpath.
        call findwhichdbpath('prs',lworking,ier)
      elseif(WORD(1:7).eq.'*stdprs')then
        ipathapres=2  ! standard folder for pressure coef database
        CALL EGETRM(LOUTSTR,K,LAPRES,'W','pressure db',IER)

      elseif(WORD(1:6).eq.'*cfcdb')then
        CALL EGETRM(LOUTSTR,K,lworking,'W','CFClayers db',IER)

C Process the file name following *cfcdb and depending on its path
C set the LCFCDB common block string via call to findwhichdbpath.
        call findwhichdbpath('cfc',lworking,ier)
      elseif(WORD(1:9).eq.'*stdcfcdb')then
        CALL EGETRM(LOUTSTR,K,LCFCDB,'W','CFClayers db',IER)
        ipathcfc=2  ! standard folder for material database

      elseif(WORD(1:4).eq.'*evn')then
        CALL EGETRM(LOUTSTR,K,lworking,'W','profiles db',IER)

C Process the file name following *evn and depending on its path
C set the LPRFDB common block string via call to findwhichdbpath.
        call findwhichdbpath('evn',lworking,ier)
      elseif(WORD(1:7).eq.'*stdevn')then
        CALL EGETRM(LOUTSTR,K,LPRFDB,'W','profile db',IER)
        ipathprodb=2  ! standard folder for profiles database

      elseif(WORD(1:6).eq.'*mould')then
        CALL EGETRM(LOUTSTR,K,lworking,'W','mould isopleths',IER)

C Process the file name following *mld and depending on its path
C set the lfmould common block string via call to findwhichdbpath.
        call findwhichdbpath('mld',lworking,ier)
      elseif(WORD(1:9).eq.'*stdmould')then
        ipathmould=2  ! standard folder for mould isopleths
        CALL EGETRM(LOUTSTR,K,lfmould,'W','mould isopleths',IER)
      elseif(WORD(1:4).eq.'*clm')then
        CALL EGETRM(LOUTSTR,K,lworking,'W','weather file',IER)

C Read name of standard weather file and set a flag indicating
C single-year weather.
        call findwhichdbpath('clm',lworking,ier)
        bSY_climate_defined = .true.        
      elseif(WORD(1:7).eq.'*stdclm')then
        ipathclim=2  ! standard folder for weather
        CALL EGETRM(LOUTSTR,K,lclim,'W','weather file',IER)
        bSY_climate_defined = .true.
        
      elseif(WORD(1:6).eq.'*myclm') then

C Read multi-year weather file and set flag indicating muli-year
C weather.     
        CALL EGETRM(LOUTSTR,K,cTemp,
     &    'W','Multi-year climate db',IER)
        cMY_climate_db_name = cTemp(1:lnblnk(cTemp))
        bMY_climates_defined = .true.

C Special materials and miscellaneous component database.
      elseif(WORD(1:7).eq.'*mscldb')then
        CALL EGETRM(LOUTSTR,K,lworking,'W','misc comp db',IER)

C Process the file name following *mscldb and depending on its path
C set the MCMPDBFL common block string.
        call findwhichdbpath('msc',lworking,ier)

C Standard folder for miscellaneous component database.
      elseif(WORD(1:10).eq.'*stdmscldb')then
        CALL EGETRM(LOUTSTR,K,MCMPDBFL,'W','misc comp db',IER)
        ipathmsc=2  

C Process the file name following *pdb and depending on its path
C set the LPCDB common block string.
      elseif(WORD(1:4).eq.'*pdb')then
        CALL EGETRM(LOUTSTR,K,lworking,'W','plant comp db',IER)
        call findwhichdbpath('pdb',lworking,ier)

      elseif(WORD(1:7).eq.'*stdpdb')then
        ipathpcdb=2  ! standard folder for plant template database
        CALL EGETRM(LOUTSTR,K,LPCDB,'W','plant template db',IER)
      elseif(WORD(1:7).eq.'*predef')then
        CALL EGETRM(LOUTSTR,K,lworking,'W','predefined obj db',IER)

C Process the file name following *predef and depending on its path
C set the LPREDEF common block string.
        call findwhichdbpath('pre',lworking,ier)
        literalfile=lworking(1:72)
        goto 44
      elseif(WORD(1:10).eq.'*stdpredef')then
        ipathpredef=2  ! standard folder for predefined objects database
        CALL EGETRM(LOUTSTR,K,LPREDEF,'W','predefined db',IER)
      elseif(WORD(1:5).eq.'*sbem')then

C Process the file name following *sbem and depending on its path
C set the LSBEM common block string.
        CALL EGETRM(LOUTSTR,K,lworking,'W','SBEM db',IER)
        call findwhichdbpath('sbm',lworking,ier)

      elseif(WORD(1:8).eq.'*stdsbem')then

C SBEM database file in standard location.
        ipathsbem=2
        CALL EGETRM(LOUTSTR,K,LSBEM,'W','standard SBEM db',IER)

      elseif(WORD(1:13).eq.'*database_end')then
        return
      else                                      ! not a known tag backspace and return.
         backspace(IUC)
         return
      endif
      goto 44                                    ! scan another line.
      end

