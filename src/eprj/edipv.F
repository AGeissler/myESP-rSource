C This file is part of the ESP-r system.
C Copyright Energy Systems Research Unit, University of
C Strathclyde, Glasgow Scotland, 2001.

C ESP-r is free software.  You can redistribute it and/or
C modify it under the terms of the GNU General Public
C License as published by the Free Software Foundation 
C (version 2 orlater).

C ESP-r is distributed in the hope that it will be useful
C but WITHOUT ANY WARRANTY; without even the implied
C warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
C PURPOSE. See the GNU General Public License for more
C details.


C Subroutine included:
C IPVDAT: define information used in generating IPV.
C ipvinitialq: auto-setup a set of classic IPV topics.
C ipvseasons: presents a list of typical seasons/assessments and
C             returns a string which sumarizes the choice.
C ipvdatinit: initializes IPV data structures based on the passed act
C IPV2SIMPAR: copy relevant data from IPV description to
C             simulation parameter sets.
C getmultip: looks for typical week in each season based on closest
C            degree days and radiation patterns.


C ************* IPVDAT
C Define information used in generating IPV.
C If act = 'i' then initialise variables.
      subroutine ipvdat(act)
#include "building.h"
#include "model.h"
#include "geometry.h"
#include "esprdbfile.h"
#include "espriou.h"  
#include "net_flow.h"
#include "net_flow_data.h"
#include "control.h"
#include "seasons.h"
#include "ipvdata.h"
#include "schedule.h"
#include "help.h"
      
      integer lnblnk  ! function definition

      COMMON/FILEP/IFIL
      common/OUTIN/IUOUT,IUIN,IEOUT
      COMMON/SPAD/MMOD,LIMIT,LIMTTY
      
      integer ncomp,ncon
      COMMON/C1/NCOMP,NCON
      COMMON/SET1/IYEAR,IBDOY,IEDOY,IFDAY,IFTIME

C IPV description.
      common/IPVF/lipvdatf
      CHARACTER path*72                ! model path
      COMMON/rpath/path

C Simulation parameter sets.
      common/spfldat/nsset,isset,isstup,isbnstep,ispnstep,issave,isavgh
      INTEGER :: nsset,isset,isstup,isbnstep,ispnstep,issave,isavgh
      common/cctlnm/ctldoc,lctlf
      character CTLDOC*248,LCTLF*72

      common/caleni/nbdaytype,nbcaldays(MDTY),icalender(365)
      INTEGER nbdaytype,nbcaldays,icalender
      common/dynamico/isdynamicocup(MCOM)
      integer isdynamicocup
      integer icascf
      common/cctl/icascf(mcom)

      COMMON/CONTM0/NCONTM,NOCNTM,CONTMNAM(MCONTM)
      CHARACTER CONTMNAM*12

C Primary energy conversions and emissions.
      common/PCONV/ipconv,pcnvht,pcnvcl,pcnvlt,pcnvfn,pcnvsp,pcnvhw
      INTEGER :: ipconv
      REAL :: pcnvht,pcnvcl,pcnvlt,pcnvfn,pcnvsp,pcnvhw
      
      integer ipvictmindex    ! which contaminate source
      character ipvctmname*12 ! name of contaminate
      character ipvndname*12  ! associated flow node name
      common/contm10/ipvictmindex,ipvctmname,ipvndname

      COMMON/AFN/IAIRN,LAPROB,ICAAS(MCOM)
      INTEGER :: iairn,icaas
      CHARACTER LAPROB*72

      real vcp  ! Guth visual comfort probablility -60 - 60 values from Radiance.
      common/guth/vcp(13)

C ITM is main IPV menu.
C ITMGET is menu for selecting metric types.
C ITME is menu for energy deman sets.
C ITMM is menu for metric sets.
      DIMENSION ITMGET(17),ITM(31),ITMM(13),ITME(18),IVALS(MCOM)
      character ITM*48,ITMM*41,ITME*33,ITMGET*33,hold*40,hold5*50
      character act*1
      character lipvdatf*72,outs*124
      character DS*7,DS1*10,DS2*8
      character*72 LTMP
      character doit*248           ! temporary hvac folder path if folder does not exist
      character t40*40,t248*248,summary*24
      character descr*7,descrst*10,descrfn*10,descr2st*8,descr2fn*8
      character t12*12,key*1
      character simact*6,ipvaction*3,zdescr*124
      CHARACTER VERTC(MCONTM)*12
      character vert(MNOD)*12
      logical MODSIT,XST,OK,have_occup
      integer nitms,INO,nitmms,INOM,nitmget,INOG,nitmes,INOE ! max items and current menu item

C Variables for handling climate file.
      character llclmdb*144
      integer lndbp
      logical unixok
      character fs*1
      character TempString*148     ! temporary msc folder path to check if folder exists
      character iformat*4,ifocus*4
      character sfile*96,snpfile*96
      character limg*72
      character vcptext*96
      logical lexist               ! true if msc folder exists
      integer ISTRW
      DIMENSION ICtmPK(MCONTM)
      DIMENSION INPK(MNOD)

      helpinsub='edipv'  ! set for subroutine

C Set folder separator (fs) to \ or / as required.
      call isunix(unixok)
      if(unixok)then
        fs = char(47)
      else
        fs = char(92)
      endif

C If act = 'i' then initialise variables. Initial assumption
C is for an annual assessment.
      t40=' '
      lr=lnblnk(cfgroot)
      simact='------'  ! to signal initial state of choice

C Rescan the `climatelist` file. Check if this climate is
C in the list. If not instantiate season and typical start and end dates.
      if(ihaveseason.gt.0)then
        continue
      else

C Setup string buffer with distribution weather folder name.
        lndbp=lnblnk(standardclmpath)
        if(ipathclim.eq.0.or.ipathclim.eq.1)then
          llclmdb=LCLIM
        elseif(ipathclim.eq.2)then
          write(llclmdb,'(3a)') standardclmpath(1:lndbp),fs,
     &      LCLIM(1:lnblnk(LCLIM))
        endif
        INQUIRE (FILE=cdblfil,EXIST=XST)
        if(XST)then
          IUF=IFIL+2
          call scancdblist(IUF,llclmdb,'p',ok,ier)
          if(ok)then
            continue
          else

C Set default early winter, spring, summer, autumn, late winter periods.
            CALL EDAY(9,1,ia1wins); CALL EDAY(15,1,ia1winf)
            CALL EDAY(6,3,ia1sprs); CALL EDAY(12,3,ia1sprf)
            CALL EDAY(11,7,iasums); CALL EDAY(17,7,iasumf)
            CALL EDAY(2,10,ia2sprs); CALL EDAY(8,10,ia2sprf)
            CALL EDAY(20,11,ia2wins); CALL EDAY(26,11,ia2winf)

C Default season definitions.
            CALL EDAY(1,1,is1wins); CALL EDAY(28,2,is1winf)
            CALL EDAY(1,11,is2wins); CALL EDAY(31,12,is2winf)
            CALL EDAY(1,3,is1sprs); CALL EDAY(30,4,is1sprf)
            CALL EDAY(1,9,is2sprs); CALL EDAY(31,10,is2sprf)
            CALL EDAY(1,5,is1sums); CALL EDAY(31,8,is1sumf)
          endif
        endif
      endif

C If creating and IPV from scratch clear, select season and
C initialise the IPV common blocks.
      if(act(1:1).eq.'i'.or.act(1:1).eq.'I')then
        call clearipvdat(act)
        call ipvseasons(simact)
        call ipvdatinit(simact)
        ipvsimu=simact   ! remember which one was selected
        MODSIT=.true.
      else
        MODSIT=.false.
      endif

C Scan any demands file.
      IUO=IFIL+1
      XST=.FALSE.
      call FINDFIL(bdmds,XST)
      IF(XST)THEN
        CALL ERPFREE(IUO,ISTAT)
        CALL ERBDMD(0,ITRU,IUO,IER)
        CALL ERPFREE(IUO,ISTAT)
      ENDIF

    3 INO=-4
      IIER=0

      WRITE(ITM(1),'(2A)')     'a title: ',ipvtitl(1:36)
      WRITE(ITM(2),'(2A)')     'b synopsis: ',ipvsynop(1:34)
      if(ipvform.eq.0)then
        WRITE(ITM(3),'(a)')    'c report format >> unspecified'
      elseif(ipvform.eq.1)then
        WRITE(ITM(3),'(a)')    'c report format >> normal'
      elseif(ipvform.eq.2)then
        WRITE(ITM(3),'(a)')    'c report format >> csv   '
      endif
      if(lipvimg(1).eq.'UNKNOWN')then
        WRITE(ITM(4),'(A)')     'd model image: -'
      else
        WRITE(ITM(4),'(2A)')     'd model image: ',lipvimg(1)(1:32)
      endif
      if(lipvimg(2).eq.'UNKNOWN')then
        WRITE(ITM(5),'(A)')     'e glare image: -'
      else
        WRITE(ITM(5),'(2A)')     'e glare image: ',lipvimg(2)(1:32)
      endif
      if(ivpfocname(1:2).eq.'  '.or.ivpfocname(1:4).eq.'UNKN')then
        WRITE(ITM(6),'(a)')      'f focus zone (-)'
      else
        lnw=lnblnk(ivpfocname)  ! Re-establish the index of the focus zone.
        WRITE(ITM(6),'(a,1x,2a)') 'f focus zone (',
     &    ivpfocname(1:lnw),') '
        if(ncomp.gt.0)then
          do iz=1,ncomp
            lnzz=lnblnk(zname(iz))  ! Find matching zone index.
            if(ivpfocname(1:lnw).eq.zname(iz)(1:lnzz))then
              ipvfoczone=iz
              cycle
            endif
          enddo
        endif
      endif
      WRITE(ITM(7),'(a,i2,a)') 'g performance metrics (',nms,') '
      WRITE(ITM(8),'(a,i2,a)') 'h demand sets & report options (',
     &  neds,') '

C Echo whether all days or typical week etc.
      write(summary,'(a)') ipvsimu
      if(ipvsimu(1:6).eq.'icwint'.or.ipvsimu(1:6).eq.'icsprt'.or.
     &   ipvsimu(1:6).eq.'icsumt'.or.ipvsimu(1:6).eq.'icautt')then
        write(summary,'(a)') 'typical week'
      elseif(ipvsimu(1:6).eq.'icwins'.or.ipvsimu(1:6).eq.'icsprs'.or.
     &   ipvsimu(1:6).eq.'icsums'.or.ipvsimu(1:6).eq.'icauts')then
        write(summary,'(a)') 'all days'
      elseif(ipvsimu(1:6).eq.'icwinf'.or.ipvsimu(1:6).eq.'icsprf'.or.
     &   ipvsimu(1:6).eq.'icsumf'.or.ipvsimu(1:6).eq.'icautf')then
        write(summary,'(a)') 'typical fortnight'
      elseif(ipvsimu(1:3).eq.'ias')then    
        write(summary,'(a)') 'full year'
      elseif(ipvsimu(1:3).eq.'i3s')then    
        write(summary,'(a)') 'full seasons'
      elseif(ipvsimu(1:3).eq.'i3t')then    
        write(summary,'(a)') 'typical weeks'
      elseif(ipvsimu(1:3).eq.'i3f')then    
        write(summary,'(a)') 'typical fortnights'
      elseif(ipvsimu(1:3).eq.'i5s')then    
        write(summary,'(a)') 'full seasons'
      elseif(ipvsimu(1:3).eq.'i5t')then    
        write(summary,'(a)') 'typical weeks'
      elseif(ipvsimu(1:3).eq.'i5f')then    
        write(summary,'(a)') 'typical fortnights'
      endif
      lnsum=lnblnk(summary)
      WRITE(ITM(9),'(a,i2,3a)') 'i ___simulations (',nipvassmt,' ',
     &  summary(1:lnsum),') ___days'
      
      if(nipvassmt.eq.0)then
        m=10
        WRITE(ITM(10),'(a)')        'j  simulation not yet defined '
      elseif(nipvassmt.eq.1)then
        m=9
        call stdate(iyear,ipvastjd(1),descr,descrst,descr2st)
        call stdate(iyear,ipvafnjd(1),descr,descrfn,descr2fn)
        jjd1=(ipvafnjd(1)-ipvastjd(1))+1
        WRITE(ITM(10),'(6a,i4)') 'j ',ipvadesc(1)(1:15),' ',descrst,
     &    ' ',descrfn,jjd1
      elseif(nipvassmt.eq.3)then
        m=12
        call stdate(iyear,ipvastjd(1),descr,descrst,descr2st)
        call stdate(iyear,ipvafnjd(1),descr,descrfn,descr2fn)
        jjd1=(ipvafnjd(1)-ipvastjd(1))+1
        WRITE(ITM(10),'(6a,i4)') 'j ',ipvadesc(1)(1:15),' ',descrst,
     &    ' ',descrfn,jjd1
        call stdate(iyear,ipvastjd(2),descr,descrst,descr2st)
        call stdate(iyear,ipvafnjd(2),descr,descrfn,descr2fn)
        jjd2=(ipvafnjd(2)-ipvastjd(2))+1
        WRITE(ITM(11),'(6a,i4)') 'k ',ipvadesc(2)(1:15),' ',descrst,
     &    ' ',descrfn,jjd2
        call stdate(iyear,ipvastjd(3),descr,descrst,descr2st)
        call stdate(iyear,ipvafnjd(3),descr,descrfn,descr2fn)
        jjd3=(ipvafnjd(3)-ipvastjd(3))+1
        WRITE(ITM(12),'(6a,i4)') 'l ',ipvadesc(3)(1:15),' ',descrst,
     &    ' ',descrfn,jjd3
      elseif(nipvassmt.eq.5)then
        m=14
        call stdate(iyear,ipvastjd(1),descr,descrst,descr2st)
        call stdate(iyear,ipvafnjd(1),descr,descrfn,descr2fn)
        jjd1=(ipvafnjd(1)-ipvastjd(1))+1
        WRITE(ITM(10),'(6a,i4)') 'j ',ipvadesc(1)(1:15),' ',descrst,
     &    ' ',descrfn,jjd1
        call stdate(iyear,ipvastjd(2),descr,descrst,descr2st)
        call stdate(iyear,ipvafnjd(2),descr,descrfn,descr2fn)
        jjd2=(ipvafnjd(2)-ipvastjd(2))+1
        WRITE(ITM(11),'(6a,i4)') 'k ',ipvadesc(2)(1:15),' ',descrst,
     &    ' ',descrfn,jjd2
        call stdate(iyear,ipvastjd(3),descr,descrst,descr2st)
        call stdate(iyear,ipvafnjd(3),descr,descrfn,descr2fn)
        jjd3=(ipvafnjd(3)-ipvastjd(3))+1
        WRITE(ITM(12),'(6a,i4)') 'l ',ipvadesc(3)(1:15),' ',descrst,
     &    ' ',descrfn,jjd3
        call stdate(iyear,ipvastjd(4),descr,descrst,descr2st)
        call stdate(iyear,ipvafnjd(4),descr,descrfn,descr2fn)
        jjd4=(ipvafnjd(4)-ipvastjd(4))+1
        WRITE(ITM(13),'(6a,i4)') 'm ',ipvadesc(4)(1:15),' ',descrst,
     &    ' ',descrfn,jjd4
        call stdate(iyear,ipvastjd(5),descr,descrst,descr2st)
        call stdate(iyear,ipvafnjd(5),descr,descrfn,descr2fn)
        jjd5=(ipvafnjd(5)-ipvastjd(5))+1
        WRITE(ITM(14),'(6a,i4)') 'n ',ipvadesc(5)(1:15),' ',descrst,
     &    ' ',descrfn,jjd5
      endif
      if(nipvdispjd.eq.1)then
        WRITE(ITM(m+1),'(a,i3,a)') 'o display days (',ipvdispjd(1),') '
      elseif(nipvdispjd.eq.3)then
        WRITE(ITM(m+1),'(a,3i4,a)') 'o display days (',ipvdispjd(1),
     &     ipvdispjd(2),ipvdispjd(3),') '
      elseif(nipvdispjd.eq.5)then
        WRITE(ITM(m+1),'(a,5i4,a)') 'o  display days (',ipvdispjd(1),
     &     ipvdispjd(2),ipvdispjd(3),ipvdispjd(4),ipvdispjd(5),') '
      endif
      if(nipvassmt.eq.0.or.nipvassmt.eq.1)then
        ITM(m+2)=                    '  scaling: all seasons...     '
        WRITE(ITM(m+3), '(a,f7.2)')  'q  heating    ',ddmheat(1)
        WRITE(ITM(m+4), '(a,f7.2)')  'r  cooling    ',ddmcool(1)
        WRITE(ITM(m+5),'(a,f7.2)')   's  time-based ',ddmtime(1)
      elseif(nipvassmt.eq.3)then
        ITM(m+2)=      '  scaling: winter transition & summer       '
        WRITE(ITM(m+3),'(a,3f7.2)')  'q  heating    ',ddmheat(1),
     &    ddmheat(2),ddmheat(3)
        WRITE(ITM(m+4),'(a,3f7.2)')  'r  cooling    ',ddmcool(1),
     &    ddmcool(2),ddmcool(3)
        WRITE(ITM(m+5),'(a,3f7.2)')  's  time-based ',ddmtime(1),
     &    ddmtime(2),ddmtime(3)
      elseif(nipvassmt.eq.5)then
        ITM(m+2)=    '  scaling:  winter spring summer autumn winter'
        WRITE(ITM(m+3),'(a,5f7.2)')  'q heating    ',ddmheat(1),
     &    ddmheat(2),ddmheat(3),ddmheat(4),ddmheat(5)
        WRITE(ITM(m+4),'(a,5f7.2)')  'r cooling    ',ddmcool(1),
     &    ddmcool(2),ddmcool(3),ddmcool(4),ddmcool(5)
        WRITE(ITM(m+5),'(a,5f7.2)')  's time-based ',ddmtime(1),
     &    ddmtime(2),ddmtime(3),ddmtime(4),ddmtime(5)
      endif
      ITM(m+6)=    'w re-scan climate for seasons or HDD ratios '
      ITM(m+7)=   '! list current directives                   '
      ITM(m+8)=   '? help                                      '
      ITM(m+9)=   '- exit menu                                 '
      nitms=m+9

C Help text for this menu.
  4   helptopic='ipv_menu_overview'
      call gethelptext(helpinsub,helptopic,nbhelp)

      if(mmod.eq.8)then
        CALL EMENU('Integrated Performance View',ITM,nitms,INO)
      else
        CALL EMENU('IPV description',ITM,nitms,INO)
      endif

      if(INO.EQ.nitms)then

C If data updated, write IPV definition file. 
C Also update or create matched simulation parameter sets.
C Check to see if msc folder exists, if not create and set mscpth. 
        WRITE(TempString,'(4a)')'..',fs,'msc',fs
        INQUIRE(FILE=TempString,EXIST=lexist)
        if (.not.lexist) then
          write(doit,'(4a)') 'mkdir ',path(1:(lnblnk(path)-4)),'msc'
          call usrmsg('Creating folder:',doit,'-')
          call runit(doit,'-')
          call pausems(400)
          write(mscpth,'(3a)')'..',fs,'msc'
        endif
        if(MODSIT)then
          lr=lnblnk(cfgroot)
          WRITE(ltmp,'(4a)') mscpth(1:lnblnk(mscpth)),fs,
     &      cfgroot(1:lr),'.ipv'
          if(lipvdatf(1:8).eq.'internal')then
            lipvdatf=ltmp
            ipvaction='ipv'
          elseif(lipvdatf(1:4).eq.'UNKN')then
            lipvdatf=ltmp
            ipvaction='ipv'
          endif
          CALL EASKS(ltmp,'IPV description file ?',
     &      ' ',72,'xxx.ipv','IPV file',IER,nbhelp)
          if(ltmp(1:2).ne.'  '.and.ltmp(1:4).ne.'UNKN') lipvdatf=ltmp

C Write out the current IPV definitions.
          ipvaction='ipv'
          call mkipvdat(ifil+1,lipvdatf,ipvaction)

          if(nipvassmt.eq.1.or.nipvassmt.eq.3.or.nipvassmt.eq.5)then

C If one of the standard number of IPV assessments transfer the relevant
C data into the simulation parameter set data structure.
            call ipv2simpar(ipvsimu)
          endif
          CALL EMKCFG('-',IER)
        endif
        return
      elseif(INO.EQ.nitms-1)then
        helptopic='ipv_menu_overview'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL PHELPD('IPV setup',nbhelp,'-',0,0,IER)
      elseif(INO.EQ.1)then
        t40=ipvtitl
        CALL EASKS(t40,'IPV tile?',' ',40,'TITLE','IPV Tile',
     &    IER,nbhelp)
        if(t40(1:2).ne.'  ')ipvtitl=t40
        MODSIT=.true.
      elseif(INO.EQ.2)then
        call edisp(iuout,'Current synopsis...')
        t248=ipvsynop
        ISTRW=72
        CALL EASKS248(t248,'Synopsis',' ',ISTRW,
     &    'No sysnopsis provided for this model.','synopsis',
     &    IER,nbhelp)
        if(t248(1:2).ne.'  ')ipvsynop=t248
        MODSIT=.true.
      elseif(INO.EQ.3)then

C Toggle between reporting formats.
        ipvform=ipvform+1
        if(ipvform.gt.2) ipvform = 0
        MODSIT=.true.
        goto 3
      elseif(INO.EQ.4.or.INO.EQ.5)then

        if(INO.EQ.4) im=1
        if(INO.EQ.5) im=2

C Model Image associated with the IPV.
        INQUIRE (FILE=lipvimg(im),EXIST=XST)
        if(XST)then
          iformat='PNG '
          ifocus='****'
          call edimage(lipvimg(im),iformat,ifocus,iier)
        endif
        if(INO.EQ.4)then
          CALL EASKmbox(' ','Project image options:',
     &      'add from ../images folder','add from another folder',
     &      'accept current','cancel',' ',' ',' ',' ',IOK,nbhelp)
        elseif(INO.EQ.5)then
          CALL EASKmbox(' ','Glare image options:',
     &      'add from ../images folder','add from another folder',
     &      'accept current','cancel',' ',' ',' ',' ',IOK,nbhelp)
        endif
        if(iok.eq.1)then

C Ask for names of files in the ../images folder.
          sfile=' '
          snpfile=' '
          call edisp(iuout,' ')
          call browsefilelist('?','img','fil',sfile,snpfile,nlist,
     &      iier)

C Try to use information gathered from file scan.
          if(nlist.gt.0)then
            sfile=' '
            snpfile=' '
            call browsefilelist('b','img','fil',sfile,snpfile,nlist,
     &        iier)
            if(snpfile(1:2).ne.'  ')then
              write(limg,'(3a)')imgpth(1:lnblnk(imgpth)),fs,
     &          snpfile(1:lnblnk(snpfile))
            else
              write(limg,'(a)')imgpth(1:lnblnk(imgpth))
            endif
          else
            write(limg,'(a)')imgpth(1:lnblnk(imgpth))
          endif
          iformat='PNG '
          ifocus='****'
          call edimage(limg,iformat,ifocus,iier)
          if(iier.eq.2)then
            continue
          elseif(iier.eq.0)then
            if(INO.EQ.4)then
              lipvimg(1)=limg
              if(nipvimg.eq.0) nipvimg=1
            elseif(INO.EQ.5)then
              lipvimg(2)=limg
              if(nipvimg.eq.1) nipvimg=2
            endif
            MODSIT=.true.
          endif
          INO=-1
          goto 3
        elseif(iok.eq.2)then
          limg=' '
          call easks(limg,' ','Image file name?',72,'xxx.gif',
     &      'image file name',ier,12)
          iformat='PNG '
          ifocus='****'
          call edimage(limg,iformat,ifocus,iier)
          if(iier.eq.0)then
            if(INO.EQ.4)then
              lipvimg(1)=limg
              if(nipvimg.eq.0) nipvimg=1
            elseif(INO.EQ.5)then
              lipvimg(2)=limg
              if(nipvimg.eq.1) nipvimg=2
            endif
            MODSIT=.true.
          endif
        elseif(iok.eq.3)then
          INO=-1
          goto 3
        elseif(iok.eq.4)then
          INO=-1
          goto 3
        endif
      elseif(INO.EQ.6)then

C Choose focus zone.
        INPIC=1
        CALL EPICKS(INPIC,IVALS,' ',
     &    'Which zone is focus for visual metrics: ',
     &    12,NCOMP,zname,' zone list',IER,nbhelp)
        if(ivals(1).ne.0)then
          ipvfoczone=ivals(1)
          write(ivpfocname,'(a)') zname(ipvfoczone)
        endif
        MODSIT=.true.
      elseif(INO.EQ.7)then

C Deal with performance metrics.
        goto 42
      elseif(INO.EQ.8)then

C Demand sets - list current sets and then display as a menu.
        call listipvdat(iuout,'d',ier)
        goto 142
      elseif(INO.EQ.9)then

C Toggle the number of assessments. Also set initial display days
C so they fall within the assessments being caried out.
        call ipvseasons(simact)  ! Ask for which season/assessment period(s).
        call ipvdatinit(simact)
        ipvsimu=simact   ! remember which one was selected

        if(nipvassmt.eq.1)then
          m=10
          MODSIT=.true.
        elseif(nipvassmt.eq.3)then
          m=12
          MODSIT=.true.
        elseif(nipvassmt.eq.5)then
          m=14
          MODSIT=.true.
        endif

C Debug.
C        write(6,*)ipvastjd
C        write(6,*)ipvafnjd
C        write(6,*)ipvdispjd
C        write(6,*) ok
        call usrmsg('After switching the number of assessments',
     &    'you should re-establish scaling factors.','W')

      elseif(INO.EQ.10.and.
     &      (nipvassmt.eq.1.or.nipvassmt.eq.3.or.nipvassmt.eq.5))then

C Present pop-up of the current period(s) and then allow editing
C in terms of day and month which is then converted to a julian
C day for storage.
   91   hold = ' '
        CALL STDATE(IYEAR,ipvastjd(1),DS,DS1,DS2)
        write(outs,'(a,i4,2x,a)') ' Start : ',ipvastjd(1),DS1
        call edisp(iuout,outs)
        CALL STDATE(IYEAR,ipvafnjd(1),DS,DS1,DS2)
        write(outs,'(a,i4,2x,a)') ' Finish: ',ipvafnjd(1),DS1
        call edisp(iuout,outs)
        call edayr(ipvastjd(1),id1,im1)
        call edayr(ipvafnjd(1),id2,im2)
        WRITE(HOLD,'(4I4)')id1,im1,id2,im2
        if(nipvassmt.eq.1)then
          CALL EASKS(HOLD,'Simulation period: start (day & month) and',
     &      'finish (day & month):',40,' 9 1 15 1','single sim period',
     &      IIER,nbhelp)
        elseif(nipvassmt.eq.3)then
          CALL EASKS(HOLD,'Winter period: start (day & month) and',
     &      'finish (day & month):',40,' 9 1 15 1','winter sim period',
     &      IIER,nbhelp)
        elseif(nipvassmt.eq.5)then
          CALL EASKS(HOLD,'1st winter period: start (day & month) and',
     &      'finish (day & month):',40,' 9 1 15 1','winter sim period',
     &      IIER,nbhelp)
        endif
        K=0
        CALL EGETWI(HOLD,K,id1,1,31,'W','sim st day',IIER)
        CALL EGETWI(HOLD,K,im1,1,12,'W','sim st month',IIER)
        CALL EGETWI(HOLD,K,id2,1,31,'W','sim fn day',IIER)
        CALL EGETWI(HOLD,K,im2,1,12,'W','sim fn month',IIER)
        CALL EDAY(id1,im1,ipvastjd(1))
        CALL EDAY(id2,im2,ipvafnjd(1))
        if(iier.ne.0)goto 91
        MODSIT=.true.
      elseif(INO.EQ.11.and.(nipvassmt.eq.3.or.nipvassmt.eq.5))then
   92   hold = ' '
        CALL STDATE(IYEAR,ipvastjd(2),DS,DS1,DS2)
        write(outs,'(a,i4,2x,a)') ' Start : ',ipvastjd(2),DS1
        call edisp(iuout,outs)
        CALL STDATE(IYEAR,ipvafnjd(2),DS,DS1,DS2)
        write(outs,'(a,i4,2x,a)') ' Finish: ',ipvafnjd(2),DS1
        call edisp(iuout,outs)
        call edayr(ipvastjd(2),id1,im1)
        call edayr(ipvafnjd(2),id2,im2)
        WRITE(HOLD,'(4I4)')id1,im1,id2,im2
        if(nipvassmt.eq.3)then
          CALL EASKS(HOLD,'Transition period: start (day & month) and',
     &      'finish (day & month):',40,' 9 1 15 1','trns sim period',
     &      IIER,nbhelp)
        elseif(nipvassmt.eq.5)then
          CALL EASKS(HOLD,'Spring period: start (day & month) and',
     &      'finish (day & month):',40,' 9 1 15 1','spring sim period',
     &      IIER,nbhelp)
        endif
        K=0
        CALL EGETWI(HOLD,K,id1,1,31,'W','sim st day',IIER)
        CALL EGETWI(HOLD,K,im1,1,12,'W','sim st month',IIER)
        CALL EGETWI(HOLD,K,id2,1,31,'W','sim fn day',IIER)
        CALL EGETWI(HOLD,K,im2,1,12,'W','sim fn month',IIER)
        CALL EDAY(id1,im1,ipvastjd(2))
        CALL EDAY(id2,im2,ipvafnjd(2))
        if(iier.ne.0)goto 92
        MODSIT=.true.
      elseif(INO.EQ.12.and.(nipvassmt.eq.3.or.nipvassmt.eq.5))then
   93   hold = ' '
        CALL STDATE(IYEAR,ipvastjd(3),DS,DS1,DS2)
        write(outs,'(a,i4,2x,a)') ' Start : ',ipvastjd(3),DS1
        call edisp(iuout,outs)
        CALL STDATE(IYEAR,ipvafnjd(3),DS,DS1,DS2)
        write(outs,'(a,i4,2x,a)') ' Finish: ',ipvafnjd(3),DS1
        call edisp(iuout,outs)
        call edayr(ipvastjd(3),id1,im1)
        call edayr(ipvafnjd(3),id2,im2)
        WRITE(HOLD,'(4I4)')id1,im1,id2,im2
        CALL EASKS(HOLD,'Summer period: start (day & month) and',
     &    'finish (day & month):',40,' 9 1 15 1','summer sim period',
     &    IIER,nbhelp)
        K=0
        CALL EGETWI(HOLD,K,id1,1,31,'W','sim st day',IIER)
        CALL EGETWI(HOLD,K,im1,1,12,'W','sim st month',IIER)
        CALL EGETWI(HOLD,K,id2,1,31,'W','sim fn day',IIER)
        CALL EGETWI(HOLD,K,im2,1,12,'W','sim fn month',IIER)
        CALL EDAY(id1,im1,ipvastjd(3))
        CALL EDAY(id2,im2,ipvafnjd(3))
        if(iier.ne.0)goto 93
        MODSIT=.true.
      elseif(INO.EQ.13.and.nipvassmt.eq.5)then
   94   hold = ' '
        CALL STDATE(IYEAR,ipvastjd(4),DS,DS1,DS2)
        write(outs,'(a,i4,2x,a)') ' Autumn start : ',ipvastjd(4),DS1
        call edisp(iuout,outs)
        CALL STDATE(IYEAR,ipvafnjd(4),DS,DS1,DS2)
        write(outs,'(a,i4,2x,a)') ' Autumn finish: ',ipvafnjd(4),DS1
        call edisp(iuout,outs)
        call edayr(ipvastjd(4),id1,im1)
        call edayr(ipvafnjd(4),id2,im2)
        WRITE(HOLD,'(4I4)')id1,im1,id2,im2
        CALL EASKS(HOLD,'Autumn period: start (day & month) and',
     &    'finish (day & month):',40,' 9 1 15 1','autm sim period',
     &    IIER,4)
        K=0
        CALL EGETWI(HOLD,K,id1,1,31,'W','sim st day',IIER)
        CALL EGETWI(HOLD,K,im1,1,12,'W','sim st month',IIER)
        CALL EGETWI(HOLD,K,id2,1,31,'W','sim fn day',IIER)
        CALL EGETWI(HOLD,K,im2,1,12,'W','sim fn month',IIER)
        CALL EDAY(id1,im1,ipvastjd(4))
        CALL EDAY(id2,im2,ipvafnjd(4))
        if(iier.ne.0)goto 94
        MODSIT=.true.
      elseif(INO.EQ.14.and.nipvassmt.eq.5)then
   95   hold = ' '
        CALL STDATE(IYEAR,ipvastjd(5),DS,DS1,DS2)
        write(outs,'(a,i4,2x,a)')' 2nd winter start : ',ipvastjd(4),DS1
        call edisp(iuout,outs)
        CALL STDATE(IYEAR,ipvafnjd(5),DS,DS1,DS2)
        write(outs,'(a,i4,2x,a)')' 2nd winter finish: ',ipvafnjd(4),DS1
        call edisp(iuout,outs)
        call edayr(ipvastjd(5),id1,im1)
        call edayr(ipvafnjd(5),id2,im2)
        WRITE(HOLD,'(4I4)')id1,im1,id2,im2
        CALL EASKS(HOLD,'2nd winter period: start (day & month) and',
     &    'finish (day & month):',40,' 9 1 15 1','2nd win sim period',
     &    IIER,nbhelp)
        K=0
        CALL EGETWI(HOLD,K,id1,1,31,'W','sim st day',IIER)
        CALL EGETWI(HOLD,K,im1,1,12,'W','sim st month',IIER)
        CALL EGETWI(HOLD,K,id2,1,31,'W','sim fn day',IIER)
        CALL EGETWI(HOLD,K,im2,1,12,'W','sim fn month',IIER)
        CALL EDAY(id1,im1,ipvastjd(5))
        CALL EDAY(id2,im2,ipvafnjd(5))
        if(iier.ne.0)goto 95
        MODSIT=.true.
      elseif(INO.EQ.m+1)then

C Display days.
C Initial values should be defined when the number of assessments
C and the assessment periods were setup (e.g. the 2nd day of each
C assessment.
        if(nipvdispjd.gt.0)then
          do 195 ij=1,nipvdispjd
  194       hold = ' '
            CALL STDATE(IYEAR,ipvdispjd(ij),DS,DS1,DS2)
            write(outs,'(a,i4,2x,a)') ' Display day: ',ipvdispjd(ij),DS1
            call edisp(iuout,outs)
            call edayr(ipvdispjd(ij),id1,im1)
            WRITE(HOLD,'(2I4)')id1,im1
            CALL EASKS(HOLD,'Display day & month: (should be',
     &        '(within a simulation period):',40,' 9 1','disp day',
     &        IIER,nbhelp)
            K=0
            CALL EGETWI(HOLD,K,id1,1,31,'W','sim st day',IIER)
            CALL EGETWI(HOLD,K,im1,1,12,'W','sim st month',IIER)
            CALL EDAY(id1,im1,ipvdispjd(ij))
            if(iier.ne.0)goto 194
            MODSIT=.true.
  195     continue
        endif
        CALL EASKMBOX(' ','Options:','add another day','cancel',
     &    ' ',' ',' ',' ',' ',' ',IW,nbhelp)
        if(IW.eq.1)then
          nipvdispjd=nipvdispjd+1
          ij=nipvdispjd
  196     HOLD='  12   12 '
          CALL EASKS(HOLD,'Display day & month: (should be',
     &      '(within a simulation period):',40,' 9 1','disp day',
     &      IIER,nbhelp)
          K=0
          CALL EGETWI(HOLD,K,id1,1,31,'W','display day',IIER)
          CALL EGETWI(HOLD,K,im1,1,12,'W','display month',IIER)
          CALL EDAY(id1,im1,ipvdispjd(ij))
          if(iier.ne.0)goto 196
          CALL STDATE(IYEAR,ipvdispjd(ij),DS,DS1,DS2)
          write(outs,'(a,i4,2x,a)') ' New day: ',ipvdispjd(ij),DS1
          call edisp(iuout,outs)
          MODSIT=.true.
        endif
      elseif(INO.EQ.m+3)then

C Heating degree day multiplier or day ratio.
   96   hold5 = ' '
        WRITE(HOLD5,'(5F7.3)') ddmheat(1),ddmheat(2),ddmheat(3),
     &    ddmheat(4),ddmheat(5)
        helptopic='ipv_season_multip'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL EASKS(HOLD5,'Seasonal multiplier for heating:',
     &    ' ',50,' 17. 17. 17. 17.  17.','heating dd',IIER,nbhelp)
        K=0
        CALL EGETWR(HOLD5,K,ddmheat(1),1.,999.,'W','win ht',IIER)
        CALL EGETWR(HOLD5,K,ddmheat(2),1.,999.,'W','spr ht',IIER)
        CALL EGETWR(HOLD5,K,ddmheat(3),1.,999.,'W','sum ht',IIER)
        CALL EGETWR(HOLD5,K,ddmheat(4),1.,999.,'W','aut ht',IIER)
        CALL EGETWR(HOLD5,K,ddmheat(5),1.,999.,'W','win ht',IIER)
        if(iier.ne.0)goto 96
        MODSIT=.true.
      elseif(INO.EQ.m+4)then

C Cooling degree day multiplier or day ratio.
   97   hold5 = ' '
        WRITE(HOLD5,'(5F7.3)') ddmcool(1),ddmcool(2),ddmcool(3),
     &    ddmcool(4),ddmcool(5)
        helptopic='ipv_season_multip'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL EASKS(HOLD5,'Seasonal multiplier for cooling:',
     &    ' ',50,' 17. 17. 17. 17.  17.','cooling dd',IIER,nbhelp)
        K=0
        CALL EGETWR(HOLD5,K,ddmcool(1),1.,999.,'W','win cl',IIER)
        CALL EGETWR(HOLD5,K,ddmcool(2),1.,999.,'W','spr cl',IIER)
        CALL EGETWR(HOLD5,K,ddmcool(3),1.,999.,'W','sum cl',IIER)
        CALL EGETWR(HOLD5,K,ddmcool(4),1.,999.,'W','aut cl',IIER)
        CALL EGETWR(HOLD5,K,ddmcool(5),1.,999.,'W','win cl',IIER)
        if(iier.ne.0)goto 97
        MODSIT=.true.
      elseif(INO.EQ.m+5)then

C Time based ratio for fans and domestic hot water.
  101   hold5 = ' '
        WRITE(HOLD5,'(5F7.3)') ddmtime(1),ddmtime(2),ddmtime(3),
     &    ddmtime(4),ddmtime(5)
        helptopic='ipv_season_multip'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL EASKS(HOLD5,'Seasonal multipliers for time based items:',
     &    '(dispersed fans & DHW)',50,' 17. 17. 17. 17.  17.',
     &     'time-based dd',IIER,nbhelp)
        K=0
        CALL EGETWR(HOLD5,K,ddmtime(1),1.,999.,'W','win time',IIER)
        CALL EGETWR(HOLD5,K,ddmtime(2),1.,999.,'W','spr time',IIER)
        CALL EGETWR(HOLD5,K,ddmtime(3),1.,999.,'W','sum time',IIER)
        CALL EGETWR(HOLD5,K,ddmtime(4),1.,999.,'W','aut time',IIER)
        CALL EGETWR(HOLD5,K,ddmtime(5),1.,999.,'W','win time',IIER)
        if(iier.ne.0)goto 101
        MODSIT=.true.
      elseif(INO.EQ.m+6)then
        helptopic='ipv_season_scan'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL EASKMBOX(' ','Options:',
     &    'scan climatelist for seasons & periods',
     &    'scan for day/DD ratios','cancel',
     &    ' ',' ',' ',' ',' ',IW,nbhelp)
        if(IW.eq.1.or.IW.eq.2)then

C Expand climate file name in case of standard location.
          lndbp=lnblnk(standardclmpath)
          if(ipathclim.eq.0.or.ipathclim.eq.1)then
            llclmdb=LCLIM
          elseif(ipathclim.eq.2)then
            write(llclmdb,'(3a)') standardclmpath(1:lndbp),fs,
     &        LCLIM(1:lnblnk(LCLIM))
          endif

C Rescan climate to acquire ratios from the climatelist file.
          INQUIRE (FILE=cdblfil,EXIST=XST)
          if(XST)then
            IUF=IFIL+2
            call scancdblist(IUF,llclmdb,'p',ok,ier)
          endif
        endif
        if(IW.eq.1)then

C Re-scan climate seasons and periods. Also reset the display days.
          ipvastjd(1)=ia1wins; ipvafnjd(1)=ia1winf
          ipvastjd(2)=ia1sprs; ipvafnjd(2)=ia1sprf
          ipvastjd(3)=iasums; ipvafnjd(3)=iasumf
          ipvastjd(4)=ia2sprs; ipvafnjd(4)=ia2sprf
          ipvastjd(5)=ia2wins; ipvafnjd(5)=ia2winf
          if(nsset.eq.1)then
            nipvdispjd=5
            ipvdispjd(1)=ia1wins+1
            ipvdispjd(2)=ia1sprs+1
            ipvdispjd(3)=iasums+1
            ipvdispjd(4)=ia2sprs+1
            ipvdispjd(5)=ia2wins+1
          elseif(nsset.eq.3)then
            nipvdispjd=3
            ipvdispjd(1)=ia1wins+1
            ipvdispjd(2)=ia1sprs+1
            ipvdispjd(3)=iasums+1
          elseif(nsset.eq.5)then
            nipvdispjd=5
            ipvdispjd(1)=ia1wins+1
            ipvdispjd(2)=ia1sprs+1
            ipvdispjd(3)=iasums+1
            ipvdispjd(4)=ia2sprs+1
            ipvdispjd(5)=ia2wins+1
          endif
          MODSIT=.true.
        elseif(IW.eq.2)then

C Local (stripped) version of a subroutine from clmper.F
          helptopic='ipv_season_scan'
          call gethelptext(helpinsub,helptopic,nbhelp)
          CALL EASKMBOX(' ','Options:',
     &      'auto-setup of seasonal scaling',
     &      'manual setup of seasonal scaling','cancel',
     &      ' ',' ',' ',' ',' ',IW,nbhelp)
          if(IW.eq.1)then
            call getmultip('a')
          elseif(IW.eq.2)then
            call getmultip('-')
          elseif(IW.eq.3)then
            goto 4
          endif

C Loop through user number of assessments seasons and set dd & time arrays.
          do 105 ij=1,5
            ddmheat(ij)=dmheat(ij)
            ddmcool(ij)=dmcool(ij)
            ddmtime(ij)=dmtime(ij)
  105     continue
          MODSIT=.true.
        endif
      elseif(INO.EQ.m+7)then

C List current metrics and demand sets.
        call listipvdat(iuout,'a',ier)
      else
        INO=-4
        GOTO 4
      endif
      INO=-4
      GOTO 3

C Deal with metrics (and then return to appropriate point in main menue).
  42  continue
      M=1
      ITMM(1) =   '___metric___ zones area'
      if(nms.eq.0)then
        lnd=lnblnk(msdoc(1))
        if(lnd.lt.12) lnd=12
        WRITE(ITMM(2),'(2A,i6,F8.1)')'a ',msdoc(1)(1:lnd),
     &    nzmg(1),emgflr(1)
        M=2
      else
        do 102 L=1,nms
          M=M+1
          CALL EMKEY(L,KEY,IER)
          lnd=lnblnk(msdoc(L))
          if(lnd.lt.12) lnd=12
          WRITE(ITMM(M),'(3A,i6,F8.1)') key,' ',
     &      msdoc(L)(1:lnd),nzmg(L),emgflr(L)
 102    continue
      endif
      ITMM(M+1)=   '+ add/delete metric set            '
      ITMM(M+2)=   '! list metrics                     '
      ITMM(M+3)=   '? help                             '
      ITMM(M+4)=   '- exit menu                        '
      nitmms=M+4

C Help text for this menu.
 43   helptopic='ipv_menu_overview'
      call gethelptext(helpinsub,helptopic,nbhelp)

      CALL EMENU('IPV metric sets',ITMM,nitmms,INOM)
      if(INOM.GE.2.and.INOM.le.nitmms-4)then
        ijj=INOM-1
        goto 44
      elseif(INOM.EQ.nitmms-3)then

C Add delete copy metric set.
        CALL EASKMBOX('Metric set options:',' ','add','delete',
     &    'do nothing',' ',' ',' ',' ',' ',IRT,nbhelp)
        if(IRT.eq.1)then
          if(nms.lt.MIPVM)then
            nms=nms+1
            ijj=nms
            goto 44
          endif
        elseif(IRT.eq.2)then
          CALL EASKI(IFOC,'Which metric set (give index)?',
     &      '(zero to skip)',0,'F',nms,'F',1,'demand set index',
     &      IERI,nbhelp)
          if(ieri.eq.-3)then
            continue
          elseif(ifoc.eq.0)then
            continue
          else

C If last metric set then clear and decrement.
            if(ifoc.eq.nms)then
              imetget(nms)=0
              imetmsc(nms,1)=0
              imetmsc(nms,2)=0
              nzmg(nms)=0
              emgflr(nms)=0
              metrglbl(nms)=' '
              msdoc(nms)=' '
              metgroup(nms)=' '
              nms=nms-1
              MODSIT=.true.
              goto 42
            else
              do 111 is=ifoc,nms-1
                imetget(is)=imetget(is+1)
                imetmsc(is,1)=imetmsc(is+1,1)
                imetmsc(is,2)=imetmsc(is+1,2)
                nzmg(is)=nzmg(is+1)
                emgflr(is)=emgflr(is+1)
                metrglbl(is)=metrglbl(is+1)
                msdoc(is)=msdoc(is+1)
                metgroup(is)=metgroup(is+1)
                do iss=1,MCOM
                  izmg(is,iss)=izmg(is+1,iss)
                enddo
 111          continue
              nms=nms-1
              MODSIT=.true.
              goto 42
            endif
          endif
        elseif(IRT.eq.3)then
          continue
        endif
      elseif(INOM.EQ.nitmms-2)then

C List current metrics.
        call listipvdat(iuout,'m',ier)
      elseif(INOM.EQ.nitmms-1)then
        helptopic='ipv_menu_overview'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL PHELPD('IPV metrics',nbhelp,'-',0,0,IER)
      elseif(INOM.EQ.nitmms)then
        INO=-4
        GOTO 3
      else
        INOM=-3
        GOTO 43
      endif
      INOM=-4
      GOTO 42

 44   continue
 
C In esrures IGET indices are used. Currently comfort for an IPV makes
C use of resultant temperature (IGET=6).
C IGET = 70 IPV demands set (heating/cooling/lights/small_power/fans/DHW/PV)
C IGET = 71 IPV distributed demands (light (unctld)/light (ctld)/ fans
C            /pumps/lifts/small_power/DHW/PV
C IGET = 73 IPV visual comfort (various IGET,4) where 1=Guth VCP, 2=glare,
C           3=daylight factors
C IGET = 74 IPV conv+radiant & latent zone injection
C Contaminates requires access to flow results.
 777  ITMGET(1) ='a zone resultant temperature'
      ITMGET(2) ='b visual comfort (Guth)     '
      ITMGET(3) ='c visual comfort (glare)    '
      ITMGET(4) ='d daylight factors          '
      ITMGET(5) ='e contaminant               '
      ITMGET(6) ='f emissions               '
      ITMGET(7) ='g zone dry bulb temperature '
      ITMGET(8) ='h zone relative humidity (%)'
      ITMGET(9) ='i zone infiltration load    '
      ITMGET(10)='j zone ventilation load     '
      ITMGET(11)='k zone casual gains (all)   '
      ITMGET(12)='l zone solar (from outside) '
      ITMGET(13)='m zone solar (absorbed in)  '
      ITMGET(14)='n plant cnv + rad & latent  '
      ITMGET(15)='* NO CHOICE AT THIS TIME    '
      ITMGET(16)='? help                      '
      ITMGET(17)='- exit                      '
      nitmget=17

C Help text for this menu.
      helptopic='ipv_metric_overview'
      call gethelptext(helpinsub,helptopic,nbhelp)

C Remind use of the current topic.
      call edisp(iuout,' ')
      lnd=lnblnk(msdoc(ijj))
      lnz=lnblnk(metgroup(ijj))
      lng=lnblnk(metrglbl(ijj))
      write(outs,'(6a)') 'The current topic is ',msdoc(ijj)(1:lnd),
     &  ' with units ',metrglbl(ijj)(1:lng),' grouped as ',
     &  metgroup(ijj)(1:lnz)
      call edisp(iuout,outs)
      call EASKOK(outs,'Change/update topic attributes?',OK,nbhelp)

      if(ok)then
        CALL EMENU('Metric options',ITMGET,nitmget,INOG)
        if(INOG.eq.nitmget)then
          goto 43
        elseif(INOG.eq.nitmget-1)then
          helptopic='ipv_metric_overview'
          call gethelptext(helpinsub,helptopic,nbhelp)
          CALL PHELPD('metric options',nbhelp,'-',0,0,IER)
          goto 777
        elseif(INOG.eq.nitmget-2)then
          goto 43
        elseif(INOG.eq.1)then  ! Comfort proxy show rooms with occupants.
          imetget(ijj)=6
          imetmsc(ijj,1)=6
          imetmsc(ijj,2)=1
          msdoc(ijj)='comfort'
          metgroup(ijj)='ocup_zones'
          metrglbl(ijj)='Resultant T (degC)'
          call edisp(iuout,' ')
          call edisp(iuout,'The following zones have occupants:')
          icount=0
          IUF=IFIL+2
          do iz = 1,NCOMP
            CALL ERPFREE(IUF,ISTAT)
            CALL EROPER(0,iuout,IUF,iz,IER)
            have_occup=.FALSE.
            if(isdynamicocup(iz).ne.0)then  ! dynamic occupants
              have_occup=.TRUE.
            endif
            do IDTY=1,NBDAYTYPE        ! If gain 1st slot and > zero.
              if(NCAS(IDTY).gt.0)then
                do I = 1,NCAS(IDTY)
                  icur=ICGT(IDTY,I)
                  if(iabs(icur).eq.1.and.CMGS(IDTY,I).gt.1.0)then
                    have_occup=.TRUE.
                    cycle
                  endif
                enddo
              endif
            enddo
            if(have_occup)then
              icount=icount+1
              IVALS(icount)=iz
            endif
          enddo
          if(icount.gt.0)then
            call ZNALIST(icount,ivals,zdescr,length,ierr)
            call edisp(iuout,zdescr)
          endif
        elseif(INOG.eq.2)then  ! Guth
          imetget(ijj)=73
          imetmsc(ijj,1)=1
          imetmsc(ijj,2)=0
          msdoc(ijj)='Guth'
          lnfoc=lnblnk(ivpfocname)
          write(metgroup(ijj),'(a)') ivpfocname(1:lnfoc)
          metrglbl(ijj)='Guth comfort (-)'
          izmg(ijj,1)=ipvfoczone
          nzmg(ijj)=1
        elseif(INOG.eq.3)then  ! Glare
          imetget(ijj)=73
          imetmsc(ijj,1)=2
          imetmsc(ijj,2)=0
          msdoc(ijj)='glare'
          lnfoc=lnblnk(ivpfocname)
          write(metgroup(ijj),'(a)') ivpfocname(1:lnfoc)
          metrglbl(ijj)='Glare (-)'
          izmg(ijj,1)=ipvfoczone
          nzmg(ijj)=1
        elseif(INOG.eq.4)then
          imetget(ijj)=73
          imetmsc(ijj,1)=3
          imetmsc(ijj,2)=0
          msdoc(ijj)='daylightfact'
          lnfoc=lnblnk(ivpfocname)
          write(metgroup(ijj),'(a)') ivpfocname(1:lnfoc)
          metrglbl(ijj)='Daylight factors (%)'
          izmg(ijj,1)=ipvfoczone
          nzmg(ijj)=1
        elseif(INOG.eq.5)then
          if(NOCNTM.eq.0)then  ! Contamintes need setting up.
            if(IAIRN.ge.1)then
              call edisp(iuout,' ')
              call edisp(iuout,
     &          'A flow network exists so contaminates can be defined.')
            else
              call edisp(iuout,' ')
              call edisp(iuout,
     &          'A flow network needs to be created before ')
              call edisp(iuout,'contaminates can be defined.')
              goto 43
            endif
            ipvictmindex=0
            ipvctmname='UNKNOWN'
            ipvndname='UNKNOWN'
            if(IAIRN.ge.1)then
              call edisp(iuout,
     &        'An auto setup option can setup ambient CO2 levels plus')
              call edisp(iuout,'occupant based CO2 emissions.')
              call ctprob_init()
              ipvictmindex=1
              ipvctmname='occupants'
              if(ipvfoczone.gt.0)then
                if(ICAAS(ipvfoczone).gt.0)then
                  i=ICAAS(ipvfoczone)
                  write(ipvndname,'(a)') NDNAM(i)(1:12)   ! Remember the associated flow node.
                endif
              else

C Select associated fow node.
                call edisp(iuout,' ')
                call edisp(iuout,
     &            'Remember to set focus zone for contaminates.')
                do I=1,NNOD
                  VERT(I)=NDNAM(I)
                enddo
                NNPK=1
                call EPICKS(NNPK,INPK,'Select associated flow node',' ',
     &            12,NNOD,VERT,'Node select',IER,NH)
                if(INPK(i).gt.0)then
                  write(ipvndname,'(a)') NDNAM(INPK(1))(1:12)   ! Remember the associated flow node.
                endif
              endif
            endif
          else
            do I=1,NCONTM
              VERTC(I)=CONTMNAM(I)
            enddo
            NCPK=1
            call EPICKS(NCPK,ICtmPK,'Select contaminate','  ',
     &        12,NCONTM,VERTC,'Contaminant select',IER,NH)
            if (NCPK.gt.0) then
              ipvictmindex=ICtmPK(1)
              ipvctmname=(CONTMNAM(ICtmPK(1)))
            endif
            if(ICAAS(ipvfoczone).gt.0)then
              i=ICAAS(ipvfoczone)
              write(ipvndname,'(a)') NDNAM(i)(1:12)   ! Remember the associated flow node.
            endif
          endif
          imetget(ijj)=14
          imetmsc(ijj,1)=0
          imetmsc(ijj,2)=0
          msdoc(ijj)='contaminant'
          write(metgroup(ijj),'(a)') ivpfocname(1:lnfoc)
          metrglbl(nms)='contaminant (g/kg)'
        elseif(INOG.eq.6)then
          if(ipconv.eq.0)then
            call edisp(iuout,' ')
            call edisp(iuout,
     &      'Please add primary energy conversions & emissions to your')
            call edisp(iuout,'model via the following dialog.')
            call peconv('i')
          endif
          imetget(ijj)=75
          imetmsc(ijj,1)=0
          imetmsc(ijj,2)=0
          msdoc(ijj)='emissions'
          write(metgroup(ijj),'(a)') 'emission_zn'
          metrglbl(nms)='kg/m^2.a'
        elseif(INOG.eq.7)then
          imetget(ijj)=1
          imetmsc(ijj,1)=1
          imetmsc(ijj,2)=0
          metgroup(ijj)='dbt_zones'
          msdoc(ijj)='ZonedbT'
          metrglbl(ijj)='Zone db T (degC)'
        elseif(INOG.eq.8)then
          imetget(ijj)=13
          imetmsc(ijj,1)=13
          imetmsc(ijj,2)=0
          msdoc(ijj)='ZoneRH'
          metgroup(ijj)='rh_zones'
          metrglbl(ijj)='Zone rel humid (%)'
        elseif(INOG.eq.9)then
          imetget(ijj)=11
          imetmsc(ijj,1)=11
          imetmsc(ijj,2)=0
          msdoc(ijj)='Infiltration'
          metgroup(ijj)='infil_zones'
          metrglbl(ijj)='Infiltration (W)'
        elseif(INOG.eq.10)then
          imetget(ijj)=12
          imetmsc(ijj,1)=12
          imetmsc(ijj,2)=0
          msdoc(ijj)='Ventilation'
          metgroup(ijj)='vent_zones'
          metrglbl(ijj)='Ventilation (W)'
        elseif(INOG.eq.11)then
          imetget(ijj)=15
          imetmsc(ijj,1)=1
          imetmsc(ijj,2)=0
          msdoc(ijj)='TotalCasGain'
          metgroup(ijj)='cas_zones'
          metrglbl(ijj)='Total casual gn (W)'
        elseif(INOG.eq.12)then
          imetget(ijj)=38
          imetmsc(ijj,1)=1
          imetmsc(ijj,2)=0
          msdoc(ijj)='SolarFacade'
          metgroup(ijj)='sol_zones'
          metrglbl(ijj)='Solar via outside(W)'
        elseif(INOG.eq.13)then
          imetget(ijj)=40
          imetmsc(ijj,1)=1
          imetmsc(ijj,2)=0
          metgroup(ijj)='solabs_zones'
          msdoc(ijj)='Solarabsorb'
          metrglbl(ijj)='Solar absorbed (W)'
        elseif(INOG.eq.14)then
          imetget(ijj)=74
          imetmsc(ijj,1)=1
          imetmsc(ijj,2)=0
          msdoc(ijj)='Plt_C+R&L'
          metgroup(ijj)='pltcrl_zones'
          metrglbl(ijj)='Plnt cnv+rad&lat (W)'
        endif
      else
        if(imetget(ijj).eq.6)then  ! Keeping comfort remind user of likely zones.
          call edisp(iuout,' ')
          call edisp(iuout,'The following zones have occupants:')
          icount=0
          IUF=IFIL+2
          do iz = 1,NCOMP
            CALL ERPFREE(IUF,ISTAT)
            CALL EROPER(0,iuout,IUF,iz,IER)
            have_occup=.FALSE.
            if(isdynamicocup(iz).ne.0)then  ! dynamic occupants
              have_occup=.TRUE.
            endif
            do IDTY=1,NBDAYTYPE        ! If gain 1st slot and > zero.
              if(NCAS(IDTY).gt.0)then
                do I = 1,NCAS(IDTY)
                  icur=ICGT(IDTY,I)
                  if(iabs(icur).eq.1.and.CMGS(IDTY,I).gt.1.0)then
                    have_occup=.TRUE.
                    cycle
                  endif
                enddo
              endif
            enddo
            if(have_occup)then
              icount=icount+1
              IVALS(icount)=iz
            endif
          enddo
          if(icount.gt.0)then
            call ZNALIST(icount,ivals,zdescr,length,ierr)
            call edisp(iuout,zdescr)
          endif
        endif
      endif

C If one of the topics associated with the focus zone insert the
C relevant association. Otherwise confirm with user.
      if(imetget(ijj).eq.73.or.imetget(ijj).eq.14)then
        nzmg(ijj)=1
        izmg(ijj,1)=ipvfoczone
        emgflr(ijj)=ZBASEA(ipvfoczone)
        TFLA=ZBASEA(ipvfoczone)
        lnfoc=lnblnk(ivpfocname)
        write(t12,'(a)') ivpfocname(1:lnfoc)
        if(imetget(ijj).eq.73.and.imetmsc(ijj,1).eq.1)then     ! Guth VCP
          call edisp(iuout,' ')
          call edisp(iuout,'Current Guth VCP values:')
          call edisp(iuout,
     &'-60  -50   -40   -30   -20   -10   0   10  20  30  40   50   60')
          write(vcptext,'(13f6.2)') (vcp(ij),ij=1,13)
          call edisp(iuout,vcptext)
          CALL EASKS(vcptext,
     &      'Guth VCP values 10 deg increments -60 to 60',
     &      ' ',96,'Guth VCP','Guth VCP',IER,nbhelp)
          if(lnblnk(vcptext).gt.30)then
            k=0
            do iv=1,13
              CALL EGETWR(vcptext,K,vcp(iv),0.,90.,'W','vcp',IER)
            enddo
            write(6,*) 'vcp is now ',vcp
          endif
        elseif(imetget(ijj).eq.73.and.imetmsc(ijj,1).eq.2)then ! Glare
        elseif(imetget(ijj).eq.73.and.imetmsc(ijj,1).eq.3)then ! DF
        endif
      else
        INPIC=NCOMP
        CALL EPICKS(INPIC,IVALS,' ',
     &    ' Which zones to associate with this metric: ',
     &    12,NCOMP,zname,' zone list',IER,nbhelp)

        t12=metgroup(ijj)
        CALL EASKS(t12,
     &    'Unique identifier for zones associated with metric?',
     &    ' ',12,'offices','metric group name',IER,nbhelp)

C Find base area of each of the zones.
        TFLA=0.
        do 24 mz=1,INPIC
          izmg(ijj,mz)=IVALS(mz)
          TFLA=TFLA+ZBASEA(ivals(mz))
  24    continue
        CALL EASKR(TFLA,' ',' Associated floor area (m^2)? ',
     &    0.0,'F',9999.0,'W',0.,'associated floor m2',IERR,nbhelp)
        nzmg(ijj)=INPIC
        if(t12(1:2).ne.'  '.and.t12(1:4).ne.'UNKN')metgroup(ijj)=t12
        emgflr(ijj)=TFLA
      endif

C Only ask about the casual gain associated for relevant metrics.
      if(imetget(ijj).eq.11.or.imetget(ijj).eq.13.or.
     &   imetget(ijj).eq.14.or.imetget(ijj).eq.15.or.
     &   imetget(ijj).eq.38.or.imetget(ijj).eq.40.or.
     &   imetget(ijj).eq.73.or.imetget(ijj).eq.74.or.
     &   imetget(ijj).eq.75)then
        MODSIT=.true.
        goto 42
      else    ! Check which casual gain to filter on.
        icgv=imetmsc(ijj,1)
        CALL EASKI(icgv,'Casual gain type associated with occupancy',
     &    '(0 = full time): ',
     &    0,'F',3,'F',1,'casual type for comfort',IERI,nbhelp)
        if(ieri.eq.-3) then
          goto 43
        else

C Ok to instantiate the data and then jump to 42 to process it.
          imetmsc(ijj,1)=icgv
          MODSIT=.true.
          goto 42
        endif
      endif

C Deal with energy sets (and then return to appropriate point in main menue).
 142  continue
      M=1
      ITME(1) =   ' name       zones area scaling'
      if(neds.eq.0)then
        WRITE(ITME(2),'(2A,i4,F7.1,f6.1)')'a ',zedsdoc(1),nzedg(1),
     &    edgflr(1),edgsca(1)
        M=2

C Demand set - show user which zones are linked to a control law.
C Build up initial list of associated zones.
        ipossible=0
        if(lnblnk(lctlf).eq.0.or.lctlf(1:7).eq.'UNKNOWN')then
          call edisp(iuout,' ')
          call edisp(iuout,
     &      'No ideal or plant controls found - you may need to')
          call edisp(iuout,
     &      'manuall select zones which have environmental controls.')
          ipossible=-1
        endif
        if(ncf.eq.0.and.ncl.eq.0)then
          call edisp(iuout,' ')
          call edisp(iuout,
     &      'No ideal or plant controls found - you may need to')
          call edisp(iuout,
     &      'manuall select zones which have environmental controls.')
          ipossible=-1
        endif
        if(ipossible.eq.0)then  ! If not already filtered out proceed.
          call edisp(iuout,' ')
          call edisp(iuout,
     &    'The following zones have environmental controls:')
          DO LC=1,NCOMP
            IF(ICASCF(LC).NE.0)THEN
              ipossible=ipossible+1
              WRITE(OUTS,85)LC,ZNAME(LC),ICASCF(LC)
  85          FORMAT(' zone (',I2,') ',A,' << control ',I2)
              call edisp(iuout,outs)
            endif
          enddo
        endif
      else
        do 10 L=1,neds
          M=M+1
          CALL EMKEY(L,KEY,IER)
          WRITE(ITME(M),'(3A,i4,F7.1,f6.1)')key,' ',zedsdoc(L),
     &      nzedg(L),edgflr(L),edgsca(L)
  10    continue
      endif
      ITME(M+1)=   '+ add/delete demand set        '
      if(ifbhits.eq.1)then
        ITME(M+2)= 'j frequency bins report hits   '
      else
        ITME(M+2)= 'j frequency bins report %      '
      endif
      ITME(M+3)=   '! list demand sets             '
      ITME(M+4)=   '? help                         '
      ITME(M+5)=   '- exit menu                    '
      nitmes=M+5

C Help text for this menu.
 143  helptopic='ipv_demand_overview'
      call gethelptext(helpinsub,helptopic,nbhelp)

      CALL EMENU('IPV energy demand sets',ITME,nitmes,INOE)
      if(INOE.GE.2.and.INOE.le.nitmes-5)then
        iset=INOE-1
        goto 144
      elseif(INOE.EQ.nitmes-4)then

C Add delete copy demand set.
        CALL EASKMBOX('Demand set options:',' ','add','delete',
     &    'do nothing',' ',' ',' ',' ',' ',IRT,nbhelp)
        if(IRT.eq.1)then
          if(neds.lt.MIPVM)then
            neds=neds+1
            iset=neds
            goto 144  ! Jump to creation & editing.
          endif
        elseif(IRT.eq.2)then
          CALL EASKI(IFOC,'Which demand set (give index) ? ',
     &      '(zero to skip)',0,'F',neds,'F',1,'demand set index',
     &      IERI,nbhelp)
          if(ieri.eq.-3) then
            continue
          elseif(ifoc.eq.0)then
            continue
          else

C If the last demand set clear and then decrement.
            if(ifoc.eq.neds)then
              zedsdoc(neds)=' '
              nzedg(neds)=0
              edgflr(neds)=1.
              edgsca(is)=1.0
              do iss=1,MCOM
                izedg(neds,iss)=0
              enddo
              neds=neds-1
              MODSIT=.true.
              goto 142  ! Redisplay demand set menu.
            else
              do 11 is=ifoc,neds-1
                zedsdoc(is)=zedsdoc(is+1)
                nzedg(is)=nzedg(is+1)
                edgflr(is)=edgflr(is+1)
                edgsca(is)=edgsca(is+1)
                do iss=1,MCOM
                  izedg(is,iss)=izedg(is+1,iss)
                enddo
  11          continue
              neds=neds-1
              MODSIT=.true.
              goto 142  ! Redisplay demand set menu.
            endif
          endif
        elseif(IRT.eq.3)then
          continue
        endif
      elseif(INOE.EQ.nitmes-3)then
        if(ifbhits.eq.1)then
          ifbhits=0
          call edisp(iuout,'Frequency bins report %.')
          MODSIT=.true.
        else
          ifbhits=1
          call edisp(iuout,'Frequency bins report hits.')
          MODSIT=.true.
        endif
      elseif(INOE.EQ.nitmes-2)then

C List current demand sets..
        call listipvdat(iuout,'d',ier)
      elseif(INOE.EQ.nitmes-1)then
        helptopic='ipv_demand_overview'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL PHELPD('IPV demand sets',nbhelp,'-',0,0,IER)
      elseif(INOE.EQ.nitmes)then
        INO=-4
        GOTO 3
      else
        INOM=-3
        GOTO 143
      endif
      INOM=-4
      GOTO 142

C Editing or creation of a deman set.
  144 continue
      t12=zedsdoc(iset)
      CALL EASKS(t12,'Identifier for this demand set ?',
     &  ' ',12,'offices','demand set name',IER,nbhelp)
      if(t12(1:2).ne.'  '.and.t12(1:4).ne.'UNKN')zedsdoc(iset)=t12

      icount=0
      call edisp(iuout,' ')
      call edisp(iuout,
     &  'The following zones have environmental controls:')
      DO LC=1,NCOMP
        IF(ICASCF(LC).NE.0)THEN
          icount=icount+1
          IVALS(icount)=LC
        endif
      enddo
      if(icount.gt.0)then
        call ZNALIST(icount,ivals,zdescr,length,ierr)
        call edisp(iuout,zdescr)
      endif

      INPIC=NCOMP
      CALL EPICKS(INPIC,IVALS,' ',
     &  ' Which zones to include in this energy demand set: ',
     &  12,NCOMP,zname,'zone list',IER,nbhelp)
      nzedg(iset)=INPIC
      TFLA=0.
      do 26 mz=1,nzedg(iset)
        izedg(iset,mz)=IVALS(mz)
        TFLA=TFLA+ZBASEA(ivals(mz))
  26  continue
      edgflr(iset)=TFLA
      CALL EASKR(edgflr(iset),' ','Energy set floor area (m^2)? ',
     &  0.0,'F',9999.0,'W',1.,'energy set m2',IER,nbhelp)
      CALL EASKR(edgsca(iset),' ','Multiplier for energy set? ',
     &  0.0,'F',9999.0,'W',1.,'1st set mult',IER,nbhelp)
      MODSIT=.true.
      goto 142

      end

C ************* ipvinitialq
C Initial questions to populate a classic IPV. An initial check
C of data structures and warnings about dependencies. Returns
C iresponse = 0 if no issues, -1 if user chose not to proceed,
C or > 0 if dependencies missing.
C act = 'd' check dependencies and return.
      subroutine ipvinitialq(act,iresponse)
#include "building.h"
#include "model.h"
#include "geometry.h"
#include "esprdbfile.h"
#include "espriou.h"  
#include "net_flow.h"
#include "net_flow_data.h"
#include "control.h"
#include "seasons.h"
#include "ipvdata.h"
#include "schedule.h"
#include "help.h"

      COMMON/FILEP/IFIL
      common/OUTIN/IUOUT,IUIN,IEOUT

C IPV description.
      character lipvdatf*72
      common/IPVF/lipvdatf
      CHARACTER path*72                ! model path
      COMMON/rpath/path
      integer ncomp,ncon
      common/c1/ncomp,ncon
      common/cctlnm/ctldoc,lctlf
      character CTLDOC*248,LCTLF*72
      integer icascf
      common/cctl/icascf(mcom)

C Primary energy conversions and emissions.
      common/PCONV/ipconv,pcnvht,pcnvcl,pcnvlt,pcnvfn,pcnvsp,pcnvhw
      INTEGER :: ipconv
      REAL :: pcnvht,pcnvcl,pcnvlt,pcnvfn,pcnvsp,pcnvhw

      common/caleni/nbdaytype,nbcaldays(MDTY),icalender(365)
      INTEGER nbdaytype,nbcaldays,icalender
      common/dynamico/isdynamicocup(MCOM)
      integer isdynamicocup

      COMMON/CONTM0/NCONTM,NOCNTM,CONTMNAM(MCONTM)
      CHARACTER CONTMNAM*12

      COMMON/AFN/IAIRN,LAPROB,ICAAS(MCOM)
      INTEGER :: iairn,icaas
      CHARACTER LAPROB*72

      integer ipvictmindex    ! which contaminate source
      character ipvctmname*12 ! name of contaminate
      character ipvndname*12  ! associated flow node name
      common/contm10/ipvictmindex,ipvctmname,ipvndname

      real vcp  ! Guth visual comfort probablility -60 - 60 values from Radiance.
      common/guth/vcp(13)

      DIMENSION ICtmPK(MCONTM),INPK(MNOD)
      CHARACTER VERTC(MCONTM)*12
      character vert(MNOD)*12
      character outs*124
      character act*1,fs*1
      character ipvaction*3
      logical unixok,ok,have_occup,XST
      integer IVALS(MCOM)  ! the array of zones to include
      character TempString*148     ! temporary msc folder path to check if folder exists
      character doit*248           ! temporary hvac folder path if folder does not exist
      character*72 LTMP,zdescr*124
      character llclmdb*144
      logical lexist               ! true if msc folder exists

      helpinsub='edipv'  ! set for subroutine
      helptopic='ipv_initial_q'
      call gethelptext(helpinsub,helptopic,nbhelp)

C Set folder separator (fs) to \ or / as required.
      call isunix(unixok)
      if(unixok)then
        fs = char(47)
      else
        fs = char(92)
      endif
  
C Initial pass on dependencies.
      if(act(1:1).eq.'d')then
        ok=.true.
      else
        CALL EASKOK(' ',
     &  'Check for common IPV dependencies in this model? ',OK,nbhelp)
      endif
      if(OK)then
        idepend=0
        if(ipconv.eq.0)then
          call edisp(iuout,' ')
          call edisp(iuout,
     &      'Please add primary energy conversions & emissions to your')
          call edisp(iuout,'model via the following dialog.')
          call peconv('i')
        endif

        if(bdmds(1:7).eq.'UNKNOWN')then
          call edisp(iuout,' ')
          call edisp(iuout,
     &      'If you want to include non-zone specific demands such as')
          call edisp(iuout,
     &      'lifts and pumps and DHW create relevant schedules via')
          call edisp(iuout,
     &      'browse->context->fan/lift/DHW dispersed demands')
          idepend=idepend +1
        endif
        if(IAIRN.ge.1)then
          continue
        else
          call edisp(iuout,' ')
          call edisp(iuout,
     &      'If you want to include IAQ a flow network is needed.')
          call edisp(iuout,
     &      'The easiest approach is to attribute surfaces leakage')
          call edisp(iuout,
     &      'characteristics and then auto-build a flow network.')
        endif

C If there are no contaminants yet but there is a flow network and occupants
C in at least some of the zones then offer to create.
        if(NOCNTM.eq.0)then  ! Contaminants.
          call edisp(iuout,' ')
          call edisp(iuout,
     &    'If you want to include IAQ contaminate sources are needed.')
          if(IAIRN.ge.1)then
            call zones_with_occupants(icount,ivals,zdescr,ierr)
            if(icount.eq.0)then
              call usrmsg(
     &          'No occupants in building so CO2 tracking is not',
     &          'possible. Perhaps your model needs more work?','W')
              goto 42
            endif
            call edisp(iuout,
     &      'Enabling ambient CO2 levels plus occupant CO2 in model.')
            call ctprob_init()     ! Setup the contaminate regime.
            ipvictmindex=1         ! Assign that contaminate.
            ipvctmname='occupants'

C Select associated fow node.
            do I=1,NNOD
              VERT(I)=NDNAM(I)
            enddo
            NNPK=1
            call EPICKS(NNPK,INPK,'Select associated flow node',' ',
     &        12,NNOD,VERT,'Node select',IER,nbhelp)
            if(INPK(1).gt.0)then
              write(ipvndname,'(a)') NDNAM(INPK(1))(1:12)   ! Remember the associated flow node.
            endif
          endif
        else

C Select associated fow node.
          do I=1,NNOD
            VERT(I)=NDNAM(I)
          enddo
          NNPK=1
          call EPICKS(NNPK,INPK,'Select associated flow node',' ',
     &      12,NNOD,VERT,'Node select',IER,nbhelp)
          if(INPK(i).gt.0)then
            write(ipvndname,'(a)') NDNAM(INPK(1))(1:12)   ! Remember the associated flow node.
          endif

C Provide list of contaminates.
          do I=1,NCONTM
            VERTC(I)=CONTMNAM(I)
          enddo
          NCPK=1
          call EPICKS(NCPK,ICtmPK,'Select contaminate','  ',
     &      12,NCONTM,VERTC,'Contaminant select',IER,nbhelp)
          if (NCPK.gt.0) then
            ipvictmindex=ICtmPK(1)
            ipvctmname=(CONTMNAM(ICtmPK(1)))
          endif
        endif

 42     if(lradcf(1:7).eq.'UNKNOWN')then
          call edisp(iuout,' ')
          call edisp(iuout,
     &      'If you want to include visual comfort e.g. glare or')
          call edisp(iuout,
     &      'daylight factors invoke the e2r module, select on of')
          call edisp(iuout,
     &      'the scene types, select a focus zone and provide the')
          call edisp(iuout,
     &      'requested view point or gridding information. Use')
          call edisp(iuout,'browse->visualation to get started.')
          idepend=idepend +1
        endif

        if(lnblnk(lctlf).eq.0.or.lctlf(1:7).eq.'UNKNOWN')then
          call edisp(iuout,' ')
          call edisp(iuout,
     &      'If you want to track energy demands please include')
          call edisp(iuout,
     &      'environmental controls via browse->controls.')
          idepend=idepend +1
        endif

C More dependencies here...
        if(act(1:1).eq.'d')then
          iresponse=idepend
          return
        endif
      endif

      CALL EASKOK('Continue setting up a classic IPV?',
     &  ' ',OK,nbhelp)
      if(OK)then

C Instanciate the usual suspects for a classic IPV.
        call clearipvdat(act)
        lr=lnblnk(cfgroot)
        ipvsimu='i3t   '  ! remember in common
        nipvassmt=3       ! three seasons (one week in each)
        ipvastjd(1)=ia1wins
        ipvafnjd(1)=ia1winf
        ipvastjd(2)=ia1sprs
        ipvafnjd(2)=ia1sprf
        ipvastjd(3)=iasums
        ipvafnjd(3)=iasumf
        nipvdispjd=3
        ipvdispjd(1)=ia1wins+1
        ipvdispjd(2)=ia1sprs+1
        ipvdispjd(3)=iasums+1
        write(ipvadesc(1),'(A)') 'winter week'
        write(ipvadesc(2),'(A)') 'transition week'
        write(ipvadesc(3),'(A)') 'summer week'
        do ij=1,5
          ddmheat(ij)=1.0
          ddmcool(ij)=1.0
          ddmtime(ij)=1.0
          
        enddo
        write(ipvtitl,'(a)') modeltitle(1:40)
        write(ipvsynop,'(2a)') 'Standard IPV topics have been setup ',
     &    'and it is up to the user to implement dependencies.'

C Choose focus zone.
        INPIC=1
        CALL EPICKS(INPIC,IVALS,' ',
     &    'Which zone is focus for visual metrics: ',
     &    12,NCOMP,zname,' zone list',IER,nbhelp)
        if(ivals(1).ne.0)then
          ipvfoczone=ivals(1)
          write(ivpfocname,'(a)') zname(ipvfoczone)
        endif
        lnfoc=lnblnk(ivpfocname)

C Check that the msc folder exists and also set the directives
C file name.
        WRITE(TempString,'(4a)')'..',fs,'msc',fs
        INQUIRE(FILE=TempString,EXIST=lexist)
        if (.not.lexist) then
          write(doit,'(4a)') 'mkdir ',path(1:(lnblnk(path)-4)),'msc'
          call usrmsg('Creating folder:',doit,'-')
          call runit(doit,'-')
          call pausems(400)
          write(mscpth,'(3a)')'..',fs,'msc'
        endif
        lr=lnblnk(cfgroot)
        WRITE(ltmp,'(4a)') mscpth(1:lnblnk(mscpth)),fs,
     &    cfgroot(1:lr),'.ipv'
        if(lipvdatf(1:8).eq.'internal')then
          lipvdatf=ltmp
          ipvaction='ipv'
        elseif(lipvdatf(1:4).eq.'UNKN')then
          lipvdatf=ltmp
          ipvaction='ipv'
        endif
        CALL EASKS(ltmp,'IPV directives file?',
     &    '(confirm name)',72,'xxx.ipv','IPV file',IER,nbhelp)
        if(ltmp(1:2).ne.'  '.and.ltmp(1:4).ne.'UNKN') lipvdatf=ltmp

C Classic IPV has one week for each season. Establish these
C weeks by scanning the climatelist file.
C Expand climate file name in case of standard location.
        lndbp=lnblnk(standardclmpath)
        if(ipathclim.eq.0.or.ipathclim.eq.1)then
          llclmdb=LCLIM
        elseif(ipathclim.eq.2)then
          write(llclmdb,'(3a)') standardclmpath(1:lndbp),fs,
     &      LCLIM(1:lnblnk(LCLIM))
        endif

C Rescan climate to acquire ratios from the climatelist file.
        INQUIRE (FILE=cdblfil,EXIST=XST)
        if(XST)then
          IUF=IFIL+2
          call scancdblist(IUF,llclmdb,'p',ok,ier)
        endif

C Re-set climate seasons and periods. Also reset the display days.
        ipvastjd(1)=ia1wins; ipvafnjd(1)=ia1winf
        ipvastjd(2)=ia1sprs; ipvafnjd(2)=ia1sprf
        ipvastjd(3)=iasums; ipvafnjd(3)=iasumf
        ipvastjd(4)=ia2sprs; ipvafnjd(4)=ia2sprf
        ipvastjd(5)=ia2wins; ipvafnjd(5)=ia2winf
        nipvdispjd=3
        ipvdispjd(1)=ia1wins+1
        ipvdispjd(2)=ia1sprs+1
        ipvdispjd(3)=iasums+1

C Establish scaling - tell the user what is about to happen.
        call edisp(iuout,' ')
        call edisp(iuout,
     &    'An IPV assess one week in the winter, spring and summer and')
        call edisp(iuout,
     &    'scales these weeks to get to whole season performance by')
        call edisp(iuout,
     &    'scanning weather data to find best fit weeks and also uses')
        call edisp(iuout,
     &    'ratios of heating and cooling degree days between weeks and')
        call edisp(iuout,
     &    'seasons for scaling. There are a sequence of tasks with')
        call edisp(iuout,'a few confirmations needed.')

        call usrmsg(
     &  'To establish best-fit weeks and scaling for assessments',
     &  'please respond to the following dialogs...','W')
        CALL EASKMBOX(' ','Options:',
     &    'auto-setup of seasonal scaling',
     &    'manual setup of seasonal scaling',' ',
     &    ' ',' ',' ',' ',' ',IW,nbhelp)
        if(IW.eq.1)then
          call getmultip('a')
        elseif(IW.eq.2)then
          call getmultip('-')
        endif

C Loop through user number of assessments seasons and set dd arrays.
        do 105 ij=1,5
          ddmheat(ij)=dmheat(ij)
          ddmcool(ij)=dmcool(ij)
          ddmtime(ij)=dmtime(ij)
  105   continue

C Add in the usual metrics for thermal and visual comfort.
C For comfort - loop through each zone with an operation file
C loking for the existance of occupants.
        nms=0
        nms=nms+1
        imetget(nms)=6
        imetmsc(nms,1)=6
        imetmsc(nms,2)=1
        msdoc(nms)='comfort'
        metgroup(nms)='ocup_zones'
        metrglbl(nms)='Resultant T (degC)'
        IUF=IFIL+2
        icount=0
        TFLA=0.
        call edisp(iuout,' ')
        call edisp(iuout,'The following zones have occupants:')
        do iz = 1,NCOMP
          CALL ERPFREE(IUF,ISTAT)
          CALL EROPER(0,iuout,IUF,iz,IER)
          have_occup=.FALSE.
          if(isdynamicocup(iz).ne.0)then  ! dynamic occupants
            have_occup=.TRUE.
          endif
          do IDTY=1,NBDAYTYPE
            if(NCAS(IDTY).gt.0)then
              do I = 1,NCAS(IDTY)
                icur=ICGT(IDTY,I)
                if(iabs(icur).eq.1.and.CMGS(IDTY,I).gt.1.0)then
                  have_occup=.TRUE.
                  cycle
                endif
              enddo
            endif
          enddo
          if(have_occup)then
            icount=icount+1
            IVALS(icount)=iz
            nzmg(nms)= nzmg(nms)+1
            TFLA=TFLA+ZBASEA(iz)
            izmg(nms,nzmg(nms))=iz
          endif
        enddo
        if(icount.gt.0)then
          call ZNALIST(icount,ivals,zdescr,length,ierr)
          call edisp(iuout,zdescr)
        endif

        nms=nms+1
        imetget(nms)=73   ! Guth is linked to focus zone.
        imetmsc(nms,1)=1
        imetmsc(nms,2)=0
        msdoc(nms)='Guth'
        write(metgroup(nms),'(a)') ivpfocname(1:lnfoc)
        metrglbl(nms)='Guth comfort (-)'
        nzmg(nms)=1
        izmg(nms,1)=ipvfoczone
        emgflr(nms)=ZBASEA(ipvfoczone)
        do ii=1,13
          vcp(ii)=0.0  ! Reset the VCP predictions.
        enddo

        nms=nms+1
        imetget(nms)=73   ! Glare is linked to focus zone.
        imetmsc(nms,1)=2
        imetmsc(nms,2)=0
        msdoc(nms)='glare'
        write(metgroup(nms),'(a)') ivpfocname(1:lnfoc)
        metrglbl(nms)='Glare (-)'
        nzmg(nms)=1
        izmg(nms,1)=ipvfoczone
        emgflr(nms)=ZBASEA(ipvfoczone)

        nms=nms+1
        imetget(nms)=73   ! Daylight factors are linked to focus zone.
        imetmsc(nms,1)=3
        imetmsc(nms,2)=0
        msdoc(nms)='daylightfact'
        write(metgroup(nms),'(a)') ivpfocname(1:lnfoc)
        metrglbl(nms)='Daylight factors (%)'
        nzmg(nms)=1
        izmg(nms,1)=ipvfoczone
        emgflr(nms)=ZBASEA(ipvfoczone)

C Demand set & emissions - discover which zones are linked to a control law.
C Build up initial list of associated zones.
        ipossible=0
        if(lnblnk(lctlf).eq.0.or.lctlf(1:7).eq.'UNKNOWN')then
          call edisp(iuout,' ')
          call edisp(iuout,
     &      'No ideal or plant controls found - you may need to')
          call edisp(iuout,
     &      'manuall select zones which have environmental controls.')
          call edisp(iuout,
     &      'Remember to check zones associated with other topics and')
          call edisp(iuout,
     &      'also re-scan climate for seasons and DD ratios.')
          ipossible=-1
        endif
        if(ncf.eq.0.and.ncl.eq.0)then
          call edisp(iuout,' ')
          call edisp(iuout,
     &      'No ideal or plant controls found - you may need to')
          call edisp(iuout,
     &      'manuall select zones which have environmental controls.')
          call edisp(iuout,
     &      'Remember to check zones associated with other topics and')
          call edisp(iuout,
     &      'also re-scan climate for seasons and DD ratios.')
          ipossible=-1
        endif

C If there are demands then there will be emissions. Zone list same
C as demands.
        if(ipossible.eq.0)then
          if(ipconv.eq.0)then
            call edisp(iuout,' ')
            call edisp(iuout,
     &      'Please add primary energy conversions & emissions to your')
            call edisp(iuout,'model via the following dialog.')
            call peconv('i')
          endif
          icount=0
          TFLA=0.
          DO LC=1,NCOMP
            IF(ICASCF(LC).NE.0)THEN
              icount=icount+1
              nzmg(5)= nzmg(5)+1
              izmg(5,nzmg(5))=LC
              TFLA=TFLA+ZBASEA(LC)
            endif
          enddo
          if(icount.gt.0)then
            nms=nms+1
            imetget(nms)=75      ! Emissions IGET is 75 res.
            imetmsc(nms,1)=0
            imetmsc(nms,2)=0
            msdoc(nms)='emissions'
            write(metgroup(nms),'(a)') 'emission_zn'
            metrglbl(nms)='kg/m^2.a'
            emgflr(nms)=TFLA
          endif
        endif

C If there are contaminates IMFGET = 14 for this topic.
        if(NOCNTM.gt.0)then  ! Contaminants.
          nms=nms+1
          imetget(nms)=14
          imetmsc(nms,1)=0
          imetmsc(nms,2)=0
          msdoc(nms)='contaminant'
          write(metgroup(nms),'(a)') ivpfocname(1:lnfoc)
          metrglbl(nms)='contaminant (g/kg)'
          nzmg(nms)=1
          izmg(nms,1)=ipvfoczone
          emgflr(nms)=ZBASEA(ipvfoczone)
        endif

C If environmental controls then create a demand entry.
        if(ipossible.eq.0)then
          icount=0
          TFLA=0.
          call edisp(iuout,' ')
          call edisp(iuout,
     &    'The following zones have environmental controls:')
          DO LC=1,NCOMP
            IF(ICASCF(LC).NE.0)THEN
              ipossible=ipossible+1
              WRITE(OUTS,85)LC,ZNAME(LC),ICASCF(LC)
  85          FORMAT(' zone (',I2,') ',A,' << control ',I2)
              call edisp(iuout,outs)
              izedg(1,ipossible)=LC
              TFLA=TFLA+ZBASEA(LC)
              icount=icount+1
              IVALS(icount)=LC
            endif
          enddo
          if(ipossible.gt.0)then
            nzedg(1)=ipossible
            edgflr(1)=TFLA
            edgsca(1)=1.0
            zedsdoc(1) = 'zones'
            neds=1
            call ZNALIST(icount,ivals,zdescr,length,ierr)
            call edisp(iuout,zdescr)
          endif
        endif
        iresponse=0

C Save the current information to the IPV directives file before exiting.
        ipvaction='ipv'
        call mkipvdat(ifil+1,lipvdatf,ipvaction)
        call edisp(iuout,' ')
        call edisp(iuout,
     &    'NEXT: confirm the details from the initial setup e.g. the')
        call edisp(iuout,
     &    'associated zones for each topic, periods & DD ratios.')
        return
      else
        iresponse=-1
        return
      endif
      return
      end


C ************* ipvseasons
C ipvseasons presents a list of typical seasons/assessments and
C returns a string which sumarizes the choice (see list of values
C in subroutine ipvdatinit.
      subroutine ipvseasons(simact)
#include "help.h"

      DIMENSION ITM(22)
      character ITM*33
      character simact*6
      integer nitms,INO ! max items and current menu item

      helpinsub='edipv'  ! set for subroutine

      simact='ias   '   ! assume annual assessment

    3 INO=-4
      ITM(1)= 'a typical week in winter         '
      ITM(2)= 'b typical fortnight in winter    '
      ITM(3)= 'c winter season (from January)   '
      ITM(4)= 'd typical week in spring         '
      ITM(5)= 'e typical fortnight in spring    '
      ITM(6)= 'f spring season                  '
      ITM(7)= 'g typical week in summer         '
      ITM(8)= 'h typical fortnight in summer    '
      ITM(9)= 'i summer season                  '
      ITM(10)='j typical week in autumn         '
      ITM(11)='k typical fortnight in autumn    '
      ITM(12)='l autumn season                  '
      ITM(13)='m annual assessment (default)    '
      ITM(14)='n three seasons (all days)       '
      ITM(15)='o three seasons (typical week)   '
      ITM(16)='p three seasons (typ fortnight)  '
      ITM(17)='q five seasons (all days)        '
      ITM(18)='r five seasons (typical week)    '
      ITM(19)='s five seasons (typ fortnight)   '
      ITM(20)='? help                           '
      ITM(21)='- exit menu                      '
      nitms=21

      helptopic='ipv_assessment_regime'
      call gethelptext(helpinsub,helptopic,nbhelp)
      CALL EMENU('Options for assessments',ITM,nitms,INO)

      if(INO.EQ.nitms)then
        return
      elseif(INO.EQ.nitms-1)then
        helptopic='ipv_assessment_regime'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL PHELPD('IPV seasons',nbhelp,'-',0,0,IER)
        goto 3
      elseif(INO.eq.1)then  ! typical week in winter
        simact='icwint'
      elseif(INO.eq.2)then  ! typical fortnight in winter
        simact='icwinf'
      elseif(INO.eq.3)then  ! winter season
        simact='icwins'
      elseif(INO.eq.4)then  ! typical week in spring 
        simact='icsprt'
      elseif(INO.eq.5)then  ! typical fortnight in spring
        simact='icsprf'
      elseif(INO.eq.6)then  ! spring season
        simact='icsprs'
      elseif(INO.eq.7)then  ! typical week in summer
        simact='icsumt'
      elseif(INO.eq.8)then  ! typical fortnight in summer
        simact='icsumf'
      elseif(INO.eq.9)then  ! summer season
        simact='icsums'
      elseif(INO.eq.10)then ! typical week in autumn 
        simact='icautt'
      elseif(INO.eq.11)then ! typical fortnight in autumn
        simact='icautf'
      elseif(INO.eq.12)then ! autumn season
        simact='icauts'
      elseif(INO.eq.13)then ! annual assessment
        simact='ias   '
      elseif(INO.eq.14)then ! three seasons (all days)
        simact='i3s   '
      elseif(INO.eq.15)then ! three seasons (typical week)
        simact='i3t   '
      elseif(INO.eq.16)then ! three seasons (typical fortnight) 
        simact='i3f   '
      elseif(INO.eq.17)then ! five seasons (all days)
        simact='i5s   '
      elseif(INO.eq.18)then ! five seasons (typical week)
        simact='i5t   '
      elseif(INO.eq.19)then ! five seasons (typical fortnight) 
        simact='i5f   '
      endif
      return

      end

C ************* ipvdatinit
C ipvdatinit initializes IPV data structures based on the passed act
C parameter where
C   act='icwint' initial values using winter climate typical week
C      ='icwinf' initial values using winter climate fortnight
C      ='icwins' initial values using winter climate season
C      ='icsprt' initial values using spring climate typical week
C      ='icsprf' initial values using spring climate fortnight
C      ='icsprs' initial values using spring climate season
C      ='icsumt' initial values using summer climate typical week
C      ='icsumf' initial values using summer climate fortnight
C      ='icsums' initial values using summer climate season
C      ='icautt' initial values using autumn climate typical week
C      ='icautf' initial values using autumn climate fortnight
C      ='icauts' initial values using autumn climate season
C      ='ias' initial annual (365 days)
C      ='i3s' initial 3 season using climate seasons
C      ='i3t' initial 3 season using climate typical week
C      ='i3f' initial 3 season using climate typical fortnight
C      ='i5s' initial 5 season using climate seasons
C      ='i5t' initial 5 season using climate typical week
C      ='i5f' initial 5 season using climate typical fortnight
C      ='i--' initial values for call from interface
      subroutine ipvdatinit(act)
#include "building.h"
#include "model.h"
#include "geometry.h"
#include "seasons.h"
C seasons.h provides typper and typsea
#include "ipvdata.h"
#include "help.h"
      
      integer lnblnk  ! function definition

C Simulation parameter sets.
      common/spfldat/nsset,isset,isstup,isbnstep,ispnstep,issave,isavgh
      INTEGER :: nsset,isset,isstup,isbnstep,ispnstep,issave,isavgh
      integer ncomp,ncon
      common/c1/ncomp,ncon

C Extended simulation parameters for each set.
      common/spfldats/isstupex(MSPS),isbnstepex(MSPS),ispnstepex(MSPS),
     &  issaveex(MSPS),isavghex(MSPS),iscfdactivate(MSPS),
     &  isicfdys(MSPS),isicfdyf(MSPS),
     &  scftims(MSPS),scftimf(MSPS)
      INTEGER :: isstupex,isbnstepex,ispnstepex,issaveex,isavghex
      INTEGER :: iscfdactivate      ! zero ignore domains
      INTEGER :: isicfdys,isicfdyf  ! CFD simulation start & finish days
      REAL :: scftims,scftimf       ! CFD simulation start & finish time
      integer IVALS(MCOM)  ! the array of zones to include

      character act*6,cr*32

      helpinsub='edipv'  ! set for subroutine
      helptopic='ipv_initial'
      call gethelptext(helpinsub,helptopic,nbhelp)

C Set root name string and length.
      cr=cfgroot
      lr=lnblnk(cfgroot)

C Choose focus zone.
      INPIC=1
      CALL EPICKS(INPIC,IVALS,' ',
     &  'Which zone is focus for visual metrics: ',
     &  12,NCOMP,zname,' zone list',IER,nbhelp)
      if(ivals(1).ne.0)then
        ipvfoczone=ivals(1)
        write(ivpfocname,'(a)') zname(ipvfoczone)
      endif

      if(act(1:3).eq.'i--'.or.act(1:3).eq.'I--')then

C If there are 1/3/5 simulation parameter sets then assume the
C same for number of assessments, otherwise set to one annual.
        if(nsset.eq.1)then
          nipvassmt=1
          CALL EDAY(1,1,ipvastjd(1))
          CALL EDAY(31,12,ipvafnjd(1))
          write(ipvadesc(1),'(A)') 'annual '
          nipvdispjd=5
          ipvdispjd(1)=ia1wins+1
          ipvdispjd(2)=ia1sprs+1
          ipvdispjd(3)=iasums+1
          ipvdispjd(4)=ia2sprs+1
          ipvdispjd(5)=ia2wins+1
        elseif(nsset.eq.3)then
          nipvassmt=3
          ipvastjd(1)=ia1wins
          ipvafnjd(1)=ia1winf
          ipvastjd(2)=ia1sprs
          ipvafnjd(2)=ia1sprf
          ipvastjd(3)=iasums
          ipvafnjd(3)=iasumf
          nipvdispjd=3
          ipvdispjd(1)=ia1wins+1
          ipvdispjd(2)=ia1sprs+1
          ipvdispjd(3)=iasums+1
          write(ipvadesc(1),'(A)') 'winter'
          write(ipvadesc(2),'(A)') 'transition'
          write(ipvadesc(3),'(A)') 'summer'
        elseif(nsset.eq.5)then
          nipvassmt=5
          ipvastjd(1)=ia1wins
          ipvafnjd(1)=ia1winf
          ipvastjd(2)=ia1sprs
          ipvafnjd(2)=ia1sprf
          ipvastjd(3)=iasums
          ipvafnjd(3)=iasumf
          ipvastjd(4)=ia2sprs
          ipvafnjd(4)=ia2sprf
          ipvastjd(5)=ia2wins
          ipvafnjd(5)=ia2winf
          nipvdispjd=5
          ipvdispjd(1)=ia1wins+1
          ipvdispjd(2)=ia1sprs+1
          ipvdispjd(3)=iasums+1
          ipvdispjd(4)=ia2sprs+1
          ipvdispjd(5)=ia2wins+1
          write(ipvadesc(1),'(A)') '1st winter'
          write(ipvadesc(2),'(A)') 'spring'
          write(ipvadesc(3),'(A)') 'summer'
          write(ipvadesc(4),'(A)') 'autumn'
          write(ipvadesc(5),'(A)') '2nd winter'
        else

C Something other than 1/3/5
          nipvassmt=1
          CALL EDAY(1,1,ipvastjd(1))
          CALL EDAY(31,12,ipvafnjd(1))
          nipvdispjd=5
          ipvdispjd(1)=ia1wins+1
          ipvdispjd(2)=ia1sprs+1
          ipvdispjd(3)=iasums+1
          ipvdispjd(4)=ia2sprs+1
          ipvdispjd(5)=ia2wins+1
          write(ipvadesc(1),'(A)') 'annual '
        endif
      elseif(act(1:6).eq.'icwint'.or.act(1:6).eq.'Icwint')then
        nipvassmt=1    ! winter typical week
        ipvastjd(1)=ia1wins
        ipvafnjd(1)=ia1winf
        nipvdispjd=1
        ipvdispjd(1)=ia1wins+1
        write(ipvadesc(1),'(A)') 'winter week'
      elseif(act(1:6).eq.'icwinf'.or.act(1:6).eq.'Icwinf')then
        nipvassmt=1    ! winter typical fortnight
        if(ia1wins.gt.3)then
          ipvastjd(1)=ia1wins-3
          ipvafnjd(1)=ia1winf+4
        else
          ipvastjd(1)=ia1wins
          ipvafnjd(1)=ia1winf+7
        endif
        nipvdispjd=1
        ipvdispjd(1)=ia1wins+1
        write(ipvadesc(1),'(A)') 'winter fortnight'
      elseif(act(1:6).eq.'icwins'.or.act(1:6).eq.'Icwins')then
        nipvassmt=1    ! winter (early year) full season
        ipvastjd(1)=is1wins
        ipvafnjd(1)=is1winf
        nipvdispjd=1
        ipvdispjd(1)=ia1wins+1
        write(ipvadesc(1),'(A)') 'winter season'
      elseif(act(1:6).eq.'icsprt'.or.act(1:6).eq.'Icsprt')then
        nipvassmt=1    ! spring typical week
        ipvastjd(1)=ia1sprs
        ipvafnjd(1)=ia1sprf
        nipvdispjd=1
        ipvdispjd(1)=ia1sprs+1
        write(ipvadesc(1),'(A)') 'spring week'
      elseif(act(1:6).eq.'icsprf'.or.act(1:6).eq.'Icsprf')then
        nipvassmt=1    ! spring typical fortnight
        ipvastjd(1)=ia1sprs-3
        ipvafnjd(1)=ia1sprf+4
        nipvdispjd=1
        ipvdispjd(1)=ia1sprs+1
        write(ipvadesc(1),'(A)') 'spring fortnight'
      elseif(act(1:6).eq.'icsprs'.or.act(1:6).eq.'Icsprs')then
        nipvassmt=1    ! spring full season
        ipvastjd(1)=is1sprs
        ipvafnjd(1)=is1sprf
        nipvdispjd=1
        ipvdispjd(1)=ia1sprs+1
        write(ipvadesc(1),'(A)') 'spring season'
      elseif(act(1:6).eq.'icsumt'.or.act(1:6).eq.'Icsumt')then
        nipvassmt=1    ! summer typical week
        ipvastjd(1)=iasums
        ipvafnjd(1)=iasumf
        nipvdispjd=1
        ipvdispjd(1)=iasums+1
        write(ipvadesc(1),'(A)') 'summer week'
      elseif(act(1:6).eq.'icsumf'.or.act(1:6).eq.'Icsumf')then
        nipvassmt=1   ! summer typical fortnight
        ipvastjd(1)=iasums-3
        ipvafnjd(1)=iasumf+4
        nipvdispjd=1
        ipvdispjd(1)=iasums+1
        write(ipvadesc(1),'(A)') 'summer fortnight'
      elseif(act(1:6).eq.'icsums'.or.act(1:6).eq.'Icsums')then
        nipvassmt=1    ! summer full season
        ipvastjd(1)=is1sums
        ipvafnjd(1)=is1sumf
        nipvdispjd=1
        ipvdispjd(1)=iasums+1
        write(ipvadesc(1),'(A)') 'summer season'
      elseif(act(1:6).eq.'icautt'.or.act(1:6).eq.'Icautt')then
        nipvassmt=1    ! autumn typical week
        ipvastjd(1)=ia2sprs
        ipvafnjd(1)=ia2sprf
        nipvdispjd=1
        ipvdispjd(1)=ia2sprs+1
        write(ipvadesc(1),'(A)') 'autumn week'
      elseif(act(1:6).eq.'icautf'.or.act(1:6).eq.'Icautf')then
        nipvassmt=1   ! autumn typical fortnight
        ipvastjd(1)=ia2sprs-3
        ipvafnjd(1)=ia2sprf+4
        nipvdispjd=1
        ipvdispjd(1)=ia2sprs+1
        write(ipvadesc(1),'(A)') 'autumn fortnight'
      elseif(act(1:6).eq.'icauts'.or.act(1:6).eq.'Icauts')then
        nipvassmt=1    ! autumn full season
        ipvastjd(1)=is2sprs
        ipvafnjd(1)=is2sprf
        nipvdispjd=1
        ipvdispjd(1)=ia2sprs+1
        write(ipvadesc(1),'(A)') 'autumn season'
      elseif(act(1:3).eq.'ias'.or.act(1:3).eq.'Ias')then
        nipvassmt=1    ! annual (full year)
        CALL EDAY(1,1,ipvastjd(1))
        CALL EDAY(31,12,ipvafnjd(1))
        nipvdispjd=5
        ipvdispjd(1)=ia1wins+1
        ipvdispjd(2)=ia1sprs+1
        ipvdispjd(3)=iasums+1
        ipvdispjd(4)=ia2sprs+1
        ipvdispjd(5)=ia2wins+1
        write(ipvadesc(1),'(A)') 'annual'
      elseif(act(1:3).eq.'i3s'.or.act(1:3).eq.'I3s')then
        nipvassmt=3    ! three seasons (all days in each)
        ipvastjd(1)=is1wins
        ipvafnjd(1)=is1winf
        ipvastjd(2)=is1sprs
        ipvafnjd(2)=is1sprf
        ipvastjd(3)=is1sums
        ipvafnjd(3)=is1sumf
        nipvdispjd=3
        ipvdispjd(1)=ia1wins+1
        ipvdispjd(2)=ia1sprs+1
        ipvdispjd(3)=iasums+1
        write(ipvadesc(1),'(A)') 'winter season'
        write(ipvadesc(2),'(A)') 'transition season'
        write(ipvadesc(3),'(A)') 'summer season'
      elseif(act(1:3).eq.'i3t'.or.act(1:3).eq.'I3t')then
        nipvassmt=3    ! three seasons (one week in each)
        ipvastjd(1)=ia1wins
        ipvafnjd(1)=ia1winf
        ipvastjd(2)=ia1sprs
        ipvafnjd(2)=ia1sprf
        ipvastjd(3)=iasums
        ipvafnjd(3)=iasumf
        nipvdispjd=3
        ipvdispjd(1)=ia1wins+1
        ipvdispjd(2)=ia1sprs+1
        ipvdispjd(3)=iasums+1
        write(ipvadesc(1),'(A)') 'winter week'
        write(ipvadesc(2),'(A)') 'transition week'
        write(ipvadesc(3),'(A)') 'summer week'
      elseif(act(1:3).eq.'i3f'.or.act(1:3).eq.'I3f')then
        nipvassmt=3    ! three seasons (fortnight in each)
        ipvastjd(1)=ia1wins
        ipvafnjd(1)=ia1winf
        ipvastjd(2)=ia1sprs
        ipvafnjd(2)=ia1sprf
        ipvastjd(3)=iasums
        ipvafnjd(3)=iasumf
        nipvdispjd=3
        ipvdispjd(1)=ia1wins+1
        ipvdispjd(2)=ia1sprs+1
        ipvdispjd(3)=iasums+1
        write(ipvadesc(1),'(A)') 'winter fortnight'
        write(ipvadesc(2),'(A)') 'transition frtnt'
        write(ipvadesc(3),'(A)') 'summer fortnight'
      elseif(act(1:3).eq.'i5s'.or.act(1:3).eq.'I5s')then
        nipvassmt=5    ! five seasons (all days in each)
        ipvastjd(1)=is1wins
        ipvafnjd(1)=is1winf
        ipvastjd(2)=is1sprs
        ipvafnjd(2)=is1sprf
        ipvastjd(3)=is1sums
        ipvafnjd(3)=is1sumf
        ipvastjd(4)=is2sprs
        ipvafnjd(4)=is2sprf
        ipvastjd(5)=is2wins
        ipvafnjd(5)=is2winf
        nipvdispjd=5
        ipvdispjd(1)=is1wins+1
        ipvdispjd(2)=is1sprs+1
        ipvdispjd(3)=is1sums+1
        ipvdispjd(4)=is2sprs+1
        ipvdispjd(5)=is2wins+1
        write(ipvadesc(1),'(A)') '1st winter season'
        write(ipvadesc(2),'(A)') 'spring season'
        write(ipvadesc(3),'(A)') 'summer season'
        write(ipvadesc(4),'(A)') 'autumn season'
        write(ipvadesc(5),'(A)') '2nd winter season'
      elseif(act(1:3).eq.'i5t'.or.act(1:3).eq.'I5t')then
        nipvassmt=5    ! five seasons (one week in each)
        ipvastjd(1)=ia1wins
        ipvafnjd(1)=ia1winf
        ipvastjd(2)=ia1sprs
        ipvafnjd(2)=ia1sprf
        ipvastjd(3)=iasums
        ipvafnjd(3)=iasumf
        ipvastjd(4)=ia2sprs
        ipvafnjd(4)=ia2sprf
        ipvastjd(5)=ia2wins
        ipvafnjd(5)=ia2winf
        nipvdispjd=5
        ipvdispjd(1)=ia1wins+1
        ipvdispjd(2)=ia1sprs+1
        ipvdispjd(3)=iasums+1
        ipvdispjd(4)=ia2sprs+1
        ipvdispjd(5)=ia2wins+1
        write(ipvadesc(1),'(A)') '1st winter week'
        write(ipvadesc(2),'(A)') 'spring week'
        write(ipvadesc(3),'(A)') 'summer week'
        write(ipvadesc(4),'(A)') 'autumn week'
        write(ipvadesc(5),'(A)') '2nd winter week'
      elseif(act(1:3).eq.'i5f'.or.act(1:3).eq.'I5f')then
        nipvassmt=5    ! five seasons (fortnight in each)
        ipvastjd(1)=ia1wins
        ipvafnjd(1)=ia1winf
        ipvastjd(2)=ia1sprs
        ipvafnjd(2)=ia1sprf
        ipvastjd(3)=iasums
        ipvafnjd(3)=iasumf
        ipvastjd(4)=ia2sprs
        ipvafnjd(4)=ia2sprf
        ipvastjd(5)=ia2wins
        ipvafnjd(5)=ia2winf
        nipvdispjd=5
        ipvdispjd(1)=ia1wins+1
        ipvdispjd(2)=ia1sprs+1
        ipvdispjd(3)=iasums+1
        ipvdispjd(4)=ia2sprs+1
        ipvdispjd(5)=ia2wins+1
        write(ipvadesc(1),'(A)') '1st winter frtnt'
        write(ipvadesc(2),'(A)') 'spring fortnight'
        write(ipvadesc(3),'(A)') 'summer fortnight'
        write(ipvadesc(4),'(A)') 'autumn fortnight'
        write(ipvadesc(5),'(A)') '2nd winter frtnt'
      endif

C Reset ddm* data structures for assessments that span a
C whole season.
      if(act(1:6).eq.'icwins'.or.act(1:6).eq.'Icwins')then
        continue
      elseif(act(1:6).eq.'icsprs'.or.act(1:6).eq.'Icsprs')then
        continue
      elseif(act(1:6).eq.'icsums'.or.act(1:6).eq.'Icsums')then
        continue
      elseif(act(1:6).eq.'icauts'.or.act(1:6).eq.'Icauts')then
        continue
      elseif(act(1:3).eq.'ias'.or.act(1:3).eq.'Ias')then
        continue
      elseif(act(1:3).eq.'i3s'.or.act(1:3).eq.'I3s')then
        continue
      elseif(act(1:3).eq.'i5s'.or.act(1:3).eq.'I5s')then
        continue
      else
        return
      endif
      do 106 ij=1,5
        ddmheat(ij)=1.0
        ddmcool(ij)=1.0
        ddmtime(ij)=1.0
  106 continue
      return
      end


C ************* IPV2SIMPAR
C Copy relevant data from IPV description to simulation parameter
C sets. See subroutine IPVDAT for common block descriptions.
C act is one of the character tags from ipvseasons or recovered
C from the IPV description file.

      subroutine ipv2simpar(act)
#include "building.h"
#include "model.h"
#include "ipvdata.h"
#include "plant.h"
#include "power.h"

      integer lnblnk  ! function definition

      COMMON/FILEP/IFIL
      common/OUTIN/IUOUT,IUIN,IEOUT
      COMMON/C6/INDCFG
      COMMON/PREC7/ITCNST

C IPV description via ipvdata.h.
      common/IPVF/lipvdatf
      COMMON/MOIST01/MSTROK,MSTRZN(MCOM)
      LOGICAL MSTROK,MSTRZN
      COMMON/AFN/IAIRN,LAPROB,ICAAS(MCOM)
      INTEGER :: iairn,icaas

C Simulation parameter sets.
      common/ACTDOM/CFDOK
      logical CFDOK
      common/spfldes/spfdescr(MSPS)
      common/spflper/isstday(MSPS),isstmon(MSPS),isfnday(MSPS),
     &               isfnmon(MSPS)
      common/spfldat/nsset,isset,isstup,isbnstep,ispnstep,issave,isavgh
      INTEGER :: nsset,isset,isstup,isbnstep,ispnstep,issave,isavgh
      common/spflres/sblres(MSPS),sflres(MSPS),splres(MSPS),
     &  smstres(MSPS),selres(MSPS),scfdres(MSPS),sipvres
      character sblres*72,sflres*72,splres*72,smstres*72,
     &  selres*72,scfdres*72,sipvres*72

C Extended simulation parameters for each set.
      common/spfldats/isstupex(MSPS),isbnstepex(MSPS),ispnstepex(MSPS),
     &  issaveex(MSPS),isavghex(MSPS),iscfdactivate(MSPS),
     &  isicfdys(MSPS),isicfdyf(MSPS),
     &  scftims(MSPS),scftimf(MSPS)
      INTEGER :: isstupex,isbnstepex,ispnstepex,issaveex,isavghex
      INTEGER :: iscfdactivate      ! zero ignore domains
      INTEGER :: isicfdys,isicfdyf  ! CFD simulation start & finish days
      REAL :: scftims,scftimf       ! CFD simulation start & finish time

      character act*6,lipvdatf*72,cr*32,tag*4
      character spfdescr*30,outs*124
      CHARACTER LAPROB*72
      logical unixok
      logical shift_needed
      integer old_retained,next_avail

      if(nipvassmt.eq.1.or.nipvassmt.eq.3.or.nipvassmt.eq.5)then
        continue
      else
        call usrmsg('Request not for one, three or five simulations',
     &    'so unsure how to setup simulation parameters.','W')
        return
      endif

C If simulation parameter set save level is still -1 then there has
C probably not been a simulation parameter set defined, so set some
C of the common variables.
      if(issave.eq.-1)then

C This block of code may not be necessary.
        call edisp(iuout,'Scanning for startup days')
        call scntcnst(TDM,istd,TCM,ISTC,ITCN)
        if(isstup.eq.0)isstup=ITCNST
        if(isstupex(1).eq.0)isstupex(1)=ITCNST
        isbnstep=4;ispnstep=10
        issave=4; isavgh=0
        isbnstepex(1)=4; isbnstepex(2)=4;isbnstepex(3)=4
        isbnstepex(4)=4; isbnstepex(5)=4
        ispnstepex(1)=10; ispnstepex(2)=10; ispnstepex(3)=10
        ispnstepex(4)=10; ispnstepex(5)=10
        issaveex(1)=4; issaveex(2)=4; issaveex(3)=4
        issaveex(4)=4; issaveex(5)=4
        isavghex(1)=0; isavghex(2)=0; isavghex(3)=0
        isavghex(4)=0; isavghex(5)=0
      else
        call edisp(iuout,'Scanning for startup days')
        call scntcnst(TDM,istd,TCM,ISTC,ITCN)
        if(isstup.eq.0)isstup=ITCNST
      endif

      cr=cfgroot
      lr=lnblnk(cfgroot)
      if(nsset.eq.0)then

C If no simulation parameter sets then instanciate for the usual slots.
        nsset=nipvassmt            ! SPS are all from IPV definitions.
        iupd_win1=-1; iupd_spr=-1; iupd_sum=-1; iupd_aut=-1
        iupd_win2=-1; iupd_ann=-1  ! Mark assuming no matches.
        if(nipvassmt.eq.1)then     ! Assign usual slots.
          iupd_win1=1; iupd_spr=1; iupd_sum=1; iupd_aut=1
          iupd_win2=1; iupd_ann=1  ! Mark for first slot.
        elseif(nipvassmt.eq.3)then
          iupd_win1=1; iupd_spr=2; iupd_sum=3
        elseif(nipvassmt.eq.5)then
          iupd_win1=1; iupd_spr=2; iupd_sum=3
          iupd_aut=4; iupd_win2=5
        endif
C        write(6,*) 'win1 spr sum aut win2 ann',iupd_win1,iupd_spr,
C     &    iupd_sum,iupd_aut,iupd_win2,iupd_ann
      elseif(nsset.gt.0)then

C If there are existing SPS loop thorough to see which ones match
C the key phrases used with an IPV.  If nipvassmt is 1 then
C it should be the first slot (so move any existing SPS
C up to make room if slot one is not alredy an IPV slot.
C If nipvassmt is 3 then the first three slots should be for IPV.
C If nipvassmt is 5 then the first five slots should be for IPV.
        iupd_win1=-1; iupd_spr=-1; iupd_sum=-1; iupd_aut=-1
        iupd_win2=-1; iupd_ann=-1  ! Mark assuming no matches.
        next_avail=nsset           ! End of current SPS.
        old_retained=0
        do ij=1,nsset              ! loop to see which SPS can be updated.
          lde=lnblnk(spfdescr(ij))
          if(spfdescr(ij)(1:lde).eq.'win1 ')then
            write(outs,'(a,i2,a)') ' Existing win1 SPS ',ij,' to update'
            call edisp(iuout,outs)
            iupd_win1=ij
          elseif(spfdescr(ij)(1:lde).eq.'spr')then
            write(outs,'(a,i2,a)') ' Existing spr SPS ',ij,' to update'
            call edisp(iuout,outs)
            iupd_spr=ij
          elseif(spfdescr(ij)(1:lde).eq.'sum ')then
            write(outs,'(a,i2,a)') ' Existing sum SPS ',ij,' to update'
            call edisp(iuout,outs)
            iupd_sum=ij
          elseif(spfdescr(ij)(1:lde).eq.'aut ')then
            write(outs,'(a,i2,a)') ' Existing aut SPS ',ij,' to update'
            call edisp(iuout,outs)
            iupd_aut=ij
          elseif(spfdescr(ij)(1:lde).eq.'win2 ')then
            write(outs,'(a,i2,a)') ' Existing win2 SPS ',ij,' to update'
            call edisp(iuout,outs)
            iupd_win2=ij
          elseif(spfdescr(ij)(1:lde).eq.'win ')then
            write(outs,'(a,i2,a)') ' Existing win SPS ',ij,' to update'
            call edisp(iuout,outs)
            iupd_win1=ij
          elseif(spfdescr(ij)(1:lde).eq.'ann ')then
            write(outs,'(a,i2,a)') ' Existing ann SPS ',ij,' to update'
            call edisp(iuout,outs)
            iupd_win2=ij
          else
            old_retained=old_retained+1
            write(outs,'(a,i2,2a)')' Existing SPS ',ij,spfdescr(ij),
     &        ' retained.'
            call edisp(iuout,outs)
          endif
        enddo
        write(6,*) 'win1 spr sum aut win2 ann ',iupd_win1,iupd_spr,
     &    iupd_sum,iupd_aut,iupd_win2,iupd_ann

C Second pass to identify where new IPV slots can be placed.
C If a single IPV assessment then it should go in slot 1.
        if(nipvassmt.eq.1)then     ! If single IPV assessment
          islot=-1                 ! SPS index initial state
          shift_needed=.false.
          if(act(1:5).eq.'icwin'.or.act(1:5).eq.'Icwin')then
            if(iupd_win1.eq.-1)then  ! Does not already exist so set islot.
              iupd_win1=1
              islot=iupd_win1; shift_needed=.true.
            else
              islot=iupd_win1        ! Already exists, set islot to match.
            endif
            tag='win'
          elseif(act(1:5).eq.'icspr'.or.act(1:5).eq.'Icspr')then
            if(iupd_spr.eq.-1)then
              iupd_spr=1
              islot=iupd_spr; shift_needed=.true.
            else
              islot=iupd_spr
            endif
            tag='spr'
          elseif(act(1:5).eq.'icsum'.or.act(1:5).eq.'Icsum')then
            if(iupd_sum.eq.-1)then
              iupd_sum=1
              islot=iupd_sum; shift_needed=.true.
            else
              islot=iupd_sum
            endif
            tag='sum'
          elseif(act(1:5).eq.'icaut'.or.act(1:5).eq.'Icaut')then
            if(iupd_aut.eq.-1)then
              iupd_aut=1
              islot=iupd_aut; shift_needed=.true.
            else
              islot=iupd_aut
            endif
            tag='aut'
          elseif(act(1:3).eq.'ias'.or.act(1:3).eq.'Ias')then
            if(iupd_ann.eq.-1)then
              iupd_ann=1
              islot=iupd_ann; shift_needed=.true.
            else
              islot=iupd_ann
            endif
            tag='ann'
          endif

C If islot is still -1 then the past-in act has a match with an existing
C SPS and no shifting of subsequent SPS are required. Otherwise we need to
C shift subsequent SPS to make space for the IPV single assessment. For
C a single IPV assessment copy into 'islot' the 1st IPV array values.
          if(islot.gt.0)then
            if(shift_needed)then
              next_avail=next_avail+1
              lower=1
              call shift_sps(lower,next_avail,ier)   ! Make room.
              write(outs,'(a,i2,a)') ' Insert single SPS @',islot
              call edisp(iuout,outs)
            else
              write(outs,'(a,i2,a)') ' Update single SPS @',islot
              call edisp(iuout,outs)
            endif
            write(spfdescr(islot),'(a)') tag
            isbnstepex(islot)=4
            ispnstepex(islot)=10
            issaveex(islot)=4
            isavghex(islot)=0
            call EDAYR(ipvastjd(1),isstday(islot),isstmon(islot))
            call EDAYR(ipvafnjd(1),isfnday(islot),isfnmon(islot))
            isstupex(islot)=ITCNST
            iscfdactivate(islot)=0 
            isicfdys(islot)=ipvastjd(1); isicfdyf(islot)=ipvafnjd(1)
            scftims(islot)=1.0; scftimf(islot)=23.0
            write(6,*) 'rev attributes ',islot,spfdescr(islot),
     &        isbnstepex(islot),ispnstepex(islot),issaveex(islot),
     &        isavghex(islot),isstday(islot),isstmon(islot),
     &        isfnday(islot),isfnmon(islot),isstupex(islot)
          endif
        elseif(nipvassmt.eq.3)then     ! If three IPV assessments
          islot=-1                     ! SPS index initial state
          shift_needed=.false.
          if(iupd_win1.eq.-1)then  ! Does not already exist so set islot.
            iupd_win1=1
            islot=iupd_win1; shift_needed=.true.
          else
            islot=iupd_win1        ! Already exists, set islot to match.
          endif
          tag='win'
          if(islot.gt.0)then
            if(shift_needed)then
              next_avail=next_avail+1
              lower=1
              call shift_sps(lower,next_avail,ier)   ! Make room.
              write(outs,'(a,i2,a)') ' Insert win SPS @',islot
              call edisp(iuout,outs)
            else
              write(outs,'(a,i2,a)') ' Update win SPS @',islot
              call edisp(iuout,outs)
            endif
            write(spfdescr(islot),'(a)') tag
            isbnstepex(islot)=4
            ispnstepex(islot)=10
            issaveex(islot)=4
            isavghex(islot)=0
            call EDAYR(ipvastjd(1),isstday(islot),isstmon(islot))
            call EDAYR(ipvafnjd(1),isfnday(islot),isfnmon(islot))
            isstupex(islot)=ITCNST
            iscfdactivate(islot)=0 
            isicfdys(islot)=ipvastjd(1); isicfdyf(islot)=ipvafnjd(1)
            scftims(islot)=1.0; scftimf(islot)=23.0
            write(6,*) 'rev attributes ',islot,spfdescr(islot),
     &        isbnstepex(islot),ispnstepex(islot),issaveex(islot),
     &        isavghex(islot),isstday(islot),isstmon(islot),
     &        isfnday(islot),isfnmon(islot),isstupex(islot)
          endif

C Similar logic for spring.
          islot=-1                     ! SPS index initial state
          shift_needed=.false.
          if(iupd_spr.eq.-1)then  ! Does not already exist so set islot.
            iupd_spr=2
            islot=iupd_spr; shift_needed=.true.
          else
            islot=iupd_spr        ! Already exists, set islot to match.
          endif
          tag='spr'
          if(islot.gt.0)then
            if(shift_needed)then
              next_avail=next_avail+1
              lower=2
              call shift_sps(lower,next_avail,ier)   ! Make room.
              write(outs,'(a,i2,a)') ' Insert spring SPS @',islot
              call edisp(iuout,outs)
            else
              write(outs,'(a,i2,a)') ' Update spring SPS @',islot
              call edisp(iuout,outs)
            endif
            write(spfdescr(islot),'(a)') tag
            isbnstepex(islot)=4
            ispnstepex(islot)=10
            issaveex(islot)=4
            isavghex(islot)=0
            call EDAYR(ipvastjd(2),isstday(islot),isstmon(islot))
            call EDAYR(ipvafnjd(2),isfnday(islot),isfnmon(islot))
            isstupex(islot)=ITCNST
            iscfdactivate(islot)=0 
            isicfdys(islot)=ipvastjd(2); isicfdyf(islot)=ipvafnjd(2)
            scftims(islot)=1.0; scftimf(islot)=23.0
C            write(6,*) 'rev attributes ',islot,spfdescr(islot),
C     &        isbnstepex(islot),ispnstepex(islot),issaveex(islot),
C     &        isavghex(islot),isstday(islot),isstmon(islot),
C     &        isfnday(islot),isfnmon(islot),isstupex(islot)
          endif

C Similar logic for summmer.
          islot=-1                     ! SPS index initial state
          shift_needed=.false.
          if(iupd_sum.eq.-1)then  ! Does not already exist so set islot.
            iupd_sum=3
            islot=iupd_sum; shift_needed=.true.
          else
            islot=iupd_sum        ! Already exists, set islot to match.
          endif
          tag='sum'
          if(islot.gt.0)then
            if(shift_needed)then
              next_avail=next_avail+1
              lower=3
              call shift_sps(lower,next_avail,ier)   ! Make room.
              write(outs,'(a,i2,a)') ' Insert summer SPS @',islot
              call edisp(iuout,outs)
            else
              write(outs,'(a,i2,a)') ' Update summer SPS @',islot
              call edisp(iuout,outs)
            endif
            write(spfdescr(islot),'(a)') tag
            isbnstepex(islot)=4
            ispnstepex(islot)=10
            issaveex(islot)=4
            isavghex(islot)=0
            call EDAYR(ipvastjd(3),isstday(islot),isstmon(islot))
            call EDAYR(ipvafnjd(3),isfnday(islot),isfnmon(islot))
            isstupex(islot)=ITCNST
            iscfdactivate(islot)=0 
            isicfdys(islot)=ipvastjd(3); isicfdyf(islot)=ipvafnjd(3)
            scftims(islot)=1.0; scftimf(islot)=23.0
C            write(6,*) 'rev attributes ',islot,spfdescr(islot),
C     &        isbnstepex(islot),ispnstepex(islot),issaveex(islot),
C     &        isavghex(islot),isstday(islot),isstmon(islot),
C     &        isfnday(islot),isfnmon(islot),isstupex(islot)
          endif

        elseif(nipvassmt.eq.5)then     ! If five IPV assessments
          islot=-1                     ! SPS index initial state
          shift_needed=.false.
          if(iupd_win1.eq.-1)then  ! Does not already exist so set islot.
            iupd_win1=1
            islot=iupd_win1; shift_needed=.true.
          else
            islot=iupd_win1        ! Already exists, set islot to match.
          endif
          tag='win1'
          if(islot.gt.0)then
            if(shift_needed)then
              next_avail=next_avail+1
              lower=1
              call shift_sps(lower,next_avail,ier)   ! Make room.
              write(outs,'(a,i2,a)') ' Insert win1 SPS @',islot
              call edisp(iuout,outs)
            else
              write(outs,'(a,i2,a)') ' Update win1 SPS @',islot
              call edisp(iuout,outs)
            endif
            write(spfdescr(islot),'(a)') tag
            isbnstepex(islot)=4
            ispnstepex(islot)=10
            issaveex(islot)=4
            isavghex(islot)=0
            call EDAYR(ipvastjd(1),isstday(islot),isstmon(islot))
            call EDAYR(ipvafnjd(1),isfnday(islot),isfnmon(islot))
            isstupex(islot)=ITCNST
            iscfdactivate(islot)=0 
            isicfdys(islot)=ipvastjd(1); isicfdyf(islot)=ipvafnjd(1)
            scftims(islot)=1.0; scftimf(islot)=23.0
C            write(6,*) 'rev attributes ',islot,spfdescr(islot),
C     &        isbnstepex(islot),ispnstepex(islot),issaveex(islot),
C     &        isavghex(islot),isstday(islot),isstmon(islot),
C     &        isfnday(islot),isfnmon(islot),isstupex(islot)
          endif

C Similar logic for spring.
          islot=-1                     ! SPS index initial state
          shift_needed=.false.
          if(iupd_spr.eq.-1)then  ! Does not already exist so set islot.
            iupd_spr=2
            islot=iupd_spr; shift_needed=.true.
          else
            islot=iupd_spr        ! Already exists, set islot to match.
          endif
          tag='spr'
          if(islot.gt.0)then
            if(shift_needed)then
              next_avail=next_avail+1
              lower=2
              call shift_sps(lower,next_avail,ier)   ! Make room.
              write(outs,'(a,i2,a)') ' Insert spring SPS @',islot
              call edisp(iuout,outs)
            else
              write(outs,'(a,i2,a)') ' Update spring SPS @',islot
              call edisp(iuout,outs)
            endif
            write(spfdescr(islot),'(a)') tag
            isbnstepex(islot)=4
            ispnstepex(islot)=10
            issaveex(islot)=4
            isavghex(islot)=0
            call EDAYR(ipvastjd(2),isstday(islot),isstmon(islot))
            call EDAYR(ipvafnjd(2),isfnday(islot),isfnmon(islot))
            isstupex(islot)=ITCNST
            iscfdactivate(islot)=0 
            isicfdys(islot)=ipvastjd(2); isicfdyf(islot)=ipvafnjd(2)
            scftims(islot)=1.0; scftimf(islot)=23.0
C            write(6,*) 'rev attributes ',islot,spfdescr(islot),
C     &        isbnstepex(islot),ispnstepex(islot),issaveex(islot),
C     &        isavghex(islot),isstday(islot),isstmon(islot),
C     &        isfnday(islot),isfnmon(islot),isstupex(islot)
          endif

C Similar logic for summmer.
          islot=-1                     ! SPS index initial state
          shift_needed=.false.
          if(iupd_sum.eq.-1)then  ! Does not already exist so set islot.
            iupd_sum=3
            islot=iupd_sum; shift_needed=.true.
          else
            islot=iupd_sum        ! Already exists, set islot to match.
          endif
          tag='sum'
          if(islot.gt.0)then
            if(shift_needed)then
              next_avail=next_avail+1
              lower=3
              call shift_sps(lower,next_avail,ier)   ! Make room.
              write(outs,'(a,i2,a)') ' Insert summer SPS @',islot
              call edisp(iuout,outs)
            else
              write(outs,'(a,i2,a)') ' Update summer SPS @',islot
              call edisp(iuout,outs)
            endif
            write(spfdescr(islot),'(a)') tag
            isbnstepex(islot)=4
            ispnstepex(islot)=10
            issaveex(islot)=4
            isavghex(islot)=0
            call EDAYR(ipvastjd(3),isstday(islot),isstmon(islot))
            call EDAYR(ipvafnjd(3),isfnday(islot),isfnmon(islot))
            isstupex(islot)=ITCNST
            iscfdactivate(islot)=0 
            isicfdys(islot)=ipvastjd(3); isicfdyf(islot)=ipvafnjd(3)
            scftims(islot)=1.0; scftimf(islot)=23.0
C            write(6,*) 'rev attributes ',islot,spfdescr(islot),
C     &        isbnstepex(islot),ispnstepex(islot),issaveex(islot),
C     &        isavghex(islot),isstday(islot),isstmon(islot),
C     &        isfnday(islot),isfnmon(islot),isstupex(islot)
          endif

C Similar logic for autumn.
          islot=-1                     ! SPS index initial state
          shift_needed=.false.
          if(iupd_aut.eq.-1)then  ! Does not already exist so set islot.
            iupd_aut=4
            islot=iupd_aut; shift_needed=.true.
          else
            islot=iupd_aut        ! Already exists, set islot to match.
          endif
          tag='aut'
          if(islot.gt.0)then
            if(shift_needed)then
              next_avail=next_avail+1
              lower=4
              call shift_sps(lower,next_avail,ier)   ! Make room.
              write(outs,'(a,i2,a)') ' Insert autumn SPS @',islot
              call edisp(iuout,outs)
            else
              write(outs,'(a,i2,a)') ' Update autumn SPS @',islot
              call edisp(iuout,outs)
            endif
            write(spfdescr(islot),'(a)') tag
            isbnstepex(islot)=4
            ispnstepex(islot)=10
            issaveex(islot)=4
            isavghex(islot)=0
            call EDAYR(ipvastjd(4),isstday(islot),isstmon(islot))
            call EDAYR(ipvafnjd(4),isfnday(islot),isfnmon(islot))
            isstupex(islot)=ITCNST
            iscfdactivate(islot)=0 
            isicfdys(islot)=ipvastjd(4); isicfdyf(islot)=ipvafnjd(4)
            scftims(islot)=1.0; scftimf(islot)=23.0
C            write(6,*) 'rev attributes ',islot,spfdescr(islot),
C     &        isbnstepex(islot),ispnstepex(islot),issaveex(islot),
C     &        isavghex(islot),isstday(islot),isstmon(islot),
C     &        isfnday(islot),isfnmon(islot),isstupex(islot)
          endif

C Similar logic for late year winter.
          islot=-1                     ! SPS index initial state
          shift_needed=.false.
          if(iupd_win2.eq.-1)then  ! Does not already exist so set islot.
            iupd_win2=5
            islot=iupd_win2; shift_needed=.true.
          else
            islot=iupd_win2        ! Already exists, set islot to match.
          endif
          tag='win2'
          if(islot.gt.0)then
            if(shift_needed)then
              next_avail=next_avail+1
              lower=5
              call shift_sps(lower,next_avail,ier)   ! Make room.
              write(outs,'(a,i2,a)') ' Insert late win SPS @',islot
              call edisp(iuout,outs)
            else
              write(outs,'(a,i2,a)') ' Update late win SPS @',islot
              call edisp(iuout,outs)
            endif
            write(spfdescr(islot),'(a)') tag
            isbnstepex(islot)=4
            ispnstepex(islot)=10
            issaveex(islot)=4
            isavghex(islot)=0
            call EDAYR(ipvastjd(5),isstday(islot),isstmon(islot))
            call EDAYR(ipvafnjd(5),isfnday(islot),isfnmon(islot))
            isstupex(islot)=ITCNST
            iscfdactivate(islot)=0 
            isicfdys(islot)=ipvastjd(5); isicfdyf(islot)=ipvafnjd(5)
            scftims(islot)=1.0; scftimf(islot)=23.0
            write(6,*) 'rev attributes ',islot,spfdescr(islot),
     &        isbnstepex(islot),ispnstepex(islot),issaveex(islot),
     &        isavghex(islot),isstday(islot),isstmon(islot),
     &        isfnday(islot),isfnmon(islot),isstupex(islot)
          endif
        endif
        
        write(6,*) 'win1 spr sum aut win2 ann islot',iupd_win1,iupd_spr,
     &    iupd_sum,iupd_aut,iupd_win2,iupd_ann,islot
      endif

C Update the associated results files to reflect the new tags.
      if(nipvassmt.eq.1)then
        call isunix(unixok)
        if(INDCFG.ne.2)then
          if(unixok)then
            WRITE(sblres(islot),'(5A)') '../tmp/',cr(1:lr),'_',tag,
     &        '.res'
          else
            WRITE(sblres(islot),'(4A)') cr(1:lr),'_',tag,'.res'
          endif
        endif
        write(6,*) islot,' ',sblres(islot)
        if(IAIRN.ge.1)then
          if(unixok)then
            WRITE(sflres(islot),'(5A)') '../tmp/',cr(1:lr),'_',tag,
     &       '.mfr'
          else
            WRITE(sflres(islot),'(4A)')cr(1:lr),'_',tag,'.mfr'
          endif
        endif
        write(6,*) islot,' ',sflres(islot)
        if(INDCFG.eq.2.or.INDCFG.eq.3)then
          if(unixok)then
            WRITE(splres(islot),'(5A)') '../tmp/',cr(1:lr),'_',tag,
     &        '.plr'
          else
            WRITE(splres(islot),'(4A)')cr(1:lr),'_',tag,'.plr'
          endif
        endif
        if(MSTROK)then
          if(unixok)then
            WRITE(smstres(islot),'(5A)') '../tmp/',cr(1:lr),'_',tag,
     &        '.msr'
          else
            WRITE(smstres(islot),'(4A)')cr(1:lr),'_',tag,'.msr'
          endif
        endif
        if(ientxist.gt.0)then
          if(unixok)then
            WRITE(selres(islot),'(5A)') '../tmp/',cr(1:lr),'_',tag,
     &        '.elr'
          else
            WRITE(selres(islot),'(4A)')cr(1:lr),'_',tag,'.elr'
          endif
        endif
        if(CFDOK)then
          if(unixok)then
            WRITE(scfdres(islot),'(5A)') '../tmp/',cr(1:lr),'_',tag,
     &        '.dfr'
          else
            WRITE(scfdres(islot),'(4A)')cr(1:lr),'_',tag,'.dfr'
          endif
        endif
      elseif(nipvassmt.eq.3)then
        call isunix(unixok)
        if(INDCFG.ne.2)then
          if(unixok)then
            WRITE(sblres(iupd_win1),'(3A)')'../tmp/',cr(1:lr),'_win.res'
            WRITE(sblres(iupd_spr),'(3A)') '../tmp/',cr(1:lr),'_trn.res'
            WRITE(sblres(iupd_sum),'(3A)') '../tmp/',cr(1:lr),'_sum.res'
          else
            WRITE(sblres(iupd_win1),'(2A)')cr(1:lr),'_win.res'
            WRITE(sblres(iupd_spr),'(2A)')cr(1:lr),'_trn.res'
            WRITE(sblres(iupd_sum),'(2A)')cr(1:lr),'_sum.res'
          endif
        endif
        if(IAIRN.ge.1)then
          if(unixok)then
            WRITE(sflres(iupd_win1),'(3A)')'../tmp/',cr(1:lr),'_win.mfr'
            WRITE(sflres(iupd_spr),'(3A)') '../tmp/',cr(1:lr),'_trn.mfr'
            WRITE(sflres(iupd_sum),'(3A)') '../tmp/',cr(1:lr),'_sum.mfr'
          else
            WRITE(sflres(iupd_win1),'(2A)')cr(1:lr),'_win.mfr'
            WRITE(sflres(iupd_spr),'(2A)')cr(1:lr),'_trn.mfr'
            WRITE(sflres(iupd_sum),'(2A)')cr(1:lr),'_sum.mfr'
          endif
        endif
        if(INDCFG.eq.2.or.INDCFG.eq.3)then
          if(unixok)then
            WRITE(splres(iupd_win1),'(3A)')'../tmp/',cr(1:lr),'_win.plr'
            WRITE(splres(iupd_spr),'(3A)') '../tmp/',cr(1:lr),'_trn.plr'
            WRITE(splres(iupd_sum),'(3A)') '../tmp/',cr(1:lr),'_sum.plr'
          else
            WRITE(splres(iupd_win1),'(2A)')cr(1:lr),'_win.plr'
            WRITE(splres(iupd_spr),'(2A)')cr(1:lr),'_trn.plr'
            WRITE(splres(iupd_sum),'(2A)')cr(1:lr),'_sum.plr'
          endif
        endif
        if(MSTROK)then
          if(unixok)then
            WRITE(smstres(iupd_win1),'(3A)')'../tmp/',cr(1:lr),
     &        '_win.msr'
            WRITE(smstres(iupd_spr),'(3A)') '../tmp/',cr(1:lr),
     &        '_trn.msr'
            WRITE(smstres(iupd_sum),'(3A)') '../tmp/',cr(1:lr),
     &         '_sum.msr'
          else
            WRITE(smstres(iupd_win1),'(2A)')cr(1:lr),'_win.msr'
            WRITE(smstres(iupd_spr),'(2A)')cr(1:lr),'_trn.msr'
            WRITE(smstres(iupd_sum),'(2A)')cr(1:lr),'_sum.msr'
          endif
        endif
        if(ientxist.gt.0)then
          if(unixok)then
            WRITE(selres(iupd_win1),'(3A)')'../tmp/',cr(1:lr),'_win.elr'
            WRITE(selres(iupd_spr),'(3A)') '../tmp/',cr(1:lr),'_trn.elr'
            WRITE(selres(iupd_sum),'(3A)') '../tmp/',cr(1:lr),'_sum.elr'
          else
            WRITE(selres(iupd_win1),'(2A)')cr(1:lr),'_win.elr'
            WRITE(selres(iupd_spr),'(2A)')cr(1:lr),'_trn.elr'
            WRITE(selres(iupd_sum),'(2A)')cr(1:lr),'_sum.elr'
          endif
        endif
        if(CFDOK)then
          if(unixok)then
            WRITE(scfdres(iupd_win1),'(3A)')'../tmp/',cr(1:lr),
     &        '_win.dfr'
            WRITE(scfdres(iupd_spr),'(3A)') '../tmp/',cr(1:lr),
     &        '_trn.dfr'
            WRITE(scfdres(iupd_sum),'(3A)') '../tmp/',cr(1:lr),
     &        '_sum.dfr'
          else
            WRITE(scfdres(iupd_win1),'(2A)')cr(1:lr),'_win.dfr'
            WRITE(scfdres(iupd_spr),'(2A)')cr(1:lr),'_trn.dfr'
            WRITE(scfdres(iupd_sum),'(2A)')cr(1:lr),'_sum.dfr'
          endif
        endif
      elseif(nipvassmt.eq.5)then
        call isunix(unixok)
        if(INDCFG.ne.2)then
          if(unixok)then
            WRITE(sblres(iupd_win1),'(3A)') '../tmp/',cr(1:lr),
     &        '_win1.res'
            WRITE(sblres(iupd_spr),'(3A)') '../tmp/',cr(1:lr),'_spr.res'
            WRITE(sblres(iupd_sum),'(3A)') '../tmp/',cr(1:lr),'_sum.res'
            WRITE(sblres(iupd_aut),'(3A)') '../tmp/',cr(1:lr),'_aut.res'
            WRITE(sblres(iupd_win2),'(3A)') '../tmp/',cr(1:lr),
     &        '_win2.res'
          else
            WRITE(sblres(iupd_win1),'(2A)')cr(1:lr),'_win1.res'
            WRITE(sblres(iupd_spr),'(2A)')cr(1:lr),'_spr.res'
            WRITE(sblres(iupd_sum),'(2A)')cr(1:lr),'_sum.res'
            WRITE(sblres(iupd_aut),'(2A)')cr(1:lr),'_aut.res'
            WRITE(sblres(iupd_win2),'(2A)')cr(1:lr),'_win2.res'
          endif
        endif
        if(IAIRN.ge.1)then
          if(unixok)then
            WRITE(sflres(iupd_win1),'(3A)') '../tmp/',cr(1:lr),
     &        '_win1.mfr'
            WRITE(sflres(iupd_spr),'(3A)') '../tmp/',cr(1:lr),'_spr.mfr'
            WRITE(sflres(iupd_sum),'(3A)') '../tmp/',cr(1:lr),'_sum.mfr'
            WRITE(sflres(iupd_aut),'(3A)') '../tmp/',cr(1:lr),'_aut.mfr'
            WRITE(sflres(iupd_win2),'(3A)') '../tmp/',cr(1:lr),
     &        '_win2.mfr'
          else
            WRITE(sflres(iupd_win1),'(2A)')cr(1:lr),'_win1.mfr'
            WRITE(sflres(iupd_spr),'(2A)')cr(1:lr),'_spr.mfr'
            WRITE(sflres(iupd_sum),'(2A)')cr(1:lr),'_sum.mfr'
            WRITE(sflres(iupd_aut),'(2A)')cr(1:lr),'_aut.mfr'
            WRITE(sflres(iupd_win2),'(2A)')cr(1:lr),'_win2.mfr'
          endif
        endif
        if(INDCFG.eq.2.or.INDCFG.eq.3)then
          if(unixok)then
            WRITE(splres(iupd_win1),'(3A)') '../tmp/',cr(1:lr),
     &        '_win1.plr'
            WRITE(splres(iupd_spr),'(3A)') '../tmp/',cr(1:lr),'_spr.plr'
            WRITE(splres(iupd_sum),'(3A)') '../tmp/',cr(1:lr),'_sum.plr'
            WRITE(splres(iupd_aut),'(3A)') '../tmp/',cr(1:lr),'_aut.plr'
            WRITE(splres(iupd_win2),'(3A)') '../tmp/',cr(1:lr),
     &        '_win2.plr'
          else
            WRITE(splres(iupd_win1),'(2A)')cr(1:lr),'_win1.plr'
            WRITE(splres(iupd_spr),'(2A)')cr(1:lr),'_spr.plr'
            WRITE(splres(iupd_sum),'(2A)')cr(1:lr),'_sum.plr'
            WRITE(splres(iupd_aut),'(2A)')cr(1:lr),'_aut.plr'
            WRITE(splres(iupd_win2),'(2A)')cr(1:lr),'_win2.plr'
          endif
        endif
        if(MSTROK)then
          if(unixok)then
            WRITE(smstres(iupd_win1),'(3A)') '../tmp/',cr(1:lr),
     &        '_win1.msr'
            WRITE(smstres(iupd_spr),'(3A)') '../tmp/',cr(1:lr),
     &        '_spr.msr'
            WRITE(smstres(iupd_sum),'(3A)') '../tmp/',cr(1:lr),
     &        '_sum.msr'
            WRITE(smstres(iupd_aut),'(3A)') '../tmp/',cr(1:lr),
     &        '_aut.msr'
            WRITE(smstres(iupd_win2),'(3A)') '../tmp/',cr(1:lr),
     &        '_win2.msr'
          else
            WRITE(smstres(iupd_win1),'(2A)')cr(1:lr),'_win1.msr'
            WRITE(smstres(iupd_spr),'(2A)')cr(1:lr),'_spr.msr'
            WRITE(smstres(iupd_sum),'(2A)')cr(1:lr),'_sum.msr'
            WRITE(smstres(iupd_aut),'(2A)')cr(1:lr),'_aut.msr'
            WRITE(smstres(iupd_win2),'(2A)')cr(1:lr),'_win2.msr'
          endif
        endif
        if(ientxist.gt.0)then
          if(unixok)then
            WRITE(selres(iupd_win1),'(3A)') '../tmp/',cr(1:lr),
     &        '_win1.elr'
            WRITE(selres(iupd_spr),'(3A)') '../tmp/',cr(1:lr),'_spr.elr'
            WRITE(selres(iupd_sum),'(3A)') '../tmp/',cr(1:lr),'_sum.elr'
            WRITE(selres(iupd_aut),'(3A)') '../tmp/',cr(1:lr),'_aut.elr'
            WRITE(selres(iupd_win2),'(3A)') '../tmp/',cr(1:lr),
     &        '_win2.elr'
          else
            WRITE(selres(iupd_win1),'(2A)')cr(1:lr),'_win1.elr'
            WRITE(selres(iupd_spr),'(2A)')cr(1:lr),'_spr.elr'
            WRITE(selres(iupd_sum),'(2A)')cr(1:lr),'_sum.elr'
            WRITE(selres(iupd_aut),'(2A)')cr(1:lr),'_aut.elr'
            WRITE(selres(iupd_win2),'(2A)')cr(1:lr),'_win2.elr'
          endif
        endif
        if(CFDOK)then
          if(unixok)then
            WRITE(scfdres(iupd_win1),'(3A)') '../tmp/',cr(1:lr),
     &        '_win1.dfr'
            WRITE(scfdres(iupd_spr),'(3A)') '../tmp/',cr(1:lr),
     &        '_trn.dfr'
            WRITE(scfdres(iupd_sum),'(3A)') '../tmp/',cr(1:lr),
     &        '_sum.dfr'
            WRITE(scfdres(iupd_aut),'(3A)') '../tmp/',cr(1:lr),
     &        '_aut.dfr'
            WRITE(scfdres(iupd_win2),'(3A)') '../tmp/',cr(1:lr),
     &        '_win2.dfr'
          else
            WRITE(scfdres(iupd_win1),'(2A)')cr(1:lr),'_win1.dfr'
            WRITE(scfdres(iupd_spr),'(2A)')cr(1:lr),'_trn.dfr'
            WRITE(scfdres(iupd_sum),'(2A)')cr(1:lr),'_sum.dfr'
            WRITE(scfdres(iupd_aut),'(2A)')cr(1:lr),'_aut.dfr'
            WRITE(scfdres(iupd_win2),'(2A)')cr(1:lr),'_win2.dfr'
          endif
        endif
      endif
      if(lnblnk(lipvdatf).eq.0)then
      elseif(lipvdatf(1:7).eq.'UNKNOWN')then
      else
        WRITE(sipvres,'(2A)')cr(1:lr),'ipv.rep'
      endif
      return
      end

C ******* shift_sps *******
C Shift existing SPS upwards to make room for IPV slots.
C ifrom is the lower range and ito is one more than nsset.

      subroutine shift_sps(ifrom,ito,ier)
#include "building.h"

      common/spfldes/spfdescr(MSPS)
      character spfdescr*30
      common/spflper/isstday(MSPS),isstmon(MSPS),isfnday(MSPS),
     &               isfnmon(MSPS)
      common/spfldat/nsset,isset,isstup,isbnstep,ispnstep,issave,isavgh
      INTEGER :: nsset,isset,isstup,isbnstep,ispnstep,issave,isavgh
      common/spflres/sblres(MSPS),sflres(MSPS),splres(MSPS),
     &  smstres(MSPS),selres(MSPS),scfdres(MSPS),sipvres
      character sblres*72,sflres*72,splres*72,smstres*72,
     &  selres*72,scfdres*72,sipvres*72

C Extended simulation parameters for each set.
      common/spfldats/isstupex(MSPS),isbnstepex(MSPS),ispnstepex(MSPS),
     &  issaveex(MSPS),isavghex(MSPS),iscfdactivate(MSPS),
     &  isicfdys(MSPS),isicfdyf(MSPS),
     &  scftims(MSPS),scftimf(MSPS)
      INTEGER :: isstupex,isbnstepex,ispnstepex,issaveex,isavghex
      INTEGER :: iscfdactivate      ! zero ignore domains
      INTEGER :: isicfdys,isicfdyf  ! CFD simulation start & finish days
      REAL :: scftims,scftimf       ! CFD simulation start & finish time

C If nsset is less than ito then increment nsset.  Loop downwards from
C ito until one more than ifrom.
      if(ito.le.MSPS)then
        continue
      else
        call usrmsg('Unable to add more simulation parameter sets.',
     &              'Please adapt model.','W')
        ier=1
        return
      endif
      lower=ifrom+1
      do ij=ito,lower,-1
        spfdescr(ij)=spfdescr(ij-1)
        isstday(ij)=isstday(ij-1)
        isstmon(ij)=isstmon(ij-1)
        isfnday(ij)=isfnday(ij-1)
        isfnmon(ij)=isfnmon(ij-1)
        sblres(ij)=sblres(ij-1)
        sflres(ij)=sflres(ij-1)
        splres(ij)=splres(ij-1)
        smstres(ij)=smstres(ij-1)
        selres(ij)=selres(ij-1)
        scfdres(ij)=scfdres(ij-1)
        isstupex(ij)=isstupex(ij-1)
        isbnstepex(ij)=isbnstepex(ij-1)
        ispnstepex(ij)=ispnstepex(ij-1)
        issaveex(ij)=issaveex(ij-1)
        isavghex(ij)=isavghex(ij-1)
        iscfdactivate(ij)=iscfdactivate(ij-1)
        isicfdys(ij)=isicfdys(ij-1)
        isicfdyf(ij)=isicfdyf(ij-1)
        scftims(ij)=scftims(ij-1)
        scftimf(ij)=scftimf(ij-1)
      enddo
      if(ito.gt.nsset) nsset=ito   ! Update number of SPS.
      return
      end


C ******* getmultip *******
C getmultip looks for typical week in each season based on closest
C degree days and radiation patterns. Note: adapted from DDRADSUM
C in clmsyn.F  Uses clmddscan to establish the initial ratios.

      SUBROUTINE getmultip(act)
#include "building.h"
#include "model.h"
#include "esprdbfile.h"
#include "seasons.h"
#include "ipvdata.h"
#include "help.h"

      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      COMMON/SET1/IYEAR,IBDOY,IEDOY,IFDAY,IFTIME

C Redirected text/graphics parameters.
      character xfile*144,tg*1,delim*1
      common/exporttg/xfile,tg,delim
      common/exporttgi/ixopen,ixloc,ixunit

C dm* variables are from ipvdata.h.
C sea* variables are from season.h

C Variables for handling climate file.
      character outs*124
      logical unixok,isfortnight
      character fs*1
      character t72*72
      character act*1  ! act = a for automatic (minimal interaction)
      logical MY       ! to signal not-muilti-year weather file.

      real HBT,CBT     ! local base temperatures for heating and cooling DD calcs

C For each of the seasons:
C  cddratio is the season / assessment ratio for cooling
C  hddratio is the season / assessment ratio for heating
      dimension cddratio(MIPVA)
      dimension hddratio(MIPVA)

      helpinsub='edipv'  ! set for subroutine

C Set folder separator (fs) to \ or / as required.
      call isunix(unixok)
      if(unixok)then
        fs = char(47)
      else
        fs = char(92)
      endif

      if(is1wins.eq.0.or.is2wins.eq.0.or.is1sprs.eq.0)then
        call usrmsg('No winter|transition|summer season definitions',
     &    'found in the climate database. Skipping.','W')
        return
      endif

C Initial weightings for heating dd, cooling dd, solar radiation (from seasons.h).
      hddw=1.0; cddw=1.0; radw=1.0

C Set delimiter in reports to standard/literal format.
      delim = '-'

C First check if climate datase exists: if so open.
      CALL ERPFREE(ICLIM,ISTAT)
      MY=.false.
      call CLMOPB(MY,0,ISTAT)
      if(ISTAT.ge.0)then
        CALL ERPFREE(ICLIM,ISTAT)
        MY=.false.
        call CLMOPB(MY,1,IER)
        call edisp(iuout,' opened climate file')
      endif

C Scan the climate file for CLAT and CLONG details.
      call CLMRDBMD(IER)

C Help text for the dialogs
      helptopic='seasonal_climate_scan'
      call gethelptext(helpinsub,helptopic,nbhelp)

      if(act.eq.'a')then
        hddbaset=17.0
      else
        if(hddbaset.lt.0.1)then
          hddbaset=17.0
        endif
        HBT=hddbaset
        CALL EASKR(HBT,' ',' Heating base temperature ? ',
     &    -10.,'W',40.,'W',17.0,'DD heating base temp',IER,nbhelp)
        hddbaset=HBT
      endif

      if(act.eq.'a')then
        cddbaset=21.0
      else
        if(cddbaset.lt.0.1)then
          cddbaset=21.0
        endif
        CBT=cddbaset
        CALL EASKR(CBT,' ',' Cooling base temperature ? ',
     &    -10.,'W',40.,'W',21.0,'DD cooling base temp',IER,nbhelp)
        cddbaset=CBT
      endif

C If output to file alter the edisp unit number.
      itru = iuout
      if(ixopen.eq.1)then
        itru = ixunit
        call usrmsg(' Output being directed to file... ',' ','-')
      endif

C Write out beginning of report.
      WRITE(outs,12) hddbaset,cddbaset
  12  FORMAT(' Degree day analysis: heating base at',F6.1,' & cooling',
     &       F6.1,' Deg C')
      call edisp(itru,outs)
      call edisp(itru,'  ')
      
C Loop through the climate file and calculate data. Do not bother
C to update the variables for typical weeks in climatelist.
      itrc=2
      call clmddscan(itrc,itru,'s','p',5,iyear)

C Recover the assessment periods from clmddscan.
      ipvastjd(1)=ia1wins; ipvafnjd(1)=ia1winf
      ipvastjd(2)=ia1sprs; ipvafnjd(2)=ia1sprf
      ipvastjd(3)=iasums; ipvafnjd(3)=iasumf
      ipvastjd(4)=ia2sprs; ipvafnjd(4)=ia2sprf
      ipvastjd(5)=ia2wins; ipvafnjd(5)=ia2winf

C If there are three assessments used in the IPV, advise the user of
C alternative ratios which combine both winters and both transition
C seasons.
      if(nipvassmt.eq.3)then
        call edisp(itru,'  ')
        call edisp(itru,
     &  'A 3 season IPV which adapts ratios to reflect whole year...')

C Check if a week or fortnight applicable.
        if(ipvsimu(1:6).eq.'icwinf'.or.ipvsimu(1:6).eq.'icsprf'.or.
     &     ipvsimu(1:6).eq.'icsumf'.or.ipvsimu(1:6).eq.'icautf'.or.
     &     ipvsimu(1:3).eq.'i3f'.or.ipvsimu(1:3).eq.'i5f')then
          isfortnight=.true.
        endif

C Convert from julian day of the start of the assessment to the end
C of the assessment to number of days in period and get ratio with
C the number of days in the two seasons.
        ijulst=iwkbstrt(1)
        ijulfn=iwkbstrt(1)+6
        if(isfortnight)then
          ijulfnst=ifortbstrt(1)
          ijulfnfn=ifortbstrt(1)+13
        endif
        if(isfortnight)then
          fortdayratio=(((is1winf-is1wins)+1)+((is2winf-is2wins)+1))/
     &      ((ijulfnfn-ijulfnst)+1)
          call edisp(iuout,' ')
          write(outs,'(2a,f6.3)') 
     &      'The ratio between the winter fortnight and all ',
     &      'winter weeks is',fortdayratio
          call edisp(iuout,outs)
        else
          dayratio=(((is1winf-is1wins)+1)+((is2winf-is2wins)+1))/
     &           ((ijulfn-ijulst)+1)
          call edisp(iuout,' ')
          write(outs,'(2a,f6.3)') 
     &      'The ratio between the winter week and all',
     &      'winter weeks is',dayratio
          call edisp(iuout,outs)
        endif

C Ratio of first and second winter total heating and cooling degree
C days and the best first winter assessment week degree days.
        if(wkheatdd(1).gt.0.01)then
          hddratio(1)=(seahddtot(1)+seahddtot(5))/wkheatdd(1)
        else
          hddratio(1)=1.0
        endif
        if(wkcooldd(1).gt.0.01)then
          cddratio(1)=(seacddtot(1)+seacddtot(5))/wkcooldd(1)
        else
          cddratio(1)=1.0
        endif
        if(isfortnight)then
          if(fortheatdd(1).gt.0.01)then
            forthddratio(1)=(seahddtot(1)+seahddtot(5))/fortheatdd(1)
          else
            forthddratio(1)=1.0
          endif
          if(fortcooldd(1).gt.0.01)then
            fortcddratio(1)=(seacddtot(1)+seacddtot(5))/fortcooldd(1)
          else
            fortcddratio(1)=1.0
          endif
        endif
        if(isfortnight)then
          write(outs,'(a,F5.2,a,F5.2,a,F5.2,a)')'is ',fortdayratio,
     &    ' days & ',forthddratio(1),' HDD & ',fortcddratio(1),' CDD.'
        else
          write(outs,'(a,F5.2,a,F5.2,a,F5.2,a)')'is ',dayratio,
     &    ' days & ',hddratio(1),' HDD & ',cddratio(1),' CDD.'
        endif
        if(act.eq.'a')then
          iok=1
        else
          call easkmbox(
     &      'The ratio of the winter assessment to all winter seasons',
     &      outs,'accept','edit',' ',' ',' ',' ',' ',' ',iok,nbhelp)
        endif
        if(iok.eq.1)then
          if(isfortnight)then
            dmheat(1)=forthddratio(1)
            dmcool(1)=fortcddratio(1)
            dmtime(1)=fortdayratio
          else
            dmheat(1)=hddratio(1)
            dmcool(1)=cddratio(1)
            dmtime(1)=dayratio
          endif
        elseif(iok.eq.2)then
          if(isfortnight)then
            write(t72,'(3f8.3)') forthddratio(1),fortcddratio(1),
     &        fortdayratio
          else
            write(t72,'(3f8.3)') hddratio(1),cddratio(1),dayratio
          endif
          CALL EASKS(t72,
     &      'Scaling winter heating cooling time-based',
     &      ' ',72,'scaling','scaling',IER,nbhelp)
          K=0
          CALL EGETWR(t72,K,dmheat(1),1.,999.,'W','win ht',IIER)
          CALL EGETWR(t72,K,dmcool(1),1.,999.,'W','win cl',IIER)
          CALL EGETWR(t72,K,dmtime(1),1.,999.,'W','win time',IIER)
        endif

C Do the transitional seasons.
        ijulst=iwkbstrt(2)
        ijulfn=iwkbstrt(2)+6
        if(isfortnight)then
          ijulfnst=ifortbstrt(1)
          ijulfnfn=ifortbstrt(1)+13
        endif
        dayratio=(((is1sprf-is1sprs)+1)+((is2sprf-is2sprs)+1))/
     &           ((ijulfn-ijulst)+1)
        fortdayratio=(((is1sprf-is1sprs)+1)+((is2sprf-is2sprs)+1))/
     &           ((ijulfnfn-ijulfnst)+1)
        call edisp(iuout,' ')
        if(isfortnight)then
          write(outs,'(2a,f6.3)') 
     &      'The ratio between the spring fortnight and all ',
     &      'spring and autmun weeks is',fortdayratio
          call edisp(iuout,outs)
        else
          write(outs,'(2a,f6.3)')
     &     'The ratio between the spring week and all ',
     &     'spring and autumn weeks is',dayratio
          call edisp(iuout,outs)
        endif
        if(wkheatdd(2).gt.0.01)then
          hddratio(2)=(seahddtot(2)+seahddtot(4))/wkheatdd(2)
        else
          hddratio(2)=1.0
        endif
        if(wkcooldd(2).gt.0.01)then
          cddratio(2)=(seacddtot(2)+seacddtot(4))/wkcooldd(2)
        else
          cddratio(2)=1.0
        endif
        if(isfortnight)then
          if(fortheatdd(2).gt.0.01)then
            forthddratio(2)=(seahddtot(2)+seahddtot(4))/fortheatdd(2)
          else
            forthddratio(2)=1.0
          endif
          if(fortcooldd(2).gt.0.01)then
            fortcddratio(2)=(seacddtot(2)+seacddtot(4))/fortcooldd(2)
          else
            fortcddratio(2)=1.0
          endif
        endif
        if(isfortnight)then
          write(outs,'(a,F5.2,a,F5.2,a,F5.2,a)')'is ',fortdayratio,
     &    ' days & ',forthddratio(2),' HDD & ',fortcddratio(2),' CDD.'
        else
          write(outs,'(a,F5.2,a,F5.2,a,F5.2,a)')'is ',dayratio,
     &      ' days & ',hddratio(2),' HDD & ',cddratio(2),' CDD.'
        endif
        if(act.eq.'a')then
          iok=1
        else
          call easkmbox(
     &      'The ratio of the trans assessment to all trans seasons',
     &      outs,'accept','cancel',
     &      ' ',' ',' ',' ',' ',' ',iok,nbhelp)
        endif
        if(iok.eq.1)then
          if(isfortnight)then
            dmheat(2)=forthddratio(2)
            dmcool(2)=fortcddratio(2)
            dmtime(2)=fortdayratio
          else
            dmheat(2)=hddratio(2)
            dmcool(2)=cddratio(2)
            dmtime(2)=dayratio
          endif
        elseif(iok.eq.2)then
          if(isfortnight)then
            write(t72,'(3f8.3)') forthddratio(2),fortcddratio(2),
     &        fortdayratio
          else
            write(t72,'(3f8.3)') hddratio(2),cddratio(2),dayratio
          endif
          CALL EASKS(t72,
     &       'Scaling trans heating cooling time-based',
     &       ' ',72,'scaling','scaling',IER,nbhelp)
          K=0
          CALL EGETWR(t72,K,dmheat(2),1.,999.,'W','trn ht',IIER)
          CALL EGETWR(t72,K,dmcool(2),1.,999.,'W','trn cl',IIER)
          CALL EGETWR(t72,K,dmtime(2),1.,999.,'W','trn time',IIER)
        endif
      endif

      RETURN
      END


