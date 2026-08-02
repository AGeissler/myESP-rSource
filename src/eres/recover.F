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

C recover.F contains the following routines concerned with creating
C performance information files (PIF) which record user selections
C within the res model to support subsequent '-act recover' sessions.
C READPIF: read performance information file and gather data.
C WRITEPIF: write (append) performance information file.

C ************************ READPIF ************************
C Read performance information file and gather data.

      SUBROUTINE READPIF(iier)
#include "building.h"
#include "espriou.h"
      
      integer lnblnk  ! function definition

      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      COMMON/FILEP/IFIL
      COMMON/OUTPCH/ICOUT

      COMMON/DEFLT/IDEFLT
      common/recov01/pifltog,lpifl
      character pifltog*4,lpifl*96
      common/recov02/lname,ldisc

      COMMON/SIMPIK/ISIM,ISTADD,ID1,IM1,ID2,IM2,ISDS,ISDF,NTS,ISAVE
      COMMON/PERO/IOD1,IOM1,IOH1,IOD2,IOM2,IOH2,IODS,IODF,NOUT,IAV
      COMMON/SIMPKA/NSIM
      common/recver/izver,ipver,iever

      COMMON/ZONPIK/NZ,NZNO(MCOM)
      COMMON/IGETFLG/IOCUPF,ialstused,IROC
      COMMON/GETPIK/NGET,IGETNO(MZS,9)

      character SLABEL*32,GLABEL*24,TABLABEL*36
      COMMON/GETLABEL/SLABEL(MZS),GLABEL(MZS),TABLABEL(MZS)
      integer LNSLABEL,LNGLABEL,LNTABLABEL  ! lengths for label strings
      COMMON/LNGETLABEL/LNSLABEL(MZS),LNGLABEL(MZS),LNTABLABEL(MZS)

      common/cmftpar/WVEL,ACTL,CLO,iocut,iocuset

      COMMON/RESLIB/RFILE,PFILE,MSTRFILE,LAFRES
      character rfile*72,PFILE*72,MSTRFILE*72,LAFRES*72
      COMMON/SETNAM/RSNAME(MNRS)
      character RSNAME*40

      logical libheading ! have we printed report heading (library name etc).
      logical headingcontext ! has set or period changed so heading needs reprint
      common/libhead/libheading,headingcontext
      

C ihflag = 0 write 13h30, ihflag = 1 write 0.5625
C idhflg = 0 no day demarcations, idhflg = 1 write demarcation
C          between tabular reporting days.
C ilflag = 0 tabular labels on multi-lines, ilflag = 1 on one line
C          ilflag = 2 do not include # header lines in file.
      COMMON/GRTOOL/IHFLAG,IDHFLG,ILFLAG
      
C exporttg common defn:
C  xfile - name of file data exported to
C  tg - data type tag (tabular, graphics etc), used in subroutine: ctlexp
C  delim - delimeter to use between data (comma, space, tab etc), used 
C       by calling eddisp (as opposed to edisp).
      common/exporttg/xfile,tg,delim
      COMMON/EXPORTI/ixopen,ixunit,ixpunit

C User preferences for Monthly table.
      integer imshtr,imtcas,imtaiv,imtctl,imtsol,imunit       ! monthly topic inclusion
      common/monthlyoptions/imshtr,imtcas,imtaiv,imtctl,imtsol,imunit
      COMMON/SNAP/SNAPSH
      logical SNAPSH
      common/getmenu/menutype,igetind(65),igetflux(65)

      integer ireportunit ! 0 is default, 1 is W, 2 is kW, 3 is J
      common/repunit/ireportunit
      logical markdown
      common/markdownflag/markdown
      character hcffpattern*12    ! heat, cool, or heat+cool plus detail
      common/hlcontrol/hcffpattern(2)     

      character xfile*144,tg*1,delim*1,outs*124
      character outfile*72
      character TSLABEL*32,TLLABEL*64,TGLABEL*24
      character OUTSTR*124,lname*24,ldisc*64,key*8
      character prompt*124,mark*12
      character PDESCR*64,SDESCR*44

      logical recal,newzn
      integer iflen  ! length of file

C IOFMT - save performance groupings in same file=1, or
C separately=0 (default), using IGRP to remember which group we're 
C recovering.
      IOFMT=0
      IGRP=0

C Open performance information file.
      ier=0
      IUNIT=IFIL+13
      CALL EFOPSEQ(IUNIT,lpifl,1,IER)
      if(IER.ne.0)THEN
        IER=1
        return
      endif
      write (currentfile,'(a)') lpifl(1:lnblnk(lpifl))

C Read header line.
      CALL STRIPC(IUNIT,OUTSTR,99,ND,1,'pif header line',IER)
      if(IER.NE.0)return
      if(OUTSTR(1:29).ne.'*Performance information file')then
        call usrmsg('Incorrect file format in:',lpifl,'W')
        return
      endif

C Read output format specification. The first tokens specify the
C output file name and if 'markdown' also found then remember.
      CALL STRIPC(IUNIT,OUTSTR,99,ND,1,'output fmt',IER)
      K=0
      call EGETWI(OUTSTR,K,IOFMT,0,99,'W','output format',IER)
      outfile='UNKNOWN'
      if (ND.eq.2) then
        call EGETW(OUTSTR,K,outfile,'W','output file name',IER)
        markdown=.false.
      elseif (ND.eq.3) then
        call EGETW(OUTSTR,K,outfile,'W','output file name',IER)
        call EGETW(OUTSTR,K,mark,'W','markdown',IER)
        if(mark(1:8).eq.'markdown')then
          markdown=.true.
        else
          markdown=.false.
        endif
      endif

C Open output redirect (if saving all outputs together).
      if (IOFMT.eq.1) then
        if (outfile(1:7).eq.'UNKNOWN') then

C Output file not specified - use pif file name but changed to *.dat
          iflen=lnblnk(lpifl)-3
          if (lpifl(iflen:iflen).eq.'.') then
            write (xfile,'(2a)') lpifl(1:iflen),'dat'
          else
            write (xfile,'(2a)') lpifl(1:lnblnk(lpifl)),'.dat'
          endif
        else
          write(xfile,'(a)') outfile(1:lnblnk(outfile))
        endif
        call ctlexp(xfile,ixopen,ixunit,ixpunit,'T','Text',IER)

C If user request cancel then return after resetting ixopen to zero.
        if(ier.eq.-3)then
          ixopen=0
          return
        endif
      endif

C Allow output period and filters to be changed on subsequent groups.
C If changes detected then db is re-scanned.
 1    IDEFLT=0
      IOCUPF=0
      ialstused=0
      IROC=0
      recal=.FALSE.

C Read group data.  Return from routine when end of file reached (goto 10).
      IER=0
      CALL STRIPC(IUNIT,OUTSTR,99,ND,0,'group data',IER)
      if (IER.ne.0) goto 10
      K=0
      call EGETW(OUTSTR,K,lname,'W','group name',IER)
      call EGETRM(OUTSTR,K,ldisc,'W','group description',IER)

C Read group output type.
      CALL STRIPC(IUNIT,OUTSTR,99,ND,1,'output style',IER)
      K=0
      call EGETWI(OUTSTR,K,IOTYPE,0,13,'F','output type',IER)

C Get additional data required for output choice.
      if (IOTYPE.eq.7.or.IOTYPE.eq.8) then
        call EGETWR(OUTSTR,K,tval,0.,0.,'-','query point',IER)
      endif

C Check for keywords.
 200  CALL STRIPC(IUNIT,OUTSTR,99,ND,1,'key words',IER)
      K=0
      call EGETW(OUTSTR,K,key,'W','key name',IER)
      if (key(1:7).eq.'*period') then

C If more items on the line then recover the period information.
        if(ND.gt.7)then
          call EGETWI(OUTSTR,K,ITSD,1,31,'F','start day',IER)
          call EGETWI(OUTSTR,K,ITSM,1,12,'F','start month',IER)
          call EGETWI(OUTSTR,K,ITSH,0,24,'F','start hour',IER)
          call EGETWI(OUTSTR,K,ITFD,1,31,'F','finish day',IER)
          call EGETWI(OUTSTR,K,ITFM,1,12,'F','finish month',IER)
          call EGETWI(OUTSTR,K,ITFH,0,24,'F','finish hour',IER)
          call EGETWI(OUTSTR,K,IOTSA,1,1440,'W','output TS',IER)
          call EGETWI(OUTSTR,K,IAFLG,0,1,'W','output aver',IER)
        else

C Read period information over several lines.
          CALL STRIPC(IUNIT,OUTSTR,99,ND,1,'start date',IER)
          K=0
          call EGETWI(OUTSTR,K,ITSD,1,31,'F','start day',IER)
          call EGETWI(OUTSTR,K,ITSM,1,12,'F','start month',IER)
          call EGETWI(OUTSTR,K,ITSH,0,24,'F','start hour',IER)

          CALL STRIPC(IUNIT,OUTSTR,99,ND,1,'finish date',IER)
          K=0
          call EGETWI(OUTSTR,K,ITFD,1,31,'F','finish day',IER)
          call EGETWI(OUTSTR,K,ITFM,1,12,'F','finish month',IER)
          call EGETWI(OUTSTR,K,ITFH,0,24,'F','finish hour',IER)

C Get averaging data.
          CALL STRIPC(IUNIT,OUTSTR,99,ND,1,'averaging data',IER)
          K=0
          call EGETWI(OUTSTR,K,IOTSA,1,1440,'W','output TS',IER)
          call EGETWI(OUTSTR,K,IAFLG,0,1,'W','output aver',IER)
        endif

C Range check on start and finish days.
        recal=.FALSE.
        if (ITFM.lt.ITSM) then
          call edisp(iuout,'Start/ finish months out of order!')
        elseif (ITFM.eq.ITSM) then
          if (ITFD.lt.ITSD) then
            call edisp(iuout,'Start/ finish days out of order!')
          elseif (ITFD.eq.ITSD) then
            if (ITFH.le.ITSH) then
              call edisp(iuout,'Start/ finish hours out of order!')
            endif
          endif
        endif

C Check output period validity.
        if (ITSM.lt.IM1.or.ITFM.gt.IM2) then
          IDEFLT=1
          call MOOPER
          IDEFLT=0
        elseif (ITSM.eq.IM1.and.ITSD.lt.ID1) then
          IDEFLT=1
          call MOOPER
          IDEFLT=0
        elseif (ITFM.eq.IM2.and.ITFD.gt.ID2) then
          IDEFLT=1
          call MOOPER
          IDEFLT=0
        else
          if (ITSM.ne.IOM1) then  ! Start of output period.
            IOM1=ITSM
            recal=.TRUE.
          endif
          if (ITSD.ne.IOD1) then
            IOD1=ITSD
            recal=.TRUE.
          endif
          if (ITSH.ne.IOH1) then
            IOH1=ITSH
            recal=.TRUE.
          endif
          if (ITFM.ne.IOM2) then  ! End of output period.
            IOM2=ITFM
            recal=.TRUE.
          endif
          if (ITFD.ne.IOD2) then
            IOD2=ITFD
            recal=.TRUE.
          endif
          if (ITFH.ne.IOH2) then
            IOH2=ITFH
            recal=.TRUE.
          endif
          if (IOTSA.ne.NOUT) then  ! Averaging.
            NOUT=IOTSA
            recal=.TRUE.
          endif
          if (NOUT.gt.1.and.IAFLG.ne.IAV) then
            IAV=IAFLG
            recal=.TRUE.
          endif
        endif
        if (recal) then    ! If output period has changed.
          CALL EDAY(IOD1,IOM1,IODS)
          CALL EDAY(IOD2,IOM2,IODF)
        endif

C Filtering directives.  If data follows *filters then parse, otherwise
C read the next line.
      elseif(key(1:8).eq.'*filters') then
        if(ND.gt.6)then
          call EGETWI(OUTSTR,K,IOCUPF,0,1,'F','occup filter',IER)
          call EGETWI(OUTSTR,K,IROC,-1,1,'W','ROC filter',IER)
          call EGETWI(OUTSTR,K,ialstused,0,1,'F','anch filt',IER)
          call EGETW(OUTSTR,K,delim,'W','delimeter',IER)
          call EGETW(OUTSTR,K,hcffpattern(1),'W','hcfpattern',IER)
          call EGETWI(OUTSTR,K,ireportunit,0,3,'W','units',IER)
          call EGETWI(OUTSTR,K,IHFLAG,0,1,'W','hour flag',IER)
          call EGETWI(OUTSTR,K,IDHFLG,0,1,'W','day flag',IER)
          call EGETWI(OUTSTR,K,ILFLAG,0,2,'W','tabular flag',IER)
        else
          CALL STRIPC(IUNIT,OUTSTR,99,ND,1,'filtering data',IER)
          K=0
          call EGETWI(OUTSTR,K,IOCUPF,0,1,'F','occup filter',IER)
          call EGETWI(OUTSTR,K,IROC,-1,1,'W','ROC filter',IER)
          if(ND.gt.2)then
            call EGETWI(OUTSTR,K,ialstused,0,1,'F','anch filt',IER)
            if(ND.gt.3)then
              call EGETW(OUTSTR,K,delim,'W','delimeter',IER)
            endif
          endif
        endif

C Based on what was read set additional attributes.
        if (IOCUPF.eq.1) then
          iocut=1
          IFILT=1
          call usrmsg('Assuming casual gain type 1=occupancy','  ','-')
        endif

C Read IGETNO data.
      elseif(key(1:8).eq.'*metlist') then
        NGET=0

C Labels read into temp variables then copy below.
        TSLABEL='  '
        TLLABEL='  '
        TGLABEL='  '

C Set inital assumptions for monthly column data.
        imshtr=1; imtcas=1; imtaiv=1; imtctl=1; imtsol=1; imunit=1
 100    CALL STRIPC(IUNIT,OUTSTR,99,ND,1,'text list data',IER)
        K=0
        call EGETW(OUTSTR,K,key,'W','key name',IER)
        if(key(1:5).eq.'*data') then
          NGET=NGET+1
          call EGETWI(OUTSTR,K,IGETNO(NGET,1),0,76,'F',
     &      'iget number',IER)
          call EGETWI(OUTSTR,K,IGETNO(NGET,2),0,MCOM,'F',
     &      'zone number',IER)
          call EGETWI(OUTSTR,K,IGETNO(NGET,3),0,MZS,'W',
     &      'misc data 1',IER)
          call EGETWI(OUTSTR,K,IGETNO(NGET,4),-MZS,MZS,'W',
     &      'misc data 2',IER)
          call EGETWI(OUTSTR,K,IGETNO(NGET,5),0,MNRS,'F',
     &      'set number',IER)
          call EGETWI(OUTSTR,K,IGETNO(NGET,6),0,0,'-',
     &      'line style',IER)
          call EGETWI(OUTSTR,K,IGETNO(NGET,7),1,15,'W',
     &      'symbol type',IER)
          call EGETWI(OUTSTR,K,IGETNO(NGET,8),1,4,'F',
     &      'axis number',IER)
          if (IOTYPE.eq.5.and.ND.gt.8) then        ! For stats read label width.
            call EGETWI(OUTSTR,K,IGETNO(NGET,9),0,0,'-',
     &      'label width',IER)
          elseif (IOTYPE.eq.13.and.ND.gt.10) then   ! For monthly get column choices.
            call EGETWI(OUTSTR,K,imshtr,0,0,'-','surf ht tr col',IER)
            call EGETWI(OUTSTR,K,imtcas,0,0,'-','cas gains col',IER)
            call EGETWI(OUTSTR,K,imtaiv,0,0,'-','imtaiv columns',IER)
            call EGETWI(OUTSTR,K,imtctl,0,0,'-','imtctl columns',IER)
            call EGETWI(OUTSTR,K,imtsol,0,0,'-','solar columns',IER)
            if(ND.gt.14)then
              call EGETWI(OUTSTR,K,imunit,0,0,'-','monthly units',IER)
            endif
          endif
        elseif(key(1:5).eq.'*glab') then
          call EGETRM(OUTSTR,K,TGLABEL,'W','graph label',IER)
        elseif(key(1:5).eq.'*slab') then
          call EGETRM(OUTSTR,K,TSLABEL,'W','short label',IER)
        elseif(key(1:5).eq.'*llab') then
          call EGETRM(OUTSTR,K,TLLABEL,'W','long label',IER)
        elseif (OUTSTR(1:4).eq.'*end') then
          goto 300
        endif

C Copy label names to common area.
        SLABEL(NGET)=TSLABEL
        LNSLABEL(NGET)=lnblnk(TSLABEL)
        GLABEL(NGET)=TGLABEL
        LNGLABEL(NGET)=lnblnk(TGLABEL)
        TABLABEL(NGET)=TGLABEL
        LNTABLABEL(NGET)=lnblnk(TGLABEL)
        goto 100
      else
        goto 200
      endif
      goto 200

C Check for unique zone numbers and add to NZNO and NZ in common,
C this is needed for generating correct header information.
 300  NZ=0
      do 20 JZ=1,NGET
        newzn=.TRUE.
        if (NZ.gt.0) then
          do 30 KZ=1,NZ
            if (NZNO(KZ).eq.IGETNO(JZ,2)) newzn=.FALSE.
 30       continue
          if (newzn) then
            NZ=NZ+1
            NZNO(NZ)=IGETNO(JZ,2)
          endif
        else
          NZ=1
          NZNO(NZ)=IGETNO(JZ,2)
        endif
 20   continue

C Open output redirect (if saving separately).
      if (IOFMT.eq.0) then
        IGRP=IGRP+1
        if (outfile(1:7).eq.'UNKNOWN') then

C Output file not specified - use pif file name but changed to *.dat
          iflen=lnblnk(lpifl)-3
          if (lpifl(iflen:iflen).eq.'.') then
            write (xfile,'(2a,i4.4)') lpifl(1:iflen),'dat',IGRP
          else
            write (xfile,'(2a,i4.4)') lpifl(1:lnblnk(lpifl)),'.dat',IGRP
          endif
        else
           write (xfile,'(a,i4.4)') outfile,IGRP
        endif
        if(IOTYPE.eq.3)then
          call ctlexp(xfile,ixopen,ixunit,ixpunit,'G','Graphics',IER)
        elseif(IOTYPE.eq.4)then
          call ctlexp(xfile,ixopen,ixunit,ixpunit,'X','Tabular',IER)
        else
          call ctlexp(xfile,ixopen,ixunit,ixpunit,'T','Text',IER)
        endif

C If user request cancel then return after resetting ixopen to zero.
        if(ier.eq.-3)then
          ixopen=0
          return
        endif
      endif

C Heading string.
      if(lname(1:4).ne.'name')then
        if(delim.eq.'X')then
          write(outs,'(2a)') '*output ',lname(1:lnblnk(lname))
        else
          write(outs,'(2a)') 'output ',lname(1:lnblnk(lname))
        endif
        call edisp(ixunit,outs)
      endif
      if(ldisc(1:4).ne.'desc')then
C        write(outs,'(2a)') 'Recovering ',ldisc(1:lnblnk(ldisc))
C        call edisp(iuout,outs)
        call edisp(ixunit,ldisc)
      else
        continue
      endif

C Now recover data for this group.
C IOTYPE generally uses same definitions as MENUTYPE.
      if (IOTYPE.eq.1) then
        continue
      elseif (IOTYPE.eq.2) then

C MODGRF Needs: SLABEL.
        ivvv=0
        call MODGRF('p',ivvv)
      elseif (IOTYPE.eq.3) then

C MOWLIN Needs: GLABEL and SLABEL and graphics mode.
        call MOWLIN('p')
      elseif (IOTYPE.eq.4) then

C Set format flags for tabular listing to file (then reset to default).
C ZONTAB Needs: GLABEL and SLABEL.
        IHFLAG=1
        IDHFLG=0
        call ZONTAB('p')
        IHFLAG=0
        IDHFLG=0
      elseif (IOTYPE.eq.5) then

C MOSMSD Needs: GLABEL (SLABEL used in prompts).
        call MOSMSD('p',VMX,VMN,VAVE,VTOT,DIVMX,DIVMN,DISTMX,DISTMN)
      elseif (IOTYPE.eq.6) then
        continue
      elseif (IOTYPE.eq.7) then

C MOSQHR Needs: GLABEL (SLABEL used in prompts).
        call MOSQHR('p',1,tval)
      elseif (IOTYPE.eq.8) then

C MOSQHR Needs: GLABEL (SLABEL used in prompts).
        call MOSQHR('p',-1,tval)
      elseif (IOTYPE.eq.9) then

C Zone energy balance. Similar code as in stats.F
        SNAPSH=.FALSE.  ! over period of assessment
        IWHICH=2        ! gains and losses
        itru = icout
        prompt = ' '
        if(ixopen.eq.1)then
          itru = ixunit
        endif
        MENUTYPE=9
        call GOMSETUP
        call GOMENU
        if(ixopen.eq.1)then
          itru = ixunit
          if(NGET.ge.1)then
            write(prompt,'(a,a)') SLABEL(1)(1:LNSLABEL(1)),' >> file.'
          else
            write(prompt,'(a)')' Output being directed to file. ' 
          endif
          call usrmsg(prompt,'Scanning for range of values...','-')
        else
          write(prompt,'(3a)') 'Scanning for ',
     &      SLABEL(1)(1:LNSLABEL(1)),'...'
          itru = icout
          call usrmsg(prompt,'  ','-')
        endif

        ISET=ISIM
        CALL HDDATE(PDESCR)
        CALL HDSTEP(SDESCR)
        lsn1=MIN0(lnblnk(RFILE),42)
        if((.NOT.libheading).and.(.NOT.headingcontext))then
          libheading=.true.; headingcontext=.true.
          if(NSIM.gt.1)then
          write(outs,'(3A,I4,3A)')'Results library: ',
     &      RFILE(1:lsn1),'; results set:',ISIM,
     &      ' (',RSNAME(ISIM)(1:lnblnk(RSNAME(ISIM))),')'
          else
            write(outs,'(5A)')'Results library: ',RFILE(1:lsn1),
     &      '; (',RSNAME(ISIM)(1:lnblnk(RSNAME(ISIM))),')'
          endif
          call edisp(itru,outs)
          write(outs,'(3A)')PDESCR(1:lnblnk(PDESCR)),' ',
     &     SDESCR(12:lnblnk(SDESCR))
          call edisp(itru,outs)
        endif

C For each selected zone generate energy balance report.
        DO 65 JZ=1,NZ
          izone=NZNO(JZ)
          if(ISAVE.EQ.4)then
            CALL MOZCAS(IZONE,IWHICH,NIT)
          elseif(ISAVE.eq.2.or.ISAVE.eq.3)then
            if(izver.ge.2)then
              CALL MOZCAS(izone,IWHICH,NIT)
            else
              call usrmsg(
     &          'The results library is an older format and a full',
     &          'zone energy balance cannot be recovered.','W')
              CALL MOZCAS(izone,IWHICH,NIT)
            endif
          endif
 65     continue
        continue
      elseif (IOTYPE.eq.10) then

C MENUTYPE 10:  Energy demands.
        call MOEDMDS('p',the,tzhtime,tce,tzctime,tlhe,tzlhtime,
     &    tlce,tzlctime)
      elseif (IOTYPE.eq.11) then

C MENUTYPE 11: Zone casual gains distribution.
        call CASDIST('p')
      elseif (IOTYPE.eq.12) then

C MENUTYPE 12: Surface energy balance flux paths. << to be included >>
        continue
      elseif (IOTYPE.eq.13) then

C Monthly gains & losses.
        CALL MGLTBL
      elseif (IOTYPE.eq.14) then

C Surface condensation. << to be included >>
        continue
      endif

C Close output redirect (if saving seperately).
      if (IOFMT.eq.0) then
        if(IOTYPE.eq.3)then
          call ctlexp(xfile,ixopen,ixunit,ixpunit,'G','Graphics',IER)
        elseif(IOTYPE.eq.4)then
          call ctlexp(xfile,ixopen,ixunit,ixpunit,'X','Tabular',IER)
        else
          call ctlexp(xfile,ixopen,ixunit,ixpunit,'T','Text',IER)
        endif

C If user request cancel then return after resetting ixopen to zero.
        if(ier.eq.-3)then
          ixopen=0
          return
        endif
      endif

C Read next group.
      goto 1

C Close output redirect (if saving all outputs together).
 10   if (IOFMT.eq.1) then
        call ctlexp(xfile,ixopen,ixunit,ixpunit,'T','Text',IER)
      endif

C Close PIF file.
      call ERPFREE(IUNIT,ISTAT)

      return
      end


C ************************ WRITEPIF ************************
C WRITEPIF - write out performance information file based on
C current menu (moget) selections.

      SUBROUTINE WRITEPIF(IOTYPE)
#include "building.h"
#include "model.h"
      
      integer lnblnk  ! function definition

      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      COMMON/FILEP/IFIL
      common/recov01/pifltog,lpifl
      character pifltog*4,lpifl*96
      common/recov02/lname,ldisc

      COMMON/PERO/IOD1,IOM1,IOH1,IOD2,IOM2,IOH2,IODS,IODF,NOUT,IAV

      COMMON/IGETFLG/IOCUPF,ialstused,IROC
      COMMON/GETPIK/NGET,IGETNO(MZS,9)

      character SLABEL*32,GLABEL*24,TABLABEL*36
      COMMON/GETLABEL/SLABEL(MZS),GLABEL(MZS),TABLABEL(MZS)
      integer LNSLABEL,LNGLABEL,LNTABLABEL  ! lengths for label strings
      COMMON/LNGETLABEL/LNSLABEL(MZS),LNGLABEL(MZS),LNTABLABEL(MZS)

      COMMON/GRTOOL/IHFLAG,IDHFLG,ILFLAG
      integer imshtr,imtcas,imtaiv,imtctl,imtsol,imunit       ! monthly topic inclusion
      common/monthlyoptions/imshtr,imtcas,imtaiv,imtctl,imtsol,imunit

      integer ireportunit ! 0 is default, 1 is W, 2 is kW, 3 is J
      common/repunit/ireportunit
      logical markdown
      common/markdownflag/markdown
      character hcffpattern*12    ! heat, cool, or heat+cool plus detail
      common/hlcontrol/hcffpattern(2)     

      common/exporttg/xfile,tg,delim

      character xfile*144,tg*1,delim*1
      character lname*24,ldisc*64
      character pifexport*96

      logical XST

C Use temp file unit.
      IUNIT=IFIL+1
      XST=.false.

      INQUIRE(FILE=lpifl,EXIST=XST)

C New file write first few lines else open in append mode.
C The name of the file res will export to when the PIF
C PIF file is used in a subsequent run is assumed to
C be the model root name + .dat.
      if (.NOT.XST) then
        call ERPFREE(IUNIT,ISTAT)
        lpf=lnblnk(lpifl)
        open(IUNIT,file=lpifl(1:lpf),status='UNKNOWN',err=2)
        write(IUNIT,'(a)') '*Performance information file'
        write(pifexport,'(2a)') cfgroot(1:lnblnk(cfgroot)),'.dat'
        if(markdown)then
          write(IUNIT,'(4a)') '1 ',pifexport(1:lnblnk(pifexport)),
     &    ' markdown  # Output file in markdown format'
        else
          write(IUNIT,'(3a)') '1 ',pifexport(1:lnblnk(pifexport)),
     &    ' # Output file'
        endif
      else
        close(IUNIT)
        lpf=lnblnk(lpifl)
        open(IUNIT,file=lpifl(1:lpf),position='APPEND',status='UNKNOWN',
     &        err=1)
      endif

C Write current user topic selection to the PIF file.
      write(IUNIT,'(a)') '# Group name and description.'
      write(IUNIT,'(a)') 'name_here description_here'
      if(IOTYPE.eq.2)then
        write(IUNIT,'(i3,a)') IOTYPE,' # histogram'
      elseif(IOTYPE.eq.3)then
        write(IUNIT,'(i3,a)') IOTYPE,' # line graph'
      elseif(IOTYPE.eq.4)then
        write(IUNIT,'(i3,a)') IOTYPE,' # tabular listing'
      elseif(IOTYPE.eq.5)then
        write(IUNIT,'(i3,a)') IOTYPE,' # statistics'
      elseif(IOTYPE.eq.7)then
        write(IUNIT,'(i3,a)') IOTYPE,' # hours greater'
      elseif(IOTYPE.eq.8)then
        write(IUNIT,'(i3,a)') IOTYPE,' # hours less than'
      elseif(IOTYPE.eq.9)then
        write(IUNIT,'(i3,a)') IOTYPE,' # zone energy balance'
      elseif(IOTYPE.eq.10)then
        write(IUNIT,'(i3,a)') IOTYPE,' # energy demands'
      elseif(IOTYPE.eq.11)then
        write(IUNIT,'(i3,a)') IOTYPE,' # casual gains'
      elseif(IOTYPE.eq.12)then
        write(IUNIT,'(i3,a)') IOTYPE,' # surface energy balance'
      elseif(IOTYPE.eq.13)then
        write(IUNIT,'(i3,a)') IOTYPE,' # monthly reporting'
      else
        write(IUNIT,'(i3,a)') IOTYPE,' # style'
      endif
      write(IUNIT,'(2a)') '# Start  day month  hour finish day',
     &  ' month hour output-ts averaging'
      write(IUNIT,'(a,8i6)') '*period',IOD1,IOM1,IOH1,IOD2,IOM2,
     &  IOH2,NOUT,IAV
      
      write(IUNIT,'(3a)') '# directives: occupancy',
     &  ' rate-of-change anchors delimeter ctl-patterns units hours',
     &  ' day-demarcation tabular-labels'
      write(IUNIT,'(a,3i6,1x,3a,4i6)') '*filters',IOCUPF,IROC,
     &   ialstused,delim,' ',hcffpattern(1)(1:lnblnk(hcffpattern(1))),
     &   ireportunit,IHFLAG,IDHFLG,ILFLAG
      write(IUNIT,'(a)') '*metlist'
      if(IOTYPE.eq.5)then
        write (IUNIT,'(2a)') '# metric zone surf/type node/layer',
     &    ' res-set line-type symbol axis label-width'
      elseif(IOTYPE.eq.13)then
        write (IUNIT,'(2a)') '# metric zone surf/type node/layer',
     &    ' res-set line-type symbol axis monthly-column-choices'
      else
        write (IUNIT,'(2a)') '# metric zone surf/type node/layer',
     &    ' res-set line-type symbol axis'
      endif
      do 10 I=1,NGET
        if(IOTYPE.eq.5)then
          write (IUNIT,'(a,9i6)') '*data ',(IGETNO(I,J),J=1,9)
        elseif(IOTYPE.eq.13)then
          write (IUNIT,'(a,8i6,6i6)') '*data ',(IGETNO(I,J),J=1,8),
     &      imshtr,imtcas,imtaiv,imtctl,imtsol,imunit
        else
          write (IUNIT,'(a,8i6)') '*data ',(IGETNO(I,J),J=1,8)
        endif
        write (IUNIT,'(2a)') '*slab ',SLABEL(I)(1:LNSLABEL(I))
        write (IUNIT,'(2a)') '*glab ',GLABEL(I)(1:LNGLABEL(I))
 10   continue
      write (IUNIT,'(a)') '*end '
      write (IUNIT,'(a)') '# --------------------------------- '

      return
   1  call edisp(iuout,'Error re-opening recovery file, continuing.')
      return
   2  call edisp(iuout,'Error opening recovery file, continuing.')
      return
      end
