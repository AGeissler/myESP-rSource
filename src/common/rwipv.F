C Subroutines to read and write IPV descriptions.

C *********** mkipvdat
C Write IPV data to an external file or within cfg file based
C on the parameter act (3 char) ='ipv' or ='cfg'
      subroutine mkipvdat(ifu,file,act)
#include "building.h"
#include "ipvdata.h"

      integer lnblnk  ! function definition
      
      integer ipvictmindex    ! which contaminate source
      character ipvctmname*12 ! name of contaminate
      character ipvndname*12  ! associated flow node name
      common/contm10/ipvictmindex,ipvctmname,ipvndname

      real vcp  ! Guth visual comfort probablility -60 - 60 values from Radiance.
      common/guth/vcp(13)

      character file*72,tab*1
      character outs*124,outsd*124,outsn*124,louts*248
      character act*3   ! indicate which file destination
      character tokens*156,comment*76,aligned_str*156
      character decode*36
      character dstmp*24
      integer list
      dimension list(MCOM)  ! zone lists currently set to MCOM
      integer loutlen,itrunc,ipos  ! used with ailist call

C      tab=CHAR(9)
      tab=','

C Compose comment for opening line.
      if(ipvsimu(1:6).eq.'icwint'.or.ipvsimu(1:6).eq.'Icwint')then
        decode='winter typical week'
      elseif(ipvsimu(1:6).eq.'icwinf'.or.ipvsimu(1:6).eq.'Icwinf')then
        decode='winter typical fortnight'
      elseif(ipvsimu(1:6).eq.'icwins'.or.ipvsimu(1:6).eq.'Icwins')then
        decode='winter full season '
      elseif(ipvsimu(1:6).eq.'icsprt'.or.ipvsimu(1:6).eq.'Icsprt')then
        decode='spring typical week'
      elseif(ipvsimu(1:6).eq.'icsprf'.or.ipvsimu(1:6).eq.'Icsprf')then
        decode='spring typical fortnight'
      elseif(ipvsimu(1:6).eq.'icsprs'.or.ipvsimu(1:6).eq.'Icsprs')then
        decode='spring full season'
      elseif(ipvsimu(1:6).eq.'icsumt'.or.ipvsimu(1:6).eq.'Icsumt')then
        decode='summer typical week'
      elseif(ipvsimu(1:6).eq.'icsumf'.or.ipvsimu(1:6).eq.'Icsumf')then
        decode='summer typical fortnight'
      elseif(ipvsimu(1:6).eq.'icsums'.or.ipvsimu(1:6).eq.'Icsums')then
        decode='summer full season'
      elseif(ipvsimu(1:6).eq.'icautt'.or.ipvsimu(1:6).eq.'Icautt')then
        decode='autumn typical week'
      elseif(ipvsimu(1:6).eq.'icautf'.or.ipvsimu(1:6).eq.'Icautf')then
        decode='autumn typical fortnight'
      elseif(ipvsimu(1:6).eq.'icauts'.or.ipvsimu(1:6).eq.'Icauts')then
        decode='autumn full season'
      elseif(ipvsimu(1:3).eq.'ias'.or.ipvsimu(1:3).eq.'Ias')then
        decode='annual'
      elseif(ipvsimu(1:3).eq.'i1d'.or.ipvsimu(1:3).eq.'I1d')then
        decode='single default'
      elseif(ipvsimu(1:3).eq.'i3s'.or.ipvsimu(1:3).eq.'I3s')then
        decode='three full seasons'
      elseif(ipvsimu(1:3).eq.'i3t'.or.ipvsimu(1:3).eq.'I3t')then
        decode='three seasons (typical weeks)'
      elseif(ipvsimu(1:3).eq.'i3f'.or.ipvsimu(1:3).eq.'I3f')then
        decode='three seasons (typical fortnight)'
      elseif(ipvsimu(1:3).eq.'i5s'.or.ipvsimu(1:3).eq.'I5s')then
        decode='five full seasons'
      elseif(ipvsimu(1:3).eq.'i5t'.or.ipvsimu(1:3).eq.'I5t')then
        decode='five seasons (typical weeks)'
      elseif(ipvsimu(1:3).eq.'i5f'.or.ipvsimu(1:3).eq.'I5f')then
        decode='five seasons (typical fortnight)'
      endif

C Open any existing file by this name or create a new file.
      if(act(1:3).eq.'ipv')then
        CALL  ERPFREE(ifu,ISTAT)
        CALL EFOPSEQ(ifu,file,3,IER)
        IF(IER.LT.0)THEN
          IER=1
          RETURN
        ENDIF
        ipvversion=5
        write(ifu,'(a,i1)') '*IPV  V',ipvversion
        call dstamp(dstmp)
        write(ifu,'(2a)',IOSTAT=IOS,ERR=1003) '*date ',dstmp
        if(ipvform.eq.0)then
          ipvform=2       ! Change to csv if not set.
        endif
        if(ipvform.eq.1)then
          continue        ! If normal format skip
        else
          write(ifu,'(a)') '*report_format csv'
        endif
      elseif(act(1:3).eq.'cfg')then
        ipvversion=5
        write(tokens,'(a,i1,a,i1,2a)') '*ipv  ',ipvversion,
     &    ' ',ipvform,' ',ipvsimu
        call align_comment(32,tokens,comment,aligned_str)
        write(ifu,'(a)') aligned_str(1:lnblnk(aligned_str))
      endif
      WRITE(ifu,'(2A)',IOSTAT=IOS,ERR=1003)'*title ',
     &  ipvtitl(1:lnblnk(ipvtitl))

C Note: ipvsynop might be a long line.
      WRITE(ifu,'(2A)',IOSTAT=IOS,ERR=1003)'*synop ',
     &  ipvsynop(1:lnblnk(ipvsynop))

C If user nominated a focus zone include its name.
      if(ivpfocname(1:2).eq.'  '.or.ivpfocname(1:4).eq.'UNKN')then
        continue
      else
        lnipvz=lnblnk(ivpfocname)
        WRITE(ifu,'(2A)',IOSTAT=IOS,ERR=1003)'*focus ',
     &    ivpfocname(1:lnipvz)
      endif

C There are two likely images - one for building and one for glare.
      if(nipvimg.gt.0)then
        if(lipvimg(1)(1:2).eq.'  '.or.lipvimg(1)(1:4).eq.'UNKN')then
          continue
        else
          WRITE(outsn,'(2A)',IOSTAT=IOS,ERR=1003) '*building_image ',
     &      lipvimg(1)(1:lnblnk(lipvimg(1)))
          WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003) outsn(1:lnblnk(outsn))
        endif
        if(lipvimg(2)(1:2).eq.'  '.or.lipvimg(2)(1:4).eq.'UNKN')then
          continue
        else
          WRITE(outsn,'(2A)',IOSTAT=IOS,ERR=1003) '*glare_image ',
     &      lipvimg(2)(1:lnblnk(lipvimg(2)))
          WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003) outsn(1:lnblnk(outsn))
        endif
      endif

C Loop through each assessment period.
      WRITE(tokens,'(A,i3,2a)',IOSTAT=IOS,ERR=1003)'*simulations',
     &  nipvassmt,' ',ipvsimu
      write(comment,'(a)') decode(1:lnblnk(decode))
      call align_comment(32,tokens,comment,aligned_str)
      write(ifu,'(a)') aligned_str(1:lnblnk(aligned_str))
      WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003) '# ' 
      WRITE(ifu,'(2A)',IOSTAT=IOS,ERR=1003)
     & '# season start  finish focus heating cooling time-based',
     & '  description'
      WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003)
     & '# season  day    day    day  scaling scaling scaling '

      if(ipvsimu(1:6).eq.'icwint'.or.ipvsimu(1:6).eq.'icwinf'.or.
     &   ipvsimu(1:6).eq.'icwins')then     ! Only one winter assessment
        WRITE(outsn,'(A,i4,2i7,3F8.3,2a)',IOSTAT=IOS,ERR=1003)
     &    '*winter ',
     &    ipvastjd(1),ipvafnjd(1),ipvdispjd(1),ddmheat(1),
     &    ddmcool(1),ddmtime(1),'    ',ipvadesc(1)
        WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003) outsn(1:lnblnk(outsn))
      elseif(ipvsimu(1:6).eq.'icsprt'.or.ipvsimu(1:6).eq.'icsprf'.or.
     &       ipvsimu(1:6).eq.'icsprs')then    
        WRITE(outsn,'(A,i4,2i7,3F8.3,2a)',IOSTAT=IOS,ERR=1003)
     &    '*spring  ',
     &    ipvastjd(2),ipvafnjd(2),ipvdispjd(2),ddmheat(2),
     &    ddmcool(2),ddmtime(2),'    ',ipvadesc(2)
        WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003) outsn(1:lnblnk(outsn))
      elseif(ipvsimu(1:6).eq.'icsumt'.or.ipvsimu(1:6).eq.'icsumf'.or.
     &       ipvsimu(1:6).eq.'icsums')then
        WRITE(outsn,'(A,i4,2i7,3F8.3,2a)',IOSTAT=IOS,ERR=1003)
     &    '*summer  ',
     &    ipvastjd(3),ipvafnjd(3),ipvdispjd(3),ddmheat(3),
     &    ddmcool(3),ddmtime(3),'    ',ipvadesc(3)
        WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003) outsn(1:lnblnk(outsn))
      elseif(ipvsimu(1:6).eq.'icautt'.or. ipvsimu(1:6).eq.'icautf'.or.
     &       ipvsimu(1:6).eq.'icauts')then
        WRITE(outsn,'(A,i4,2i7,3F8.3,2a)',IOSTAT=IOS,ERR=1003)
     &    '*autumn  ',
     &    ipvastjd(4),ipvafnjd(4),ipvdispjd(4),ddmheat(4),
     &    ddmcool(4),ddmtime(4),'    ',ipvadesc(4)
        WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003) outsn(1:lnblnk(outsn))
      elseif(ipvsimu(1:3).eq.'ias')then   ! Annual simulation
        WRITE(outsn,'(A,i4,2i7,3F8.3,2a)',IOSTAT=IOS,ERR=1003)
     &    '*annual  ',
     &    ipvastjd(1),ipvafnjd(1),ipvdispjd(1),ddmheat(1),
     &    ddmcool(1),ddmtime(1),'    ',ipvadesc(1)
        WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003) outsn(1:lnblnk(outsn))
      elseif(ipvsimu(1:3).eq.'i3s'.or.ipvsimu(1:3).eq.'i3t'.or.
     &       ipvsimu(1:3).eq.'i3f')then  ! Three season variant
        WRITE(outsn,'(A,i4,2i7,3F8.3,2a)',IOSTAT=IOS,ERR=1003)
     &    '*winter  ',
     &    ipvastjd(1),ipvafnjd(1),ipvdispjd(1),ddmheat(1),
     &    ddmcool(1),ddmtime(1),'    ',ipvadesc(1)
        WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003) outsn(1:lnblnk(outsn))
        WRITE(outsn,'(A,i4,2i7,3F8.3,2a)',IOSTAT=IOS,ERR=1003)
     &    '*spring  ',
     &    ipvastjd(2),ipvafnjd(2),ipvdispjd(2),ddmheat(2),
     &    ddmcool(2),ddmtime(2),'    ',ipvadesc(2)
        WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003) outsn(1:lnblnk(outsn))
        WRITE(outsn,'(A,i4,2i7,3F8.3,2a)',IOSTAT=IOS,ERR=1003)
     &    '*summer  ',
     &    ipvastjd(3),ipvafnjd(3),ipvdispjd(3),ddmheat(3),
     &    ddmcool(3),ddmtime(3),'    ',ipvadesc(3)
        WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003) outsn(1:lnblnk(outsn))
      elseif(ipvsimu(1:3).eq.'i5s'.or.ipvsimu(1:3).eq.'i5t'.or.
     &       ipvsimu(1:3).eq.'i5f')then  ! Five season variant
        WRITE(outsn,'(A,i4,2i7,3F8.3,2a)',IOSTAT=IOS,ERR=1003)
     &    '*winter1 ',
     &    ipvastjd(1),ipvafnjd(1),ipvdispjd(1),ddmheat(1),
     &    ddmcool(1),ddmtime(1),'    ',ipvadesc(1)
        WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003) outsn(1:lnblnk(outsn))
        WRITE(outsn,'(A,i4,2i7,3F8.3,2a)',IOSTAT=IOS,ERR=1003)
     &    '*spring  ',
     &    ipvastjd(2),ipvafnjd(2),ipvdispjd(2),ddmheat(2),
     &    ddmcool(2),ddmtime(2),'    ',ipvadesc(2)
        WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003) outsn(1:lnblnk(outsn))
        WRITE(outsn,'(A,i4,2i7,3F8.3,2a)',IOSTAT=IOS,ERR=1003)
     &    '*summer  ',
     &    ipvastjd(3),ipvafnjd(3),ipvdispjd(3),ddmheat(3),
     &    ddmcool(3),ddmtime(3),'    ',ipvadesc(3)
        WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003) outsn(1:lnblnk(outsn))
        WRITE(outsn,'(A,i4,2i7,3F8.3,2a)',IOSTAT=IOS,ERR=1003)
     &    '*autumn  ',
     &    ipvastjd(4),ipvafnjd(4),ipvdispjd(4),ddmheat(4),
     &    ddmcool(4),ddmtime(4),'    ',ipvadesc(4)
        WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003) outsn(1:lnblnk(outsn))
        WRITE(outsn,'(A,i4,2i7,3F8.3,2a)',IOSTAT=IOS,ERR=1003)
     &    '*winter2 ',
     &    ipvastjd(5),ipvafnjd(5),ipvdispjd(1),ddmheat(5),
     &    ddmcool(5),ddmtime(5),'    ',ipvadesc(5)
        WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003) outsn(1:lnblnk(outsn))
      endif

C Only in the case of an annual assessment write out the display days.
      if(ipvsimu(1:3).eq.'ias')then   ! Annual simulation
        WRITE(outsn,'(A,i3,10I4)',IOSTAT=IOS,ERR=1003)'*display_days',
     &    nipvdispjd,(ipvdispjd(L),L=1,nipvdispjd)
        call SDELIM(outsn,outsd,'C',IW)
        WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003) outsd(1:lnblnk(outsd))
      endif

C If no metrics or demands close off the file.
      if(nms.eq.0.and.neds.eq.0)then
        WRITE(outsn,'(A)',IOSTAT=IOS,ERR=1003)'*metrics none'
        write(ifu,'(a)')  outsn(1:lnblnk(outsn))
        WRITE(outsn,'(A)',IOSTAT=IOS,ERR=1003)'*demands  none'
        write(ifu,'(a)') outsn(1:lnblnk(outsn))
        if(act(1:3).eq.'ipv')then
          WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003)'*end'
          CALL  ERPFREE(ifu,ISTAT)
        elseif(act(1:3).eq.'cfg')then
          WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003)'*end_ipv'
        endif
        return
      endif

C If demands only jump.
      if(nms.eq.0.and.neds.gt.0)then
        WRITE(outsn,'(A)',IOSTAT=IOS,ERR=1003)'*metrics none'
        write(ifu,'(a)')  outsn(1:lnblnk(outsn))
        goto 42
      endif


C Write metrics column labels and data depending on the type. Some topics
C are restricted to the focus zone and adapt for these.
C [ idea ] set metgroup for glare & DF to the focus zone name]
      do 28 ijj=1,nms
        lnd=lnblnk(msdoc(ijj))
        lnz=lnblnk(metgroup(ijj))
        lng=lnblnk(metrglbl(ijj))
        lnipvz=lnblnk(ivpfocname)
        if(msdoc(ijj)(1:4).eq.'Guth')then
          WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003) '# ' 
          WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003) '# metrics' 
          WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003)
     &    '# metric      focus     units    VCP -60 to 60 deg'
          WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003)
     &    '# name        zone           '
          WRITE(outsn,'(4a)',IOSTAT=IOS,ERR=1003)
     &      '*',msdoc(ijj),' ',metgroup(ijj)(1:lnz)
          WRITE(ifu,'(4A,13f6.2)',IOSTAT=IOS,ERR=1003)
     &      outsn(1:lnblnk(outsn)),tab,metrglbl(ijj)(1:lng),'  ',
     &      (vcp(ij),ij=1,13)
        elseif(msdoc(ijj)(1:5).eq.'glare')then
          WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003) '# ' 
          WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003)
     &    '# metric      focus     units'
          WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003)
     &    '# name        zone           '
          WRITE(outsn,'(4a)',IOSTAT=IOS,ERR=1003)
     &      '*',msdoc(ijj),' ',metgroup(ijj)(1:lnz)
          WRITE(ifu,'(4A)',IOSTAT=IOS,ERR=1003)
     &      outsn(1:lnblnk(outsn)),tab,metrglbl(ijj)(1:lng)
        elseif(msdoc(ijj)(1:12).eq.'daylightfact')then
          WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003) '# ' 
          WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003)
     &    '# metric      focus      units'
          WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003)
     &    '# name        zone            '
          WRITE(outsn,'(4a)',IOSTAT=IOS,ERR=1003)
     &      '*',msdoc(ijj),' ',metgroup(ijj)(1:lnz)
          WRITE(ifu,'(4A)',IOSTAT=IOS,ERR=1003)
     &      outsn(1:lnblnk(outsn)),tab,metrglbl(ijj)(1:lng)
        elseif(msdoc(ijj)(1:12).eq.'contaminant')then
          WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003) '# ' 
          WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003)
     &    '# metric      focus      flow      contam-      units'
          WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003)
     &    '# name        zone       node      inate'
          WRITE(outsn,'(7a)',IOSTAT=IOS,ERR=1003)
     &      '*',msdoc(ijj),',',metgroup(ijj)(1:lnz),',',
     &      ipvndname(1:lnblnk(ipvndname)),',',
     &      ipvctmname(1:lnblnk(ipvctmname))
          WRITE(ifu,'(4A)',IOSTAT=IOS,ERR=1003)
     &      outsn(1:lnblnk(outsn)),tab,metrglbl(ijj)(1:lng)
        else
          WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003) '# ' 
          WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003)
     &  '# metric     number floor  group    units'
          WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003)
     &  '# name       zones  area   name          '
          WRITE(outsn,'(2A,i4,F8.1,2a)',
     &      IOSTAT=IOS,ERR=1003)'*',msdoc(ijj),
     &      nzmg(ijj),emgflr(ijj),',',metgroup(ijj)(1:lnz)
          WRITE(ifu,'(3A)',IOSTAT=IOS,ERR=1003)
     &      outsn(1:lnblnk(outsn)),tab,metrglbl(ijj)(1:lng)

C Copy izmg to single array list and then write via ailist call.
          if(nzmg(ijj).gt.0)then
            do ijjj=1,nzmg(ijj)
              list(ijjj)=izmg(ijj,ijjj)
            enddo
            itrunc=1
            ipos=1
            do while (itrunc.ne.0)
              call ailist(ipos,nzmg(ijj),list,MCOM,'S',louts,loutlen,
     &          itrunc)
              write(ifu,'(1x,a)',IOSTAT=ios,ERR=1003) louts(1:loutlen)
              ipos=itrunc+1
            end do
          endif
        endif
  28  continue

C Loop through each demand set and write out associated data.
  42  if(neds.gt.0)then
C        WRITE(outsn,'(A,i4)',IOSTAT=IOS,ERR=1003)'*demand_sets',neds
C        write(ifu,'(a)') outsn(1:lnblnk(outsn))
        WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003) '# ' 
        WRITE(ifu,'(2a)',IOSTAT=IOS,ERR=1003) 
     &    '# metric  number  floor  scal- group'
        WRITE(ifu,'(2a)',IOSTAT=IOS,ERR=1003)
     &    '# name    zones   area   ing   name'
        do 43 ij=1,neds

C Write set of energy (demand/capacity) associated data.
          if(ij.le.9)then
            WRITE(outsn,'(A,i1,I6,F8.1,F7.3,2a)',IOSTAT=IOS,ERR=1003)
     &      '*demand_',ij,nzedg(ij),
     &      edgflr(ij),edgsca(ij),',',zedsdoc(ij)(1:lnblnk(zedsdoc(ij)))
            WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003)outsn(1:lnblnk(outsn))
          else
            WRITE(outsn,'(A,i2,I6,F8.1,F7.3,2a)',IOSTAT=IOS,ERR=1003)
     &      '*demand_',ij,nzedg(ij),
     &      edgflr(ij),edgsca(ij),',',zedsdoc(ij)(1:lnblnk(zedsdoc(ij)))
            WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003)outsn(1:lnblnk(outsn))
          endif

C Copy izedg into list and the write out via ailist call.
          if(nzedg(ij).gt.0)then
            do ijj=1,nzedg(ij)
              list(ijj)=izedg(ij,ijj)
            enddo
            itrunc=1
            ipos=1
            do while (itrunc.ne.0)
              call ailist(ipos,nzedg(ij),list,MCOM,'S',louts,loutlen,
     &          itrunc)
              write(ifu,'(1x,a)',IOSTAT=ios,ERR=1003) louts(1:loutlen)
              ipos=itrunc+1
            end do
          endif
  43    continue
      endif

C Write nonspecific demands data.
      WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003)
     &'# non-specific demands'
      WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003)
     &'# tag occup/light/sm pwr/fans/pumps/lifts/dhw'
      WRITE(outsn,'(A,7I6)',IOSTAT=IOS,ERR=1003)'*dmds',idmdinc(1),
     & idmdinc(2),idmdinc(3),idmdinc(4),idmdinc(5),idmdinc(6),idmdinc(7)
      WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003)outsn(1:lnblnk(outsn))

C If writing to external file use a different syntax and free the file.
C If writing to cfg file use *end_ipv and do not close file.
      if(act(1:3).eq.'ipv')then
        WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003)'*end'
        CALL  ERPFREE(ifu,ISTAT)
      elseif(act(1:3).eq.'cfg')then
        WRITE(ifu,'(A)',IOSTAT=IOS,ERR=1003)'*end_ipv'
      endif
      return

  999 CALL  ERPFREE(ifu,ISTAT)
      return

 1003 if(ios.eq.2)then
        CALL USRMSG('No premissopm to write IPV data.',' ','W')
        GOTO 999
      else
        CALL USRMSG('Problem writing IPV data.',' ','W')
        GOTO 999
      endif

      end

C *********** clearipvdat
C Clear IPV common blocks. If act = 'b' blank all the data, if act = 'i'
C then include initial data for energy demand sets and metrics.
      subroutine clearipvdat(act)
#include "building.h"
#include "model.h"
#include "ipvdata.h"
      
      integer lnblnk  ! function definition

C Comments about common blocks are in ripvdat.
      character act*1

C Days to display (default is three (early winter, spring, summer).
      nipvdispjd=3
      ipvdispjd(1)=11
      ipvdispjd(2)=67
      ipvdispjd(3)=186
      ipvfoczone=0
      ivpfocname='UNKNOWN'
      do 342 ij=1,MIPVA
        ddmheat(ij)=1.0
        ddmcool(ij)=1.0
        ddmtime(ij)=1.0
        ipvadesc(ij)=' '
  342 continue

C Images.
      nipvimg=0
      lipvimg(1)='UNKNOWN'
      lipvimg(2)='UNKNOWN'
      lipvimg(3)='UNKNOWN'
      lipvimg(4)='UNKNOWN'

C Items to clear in any case.
      do 50 ij=1,MIPVM
        nzmg(ij)=0
        imetmsc(ij,1)=0
        imetmsc(ij,2)=0
        emgflr(ij)=1.0
        izmg(ij,1)=0
        nzedg(ij)=0
        edgflr(ij)=1.0
        edgsca(ij)=1.0
        izedg(ij,1)=0
 50   continue

C Clear commons with blanking.
      if(act.eq.'b')then
        ipvsynop=' '
        nms=0
        do 65 ij=1,MIPVM
          imetget(ij)=0
          msdoc(ij)=' '
          metrglbl(ij)='- '
          metgroup(ij)='- '
  65    continue

C Energy demand sets. Assume aggregate demands included, timestep
C reporting and the display day only for tabular reports. 
        neds=0
        iaggr=1     ! Always include for focus day or week or fortnight.
        ifbhits=1   ! Assume frequency bins report hits
        do 67 ij=1,MIPVM
          zedsdoc(ij)='- '
  67    continue
      elseif(act.eq.'i'.or.act.eq.'I')then
        ipvsimu=' '
        lr=lnblnk(cfgroot)
        write(ipvtitl,'(a)') modeltitle(1:40)
        write(ipvsynop,'(2a)') ':The ',cfgroot(1:lr)

C Metric sets.
        nms=1
        do 165 ij=1,MIPVM
          if(ij.eq.1)then
            imetget(ij)=6
            msdoc(ij)='zones'
            metgroup(ij)='ocup_zones'
            metrglbl(ij)='Resultant T (degC)'
          else
            imetget(ij)=0
            msdoc(ij)=' '
            metgroup(ij)='- '
            metrglbl(ij)='- '
          endif
 165    continue

C Energy demand sets.
        neds=1
        do 167 ij=1,MIPVM
          if(ij.eq.1)then
            zedsdoc(ij)='zones'
          else
            zedsdoc(ij)='- '
          endif
 167    continue
      endif

      return
      end

C *********** ripvdat
C Read IPV data from an external file or within cfg file based
C on the parameter act (3 char) ='ipv' or ='cfg'.
      subroutine ripvdat(ifu,file,act,ier)
#include "building.h"
#include "model.h"
#include "geometry.h"
#include "espriou.h"
#include "ipvdata.h"

      integer lnblnk  ! function definition
      integer ncomp,ncon
      common/c1/ncomp,ncon
      
      integer ipvictmindex    ! which contaminate source
      character ipvctmname*12 ! name of contaminate
      character ipvndname*12  ! associated flow node name
      common/contm10/ipvictmindex,ipvctmname,ipvndname

      real vcp  ! Guth visual comfort probablility -60 - 60 values from Radiance.
      common/guth/vcp(13)

      character file*72,word*20,vword*20
      character OUTSTR*124,LOUTSTR*248
      character act*3   ! indicate which file destination
      integer list
      dimension list(MCOM)  ! zone lists currently set to MCOM
      character tipvadesc*40
      character vcptext*96

C Comments about common blocks are in ripvdat.

      iset = 0
      imset= 0

C Open any existing file by this name.
      if(act(1:3).eq.'ipv')then
        CALL EFOPSEQ(ifu,file,1,IER)
        IF(IER.LT.0)THEN
          IER=1
          RETURN
        ENDIF
        write(currentfile,'(a)') file(1:lnblnk(file))

C Check identity and if that is successful then get the version
C number and the format indicator.
        CALL STRIPC(ifu,OUTSTR,99,ND,1,'line 1',IER)
        IF(IER.NE.0)RETURN
        if(OUTSTR(1:5).ne.'* IPV'.and.OUTSTR(1:4).ne.'*IPV')then
          call usrmsg('Is not an IPV description file.',' ','W')
          IER=1
          return
        endif
        K=0
        CALL EGETW(OUTSTR,K,WORD,'W','*IPV tag',IFLAG)
        if(ND.eq.2)then
          CALL EGETW(OUTSTR,K,VWORD,'W','IPV V',IFLAG)  ! Compact header
          if(VWORD(1:2).eq.'V4')then
            ipvversion=4
          elseif(VWORD(1:2).eq.'V5')then
            ipvversion=5
          elseif(VWORD(1:1).eq.'4')then
            ipvversion=4
          elseif(VWORD(1:1).eq.'5')then
            ipvversion=5
          endif
        elseif(ND.gt.2)then
          CALL EGETWI(OUTSTR,K,ipvversion,0,5,'W','version index',IER)
          CALL EGETWI(OUTSTR,K,ipvform,0,2,'W','ipv rep format',IER)
          if(ND.gt.3)then
            CALL EGETW(OUTSTR,K,ipvsimu,'W','*IPV assessment',IFLAG)
          endif
        else
          ipvversion=3
          ipvform=0
          ipvsimu='------'
        endif
      elseif(act(1:3).eq.'cfg')then
        continue  ! the header information has already been scanned.
      endif

C Clear (blank) IPV commons before reading.
      call clearipvdat('b')
      iper=0
      imset=0
      iset=0
      ipvform=1       ! Assume normal report format unless reset.

C Scan lines.
  43  CALL LSTRIPC(ifu,LOUTSTR,99,ND,1,'IPV data',IER)
      K=0
      CALL EGETW(LOUTSTR,K,WORD,'W','tag',IFLAG)
      if(WORD(1:6).eq.'*title')then
        CALL EGETRM(LOUTSTR,K,ipvtitl,'W','title string',IER)
      elseif(WORD(1:8).eq.'*version')then
        continue   ! optional ignore if tag found
      elseif(WORD(1:5).eq.'*date')then
        continue   ! optional ignore if tag found
      elseif(WORD(1:6).eq.'*synop')then
        write(ipvsynop,'(a)') LOUTSTR(8:lnblnk(LOUTSTR))
        goto 43
      elseif(WORD(1:14).eq.'*report_format')then
        CALL EGETW(LOUTSTR,K,VWORD,'W','report format',IFLAG)
        if(VWORD(1:3).eq.'csv')then
          ipvform=2
        else
          ipvform=1
        endif
        goto 43
      elseif(WORD(1:7).eq.'*images')then
        CALL EGETWI(LOUTSTR,K,nipvimg,0,4,'W','no of images',IER)
        if(nipvimg.gt.0)then
          do i=1,nipvimg
            CALL LSTRIPC(ifu,LOUTSTR,0,ND,1,'image file',IER)
            lipvimg(i)=LOUTSTR(1:72)
          enddo
        endif
      elseif(WORD(1:15).eq.'*building_image')then
        CALL EGETRM(LOUTSTR,K,lipvimg(1),'W','building image',IER)
        nipvimg=1
        goto 43
      elseif(WORD(1:12).eq.'*glare_image')then
        CALL EGETRM(LOUTSTR,K,lipvimg(2),'W','glare image',IER)
        nipvimg=1
        goto 43
      elseif(WORD(1:6).eq.'*focus')then
        CALL EGETRM(LOUTSTR,K,ivpfocname,'W','focus zone',IER)
        lnw=lnblnk(ivpfocname)
        if(ncomp.gt.0)then
          do iz=1,ncomp
            lnzz=lnblnk(zname(iz))  ! Find matching zone index.
            if(VWORD(1:lnw).eq.zname(iz)(1:lnzz))then
              ipvfoczone=iz
              cycle
            endif
          enddo
        endif
        goto 43
      elseif(WORD(1:12).eq.'*simulations')then
        CALL EGETWI(LOUTSTR,K,nipvassmt,1,5,'W','simulations',IER)
        if(ND.gt.2)then
          CALL EGETW(LOUTSTR,K,ipvsimu,'W','*IPV assessment',IFLAG)
        endif
        iper=0
        if(ipvsimu.eq.'------')then

C We have an older file which did not contain ipvsimu so make a guess
C at its value by looking at the number of assessments.
          if(nipvassmt.eq.1)then
            ipvsimu='ias   '
          elseif(nipvassmt.eq.3)then
            ipvsimu='i3t   '
          elseif(nipvassmt.eq.5)then
            ipvsimu='i5t   '
          else
            ipvsimu='i--   '
          endif
        endif
      elseif(WORD(1:7).eq.'*period')then

C If the initial portion of ipvadesc is the same as the model root
C name then discard the initial portion.
        iper=iper+1
        CALL EGETWI(LOUTSTR,K,ipvastjd(iper),1,365,'W','assmt st',IER)
        CALL EGETWI(LOUTSTR,K,ipvafnjd(iper),1,365,'W','assmt fn',IER)
        tipvadesc=' '
        CALL EGETRM(LOUTSTR,K,tipvadesc,'W','assmt descr',IER)
        lncf=lnblnk(cfgroot)
        lnipva=lnblnk(tipvadesc)
        if(tipvadesc(1:lncf).eq.cfgroot(1:lncf))then
          write(ipvadesc(iper),'(a)') tipvadesc(lncf+1:lnipva)
        else
          write(ipvadesc(iper),'(a)') tipvadesc(1:lnipva)
        endif

      elseif(WORD(1:13).eq.'*display_days')then
        CALL EGETWI(LOUTSTR,K,nipvdispjd,1,10,'W','displays',IER)
        do 343 ij=1,nipvdispjd
          CALL EGETWI(LOUTSTR,K,ipvdispjd(ij),1,365,'W','display',IER)
  343   continue

      elseif(WORD(1:7).eq.'*annual')then
        if(ipvversion.le.4)then
          CALL EGETWR(LOUTSTR,K,ddmheat(1),0.,100.,'W','ann hdd',IER)
          CALL EGETWR(LOUTSTR,K,ddmcool(1),0.,100.,'W','ann cdd',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(1),0.,100.,'W','ann time',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(1),0.,100.,'W','ann time',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(1),0.,100.,'W','ann time',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(1),0.,100.,'W','ann time',IER)
        else
          CALL EGETWR(LOUTSTR,K,ddmheat(1),0.,100.,'W','ann hdd',IER)
          CALL EGETWR(LOUTSTR,K,ddmcool(1),0.,100.,'W','ann cdd',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(1),0.,100.,'W','ann time',IER)
        endif
      elseif(WORD(1:5).eq.'*win1')then
        if(ipvversion.le.4)then
          CALL EGETWR(LOUTSTR,K,ddmheat(1),0.,100.,'W','win h dd',IER)
          CALL EGETWR(LOUTSTR,K,ddmcool(1),0.,100.,'W','win c dd',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(1),0.,100.,'W','win fn tim',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(1),0.,100.,'W','win fn tim',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(1),0.,100.,'W','win fn tim',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(1),0.,100.,'W','win hw tim',IER)
        else
          CALL EGETWR(LOUTSTR,K,ddmheat(1),0.,100.,'W','win h dd',IER)
          CALL EGETWR(LOUTSTR,K,ddmcool(1),0.,100.,'W','win c dd',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(1),0.,100.,'W','win time',IER)
        endif
      elseif(WORD(1:8).eq.'*winter2')then  ! Trap first for 5 seasons
        iper=iper+1
        CALL EGETWI(LOUTSTR,K,ipvastjd(iper),1,365,'W','assmt st',IER)
        CALL EGETWI(LOUTSTR,K,ipvafnjd(iper),1,365,'W','assmt fn',IER)
        CALL EGETWI(LOUTSTR,K,ipvdispjd(iper),1,365,'W','display',IER)
        CALL EGETWR(LOUTSTR,K,ddmheat(iper),0.,100.,'W','win2 h dd',IER)
        CALL EGETWR(LOUTSTR,K,ddmcool(iper),0.,100.,'W','win2 c dd',IER)
        CALL EGETWR(LOUTSTR,K,ddmtime(iper),0.,100.,'W','win2 time',IER)
        tipvadesc=' '
        CALL EGETRM(LOUTSTR,K,tipvadesc,'W','assmt descr',IER)
        lncf=lnblnk(cfgroot)
        lnipva=lnblnk(tipvadesc)
        if(tipvadesc(1:lncf).eq.cfgroot(1:lncf))then
          write(ipvadesc(iper),'(a)') tipvadesc(lncf+1:lnipva)
        else
          write(ipvadesc(iper),'(a)') tipvadesc(1:lnipva)
        endif
      elseif(WORD(1:8).eq.'*winter1')then  ! Trap for 5 seasons
        iper=iper+1
        CALL EGETWI(LOUTSTR,K,ipvastjd(iper),1,365,'W','assmt st',IER)
        CALL EGETWI(LOUTSTR,K,ipvafnjd(iper),1,365,'W','assmt fn',IER)
        CALL EGETWI(LOUTSTR,K,ipvdispjd(iper),1,365,'W','display',IER)
        CALL EGETWR(LOUTSTR,K,ddmheat(iper),0.,100.,'W','win h dd',IER)
        CALL EGETWR(LOUTSTR,K,ddmcool(iper),0.,100.,'W','win c dd',IER)
        CALL EGETWR(LOUTSTR,K,ddmtime(iper),0.,100.,'W','win time',IER)
        tipvadesc=' '
        CALL EGETRM(LOUTSTR,K,tipvadesc,'W','assmt descr',IER)
        lncf=lnblnk(cfgroot)
        lnipva=lnblnk(tipvadesc)
        if(tipvadesc(1:lncf).eq.cfgroot(1:lncf))then
          write(ipvadesc(iper),'(a)') tipvadesc(lncf+1:lnipva)
        else
          write(ipvadesc(iper),'(a)') tipvadesc(1:lnipva)
        endif
      elseif(WORD(1:7).eq.'*winter')then  ! Trap for 1 or 3 seasons
        iper=iper+1
        CALL EGETWI(LOUTSTR,K,ipvastjd(iper),1,365,'W','assmt st',IER)
        CALL EGETWI(LOUTSTR,K,ipvafnjd(iper),1,365,'W','assmt fn',IER)
        CALL EGETWI(LOUTSTR,K,ipvdispjd(iper),1,365,'W','display',IER)
        CALL EGETWR(LOUTSTR,K,ddmheat(iper),0.,100.,'W','win h dd',IER)
        CALL EGETWR(LOUTSTR,K,ddmcool(iper),0.,100.,'W','win c dd',IER)
        CALL EGETWR(LOUTSTR,K,ddmtime(iper),0.,100.,'W','win time',IER)
        tipvadesc=' '
        CALL EGETRM(LOUTSTR,K,tipvadesc,'W','assmt descr',IER)
        lncf=lnblnk(cfgroot)
        lnipva=lnblnk(tipvadesc)
        if(tipvadesc(1:lncf).eq.cfgroot(1:lncf))then
          write(ipvadesc(iper),'(a)') tipvadesc(lncf+1:lnipva)
        else
          write(ipvadesc(iper),'(a)') tipvadesc(1:lnipva)
        endif
      elseif(WORD(1:7).eq.'*spring')then
        iper=iper+1
        CALL EGETWI(LOUTSTR,K,ipvastjd(iper),1,365,'W','assmt st',IER)
        CALL EGETWI(LOUTSTR,K,ipvafnjd(iper),1,365,'W','assmt fn',IER)
        CALL EGETWI(LOUTSTR,K,ipvdispjd(iper),1,365,'W','display',IER)
        CALL EGETWR(LOUTSTR,K,ddmheat(iper),0.,100.,'W','spr h dd',IER)
        CALL EGETWR(LOUTSTR,K,ddmcool(iper),0.,100.,'W','spr c dd',IER)
        CALL EGETWR(LOUTSTR,K,ddmtime(iper),0.,100.,'W','pr time',IER)
        tipvadesc=' '
        CALL EGETRM(LOUTSTR,K,tipvadesc,'W','assmt descr',IER)
        lncf=lnblnk(cfgroot)
        lnipva=lnblnk(tipvadesc)
        if(tipvadesc(1:lncf).eq.cfgroot(1:lncf))then
          write(ipvadesc(iper),'(a)') tipvadesc(lncf+1:lnipva)
        else
          write(ipvadesc(iper),'(a)') tipvadesc(1:lnipva)
        endif
      elseif(WORD(1:7).eq.'*summer')then
        iper=iper+1
        CALL EGETWI(LOUTSTR,K,ipvastjd(iper),1,365,'W','assmt st',IER)
        CALL EGETWI(LOUTSTR,K,ipvafnjd(iper),1,365,'W','assmt fn',IER)
        CALL EGETWI(LOUTSTR,K,ipvdispjd(iper),1,365,'W','display',IER)
        CALL EGETWR(LOUTSTR,K,ddmheat(iper),0.,100.,'W','sum h dd',IER)
        CALL EGETWR(LOUTSTR,K,ddmcool(iper),0.,100.,'W','sum c dd',IER)
        CALL EGETWR(LOUTSTR,K,ddmtime(iper),0.,100.,'W','pr time',IER)
        tipvadesc=' '
        CALL EGETRM(LOUTSTR,K,tipvadesc,'W','assmt descr',IER)
        lncf=lnblnk(cfgroot)
        lnipva=lnblnk(tipvadesc)
        if(tipvadesc(1:lncf).eq.cfgroot(1:lncf))then
          write(ipvadesc(iper),'(a)') tipvadesc(lncf+1:lnipva)
        else
          write(ipvadesc(iper),'(a)') tipvadesc(1:lnipva)
        endif
      elseif(WORD(1:7).eq.'*autumn')then
        iper=iper+1
        CALL EGETWI(LOUTSTR,K,ipvastjd(iper),1,365,'W','assmt st',IER)
        CALL EGETWI(LOUTSTR,K,ipvafnjd(iper),1,365,'W','assmt fn',IER)
        CALL EGETWI(LOUTSTR,K,ipvdispjd(iper),1,365,'W','display',IER)
        CALL EGETWR(LOUTSTR,K,ddmheat(iper),0.,100.,'W','aut h dd',IER)
        CALL EGETWR(LOUTSTR,K,ddmcool(iper),0.,100.,'W','aut c dd',IER)
        CALL EGETWR(LOUTSTR,K,ddmtime(iper),0.,100.,'W','pr time',IER)
        tipvadesc=' '
        CALL EGETRM(LOUTSTR,K,tipvadesc,'W','assmt descr',IER)
        lncf=lnblnk(cfgroot)
        lnipva=lnblnk(tipvadesc)
        if(tipvadesc(1:lncf).eq.cfgroot(1:lncf))then
          write(ipvadesc(iper),'(a)') tipvadesc(lncf+1:lnipva)
        else
          write(ipvadesc(iper),'(a)') tipvadesc(1:lnipva)
        endif
      elseif(WORD(1:4).eq.'*trn'.or.WORD(1:4).eq.'*spr')then
        if(ipvversion.le.4)then
          CALL EGETWR(LOUTSTR,K,ddmheat(2),0.,100.,'W','spr h dd',IER)
          CALL EGETWR(LOUTSTR,K,ddmcool(2),0.,100.,'W','spr c dd',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(2),0.,100.,'W','spr fn tim',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(2),0.,100.,'W','spr fn tim',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(2),0.,100.,'W','spr fn tim',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(2),0.,100.,'W','spr hw tim',IER)
        else
          CALL EGETWR(LOUTSTR,K,ddmheat(2),0.,100.,'W','spr h dd',IER)
          CALL EGETWR(LOUTSTR,K,ddmcool(2),0.,100.,'W','spr c dd',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(2),0.,100.,'W','spr fn tim',IER)
        endif
      elseif(WORD(1:4).eq.'*sum')then
        if(ipvversion.le.4)then
          CALL EGETWR(LOUTSTR,K,ddmheat(3),0.,100.,'W','sum h dd',IER)
          CALL EGETWR(LOUTSTR,K,ddmcool(3),0.,100.,'W','sum c dd',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(3),0.,100.,'W','sum fn tim',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(3),0.,100.,'W','sum fn tim',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(3),0.,100.,'W','sum fn tim',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(3),0.,100.,'W','sum hw tim',IER)
        else
          CALL EGETWR(LOUTSTR,K,ddmheat(3),0.,100.,'W','sum h dd',IER)
          CALL EGETWR(LOUTSTR,K,ddmcool(3),0.,100.,'W','sum c dd',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(3),0.,100.,'W','sum fn tim',IER)
        endif
      elseif(WORD(1:4).eq.'*aut')then
        if(ipvversion.le.4)then
          CALL EGETWR(LOUTSTR,K,ddmheat(4),0.,100.,'W','aut h dd',IER)
          CALL EGETWR(LOUTSTR,K,ddmcool(4),0.,100.,'W','aut c dd',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(4),0.,100.,'W','aut fn tim',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(4),0.,100.,'W','aut fn tim',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(4),0.,100.,'W','aut fn tim',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(4),0.,100.,'W','aut hw tim',IER)
        else
          CALL EGETWR(LOUTSTR,K,ddmheat(4),0.,100.,'W','aut h dd',IER)
          CALL EGETWR(LOUTSTR,K,ddmcool(4),0.,100.,'W','aut c dd',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(4),0.,100.,'W','aut fn tim',IER)
        endif
      elseif(WORD(1:5).eq.'*win2')then
        if(ipvversion.le.4)then
          CALL EGETWR(LOUTSTR,K,ddmheat(5),0.,100.,'W','win2 h dd',IER)
          CALL EGETWR(LOUTSTR,K,ddmcool(5),0.,100.,'W','win2 c dd',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(5),0.,100.,'W','win2 fn tm',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(5),0.,100.,'W','win2 fn tm',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(5),0.,100.,'W','win2 fn tm',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(5),0.,100.,'W','win2 hw tm',IER)
        else
          CALL EGETWR(LOUTSTR,K,ddmheat(5),0.,100.,'W','win2 h dd',IER)
          CALL EGETWR(LOUTSTR,K,ddmcool(5),0.,100.,'W','win2 c dd',IER)
          CALL EGETWR(LOUTSTR,K,ddmtime(5),0.,100.,'W','win2 fn tm',IER)
        endif

C One of the sets of performance metrics and associated data.
C First read out how many metrics, then for each metric write
C out its zones, floor area, weighting and set description.

      elseif(WORD(1:8).eq.'*metrics')then
        CALL EGETWI(LOUTSTR,K,nms,0,12,'W','nb of metric sets',IER)
        imset=0
      elseif(WORD(1:8).eq.'*metric_')then
        imset=imset+1
        CALL EGETWI(LOUTSTR,K,imetget(imset),0,99,'W','metric iget',IER)
        CALL EGETWI(LOUTSTR,K,imetmsc(imset,1),0,32,'W','metric msc 1',
     &    IER)
        CALL EGETWI(LOUTSTR,K,imetmsc(imset,2),0,32,'W','metric msc 2',
     &    IER)
        CALL EGETWI(LOUTSTR,K,nzmg(imset),0,MCOM,'W','nb of zones',IER)
        CALL EGETWR(LOUTSTR,K,emgflr(imset),1.,9999.,'W','metric flr',
     &    IER)
        CALL EGETWR(LOUTSTR,K,dummy,1.,9999.,'W','metric scale',IER)
        CALL EGETWR(LOUTSTR,K,dummy,0.,99.,'W','metric wtg',IER)

C Read the user defined set name and the short metric label. If version
C 4 also read in the group name for the zones associated with the metric.
        CALL EGETP(LOUTSTR,K,msdoc(imset),'W','metric set descr',IER)
        if(ipvversion.ge.4)then
          CALL EGETP(LOUTSTR,K,metgroup(imset),'W','metric group',IER)
        endif
        CALL EGETP(LOUTSTR,K,metrglbl(imset),'W','metric 20 descr',IER)

C Read in the associated zones (on multiple lines is needed). Be sure
C to update warning range if model complexity changes in building.h
        if(nzmg(imset).gt.0)then
          CALL EGETWIA(ifu,list,nzmg(imset),1,82,'W',
     &      'IPV metrics associated zones',IER)
          DO KV=1,nzmg(imset)
            izmg(imset,KV)=list(KV)
          ENDDO
        endif
        nms=imset
      elseif(WORD(1:8).eq.'*comfort'.or.WORD(1:8).eq.'*Comfort')then
        if(ipvversion.le.4)then
          imset=imset+1
          write(msdoc(imset),'(a)') 'comfort'
          CALL EGETWI(LOUTSTR,K,imetget(imset),0,99,'W',
     &       'comfort iget',IER)
          CALL EGETWI(LOUTSTR,K,imetmsc(imset,1),0,32,'W',
     &       'metric msc 1',IER)
          CALL EGETWI(LOUTSTR,K,imetmsc(imset,2),0,32,'W',
     &       'metric msc 2',IER)
          CALL EGETWI(LOUTSTR,K,nzmg(imset),0,MCOM,'W','nb zones',IER)
          CALL EGETWR(LOUTSTR,K,emgflr(imset),1.,9999.,'W',
     &       'comfort flr',IER)
          CALL EGETP(LOUTSTR,K,metgroup(imset),'W','comfort group',IER)
          CALL EGETP(LOUTSTR,K,metrglbl(imset),'W','comfort descr',IER)
        else
          imset=imset+1
          write(msdoc(imset),'(a)') 'comfort'
          CALL EGETWI(LOUTSTR,K,nzmg(imset),0,MCOM,'W','nb zones',IER)
          CALL EGETWR(LOUTSTR,K,emgflr(imset),1.,9999.,'W',
     &       'comfort flr',IER)
          CALL EGETP(LOUTSTR,K,metgroup(imset),'W','comfort group',IER)
          CALL EGETP(LOUTSTR,K,metrglbl(imset),'W','comfort descr',IER)
        endif

C Read in the associated zones (on multiple lines is needed). Be sure
C to update warning range if model complexity changes in building.h
        if(nzmg(imset).gt.0)then
          CALL EGETWIA(ifu,list,nzmg(imset),1,82,'W',
     &      'IPV metrics associated zones',IER)
          DO KV=1,nzmg(imset)
            izmg(imset,KV)=list(KV)
          ENDDO
        endif
        nms=imset
      elseif(WORD(1:10).eq.'*emissions')then
        imset=imset+1
        write(msdoc(imset),'(a)') 'emissions'
        if(ipvversion.le.4)then
          CALL EGETWI(LOUTSTR,K,imetget(imset),0,99,'W','emis iget',IER)
          CALL EGETWI(LOUTSTR,K,imetmsc(imset,1),0,32,'W',
     &      'metric msc 1',IER)
          CALL EGETWI(LOUTSTR,K,imetmsc(imset,2),0,32,'W',
     &      'metric msc 2',IER)
          CALL EGETWI(LOUTSTR,K,nzmg(imset),0,MCOM,'W','nb zones',IER)
          CALL EGETWR(LOUTSTR,K,emgflr(imset),1.,9999.,'W','emis flr',
     &      IER)
          CALL EGETP(LOUTSTR,K,metgroup(imset),'W','emis group',IER)
          CALL EGETP(LOUTSTR,K,metrglbl(imset),'W','emis descr',IER)
        else
          CALL EGETWI(LOUTSTR,K,nzmg(imset),0,MCOM,'W','nb zones',IER)
          CALL EGETWR(LOUTSTR,K,emgflr(imset),1.,9999.,'W','emis flr',
     &      IER)
          CALL EGETP(LOUTSTR,K,metgroup(imset),'W','emis group',IER)
          CALL EGETP(LOUTSTR,K,metrglbl(imset),'W','emis descr',IER)
        endif

C Read in the associated zones (on multiple lines is needed). Be sure
C to update warning range if model complexity changes in building.h
        if(nzmg(imset).gt.0)then
          CALL EGETWIA(ifu,list,nzmg(imset),1,82,'W',
     &      'IPV metrics associated zones',IER)
          DO KV=1,nzmg(imset)
            izmg(imset,KV)=list(KV)
          ENDDO
        endif
        nms=imset
      elseif(WORD(1:13).eq.'*Infiltration'.or.
     &       WORD(1:12).eq.'*Ventilation'.or.
     &       WORD(1:13).eq.'*TotalCasGain'.or.
     &       WORD(1:12).eq.'*SolarFacade'.or.
     &       WORD(1:12).eq.'*Solarabsorb'.or.
     &       WORD(1:7).eq.'*ZoneRH'.or.WORD(1:8).eq.'*ZonedbT')then
        imset=imset+1
        write(msdoc(imset),'(a)') WORD(2:lnblnk(WORD))
        if(ipvversion.le.4)then
          CALL EGETWI(LOUTSTR,K,imetget(imset),0,99,'W','met iget',IER)
          CALL EGETWI(LOUTSTR,K,imetmsc(imset,1),0,32,'W',
     &      'metric msc 1',IER)
          CALL EGETWI(LOUTSTR,K,imetmsc(imset,2),0,32,'W',
     &      'metric msc 2',IER)
          CALL EGETWI(LOUTSTR,K,nzmg(imset),0,MCOM,'W','nb zones',IER)
          CALL EGETWR(LOUTSTR,K,emgflr(imset),1.,9999.,'W','met flr',
     &      IER)
          CALL EGETP(LOUTSTR,K,metgroup(imset),'W','met group',IER)
          CALL EGETP(LOUTSTR,K,metrglbl(imset),'W','met descr',IER)
        else
          CALL EGETWI(LOUTSTR,K,nzmg(imset),0,MCOM,'W','nb zones',IER)
          CALL EGETWR(LOUTSTR,K,emgflr(imset),1.,9999.,'W','met flr',
     &      IER)
          CALL EGETP(LOUTSTR,K,metgroup(imset),'W','met group',IER)
          CALL EGETP(LOUTSTR,K,metrglbl(imset),'W','met descr',IER)
        endif

C Read in the associated zones (on multiple lines is needed). Be sure
C to update warning range if model complexity changes in building.h
        if(nzmg(imset).gt.0)then
          CALL EGETWIA(ifu,list,nzmg(imset),1,82,'W',
     &      'IPV metrics associated zones',IER)
          DO KV=1,nzmg(imset)
            izmg(imset,KV)=list(KV)
          ENDDO
        endif
        nms=imset
      elseif(WORD(1:10).eq.'*Plt_C+R&L')then
        imset=imset+1
        write(msdoc(imset),'(a)') 'Plt_C+R&L'
        if(ipvversion.le.4)then
          CALL EGETWI(LOUTSTR,K,imetget(imset),0,99,'W','plt iget',IER)
          CALL EGETWI(LOUTSTR,K,imetmsc(imset,1),0,32,'W',
     &      'metric msc 1',IER)
          CALL EGETWI(LOUTSTR,K,imetmsc(imset,2),0,32,'W',
     &      'metric msc 2',IER)
          CALL EGETWI(LOUTSTR,K,nzmg(imset),0,MCOM,'W','nb zones',IER)
          CALL EGETWR(LOUTSTR,K,emgflr(imset),1.,9999.,'W','plt flr',
     &      IER)
          CALL EGETP(LOUTSTR,K,metgroup(imset),'W','plt group',IER)
          CALL EGETP(LOUTSTR,K,metrglbl(imset),'W','plt descr',IER)
        else
          CALL EGETWI(LOUTSTR,K,nzmg(imset),0,MCOM,'W','nb zones',IER)
          CALL EGETWR(LOUTSTR,K,emgflr(imset),1.,9999.,'W','plt flr',
     &      IER)
          CALL EGETP(LOUTSTR,K,metgroup(imset),'W','plt group',IER)
          CALL EGETP(LOUTSTR,K,metrglbl(imset),'W','plt descr',IER)
        endif

C Read in the associated zones (on multiple lines is needed). Be sure
C to update warning range if model complexity changes in building.h
        if(nzmg(imset).gt.0)then
          CALL EGETWIA(ifu,list,nzmg(imset),1,82,'W',
     &      'IPV metrics associated zones',IER)
          DO KV=1,nzmg(imset)
            izmg(imset,KV)=list(KV)
          ENDDO
        endif
        nms=imset
      elseif(WORD(1:6).eq.'*glare'.or.WORD(1:6).eq.'*Glare')then
        imset=imset+1
        write(msdoc(imset),'(a)') 'glare'
        if(ipvversion.le.4)then
          CALL EGETWI(LOUTSTR,K,imetget(imset),0,99,'W','glare igt',IER)
          CALL EGETWI(LOUTSTR,K,imetmsc(imset,1),0,32,'W',
     &      'metric msc 1',IER)
          CALL EGETWI(LOUTSTR,K,imetmsc(imset,2),0,32,'W',
     &      'metric msc 2',IER)
          CALL EGETP(LOUTSTR,K,metgroup(imset),'W','glare group',IER)
          CALL EGETP(LOUTSTR,K,metrglbl(imset),'W','glare descr',IER)
        else
          CALL EGETP(LOUTSTR,K,metgroup(imset),'W','glare group',IER)
          CALL EGETP(LOUTSTR,K,metrglbl(imset),'W','glare descr',IER)
        endif
        nzmg(imset)=1; izmg(imset,1)=ipvfoczone
        nms=imset
      elseif(WORD(1:5).eq.'*DayF'.or.WORD(1:13).eq.'*daylightfact')then
        imset=imset+1
        write(msdoc(imset),'(a)') 'daylightfact'
        if(ipvversion.le.4)then
          CALL EGETWI(LOUTSTR,K,imetget(imset),0,99,'W','DF iget',IER)
          CALL EGETWI(LOUTSTR,K,imetmsc(imset,1),0,32,'W',
     &      'metric msc 1',IER)
          CALL EGETWI(LOUTSTR,K,imetmsc(imset,2),0,32,'W',
     &      'metric msc 2',IER)
          CALL EGETP(LOUTSTR,K,metgroup(imset),'W','DF group',IER)
          CALL EGETP(LOUTSTR,K,metrglbl(imset),'W','DF descr',IER)
        else
          CALL EGETP(LOUTSTR,K,metgroup(imset),'W','DF group',IER)
          CALL EGETP(LOUTSTR,K,metrglbl(imset),'W','DF descr',IER)
        endif
        nzmg(imset)=1; izmg(imset,1)=ipvfoczone
        nms=imset
      elseif(WORD(1:5).eq.'*Guth')then
        imset=imset+1
        write(msdoc(imset),'(a)') 'Guth'
        if(ipvversion.le.4)then
          CALL EGETWI(LOUTSTR,K,imetget(imset),0,99,'W','Guth igt',IER)
          CALL EGETWI(LOUTSTR,K,imetmsc(imset,1),0,32,'W',
     &      'metric msc 1',IER)
          CALL EGETWI(LOUTSTR,K,imetmsc(imset,2),0,32,'W',
     &      'metric msc 2',IER)
          CALL EGETP(LOUTSTR,K,metgroup(imset),'W','Guth group',IER)
          CALL EGETP(LOUTSTR,K,metrglbl(imset),'W','Guth descr',IER)
        else

C In newer files there will be 13 Guth VCP values following metrglbl.
          CALL EGETP(LOUTSTR,K,metgroup(imset),'W','Guth group',IER)
          CALL EGETP(LOUTSTR,K,metrglbl(imset),'W','Guth descr',IER)
          write(6,*) 'vcp? ',lnblnk(LOUTSTR),K
          vcptext='  '
          CALL EGETRM(LOUTSTR,K,vcptext,'W','vcp remainder',IER)
          if(lnblnk(vcptext).gt.30)then
            kk=0
            do iv=1,13
              CALL EGETWR(vcptext,KK,vcp(iv),0.,90.,'W','vcp',IER)
            enddo
            write(6,*) 'vcp is now ',vcp
          endif
        endif
        nzmg(imset)=1; izmg(imset,1)=ipvfoczone
        nms=imset
      elseif(WORD(1:12).eq.'*contaminant')then
        imset=imset+1
        write(msdoc(imset),'(a)') 'contaminant'
        if(ipvversion.le.4)then
          CALL EGETWI(LOUTSTR,K,imetget(imset),0,99,'W','cntm iget',IER)
          CALL EGETWI(LOUTSTR,K,imetmsc(imset,1),0,32,'W',
     &      'metric msc 1',IER)
          CALL EGETWI(LOUTSTR,K,imetmsc(imset,2),0,32,'W',
     &      'metric msc 2',IER)
          CALL EGETP(LOUTSTR,K,metgroup(imset),'W','cntm group',IER)
          CALL EGETP(LOUTSTR,K,metrglbl(imset),'W','cntm descr',IER)
        else
          CALL EGETP(LOUTSTR,K,metgroup(imset),'W','cntm group',IER)
          CALL EGETP(LOUTSTR,K,ipvndname,'W','cntm node',IER)
          CALL EGETP(LOUTSTR,K,ipvctmname,'W','source name',IER)
          CALL EGETP(LOUTSTR,K,metrglbl(imset),'W','cntm descr',IER)
        endif
        nzmg(imset)=1; izmg(imset,1)=ipvfoczone
        nms=imset

C One of the sets of energy (demand/capacity) associated data.
      elseif(WORD(1:12).eq.'*demand_sets')then
        CALL EGETWI(LOUTSTR,K,neds,0,MIPVM,'W','no of demand sets',IER)
        iset=0
      elseif(WORD(1:8).eq.'*energy_'.or.WORD(1:8).eq.'*demand_')then
        iset=iset+1
        if(ipvversion.le.4)then
          CALL EGETWI(LOUTSTR,K,nzedg(iset),0,MCOM,'W','nb zones',IER)
          CALL EGETWI(LOUTSTR,K,idummy,0,32,'W','demand msc',IER)
          CALL EGETWI(LOUTSTR,K,idummy,0,32,'W','demand msc',IER)
          CALL EGETWI(LOUTSTR,K,iaggr,0,1,'W','aggregate toggle',IER)
          CALL EGETWR(LOUTSTR,K,edgflr(iset),1.,9999.,'W','nrg flr',IER)
          CALL EGETWR(LOUTSTR,K,edgsca(iset),0.,99.,'W','nrg wtg',IER)
          CALL EGETRM(LOUTSTR,K,zedsdoc(iset),'W','enrg descr',IER)
        else
          CALL EGETWI(LOUTSTR,K,nzedg(iset),0,MCOM,'W','nb zones',IER)
          CALL EGETWR(LOUTSTR,K,edgflr(iset),1.,9999.,'W','nrg flr',IER)
          CALL EGETWR(LOUTSTR,K,edgsca(iset),0.,99.,'W','enrg wtg',IER)
          CALL EGETRM(LOUTSTR,K,zedsdoc(iset),'W','enrg descr',IER)
        endif

C Fill izedg array via call to egetwia returning list.
        if(nzedg(iset).gt.0)then
          CALL EGETWIA(ifu,list,nzedg(iset),1,MCOM,'W',
     &      'IPV dmds associated zones',IER)
          DO 119 KV=1,nzedg(iset)
            izedg(iset,KV)=list(KV)
  119     CONTINUE
        endif
        neds = iset   ! Update neds to reflect this demand.
C        write(6,*) 'neds is now ',neds

C Non-specific demands toggles for inclusion.
      elseif(WORD(1:5).eq.'*dmds')then
        CALL EGETWI(LOUTSTR,K,idmdinc(1),0,1,'W','occup dmd',IER)
        CALL EGETWI(LOUTSTR,K,idmdinc(2),0,1,'W','light dmd',IER)
        CALL EGETWI(LOUTSTR,K,idmdinc(3),0,1,'W','smpwr dmd',IER)
        CALL EGETWI(LOUTSTR,K,idmdinc(4),0,1,'W','fans dmd',IER)
        CALL EGETWI(LOUTSTR,K,idmdinc(5),0,1,'W','pumps dmd',IER)
        CALL EGETWI(LOUTSTR,K,idmdinc(6),0,1,'W','lifts dmd',IER)
        CALL EGETWI(LOUTSTR,K,idmdinc(7),0,1,'W','DHW dmd',IER)
      elseif(WORD(1:8).eq.'*end_ipv')then
        return  ! for ipv data within cfg file
      elseif(WORD(1:5).eq.'*end ')then
        goto 999  ! for separate ipv file
      else
        goto 43
      endif
      goto 43

  999 CALL  ERPFREE(ifu,ISTAT)
      return

      end


C *********** listipvdat ************
C List current IPV data common blocks.
C Where act='a' list all, act='m' list metrics, act='d' list
C demand sets. Itru is the channel for listing.
      subroutine listipvdat(itru,act,ier)
#include "building.h"
#include "seasons.h"
C seasons.h provides typper and typsea
#include "ipvdata.h"
      
      integer lnblnk  ! function definition

      COMMON/SET1/IYEAR,IBDOY,IEDOY,IFDAY,IFTIME

C Markdown flag.
      logical markdown
      common/markdownflag/markdown

C See ipvdata.h and ripvdat for IPV common block information.
      integer list
      dimension list(MCOM)  ! zone lists set to MCOM

      character outs*124,outs2*248
      character act*1
      character descr*7,descrst*10,descrfn*10,descr2st*8,descr2fn*8
      character zdescr*248   ! sized to hold 72 zone names
      character descra*7,descrb*7,descrc*7,descrd*7,descre*7

      if(act.eq.'a')then
        if(markdown)then
          write(outs,'(2A)')     'Title ',ipvtitl
          call edisp2tr(itru,outs)
          call edisp(itru,'synopsis ')
          call edisp248(itru,ipvsynop,120)
          write(outs,'(A,I3)')   'Images  ',nipvimg
          call edisp2tr(itru,outs)
        else
          write(outs,'(2A)')     ' Title :',ipvtitl
          call edisp(itru,outs)
          call edisp(itru,'synopsis:')
          call edisp248(itru,ipvsynop,80)
          write(outs,'(A,I3)')   'Images    : ',nipvimg
          call edisp(itru,outs)
        endif
        if(nipvimg.gt.0)then
          do i=1,nipvimg
            call edisp(itru,lipvimg(i))
          enddo
        endif
        call edisp(itru,'  ')
        if(markdown)then
          call edisp2tr(itru,'Simulations to be undertaken  ')
          call edisp(itru,' ')
          write(outs,'(a)') 
     &    'name    simulations                days  descripion'
          call edisp(itru,outs)
          write(outs,'(a)') 
     &    '------  -------------------------  ----  ----------'
          call edisp(itru,outs)
        else
          call edisp(itru,'Simulations to be undertaken:')
          write(outs,'(a,i2,a)') '____simulations (',nipvassmt,
     &    ')____________days descrip___'
          call edisp(itru,outs)
        endif
        if(nipvassmt.eq.0)then
          call edisp(itru,'Simulations not yet defined.')
        elseif(nipvassmt.eq.1)then
          call stdate(iyear,ipvastjd(1),descr,descrst,descr2st)
          call stdate(iyear,ipvafnjd(1),descr,descrfn,descr2fn)
          jjd1=(ipvafnjd(1)-ipvastjd(1))+1
          if(markdown)then
            write(outs,'(4a,i4,2a)') 'annual  ',descrst,' to ',
     &      descrfn,jjd1,'      ',ipvadesc(1)
          else
            write(outs,'(4a,i4,2a)') 'annual ',descrst,' to ',
     &      descrfn,jjd1,' ',ipvadesc(1)
          endif
          call edisp(itru,outs)
        elseif(nipvassmt.eq.3)then
          call stdate(iyear,ipvastjd(1),descr,descrst,descr2st)
          call stdate(iyear,ipvafnjd(1),descr,descrfn,descr2fn)
          jjd1=(ipvafnjd(1)-ipvastjd(1))+1
          if(markdown)then
            write(outs,'(4a,i4,2a)') 'winter  ',descrst,' to ',
     &      descrfn,jjd1,'     ',ipvadesc(1)
          else
            write(outs,'(4a,i4,2a)') 'winter ',descrst,' to ',
     &      descrfn,jjd1,' ',ipvadesc(1)
          endif
          call edisp(itru,outs)
          call stdate(iyear,ipvastjd(2),descr,descrst,descr2st)
          call stdate(iyear,ipvafnjd(2),descr,descrfn,descr2fn)
          jjd2=(ipvafnjd(2)-ipvastjd(2))+1
          if(markdown)then
            write(outs,'(4a,i4,2a)') 'trans   ',descrst,' to ',
     &      descrfn,jjd2,'     ',ipvadesc(2)
          else
            write(outs,'(4a,i4,2a)') 'trans  ',descrst,' to ',
     &      descrfn,jjd2,' ',ipvadesc(2)
          endif
          call edisp(itru,outs)
          call stdate(iyear,ipvastjd(3),descr,descrst,descr2st)
          call stdate(iyear,ipvafnjd(3),descr,descrfn,descr2fn)
          jjd3=(ipvafnjd(3)-ipvastjd(3))+1
          if(markdown)then
            write(outs,'(4a,i4,2a)')  'summer  ',descrst,' to ',
     &      descrfn,jjd3,'     ',ipvadesc(3)
          else
            write(outs,'(4a,i4,2a)')  'summer ',descrst,' to ',
     &      descrfn,jjd3,' ',ipvadesc(3)
          endif
          call edisp(itru,outs)
        elseif(nipvassmt.eq.5)then
          call stdate(iyear,ipvastjd(1),descr,descrst,descr2st)
          call stdate(iyear,ipvafnjd(1),descr,descrfn,descr2fn)
          jjd1=(ipvafnjd(1)-ipvastjd(1))+1
          if(markdown)then
            write(outs,'(4a,i4,2a)')  '1st win  ',descrst,' to ',
     &        descrfn,jjd1,'     ',ipvadesc(1)
          else
            write(outs,'(4a,i4,2a)')  '1st win ',descrst,' to ',
     &        descrfn,jjd1,' ',ipvadesc(1)
          endif
          call edisp(itru,outs)
          call stdate(iyear,ipvastjd(2),descr,descrst,descr2st)
          call stdate(iyear,ipvafnjd(2),descr,descrfn,descr2fn)
          jjd2=(ipvafnjd(2)-ipvastjd(2))+1
          if(markdown)then
            write(outs,'(4a,i4,2a)')  'spring   ',descrst,' to ',
     &      descrfn,jjd2,'     ',ipvadesc(2)
          else
            write(outs,'(4a,i4,2a)')  'spring  ',descrst,' to ',
     &      descrfn,jjd2,' ',ipvadesc(2)
          endif
          call edisp(itru,outs)
          call stdate(iyear,ipvastjd(3),descr,descrst,descr2st)
          call stdate(iyear,ipvafnjd(3),descr,descrfn,descr2fn)
          jjd3=(ipvafnjd(3)-ipvastjd(3))+1
          if(markdown)then
            write(outs,'(4a,i4,2a)')  'summer   ',descrst,' to ',
     &      descrfn,jjd3,'     ',ipvadesc(3)
          else
            write(outs,'(4a,i4,2a)')  'summer  ',descrst,' to ',
     &      descrfn,jjd3,' ',ipvadesc(3)
          endif
          call edisp(itru,outs)
          call stdate(iyear,ipvastjd(4),descr,descrst,descr2st)
          call stdate(iyear,ipvafnjd(4),descr,descrfn,descr2fn)
          jjd4=(ipvafnjd(4)-ipvastjd(4))+1
          if(markdown)then
            write(outs,'(4a,i4,2a)')  'autumn   ',descrst,' to ',
     &      descrfn,jjd4,'     ',ipvadesc(4)
          else
            write(outs,'(4a,i4,2a)')  'autumn  ',descrst,' to ',
     &      descrfn,jjd4,' ',ipvadesc(4)
          endif
          call edisp(itru,outs)
          call stdate(iyear,ipvastjd(5),descr,descrst,descr2st)
          call stdate(iyear,ipvafnjd(5),descr,descrfn,descr2fn)
          jjd5=(ipvafnjd(5)-ipvastjd(5))+1
          if(markdown)then
            write(outs,'(4a,i4,2a)')  '2nd win  ',descrst,' to ',
     &      descrfn,jjd5,'     ',ipvadesc(5)
          else
            write(outs,'(4a,i4,2a)')  '2nd win ',descrst,' to ',
     &      descrfn,jjd5,' ',ipvadesc(5)
          endif
          call edisp(itru,outs)
        endif

C Display days.
        if(nipvdispjd.gt.0)then
          call edisp(itru,'  ')
          do 155 ix=1,nipvdispjd
            call stdate(iyear,ipvdispjd(ix),descr,descrst,descr2st)
            write(outs,'(3a)')  'display day ',descrst,'.'
            call edisp(itru,outs)
  155     continue
        endif

        call edisp(itru,'  ')
        if(markdown)then
          call edisp2tr(itru,'Seasons (as found in climatelist)  ')
          call edisp(itru,' ')
          call edisp(itru,
     &    '         winter  spring  summer  autumn  winter')
          call edisp(itru,
     &    '-------  ------- ------  ------  ------  ------')
        else
          call edisp(itru,'Seasons (as found in climatelist file):')
          call edisp(itru,
     &    '         winter  spring  summer  autumn  winter')
        endif
        call stdate(iyear,is1wins,descra,descrst,descr2st)
        call stdate(iyear,is1sprs,descrb,descrst,descr2st)
        call stdate(iyear,is1sums,descrc,descrst,descr2st)
        call stdate(iyear,is2sprs,descrd,descrst,descr2st)
        call stdate(iyear,is2wins,descre,descrst,descr2st)
        write(outs,'(10a)') '  start  ',descra,' ',descrb,' ',descrc,
     &    ' ',descrd,' ',descre
        call edisp(itru,outs)
        call stdate(iyear,is1winf,descra,descrfn,descr2fn)
        call stdate(iyear,is1sprf,descrb,descrfn,descr2fn)
        call stdate(iyear,is1sumf,descrc,descrfn,descr2fn)
        call stdate(iyear,is2sprf,descrd,descrfn,descr2fn)
        call stdate(iyear,is2winf,descre,descrfn,descr2fn)
        write(outs,'(10a)') '  finish ',descra,' ',descrb,' ',descrc,
     &    ' ',descrd,' ',descre
        call edisp(itru,outs)
        jd1=(is1winf-is1wins)+1
        jd2=(is1sprf-is1sprs)+1
        jd3=(is1sumf-is1sums)+1
        jd4=(is2sprf-is2sprs)+1
        jd5=(is2winf-is2wins)+1
        write(outs,'(a,i4,1x,4i8)')'  days   ',jd1,jd2,jd3,jd4,jd5
        call edisp(itru,outs)

        call edisp(itru,'  ')
        if(nipvassmt.eq.1)then
          if(markdown)then
            call edisp(itru,'Ratios between assessments & all seasons')
            call edisp(itru,'  ')
            call edisp(itru,      '-------------  ---------  ')
          else
            call edisp(itru,'Ratios between assessments & all seasons:')
          endif 
          if(markdown)then
            WRITE(outs,'(a,f7.2)')  ' heating      ',ddmheat(1)
          else
            WRITE(outs,'(a,f7.2)')  ' heating    ',ddmheat(1)
          endif
          call edisp(itru,outs)
          if(markdown)then
            WRITE(outs,'(a,f7.2)')  ' cooling      ',ddmcool(1)
          else
            WRITE(outs,'(a,f7.2)')  ' cooling    ',ddmcool(1)
          endif
          call edisp(itru,outs)
          if(markdown)then
            WRITE(outs,'(a,f7.2)')  ' time-based     ',ddmtime(1)
          else
            WRITE(outs,'(a,f7.2)')  ' time-based ',ddmtime(1)
          endif
          call edisp(itru,outs)
        elseif(nipvassmt.eq.3)then
          if(markdown)then
            call edisp(itru,'Ratios between assessments')
            call edisp(itru,'')
            call edisp(itru,
     &      '              winter  transition  summer')
            call edisp(itru,
     &      '------------  ------  ----------  ------')
          else
            call edisp(itru,
     &      'Ratios between assessments and winter transition summer:')
          endif
          if(markdown)then
            WRITE(outs,'(a,3f8.2)') ' heating      ',ddmheat(1),
     &      ddmheat(2),ddmheat(3)
          else
            WRITE(outs,'(a,3f7.2)') ' heating    ',ddmheat(1),
     &      ddmheat(2),ddmheat(3)
          endif
          call edisp(itru,outs)
          if(markdown)then
            WRITE(outs,'(a,3f8.2)') ' cooling      ',ddmcool(1),
     &      ddmcool(2),ddmcool(3)
          else
            WRITE(outs,'(a,3f7.2)') ' cooling    ',ddmcool(1),
     &      ddmcool(2),ddmcool(3)
          endif
          call edisp(itru,outs)
          if(markdown)then
            WRITE(outs,'(a,3f8.2)') ' time-based   ',ddmtime(1),
     &      ddmtime(2),ddmtime(3)
          else
            WRITE(outs,'(a,3f7.2)') ' time-based ',ddmtime(1),
     &      ddmtime(2),ddmtime(3)
          endif
          call edisp(itru,outs)
        elseif(nipvassmt.eq.5)then
          if(markdown)then
            call edisp(itru,
     &      'Ratios for:   winter  spring summer autumn winter ')
            call edisp(itru,
     &      '------------  ------- ------ ------ ------ -------')
          else
            call edisp(itru,
     &      ' Ratios for:  winter spring summer autumn winter')
          endif
          if(markdown)then
            WRITE(outs,'(a,5f7.2)') ' heating      ',ddmheat(1),
     &      ddmheat(2),ddmheat(3),ddmheat(4),ddmheat(5)
          else
            WRITE(outs,'(a,5f7.2)') ' heating    ',ddmheat(1),
     &      ddmheat(2),ddmheat(3),ddmheat(4),ddmheat(5)
          endif
          call edisp(itru,outs)
          if(markdown)then
            WRITE(outs,'(a,5f7.2)') ' cooling      ',ddmcool(1),
     &      ddmcool(2),ddmcool(3),ddmcool(4),ddmcool(5)
          else
            WRITE(outs,'(a,5f7.2)') ' cooling    ',ddmcool(1),
     &      ddmcool(2),ddmcool(3),ddmcool(4),ddmcool(5)
          endif
          call edisp(itru,outs)
          if(markdown)then
            WRITE(outs,'(a,5f7.2)') ' time-based   ',ddmtime(1),
     &      ddmtime(2),ddmtime(3),ddmtime(4),ddmtime(5)
          else
            WRITE(outs,'(a,5f7.2)') ' time-based ',ddmtime(1),
     &      ddmtime(2),ddmtime(3),ddmtime(4),ddmtime(5)
          endif
          call edisp(itru,outs)
        endif
      endif
      if(act.eq.'a'.or.act.eq.'m')then

C List current metrics and demand sets.
        if(nms.ge.1)then
          call edisp(itru,' ')
          if(markdown)then
            write(outs2,'(2a)')
     &      'ID  Metric                zones  area    ',
     &      'type      name          zones'
            call edisp(itru,outs2)
            write(outs,'(2a)')
     &      '--  --------------------  -----  ------  ',
     &      '--------- ------------  ------------------------'
            call edisp(itru,outs)
          else
            call edisp(itru,
     &'Metric                 zones area  type   name')
          endif
          do 269 ij=1,nms
            lnd=lnblnk(msdoc(ij))
            lnz=lnblnk(metgroup(ij))
            lng=lnblnk(metrglbl(ij))

C If one of the metrics focused on a singe zone do not bother with
C the znarlist line. 
            if(msdoc(ij)(1:4).eq.'Guth'.or.
     &         msdoc(ij)(1:5).eq.'glare'.or.
     &         msdoc(ij)(1:12).eq.'daylightfact'.or.
     &         msdoc(ij)(1:4).eq.'DayF'.or.
     &         msdoc(ij)(1:11).eq.'contaminant')then
              if(markdown)then
                if(lng.lt.20) lng=20  ! to assist alignment
                if(lnd.lt.8) lnd=8
                if(lnz.lt.10) lnz=10
                write(outs2,'(i2,3a,i5,f9.2,4a)') ij,'  ',
     &            metrglbl(ij)(1:lng),' ',nzmg(ij),emgflr(ij),
     &            '  ',msdoc(ij)(1:lnd),'  ',
     &            metgroup(ij)(1:lnz)
                call edisp(itru,outs2)
              else
                if(lng.lt.20) lng=20  ! to assist alignment
                write(outs2,'(i2,3a,i4,f8.2,4a)') ij,' ',
     &            metrglbl(ij)(1:lng),' ',nzmg(ij),emgflr(ij),
     &            ' ',msdoc(ij)(1:lnd),' ',metgroup(ij)(1:lnz)
                call edisp(itru,outs2)
              endif
            else
              inlist=nzmg(ij)
              do 829 ijj=1,nzmg(ij)
                list(ijj)=izmg(ij,ijj)
 829          continue

C 3rd parameter needs to be MCOM.
              call znarlist(inlist,list,MCOM,zdescr,length,ier)
              if(markdown)then
                if(lng.lt.20) lng=20  ! to assist alignment
                if(lnd.lt.8) lnd=8
                if(lnz.lt.10) lnz=10
                write(outs2,'(i2,3a,i5,f9.2,6a)') ij,'  ',
     &            metrglbl(ij)(1:lng),' ',nzmg(ij),emgflr(ij),
     &            '  ',msdoc(ij)(1:lnd),'  ',
     &            metgroup(ij)(1:lnz),'  ',zdescr(7:lnblnk(zdescr))
                call edisp(itru,outs2)
              else
                if(lng.lt.20) lng=20  ! to assist alignment
                write(outs2,'(i2,3a,i4,f8.2,4a)') ij,' ',
     &            metrglbl(ij)(1:lng),' ',nzmg(ij),emgflr(ij),
     &            ' ',msdoc(ij)(1:lnd),' ',metgroup(ij)(1:lnz)
                call edisp(itru,outs2)
              endif
              if(.NOT.markdown) call edisp248(itru,zdescr,100)
            endif
 269      continue
        else
          call edisp(itru,'No metrics have been defined.')
        endif
      endif
      if(act.eq.'a'.or.act.eq.'d')then

C List current demand sets.
        if(neds.ge.1)then
          call edisp(itru,' ')
          if(markdown)then
            call edisp(itru,' ')
            write(outs,'(a)')
     &      'ID Name          Zones  Area     Scale   Notes   '
            call edisp(itru,outs)
            write(outs,'(a)')
     &      '-- ------------  ------ -------  ------  --------'
            call edisp(itru,outs)
          else
            call edisp(itru,
     &      'Demand set   zones  floor area  scaling factor')
          endif
          do 268 ij=1,neds
            inlist=nzedg(ij)
            do 828 ijj=1,nzedg(ij)
              list(ijj)=izedg(ij,ijj)
 828        continue
            call znarlist(inlist,list,MCOM,zdescr,length,ier)
            if(markdown)then
              write(outs2,'(i2,2a,i8,f8.2,f8.3,2a)') ij,' ',zedsdoc(ij),
     &          nzedg(ij),edgflr(ij),edgsca(ij),'  ',
     &          zdescr(7:lnblnk(zdescr))
              call edisp(itru,outs2)
            else
              write(outs,'(i2,2a,i4,f8.2,f6.3)') ij,' ',zedsdoc(ij),
     &          nzedg(ij),edgflr(ij),edgsca(ij)
              call edisp(itru,outs)
            endif
            if(.NOT.markdown) call edisp248(itru,zdescr,100)
 268      continue
        else
          call edisp(itru,'No demand sets have been defined.')
        endif
        call edisp(itru,' ')
        if(idmdinc(1).eq.1)then
          call edisp(itru,'Dispersed occupant demands included.')
        endif
        if(idmdinc(2).eq.1)then
          call edisp(itru,'Dispersed lighting demands included.')
        endif
        if(idmdinc(3).eq.1)then
          call edisp(itru,'Dispersed small power demands included.')
        endif
        if(idmdinc(4).eq.1)then
          call edisp(itru,'Dispersed fan demands included.')
        endif
        if(idmdinc(5).eq.1)then
          call edisp(itru,'Dispersed pump demands included.')
        endif
        if(idmdinc(6).eq.1)then
          call edisp(itru,'Dispersed lift demands included.')
        endif
        if(idmdinc(7).eq.1)then
          call edisp(itru,'Dispersed DHW demands included.')
        endif
      endif
      return
      end
