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


C RADCFG.F contains routines to control the exectution of RADIANCE from E2R.
C This is mainl done via the RADCFG file which contains data on the RADIANCE 
C rif file to be used and its intended purpose.

C RADCFGIN:  reads the contents of the radiance configuration file.
C RADCFGOUT: writes the contents of the radiance configuration file.
C whichrif:  presents a list of existing rif files.
C newscene:  asks questions supporting creation of a new rif file.
C ADDSCENE:  allows the creation of a default translation from ESP-r to radiance.
C setdfgrid: allows the creation of the daylight factor grid points.
C DFGRIDLBL: writes out an index of daylight grid point (act = g|G) or
C            the daylight factor (act = d|D).
C SENLOC:    creates a file containing the sensor locations and orientations
C checkrif:  checks for existing scenes or creates a new scene with the given 
C            parameters if none found.
C RDCFTMP:   Read in daylight coefficients to common and rewrites to casual 
C            gain control file.

C Logic that follows establishes whether both_esp_wave is true when
C reading in an existing file. Check ipairedscene & SCENEGSRC for 
C adjacent scenes.

C ********************* RADCFGIN *********************
C Reads contents of the rcf (Radiance configuration) file
C which holds meta data that cannot be held in the Radiance rif file.
      SUBROUTINE RADCFGIN(tfile)
#include "building.h"
#include "model.h"
#include "e2r_common.h"
#include "espriou.h"
      
      integer lnblnk  ! function definition
      
      character tfile*96  ! the radiance configuration file

      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      common/FILEP/IFIL
      common/rpath/path

      logical both_esp_wave      ! If true manage pair of scenes.
      logical focus_espg         ! If true editing in context of ESP-r geometry
      integer ipairedscene       ! Index of paired scene.
      COMMON/ESPWAVE/both_esp_wave,focus_espg,ipairedscene(16)

      common/user/browse

      character WORD*124, OUTS*124, OUTSTR*124
      character pfile*72,path*72,fs*1,datestr*24
      character oroot*28   ! Paired scene root name.
      dimension oroot(16)

      logical unixok
      logical browse
      logical newerfile
      
      IRCFG = IFIL+17

C Set folder separator (fs) to \ or / as required.
      call isunix(unixok)
      if(unixok)then
        fs = char(47)
      else
        fs = char(92)
      endif
      newerfile=.false.

C Clear oroot.
      do ij=1,16
        oroot(ij)='  '
      enddo

C Open the ESP-r radiance cfg (*.rcf) file.
      write(currentfile,'(a)') tfile(1:lnblnk(tfile))
      call ERPFREE(IRCFG,ISTAT)
      call EFOPSEQ(IRCFG,tfile,1,IER)
C      write(6,*) ' rcf status ',ier
C      write(6,*) ' rcf file ',LRADCF(1:lnblnk(LRADCF))
C      write(6,*) ' tfile file ',tfile(1:lnblnk(tfile))

C Read lines from file, header first.
      CALL STRIPC(IRCFG,OUTSTR,99,ND,1,'line 1',IER)
      IF(IER.NE.0)RETURN
      if (OUTSTR(1:13).ne.'*ESP-r visual') then
        call usrmsg('The specified file is not an ESP-r - ',
     &              'RADIANCE configuration file','F')
      endif

      CALL STRIPC(IRCFG,OUTSTR,99,ND,1,'date',IER)
      if(OUTSTR(1:5).eq.'*date')then
        K=6
        CALL EGETRM(OUTSTR,K,datestr,'W','date',IFLAG)
        write (outs,'(a,a)') 'Scenes last modified on ',
     &  datestr(1:lnblnk(datestr))
        newerfile=.true.
      else
        write (outs,'(a,a)') 'Scenes last modified on ',
     &    OUTSTR(1:lnblnk(OUTSTR))
      endif
      call edisp (iuout,' ')
      call edisp (iuout,outs)

C Read number of scenes contained if older file format.
      if(.NOT.newerfile)then
        CALL STRIPC(IRCFG,OUTSTR,99,ND,1,'scene nb line',IER)
        K=0
        CALL EGETWI(OUTSTR,K,NSCENE,1,16,'W','no of scenes',IER)
      
C For each scene read at least 5 items of data.
        IER=0
        do 10 I=1,NSCENE

C Reset the number of modelled blind states to zero.
          NBSRIF(I)=-1
 30       CALL STRIPC(IRCFG,OUTSTR,99,ND,1,'rif data',IER)
          if (OUTSTR(1:9).eq.'*sceneEND') goto 20
          if (IER.eq.0) then
            K=0
            CALL EGETW(OUTSTR,K,WORD,'W','data id',IFLAG)
            if (WORD(1:4).eq.'*snm') then

C SCENE is depreciated if found just continue
              continue
            elseif (WORD(1:4).eq.'*rnm') then

C First rif name corresponds to default blind state, subsquent rifs within 
C this scene represent alternate blind states.
              NBSRIF(I)=NBSRIF(I)+1
              if (NBSRIF(I).gt.NABS) then
                call usrmsg('Scene has too many rif files.','  ','W')
              endif
              if (NBSRIF(I).eq.0) then
                CALL EGETW(OUTSTR,K,RIFNAME(I),'W','rif name',IFLAG)
              else
                CALL EGETW(OUTSTR,K,LBSRIF(I),'W','bl st',IFLAG)
              endif
            elseif (WORD(1:4).eq.'*srt') then
              CALL EGETW(OUTSTR,K,SCENERT(I),'W','scene root',IFLAG)
            elseif (WORD(1:4).eq.'*rde') then
              CALL EGETRM(OUTSTR,K,SCENEDESC(I),'W','scene desc',IFLAG)
            elseif (WORD(1:4).eq.'*rpu') then
              CALL EGETW(OUTSTR,K,SCENEPURP(I),'W','scene purp',IFLAG)
            elseif (WORD(1:4).eq.'*rgd') then

C Note: there is one grid or sensor data file name per scene.
              CALL EGETW(OUTSTR,K,LDFGRID(I),'W','grid data',IFLAG)
            elseif (WORD(1:4).eq.'*slc') then
              CALL EGETW(OUTSTR,K,LDFGRID(I),'W','sensor data',IFLAG)
            elseif (WORD(1:5).eq.'*gsrc') then
              CALL EGETW(OUTSTR,K,SCENEGSRC(I),'W','geom source',IFLAG)
            elseif (WORD(1:5).eq.'*wave') then
              CALL EGETW(OUTSTR,K,SCENEWAVE(I),'W','wavefront',IFLAG)
            elseif (WORD(1:7).eq.'*wglaze') then
              CALL EGETW(OUTSTR,K,SCENEWGLZ(I),'W','wave glazing',IFLAG)
            elseif (WORD(1:6).eq.'*focus') then
              CALL EGETW(OUTSTR,K,SCENEFZNM(I),'W','focus zone',IFLAG)
            elseif (WORD(1:7).eq.'*mottle') then
              CALL EGETW(OUTSTR,K,SCENEMOTL(I),'W','mottle texture',
     &          IFLAG)
            else
              call usrmsg(' Unrecognised data prefix.',' ','W')
            endif
          endif
          goto 30
 20       continue
 10     continue
      else
      
C Look for *scene and then various tokens.  Loop till '*end'.
        IER=0
        NSCENE=0; I=0
 42     CALL STRIPC(IRCFG,OUTSTR,99,ND,1,'rcf data',IER)
        if (IER.ne.0) then
          call edisp(iuout,'Problem reading rcf line:')
          call edisp(iuout,outstr)
          return
        endif
        K=0
        CALL EGETW(OUTSTR,K,WORD,'W','data id',IFLAG)
        if (OUTSTR(1:6).eq.'*scene')then
          NSCENE=NSCENE+1
          I=NSCENE
          NBSRIF(I)=-1
          goto 42
        elseif (WORD(1:4).eq.'*snm') then

C SCENE is depreciated if found just continue
          continue
        elseif (WORD(1:4).eq.'*rnm') then

C First rif name corresponds to default blind state, subsquent rifs within 
C this scene represent alternate blind states.
          NBSRIF(I)=NBSRIF(I)+1
          if (NBSRIF(I).gt.NABS) then
            call usrmsg('Scene has too many rif files.','  ','W')
          endif
          if (NBSRIF(I).eq.0) then
            CALL EGETW(OUTSTR,K,RIFNAME(I),'W','rif name',IFLAG)
          else
            CALL EGETW(OUTSTR,K,LBSRIF(I),'W','bl st',IFLAG)
          endif
        elseif (WORD(1:4).eq.'*srt') then
          CALL EGETW(OUTSTR,K,SCENERT(I),'W','scene root',IFLAG)
        elseif (WORD(1:4).eq.'*rde') then
          CALL EGETRM(OUTSTR,K,SCENEDESC(I),'W','scene desc',IFLAG)
        elseif (WORD(1:4).eq.'*rpu') then
          CALL EGETW(OUTSTR,K,SCENEPURP(I),'W','scene purp',IFLAG)
        elseif (WORD(1:4).eq.'*rgd') then

C Note: there is one grid or sensor data file name per scene.
          CALL EGETW(OUTSTR,K,LDFGRID(I),'W','grid data',IFLAG)
        elseif (WORD(1:4).eq.'*gda') then
          CALL EGETW(OUTSTR,K,DFSURFNAME(I),'W','grid surf',IFLAG)
          CALL EGETWR(OUTSTR,K,WPDIST(I),0.,0.,'-','grid offset',IER)
          CALL EGETWI(OUTSTR,K,DFEDGE(I,1),1,99,'W','edge v start',IER)
          CALL EGETWI(OUTSTR,K,DFEDGE(I,2),1,99,'W','edge v end ',IER)
          CALL EGETWI(OUTSTR,K,igmajor(i),1,20,'W','major grid',IER)
          CALL EGETWI(OUTSTR,K,igminor(i),1,20,'W','minor grid',IER)
          CALL EGETWR(OUTSTR,K, CONV(I),0.,0.,'-','converg',IER)
          write(outs,'(2a,f6.3,4i4,f6.3)')'*gda ',
     &      DFSURFNAME(I)(1:lnblnk(DFSURFNAME(I))),WPDIST(i),
     &      DFEDGE(i,1),DFEDGE(i,2),igmajor(i),igminor(i),CONV(i)
          call edisp(iuout,outs)
          
        elseif (WORD(1:4).eq.'*slc') then
          CALL EGETW(OUTSTR,K,LDFGRID(I),'W','sensor data',IFLAG)
        elseif (WORD(1:5).eq.'*gsrc') then
          CALL EGETW(OUTSTR,K,SCENEGSRC(I),'W','geom source',IFLAG)
        elseif (WORD(1:5).eq.'*wave') then
          CALL EGETW(OUTSTR,K,SCENEWAVE(I),'W','wavefront',IFLAG)
        elseif (WORD(1:7).eq.'*wglaze') then
          CALL EGETW(OUTSTR,K,SCENEWGLZ(I),'W','wave glazing',IFLAG)
        elseif (WORD(1:6).eq.'*focus') then
          CALL EGETW(OUTSTR,K,SCENEFZNM(I),'W','focus zone',IFLAG)
        elseif (WORD(1:7).eq.'*mottle') then
          CALL EGETW(OUTSTR,K,SCENEMOTL(I),'W','mottle texture',IFLAG)
        elseif (WORD(1:7).eq.'*paired')then 
          CALL EGETW(OUTSTR,K,oroot(i),'W','other root name',IFLAG)
        elseif (WORD(1:4).eq.'*end') then
          goto 43
        else
          call usrmsg(' Unrecognised data prefix.',' ','W')
        endif
        goto 42  ! Read another line.
      endif

C Generate the ipairedscene array by comparing scene root names.
 43   do is=1,NSCENE
        lno=lnblnk(oroot(is))
        if(lno.gt.0)then
          do iso=1,NSCENE
            if(is.eq.iso)then
              continue    ! Skip if equal.
            else
              lnoo=lnblnk(oroot(iso))
              if(lnoo.gt.0)then
                if(oroot(is)(1:lno).eq.oroot(iso)(1:lnoo))then
                  ipairedscene(is)=iso
                endif
              endif
            endif
          enddo
        endif
      enddo
C      write(6,*) 'ipairedscene ',ipairedscene

      CALL ERPFREE(IRCFG,ISTAT)
      RETURN
      END

C ********************* RADCFGOUT *********************
C Writes the contents of the ESP-r radiance configuration file.
C Passed ICHK =3 update, =4 create. If focus zone needed user
C asked to confirm.
      SUBROUTINE RADCFGOUT(ICHK)
#include "building.h"
#include "model.h"
#include "geometry.h"
#include "e2r_common.h"
#include "help.h"
      
      integer lnblnk  ! function definition

      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      common/FILEP/IFIL

      integer ncomp,ncon
      COMMON/C1/NCOMP,NCON
      dimension IVALS(MCOM)

      logical both_esp_wave      ! If true manage pair of scenes.
      logical focus_espg         ! If true editing in context of ESP-r geometry
      integer ipairedscene       ! Index of paired scene.
      COMMON/ESPWAVE/both_esp_wave,focus_espg,ipairedscene(16)

      character OUTS*72
      character dstmp*24,pfile*72,msg*72
      character tokens*124,comment*72,aligned_str*124

      helpinsub='radcfg'  ! set for subroutine

C Ask if a default image should be created.
      helptopic='select_focus_zone'
      call gethelptext(helpinsub,helptopic,nbhelp)

C Check that ICHK is a sensible value. =3 updates =4 creates new.
      IRCFG = IFIL+17
      if (ICHK.ne.3.and.ICHK.ne.4) ICHK=4

C Get current time.
      call dstamp(dstmp)

      IML=10
      do I=1,NSCENE   ! Find the longest string.
        ISLENG=LNBLNK(RIFNAME(I))+5
        IF(ISLENG.GT.IML)IML=ISLENG
        ISLENG=LNBLNK(LBSRIF(I))+5
        IF(ISLENG.GT.IML)IML=ISLENG
        ISLENG=LNBLNK(SCENERT(I))+5
        IF(ISLENG.GT.IML)IML=ISLENG
        ISLENG=LNBLNK(SCENEDESC(I))+5
        IF(ISLENG.GT.IML)IML=ISLENG
        ISLENG=LNBLNK(SCENEPURP(I))+5
        IF(ISLENG.GT.IML)IML=ISLENG
        ISLENG=LNBLNK(SCENEGSRC(I))+6
        IF(ISLENG.GT.IML)IML=ISLENG
        ISLENG=LNBLNK(SCENEWAVE(I))+6
        IF(ISLENG.GT.IML)IML=ISLENG
        ISLENG=LNBLNK(SCENEWGLZ(I))+7
        IF(ISLENG.GT.IML)IML=ISLENG
      enddo

C Set file name including path and open.
      write (pfile,'(a,a)') runpath(1:lnrp),LRADCF(1:lnblnk(LRADCF))
      call ERPFREE(IRCFG,ISTAT)
      call FPOPEN(IRCFG,ISTAT,1,ICHK,pfile)

C Write header first. *date signals a newer file format.
      write(IRCFG,'(a)')'*ESP-r visualisation scene configuration file.'
      write(tokens,'(2a)') '*date ',dstmp
      write(IRCFG,'(a)') tokens(1:lnblnk(tokens))

C For each scene write 1) root identifier, 2) short description 
C and 3) intended use/ purpose.
      do 10 I=1,NSCENE
        write(IRCFG,'(a)') '*scene'
        write(tokens,'(2a)') '*rnm ',RIFNAME(I)(1:lnblnk(RIFNAME(I)))
        write(comment,'(a)') 'Radiance input file.'
        call align_comment(IML+3,tokens,comment,aligned_str)
        if(I.eq.1)then
          write(IRCFG,'(a)') aligned_str(1:lnblnk(aligned_str))
        else
          write(IRCFG,'(a)') tokens(1:lnblnk(tokens))
        endif
        if (NBSRIF(I).gt.0) then
          write(tokens,'(2a)') '*rnm ',LBSRIF(I)(1:lnblnk(LBSRIF(I)))
          write(comment,'(a)') 'Alternative Radiance input file.'
          call align_comment(IML+3,tokens,comment,aligned_str)
          write(IRCFG,'(a)') aligned_str(1:lnblnk(aligned_str))
        endif
        write(tokens,'(2a)') '*srt ',SCENERT(I)(1:lnblnk(SCENERT(I)))
        write(comment,'(a)') 'Scene root name.'
        call align_comment(IML+3,tokens,comment,aligned_str)
        if(I.eq.1)then
          write(IRCFG,'(a)') aligned_str(1:lnblnk(aligned_str))
        else
          write(IRCFG,'(a)') tokens(1:lnblnk(tokens))
        endif

        write(tokens,'(2a)')'*rde ',SCENEDESC(I)(1:lnblnk(SCENEDESC(I)))
        write(comment,'(a)') 'Scene description.'
        call align_comment(IML+3,tokens,comment,aligned_str)
        if(I.eq.1)then
          write(IRCFG,'(a)') aligned_str(1:lnblnk(aligned_str))
        else
          write(IRCFG,'(a)') tokens(1:lnblnk(tokens))
        endif

        write(tokens,'(2a)')'*rpu ',SCENEPURP(I)(1:lnblnk(SCENEPURP(I)))
        write(comment,'(a)') 'Intended use.'
        call align_comment(IML+3,tokens,comment,aligned_str)
        if(I.eq.1)then
          write(IRCFG,'(a)') aligned_str(1:lnblnk(aligned_str))
        else
          write(IRCFG,'(a)') tokens(1:lnblnk(tokens))
        endif

        write(tokens,'(2a)')'*gsrc ',
     &    SCENEGSRC(I)(1:lnblnk(SCENEGSRC(I)))
        write(comment,'(a)') 'Source of geometry.'
        call align_comment(IML+3,tokens,comment,aligned_str)
        if(I.eq.1)then
          write(IRCFG,'(a)') aligned_str(1:lnblnk(aligned_str))
        else
          write(IRCFG,'(a)') tokens(1:lnblnk(tokens))
        endif

        if(SCENEGSRC(I)(1:5).eq.'ESP-r')then
          write(tokens,'(a)')'*wave UNKNOWN'
        else
          write(tokens,'(2a)')'*wave ',
     &      SCENEWAVE(I)(1:lnblnk(SCENEWAVE(I)))
        endif
        write(comment,'(a)') 'Wavefront source file.'
        call align_comment(IML+3,tokens,comment,aligned_str)
        if(I.eq.1)then
          write(IRCFG,'(a)') aligned_str(1:lnblnk(aligned_str))
        else
          write(IRCFG,'(a)') tokens(1:lnblnk(tokens))
        endif

        write(tokens,'(2a)')'*wglaze ',
     &    SCENEWGLZ(I)(1:lnblnk(SCENEWGLZ(I)))
        write(comment,'(a)') 'Wavefront glazing entities.'
        call align_comment(IML+3,tokens,comment,aligned_str)
        if(I.eq.1)then
          write(IRCFG,'(a)') aligned_str(1:lnblnk(aligned_str))
        else
          write(IRCFG,'(a)') tokens(1:lnblnk(tokens))
        endif

        write(tokens,'(2a)')'*mottle ',
     &    SCENEMOTL(I)(1:lnblnk(SCENEMOTL(I)))
        write(comment,'(a)') 'Surface texture.'
        call align_comment(IML+3,tokens,comment,aligned_str)
        if(I.eq.1)then
          write(IRCFG,'(a)') aligned_str(1:lnblnk(aligned_str))
        else
          write(IRCFG,'(a)') tokens(1:lnblnk(tokens))
        endif

C If one of the views needing a focus zone and SCENEFZNM is UNKNOWN
C and ifocz is not a positive number confirm with user.
        if(SCENEFZNM(I)(1:7).eq.'UNKNOWN')then
          if(SCENEPURP(I)(1:8).eq.'External')then
            SCENEFZNM(I)='all'
        elseif(SCENEPURP(I)(1:8).eq.'Internal'.or.
     &         SCENEPURP(I)(1:5).eq.'Glare'.or.
     &         SCENEPURP(I)(1:8).eq.'Day_fact')then
          write(msg,'(a,i2,3a,i4)') 'For scene ',i,' ',
     &      SCENEPURP(I),' the current focus zone is ',ifocz
          call edisp(iuout,msg)
 30       INPIC=1
          CALL EPICKS(INPIC,IVALS,'Expecting a focus zone',
     &      'please confirm: ',
     &      12,NCOMP,zname,' zone list',IER,nbhelp)
          if (INPIC.eq.0) goto 30
          ifocz=IVALS(1)
          write(SCENEFZNM(I),'(a)') zname(ifocz)
        endif

        endif
        write(tokens,'(2a)')'*focus ',
     &    SCENEFZNM(I)(1:lnblnk(SCENEFZNM(I)))
        write(comment,'(a)') 'Focus zone.'
        call align_comment(IML+3,tokens,comment,aligned_str)
        if(I.eq.1)then
          write(IRCFG,'(a)') aligned_str(1:lnblnk(aligned_str))
        else
          write(IRCFG,'(a)') tokens(1:lnblnk(tokens))
        endif

        if (SCENEPURP(I)(1:8).eq.'Day_fact') then
          write(tokens,'(2a)')'*rgd ',LDFGRID(I)(1:lnblnk(LDFGRID(I)))
          write(comment,'(a)') 'Daylight factor grid.'
          call align_comment(IML+3,tokens,comment,aligned_str)
          write(IRCFG,'(a)') aligned_str(1:lnblnk(aligned_str))

          write(tokens,'(2a,f6.3,4i4,f6.3)')'*gda ',
     &      DFSURFNAME(I)(1:lnblnk(DFSURFNAME(I))),WPDIST(i),
     &      DFEDGE(i,1),DFEDGE(i,2),igmajor(i),igminor(i),CONV(i)
          write(comment,'(a)') 'grid attributes'
          call align_comment(IML+3,tokens,comment,aligned_str)
          write(IRCFG,'(a)') aligned_str(1:lnblnk(aligned_str))

        elseif (SCENEPURP(I)(1:7).eq.'Day_lux'.or.
     &          SCENEPURP(I)(1:5).eq.'Illum') then
          write(tokens,'(2a)')'*rgd ',LDFGRID(I)(1:lnblnk(LDFGRID(I)))
          write(comment,'(a)') 'Lux grid.'
          call align_comment(IML+3,tokens,comment,aligned_str)
          write(IRCFG,'(a)') aligned_str(1:lnblnk(aligned_str))

          write(tokens,'(2a,f6.3,4i4,f6.3)')'*gda ',
     &      DFSURFNAME(I)(1:lnblnk(DFSURFNAME(I))),WPDIST(i),
     &      DFEDGE(i,1),DFEDGE(i,2),igmajor(i),igminor(i),CONV(i)
          write(comment,'(a)') 'grid attributes'
          call align_comment(IML+3,tokens,comment,aligned_str)
          write(IRCFG,'(a)') aligned_str(1:lnblnk(aligned_str))

        elseif (SCENEPURP(I)(1:8).eq.'Coupling') then
          write(tokens,'(2a)')'*slc ',LDFGRID(I)(1:lnblnk(LDFGRID(I)))
          write(comment,'(a)') 'Sensor loactions.'
          call align_comment(IML+3,tokens,comment,aligned_str)
          write(IRCFG,'(a)') aligned_str(1:lnblnk(aligned_str))
        elseif (SCENEPURP(I)(1:8).eq.'Day_coef') then
          write(tokens,'(2a)')'*slc ',LDFGRID(I)(1:lnblnk(LDFGRID(I)))
          write(comment,'(a)') 'Sensor loactions.'
          call align_comment(IML+3,tokens,comment,aligned_str)
          write(IRCFG,'(a)') aligned_str(1:lnblnk(aligned_str))
        elseif(ipairedscene(I).gt.0)then
          iother=ipairedscene(I)
          write(tokens,'(2a)')'*paired ',
     &      SCENERT(iother)(1:lnblnk(SCENERT(iother)))
          write(comment,'(a)') 'Paired with.'
          call align_comment(IML+3,tokens,comment,aligned_str)
          write(IRCFG,'(a)') aligned_str(1:lnblnk(aligned_str))
        endif
C        write(IRCFG,'(a)')'*sceneEND'
 10   continue

      write(IRCFG,'(a)') '*end'
      CALL ERPFREE(IRCFG,ISTAT)

      RETURN
      END

C ********************* whichrif *********************
C Presents a list of existing scenes pointing to specific
C rif files.
C If MODE is 's' then create a static list, 'm' then allow modification.

      SUBROUTINE whichrif(MODE)
#include "building.h"
#include "model.h"
#include "geometry.h"
#include "epara.h"
#include "e2r_common.h"
#include "help.h"

      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      COMMON/FILEP/IFIL

      integer ncomp,ncon
      COMMON/C1/NCOMP,NCON

      logical both_esp_wave      ! If true manage pair of scenes.
      logical focus_espg         ! If true editing in context of ESP-r geometry
      integer ipairedscene       ! Index of paired scene.
      COMMON/ESPWAVE/both_esp_wave,focus_espg,ipairedscene(16)
  
      character priorpurpose*12
      common/prior/priorpurpose
      
      character MODE*1, VERT*30, KEY*1
      character ltmp*72,outs*124,outse*124,outsw*124,APP*6
      character ldfg*72   ! local file name
      character doit*248
      dimension VERT(35)
      dimension IVALS(MCOM)

      logical SELECT,delzone,cpyzone,modscene,clkok
      logical unixok
      integer IW  ! for radio button
      integer MVERT,IVERT  ! max items and current menu item

      helpinsub='radcfg'  ! set for subroutine

 1    IER=0
      SELECT=.false.
      delzone=.false.
      cpyzone=.false.
      modscene=.false.
      ITA2=IFIL+7
      ITA3=IFIL+10

C Initialise scene menu variables based on window size. 
C IVERT is the menu position, MVERT the current number of menu lines.
C Use shorter menu depending on what MOD is.
      SELECT=.FALSE.
      MHEAD=0
      IF(MODE.EQ.'M'.OR.MODE.EQ.'m')THEN
        MCTL=4
      ELSE
        MCTL=3
      endif
      ILEN=NSCENE
      IPACT=CREATE
      CALL EKPAGE(IPACT)

C Initial menu entry setup.
   92 IER=0
      IVERT=-3

C Loop through the items until the page to be displayed. M is the 
C current menu line index. Build up text strings for the menu. 
    3 M=MHEAD
      DO 10 L=1,ILEN
        IF(L.GE.IST.AND.(L.LE.(IST+MIFULL)))THEN
          M=M+1
          CALL EMKEY(L,KEY,IER)
          WRITE(VERT(M),14)KEY,SCENERT(L)
   14     FORMAT(A1,1X,A)
        ENDIF
   10 CONTINUE

C Number of actual items displayed.
      MVERT=M+MCTL

C If a long list include page facility text and info on portion seen.      
      IF(IPFLG.EQ.0)THEN
        VERT(M+1)=  '  ----------------------'
      ELSE
        WRITE(VERT(M+1),15)IPM,MPM 
   15   FORMAT ('0 page part: ',I1,' -- of:',I1)
      ENDIF

C If MODE has been passed as a M then add an item at the end allowing
C list to be modified otherwise blank choice.  If user is browsing
C no need to show additional options.
      IF(MODE.EQ.'M'.OR.MODE.EQ.'m')THEN
        VERT(M+2)  ='* add/delete/variant scene'
        VERT(M+3)  ='? help                '
        VERT(M+4)  ='- exit menu           '
      ELSE
        VERT(M+2)  ='? help                '
        VERT(M+3)  ='- exit menu           '
      ENDIF

C Display the menu.
      CALL EMENU('Scenes',VERT,MVERT,IVERT)
      IF(IVERT.LE.MHEAD)THEN
        IVERT=-1
        goto 3
      ELSEIF(IVERT.EQ.MVERT)THEN

C If no selection has been made before exit then return.
        IF(.NOT.SELECT)then
          ICHK=4
          if(modscene)call RADCFGOUT(ICHK)
          RETURN
        endif
      ELSEIF(IVERT.EQ.(MVERT-1))THEN

C Produce help text.
        helptopic='scene_selection_list'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL PHELPD('scene section',nbhelp,'-',0,0,IER)
      ELSEIF(IVERT.GT.MHEAD.AND.IVERT.LT.(MVERT-MCTL+1))THEN

C Decode from the potentially long list of scenes via KEYIND.
C If delete or copy scene previously selected then new selection
C should be deleted or copied, respectively.
        CALL KEYIND(MVERT,IVERT,IFOC,IO)
        SELECT=.TRUE.
        ISCENE=IFOC
        if(delzone)then   ! act on prior request to delete
          do 1110 I=IFOC,NSCENE
            if (I.eq.NSCENE) then
              RIFNAME(I)='UNKNOWN'
              SCENERT(I)='UNKNOWN'
              SCENEDESC(I)='UNKNOWN'
              SCENEPURP(I)='UNKNOWN'
              LBSRIF(I)='NONE'
              SCENEGSRC(I)='UNKNOWN'
              SCENEWAVE(I)='UNKNOWN'
              SCENEWGLZ(I)='UNKNOWN'
              SCENEFZNM(I)='UNKNOWN'
              SCENEMOTL(I)='medium'
              LDFGRID(I)='UNKNOWN'
              DFSURFNAME(I)='UNKNOWN'
              WPDIST(I)=0.0
              DFEDGE(I,1)=-1; DFEDGE(I,2)=-1
              igmajor(I)=2
              igminor(I)=2
              CONV(I)=0.250
            else
              RIFNAME(I)=RIFNAME(I+1)
              SCENERT(I)=SCENERT(I+1)
              SCENEDESC(I)=SCENEDESC(I+1)
              SCENEPURP(I)=SCENEPURP(I+1)
              LBSRIF(I)=LBSRIF(I+1)
              SCENEGSRC(I)=SCENEGSRC(I+1)
              SCENEWAVE(I)=SCENEWAVE(I+1)
              SCENEWGLZ(I)=SCENEWGLZ(I+1)
              SCENEFZNM(I)=SCENEFZNM(I+1)
              SCENEMOTL(I)=SCENEMOTL(I+1)
              LDFGRID(I)=LDFGRID(I+1)
              DFSURFNAME(I)=DFSURFNAME(I+1)
              WPDIST(I)=WPDIST(I+1)
              DFEDGE(I,1)=DFEDGE(I+1,1)
              DFEDGE(I,2)=DFEDGE(I+1,2)
              igmajor(I)=igmajor(I+1)
              igminor(I)=igminor(I+1)
              CONV(I)=CONV(I+1)
            endif
 1110     continue
          NSCENE=NSCENE-1
          ISCENE=1
          modscene=.true.
        elseif(cpyzone)then  ! act on previous request to make variant scene

C Create a variant of an existing scene. 
          if(NSCENE+1.gt.14)then
            call usrmsg('Cannot add another scene (limit is 16).',
     &        ' ','W')
            goto 1
          endif

C Rescan the source rif fie. Reset both_esp_wave and focus_espg
C as required, also reset logicals about scene topics prior to
C scanning the source rif file.
          rofil='UNKNOWN'; matfil='UNKNOWN'
          outdone=.false.; rzfil='UNKNOWN'
          rmfil='UNKNOWN'; rmmfil='UNKNOWN'
          glzfil='UNKNOWN'; aglzfil='UNKNOWN'
          indone=.false.; rskyfil='UNKNOWN'
          skydone=.false.; misdone=.false.; vewdone=.false.
          glzdone=.false.; iesdone=.false.
          octfil='UNKNOWN'; picfil='UNKNOWN'
          outdone=.false.;  indone=.false.

          call RRIF(ita2,ita3,RIFNAME(IFOC),'s',ier)
          if(SCENEGSRC(IFOC)(1:5).eq.'ESP-r')then
            both_esp_wave=.false.
            focus_espg=.true.
          else
            both_esp_wave=.false.
            focus_espg=.false.
          endif

C Request a 6 char phrase to apply to the root name and file names.
C Depending on what the user wants to alter in the variant a number of
C the scene attributes e.g. SCENEPURP, SCENEGSRC, SCENEWAVE, SCENEMOTL
C will not be altered.
C Possible variants:
C - if External can't think of any variants to implement.
C - if Internal then could focus on a different zone & viewpoint.
C - if DF then might focus on a different zone or different grid density
C - if Lux or Illuminance might focus on a different zone or grid density.
C - if Glare then might focus on a different zone and viewpoint.
          if(SCENEPURP(IFOC)(1:8).eq.'External')then
            call usrmsg('Source scene is external no variant options',
     &        'are available.','W')
            goto 1
          elseif(SCENEPURP(IFOC)(1:8).eq.'Internal'.or.
     &           SCENEPURP(IFOC)(1:5).eq.'Glare')then
            helptopic='pick_focus_zone'
            call gethelptext(helpinsub,helptopic,nbhelp)
 30         INPIC=1
            CALL EPICKS(INPIC,IVALS,
     &        'Different focus for an internal or glare study',
     &        'Pick a focus zone: ',
     &        12,NCOMP,zname,' zone list',IER,nbhelp)
            if (INPIC.eq.0) goto 30
            ifocz=IVALS(1)

          elseif(SCENEPURP(IFOC)(1:8).eq.'Day_fact'.or.
     &           SCENEPURP(IFOC)(1:7).eq.'Day_lux'.or.
     &           SCENEPURP(IFOC)(1:5).eq.'Illum')then
            helptopic='pick_focus_zone'
            call gethelptext(helpinsub,helptopic,nbhelp)
 31         INPIC=1
            CALL EPICKS(INPIC,IVALS,
     &        'Different focus zone for grid-based assessments.',
     &        'Pick a focus zone: ',
     &        12,NCOMP,zname,' zone list',IER,nbhelp)
            if (INPIC.eq.0) goto 31
            ifocz=IVALS(1)
          else
            call usrmsg('Sorry no variant options available for the',
     &        'source scene current purpose.','W')
            goto 1
          endif

C Ask for short phrase to mark this particular variant.
          APP='_xyz  '
          ISTRW=6
          CALL EASKSCMD(APP,'  ',
     &    'Specify up to 6 characters to append to names?',
     &    'cancel',clkok,ISTRW,'_xyz','root name modifier',IER,nbhelp)
          call usrmsg(' ',' ','-')
          if(clkok) goto 1         ! User selected cancel.

          NSCENE=NSCENE+1          ! Increment scene index.

C Compose a new root name for the scene depending on whether ESP-r geo
C or Wavefront geometry. outse & outsw are then used to build up the
C various file names.
          IL=min(lnblnk(cfgroot),21) 
          isp=lnblnk(SCENEPURP(IFOC))
          lna=lnblnk(APP)
          if(SCENEPURP(IFOC)(1:8).eq.'Internal'.or.
     &       SCENEPURP(IFOC)(1:5).eq.'Glare'.or.
     &       SCENEPURP(IFOC)(1:8).eq.'Day_fact'.or.
     &       SCENEPURP(IFOC)(1:7).eq.'Day_lux'.or.
     &       SCENEPURP(IFOC)(1:5).eq.'Illum')then
            lnzn=lnblnk(zname(ifocz))
            write (outse,'(4a)') zname(ifocz)(1:lnzn),'_',
     &       SCENEPURP(IFOC)(1:isp),APP(1:lna)

C As Wavefront files tend to be named obj_... append APP.
            if(SCENEPURP(ifoc)(1:8).eq.'External')then
              write(outsw,'(2a)') 'obj_ex',APP(1:lna)
            elseif(SCENEPURP(ifoc)(1:8).eq.'Internal')then
              write(outsw,'(2a)') 'obj_in',APP(1:lna)
            elseif(SCENEPURP(ifoc)(1:5).eq.'Glare')then
              write(outsw,'(2a)') 'obj_gl',APP(1:lna)
            elseif(SCENEPURP(ifoc)(1:8).eq.'Day_fact') then
              write(outsw,'(2a)') 'obj_df',APP(1:lna)
            else
              write(outsw,'(2a)') 'objn',APP(1:lna)
            endif
            if(focus_espg)then
              write (SCENERT(NSCENE),'(a)') outse(1:28)
            elseif(.NOT.focus_espg)then
              write (SCENERT(NSCENE),'(a)') outsw(1:28)
            endif
          else
            write (outse,'(4a)') cfgroot(1:IL),'_',
     &       SCENEPURP(IFOC)(1:isp),APP(1:lna)
            if(SCENEPURP(ifoc)(1:8).eq.'External')then
              write(outsw,'(3a)') 'obj_ex_',APP(1:lna),'_W'
            elseif(SCENEPURP(ifoc)(1:8).eq.'Internal')then
              write(outsw,'(3a)') 'obj_in_',APP(1:lna),'_W'
            elseif(SCENEPURP(ifoc)(1:5).eq.'Glare')then
              write(outsw,'(3a)') 'obj_gl_',APP(1:lna),'_W'
            elseif(SCENEPURP(ifoc)(1:8).eq.'Day_fact') then
              write(outsw,'(3a)') 'obj_df_',APP(1:lna),'_W'
            else
              write(outsw,'(3a)') 'obj_',APP(1:lna),'_W'
            endif
            if(focus_espg)then
              write (SCENERT(NSCENE),'(a)') outse(1:28)
            elseif(.NOT.focus_espg)then
              write (SCENERT(NSCENE),'(a)') outsw(1:28)
            endif
          endif

          call isunix(unixok)

C The variant rif file will be generated once all of the internal
C attributes have been sorted.
          write (RIFNAME(NSCENE),'(2a)') 
     &       SCENERT(NSCENE)(1:lnblnk(SCENERT(NSCENE))),'.rif'

          lnr=lnblnk(SCENEDESC(IFOC))
          write(SCENEDESC(NSCENE),'(2a)') 
     &      SCENEDESC(IFOC)(1:lnr),APP(1:lna)
          SCENEPURP(NSCENE)=SCENEPURP(IFOC)
          LBSRIF(NSCENE)=LBSRIF(IFOC)
          SCENEGSRC(NSCENE)=SCENEGSRC(IFOC)
          SCENEWAVE(NSCENE)=SCENEWAVE(IFOC)
          SCENEWGLZ(NSCENE)=SCENEWGLZ(IFOC)
          if(SCENEPURP(IFOC)(1:8).eq.'Internal'.or.
     &       SCENEPURP(IFOC)(1:5).eq.'Glare'.or.
     &       SCENEPURP(IFOC)(1:8).eq.'Day_fact'.or.
     &       SCENEPURP(IFOC)(1:7).eq.'Day_lux'.or.
     &       SCENEPURP(IFOC)(1:5).eq.'Illum')then
            if(ifocz.gt.0)then
              write(SCENEFZNM(NSCENE),'(a)') zname(ifocz)
            else
              SCENEFZNM(NSCENE)=SCENEFZNM(IFOC)
            endif
	  else
            SCENEFZNM(NSCENE)=SCENEFZNM(IFOC)
          endif
          SCENEMOTL(NSCENE)=SCENEMOTL(IFOC)
	  
          lnrne=lnblnk(outse)
          lnrnw=lnblnk(outsw)

C For scenes which rely on a grid file setup and then copy.
          if(SCENEPURP(IFOC)(1:8).eq.'Day_fact'.or.
     &       SCENEPURP(IFOC)(1:7).eq.'Day_lux'.or.
     &       SCENEPURP(IFOC)(1:5).eq.'Illum')then
            if(focus_espg)then
              write (ldfg,'(2a)') outse(1:lnrne),'.grd'
            elseif(.NOT.focus_espg)then
              write (ldfg,'(2a)') outsw(1:lnrnw),'.grd'
            endif
            LDFGRID(NSCENE)=ldfg
            if(unixok)then
              WRITE(doit,'(6A)') 'cd ',runpath(1:lnblnk(runpath)),
     &          '; cp ',
     &          LDFGRID(IFOC)(1:LNBLNK(LDFGRID(IFOC))),' ',
     &          LDFGRID(NSCENE)(1:LNBLNK(LDFGRID(NSCENE)))
            else
              WRITE(doit,'(6A)')  'cd ',runpath(1:lnblnk(runpath)),
     &          '; copy /y ',
     &          LDFGRID(IFOC)(1:LNBLNK(LDFGRID(IFOC))),' ',
     &          LDFGRID(NSCENE)(1:LNBLNK(LDFGRID(NSCENE)))
            endif
            CALL USRMSG('Copying grd file via:',doit,'-')
            CALL RUNIT(doit,'-')
            call sleep(1)
            call usrmsg('Remember to re-define grid for new focus',
     &        'zone e.g. pick reference surface and edge.','W')
          else
            LDFGRID(NSCENE)=LDFGRID(IFOC)
          endif
          DFSURFNAME(NSCENE)=DFSURFNAME(IFOC)
          WPDIST(NSCENE)=WPDIST(IFOC)
          DFEDGE(NSCENE,1)=DFEDGE(IFOC,1)
          DFEDGE(NSCENE,2)=DFEDGE(IFOC,2)
          igmajor(NSCENE)=igmajor(IFOC)
          igminor(NSCENE)=igminor(IFOC)
          CONV(NSCENE)=CONV(IFOC)
          ISCENE=NSCENE

C Adapt the scene names within the rif file. Some will then be
C filled by Radiance and others will need to be copied. Start
C with the sky file.
          if(focus_espg)then
            write (ldfg,'(2a)') outse(1:lnrne),'.sky'
          elseif(.NOT.focus_espg)then
            write (ldfg,'(2a)') outsw(1:lnrnw),'.sky'
          endif
          if(unixok)then
            WRITE(doit,'(6A)') 'cd ',runpath(1:lnblnk(runpath)),
     &        '; cp ',
     &        rskyfil(1:LNBLNK(rskyfil)),' ',ldfg(1:LNBLNK(ldfg))
          else
            WRITE(doit,'(6A)')  'cd ',runpath(1:lnblnk(runpath)),
     &        '; copy /y ',
     &        rskyfil(1:LNBLNK(rskyfil)),' ',ldfg(1:LNBLNK(ldfg))
          endif
          CALL USRMSG('Copying sky file via:',doit,'-')
          CALL RUNIT(doit,'-')
          call sleep(1)
          rskyfil=ldfg  ! update the sky file with variant name.

C Radiance glazing file.  Note: for Wavefront geometry it also includes
C the glazing materials so a separate glazing file is not needed.
          if(glzfil(1:7).eq.'UNKNOWN')then
            continue
          else
            if(focus_espg)then
              write (ldfg,'(2a)') outse(1:lnrne),'-glz.rad'
              if(unixok)then
                WRITE(doit,'(6A)') 'cd ',runpath(1:lnblnk(runpath)),
     &          '; cp ',
     &          glzfil(1:LNBLNK(glzfil)),' ',ldfg(1:LNBLNK(ldfg))
              else
                WRITE(doit,'(6A)')  'cd ',runpath(1:lnblnk(runpath)),
     &          '; copy /y ',
     &          glzfil(1:LNBLNK(glzfil)),' ',ldfg(1:LNBLNK(ldfg))
              endif
              CALL USRMSG('Copying glazing file via:',doit,'-')
              CALL RUNIT(doit,'-')
              call sleep(1)
              glzfil=ldfg  ! update the glazing file with variant name.
            elseif(.NOT.focus_espg)then
              write (ldfg,'(2a)') outsw(1:lnrnw),'-glz.rad'
              if(unixok)then
                WRITE(doit,'(6A)') 'cd ',runpath(1:lnblnk(runpath)),
     &          '; cp ',
     &          glzfil(1:LNBLNK(glzfil)),' ',ldfg(1:LNBLNK(ldfg))
              else
                WRITE(doit,'(6A)')  'cd ',runpath(1:lnblnk(runpath)),
     &          '; copy /y ',
     &          glzfil(1:LNBLNK(glzfil)),' ',ldfg(1:LNBLNK(ldfg))
              endif
              CALL USRMSG('Copying WF glazing file via:',doit,'-')
              CALL RUNIT(doit,'-')
              call sleep(1)
              glzfil=ldfg  ! update the glazing file with variant name.
              continue
            endif
          endif

C Radiance material file (which also includes glazing properties 
C in the case of a Wavefront source geometry).
          if(matfil(1:7).eq.'UNKNOWN')then
            continue
          else
            if(focus_espg)then
              write (ldfg,'(2a)') outse(1:lnrne),'.mat'
            elseif(.NOT.focus_espg)then
              write (ldfg,'(2a)') outsw(1:lnrnw),'.mat'
            endif
            if(unixok)then
              WRITE(doit,'(6A)') 'cd ',runpath(1:lnblnk(runpath)),
     &          '; cp ',
     &          matfil(1:LNBLNK(matfil)),' ',ldfg(1:LNBLNK(ldfg))
            else
              WRITE(doit,'(6A)')  'cd ',runpath(1:lnblnk(runpath)),
     &          '; copy /y ',
     &          matfil(1:LNBLNK(matfil)),' ',ldfg(1:LNBLNK(ldfg))
            endif
            CALL USRMSG('Copying material file via:',doit,'-')
            CALL RUNIT(doit,'-')
            call sleep(1)
            matfil=ldfg  ! update the material file with variant name.
          endif

C Radiance inside composition file for ESP-r geometry sources.
          if(rzfil(1:7).eq.'UNKNOWN')then
            continue
          else
            if(focus_espg)then
              write (ldfg,'(2a)') outse(1:lnrne),'-in.rad'
              if(unixok)then
                WRITE(doit,'(6A)') 'cd ',runpath(1:lnblnk(runpath)),
     &            '; cp ',
     &            rzfil(1:LNBLNK(rzfil)),' ',ldfg(1:LNBLNK(ldfg))
              else
                WRITE(doit,'(6A)')  'cd ',runpath(1:lnblnk(runpath)),
     &            '; copy /y ',
     &            rzfil(1:LNBLNK(rzfil)),' ',ldfg(1:LNBLNK(ldfg))
              endif
              CALL USRMSG('Copying Radiance inside comp file via:',
     &          doit,'-')
              CALL RUNIT(doit,'-')
              call sleep(1)
              rzfil=ldfg  ! update the inside comp file with variant name.
            elseif(.NOT.focus_espg)then
              continue
            endif
          endif

C Radiance outside composition file or Wavefront composition file.
          if(rofil(1:7).eq.'UNKNOWN')then
            continue
          else
            if(focus_espg)then
              write (ldfg,'(2a)') outse(1:lnrne),'-out.rad'
            elseif(.NOT.focus_espg)then
              write (ldfg,'(2a)') outsw(1:lnrnw),'.rad'
            endif
            if(unixok)then
              WRITE(doit,'(6A)') 'cd ',runpath(1:lnblnk(runpath)),
     &          '; cp ',
     &          rofil(1:LNBLNK(rofil)),' ',ldfg(1:LNBLNK(ldfg))
            else
              WRITE(doit,'(6A)')  'cd ',runpath(1:lnblnk(runpath)),
     &          '; copy /y ',
     &          rofil(1:LNBLNK(rofil)),' ',ldfg(1:LNBLNK(ldfg))
            endif
            CALL USRMSG('Copying Radiance outside comp file via:',
     &        doit,'-')
            CALL RUNIT(doit,'-')
            call sleep(1)
            rofil=ldfg  ! update the outside comp file with variant name.
          endif

C The empty files which Radiance will fill.
          if(focus_espg)then
            write(picfil,'(a)') outse(1:lnblnk(outse)) 
	    write(optnfil,'(2a)')outse(1:lnblnk(outse)),'.opt'
            write(octfil,'(2a)') outse(1:lnblnk(outse)),'.oct'

          elseif(.NOT.focus_espg)then
            write (picfil,'(a)') outsw(1:lnblnk(outsw))
	    write(optnfil,'(2a)')outsw(1:lnblnk(outsw)),'.opt'
            write(octfil,'(2a)') outsw(1:lnblnk(outsw)),'.oct'
          endif
          CALL USRMSG('  ','  ','-')

C Write updates to the new rif file.
          ITA2 = IFIL+7
          call mkrif(ita2,RIFNAME(NSCENE),4,0,'-',ier)

C Update the list of scenes.
          ICHK=4
          call radcfgout(ICHK)

C Check new scenes data by calling newscene.  As newscene adds one to NSCENE 
C subtract one just now, no other data will be effected. It assumes that
C the user will provide appropriate alternative names.
          call edisp(iuout,' ')
          write(outs,'(2a)') 
     &      'Next - change to the variant scene and adapt as required.'
          call edisp(iuout,outs)
          modscene=.true.
        endif

C If a zone has been deleted then redisplay menu, 
C foremost to reset ISCENE.
        if (delzone) then
          goto 1
        else

C Re-scan rif file if not a new zone.
          ITA2 = IFIL+7
          ITA3 = IFIL+10
          call RRIF(ita2,ita3,RIFNAME(ISCENE),'r',ier)
          RETURN
        endif
      endif

C Scene list manipulation.
      IF(IVERT.EQ.(MVERT-2))THEN

C Allow a scene to be added or deleted.
        IF(MODE.EQ.'M'.OR.MODE.EQ.'m')THEN
          helptopic='scene_selection_list'
          call gethelptext(helpinsub,helptopic,nbhelp)
          IW=1
          CALL EASKMBOX(' ','Action: ','add scene','delete scene',
     &      'make varant of scene','via loading a rif file','cancel',
     &      ' ',' ',' ',IW,nbhelp)
          if(IW.EQ.1)then

C Add a new scene (new rif file etc).
            delzone = .false.

            CALL EASKMBOX('Create scene:',' ','via ESP-r geometry)',
     &        'via both ESP-r & Wavefront geometry','via a rif file',
     &        'cancel',' ',' ',' ',' ',IWDV,nbhelp)
            if (IWDV.eq.1) then
              both_esp_wave=.false.; focus_espg=.true.
              call newscene('-')            ! interactive creation
              cmdact='Create'
              call autorad
              SCENEGSRC(ISCENE)='ESP-r'
              ipairedscene(ISCENE)=0        ! Not paired.
              ICHK=4
              call radcfgout(ICHK)
            elseif (IWDV.eq.2) then
              both_esp_wave=.true.; focus_espg=.true.
              call newscene('-')            ! Interactive creation.
              cmdact='Create'
              call autorad
              SCENEGSRC(ISCENE)='ESP-r'
              ipairedscene(ISCENE)=ISCENE+1 ! Paired with +1.
              ipairedscene(ISCENE+1)=ISCENE ! Reciprical
              SCENEGSRC(ipairedscene(ISCENE))='Wavefront'
              ICHK=4
              call radcfgout(ICHK)
            elseif (IWDV.eq.3) then

C Add scene from an existing rif file (i.e. a wavefront related rif).
              goto 42
            else
              goto 92
            endif
            RETURN
          elseif(IW.EQ.2)then

C Delete a scene.
            delzone = .true.
            call usrmsg(' ','Select scene to delete from list.','-')
            IVERT=-2
            goto 3
          elseif(IW.EQ.3)then

C Copy a scene
            cpyzone = .true.
            call usrmsg(' ',
     &        'Select scene to make variant of from list.','-')
            IVERT=-2
            goto 3
          elseif(IW.EQ.4)then
            goto 42
          endif

        ELSE

C If there are enough items allow paging control via EKPAGE.
          IF(IPFLG.EQ.1)THEN
            IPACT=EDIT
            CALL EKPAGE(IPACT)
          ENDIF
        ENDIF
      
      else

C Not one of the legal menu choices.
        IVERT=-1
        goto 92
      ENDIF
      IVERT=-2
      goto 3


C Add scene from an existing rif file (i.e. a wavefront related rif).
  42  if(NSCENE+1.gt.14)then
        call usrmsg('Cannot add another scene (limit is 16).',
     &    ' ','W')
        return
      endif
 642  helptopic='e2r_RIF_file_name'
      call gethelptext(helpinsub,helptopic,nbhelp)
      lnrp=lnblnk(runpath)  ! remember this 
      write(ltmp,'(a)') runpath(1:lnrp)
      CALL EASKS(ltmp,'Radiance RIF file?',
     &  '  ',72,'xxx.rif','rif file name',IER,nbhelp)
      if(ltmp.eq.' ')goto 642

C Set new scene to focus scene.
      NSCENE = NSCENE+1
      ISCENE=NSCENE

      RIFNAME(ISCENE) = ltmp
      ITA2 = IFIL+7
      ITA3 = IFIL+10
      call RRIF(ita2,ita3,RIFNAME(ISCENE),'s',ier)

C Some details are sensitive to the geometric source.
      CALL EASKMBOX('Source of geometry:',' ',
     &  'Generated by ESP-r','Based on wavefront file',' ',
     &  ' ',' ',' ',' ',' ',IWS,nbhelp)
      if(IWS.eq.1)then

C When importing from a user supplied rif file which uses ESP-r
C geometry there may or may not be a Wavefront version scene.
C Ask the user.
        write(SCENEGSRC(ISCENE),'(a)') 'ESP-r'
        focus_espg=.true.
        IW=1
        CALL EASKMBOX('Options for equivalent wavefront scene:',
     &    ' ','no wavefront scene',
     &    'one of the current scenes','will add later',
     &    ' ',' ',' ',' ',' ',IW,nbhelp)
        if(IW.eq.1)then
          both_esp_wave=.false.
        elseif(IW.eq.2)then
          both_esp_wave=.true.
        elseif(IW.eq.3)then
          both_esp_wave=.false.
        endif
      elseif(IWS.eq.2)then

C When importing from a user supplied rif file which is based
C on Wavefront geometry there may or may not be an ESP-r version
C of the scene.  Ask the user....
        write(SCENEGSRC(ISCENE),'(a)') 'Wavefront'
        focus_espg=.false.
        IW=1
        CALL EASKMBOX('Options for equivalent ESP-r geo scene:',
     &    ' ','no ESP-r geo scene',
     &    'one of the current scenes','will add later',
     &    ' ',' ',' ',' ',' ',IW,nbhelp)
        if(IW.eq.1)then
          both_esp_wave=.false.
        elseif(IW.eq.2)then
          both_esp_wave=.true.
        elseif(IW.eq.3)then
          both_esp_wave=.false.
        endif
      endif

C Ask user to confirm use case hint from the name of the file.
      call edisp(iuout,' ')      
      call edisp(iuout,
     &  'Please re-confirm the original purpose of that rif file.')      
      helptopic='purpose_of_rif_scene'
      call gethelptext(helpinsub,helptopic,nbhelp)
      ino=0
      idno=1
      call MENUATOL('Select a scene purpose.','Scene purpose',
     &  'a External images','b Internal images',
     &  'c Glare sources','d Daylight factors','e Coupling',
     &  'f Daylight coefficients','g Daylight lux @points',
     &  'h Illuminance @points','i Night external',' ',' ',
     &  ' ',ino,idno,nbhelp)
      IRP=INO

C Ask focus zone if internal, glare, coupling, daylight coefficients,
C lux @ points or illuminance @ points.
      if (IRP.eq.2.or.IRP.eq.3.or.IRP.eq.5.or.IRP.eq.6.or.
     &    IRP.eq.7.or.IRP.eq.8) then
        call edisp(iuout,' ')      
        call edisp(iuout,
     &  'Please re-confirm the original focus zone of that rif file.')      
        helptopic='pick_focus_zone'
        call gethelptext(helpinsub,helptopic,nbhelp)
 33     INPIC=1
        CALL EPICKS(INPIC,IVALS,' ','Pick a focus zone: ',
     &    12,NCOMP,zname,' zone list',IER,nbhelp)
        if (INPIC.eq.0) goto 33
        ifocz=IVALS(1)

C If focused on ESP-r geometry and also managing Wavefront then update
C the appropriate SCENEFZNM.
        if(focus_espg.and.both_esp_wave)then
          write(SCENEFZNM(ISCENE),'(a)') zname(ifocz)   ! ESP-r scene
          write(SCENEFZNM(ISCENE+1),'(a)') zname(ifocz) ! Wavefront scene
        elseif(.NOT.focus_espg.and.both_esp_wave)then
          write(SCENEFZNM(ISCENE),'(a)') zname(ifocz)   ! Wavefront scene
          write(SCENEFZNM(ISCENE-1),'(a)') zname(ifocz) ! ESP-r scene
        elseif(focus_espg.and.(.NOT.both_esp_wave))then
          write(SCENEFZNM(ISCENE),'(a)') zname(ifocz)   ! ESP-r scene
        elseif(.NOT.focus_espg.and.(.NOT.both_esp_wave))then
          write(SCENEFZNM(ISCENE),'(a)') zname(ifocz)  ! ISCENE is Wavefront
        endif
      endif
      if (IRP.eq.1) then 
        if(focus_espg.and.both_esp_wave)then
          SCENEPURP(ISCENE)='External'
          SCENEPURP(ISCENE+1)='External'
        elseif(.NOT.focus_espg.and.both_esp_wave)then
          SCENEPURP(ISCENE)='External'
          SCENEPURP(ISCENE-1)='External'
        elseif(focus_espg.and.(.NOT.both_esp_wave))then
          SCENEPURP(ISCENE)='External'
        elseif(.NOT.focus_espg.and.(.NOT.both_esp_wave))then
          SCENEPURP(ISCENE)='External'
        endif
      elseif (IRP.eq.2) then 
        if(focus_espg.and.both_esp_wave)then
          SCENEPURP(ISCENE)='Internal'
          SCENEPURP(ISCENE+1)='Internal'
        elseif(.NOT.focus_espg.and.both_esp_wave)then
          SCENEPURP(ISCENE)='Internal'
          SCENEPURP(ISCENE-1)='Internal'
        elseif(focus_espg.and.(.NOT.both_esp_wave))then
          SCENEPURP(ISCENE)='Internal'
        elseif(.NOT.focus_espg.and.(.NOT.both_esp_wave))then
          SCENEPURP(ISCENE)='Internal'
        endif
      elseif (IRP.eq.3) then 
        if(focus_espg.and.both_esp_wave)then
          SCENEPURP(ISCENE)='Glare'
          SCENEPURP(ISCENE+1)='Glare'
        elseif(.NOT.focus_espg.and.both_esp_wave)then
          SCENEPURP(ISCENE)='Glare'
          SCENEPURP(ISCENE-1)='Glare'
        elseif(focus_espg.and.(.NOT.both_esp_wave))then
          SCENEPURP(ISCENE)='Glare'
        elseif(.NOT.focus_espg.and.(.NOT.both_esp_wave))then
          SCENEPURP(ISCENE)='Glare'
        endif
      elseif (IRP.eq.4) then 
        if(focus_espg.and.both_esp_wave)then
          SCENEPURP(ISCENE)='Day_fact'
          SCENEPURP(ISCENE+1)='Day_fact'
        elseif(.NOT.focus_espg.and.both_esp_wave)then
          SCENEPURP(ISCENE)='Day_fact'
          SCENEPURP(ISCENE-1)='Day_fact'
        elseif(focus_espg.and.(.NOT.both_esp_wave))then
          SCENEPURP(ISCENE)='Day_fact'
        elseif(.NOT.focus_espg.and.(.NOT.both_esp_wave))then
          SCENEPURP(ISCENE)='Day_fact'
        endif
      elseif (IRP.eq.5) then 
        if(focus_espg.and.both_esp_wave)then
          SCENEPURP(ISCENE)='Coupling'
          SCENEPURP(ISCENE+1)='Coupling'
        elseif(.NOT.focus_espg.and.both_esp_wave)then
          SCENEPURP(ISCENE)='Coupling'
          SCENEPURP(ISCENE-1)='Coupling'
        elseif(focus_espg.and.(.NOT.both_esp_wave))then
          SCENEPURP(ISCENE)='Coupling'
        elseif(.NOT.focus_espg.and.(.NOT.both_esp_wave))then
          SCENEPURP(ISCENE)='Coupling'
        endif
      elseif (IRP.eq.6) then 
        if(focus_espg.and.both_esp_wave)then
          SCENEPURP(ISCENE)='Day_coef'
          SCENEPURP(ISCENE+1)='Day_coef'
        elseif(.NOT.focus_espg.and.both_esp_wave)then
          SCENEPURP(ISCENE)='Day_coef'
          SCENEPURP(ISCENE-1)='Day_coef'
        elseif(focus_espg.and.(.NOT.both_esp_wave))then
          SCENEPURP(ISCENE)='Day_coef'
        elseif(.NOT.focus_espg.and.(.NOT.both_esp_wave))then
          SCENEPURP(ISCENE)='Day_coef'
        endif
      elseif (IRP.eq.7) then 
        if(focus_espg.and.both_esp_wave)then
          SCENEPURP(ISCENE)='Day_lux'
          SCENEPURP(ISCENE+1)='Day_lux'
        elseif(.NOT.focus_espg.and.both_esp_wave)then
          SCENEPURP(ISCENE)='Day_lux'
          SCENEPURP(ISCENE-1)='Day_lux'
        elseif(focus_espg.and.(.NOT.both_esp_wave))then
          SCENEPURP(ISCENE)='Day_lux'
        elseif(.NOT.focus_espg.and.(.NOT.both_esp_wave))then
          SCENEPURP(ISCENE)='Day_lux'
        endif
      elseif (IRP.eq.8) then 
        if(focus_espg.and.both_esp_wave)then
          SCENEPURP(ISCENE)='Illum'
          SCENEPURP(ISCENE+1)='Illum'
        elseif(.NOT.focus_espg.and.both_esp_wave)then
          SCENEPURP(ISCENE)='Illum'
          SCENEPURP(ISCENE-1)='Illum'
        elseif(focus_espg.and.(.NOT.both_esp_wave))then
          SCENEPURP(ISCENE)='Illum'
        elseif(.NOT.focus_espg.and.(.NOT.both_esp_wave))then
          SCENEPURP(ISCENE)='Illum'
        endif
      elseif (IRP.eq.9) then 
        if(focus_espg.and.both_esp_wave)then
          SCENEPURP(ISCENE)='Night_ext'
          SCENEPURP(ISCENE+1)='Night_ext'
        elseif(.NOT.focus_espg.and.both_esp_wave)then
          SCENEPURP(ISCENE)='Night_ext'
          SCENEPURP(ISCENE-1)='Night_ext'
        elseif(focus_espg.and.(.NOT.both_esp_wave))then
          SCENEPURP(ISCENE)='Night_ext'
        elseif(.NOT.focus_espg.and.(.NOT.both_esp_wave))then
          SCENEPURP(ISCENE)='Night_ext'
        endif
      endif

C Set default names based on above.  If wavefront source
C adjust the root name to avoid duplication.
      if (IRP.eq.2.or.IRP.eq.3.or.IRP.eq.5.or.IRP.eq.6.or.
     &    IRP.eq.7.or.IRP.eq.8) then
        lnfzn=lnblnk(zname(ifocz))
        write (outse,'(3a)') zname(ifocz)(1:lnfzn),
     &    '_',SCENEPURP(ISCENE)(1:6)

        write (outsw,'(4a)') zname(ifocz)(1:lnfzn),
     &    '_',SCENEPURP(ISCENE)(1:6),'_W'
        if(focus_espg.and.both_esp_wave)then
          write (SCENEDESC(ISCENE),'(5a)') 'Scene of ',
     &      zname(ifocz)(1:lnfzn),' for ',
     &      SCENEPURP(ISCENE)(1:lnblnk(SCENEPURP(ISCENE))),
     &      ' analysis.'
          write (SCENEDESC(ISCENE+1),'(5a)') 'Scene of ',
     &      zname(ifocz)(1:lnfzn),' for ',
     &      SCENEPURP(ISCENE)(1:lnblnk(SCENEPURP(ISCENE))),
     &      ' analysis _W.'
        elseif(.NOT.focus_espg.and.both_esp_wave)then
          write (SCENEDESC(ISCENE),'(5a)') 'Scene of ',
     &      zname(ifocz)(1:lnfzn),' for ',
     &      SCENEPURP(ISCENE)(1:lnblnk(SCENEPURP(ISCENE))),
     &      ' analysis _W.'
          write (SCENEDESC(ISCENE-1),'(5a)') 'Scene of ',
     &      zname(ifocz)(1:lnfzn),' for ',
     &      SCENEPURP(ISCENE)(1:lnblnk(SCENEPURP(ISCENE))),
     &      ' analysis.'
        elseif(focus_espg.and.(.NOT.both_esp_wave))then
          write (SCENEDESC(ISCENE),'(5a)') 'Scene of ',
     &      zname(ifocz)(1:lnfzn),' for ',
     &      SCENEPURP(ISCENE)(1:lnblnk(SCENEPURP(ISCENE))),
     &      ' analysis.'
        elseif(.NOT.focus_espg.and.(.NOT.both_esp_wave))then
          write (SCENEDESC(ISCENE),'(5a)') 'Scene of ',
     &      zname(ifocz)(1:lnfzn),' for ',
     &      SCENEPURP(ISCENE)(1:lnblnk(SCENEPURP(ISCENE))),
     &      ' analysis _W.'
        endif
      else
        IL=min(lnblnk(cfgroot),21) 
        write (outse,'(3a)') cfgroot(1:IL),'_',
     &    SCENEPURP(ISCENE)(1:6)
        write (outsw,'(4a)') cfgroot(1:IL),'_',
     &    SCENEPURP(ISCENE)(1:6),'_W'
        if(focus_espg.and.both_esp_wave)then
          write (SCENEDESC(ISCENE),'(5a)') 'Scene of ',
     &      cfgroot(1:lnblnk(cfgroot)),' for ',
     &      SCENEPURP(ISCENE)(1:lnblnk(SCENEPURP(ISCENE))),
     &      ' analysis.'
          write (SCENEDESC(ISCENE+1),'(5a)') 'Scene of ',
     &      cfgroot(1:lnblnk(cfgroot)),' for ',
     &      SCENEPURP(ISCENE)(1:lnblnk(SCENEPURP(ISCENE))),
     &      ' analysis _W.'
        elseif(.NOT.focus_espg.and.both_esp_wave)then
          write (SCENEDESC(ISCENE),'(5a)') 'Scene of ',
     &      cfgroot(1:lnblnk(cfgroot)),' for ',
     &      SCENEPURP(ISCENE)(1:lnblnk(SCENEPURP(ISCENE))),
     &      ' analysis _W.'
          write (SCENEDESC(ISCENE-1),'(5a)') 'Scene of ',
     &      cfgroot(1:lnblnk(cfgroot)),' for ',
     &      SCENEPURP(ISCENE)(1:lnblnk(SCENEPURP(ISCENE))),
     &      ' analysis.'
        elseif(focus_espg.and.(.NOT.both_esp_wave))then
          write (SCENEDESC(ISCENE),'(5a)') 'Scene of ',
     &      cfgroot(1:lnblnk(cfgroot)),' for ',
     &      SCENEPURP(ISCENE)(1:lnblnk(SCENEPURP(ISCENE))),
     &      ' analysis.'
        elseif(.NOT.focus_espg.and.(.NOT.both_esp_wave))then
          write (SCENEDESC(ISCENE),'(5a)') 'Scene of ',
     &      cfgroot(1:lnblnk(cfgroot)),' for ',
     &      SCENEPURP(ISCENE)(1:lnblnk(SCENEPURP(ISCENE))),
     &      ' analysis _W.'
        endif
      endif
      if(focus_espg.and.both_esp_wave)then
        write (SCENERT(ISCENE),'(a)') outse(1:28)
        write (SCENERT(ISCENE+1),'(a)') outsw(1:28)
      elseif(.NOT.focus_espg.and.both_esp_wave)then
        write (SCENERT(ISCENE),'(a)') outsw(1:28)
        write (SCENERT(ISCENE-1),'(a)') outse(1:28)
      elseif(focus_espg.and.(.NOT.both_esp_wave))then
        write (SCENERT(ISCENE),'(a)') outse(1:28)
      elseif(.NOT.focus_espg.and.(.NOT.both_esp_wave))then
        write (SCENERT(ISCENE),'(a)') outsw(1:28)
      endif

C Additional details for specific purposes.
      if (IRP.eq.4.or.IRP.eq.7.or.IRP.eq.8) then

C Set up daylight factor or illuminance calc points.
        if(SCENEFZNM(ISCENE)(1:7).eq.'UNKNOWN')then
          continue
        else
          write(outs,'(2a)') 'Current focus is ',SCENEFZNM(ISCENE)
          call edisp(iuout,' ')
          call edisp(iuout,outs)
        endif
        call setdfgrid(izone,ier)

      elseif (IRP.eq.5.or.IRP.eq.6) then 

C Get sensor information from cgc file and create a slc file.
         call SENLOC
      endif

      CALL EASKMBOX('Source of geometry:',' ',
     &  'Generated by ESP-r','Based on wavefront file',' ',
     &  ' ',' ',' ',' ',' ',IWS,nbhelp)
      if(IWS.eq.1)then
        focus_espg=.true.
        write(SCENEGSRC(ISCENE),'(a)') 'ESP-r'
        ICHK=4
        call radcfgout(ICHK)
      elseif(IWS.eq.2)then
        focus_espg=.false.
        write(SCENEGSRC(ISCENE),'(a)') 'Wavefront'
 643    helptopic='e2r_RIF_file_name'
        call gethelptext(helpinsub,helptopic,nbhelp)
        lnrp=lnblnk(runpath)  ! remember this 
        write(ltmp,'(a)') runpath(1:lnrp)
        CALL EASKS(ltmp,'Wavefront file?',
     &    '  ',72,'xxx.rif','wavefront file name',IER,nbhelp)
        if(ltmp.eq.' ')goto 643
        SCENEWAVE(ISCENE)=ltmp
        ICHK=4
        call radcfgout(ICHK)
      endif
      return

      END


C ********************* newscene *********************
C Presents choices for a new rif file or sets
C values if MODE='s' (silent) from command line.
C if MODE='-' ask purposes and confirm root and file names.
C if MODE='s' command line purpose not matched so make one to match.
C Parent subroutine will already have set both_esp_wave and focus_espg.

      SUBROUTINE newscene(MODE)
#include "building.h"
#include "model.h"
#include "geometry.h"
#include "e2r_common.h"
#include "help.h"

      integer lnblnk  ! function definition

      integer ncomp,ncon
      common/C1/NCOMP,NCON

      logical both_esp_wave      ! If true manage pair of scenes.
      logical focus_espg         ! If true editing in context of ESP-r geometry
      integer ipairedscene       ! Index of paired scene.
      COMMON/ESPWAVE/both_esp_wave,focus_espg,ipairedscene(16)

      COMMON/gzonpik/izgfoc,nzg,nznog(mcom)

      character priorpurpose*12   ! use when copying a scene.
      common/prior/priorpurpose

      dimension IVALS(MCOM)

      character MODE*1,outs*124,outse*124,outsw*124
      character dname*28,tname*28
      character dname72*72,tname72*72
      character rifwave*72
      logical paired,fict,purpchanged

      helpinsub='radcfg'  ! set for subroutine
      outse='  '; outsw='  '

C If e2r passed a purpose but not found then checkrif will
C call newscene with mode = 's' so as to use the passed
C purpose when building the files.
      if (MODE.eq.'s') then
        IZONE=0
        IRP=0
        if (aim(1:8).eq.'External') then
          IRP=1
          izone=-1
          zone='FRED'
        elseif (aim(1:8).eq.'Internal') then
          IRP=2
          call zindex(zone,izone)
        elseif (aim(1:5).eq.'Glare') then
          IRP=3
          call zindex(zone,izone)
        elseif (aim(1:8).eq.'Day_fact') then
          IRP=4
          call zindex(zone,izone)
        elseif (aim(1:8).eq.'Coupling') then
          IRP=5
          call zindex(zone,izone)
        elseif (aim(1:8).eq.'Day_coef') then
          IRP=6
          call zindex(zone,izone)
        elseif (aim(1:7).eq.'Day_lux') then
          IRP=7
          call zindex(zone,izone)
        elseif (aim(1:5).eq.'Illum') then
          IRP=8
          call zindex(zone,izone)
        elseif (aim(1:9).eq.'Night_ext') then
          IRP=9
          izone=-1
          zone='FRED'
        endif
      else

C Mode is not s then confirm purpose.
        helptopic='purpose_of_rif_scene'
        call gethelptext(helpinsub,helptopic,nbhelp)
        ino=0
        idno=1
        call MENUATOL('Select a scene purpose.',' Scene purpose',
     &    'a External images','b Internal images','c Glare sources',
     &    'd Daylight factors','e Coupling','f Daylight coefficients',
     &    'g Daylight lux @points','h Illuminance @points',
     &    'i Night external',' ',' ',' ',ino,idno,nbhelp)
        IRP=INO
      endif

C If no purpose then exit.
      if(IRP.eq.0)return

C Check for tmc file
      call inqu_tmc

C Reset glazing treatment to default.
      iglzty=1

      NSCENE=NSCENE+1      ! Increment scene counter.

C Set new scene to focus scene.
      ISCENE=NSCENE        ! Set new scene as current.

C If user creating both ESP-r and Wavefront variants set
C additional pointers.
      paired=.false.
      iother=0
      if(both_esp_wave)then 
        if(focus_espg)then
          ipairedscene(ISCENE)=ISCENE+1 ! Pair with Wavefront variant.
          NSCENE = NSCENE+1             ! Increment for Wavefront variant.
          ipairedscene(NSCENE)=ISCENE   ! Link back to ESP-r variant.
          paired=.true.; iother= ipairedscene(ISCENE)
        else

C << Ask user which ESP-r geometric scene.
          call usrmsg('Logic to be added to link Wavefront to',
     &      'an ESP-r scene.','W')
        endif
      endif

C Reset the number of chosen zones (nzg) to NCOMP.
      nzg=NCOMP
      izgfoc=0
      do 3 I=1,NCOMP
        nznog(I)=I
 3    continue

C Now ask for zone for internal and glare purposes.
      if (IRP.eq.2.or.IRP.eq.3.or.IRP.eq.5.or.IRP.eq.6.or.
     &    IRP.eq.7.or.IRP.eq.8) then
        if (IZONE.le.0) then
          helptopic='pick_focus_zone'
          call gethelptext(helpinsub,helptopic,nbhelp)
 30       INPIC=1
          CALL EPICKS(INPIC,IVALS,' ','Pick a focus zone: ',
     &      12,NCOMP,zname,' zone list',IER,nbhelp)
          if (INPIC.eq.0) goto 30
          ifocz=IVALS(1)
        else
          ifocz=IZONE
        endif
        write(SCENEFZNM(ISCENE),'(a)') zname(ifocz)
        if(paired) write(SCENEFZNM(iother),'(a)') zname(ifocz)
      endif
      
C Set assumptions for each scene purpose and gather any additional 
C information.
      if (IRP.eq.1) then 
        SCENEPURP(ISCENE)='External'
        intext = 0   ! external view
        iadobe = 2   ! medium adobe texture
        SCENEMOTL(ISCENE)='medium'
        imgqua = 'Medium'
        detlvl = 'Medium'
        llvar  = 'High'
        indrcb = 1
        ipicx  = 500
        penumb = 'False'
        crenrp = '-av 10 10 10 '
        coconv = '  '
        ifocz = -1
        write(SCENEFZNM(ISCENE),'(a)') 'all'
	if(paired)then
          SCENEPURP(iother)='External'
          SCENEMOTL(iother)='medium'
          write(SCENEFZNM(iother),'(a)') 'all'
        endif
      elseif (IRP.eq.2) then 

C Internal image selected.
        SCENEPURP(ISCENE)='Internal'
        intext = 1  ! internal view
        iadobe = 3  ! fine adobe texture
        SCENEMOTL(ISCENE)='fine'
        imgqua = 'Medium'
        detlvl = 'Medium'
        llvar  = 'High'
        indrcb = 1
        ipicx  = 500
        penumb = 'True'
        crenrp = '-av 1 1 1 '
        coconv = '  '
	if(paired)then
          SCENEPURP(iother)='Internal'
          SCENEMOTL(iother)='fine'
        endif
      elseif (IRP.eq.3) then 

C Glare image selected.
        SCENEPURP(ISCENE)='Glare'
        intext = 1  ! internal view
        iadobe = 2  ! medium adobe texture
        SCENEMOTL(ISCENE)='medium'
        imgqua = 'Medium'
        detlvl = 'Medium'
        llvar  = 'High'
        indrcb = 2
        ipicx  = 500
        penumb = 'True'
        crenrp = '-av 0 0 0 '
        coconv = '  '
        if(paired)then
          SCENEPURP(iother)='Glare'
          SCENEMOTL(iother)='medium'
        endif
      elseif (IRP.eq.4) then 

C Daylight factor calculation selected.
        SCENEPURP(ISCENE)='Day_fact'
        intext = 1  ! internal view
        iadobe = 2  ! medium adobe texture
        SCENEMOTL(ISCENE)='medium'
        imgqua = 'Medium'
        detlvl = 'Medium'
        llvar  = 'High'
        indrcb = 2
        ipicx  = 500
        penumb = 'True'
        crenrp = '-av 0 0 0 -ab 1'
        coconv = '  '
        if(paired)then
          SCENEPURP(iother)='Day_fact'
          SCENEMOTL(iother)='medium'
        endif
      elseif (IRP.eq.5) then 

C Coupling calculation selected. (thermal/lighting simulation)
        SCENEPURP(ISCENE)='Coupling'
        intext = 1  ! internal view
        iadobe = 2  ! medium adobe texture
        SCENEMOTL(ISCENE)='medium'
        imgqua = 'Medium'
        detlvl = 'Medium'
        llvar  = 'High'
        indrcb = 2
        ipicx  = 500
        penumb = 'False'
        crenrp = '-av 0 0 0  '
        coconv = '  '
        vewcmds(ISCENE)='  '
        indxvew = 0
        if(paired)then
          SCENEPURP(iother)='Coupling'
          SCENEMOTL(iother)='medium'
        endif
      elseif (IRP.eq.6) then 

C Daylight coefficient calculation selected. (thermal/lighting simulation)
        SCENEPURP(ISCENE)='Day_coef'
        intext = 1  ! internal view
        iadobe = 2  ! medium adobe texture
        SCENEMOTL(ISCENE)='medium'
        imgqua = 'Medium'
        detlvl = 'Medium'
        llvar  = 'High'
        indrcb = 2
        ipicx  = 500
        penumb = 'False'
        crenrp = '-av 0 0 0 '
        coconv = '  '
        if(paired)then
          SCENEPURP(iother)='Day_coef'
          SCENEMOTL(iother)='medium'
        endif
      elseif (IRP.eq.7) then 

C Daylight lux (illuminance at grid points with CIE sky) calculation selected.
        SCENEPURP(ISCENE)='Day_lux'
        intext = 1  ! internal view
        iadobe = 2  ! medium adobe texture
        SCENEMOTL(ISCENE)='medium'
        imgqua = 'Medium'
        detlvl = 'Medium'
        llvar  = 'High'
        indrcb = 2
        ipicx  = 500
        penumb = 'True'
        crenrp = '-av 0 0 0 -ab 1'
        coconv = '  '
        if(paired)then
          SCENEPURP(iother)='Day_lux'
          SCENEMOTL(iother)='medium'
        endif
      elseif (IRP.eq.8) then 

C Illuminance at grid points with std sky calculation selected.
        SCENEPURP(ISCENE)='Illum'
        intext = 1  ! internal view
        iadobe = 2  ! medium adobe texture
        SCENEMOTL(ISCENE)='medium'
        imgqua = 'Medium'
        detlvl = 'Medium'
        llvar  = 'High'
        indrcb = 2
        ipicx  = 500
        penumb = 'True'
        crenrp = '-av 0 0 0 -ab 1'
        coconv = '  '
        if(paired)then
          SCENEPURP(iother)='Illum'
          SCENEMOTL(iother)='medium'
        endif
      elseif (IRP.eq.9) then 

C Sun positions.
        SCENEPURP(ISCENE)='Night_ext'
        intext = 0  ! external view
        iadobe = 2  ! medium adobe texture
        SCENEMOTL(ISCENE)='medium'
        imgqua = 'Medium'
        detlvl = 'Medium'
        llvar  = 'High'
        indrcb = 2
        ipicx  = 500
        penumb = 'True'
        crenrp = '-av 0.0005 0.0005 0.0005 '
        coconv = '  '
        ifocz = -1
        write(SCENEFZNM(ISCENE),'(a)') 'all'
        if(paired)then
          SCENEPURP(iother)='Night_ext'
          SCENEMOTL(iother)='medium'
        endif
      endif

C Set default names based on above.
      if (IRP.eq.2.or.IRP.eq.3.or.IRP.eq.5.or.IRP.eq.6.or.
     &    IRP.eq.7.or.IRP.eq.8) then
        write(outse,'(3a)') zname(ifocz)(1:lnzname(ifocz)),
     &     '_',SCENEPURP(ISCENE)(1:6)
        write(outsw,'(4a)') zname(ifocz)(1:lnzname(ifocz)),
     &     '_',SCENEPURP(ISCENE)(1:6),'_W'
        write (SCENEDESC(ISCENE),'(5a)') 'Scene of ',
     &     zname(ifocz)(1:lnzname(ifocz)),' for ',
     &     SCENEPURP(ISCENE)(1:lnblnk(SCENEPURP(ISCENE))),' analysis.'
        if(paired)then
          write (SCENEDESC(iother),'(5a)') 'Scene of ',
     &      zname(ifocz)(1:lnzname(ifocz)),' for ',
     &      SCENEPURP(iother)(1:lnblnk(SCENEPURP(iother))),' analysis.'
        endif
      else
        IL=min(lnblnk(cfgroot),21)
        write (outse,'(3a)') cfgroot(1:IL),'_',SCENEPURP(ISCENE)(1:6)
        write (SCENEDESC(ISCENE),'(5a)') 'Scene of ',
     &     cfgroot(1:lnblnk(cfgroot)),' for ',
     &     SCENEPURP(ISCENE)(1:lnblnk(SCENEPURP(ISCENE))),' analysis.'
        if(paired)then
          write (SCENEDESC(iother),'(5a)') 'Scene of ',
     &      cfgroot(1:lnblnk(cfgroot)),' for ',
     &      SCENEPURP(iother)(1:lnblnk(SCENEPURP(iother))),' analysis.'
        endif
      endif
      write (SCENERT(ISCENE),'(a)') outse(1:28)
      if(paired)then
        write (SCENERT(iother),'(a)') outsw(1:28)
      endif
      write (RIFNAME(ISCENE),'(a,a)') 
     &       SCENERT(ISCENE)(1:lnblnk(SCENERT(ISCENE))),'.rif'
      write (LBSRIF(ISCENE),'(a,a)') 
     &       SCENERT(ISCENE)(1:lnblnk(SCENERT(ISCENE))),'_a.rif'
      if(paired)then
        if (SCENEPURP(iother)(1:8).eq.'External')then
          rifwave='obj_ex.rif'
        elseif(SCENEPURP(iother)(1:8).eq.'Internal')then
          rifwave='obj_in.rif'
        elseif(SCENEPURP(iother)(1:5).eq.'Glare')then
          rifwave='obj_gl.rif'
        elseif (SCENEPURP(iother)(1:8).eq.'Day_fact') then
          rifwave='obj_df.rif'
        else
          rifwave='obj.rif'
        endif
        write (RIFNAME(iother),'(a)') rifwave(1:lnblnk(rifwave))
      endif

C Additional details for specific purposes.
      if (IRP.eq.4.or.IRP.eq.7.or.IRP.eq.8) then

C Set up daylight factor or illuminance calc points.
        call setdfgrid(izone,ier)

      elseif (IRP.eq.5.or.IRP.eq.6) then 

C Get sensor information from cgc file and create a slc file.
        call SENLOC
      endif

C Clear command line variables. << why do we do this? >>
      zone='NONE'
C      aim='NONE'
      cmdact='NONE'

C If not in silent mode then confirm names with user.
      if (MODE.eq.'-') then

C Ask for root name and use this for rif name.
 15     write(tname,'(a)') SCENERT(ISCENE)
        write(dname,'(a)') SCENERT(ISCENE)
        CALL EASKS(tname,'Scene root name?',' ',28,dname,
     &    'root name',IER,nbhelp)
        IF(tname(1:2).EQ.'  ') GOTO 15
        write(SCENERT(ISCENE),'(a)') tname(1:lnblnk(tname))
        if(paired)then
 16       write(tname,'(2a)') 
     &      SCENERT(ISCENE)(1:lnblnk(SCENERT(ISCENE))),'_W'
          write(dname,'(a)') tname(1:lnblnk(tname))
          CALL EASKS(tname,'Scene root name (Wavefront)?',' ',28,dname,
     &      'root name',IER,nbhelp)
          IF(tname(1:2).EQ.'  ') GOTO 16
          write(SCENERT(iother),'(a)') tname(1:lnblnk(tname))
        endif
        write (RIFNAME(ISCENE),'(a,a)') 
     &    SCENERT(ISCENE)(1:lnblnk(SCENERT(ISCENE))),'.rif'
        write (LBSRIF(ISCENE),'(a,a)') 
     &    SCENERT(ISCENE)(1:lnblnk(SCENERT(ISCENE))),'_a.rif'
        if(paired)then
          write(SCENEGSRC(iother),'(a)') 'Wavefront'
        endif

C Ask for additional description.
 25     write(tname72,'(a)') SCENEDESC(ISCENE)
        write(dname72,'(a)') SCENEDESC(ISCENE)
        CALL EASKS(tname72,'Description of the scene?',' ',72,
     &    dname72,'scene desc',IER,nbhelp)
        IF(tname72(1:2).EQ.'  ') GOTO 25
        write(SCENEDESC(ISCENE),'(a)') tname72(1:lnblnk(tname72))
        if(paired)then
 26       write(tname72,'(2a)') 
     &      SCENEDESC(ISCENE)(1:lnblnk(SCENEDESC(ISCENE))),'_W'
          write(dname72,'(2a)') 
     &      SCENEDESC(ISCENE)(1:lnblnk(SCENEDESC(ISCENE))),'_W'
          CALL EASKS(tname72,'Description of the scene (Wavefront)?',
     &      ' ',72,dname72,'scene desc',IER,nbhelp)
          IF(tname72(1:2).EQ.'  ') GOTO 26
          write(SCENEDESC(iother),'(a)') tname72(1:lnblnk(tname72))
        endif

C Debug
C        write(6,*) 'scene root ',SCENERT(ISCENE)
C        write(6,*) 'scene desc ',SCENEDESC(ISCENE)
C        if(paired)then
C          write(6,*) 'W scene root ',SCENERT(iother)
C          write(6,*) 'W scene desc ',SCENEDESC(iother)
C        endif
      endif

C Call rifset to set up all remaining data.
      call rifset(ifocz,'f')

      RETURN
      END

C ********************** ADDSCENE **********************
C Manages the creation of a default translation from ESP-r to
C Radiance or loaded via scanning an existing rif file.
      subroutine ADDSCENE(IRCFG)
#include "building.h"
#include "model.h"
#include "e2r_common.h"
#include "help.h"

      integer lnblnk  ! function definition
      COMMON/OUTIN/IUOUT,IUIN,IEOUT

      logical both_esp_wave      ! If true manage pair of scenes.
      logical focus_espg         ! If true editing in context of ESP-r geometry
      integer ipairedscene       ! Index of paired scene.
      COMMON/ESPWAVE/both_esp_wave,focus_espg,ipairedscene(16)

      CHARACTER pfile*72
      character ltmp*72,ltmpsw*72
      logical WFOK

      helpinsub='radcfg'  ! set for subroutine

C Ask if a default image should be created.
      helptopic='create_scene_interactive'
      call gethelptext(helpinsub,helptopic,nbhelp)

      CALL EASKMBOX('This is a new visualisation.',
     &  'Create scene:','via ESP-r geometry',
     &  'via both ESP-r & Wavefront geometry','via loading a rif file',
     &  'cancel',' ',' ',' ',' ',IWDV,nbhelp)
      if(IWDV.eq.1)then

C Create new scene file asking the user what the purpose is.
C Signal geometric source will only be ESP-r and no paired scene is to
C be created. Process with minimal user interactions.
        write (pfile,'(a,a)') runpath(1:lnrp),LRADCF(1:lnblnk(LRADCF))
        call ERPFREE(IRCFG,ISTAT)
        call FPOPEN(IRCFG,ISTAT,1,3,pfile)
        both_esp_wave=.false.; focus_espg=.true.
        call newscene('-')     ! As 'a' but dialog without default.
        cmdact='Create'
        call autorad
        ICHK=4
        SCENEGSRC(ISCENE)='ESP-r'
        call radcfgout(ICHK)
      elseif(IWDV.eq.2)then

C Generate both ESP-r and Wavefront based scenes.
        write (pfile,'(a,a)') runpath(1:lnrp),LRADCF(1:lnblnk(LRADCF))
        call ERPFREE(IRCFG,ISTAT)
        call FPOPEN(IRCFG,ISTAT,1,3,pfile)
        both_esp_wave=.true.; focus_espg=.true.
        call newscene('-')     ! As 'a' but dialog without default.
        cmdact='Create'
        call autorad
        ICHK=4
        SCENEGSRC(ISCENE)='ESP-r'
        SCENEGSRC(ISCENE+1)='Wavefront'

C Prior to calling mkriofil if a Wavefront variant is being processed
C then the name of the Wavefront object file would be good to know. 
        WFOK=.false.
        iother=ipairedscene(ISCENE)
        if(SCENEGSRC(ISCENE)(1:9).eq.'Wavefront')then
          if(SCENEWAVE(ISCENE)(1:7).eq.'UNKNOWN')then
            CALL EASKOK(' ','Supply name of Wavefront object file?',
     &        WFOK,nbhelp)
          endif
        endif
        if(iother.gt.0)then
          if(SCENEGSRC(iother)(1:9).eq.'Wavefront')then
            if(SCENEWAVE(iother)(1:7).eq.'UNKNOWN')then
              CALL EASKOK(' ','Supply name of Wavefront object file?',
     &          WFOK,nbhelp)
            endif
          endif
        endif
        if(WFOK)then
          call edisp(iuout,' ')
          call edisp(iuout,
     &      'Place the Wavefront obj and mtl files in the model rad')
          call edisp(iuout,
     &      'folder. The conversion to a Radiance mesh file will be')
          call edisp(iuout,'done in that folder.')
 642      helptopic='e2r_RIF_file_name'
          call gethelptext(helpinsub,helptopic,nbhelp)
C          lnrp=lnblnk(runpath)  ! remember this 
C          write(ltmpsw,'(2a)') runpath(1:lnrp),'xx.obj'
          write(ltmpsw,'(a)') 'xx.obj'
          CALL EASKS(ltmpsw,
     &      'Wavefront object file (assumed to be in rad folder)?',
     &      '  ',72,'xxx.obj','wavefront file',IER,nbhelp)
          if(ltmpsw.eq.' ')goto 642
          write(SCENEWAVE(ISCENE+1),'(a)') ltmpsw(1:lnblnk(ltmpsw))
        endif
        ICHK=4
        call radcfgout(ICHK)
      elseif(IWDV.eq.3)then

C Create scene from an existing rif file (i.e. a wavefront related rif).
        if(NSCENE+1.gt.10)then
          call usrmsg('Cannot add another scene (limit is 10).',
     &      ' ','W')
          return
        endif
 643    helptopic='e2r_RIF_file_name'
        call gethelptext(helpinsub,helptopic,nbhelp)
        write(ltmp,'(a)') runpath(1:lnrp)
        CALL EASKS(ltmp,'Radiance RIF file?',
     &    '  ',72,'xxx.rif','rif file name',IER,nbhelp)
        if(ltmp.eq.' ')goto 643

C Set new scene to focus scene.
        NSCENE = NSCENE+1
        ISCENE=NSCENE
        RIFNAME(ISCENE) = ltmp
        call RRIF(ita2,ita3,RIFNAME(ISCENE),'s',ier)
        CALL EASKMBOX('Source of geometry:',' ',
     &     'Generated by ESP-r','Based on wavefront file',' ',
     &     ' ',' ',' ',' ',' ',IW,nbhelp)
        if(IW.eq.1)then
          write(SCENEGSRC(ISCENE),'(a)') 'ESP-r'
        elseif(IW.eq.2)then
          write(SCENEGSRC(ISCENE),'(a)') 'Wavefront'
        endif
        ICHK=4
        call radcfgout(ICHK)
      endif

      RETURN
      END

C ********************** setdfgrid **********************
C Manages creation of the daylight factor grid points or
C illuminance grid points. Passed in izone (indicator of
C whether the focus zone has already been set zero = unset).
      subroutine setdfgrid(izone,ier)
#include "building.h"
#include "model.h"
#include "geometry.h"
#include "prj3dv.h"
#include "e2r_common.h"
#include "help.h"

      integer lnblnk  ! function definition

      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      COMMON/SPAD/MMOD,LIMIT,LIMTTY
      COMMON/FILEP/IFIL
      
      integer ncomp,ncon
      common/C1/NCOMP,NCON

      COMMON/gzonpik/izgfoc,nzg,nznog(mcom)
      logical keep
      COMMON/gridfoc/keep

      character DESCRC*25
      character outs*124,tfile*72
      character ITEM(MV)*15, t12*12
      character ldfg*72   ! local file name

      dimension IVALS(MCOM), PEQN(4), ITEMV(MV)

      helpinsub='radcfg'  ! set for subroutine

C New grid so set number of grid points to zero.
      NDFP=0

C Remember current graphics display options.
      call SVDSOPT

C Ask for focus zone (if not already defined), surface and working plane height.
      if (IZONE.eq.0) then
        helptopic='select_focus_zone'
        call gethelptext(helpinsub,helptopic,nbhelp)
 32     INPIC=1
        CALL EPICKS(INPIC,IVALS,
     &    'Generate DF or illuminance for which zone. ',
     &    '(You will then be asked to grid a surface.) ',
     &    12,NCOMP,zname,'Select zone',IER,nbhelp)
        if (INPIC.eq.0) goto 32
        ifocz=IVALS(1)
        write(SCENEFZNM(ISCENE),'(a)') zname(ifocz)
        nzg=1          ! one zone while defining grid
        nznog(1)=ifocz
        izgfoc=0
      else
        ifocz=izone
        write(SCENEFZNM(ISCENE),'(a)') zname(ifocz)
        nzg=1          ! one zone while defining grid
        nznog(1)=izone
        izgfoc=0
      endif
      MODIFYVIEW=.TRUE.
      MODLEN=.TRUE.
      MODBND=.TRUE.
      ITSNM=0
      CALL redraw(IER)

C Ask for surface and highlight.
      call EASKSUR(ifocz,IS,'-',
     &  'Select a surface to grid. (You can then transform ',
     &  'the grid perpendicular to this surface.)',IER)
      CALL INLNST(1)
      CALL SURADJ(ifocz,IS,IE,TMP,IZC,ISC,IC,DESCRC)
      LINSTY(IC)=2
      MODIFYVIEW=.TRUE.
      ITVNO=0
      CALL redraw(IER)
      CALL INLNST(1)
      write(DFSURFNAME(ISCENE),'(a)') sname(ifocz,is)

C Setup help for the next few dialogs.
      helptopic='distance_to_df_grid'
      call gethelptext(helpinsub,helptopic,nbhelp)
      if(WPDIST(ISCENE).lt.0.01) WPDIST(ISCENE)=0.9
      CALL EASKR(WPDIST(ISCENE),'Distance of grid from surface?',' ',
     &   0.001,'F',2.5,'W',0.9,'grid transform',IER,nbhelp)
      call EASKMBOX('Transform grid to: ',' ','Inside','Outside',
     &  ' ',' ',' ',' ',' ',' ',IW,nbhelp)
      if (IW.eq.1) then
        dirf=-1.
      else
        dirf=1.
      endif

C Now read zone geometry file to enable edge selection.
      ITA1 = IFIL+6
      call georead(ITA1,LGEOM(ifocz),ifocz,0,IUOUT,IER)

C Generate a list of edges. JVN(MS,MV)
      if(DFEDGE(ISCENE,1).gt.0.and.DFEDGE(ISCENE,2).gt.0)then
         write(outs,'(a,i4,a,i4)') 'Current edge is ',
     &     DFEDGE(ISCENE,1),' -> ',DFEDGE(ISCENE,2)
         call edisp(iuout,'  ')
         call edisp(iuout,outs)
      endif
      NI=0
      do 10 I=1,NVER(IS)-1
        write(ITEM(I),'(2(a,i3))')'edge v',JVN(IS,I),'->v',JVN(IS,I+1)
        NI=NI+1
 10   continue
      NI=NI+1
      write(ITEM(NI),'(2(a,i3))')'edge v',JVN(IS,NI),'->v',JVN(IS,1)
 12   INPIC=1
      CALL EPICKS(INPIC,ITEMV,
     &  'Select an edge (defined by vertex numbers). ',
     &  'The grid will be parallel to this edge. ',
     &  15,NI,ITEM,'Select edge',IER,nbhelp)
      if (INPIC.eq.0) goto 12

C Translate edge selection into vertex numbers.
      if (ITEMV(1).eq.NI) then
        IVS=JVN(IS,NI)
        IVE=JVN(IS,1)
      else
        IVS=JVN(IS,ITEMV(1))
        IVE=JVN(IS,ITEMV(1)+1)
      endif
      DFEDGE(ISCENE,1)=IVS
      DFEDGE(ISCENE,2)=IVE

C Grid density.
      helptopic='df_grid_density'
      call gethelptext(helpinsub,helptopic,nbhelp)

C Major grid lines are parallel with the selected edge.
      if(igmajor(ISCENE).eq.0) igmajor(ISCENE)=4
      if(igminor(ISCENE).eq.0) igminor(ISCENE)=6
      write (t12,'(2i4)') igmajor(ISCENE),igminor(ISCENE)
 25   CALL EASKS(t12,'Grid density? (First item is number of grid',
     &   'points in the rows parallel to the chosen edge.) ',
     &   12,' 2  2','grid density',IER,nbhelp)
      K=0
      CALL EGETWI(t12,K,igma,1,25,'W','grid lines',IER)
      CALL EGETWI(t12,K,igmi,1,25,'W','grid lines',IER)
      igmajor(ISCENE)=igma
      igminor(ISCENE)=igmi
      if (igma.lt.1.or.igmi.lt.1) then
        call usrmsg('There must be at least 1 grid',
     &              'line in each direction.','W')
        goto 25
      endif

C Calculate grid points.
      WPD=WPDIST(ISCENE)*dirf
      call GRID(IS,igma,igmi,IVS,IVE,WPD,DFPTS,NDFP,PEQN,IER)
      call ERPFREE(ITA1,ISTAT)

C Get file name for grid points.
      write (ldfg,'(a,a4)') 
     &  SCENERT(ISCENE)(1:lnblnk(SCENERT(ISCENE))),'.grd'
 55   CALL EASKS(LDFG,'DF or illuminance grid points file name?',
     &  ' ',72,'points.grd','df points file name',IER,nbhelp)
      IF(LDFG.EQ.' ') GOTO 55
      LDFGRID(ISCENE)=LDFG

C Open df points file in working directory.
      write (tfile,'(a,a)') runpath(1:lnrp),
     &  LDFGRID(ISCENE)(1:lnblnk(LDFGRID(ISCENE)))
      call FPOPEN(ITA1,ISTAT,1,3,tfile)
      do 552 IC=1,3
        PEQN(IC)=PEQN(IC)*dirf
 552  continue

C Write points and sensor direction (plane normal).
      do 555 IV=1,NDFP
       write(ITA1,'(6(f8.4,1x))')(DFPTS(IV,IX),IX=1,3),(PEQN(IY),IY=1,3)
 555  continue
      DFDIR(1)=PEQN(1)
      DFDIR(2)=PEQN(2)
      DFDIR(3)=PEQN(3)
      call ERPFREE(ITA1,ISTAT)

C Display grid on the zone... (if in graphic mode)
      if (MMOD.EQ.8) then
        call edisp(iuout,'  ')
        call edisp(iuout,'The proposed calculation points are shown')
        call edisp(iuout,'as dots in the graphic feedback.')
        MODIFYVIEW=.TRUE.
        ITZNM=1
        ITSNM=1
        ITVNO=1
        CALL redraw(IER)
        CALL DFGRIDLBL('g')
      endif

C Ask user if focus should remain on this zone or whole model
      helptopic='df_stay_focused'
      call gethelptext(helpinsub,helptopic,nbhelp)
      CALL EASKOK(' ','Keep graphic feedback focused on this zone?',
     &  KEEP,nbhelp)

C Note that DF calculations often are dependent on adjacent spaces so
C reset nzg to include all zones.
      nzg=NCOMP 
      izgfoc=0
      do 300 I=1,NCOMP
        nznog(I)=I
 300  continue

C Reset display options.
      call RCDSOPT

      return
      end

C ********************** DFGRIDLBL **********************
C Writes an index of daylight grid points (act = g|G) or
C daylight factors (act = d|D) or lux (act = l|L).
      SUBROUTINE DFGRIDLBL(act)
#include "building.h"
#include "prj3dv.h"
#include "e2r_common.h"
      
      integer lnblnk  ! function definition

C Parameters
      character act*1  ! action

C      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      COMMON/FILEP/IFIL
      COMMON/SPAD/MMOD,LIMIT,LIMTTY
      integer ifs,itfs,imfs
      COMMON/GFONT/IFS,ITFS,IMFS
      COMMON/gzonpik/izgfoc,nzg,nznog(mcom)
      logical keep
      COMMON/gridfoc/keep

      DIMENSION  COG1(3),COG2(3)

      CHARACTER VT*6
      character OUTSTR*124,tfile*72

#ifdef OSI
      integer iix,iiy,iid1,iid2,iicol
#else
      integer*8 iix,iiy,iid1,iid2,iicol
#endif
      
C If not in graphical mode return.
      if (MMOD.ne.8) return

C If there are no grid points warn the user and check if
C there is a *.df file or a *.df- file that can be scanned again.
      if(NDFP.eq.0)then
        ITA1 = IFIL+6

        call ERPFREE(ITA1,ISTAT)
        write(tfile,'(a,a)')runpath(1:lnrp),
     &    LDFGRID(ISCENE)(1:lnblnk(LDFGRID(ISCENE)))
        call FPOPEN(ITA1,ISTAT,1,1,tfile)
        if(ISTAT.eq.0)then
          continue
        else
          call ERPFREE(ITA1,ISTAT)
          write(tfile,'(3a)')runpath(1:lnrp),
     &      LDFGRID(ISCENE)(1:lnblnk(LDFGRID(ISCENE))),'-'
          call FPOPEN(ITA1,ISTAT,1,1,tfile)
          if(ISTAT.eq.0)then
            continue
          else
            CALL USRMSG('Could not read grid points file.',
     &       'Cannot report without this information. Returning.','W')
            return
          endif
        endif
  42    CALL STRIPC(ITA1,OUTSTR,0,ND,0,'grid points',IFER)
        if (IFER.eq.0) then
          NDFP=NDFP+1
          K=0
          CALL EGETWR(OUTSTR,K,XTMP,0.,0.,'-','grid x',IER)
          DFPTS(NDFP,1)=XTMP
          CALL EGETWR(OUTSTR,K,XTMP,0.,0.,'-','grid y',IER)
          DFPTS(NDFP,2)=XTMP
          CALL EGETWR(OUTSTR,K,XTMP,0.,0.,'-','grid z',IER)
          DFPTS(NDFP,3)=XTMP
          CALL EGETWR(OUTSTR,K,DFDIR(1),0.,0.,'-','vec x',IER)
          CALL EGETWR(OUTSTR,K,DFDIR(2),0.,0.,'-','vec y',IER)
          CALL EGETWR(OUTSTR,K,DFDIR(3),0.,0.,'-','vec z',IER)
          goto 42
        endif
        call ERPFREE(ITA1,ISTAT)
        if (NDFP.eq.0) then
          CALL USRMSG('Could not read grid points file.',
     &      'Cannot report without this information. Returning.','W')
          return
        endif
      endif

C If vertices are to be labeled do it in the smallest font.
      call SVDSOPT
      IF(IFS.GT.0)CALL winfnt(0)
      MODIFYVIEW=.TRUE.
      MODLEN=.TRUE.
      MODBND=.TRUE.
      ITZNM=1
      ITSNM=1
      ITVNO=1

C If viewing DF and user asked to keep focus on the zone
C reset before calling redraw.
      if(keep)then
        nzg=1
        nznog(1)=ifocz
        izgfoc=ifocz
      endif
      CALL redraw(IER)

C Loop through grid points.
      do 10 I=1,NDFP
        COG1(1)=DFPTS(I,1)
        COG1(2)=DFPTS(I,2)
        COG1(3)=DFPTS(I,3)
        CALL VECTRN(COG1,TSMAT,COG2,IER)
        call u2pixel(COG2(1),COG2(2),iix,iiy)
        CALL ecirc(iix,iiy,3,1)
        if(act.eq.'g'.or.act.eq.'G')then
          IF(I.le.9)THEN
            WRITE(VT,'(i1)') I
          ELSEIF(I.gt.9.and.I.le.99)THEN
            WRITE(VT,'(i2)') I
          ELSEIF(I.gt.99)THEN
            WRITE(VT,'(i3)') I
          ENDIF
        elseif(act.eq.'d'.or.act.eq.'D')then
          IF(DFVALS(I).lt.10.0)THEN
            WRITE(VT,'(f4.2)') DFVALS(I)
          ELSEIF(DFVALS(I).ge.10.0)THEN
            WRITE(VT,'(f4.1)') DFVALS(I)
          ENDIF
        elseif(act.eq.'l'.or.act.eq.'L')then
          IF(DFVALS(I).lt.1000.0)THEN
            WRITE(VT,'(f5.0)') DFVALS(I)
          ELSEIF(DFVALS(I).ge.1000.0)THEN
            WRITE(VT,'(f6.0)') DFVALS(I)
          ENDIF
        endif
        iid1=iix+1; iid2=iiy-1; iicol=0
        call textatxy(iid1,iid2,VT,'-',iicol)
 10   continue

C Reset display options.
      call RCDSOPT
      IF(IFS.GT.0)CALL winfnt(IFS)

      RETURN
      END


C ********************** SENLOC **********************
C Creates a file of sensor locations and orientations
C for use with a 'Coupling' scene purpose.

      SUBROUTINE SENLOC
#include "building.h"
#include "model.h"
#include "e2r_common.h"
      
      integer lnblnk  ! function definition

      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      COMMON/FILEP/IFIL

      COMMON/CGCIN2N/NCGTC(MCOM,MDTY),NLITZ(MCOM),IDFST(MCOM,MLCOM),
     & CGX(MCOM,MLCOM,MDF),CGY(MCOM,MLCOM,MDF),CGH(MCOM,MLCOM,MDF),
     & UX(MCOM,MLCOM,MDF),UY(MCOM,MLCOM,MDF),UH(MCOM,MLCOM,MDF),
     & SETPT(MCOM,MLCOM),SYSPER(MCOM,MLCOM),SOFFLL(MCOM,MLCOM),
     & IOFFDT(MCOM,MLCOM),SMLOUT(MCOM,MLCOM),SMEOUT(MCOM,MLCOM)

      COMMON/DAYF/NDF(MCOM,MLCOM),DFDAT(MCOM,MLCOM,MGT,MDF),
     & AZIJL(MCOM,MGT),SUMDF(MCOM,MLCOM,MDF),NUDWIN(MCOM)

      CHARACTER LCGCINT*72
      character outs*124

      logical XST

      write(outs,'(a)')'Reading casual gain control file.'
      call edisp (iuout,outs)
      write(outs,'(a,a)')'Searching for sensor locations ',
     &                                        'and orientations.'
      call edisp (iuout,outs)

C Define temporary unit numbers.
      ITA1 = IFIL+6
      call ERPFREE(ITA1,ISTAT)

C Set local file name.
      LCGCINT=LCGCIN(ifocz)
      ICOMP=ifocz
      ier=0
      if (LCGCINT(1:2).eq.'  ') then
        write(outs,'(a)')'Error: casual gain control file undefined!'
        call edisp (iuout,outs)
        return
      endif
      itrc=0
      call ercgcf(itrc,itru,LCGCINT,ICOMP,ier)

C Open sensor info file (do not add path as radiance is run locally).
      write(LDFGRID(ISCENE),'(a,a)')
     &  SCENERT(ISCENE)(1:lnblnk(SCENERT(ISCENE))),'.slc'

C Check for file (add working dir path)
      write (outs,'(a,a)') runpath(1:lnblnk(runpath)),
     &  LDFGRID(ISCENE)(1:lnblnk(LDFGRID(ISCENE)))
      call FINDFIL(outs,XST)
      if(XST)then
        call FPOPEN(ITA1,ISTAT,1,0,outs)
        call EFDELET(ITA1,ISTAT)
      endif
      call FPOPEN(ITA1,ISTAT,1,2,outs)

C Write sensor information to radiance grid file.
      do 10 I=1,NLITZ(ICOMP)
        if ((IDFST(ICOMP,I).EQ.4).or.(IDFST(ICOMP,I).EQ.5)) then
          do 20 J=1,NDF(ICOMP,I)
            write (ITA1,'(6f10.3)') CGX(ICOMP,I,J),CGY(ICOMP,I,J),
     &         CGH(ICOMP,I,J),UX(ICOMP,I,J),UY(ICOMP,I,J),UH(ICOMP,I,J)
 20       continue
        endif
 10   continue
      call ERPFREE(ITA1,ISTAT)

      RETURN
      END

C ********************* checkrif *********************
C Checks for existing scenes or creates a new scene  
C the given parameters if none found. If filexst
C is 'o' then old file, 'n' then create a new file.
      SUBROUTINE checkrif(filexst)
#include "building.h"
#include "model.h"
#include "e2r_common.h"
      
      integer lnblnk  ! function definition

      COMMON/FILEP/IFIL

      logical both_esp_wave      ! If true manage pair of scenes.
      logical focus_espg         ! If true editing in context of ESP-r geometry
      integer ipairedscene       ! Index of paired scene.
      COMMON/ESPWAVE/both_esp_wave,focus_espg,ipairedscene(16)

      character filexst*1
      character pfile*72
      character rifwave*72

      LOGICAL modscene

      IER=0
      ISCENE=0
      modscene=.false.

C Clear rcf file unit number/ channel.
      CALL ERPFREE(IRCFG,ISTAT)

C Recover zone index from command line variable zone (used for
C coupling and daylight coefficients).
      if (zone(1:4).eq.'Def_') then
        izone=-1
      else
        if (aim(1:8).eq.'External') then
          izone=-1
          zone='  '
        else
          call zindex(zone,izone)
        endif
      endif

C If rcf file does exist then search for existing rif file.
C Loop over existing scenes to see if one matches the command
C line argument passed to e2r.
      if (filexst.eq.'o') then
        do 10 i=1,NSCENE
          if (SCENEPURP(I)(1:5).eq.aim(1:5)) then

C Possible match, further check - look at rif file.
            ifocz=-1
            call RRIF(ITA2,ITA3,RIFNAME(I),'s',IER)
            if (ifocz.eq.IZONE) then

C Got match.
              ISCENE=I

C Update sensor information from cgc file.
              if(aim(1:8).eq.'Coupling'.or.aim(1:8).eq.'Day_coef')then
                call SENLOC
              endif

C If daylight coefficients and cmdact=create and model exists then issue warning.
              if(aim(1:8).eq.'Day_coef'.and.
     &           cmdact(1:6).eq.'Create')then
                call usrmsg('Model already exists!',
     &                      'Using existing model.','W')
                cmdact='Calculate'
              endif
              modscene=.true.
            endif
          endif
 10     continue
      else

C Create a new rcf file.
        write (LRADCF,'(a,a)') cfgroot(1:lnblnk(cfgroot)),'.rcf'

C Set file name including path and open.
        write(pfile,'(a,a)') runpath(1:lnblnk(runpath)),
     &    LRADCF(1:lnblnk(LRADCF))
        call ERPFREE(IRCFG,ISTAT)
        call FPOPEN(IRCFG,ISTAT,1,2,pfile)
      endif

C If not an existing scene then use newscene to configure
C it and then autorad to build the matching descriptive
C files and save to rif file and lastly the esp-r rcf file.
      if (.not.modscene) then
        call newscene('s')
        cmdact='Create'
        call autorad
        IFC=3
        SCENEGSRC(ISCENE)='ESP-r'
        if(both_esp_wave)then

C Ensure both scenes have same purpose.
          write(SCENEGSRC(ISCENE+1),'(a)') 'Wavefront'
          SCENEPURP(ISCENE+1)=SCENEPURP(ISCENE)
          call mkrif(ITA2,RIFNAME(ISCENE),IFC,0,'-',IER)
          if (SCENEPURP(ISCENE)(1:8).eq.'External')then
            rifwave='obj_ex.rif'
          elseif(SCENEPURP(ISCENE)(1:8).eq.'Internal')then
            rifwave='obj_in.rif'
          elseif(SCENEPURP(ISCENE)(1:5).eq.'Glare')then
            rifwave='obj_gl.rif'
          elseif (SCENEPURP(ISCENE)(1:8).eq.'Day_fact') then
            rifwave='obj_df.rif'
          else
            rifwave='obj.rif'
          endif
          call mkrif(ITA2,rifwave,IFC,0,'w',IER)
        endif
        if (NBSRIF(ISCENE).eq.1) then
          call mkrif(ita2,LBSRIF(ISCENE),IFC,1,'-',ier)
        endif

C Update rcf file.
        ICHK=3
        call RADCFGOUT(ICHK)
        cmdact='NONE'   ! reset because we have done it
      endif

      return
      end


C ********************** RDCFTMP **********************
C Read daylight coefficients to common and updates casual 
C gain control file.

      SUBROUTINE RDCFTMP(LFTMP,ISET)
#include "building.h"
#include "model.h"
#include "e2r_common.h"

      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      COMMON/FILEP/IFIL

      COMMON/CGCIN2N/NCGTC(MCOM,MDTY),NLITZ(MCOM),IDFST(MCOM,MLCOM),
     & CGX(MCOM,MLCOM,MDF),CGY(MCOM,MLCOM,MDF),CGH(MCOM,MLCOM,MDF),
     & UX(MCOM,MLCOM,MDF),UY(MCOM,MLCOM,MDF),UH(MCOM,MLCOM,MDF),
     & SETPT(MCOM,MLCOM),SYSPER(MCOM,MLCOM),SOFFLL(MCOM,MLCOM),
     & IOFFDT(MCOM,MLCOM),SMLOUT(MCOM,MLCOM),SMEOUT(MCOM,MLCOM)

      COMMON/DAYF/NDF(MCOM,MLCOM),DFDAT(MCOM,MLCOM,MGT,MDF),
     & AZIJL(MCOM,MGT),SUMDF(MCOM,MLCOM,MDF),NUDWIN(MCOM)

C Daylight coefficient common block data definition:
C MCOM  - max no of thermal zones.
C MLCOM - max no of lighting zones per thermal zone.
C MDF   - max no of sensors per lighting zone.
C MDCFP - max no of daylight coefficient data sets.
C MDCFS - max no of alternative daylight coefficient data sets.
C MDCFV - max no of sky patches.
      COMMON/DCOEF/NDCFP,NDCFS(MDCFP),DCOEF(MDCFP,MDCFS,MDCFV),
     & IDCFID(MCOM,MLCOM,MDF)

      CHARACTER*72 LCGCINT,LFTMP
      character outs*124

C Define temporary unit numbers.
      ITA1 = IFIL+6
      call ERPFREE(ITA1,ISTAT)

C Set local file name.
      LCGCINT=LCGCIN(ifocz)
      ICOMP=ifocz
      ier=0
      if (LCGCINT(1:2).eq.'  ') then
        write(outs,'(a)')'Error: casual gain control file undefined!'
        call edisp (iuout,outs)
        return
      endif
      itrc=0
      call ercgcf(itrc,itru,LCGCINT,ICOMP,ier)

C Read illuminance information from temp file.
      call FPOPEN(ITA1,ISTAT,1,0,LFTMP)
      do 5 IP=1,145
       do 10 ILZ=1,NLITZ(ICOMP)
        if (IDFST(ICOMP,ILZ).EQ.5) then
         do 20 ISE=1,NDF(ICOMP,ILZ)
          CALL STRIPC(ITA1,OUTS,0,ND,0,'transfer file',IER)
          if (ier.ne.0) then 
           call usrmsg('Cannot read transfer file:',LFTMP,'F')
           return
          endif
          K=0
          call EGETWR(OUTS,K,SIL,0.,10000.,'W','illuminance',IER)
          JSEN=IDCFID(ICOMP,ILZ,ISE)
          DCOEF(JSEN,ISET,IP)=SIL
 20      continue
        endif
 10    continue
 5    continue
      call ERPFREE(ITA1,ISTAT)

C Write data to casual gain control file.
      call CASCTMK(LCGCINT,ICOMP,'s',IER)

      RETURN
      END
