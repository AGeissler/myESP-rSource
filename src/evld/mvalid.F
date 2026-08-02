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

C ******
C This file contains the following subroutines:
C  MValid   - Menu for conduction, BESTEST and Selftests.
C  selftest - Scan tester_list file & user selections for self-test assessments.
C  RUNSELF  - Read model details and directives and and invoke self-test runs.
C  CMPASCI  - Line by line compare of two ASCII files.
C  CMPSUMRY - Line by line compare of values in .summary files.
C  exemplartest - Scan exemplar_list file & user selections for tests. [work-in-progress]
C  runexemplars - read model and directives and invoke tests runs. [work-in-progress]

C *********** MValid ***************
C The main menu for conduction, BESTEST and self-testing.
      subroutine MValid
#include "help.h"

      common/OUTIN/IUOUT,IUIN,IEOUT
      common/GR1D04/GAM,RGAM
      common/Vld04/DTime,NTime,MonNod,NRoot,NIter
      common/Vld10/VlDirR,VlDirL
      common/FILEP/IFIL

      character*25 item(8)
      character VlDirR*96,VlDirL*72
      integer MITEM,INO ! max items and current menu item

      helpinsub='mvalid'  ! set for subroutine

C  VlDirR='/opt/esp-r/validation'
      DTime=3600.0
      NTime=24
      MonNod=1
      NRoot=500
      NIter=1000
      Gam=0.5

   10 INO=-3
      ITEM(1)='a conduction'
      ITEM(2)='b BESTEST'
      ITEM(3)='c self-test'
      ITEM(4)='d exemplar-test'
      ITEM(5)='? help'
      ITEM(6)='- exit menu'
      MITEM=6

      helptopic='validation-overview'
      call gethelptext(helpinsub,helptopic,nbhelp)
      CALL EMENU('Self testing',ITEM,MITEM,INO)

C Help.
      IF(INO.EQ.MITEM-1)THEN
        CALL PHELPD('vld opening',5,'-',0,0,IER)

C Analytical solution for the transient conduction within a
C multi-layered construction.
      ELSEIF(INO.EQ.1)THEN
        CALL ANLCND

      ELSEIF(INO.EQ.2)THEN  ! BESTEST tests.
        call STNDRDS
          
      ELSEIF(INO.EQ.3)THEN  ! self testing.
        call selftest
          
      ELSEIF(INO.EQ.4)THEN  ! exemplar testing.
        call exemplartest
          
      ELSEIF(INO.EQ.MITEM)THEN  ! Exit menu.
        return
 
      ENDIF
      
      GOTO 10
      END

C ************** selftest ***************
C Scan the tester_list file [selftest_list], parse information and offer
C the user a selection of or information about one or more of these.
C iexfil is the file unit number.
C exfil is the tester_list file name. 
      subroutine selftest
#include "espriou.h"
#include "help.h"

      COMMON/SPAD/MMOD,LIMIT,LIMTTY
      common/OUTIN/IUOUT,IUIN,IEOUT
      common/FILEP/IFIL
      common/selflst/selflbl,selffl
      character selflbl*20,selffl*96
      common/selflstopen/summary_opnd,sum_hdr,adiffs,tdiffs,anots,mtest
      logical summary_opnd,sum_hdr
      real adiffs,tdiffs,anots

C name is array of string holding *name from the file.
C inl2ngr is a pointer from the menu item back to the group index.
C group is array of strings selname(i)='UNKNOWN' holding *groups in the file.
C gitem is menu string array for high level group selection.
C gh() is array of help entries.

      dimension group(25),gitem(30)
      dimension inl2ngr(50)
      character gitem*44,outstr*124,word*20
      character exfil*96,group*42,tgroup*42,key*1
      character separ*42
      logical XST
      integer irpt  ! for radio button
      integer NCOG,ICOG ! max items and current menu item

#ifdef OSI
      integer impx,impy,iwe
#else
      integer*8 impx,impy,iwe
#endif

      helpinsub='mvalid'  ! set for subroutine
      helptopic='self-checks-overview'
      call gethelptext(helpinsub,helptopic,nbhelp)

C selffl is located in the source distribution in the file
C models/validation/selftest/selftest_list.
      currentfile='  '; selffl='  '
      summary_opnd = .false.; sum_hdr= .false.
      adiffs = 0; tdiffs = 0; anots = 0
      XST=.false.
      ISTRW=96
      CALL EASKF(selffl,' ','Self test list file?',
     &   ISTRW,'selftest_list','Self test list',IER,nbhelp)
      
      if(selffl(1:2).ne.'  '.and.selffl(1:4).ne.'UNKN')then
        INQUIRE (FILE=selffl,EXIST=XST)
        if (.not.XST) then
          call edisp(iuout,'No selftest list available (file missing)')
          call edisp(iuout,'returning to main menu.')
          return
        endif
      else
        call edisp(iuout,'No selftest list available (undefined)')
        call edisp(iuout,'returning to main menu.')
        return
      endif
      exfil=selffl
      write(currentfile,'(a)') selffl(1:lnblnk(selffl))
      
C Open the self test models file.
      IER=0
      iexfil=IFIL+1
      CALL EFOPSEQ(iexfil,exfil,1,IER)
      if(ier.ne.0) then
        call edisp(iuout,'Cannot open selftest description file')
        call edisp(iuout,'returning to main menu.')
        ier=1
        return
      endif
      
C Read the header and see if '*SELFTEST'.
      CALL STRIPC(iexfil,outstr,0,ND,0,'selftest header',IER)
      if(outstr(1:9).ne.'*SELFTEST')then
        call edisp(iuout,' Not a selftest models file!')
        ier=1
        return
      endif
 
C ngr is the number of groups found in the file (max 15).
C ing is index of menu string to be written to. This is separate from
C     inl so that labels do not count in the generation of KEY characters.
C inl is the total number of groups and labels to display.
      ngr=0
      ing=0
      inl=0
      inl2ngr(1)=0
      igrsel=0

 12   CALL STRIPC(iexfil,outstr,0,ND,0,'group id',IER)
      if(ier.ne.0)goto 901
      K=0
      CALL EGETW(outstr,K,word,'W','group id',IER)
      if(word(1:6).eq.'*group')then

C If we have less than 25 groups, increment ngr, inl, ing. Set
C KEY character using ing value, write menu string gitem and set
C pointer inl2ngr for menu line inl to ngr value.
        if(ngr.lt.15)then
          ngr=ngr+1
          inl=inl+1
          ing=ing+1
          CALL EMKEY(ing,KEY,IER)
          call egetrm(outstr,K,group(ngr),'W','group id',IER)
          write(gitem(inl),'(a1,1x,a42)') key, group(ngr)(1:42)
          inl2ngr(inl)=ngr
        else
          call edisp(iuout,' Too many groups!')
        endif
      elseif(word(1:6).eq.'*label')then

C A label included. Only increment inl and set pointer to zero.
        inl=inl+1
        separ=' '
        call egetrm(outstr,K,separ,'W','label',IER)
        ln=max(1,lnblnk(separ))
        ln=min0(ln,42)
        write(gitem(inl),'(2x,a)') separ(1:ln)
        inl2ngr(inl)=0
      endif
      
C Loop unti EOF error detected above (jumps to 901 and closes file).
      goto 12
 901  CALL ERPFREE(iexfil,ISTAT)

C Present the display, check the maximum width of gitem and
C use vwmenu or emenu depeding on graphic or text display.
 444  ICOG= -2
      gitem(inl+1)=  '  _____________________ '
      gitem(inl+2)=  '* all available models  '
      gitem(inl+3)=  '  _____________________ '
      gitem(inl+4)=  '? help                  '
      gitem(inl+5)=  '- exit menu'
      NCOG=inl+5

C Help @ validation list file level.
      helptopic='selftest-overview'
      call gethelptext(helpinsub,helptopic,nbhelp)
      IW=0
      do 141 ij=1,NCOG
        IW=max(lnblnk(gitem(ij)),IW)
 141  continue
 
C If in graphics mode use a variable width menu.
      if(MMOD.EQ.8)then
        impx=0
        impy=0
        iwe=iw
        CALL VWMENU('Selftest sets',gitem,NCOG,impx,impy,iwe,
     &     irpx,irpy,ICOG)
      else
        CALL EMENU('Selftest sets',gitem,NCOG,ICOG)
      endif
      
C Deal with the user selection.
      IF(ICOG.EQ.NCOG)THEN

C Exit.  Close file and return to calling menu.
        CALL ERPFREE(iexfil,ISTAT)
        return
      ELSEIF(ICOG.EQ.(NCOG-1))THEN

C Display help messages.
        CALL PHELPD('selftest list menu',nbhelp,'-',0,0,IER)
      ELSEIF(ICOG.EQ.(NCOG-3))THEN

C Automatically run all of the available models.
        irpt=0
        if(irpt.ne.3)then
          do 200 IGRP=1,NGR
            TGROUP=group(IGRP)
            IFLG=1
            if (IGRP.eq.NGR) IFLG=-1
            call runself(exfil,TGROUP,IFLG,IRPT)
 200      continue
        endif

C Select specific models from a group.
      ELSEIF(ICOG.GE.1.AND.ICOG.LT.(NCOG-4))THEN
        if(gitem(ICOG)(1:2).eq.'  ')goto 444
        if(inl2ngr(ICOG).gt.0.and.inl2ngr(ICOG).le.25)then
          igrsel=inl2ngr(ICOG)
          IFLG=0
          IRPT=0
          TGROUP=group(igrsel)
          call runself(exfil,TGROUP,IFLG,IRPT)
        endif
      endif

C Redisplay menu.
      goto 444
      end


C *********** RUNSELF ***************
C Read model details from the specified file and group and run it.
C LVFILE - name of the file containing model details
C LGROUP - name of group
C IFLG   - =1 select all models automatically
C          =0 allow user selection
C IRPT   - report to file (=1) or screen (=2) or unknown (=0, will ask)

      subroutine runself(LVFILE,LGROUP,IFLG,IRPT)
#include "building.h"
#include "model.h"
#include "esprdbfile.h"
#include "help.h"

      common/OUTIN/IUOUT,IUIN,IEOUT
      common/FILEP/IFIL
      COMMON/exporttg/xfile,tg,delim
      COMMON/exporttgi/ixopen,ixloc,ixunit
      common/selflstopen/summary_opnd,sum_hdr,adiffs,tdiffs,anots,mtest
      logical summary_opnd,sum_hdr
      real adiffs,tdiffs,anots
      
      parameter (MNVCASE=124)

      character name(MNVCASE)*36
      character cfgitem(MNVCASE)*72
      character rootitem(MNVCASE)*72
      character spsitem(MNVCASE)*24   ! simulation parameter set name
      character pifitem(MNVCASE)*72   ! PIF directives
      character resitem(MNVCASE,2)*72
      character repitem(MNVCASE,2)*72
      character h3kitem(MNVCASE,2)*72
      character smryitem(MNVCASE,2)*72
      character stdbin*72,tstbin*72,models_path*72
      character stdrep_path*72,tstrep_path*72
      character MODE*4
      
      dimension itemval(MNVCASE)    ! array of items to test
      dimension itemadiffs(MNVCASE) ! ascii difference found [??]
      dimension itemtdiffs(MNVCASE) ! tollerance difference found [??]
      dimension itemnot(MNVCASE)    ! item report missing [??]
      real agdiffs,tgdiffs          ! group differences

      character outstr*124,word*20,t72*72
      character LVFILE*72,LGROUP*42,tgroup*42
      character files*72,filesA*96,filesB*96,act*1
      character filesAs*96,filesBs*96
      character outs*256,doit*360
      character xfile*144,tg*1,delim*1
      character odir*84,pwd*84,the_cfg_folder*96
      character dstmp*24,vmsg*24
      
      logical CLOSE,unixok,same
      logical xstm,xstA,xstB
      logical passed_report_name
      integer irpt  ! for radio button

      helpinsub='mvalid'  ! set for subroutine

      call isunix(unixok)  ! Detect if linux and find user folder.
      odir=' '
      call usrdir(odir)
C      write(6,*) 'usrdir is ',odir(1:lnblnk(odir))

C Initialise variables.
      IVFILE=IFIL+1
      iunitA=IFIL+2
      iunitB=IFIL+3

C Values set for development testing. Normally commented out.
C      stdbin='/home/jon/esp-r_v16/bin/' 
C      tstbin='/home/jon/esp-r_jwh/bin/'
C      models_path='/home/jon/alt_tester/'
C      stdrep_path='/home/jon/alt_tester/std_reports/' 
C      tstrep_path='/home/jon/alt_tester/tst_reports/ '
      agdiffs=0.0; tgdiffs=0.0; tol = 0; iverb =0; irpt=0
      passed_report_name = .false.
      do 10 I=1,MNVCASE
        cfgitem(I)=' '; rootitem(I)=' '
        spsitem(I)=' '; pifitem(I)=' '
        resitem(I,1)=' '; resitem(I,2)=' '
        repitem(I,1)=' '; repitem(I,2)=' '
        h3kitem(I,1)=' '; h3kitem(I,2)=' '
        smryitem(I,1)=' '; smryitem(I,2)=' '
        name(I)=' '
        itemadiffs(I)=0; itemtdiffs(I)=0; itemnot(I)=0
 10   continue
 
C Open file.
      call EFOPSEQ(IVFILE,LVFILE,1,IER)
      call STRIPC(IVFILE,outstr,0,ND,0,'Standards 1st line',IER)
      if (outstr(1:9).ne.'*SELFTEST') then
        call edisp(iuout,' Not a selftest models file!')
        ier=1
        return
      endif

C File opened successfully.  Read contents line by line.
  7   call STRIPC(IVFILE,outstr,0,ND,0,'group id',IER)
      if(ier.ne.0)goto 100
      K=0
      call EGETW(outstr,K,word,'W','group id',IER)
      if (word(1:13).eq.'*bin_path_std') then
        call egetrm(outstr,K,stdbin,'W','path to std bin',IER)
      elseif (word(1:13).eq.'*bin_path_tst') then
        call egetrm(outstr,K,tstbin,'W','path to tst bin',IER)
      elseif (word(1:11).eq.'*model_path') then
        call egetrm(outstr,K,models_path,'W','path selftest models',IER)
      elseif (word(1:16).eq.'*report_path_std') then
        call egetrm(outstr,K,stdrep_path,'W','path std reports',IER)
      elseif (word(1:16).eq.'*report_path_tst') then
        call egetrm(outstr,K,tstrep_path,'W','path tst reports',IER)
      elseif (word(1:11).eq.'*tollerance') then
        call egetwr(outstr,K,tol,0.,0.,'-','tollerance',IER)
      elseif (word(1:8).eq.'*verbose') then
        call egetw(outstr,K,vmsg,'W','verbose',IER)
        if(vmsg(1:2).eq.'NO'.or.vmsg(1:2).eq.'no')then
          iverb=1
        elseif(vmsg(1:3).eq.'YES'.or.vmsg(1:3).eq.'yes')then
          iverb=2
        endif
      elseif (word(1:7).eq.'*output') then
        call egetw(outstr,K,vmsg,'W','output chnl',IER)
        if(vmsg(1:4).eq.'FILE'.or.vmsg(1:4).eq.'file')then
          irpt=1
          call egetrm(outstr,K,xfile,'W','output report',IER)
          passed_report_name = .true.
          write(6,*) 'output ',xfile(1:lnblnk(xfile))
        elseif(vmsg(1:6).eq.'SCREEN'.or.vmsg(1:6).eq.'screen')then
          irpt=2
        endif
      elseif (word(1:6).eq.'*group') then
        call egetrm(outstr,K,TGROUP,'W','group name',IER)
        IL=max(1,lnblnk(TGROUP))
        if (TGROUP(1:IL).eq.LGROUP(1:IL)) then
        
C This is our group (as selected in previous menu).  Read model details.
C Reset number of models and minimum menu width.
          inm=0
          IW=12
  39      call STRIPC(IVFILE,outstr,0,ND,0,'item id',IER)
          K=0
          call EGETW(outstr,K,word,'W','member tags',IER)
          if (word(1:5).eq.'*name') then
            inm=inm+1
            call egetrm(outstr,K,name(inm),'W','test name',IER)
            IW=max(IW,lnblnk(name(inm)))
          elseif(word(1:4).eq.'*cfg')then
            call egetrm(outstr,K,cfgitem(inm),'W','cfg name',IER)
          elseif(word(1:5).eq.'*root')then
            call egetrm(outstr,K,rootitem(inm),'W','model folder',IER)
          elseif(word(1:4).eq.'*SPS ')then
            call egetrm(outstr,K,spsitem(inm),'W','sps tag',IER)
          elseif(word(1:4).eq.'*PIF ')then
            call egetrm(outstr,K,pifitem(inm),'W','pif file',IER)
          elseif(word(1:7).eq.'*stdres')then
            call egetrm(outstr,K,resitem(inm,1),'W','std res file',IER)
          elseif(word(1:7).eq.'*tstres')then
            call egetrm(outstr,K,resitem(inm,2),'W','tst res file',IER)
          elseif(word(1:4).eq.'*res')then
            call egetrm(outstr,K,resitem(inm,1),'W','res sl4 file',IER)
            resitem(inm,2) = resitem(inm,1)  ! both use same file name
          elseif(word(1:7).eq.'*stdrep')then
            call egetrm(outstr,K,repitem(inm,1),'W','std report',IER)
          elseif(word(1:7).eq.'*tstrep')then
            call egetrm(outstr,K,repitem(inm,2),'W','tst report',IER)
          elseif(word(1:4).eq.'*rep')then
            call egetrm(outstr,K,repitem(inm,1),'W','sl4 report',IER)
            repitem(inm,2) = repitem(inm,1)  ! both use same file name
          elseif(word(1:7).eq.'*stdsum')then
            call egetrm(outstr,K,smryitem(inm,1),'W','std summary',IER)
          elseif(word(1:7).eq.'*tstsum')then
            call egetrm(outstr,K,smryitem(inm,2),'W','tst summary',IER)
          elseif(word(1:4).eq.'*sum')then
            call egetrm(outstr,K,smryitem(inm,1),'W','summary rep',IER)
            smryitem(inm,2) = smryitem(inm,1) ! both use same file name
          elseif(word(1:7).eq.'*stdh3k')then
            call egetrm(outstr,K,h3kitem(inm,1),'W','std h3k rpt',IER)
          elseif(word(1:7).eq.'*tsth3k')then
            call egetrm(outstr,K,h3kitem(inm,2),'W','tst h3k rpt',IER)
          elseif(word(1:4).eq.'*h3k')then
            call egetrm(outstr,K,h3kitem(inm,1),'W','h3k rpt',IER)
            h3kitem(inm,2) = h3kitem(inm,1)  ! both use same file name
          elseif(word(1:6).eq.'*group')then
          
C End of group - stop reading.
            goto 100
          elseif(word(1:6).eq.'*end')then
          
C End of file (hence group) - stop reading.
            goto 100
          endif
          
C Continue reading group details from file.
          goto 39
        endif
      elseif(word(1:6).eq.'*end')then
          
C End of file (hence group) - stop reading.
        goto 100
      endif
          
C Continue reading details from file - still looking for selected group.
      goto 7
      
C If we have reached the next group it is ok to present the list.
C First close file.
 100  CALL ERPFREE(IVFILE,ISTAT)
      if (INM.eq.0) then
        call usrmsg('No models found',' ','W')
        return
      endif
      
C Allow model selection if IFLG=0, otherwise all models in group will be 
C processed.
      INPIC=INM 

      if (IFLG.eq.0) then
        helptopic='self-test-groups'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL EPICKS(INPIC,itemval,' ',' Which tests?',
     &    IW,INM,name,'test list',IER,nbhelp)
      endif

C Check models have been selected.
      if (INPIC.eq.0) then
        return
      endif

C Process selected models.
      write (outs,'(2a)') 'Selected tests from group: ',
     &  LGROUP(1:lnblnk(LGROUP))
      call edisp(iuout,outs)
      do 999 I=1,inpic
        if (IFLG.eq.0) then
          imod=itemval(I)
        else
          imod=I
        endif
        write (outs,'(2a)') ' ',name(imod)(1:lnblnk(name(imod)))
        call edisp(iuout,outs)
 999  continue
      if(passed_report_name)then
        continue
      else
        if (IRPT.eq.0) then
          helptopic='selftest-reporting'
          call gethelptext(helpinsub,helptopic,nbhelp)
          irpt=1
          CALL EASKMBOX(' Reporting to:',' ','file',
     &      'screen','cancel',' ',' ',' ',' ',' ',irpt,nbhelp)
        endif
      endif
      if (irpt.eq.1) then

C If not yet opened ixunit do so. Set summary_opnd so
C subsequent use keeps appending more text to it.
 9      itu=ixunit
        if(passed_report_name)then  ! Open only on first pass.
          if(.NOT.summary_opnd)then
            call ctldumpt(xfile,ixopen,ixloc,ixunit,'exported',IER)
            summary_opnd=.true.
          endif
        else
          if(unixok)then
            write(xfile,'(a)')'/tmp/selftesta.txt'
          else
            write(xfile,'(a)')'selftests.txt'
          endif
          if(.NOT.summary_opnd)then  ! Open only on first pass.
            call ctlexp(xfile,ixopen,ixloc,ixunit,'T','exported',IER)
            summary_opnd=.true.
          endif
        endif
        
C Set irpt to -1 so that if running multiple groups this question will 
C not be asked again.
        irpt=-1
      elseif (irpt.eq.2) then
        itu=iuout                     
      elseif (irpt.eq.3) then
        return
      endif

C Echo context to the report if initial pass.
      if(.NOT.sum_hdr)then
        call dstamp(dstmp)
        write (outs,'(2a)') 'Selftest: ',dstmp(1:lnblnk(dstmp))
        call edisp(itu,outs)
        write (outs,'(2a)') 'Standard ESP-r: ',stdbin(1:lnblnk(stdbin))
        call edisp(itu,outs)
        write (outs,'(2a)') 'Test ESP-r: ',tstbin(1:lnblnk(tstbin))
        call edisp(itu,outs)
        write (outs,'(2a)') 'Standard reports: ',
     &    stdrep_path(1:lnblnk(stdrep_path))
        call edisp(itu,outs)
        write (outs,'(2a)') 'Test reports: ',
     &    tstrep_path(1:lnblnk(tstrep_path))
        call edisp(itu,outs)
        write (outs,'(2a)') 'Model folders: ',
     &    models_path(1:lnblnk(models_path))
        call edisp(itu,outs)
        if(tol.lt.0.1) tol=1.0  ! If not in selftest list file set.
        write (outs,'(a,f5.2,a,f5.2,a)') 'Tollerance: ',tol,
     &    ' W or ',tol*0.1,' GJ oC W/m2.'
        call edisp(itu,outs)
        sum_hdr = .true.
      endif
        
C The models selected are now processed. bps (simulations) start for each 
C one of them. It reads the models configuration name and the name of the 
C simulation settings from the text file.
      do 200 ipick=1,INPIC

C If user selected models then selection will be held in itemval array,
C otherwise run all models.
        if (IFLG.eq.0) then
          imod=itemval(ipick)
        else
          imod=ipick
        endif

C This will write out the titles for the results.
        write (outs,'(4a)') 'Test: ',LGROUP(1:lnblnk(LGROUP)),
     &    ' - ',name(imod)(1:lnblnk(name(imod)))
        call edisp(itu,outs)

C Change to the directory of the model
        lnmp = lnblnk(models_path)
        lnroot = lnblnk(rootitem(imod))
        write(the_cfg_folder,'(2a)') models_path(1:lnmp),
     &    rootitem(imod)(1:lnroot)
        if(iverb.eq.2)then
          write(iuout,*) '  '
          write(iuout,*) 'Path to cfg ',
     &      the_cfg_folder(1:lnblnk(the_cfg_folder))
        endif
        call chdir(the_cfg_folder(1:lnblnk(the_cfg_folder)))
        call usrdir(pwd)
        if(iverb.eq.2)then
          write(iuout,*) 'Working in folder ',pwd(1:lnblnk(pwd))
        endif

        mtest = mtest + 1  ! Increment count of tests.

C Run the std version of bps.
        lncfg = lnblnk(cfgitem(imod))
        lnsps = lnblnk(spsitem(imod))
        write (outs,'(4a)') cfgitem(imod)(1:lncfg),
     &    ' -p ',spsitem(imod)(1:lnsps),' silent'
        if(iverb.eq.2)then
          write(doit,'(3a)') stdbin(1:lnblnk(stdbin)),
     &      'bps -mode text -file ',outs(1:lnblnk(outs))
        else
          write(doit,'(4a)') stdbin(1:lnblnk(stdbin)),
     &      'bps -mode text -file ',outs(1:lnblnk(outs)),
     &      ' >/dev/null 2>&1'
        endif
        if(iverb.eq.2)then
          write(iuout,*) '  '
          write(iuout,*) 'Running std: ',doit(1:lnblnk(doit))
        else
          write(iuout,*) '  '
          write(iuout,*) 'Assessing std: ',cfgitem(imod)(1:lncfg)
        endif
        call runit(doit,'-')

C After the bps (simulation), the results recovery starts.
C It uses PIF file for save level 4 runs. If pif is 'SL5' then
C no need to run res.
        if(pifitem(imod)(1:3).eq.'SL5')then

C SL5 run post tasks: delete out.xml & out.csv
          t72 = 'out.xml'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'out.csv'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = resitem(imod,1)  ! If plant only no useful results file.
          ltf=LNBLNK(t72); xstm=.false.
          INQUIRE (FILE=t72(1:ltf),EXIST=xstm)
          if(xstm)then
            write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
            call runit(doit,'-')
          endif

C For some models there will be files with the model root name with .fcts1
C or .fcts3 ending. Clean these if found.
          lncfg = lnblnk(cfgitem(imod))
          write(files,'(2a)') cfgitem(imod)(1:lncfg),'.fcts1'
          ltf=max(1,LNBLNK(files)); xstm=.false.
          INQUIRE (FILE=files(1:ltf),EXIST=xstm)
          if(xstm)then
            write(doit,'(2a)') 'rm -f ',files(1:ltf)
            call runit(doit,'-')
          endif
          write(files,'(2a)') cfgitem(imod)(1:lncfg),'.fcts3'
          ltf=max(1,LNBLNK(files)); xstm=.false.
          INQUIRE (FILE=files(1:ltf),EXIST=xstm)
          if(xstm)then
            write(doit,'(2a)') 'rm -f ',files(1:ltf)
            call runit(doit,'-')
          endif

C Clean up flow, electrical and pant results if results.* found.
          t72 = 'results.mfr'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'results.elr'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'results.eres'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'results.plr'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'results.pres'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')

C Clean other extraneous files.
          t72 = 'fort.11'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'JAN'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'JUL'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')

C Move the .h3k file to h3kitem(imod,1) and then the out.summary to repitem.
C If plant only there will not be a h3k file so check existance.
          lncfg = lnblnk(cfgitem(imod))
          lnh3k = lnblnk(h3kitem(imod,1))
          lnstdrp = lnblnk(stdrep_path)
          write(files,'(2a)') cfgitem(imod)(1:lncfg),'.h3k'
          ltf=max(1,LNBLNK(files)); xstm=.false.
          INQUIRE (FILE=files(1:ltf),EXIST=xstm)

          write(filesA,'(2a)') stdrep_path(1:lnstdrp),
     &      h3kitem(imod,1)(1:lnh3k)
          if(xstm)then
            write(doit,'(4a)') 'mv ',files(1:lnblnk(files)),' ',
     &        filesA(1:lnblnk(filesA))
            call runit(doit,'-')
          endif
          lnsum= lnblnk(smryitem(imod,1))
          write(filesAs,'(2a)') stdrep_path(1:lnstdrp),
     &      smryitem(imod,1)(1:lnsum)
          t72 = 'out.summary'
          ltf=LNBLNK(t72); xstm=.false.
          INQUIRE (FILE=t72(1:ltf),EXIST=xstm)
          if(xstm)then
            write(doit,'(2a)') 'mv out.summary ',
     &        filesAs(1:lnblnk(filesAs))
            call runit(doit,'-')
          endif
        else

C Invoke the std version of res with PIF file recovery instructions.
          lnres = lnblnk(resitem(imod,1))
          lnpif = lnblnk(pifitem(imod))
          if(iverb.eq.2)then
            write(doit,'(5a)') stdbin(1:lnblnk(stdbin)),
     &        'res -mode text -file ',resitem(imod,1)(1:lnres),
     &        ' -act recover silent -actf ',pifitem(imod)(1:lnpif)
          else
            write(doit,'(6a)') stdbin(1:lnblnk(stdbin)),
     &        'res -mode text -file ',resitem(imod,1)(1:lnres),
     &        ' -act recover silent -actf ',pifitem(imod)(1:lnpif),
     &        ' >/dev/null 2>&1'
          endif
          if(iverb.eq.2)then
            write(iuout,*) '  '
            write(iuout,*) 'Running: ',doit(1:lnblnk(doit))
          endif
          call runit(doit,'-')

C Rename the standard output file from the PIF session (results.dat) to
C the name in repitem
          lnstdrp=lnblnk(stdrep_path)
          lnrep= lnblnk(repitem(imod,1))
          lnsum= lnblnk(smryitem(imod,1))
          write(filesA,'(2a)') stdrep_path(1:lnstdrp),
     &      repitem(imod,1)(1:lnrep)
          t72 = 'results.dat'
          ltf=LNBLNK(t72); xstm=.false.
          INQUIRE (FILE=t72(1:ltf),EXIST=xstm)
          if(xstm)then
            write(doit,'(2a)') 'mv results.dat ',
     &        filesA(1:lnblnk(filesA))     
            call runit(doit,'-')
          endif

          lnsum= lnblnk(smryitem(imod,1))
          write(filesAs,'(2a)') stdrep_path(1:lnstdrp),
     &      smryitem(imod,1)(1:lnsum)
          t72 = 'out.summary'
          ltf=LNBLNK(t72); xstm=.false.
          INQUIRE (FILE=t72(1:ltf),EXIST=xstm)
          if(xstm)then
            write(doit,'(2a)') 'mv out.summary ',
     &        filesAs(1:lnblnk(filesAs))
            call runit(doit,'-')
          endif

C If there is an input.xml file in the folder running bps will likely
C generate an out.xml, out.csv file and out.summary files.
          t72 = 'out.xml'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'out.csv'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = resitem(imod,1)   ! Delete the standard results library.
          ltf=LNBLNK(t72); xstm=.false.
          INQUIRE (FILE=t72(1:ltf),EXIST=xstm)
          if(xstm)then
            write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
            call runit(doit,'-')
          endif

C For some models there will be files with the model root name with .fcts1
C or .fcts3 ending. Clean these if found.
          lncfg = lnblnk(cfgitem(imod))
          write(files,'(2a)') cfgitem(imod)(1:lncfg),'.fcts1'
          ltf=max(1,LNBLNK(files)); xstm=.false.
          INQUIRE (FILE=files(1:ltf),EXIST=xstm)
          if(xstm)then
            write(doit,'(2a)') 'rm -f ',files(1:ltf)
            call runit(doit,'-')
          endif
          write(files,'(2a)') cfgitem(imod)(1:lncfg),'.fcts3'
          ltf=max(1,LNBLNK(files)); xstm=.false.
          INQUIRE (FILE=files(1:ltf),EXIST=xstm)
          if(xstm)then
            write(doit,'(2a)') 'rm -f ',files(1:ltf)
            call runit(doit,'-')
          endif

C Also clean up flow / electrical / plant results if found.
          t72 = 'results.mfr'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'results.elr'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'results.eres'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'results.plr'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'results.pres'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')

C Clean other extraneous files.
          t72 = 'fort.11'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'JAN'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'JUL'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
        endif

C Run the tst version of bps.
        lncfg = lnblnk(cfgitem(imod))
        lnsps = lnblnk(spsitem(imod))
        write (outs,'(4a)') cfgitem(imod)(1:lncfg),
     &    ' -p ',spsitem(imod)(1:lnsps),' silent'
        if(iverb.eq.2)then
          write(doit,'(3a)') tstbin(1:lnblnk(tstbin)),
     &      'bps -mode text -file ',outs(1:lnblnk(outs))
        else
          write(doit,'(4a)') tstbin(1:lnblnk(tstbin)),
     &      'bps -mode text -file ',outs(1:lnblnk(outs)),
     &      ' >/dev/null 2>&1'
        endif
        if(iverb.eq.2)then
          write(iuout,*) '  '
          write(iuout,*) 'Running tst: ',doit(1:lnblnk(doit))
        else
C          write(iuout,*) '  '
          write(iuout,*) 'Assessing tst: ',cfgitem(imod)(1:lncfg)
        endif
        call runit(doit,'-')

C After the bps (simulation), the results recovery starts.
C It uses PIF file for save level 4 runs. If pif is 'SL5' then
C no need to run res.
        if(pifitem(imod)(1:3).eq.'SL5')then

C SL5 run post tasks: delete out.xml & out.csv
          t72 = 'out.xml'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'out.csv'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = resitem(imod,2)
          ltf=LNBLNK(t72); xstm=.false.
          INQUIRE (FILE=t72(1:ltf),EXIST=xstm)
          if(xstm)then
            write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
            call runit(doit,'-')
          endif

C For some models there will be files with the model root name with .fcts1
C or .fcts3 ending. Clean these if found.
          lncfg = lnblnk(cfgitem(imod))
          write(files,'(2a)') cfgitem(imod)(1:lncfg),'.fcts1'
          ltf=max(1,LNBLNK(files)); xstm=.false.
          INQUIRE (FILE=files(1:ltf),EXIST=xstm)
          if(xstm)then
            write(doit,'(2a)') 'rm -f ',files(1:ltf)
            call runit(doit,'-')
          endif
          write(files,'(2a)') cfgitem(imod)(1:lncfg),'.fcts3'
          ltf=max(1,LNBLNK(files)); xstm=.false.
          INQUIRE (FILE=files(1:ltf),EXIST=xstm)
          if(xstm)then
            write(doit,'(2a)') 'rm -f ',files(1:ltf)
            call runit(doit,'-')
          endif

C Also clean up flow / electrical / plant results if found.
          t72 = 'results.mfr'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'results.elr'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'results.eres'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'results.plr'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'results.pres'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')

C Clean other extraneous files.
          t72 = 'fort.11'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'JAN'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'JUL'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')

C Move the .h3k file to h3kitem(imod,1) and then the out.summary to repitem
          lncfg = lnblnk(cfgitem(imod))
          lnh3k = lnblnk(h3kitem(imod,2))
          lntstrp = lnblnk(tstrep_path)
          write(files,'(2a)') cfgitem(imod)(1:lncfg),'.h3k'
          ltf=max(1,LNBLNK(files)); xstm=.false.
          INQUIRE (FILE=files(1:ltf),EXIST=xstm)
          write(filesB,'(2a)') tstrep_path(1:lntstrp),
     &      h3kitem(imod,2)(1:lnh3k)
          if(xstm)then
            write(doit,'(4a)') 'mv ',files(1:lnblnk(files)),' ',
     &        filesB(1:lnblnk(filesB))
            call runit(doit,'-')
          endif
          lnsum= lnblnk(smryitem(imod,2))
          write(filesBs,'(2a)') tstrep_path(1:lntstrp),
     &      smryitem(imod,2)(1:lnsum)
          t72 = 'out.summary'
          ltf=LNBLNK(t72); xstm=.false.
          INQUIRE (FILE=t72(1:ltf),EXIST=xstm)
          if(xstm)then
            write(doit,'(2a)') 'mv out.summary ',
     &        filesBs(1:lnblnk(filesBs))
            call runit(doit,'-')
          endif
        else

C Invoke the tst version of res with PIF file recovery instructions.
          lnres = lnblnk(resitem(imod,1))
          lnpif = lnblnk(pifitem(imod))
          if(iverb.eq.2)then
            write(doit,'(5a)') tstbin(1:lnblnk(tstbin)),
     &        'res -mode text -file ',resitem(imod,1)(1:lnres),
     &        ' -act recover silent -actf ',pifitem(imod)(1:lnpif)
          else
            write(doit,'(6a)') tstbin(1:lnblnk(tstbin)),
     &        'res -mode text -file ',resitem(imod,1)(1:lnres),
     &        ' -act recover silent -actf ',pifitem(imod)(1:lnpif),
     &        ' >/dev/null 2>&1'
          endif
          if(iverb.eq.2)then
            write(iuout,*) '  '
            write(iuout,*) 'Running: ',doit(1:lnblnk(doit))
          endif
          call runit(doit,'-')

C Rename the standard output file from the PIF session (results.dat) to
C the name in repitem
          lntstrep=lnblnk(tstrep_path)
          lnrep= lnblnk(repitem(imod,2))
          lnsum= lnblnk(smryitem(imod,2))
          write(filesB,'(2a)') tstrep_path(1:lntstrep),
     &      repitem(imod,2)(1:lnrep)
          t72 = 'results.dat'
          ltf=LNBLNK(t72); xstm=.false.
          INQUIRE (FILE=t72(1:ltf),EXIST=xstm)
          if(xstm)then
            write(doit,'(2a)') 'mv results.dat ',
     &        filesB(1:lnblnk(filesB))      
            call runit(doit,'-')
          endif

C If an input.xml file in the folder then an out.summary will be created.
          write(filesBs,'(2a)') tstrep_path(1:lntstrep),
     &      smryitem(imod,2)(1:lnsum)
          t72 = 'out.summary'
          ltf=LNBLNK(t72); xstm=.false.
          INQUIRE (FILE=t72(1:ltf),EXIST=xstm)
          if(xstm)then
            write(doit,'(2a)') 'mv out.summary ',
     &        filesBs(1:lnblnk(filesBs))
            call runit(doit,'-')
          endif

C If there is an input.xml file in the folder running bps will likely
C generate an out.xml, out.csv file and out.summary files.
          t72 = 'out.xml'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'out.csv'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = resitem(imod,2)  ! Delete the tst results library.
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')

C For some models there will be files with the model root name with .fcts1
C or .fcts3 ending. Clean these if found.
          lncfg = lnblnk(cfgitem(imod))
          write(files,'(2a)') cfgitem(imod)(1:lncfg),'.fcts1'
          ltf=max(1,LNBLNK(files)); xstm=.false.
          INQUIRE (FILE=files(1:ltf),EXIST=xstm)
          if(xstm)then
            write(doit,'(2a)') 'rm -f ',files(1:ltf)
            call runit(doit,'-')
          endif
          write(files,'(2a)') cfgitem(imod)(1:lncfg),'.fcts3'
          ltf=max(1,LNBLNK(files)); xstm=.false.
          INQUIRE (FILE=files(1:ltf),EXIST=xstm)
          if(xstm)then
            write(doit,'(2a)') 'rm -f ',files(1:ltf)
            call runit(doit,'-')
          endif

C Also clean up flow / electrical / plant results if found.
          t72 = 'results.mfr'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'results.elr'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'results.eres'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'results.plr'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'results.pres'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')

C Clean other extraneous files.
          t72 = 'fort.11'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'JAN'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
          t72 = 'JUL'
          write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
          call runit(doit,'-')
        endif

C If save 5 do asci comparison of h3k files if both files exist.
        if(pifitem(imod)(1:3).eq.'SL5')then
          ltf=max(1,LNBLNK(filesA)); xstA=.false.
          INQUIRE (FILE=filesA(1:ltf),EXIST=xstA)
          ltf=max(1,LNBLNK(filesB)); xstB=.false.
          INQUIRE (FILE=filesB(1:ltf),EXIST=xstB)
          if(xstA.and.xstB)then
            write(iuout,*) '  '
            act = 'v'
            call CMPASCI(itu,iunitA,iunitB,filesA,filesB,act,same,
     &        icount,IER)
            if(same)then
              lnh3k = lnblnk(h3kitem(imod,2))
              write(outs,'(3a)') 'Comparing ',h3kitem(imod,1)(1:lnh3k),
     &          ' identical'
              lnout=lnblnk(outs)
              write(iuout,*) outs(1:lnout)
              if(itu.ne.iuout)then
                call edisp(itu,outs)
              endif
            else
              lnh3k = lnblnk(h3kitem(imod,2))
              write(outs,'(3a,i5,a)')'Comparing ',
     &          h3kitem(imod,1)(1:lnh3k),
     &          ' differs by ',icount,' lines.'
              lnout=lnblnk(outs)
              write(iuout,*) outs(1:lnout)
              if(itu.ne.iuout)then
                call edisp(itu,outs)
              endif
              itemadiffs(imod) = itemadiffs(imod) +icount
              adiffs = adiffs + icount
            endif
          else  ! No h3k file produced (could be plant only)
            lnh3k = lnblnk(h3kitem(imod,2))
            write(outs,'(3a)')'Report ',
     &        h3kitem(imod,1)(1:lnh3k),' not produced.'
            lnout=lnblnk(outs)
            write(iuout,*) outs(1:lnout)
            if(itu.ne.iuout)then
              call edisp(itu,outs)
            endif
            itemadiffs(imod) = itemadiffs(imod) + 1
            itemnot(imod) = itemnot(imod) + 1
            anots = anots + 1   ! Remember missing report
            adiffs = adiffs + 1
          endif

C Now compare the report files from save level 4 if both files exist. 
C Open fileA and fileB via cmpasci to echo the differences found to channel itu.
        else
          ltf=max(1,LNBLNK(filesA)); xstA=.false.
          INQUIRE (FILE=filesA(1:ltf),EXIST=xstA)
          ltf=max(1,LNBLNK(filesB)); xstB=.false.
          INQUIRE (FILE=filesB(1:ltf),EXIST=xstB)
          if(xstA.and.xstB)then
            write(iuout,*) '  '
            act = 'v'
            call CMPASCI(itu,iunitA,iunitB,filesA,filesB,act,same,
     &        icount,IER)
            if(same)then
              lnrep= lnblnk(repitem(imod,1))
              write(outs,'(3a)') 'Comparing ',repitem(imod,1)(1:lnrep),
     &        ' identical'
              lnout=lnblnk(outs)
              write(iuout,*) outs(1:lnout)
              if(itu.ne.iuout)then
                call edisp(itu,outs)
              endif
            else
              lnrep= lnblnk(repitem(imod,1))
              write(outs,'(3a,i5,a)')'Comparing ',
     &          repitem(imod,1)(1:lnrep),
     &          ' differs by ',icount,' lines.'
              lnout=lnblnk(outs)
              write(iuout,*) outs(1:lnout)
              if(itu.ne.iuout)then
                call edisp(itu,outs)
              endif
              itemadiffs(imod) = itemadiffs(imod) +icount
              adiffs = adiffs + icount   ! Overall differences
              agdiffs = agdiffs + icount ! Group differences
            endif
          else  ! No report file produced
            lnrep= lnblnk(repitem(imod,1))
            write(outs,'(3a)')'Report ',
     &          repitem(imod,1)(1:lnrep),' not produced.'
            lnout=lnblnk(outs)
            write(iuout,*) outs(1:lnout)
            if(itu.ne.iuout)then
              call edisp(itu,outs)
            endif
            itemadiffs(imod) = itemadiffs(imod) + 1
            itemnot(imod) = itemnot(imod) + 1
            anots = anots + 1   ! Remember missing report
            adiffs = adiffs + 1
          endif
        endif

C Now compare the .summary files.  Open filesAs and filesBs.
C via cmpasci to echo the differences found to channel itu.
        ltf=max(1,LNBLNK(filesAs)); xstA=.false.
        INQUIRE (FILE=filesAs(1:ltf),EXIST=xstA)
        ltf=max(1,LNBLNK(filesBs)); xstB=.false.
        INQUIRE (FILE=filesBs(1:ltf),EXIST=xstB)
        if(xstA.and.xstB)then
          act = 'v'
          call CMPSUMRY(itu,iunitA,iunitB,filesAs,filesBs,tol,act,same,
     &      icount,IER)
          if(same)then
            lnsum= lnblnk(smryitem(imod,1))
            write(outs,'(3a)') 'Comparing ',smryitem(imod,1)(1:lnsum),
     &        ' within tollerance.'
            lnout=lnblnk(outs)
            write(iuout,*) outs(1:lnout)
            if(itu.ne.iuout)then
              call edisp(itu,outs)
            endif
          else
            lnsum= lnblnk(smryitem(imod,1))
            write(outs,'(2a,i5,a)') 'Comparing ',
     &        smryitem(imod,1)(1:lnsum),icount,
     &        ' outside tollerance.'
            lnout=lnblnk(outs)
            write(iuout,*) outs(1:lnout)
            if(itu.ne.iuout)then
              call edisp(itu,outs)
            endif
            itemtdiffs(imod) = itemtdiffs(imod) +icount
            tdiffs = tdiffs + icount
            tgdiffs = tgdiffs + icount
          endif
        else   ! summary file not created.
          lnsum= lnblnk(smryitem(imod,1))
          write(outs,'(3a)') 'Summary ',
     &        smryitem(imod,1)(1:lnsum),' not produced.'
          lnout=lnblnk(outs)
          write(iuout,*) outs(1:lnout)
          if(itu.ne.iuout)then
            call edisp(itu,outs)
          endif
          itemadiffs(imod) = itemadiffs(imod) + 1
          itemnot(imod) = itemnot(imod) + 1
          anots = anots + 1   ! Remember missing report
          adiffs = adiffs + 1
        endif

        call chdir(odir)
 
C Add overall differences to report.
 200  continue
      write(outs,'(a)') ' '
      write(iuout,*) outs(1:lnblnk(outs))
      if(itu.ne.iuout)then
        call edisp(itu,outs)
      endif
      write(outs,'(a,i3)') 'Test cases run: ',mtest
      write(iuout,*) outs(1:lnblnk(outs))
      if(itu.ne.iuout)then
        call edisp(itu,outs)
      endif
      write(outs,'(a,i5)') 'ASCII line differences (group)',
     &  int(agdiffs)
      write(iuout,*) outs(1:lnblnk(outs))
      if(itu.ne.iuout)then
        call edisp(itu,outs)
      endif
      write(outs,'(a,i5)') 'Tollerance line differences (group)',
     &  int(tgdiffs)
      write(iuout,*) outs(1:lnblnk(outs))
      if(itu.ne.iuout)then
        call edisp(itu,outs)
      endif
      write(outs,'(a,i5)') 'ASCII line differences (overall)',
     &  int(adiffs)
      write(iuout,*) outs(1:lnblnk(outs))
      if(itu.ne.iuout)then
        call edisp(itu,outs)
      endif
      write(outs,'(a,i5)') 'Tollerance line differences (overall)',
     &  int(tdiffs)
      write(iuout,*) outs(1:lnblnk(outs))
      if(itu.ne.iuout)then
        call edisp(itu,outs)
      endif
      write(outs,'(a,i5)') 'Missing reports (overall)',
     &  int(anots)
      write(iuout,*) outs(1:lnblnk(outs))
      if(itu.ne.iuout)then
        call edisp(itu,outs)
      endif

C If writing to file close it on last model group.
      if (irpt.eq.1.and.iflg.eq.-1)then
        CALL ERPFREE(ixunit,ISTAT)
      elseif (irpt.eq.-1.and.iflg.eq.-1)then
        CALL ERPFREE(ixunit,ISTAT)
      endif

C If writing to file close it if only running one group.
      if (irpt.eq.1.and.iflg.eq.0)then
        CALL ERPFREE(ixunit,ISTAT)
      elseif (irpt.eq.-1.and.iflg.eq.0)then
        CALL ERPFREE(ixunit,ISTAT)
      endif
      RETURN
      END

C ******************** CMPASCI ********************
C Line by line compare of two ASCII files. act = '-' do silents returning logical same.
C act = 'v' display differences, act = 'c' return icount. If itu is the same as
C iuout then display will be to text feedback.

      SUBROUTINE CMPASCI(itu,iunitA,iunitB,LFILA,LFILB,act,same,
     &  icount,IER)
      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      character LFILA*96,LFILB*96
      character outstrA*248,outstrB*248,act*1
      character msgA*320,msgB*320
      logical :: XSTA,XSTB
      logical :: unixok,same

C Set return values.
      same = .true.; icount = 0

C Open flila on unit iunita.
      call isunix(unixok)
      if(LFILA(1:1).eq.'/')then
        ltf=max(1,LNBLNK(LFILA))
        INQUIRE (FILE=LFILA(1:ltf),EXIST=XSTA)
        IF (XSTA) THEN
          LN=max(1,lnblnk(LFILA))
          OPEN (IUNITA,FILE=LFILA(1:LN),ACCESS='SEQUENTIAL',
     &              STATUS='OLD',IOSTAT=ISTAT)
        else
          call usrmsg('CMPASCI: problem opening file!',LFILA,'W')
          same = .false.; icount = icount + 1
          return
        endif
      else
        if(.NOT.unixok)then   ! For W10 use FPOPEN
          CALL ERPFREE(IUNITA,ISTAT)
          call FPOPEN(IUNITA,ISTAT,1,0,LFILA)
          if(ISTAT.ge.0)XSTA=.true.
          if(ISTAT.eq.-301)XSTA=.false.
        else
          CALL EFOPSEQ(IUNITA,LFILA,1,IER)
        endif
      endif
      IF(IER.LT.0)THEN
        call usrmsg(' Problem detected while trying to open',LFILA,'W')
        same = .false.; icount = icount + 1
        return
      ENDIF

C Open flilB on unit iunitB.
      if(LFILB(1:1).eq.'/')then
        ltf=max(1,LNBLNK(LFILB))
        INQUIRE (FILE=LFILB(1:ltf),EXIST=XSTB)
        IF (XSTB) THEN
          LN=max(1,lnblnk(LFILB))
          OPEN (IUNITB,FILE=LFILB(1:LN),ACCESS='SEQUENTIAL',
     &              STATUS='OLD',IOSTAT=ISTAT)
        else
          call usrmsg('CMPASCI: problem opening file!',LFILB,'W')
          same = .false.; icount = icount + 1
          return
        endif
      else
        if(.NOT.unixok)then   ! For W10 use FPOPEN
          CALL ERPFREE(IUNITB,ISTAT)
          call FPOPEN(IUNITB,ISTAT,1,0,LFILB)
          if(ISTAT.ge.0)XSTB=.true.
          if(ISTAT.eq.-301)XSTB=.false.
        else
          CALL EFOPSEQ(IUNITB,LFILB,1,IER)
        endif
      endif
      IF(IER.LT.0)THEN
        call usrmsg(' Problem detected while trying to open',LFILB,'W')
        return
      ENDIF

C Loop line by line through two ascii files. 

C Read a line from each ASCII file.
    7 READ(IUNITA,10,IOSTAT=IOS,END=102)outstrA
      READ(IUNITB,10,IOSTAT=IOS,END=102)outstrB
   10 FORMAT(A248)
      lnA = lnblnk(outstrA)    ! Check lengths
      lnB = lnblnk(outstrB)
      lnAA = lnA; lnBB = lnB
      if(lnAA.gt.315) lnAA=315 ! Avoid buffer overflow
      if(lnBB.gt.315) lnBB=315
      if(lnA.ne.lnB)then       ! If lengths differ and act = -v display
        icount = icount +1
        write(msgA,'(2a)') 'A: ',outstrA(1:lnAA)
        write(msgB,'(2a)') 'B: ',outstrB(1:lnBB)
        same = .false.
      else
        if(outstrA(1:lnA).eq.outstrB(1:lnB))then
          goto 7
        else
          icount = icount +1
          write(msgA,'(2a)') 'A: ',outstrA(1:lnAA)
          write(msgB,'(2a)') 'B: ',outstrB(1:lnBB)
          same = .false.
        endif
      endif
      if(act.eq.'v')then
        if(itu.eq.iuout) then
          call edisp248(iuout,msgA,100)
          call edisp248(iuout,msgB,100)
        else
          write(itu,'(a)') msgA(1:lnblnk(msgA))
          write(itu,'(a)') msgB(1:lnblnk(msgB))
        endif
      endif
      goto 7

C End of file encountered, close the files and return to menu.
  102 if(IOS.eq.2)then
        CALL USRMSG(' ','CMPASCI: end of file problem!','-')
      else
        continue
      endif
      CALL ERPFREE(IUNITA,ISTAT)
      CALL ERPFREE(IUNITB,ISTAT)
      return
      end


C ******************** CMPSUMRY ********************
C Line by line comparison of two .summary files (produced by bps when there is
C an input.xml file in the folder). act = '-' do silent scan returning logical
C same and icount. act = 'v' display differences. If itu is the same as
C iuout then display will be to text feedback.
C tol is tollerance for differences in value e.g. 1.0 will treate as identical
C differences less than 1W for lines with (W) units or 0.1 GJ or 0.1 degC.
C If tol = 10 then 10W for lines with (W) or 1GJ or 1 degC.

C Example of a .summary file line: 
C building/all_zones/energy_balance/net::Total_Average 388.777924 (W)


      SUBROUTINE CMPSUMRY(itu,iunitA,iunitB,LFILA,LFILB,tol,act,same,
     &  icount,IER)
      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      character LFILA*96,LFILB*96
      character outstrA*248,outstrB*248,act*1
      character msgA*248,msgB*248,unitA*8,unitB*8
      logical :: XSTA,XSTB
      logical :: unixok,same

C Set return values.
      same = .true.; icount = 0

C Open flila on unit iunita.
      call isunix(unixok)
      if(LFILA(1:1).eq.'/')then
        ltf=max(1,LNBLNK(LFILA))
        INQUIRE (FILE=LFILA(1:ltf),EXIST=XSTA)
        IF (XSTA) THEN
          LN=max(1,lnblnk(LFILA))
          OPEN (IUNITA,FILE=LFILA(1:LN),ACCESS='SEQUENTIAL',
     &              STATUS='OLD',IOSTAT=ISTAT)
        else
          call usrmsg('CMPSUMRY: problem opening file!',LFILA,'W')
          same = .false.
          icount = icount +1
          return
        endif
      else
        if(.NOT.unixok)then   ! For W10 use FPOPEN
          CALL ERPFREE(IUNITA,ISTAT)
          call FPOPEN(IUNITA,ISTAT,1,0,LFILA)
          if(ISTAT.ge.0)XSTA=.true.
          if(ISTAT.eq.-301)XSTA=.false.
        else
          CALL EFOPSEQ(IUNITA,LFILA,1,IER)
        endif
      endif
      IF(IER.LT.0)THEN
        call usrmsg(' Problem detected while trying to open',LFILA,'W')
        same = .false.
        icount = icount +1
        return
      ENDIF

C Open flilB on unit iunitB.
      if(LFILB(1:1).eq.'/')then
        ltf=max(1,LNBLNK(LFILB))
        INQUIRE (FILE=LFILB(1:ltf),EXIST=XSTB)
        IF (XSTB) THEN
          LN=max(1,lnblnk(LFILB))
          OPEN (IUNITB,FILE=LFILB(1:LN),ACCESS='SEQUENTIAL',
     &              STATUS='OLD',IOSTAT=ISTAT)
        else
          call usrmsg('CMPSUMRY: problem opening file!',LFILB,'W')
          same = .false.
          icount = icount +1
          return
        endif
      else
        if(.NOT.unixok)then   ! For W10 use FPOPEN
          CALL ERPFREE(IUNITB,ISTAT)
          call FPOPEN(IUNITB,ISTAT,1,0,LFILB)
          if(ISTAT.ge.0)XSTB=.true.
          if(ISTAT.eq.-301)XSTB=.false.
        else
          CALL EFOPSEQ(IUNITB,LFILB,1,IER)
        endif
      endif
      IF(IER.LT.0)THEN
        call usrmsg(' Problem detected while trying to open',LFILB,'W')
        same = .false.
        icount = icount +1
        return
      ENDIF

C Loop line by line through two .summary files. 

C Read a line from each ASCII file.
    7 READ(IUNITA,10,IOSTAT=IOS,END=102)outstrA
      READ(IUNITB,10,IOSTAT=IOS,END=102)outstrB
   10 FORMAT(A248)
      lnA = lnblnk(outstrA)  ! Check lengths
      lnB = lnblnk(outstrB)
      if(outstrA(1:lnA).eq.outstrB(1:lnB))then
        goto 7     ! Line are identical read another.
      else
        write(msgA,'(2a)') 'A: ',outstrA(1:lnA)
        write(msgB,'(2a)') 'B: ',outstrB(1:lnB)
        K=0
        call findspace(outstrA,K,ispace,ACT,MSGA,ier)
        K=ispace
        CALL EGETWR(outstrA,K,VALA,0.,0.,'-','A value',IER)
        CALL EGETW(outstrA,K,unitA,'W','unit for A',IER)
        K=0
        call findspace(outstrB,K,ispace,ACT,MSGB,ier)
        K=ispace
        CALL EGETWR(outstrB,K,VALB,0.,0.,'-','B value',IER)
        CALL EGETW(outstrB,K,unitB,'W','unit for B',IER)
C        write(6,*) 'A comparison',VALA,unitA,B comparison',VALB,unitB
        if(unitA(1:3).eq.'(W)')then
          difference = abs(VALB-VALA)
          if(difference.gt.tol)then
            icount = icount +1
            same = .false.
          else
            goto 7    ! Within tol read another line.
          endif
        elseif(unitA(1:4).eq.'(GJ)')then
          difference = abs(VALB-VALA)
          if(difference.gt.(tol*0.001))then  ! adjust scale for (GJ)
            icount = icount +1
            same = .false.
          else
            goto 7    ! Within tol read another line.
          endif
        elseif(unitA(1:4).eq.'(oC)')then
          difference = abs(VALB-VALA)
          if(difference.gt.(tol*0.1))then  ! adjust scale for (oC)
            icount = icount +1
            same = .false.
          else
            goto 7    ! Within tol read another line.
          endif
        elseif(unitA(1:6).eq.'(W/m2)')then
          difference = abs(VALB-VALA)
          if(difference.gt.(tol*0.1))then  ! adjust scale for (W/m2)
            icount = icount +1
            same = .false.
          else
            goto 7    ! Within tol read another line.
          endif
        elseif(unitA(1:6).eq.'(W/(m2')then
          difference = abs(VALB-VALA)
          if(difference.gt.(tol*0.1))then  ! adjust scale for (W/(m2) K))
            icount = icount +1
            same = .false.
          else
            goto 7    ! Within tol read another line.
          endif
        elseif(unitA(1:3).eq.'(-)')then
          difference = abs(VALB-VALA)
          if(difference.gt.tol)then
            icount = icount +1
            same = .false.
          else
            goto 7    ! Within tol read another line.
          endif
        endif
      endif
      if(act.eq.'v')then
        if(itu.eq.iuout) then
          call edisp248(iuout,msgA,100)
          call edisp248(iuout,msgB,100)
        else
          write(itu,'(a)') msgA(1:lnblnk(msgA))
          write(itu,'(a)') msgB(1:lnblnk(msgB))
        endif
      endif
      goto 7

C End of file encountered, close the files and return to menu.
  102 if(IOS.eq.2)then
        CALL USRMSG(' ','CMPSUMRY: end of file problem!','-')
      else
        continue
      endif
      CALL ERPFREE(IUNITA,ISTAT)
      CALL ERPFREE(IUNITB,ISTAT)
      return
      end


C ******************** Findspace ********************
C returns position of the first space or tab or comma in a string.
C Provides a warning message if ACT='W', a failure message if ACT='F' and does
C no message if ACT='-'.

      SUBROUTINE findspace(STRING,K,ispace,ACT,MSG,ier)
#include "espriou.h"
      common/OUTIN/IUOUT,IUIN,IEOUT

      CHARACTER*(*) STRING, MSG
      CHARACTER ACT*1,A*1,loutstr*248,outs*124
      logical unixok

C LS is th maximum length of STRING, L the current position,
C LW the maximum length of PHRASE.
      ier=0
      LS=LEN(STRING)
      L=0
      ispace = 0

C Traverse characters.
   10 K=K+1
      IF(K.GT.LS) GOTO 999
      A=STRING(K:K)
      if(A.EQ.' '.OR.A.EQ.','.OR.(ICHAR(A).EQ.9))then
        ispace=K
      else
        GOTO 10
      endif

  100 RETURN

  999 if(ACT.EQ.'-')then
        RETURN
      elseif(ACT.EQ.'W')then
        ier=1
        LN=max(1,lnblnk(currentfile))
        LNM=max(1,lnblnk(MSG))
        if(currentfile(1:2).ne.'  ')then
          WRITE(loutstr,'(5a)',IOSTAT=IOS,ERR=1)
     &      'WARNING: in ',currentfile(1:LN),
     &      ' past end of line for ',MSG(1:LNM),' (phrase) in...'
        else
          WRITE(loutstr,'(3a)',IOSTAT=IOS,ERR=1)
     &      'WARNING: Past end of line for ',
     &      MSG(1:LNM),' (phrase) in...'
        endif
      elseif(ACT.EQ.'F')then
        LN=max(1,lnblnk(currentfile))
        LNM=max(1,lnblnk(MSG))
        if(currentfile(1:2).ne.'  ')then
          WRITE(loutstr,'(5a)',IOSTAT=IOS,ERR=1)
     &      'FAILURE: in ',currentfile(1:LN),
     &      ' past end of line for ',MSG(1:LNM),' (phrase) in...'
        else
          WRITE(loutstr,'(3a)',IOSTAT=IOS,ERR=1)
     &      'FAILURE: Past end of line for ',
     &      MSG(1:LNM),' (phrase) in...'
        endif
      endif
      CALL EDISP248(iuout,loutstr,100)
      if(lnblnk(STRING).gt.123)then
        CALL EDISP248(iuout,STRING,100)
      else
        CALL EDISP(iuout,STRING)
      endif
      GOTO 100

   1  call isunix(unixok)
      if(.NOT.unixok) return  ! if DOS return because of lack of I/O channel
      LNM=max(1,lnblnk(STRING))
      if(IOS.eq.2)then
        write(outs,*) 'findspace: permission error writing warning re:'
        call edisp(iuout,outs)
        call edisp(iuout,STRING(1:LNM))
      else
        write(outs,*) 'findspace: error writing warning re: '
        call edisp(iuout,outs)
        call edisp(iuout,STRING(1:LNM))
      endif
      return

      END


C ************** exemplartest ***************
C Scan the [exemplars_selftest_list] file, parse information and offer
C the user selection of or information about one or more models.
C iexfil is the file unit number.
C exfil is the exemplars_selftest_list file name. 
      subroutine exemplartest
#include "espriou.h"
#include "help.h"

C espriou.h provides currentfile.
      COMMON/SPAD/MMOD,LIMIT,LIMTTY
      common/OUTIN/IUOUT,IUIN,IEOUT
      common/FILEP/IFIL
      common/exemtlst/exemtlbl,exemtfl
      character exemtlbl*20,exemtfl*96
      common/selflstopen/summary_opnd,sum_hdr,adiffs,tdiffs,anots,mtest
      logical summary_opnd,sum_hdr
      real adiffs,tdiffs,anots

C name is array of string holding *name from the file.
C inl2ngr is a pointer from the menu item back to the group index.
C group is array of strings selname(i)='UNKNOWN' holding *groups in the file.
C gitem is menu string array for high level group selection.
C gh() is array of help entries.  Variable names same as in selftest.

      dimension group(25),gitem(30)
      dimension inl2ngr(50)
      character gitem*44,outstr*124,word*20
      character exfil*96,group*42,tgroup*42,key*1
      character separ*42
      logical XST
      integer irpt  ! for radio button
      integer NCOG,ICOG ! max items and current menu item

#ifdef OSI
      integer impx,impy,iwe
#else
      integer*8 impx,impy,iwe
#endif

      helpinsub='mvalid'  ! set for subroutine
      helptopic='exemplar-checks-overview'
      call gethelptext(helpinsub,helptopic,nbhelp)

      currentfile='  '; exemtfl='  '
      summary_opnd = .false.; sum_hdr= .false.
      adiffs = 0; tdiffs = 0; anots = 0; mtest = 0
      XST=.false.
      ISTRW=96
      CALL EASKF(exemtfl,' ','Exemplar list file?',
     &   ISTRW,'exemplar_selftest_list','Exemplar list',IER,nbhelp)
      
      if(exemfl(1:2).ne.'  '.and.exemfl(1:4).ne.'UNKN')then
        INQUIRE (FILE=exemfl,EXIST=XST)
        if (.not.XST) then
          call edisp(iuout,'No exemplar list available (file missing)')
          call edisp(iuout,'returning to main menu.')
          return
        endif
      else
        call edisp(iuout,'No exemplar list available (undefined)')
        call edisp(iuout,'returning to main menu.')
        return
      endif
      exfil=exemtfl
      write(currentfile,'(a)') exemtfl(1:lnblnk(exemtfl))
      
C Open the exemplar models list file.
      IER=0
      iexfil=IFIL+1
      CALL EFOPSEQ(iexfil,exfil,1,IER)
      if(ier.ne.0) then
        call edisp(iuout,'Cannot open exemplars description file')
        call edisp(iuout,'[exemplars_selftest_list] returning...')
        ier=1
        return
      endif
      
C Read the header and see if '*EXEMPLARTEST'.
      CALL STRIPC(iexfil,outstr,0,ND,0,'selftest header',IER)
      if(outstr(1:13).ne.'*EXEMPLARTEST')then
        call edisp(iuout,' Not a exemplars models list file!')
        ier=1
        return
      endif
 
C ngr is the number of groups found in the file (max 15).
C ing is index of menu string to be written to. This is separate from
C     inl so that labels do not count in the generation of KEY characters.
C inl is the total number of groups and labels to display.
      ngr=0
      ing=0
      inl=0
      inl2ngr(1)=0
      igrsel=0

 12   CALL STRIPC(iexfil,outstr,0,ND,0,'group id',IER)
      if(ier.ne.0)goto 901
      K=0
      CALL EGETW(outstr,K,word,'W','group id',IER)
      if(word(1:6).eq.'*group')then

C If we have less than 15 groups, increment ngr, inl, ing. Set
C KEY character using ing value, write menu string gitem and set
C pointer inl2ngr for menu line inl to ngr value.
        if(ngr.lt.15)then
          ngr=ngr+1
          inl=inl+1
          ing=ing+1
          CALL EMKEY(ing,KEY,IER)
          call egetrm(outstr,K,group(ngr),'W','group id',IER)
          write(gitem(inl),'(a1,1x,a42)') key, group(ngr)(1:42)
          inl2ngr(inl)=ngr
        else
          call edisp(iuout,' Too many groups!')
        endif
      elseif(word(1:6).eq.'*label')then

C A label included. Only increment inl and set pointer to zero.
        inl=inl+1
        separ=' '
        call egetrm(outstr,K,separ,'W','label',IER)
        ln=max(1,lnblnk(separ))
        ln=min0(ln,42)
        write(gitem(inl),'(2x,a)') separ(1:ln)
        inl2ngr(inl)=0
      endif
      
C Loop unti EOF error detected above (jumps to 901 and closes file).
      goto 12
 901  CALL ERPFREE(iexfil,ISTAT)

C Present the display, check the maximum width of gitem and
C use vwmenu or emenu depeding on graphic or text display.
 444  ICOG= -2
      gitem(inl+1)=  '  _____________________ '
      gitem(inl+2)=  '* all available models  '
      gitem(inl+3)=  '  _____________________ '
      gitem(inl+4)=  '? help                  '
      gitem(inl+5)=  '- exit menu'
      NCOG=inl+5

C Help @ validation list file level.
      helptopic='selftest-overview'
      call gethelptext(helpinsub,helptopic,nbhelp)
      IW=0
      do 141 ij=1,NCOG
        IW=max(lnblnk(gitem(ij)),IW)
 141  continue
 
C If in graphics mode use a variable width menu.
      if(MMOD.EQ.8)then
        impx=0
        impy=0
        iwe=iw
        CALL VWMENU('Exemplar sets',gitem,NCOG,impx,impy,iwe,
     &     irpx,irpy,ICOG)
      else
        CALL EMENU('Exemplar sets',gitem,NCOG,ICOG)
      endif
      
C Deal with the user selection.
      IF(ICOG.EQ.NCOG)THEN

C Exit.  Close file and return to calling menu.
        CALL ERPFREE(iexfil,ISTAT)
        return
      ELSEIF(ICOG.EQ.(NCOG-1))THEN

C Display help messages.
        CALL PHELPD('exemplar list menu',nbhelp,'-',0,0,IER)
      ELSEIF(ICOG.EQ.(NCOG-3))THEN

C Automatically run all of the available models. Loop thrugh each of
C the groups in the menu. 
        mtest = 0
        irpt=0
        if(irpt.ne.3)then
          do 200 IGRP=1,NGR
            TGROUP=group(IGRP)
            IFLG=1
            if (IGRP.eq.NGR) IFLG=-1
            write(6,*) 'ngr tgroup iflg',ngr,' ',tgroup,iflg
            call runself(exfil,TGROUP,IFLG,IRPT)  ! [??]
 200      continue
        endif

C Select specific models from a group.
      ELSEIF(ICOG.GE.1.AND.ICOG.LT.(NCOG-4))THEN
        if(gitem(ICOG)(1:2).eq.'  ')goto 444
        if(inl2ngr(ICOG).gt.0.and.inl2ngr(ICOG).le.25)then
          mtest = 0
          igrsel=inl2ngr(ICOG)
          IFLG=0
          IRPT=0
          TGROUP=group(igrsel)
          call runexemplars(exfil,TGROUP,IFLG,IRPT)  ! [??]
        endif
      endif

C Redisplay menu.
      goto 444
      end


C *********** RUNEXEMPLARS ***************
C Read model details from the specified exemplars file and group and run it.
C LVFILE - name of the file containing list of model details
C LGROUP - name of group
C IFLG   - =1 select all models automatically
C          =0 allow user selection
C IRPT   - report to file (=1) or screen (=2) or unknown (=0, will ask)

      subroutine runexemplars(LVFILE,LGROUP,IFLG,IRPT)
#include "building.h"
#include "model.h"
#include "esprdbfile.h"
#include "help.h"

      common/OUTIN/IUOUT,IUIN,IEOUT
      common/FILEP/IFIL
      COMMON/exporttg/xfile,tg,delim
      COMMON/exporttgi/ixopen,ixloc,ixunit
      common/selflstopen/summary_opnd,sum_hdr,adiffs,tdiffs,anots,mtest
      logical summary_opnd,sum_hdr
      real adiffs,tdiffs,anots
      
      parameter (MNVCASE=124)

      character name(MNVCASE)*36
      character cfgitem(MNVCASE)*72
      character rootitem(MNVCASE)*72
      character spsitem(MNVCASE)*24   ! simulation parameter set name
      character pifitem(MNVCASE)*72   ! PIF directives
      character resitem(MNVCASE)*72
      character repitem(MNVCASE)*72
      character stdbin*72,tstbin*72,models_path*72
      character stdrep_path*72
      
      dimension itemval(MNVCASE)    ! array of items to test
      dimension itemadiffs(MNVCASE) ! ascii difference found [??]
      dimension itemtdiffs(MNVCASE) ! tollerance difference found [??]
      dimension itemnot(MNVCASE)    ! item report missing [??]
      real agdiffs,tgdiffs          ! group differences and missing reports

      character outstr*124,word*20,t72*72
      character LVFILE*72,LGROUP*42,tgroup*42
      character files*72,filesA*96,filesB*96,act*1
      character outs*256,doit*360
      character xfile*144,tg*1,delim*1
      character odir*84,pwd*84,the_cfg_folder*96
      character dstmp*24,vmsg*24
      
      logical XST,unixok,same
      logical passed_report_name
      integer irpt  ! for radio button

      helpinsub='mvalid'  ! set for subroutine

      call isunix(unixok)  ! Detect if linux and find user folder.
      odir=' '
      call usrdir(odir)
      write(6,*) 'usrdir is ',odir(1:lnblnk(odir))

C Initialise variables. Edit stdbin models_path etc for development testing.
      IVFILE=IFIL+1
      iunitA=IFIL+2
      iunitB=IFIL+3
C      stdbin='/home/jon/esp-r_jwh/bin/' 
C      models_path='/home/jon/alt_exemplars/'
C      stdrep_path='/home/jon/alt_exemplars/std_reports/' 
      agdiffs=0.0; tgdiffs=0.0; tol = 0; iverb =0; irpt=0
      passed_report_name = .false.
      filesA=' '; filesB=' '
      do 10 I=1,MNVCASE
        cfgitem(I)=' '; rootitem(I)=' '
        spsitem(I)=' '; pifitem(I)=' '
        resitem(I)=' '; repitem(I)=' '
        name(I)=' '
        itemadiffs(I)=0; itemtdiffs(I)=0; itemnot(I)=0
 10   continue
 
C Open file.
      call EFOPSEQ(IVFILE,LVFILE,1,IER)
      call STRIPC(IVFILE,outstr,0,ND,0,'Exemplars 1st line',IER)
      if (outstr(1:13).ne.'*EXEMPLARTEST') then
        call edisp(iuout,' Not an exemplars models list file!')
        ier=1
        return
      endif

C File opened successfully.  Read contents line by line.
  7   call STRIPC(IVFILE,outstr,0,ND,0,'group id',IER)
      if(ier.ne.0)goto 100
      K=0
      call EGETW(outstr,K,word,'W','group id',IER)
      if (word(1:9).eq.'*bin_path') then
        call egetrm(outstr,K,stdbin,'W','path to std bin',IER)
      elseif (word(1:11).eq.'*model_path') then
        call egetrm(outstr,K,models_path,'W','path exemplar models',IER)
      elseif (word(1:16).eq.'*report_path') then
        call egetrm(outstr,K,stdrep_path,'W','path reports',IER)
      elseif (word(1:11).eq.'*tollerance') then
        call egetwr(outstr,K,tol,0.,0.,'-','tollerance',IER)
      elseif (word(1:8).eq.'*verbose') then
        call egetw(outstr,K,vmsg,'W','verbose',IER)
        if(vmsg(1:2).eq.'NO'.or.vmsg(1:2).eq.'no')then
          iverb=1
        elseif(vmsg(1:3).eq.'YES'.or.vmsg(1:3).eq.'yes')then
          iverb=2
        endif
      elseif (word(1:7).eq.'*output') then
        call egetw(outstr,K,vmsg,'W','output chnl',IER)
        if(vmsg(1:4).eq.'FILE'.or.vmsg(1:4).eq.'file')then
          irpt=1
          call egetrm(outstr,K,xfile,'W','output report',IER)
          passed_report_name = .true.
          write(6,*) 'output ',xfile(1:lnblnk(xfile))
        elseif(vmsg(1:6).eq.'SCREEN'.or.vmsg(1:6).eq.'screen')then
          irpt=2
        endif
      elseif (word(1:6).eq.'*group') then
        call egetrm(outstr,K,TGROUP,'W','group name',IER)
        IL=max(1,lnblnk(TGROUP))
        if (TGROUP(1:IL).eq.LGROUP(1:IL)) then
        
C This is our group (selected in previous menu).  Read model details.
C Reset number of models and minimum menu width.
          inm=0
          IW=12
  39      call STRIPC(IVFILE,outstr,0,ND,0,'item id',IER)
          K=0
          call EGETW(outstr,K,word,'W','member tags',IER)
          if (word(1:5).eq.'*name') then
            inm=inm+1
            call egetrm(outstr,K,name(inm),'W','test name',IER)
            IW=max(IW,lnblnk(name(inm)))
          elseif(word(1:4).eq.'*cfg')then
            call egetrm(outstr,K,cfgitem(inm),'W','cfg name',IER)
          elseif(word(1:5).eq.'*root')then
            call egetrm(outstr,K,rootitem(inm),'W','model folder',IER)
          elseif(word(1:4).eq.'*SPS ')then
            call egetrm(outstr,K,spsitem(inm),'W','sps tag',IER)
          elseif(word(1:4).eq.'*PIF ')then
            call egetrm(outstr,K,pifitem(inm),'W','pif file',IER)
          elseif(word(1:4).eq.'*res')then
            call egetrm(outstr,K,resitem(inm),'W','res sl4 file',IER)
          elseif(word(1:4).eq.'*rep')then
            call egetrm(outstr,K,repitem(inm),'W','sl4 report',IER)
          elseif(word(1:6).eq.'*group')then
          
C End of group - stop reading.
            goto 100
          elseif(word(1:6).eq.'*end')then
          
C End of file (hence group) - stop reading.
            goto 100
          endif
          
C Continue reading group details from file.
          goto 39
        endif
      elseif(word(1:6).eq.'*end')then
          
C End of file (hence group) - stop reading.
        goto 100
      endif
          
C Continue reading details from file - still looking for selected group.
      goto 7
      
C If we have reached the next group it is ok to present the list.
C First close file.
 100  CALL ERPFREE(IVFILE,ISTAT)
      if (INM.eq.0) then
        call usrmsg('No models found',' ','W')
        return
      endif
      
C Allow model selection if IFLG=0, otherwise all models in group will be 
C processed.
      INPIC=INM 

      if (IFLG.eq.0) then
        helptopic='exemplar-test-groups'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL EPICKS(INPIC,itemval,' ',' Which tests?',
     &    IW,INM,name,'test list',IER,nbhelp)
      endif

C Check models have been selected.
      if (INPIC.eq.0) then
        return
      endif

C Process selected models.
      write (outs,'(2a)') 'Selected tests from group: ',
     &  LGROUP(1:lnblnk(LGROUP))
      call edisp(iuout,outs)
      do 999 I=1,inpic
        if (IFLG.eq.0) then
          imod=itemval(I)
        else
          imod=I
        endif
        write (outs,'(2a)') ' ',name(imod)(1:lnblnk(name(imod)))
        call edisp(iuout,outs)
 999  continue
      if(passed_report_name)then
        continue
      else
        if (IRPT.eq.0) then
          helptopic='exemplars-test-reporting'
          call gethelptext(helpinsub,helptopic,nbhelp)
          irpt=1
          CALL EASKMBOX(' Reporting to:',' ','file',
     &      'screen','cancel',' ',' ',' ',' ',' ',irpt,nbhelp)
        endif
      endif
      if (irpt.eq.1) then

C If not yet opened ixunit do so. Set summary_opnd so
C subsequent use keeps appending more text to it.
        itu=ixunit
        if(passed_report_name)then  ! Open only on first pass.
          if(.NOT.summary_opnd)then
            call ctldumpt(xfile,ixopen,ixloc,ixunit,'exported',IER)
            summary_opnd=.true.
          endif
        else
          if(unixok)then
            write(xfile,'(a)')'/tmp/all_exemplars.txt'
          else
            write(xfile,'(a)')'all_exemplars.txt'
          endif
          if(.NOT.summary_opnd)then  ! Open only on first pass.
            call ctlexp(xfile,ixopen,ixloc,ixunit,'T','exported',IER)
            summary_opnd=.true.
          endif
        endif
        
C Set irpt to -1 so that if running multiple groups this question will 
C not be asked again.
        irpt=-1
      elseif (irpt.eq.2) then
        itu=iuout                     
      elseif (irpt.eq.3) then
        return
      endif

C Echo context to the report if initial pass.
      if(.NOT.sum_hdr)then
        call dstamp(dstmp)
        write (outs,'(2a)') 'Exemplars: ',dstmp(1:lnblnk(dstmp))
        call edisp(itu,outs)
        write (outs,'(2a)') 'ESP-r: ',stdbin(1:lnblnk(stdbin))
        call edisp(itu,outs)
        write (outs,'(2a)') 'Reports: ',
     &    stdrep_path(1:lnblnk(stdrep_path))
        call edisp(itu,outs)
        write (outs,'(2a)') 'Model folders: ',
     &    models_path(1:lnblnk(models_path))
        call edisp(itu,outs)
        if(tol.lt.0.1) tol=1.0  ! If not in exemplars list file set.
        write (outs,'(a,f5.2,a,f5.2,a)') 'Tollerance: ',tol,
     &    ' W or ',tol*0.1,' GJ oC W/m2.'
        call edisp(itu,outs)
        sum_hdr = .true.
      endif
        
C The models selected are now processed. bps (simulations) start for each 
C one of them. It reads the models configuration name and the name of the 
C simulation settings from the text file.
      do 200 ipick=1,INPIC

C If user selected models then selection will be held in itemval array,
C otherwise run all models.
        if (IFLG.eq.0) then
          imod=itemval(ipick)
        else
          imod=ipick
        endif

C This will write out the titles for the results.
        write (outs,'(4a)') 'Test: ',LGROUP(1:lnblnk(LGROUP)),
     &    ' - ',name(imod)(1:lnblnk(name(imod)))
        call edisp(itu,outs)
        call edisp(iuout,outs)

C Change to the directory of the model
        lnmp = lnblnk(models_path)
        lnroot = lnblnk(rootitem(imod))
        write(the_cfg_folder,'(2a)') models_path(1:lnmp),
     &    rootitem(imod)(1:lnroot)
        if(iverb.eq.2)then
          write(iuout,*) '  '
          write(iuout,*) 'Path to cfg ',
     &      the_cfg_folder(1:lnblnk(the_cfg_folder))
        endif
        call chdir(the_cfg_folder(1:lnblnk(the_cfg_folder)))
        call usrdir(pwd)
        if(iverb.eq.2)then
          write(iuout,*) 'Working in folder ',pwd(1:lnblnk(pwd))
        endif

C Scan the cfg file and then loop through the zones to see if shading
C files exist via checkmodel. Run the std version of bps.
        lncfg = lnblnk(cfgitem(imod))
        lnsps = lnblnk(spsitem(imod))
        write (outs,'(4a)') cfgitem(imod)(1:lncfg),
     &    ' -p ',spsitem(imod)(1:lnsps),' silent'
        if(iverb.eq.2)then
          write(doit,'(3a)') stdbin(1:lnblnk(stdbin)),
     &      'bps -mode text -file ',outs(1:lnblnk(outs))
        else
          write(doit,'(4a)') stdbin(1:lnblnk(stdbin)),
     &      'bps -mode text -file ',outs(1:lnblnk(outs)),
     &      ' >/dev/null 2>&1'
        endif
        write(outs,*) 'Testing: ',name(imod)(1:lnblnk(name(imod))),
     &    '  ',cfgitem(imod)(1:lncfg)
        call edisp(iuout,outs)
        if(iverb.eq.2)then
          write(iuout,*) '  '
          write(iuout,*) 'Running: ',doit(1:lnblnk(doit))
        else
          write(iuout,*) '  '
          write(iuout,*) 'Assessing: ',cfgitem(imod)(1:lncfg)
        endif
        call runit(doit,'-')

C After the bps (simulation), the results recovery starts.
C It uses PIF file for save level 4 runs. 

C Invoke the std version of res with PIF file recovery instructions.
        lnres = lnblnk(resitem(imod))
        lnpif = lnblnk(pifitem(imod))
        if(iverb.eq.2)then
          write(doit,'(5a)') stdbin(1:lnblnk(stdbin)),
     &      'res -mode text -file ',resitem(imod)(1:lnres),
     &      ' -act recover silent -actf ',pifitem(imod)(1:lnpif)
        else
          write(doit,'(6a)') stdbin(1:lnblnk(stdbin)),
     &      'res -mode text -file ',resitem(imod)(1:lnres),
     &      ' -act recover silent -actf ',pifitem(imod)(1:lnpif),
     &      ' >/dev/null 2>&1'
        endif
        if(iverb.eq.2)then
          write(iuout,*) '  '
          write(iuout,*) 'Running: ',doit(1:lnblnk(doit))
        endif
        call runit(doit,'-')

C Setup comparison of the output file from the PIF session.
        lnstdrp=lnblnk(stdrep_path)
        lnrep= lnblnk(repitem(imod))
        write(filesA,'(a)') repitem(imod)(1:lnrep)
        write(6,*) 'filesA ',filesB(1:lnblnk(filesA))

C This same .dat file should be found in model test folder.
        lnrep= lnblnk(repitem(imod))
        write(filesB,'(2a)') '../test/',repitem(imod)(1:lnrep)
        write(6,*) 'filesB ',filesB(1:lnblnk(filesB))

C If there is an input.xml file in the folder running bps will likely
C generate an out.xml, out.csv file and out.summary files.
        t72 = resitem(imod)   ! Delete the standard results library.
        write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
        call runit(doit,'-')

C Also clean up flow results if results.mfr found.
        t72 = 'results.mfr'
        write(doit,'(2a)') 'rm -f ',t72(1:lnblnk(t72))
        call runit(doit,'-')

C Now compare the report files from save level 4.  Open fileA and fileB.
C via cmpasci to echo the differences found to channel itu.
        write(iuout,*) '  '
        act = 'v'
        call CMPASCI(itu,iunitA,iunitB,filesA,filesB,act,same,
     &    icount,IER)
        if(same)then
          lnrep= lnblnk(repitem(imod))
          write(outs,'(3a)') 'Comparing ',repitem(imod)(1:lnrep),
     &      ' identical'
          call edisp(iuout,outs)
          lnout=lnblnk(outs)
          write(iuout,*) outs(1:lnout)
          if(itu.ne.iuout)then
            call edisp(itu,outs)
          endif
        else
          lnrep= lnblnk(repitem(imod))
          write(outs,'(3a,i5,a)')'Comparing ',
     &      repitem(imod)(1:lnrep),
     &        ' differs by ',icount,' lines.'
          lnout=lnblnk(outs)
          call edisp(iuout,outs)
          write(iuout,*) outs(1:lnout)
          if(itu.ne.iuout)then
            call edisp(itu,outs)
          endif
          itemadiffs(imod) = itemadiffs(imod) +icount
          adiffs = adiffs + icount
          agdiffs = agdiffs + icount
        endif

        call chdir(odir)
 
C Add overall differences to report.
 200  continue
      write(outs,'(a,i5)') 'ASCII line differences (group)',
     &  int(agdiffs)
      call edisp(iuout,outs)
      write(iuout,*) outs(1:lnblnk(outs))
      if(itu.ne.iuout)then
        call edisp(itu,outs)
      endif
      write(outs,'(a,i5)') 'ASCII line differences (overall)',
     &  int(adiffs)
      call edisp(iuout,outs)
      if(itu.ne.iuout)then
        call edisp(itu,outs)
      endif


C If writing to file close it on last model group.
      if (irpt.eq.1.and.iflg.eq.-1)then
        CALL ERPFREE(ixunit,ISTAT)
      elseif (irpt.eq.-1.and.iflg.eq.-1)then
        CALL ERPFREE(ixunit,ISTAT)
      endif

C If writing to file close it if only running one group.
      if (irpt.eq.1.and.iflg.eq.0)then
        CALL ERPFREE(ixunit,ISTAT)
      elseif (irpt.eq.-1.and.iflg.eq.0)then
        CALL ERPFREE(ixunit,ISTAT)
      endif
      RETURN
      END


      
