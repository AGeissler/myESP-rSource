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


C This file comprises:
C exportcad - model export facility.
C cadin     - conversion of CAD model to ESP-r model.
C rexmpl    - scan the exemplars file and offers a list of models.
C wtarentry - write entry in tar names file summarising ../

C ******************** exportcad ********************
C Controls model export facilities.

      subroutine exportcad(itrc,IER)
#include "building.h"
#include "model.h"
#include "prj3dv.h"

      common/SPAD/MMOD,LIMIT,LIMTTY
      common/OUTIN/IUOUT,IUIN,IEOUT
      integer childterminal  ! picks up mmod from prj invocation
      common/childt/childterminal
      common/FILEP/IFIL
      integer ncomp,ncon
      common/C1/NCOMP,NCON
      common/rpath/path
      common/appw/iappw,iappx,iappy
      integer iappw,iappx,iappy
      common/FOPENED/CFGOK,MLDBOK,MATDBOK,CTLOK,OPTKOK,CFCDBOK
      common/exporttg/xfile,tg,delim
      common/exporttgi/ixopen,ixloc,ixunit
      common/gzonpik/izgfoc,nzg,nznog(mcom)

      integer lnblnk  ! function
      logical OK,CFGOK,MLDBOK,MATDBOK,CTLOK,OPTKOK,concat,CFCDBOK
      logical unixok,found
      integer iappwpc ! application height as % of nominal size

      DIMENSION ITEMS1(13)
      DIMENSION IVALS(MCOM)

      character LFIL*72,path*72,ltmp*72
      character dradcf*72   ! default name
      character ITEMS1*34,doit*248,outs*124,outs248*248
      character tmode*12
      character xfile*144,tg*1,delim*1,ffile*72,dvfil*72
      character longtfile*144,longtfiledos*144,L144*144,lltmp*144
      character brw*8   ! to signal model is being browsed
      character fs*1
      character hold32*32
      logical oke,modcfg
      integer NITMS,INO1 ! max items and current menu item
      integer lnradcf,lnltmp  ! lengths of file names
      character helpinsub*24 ! subroutine name
      character helptopic*24 ! string (unique) for topic
      integer nbhelp     ! number of help lines found

      helpinsub='exportcad'  ! set for subroutine

      if(.NOT.cfgok)then
        call usrmsg('A model must be loaded before',
     &              'it can be exported!','W')
         return
      endif

C Set folder separator (fs) to \ or / as required.
      call isunix(unixok)
      if(unixok)then
        fs = char(47)
      else
        fs = char(92)
      endif
      ioout=2  ! set initial value
      lfil='  '

C Check if xfig exists.
      if(unixok)then
        found_xfig=.false.; hold32='xfig'
        call isinstalled(hold32,found_xfig)
        if(.NOT.found_xfig)then
          call usrmsg('You have not installed xfig so vector',
     &                'plots cannot be viewed.','W')
        endif
      else
        found_xfig=.false.  ! not in Windows
      endif

   33 INO1=-4
      ITEMS1(1)  ='a to VIEWER                '
      ITEMS1(2)  ='b to DXF (V12-V14)         '
      ITEMS1(3)  ='c to RADIANCE              '
      if(found_xfig)then
        ITEMS1(4)='d to xfig (vector ploting) '
      else
        ITEMS1(4)='d to xfig (not installed)  '
      endif
      ITEMS1(5)  ='e to MicroGDS (THINGS file)'
      ITEMS1(6)  ='f to EnergyPlus            '
      ITEMS1(7)  ='g to VRML #2.0             '
      ITEMS1(8)  ='h to ESP-r meta file       '
      ITEMS1(9)  ='i to VEF file              '
      ITEMS1(10) ='! select zones to include  '
      ITEMS1(11) ='* adjust viewing parameters'
      ITEMS1(12) ='? help                     '
      ITEMS1(13) ='- exit menu                '

C Help text for this menu.
      helptopic='cad_export_choices'
      call gethelptext(helpinsub,helptopic,nbhelp)
      NITMS=13
      CALL EMENU('Export model',ITEMS1,NITMS,INO1)

C Prepare filter command lines with path.
      if(INO1.EQ.1.or.INO1.EQ.2.or.INO1.EQ.3.or.INO1.EQ.5.or.
     &   INO1.eq.6.or.INO1.eq.7.or.INO1.eq.8.or.INO1.eq.9.or.
     &   INO1.eq.10.or.INO1.eq.11.or.INO1.eq.12.or.
     &   INO1.eq.13)then
        continue
      else

C Set up initial part of LFIL.
C << Later code might overwrite this and thus the user choice is lost. >>
        iw=0
        if(path(1:2).ne.'./'.and.path(1:2).ne.'  ')then
          call edisp(iuout,' ')
          write(outs,'(A,A)') 'The current path is: ',path
          call edisp(iuout,outs)
          CALL EASKAB('The model is in a remote folder.',
     &    'Export to:','current folder','home folder',
     &    IW,nbhelp)
        endif
        if(iw.eq.1)then
          write(LFIL,'(a)') path(1:lnblnk(path))
        else
          LFIL='./'     ! environment call required here
        endif
      endif

C Construct a 'VIEWER' format file. First compose file name taking
C into account pwdinitial and pwdtocfg.
 11   IF(INO1.EQ.1)THEN
        if(pwdtocfg(1:1).eq.'!')then
          write(L144,'(4a)')path(1:lnblnk(path)),fs,
     &      cfgroot(1:lnblnk(cfgroot)),'.vew'
        else
          write(L144,'(6a)') pwdinitial(1:lnpwdi),fs,
     &      pwdtocfg(1:lnblnk(pwdtocfg)),fs,
     &      cfgroot(1:lnblnk(cfgroot)),'.vew'
        endif
        call tstamp('>','PRJ: export to viewer')
        dvfil='std.vew'  ! default file name
        doit = ' '

C The X11 version will return only the name of the
C file, while the GTK version will return the
C name and full path. L144 is the suggested name and
C lltmp is the file name returned by the user selection.
        lltmp='  '
        call edisp(iuout,' ')
        CALL EASKXORGTKF(L144,
     &    'Viewer export file name?',' ',
     &    dvfil,lltmp,'viewer file',IER,nbhelp)

C On user request, jump back and re-display the menu.
        if(ier.eq.-3)then
          return  ! cancel detected, return.
        endif
        if(ier.ne.0.or.lltmp(1:2).eq.'  ')then
          goto 11
        endif

C Ensure full explicit path for unix or non-unix.
        if(unixok)then
          call addpath(LCFGF,longtfile,concat)  ! for the model cfg file
        else

C If running on a non-unix machine see if there are spaces in the name
C and change any / to \.
          call addpath(LCFGF,longtfile,concat)  ! for the model cfg file
          call cmdfiledos(longtfile,longtfiledos,ier)
          longtfile=' '
          longtfile=longtfiledos
        endif
        if(itrc.le.1)then
          write(doit,'(4a)') 'ecnv -obs -if esp -in ',
     &      longtfile(1:lnblnk(longtfile)),' -of viewer -out ',
     &      lltmp(1:lnblnk(lltmp))
          call usrmsg('Starting conversion via',doit,'-')
          call runit(doit,'-')
        elseif(itrc.gt.1)then
          write(doit,'(4a)') 'ecnv -v -obs -if esp -in ',
     &      longtfile(1:lnblnk(longtfile)),' -of viewer -out ',
     &      lltmp(1:lnblnk(lltmp))
          call usrmsg('Starting conversion via',doit,'-')
        endif
        tmode = 'text'
        call runit(doit,tmode)
        write(outs248,'(3a)')'The file ',lltmp(1:lnblnk(lltmp)),
     &    ' has been created.'
        call edisp248(iuout,outs248,90)
        goto 33  ! display the menu again

C DXF creation. Process similar to ecnv filter for VIEWER.
      ELSEIF(INO1.EQ.2)THEN
        if(pwdtocfg(1:1).eq.'!')then
          write(L144,'(4a)')path(1:lnblnk(path)),fs,
     &      cfgroot(1:lnblnk(cfgroot)),'.DXF'
        else
          write(L144,'(6a)') pwdinitial(1:lnpwdi),fs,
     &      pwdtocfg(1:lnblnk(pwdtocfg)),fs,
     &      cfgroot(1:lnblnk(cfgroot)),'.DXF'
        endif
        call tstamp('>','PRJ: export to DXF')
        dvfil='std.DXF'  ! default file name
        doit = ' '

C L144 is the suggested name and lltmp is the file name returned
C by the user selection.
        lltmp='  '
        call edisp(iuout,' ')
        CALL EASKXORGTKF(L144,'DXF export file name?',' ',
     &    dvfil,lltmp,'dxf file',IER,nbhelp)

C On user request, jump back and re-display the menu.
        if(ier.eq.-3)then
          return  ! cancel detected, return.
        endif
        if(ier.ne.0.or.lltmp(1:2).eq.'  ')then
          goto 11
        endif
        call edisp(iuout,'Performing model conversion.')
        if(unixok)then
          call addpath(LCFGF,longtfile,concat)
        else

C If running on a non-unix machine see if there are spaces in the name
C and change any / to \.
          call addpath(LCFGF,longtfile,concat)
          call cmdfiledos(longtfile,longtfiledos,ier)
          longtfile=' '
          longtfile=longtfiledos
        endif
        if(itrc.le.1)then
          write(doit,'(4a)') 'ecnv -obs -if esp -in ',
     &       longtfile(1:lnblnk(longtfile)),
     &      ' -of dxf -out ',lltmp(1:lnblnk(lltmp))
        elseif(itrc.gt.1)then
          write(doit,'(4a)') 'ecnv -v -obs -if esp -in ',
     &       longtfile(1:lnblnk(longtfile)),
     &      ' -of dxf -out ',lltmp(1:lnblnk(lltmp))
        endif
        call usrmsg('Starting conversion via',doit,'-')
        tmode = 'text'
        call runit(doit,tmode)
        write(outs248,'(3a)')'The file ',lltmp(1:lnblnk(lltmp)),
     &    ' has been created.'
        call edisp248(iuout,outs248,90)
        goto 33

C Check if Radiance has been installed. If not, do not start e2r.
      ELSEIF(INO1.EQ.3)THEN
        if(unixok)then
          found_ximage=.false.; hold32='ximage'
          call isinstalled(hold32,found_ximage)
          if(.NOT.found_ximage)then
            call usrmsg(
     &      'Radiance noy installed so images',
     &      'cannot be generated.','W')
            goto 33
          endif
        else
          found_ximage=.false.  ! not in Windows
        endif

C Check if model path for radiance has been set, otherwise create a
C ../rad folder. If there is no Radiance Scene file make it based
C on the model root name. Take into account where prj was started.
        modcfg=.false.
        if(pwdtocfg(1:1).eq.'!')then
          if(radpth(1:3).eq.'./ '.or.radpth(1:3).eq.'.\\ ')then
            CALL PHELPD('visu-rad-warning',7,'-',0,0,IER)
            write(radpth,'(3a)')'..',fs,'rad'
            write(doit,'(4a)') 'mkdir ',
     &        cfgroot(1:lnblnk(cfgroot)),fs,'rad'
            call usrmsg('Creating folder for radiance model:',doit,'P')
            call runit(doit,'-')

C And update radpath. The L144 is for future use of EASKXORGTKF.
            write(radpth,'(3a)') '..',fs,'rad'
            modcfg=.true.
          endif

          if(lradcf(1:7).eq.'UNKNOWN')then  ! setup in form of ../rad/xx
            write(lradcf,'(6a)') '..',fs,'rad',fs,
     &        cfgroot(1:lnblnk(cfgroot)),'.rcf'
            modcfg=.true.
          endif
        else

C Prj has not been invoked from within the model cfg folder.
          if(radpth(1:3).eq.'./ '.or.radpth(1:3).eq.'.\\ ')then
            CALL PHELPD('visu-rad-warning',7,'-',0,0,IER)
            write(radpth,'(3a)')'..',fs,'rad'
            write(doit,'(8a)') 'mkdir ',pwdinitial(1:lnpwdi),fs,
     &        pwdtocfg(1:lnblnk(pwdtocfg)),fs,'..',fs,'rad'
            call usrmsg('Creating folder for radiance model:',doit,'P')
            call runit(doit,'-')

C And update radpath. The L144 is for future use of EASKXORGTKF.
            write(radpth,'(3a)') '..',fs,'rad'
            modcfg=.true.
          endif

          if(lradcf(1:7).eq.'UNKNOWN')then
            write(lradcf,'(2a)') cfgroot(1:lnblnk(cfgroot)),'.rcf'
            modcfg=.true.
          endif
        endif

C Check rcf name and contents of file on return from e2r.
 319    ltmp=lradcf
        dradcf='scene.rcf'
        CALL EASKS(ltmp,'Radiance scene file name?',
     &    '  ',72,dradcf,'Radiance scene file name',IER,nbhelp)
        if(ltmp.eq.' ')goto 319 

C rcf name OK and differs from prior name, update cfg file and start e2r.
        lnradcf=lnblnk(lradcf)
        lnltmp=lnblnk(ltmp)
        if(ltmp(1:lnltmp).eq.lradcf(1:lnradcf))then
          continue
        else
          modcfg=.true.
          lradcf = ltmp
        endif

C Start e2r with current configuration file.
        CALL EASKOK(' ','Run Radiance controller?',OK,nbhelp)
        if(OK)then

C Save cfg file including ../rad folder and lradcf file name and start e2r.
          if(modcfg)then
            call tstamp('>','PRJ: save configuration (with rcf file)')
            CALL EMKCFG('-',IER)
          endif

C Set browsing to no.
          brw = ' -b no '
          call tstamp('>','PRJ: export to RADIANCE')
          doit = ' '
          call terminalmode(childterminal,tmode)
          if(unixok)then
            call addpath(LCFGF,longtfile,concat)
          else
            call addpath(LCFGF,longtfile,concat)
            call cmdfiledos(longtfile,longtfiledos,ier)
            longtfile=' '
            longtfile=longtfiledos
          endif
          if(iappw.eq.690)then
            iappwpc=100
          else
            iappwpc=nint(100.0*(real(iappw)/690.0))  ! reconstitute %
          endif
          if(iappwpc.gt.0.and.iappwpc.le.200)then
            write(doit,'(a,3i4,4a)') 'e2r -s ',
     &        iappwpc,iappx+35,iappy+40,' -file ',
     &        longtfile(1:lnblnk(longtfile)),brw,' & '
          else
            write(doit,'(4a)') 'e2r -s 0 0 0 -file ',
     &        longtfile(1:lnblnk(longtfile)),brw,' & '
          endif
          call usrmsg('Starting Radiance conversion: ',doit,'-')
          call runit(doit,'-')
        endif

      ELSEIF(INO1.EQ.4)THEN

C Toggle output redirection to xfig.
        call easkok(' ','Accept wireframe composition?',oke,nbhelp)
        if(.NOT.oke)then
          call usrmsg('Adjust wireframe image and',
     &      'then re-select this option.','W')
          goto 33
        endif

        if (MMOD.ne.8) then
          call edisp(iuout,' ')
          call edisp(iuout,'Drawing to virtual screen.')
        endif

C Open the export file and then request that the wireframe
C or virtual wireframe is redrawn (in the C code the low level calls
C will be writing out the vector information to the export file.
        write(xfile,'(2a)') cfgroot(1:lnblnk(cfgroot)),'_wire'
        call ctlexp(xfile,ixopen,ixloc,ixunit,'G','wire frame',IER)

C If user canceled the wireframe export loop back.
        if(ier.eq.-3)then
          goto 33
        endif
        call usrmsg('Sending wireframe ...',' ','-')
        MODIFYVIEW=.TRUE.
        CALL redraw(IER)

C Toggle the wirefram export file closed.
        call usrmsg('Sending wireframe ... done.',' ','-')
        call ctlexp(xfile,ixopen,ixloc,ixunit,'G','Wire frame',IER)
        MODIFYVIEW=.TRUE.

C Check if conversion to xfig wanted.
        call easkok(' ','Convert output to xfig format?',ok,nbhelp)
        if(ok)then
          write(ffile,'(a,a4)') xfile(1:lnblnk(xfile)),'.fig'
          doit = ' '
          write(doit,'(4a)') 'ecnv -if ww -in ',
     &      xfile(1:lnblnk(xfile)),' -of xfig -out ',
     &      ffile(1:lnblnk(ffile))
          call usrmsg('Starting xfig conversion: ',doit,'-')
          call runit(doit,'-')

          iuj=ifil+1  ! remove the intermediate _wire file
          close(iuj)
          CALL EFOPSEQ(iuj,xfile,1,IER)
          CALL EFDELET(iuj,ISTAT)

          if(found_xfig)then
            call easkok(' ','Invoke xfig?',ok,nbhelp)
            if(ok)then
              call edisp(iuout,' ')
              call edisp(iuout,'Optional: use fig2dev to create a png')
              call edisp(iuout,'image via:')
              call edisp(iuout,'fig2dev -L png -m 1.5 -S 4 -b 10 ')
              call edisp(iuout,'  xx.fig > xx.png')
              doit = ' '
              write(doit,'(3a)') 'xfig ',ffile(1:lnblnk(ffile)),' & '
              call usrmsg('Starting xfig: ',doit,'-')
              call runit(doit,'-')
            endif
          endif
        endif

C Export to THINGS (MicroGDS THF) file. Treat obstructions as clusters.
      ELSEIF(INO1.EQ.5)THEN
        if(pwdtocfg(1:1).eq.'!')then
          write(L144,'(4a)')path(1:lnblnk(path)),fs,
     &      cfgroot(1:lnblnk(cfgroot)),'.THF'
        else
          write(L144,'(6a)') pwdinitial(1:lnpwdi),fs,
     &      pwdtocfg(1:lnblnk(pwdtocfg)),fs,
     &      cfgroot(1:lnblnk(cfgroot)),'.THF'
        endif
        call tstamp('>','PRJ: export to DXF')
        dvfil='std.THF'  ! default file name
        doit = ' '

C L144 is the suggested name and lltmp is the file name returned
C by the user selection.
        lltmp='  '
        call edisp(iuout,' ')
        CALL EASKXORGTKF(L144,'MicroGDS THINGS file name?',' ',
     &    dvfil,lltmp,'THF file',IER,nbhelp)

C If user request jump back and re-display the menu.
        if(ier.eq.-3)then
          return  ! cancel detected, return.
        endif
        if(ier.ne.0.or.lltmp(1:2).eq.'  ')then
          goto 11
        endif
        call tstamp('>','PRJ: export to THINGS')
        doit = ' '
        call edisp(iuout,'Performing model conversion.')

C << Why is the addpath not used? >>
        if(unixok)then
          call addpath(LCFGF,longtfile,concat)
        else
          call addpath(LCFGF,longtfile,concat)
          call cmdfiledos(longtfile,longtfiledos,ier)
          longtfile=' '
          longtfile=longtfiledos
        endif
        if(itrc.le.1)then
          write(doit,'(4a)') 'ecnv -obs -if esp -in ',
     &       longtfile(1:lnblnk(longtfile)),
     &      ' -of thf -out ',lltmp(:lnblnk(lltmp))
        elseif(itrc.gt.1)then
          write(doit,'(4a)') 'ecnv -v -obs -if esp -in ',
     &       longtfile(1:lnblnk(longtfile)),
     &      ' -of thf -out ',lltmp(1:lnblnk(lltmp))
        endif
        call usrmsg('Starting THF conversion: ',doit,'-')
        tmode = 'text'
        call runit(doit,tmode)
        write(outs248,'(3a)')'The file ',lltmp(1:lnblnk(lltmp)),
     &    ' has been created.'
        call edisp248(iuout,outs248,90)
        goto 33

C Export to EnergyPlus IDF file (work in progress).
      ELSEIF(INO1.EQ.6)THEN
        if(pwdtocfg(1:1).eq.'!')then
          write(L144,'(4a)')path(1:lnblnk(path)),fs,
     &      cfgroot(1:lnblnk(cfgroot)),'.idf'
        else
          write(L144,'(6a)') pwdinitial(1:lnpwdi),fs,
     &      pwdtocfg(1:lnblnk(pwdtocfg)),fs,
     &      cfgroot(1:lnblnk(cfgroot)),'.idf'
        endif
        call tstamp('>','PRJ: export to energyplus')
        dvfil='std.idf'  ! default file name
        doit = ' '

C L144 is the suggested name and lltmp is the file name returned
C by the user selection.
        lltmp='  '
        call edisp(iuout,' ')
        CALL EASKXORGTKF(L144,'EnergyPlus input file name?',' ',
     &    dvfil,lltmp,'E+ idf file',IER,nbhelp)

C If user request jump back and re-display the menu.
        if(ier.eq.-3)then
          return  ! cancel detected, return.
        endif
        if(ier.ne.0.or.lltmp(1:2).eq.'  ')then
          goto 11
        endif

C Gather information about IDF file version.
        ever=7.2
        CALL EASKATOG('EnergyPlus IDF file version',' ',
     &    '7.2','8.0','8.5','8.8','9.0','cancel','-',IW,nbhelp)
        if(IW.eq.1)then
          ever=7.2
        elseif(IW.eq.2)then
          ever=8.0
        elseif(IW.eq.3)then
          ever=8.5
        elseif(IW.eq.4)then
          ever=8.8
        elseif(IW.eq.5)then
          ever=9.0
        elseif(IW.eq.6)then
          continue
        endif
        if(IW.eq.6)then
          continue
        else
          call edisp(iuout,'Performing model conversion.')
          itrcc=2
          ioout=2
          open(ioout,file=lltmp,status='UNKNOWN',err=901)
          call e2eplus(itrcc,ioout,ever,ier)
          close(ioout)
          write(outs248,'(3a)')'The file ',lltmp(1:lnblnk(lltmp)),
     &      ' has been created.'
          call edisp248(iuout,outs248,90)
        endif
        goto 33

C Export to VRML (2.0) file.
      ELSEIF(INO1.EQ.7)THEN
        if(pwdtocfg(1:1).eq.'!')then
          write(L144,'(4a)')path(1:lnblnk(path)),fs,
     &      cfgroot(1:lnblnk(cfgroot)),'.vrml'
        else
          write(L144,'(6a)') pwdinitial(1:lnpwdi),fs,
     &      pwdtocfg(1:lnblnk(pwdtocfg)),fs,
     &      cfgroot(1:lnblnk(cfgroot)),'.vrml'
        endif
        call tstamp('>','PRJ: export to VRML')
        dvfil='std.vrml'  ! default file name
        doit = ' '

C The X11 version will return only the name of the
C file, while the GTK version will return the
C name and full path. L144 is the suggested name and
C lltmp is the file name returned by the user selection.
        lltmp='  '
        call edisp(iuout,' ')
        CALL EASKXORGTKF(L144,'VRML file name?',' ',
     &    dvfil,lltmp,'vrml file',IER,nbhelp)

C On user request, jump back and re-display the menu.
        if(ier.eq.-3)then
          return  ! cancel detected, return.
        endif
        if(ier.ne.0.or.lltmp(1:2).eq.'  ')then
          goto 11
        endif

        if(unixok)then
          call addpath(LCFGF,longtfile,concat)
        else

C If running on a non-unix machine see if there are spaces in the name
C and change any / to \.
          call addpath(LCFGF,longtfile,concat)
          call cmdfiledos(longtfile,longtfiledos,ier)
          longtfile=' '
          longtfile=longtfiledos
        endif
        if(itrc.le.1)then
          write(doit,'(4a)') 'ecnv -if esp -in ',
     &      longtfile(1:lnblnk(longtfile)),
     &      ' -of vrml -out ',lltmp(1:lnblnk(lltmp))
        elseif(itrc.gt.1)then
          write(doit,'(4a)') 'ecnv -v -if esp -in ',
     &      longtfile(1:lnblnk(longtfile)),
     &      ' -of vrml -out ',lltmp(1:lnblnk(lltmp))
        endif
C        call usrmsg('Starting VRML conversion: ',doit,'-')
        tmode = 'text'
        call runit(doit,tmode)
        write(outs248,'(3a)')'The file ',lltmp(1:lnblnk(lltmp)),
     &    ' has been created.'
        call edisp248(iuout,outs248,90)
        call edisp(iuout,'Use a VRML viewer to explore it.')
        goto 33  ! display the menu again

C Export to ESP-r meta file.
      ELSEIF(INO1.EQ.8)THEN
        if(pwdtocfg(1:1).eq.'!')then
          write(L144,'(4a)')path(1:lnblnk(path)),fs,
     &      cfgroot(1:lnblnk(cfgroot)),'.mesp'
        else
          write(L144,'(6a)') pwdinitial(1:lnpwdi),fs,
     &      pwdtocfg(1:lnblnk(pwdtocfg)),fs,
     &      cfgroot(1:lnblnk(cfgroot)),'.mesp'
        endif

        call tstamp('>','PRJ: export to ESP-r meta')
        dvfil='std.mesp'  ! default file name
        doit = ' '

C L144 is the suggested name and lltmp is the file name returned
C by the user selection.
        lltmp='  '
        call edisp(iuout,'  ')  ! echo a blank line
        CALL EASKXORGTKF(L144,' ','ESP-r meta file name?',
     &    dvfil,lltmp,'ESP-r Meta file',IER,nbhelp)

C On user request, jump back and re-display the menu.
        if(ier.eq.-3)then
          return  ! cancel detected, return.
        endif
        if(ier.ne.0.or.lltmp(1:2).eq.'  ')then
          goto 11
        endif
        ioout=2
        call METAWRITE(ioout,lltmp,'-',IER)
        write(outs248,'(3a)')'The file ',lltmp(1:lnblnk(lltmp)),
     &    ' has been created.'
        call edisp248(iuout,outs248,90)
        goto 33

C Export to acrobat VEF file.
      ELSEIF(INO1.EQ.9)THEN
        if(pwdtocfg(1:1).eq.'!')then
          write(L144,'(4a)')path(1:lnblnk(path)),fs,
     &      cfgroot(1:lnblnk(cfgroot)),'.vef'
        else
          write(L144,'(6a)') pwdinitial(1:lnpwdi),fs,
     &      pwdtocfg(1:lnblnk(pwdtocfg)),fs,
     &      cfgroot(1:lnblnk(cfgroot)),'.vef'
        endif
        call tstamp('>','PRJ: export to VEF')
        dvfil='std.vef'  ! default file name
        doit = ' '

C L144 is the suggested name and lltmp is the file name returned
C by the user selection.
        lltmp='  '
        call edisp(iuout,'  ')  ! echo a blank line
        CALL EASKXORGTKF(L144,'acrobat VEF file name?',' ',
     &    dvfil,lltmp,'acrobat vef file',IER,nbhelp)

C On user request, jump back and re-display the menu.
        if(ier.eq.-3)then
          return  ! cancel detected, return.
        endif
        if(ier.ne.0.or.lltmp(1:2).eq.'  ')then
          goto 11
        endif

        ioout=2
        CALL EASKABC('Specify GEN bodies or separate polygons with',
     &    'glazing clear:','GEN bodies','with transparent windows',
     &    'cancel',IRT,nbhelp)
        if(IRT.eq.3)then
          continue
        else
          CALL EASKABC(' ','Colour options:','grey zones',
     &                          'coloured zones','cancel',IRT2,nbhelp)
          if (IRT2.eq.3) then
            continue
          else
            IRT2=IRT2-1
            CALL EFOPSEQ(ioout,lltmp,3,IER)
            if(IRT.eq.1)then
              call e2vef(itrc,ioout,'g',IRT2)  ! see end of this file
            elseif(IRT.eq.2)then
              call e2vef(itrc,ioout,'p',IRT2)
            endif
            CALL ERPFREE(ioout,ios)
            write(outs248,'(3a)')'The file ',lltmp(1:lnblnk(lltmp)),
     &      ' has been created.'
            call edisp248(iuout,outs248,90)
          endif
        endif
        goto 33

C Zones selection.
      ELSEIF(INO1.EQ.10)THEN
        INPIC=NCOMP
        call ASKMULTIZON(INPIC,IVALS,'Zones to include:',
     &    'zone list','-',IER) 
        nzg = inpic
        if(izgfoc.ne.0)then
          found=.false.
          do mz=1,nzg
            nznog(mz)=IVALS(mz)
            if(ivals(mz).eq.izgfoc)found=.true.
          enddo
          if(.NOT.found)then
            nzg=nzg+1
            nznog(nzg)=izgfoc
          endif
        elseif(izgfoc.eq.0)then
          do mz=1,nzg
            nznog(mz)=IVALS(mz)
          enddo
        endif
        MODBND=.TRUE.
        MODLEN=.TRUE.
        MODIFYVIEW=.TRUE.

C Viewing parameters.
      ELSEIF(INO1.EQ.11)THEN
        call evset('-',IER)                          ! adjust viewing parameters
      ELSEIF(INO1.EQ.12)THEN
       CALL PHELPD('export menu',nbhelp,'-',0,0,IER) ! pop up help text
      ELSEIF(INO1.EQ.13)THEN
        RETURN
      else
        goto 33
      endif
      goto 33

C Error.
 901  call isunix(unixok)
      if(unixok)write(6,*) 'Unable to open IDF file ',LFIL,' on ',ioout
      return
     
      END

C ******************** cadin ********************
C Converts a CAD input model to an ESP-r model.

      subroutine cadin(itrc,ier)
#include "building.h"
#include "model.h"
      
      integer lnblnk  ! function definition

      common/FILEP/IFIL
      INTEGER :: ifil
      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      INTEGER :: iuout,iuin,ieout
      common/warncounter/iecount
      character iefile*72
      common/logs/iefile
      
      integer ncomp,ncon
      common/C1/NCOMP,NCON
      common/rcmd/LCMDFL
      common/user/browse

C External cad package.
      common/cad3rd/cadlbl,cadexe,cadfmt

C The model root and mpath passed to silentmodel call.
      character metaroot*32,metampath*72
      common/cadmeta/metaroot,metampath

      LOGICAL OK,concat,browse,confirm,ckpath,unixok
      logical itisanexemplar  ! passed to newprb to ensure correct browse mode

      character LFIL*72,LCMDFL*144
      character outs*124,OUTSTR*124,doit*248
      character cadlbl*20,cadexe*20,cadfmt*5
      character longtfile*144,longtfiledos*144
      character placemodelin*96
      character exemplar*144 ! full path to new cfg file
      character fs*1
      character thecfgis*72  ! cfg file only
      character cfgpath*144  ! full path to cfg folder for refocus
      character theprimedirective*8
      logical XST
      character hold32*32,msg*72,msg2*30
      logical found_cad      ! to test application existence

C For help messages
      character helpinsub*24 ! subroutine name
      character helptopic*24 ! string (unique) for topic
      integer nbhelp         ! number of help lines found

      helpinsub='cadin'      ! set for subroutine

C Check if running Unix-like OS.
      placemodelin='  '
      call isunix(unixok)
      if(unixok)then
        fs = char(47)
      else
        fs = char(92)
      endif

C CAD definition and importing.
      helptopic='cad_import_choices'
      call gethelptext(helpinsub,helptopic,nbhelp)

C Tell the user to create a model configuration file.
      CALL PHELPD('CAD warning',nbhelp,'-',0,0,IER)

C Determine which format of CAD file will be scanned in.
      write(msg2,'(2a)') 'invoke ',cadexe(1:lnblnk(cadexe))
      CALL EASKATOG(' ','Options:',
     &  'parse DXF (V12-14)','parse gbXML',msg2,
     &  'parse ESP-r META','cancel',' ',' ',IRT,nbhelp)
      if(irt.eq.1)then
        cadfmt='DXF '
      elseif(irt.eq.2)then
        cadfmt='gbXML'
      elseif(irt.eq.3)then

C Check if CAD tool has been installed. If not, do not start.
        if(unixok)then
          found_cad=.false. 
          write(hold32,'(a)') cadexe(1:lnblnk(cadexe))
          call isinstalled(hold32,found_cad)
          if(.NOT.found_cad)then
            write(msg,'(3a)') 'The CAD tool ',
     &        cadexe(1:lnblnk(cadexe)),' is not installed.'
            call usrmsg(msg,'Exiting.','W')
            return
          endif
        else
          found_cad=.false.  ! not in Windows
        endif

C Create one or more dummy zones for gtool to work with.
        call edisp(iuout,'CAD tool requires dummy zones to exist')
        call edisp(iuout,'prior to starting. Specify how many to')
        call edisp(iuout,'create. Note that zone names must follow the')
        call edisp(iuout,'same root names as the ESP-r model.')
        call create_dummy_zones()
        write(outs,'(3a)') 'Proceed with ',
     &    cadexe(1:lnblnk(cadexe)),'?'

        call usrmsg(
     &    'After the CAD tool starts, Project Manager will exit.',
     &    'Please restart it when the CAD session is finished.','W')
        CALL EASKOK(' ',outs,OK,nbhelp)
        IF(OK)then
          call pausems(1000)
          write(doit,'(2a)') cadexe(1:lnblnk(cadexe)),' & '
          call runit(doit,'-')
          call pausems(1000)
          call usrmsg(
     &    'Shutting down Project Manager.','...','-')
          call pausems(1000)

C Clear arrays.
          CALL DeallocateAllArrays
          close(ieout)
          if(iecount.eq.0)then     ! if there were no warnings delete
            open(ieout,file=iefile,status='UNKNOWN')
            call efdelet(ieout,istat)
          endif
          CALL EPAGEND
          STOP
        endif
        return
      elseif(irt.eq.4)then

C Attempt to scan ESP-r META file.
        LFIL='  '
  29    CALL EASKS(LFIL,' ','ESP-r META file to import?',
     &    72,' ','META file',IER,nbhelp)

C Read the first couple of lines of the import file to see if it
C is a DXF file.
        if(LFIL(1:2).ne.'  '.and.LFIL(1:4).ne.'UNKN')then
          INQUIRE (FILE=LFIL,EXIST=XST)
          if(.NOT.XST)then
            write(outs,'(a,a)')'Could not locate ',LFIL(1:lnblnk(LFIL))
            call easkok(outs,'Retry?',OK,nbhelp)
            if(OK)then
              goto 29
            else
              return
            endif
          endif
          call clrprb  ! clear the model common blocks prior to scan
          write(longtfile,'(a)') LFIL(1:lnblnk(LFIL))
          IUNIT=IFIL+1
          modeltitle='UNKNOWN'

C If iverb passed as non-zero then change directive.
          call silentread(iunit,longtfile,'-',theprimedirective,ier)
          call edisp(iuout,theprimedirective)
          call edisp(iuout,modeltitle)
          call edisp248(iuout,modeldocblock,80)

C Figure out where the new model cfg is going to be placed.
C          call usrdir(pwdinitial)
          write(6,*) 'the pwd command is ',
     &      pwdinitial(1:lnblnk(pwdinitial))
          thecfgis='  '
          write(thecfgis,'(2a)') metaroot(1:lnblnk(metaroot)),'.cfg'
          cfgpath='  '
          write(cfgpath,'(5a)') pwdinitial(1:lnblnk(pwdinitial)),fs,
     &      metampath(1:lnblnk(metampath)),fs,'cfg'
          call edisp(iuout,'  ')
          call edisp(iuout,'The model cfg folder is:')
          call edisp(iuout,cfgpath)
          call edisp(iuout,'The cfg file is:')
          call edisp(iuout,thecfgis)
          call edisp(iuout,'It is now OK to [YES RESTART].')

C Suggest to user to refocus prj within the model cfg folder at
C this specific cfg file.
          call refocus(cfgpath,thecfgis,iret)
          if(iret.eq.0)then
            continue
          elseif(iret.eq.-2)then
            call edisp(iuout,
     &      'Features using Radiance calls may not work as expected.')
          endif
          return
        endif
      elseif(irt.eq.5)then
        return
      endif

C If there is a currently loaded model then convert it.
      if(NCOMP.gt.0.and.cadfmt(1:3).eq.'DXF')then
        call tstamp('>','prj: start CAD')
        doit = ' '
        CALL EASKS(LFIL,'Putting current model into DXF file named:',
     &    ' ',72,' ','dxf export file',IER,nbhelp)
        call edisp(iuout,' Performing model conversion.')
        if(unixok)then
          call addpath(LCFGF,longtfile,concat)
        else

C If running on a non-unix machine see if there are spaces in the name
C and change any / to \.
          call addpath(LCFGF,longtfile,concat)
          call cmdfiledos(longtfile,longtfiledos,ier)
          longtfile=' '
          longtfile=longtfiledos
        endif
        if(itrc.le.1)then
          write(doit,'(4a)') 'ecnv -obs -if esp -in ',
     &      longtfile(1:lnblnk(longtfile)),
     &      ' -of dxf -out ',LFIL(1:lnblnk(LFIL))
        elseif(itrc.gt.1)then
          write(doit,'(4a)') 'ecnv -v -obs -if esp -in ',
     &      longtfile(1:lnblnk(longtfile)),
     &      ' -of dxf -out ',LFIL(1:lnblnk(LFIL))
        endif

C Invoke export of current ESP-r model to DXF.
        call runit(doit,'-')
        write(outs,'(3a)')' The file ',LFIL(1:lnblnk(LFIL)),
     &    ' has been created.'
        call edisp(iuout,outs)
        write(outs,'(3a)') ' Proceed with ',
     &    cadlbl(1:lnblnk(cadlbl)),' ?'
        CALL EASKOK(' ',outs,OK,nbhelp)
        IF(OK)then
          call runit(cadexe,'graph')
        endif
      endif

 28   if(cadfmt(1:3).eq.'DXF')then
        CALL EASKOK(' ','Import DXF file?',OK,nbhelp)
      elseif(cadfmt(1:5).eq.'gbXML')then
        OK=.true.
      endif
      IF(OK)then
        CALL EASKS(LFIL,'File to import?',
     &    ' ',72,' ','cad import file',IER,nbhelp)

C Read the first couple of lines of the import file to see if it
C is a DXF file.
        CALL ERPFREE(ifun,ISTAT)
        call FPOPEN(ifun,ISTAT,1,0,LFIL)
        IF(ISTAT.LT.0)THEN
          write(outs,'(a,a)')'Could not locate ',LFIL(1:lnblnk(LFIL))
          call easkok(outs,'Retry?',OK,nbhelp)
          if(OK)then
            goto 28
          else
            return
          endif
        endif
        CALL STRIPC(ifun,OUTSTR,0,ND,1,'CAD first line',IER)
        if(OUTSTR(1:13).eq.'<?xml version')then
          call edisp(iuout,'Based on 1st line of file:')
          call edisp(iuout,outstr)
          call edisp(iuout,'this file is assumed to be a gbXML file.')
          call edisp(iuout,' ')
          call edisp(iuout,'An example of [full path to model folder]:')
          call edisp(iuout,'/home/fred/Models/office')
          cadfmt='gbXML'
          CALL EASKS(placemodelin,'Folder for model (give full path):',
     &      ' ',96,' ','model folder',IER,nbhelp)

C Setup various file names needed.
          exemplar='  '
          write(exemplar,'(5a)') placemodelin,fs,'cfg',fs,'bldg-1.cfg'
          cfgpath='  '
          write(cfgpath,'(3a)')  placemodelin,fs,'cfg'
          thecfgis='  '
          write(thecfgis,'(a)') 'bldg-1.cfg'

        else
          CALL STRIPC(ifun,OUTSTR,0,ND,1,'CAD 2nd line',IER)
          if(OUTSTR(1:7).eq.'SECTION')then
            call edisp(iuout,'Based on 2nd line of file:')
            call edisp(iuout,outstr)
            call edisp(iuout,'this file is assumed to be a DXF file.')
            cadfmt='DXF '
          else
            call edisp(iuout,
     &      'Initial lines of file do not indicate the file is either')
            call edisp(iuout,
     &      'a DXF file or an ASCII gbXML file.')
            cadfmt='UNKN'
          endif
        endif
        CALL ERPFREE(ifun,ISTAT)

        if(cadfmt(1:3).eq.'DXF')then
          if(unixok)then
            call addpath(LCFGF,longtfile,concat)
          else

C If running on a non-unix machine see if there are spaces in the name
C and change any / to \.
            call addpath(LCFGF,longtfile,concat)
            call cmdfiledos(longtfile,longtfiledos,ier)
            longtfile=' '
            longtfile=longtfiledos
          endif
          call edisp(iuout,' Performing model conversion.')
          call tstamp('>','PRJ: import from CAD')
          if(itrc.le.1)then
            write(doit,'(4a)') 'ecnv -obs -mm -if dxf -in ',
     &        LFIL(1:lnblnk(LFIL)),' -of esp -u upgrade -out ',
     &        longtfile(1:lnblnk(longtfile))
          elseif(itrc.gt.1)then
            write(doit,'(4a)')'ecnv -v -obs -mm -if dxf -in ',
     &        LFIL(1:lnblnk(LFIL)),' -of esp -u upgrade -out ',
     &        longtfile(1:lnblnk(longtfile))
          endif
          call usrmsg('doing conversion via',doit,'-')
          call runit(doit,'-')

C Re-read the upgraded configuration file. Reconstruct the path.
          itisanexemplar=.false.
          browse=.false.
          confirm=.true.
          ckpath=.false.
          LCMDFL=longtfile
          call NEWPRB(ITRC,confirm,ckpath,itisanexemplar,IER)
        elseif(cadfmt(1:5).eq.'gbXML')then
          call edisp(iuout,' Performing model conversion.')
          call tstamp('>','PRJ: parsing gbXML')
          if(itrc.le.1)then
            write(doit,'(4a)') 'ecnv -if gbxml -in ',
     &        LFIL(1:lnblnk(LFIL)),' -of esp -out ',
     &        placemodelin(1:lnblnk(placemodelin))
          elseif(itrc.gt.1)then
            write(doit,'(4a)')'ecnv -v -if gbxml -in ',
     &        LFIL(1:lnblnk(LFIL)),' -of esp -out ',
     &        placemodelin(1:lnblnk(placemodelin))
          endif
          write(6,*) doit
          call usrmsg('doing conversion via',doit,'-')
          call pausems(1000)
          call runit(doit,'-')
          call pausems(1000)

C Use call to refocus to start a new prj within the model cfg folder.
          call refocus(cfgpath,thecfgis,iret)
          if(iret.eq.0)then
            continue
          elseif(iret.eq.-2)then
            call edisp(iuout,
     &      'Features using Radiance calls may not work as expected.')
          endif
        endif
      endif
      return
      end


C ******** gtoolin
C gtoolin sets up files and invokes the CAD tool gtool.
      subroutine gtoolin(itrc,ier)
#include "building.h"
#include "model.h"
      
      integer lnblnk  ! function definition

      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      INTEGER :: iuout,iuin,IEOUT
      common/warncounter/iecount

C Session log file (for errors etc.)
      character iefile*72
      common/logs/iefile

C External cad package.
      common/cad3rd/cadlbl,cadexe,cadfmt

      LOGICAL OK,unixok

      character outs*124,OUTSTR*124,doit*248
      character cadlbl*20,cadexe*20,cadfmt*5
      character hold32*32,msg*72,msg2*30
      logical found_cad ! to test application existance

C For help messages
      character helpinsub*24 ! subroutine name
      character helptopic*24 ! string (unique) for topic
      integer nbhelp         ! number of help lines found

      helpinsub='cadin'      ! set for subroutine

C Check if running Unix-like OS.
      call isunix(unixok)

C CAD definition and importing.
      helptopic='cad_import_choices'
      call gethelptext(helpinsub,helptopic,nbhelp)

C Tell the user to create a model configuration file.
      CALL PHELPD('CAD warning',nbhelp,'-',0,0,IER)

C Ask user which format of CAD file will be scanned in.
      write(msg2,'(3a)') 'Invoke ',cadexe(1:lnblnk(cadexe)),'?'
      CALL EASKOK(' ',msg2,OK,nbhelp)
      if(OK)then

C Check if CAD tool has been installed. If not do not start.
        if(unixok)then
          found_cad=.false. 
          write(hold32,'(a)') cadexe(1:lnblnk(cadexe))
          call isinstalled(hold32,found_cad)
          if(.NOT.found_cad)then
            write(msg,'(3a)') 'CAD tool ',
     &        cadexe(1:lnblnk(cadexe)),' is not installed.'
            call usrmsg(msg,'Exiting.','W')
            return
          endif
        else
          found_cad=.false.  ! not in Windows
        endif

C Create one or more dummy zones for gtool to work with.
        call edisp(iuout,'CAD tool requires dummy zones to exist')
        call edisp(iuout,'prior to starting. Specify how many to')
        call edisp(iuout,'create. Note: zone names must follow the')
        call edisp(iuout,'same root name pattern as an ESP-r model.')
        call create_dummy_zones()
        write(outs,'(3a)') 'Proceed with ',
     &    cadexe(1:lnblnk(cadexe)),'?'

        call usrmsg(
     &    'After the CAD tool starts the Project Manager will exit.',
     &    'Please restart it when you finish the CAD tool session.',
     &    'W')
        CALL EASKOK(' ',outs,OK,nbhelp)
        IF(OK)then
          call pausems(1000)
          write(doit,'(2a)') cadexe(1:lnblnk(cadexe)),' & '
          call runit(doit,'-')
          call pausems(1000)
          call usrmsg(
     &    'Shutting down Project Manager','...','-')
          call pausems(1000)
          close(ieout)
          if(iecount.eq.0)then     ! if there were no warnings delete
            open(ieout,file=iefile,status='UNKNOWN')
            call efdelet(ieout,istat)
          endif

C Clear allocatable arrays
          CALL DeallocateAllArrays
          CALL EPAGEND
          STOP
        endif
        return
      else
        return
      endif

      end  ! of gtoolin


C ********* rexmpl
C rexmpl: scan the exemplars file, parse information and offer the user
C selection of, or information about, exemplars. 
C useraction is returned - 'ownit' or 'continue'
C Currently supports up to 30 items in a category.
      subroutine rexmpl(itrc,iexfil,exemplar,useraction,ier)
#include "building.h"
#include "model.h"
#include "espriou.h"
C espriou.h provides currentfile, defdbfl exemfl.
      
      integer lnblnk  ! function definition

C Passed parameters.
      integer itrc    ! indicates feedback level
      integer iexfil  ! unit number for the exemplar file
      character exemplar*144 ! full path to selected cfg file
      character useraction*8 ! passed in to specify the task
      integer ier     ! returned as 1 if no exemplars file found

      COMMON/SPAD/MMOD,LIMIT,LIMTTY
      common/OUTIN/IUOUT,IUIN,IEOUT
      integer childterminal  ! picks up mmod from starting of prj
      common/childt/childterminal
      
      integer ncomp,ncon
      common/C1/NCOMP,NCON
      INTEGER :: INDFCG
      common/C6/INDCFG

C Images - included to warn user on subsequent browsing sessions
C to quit the image browser.
C      character imgfmt*4  ! GIF XBMP TIF JPG
C      character imgfoc*4  ! FZON FNET FCTL FDFS
C      character limgfil*72  ! file name (extend to 144 char)
C      character imgdoc*248  ! text associated with image
C      common/imagf/imgfmt(MIMG),imgfoc(MIMG),limgfil(MIMG),imgdoc(MIMG)

      integer noimg  ! number of images
      integer iton   ! zero if images not yet shown, one if yes
      common/imagfi/noimg,iton
      integer IVALS
      dimension IVALS(30)
      character name*70,subitem*70,label*42
      dimension name(30),label(30),inl2in(30),subitem(30)
      character group*42
      dimension group(15)
      character gitem*44,gh*72,cfgitem*144,itemlog*144
      dimension gitem(18),gh(13),cfgitem(30),inl2ngr(30),itemlog(30)
      character OUTSTR*124,WORD*20
      character tgroup*42,key*1,look*70,tlook*70
      character outs*124,rootp*84,separ*42
      character doit*248
      character exdescr*81,doc*72,head*50
      character cfg*144,cfg_path*84,cfg_root*72,doc_file*96
      character ilog*144
      character tmode*12

C Strings for dos version of rootp and upath and a double quote.
      character rootpdos*84,upathdos*72,rootpath*72,modelroot*72
      character dq*1,fs*1
      character lltmp*144,lguess*144,dguess*144
      character thecfgis*72,cfgcfg*42 ! cfg file only & cfg/cfg_file
      character cfgpath*144  ! full path to cfg folder for refocus
      logical XST,unixok
      integer NCOG,ICOG,ICO ! max items and current menu item
      integer loop  ! for looping
      integer iglib   ! if 1 then X11, if 2 then GTK, if 3 then text only.

#ifdef OSI
      integer impx,impy,iwe
#else
      integer*8 impx,impy,iwe
#endif

C For help messages
      common/pophelp/h(60)
      character h*72
      character helpinsub*24 ! subroutine name
      character helptopic*24 ! string (unique) for topic
      integer nbhelp         ! number of help lines found

      helpinsub='rexmpl'     ! set for subroutine

C Make up a double quote and file separator for xcopy.
      useraction ='-'
      call isunix(unixok)
      if(unixok)then
        fs = char(47)
        dq = char(34)
      else
        fs = char(92)
        dq = char(34)
      endif
      iglib = igraphiclib()  ! find out if GTK.
      call terminalmode(childterminal,tmode)

C Open the exemplar file. This is also a jump back point in case
C the user decides to cancel the selection.
  41  IER=0
      if(exemfl(1:1).eq.'/')then
        call FPOPEN(iexfil,ISTAT,1,0,exemfl)
      else
        CALL EFOPSEQ(iexfil,exemfl,1,IER)
      endif
      if(ier.ne.0) return

C Clear some of the local strings.
      do 38 loop=1,30
        IVALS(loop)=0
        name(loop)='  '
        subitem(loop)='  '
        label(loop)='  '
        cfgitem(loop)='  '
        itemlog(loop)='  '
        inl2in(loop)=0
        inl2ngr(loop)=0
        if(loop.le.15) group(loop)='  '
  38  continue

C Read the header and see if '*EXEMPLARS'
      CALL STRIPC(iexfil,OUTSTR,0,ND,0,'exemplars header',IER)
      if(outstr(1:10).ne.'*EXEMPLARS')then
        call edisp(iuout,' ')
        call edisp(iuout,'This is not an exemplars file!')
        ier=1
        return
      endif

C In the first pass make up the menu, beginning with
C an organisational header and the `names' of each exemplar.
C Clear group related help and length of menu.
      do i=1,13
        gh(i)=' '
      enddo
      exemplar='UNKNOWN'; nghelp=0
      ilog='  '; ngr=0; ing=0; inl=0; igrsel=0
      cfgcfg='  '; thecfgis='  '; cfgpath='  '

 12   CALL STRIPC(iexfil,OUTSTR,0,ND,0,'group id',IER)
      if(ier.ne.0)goto 901
      K=0
      CALL EGETW(outstr,K,WORD,'W','group id',IER)
      if(WORD(1:6).eq.'*group')then
        if(ngr.lt.15)then
          ngr=ngr+1; inl=inl+1; ing=ing+1
          CALL EMKEY(ing,KEY,IER)
          call egetrm(outstr,K,group(ngr),'W','group id',IER)
          write(gitem(inl),'(a1,1x,a42)') key, group(ngr)(1:42)
          inl2ngr(inl)=ngr
        endif
        goto 12
      elseif(WORD(1:6).eq.'*label')then
        inl=inl+1
        separ=' '
        call egetrm(outstr,K,separ,'W','label',IER)
        write(gitem(inl),'(2x,a42)') separ(1:42)
        inl2ngr(inl)=0
        goto 12
      else
        goto 12
      endif
 901  CALL ERPFREE(iexfil,ISTAT)

C Present the display.
 142  ICOG= -2
      gitem(inl+1)=  '  _______________________________'
      gitem(inl+2)=  '? help                           '
      gitem(inl+3)=  '- exit menu                      '
      NCOG=inl+3
      IW=0
      do 141 ij=1,NCOG
        if(lnblnk(gitem(ij)).gt.IW)IW=lnblnk(gitem(ij))
 141  continue

C Help text for this menu.
      helptopic='exemplar_menu_help'
      call gethelptext(helpinsub,helptopic,nbhelp)

      if(MMOD.EQ.8)then
        impx=0; impy=0; iwe=iw
        CALL VWMENU('Exemplar categories',gitem,NCOG,
     &    impx,impy,iwe,irpx,irpy,ICOG)
      else
        CALL EMENU('Exemplar categories',gitem,NCOG,ICOG)
      endif

      IF(ICOG.EQ.NCOG)THEN

C Exit from the menu.
        CALL ERPFREE(iexfil,ISTAT)
        return

      ELSEIF(ICOG.EQ.(NCOG-1))THEN

C Help @ db file level.
        helptopic='exemplar_menu_help'
        call gethelptext(helpinsub,helptopic,nbhelp)
        CALL PHELPD('exemplar menu',nbhelp,'-',0,0,IER)

      ELSEIF(ICOG.GT.0.AND.ICOG.LT.(NCOG-2))THEN

C User selected an exemplar category from the menu.
        if(gitem(ICOG)(1:2).eq.'  ')goto 142
        if(inl2ngr(ICOG).gt.0.and.inl2ngr(ICOG).le.15)then
          igrsel=inl2ngr(ICOG)
        else
          goto 142
        endif

C Re-open the exemplar file.
        if(exemfl(1:1).eq.'/')then
          call FPOPEN(iexfil,ISTAT,1,0,exemfl)
        else
          CALL EFOPSEQ(iexfil,exemfl,1,IER)
        endif
 14     CALL STRIPC(iexfil,OUTSTR,0,ND,0,'group id',IER)
        if(ier.ne.0)goto 902
        K=0
        CALL EGETW(outstr,K,WORD,'W','group id',IER)
        if(WORD(1:6).eq.'*group')then
          call egetrm(outstr,K,tgroup,'W','group id',IER)
          if(tgroup(1:42).eq.group(igrsel)(1:42))then

C We have a match with the selected group so now fill the menu
C strings for presentation to the user.
            nghelp=0; in=0; inl=0
 16         CALL STRIPC(iexfil,OUTSTR,0,ND,0,'group',IER)
            if(ier.ne.0)goto 90
            K=0
            CALL EGETW(outstr,K,WORD,'W','group id',IER)
            if(WORD(1:6).eq.'*group')then

C If we have reached the next group then we have scanned far enough.
              goto 90
            elseif(WORD(1:5).eq.'*help')then
              nghelp=nghelp+1
              call egetrm(outstr,K,gh(nghelp),'W','help',IER)
              goto 16
            elseif(WORD(1:5).eq.'*name')then

C Grab the name and assign a key to each item in the category.
C if there are more than 26 items then present as a multi-page
C selection list.
              in=in+1; inl=inl+1
              call egetrm(outstr,K,name(in),'W','name',IER)
              write(label(inl),'(a)')'-'
              CALL EMKEY(in,KEY,IER)
              write(subitem(inl),'(a)') name(in)
              inl2in(inl)=in
              goto 16
            elseif(WORD(1:4).eq.'*cfg')then
              call egetrm(outstr,K,cfgitem(in),'W','cfgitem',IER)
              goto 16
            elseif(WORD(1:4).eq.'*log')then

C In case there is a log file mentioned.
              call egetrm(outstr,K,itemlog(in),'W','item log',IER)
              goto 16
            elseif(WORD(1:6).eq.'*label')then
              inl=inl+1
              call egetrm(outstr,K,label(inl),'W','label',IER)
              inl2in(inl)=0
              goto 16
            else
              goto 16
            endif
  90        CALL ERPFREE(iexfil,ISTAT)
            goto 42
          else
            goto 14
          endif
          goto 14
        else
          goto 14
        endif
      endif
 902  CALL ERPFREE(iexfil,ISTAT)
      goto 142

C Present the list of this category of exemplars. If there
C are more than a-z to choose from use alternative.
  42  continue
      if(in.eq.0)then
        CALL ERPFREE(iexfil,ISTAT)
        goto 41
      endif

C Help text for this menu. Append the help text read from the
C exemplars file to the existing help text found in the esruprj.help
C file for this topic.
      iwsub=0
      do ij=1,in
        if(lnblnk(subitem(ij)).gt.iwsub)iwsub=lnblnk(subitem(ij))
      enddo
      helptopic='exemplar_set_browse'
      call gethelptext(helpinsub,helptopic,nbhelp)
      if(nghelp.gt.0)then
        do ih=1,nghelp
          lex=MIN0(lnblnk(gh(ih)),72)
          write(h(nbhelp+ih),'(a)')gh(ih)(1:lex)
        enddo
      endif
      write(head,'(2a)') 'Exemplar group: ',
     &           group(igrsel)(1:lnblnk(group(igrsel)))
      IX=1
      CALL EPICKS(IX,IVALS,' ',' ',
     &  iwsub,in,subitem,head,IER,nbhelp+nghelp)
      IF(IX.EQ.0)then

C Exit from detail level back to group level.
        CALL ERPFREE(iexfil,ISTAT)
        goto 41
      endif
      ICO=IVALS(1)
      if(ICO.eq.0)then

C Exit from detail level back to group level.
        CALL ERPFREE(iexfil,ISTAT)
        goto 41
      endif
      if(inl2in(ICO).eq.0)goto 42
      look=name(inl2in(ICO))

C Open the related configuration file (if it exists), find its log file.
      CALL ERPFREE(iexfil+1,ISTAT)
      cfg=cfgitem(inl2in(ICO))
      INQUIRE (FILE=cfg,EXIST=XST)
      if(XST)then
        call fdroot(cfg,cfg_path,cfg_root)
        call FPOPEN(iexfil+1,ISTAT,1,0,cfg)
  242   CALL STRIPC(iexfil+1,OUTSTR,0,ND,0,'cfg scan',IER)
        if(ier.ne.0)goto 902
        if(OUTSTR(1:10).eq.'* PROJ LOG')then
          CALL STRIPC(iexfil+1,OUTSTR,0,ND,1,'log file',IER)
          doc=OUTSTR(1:LNBLNK(OUTSTR))
          write(doc_file,'(a,a)')cfg_path(1:lnblnk(cfg_path)),
     &      doc(1:lnblnk(doc))
          CALL ERPFREE(iexfil+1,ISTAT)
          call tstamp('>','PRJ: list exemplar documentation')
          call edisp(iuout,' ')
          CALL LISTAS(iexfil+1,doc_file,IER)
        else
          goto 242  ! read another line
        endif

C Attempt to access the exempler.
        continue

      else
         call usrmsg('This exemplar is currently off-line.',
     &     'Arrange access via your ESP-r administrator.','W')
         exemplar='UNKNOWN'
         CALL ERPFREE(iexfil+1,ISTAT)
         rewind(iexfil,ERR=92)
         goto 41
      endif
      
C Process the user selection.
      IF(ICO.GE.1.AND.ICO.LE.(IN))THEN

C User selected an exemplar, attempt to access it.
        if(inl2in(ICO).eq.0)goto 42
        look=name(inl2in(ICO))
        if(exemfl(1:1).eq.'/')then
          call FPOPEN(iexfil,ISTAT,1,0,exemfl)
        else
          CALL EFOPSEQ(iexfil,exemfl,1,IER)
        endif
        CALL STRIPC(iexfil,OUTSTR,0,ND,0,'exemplar 1st line',IER)
        CALL STRIPC(iexfil,OUTSTR,0,ND,0,'group id',IER)
  3     CALL STRIPC(iexfil,OUTSTR,0,ND,0,'item lines',IER)
        if(ier.ne.0)goto 92

        exemplar=' '; rootp=' '
        K=0
        CALL EGETW(outstr,K,WORD,'W','*name tag',IER)
        if(WORD(1:5).eq.'*name')then
          call egetrm(outstr,K,tlook,'W','name',IER)
          if(tlook.eq.look)then
            lex=MIN0(lnblnk(tlook),70)
            write(exdescr,'(a,a)') 'Model: ',tlook(1:lex)
  4         CALL STRIPC(iexfil,OUTSTR,0,ND,0,'item line',IER)
            if(ier.ne.0)goto 92
            K=0
            CALL EGETW(outstr,K,WORD,'W','*cfg tag',IER)
            if(WORD(1:4).eq.'*end')then
              goto 92
            elseif(WORD(1:4).eq.'*cfg')then
              call egetrm(outstr,K,exemplar,'W','name',IER)
              if(iglib.eq.2.or.iglib.eq.3)then
                call edisp(iuout,'  ')
                call edisp(iuout,'Exemplar model is:')  ! echo if gtk
                call edisp(iuout,exemplar)
              endif

C If running on a Windows machine get the last token within the
C string exemplar, which is the name of the configuration file.
C Work backwards from the end of the string till a file separator
C of either DOS or Unix type is reached.
              call isunix(unixok)
              if(.NOT.unixok)then
                thecfgis=' '
                ilexem = lnblnk(exemplar)  ! last character
                ipos = ilexem
 344            continue
                ipos=ipos-1
                if(exemplar(ipos:ipos).eq.fs)then
                  write(thecfgis,'(a)') exemplar(ipos+1:ilexem)
                  write(cfgpath,'(a)') exemplar(1:ipos)
                  if(itrc.ge.1)then
                    call edisp(iuout,'  ')
                    call edisp(iuout,'The model configuration file is:')
                    call edisp(iuout,thecfgis)
                    call edisp(iuout,'The model cfg folder is:')
                    call edisp(iuout,cfgpath)
                  endif
                elseif(exemplar(ipos:ipos).eq.'/')then
                  write(thecfgis,'(a)') exemplar(ipos+1:ilexem)
                  write(cfgpath,'(a)') exemplar(1:ipos)
                  if(itrc.ge.1)then
                    call edisp(iuout,'  ')
                    call edisp(iuout,'The model configuration file is:')
                    call edisp(iuout,thecfgis)
                    call edisp(iuout,'The model cfg folder is:')
                    call edisp(iuout,cfgpath)
                  endif
                else
                  if(ipos.eq.1)then
                    call edisp(iuout,
     &                'Warning: No folder separator found.')
                  elseif(ipos.gt.1)then
                    goto 344
                  endif
                endif
              endif

C Check to see if the configuration file exists. If it does not then
C warn the user that the model is off-line and return
C to the list of exemplars.
              CALL ERPFREE(iexfil+1,ISTAT)
              INQUIRE (FILE=exemplar,EXIST=XST)
              if(.NOT.XST)then
                call tstamp('>','prj: found off-line exemplar')
                call usrmsg('The selected exemplar is off-line.',
     &             'Arrange access via your ESP-r administrator.','W')
                exemplar='UNKNOWN'
                CALL ERPFREE(iexfil+1,ISTAT)
                rewind(iexfil,ERR=92)
                goto 41
              endif
              goto 4
            elseif(WORD(1:5).eq.'*root')then
              call egetrm(outstr,K,rootp,'W','exem root dir',IER)
              goto 4
            elseif(WORD(1:4).eq.'*log')then

C In case there is a log file mentioned.
              call egetrm(outstr,K,ilog,'W','item log',IER)
              goto 4
            elseif(WORD(1:5).eq.'*item')then

C If a log file mentioned in exemplars file show it.
              if(ilog(1:2).ne.'  ')then
                call edisp(iuout,' ')
                CALL LISTAS(iexfil+1,ilog,IER)
              endif

C Reached the end of the item, check if it is to be chosen.
C If there are zones existing, remove them from memory.
              CALL EASKAB(exdescr,'Options:',
     &          'proceed','cancel',IW,0)
              if(IW.eq.1)then
                if(upath.ne.' ')then

C See if can derive the local configuration name. cfg_root is the last
C word in the path to the exemplar root folder.
                  call fdroot(rootp,cfg_path,cfg_root)
                  lex=lnblnk(exemplar)
                  lcfgr=lnblnk(cfg_root)

C Loop along the exemplar string to locate the position of the root
C folder in the path of the exemplar.
                  ilex=1
                  do 342 ilex=1,lex-lcfgr
                    if(cfg_root(1:lcfgr).eq.
     &                exemplar(ilex:ilex+lcfgr-1))then
                      goto 343
                    endif
  342             continue
  343             continue

C The balance of the exemplar string after the cfg_root.
                  write(cfgcfg,'(a)') exemplar(ilex+lcfgr:lex)
                  write(thecfgis,'(a)') cfgcfg(6:lnblnk(cfgcfg))
                  call isunix(unixok)
                  if(unixok)then

C Offer user choice of home folder plus the exemplar root name. In GTK
C there is a file browser returning folder name via lltmp. If lltmp
C folder does not exist then need to create it.
C If prj was started in a Models folder prepend to path.
                    if(pwdinitial(lnpwdi-5:lnpwdi).eq.'Models')then
                      write(lguess,'(3a)') pwdinitial(1:lnpwdi),fs,
     &                  cfg_root(1:lcfgr)
                    else
                      write(lguess,'(3a)') upath(1:lnblnk(upath)),fs,
     &                  cfg_root(1:lcfgr)
                    endif
                    write(dguess,'(3a)') upath(1:lnblnk(upath)),fs,
     &                cfg_root(1:lcfgr)
                    call edisp(iuout,' ')  ! echo a blank line
                    CALL EASKXORGTKF(lguess,' ','Copy model to:',               
     &                dguess,lltmp,'models folder',IER,0)

                    write(doit,'(2a)') 'mkdir ',lltmp(1:lnblnk(lltmp))
                    call runit(doit,'-')
                    call pausems(100)

C Add /cfg to lltmp to get the configuration folder
                    write(cfgpath,'(3a)')lltmp(1:lnblnk(lltmp)),fs,'cfg'

C Make up the command to copy from the exemplar source folder to the
C users defined destination. If we add cfgcfg to lltmp we should have
C the folder where the model cfg file is located.
                    write(doit,'(5a)') 'cp -r ',rootp(1:lnblnk(rootp)),
     &                fs,'*  ',lltmp(1:lnblnk(lltmp))
                  else

C Non-unix, if lnblankrp or lnblankup is one then there are no spaces in the 
C string so double quotes are not required.  Change any '/'
C in rootp and upath to '\'. If prj started in folder Models then insert
C this into the suggested path. 
                    rootpdos=' '
                    upathdos=' '
                    lnrp=lnblnk(rootp)
                    lnblankrp=iprevblnk(rootp,lnrp)

                    call backslashit(rootp,rootpdos)
                    if(itrc.gt.0)then
                      write(outs,'(2a)') 'Getting model from...',
     &                  rootpdos(1:lnrp)
                      call edisp(iuout,outs)
                    endif

C For an xcopy command we need the destination to also include the
C model folder so use fdroot to extract this from rootp.
                    call fdroot(rootpdos,rootpath,modelroot)

                    call edisp(iuout,' ')

C Offer user choice of home folder plus the exemplar root name. In GTK
C there is a file browser returning folder name via lltmp.
                    if(pwdinitial(lnpwdi-5:lnpwdi).eq.'Models')then
                      write(lguess,'(3a)') pwdinitial(1:lnpwdi),fs,
     &                  cfg_root(1:lcfgr)
                    else
                      write(lguess,'(3a)') upath(1:lnblnk(upath)),fs,
     &                  cfg_root(1:lcfgr)
                    endif
                    write(dguess,'(3a)') upath(1:lnblnk(upath)),fs,
     &                cfg_root(1:lcfgr)
                    call edisp(iuout,' ')  ! echo a blank line
                    CALL EASKXORGTKF(lguess,' ','Copy model to:',
     &                dguess,lltmp,'models folder',IER,0)

                    lnup=lnblnk(lltmp)
                    lnblankup=iprevblnk(lltmp,lnup)
                    call backslashit(lltmp,upathdos)

C If upathdos folder does not exist then need to create it.
                    if(lnblankup.gt.1)then

C There are blanks in destination strings.
                      write(doit,'(4a)') 'mkdir ',dq,upathdos(1:lnup),dq
                    elseif(lnblankup.eq.1)then

C There are no blanks in destination string.
                      write(doit,'(2a)') 'mkdir ',upathdos(1:lnup)
                    endif
                    call runit(doit,'-')

                    if(lnblankrp.gt.1.and.lnblankup.gt.1)then

C There are blanks in both source or destination strings.
                      write(doit,'(11a)') 'xcopy /e /i /y /c ',dq,
     &                  rootpdos(1:lnrp),fs,'*',dq,'  ',dq,
     &                  upathdos(1:lnup),fs,dq
                    elseif(lnblankrp.eq.1.and.lnblankup.gt.1)then

C There are blanks in destination string.
                      write(doit,'(9a)') 'xcopy /e /i /y /c ',
     &                  rootpdos(1:lnrp),fs,'*  ',dq,
     &                  upathdos(1:lnup),fs,dq
                    elseif(lnblankrp.gt.1.and.lnblankup.eq.1)then

C There are blanks in source string.
                      write(doit,'(9a)') 'xcopy /e /i /y /c ',dq,
     &                  rootpdos(1:lnrp),fs,'*',dq,'  ',
     &                  upathdos(1:lnup),fs
                    elseif(lnblankrp.eq.1.and.lnblankup.eq.1)then

C There are no blanks in either source or destination strings.
                      write(doit,'(6a)') 'xcopy /e /i /y /c ',
     &                  rootpdos(1:lnrp),fs,'*  ',
     &                  upathdos(1:lnup),fs
                    endif
                  endif  ! of non-unix

C If trace is verbose echo the command.
C                  if(itrc.ge.0)then
C                    call edisp(iuout,
C     &                'Copying exemplar to destination folder via:')
C                    call edisp248(iuout,doit,100)
C                    call pausems(5000)
C                  endif
                  call runit(doit,'-')
                  if(unixok)then

C If Unix-like then append tail of exemplar string (cfgcfg) to the users selected
C path to get the new exemplar name.
                    write(exemplar,'(2a)') lltmp(1:lnblnk(lltmp)),
     &                cfgcfg(1:lnblnk(cfgcfg))
                  else

C If Windows based then append cfg & cfg file string to the path to make
C up the exemplar name (this assumes a standard layout of model).
                    write(exemplar,'(5a)')upathdos(1:lnblnk(upathdos)),
     &                fs,'cfg',fs,thecfgis(1:lnblnk(thecfgis))
                    write(cfgpath,'(3a)')upathdos(1:lnblnk(upathdos)),
     &                fs,'cfg'
                  endif

                  call tstamp('>','prj: owning exemplar')
                  call tstamp('>',exemplar)
                  useraction='ownit   '
                endif

C Suggest to user to refocus prj within the model cfg folder at
C this specific cfg file.
                call refocus(cfgpath,thecfgis,iret)
                if(iret.eq.0)then
                  continue
                elseif(iret.eq.-2)then
                  call edisp(iuout,
     &       'Features using Radiance calls may not work as expected.')
                endif

              elseif(IW.eq.2)then

C User asked to cancel current selection, rewind the file and jump back to
C the top menu to re-resent groups.
                exemplar='UNKNOWN'
                REWIND(iexfil,ERR=92)
                useraction='continue'
                goto 41
              endif

C If there was a prior model in memory (e.g. some zones or a plant
C only model) clear it.
              CALL ERPFREE(iexfil,ISTAT)
              if(ncomp.gt.0.or.INDCFG.eq.2)then
                call usrmsg(' ',
     &                      'Clearing previous model.','-')
                call startbuffer()
                call clrprb
              endif

C If the image browser was started in a previous session warn
C the user to kill it.
              if(iton.eq.1)then
                noimg=0
                iton=0
              endif
              return  ! we now know the model we want to load
            else
              goto 4
            endif  ! of finding a matching item token
          else
            goto 3
          endif
        else
          goto 3
        endif      ! of processing a name token
      else
        goto 42
      endif        ! of processing user selection
      goto 42

  92  CALL ERPFREE(iexfil,ISTAT)
      goto 42

      end

C ***** wtarentry
C Write entry in tar names file taking into account case where file
C name begins with ../ .
      subroutine wtarentry(iunit,pathome,tr,ltr,tp,ltp,lf)
      integer lnblnk  ! function definition

C Parameters
      integer iunit   ! file unit for written entries
      logical pathome ! true uses tp, false uses tr
      character tp*32 ! path including final file separator
      character tr*32 ! path without 
      character*(*) lf ! file name for entries
      integer ltr,ltp  ! length of strings tp and tr
      
      character fs*1
      logical unixok

      call isunix(unixok)
      if(unixok)then
        fs = char(47)
      else
        fs = char(92)
      endif

      llf=lnblnk(lf)
      lltr=ltr
      lltp=ltp
      if(llf.eq.0)then
        return
      elseif(lf(1:7).eq.'UNKNOWN')then
        return
      elseif(lf(1:3).eq.'../')then
        if(.NOT.pathome)then
          write(IUNIT,'(3a)',IOSTAT=IOS,ERR=2) tr(1:lltr),fs,lf(4:llf)
        else
          write(IUNIT,'(3a)',IOSTAT=IOS,ERR=2) tp(1:lltp),fs,lf(1:llf)
        endif
      elseif(lf(1:2).eq.'./')then
        write(IUNIT,'(3a)',IOSTAT=IOS,ERR=2) tp(1:lltp),fs,lf(3:llf)
      else
        write(IUNIT,'(3a)',IOSTAT=IOS,ERR=2) tp(1:lltp),fs,lf(1:llf)
      endif
      return

C Error.
    2 if(IOS.eq.2)then
        CALL USRMSG(
     &  'Permissions problem when writing tar names entry',lf,'W')
      else
        CALL USRMSG('Problem when writing tar names entry',lf,'W')
      endif
      call ERPFREE(IUNIT,ISTAT)
      return

      end

C e2vef takes an esp-r problem definition and creates a vef
C vertex-edge-face file in file unit ioout (assumed to already be opened).
      subroutine e2vef(itrc,ioout,act,icol)
#include "building.h"
#include "model.h"
#include "geometry.h"

C Parameter
      character act*1  ! g for GEN p for PLA
      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      COMMON/FILEP/IFIL

      COMMON/C1/NCOMP,NCON
      COMMON/GB1/XB(12),YB(12),ZB(12),JVNB(6,4)
      integer IZSTOCN
      COMMON/C24/IZSTOCN(MCOM,MS)

      dimension ndum(MS)

      character outs*124,zsn*28

C Setup standard assumptions.
      IFIL=11

C Assume geometry or obstructions on IFIL+8 
      ITA1 = IFIL+8

C If VEF output required.
      INPIC=NCOMP
      do 4 mz=1,inpic
        newfoc=mz

C Read in the zone geometry.
        WRITE(outs,'(a,a)')' Scanning : ',LGEOM(newfoc)
        CALL edisp(iuout,outs)

C Depending on version of geometry file scan it.
        if(gversion(mz).lt.1.1) then
          CALL EGOMIN(ITA1,LGEOM(newfoc),newfoc,1,ITRC,IUOUT,IER)
          if(iobs(mz).eq.0)then
            continue   ! no obstructions
          elseif(iobs(mz).eq.1)then
            CALL ERPFREE(ITA1,ISTAT)
            CALL EGOMST(ITA1,newfoc,ZOBS(newfoc),0,ITRC,ITRU,IER)
          endif
        elseif(gversion(mz).ge.1.1) then
          call georead(ITA1,LGEOM(newfoc),newfoc,1,iuout,ier)
          if(iobs(mz).eq.0)then
            continue   ! no obstructions
          elseif(iobs(mz).eq.1)then
            CALL ERPFREE(ITA1,ISTAT)
            CALL EGOMST(ITA1,newfoc,ZOBS(newfoc),0,ITRC,ITRU,IER)
          elseif(iobs(mz).eq.2)then
            continue   ! obs within geo file
          endif
        endif

C Write colour directive to VEF file.
C << Currently hard-coded for 1 colour per zone.
C << Currently limited to default 15 colours in Akrobat.
        if (icol.eq.1) then
          WRITE(ioout,'(a)')' COL'
C Calculate colour reference. Repeat colours every 15 zones.
          icolref=newfoc-15*(((newfoc+14)/15)-1)
          WRITE(ioout,'(i2.1)')icolref
        endif

C Write surface information to the VEF file. (Rotation set to 0 deg).
C If act is g then use GEN bodies.
        if(act.eq.'g'.or.act.eq.'G')then
          WRITE(ioout,'(a,a12)')'GEN ',zname(newfoc)
          WRITE(ioout,'(2I4,a)')NTV,NSUR,' 0'
          DO J=1,NTV
            WRITE(ioout,'(3F10.5)')X(J),Y(J),Z(J)
          ENDDO

          DO J=1,NSUR
            WRITE(ioout,'(I5,102I4)')NVER(J),(JVN(J,K),K=1,NVER(J))
            ndum(J)=0
          ENDDO

          write(ioout,5611)(ndum(I),I=1,NSUR)
          write(ioout,5611)(ndum(I),I=1,NSUR)
5611      FORMAT(1X,32(I2))
          write(ioout,'(a)')' 3 0 0 0 '
        elseif(act.eq.'p'.or.act.eq.'P')then

C Write each surface as a PLA type. Proceed with VIS 0 or 2 depending
C on if it is transparent. Also mark each with an identity zone:name.
C Note: does not know quite how to treat CFC and CFC2 so just treats
C as opaque.
          do J=1,NSUR
            icn=IZSTOCN(newfoc,j)  ! get connection
            lsn=lnblnk(sname(newfoc,j))
            write(zsn,'(3a)') zname(newfoc)(1:lnzname(newfoc)),
     &        ':',sname(newfoc,j)(1:lsn)
            if (SOTF(newfoc,j)(1:4).ne.'OPAQ'.and.
     &          SOTF(newfoc,j)(1:4).ne.'CFC '.and.
     &          SOTF(newfoc,j)(1:4).ne.'CFC2') then
              write(ioout,'(a)')' VIS '
              write(ioout,'(a)')' 0 '
            else
              write(ioout,'(a)')' VIS '
              write(ioout,'(a)')' 2 '
            endif
            WRITE(ioout,'(2a)')'PLA ',zsn(1:lnblnk(zsn))
            WRITE(ioout,'(I4)') NVER(J)
            DO I=1,NVER(J)
              WRITE(ioout,'(3F11.5)')
     &          X(JVN(J,I)),Y(JVN(J,I)),Z(JVN(J,I))
            enddo
          enddo
        endif

C If there are obstructions convert them to general polygons
C and write them out.
        if(iobs(mz).eq.0)then
          continue
        elseif(iobs(mz).eq.1.or.iobs(mz).eq.2)then

          DO 301 IB=1,nbobs(mz)
            if(BLOCKTYP(mz,IB)(1:4).eq.'obs ')then
              CALL CNVBLK(XOB(mz,IB),YOB(mz,IB),ZOB(mz,IB),
     &          DXOB(mz,IB),DYOB(mz,IB),DZOB(mz,IB),
     &          BANGOB(mz,IB,1))
            elseif(BLOCKTYP(mz,IB)(1:4).eq.'obs3')then
              CALL CNVBLK3A(XOB(mz,IB),YOB(mz,IB),ZOB(mz,IB),
     &          DXOB(mz,IB),DYOB(mz,IB),DZOB(mz,IB),
     &          BANGOB(mz,IB,1),BANGOB(mz,IB,2),BANGOB(mz,IB,3))
            elseif(BLOCKTYP(mz,IB)(1:4).eq.'obsp')then
              call CNVBLKP(mz,IB) ! convert obsp type.
            endif
C << and then what... look in e2rform.F near line 1120
C << there are 8 vertices and 6 surfaces
            lnz=lnblnk(zname(mz))
            WRITE(ioout,'(4a)') 'GEN OBS:',zname(mz)(1:lnz),':',
     &        BLOCKNAME(mz,IB)(1:LNBLOCKNAME(mz,IB))
            WRITE(ioout,'(a)')' 8  6  0'
            DO 41 J=1,8
              if (BLOCKTYP(mz,ib).eq.'obs ') then
                WRITE(ioout,'(3F10.5)')XB(J),YB(J),ZB(J)
              elseif (BLOCKTYP(mz,ib).eq.'obs3') then
                WRITE(ioout,'(3F10.5)')XB(J),YB(J),ZB(J)
              elseif (BLOCKTYP(mz,ib).eq.'obsp') then
                WRITE(ioout,'(3F10.5)')XBP(mz,IB,J),YBP(mz,IB,J),
     &            ZBP(mz,IB,J)
              endif
   41       continue
            DO 51 J=1,6
              WRITE(ioout,'(a,4I3)') ' 4 ',(JVNB(J,K),K=1,4)
              ndum(J)=0
   51       continue

            write(ioout,5611)(ndum(I),I=1,6)
            write(ioout,5611)(ndum(I),I=1,6)
            write(ioout,'(a)')' 3 0 0 0 '

  301     CONTINUE
        endif

   4  continue

      return
      END

C ************ create_dummy_zones
      subroutine create_dummy_zones()
#include "building.h"
#include "model.h"
#include "geometry.h"
#include "prj3dv.h"
      
      integer lnblnk  ! function definition

      common/OUTIN/IUOUT,IUIN,IEOUT
      common/FILEP/IFIL

      integer ncomp,ncon
      COMMON/C1/NCOMP,NCON
      COMMON/C3/IC1(MCON),IE1(MCON),ICT(MCON),IC2(MCON),IE2(MCON)
      integer IZSTOCN
      COMMON/C24/IZSTOCN(MCOM,MS)

      logical unixok
      character fs*1
      character ZN*12,DFILE*72,CFILE*72,OFILE*72

C For help messages
      character helpinsub*24 ! subroutine name
      character helptopic*24 ! string (unique) for topic
      integer nbhelp         ! number of help lines found

      helpinsub='SETUP'      ! set for subroutine

C Check if Unix-based or DOS based.
      call isunix(unixok)
      if(unixok)then
        fs = char(47)
      else
        fs = char(92)
      endif

C Ask how many dummy zones to create.
      itocreate=3
      CALL EASKI(itocreate,' ',
     &  'Number of zones to create in CAD tool?',
     &  1,'F',10,'W',4,'nb to create',IERI,nbhelp)
      if(NCOMP.eq.0)then
        loopstart=1; loopend=itocreate
      else
        loopstart=NCOMP+1; loopend=NCOMP+itocreate+1
      endif

C All dummy zones are 1m cubes.
      XO1=0.0; YO1=0.0; ZO1=0.0
      DX1=1.0; DY1=1.0; DZ1=1.0
      ANGR=0.
      ANGE=0.

      do icomp=loopstart,loopend

C Set default name of the zone and zone geometry file.
        write(ZN,'(2a,i2.2)') cfgroot(1:lnblnk(cfgroot)),'_',icomp
        call st2name(ZN,zname(ICOMP))
        lnzname(icomp)=lnblnk(zname(ICOMP))  ! update this string length.

C << Unix vs NT needs to be updated >>
        if(zonepth(1:2).eq.'  '.or.zonepth(1:2).eq.'./')then
          WRITE(DFILE,'(2a)') zname(ICOMP)(1:lnzname(ICOMP)),'.geo'
          WRITE(CFILE,'(2a)') zname(ICOMP)(1:lnzname(ICOMP)),'.con'
          WRITE(OFILE,'(2a)') zname(ICOMP)(1:lnzname(ICOMP)),'.opr'
        else
          WRITE(DFILE,'(4a)') zonepth(1:lnblnk(zonepth)),fs,
     &     zname(ICOMP)(1:lnzname(ICOMP)),'.geo'
          WRITE(CFILE,'(4a)') zonepth(1:lnblnk(zonepth)),fs,
     &     zname(ICOMP)(1:lnzname(ICOMP)),'.con'
          WRITE(OFILE,'(4a)') zonepth(1:lnblnk(zonepth)),fs,
     &     zname(ICOMP)(1:lnzname(ICOMP)),'.opr'
        endif
        LGEOM(ICOMP)=DFILE
        LTHRM(ICOMP)=CFILE
        LPROJ(ICOMP)=OFILE

        XO1=XO1+1.0  ! shift each dummy zone by 1m in X.
          
C Convert box into a gen description using two rotation angles.
C Pass in critical dimensions and expect the data back via common G1.
        CALL ERECC3A(XO1,YO1,ZO1,DX1,DY1,DZ1,ANGR,ANGE,0.0)
        NSUR=6
        IUZBASEA(icomp)=0
        IBASES(ICOMP,1)=6
        MODBND=.TRUE.; MODIFYVIEW=.TRUE.; MODLEN=.TRUE.
        iZBFLG(ICOMP)=0

C Use subroutine filsur to instantiate G5 common block.
        CALL FILSUR(ICOMP,0)

C Update the connection list and hash.
        ICCC=NCON
        DO ICC=1,NSUR
          ICCC=ICCC+1
          IC1(ICCC)=ICOMP; IE1(ICCC)=ICC; ICT(ICCC)=-1   ! mark as unknown at this time.
          IC2(ICCC)=0; IE2(ICCC)=0
          IZSTOCN(icomp,icc)=iccc
        ENDDO
        NCON=ICCC

        write(zdesc(ICOMP),'(2a)') zname(ICOMP)(1:lnzname(ICOMP)),
     &    ' describes a '
        lnzdesc(ICOMP)=lnblnk(zdesc(ICOMP))  ! update the length of this string.
        CTYPE(icomp)='GEN '
        NDP(ICOMP)=3
        IDPN(ICOMP,1)=0; IDPN(ICOMP,2)=0; IDPN(ICOMP,3)=0
        NZSUR(ICOMP)=NSUR  ! update nzsur() it is needed by zgupdate.
        NZTV(ICOMP)=NTV

C Update the G7 common blocks and then assign ZBASEA.
        call zgupdate(1,ICOMP,ier)

        ZBASEA(icomp)= SNA(icomp,IBASES(ICOMP,1))

C Update the global coordinates for this zones surfaces so that
C the subsequent wireframe image can be drawn and the bounds of
C the zone can be calculated. Also instantiate isznver and iszjvn
C for the new zone (as is done in subroutine zdata).
        DO 40 J=1,NZTV(ICOMP)
          szcoords(ICOMP,J,1)=X(J)
          szcoords(ICOMP,J,2)=Y(J)
          szcoords(ICOMP,J,3)=Z(J)
   40   CONTINUE
        DO 50 J=1,nzsur(ICOMP)
          icc=IZSTOCN(icomp,j)
          if(icc.ne.0)then
            isznver(icomp,j)=NVER(J)
            N = isznver(icomp,j)
            DO 60 K=1,N
              iszjvn(icomp,j,K)=JVN(J,K)
   60       CONTINUE
          endif
   50   continue

C Save this to file before passing into the geometry editing facility.
        gversion(icomp) =1.1
        call geowrite(IFIL+2,LGEOM(ICOMP),ICOMP,iuout,3,IER)
        NCOMP=NCOMP+1
        NCCODE(ICOMP)=NCOMP                                                                                                                
        CALL EMKCFG('-',IER)                                                       
      enddo

      return
      end  ! of create_dummy_zones

