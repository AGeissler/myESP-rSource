C This file is part of the ESP-r system.
C Copyright Energy Systems Research Unit, University of
C Strathclyde, Glasgow, Scotland, 2001-7.

C ESP-r is free software.  You can redistribute it and/or
C modify it under the terms of the GNU General Public
C Licence as published by the Free Software Foundation 
C (version 2 or later).

C ESP-r is distributed in the hope that it will be useful
C but WITHOUT ANY WARRANTY; without even the implied
C warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
C PURPOSE. See the GNU General Public Licence for more
C details.

C You should have received a copy of the GNU General Public
C Licence along with ESP-r. If not, write to the Free
C Software Foundation, Inc., 59 Temple Place, Suite 330,
C Boston, MA 02111-1307 USA.

C esru_ask.F includes notices used by ESP-r modules.

C   cpwpk()    Displays a copyright notice via help text data structure.
C   askabout() Creates text for a popup notice for each module.

C ***** copyright notice.
C cpwpk displays a copyright notice via help text data structure.
C It saves and recover any current contextual help messages.
C The integer parameter il is not used.
      subroutine cpwpk(il)
      common/pophelp/h(60)
      dimension dh(60)
      character H*72,dh*72

C First copy any existing help into dh array.
      do 50 i=1,60
        dh(i)=h(i)
   50 continue

      H(1)='ESP-r is distributed as free software under the terms of'
      H(2)='the GNU General Public License as published by the Free'
      H(3)='Software Foundation (Version 2 or later), Details available'
      H(4)='at https://opensource.org/licenses/gpl-license.'
      CALL PHELPD('ESP-r license',4,'-',0,0,IER)

C Recover previous help via dh array.
      do 51 i=1,60
        h(i)= dh(i)
   51 continue

      return
      end

C **** askabout() Creates text for a popup notice for each module. ***
C The intent of askabout is to have one place for the opening help
C message of each applicaton (and thus support an application independent
C static help menu under GTK).
C If initonly is 1 then only assign the h() array and do not display
C If initonly is 0 save prior h(), update h(), display and re-establish prior.
C module (char *4) identifies the ESP-r module.
      subroutine askabout(module,initonly)
      common/pophelp/h(60)
      dimension dh(60)
      character H*72,dh*72,module*4

C First copy any existing help into dh array.
      if(initonly.eq.0)then

        do 50 i=1,60
          dh(i)=h(i)
   50   continue
      else

      endif

      if(module(1:3).eq.'aco')then
        H(1)='This module supports an appraisal of the acoustic'
        H(2)='reverberation time within a zone.'
        if(initonly.eq.0)then
          CALL PHELPD('aco opening',2,'-',0,0,IER)
        endif
      elseif(module(1:3).eq.'b2e')then
        H(1)='This module manages the process of updating measured'
        H(2)='data when running ESP-r in BEMS emulation mode.'
        if(initonly.eq.0)then
          CALL PHELPD('b2e opening',2,'-',0,0,IER)
        endif
      elseif(module(1:3).eq.'bps')then
        H(1)='This module subjects a model to an integrated simulation'
        H(2)='depending on the available domain information relating'
        H(3)='to the building(s), plant and/or control elements.'
        if(initonly.eq.0)then
          CALL PHELPD('bps opening',3,'-',0,0,IER)
        endif
      elseif(module(1:3).eq.'cfg')then
        H(1)='This module checks the building topology by matching'
        H(2)='zone surfaces. Where a probable match is found, it is'
        H(3)='presented for user confirmation or association with a'
        H(4)='standard boundary condition (external, ground, adiabetic'
        H(5)='etc.). The test tolerance may be set to accommodate'
        H(6)='minor differences between zone descriptions.'
        if(initonly.eq.0)then
          CALL PHELPD('cfg opening',6,'-',0,0,IER)
        endif
      elseif(module(1:3).eq.'cdb')then
        H(1)='This module manages icons relating to flow, plant'
        H(2)='and electrical components as used in graphical'
        H(3)='network definition.'
        if(initonly.eq.0)then
          CALL PHELPD('cdb opening',3,'-',0,0,IER)
        endif
C      elseif(module(1:3).eq.'clm')then
C Help moved to esruclm.help
C        if(initonly.eq.0)then
C          CALL PHELPD('clm opening',24,'clm_dbitems ',0,0,IER)
C        endif
      elseif(module(1:3).eq.'dbm')then
        H(1)='This module provides information on the various'
        H(2)='databases used by ESP-r modules.'
        if(initonly.eq.0)then
          CALL PHELPD('dbm opening',2,'-',0,0,IER)
        endif
      elseif(module(1:3).eq.'dfs')then
        H(1)='This module allows an independent solution of a '
        H(2)='CFD domains (i.e. with no thermal zone conflation).'
        if(initonly.eq.0)then
          CALL PHELPD('dfs opening',2,'-',0,0,IER)
        endif
      elseif(module(1:3).eq.'dfv')then
        H(1)='This module allows CFD flow results to be viewed as'
        H(2)='a simulation progresses.'
        if(initonly.eq.0)then
          CALL PHELPD('dfv opening',2,'-',0,0,IER)
        endif
      elseif(module(1:4).eq.'ecnv')then
        h(1)='This module controls various data conversions, e.g. from'
        h(2)='an ESP-r geometry format to DXF.'
        if(initonly.eq.0)then
          CALL PHELPD('ecnv opening',2,'-',0,0,IER)
        endif
      elseif(module(1:3).eq.'e2r')then
        H(1)='This module converts an ESP-r model definition into'
        H(2)='a form that can be interpreted by the RADIANCE lighting'
        H(3)='simulation package. It then executes RADIANCE as required'
        H(4)='by the ESP-r integrated simulation that is underway.'
        if(initonly.eq.0)then
          CALL PHELPD('e2r opening',4,'-',0,0,IER)
        endif
      elseif(module(1:3).eq.'eco')then
        H(1)='This module support a life cycle analysis of the'
        H(2)='of the current model.'
        if(initonly.eq.0)then
          CALL PHELPD('eco opening',2,'-',0,0,IER)
        endif
      elseif(module(1:3).eq.'grd')then
        H(1)='This module supports the definition of multi-dimensional'
        H(2)='grids for building constructions and the ground.'
        if(initonly.eq.0)then
          CALL PHELPD('grd opening',2,'-',0,0,IER)
        endif
      elseif(module(1:3).eq.'ish')then
        H(1) ='This module pre-calculates the shading of external'
        H(2) ='surfaces and the insolation of internal surfaces for'
        H(3) ='a typical day in each month and passes these data to'
        H(4) ='the ESP-r Integrated Simulator.'
        if(initonly.eq.0)then
          CALL PHELPD('ish opening',4,'-',0,0,IER)
        endif
      elseif(module(1:3).eq.'mfs')then
        H(1)='This module is a standalone solver for network flow'
        H(2)='problems (i.e. with no building/ plant conflation).'
        if(initonly.eq.0)then
          CALL PHELPD('mfs opening',2,'-',0,0,IER)
        endif
      elseif(module(1:3).eq.'mld')then
        h(1)='This is the main menu for mould growth results'
        h(2)='analysis. The user can examine the growing'
        h(3)='conditions for a particular mould from the results'
        h(4)='of a building simulation (with combined heat and'
        h(5)='moisture transfer).'  
        h(6)='Statistics - reports the time conditions are'
        h(7)='suitable for the growth of a particular mould.'
        if(initonly.eq.0)then
          CALL PHELPD('mld opening',7,'-',0,0,IER)
        endif
      elseif(module(1:3).eq.'mrt'.or.module(1:3).eq.'vwf')then
        H(1)='MRT attrib & calc: Reads a viewfactor file which'
        H(2)='  includes the definition of MRT sensors within a'
        H(3)='  zone.  The user may calculate viewfactors for the'
        H(4)='  sensors or use previously computed viewfactors to'
        H(5)='  be used for re-evaluation of position dependent'
        H(6)='  MRT with different sets of surface temperatures.'
        H(7)='  '
        H(8)='View factor calc: Calculate viewfactors for all'
        H(9)='  surfaces of the input geometry. '
        H(10)=' '
        H(11)='Merge: Write zone viewfactors to file. '
        H(12)=' '
        H(13)='Calc Parameters: Change number of bubbles on a'
        H(14)='  surface, or number of patches on bubble. '
        H(15)=' '
        H(16)='Display surface info: list out the attributes of'
        H(17)='  the surfaces [depending on trace verbosity] '
        if(initonly.eq.0)then
          CALL PHELPD('vwf opening',17,'-',0,0,IER)
        endif
      elseif(module(1:3).eq.'net'.or.module(1:4).eq.'enet')then
        H(1)='Welcome to the ESP-r graphic network tool. '
        H(2)=' '
        h(3)='You can access a description of the facilities and hints'
        h(4)='for use in the menu ? options. '
        if(initonly.eq.0)then
          CALL PHELPD('enet opening',4,'-',0,0,IER)
        endif
      elseif(module(1:3).eq.'pdb')then
        h(1) =' The MENU commands have the following meaning :'
        h(2) =' '
        h(3) =' HELP    this facility.'
        h(4) =' '
        h(5) =' LIST    chronological, sequential, alphabetical or'
        h(6) ='         user-defined-order list of the database'
        h(7) =' '
        h(8) =' INSERT  allows plant components to be entered into'
        h(9) ='         the database'
        h(10)=' '
        h(11)=' DELETE  allows components to be deleted from'
        h(12)='         the database'
        h(13)=' '
        h(14)=' EDIT    enables component attributes to be'
        h(15)='         individually edited'
        h(16)=' '
        h(17)=' FINISH  terminates the current program run'
        if(initonly.eq.0)then
          CALL PHELPD('pdb opening',17,'-',0,0,IER)
        endif
      elseif(module(1:3).eq.'plt')then
        h(1) =' Menu "ESRUplt" commands are as follows :'
        h(2) =' 1:  initiates the input section of the program'
        h(3) ='     reading all data files and performing all'
        h(4) ='     time independant pre-simulation calculations.'
        h(5) ='     Menu picks "2" and "3" cannot be selected'
        h(6) ='     until the data has been successfully entered.'
        h(7) =' 2:  invokes the simulation section of the program'
        h(8) =' 3:  permits the re-specification of a climate file'
        h(9) =' t:  allows the user to retrieve selected results'
        h(10)='     from a subroutine trace facility. This facility'
        h(11)='     is for use primarily by the experienced user.'
        h(12)=' '
        h(13)=' w:  stops non-fatal warning messages being'
        h(14)='     printed. However, the total number of warning'
        h(15)='     messages will be output in the simulation'
        h(16)='     description prior to commencement of simulation.'
        h(17)=' ?:  invokes this and other descriptive messages'
        h(18)='     as an aid to inexperienced users. These '
        h(19)='     messages are a synopsis of, and should be used'
        h(20)='     in conjunction with, the ESP documentation set'
        h(21)='     "ESP-r A Building and Plant Energy Simulation'
        h(22)='     System".'
        if(initonly.eq.0)then
          CALL PHELPD('plt opening',22,'-',0,0,IER)
        endif
C      elseif(module(1:3).eq.'prj')then
C        H(1)=' '
C        if(initonly.eq.0)then
C          CALL PHELPD('prj opening',50,'-',0,0,IER)
C        endif
      elseif(module(1:3).eq.'pro')then
        H(1)= '     ESP-r Standard Profiles Database '
        H(2)= '   ================================ '
        H(3)= ' '
        H(4)= ' The data requirements are summarised below :- '
        H(5)= ' '
        H(6)= ' for each profile - number of events.'
        H(7)= ' '
        H(8)= ' for each event - start and finish hours'
        H(9)= '                  percentage of first magnitude '
        H(10)='                  percentage of second magnitude '
        H(11)=' '
        H(12)='Currently you may add new profiles at the end of the'
        H(13)='database or, for existing profiles alter event con- '
        H(14)='tents but not the number of events. '
        H(15)=' '
        H(16)='Profiles in this database are currently used to  '
        H(17)='assist in the defintion of zone operations.  The'
        H(18)='standard is found in /usr/esru/esp-r/databases '
        if(initonly.eq.0)then
          CALL PHELPD('pro opening',18,'-',0,0,IER)
        endif
      elseif(module(1:3).eq.'res')then
        H(1)='Results analysis is the simulation results '
        H(2)='recovery facility for ESP-r. '
        H(4)=' '
        H(5)='It always begins with the specification of a results'
        H(6)='library which may contain results for one or more '
        H(7)='simulations performed by the Simulator. '
        H(8)=' '
        H(9) ='The symbols >> in a menu denotes a toggle selection,'
        H(10)='the inclusion of a [?] in a prompt means that a'
        H(11)='default and explanation of the item is available.'
        H(12)=' '
        H(13)='Current menu options are: '
        H(14)=' 1 select result file - select this if the results'
        H(15)='   file needs to be changed. Not used very often'
        H(16)='   as res is usually started with a specified file'
        H(17)='   name.'
        H(18)=' 2 select result set - if more than one assessment'
        H(19)='   is stored in the file (perhaps for a different'
        H(20)='   period) you can re-select using this facility.'
        H(21)=' 3 define output period - allows you to `zoom in`'
        H(22)='   on a portion of an assessment.'
        H(23)=' 4 select zones - allows you to select one or more'
        H(24)='   zones to include in a graph or report.'
        H(25)=' a graphs - includes a range of graph types e.g.'
        H(26)='   variable(s) vs time, frequency binns.'
        H(27)=' c timestep reports - tabular listings of data at'
        H(28)='   each timestep over the selected period.'
        H(29)=' d enquire about - e.g. statistics, energy balances,'
        H(30)='   hours over/under, frequency binn tables'
        H(31)=' e plant results - anything to do with detailed'
        H(32)='   system components and their performance metrics.'
        H(33)=' f indoor environmental quality - comfort and air'
        H(34)='   quality issues'
        H(35)=' g electrical results - if an electical domain is'
        H(36)='   included in the model all reports/graphs are here.'
        H(37)=' h CFD - if there is a CFD domain, 2D and 3D graphs'
        H(38)='   and images at slices can be produced here.'
        if(initonly.eq.0)then
          CALL PHELPD('res opening',38,'-',0,0,IER)
        endif
      elseif(module(1:3).eq.'tdf')then
        continue  ! tdf functionality has been moved into esruprj
      elseif(module(1:3).eq.'vew')then
        H(1)= ' ? HELP     to list command options available'
        H(2)= ' Byebye     to exit'
        H(3)= ' /Input     to TYPE IN numerical values of E F or M '
        H(4)= ' Eye point  to pick the Eye-point with cursor'
        H(5)= ' Foc &Midp  to pick the Focus and Mid-points together'
        H(6)= ' Mid point  to pick the Mid-point alone'
        H(7)= ' Hidden     to remove hidden lines'
        H(8)= ' Dashed     to draw them dashed'
        H(9)= ' Visible    to draw them solid'
        H(10)=' Lens mm    to specify Lens focal length (in mm)'
        H(11)=' *Enlarg %  to specify Enlargement factor'
        H(12)=' Angle vew  to specify viewcone angle (degrees)'
        H(13)='       (negative for biggest picture that will fit)'
        H(14)='  '
        H(15)= ' 4 Orthogr  to select orthogonal projection'
        H(16)= ' 3 Perspec  to select perspective projection'
        H(17)= ' 2 || horz  to select parallel proj. onto horiz. plane'
        H(18)= ' 1 || vert  to select parallel proj. onto vert. plane'
        H(19)= ' Type all   to list all current view parameters'
        H(20)= ' Input all  to type in all view parameters'
        H(21)= ' Screen     to draw the picture'
        H(22)= ' Out file   to file it in a picture-file'
        H(23)= ' New file   to select a new input file'
        H(24)=' +Def vews  to select default views'
        if(initonly.eq.0)then
          CALL PHELPD('vew opening',24,'-',0,0,IER)
        endif
      elseif(module(1:4).eq.'wire')then
        h(1)= 'This interface provides control of the wireframe image'
        h(2)= 'of the model geometry. There are several sections as'
        h(3)= 'follows: '
        h(4)= ' View orientation - allows you to defined the eye point'
        h(5)= '   (where you are looking from), the view point (what'
        h(6)= '   you are looking at). Both of these can be edited by'
        h(7)= '   typing in new values or by using the up/down arrows'
        h(8)= '   to the right of each number.'
        h(9)= '   The view bounds can be toggled from optimum to use'
        h(10)='   a user specified angle of view. Reducing the angle'
        h(11)='   of view will act to zoom into the view.'
        h(12)=' '
        h(13)=' Display information toggles - this toggles off various'
        h(14)='   attributes (e.g. zone names, site grid).'
        h(15)=' '
        h(16)=' Display grid - controls whether a dotted grid is drawn'
        h(17)='   and the density of the grid.'
        h(18)=' '
        h(19)=' Zones to display - one box for each zone in the model'
        h(20)='   (unticking removes that zone from the view)'
        h(21)=' '
        h(22)=' Display options - three sets of options:'
        h(23)='   View type: pick one of these'
        h(24)='   Include: allows you to draw all or a subset of the'
        h(25)='     model geometry e.g. only surfaces facing outside.'
        h(26)='   Hilight: allows you to identify specific attributes'
        h(27)='     of entities by drawing them in a bold line.'
        h(28)=' '
        h(29)='To apply one or more changes click on the Apply button.'
        h(30)=' '
        h(31)='To apply one or more changes and close the wireframe '
        h(32)='control interface click on the OK button. '
        h(33)=' '
        h(34)='Note: a few of the selections are not yet operational. '
        if(initonly.eq.0)then
          CALL PHELPD('wirefram ctl',34,'-',0,0,IER)
        endif
      else
        H(1)='Did not find a matching ESP-r application name.'
        H(2)='So no one has linked this facility yet. '
        if(initonly.eq.0)then
          CALL PHELPD('ask about mismatch',2,'-',0,0,IER)
        endif
      endif

C Recover previous help via dh array.
      if(initonly.eq.0)then
        do 51 i=1,60
          h(i)= dh(i)
   51   continue
      endif
      return
      end

C test subroutinhe to capture help text from external file (to be
C found in the ESP-r distribution lib folder) with file ending
C of *.help
      subroutine gethelptext(sub,topic,nbhelp)
#include "espriou.h"
      
      integer lnblnk  ! function definition

C Parameters.
      character sub*24   ! subroutine name
      character topic*24 ! string (unique) for topic
      integer nbhelp     ! number of help lines found (will be non-zero
                         ! because default message will be created

C Commons.
      integer ifil
      COMMON/FILEP/IFIL
      character h*72
      common/pophelp/h(60)
      integer iuout,iuin,ieout
      COMMON/OUTIN/IUOUT,IUIN,IEOUT

C Local variables.
      logical unixok
      integer ln,lna  ! length of install path, app name app folder
      integer lns,lnt,lnhl ! length of subroutine, topic, helpline
      integer lno          ! length of outstr
      integer K            ! position
      integer trywhich     ! if one then primary help file if two use secondary

      character helpfil*96 ! the help file to load
      dimension folder(4),lnf(4)  ! up to 3 subsequent folders to look in and their lengths
      character folder*8   ! the source folder for the application
      integer lnf
      character WORD*32    ! a token to read
      character outstr*124 ! for help phrases
      character outs*124 ! for messages
      character helpline*72 ! line of help text
      character sq*1        ! single quote
      character fs*1
      logical havesub       ! set true if subroutine matches
      logical havetopic     ! set true if topic matches
      character app*24      ! application name (local)
#include "espinstalldir.h"

C Get to work.
      app = helpinapp      ! set local app
      trywhich = 1         ! assume primary help file
      call isunix(unixok)  ! which folder separator to use
      if(unixok)then
        fs = char(47)
      else
        fs = char(92)
      endif
      ln=lnblnk(instpath)  ! where is ESP-r installed
      lna=lnblnk(app)      ! length of app name
      folder(1)='UNKNOWN'
      folder(2)='UNKNOWN'
      folder(3)='UNKNOWN'
      folder(4)='UNKNOWN'
      if(app(1:3).eq.'aco')then
        folder='esruaco'
        folder(2)='common'
      elseif(app(1:3).eq.'cdb')then
        folder(1)='esrucdb'
        folder(2)='common'
      elseif(app(1:3).eq.'clm')then
        folder(1)='esruclm'
        folder(2)='common'
      elseif(app(1:3).eq.'bps')then
        folder(1)='esrubps'
        folder(2)='common'
      elseif(app(1:3).eq.'dbm')then
        folder(1)='esrudbm'
        folder(2)='common'
      elseif(app(1:3).eq.'dfs')then
        folder(1)='esrudfs'
        folder(2)='common'
      elseif(app(1:3).eq.'dfv')then
        folder(1)='esrudfs'
        folder(2)='common'
      elseif(app(1:3).eq.'e2r')then
        folder(1)='esrue2r'
        folder(2)='common'
      elseif(app(1:3).eq.'eco')then
        folder(1)='esrueco'
        folder(2)='common'
      elseif(app(1:4).eq.'ecnv')then
        folder(1)='esrucnv'
        folder(2)='common'
      elseif(app(1:3).eq.'grd')then
        folder(1)='esrugrd'
        folder(2)='common'
      elseif(app(1:3).eq.'ish')then
        folder(1)='esruish'
        folder(2)='common'
      elseif(app(1:3).eq.'mfs')then
        folder(1)='esrumfs'
        folder(2)='common'
      elseif(app(1:3).eq.'mld')then
        folder(1)='esrumld'
        folder(2)='common'
      elseif(app(1:3).eq.'mrt')then
        folder(1)='esrumrt'
        folder(2)='common'
      elseif(app(1:3).eq.'net')then
        folder(1)='esrunet'
        folder(2)='common'
      elseif(app(1:3).eq.'pdb')then
        folder(1)='esrupdb'
        folder(2)='common'
      elseif(app(1:3).eq.'prj')then
        folder(1)='esruprj'
        folder(2)='esrumfs'
        folder(3)='esrudfs'
        folder(4)='common'
      elseif(app(1:3).eq.'res')then
        folder(1)='esrures'
        folder(2)='esrumfs'
        folder(3)='esrubps'
        folder(4)='common'
      elseif(app(1:3).eq.'vew')then
        folder(1)='esruvew'
        folder(2)='common'
      endif
      lnf(1)=lnblnk(folder(1))   ! length of folder name
      lnf(2)=lnblnk(folder(2))   ! length of folder name
      lnf(3)=lnblnk(folder(3))   ! length of folder name
      lnf(4)=lnblnk(folder(4))   ! length of folder name
      lns=lnblnk(sub)      ! length of subroutine
      lnt=lnblnk(topic)    ! length of topic
      sq=char(39)   ! assign single quote

C Assume that the installer has placed the help text file in the
C distribution lib folder.
  42  if(trywhich.eq.1)then
        write(helpfil,'(6a)') instpath(1:ln),fs,'lib',fs,
     &    folder(1)(1:lnf(1)),'.help'
C        call edisp(iuout,helpfil)   ! debug
      elseif(trywhich.ge.2)then
        write(helpfil,'(6a)') instpath(1:ln),fs,'lib',fs,
     &    folder(trywhich)(1:lnf(trywhich)),'.help'
C        call edisp(iuout,helpfil)   ! debug
      endif
      if(folder(1)(1:7).eq.'esrubps')then
        ihfil=ifil+8              ! set file unit number for bps
      elseif(folder(1)(1:7).eq.'esruprj')then
        ihfil=ifil+11              ! set file unit number for prj
      elseif(folder(1)(1:7).eq.'esrures')then
        ihfil=ifil+8              ! set file unit number for res
      elseif(folder(1)(1:7).eq.'esrumrt')then
        ihfil=ifil+8              ! set file unit number for mrt
      else
        ihfil=ifil+2              ! set file unit number
      endif
      CALL ERPFREE(ihfil,ISTAT) ! free it
      CALL EFOPSEQ(ihfil,helpfil,1,IER)
      if(ier.ne.0)then
        call edisp(iuout,'Could not find help text file.')
        call edisp(iuout,helpfil)
        nbhelp=-1
        return
      endif

      nbhelp=0  ! set initial assumption
      havesub= .false.
      havetopic= .false.

C Check the header to see if this is a help text file.
      read(ihfil,'(a)',IOSTAT=IOS,ERR=1) outstr   ! scan header
      if(ios.ne.0) goto 1   ! jump if EOF
      K=0
      CALL EGETW(outstr,K,WORD,'W','*help_text',IER)
      if(WORD(1:10).ne.'*help_text')then
        call edisp(iuout,'Not a help text file.')
        call edisp(iuout,helpfil)
        goto 1
      endif
      CALL EGETW(outstr,K,WORD,'W','application',IER)

C If we are looking in the primary help file check the app name.
C Subsequent help texts will not have a matching app line.
      if(trywhich.eq.1)then
        if(WORD(1:3).eq.'dfv'.and.app(1:3).eq.'dfs')then
          continue  ! dfv uses dfs help text
        elseif(WORD(1:3).eq.'dfs'.and.app(1:3).eq.'dfv')then
          continue  ! dfv uses dfs help text
        elseif(WORD(1:lna).ne.app(1:lna))then
          write(outs,'(3a)') WORD(1:lna),
     &      ' is not a help text file for ',app(1:lna)
          call edisp(iuout,outs)
          goto 1
        endif
      endif

C If we reached this point then we can go looking for a
C subroutine and topic match on each new line.
      scan: do
        read(ihfil,'(a)',IOSTAT=IOS,ERR=1) outstr   ! scan line
        if(ios.ne.0)  exit scan   ! jump if EOF
        lno=lnblnk(outstr)
        if(lno.le.1) cycle  ! jump if a blank line

        if(outstr(1:1).eq.'*')then  ! we had a key phrase
          K=0
          CALL EGETW(outstr,K,WORD,'W','*item',IER)
          if(WORD(1:5).eq.'*item')then
            CALL EGETW(outstr,K,WORD,'W','subroutine',IER)
            if(WORD(1:lns).eq.sub(1:lns))then

C We have a matching subroutine name remember this and if
C there are still more characters in the line then it
C probably is a topic.
              havesub = .true.  ! set true subroutine matches
              if(K.le.lno)then
                CALL EGETW(outstr,K,WORD,'W','topic',IER)
                if(WORD(1:lnt).eq.topic(1:lnt))then
                  havetopic = .true.     ! set true topic matches
                  cycle   ! keep scanning
                else
                  cycle   ! keep scanning
                endif
              else
                cycle   ! keep scanning
              endif
            else
              cycle   ! keep scanning
            endif
          elseif(WORD(1:8).eq.'*enditem')then

C Reached end of item exit scan, if this was the correct sub and topic then
C set number of help lines.
            if(havesub.and.havetopic)then
              CALL ERPFREE(ihfil,ISTAT) ! free it
              return
            else
              cycle  ! keep scanning
            endif
          elseif(WORD(1:6).eq.'*topic')then

C In case there are some old *topic lines in the help f
            if(lno.le.7) cycle  ! jump if no other token
            K=7
            CALL EGETW(outstr,K,WORD,'W','topic',IER)
            if(WORD(1:lnt).eq.topic(1:lnt))then
              havetopic = .true.     ! set true topic matches
              cycle   ! keep scanning
            else
              cycle   ! keep scanning
            endif
          elseif(WORD(1:11).eq.'*subroutine')then
            cycle   ! keep scanning
          endif

        elseif(outstr(1:1).eq.sq)then

C We have a line that starts with a quote. If we have reached the
C correct topic then read it otherwise scan another line.
          if(havesub.and.havetopic)then
            K=0
            call EGETDQ(outstr,K,helpline,'W','phrase',IER)
            lnhl=lnblnk(helpline)
C Debug (to text feedback).
C            call edisp(iuout,helpline)
            nbhelp=nbhelp+1
            if(lnhl.le.1)then
              h(nbhelp)='  '
            else
              write(h(nbhelp),'(a)') helpline(1:lnhl)
            endif
            cycle   ! keep scanning
          else
            cycle   ! keep scanning
          endif
        else
 
C Does not match expectations.
          cycle   ! keep scanning
        endif
      end do scan

C The end of file reached. If there has not been a match try other
C likely help files for the current application. If that fails
C make up a standard message.
  1   close(ihfil)
      if(havesub.and.havetopic)then
        return
      else

C Check if there is an alternative help text file. If not give up.
        if(trywhich.eq.1)then
          if(folder(2)(1:7).eq.'UNKNOWN')then
C            call edisp(iuout,'alt folder was UNKNOWN')
            h(1)='No alt contextual help text was found for the topic'
            write(h(2),'(5a)') topic(1:lnt),' in ',sub(1:lns),
     &        ' of application ',app(1:lna)
            nbhelp=2
            return
          endif
          trywhich=trywhich+1
          goto 42  ! try scaning the alternative file
        elseif(trywhich.eq.2)then
          if(folder(3)(1:7).eq.'UNKNOWN')then
C            call edisp(iuout,'2nd alt folder was UNKNOWN')
            h(1)='No alt contextual help text was found for the topic'
            write(h(2),'(5a)') topic(1:lnt),' in ',sub(1:lns),
     &        ' of application ',app(1:lna)
            nbhelp=2
            return
          endif
          trywhich=trywhich+1
          goto 42  ! try scaning the alternative file
        elseif(trywhich.eq.3)then
          if(folder(4)(1:7).eq.'UNKNOWN')then
C            call edisp(iuout,'3rd alt folder was UNKNOWN')
            h(1)='No alt contextual help text was found for the topic'
            write(h(2),'(5a)') topic(1:lnt),' in ',sub(1:lns),
     &        ' of application ',app(1:lna)
            nbhelp=2
            return
          endif
          trywhich=trywhich+1
          goto 42  ! try scaning the alternative file
        elseif(trywhich.eq.4)then
          continue
        endif
        h(1)='No help found for this topic!'
        write(h(2),'(5a)') topic(1:lnt),' in ',sub(1:lns),
     &    ' of application ',app(1:lna)
        nbhelp=2
        return
      endif

      end

C subroutine to generate line of help text for an integer dialog
C in terms of default and current value. 
      subroutine helpcurrentint(idef,icur,ermsg,nbhelp,newnbhelp,ier)
      
      integer lnblnk  ! function definition

C Parameters.
      integer idef         ! default value for this number
      integer icur         ! current value
      CHARACTER*(*) ERMSG  ! error message from dialog
      integer nbhelp       ! original nb of help lines
      integer newnbhelp    ! number of help lines after additional line(s)
      integer ier          ! zero is ok

C Commons.
      character h*72
      common/pophelp/h(60)
      integer iuout,iuin,ieout
      COMMON/OUTIN/IUOUT,IUIN,IEOUT

C Local variables.
      integer ln           ! length of error message
      integer iw2,iw1      ! length of two integers
      character outs*124 ! for messages
      CHARACTER STR1*10,STR2*10

C Get to work.
      CALL INTSTR(icur,STR2,IW2,IER)
      CALL INTSTR(idef,STR1,IW1,IER)
      IF(nbhelp.GT.0)THEN
        newnbhelp=nbhelp+1
        WRITE(H(nbhelp+1),65,iostat=ios,err=1)STR2(1:IW2),STR1(1:IW1)
      ELSE
        LN=max(1,LNBLNK(ERMSG))
        newnbhelp=2
        WRITE(H(1),64,iostat=ios,err=1)ERMSG(1:LN)
   64   FORMAT('No help available for ',a,'!')
        WRITE(H(2),65,iostat=ios,err=1)STR2(1:IW2),STR1(1:IW1)
   65   FORMAT('Current value is (',a,'); Default value is (',a,').')
      ENDIF

      return

   1  if(IOS.eq.2)then
        write(outs,*) 'helpcurretint: permission error composing help.'
        call edisp(iuout,outs)
      else
        write(outs,*) 'helpcurretint: error composing help.'
        call edisp(iuout,outs)
      endif
      return
      end


C subroutine to generate line of help text for an real dialog
C in terms of default and current value. 
      subroutine helpcurrentreal(def,cur,ermsg,nbhelp,newnbhelp,ier)
      
      integer lnblnk  ! function definition

C Parameters.
      real def         ! default value for this number
      real cur         ! current value
      CHARACTER*(*) ERMSG  ! error message from dialog
      integer nbhelp       ! original nb of help lines
      integer newnbhelp    ! number of help lines after additional line(s)
      integer ier          ! zero is ok

C Commons.
      character h*72
      common/pophelp/h(60)
      integer iuout,iuin,ieout
      COMMON/OUTIN/IUOUT,IUIN,IEOUT

C Local variables.
      integer ln           ! length of error message
      integer iw2,iw1      ! length of two integers
      character outs*124 ! for messages
      CHARACTER STR1*16,STR2*16

C Get to work.
      CALL REL16STR(def,STR1,IW1,IER)
      CALL REL16STR(cur,STR2,IW2,IER)
      IF(nbhelp.GT.0)THEN
        newnbhelp=nbhelp+1
        WRITE(H(nbhelp+1),65,iostat=ios,err=1)STR2(1:IW2),STR1(1:IW1)
      ELSE
        LN=max(1,LNBLNK(ERMSG))
        newnbhelp=2
        WRITE(H(1),64,iostat=ios,err=1)ERMSG(1:LN)
   64   FORMAT('No help available for ',a,'!')
        WRITE(H(2),65,iostat=ios,err=1)STR2(1:IW2),STR1(1:IW1)
   65   FORMAT('Current value is (',a,'); Default value is (',a,').')
      ENDIF

      return

   1  if(IOS.eq.2)then
        write(outs,*) 
     &    'helpcurrentreal: permission error composing help.'
        call edisp(iuout,outs)
      else
        write(outs,*) 'helpcurrentreal: error composing help.'
        call edisp(iuout,outs)
      endif
      return
      end

C subroutine to generate line of help text for an string dialog
C if no lines of help have been defined in the calling code. 
      subroutine helpwithblank(ermsg,nbhelp,newnbhelp,ier)
      
      integer lnblnk  ! function definition

C Parameters.
      CHARACTER*(*) ERMSG  ! error message from dialog
      integer nbhelp       ! original nb of help lines
      integer newnbhelp    ! number of help lines after additional line(s)
      integer ier          ! zero is ok

C Commons.
      character h*72
      common/pophelp/h(60)
      integer iuout,iuin,ieout
      COMMON/OUTIN/IUOUT,IUIN,IEOUT

C Local variables.
      integer ln           ! length of error message
      character outs*124 ! for messages

C Get to work.
      IF(nbhelp.GT.0)THEN
        newnbhelp=nbhelp  ! nothing to do
      ELSE
        LN=max(1,LNBLNK(ERMSG))
        newnbhelp=1
        WRITE(H(1),64,iostat=ios,err=1)ERMSG(1:LN)
   64   FORMAT('No help available for ',a,'!')
      ENDIF

      return

   1  if(IOS.eq.2)then
        write(outs,*) 
     &    'helpwithblank: permission error composing help.'
        call edisp(iuout,outs)
        ier=2
      else
        write(outs,*) 'helpwithblank: error composing help.'
        call edisp(iuout,outs)
        ier=1
      endif
      return
      end
