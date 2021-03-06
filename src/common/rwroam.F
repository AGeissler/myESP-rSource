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

C You should have received a copy of the GNU General Public
C License along with ESP-r. If not, write to the Free
C Software Foundation, Inc., 59 Temple Place, Suite 330,
C Boston, MA 02111-1307 USA.

C This file includes:
C  RROAM   reads roaming file
C  WROAM   writes roaming file

C ********************* RROAM
C Subroutine to read roaming file. This file supports changing location
C in terms of latitude and longitude, and orientation of a model
C temporally. Furthermore day types can be associated with each period
C as defined in this file and can later be used to link casual gains as
C defined in the operations files.

      SUBROUTINE RROAM
      IMPLICIT NONE

C espriou.h provides currentfile.
#include "espriou.h"
#include "roam.h"
      
      integer lnblnk  ! function definition

      integer iuout,iuin,ieout
      COMMON/OUTIN/IUOUT,IUIN,IEOUT
      COMMON/FILEP/IFIL
      CHARACTER OUTSTR*124,outs*124
      INTEGER IDEBUG,IER,ILOC,IUF,IDEL,IPDR,IFIL,ND,K,ISTAT
      LOGICAL GOON

C Read in file and initialise variables
      IPDR=0
      IDEL=0
      ILOC=0
      IDEBUG=0
      IUF=IFIL+1
      GOON=.TRUE.
      CALL EFOPSEQ(IUF,LROAM,1,IER)
      write(currentfile,'(A)') LROAM(1:lnblnk(LROAM))
      CALL STRIPC(IUF,OUTSTR,99,ND,1,'Roaming file',IER)

C Read version number (currently 0.1)
      IF(OUTSTR(1:8).EQ.'*version')THEN      
        K=8
        CALL EGETWR(OUTSTR,K,VERS,0.01,10.,'W','version #',IER)
      ELSE
        IER=1
      ENDIF

C Read location information
      CALL STRIPC(IUF,OUTSTR,99,ND,1,'Roaming file',IER)
      IF(OUTSTR(1:9).NE.'*location')IER=1

      DO WHILE(OUTSTR(1:13).NE.'*end_location')
        IF(GOON)THEN
          CALL STRIPC(IUF,OUTSTR,99,ND,1,'Roaming file',IER)
          GOON=.FALSE.
        ENDIF
        ILOC=ILOC+1
        K=0
        CALL EGETWR(OUTSTR,K,XLOCLT(ILOC),-90.,90.,'W','latitude',IER)
        CALL EGETWR(OUTSTR,K,XLOCLG(ILOC),-180.,180.,'W','longtude',IER)
        CALL EGETWI(OUTSTR,K,IDEL,0,0,'-','location #',IER)
        CALL EGETRM(OUTSTR,K,LOCRM(ILOC),'W','location name',IER)
        NLOC=ILOC
        CALL STRIPC(IUF,OUTSTR,99,ND,1,'Roaming file',IER)
      END DO

C Read start time of each period, location and day type
      GOON=.TRUE.
      CALL STRIPC(IUF,OUTSTR,99,ND,1,'Roaming file',IER)
      IF(OUTSTR(1:6).NE.'*start')IER=1
      DO WHILE(OUTSTR(1:13).NE.'*end_temporal')
        IF(GOON)THEN
          CALL STRIPC(IUF,OUTSTR,99,ND,1,'Roaming file',IER)
          GOON=.FALSE.
        ENDIF
        IPDR=IPDR+1
        K=0
        CALL EGETWR(OUTSTR,K,XLHRS(IPDR),0.,0.,'-','start hour',IER)
        CALL EGETWI(OUTSTR,K,LDYS(IPDR),0,0,'-','start hour',IER)
        CALL EGETWI(OUTSTR,K,LMTS(IPDR),0,0,'-','start hour',IER)
        CALL EGETWI(OUTSTR,K,LORNT(IPDR),0,0,'-','start hour',IER)
        CALL EGETWI(OUTSTR,K,LDTP(IPDR),0,0,'-','start hour',IER)
        CALL EGETWI(OUTSTR,K,LLOCT(IPDR),0,0,'-','start hour',IER)
        NPDR=IPDR
        CALL STRIPC(IUF,OUTSTR,99,ND,1,'Roaming file',IER)
      END DO

C Debug...
      IF(IDEBUG.EQ.1)THEN
        WRITE(outs,*)'version number = ',VERS
        call edisp(iuout,outs)
        WRITE(outs,*)'Location Latitude Longitude'
        call edisp(iuout,outs)
        DO 1000 ILOC=1,NLOC
          WRITE(outs,*)LOCRM(ILOC),XLOCLT(ILOC),XLOCLG(ILOC)
          call edisp(iuout,outs)
 1000   CONTINUE
        WRITE(outs,*)
     &    'start hr, dy, mth, orientation, day type, location'
        call edisp(iuout,outs)
        DO 1100 IPDR=1,NPDR
          WRITE(outs,*)XLHRS(IPDR),LDYS(IPDR),LMTS(IPDR),LORNT(IPDR),
     &    LDTP(IPDR),LLOCT(IPDR)
          call edisp(iuout,outs)
 1100   CONTINUE
      ENDIF

C Error handling (still to be done)
      IF(IER.NE.0)THEN
        write(outs,*)'some error has occured while reading',LROAM
        call edisp(iuout,outs)
        write(outs,*)'returning ...'
        call edisp(iuout,outs)
      ENDIF
      CALL ERPFREE(IUF,ISTAT)
      RETURN
      END

C ********************* WROAM
C Subroutine to write roaming file template. The intention is to write
C actual data using a text editor
      SUBROUTINE WROAM
      IMPLICIT NONE

C espriou.h provides currentfile.
#include "espriou.h"
#include "roam.h"

      COMMON/FILEP/IFIL
      INTEGER IUF,IFIL,IOS,IER
      integer lnblnk  ! function type

C Set output channel
      IUF=IFIL+1
      write(currentfile,'(A)',IOSTAT=IOS,ERR=1) LROAM(1:lnblnk(LROAM))
      CALL EFOPSEQ(IUF,LROAM,3,IER)
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)'*version 0.1'
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)
     &'*location information # '
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)
     &'# '
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)
     &'# '
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)
     &'# '
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)
     &'# '
      WRITE(IUF,'(A,1X,I4,1X,A)',IOSTAT=IOS,ERR=1)
     &'# , ',MLOC,' '
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)'0.0 -0.0 1 AAA'
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)'0.0 -0.0 2 BBB'
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)
     &'0.0 -0.0 3 CCC'
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)
     &'0.0 -0.0 4 DDD'
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)'0.0 -0.0 5 EEE'
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)
     &'0.0 -0.0 6 FFF'
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)
     &'0.0 -0.0 7 GGG'
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)'0.0 -0.0 8 HHH'
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)'*end_location information'
      WRITE(IUF,'(2A)',IOSTAT=IOS,ERR=1)'*start hour, day, month, ',
     &'orientation, day type, '
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)
     &'# '
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)
     &'# '
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)
     &'# '
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)
     &'# '
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)
     &'# '
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)
     &'# '
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)
     &'# '
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)
     &'# '
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)
     &'# '
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)
     &'# '
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)
     &'# '
      WRITE(IUF,'(A,1X,I4,1X,A)',IOSTAT=IOS,ERR=1)
     &'# , ',MPDR,' '
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)
     &'# Finally you may delete all lines begining with # (optional)'
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)'1 1 3 0 1 3'
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)'18 1 3 0 1 0'
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)'8 2 3 0 1 4'
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)'16 2 3 0 1 0'
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)'10 3 3 0 1 5'
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)'19 3 3 0 1 0'
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)'7 5 3 0 1 1'
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)'*end_temporal information'
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)
     &'##################### IMPORTANT #######################'
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)
     &'# '
      WRITE(IUF,'(A)',IOSTAT=IOS,ERR=1)
     &'# '

      RETURN

C Error messages.
 1    if(IOS.eq.2)then
        call usrmsg(' No permission to write roam file!',
     &            ' returning to menu...','W')
        ier=1
      else
        call usrmsg(' Roam file transfer error !',
     &            ' returning to menu...','W')
        ier=1
      endif
      return
      
      END
