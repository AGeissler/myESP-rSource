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

C You should have received a copy of the GNU General Public
C License along with ESP-r. If not, write to the Free
C Software Foundation, Inc., 59 Temple Place, Suite 330,
C Boston, MA 02111-1307 USA.

C This file contains subroutines for imposing momentum and turbulence 
C boundary conditions on blockages when k-e model is active and thermal 
C boundary conditions for all turbulence models
C The subroutines are named BLKBND* with the final letters representing 
C the equations being set up/ solved.

C ********************* BLKBNDU *********************
C BLKBNDU - set boundary conditions for U velocity for blockages.
      SUBROUTINE BLKBNDU(IBLK)
#include "building.h"
#include "cfd.h"

      COMMON/ICFNOD/ICFD,ICP
      COMMON/GEOM/XP(ntcelx),YP(ntcely),ZP(ntcelz),
     1            DXEP(ntcelx),DXPW(ntcelx),DYNP(ntcely),DYPS(ntcely),
     2            DZHP(ntcelz),DZPL(ntcelz),
     3            SEW(ntcelx),SNS(ntcely),SHL(ntcelz),
     4            XU(ntcelx),YV(ntcely),ZW(ntcelz)
      COMMON/GEOM3/AreaHLP(ntcelx,ntcely),AreaHLU(ntcelx,ntcely),
     &             AreaHLV(ntcelx,ntcely),AreaEWP(ntcely,ntcelz),
     &             AreaEWV(ntcely,ntcelz),AreaEWW(ntcely,ntcelz),
     &             AreaNSP(ntcelx,ntcelz),AreaNSU(ntcelx,ntcelz),
     &             AreaNSW(ntcelx,ntcelz)
      COMMON/ALL/NI,NJ,NK,NIM1,NJM1,NKM1,NIM2,NJM2,NKM2
      COMMON/VARf/Uf(ntcelx,ntcely,ntcelz),Vf(ntcelx,ntcely,ntcelz),
     1            Wf(ntcelx,ntcely,ntcelz),
     2            P(ntcelx,ntcely,ntcelz),PP(ntcelx,ntcely,ntcelz),
     3            TEf(ntcelx,ntcely,ntcelz),EDf(ntcelx,ntcely,ntcelz)
      COMMON/WALLF/TAUEV(ntcelx,ntcely,ntcelz),
     1          TAUWV(ntcelx,ntcely,ntcelz),TAUEW(ntcelx,ntcely,ntcelz),
     2          TAUWW(ntcelx,ntcely,ntcelz),TAUNU(ntcelx,ntcely,ntcelz),
     3          TAUSU(ntcelx,ntcely,ntcelz),TAUNW(ntcelx,ntcely,ntcelz),
     4          TAUSW(ntcelx,ntcely,ntcelz),TAUHU(ntcelx,ntcely,ntcelz),
     5          TAULU(ntcelx,ntcely,ntcelz),TAUHV(ntcelx,ntcely,ntcelz),
     6          TAULV(ntcelx,ntcely,ntcelz)
      COMMON/COEF/AP(ntcelx,ntcely,ntcelz),AE(ntcelx,ntcely,ntcelz),
     1            AW(ntcelx,ntcely,ntcelz),AN(ntcelx,ntcely,ntcelz),
     2            AS(ntcelx,ntcely,ntcelz),AH(ntcelx,ntcely,ntcelz),
     3            AL(ntcelx,ntcely,ntcelz),SU(ntcelx,ntcely,ntcelz),
     4            SP(ntcelx,ntcely,ntcelz)

      common/KEYVOLS/NVOL(MNZ),IVOLF(MNVLS,MNZ),IVCELLS(MNVLS,MNZ,2),
     &               JVCELLS(MNVLS,MNZ,2),KVCELLS(MNVLS,MNZ,2)

C Set values for all faces.
      do 40 IFACE=3,6
        IST=max(3,IVCELLS(IBLK,ICFD,1))
        IFN=IVCELLS(IBLK,ICFD,2)
        JST=JVCELLS(IBLK,ICFD,1)
        JFN=JVCELLS(IBLK,ICFD,2)
        KST=KVCELLS(IBLK,ICFD,1)
        KFN=KVCELLS(IBLK,ICFD,2)

C Values are set for cell adjacent to blockage.  Adjust values of I,J,K 
C to adjacent cell in all directions and set start and end cell to same 
C cell for the current face. 
        if (IFACE.eq.3) then
          JST=JVCELLS(IBLK,ICFD,2)+1
          JFN=JVCELLS(IBLK,ICFD,2)+1
        elseif (IFACE.eq.4) then
          JST=JVCELLS(IBLK,ICFD,1)-1
          JFN=JVCELLS(IBLK,ICFD,1)-1
        elseif (IFACE.eq.5) then
          KST=KVCELLS(IBLK,ICFD,2)+1
          KFN=KVCELLS(IBLK,ICFD,2)+1
        elseif (IFACE.eq.6) then
          KST=KVCELLS(IBLK,ICFD,1)-1
          KFN=KVCELLS(IBLK,ICFD,1)-1
        endif

C Skip if outside domain.
        if (IST.eq.1.or.IFN.eq.NI.or.JST.eq.1.or.JFN.eq.NJ.or.
     &      KST.eq.1.or.KFN.eq.NK) goto 40

C Loop through whole blockage.
        do 10 II=IST,IFN
          do 20 JJ=JST,JFN
            do 30 KK=KST,KFN

C Calculate distance from wall to next to wall grid point ($DELTA y sub p$).
C Then set source terms in adjacent cells. 
C TAU??: wall shear stress (+ive value means stress operates in axis direction)
C SP: source term, units are force/velocity
C A?=0. : disconnect cell from blocked cell.
C Blockage on South face of adjacent cell.
C Blockage on South face of adjacent cell.
              if (IFACE.eq.3) then
                Dface=YP(JJ)-YV(JJ)
                TMULT=TRBSRCU(II,JJ,KK,Dface)
                TAUSU(II,JJ,KK)=-TMULT*Uf(II,JJ,KK)
                SP(II,JJ,KK)=SP(II,JJ,KK)-TMULT*AreaNSU(II,KK)
                AS(II,JJ,KK)=0.

C Blockage on North face of adjacent cell.
              elseif (IFACE.eq.4) then
                Dface=YV(JJ+1)-YP(JJ)
                TMULT=TRBSRCU(II,JJ,KK,Dface)
                TAUNU(II,JJ,KK)=-TMULT*Uf(II,JJ,KK)
                SP(II,JJ,KK)=SP(II,JJ,KK)-TMULT*AreaNSU(II,KK)
                AN(II,JJ,KK)=0.

C Blockage on Low face of adjacent cell.
              elseif (IFACE.eq.5) then
                Dface=ZP(KK)-ZW(KK)
                TMULT=TRBSRCU(II,JJ,KK,Dface)
                TAULU(II,JJ,KK)=-TMULT*Uf(II,JJ,KK)
                SP(II,JJ,KK)=SP(II,JJ,KK)-TMULT*AreaHLU(II,JJ)
                AL(II,JJ,KK)=0.

C Blockage on High face of adjacent cell.
              elseif (IFACE.eq.6) then
                Dface=ZW(KK+1)-ZP(KK)
                TMULT=TRBSRCU(II,JJ,KK,Dface)
                TAUHU(II,JJ,KK)=-TMULT*Uf(II,JJ,KK)
                SP(II,JJ,KK)=SP(II,JJ,KK)-TMULT*AreaHLU(II,JJ)
                AH(II,JJ,KK)=0.
              endif
 30         continue
 20       continue
 10     continue
 40   continue

      RETURN
      END



C ********************* BLKBNDV *********************
C BLKBNDV - set boundary conditions for V velocity for blockages.
      SUBROUTINE BLKBNDV(IBLK)
#include "building.h"
#include "cfd.h"

      COMMON/ICFNOD/ICFD,ICP
      COMMON/GEOM/XP(ntcelx),YP(ntcely),ZP(ntcelz),
     1            DXEP(ntcelx),DXPW(ntcelx),DYNP(ntcely),DYPS(ntcely),
     2            DZHP(ntcelz),DZPL(ntcelz),
     3            SEW(ntcelx),SNS(ntcely),SHL(ntcelz),
     4            XU(ntcelx),YV(ntcely),ZW(ntcelz)
      COMMON/GEOM3/AreaHLP(ntcelx,ntcely),AreaHLU(ntcelx,ntcely),
     &             AreaHLV(ntcelx,ntcely),AreaEWP(ntcely,ntcelz),
     &             AreaEWV(ntcely,ntcelz),AreaEWW(ntcely,ntcelz),
     &             AreaNSP(ntcelx,ntcelz),AreaNSU(ntcelx,ntcelz),
     &             AreaNSW(ntcelx,ntcelz)
      COMMON/ALL/NI,NJ,NK,NIM1,NJM1,NKM1,NIM2,NJM2,NKM2
      COMMON/VARf/Uf(ntcelx,ntcely,ntcelz),Vf(ntcelx,ntcely,ntcelz),
     1            Wf(ntcelx,ntcely,ntcelz),
     2            P(ntcelx,ntcely,ntcelz),PP(ntcelx,ntcely,ntcelz),
     3            TEf(ntcelx,ntcely,ntcelz),EDf(ntcelx,ntcely,ntcelz)
      COMMON/WALLF/TAUEV(ntcelx,ntcely,ntcelz),
     1          TAUWV(ntcelx,ntcely,ntcelz),TAUEW(ntcelx,ntcely,ntcelz),
     2          TAUWW(ntcelx,ntcely,ntcelz),TAUNU(ntcelx,ntcely,ntcelz),
     3          TAUSU(ntcelx,ntcely,ntcelz),TAUNW(ntcelx,ntcely,ntcelz),
     4          TAUSW(ntcelx,ntcely,ntcelz),TAUHU(ntcelx,ntcely,ntcelz),
     5          TAULU(ntcelx,ntcely,ntcelz),TAUHV(ntcelx,ntcely,ntcelz),
     6          TAULV(ntcelx,ntcely,ntcelz)
      COMMON/COEF/AP(ntcelx,ntcely,ntcelz),AE(ntcelx,ntcely,ntcelz),
     1            AW(ntcelx,ntcely,ntcelz),AN(ntcelx,ntcely,ntcelz),
     2            AS(ntcelx,ntcely,ntcelz),AH(ntcelx,ntcely,ntcelz),
     3            AL(ntcelx,ntcely,ntcelz),SU(ntcelx,ntcely,ntcelz),
     4            SP(ntcelx,ntcely,ntcelz)

      common/KEYVOLS/NVOL(MNZ),IVOLF(MNVLS,MNZ),IVCELLS(MNVLS,MNZ,2),
     &               JVCELLS(MNVLS,MNZ,2),KVCELLS(MNVLS,MNZ,2)

C Set values for all faces.
      do 40 IFACE=1,6
        if (IFACE.lt.3.or.IFACE.gt.4) then
          IST=IVCELLS(IBLK,ICFD,1)
          IFN=IVCELLS(IBLK,ICFD,2)
          JST=max(3,JVCELLS(IBLK,ICFD,1))
          JFN=JVCELLS(IBLK,ICFD,2)
          KST=KVCELLS(IBLK,ICFD,1)
          KFN=KVCELLS(IBLK,ICFD,2)

C Values are set for cell adjacent to blockage.  Adjust values of I,J,K 
C to adjacent cell in all directions and set start and end cell to same 
C cell for the current face. 
          if (IFACE.eq.1) then
            IST=IVCELLS(IBLK,ICFD,2)+1
            IFN=IVCELLS(IBLK,ICFD,2)+1
          elseif (IFACE.eq.2) then
            IST=IVCELLS(IBLK,ICFD,1)-1
            IFN=IVCELLS(IBLK,ICFD,1)-1
          elseif (IFACE.eq.5) then
            KST=KVCELLS(IBLK,ICFD,2)+1
            KFN=KVCELLS(IBLK,ICFD,2)+1
          elseif (IFACE.eq.6) then
            KST=KVCELLS(IBLK,ICFD,1)-1
            KFN=KVCELLS(IBLK,ICFD,1)-1
          endif

C Skip if outside domain.
        if (IST.eq.1.or.IFN.eq.NI.or.JST.eq.1.or.JFN.eq.NJ.or.
     &      KST.eq.1.or.KFN.eq.NK) goto 40

C Loop through whole blockage.
          do 10 II=IST,IFN
            do 20 JJ=JST,JFN
              do 30 KK=KST,KFN

C Calculate distance from wall to next to wall grid point ($DELTA y sub p$).
C Then set source terms in adjacent cells. 
C TAU??: wall shear stress (+ive value means stress operates in axis direction)
C SP: source term, units are force/velocity
C A?=0. : disconnect cell from blocked cell.
C Blockage on West face of adjacent cell.
                if (IFACE.eq.1) then
                  Dface=XP(II)-XU(II)
                  TMULT=TRBSRCV(II,JJ,KK,Dface)
                  TAUWV(II,JJ,KK)=-TMULT*Vf(II,JJ,KK)
                  SP(II,JJ,KK)=SP(II,JJ,KK)-TMULT*AreaEWV(JJ,KK)
                  AW(II,JJ,KK)=0.

C Blockage on East face of adjacent cell.
                elseif (IFACE.eq.2) then
                  Dface=XU(II+1)-XP(II)
                  TMULT=TRBSRCV(II,JJ,KK,Dface)
                  TAUEV(II,JJ,KK)=-TMULT*Vf(II,JJ,KK)
                  SP(II,JJ,KK)=SP(II,JJ,KK)-TMULT*AreaEWV(JJ,KK)
                  AE(II,JJ,KK)=0.

C Blockage on Low face of adjacent cell.
                elseif (IFACE.eq.5) then
                  Dface=ZP(KK)-ZW(KK)
                  TMULT=TRBSRCV(II,JJ,KK,Dface)
                  TAULV(II,JJ,KK)=-TMULT*Vf(II,JJ,KK)
                  SP(II,JJ,KK)=SP(II,JJ,KK)-TMULT*AreaHLV(II,JJ)
                  AL(II,JJ,KK)=0.

C Blockage on High face of adjacent cell.
                elseif (IFACE.eq.6) then
                  Dface=ZW(KK+1)-ZP(KK)
                  TMULT=TRBSRCV(II,JJ,KK,Dface)
                  TAUHV(II,JJ,KK)=-TMULT*Vf(II,JJ,KK)
                  SP(II,JJ,KK)=SP(II,JJ,KK)-TMULT*AreaHLV(II,JJ)
                  AH(II,JJ,KK)=0.
                endif
 30           continue
 20         continue
 10       continue
        endif
 40   continue

      RETURN
      END


C ********************* BLKBNDW *********************
C BLKBNDW - set boundary conditions for W velocity for blockages.
      SUBROUTINE BLKBNDW(IBLK)
#include "building.h"
#include "cfd.h"

      COMMON/ICFNOD/ICFD,ICP
      COMMON/GEOM/XP(ntcelx),YP(ntcely),ZP(ntcelz),
     1            DXEP(ntcelx),DXPW(ntcelx),DYNP(ntcely),DYPS(ntcely),
     2            DZHP(ntcelz),DZPL(ntcelz),
     3            SEW(ntcelx),SNS(ntcely),SHL(ntcelz),
     4            XU(ntcelx),YV(ntcely),ZW(ntcelz)
      COMMON/GEOM3/AreaHLP(ntcelx,ntcely),AreaHLU(ntcelx,ntcely),
     &             AreaHLV(ntcelx,ntcely),AreaEWP(ntcely,ntcelz),
     &             AreaEWV(ntcely,ntcelz),AreaEWW(ntcely,ntcelz),
     &             AreaNSP(ntcelx,ntcelz),AreaNSU(ntcelx,ntcelz),
     &             AreaNSW(ntcelx,ntcelz)
      COMMON/ALL/NI,NJ,NK,NIM1,NJM1,NKM1,NIM2,NJM2,NKM2
      COMMON/VARf/Uf(ntcelx,ntcely,ntcelz),Vf(ntcelx,ntcely,ntcelz),
     1            Wf(ntcelx,ntcely,ntcelz),
     2            P(ntcelx,ntcely,ntcelz),PP(ntcelx,ntcely,ntcelz),
     3            TEf(ntcelx,ntcely,ntcelz),EDf(ntcelx,ntcely,ntcelz)
      COMMON/WALLF/TAUEV(ntcelx,ntcely,ntcelz),
     1          TAUWV(ntcelx,ntcely,ntcelz),TAUEW(ntcelx,ntcely,ntcelz),
     2          TAUWW(ntcelx,ntcely,ntcelz),TAUNU(ntcelx,ntcely,ntcelz),
     3          TAUSU(ntcelx,ntcely,ntcelz),TAUNW(ntcelx,ntcely,ntcelz),
     4          TAUSW(ntcelx,ntcely,ntcelz),TAUHU(ntcelx,ntcely,ntcelz),
     5          TAULU(ntcelx,ntcely,ntcelz),TAUHV(ntcelx,ntcely,ntcelz),
     6          TAULV(ntcelx,ntcely,ntcelz)
      COMMON/COEF/AP(ntcelx,ntcely,ntcelz),AE(ntcelx,ntcely,ntcelz),
     1            AW(ntcelx,ntcely,ntcelz),AN(ntcelx,ntcely,ntcelz),
     2            AS(ntcelx,ntcely,ntcelz),AH(ntcelx,ntcely,ntcelz),
     3            AL(ntcelx,ntcely,ntcelz),SU(ntcelx,ntcely,ntcelz),
     4            SP(ntcelx,ntcely,ntcelz)

      common/KEYVOLS/NVOL(MNZ),IVOLF(MNVLS,MNZ),IVCELLS(MNVLS,MNZ,2),
     &               JVCELLS(MNVLS,MNZ,2),KVCELLS(MNVLS,MNZ,2)

C Set values for all faces.
      do 40 IFACE=1,4
        IST=IVCELLS(IBLK,ICFD,1)
        IFN=IVCELLS(IBLK,ICFD,2)
        JST=JVCELLS(IBLK,ICFD,1)
        JFN=JVCELLS(IBLK,ICFD,2)
        KST=max(3,KVCELLS(IBLK,ICFD,1))
        KFN=KVCELLS(IBLK,ICFD,2)

C Values are set for cell adjacent to blockage.  Adjust values of I,J,K 
C to adjacent cell in all directions and set start and end cell to same 
C cell for the current face. 
        if (IFACE.eq.1) then
          IST=IVCELLS(IBLK,ICFD,2)+1
          IFN=IVCELLS(IBLK,ICFD,2)+1
        elseif (IFACE.eq.2) then
          IST=IVCELLS(IBLK,ICFD,1)-1
          IFN=IVCELLS(IBLK,ICFD,1)-1
        elseif (IFACE.eq.3) then
          JST=JVCELLS(IBLK,ICFD,2)+1
          JFN=JVCELLS(IBLK,ICFD,2)+1
        elseif (IFACE.eq.4) then
          JST=JVCELLS(IBLK,ICFD,1)-1
          JFN=JVCELLS(IBLK,ICFD,1)-1
        endif

C Skip if outside domain.
        if (IST.eq.1.or.IFN.eq.NI.or.JST.eq.1.or.JFN.eq.NJ.or.
     &      KST.eq.1.or.KFN.eq.NK) goto 40

C Loop through whole blockage.
        do 10 II=IST,IFN
          do 20 JJ=JST,JFN
            do 30 KK=KST,KFN

C Calculate distance from wall to next to wall grid point ($DELTA y sub p$).
C Then set source terms in adjacent cells. 
C TAU??: wall shear stress (+ive value means stress operates in axis direction)
C SP: source term, units are force/velocity
C A?=0. : disconnect cell from blocked cell.
C Blockage on West face of adjacent cell.
              if (IFACE.eq.1) then
                Dface=XP(II)-XU(II)
                TMULT=TRBSRCW(II,JJ,KK,Dface)
                TAUWW(II,JJ,KK)=-TMULT*Wf(II,JJ,KK)
                SP(II,JJ,KK)=SP(II,JJ,KK)-TMULT*AreaEWW(JJ,KK)
                AW(II,JJ,KK)=0.

C Blockage on East face of adjacent cell.
              elseif (IFACE.eq.2) then
                Dface=XU(II+1)-XP(II)
                TMULT=TRBSRCW(II,JJ,KK,Dface)
                TAUEW(II,JJ,KK)=-TMULT*Wf(II,JJ,KK)
                SP(II,JJ,KK)=SP(II,JJ,KK)-TMULT*AreaEWW(JJ,KK)
                AE(II,JJ,KK)=0.

C Blockage on South face of adjacent cell.
              elseif (IFACE.eq.3) then
                Dface=YP(JJ)-YV(JJ)
                TMULT=TRBSRCW(II,JJ,KK,Dface)
                TAUSW(II,JJ,KK)=-TMULT*Wf(II,JJ,KK)
                SP(II,JJ,KK)=SP(II,JJ,KK)-TMULT*AreaNSW(II,KK)
                AS(II,JJ,KK)=0.

C Blockage on North face of adjacent cell.
              elseif (IFACE.eq.4) then
                Dface=YV(JJ+1)-YP(JJ)
                TMULT=TRBSRCW(II,JJ,KK,Dface)
                TAUNW(II,JJ,KK)=-TMULT*Wf(II,JJ,KK)
                SP(II,JJ,KK)=SP(II,JJ,KK)-TMULT*AreaNSW(II,KK)
                AN(II,JJ,KK)=0.
              endif
 30         continue
 20       continue
 10     continue
 40   continue

      RETURN
      END


C ********************* BLKBNDKE *********************
C BLKBNDKE - set boundary conditions for turbulent kinetic 
C   energy for blockages.
      SUBROUTINE BLKBNDKE(IBLK)
#include "building.h"
#include "cfd.h"

      COMMON/ICFNOD/ICFD,ICP
      COMMON/GEOM/XP(ntcelx),YP(ntcely),ZP(ntcelz),
     1            DXEP(ntcelx),DXPW(ntcelx),DYNP(ntcely),DYPS(ntcely),
     2            DZHP(ntcelz),DZPL(ntcelz),
     3            SEW(ntcelx),SNS(ntcely),SHL(ntcelz),
     4            XU(ntcelx),YV(ntcely),ZW(ntcelz)
      COMMON/GEOM2/VolP(ntcelx,ntcely,ntcelz),
     &             VolU(ntcelx,ntcely,ntcelz),
     &             VolV(ntcelx,ntcely,ntcelz),
     &             VolW(ntcelx,ntcely,ntcelz)
      COMMON/ALL/NI,NJ,NK,NIM1,NJM1,NKM1,NIM2,NJM2,NKM2
      COMMON/FLUPRf/URFVIS,VISCOS,PRANDT,SH,
     1            DENf(ntcelx,ntcely,ntcelz),VIS(ntcelx,ntcely,ntcelz),
     2            BETA(ntcelx,ntcely,ntcelz)
      COMMON/VARf/Uf(ntcelx,ntcely,ntcelz),Vf(ntcelx,ntcely,ntcelz),
     1            Wf(ntcelx,ntcely,ntcelz),
     2            P(ntcelx,ntcely,ntcelz),PP(ntcelx,ntcely,ntcelz),
     3            TEf(ntcelx,ntcely,ntcelz),EDf(ntcelx,ntcely,ntcelz)
      COMMON/TURB/GEN(ntcelx,ntcely,ntcelz),CD,CMU,C1,C2,C3,CAPPA,ELOG,
     &            TURBIN,ALAMDA,PRTE,PRED
      COMMON/WALLF/TAUEV(ntcelx,ntcely,ntcelz),
     1          TAUWV(ntcelx,ntcely,ntcelz),TAUEW(ntcelx,ntcely,ntcelz),
     2          TAUWW(ntcelx,ntcely,ntcelz),TAUNU(ntcelx,ntcely,ntcelz),
     3          TAUSU(ntcelx,ntcely,ntcelz),TAUNW(ntcelx,ntcely,ntcelz),
     4          TAUSW(ntcelx,ntcely,ntcelz),TAUHU(ntcelx,ntcely,ntcelz),
     5          TAULU(ntcelx,ntcely,ntcelz),TAUHV(ntcelx,ntcely,ntcelz),
     6          TAULV(ntcelx,ntcely,ntcelz)
      COMMON/COEF/AP(ntcelx,ntcely,ntcelz),AE(ntcelx,ntcely,ntcelz),
     1            AW(ntcelx,ntcely,ntcelz),AN(ntcelx,ntcely,ntcelz),
     2            AS(ntcelx,ntcely,ntcelz),AH(ntcelx,ntcely,ntcelz),
     3            AL(ntcelx,ntcely,ntcelz),SU(ntcelx,ntcely,ntcelz),
     4            SP(ntcelx,ntcely,ntcelz)
      COMMON/SUSP/SUKD(ntcelx,ntcely,ntcelz),SPKD(ntcelx,ntcely,ntcelz)

      COMMON/CONST/GREAT,small,GRAV

      common/KEYVOLS/NVOL(MNZ),IVOLF(MNVLS,MNZ),IVCELLS(MNVLS,MNZ,2),
     &               JVCELLS(MNVLS,MNZ,2),KVCELLS(MNVLS,MNZ,2)

C Set values for all faces.
      do 40 IFACE=1,6
        IST=IVCELLS(IBLK,ICFD,1)
        IFN=IVCELLS(IBLK,ICFD,2)
        JST=JVCELLS(IBLK,ICFD,1)
        JFN=JVCELLS(IBLK,ICFD,2)
        KST=KVCELLS(IBLK,ICFD,1)
        KFN=KVCELLS(IBLK,ICFD,2)

C Values are set for cell adjacent to blockage.  Adjust values of I,J,K 
C to adjacent cell in all directions and set start and end cell to same 
C cell for the current face. 
        if (IFACE.eq.1) then
          IST=IVCELLS(IBLK,ICFD,2)+1
          IFN=IVCELLS(IBLK,ICFD,2)+1
        elseif (IFACE.eq.2) then
          IST=IVCELLS(IBLK,ICFD,1)-1
          IFN=IVCELLS(IBLK,ICFD,1)-1
        elseif (IFACE.eq.3) then
          JST=JVCELLS(IBLK,ICFD,2)+1
          JFN=JVCELLS(IBLK,ICFD,2)+1
        elseif (IFACE.eq.4) then
          JST=JVCELLS(IBLK,ICFD,1)-1
          JFN=JVCELLS(IBLK,ICFD,1)-1
        elseif (IFACE.eq.5) then
          KST=KVCELLS(IBLK,ICFD,2)+1
          KFN=KVCELLS(IBLK,ICFD,2)+1
        elseif (IFACE.eq.6) then
          KST=KVCELLS(IBLK,ICFD,1)-1
          KFN=KVCELLS(IBLK,ICFD,1)-1
        endif

C Skip if outside domain.
        if (IST.eq.1.or.IFN.eq.NI.or.JST.eq.1.or.JFN.eq.NJ.or.
     &      KST.eq.1.or.KFN.eq.NK) goto 40

C Loop through whole blockage.
        do 10 II=IST,IFN
          do 20 JJ=JST,JFN
            do 30 KK=KST,KFN

C Calculate distance from wall to next to wall grid point.
C Blockage on West face of adjacent cell.
              if (IFACE.eq.1) then
                Dface=XP(II)-XU(II)
                GENCOU=(ABS(TAUWV(II,JJ+1,KK)*Vf(II,JJ+1,KK))+
     &                  ABS(TAUWV(II,JJ,KK)*Vf(II,JJ,KK))+
     &                  ABS(TAUWW(II,JJ,KK+1)*Wf(II,JJ,KK+1))+
     &                  ABS(TAUWW(II,JJ,KK)*Wf(II,JJ,KK)))/(2.0*Dface)
                AW(II,JJ,KK)=0.

C Blockage on East face of adjacent cell.
              elseif (IFACE.eq.2) then
                Dface=XU(II+1)-XP(II)
                GENCOU=(ABS(TAUEV(II,JJ+1,KK)*Vf(II,JJ+1,KK))+
     &                  ABS(TAUEV(II,JJ,KK)*Vf(II,JJ,KK))+
     &                  ABS(TAUEW(II,JJ,KK+1)*Wf(II,JJ,KK+1))+
     &                  ABS(TAUEW(II,JJ,KK)*Wf(II,JJ,KK)))/(2.0*Dface)
                AE(II,JJ,KK)=0.

C Blockage on South face of adjacent cell.
              elseif (IFACE.eq.3) then
                Dface=YP(JJ)-YV(JJ)
                GENCOU=(ABS(TAUSU(II+1,JJ,KK)*Uf(II+1,JJ,KK))+
     &                  ABS(TAUSU(II,JJ,KK)*Uf(II,JJ,KK))+
     &                  ABS(TAUSW(II,JJ,KK+1)*Wf(II,JJ,KK+1))+
     &                  ABS(TAUSW(II,JJ,KK)*Wf(II,JJ,KK)))/(2.0*Dface)
                AS(II,JJ,KK)=0.

C Blockage on North face of adjacent cell.
              elseif (IFACE.eq.4) then
                Dface=YV(JJ+1)-YP(JJ)
                GENCOU=(ABS(TAUNU(II+1,JJ,KK)*Uf(II+1,JJ,KK))+
     &                  ABS(TAUNU(II,JJ,KK)*Uf(II,JJ,KK))+
     &                  ABS(TAUNW(II,JJ,KK+1)*Wf(II,JJ,KK+1))+
     &                  ABS(TAUNW(II,JJ,KK)*Wf(II,JJ,KK)))/(2.0*Dface)
                AN(II,JJ,KK)=0.

C Blockage on Low face of adjacent cell.
              elseif (IFACE.eq.5) then
                Dface=ZP(KK)-ZW(KK)
                GENCOU=(ABS(TAULU(II+1,JJ,KK)*Uf(II+1,JJ,KK))+
     &                  ABS(TAULU(II,JJ,KK)*Uf(II,JJ,KK))+
     &                  ABS(TAULV(II,JJ+1,KK)*Vf(II,JJ+1,KK))+
     &                  ABS(TAULV(II,JJ,KK)*Vf(II,JJ,KK)))/(2.0*Dface)
                AL(II,JJ,KK)=0.

C Blockage on High face of adjacent cell.
              elseif (IFACE.eq.6) then
                Dface=ZW(KK+1)-ZP(KK)
                GENCOU=(ABS(TAUHU(II+1,JJ,KK)*Uf(II+1,JJ,KK))+
     &                  ABS(TAUHU(II,JJ,KK)*Uf(II,JJ,KK))+
     &                  ABS(TAUHV(II,JJ+1,KK)*Vf(II,JJ+1,KK))+
     &                  ABS(TAUHV(II,JJ,KK)*Vf(II,JJ,KK)))/(2.0*Dface)
                AH(II,JJ,KK)=0.
              endif

C Calculate and set source terms in adjacent cells if inside domain.
              if (II.gt.1.and.JJ.gt.1.and.KK.gt.1.and.
     &            II.lt.NI.and.JJ.lt.NJ.and.KK.lt.NK) then
                SQRTK=sqrt(abs(TEf(II,JJ,KK)))
                Tvis=VIS(II,JJ,KK)-VISCOS+small
                Dplus=DENf(II,JJ,KK)*SQRTK*(CMU**0.25)*Dface/VISCOS
                call VELDERIV(II,JJ,KK,DUDY,DUDZ,DVDX,DVDZ,DWDX,DWDY)
                GENRES=GEN(II,JJ,KK)-Tvis*((DUDY+DVDX)**2+
     &                                   (DUDZ+DWDX)**2+(DVDZ+DWDY)**2)
                GEN(II,JJ,KK)=GENRES+GENCOU
                if (Dplus.gt.11.63) then
                  DITERM=CD*DENf(II,JJ,KK)*(CMU**.75)*SQRTK*
     &                                   ALOG(ELOG*Dplus)/(CAPPA*Dface)
                else
                  DITERM=CD*DENf(II,JJ,KK)*(CMU**.75)*SQRTK*Dplus/Dface
                endif
                SU(II,JJ,KK)=GEN(II,JJ,KK)*VolP(II,JJ,KK)+SUKD(II,JJ,KK)
                SP(II,JJ,KK)=-DITERM*VolP(II,JJ,KK)+SPKD(II,JJ,KK)
              endif
 30         continue
 20       continue
 10     continue
 40   continue

      RETURN
      END


C ********************* BLKBNDED *********************
C BLKBNDED - set boundary conditions for the dissipation of turbulent kinetic 
C   energy due to solid face for blockages.
      SUBROUTINE BLKBNDED(IBLK)
#include "building.h"
#include "cfd.h"

      COMMON/ICFNOD/ICFD,ICP
      COMMON/GEOM/XP(ntcelx),YP(ntcely),ZP(ntcelz),
     1            DXEP(ntcelx),DXPW(ntcelx),DYNP(ntcely),DYPS(ntcely),
     2            DZHP(ntcelz),DZPL(ntcelz),
     3            SEW(ntcelx),SNS(ntcely),SHL(ntcelz),
     4            XU(ntcelx),YV(ntcely),ZW(ntcelz)
      COMMON/ALL/NI,NJ,NK,NIM1,NJM1,NKM1,NIM2,NJM2,NKM2
      COMMON/VARf/Uf(ntcelx,ntcely,ntcelz),Vf(ntcelx,ntcely,ntcelz),
     1            Wf(ntcelx,ntcely,ntcelz),
     2            P(ntcelx,ntcely,ntcelz),PP(ntcelx,ntcely,ntcelz),
     3            TEf(ntcelx,ntcely,ntcelz),EDf(ntcelx,ntcely,ntcelz)
      COMMON/TURB/GEN(ntcelx,ntcely,ntcelz),CD,CMU,C1,C2,C3,CAPPA,ELOG,
     &            TURBIN,ALAMDA,PRTE,PRED
      COMMON/COEF/AP(ntcelx,ntcely,ntcelz),AE(ntcelx,ntcely,ntcelz),
     1            AW(ntcelx,ntcely,ntcelz),AN(ntcelx,ntcely,ntcelz),
     2            AS(ntcelx,ntcely,ntcelz),AH(ntcelx,ntcely,ntcelz),
     3            AL(ntcelx,ntcely,ntcelz),SU(ntcelx,ntcely,ntcelz),
     4            SP(ntcelx,ntcely,ntcelz)

      COMMON/CONST/GREAT,small,GRAV

      common/KEYVOLS/NVOL(MNZ),IVOLF(MNVLS,MNZ),IVCELLS(MNVLS,MNZ,2),
     &               JVCELLS(MNVLS,MNZ,2),KVCELLS(MNVLS,MNZ,2)

C Set values for all faces.
      do 40 IFACE=1,6
        IST=IVCELLS(IBLK,ICFD,1)
        IFN=IVCELLS(IBLK,ICFD,2)
        JST=JVCELLS(IBLK,ICFD,1)
        JFN=JVCELLS(IBLK,ICFD,2)
        KST=KVCELLS(IBLK,ICFD,1)
        KFN=KVCELLS(IBLK,ICFD,2)

C Values are set for cell adjacent to blockage.  Adjust values of I,J,K 
C to adjacent cell in all directions and set start and end cell to same 
C cell for the current face. 
        if (IFACE.eq.1) then
          IST=IVCELLS(IBLK,ICFD,2)+1
          IFN=IVCELLS(IBLK,ICFD,2)+1
        elseif (IFACE.eq.2) then
          IST=IVCELLS(IBLK,ICFD,1)-1
          IFN=IVCELLS(IBLK,ICFD,1)-1
        elseif (IFACE.eq.3) then
          JST=JVCELLS(IBLK,ICFD,2)+1
          JFN=JVCELLS(IBLK,ICFD,2)+1
        elseif (IFACE.eq.4) then
          JST=JVCELLS(IBLK,ICFD,1)-1
          JFN=JVCELLS(IBLK,ICFD,1)-1
        elseif (IFACE.eq.5) then
          KST=KVCELLS(IBLK,ICFD,2)+1
          KFN=KVCELLS(IBLK,ICFD,2)+1
        elseif (IFACE.eq.6) then
          KST=KVCELLS(IBLK,ICFD,1)-1
          KFN=KVCELLS(IBLK,ICFD,1)-1
        endif

C Skip if outside domain.
        if (IST.eq.1.or.IFN.eq.NI.or.JST.eq.1.or.JFN.eq.NJ.or.
     &      KST.eq.1.or.KFN.eq.NK) goto 40

C Loop through whole blockage.
        do 10 II=IST,IFN
          do 20 JJ=JST,JFN
            do 30 KK=KST,KFN

C Calculate distance from 
C next to wall grid point and disconnect cell from blockage (A?=0.)
C Blockage on West face of adjacent cell.
              if (IFACE.eq.1) then
                Dface=XP(II)-XU(II)
                AW(II,JJ,KK)=0.

C Blockage on East face of adjacent cell.
              elseif (IFACE.eq.2) then
                Dface=XU(II+1)-XP(II)
                AE(II,JJ,KK)=0.

C Blockage on South face of adjacent cell.
              elseif (IFACE.eq.3) then
                Dface=YP(JJ)-YV(JJ)
                AS(II,JJ,KK)=0.

C Blockage on North face of adjacent cell.
              elseif (IFACE.eq.4) then
                Dface=YV(JJ+1)-YP(JJ)
                AN(II,JJ,KK)=0.

C Blockage on Low face of adjacent cell.
              elseif (IFACE.eq.5) then
                Dface=ZP(KK)-ZW(KK)
                AL(II,JJ,KK)=0.

C Blockage on High face of adjacent cell.
              elseif (IFACE.eq.6) then
                Dface=ZW(KK+1)-ZP(KK)
                AH(II,JJ,KK)=0.
              endif

C Calculate and set source terms.
              TERM=(CMU**.75)/(CAPPA*Dface)
              SU(II,JJ,KK)=GREAT*TERM*ABS(TEf(II,JJ,KK))**1.5
              SP(II,JJ,KK)=-GREAT
 30         continue
 20       continue
 10     continue
 40   continue

      RETURN
      END


C ********************* BLKBNDT *********************
C BLKBNDT - set boundary conditions for the temperature/ heat flux 
C   due to solid face for blockages.
C   VOLHeat is the heat flux in [W/m2], for the  whole blockage, 
C   prescribed by the user.
      SUBROUTINE BLKBNDT(IBLK)
#include "building.h"
#include "cfd.h"

      COMMON/ICFNOD/ICFD,ICP
      COMMON/GEOM3/AreaHLP(ntcelx,ntcely),AreaHLU(ntcelx,ntcely),
     &             AreaHLV(ntcelx,ntcely),AreaEWP(ntcely,ntcelz),
     &             AreaEWV(ntcely,ntcelz),AreaEWW(ntcely,ntcelz),
     &             AreaNSP(ntcelx,ntcelz),AreaNSU(ntcelx,ntcelz),
     &             AreaNSW(ntcelx,ntcelz)
      COMMON/ALL/NI,NJ,NK,NIM1,NJM1,NKM1,NIM2,NJM2,NKM2
      COMMON/FLUPRf/URFVIS,VISCOS,PRANDT,SH,
     1            DENf(ntcelx,ntcely,ntcelz),VIS(ntcelx,ntcely,ntcelz),
     2            BETA(ntcelx,ntcely,ntcelz)
      COMMON/COEF/AP(ntcelx,ntcely,ntcelz),AE(ntcelx,ntcely,ntcelz),
     1            AW(ntcelx,ntcely,ntcelz),AN(ntcelx,ntcely,ntcelz),
     2            AS(ntcelx,ntcely,ntcelz),AH(ntcelx,ntcely,ntcelz),
     3            AL(ntcelx,ntcely,ntcelz),SU(ntcelx,ntcely,ntcelz),
     4            SP(ntcelx,ntcely,ntcelz)
      COMMON/ITERAT/NITER
      common/KEYVOLS/NVOL(MNZ),IVOLF(MNVLS,MNZ),IVCELLS(MNVLS,MNZ,2),
     &               JVCELLS(MNVLS,MNZ,2),KVCELLS(MNVLS,MNZ,2)
      common/KEYVDAT/IVTYPE(MNVLS,MNZ),VOLTemp(MNVLS,MNZ),
     &          VOLHeat(MNVLS,MNZ),IVConfl(MNVLS,MNZ),VOLHum(MNVLS,MNZ),
     &          VOLCO2(MNVLS,MNZ),VOLVel(MNVLS,MNZ),VOLDir(MNVLS,MNZ,2),
     &          VOLArea(MNVLS,MNZ),VOLPres(MNVLS,MNZ),
     &          VOLPol(MCTM,MNVLS,MNZ)
      COMMON/HUMIDT/HUMf(ntcelx,ntcely,ntcelz) 
      COMMON/TEMPf/Tf(ntcelx,ntcely,ntcelz),GAMH(ntcelx,ntcely,ntcelz),
     1             RESORT,NSWPT,URFT,FSDTT,PRANDL,PFUN
      COMMON/VARf/Uf(ntcelx,ntcely,ntcelz),Vf(ntcelx,ntcely,ntcelz),
     1            Wf(ntcelx,ntcely,ntcelz),
     2            P(ntcelx,ntcely,ntcelz),PP(ntcelx,ntcely,ntcelz),
     3            TEf(ntcelx,ntcely,ntcelz),EDf(ntcelx,ntcely,ntcelz)

!!! These common passes the environmnetal parameters arround each body segment to the physiology model 
      REAL headflux1,trunkflux1,larmflux1,rarmflux1,llegflux1,
     &     rlegflux1,H2Oexp1,H2Oswt1,Qoccsens1,Qocclat1,COflux,
     &     Tsko_av1,Tco_av1,Tsk_av1
      COMMON/ocflux1/headflux1,trunkflux1,larmflux1,rarmflux1,llegflux1,
     &        rlegflux1,H2Oexp1,H2Oswt1,Qoccsens1,Qocclat1,COflux,
     &        Tsko_av1, Tco_av1, Tsk_av1
      REAL tavhead1,tavtrunk1,tavarml1,tavarmr1,tavlegl1,tavlegr1
      REAL hrhead,hrtrunk,hrarml,hrarmr,hrlegl,hrlegr
      integer cfdcheck, msegflag
      COMMON/segtemp1/tavhead1,tavtrunk1,tavarml1,tavarmr1,tavlegl1,
     &       tavlegr1,hrhead,hrtrunk,hrarml,hrarmr,hrlegl,hrlegr,
     &       cfdcheck, msegflag,HMETA2
      REAL headflux2,trunkflux2,larmflux2,rarmflux2,llegflux2,
     &     rlegflux2,H2Oexp2,H2Oswt2,Qoccsens2,Qocclat2,
     &     Tsko_av2,Tco_av2,Tsk_av2
      COMMON/ocflux2/headflux2,trunkflux2,larmflux2,rarmflux2,llegflux2,
     &        rlegflux2,H2Oexp2,H2Oswt2,Qoccsens2,Qocclat2,
     &        Tsko_av2, Tco_av2, Tsk_av2
      REAL tavhead2,tavtrunk2,tavarml2,tavarmr2,tavlegl2,tavlegr2
      COMMON/segtemp2/tavhead2,tavtrunk2,tavarml2,tavarmr2,tavlegl2,
     &       tavlegr2
      COMMON/ocfluxo/H2Oexp,H2Oswt,COnflux
      common/KEYVOLN/VOLNAME(MNVLS,MNZ),VCsurf(MNVLS,MNZ)
      character VOLNAME*12,VCsurf*12

      REAL SAT1X,SAT2X,SAT3X,SAT4X,FAT1X,FAT2X,FAT3X,FAT4X
      REAL SAT1Y,SAT2Y,SAT3Y,SAT4Y,FAT1Y,FAT2Y,FAT3Y,FAT4Y
      REAL SARH1,SARH2,SARH3,SARH4,FARH1,FARH2,FARH3,FARH4
      REAL Svel1,Svel2,Svel3,Svel4,Fvel1,Fvel2,Fvel3,Fvel4
      REAL Velhead1,Veltrunk1,Velarml1,Velarmr1,Vellegl1,Vellegr1
      SAT1X=0.0
      SAT2X=0.0
      SAT3X=0.0
      SAT4X=0.0
      FAT1X=0.0
      FAT2X=0.0
      FAT3X=0.0
      FAT4X=0.0
      SAT1Y=0.0
      SAT2Y=0.0
      SAT3Y=0.0
      SAT4Y=0.0
      FAT1Y=0.0
      FAT2Y=0.0
      FAT3Y=0.0
      FAT4Y=0.0
C Set values for all faces.
      do 40 IFACE=1,6
        IST=IVCELLS(IBLK,ICFD,1)
        IFN=IVCELLS(IBLK,ICFD,2)
        JST=JVCELLS(IBLK,ICFD,1)
        JFN=JVCELLS(IBLK,ICFD,2)
        KST=KVCELLS(IBLK,ICFD,1)
        KFN=KVCELLS(IBLK,ICFD,2)

c      write(*,*)'UP_istfn,jstfn,kstfn=',IST,IFN,JST,JFN,KST,KFN

C Values are set for cell adjacent to blockage.  Adjust values of I,J,K 
C to adjacent cell in all directions and set start and end cell to same 
C cell for the current face. 
        if (IFACE.eq.1) then
          IST=IVCELLS(IBLK,ICFD,2)+1
          IFN=IVCELLS(IBLK,ICFD,2)+1
        elseif (IFACE.eq.2) then
          IST=IVCELLS(IBLK,ICFD,1)-1
          IFN=IVCELLS(IBLK,ICFD,1)-1
        elseif (IFACE.eq.3) then
          JST=JVCELLS(IBLK,ICFD,2)+1
          JFN=JVCELLS(IBLK,ICFD,2)+1
        elseif (IFACE.eq.4) then
          JST=JVCELLS(IBLK,ICFD,1)-1
          JFN=JVCELLS(IBLK,ICFD,1)-1
        elseif (IFACE.eq.5) then
          KST=KVCELLS(IBLK,ICFD,2)+1
          KFN=KVCELLS(IBLK,ICFD,2)+1
        elseif (IFACE.eq.6) then
          KST=KVCELLS(IBLK,ICFD,1)-1
          KFN=KVCELLS(IBLK,ICFD,1)-1
        endif

c      write(*,*)'LO_istfn,jstfn,kstfn=',IST,IFN,JST,JFN,KST,KFN

C Skip if outside domain.
        if (IST.eq.1.or.IFN.eq.NI.or.JST.eq.1.or.JFN.eq.NJ.or.
     &      KST.eq.1.or.KFN.eq.NK) goto 40

C Loop through whole blockage.
        do 10 II=IST,IFN
          do 20 JJ=JST,JFN
            do 30 KK=KST,KFN

c      write(*,*)'L2=i.,j.,k.,IFACE',IST,IFN,JST,JFN,KST,KFN,IFACE 

C Calculate distance from wall to next to wall grid point.
C Get surface area of face and disconnect cell from blockage (A?=0.)
C Blockage on West face of adjacent cell.
              if (IFACE.eq.1) then
                AREA=AreaEWP(JJ,KK)
                AW(II,JJ,KK)=0.

C Blockage on East face of adjacent cell.
              elseif (IFACE.eq.2) then
                AREA=AreaEWP(JJ,KK)
                AE(II,JJ,KK)=0.

C Blockage on South face of adjacent cell.
              elseif (IFACE.eq.3) then
                AREA=AreaNSP(II,KK)
                AS(II,JJ,KK)=0.

C Blockage on North face of adjacent cell.
              elseif (IFACE.eq.4) then
                AREA=AreaNSP(II,KK)
                AN(II,JJ,KK)=0.

C Blockage on Low face of adjacent cell.
              elseif (IFACE.eq.5) then
                AREA=AreaHLP(II,JJ)
                AL(II,JJ,KK)=0.

C Blockage on High face of adjacent cell.
              elseif (IFACE.eq.6) then
                AREA=AreaHLP(II,JJ)
                AH(II,JJ,KK)=0.
              endif

C Set combined convective/ diffusive flux terms.
             if (VOLNAME(IBLK,ICFD).eq.'head1')then
                SU(II,JJ,KK)=SU(II,JJ,KK)+headflux1*AREA/SH
              Elseif(VOLNAME(IBLK,ICFD).eq.'head2')then
                SU(II,JJ,KK)=SU(II,JJ,KK)+headflux2*AREA/SH
              Elseif (VOLNAME(IBLK,ICFD).eq.'trunk1')then
                SU(II,JJ,KK)=SU(II,JJ,KK)+trunkflux1*AREA/SH
              Elseif (VOLNAME(IBLK,ICFD).eq.'trunk2')then
                SU(II,JJ,KK)=SU(II,JJ,KK)+trunkflux2*AREA/SH
              Elseif(VOLNAME(IBLK,ICFD).eq.'arml1')then
                SU(II,JJ,KK)=SU(II,JJ,KK)+larmflux1*AREA/SH
              Elseif(VOLNAME(IBLK,ICFD).eq.'arml2')then
                SU(II,JJ,KK)=SU(II,JJ,KK)+larmflux2*AREA/SH
              Elseif(VOLNAME(IBLK,ICFD).eq.'armr1')then
                SU(II,JJ,KK)=SU(II,JJ,KK)+rarmflux1*AREA/SH
              Elseif(VOLNAME(IBLK,ICFD).eq.'armr2')then
                SU(II,JJ,KK)=SU(II,JJ,KK)+rarmflux2*AREA/SH
              Elseif(VOLNAME(IBLK,ICFD).eq.'legl1')then
                SU(II,JJ,KK)=SU(II,JJ,KK)+llegflux1*AREA/SH
              Elseif(VOLNAME(IBLK,ICFD).eq.'legl2')then
                SU(II,JJ,KK)=SU(II,JJ,KK)+llegflux2*AREA/SH
              Elseif(VOLNAME(IBLK,ICFD).eq.'legr1')then
                SU(II,JJ,KK)=SU(II,JJ,KK)+rlegflux1*AREA/SH
              Elseif(VOLNAME(IBLK,ICFD).eq.'legr2')then
                SU(II,JJ,KK)=SU(II,JJ,KK)+rlegflux2*AREA/SH
              ElseIf (VOLNAME(IBLK,ICFD).eq.'Occ'.or.
     &             VOLNAME(IBLK,ICFD).eq.'Occ1'.or.
     &             VOLNAME(IBLK,ICFD).eq.'Occ2'.or.
     &             VOLNAME(IBLK,ICFD).eq.'Occ3'.or.
     &             VOLNAME(IBLK,ICFD).eq.'Occ4'.or.
     &             VOLNAME(IBLK,ICFD).eq.'Occ5'.or.
     &             VOLNAME(IBLK,ICFD).eq.'Occ6'.or.
     &             VOLNAME(IBLK,ICFD).eq.'Occ7')then
                SU(II,JJ,KK)=SU(II,JJ,KK)+COnflux*AREA/SH
              Else
               SU(II,JJ,KK)=SU(II,JJ,KK)+VOLHeat(IBLK,ICFD)*AREA/SH
               Endif
c       write(*,*)'SU_contrib=',VOLHeat(IBLK,ICFD) 
c      pause'heat blockage'
 30         continue
 20       continue
 10     continue
 40   continue

      If ((VOLNAME(IBLK,ICFD).eq.'head1').or.
     &     (VOLNAME(IBLK,ICFD).eq.'trunk1').or.
     &     (VOLNAME(IBLK,ICFD).eq.'armr1').or.
     &     (VOLNAME(IBLK,ICFD).eq.'arml1').or.
     &     (VOLNAME(IBLK,ICFD).eq.'legr1').or.
     &     (VOLNAME(IBLK,ICFD).eq.'legl1').or.
     &     (VOLNAME(IBLK,ICFD).eq.'head2').or.
     &     (VOLNAME(IBLK,ICFD).eq.'trunk2').or.
     &     (VOLNAME(IBLK,ICFD).eq.'armr2').or.
     &     (VOLNAME(IBLK,ICFD).eq.'arml2').or.
     &     (VOLNAME(IBLK,ICFD).eq.'legr2').or.
     &     (VOLNAME(IBLK,ICFD).eq.'legl2'))then

        IST=IVCELLS(IBLK,ICFD,1)
        IFN=IVCELLS(IBLK,ICFD,2)
        JST=JVCELLS(IBLK,ICFD,1)
        JFN=JVCELLS(IBLK,ICFD,2)
        KST=KVCELLS(IBLK,ICFD,1)
        KFN=KVCELLS(IBLK,ICFD,2)

C Loop through whole blockage.
        do II=IST,IFN
          do JJ=JST,JFN
            do KK=KST,KFN


C<< Calculating the average aire temperature and relative humidity of blockage nearest air cell
C<< four cell points at Zstart and four at Zfinish 
                if(II.eq.IST.and.JJ.eq.JST.and.KK.eq.KST)then
                SAT1X=Tf(II,JJ-1,KK)    
                SAT1Y=Tf(II-1,JJ,KK)                  
                SARH1=HUMF(II,JJ-1,KK)
                Svel1=SQRT((SQRT(abs(Uf(II,JJ-1,KK))**2+
     &               abs(Vf(II,JJ-1,KK))**2)))**2+abs(Wf(II,JJ-1,KK))**2)            
                elseif(II.eq.IST.and.JJ.eq.JFN.and.KK.eq.KST)then
                SAT2X=Tf(II,JJ+1,KK)
                SAT2Y=Tf(II-1,JJ,KK)
                SARH2=HUMF(II,JJ+1,KK)
                Svel2=SQRT((SQRT(abs(Uf(II,JJ+1,KK))**2+
     &               abs(Vf(II,JJ+1,KK))**2)))**2+abs(Wf(II,JJ+1,KK))**2)
                elseif(II.eq.IFN.and.JJ.eq.JST.and.KK.eq.KST)then
                SAT3X=Tf(II,JJ-1,KK)
                SAT3Y=Tf(II+1,JJ,KK)
                SARH3=HUMF(II,JJ-1,KK)
                Svel3=SQRT((SQRT(abs(Uf(II,JJ-1,KK))**2+
     &               abs(Vf(II,JJ-1,KK))**2)))**2+abs(Wf(II,JJ-1,KK))**2)
                elseif(II.eq.IFN.and.JJ.eq.JFN.and.KK.eq.KST)then
                SAT4X=Tf(II,JJ+1,KK)
                SAT4Y=Tf(II+1,JJ,KK)
                SARH4=HUMF(II,JJ+1,KK)
                Svel4=SQRT((SQRT(abs(Uf(II,JJ+1,KK))**2+
     &               abs(Vf(II,JJ+1,KK))**2)))**2+abs(Wf(II,JJ+1,KK))**2)
                elseif(II.eq.IST.and.JJ.eq.JST.and.KK.eq.KFN)then
                FAT1X=Tf(II,JJ-1,KK)
                FAT1Y=Tf(II-1,JJ,KK)                    
                FARH1=HUMF(II,JJ-1,KK)
                Fvel1=SQRT((SQRT(abs(Uf(II,JJ-1,KK))**2+
     &               abs(Vf(II,JJ-1,KK))**2)))**2+abs(Wf(II,JJ-1,KK))**2)
                elseif(II.eq.IST.and.JJ.eq.JFN.and.KK.eq.KFN)then
                FAT2X=Tf(II,JJ+1,KK)
                FAT2Y=Tf(II-1,JJ,KK)
                FARH2=HUMF(II,JJ+1,KK)
                Fvel2=SQRT((SQRT(abs(Uf(II,JJ+1,KK))**2+
     &               abs(Vf(II,JJ+1,KK))**2)))**2+abs(Wf(II,JJ+1,KK))**2)
                elseif(II.eq.IFN.and.JJ.eq.JST.and.KK.eq.KFN)then
                FAT3X=Tf(II,JJ-1,KK)
                FAT3Y=Tf(II+1,JJ,KK)
                FARH3=HUMF(II,JJ-1,KK)
                Fvel3=SQRT((SQRT(abs(Uf(II,JJ-1,KK))**2+
     &               abs(Vf(II,JJ-1,KK))**2)))**2+abs(Wf(II,JJ-1,KK))**2)
                elseif(II.eq.IFN.and.JJ.eq.JFN.and.KK.eq.KFN)then
                FAT4X=Tf(II,JJ+1,KK)
                FAT4Y=Tf(II+1,JJ,KK)
                FARH4=HUMF(II,JJ+1,KK)
                Fvel4=SQRT((SQRT(abs(Uf(II,JJ-1,KK))**2+
     &               abs(Vf(II,JJ-1,KK))**2)))**2+abs(Wf(II,JJ-1,KK))**2)
                endif

          enddo
        enddo
      enddo

      If (VOLNAME(IBLK,ICFD).eq.'head1')then
                tavhead1=(FAT1X+FAT2X+FAT3X+FAT4X+FAT1Y+FAT2Y+FAT3Y
     &                +FAT4Y)/8.0
                hrhead=(SARH1+SARH2+SARH3+SARH4+FARH1+FARH2+FARH3+
     &           FARH4)/8.0
                Velhead1=(Svel1+Svel2+Svel3+Svel4+Fvel1+Fvel2+Fvel3+
     &            Fvel4)/8
      elseif (VOLNAME(IBLK,ICFD).eq.'trunk1')then
          tavtrunk1=(SAT1X+SAT2X+SAT3X+SAT4X+FAT1X+FAT2X+FAT3X+
     &     FAT4X+SAT1Y+SAT2Y+SAT3Y+SAT4Y+FAT1Y+FAT2Y+FAT3Y+FAT4Y)/16.0
                hrtrunk=(SARH1+SARH2+SARH3+SARH4+FARH1+FARH2+FARH3+
     &           FARH4)/8.0
                Veltrunk1=(Svel1+Svel2+Svel3+Svel4+Fvel1+Fvel2+Fvel3+
     &            Fvel4)/8
      elseif (VOLNAME(IBLK,ICFD).eq.'armr1')then
                tavarmr1=(SAT1X+SAT2X+FAT1X+FAT2X+SAT1Y+SAT2Y+
     &            FAT1Y+FAT2Y)/8.0
                hrarmr=(SARH1+SARH2+SARH3+SARH4+FARH1+FARH2+FARH3+
     &           FARH4)/8.0
                Velarmr1=(Svel1+Svel2+Svel3+Svel4+Fvel1+Fvel2+Fvel3+
     &            Fvel4)/8
      elseif (VOLNAME(IBLK,ICFD).eq.'arml1')then
                tavarml1=(SAT3X+SAT4X+FAT3X+FAT4X+SAT3Y+SAT4Y+
     &           FAT3Y+FAT4Y)/8.0
                hrarml=(SARH1+SARH2+SARH3+SARH4+FARH1+FARH2+FARH3+
     &           FARH4)/8.0
                Velarml1=(Svel1+Svel2+Svel3+Svel4+Fvel1+Fvel2+Fvel3+
     &            Fvel4)/8
      elseif (VOLNAME(IBLK,ICFD).eq.'legr1')then
                tavlegr1=(SAT1X+SAT2X+SAT3X+SAT4X+SAT1Y+SAT2Y+SAT3Y+
     &           SAT4Y+FAT1Y+FAT2Y)/10.0
                hrlegr=(SARH1+SARH2+SARH3+SARH4+FARH1+FARH2+FARH3+
     &           FARH4)/8.0
                Vellegr1=(Svel1+Svel2+Svel3+Svel4+Fvel1+Fvel2+Fvel3+
     &            Fvel4)/8
      elseif (VOLNAME(IBLK,ICFD).eq.'legl1')then
                tavlegl1=(SAT1X+SAT2X+SAT3X+SAT4X+SAT1Y+SAT2Y+SAT3Y+
     &           SAT4Y+FAT3Y+FAT4Y)/10.0
                hrlegl=(SARH1+SARH2+SARH3+SARH4+FARH1+FARH2+FARH3+
     &           FARH4)/8.0
                Vellegl1=(Svel1+Svel2+Svel3+Svel4+Fvel1+Fvel2+Fvel3+
     &            Fvel4)/8
                cfdcheck=2
C********************************
      elseif (VOLNAME(IBLK,ICFD).eq.'head2')then
                tavhead2=(FAT1X+FAT2X+FAT3X+FAT4X+FAT1Y+FAT2Y+FAT3Y
     &                +FAT4Y)/8.0
      elseif (VOLNAME(IBLK,ICFD).eq.'trunk2')then
          tavtrunk2=(SAT1X+SAT2X+SAT3X+SAT4X+FAT1X+FAT2X+FAT3X+
     &     FAT4X+SAT1Y+SAT2Y+SAT3Y+SAT4Y+FAT1Y+FAT2Y+FAT3Y+FAT4Y)/16.0
      elseif (VOLNAME(IBLK,ICFD).eq.'armr2')then
                tavarmr2=(SAT1X+SAT2X+FAT1X+FAT2X+SAT1Y+SAT2Y+
     &            FAT1Y+FAT2Y)/8.0
      elseif (VOLNAME(IBLK,ICFD).eq.'arml2')then
                tavarml2=(SAT3X+SAT4X+FAT3X+FAT4X+SAT3Y+SAT4Y+
     &           FAT3Y+FAT4Y)/8.0
      elseif (VOLNAME(IBLK,ICFD).eq.'legr2')then
                tavlegr2=(SAT1X+SAT2X+SAT3X+SAT4X+SAT1Y+SAT2Y+SAT3Y+
     &           SAT4Y+FAT1Y+FAT2Y)/10.0
      elseif (VOLNAME(IBLK,ICFD).eq.'legl2')then
                tavlegl2=(SAT1X+SAT2X+SAT3X+SAT4X+SAT1Y+SAT2Y+SAT3Y+
     &           SAT4Y+FAT3Y+FAT4Y)/10.0
                cfdcheck=2

      endif
!               print *,tavhead1, Velhead1, tavlegl1, vellegl1 
      ENDIF
      RETURN
      END

C ********************* TRBSRCU *********************
C TRBSRCU: stress/velocity at next-to-wall point for solid BC's, U-vel dir.

      FUNCTION TRBSRCU(I,J,K,Dface)
#include "building.h"
#include "cfd.h"
      
      COMMON/INTERP/SIFE(ntcelx),SIFW(ntcelx),SIFN(ntcely),SIFS(ntcely),
     &              SIFH(ntcelz),SIFL(ntcelz)
      COMMON/FLUPRf/URFVIS,VISCOS,PRANDT,SH,
     1            DENf(ntcelx,ntcely,ntcelz),VIS(ntcelx,ntcely,ntcelz),
     2            BETA(ntcelx,ntcely,ntcelz)
      COMMON/VARf/Uf(ntcelx,ntcely,ntcelz),Vf(ntcelx,ntcely,ntcelz),
     1            Wf(ntcelx,ntcely,ntcelz),
     2            P(ntcelx,ntcely,ntcelz),PP(ntcelx,ntcely,ntcelz),
     3            TEf(ntcelx,ntcely,ntcelz),EDf(ntcelx,ntcely,ntcelz)
      COMMON/TURB/GEN(ntcelx,ntcely,ntcelz),CD,CMU,C1,C2,C3,CAPPA,ELOG,
     &            TURBIN,ALAMDA,PRTE,PRED
     
      CDTERM=CMU**0.25

C Calculate y+ by interpolation.
      YPLUSm1=DENf(I-1,J,K)*SQRT(ABS(TEf(I-1,J,K)))*CDTERM*Dface/VISCOS
      YPLUS  =DENf(I,J,K) * SQRT(ABS(TEf(I,J,K))) * CDTERM*Dface/VISCOS
      YPLUSA=SIFW(I)*YPLUS+SIFE(I-1)*YPLUSm1
      
C Calculate stress/velocity at next-to-wall point. Calculation depends upon
C whether we are in laminar sub-layer or log-law region.
      if (YPLUSA.GT.11.63) then

C Interpolate values of k and density.
        SQRTK=SQRT(ABS(SIFW(I)*TEf(I,J,K)+SIFE(I-1)*TEf(I-1,J,K)))
        DEN=SIFW(I)*DENf(I,J,K)+SIFE(I-1)*DENf(I-1,J,K)
        TRBSRCU=(DEN*CDTERM*SQRTK*CAPPA)/ALOG(ELOG*YPLUSA)
      else
        TRBSRCU=VISCOS/Dface
      endif
      
      return
      end
      
C ********************* TRBSRCV *********************
C TRBSRCV: stress/velocity at next-to-wall point for solid BC's, V-vel dir.

      FUNCTION TRBSRCV(I,J,K,Dface)
#include "building.h"
#include "cfd.h"
      
      COMMON/INTERP/SIFE(ntcelx),SIFW(ntcelx),SIFN(ntcely),SIFS(ntcely),
     &              SIFH(ntcelz),SIFL(ntcelz)
      COMMON/FLUPRf/URFVIS,VISCOS,PRANDT,SH,
     1            DENf(ntcelx,ntcely,ntcelz),VIS(ntcelx,ntcely,ntcelz),
     2            BETA(ntcelx,ntcely,ntcelz)
      COMMON/VARf/Uf(ntcelx,ntcely,ntcelz),Vf(ntcelx,ntcely,ntcelz),
     1            Wf(ntcelx,ntcely,ntcelz),
     2            P(ntcelx,ntcely,ntcelz),PP(ntcelx,ntcely,ntcelz),
     3            TEf(ntcelx,ntcely,ntcelz),EDf(ntcelx,ntcely,ntcelz)
      COMMON/TURB/GEN(ntcelx,ntcely,ntcelz),CD,CMU,C1,C2,C3,CAPPA,ELOG,
     &            TURBIN,ALAMDA,PRTE,PRED
     
      CDTERM=CMU**0.25

C Calculate y+ by interpolation.
      YPLUSm1=DENf(I,J-1,K)*SQRT(ABS(TEf(I,J-1,K)))*CDTERM*Dface/VISCOS
      YPLUS  =DENf(I,J,K) * SQRT(ABS(TEf(I,J,K))) * CDTERM*Dface/VISCOS
      YPLUSA=SIFS(J)*YPLUS+SIFN(J-1)*YPLUSm1
      
C Calculate stress/velocity at next-to-wall point. Calculation depends upon
C whether we are in laminar sub-layer or log-law region.
      if (YPLUSA.GT.11.63) then

C Interpolate values of k and density.
        SQRTK=SQRT(ABS(SIFS(J)*TEf(I,J,K)+SIFN(J-1)*TEf(I,J-1,K)))
        DEN=SIFS(J)*DENf(I,J,K)+SIFN(J-1)*DENf(I,J-1,K)
        TRBSRCV=(DEN*CDTERM*SQRTK*CAPPA)/ALOG(ELOG*YPLUSA)
      else
        TRBSRCV=VISCOS/Dface
      endif
      
      return
      end
      
C ********************* TRBSRCW *********************
C TRBSRCW: stress/velocity at next-to-wall point for solid BC's, W-vel dir.

      FUNCTION TRBSRCW(I,J,K,Dface)
#include "building.h"
#include "cfd.h"
      
      COMMON/INTERP/SIFE(ntcelx),SIFW(ntcelx),SIFN(ntcely),SIFS(ntcely),
     &              SIFH(ntcelz),SIFL(ntcelz)
      COMMON/FLUPRf/URFVIS,VISCOS,PRANDT,SH,
     1            DENf(ntcelx,ntcely,ntcelz),VIS(ntcelx,ntcely,ntcelz),
     2            BETA(ntcelx,ntcely,ntcelz)
      COMMON/VARf/Uf(ntcelx,ntcely,ntcelz),Vf(ntcelx,ntcely,ntcelz),
     1            Wf(ntcelx,ntcely,ntcelz),
     2            P(ntcelx,ntcely,ntcelz),PP(ntcelx,ntcely,ntcelz),
     3            TEf(ntcelx,ntcely,ntcelz),EDf(ntcelx,ntcely,ntcelz)
      COMMON/TURB/GEN(ntcelx,ntcely,ntcelz),CD,CMU,C1,C2,C3,CAPPA,ELOG,
     &            TURBIN,ALAMDA,PRTE,PRED
     
      CDTERM=CMU**0.25

C Calculate y+ by interpolation.
      YPLUSm1=DENf(I,J,K-1)*SQRT(ABS(TEf(I,J,K-1)))*CDTERM*Dface/VISCOS
      YPLUS  =DENf(I,J,K) * SQRT(ABS(TEf(I,J,K))) * CDTERM*Dface/VISCOS
      YPLUSA=SIFL(K)*YPLUS+SIFH(K-1)*YPLUSm1
      
C Calculate stress/velocity at next-to-wall point. Calculation depends upon
C whether we are in laminar sub-layer or log-law region.
      if (YPLUSA.GT.11.63) then

C Interpolate values of k and density.
        SQRTK=SQRT(ABS(SIFL(K)*TEf(I,J,K)+SIFH(K-1)*TEf(I,J,K-1)))
        DEN=SIFL(K)*DENf(I,J,K)+SIFH(K-1)*DENf(I,J,K-1)
        TRBSRCW=(DEN*CDTERM*SQRTK*CAPPA)/ALOG(ELOG*YPLUSA)
      else
        TRBSRCW=VISCOS/Dface
      endif
      
      return
      end
      
      
C ********************* YUANBC *********************
C YUANBCW: boundary conditions on solid wall for Yuan wall functions.

      FUNCTION YUANBC(I,J,K,Dface)
#include "building.h"
#include "cfd.h"

      COMMON/FLUPRf/URFVIS,VISCOS,PRANDT,SH,
     1            DENf(ntcelx,ntcely,ntcelz),VIS(ntcelx,ntcely,ntcelz),
     2            BETA(ntcelx,ntcely,ntcelz)
      COMMON/VARf/Uf(ntcelx,ntcely,ntcelz),Vf(ntcelx,ntcely,ntcelz),
     1            Wf(ntcelx,ntcely,ntcelz),
     2            P(ntcelx,ntcely,ntcelz),PP(ntcelx,ntcely,ntcelz),
     3            TEf(ntcelx,ntcely,ntcelz),EDf(ntcelx,ntcely,ntcelz)
      COMMON/TEMPf/Tf(ntcelx,ntcely,ntcelz),GAMH(ntcelx,ntcely,ntcelz),
     1             RESORT,NSWPT,URFT,FSDTT,PRANDL,PFUN
      COMMON/YUANcm/Uqstor(ntcelx,ntcely,ntcelz),Uqinit
     
      LOGICAL Uqinit
      
C Set the friction velocity ($u sub tau$) to initial guess.
      iter = 0
      Utau = 0.06

C Calculate the thermal diffusivity.
      therdif = VISCOS / (PRANDL*DENf(i,j,k))

C Iterate until solution converges.
 10   Utauold = Utau

C Calculate dimensionless distance for inner sublayer ($y sup i sup **$) then
C evaluate wall function ($f sub i$) for inner sublayer. The velocity scale
C based on the heat flux ($U sub q$) was calculated and stored when the
C thermal wall function was evaluated.
      yistst = Dface*(Uqstor(i,j,k)**3) / (therdif*Utauold**2)
      if(yistst .gt. 0.53)then
        fi = 0.228
      else
        fi = yistst*(1.41 - 3.11*yistst + 2.38*yistst*yistst)
      endif

C Calculate dimensionless distance for outer sublayer ($y sup o sup **$) then
C evaluate wall function ($f sub o$) for outer sublayer.
      yostst = Dface*(Uqstor(i,j,k)**7) / (therdif*Utauold**6)
      if(yostst.gt.0.09)then
        fo = 0.024776*EXP(-yostst)
      elseif(yostst.gt.0.005)then
        fo =-0.458-0.258*ALOG(yostst)-0.02425*(ALOG(yostst))**2
      else
        fo = 0.228
      endif

C Set the dimensionless velocity ($u sup **$) according to the wall function.
      ustst = AMIN1(fi,fo)
      ustst = AMAX1(ustst,0.00001)

C Calculate the updated friction velocity. Use relaxation to stabilize solution.
      Utau = ( ((Uqstor(i,j,k))**3) *abs(Wf(i,j,K))/ustst )**0.25
      Utau = 0.5*( Utauold+Utau )

C Check for convergence. Tolerence based on recommendation on p90 of Yuan thesis.
C Allow 15 iterations at least.
      iter = iter + 1
      if( (ABS(Utau-Utauold)/Utau).gt.0.01 .and. iter.le.15 )goto 10
      YUANBC=Utau

      return
      end
      
