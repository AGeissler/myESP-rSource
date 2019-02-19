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

C This file defines ESP-r climate metrics. The idea is that each climate
C metric has an integer index. Data structures of values, names, limits, 
C etc. are all referenced with this integer. Adding a new climate metric
C is then a simple procedure of defining it in this header, and adding 
C code into ESP-r that uses it.

C Index | Metric
C 1     | Dry bulb temperature
C 2     | Diffuse horizontal solar intensity
C 3     | Direct normal solar intensity
C 4     | Global horizonal solar intensity
C 5     | Wind speed
C 6     | Wind direction
C 7     | Relative humidity
C 8     | Total cloud cover
C 9     | Opaque cloud cover
C 10    | Atmospheric pressure

C Parameter definitions. 

C MCM    - Current number of available climate metrics (*24 = record
C          width of climate database)
C NMCM   - Number of values in MCMALL (also the currect climate database
C          version number)
C MCMALL - Values of MCM from all previous versions
C >>NOTE - Whenever the value of MCM is increased in the public release,
C the array MCMALL must be appended with the new value and NMCM must be
C incremented by 1. This maintains back-compatibility for old climate
C databases. The idea is that the ith element in MCMALL is the correct
C value of MCM to use for a version i database.
C CMNAMA - Abreviations (6 characters)
C CMNAMF - Full names (24 characters maximum)
C CMUNIT - Units for values (6 characters maximum)
C CFUNIT - Units for values in climate file (24 characters maximum)
C CFMAX  - Maximum values (in climate file units)
C CFMIN  - Minimum values (in climate file units)

      INTEGER,PARAMETER :: MCM=10

      INTEGER,PARAMETER :: NMCM=2

      INTEGER,PARAMETER :: MCMALL(NMCM)=
     &(/6,10/)

      CHARACTER*6,PARAMETER :: CMNAMA(MCM)=
     &(/'DBTEMP',
     &  'DFHSOL',
     &  'DRNSOL',
     &  'GLHSOL',
     &  'WNDSPD',
     &  'WNDDIR',
     &  'RELHUM',
     &  'TOTCLD',
     &  'OPQCLD',
     &  'ATMPRS'/)

      CHARACTER*24,PARAMETER :: CMNAMF(MCM)=
     &(/'Dry bulb temperature    ',
     &  'Diffuse horizontal solar',
     &  'Direct normal solar     ',
     &  'Global horizontal solar ',
     &  'Wind speed              ',
     &  'Wind direction          ',
     &  'Relative humidity       ',
     &  'Total cloud cover       ',
     &  'Opaque cloud cover      ',
     &  'Atmospheric pressure    '/)

      CHARACTER*6,PARAMETER :: CMUNIT(MCM)=
     &(/'deg C ',
     &  'W/m^2 ',
     &  'W/m^2 ',
     &  'W/m^2 ',
     &  'm/s   ',
     &  'deg CW',
     &  '%     ',
     &  'tenths',
     &  'tenths',
     &  'Pa    '/)

      CHARACTER*24,PARAMETER :: CFUNIT(MCM)=
     &(/'tenths deg C            ',
     &  'W/m**2                  ',
     &  'W/m**2                  ',
     &  'W/m**2                  ',
     &  'tenths m/s              ',
     &  'deg clockwise from north',
     &  'percent                 ',
     &  'tenths of sky covered   ',
     &  'tenths of sky covered   ',
     &  'Pa                      '/)

      INTEGER,PARAMETER :: CFMAX(MCM)=
     &(/570,     ! highest ever recorded
     &  1362,    ! solar constant
     &  1362,    ! solar constant
     &  1362,    ! solar constant
     &  1130,    ! highest ever recorded
     &  360,     ! maximum degrees in a circle
     &  100,     ! maximum percentage
     &  10,      ! maximum tenths
     &  10,      ! maximum tenths
     &  120000/) ! matches EPW file limits

      INTEGER,PARAMETER :: CFMIN(MCM)=
     &(/-900,   ! lowest ever recorded
     &  0,      ! cannot be below zero
     &  0,      ! cannot be below zero
     &  0,      ! cannot be below zero
     &  0,      ! cannot be below zero
     &  0,      ! cannot be below zero
     &  0,      ! cannot be below zero
     &  0,      ! cannot be below zero
     &  0,      ! cannot be below zero
     &  31000/) ! matches EPW file limits

C Common block definitions.

C CLMFIL
C CFVER  - Climate file/database version
C CFYEAR - Climate year
C CFLOC  - Location
C CFLAT  - Latitude
C CFLONG - Longitude
C CFMCM  - Maximum number of metrics in this climate file (*24 = record width)

C CLMMET
C NCM    - Number of climate metrics
C CMCOL  - Columns in the climate file of each metric (not present if 0)
C CMXST  - Flags indicating which metrics are present

C CLMVAL
C CMIVAL - One day of integer values from climate file, plus an extra
C time step for future values of the last hour
C >>Note - Currently this just has intantaneous values read from the
C climate file, and subroutines are defined, marked with tag 
C <CLMNEW2OLD>, to map this to existing old data structures. In future,
C the old data structures should be depreciated.

      COMMON/CLMFIL/CFVER,CFYEAR,CFLAT,CFLONG,CFMCM,CFLOC
      integer CFVER,CFYEAR,CFMCM
      real CFLAT,CFLONG
      character CFLOC*30

      COMMON/CLMMET/MNCM,NCM,CMCOL(MCM),CMXST(MCM)
      integer MNCM,NCM,CMCOL
      logical CMXST

      COMMON/CLMVAL/CMIVAL(MCM,25)
      integer CMIVAL








