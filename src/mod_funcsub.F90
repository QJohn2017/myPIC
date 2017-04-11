MODULE FUNCSUB

USE CONSTANTS

IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------------------------------

FUNCTION GAUSS(X,X0,FWHM) RESULT(A)

!Input numbers
REAL,INTENT(IN)							:: X0,FWHM
INTEGER									:: NPTS
REAL										:: W

!
REAL,ALLOCATABLE,DIMENSION(:),INTENT(IN)		:: X
REAL,ALLOCATABLE,DIMENSION(:)					:: A

NPTS = SIZE(X)

ALLOCATE(A(NPTS))

W = FWHM / SQRT(2.0 * LOG(2.0))

A = EXP(-2 * ((X - X0)**2) / (W**2))

END FUNCTION GAUSS

!Time grid generator
!----------------------------------------------------------------------------------------------------
SUBROUTINE TGRIDDER(TPTS,DT,TGRID)

INTEGER								:: i
INTEGER,INTENT(IN)						:: TPTS
REAL,INTENT(IN)						:: DT
REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)	:: TGRID

DO i = 0,TPTS-1

	TGRID(i+1) = i*DT 

ENDDO

END SUBROUTINE TGRIDDER

!Spatial grid generator
!----------------------------------------------------------------------------------------------------
SUBROUTINE XYZGRIDDER(ZPTS,DZ,ZGRID)

INTEGER								:: i,j
INTEGER,INTENT(IN)						:: ZPTS
REAL,INTENT(IN)						:: DZ
REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)	:: ZGRID

j = -(ZPTS-1)/2

DO i = 1,ZPTS
	
	ZGRID(i) = j*DZ 
	j = j + 1

ENDDO

END SUBROUTINE XYZGRIDDER

END MODULE FUNCSUB
