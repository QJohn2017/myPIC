MODULE CONSTANTS

IMPLICIT NONE

!Internal units:
!Mass:		m_e
!Time: 		fs
!Distance:	um
!Charge:		e
!Current:		e/fs
!Voltage:		V
!Energy:		m_e um^2 /fs^2
!Velocity:	um/fs
!Momentum:	m_e um/fs

REAL,PARAMETER					:: me = 1.0	 		!Mass of electron (in units of m_e)
REAL,PARAMETER					:: mp = 1836.152673		!Mass of proton
REAL,PARAMETER					:: c = 2.998792458d-1 	!Speed of light um/fs
REAL,PARAMETER					:: eps0 = 3.14208d8		!Permittivity of free space e^2 fs^2 / m_e um^3
REAL,PARAMETER					:: PI = 3.1415926535898
REAL,PARAMETER					:: kB = 1.51563d-11		!Boltzmann constant in m_e um^2/fs^2/K
REAL,PARAMETER					:: nref = 1.0
!Define particle derived type
!----------------------------------------------------------------------------------------------------
TYPE PARTICLE

INTEGER,ALLOCATABLE,DIMENSION(:)		:: FLAG
REAL								:: MASS
REAL								:: CHARGE
REAL,ALLOCATABLE,DIMENSION(:)			:: POS
REAL,ALLOCATABLE,DIMENSION(:)			:: VEL

END TYPE PARTICLE

END MODULE CONSTANTS
