MODULE ROUTINES

USE CONSTANTS
USE FUNCSUB 

IMPLICIT NONE

CONTAINS

!Input file reader
!----------------------------------------------------------------------------------------------------
SUBROUTINE READER(NUM,PARAMS)

INTEGER								:: i
INTEGER,INTENT(IN)						:: NUM
REAL,DIMENSION(21),INTENT(OUT)			:: PARAMS
CHARACTER,DIMENSION(40)					:: BUFFER

!Read data from file
!----------------------------------------------------------------------------------------------------
READ(NUM,*)
READ(NUM,*)

DO i = 1,5
	READ(NUM,*) PARAMS(i)
ENDDO

READ(NUM,*)
READ(NUM,*)
READ(NUM,*)

DO i = 6,14
	READ(NUM,*) PARAMS(i)
ENDDO

READ(NUM,*)
READ(NUM,*)
READ(NUM,*)

DO i = 15,21
	READ(NUM,*) PARAMS(i)
ENDDO

END SUBROUTINE READER

!Gaussian Particle Distribution Generator
!----------------------------------------------------------------------------------------------------
SUBROUTINE PGEN(NUMP,PARTS,ZMAX,RAD,MEAN)

INTEGER								:: i,j
INTEGER,INTENT(IN)						:: NUMP
REAL									:: U,V,R,C,G
REAL,INTENT(IN)						:: RAD,ZMAX,MEAN
TYPE(PARTICLE),INTENT(INOUT)				:: PARTS

i = 1

CALL RANDOM_NUMBER(U)
U = U * 2.0 - 1.0

DO WHILE(i <= NUMP)
	
	CALL RANDOM_NUMBER(V)
	V = V * 2.0 - 1.0
	R = U**2.0 + V**2.0
	
	IF(R > 1.0)THEN
	R = R-1.0
	ENDIF

		C = SQRT(-2.0*LOG(R) / R)

		G = U * C
		U = V
						
		PARTS%POS(i) = V * C * RAD + MEAN
		
		IF(ABS(PARTS%POS(i)) > ZMAX)THEN
		
			i = i
		
		ELSE
		
			i = i + 1
			PARTS%FLAG(i) = 1		
		ENDIF

ENDDO


END SUBROUTINE PGEN

!Max-Bolt Particle Velocity Generator
!---------------------------------------------------------------------------------------------------
SUBROUTINE VGEN(NUMP,PARTS,TEMP)

INTEGER								:: i,j
INTEGER,INTENT(IN)						:: NUMP
REAL									:: U,V,R,C,G,VEL
REAL,INTENT(IN)						:: TEMP
TYPE(PARTICLE),INTENT(INOUT)				:: PARTS

VEL = SQRT(kB * TEMP / PARTS%MASS)

i = 1

CALL RANDOM_NUMBER(U)
U = U * 2.0 - 1.0

DO WHILE(i <= NUMP)
	
	CALL RANDOM_NUMBER(V)
	V = V * 2.0 - 1.0
	R = U**2.0 + V**2.0
	
	IF(R > 1.0)THEN
	R = R-1.0
	ENDIF

		C = SQRT(-2.0*LOG(R) / R)

		G = U * C
		U = V
						
		PARTS%VEL(i) = V * C * VEL
		
		i = i+1

ENDDO



END SUBROUTINE VGEN

!Number Density Calculator
!----------------------------------------------------------------------------------------------------
SUBROUTINE NDENS(PARTS,GRID,DENS,NPARTS,DZ)

!0th order interpolation of particle positions to grid locations
!Weights particle position based on distance from adjacent grid points

INTEGER								:: i,j,COUNTER
INTEGER,INTENT(IN)						:: NPARTS
REAL									:: W1,W2,D1,D2
REAL,INTENT(IN)						:: DZ
REAL,ALLOCATABLE,DIMENSION(:)				:: PACKED
REAL,ALLOCATABLE,DIMENSION(:),INTENT(IN)	:: GRID
REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)	:: DENS
TYPE(PARTICLE),INTENT(IN)				:: PARTS

DENS = 0

DO i = 1,SIZE(GRID)-1
	
	COUNTER = COUNT(PARTS%POS >= GRID(i) .AND. PARTS%POS <= GRID(i+1))
	ALLOCATE(PACKED(COUNTER))
	PACKED = PACK(PARTS%POS,PARTS%POS >= GRID(i) .AND. PARTS%POS <= GRID(i+1))
	
	DO j = 1,SIZE(PACKED)
	
		D1 = PACKED(j) - GRID(i)
		D2 = GRID(i+1) - PACKED(j)
		
		W1 = 1.0 - D1
		W2 = 1.0 - D2
	
		DENS(i) = DENS(i) + W1
		DENS(i+1) = DENS(i+1) + W1
	
	ENDDO	
	DEALLOCATE(PACKED)
ENDDO

DENS = DENS / DZ

END SUBROUTINE NDENS

!Laser Generator
!---------------------------------------------------------------------------------------------------
SUBROUTINE LASER(LFIELD,THIS_T,ZGRID,T0,FWHM,I0)

REAL									:: E0
REAL,INTENT(IN)						:: THIS_T,T0,FWHM,I0
REAL,ALLOCATABLE,DIMENSION(:)				:: THIS_GRID
REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)	:: LFIELD
REAL,ALLOCATABLE,DIMENSION(:),INTENT(IN)	:: ZGRID

ALLOCATE(THIS_GRID(SIZE(ZGRID)))

THIS_GRID = THIS_T - (ZGRID / (c/10.0))

E0 = SQRT(2.0 * I0 / ((c/10.0) * eps0 * nref))

LFIELD = E0 * GAUSS(THIS_GRID,T0,FWHM) 

END SUBROUTINE LASER

!Potential Solver
!---------------------------------------------------------------------------------------------------
SUBROUTINE PHI_SOLVE(POT,DENS,ZPTS,DZ)

!Solves Poisson's equation using finite difference method
!Matrix method is used as it is more efficient than looping
!A.Phi = -(dz^2/eps_0) rho

INTEGER								:: i,INFO
INTEGER,INTENT(IN)						:: ZPTS
INTEGER,ALLOCATABLE,DIMENSION(:)			:: IPIV
REAL,INTENT(IN)						:: DZ	
REAL,ALLOCATABLE,DIMENSION(:)				:: B	
REAL,ALLOCATABLE,DIMENSION(:,:)			:: A
REAL,ALLOCATABLE,DIMENSION(:),INTENT(IN)	:: DENS
REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)	:: POT

ALLOCATE(A(ZPTS,ZPTS),B(ZPTS),IPIV(ZPTS))

!Make tridiagonal A-matrix
A = 0
DO i = 1,ZPTS-1

	A(i,i) = -2
	A(i+1,i) = 1
	A(i,i+1) = 1

ENDDO
A(ZPTS,ZPTS) = -2

B = -(DZ**2.0 / eps0) * DENS

INFO = 0

!DGETRF(# Row of A, # Col of A, A, Lead Dim of A, Pivot positions, Err)
!CALL DGETRF(ZPTS,ZPTS,A,ZPTS,IPIV,INFO)

!DGETRS(Transpose Flag, Dimension of A, # of Col in B, A, 1st Dim of A, Pivot result from LU decomp, B, Lead Dim of B, Err)
!CALL DGETRS(TFLAG,ZPTS,1,A,ZPTS,IPIV,B,ZPTS,INFO)

CALL SGESV(ZPTS,1,A,ZPTS,IPIV,B,ZPTS,INFO)

POT = B

DEALLOCATE(A,B,IPIV)

END SUBROUTINE PHI_SOLVE


!Efield Solver
!----------------------------------------------------------------------------------------------------
SUBROUTINE E_SOLVER(EFIELD,POT,ZPTS,DZ)
!Solves for the E-field from the potential

INTEGER									:: i
INTEGER,INTENT(IN)							:: ZPTS
REAL,INTENT(IN)							:: DZ
REAL,ALLOCATABLE,DIMENSION(:),INTENT(IN)		:: POT
REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)		:: EFIELD

EFIELD = 0
EFIELD(1) = (POT(2)-POT(1)) / DZ

DO i = 2,ZPTS-1

	EFIELD(i) = (POT(i+1) - POT(i-1)) / (2.0 * DZ)

ENDDO

EFIELD(ZPTS) = (POT(ZPTS)-POT(ZPTS-1)) / DZ

EFIELD = -EFIELD

END SUBROUTINE E_SOLVER


!Velocity Half-Step Pusher
!----------------------------------------------------------------------------------------------------
SUBROUTINE HALF_PUSHER(IONS,EONS,EFIELD,NE,NI,TGRID,DT,ZGRID,ZPTS,DZ)

INTEGER									:: i,LOC1,LOC2
INTEGER,DIMENSION(1)						:: LOC
INTEGER,INTENT(IN)							:: ZPTS,NE,NI
REAL										:: E1,E2,EAVG,POS 
REAL,INTENT(IN)							:: DT,DZ
REAL,ALLOCATABLE,DIMENSION(:),INTENT(IN)		:: EFIELD,TGRID,ZGRID
TYPE(PARTICLE),INTENT(INOUT)					:: IONS,EONS

!Push electrons
DO i = 1,NE
	
	IF(EONS%FLAG(i) == 1)THEN	

	POS = EONS%POS(i)
	LOC = MINLOC(ABS(ZGRID-POS))
	LOC1 = LOC(1)

	IF(ZGRID(LOC1) >= POS)THEN
	
		LOC2 = LOC1
		LOC1 = LOC1 - 1

	ELSE
	
		LOC2 = LOC1 + 1

	ENDIF

	E1 = EFIELD(LOC1)
	E2 = EFIELD(LOC2)
	EAVG = (E1 + (POS - ZGRID(LOC1))* ((E2 - E1) / DZ)) 

	EONS%VEL(i) = EONS%VEL(i) + (EONS%CHARGE / EONS%MASS) * EAVG * (DT/2.0)
	
	ENDIF

ENDDO

!Push ions
DO i = 1,NI
	
	IF(IONS%FLAG(i) == 1)THEN
	
	POS = IONS%POS(i)
	LOC = MINLOC(ABS(ZGRID-POS))
	LOC1 = LOC(1)	

	IF(ZGRID(LOC1) >= POS)THEN
	
		LOC2 = LOC1
		LOC1 = LOC1 - 1

	ELSE
	
		LOC2 = LOC1 + 1

	ENDIF

	E1 = EFIELD(LOC1)
	E2 = EFIELD(LOC2)
	EAVG = E1 + (POS - ZGRID(LOC1))* ((E2 - E1) / DZ) 

	IONS%VEL(i) = IONS%VEL(i) + (IONS%CHARGE / IONS%MASS) * EAVG * (DT/2.0)

	ENDIF

ENDDO


END SUBROUTINE HALF_PUSHER

!Particle Pusher
!----------------------------------------------------------------------------------------------------
SUBROUTINE PUSHER(IONS,EONS,EFIELD,NE,NI,TGRID,DT,ZGRID,ZPTS,DZ)

INTEGER									:: i,j,LOC1,LOC2,NEP,NEI
INTEGER,DIMENSION(1)						:: LOC
INTEGER,ALLOCATABLE,DIMENSION(:)				:: PARTS
INTEGER,INTENT(INOUT)						:: ZPTS,NE,NI
REAL										:: E1,E2,EAVG,POS 
REAL,INTENT(IN)							:: DT,DZ
REAL,ALLOCATABLE,DIMENSION(:),INTENT(IN)		:: EFIELD,TGRID,ZGRID
TYPE(PARTICLE),INTENT(INOUT)					:: IONS,EONS

!NEP = COUNT(ABS(EONS%POS) .LT. ZGRID(ZPTS))
!ALLOCATE(PARTS(NEP))

!DO i = 1,NE
!	IF(EONS%POS(i) .LT. ZGRID(ZPTS))THEN
!		PARTS(i) = i
!	ENDIF
!ENDDO

!EONS%VEL = PACK(EONS%VEL,EONS%VEL == EONS%VEL(PARTS))
!EONS%POS = PACK(EONS%POS,ABS(EONS%POS) .LT. ZGRID(ZPTS))

!Push electrons
DO i = 1,NE

	IF(EONS%FLAG(i) == 1)THEN
	
	POS = EONS%POS(i)
	LOC = MINLOC(ABS(ZGRID-POS))
	LOC1 = LOC(1)

	IF(ZGRID(LOC1) >= POS)THEN
	
		LOC2 = LOC1
		LOC1 = LOC1 - 1

	ELSE
	
		LOC2 = LOC1 + 1

	ENDIF

	E1 = EFIELD(LOC1)
	E2 = EFIELD(LOC2)
	EAVG = E1 + (POS - ZGRID(LOC1))* ((E2 - E1) / DZ) 

	EONS%VEL(i) = EONS%VEL(i) + (EONS%CHARGE / EONS%MASS) * EAVG * DT
	EONS%POS(i) = EONS%POS(i) + EONS%VEL(i) * DT

	IF(ABS(EONS%POS(i)) .GT. ZGRID(ZPTS))THEN
		EONS%FLAG(i) = 0
	ELSE
		EONS%FLAG(i) = 1
	ENDIF

	ENDIF

ENDDO

!DEALLOCATE(PARTS)


!Push ions
DO i = 1,NI
	
	IF(IONS%FLAG(i) == 1)THEN

	POS = IONS%POS(i)
	LOC = MINLOC(ABS(ZGRID-POS))
	LOC1 = LOC(1)	

	IF(ZGRID(LOC1) >= POS)THEN
	
		LOC2 = LOC1
		LOC1 = LOC1 - 1

	ELSE
	
		LOC2 = LOC1 + 1

	ENDIF

	E1 = EFIELD(LOC1)
	E2 = EFIELD(LOC2)
	EAVG = E1 + (POS - ZGRID(LOC1))* ((E2 - E1) / DZ) 

	IONS%VEL(i) = IONS%VEL(i) + (IONS%CHARGE / IONS%MASS) * EAVG * DT
	IONS%POS(i) = IONS%POS(i) + IONS%VEL(i) * DT
	
	IF(ABS(IONS%POS(i)) .GT. ZGRID(ZPTS))THEN
		IONS%FLAG(i) = 0
	ELSE
		IONS%FLAG(i) = 1
	ENDIF

	ENDIF

ENDDO

END SUBROUTINE PUSHER

END MODULE ROUTINES
