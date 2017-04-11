PROGRAM MYPIC

USE MPI
USE CONSTANTS
USE FUNCSUB
USE ROUTINES

IMPLICIT NONE

!Simple Variables
LOGICAL								:: EXISTER,OPENER
INTEGER								:: NUM,i,j

!MPI Variables
INTEGER								:: npes,rank,ierr


!Global grids
REAL									:: DT,DZ,ZMAX,ZMIN
INTEGER								:: ZPTS,TPTS
REAL,ALLOCATABLE,DIMENSION(:)				:: ZGRID,TGRID

!Input parameters
REAL									:: TFWHM,FRAD,BETA,INTENS,T0,PHI2,PHI3,PHI4,MASS_I,PRAD,LAMBDA_0,MEAN1,MEAN2,TEMP
INTEGER								:: NUM_I,I_STATE,NUM_E,TMOD
REAL,DIMENSION(21)						:: PARAMS

!Particles and fields
TYPE(PARTICLE)							:: IONS,EONS
REAL,ALLOCATABLE,DIMENSION(:)				:: NUMDENSI,NUMDENSE,CUR,CDENS,POT,VSYNCE,VSYNCI
REAL,ALLOCATABLE,DIMENSION(:)				:: EFIELD,LFIELD,BFIELD !These need to be made complex

!Laser Stuff
REAL									:: OMEGA_0
CHARACTER								:: POL

!Begin MPI
!---------------------------------------------------------------------------------------------------
CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,npes,ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

!Startup File Check
!---------------------------------------------------------------------------------------------------
IF(rank == 0)THEN

NUM = 10

INQUIRE(FILE = 'bin/../startup.txt',EXIST = EXISTER)

IF(EXISTER)THEN

	OPEN(UNIT=NUM,FILE='bin/../startup.txt')

	INQUIRE(NUM,OPENED=OPENER)
	
		IF(OPENER)THEN
	
			WRITE(*,*) 'Startup file succesfully opened!'
	
		ELSE
	
			WRITE(*,*) 'Startup file not opened! Exiting now!'
			STOP
		
		ENDIF

ELSE
	
	WRITE(*,*) 'No startup file detected. Exiting'
	STOP
	
ENDIF

ENDIF

!Read Startup File
!---------------------------------------------------------------------------------------------------
IF(rank == 0)THEN
	CALL READER(NUM,PARAMS)
	CLOSE(NUM)
ENDIF
CALL MPI_BCAST(PARAMS,SIZE(PARAMS),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!Assign grids and startup values
!---------------------------------------------------------------------------------------------------
TPTS = 2**INT(PARAMS(1))
ZPTS = INT(PARAMS(2))
DT = PARAMS(3)
DZ = PARAMS(4)
TMOD = INT(PARAMS(5))
TFWHM = PARAMS(6)
LAMBDA_0 = PARAMS(7)
FRAD = PARAMS(8)
BETA = PARAMS(9)
INTENS = PARAMS(10)
T0 = PARAMS(11)
PHI2 = PARAMS(12)
PHI3 = PARAMS(13)
PHI4 = PARAMS(14)
NUM_I = PARAMS(15)
I_STATE = PARAMS(16)
MASS_I = PARAMS(17)
PRAD = PARAMS(18)
MEAN1 = PARAMS(19)
MEAN2 = PARAMS(20)
TEMP = PARAMS(21)
NUM_E = NUM_I * I_STATE
OMEGA_0 = 2 * PI * C / LAMBDA_0

IF(rank == 0)THEN
	IF(DT > (2.0 * PI / OMEGA_0))THEN

		WRITE(*,*) 'Time step too large. Scaling to w0...'
		WRITE(*,*) 'DT is now: ',(2.0 * PI / OMEGA_0)

		DT = (2.0 * PI / OMEGA_0)

	ENDIF
	WRITE(*,*)'Inputs read'
ENDIF

!Grid Generation
!---------------------------------------------------------------------------------------------------
IF(rank == 0)THEN
	WRITE(*,*) 'The z-grid size is: ',DZ * ZPTS,'um'
	WRITE(*,*) 'The t-grid size is: ',DT * TPTS,'fs'

	WRITE(*,*) 'Number of Output Writes: ',FLOOR(REAL(TPTS) / REAL(TMOD)) + 2
	WRITE(*,*) 'Pulse Center: ',T0 * c,'(um)'
ENDIF

ALLOCATE(TGRID(TPTS),ZGRID(ZPTS))

CALL TGRIDDER(TPTS,DT,TGRID)
CALL XYZGRIDDER(ZPTS,DZ,ZGRID)

IF(rank == 0)THEN
	OPEN(UNIT=12,FILE='res/Tgrid.txt')
	OPEN(UNIT=13,FILE='res/Zgrid.txt')

	WRITE(12,*) TGRID
	WRITE(13,*) ZGRID

	CLOSE(12)
	CLOSE(13)

	WRITE(*,*)'Grids generated'
ENDIF

ZMAX = ZGRID(ZPTS)
ZMIN = ZGRID(1)
!Particle Generation
!---------------------------------------------------------------------------------------------------
ALLOCATE(IONS%FLAG(NUM_I),IONS%POS(NUM_I),IONS%VEL(NUM_I))
ALLOCATE(EONS%FLAG(NUM_E),EONS%POS(NUM_E),EONS%VEL(NUM_E))

ALLOCATE(NUMDENSI(ZPTS),NUMDENSE(ZPTS))

IONS%CHARGE = I_STATE
EONS%CHARGE = -1.0

IONS%MASS = MASS_I
EONS%MASS = me

OPEN(UNIT=100,FILE='res/Times.txt')

WRITE(100,*) TGRID(1)

CALL PGEN(NUM_I,IONS,ZGRID(ZPTS),PRAD/2.0,MEAN1)
WRITE(*,*)'Ions Generated'

CALL PGEN(NUM_E,EONS,ZGRID(ZPTS),PRAD/2.0,MEAN2)
WRITE(*,*)'Electrons Generated'

OPEN(UNIT=14,FILE='res/Pos_I.txt')
OPEN(UNIT=15,FILE='res/Pos_E.txt')

WRITE(14,*)IONS%POS
WRITE(15,*)EONS%POS

CALL VGEN(NUM_I,IONS,TEMP)
CALL VGEN(NUM_E,EONS,TEMP)

OPEN(UNIT=21,FILE='res/Vel_I.txt')
OPEN(UNIT=22,FILE='res/Vel_E.txt')

WRITE(21,*) IONS%VEL
WRITE(22,*) EONS%VEL

CALL NDENS(IONS,ZGRID,NUMDENSI,NUM_I,DZ)
CALL NDENS(EONS,ZGRID,NUMDENSE,NUM_E,DZ)

OPEN(UNIT=16,FILE='res/IDens.txt')
OPEN(UNIT=17,FILE='res/EDens.txt')

WRITE(16,*) NUMDENSI
WRITE(17,*) NUMDENSE

!CLOSE(16)
!CLOSE(17)

CDENS = IONS%CHARGE*NUMDENSI+EONS%CHARGE*NUMDENSE

OPEN(UNIT=20,FILE='res/CDens.txt')
WRITE(20,*) CDENS


PRINT*,'Initial charge density calculated'
!Generate initial Electric Potential
!---------------------------------------------------------------------------------------------------
ALLOCATE(POT(ZPTS))
CALL PHI_SOLVE(POT,CDENS,ZPTS,DZ)

OPEN(UNIT=18,FILE='res/EPot.txt')

WRITE(18,*) POT
!CLOSE(18)
PRINT*,'Initial Electric Potential calculated'

!Generate initial Electric Field
!---------------------------------------------------------------------------------------------------
ALLOCATE(EFIELD(ZPTS))

CALL E_SOLVER(EFIELD,POT,ZPTS,DZ)

OPEN(UNIT=19,FILE='res/EField.txt')
!WRITE(19,*) EFIELD
!CLOSE(19)
PRINT*, 'Initial Electric Field calculated'

!Generate Laser Electric Field
!---------------------------------------------------------------------------------------------------
ALLOCATE(LFIELD(ZPTS))

INTENS = INTENS * 1.0977d-11

CALL LASER(LFIELD,TGRID(1),ZGRID,T0,TFWHM,INTENS)
 
OPEN(UNIT=23,FILE='res/LField.txt')
WRITE(23,*) LFIELD
PRINT*, 'Laser Generated'

EFIELD = EFIELD + LFIELD

WRITE(19,*) EFIELD

!First Particle Push
!---------------------------------------------------------------------------------------------------
ALLOCATE(VSYNCE(NUM_E),VSYNCI(NUM_I))

PRINT*,'Starting the push'
CALL HALF_PUSHER(IONS,EONS,EFIELD,NUM_E,NUM_I,TGRID,DT,ZGRID,ZPTS,DZ)

!Main Loop
!---------------------------------------------------------------------------------------------------
PRINT*,'Looping...'
DO i = 2,TPTS

CALL PUSHER(IONS,EONS,EFIELD,NUM_E,NUM_I,TGRID,DT,ZGRID,ZPTS,DZ) 
CALL NDENS(IONS,ZGRID,NUMDENSI,NUM_I,DZ)
CALL NDENS(EONS,ZGRID,NUMDENSE,NUM_E,DZ)
CDENS = IONS%CHARGE*NUMDENSI+EONS%CHARGE*NUMDENSE
CALL PHI_SOLVE(POT,CDENS,ZPTS,DZ)
CALL E_SOLVER(EFIELD,POT,ZPTS,DZ)
CALL LASER(LFIELD,TGRID(i),ZGRID,T0,TFWHM,INTENS)
EFIELD = EFIELD + LFIELD
	IF(MOD(i,TMOD) == 0 .OR. i == TPTS)THEN
		
		WRITE(*,*) 'Syncing and Writing... ',i
		!WRITE(*,*) SIZE(EONS%POS)	
		VSYNCE = EONS%VEL
		VSYNCI = IONS%VEL

		CALL HALF_PUSHER(IONS,EONS,EFIELD,NUM_E,NUM_I,TGRID,DT,ZGRID,ZPTS,DZ)
		WRITE(100,*) TGRID(i)
		WRITE(14,*) IONS%POS
		WRITE(15,*) EONS%POS
		WRITE(16,*) NUMDENSI
		WRITE(17,*) NUMDENSE
		WRITE(18,*) POT
		WRITE(19,*) EFIELD
		WRITE(20,*) CDENS
		WRITE(21,*) IONS%VEL
		WRITE(22,*) EONS%VEL
		WRITE(23,*) LFIELD
	
		EONS%VEL = VSYNCE
		IONS%VEL = VSYNCI

	ENDIF

ENDDO

PRINT*,'Done Looping...'

!OPEN(UNIT=11,FILE='res/E1.txt')
!OPEN(UNIT=12,FILE='res/Pot1.txt')
!OPEN(UNIT=13,FILE='res/CDens.txt')

!WRITE(11,*)EFIELD
!WRITE(12,*)POT
!WRITE(13,*)CDENS

CALL MPI_FINALIZE(ierr)

DEALLOCATE(TGRID,ZGRID)

END PROGRAM MYPIC
