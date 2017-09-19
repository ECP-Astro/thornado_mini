MODULE TimeSteppingModule

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    UnitsDisplay
  USE ProgramHeaderModule, ONLY: &
    nX, nDOFX, nNodes, &
    nE, nDOF
  USE MeshModule, ONLY: &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    a, b, c
  USE FluidFieldsModule, ONLY: &
    uCF, nCF
  USE RadiationFieldsModule, ONLY: &
    uCR, rhsCR, nCR, nSpecies
  USE InputOutputModule, ONLY: &
    WriteFields1D, &
    WriteFieldsRestart1D
  USE TallyModule, ONLY: &
    InitializeGlobalTally, &
    ComputeGlobalTally
  USE RadiationEvolutionModule, ONLY: &
    ComputeRHS_Radiation, &
    ApplySlopeLimiter_Radiation, &
    ApplyPositivityLimiter_Radiation
  USE FluidRadiationCouplingModule, ONLY: &
    CoupleFluidRadiation
  USE BoundaryConditionsModule, ONLY: &
    ApplyBoundaryConditions_Radiation

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  LOGICAL :: DisplayTimers = .FALSE.

  LOGICAL :: EvolveGravity
  LOGICAL :: EvolveFluid
  LOGICAL :: EvolveRadiation
  INTEGER :: nStages_SSP_RK
  INTEGER :: nStages_SI_RK
  REAL(DP), DIMENSION(:,:,:,:,:),     ALLOCATABLE :: uCF_0
  REAL(DP), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: uCR_0

  PROCEDURE (TimeStepper), POINTER, PUBLIC :: &
    SSP_RK => NULL(), &
    SI_RK  => NULL()

  INTERFACE
    SUBROUTINE TimeStepper( t, dt )
      USE KindModule, ONLY: DP
      REAL(DP), INTENT(in) :: t, dt
    END SUBROUTINE TimeStepper
  END INTERFACE

  PUBLIC :: InitializeTimeStepping
  PUBLIC :: FinalizeTimeStepping
  PUBLIC :: EvolveFields
  PUBLIC :: BackwardEuler

CONTAINS


  SUBROUTINE InitializeTimeStepping &
               ( EvolveGravity_Option, EvolveFluid_Option, &
                 EvolveRadiation_Option, nStages_SSP_RK_Option, &
                 nStages_SI_RK_Option )

    LOGICAL, INTENT(in), OPTIONAL :: EvolveGravity_Option
    LOGICAL, INTENT(in), OPTIONAL :: EvolveFluid_Option
    LOGICAL, INTENT(in), OPTIONAL :: EvolveRadiation_Option
    INTEGER, INTENT(in), OPTIONAL :: nStages_SSP_RK_Option
    INTEGER, INTENT(in), OPTIONAL :: nStages_SI_RK_Option

    EvolveGravity = .FALSE.
    IF( PRESENT( EvolveGravity_Option ) ) &
      EvolveGravity = EvolveGravity_Option
    WRITE(*,*)
    WRITE(*,'(A5,A19,L1)') '', &
      'Evolve Gravity = ', EvolveGravity

    EvolveFluid = .FALSE.
    IF( PRESENT( EvolveFluid_Option ) ) &
      EvolveFluid = EvolveFluid_Option
    WRITE(*,'(A5,A19,L1)') '', &
      'Evolve Fluid = ', EvolveFluid

    EvolveRadiation = .FALSE.
    IF( PRESENT( EvolveRadiation_Option ) ) &
      EvolveRadiation = EvolveRadiation_Option
    WRITE(*,'(A5,A19,L1)') '', &
      'Evolve Radiation = ', EvolveRadiation
    WRITE(*,*)

    nStages_SSP_RK = 1
    IF( PRESENT( nStages_SSP_RK_Option ) ) &
      nStages_SSP_RK = nStages_SSP_RK_Option

    SELECT CASE ( nStages_SSP_RK )
      CASE ( 1 )

        SSP_RK => SSP_RK1

      CASE ( 2 )

        SSP_RK => SSP_RK2

      CASE ( 3 )

        SSP_RK => SSP_RK3

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A4,A43,I2.2)') &
          '', 'SSP_RK not implemented for nStatgesSSPRK = ', nStages_SSP_RK
        STOP

    END SELECT

    nStages_SI_RK = 1
    IF( PRESENT( nStages_SI_RK_Option ) ) &
      nStages_SI_RK = nStages_SI_RK_Option

    SELECT CASE ( nStages_SI_RK )
      CASE ( 1 )

        SI_RK => SI_RK1

      CASE ( 2 )

        SI_RK => SI_RK2

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A4,A43,I2.2)') &
          '', 'SI_RK not implemented for nStatges_SI_RK = ', nStages_SI_RK
        STOP

    END SELECT

  END SUBROUTINE InitializeTimeStepping


  SUBROUTINE FinalizeTimeStepping

    NULLIFY( SSP_RK, SI_RK )

  END SUBROUTINE FinalizeTimeStepping


  SUBROUTINE EvolveFields &
               ( t_begin, t_end, dt_write, UpdateFields, dt_fixed_Option )

    REAL(DP), INTENT(in) :: t_begin, t_end, dt_write
    PROCEDURE (TimeStepper) :: UpdateFields
    REAL(DP), INTENT(in), OPTIONAL :: dt_fixed_Option

    LOGICAL  :: WriteOutput = .FALSE.
    LOGICAL  :: FixedTimeStep
    INTEGER  :: iCycle
    REAL(DP) :: t, t_write, dt, dt_fixed
    REAL(DP), DIMENSION(0:1) :: WallTime

    FixedTimeStep = .FALSE.
    IF( PRESENT( dt_fixed_Option ) )THEN
      FixedTimeStep = .TRUE.
      dt_fixed = dt_fixed_Option
    END IF

    ASSOCIATE( U => UnitsDisplay )

    WRITE(*,*)
    WRITE(*,'(A4,A21)') '', 'INFO: Evolving Fields'
    WRITE(*,*)
    WRITE(*,'(A6,A10,ES10.4E2,A1,2A2,A8,ES10.4E2,&
              &A1,2A2,A11,ES10.4E2,A1,A2)') &
      '', 't_begin = ',  t_begin  / U % TimeUnit, '', TRIM( U % TimeLabel ), &
      '', 't_end = ',    t_end    / U % TimeUnit, '', TRIM( U % TimeLabel ), &
      '', 'dt_write = ', dt_write / U % TimeUnit, '', TRIM( U % TimeLabel )
    WRITE(*,*)

    iCycle  = 0
    t       = t_begin
    t_write = t_begin + dt_write

    IF( EvolveRadiation )THEN

      CALL ApplyPositivityLimiter_Radiation

    END IF

    CALL InitializeGlobalTally &
           ( Time = t, &
             TallyGravity_Option = EvolveGravity, &
             TallyFluid_Option = EvolveFluid,  &
             TallyRadiation_Option = EvolveRadiation )

    CALL WriteFields1D &
           ( Time = t, &
             WriteGeometryFields_Option = EvolveGravity, &
             WriteFluidFields_Option = EvolveFluid, &
             WriteRadiationFields_Option = EvolveRadiation )

    CALL WriteFieldsRestart1D &
           ( Time = t, &
             WriteGeometryFields_Option = EvolveGravity, &
             WriteFluidFields_Option = EvolveFluid, &
             WriteRadiationFields_Option = EvolveRadiation )

    WallTime(0) = MPI_WTIME( )

    DO WHILE( t < t_end )

      iCycle = iCycle + 1

      IF( FixedTimeStep )THEN
        dt = dt_fixed
      ELSE
        CALL ComputeTimeStep( dt )
      END IF

      IF( t + dt > t_end )THEN

        dt = t_end - t

      END IF

      IF( t + dt > t_write )THEN

        dt          = t_write - t
        t_write     = t_write + dt_write
        WriteOutput = .TRUE.

      END IF

      IF( MOD( iCycle, 100 ) == 0 )THEN

        WRITE(*,'(A8,A8,I8.8,A2,A4,ES12.6E2,A1,A2,A2,A5,ES12.6E2,A1,A2)') &
          '', 'Cycle = ', iCycle, &
          '', 't = ', t / U % TimeUnit, '', TRIM( U % TimeLabel ), &
          '', 'dt = ', dt / U % TimeUnit, '', TRIM( U % TimeLabel )

        CALL ComputeGlobalTally( Time = t )

      END IF

      CALL UpdateFields( t, dt )

      t = t + dt

      IF( WriteOutput )THEN

        CALL WriteFields1D &
               ( Time = t, &
                 WriteGeometryFields_Option = EvolveGravity, &
                 WriteFluidFields_Option = EvolveFluid, &
                 WriteRadiationFields_Option = EvolveRadiation )

        CALL WriteFieldsRestart1D &
               ( Time = t, &
                 WriteGeometryFields_Option = EvolveGravity, &
                 WriteFluidFields_Option = EvolveFluid, &
                 WriteRadiationFields_Option = EvolveRadiation )

        WriteOutput = .FALSE.

      END IF

    END DO

    WallTime(1) = MPI_WTIME( )

    WRITE(*,*)
    WRITE(*,'(A6,A15,ES10.4E2,A1,A2,A6,I8.8,A7,A4,ES10.4E2,A2)') &
      '', 'Evolved to t = ', t / U % TimeUnit, '', TRIM( U % TimeLabel ), &
      ' with ', iCycle, ' cycles', &
      ' in ', WallTime(1)-WallTime(0), ' s'
    WRITE(*,*)

    CALL WriteFields1D &
           ( Time = t, &
             WriteGeometryFields_Option = EvolveGravity, &
             WriteFluidFields_Option = EvolveFluid, &
             WriteRadiationFields_Option = EvolveRadiation )

    CALL WriteFieldsRestart1D &
           ( Time = t, &
             WriteGeometryFields_Option = EvolveGravity, &
             WriteFluidFields_Option = EvolveFluid, &
             WriteRadiationFields_Option = EvolveRadiation )

    END ASSOCIATE ! U

  END SUBROUTINE EvolveFields


  SUBROUTINE ComputeTimeStep( dt )

    REAL(DP), INTENT(out) :: dt

    REAL(DP) :: dt_Fluid, dt_Radiation

    dt_Fluid     = HUGE( 1.0_DP )
    dt_Radiation = HUGE( 1.0_DP )

    IF( EvolveFluid )THEN

      CALL ComputeTimeStep_Fluid( dt_Fluid )

    END IF

    IF( EvolveRadiation )THEN

      CALL ComputeTimeStep_Radiation( dt_Radiation )

    END IF

    dt = MIN( dt_Fluid, dt_Radiation )

  END SUBROUTINE ComputeTimeStep


  SUBROUTINE ComputeTimeStep_Fluid( dt )

    REAL(DP), INTENT(out) :: dt

    dt = HUGE( 1.0_DP )

  END SUBROUTINE ComputeTimeStep_Fluid


  SUBROUTINE ComputeTimeStep_Radiation( dt )

    REAL(DP), INTENT(out) :: dt

    INTEGER  :: iX1, iX2, iX3, iNodeX
    REAL(DP) :: CFL, dt_X1, dt_X2, dt_X3

    dt = HUGE( 1.0_DP )

    ASSOCIATE( dX1 => MeshX(1) % Width (1:nX(1)), &
               dX2 => MeshX(2) % Width (1:nX(2)), &
               dX3 => MeshX(3) % Width (1:nX(3)), &
               X1  => MeshX(1) % Center(1:nX(1)), &
               X2  => MeshX(2) % Center(1:nX(2)), &
               X3  => MeshX(3) % Center(1:nX(3)) )

    CFL = 0.45_DP / ( 2.0_DP * DBLE( nNodes - 1 ) + 1.0_DP )

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)
        DO iX1 = 1, nX(1)

          dt_X1 = CFL * dX1(iX1)
          dt_X2 = CFL * a( [ X1(iX1), X2(iX2), X3(iX3) ] ) &
                      * dX2(iX2)
          dt_X3 = CFL * b( [ X1(iX1), X2(iX2), X3(iX3) ] ) &
                      * c( [ X1(iX1), X2(iX2), X3(iX3) ] ) &
                      * dX3(iX3)

          dt = MIN( dt, dt_x1, dt_X2, dt_X3 )

        END DO
      END DO
    END DO

    END ASSOCIATE ! dX1, etc.

  END SUBROUTINE ComputeTimeStep_Radiation


  SUBROUTINE SSP_RK1( t, dt )

    REAL(DP), INTENT(in) :: t, dt

    CALL Initialize_SSP_RK

    IF( EvolveRadiation )THEN

      CALL ApplyBoundaryConditions_Radiation( t )

      CALL ComputeRHS_Radiation &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ] )

      CALL ApplyRHS_Radiation &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, alpha = 0.0_DP, beta = 1.0_DP )

      CALL ApplyBoundaryConditions_Radiation &
             ( t + dt, LimiterBC_Option = .TRUE. )

      CALL ApplySlopeLimiter_Radiation

      CALL ApplyPositivityLimiter_Radiation

    END IF

    CALL Finalize_SSP_RK

  END SUBROUTINE SSP_RK1


  SUBROUTINE SSP_RK2( t, dt )

    REAL(DP), INTENT(in) :: t, dt

    CALL Initialize_SSP_RK

    ! -- RK Stage 1 --

    IF( EvolveRadiation )THEN

      CALL ApplyBoundaryConditions_Radiation( t )

      CALL ComputeRHS_Radiation &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ] )

      CALL ApplyRHS_Radiation &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, alpha = 0.0_DP, beta = 1.0_DP )

      CALL ApplyBoundaryConditions_Radiation &
             ( t + dt, LimiterBC_Option = .TRUE. )

      CALL ApplySlopeLimiter_Radiation

      CALL ApplyPositivityLimiter_Radiation

    END IF

    ! -- RK Stage 2 --

    IF( EvolveRadiation )THEN

      CALL ApplyBoundaryConditions_Radiation( t + dt )

      CALL ComputeRHS_Radiation &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ] )

      CALL ApplyRHS_Radiation &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, alpha = 0.5_DP, beta = 0.5_DP )

      CALL ApplyBoundaryConditions_Radiation &
             ( t + dt, LimiterBC_Option = .TRUE. )

      CALL ApplySlopeLimiter_Radiation

      CALL ApplyPositivityLimiter_Radiation

    END IF

    CALL Finalize_SSP_RK

  END SUBROUTINE SSP_RK2


  SUBROUTINE SSP_RK3( t, dt )

    REAL(DP), INTENT(in) :: t, dt

    CALL Initialize_SSP_RK

    ! -- RK Stage 1 -- 

    IF( EvolveRadiation )THEN

      CALL ApplyBoundaryConditions_Radiation( t )

      CALL ComputeRHS_Radiation &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ] )

      CALL ApplyRHS_Radiation &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, alpha = 0.0_DP, beta = 1.0_DP )

      CALL ApplyBoundaryConditions_Radiation &
             ( t + dt, LimiterBC_Option = .TRUE. )

      CALL ApplySlopeLimiter_Radiation

      CALL ApplyPositivityLimiter_Radiation

    END IF

    ! -- RK Stage 2 --

    IF( EvolveRadiation )THEN

      CALL ApplyBoundaryConditions_Radiation( t + dt )

      CALL ComputeRHS_Radiation &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ] )

      CALL ApplyRHS_Radiation &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, alpha = 0.75_DP, beta = 0.25_DP )

      CALL ApplyBoundaryConditions_Radiation &
             ( t + 0.5_DP * dt, LimiterBC_Option = .TRUE. )

      CALL ApplySlopeLimiter_Radiation

      CALL ApplyPositivityLimiter_Radiation

    END IF

    ! -- RK Stage 3 --

    IF( EvolveRadiation )THEN

      CALL ApplyBoundaryConditions_Radiation( t + 0.5_DP * dt )

      CALL ComputeRHS_Radiation &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ] )

      CALL ApplyRHS_Radiation &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, alpha = 1.0_DP / 3.0_DP , beta = 2.0_DP / 3.0_DP )

      CALL ApplyBoundaryConditions_Radiation &
             ( t + dt, LimiterBC_Option = .TRUE. )

      CALL ApplySlopeLimiter_Radiation

      CALL ApplyPositivityLimiter_Radiation

    END IF

    CALL Finalize_SSP_RK

  END SUBROUTINE SSP_RK3


  SUBROUTINE Initialize_SSP_RK

    IF( EvolveRadiation )THEN

      ALLOCATE( uCR_0(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),1:nCR,1:nSpecies) )

      uCR_0(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),1:nCR,1:nSpecies) &
        = uCR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),1:nCR,1:nSpecies)

    END IF

  END SUBROUTINE Initialize_SSP_RK


  SUBROUTINE Finalize_SSP_RK

    IF( EvolveRadiation )THEN

      DEALLOCATE( uCR_0 )

    END IF

  END SUBROUTINE Finalize_SSP_RK


  SUBROUTINE SI_RK1( t, dt )

    REAL(DP), INTENT(in) :: t, dt

    CALL Initialize_SI_RK

    IF( EvolveRadiation )THEN

      CALL ApplyBoundaryConditions_Radiation( t )

      CALL ComputeRHS_Radiation &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ] )

      CALL ApplyRHS_Radiation &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, alpha = 0.0_DP, beta = 1.0_DP )

      CALL ApplyBoundaryConditions_Radiation &
             ( t + dt, LimiterBC_Option = .TRUE. )

      CALL ApplySlopeLimiter_Radiation

      CALL ApplyPositivityLimiter_Radiation

    END IF

    CALL CoupleFluidRadiation &
           ( dt, iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
             EvolveFluid_Option = EvolveFluid )

    CALL Finalize_SI_RK

  END SUBROUTINE SI_RK1


  SUBROUTINE SI_RK2( t, dt )

    REAL(DP), INTENT(in) :: t, dt

    INTEGER :: iE, iX1, iX2, iX3, iCF, iCR, iS
    REAL(DP), DIMENSION(2) :: wTime

    CALL Initialize_SI_RK

    ! -- SI-RK Stage 1 --

    wTime(1) = MPI_WTIME( )

    IF( EvolveRadiation )THEN

      CALL ApplyBoundaryConditions_Radiation( t )

      CALL ComputeRHS_Radiation &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ] )

      CALL ApplyRHS_Radiation &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, alpha = 0.0_DP, beta = 1.0_DP )

      CALL ApplyBoundaryConditions_Radiation &
             ( t + dt, LimiterBC_Option = .TRUE. )

      CALL ApplySlopeLimiter_Radiation

      CALL ApplyPositivityLimiter_Radiation

    END IF

    wTime(2) = MPI_WTIME( )
    IF( DisplayTimers )THEN
      WRITE(*,*)
      WRITE(*,'(A4,A13,ES12.4E2)') '', 'Advection 1: ', wTime(2) - wTime(1)
    END IF

    wTime(1) = MPI_WTIME( )

    CALL CoupleFluidRadiation &
           ( dt, iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
             EvolveFluid_Option = EvolveFluid )

    wTime(2) = MPI_WTIME( )
    IF( DisplayTimers ) &
      WRITE(*,'(A4,A13,ES12.4E2)') '', 'Collision 1: ', wTime(2) - wTime(1)

    IF( EvolveRadiation )THEN

      CALL ApplyPositivityLimiter_Radiation

    END IF

    ! -- SI-RK Stage 2 --

    wTime(1) = MPI_WTIME( )

    IF( EvolveRadiation )THEN

      CALL ApplyBoundaryConditions_Radiation( t + dt )

      CALL ComputeRHS_Radiation &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ] )

      CALL ApplyRHS_Radiation &
             ( iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
               dt = dt, alpha = 0.0_DP, beta = 1.0_DP )

      CALL ApplyBoundaryConditions_Radiation &
             ( t + dt, LimiterBC_Option = .TRUE. )

      CALL ApplySlopeLimiter_Radiation

      CALL ApplyPositivityLimiter_Radiation

    END IF

    wTime(2) = MPI_WTIME( )
    IF( DisplayTimers ) &
      WRITE(*,'(A4,A13,ES12.4E2)') '', 'Advection 2: ', wTime(2) - wTime(1)

    wTime(1) = MPI_WTIME( )

    CALL CoupleFluidRadiation &
           ( dt, iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
             EvolveFluid_Option = EvolveFluid )

    wTime(2) = MPI_WTIME( )
    IF( DisplayTimers ) &
      WRITE(*,'(A4,A13,ES12.4E2)') '', 'Collision 2: ', wTime(2) - wTime(1)

    IF( EvolveRadiation )THEN

      CALL ApplyPositivityLimiter_Radiation

    END IF

    ! -- Combine Steps --

    IF( EvolveFluid )THEN

      DO iCF = 1, nCF
        DO iX3 = 1, nX(3)
          DO iX2 = 1, nX(2)
            DO iX1 = 1, nX(1)

              uCF(:,iX1,iX2,iX3,iCF) &
                = 0.5_DP * ( uCF_0(:,iX1,iX2,iX3,iCF) &
                             + uCF(:,iX1,iX2,iX3,iCF) )

            END DO
          END DO
        END DO
      END DO

    END IF

    IF( EvolveRadiation )THEN

      DO iS = 1, nSpecies
        DO iCR = 1, nCR
          DO iX3 = 1, nX(3)
            DO iX2 = 1, nX(2)
              DO iX1 = 1, nX(1)
                DO iE = 1, nE

                  uCR(:,iE,iX1,iX2,iX3,iCR,iS) &
                    = 0.5_DP * ( uCR_0(:,iE,iX1,iX2,iX3,iCR,iS) &
                                 + uCR(:,iE,iX1,iX2,iX3,iCR,iS) )

                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

    END IF

    CALL Finalize_SI_RK

  END SUBROUTINE SI_RK2


  SUBROUTINE Initialize_SI_RK

    IF( EvolveFluid )THEN

      ALLOCATE( uCF_0(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),1:nCF) )

      uCF_0(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),1:nCF) &
        = uCF(1:nDOFX,1:nX(1),1:nX(2),1:nX(3),1:nCF)

    END IF

    IF( EvolveRadiation )THEN

      ALLOCATE( uCR_0(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),1:nCR,1:nSpecies) )

      uCR_0(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),1:nCR,1:nSpecies) &
        = uCR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),1:nCR,1:nSpecies)

    END IF

  END SUBROUTINE Initialize_SI_RK


  SUBROUTINE Finalize_SI_RK

    IF( EvolveFluid )THEN

      DEALLOCATE( uCF_0 )

    END IF

    IF( EvolveRadiation )THEN

      DEALLOCATE( uCR_0 )

    END IF

  END SUBROUTINE Finalize_SI_RK


  SUBROUTINE ApplyRHS_Radiation( iX_Begin, iX_End, dt, alpha, beta )

    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End
    REAL(DP),              INTENT(in) :: dt, alpha, beta

    INTEGER :: iE, iX1, iX2, iX3, iCR, iS

    DO iS = 1, nSpecies
      DO iCR = 1, nCR
        DO iX3 = iX_Begin(3), iX_End(3)
          DO iX2 = iX_Begin(2), iX_End(2)
            DO iX1 = iX_Begin(1), iX_End(1)
              DO iE = 1, nE

                uCR(:,iE,iX1,iX2,iX3,iCR,iS) &
                  = alpha * uCR_0(:,iE,iX1,iX2,iX3,iCR,iS) &
                      + beta * ( uCR(:,iE,iX1,iX2,iX3,iCR,iS) &
                                 + dt * rhsCR(:,iE,iX1,iX2,iX3,iCR,iS) )

              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE ApplyRHS_Radiation


  SUBROUTINE BackwardEuler( t, dt )

    REAL(DP), INTENT(in) :: t, dt

    CALL CoupleFluidRadiation &
           ( dt, iX_Begin = [ 1, 1, 1 ], iX_End = [ nX(1), nX(2), nX(3) ], &
             EvolveFluid_Option = EvolveFluid )

  END SUBROUTINE BackwardEuler


END MODULE TimeSteppingModule
