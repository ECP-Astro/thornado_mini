MODULE FluidRadiationCouplingSolutionModule_EmAbScatt

  USE KindModule, ONLY: &
    DP, FourPi
  USE UnitsModule, ONLY: &
    SpeedOfLight, &
    BoltzmannConstant, &
    PlanckConstant, &
    Gram, &
    Centimeter, &
    Kelvin, &
    MeV
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, nDOFX, &
    nE, nNodesE, nDOFE
  USE FluidFieldsModule, ONLY: &
    uCF, nCF, &
    uPF, nPF, iPF_D, iPF_E, iPF_Ne, &
    uAF, nAF, iAF_P, iAF_T, iAF_Ye, iAF_S, iAF_E, &
    iAF_Me, iAF_Mp, iAF_Mn, iAF_Xp, iAF_Xn, iAF_Xa, &
    iAF_Xh, iAF_Gm, iAF_Cs
  USE EulerEquationsUtilitiesModule, ONLY: &
    ComputeConserved, &
    ComputePrimitive
  USE RadiationFieldsModule, ONLY: &
    nPR, iPR_D, iPR_I1, iPR_I2, iPR_I3
  USE MomentEquationsUtilitiesModule, ONLY: &
    ComputeConservedMoments, &
    ComputePrimitiveMoments
  USE FluidRadiationCouplingUtilitiesModule, ONLY: &
    InitializeNodes, &
    InitializeNodesX, &
    InitializeWeights, &
    InitializeFluidFields, &
    InitializeRadiationFields, &
    FinalizeFluidFields, &
    FinalizeRadiationFields, &
    FermiDirac, &
    dFermiDiracdT, &
    dFermiDiracdY, &
    ENORM
  USE EquationOfStateModule, ONLY: &
    BaryonMass, &
    ComputeTemperatureFromSpecificInternalEnergy, &
    ComputeSpecificInternalEnergy, &
    ComputeElectronChemicalPotential, &
    ComputeProtonChemicalPotential, &
    ComputeNeutronChemicalPotential, &
    ComputeAuxiliary_Fluid, &
    ApplyEquationOfState, &
    ComputeThermodynamicStates_Primitive
  USE OpacityModule, ONLY: &
    ComputeAbsorptionOpacity, &
    ComputeScatteringOpacity_ES

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  LOGICAL            :: EvolveFluid
  INTEGER            :: nNodesX_G, nNodesE_G
  INTEGER, PARAMETER :: i_Y = 1, i_E = 2
  INTEGER, PARAMETER :: iOld = 0, iNew = 1
  REAL(DP), DIMENSION(:),     ALLOCATABLE :: E_N, W2_N, W3_N
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: X_N, uPF_N, uAF_N
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: duAFdY_N, duAFdT_N
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: Chi, FD, Sigma
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: dFDdT_Y, dFDdY_T
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: dFDdE_Y, dFDdY_E
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: FVEC, dUVEC
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: uPR_N
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: UVEC, FJAC

  PUBLIC :: CoupleFluidRadiation_EmAbScatt

  LOGICAL, PARAMETER :: DisplayTimers = .FALSE.
  REAL(DP) :: Timer_Total
  REAL(DP) :: Timer_Primitive
  REAL(DP) :: Timer_Initialize
  REAL(DP) :: Timer_Initialize_Fluid
  REAL(DP) :: Timer_Finalize
  REAL(DP) :: Timer_Finalize_Fluid
  REAL(DP) :: Timer_Conserved
  REAL(DP) :: Timer_Compute
  REAL(DP) :: Timer_Couple
  REAL(DP) :: Timer_Rates
  REAL(DP) :: Timer_Linear, dT_Linear
  REAL(DP) :: Timer_Specific, dT_Specific
  REAL(DP) :: Timer_Temp, dT_Temp
  REAL(DP) :: Timer_Eq, dT_Eq
  REAL(DP) :: Timer_FJAC, dT_FJAC
  REAL(DP) :: Timer_FD, dT_FD

CONTAINS


  SUBROUTINE CoupleFluidRadiation_EmAbScatt &
               ( dt, iX_Begin, iX_End, EvolveFluid_Option )

    REAL(DP),              INTENT(in) :: dt
    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End
    LOGICAL,               INTENT(in), OPTIONAL :: EvolveFluid_Option

    INTEGER :: iX1, iX2, iX3

    CALL Timer_Start( Timer_Total )

    EvolveFluid = .TRUE.
    IF( PRESENT( EvolveFluid_Option ) ) &
      EvolveFluid = EvolveFluid_Option

    CALL Timer_Start( Timer_Primitive )

    CALL ComputePrimitiveMoments &
           ( iX_Begin = iX_Begin, iX_End = iX_End )

    CALL Timer_Stop( Timer_Primitive )

    CALL Timer_Start( Timer_Initialize_Fluid )

    DO iX3 = iX_Begin(3), iX_End(3)
      DO iX2 = iX_Begin(2), iX_End(2)
        DO iX1 = iX_Begin(1), iX_End(1)

          CALL ComputePrimitive &
                 ( uCF(:,iX1,iX2,iX3,1:nCF), uPF(:,iX1,iX2,iX3,1:nPF) )

          CALL ComputeAuxiliary_Fluid &
                 ( uPF(:,iX1,iX2,iX3,iPF_D ), uPF(:,iX1,iX2,iX3,iPF_E ), &
                   uPF(:,iX1,iX2,iX3,iPF_Ne), uAF(:,iX1,iX2,iX3,iAF_P ), &
                   uAF(:,iX1,iX2,iX3,iAF_T ), uAF(:,iX1,iX2,iX3,iAF_Ye), &
                   uAF(:,iX1,iX2,iX3,iAF_S ), uAF(:,iX1,iX2,iX3,iAF_E ), &
                   uAF(:,iX1,iX2,iX3,iAF_Gm), uAF(:,iX1,iX2,iX3,iAF_Cs) )

        END DO
      END DO
    END DO

    CALL Timer_Stop( Timer_Initialize_Fluid )

    CALL Timer_Start( Timer_Initialize )

    CALL InitializeFluidRadiationCoupling

    CALL Timer_Stop( Timer_Initialize )

    CALL Timer_Start( Timer_Compute )

    CALL CoupleFluidRadiation( dt )

    CALL Timer_Stop( Timer_Compute )

    CALL Timer_Start( Timer_Finalize )

    CALL FinalizeFluidRadiationCoupling

    CALL Timer_Stop( Timer_Finalize )

    CALL Timer_Start( Timer_Finalize_Fluid )

    DO iX3 = iX_Begin(3), iX_End(3)
      DO iX2 = iX_Begin(2), iX_End(2)
        DO iX1 = iX_Begin(1), iX_End(1)

          CALL ApplyEquationOfState &
                 ( uPF(1:nDOFX,iX1,iX2,iX3,iPF_D ), &
                   uAF(1:nDOFX,iX1,iX2,iX3,iAF_T ), &
                   uAF(1:nDOFX,iX1,iX2,iX3,iAF_Ye), &
                   uAF(1:nDOFX,iX1,iX2,iX3,iAF_P ), &
                   uAF(1:nDOFX,iX1,iX2,iX3,iAF_S ), &
                   uAF(1:nDOFX,iX1,iX2,iX3,iAF_E ), &
                   uAF(1:nDOFX,iX1,iX2,iX3,iAF_Me), &
                   uAF(1:nDOFX,iX1,iX2,iX3,iAF_Mp), &
                   uAF(1:nDOFX,iX1,iX2,iX3,iAF_Mn), &
                   uAF(1:nDOFX,iX1,iX2,iX3,iAF_Xp), &
                   uAF(1:nDOFX,iX1,iX2,iX3,iAF_Xn), &
                   uAF(1:nDOFX,iX1,iX2,iX3,iAF_Xa), &
                   uAF(1:nDOFX,iX1,iX2,iX3,iAF_Xh), &
                   uAF(1:nDOFX,iX1,iX2,iX3,iAF_Gm) )

          CALL ComputeThermodynamicStates_Primitive &
                 ( uPF(1:nDOFX,iX1,iX2,iX3,iPF_D ), &
                   uAF(1:nDOFX,iX1,iX2,iX3,iAF_T ), &
                   uAF(1:nDOFX,iX1,iX2,iX3,iAF_Ye), &
                   uPF(1:nDOFX,iX1,iX2,iX3,iPF_E ), &
                   uAF(1:nDOFX,iX1,iX2,iX3,iAF_E ), &
                   uPF(1:nDOFX,iX1,iX2,iX3,iPF_Ne) )

          CALL ComputeConserved &
                 ( uPF(:,iX1,iX2,iX3,1:nPF), uCF(:,iX1,iX2,iX3,1:nCF) )

        END DO
      END DO
    END DO

    CALL Timer_Stop( Timer_Finalize_Fluid )

    CALL Timer_Start( Timer_Conserved )

    CALL ComputeConservedMoments &
           ( iX_Begin = iX_Begin, iX_End = iX_End )

    CALL Timer_Stop( Timer_Conserved )

    CALL Timer_Stop( Timer_Total )

    IF( DisplayTimers )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A30,ES10.4E2)') &
        '', 'Total: ', Timer_Total
      WRITE(*,*)
      WRITE(*,'(A4,A30,ES10.4E2)') &
        '', 'Primitive: ', Timer_Primitive
      WRITE(*,'(A4,A30,ES10.4E2)') &
        '', 'Initialize: ', Timer_Initialize
      WRITE(*,'(A4,A30,ES10.4E2)') &
        '', 'Initialize Fluid: ', Timer_Initialize_Fluid
      WRITE(*,'(A4,A30,ES10.4E2)') &
        '', 'Compute: ', Timer_Compute
      WRITE(*,'(A4,A30,ES10.4E2)') &
        '', 'Finalize: ', Timer_Finalize
      WRITE(*,'(A4,A30,ES10.4E2)') &
        '', 'Finalize Fluid: ', Timer_Finalize_Fluid
      WRITE(*,'(A4,A30,ES10.4E2)') &
        '', 'Conserved: ', Timer_Conserved
      WRITE(*,*)
      WRITE(*,'(A4,A30,ES10.4E2)') &
        '', 'Sum: ', Timer_Primitive + Timer_Initialize &
        + Timer_Initialize_Fluid + Timer_Compute + Timer_Finalize &
        + Timer_Finalize_Fluid + Timer_Conserved
      WRITE(*,*)

    END IF

  END SUBROUTINE CoupleFluidRadiation_EmAbScatt


  SUBROUTINE InitializeFluidRadiationCoupling

    nNodesX_G = PRODUCT(nX) * nDOFX
    nNodesE_G =         nE  * nDOFE

    ALLOCATE( E_N(nNodesE_G) )
    CALL InitializeNodes( E_N )

    ALLOCATE( W2_N(nNodesE_G), W3_N(nNodesE_G) )
    CALL InitializeWeights( W2_N, W3_N )

    ALLOCATE( X_N(nNodesX_G,3) )
    CALL InitializeNodesX( X_N )

    ALLOCATE( uPF_N(nPF, nNodesX_G) )
    ALLOCATE( uAF_N(nAF, nNodesX_G) )
    CALL InitializeFluidFields( uPF_N, uAF_N )

    ALLOCATE( duAFdT_N(nAF, nNodesX_G) )
    ALLOCATE( duAFdY_N(nAF, nNodesX_G) )

    ALLOCATE( uPR_N(nNodesE_G, nPR, nNodesX_G) )
    CALL InitializeRadiationFields( uPR_N )

    ALLOCATE &
      ( Chi  (nNodesE_G, nNodesX_G), &
        Sigma(nNodesE_G, nNodesX_G) )

    ALLOCATE &
      ( FD     (nNodesE_G, nNodesX_G), & ! --- Fermi-Dirac (FD) Distribution
        dFDdT_Y(nNodesE_G, nNodesX_G), & ! --- FD Deriv. wrt. T  (Const. Ye)
        dFDdY_T(nNodesE_G, nNodesX_G), & ! --- FD Deriv. wrt. Ye (Const. T )
        dFDdE_Y(nNodesE_G, nNodesX_G), & ! --- FD Deriv. wrt. E  (Const. Ye)
        dFDdY_E(nNodesE_G, nNodesX_G) )  ! --- FD Deriv. wrt. Ye (Const. E )

    ALLOCATE( UVEC (2,    nNodesX_G, 0:1) )
    ALLOCATE( FVEC (2,    nNodesX_G) )
    ALLOCATE( dUVEC(2,    nNodesX_G) )
    ALLOCATE( FJAC (2, 2, nNodesX_G) )

  END SUBROUTINE InitializeFluidRadiationCoupling


  SUBROUTINE FinalizeFluidRadiationCoupling

    CALL FinalizeFluidFields( uPF_N, uAF_N )

    CALL FinalizeRadiationFields( uPR_N )

    DEALLOCATE( E_N, W2_N, W3_N )
    DEALLOCATE( X_N, uPF_N, uPR_N )
    DEALLOCATE( uAF_N, duAFdT_N, duAFdY_N )
    DEALLOCATE( Chi, FD, Sigma )
    DEALLOCATE( dFDdT_Y, dFDdY_T )
    DEALLOCATE( dFDdE_Y, dFDdY_E )

    DEALLOCATE( UVEC, FVEC, dUVEC, FJAC )

  END SUBROUTINE FinalizeFluidRadiationCoupling


  SUBROUTINE CoupleFluidRadiation( dt )

    REAL(DP), INTENT(in) :: dt

    LOGICAL,  DIMENSION(1:nNodesX_G) :: Converged
    INTEGER                          :: iter, iX, iX_MAX
    REAL(DP)                         :: Norm, MaxNorm
    REAL(DP),              PARAMETER :: NewtonTol = 1.0d-10
    REAL(DP), DIMENSION(1:nNodesX_G) :: Work1_N, Work2_N, Work3_N

    Timer_Linear   = 0.0_DP
    Timer_Specific = 0.0_DP
    Timer_Temp     = 0.0_DP
    Timer_Eq       = 0.0_DP
    Timer_FJAC     = 0.0_DP
    Timer_FD       = 0.0_DP

    CALL Timer_Start( Timer_Couple )

    CALL SetUVEC( uAF_N(iAF_Ye,:), uAF_N(iAF_E,:), iOld )

    ! --- Compute Absorption Opacities ---

    CALL Timer_Start( Timer_Rates )

    CALL SetRates

    CALL Timer_Stop( Timer_Rates )

    Converged  = .FALSE.
    iter       = 0

    DO WHILE( .NOT. ALL( Converged ) )

      iter = iter + 1

      CALL SetUVEC( uAF_N(iAF_Ye,:), uAF_N(iAF_E,:), iNew )

      CALL Timer_Start( dT_Specific )

      CALL ComputeSpecificInternalEnergy & ! --- Computes Derivatives
             ( uPF_N(iPF_D, :), uAF_N(iAF_T,:), uAF_N(iAF_Ye,:), &
               Work1_N, dVdY = Work2_N, dVdZ = Work3_N, &
               Mask_Option = .NOT. Converged )

      WHERE( .NOT. Converged )

        uAF_N   (iAF_E,:) = Work1_N
        duAFdT_N(iAF_E,:) = Work2_N
        duAFdY_N(iAF_E,:) = Work3_N

      END WHERE

      CALL Timer_Stop( dT_Specific )

      CALL Timer_Add( Timer_Specific, dT_Specific )

      ! --- Compute FD Distribution and Derivatives ---

      CALL Timer_Start( dT_Eq )

      CALL SetEquilibrium( Converged )

      CALL Timer_Stop( dT_Eq )

      CALL Timer_Add( Timer_Eq, dT_Eq )

      ! --- Set Equations Vector ---

      CALL SetFVEC( dt, Converged )

      ! --- Set Jacobian Matrix ---

      CALL Timer_Start( dT_FJAC )

      CALL SetFJAC( dt, Converged )

      CALL Timer_Stop( dT_FJAC )

      CALL Timer_Add( Timer_FJAC, dT_FJAC )

      ! --- Invert for Correction ---

      CALL Timer_Start( dT_Linear )

      CALL SolveLinearSystems( Converged )

      CALL Timer_Stop( dT_Linear )

      CALL Timer_Add( Timer_Linear, dT_Linear )

      UVEC(:,:,iNew) = UVEC(:,:,iNew) + dUVEC(:,:)

      CALL GetUVEC( uAF_N(iAF_Ye,:), uAF_N(iAF_E,:), iNew )

      CALL Timer_Start( dT_Temp )

      CALL ComputeTemperatureFromSpecificInternalEnergy &
             ( uPF_N(iPF_D,:), uAF_N(iAF_E,:), uAF_N(iAF_Ye,:), &
               Work1_N, Mask_Option = .NOT. Converged )

      WHERE( .NOT. Converged )

        uAF_N(iAF_T,:) = Work1_N

      END WHERE

      CALL Timer_Stop( dT_Temp )

      CALL Timer_Add( Timer_Temp, dT_Temp )

      MaxNorm = 0.0_DP
      DO iX = 1, nNodesX_G

        Norm = ENORM( dUVEC(:,iX) / UVEC(:,iX,iNew) )
        IF( Norm <= NewtonTol ) Converged(iX) = .TRUE.

        IF( Norm >= MaxNorm )THEN

          MaxNorm = Norm
          iX_MAX  = iX

        END IF

      END DO

      IF( iter > 20 )THEN

        WRITE(*,*)
        WRITE(*,'(A8,A)') ' ', 'Emission/Absorption'
        WRITE(*,*)
        WRITE(*,'(A10,A12,I6.6,A2,A11,ES10.4E2,A2,A9,I6.6)') &
          ' ', 'Iteration = ', iter, &
          ', ', '||dU/U|| = ', MaxNorm, &
          ', ', 'iX_MAX = ', iX_MAX
        WRITE(*,*)
        WRITE(*,'(A12,I4.4,A19,I4.4)') &
          '', COUNT( Converged .EQV. .TRUE. ), &
          ' converged, out of ', SIZE( Converged )
        WRITE(*,*)
        WRITE(*,'(A12,A22,ES10.4E2,A3,ES10.4E2)') &
          '', 'MIN/MAX    ||dY/Y|| = ', &
          MINVAL( ABS( dUVEC(i_Y,:) / UVEC(i_Y,:,iNew) ) ),' / ', &
          MAXVAL( ABS( dUVEC(i_Y,:) / UVEC(i_Y,:,iNew) ) )
        WRITE(*,'(A12,A22,ES10.4E2,A3,ES10.4E2)') &
          '', 'MIN/MAX ||FVEC(Y)|| = ', &
          MINVAL( ABS( FVEC(i_Y,:) ) ), ' / ', &
          MAXVAL( ABS( FVEC(i_Y,:) ) )
        WRITE(*,*)
        WRITE(*,'(A12,A22,ES10.4E2,A3,ES10.4E2)') &
          '', 'MIN/MAX    ||dE/E|| = ', &
          MINVAL( ABS( dUVEC(i_E,:) / UVEC(i_E,:,iNew) ) ),' / ', &
          MAXVAL( ABS( dUVEC(i_E,:) / UVEC(i_E,:,iNew) ) )
        WRITE(*,'(A12,A22,ES10.4E2,A3,ES10.4E2)') &
          '', 'MIN/MAX ||FVEC(E)|| = ', &
          MINVAL( ABS( FVEC(i_E,:) ) ), ' / ', &
          MAXVAL( ABS( FVEC(i_E,:) ) )
        WRITE(*,*)
        WRITE(*,'(A12,A4,ES10.4E2,A2,A4,ES10.4E2,A2,A4,ES10.4E2)') &
          '', 'D = ', uPF_N(iPF_D, iX_MAX) / ( Gram / Centimeter**3 ), &
          '', 'T = ', uAF_N(iAF_T, iX_MAX) / Kelvin, &
          '', 'Y = ', uAF_N(iAF_Ye,iX_MAX)
        WRITE(*,*)

      END IF

    END DO

    CALL UpdateRadiationFields( dt )

    CALL Timer_Stop( Timer_Couple )

    IF( DisplayTimers )THEN

      WRITE(*,*)
      WRITE(*,'(A4,A)') &
        '', 'Timers:'
      WRITE(*,*)
      WRITE(*,'(A4,A30,ES10.4E2)') &
        '', 'CoupleFluidRadiation: ', Timer_Couple
      WRITE(*,*)
      WRITE(*,'(A4,A30,ES10.4E2)') &
        '', 'SetRates: ', Timer_Rates
      WRITE(*,'(A4,A30,ES10.4E2)') &
        '', 'SolveLinearSystem: ', Timer_Linear
      WRITE(*,'(A4,A30,ES10.4E2)') &
        '', 'ComputeSpecificEnergy: ', Timer_Specific
      WRITE(*,'(A4,A30,ES10.4E2)') &
        '', 'ComputeTemperature: ', Timer_Temp
      WRITE(*,'(A4,A30,ES10.4E2)') &
        '', 'SetEquilibrium: ', Timer_Eq
      WRITE(*,'(A4,A30,ES10.4E2)') &
        '', 'SetFJAC: ', Timer_FJAC
      WRITE(*,'(A4,A30,ES10.4E2)') &
        '', 'FD: ', Timer_FD
      WRITE(*,*)
      WRITE(*,'(A4,A30,ES10.4E2)') &
        '', 'Sum: ', Timer_Rates+Timer_Linear+Timer_Specific+Timer_Temp &
        +Timer_Eq+Timer_FJAC
      WRITE(*,*)

    END IF

  END SUBROUTINE CoupleFluidRadiation


  SUBROUTINE SetUVEC( Y, E, iState )

    REAL(DP), DIMENSION(nNodesX_G), INTENT(in) :: Y, E
    INTEGER,                        INTENT(in) :: iState

    INTEGER :: iX

    DO iX = 1, nNodesX_G
      UVEC(i_Y,iX,iState) = Y(iX)
      UVEC(i_E,iX,iState) = E(iX)
    END DO

  END SUBROUTINE SetUVEC


  SUBROUTINE GetUVEC( Y, E, iState )

    REAL(DP), DIMENSION(nNodesX_G), INTENT(out) :: Y, E
    INTEGER,                        INTENT(in)  :: iState

    INTEGER :: iX, iE

    DO iX = 1, nNodesX_G
      Y(iX) = UVEC(i_Y,iX,iState)
      E(iX) = UVEC(i_E,iX,iState)
    END DO

  END SUBROUTINE GetUVEC


  SUBROUTINE SetRates

    ASSOCIATE &
      ( D_N => uPF_N(iPF_D, 1:nNodesX_G), &
        T_N => uAF_N(iAF_T, 1:nNodesX_G), &
        Y_N => uAF_N(iAF_Ye,1:nNodesX_G) )

    CALL ComputeAbsorptionOpacity &
           ( E_N, D_N, T_N, Y_N, X_N(:,1), X_N(:,2), X_N(:,3), Chi )

    CALL ComputeScatteringOpacity_ES &
           ( E_N, D_N, T_N, Y_N, X_N(:,1), X_N(:,2), X_N(:,3), Sigma )

    END ASSOCIATE ! D_N, etc.

  END SUBROUTINE SetRates


  SUBROUTINE SetEquilibrium( Converged )

    LOGICAL, DIMENSION(:), INTENT(in) :: Converged

    INTEGER  :: iX
    REAL(DP) :: Mnu, dMnudT, dMnudY, kB = BoltzmannConstant
    REAL(DP), DIMENSION(nNodesX_G) :: Work1_N, Work2_N, Work3_N

    CALL ComputeElectronChemicalPotential &
           ( uPF_N(iPF_D, :), uAF_N(iAF_T, :), uAF_N(iAF_Ye,:), &
             Work1_N, dVdY = Work2_N, dVdZ = Work3_N, &
             Mask_Option = .NOT. Converged )

    WHERE( .NOT. Converged )

      uAF_N   (iAF_Me,:) = Work1_N
      duAFdT_N(iAF_Me,:) = Work2_N
      duAFdY_N(iAF_Me,:) = Work3_N

    END WHERE

    CALL ComputeProtonChemicalPotential &
           ( uPF_N(iPF_D, :), uAF_N(iAF_T, :), uAF_N(iAF_Ye,:), &
             Work1_N, dVdY = Work2_N, dVdZ = Work3_N, &
             Mask_Option = .NOT. Converged )

    WHERE( .NOT. Converged )

      uAF_N   (iAF_Mp,:) = Work1_N
      duAFdT_N(iAF_Mp,:) = Work2_N
      duAFdY_N(iAF_Mp,:) = Work3_N

    END WHERE

    CALL ComputeNeutronChemicalPotential &
           ( uPF_N(iPF_D, :), uAF_N(iAF_T, :), uAF_N(iAF_Ye,:), &
             Work1_N, dVdY = Work2_N, dVdZ = Work3_N, &
             Mask_Option = .NOT. Converged )

    WHERE( .NOT. Converged )

      uAF_N   (iAF_Mn,:) = Work1_N
      duAFdT_N(iAF_Mn,:) = Work2_N
      duAFdY_N(iAF_Mn,:) = Work3_N

    END WHERE

    CALL Timer_Start( dT_FD )

    !$OMP PARALLEL DO PRIVATE &
    !$OMP&              ( iX, Mnu, dMnudT, dMnudY )
    DO iX = 1, nNodesX_G

      IF( Converged(iX) ) CYCLE

      Mnu    = uAF_N   (iAF_Me,iX) + uAF_N   (iAF_Mp,iX) - uAF_N   (iAF_Mn,iX)
      dMnudT = duAFdT_N(iAF_Me,iX) + duAFdT_N(iAF_Mp,iX) - duAFdT_N(iAF_Mn,iX)
      dMnudY = duAFdY_N(iAF_Me,iX) + duAFdY_N(iAF_Mp,iX) - duAFdY_N(iAF_Mn,iX)

      FD(:,iX) &
        = FermiDirac &
            ( E_N, Mnu, kB * uAF_N(iAF_T,iX) )

      ! --- Derivatives w.r.t. T and E ---

      dFDdT_Y(:,iX) &
        = dFermiDiracdT &
            ( E_N, Mnu, kB * uAF_N(iAF_T,iX), dMnudT, uAF_N(iAF_T,iX) )

      dFDdE_Y(:,iX) &
        = dFDdT_Y(:,iX) / duAFdT_N(iAF_E,iX)

      ! --- Derivatives w.r.t. Ye ---

      dFDdY_T(:,iX) &
        = dFermiDiracdY &
            ( E_N, Mnu, kB * uAF_N(iAF_T,iX), dMnudY, uAF_N(iAF_T,iX) )

      dFDdY_E(:,iX) &
        = dFDdY_T(:,iX) &
            - dFDdT_Y(:,iX) * duAFdY_N(iAF_E,iX) / duAFdT_N(iAF_E,iX)

    END DO
    !$OMP END PARALLEL DO

    CALL Timer_Stop( dT_FD )

    CALL Timer_Add( Timer_FD, dT_FD )

  END SUBROUTINE SetEquilibrium


  SUBROUTINE SetFVEC( dt, Converged )

    REAL(DP),              INTENT(in) :: dt
    LOGICAL, DIMENSION(:), INTENT(in) :: Converged

    INTEGER :: iX, iE
    REAL(DP), DIMENSION(1:nNodesE_G) :: GammaT

    ASSOCIATE &
      ( hc3 => ( PlanckConstant * SpeedOfLight )**3 )

    ASSOCIATE &
      ( mB  => BaryonMass, &
        D_N => uPF_N(iPF_D,:) )

    FVEC = 0.0_DP

    DO iX = 1, nNodesX_G

      IF( Converged(iX) ) CYCLE

      GammaT(1:nNodesE_G) &
        = dt * Chi(1:nNodesE_G,iX) &
            / ( 1.0_DP + dt * Chi(1:nNodesE_G,iX) )

      ! --- Electron Fraction Equation ---

      FVEC(i_Y,iX) &
        = ( UVEC(i_Y,iX,iNew) - UVEC(i_Y,iX,iOld) )

      IF( EvolveFluid )THEN

        FVEC(i_Y,iX) &
          = FVEC(i_Y,iX) &
              + ( mB / D_N(iX) ) / hc3 &
                  * SUM( W2_N(:) * GammaT(:) &
                           * ( FourPi * FD(:,iX) - uPR_N(:,iPR_D,iX) ) )

      END IF

      ! --- Energy Equation ---

      FVEC(i_E,iX) &
        = ( UVEC(i_E,iX,iNew) - UVEC(i_E,iX,iOld) )

      IF( EvolveFluid )THEN

        FVEC(i_E,iX) &
          = FVEC(i_E,iX) &
              + ( 1.0_DP / D_N(iX) ) / hc3 &
                  * SUM( W3_N(:) * GammaT(:) &
                           * ( FourPi * FD(:,iX) - uPR_N(:,iPR_D,iX) ) )

      END IF

    END DO

    END ASSOCIATE ! mB, D_N

    END ASSOCIATE ! hc3

  END SUBROUTINE SetFVEC


  SUBROUTINE SetFJAC( dt, Converged )

    REAL(DP),              INTENT(in) :: dt
    LOGICAL, DIMENSION(:), INTENT(in) :: Converged

    INTEGER :: iX, iE
    REAL(DP), DIMENSION(1:nNodesE_G) :: GammaT

    ASSOCIATE &
      ( hc3 => ( PlanckConstant * SpeedOfLight )**3 )

    ASSOCIATE &
      ( mB   => BaryonMass, &
        D    => uPF_N(iPF_D,:) )

    FJAC = 0.0_DP

    DO iX = 1, nNodesX_G

      IF( Converged(iX) ) CYCLE

      GammaT(1:nNodesE_G) &
        = dt * Chi(1:nNodesE_G,iX) &
            / ( 1.0_DP + dt * Chi(1:nNodesE_G,iX) )

      ! --- Electron Density Equation ---

      FJAC(i_Y,i_Y,iX) &
        = 1.0_DP

      IF( EvolveFluid )THEN

        FJAC(i_Y,i_Y,iX) &
          = FJAC(i_Y,i_Y,iX) &
              + ( mB / D(iX) ) / hc3 &
                  * SUM( W2_N(:) * GammaT(:) * FourPi * dFDdY_E(:,iX) )

        FJAC(i_Y,i_E,iX) &
          =     ( mB / D(iX) ) / hc3 &
                  * SUM( W2_N(:) * GammaT(:) * FourPi * dFDdE_Y(:,iX) )

      END IF

      ! --- Internal Energy Equation ---

      FJAC(i_E,i_E,iX) &
        = 1.0_DP

      IF( EvolveFluid )THEN

        FJAC(i_E,i_Y,iX) &
          =     ( 1.0_DP / D(iX) ) / hc3 &
                  * SUM( W3_N(:) * GammaT(:) * FourPi * dFDdY_E(:,iX) )

        FJAC(i_E,i_E,iX) &
          = FJAC(i_E,i_E,iX) &
              + ( 1.0_DP / D(iX) ) / hc3 &
                  * SUM( W3_N(:) * GammaT(:) * FourPi * dFDdE_Y(:,iX) )

      END IF

    END DO

    END ASSOCIATE ! mB, etc.

    END ASSOCIATE ! hc3

  END SUBROUTINE SetFJAC


  SUBROUTINE SolveLinearSystems( Converged )

    LOGICAL, DIMENSION(:), INTENT(in) :: Converged

    INTEGER                  :: iX
    REAL(DP)                 :: detFJAC
    REAL(DP), DIMENSION(2,2) :: invFJAC

    dUVEC = 0.0_DP

    DO iX = 1, nNodesX_G

      IF( Converged(iX) ) CYCLE

      detFJAC = FJAC(i_Y,i_Y,iX) * FJAC(i_E,i_E,iX) &
                  - FJAC(i_E,i_Y,iX) * FJAC(i_Y,i_E,iX)

      invFJAC(i_Y,i_Y) =   FJAC(i_E,i_E,iX)
      invFJAC(i_Y,i_E) = - FJAC(i_Y,i_E,iX)
      invFJAC(i_E,i_Y) = - FJAC(i_E,i_Y,iX)
      invFJAC(i_E,i_E) =   FJAC(i_Y,i_Y,iX)

      invFJAC = invFJAC / detFJAC

      dUVEC(:,iX) = - MATMUL( invFJAC, FVEC(:,iX) )

    END DO

  END SUBROUTINE SolveLinearSystems


  SUBROUTINE UpdateRadiationFields( dt )

    REAL(DP), INTENT(in) :: dt

    INTEGER  :: iX, iE
    REAL(DP) :: Gamma, Gamma_T

    DO iX = 1, nNodesX_G
      DO iE = 1, nNodesE_G

        Gamma   = dt * Chi(iE,iX)
        Gamma_T = Gamma + dt * Sigma(iE,iX)

        uPR_N(iE,iPR_D,iX) &
          = ( uPR_N(iE,iPR_D,iX) &
              + Gamma * FourPi * FD(iE,iX) ) / ( 1.0_DP + Gamma )

        uPR_N(iE,iPR_I1,iX) &
          = uPR_N(iE,iPR_I1,iX) / ( 1.0_DP + Gamma_T )

        uPR_N(iE,iPR_I2,iX) &
          = uPR_N(iE,iPR_I2,iX) / ( 1.0_DP + Gamma_T )

        uPR_N(iE,iPR_I3,iX) &
          = uPR_N(iE,iPR_I3,iX) / ( 1.0_DP + Gamma_T )

      END DO
    END DO

  END SUBROUTINE UpdateRadiationFields


  SUBROUTINE Timer_Start( Timer )

    REAL(DP) :: Timer

    IF( .NOT. DisplayTimers ) RETURN

    Timer = MPI_WTIME( )

  END SUBROUTINE Timer_Start


  SUBROUTINE Timer_Stop( Timer )

    REAL(DP) :: Timer

    IF( .NOT. DisplayTimers ) RETURN

    Timer = MPI_WTIME( ) - Timer

  END SUBROUTINE Timer_Stop


  SUBROUTINE Timer_Add( Timer, dT )

    REAL(DP) :: Timer, dT

    IF( .NOT. DisplayTimers ) RETURN

    Timer = Timer + dT

  END SUBROUTINE Timer_Add


END MODULE FluidRadiationCouplingSolutionModule_EmAbScatt
