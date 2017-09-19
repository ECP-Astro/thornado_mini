MODULE OpacityModule_TABLE

  ! --- weaklib modules --------------------------

  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF
  USE wlOpacityTableIOModuleHDF, ONLY: &
    ReadOpacityTableHDF
  USE wlOpacityTableModule, ONLY: &
    OpacityTableType
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable_1D3D, &
    LogInterpolateSingleVariable_2D2D

  ! ----------------------------------------------

  USE KindModule, ONLY: &
    DP
  USE UnitsModule, ONLY: &
    Gram, &
    Centimeter, &
    Kelvin, &
    MeV

  IMPLICIT NONE
  PRIVATE

  CHARACTER(256) :: &
    OpacityTableName
  INTEGER :: &
    iD_T, iT_T, iY_T
  REAL(DP), DIMENSION(:), ALLOCATABLE :: &
    Es_T, Ds_T, Ts_T, Ys_T, Etas_T, &
    LogEs_T, LogDs_T, LogTs_T, LogEtas_T
  TYPE(OpacityTableType) :: &
    OPACITIES

  PUBLIC :: InitializeOpacities_TABLE
  PUBLIC :: FinalizeOpacities_TABLE
  PUBLIC :: ComputeAbsorptionOpacity_TABLE
  PUBLIC :: ComputeScatteringOpacity_ES_TABLE
  PUBLIC :: ComputeScatteringOpacity_NES_TABLE

CONTAINS


  SUBROUTINE InitializeOpacities_TABLE( OpacityTableName_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: OpacityTableName_Option

    OpacityTableName = 'OpacityTable.h5'
    IF( PRESENT( OpacityTableName_Option ) ) &
      OpacityTableName = TRIM( OpacityTableName_Option )

    WRITE(*,*)
    WRITE(*,'(A7,A12,A)') &
      '', 'Table Name: ', TRIM( OpacityTableName )

    CALL InitializeHDF( )

    CALL ReadOpacityTableHDF &
           ( OPACITIES, TRIM( OpacityTableName ) )

    CALL FinalizeHDF( )

    ! --- Thermodynamic State Indices ---

    iD_T = OPACITIES % TS % Indices % iRho
    iT_T = OPACITIES % TS % Indices % iT
    iY_T = OPACITIES % TS % Indices % iYe

    ! --- Thermodynamic States ---

    ! --- Density ---

    ALLOCATE( Ds_T(OPACITIES % TS % nPoints(iD_T)) )
    Ds_T = OPACITIES % TS % States(iD_T) % Values

    ALLOCATE( LogDs_T(SIZE( Ds_T )) )
    LogDs_T = LOG10( Ds_T )

    ! --- Temperature ---

    ALLOCATE( Ts_T(OPACITIES % TS % nPoints(iT_T)) )
    Ts_T = OPACITIES % TS % States(iT_T) % Values

    ALLOCATE( LogTs_T(SIZE( Ts_T )) )
    LogTs_T = LOG10( Ts_T )

    ! --- Electron Fraction ---

    ALLOCATE( Ys_T(OPACITIES % TS % nPoints(iY_T)) )
    Ys_T = OPACITIES % TS % States(iY_T) % Values

    ! --- Energy Grid ---

    ALLOCATE( Es_T(OPACITIES % EnergyGrid % nPoints) )
    Es_T = OPACITIES % EnergyGrid  % Values

    ALLOCATE( LogEs_T(SIZE( Es_T )) )
    LogEs_T = LOG10( Es_T )

    ! --- Eta Grid ---

    ALLOCATE( Etas_T(OPACITIES % EtaGrid % nPoints) )
    Etas_T = OPACITIES % EtaGrid  % Values

    ALLOCATE( LogEtas_T(SIZE( Etas_T )) )
    LogEtas_T = LOG10( Etas_T )

  END SUBROUTINE InitializeOpacities_TABLE


  SUBROUTINE FinalizeOpacities_TABLE

    DEALLOCATE( Es_T, Ds_T, Ts_T, Ys_T, Etas_T )
    DEALLOCATE( LogEs_T, LogDs_T, LogTs_T, LogEtas_T )

  END SUBROUTINE FinalizeOpacities_TABLE


  SUBROUTINE ComputeAbsorptionOpacity_TABLE &
               ( E, D, T, Y, X1, X2, X3, Chi )

    REAL(DP), DIMENSION(:),   INTENT(in)  :: E, D, T, Y, X1, X2, X3
    REAL(DP), DIMENSION(:,:), INTENT(out) :: Chi

    REAL(DP), DIMENSION(SIZE(E)) :: LogE

    ASSOCIATE &
      ( Chi_T => OPACITIES % ThermEmAb % Absorptivity(1) % Values, &
        OS    => OPACITIES % ThermEmAb % Offsets(1) )

    LogE = LOG10( E / MeV )

    CALL LogInterpolateSingleVariable_1D3D &
           ( LogE, LOG10( D / ( Gram / Centimeter**3 ) ), &
             LOG10( T / Kelvin ), Y, &
             LogEs_T, LogDs_T, LogTs_T, Ys_T, OS, &
             Chi_T, Chi )

    Chi(:,:) = Chi(:,:) / Centimeter

    END ASSOCIATE ! Chi_T, etc.

  END SUBROUTINE ComputeAbsorptionOpacity_TABLE


  SUBROUTINE ComputeScatteringOpacity_ES_TABLE &
               ( E, D, T, Y, X1, X2, X3, Sigma )

    REAL(DP), DIMENSION(:),   INTENT(in)  :: E, D, T, Y, X1, X2, X3
    REAL(DP), DIMENSION(:,:), INTENT(out) :: Sigma

    REAL(DP), DIMENSION(SIZE(E)) :: LogE

    ASSOCIATE &
      ( Sigma_T => OPACITIES % Scatt_Iso % Kernel(1) % Values(:,:,:,:,1), &
        OS      => OPACITIES % Scatt_Iso % Offsets(1,1) )

    LogE = LOG10( E / MeV )

    CALL LogInterpolateSingleVariable_1D3D &
           ( LogE, LOG10( D / ( Gram / Centimeter**3 ) ), &
             LOG10( T / Kelvin ), Y, &
             LogEs_T, LogDs_T, LogTs_T, Ys_T, OS, &
             Sigma_T, Sigma )

    Sigma(:,:) = Sigma(:,:) / Centimeter

    END ASSOCIATE ! Sigma_T, etc.

  END SUBROUTINE ComputeScatteringOpacity_ES_TABLE


  SUBROUTINE ComputeScatteringOpacity_NES_TABLE( E, T, Eta, R0_In, R0_Out )

    REAL(DP), DIMENSION(:),     INTENT(in)  :: E, T, Eta
    REAL(DP), DIMENSION(:,:,:), INTENT(out) :: R0_In, R0_Out

    INTEGER :: iX
    REAL(DP), DIMENSION(SIZE(E)) :: LogE

    LogE = LOG10( E / MeV )

    CALL LogInterpolateSingleVariable_2D2D &
           ( LogE, LogE, LOG10( T / Kelvin ), LOG10( Eta ), &
             LogEs_T, LogEs_T, LogTs_T, LogEtas_T, &
             OPACITIES % Scatt_NES % Offsets(1,1), &
             OPACITIES % Scatt_NES % Kernel(1) % Values(:,:,:,:,1), &
             R0_Out )

    R0_Out = R0_Out / ( Centimeter * MeV**3 )

    DO iX = 1, SIZE( T )

      R0_In(:,:,iX) = TRANSPOSE( R0_Out(:,:,iX) )

    END DO

  END SUBROUTINE ComputeScatteringOpacity_NES_TABLE


END MODULE OpacityModule_TABLE
