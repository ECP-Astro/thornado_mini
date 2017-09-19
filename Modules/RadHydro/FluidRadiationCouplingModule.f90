MODULE FluidRadiationCouplingModule

  USE KindModule, ONLY: &
    DP
  USE FluidRadiationCouplingSolutionModule_EmAbScatt, ONLY: &
    CoupleFluidRadiation_EmAbScatt

  IMPLICIT NONE
  PRIVATE

  CHARACTER(32) :: FluidRadiationCoupling = 'Dummy'

  PROCEDURE (CouplingProcedure), POINTER, PUBLIC :: &
    CoupleFluidRadiation => NULL()

  INTERFACE
    SUBROUTINE CouplingProcedure( dt, iX_Begin, iX_End, EvolveFluid_Option )
      USE KindModule, ONLY: DP
      REAL(DP),              INTENT(in) :: dt
      INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End
      LOGICAL,               INTENT(in), OPTIONAL :: EvolveFluid_Option
    END SUBROUTINE CouplingProcedure
  END INTERFACE

  PUBLIC :: InitializeFluidRadiationCoupling
  PUBLIC :: FinalizeFluidRadiationCoupling

CONTAINS


  SUBROUTINE InitializeFluidRadiationCoupling &
               ( FluidRadiationCoupling_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      FluidRadiationCoupling_Option

    IF( PRESENT( FluidRadiationCoupling_Option ) )THEN
      FluidRadiationCoupling = FluidRadiationCoupling_Option
    END IF

    SELECT CASE ( TRIM( FluidRadiationCoupling ) )

      CASE( 'EmAbScatt' )

        CoupleFluidRadiation &
          => CoupleFluidRadiation_EmAbScatt

      CASE DEFAULT

        CoupleFluidRadiation &
          => CoupleFluidRadiation_Dummy

    END SELECT

  END SUBROUTINE InitializeFluidRadiationCoupling


  SUBROUTINE FinalizeFluidRadiationCoupling

    NULLIFY( CoupleFluidRadiation )

  END SUBROUTINE FinalizeFluidRadiationCoupling


  SUBROUTINE CoupleFluidRadiation_Dummy &
               ( dt, iX_Begin, iX_End, EvolveFluid_Option )

    REAL(DP),              INTENT(in) :: dt
    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End
    LOGICAL,               INTENT(in), OPTIONAL :: EvolveFluid_Option

    WRITE(*,*)
    WRITE(*,'(A4,A)') &
      '', 'FluidRadiationCouplingModule: CoupleFluidRadiation_Dummy'
    WRITE(*,*)

  END SUBROUTINE CoupleFluidRadiation_Dummy


END MODULE FluidRadiationCouplingModule
