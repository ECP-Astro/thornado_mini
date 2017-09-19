MODULE FluidEvolutionModule

  USE KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  CHARACTER(32) :: FluidSolver

  PROCEDURE (ComputeRHS), POINTER, PUBLIC :: &
    ComputeRHS_Fluid => NULL()

  INTERFACE
    SUBROUTINE ComputeRHS( iX_Begin, iX_End )
      INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End
    END SUBROUTINE ComputeRHS
  END INTERFACE

  PROCEDURE (ApplyLimiter), POINTER, PUBLIC :: &
    ApplySlopeLimiter_Fluid      => NULL(), &
    ApplyPositivityLimiter_Fluid => NULL()

  INTERFACE
    SUBROUTINE ApplyLimiter
    END SUBROUTINE ApplyLimiter
  END INTERFACE

  PUBLIC :: InitializeFluidEvolution
  PUBLIC :: FinalizeFluidEvolution

CONTAINS


  SUBROUTINE InitializeFluidEvolution &
               ( FluidSolver_Option, &
                 ApplySlopeLimiter_Option, BetaTVB_Option, BetaTVD_Option, &
                 ApplyPositivityLimiter_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: FluidSolver_Option
    LOGICAL,          INTENT(in), OPTIONAL :: ApplySlopeLimiter_Option
    REAL(DP),         INTENT(in), OPTIONAL :: BetaTVB_Option
    REAL(DP),         INTENT(in), OPTIONAL :: BetaTVD_Option
    LOGICAL,          INTENT(in), OPTIONAL :: ApplyPositivityLimiter_Option

    FluidSolver = 'Dummy'
    IF( PRESENT( FluidSolver_Option ) ) &
      FluidSolver = FluidSolver_Option

    SELECT CASE ( TRIM( FluidSolver ) )

      CASE DEFAULT

        ComputeRHS_Fluid &
          => ComputeRHS_Dummy
        ApplySlopeLimiter_Fluid &
          => ApplyLimiter_Dummy
        ApplyPositivityLimiter_Fluid &
          => ApplyLimiter_Dummy

    END SELECT

  END SUBROUTINE InitializeFluidEvolution


  SUBROUTINE FinalizeFluidEvolution

    NULLIFY( ComputeRHS_Fluid )
    NULLIFY( ApplySlopeLimiter_Fluid )
    NULLIFY( ApplyPositivityLimiter_Fluid )

  END SUBROUTINE FinalizeFluidEvolution


  SUBROUTINE ComputeRHS_Dummy( iX_Begin, iX_End )

    INTEGER, DIMENSION(3), INTENT(in) :: iX_Begin, iX_End

    WRITE(*,*)
    WRITE(*,'(A4,A)') &
      '', 'FluidEvolutionModule: ComputeRHS_Dummy'
    WRITE(*,*)

    RETURN

  END SUBROUTINE ComputeRHS_Dummy


  SUBROUTINE ApplyLimiter_Dummy

    WRITE(*,*)
    WRITE(*,'(A4,A)') &
      '', 'FluidEvolutionModule: ApplyLimiter_Dummy'
    WRITE(*,*)

    RETURN

  END SUBROUTINE ApplyLimiter_Dummy


END MODULE FluidEvolutionModule
