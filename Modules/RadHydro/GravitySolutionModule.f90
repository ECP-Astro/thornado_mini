MODULE GravitySolutionModule

  USE KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  CHARACTER(32) :: GravitySolver

  PROCEDURE (GS), POINTER, PUBLIC :: &
    SolveGravity => NULL()

  INTERFACE
    SUBROUTINE GS
    END SUBROUTINE GS
  END INTERFACE

  PUBLIC :: InitializeGravitySolver
  PUBLIC :: FinalizeGravitySolver

CONTAINS


  SUBROUTINE InitializeGravitySolver( GravitySolver_Option, PointMass_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: GravitySolver_Option
    REAL(DP),         INTENT(in), OPTIONAL :: PointMass_Option

    GravitySolver = 'Dummy'
    IF( PRESENT( GravitySolver_Option ) ) &
      GravitySolver = GravitySolver_Option

    WRITE(*,*)
    WRITE(*,'(A5,A16,A)') &
      '', 'Gravity Solver: ', TRIM( GravitySolver )
    WRITE(*,'(A5,A16)') &
      '', '--------------- '

    SELECT CASE ( TRIM( GravitySolver ) )

      CASE DEFAULT

        SolveGravity &
          => SolveGravity_Dummy

    END SELECT

  END SUBROUTINE InitializeGravitySolver


  SUBROUTINE FinalizeGravitySolver

    NULLIFY( SolveGravity )

  END SUBROUTINE FinalizeGravitySolver


  SUBROUTINE SolveGravity_Dummy

    WRITE(*,*)
    WRITE(*,'(A4,A)') &
      '', 'GravitySolutionModule: SolveGravity_Dummy'
    WRITE(*,*)

  END SUBROUTINE SolveGravity_Dummy


END MODULE GravitySolutionModule
