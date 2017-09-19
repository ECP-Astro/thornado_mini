MODULE RiemannSolverModule

  USE KindModule, ONLY: &
    DP

  IMPLICIT NONE
  PRIVATE

  CHARACTER(32) :: FluidRiemannSolver
  CHARACTER(32) :: RadiationRiemannSolver

  PROCEDURE (RiemannSolver), POINTER, PUBLIC :: &
    NumericalFlux_Fluid     => NULL(), &
    NumericalFlux_Radiation => NULL()

  INTERFACE
    PURE FUNCTION RiemannSolver &
      ( u_L, u_R, Flux_L, Flux_R, a, aP, aM, nF )
      USE KindModule, ONLY: DP
      INTEGER,                   INTENT(in) :: nF
      REAL(DP), DIMENSION(1:nF), INTENT(in) :: u_L, u_R, flux_L, flux_R
      REAL(DP),                  INTENT(in) :: a, aP, aM
      REAL(DP), DIMENSION(1:nF)             :: RiemannSolver
    END FUNCTION RiemannSolver
  END INTERFACE

  PUBLIC :: InitializeRiemannSolvers
  PUBLIC :: FinalizeRiemannSolvers

CONTAINS


  SUBROUTINE InitializeRiemannSolvers &
               ( FluidRiemannSolver_Option, RadiationRiemannSolver_Option )

    CHARACTER(LEN=*), INTENT(in), OPTIONAL :: &
      FluidRiemannSolver_Option, &
      RadiationRiemannSolver_Option

    FluidRiemannSolver = 'LLF'
    IF( PRESENT( FluidRiemannSolver_Option ) ) &
      FluidRiemannSolver &
        = FluidRiemannSolver_Option

    WRITE(*,*)
    WRITE(*,'(A5,A22,A)') &
      '', 'Fluid Riemann Solver: ', TRIM( FluidRiemannSolver )
    WRITE(*,'(A5,A21)') &
      '', '---------------------'

    SELECT CASE ( TRIM( FluidRiemannSolver ) )
      CASE ( 'LLF' )
        NumericalFlux_Fluid => NumericalFlux_LLF
      CASE ( 'HLL' )
        NumericalFlux_Fluid => NumericalFlux_HLL
      CASE DEFAULT
        NumericalFlux_Fluid => NumericalFlux_LLF
    END SELECT

    RadiationRiemannSolver = 'LLF'
    IF( PRESENT( RadiationRiemannSolver_Option ) ) &
      RadiationRiemannSolver &
        = RadiationRiemannSolver_Option

    WRITE(*,*)
    WRITE(*,'(A5,A26,A)') &
      '', 'Radiation Riemann Solver: ', TRIM( RadiationRiemannSolver )
    WRITE(*,'(A5,A25)') &
      '', '-------------------------'

    SELECT CASE( TRIM( RadiationRiemannSolver ) )
      CASE ( 'LLF' )
        NumericalFlux_Radiation => NumericalFlux_LLF
      CASE ( 'HLL' )
        NumericalFlux_Radiation => NumericalFlux_HLL
      CASE DEFAULT
        NumericalFlux_Radiation => NumericalFlux_LLF
    END SELECT

  END SUBROUTINE InitializeRiemannSolvers


  SUBROUTINE FinalizeRiemannSolvers

    NULLIFY( NumericalFlux_Fluid )
    NULLIFY( NumericalFlux_Radiation )

  END SUBROUTINE FinalizeRiemannSolvers


  PURE FUNCTION NumericalFlux_LLF &
    ( u_L, u_R, Flux_L, Flux_R, alpha, alpha_P, alpha_M, nF )

    ! --- Local Lax-Friedrichs Flux ---

    INTEGER,                   INTENT(in) :: nF
    REAL(DP), DIMENSION(1:nF), INTENT(in) :: u_L, u_R, flux_L, flux_R
    REAL(DP),                  INTENT(in) :: alpha, alpha_P, alpha_M
    REAL(DP), DIMENSION(1:nF)             :: NumericalFlux_LLF

    NumericalFlux_LLF &
      = 0.5_DP * ( flux_L + flux_R - alpha * ( u_R - u_L ) )

    RETURN
  END FUNCTION NumericalFlux_LLF


  PURE FUNCTION NumericalFlux_HLL &
    ( u_L, u_R, Flux_L, Flux_R, alpha, alpha_P, alpha_M, nF )

    ! --- Harten-Lax-van Leer Flux ---

    INTEGER,                   INTENT(in) :: nF
    REAL(DP), DIMENSION(1:nF), INTENT(in) :: u_L, u_R, flux_L, flux_R
    REAL(DP),                  INTENT(in) :: alpha, alpha_P, alpha_M
    REAL(DP), DIMENSION(1:nF)             :: NumericalFlux_HLL

    NumericalFlux_HLL &
      = ( alpha_P * flux_L + alpha_M * flux_R &
            - alpha_P * alpha_M * ( u_R - u_L ) ) / ( alpha_P + alpha_M )

    RETURN
  END FUNCTION NumericalFlux_HLL


END MODULE RiemannSolverModule

