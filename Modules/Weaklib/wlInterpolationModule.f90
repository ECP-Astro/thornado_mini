MODULE wlInterpolationModule

  USE wlKindModule, ONLY: &
    dp
  USE wlThermoStateModule, ONLY: &
    ThermoStateType
  USE wlDependentVariablesModule, ONLY: &
    DependentVariablesType

  implicit none
  private

  INCLUDE 'mpif.h'

  PUBLIC :: LogInterpolateSingleVariable_3D
  PUBLIC :: LogInterpolateDifferentiateSingleVariable_3D
  PUBLIC :: LogInterpolateSingleVariable_1D3D
  PUBLIC :: LogInterpolateSingleVariable_2D2D
  PUBLIC :: ComputeTempFromIntEnergy_Lookup

  REAL(dp), PARAMETER :: One  = 1.00_dp
  REAL(dp), PARAMETER :: Ten  = 10.0_dp
  REAL(dp), PARAMETER :: ln10 = LOG(Ten)

CONTAINS


  PURE INTEGER FUNCTION Index1D( x, xx, n )

    REAL(dp), INTENT(in) :: x, xx(n)
    INTEGER,  INTENT(in) :: n

    INTEGER :: il, im, iu

    il = 0
    iu = n+1
    DO WHILE ( iu - il > 1 )
      im = (iu+il)/2
      IF ((xx(n).ge.xx(1)).eqv.(x.ge.xx(im))) THEN
        il = im
      ELSE
        iu = im
      END IF
    END DO

    IF (x.eq.xx(1)) THEN
      Index1D = 1
    ELSEIF (x.eq.xx(n)) THEN
      Index1D = n-1
    ELSE
      Index1D = il
    END IF

    RETURN
  END FUNCTION Index1D


  PURE INTEGER FUNCTION Index1D_Lin( x, xx, n )

    REAL(dp), INTENT(in) :: x, xx(n)
    INTEGER,  INTENT(in) :: n

    Index1D_Lin &
      = MAX( 1, &
             MIN( FLOOR( 1 + (n-1)*(x-xx(1))/(xx(n)-xx(1)) + 1.d-12 ), &
                  n - 1 ) )

    RETURN
  END FUNCTION Index1D_Lin


  PURE INTEGER FUNCTION Index1D_Log( x, xx, n )

    REAL(dp), INTENT(in) :: x, xx(n)
    INTEGER,  INTENT(in) :: n

    Index1D_Log &
      = FLOOR( 1 + (n-1)*LOG10(x/xx(1))/LOG10(xx(n)/xx(1)) + 1.d-12 )

    RETURN
  END FUNCTION Index1D_Log


  PURE REAL(dp) FUNCTION TriLinear &
    ( p000, p100, p010, p110, p001, p101, p011, p111, dX1, dX2, dX3 )

    REAL(dp), INTENT(in) :: &
      p000, p100, p010, p110, &
      p001, p101, p011, p111, &
      dX1, dX2, dX3

    REAL(dp) :: ddX1, ddX2, ddX3

    ddX1 = One - dX1
    ddX2 = One - dX2
    ddX3 = One - dX3

    TriLinear                                        &
      = ddX3                                         &
         * (   ddX2 * ( ddX1 * p000 + dX1 * p100 )   &
             +  dX2 * ( ddX1 * p010 + dX1 * p110 ) ) &
      +  dX3                                         &
         * (   ddX2 * ( ddX1 * p001 + dX1 * p101 )   &
             +  dX2 * ( ddX1 * p011 + dX1 * p111 ) )

    RETURN
  END FUNCTION TriLinear


  PURE REAL(dp) FUNCTION dTrilineardX1 &
    ( p000, p100, p010, p110, p001, p101, p011, p111, dX2, dX3 )

    REAL(dp), INTENT(in) :: &
      p000, p100, p010, p110, &
      p001, p101, p011, p111, &
      dX2, dX3

    REAL(dp) :: ddX2, ddX3

    ddX2 = One - dX2
    ddX3 = One - dX3

    dTrilineardX1 &
      = ddX3  * ( ddX2 * ( p100 - p000 ) + dX2 * ( p110 - p010 ) ) &
        + dX3 * ( ddX2 * ( p101 - p001 ) + dX2 * ( p111 - p011 ) )

    RETURN
  END FUNCTION dTrilineardX1


  PURE REAL(dp) FUNCTION dTrilineardX2 &
    ( p000, p100, p010, p110, p001, p101, p011, p111, dX1, dX3 )

    REAL(dp), INTENT(in) :: &
      p000, p100, p010, p110, &
      p001, p101, p011, p111, &
      dX1, dX3

    REAL(dp) :: ddX1, ddX3

    ddX1 = One - dX1
    ddX3 = One - dX3

    dTrilineardX2 &
      = ddX3  * ( ddX1 * ( p010 - p000 ) + dX1 * ( p110 - p100 ) ) &
        + dX3 * ( ddX1 * ( p011 - p001 ) + dX1 * ( p111 - p101 ) )

    RETURN
  END FUNCTION dTrilineardX2


  PURE REAL(dp) FUNCTION dTrilineardX3 &
    ( p000, p100, p010, p110, p001, p101, p011, p111, dX1, dX2 )

    REAL(dp), INTENT(in) :: &
      p000, p100, p010, p110, &
      p001, p101, p011, p111, &
      dX1, dX2

    REAL(dp) :: &
      ddX1, ddX2

    ddX1 = One - dX1
    ddX2 = One - dX2

    dTrilineardX3 &
      = ddX1  * ( ddX2 * ( p001 - p000 ) + dX2 * ( p011 - p010 ) ) &
        + dX1 * ( ddX2 * ( p101 - p100 ) + dX2 * ( p111 - p110 ) )

    RETURN
  END FUNCTION dTrilineardX3


  PURE REAL(dp) FUNCTION TetraLinear &
    ( p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
      p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
      dX1, dX2, dX3, dX4 )

    REAL(dp), INTENT(in) :: &
      p0000, p1000, p0100, p1100, p0010, p1010, p0110, p1110, &
      p0001, p1001, p0101, p1101, p0011, p1011, p0111, p1111, &
      dX1, dX2, dX3, dX4

    REAL(dp) :: ddX1, ddX2, ddX3, ddX4

    ddX1 = One - dX1
    ddX2 = One - dX2
    ddX3 = One - dX3
    ddX4 = One - dX4

    TetraLinear                                                    &
      = ddX4                                                       &
         * (   ddX3 * (  ddX2 * ( ddX1 * p0000 + dX1 * p1000 )     &
                        + dX2 * ( ddX1 * p0100 + dX1 * p1100 ) )   &
             +  dX3 * (  ddX2 * ( ddX1 * p0010 + dX1 * p1010 )     &
                        + dX2 * ( ddX1 * p0110 + dX1 * p1110 ) ) ) &
      +  dX4                                                       &
         * (   ddX3 * (  ddX2 * ( ddX1 * p0001 + dX1 * p1001 )     &
                        + dX2 * ( ddX1 * p0101 + dX1 * p1101 ) )   &
             +  dX3 * (  ddX2 * ( ddX1 * p0011 + dX1 * p1011 )     &
                        + dX2 * ( ddX1 * p0111 + dX1 * p1111 ) ) )

    RETURN
  END FUNCTION TetraLinear


  SUBROUTINE LogInterpolateSingleVariable_3D &
               ( LogD, LogT, LinY, LogDs, LogTs, LinYs, OS, Table, &
                 Interpolant, ReturnLog_Option, Mask_Option )

    REAL(dp), DIMENSION(:),     INTENT(in)  :: LogD
    REAL(dp), DIMENSION(:),     INTENT(in)  :: LogT
    REAL(dp), DIMENSION(:),     INTENT(in)  :: LinY
    REAL(dp), DIMENSION(:),     INTENT(in)  :: LogDs
    REAL(dp), DIMENSION(:),     INTENT(in)  :: LogTs
    REAL(dp), DIMENSION(:),     INTENT(in)  :: LinYs
    REAL(dp),                   INTENT(in)  :: OS
    REAL(dp), DIMENSION(:,:,:), INTENT(in)  :: Table
    REAL(dp), DIMENSION(:),     INTENT(out) :: Interpolant
    LOGICAL,                    INTENT(in), OPTIONAL :: ReturnLog_Option
    LOGICAL,  DIMENSION(:),     INTENT(in), OPTIONAL :: Mask_Option

    LOGICAL, DIMENSION(SIZE(LogD)) :: &
      Mask
    LOGICAL :: &
      ReturnLog
    INTEGER :: &
      iP, &
      iD, iT, iY, &
      pD, pT, pY
    REAL(dp) :: &
      dD, dT, dY
    REAL(dp), DIMENSION(0:1,0:1,0:1) :: &
      p

    Mask = .TRUE.
    IF( PRESENT( Mask_Option ) ) &
      Mask = Mask_Option

    ReturnLog = .FALSE.
    IF( PRESENT( ReturnLog_Option ) ) &
      ReturnLog = ReturnLog_Option

    DO iP = 1, SIZE( LogD )

      IF( .NOT. Mask(iP) ) CYCLE

      iD = Index1D_Lin( LogD(iP), LogDs, SIZE( LogDs ) )
      iT = Index1D_Lin( LogT(iP), LogTs, SIZE( LogTs ) )
      iY = Index1D_Lin( LinY(iP), LinYs, SIZE( LinYs ) )

      dD = ( LogD(iP) - LogDs(iD) ) / ( LogDs(iD+1) - LogDs(iD) )
      dT = ( LogT(iP) - LogTs(iT) ) / ( LogTs(iT+1) - LogTs(iT) )
      dY = ( LinY(iP) - LinYs(iY) ) / ( LinYs(iY+1) - LinYs(iY) )

      DO pY = 0, 1
        DO pT = 0, 1
          DO pD = 0, 1

            p(pD,pT,pY) &
              = TABLE(iD+pD,iT+pT,iY+pY)

          END DO
        END DO
      END DO

      Interpolant(iP) &
        = TriLinear &
            ( p(0,0,0), p(1,0,0), p(0,1,0), p(1,1,0), &
              p(0,0,1), p(1,0,1), p(0,1,1), p(1,1,1), &
              dD, dT, dY )

    END DO

    IF( .NOT. ReturnLog ) &
      Interpolant = Ten**( Interpolant(:) ) - OS

  END SUBROUTINE LogInterpolateSingleVariable_3D


  SUBROUTINE LogInterpolateDifferentiateSingleVariable_3D &
               ( LogD, LogT, LinY, LogDs, LogTs, LinYs, OS, Table, &
                 Interpolant, Derivative, Mask_Option )

    REAL(dp), DIMENSION(:),     INTENT(in)  :: LogD
    REAL(dp), DIMENSION(:),     INTENT(in)  :: LogT
    REAL(dp), DIMENSION(:),     INTENT(in)  :: LinY
    REAL(dp), DIMENSION(:),     INTENT(in)  :: LogDs
    REAL(dp), DIMENSION(:),     INTENT(in)  :: LogTs
    REAL(dp), DIMENSION(:),     INTENT(in)  :: LinYs
    REAL(dp),                   INTENT(in)  :: OS
    REAL(dp), DIMENSION(:,:,:), INTENT(in)  :: Table
    REAL(dp), DIMENSION(:),     INTENT(out) :: Interpolant
    REAL(dp), DIMENSION(:,:),   INTENT(out) :: Derivative
    LOGICAL,  DIMENSION(:),     INTENT(in), OPTIONAL :: Mask_Option

    LOGICAL, DIMENSION(SIZE(LogD)) :: &
      Mask
    INTEGER :: &
      iP, &
      iD, iT, iY, &
      pD, pT, pY
    REAL(dp) :: &
      dD, dT, dY, &
      aD, aT, aY
    REAL(dp), DIMENSION(0:1,0:1,0:1) :: &
      p

    Mask = .TRUE.
    IF( PRESENT( Mask_Option ) ) &
      Mask = Mask_Option

    DO iP = 1, SIZE( LogD )

      IF( .NOT. Mask(iP) ) CYCLE

      iD = Index1D_Lin( LogD(iP), LogDs, SIZE( LogDs ) )
      iT = Index1D_Lin( LogT(iP), LogTs, SIZE( LogTs ) )
      iY = Index1D_Lin( LinY(iP), LinYs, SIZE( LinYs ) )

      dD = ( LogD(iP) - LogDs(iD) ) / ( LogDs(iD+1) - LogDs(iD) )
      dT = ( LogT(iP) - LogTs(iT) ) / ( LogTs(iT+1) - LogTs(iT) )
      dY = ( LinY(iP) - LinYs(iY) ) / ( LinYs(iY+1) - LinYs(iY) )

      aD = One  / ( Ten**LogD(iP) * ( LogDs(iD+1) - LogDs(iD) ) )
      aT = One  / ( Ten**LogT(iP) * ( LogTs(iT+1) - LogTs(iT) ) )
      aY = ln10 /                   ( LinYs(iY+1) - LinYs(iY) )

      DO pY = 0, 1
        DO pT = 0, 1
          DO pD = 0, 1

            p(pD,pT,pY) &
              = TABLE(iD+pD,iT+pT,iY+pY)

          END DO
        END DO
      END DO

      Interpolant(iP) &
        = Ten**( TriLinear &
                   ( p(0,0,0), p(1,0,0), p(0,1,0), p(1,1,0), &
                     p(0,0,1), p(1,0,1), p(0,1,1), p(1,1,1), &
                     dD, dT, dY ) ) - OS

      Derivative(iP,1) &
        = Interpolant(iP) * aD &
            * dTriLineardX1 &
                ( p(0,0,0), p(1,0,0), p(0,1,0), p(1,1,0), &
                  p(0,0,1), p(1,0,1), p(0,1,1), p(1,1,1), &
                  dT, dY )

      Derivative(iP,2) &
        = Interpolant(iP) * aT &
            * dTriLineardX2 &
                ( p(0,0,0), p(1,0,0), p(0,1,0), p(1,1,0), &
                  p(0,0,1), p(1,0,1), p(0,1,1), p(1,1,1), &
                  dD, dY )

      Derivative(iP,3) &
        = Interpolant(iP) * aY &
            * dTriLineardX3 &
                ( p(0,0,0), p(1,0,0), p(0,1,0), p(1,1,0), &
                  p(0,0,1), p(1,0,1), p(0,1,1), p(1,1,1), &
                  dD, dT )

    END DO

  END SUBROUTINE LogInterpolateDifferentiateSingleVariable_3D

  
  SUBROUTINE LogInterpolateSingleVariable_1D3D &
               ( LogX1, LogX2, LogX3, LinX4, LogCoordsX1, LogCoordsX2, &
                 LogCoordsX3, LinCoordsX4, Offset, Table, Interpolant )

    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogX1
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogX2
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogX3
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LinX4
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogCoordsX1
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogCoordsX2
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogCoordsX3
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LinCoordsX4
    REAL(dp),                     INTENT(in)  :: Offset
    REAL(dp), DIMENSION(:,:,:,:), INTENT(in)  :: Table
    REAL(dp), DIMENSION(:,:),     INTENT(out) :: Interpolant

    INTEGER :: i, j
    INTEGER :: iX1, iX2, iX3, iX4
    INTEGER :: p1, p2, p3, p4
    REAL(dp), DIMENSION(4) :: dX
    REAL(dp), DIMENSION(0:1,0:1,0:1,0:1) :: p

    !$OMP PARALLEL DO PRIVATE &
    !$OMP&              ( i, j, iX1, iX2, iX3, iX4, dX, p1, p2, p3, p4, p )
    DO j = 1, SIZE( LogX2 )

      iX4   = Index1D_Lin( LinX4(j), LinCoordsX4, SIZE( LinCoordsX4 ) )
      dX(4) = ( LinX4(j) - LinCoordsX4(iX4) ) &
              / ( LinCoordsX4(iX4+1) - LinCoordsX4(iX4) )

      iX3   = Index1D_Lin( LogX3(j), LogCoordsX3, SIZE( LogCoordsX3 ) )
      dX(3) = ( LogX3(j) - LogCoordsX3(iX3) ) &
              / ( LogCoordsX3(iX3+1) - LogCoordsX3(iX3) )

      iX2   = Index1D_Lin( LogX2(j), LogCoordsX2, SIZE( LogCoordsX2 ) )
      dX(2) = ( LogX2(j) - LogCoordsX2(iX2) ) &
              / ( LogCoordsX2(iX2+1) - LogCoordsX2(iX2) )

      DO i = 1, SIZE( LogX1 )

        iX1   = Index1D_Lin( LogX1(i), LogCoordsX1, SIZE( LogCoordsX1 ) )
        dX(1) = ( LogX1(i) - LogCoordsX1(iX1) ) &
                / ( LogCoordsX1(iX1+1) - LogCoordsX1(iX1) )

        DO p4 = 0, 1
          DO p3 = 0, 1
            DO p2 = 0, 1
              DO p1 = 0, 1

                p(p1,p2,p3,p4) &
                  = TABLE(iX1+p1,iX2+p2,iX3+p3,iX4+p4)

              END DO
            END DO
          END DO
        END DO

        Interpolant(i,j) &
          = TetraLinear &
              ( p(0,0,0,0), p(1,0,0,0), p(0,1,0,0), p(1,1,0,0), &
                p(0,0,1,0), p(1,0,1,0), p(0,1,1,0), p(1,1,1,0), &
                p(0,0,0,1), p(1,0,0,1), p(0,1,0,1), p(1,1,0,1), &
                p(0,0,1,1), p(1,0,1,1), p(0,1,1,1), p(1,1,1,1), &
                dX(1), dX(2), dX(3), dX(4) )

        Interpolant(i,j) &
          = Ten**( Interpolant(i,j) ) - Offset

      END DO ! i
    END DO ! j
    !$OMP END PARALLEL DO

  END SUBROUTINE LogInterpolateSingleVariable_1D3D


  SUBROUTINE LogInterpolateSingleVariable_2D2D &
               ( LogX1, LogX2, LogX3, LogX4, LogCoordsX1, LogCoordsX2, &
                 LogCoordsX3, LogCoordsX4, Offset, Table, Interpolant )

    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogX1
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogX2
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogX3
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogX4
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogCoordsX1
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogCoordsX2
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogCoordsX3
    REAL(dp), DIMENSION(:),       INTENT(in)  :: LogCoordsX4
    REAL(dp),                     INTENT(in)  :: Offset
    REAL(dp), DIMENSION(:,:,:,:), INTENT(in)  :: Table
    REAL(dp), DIMENSION(:,:,:),   INTENT(out) :: Interpolant

    INTEGER :: i, j, k
    INTEGER :: iX1, iX2, iX3, iX4
    INTEGER :: p1, p2, p3, p4
    REAL(dp), DIMENSION(4) :: dX
    REAL(dp), DIMENSION(0:1,0:1,0:1,0:1) :: p

    DO k = 1, SIZE( LogX3 )

      iX4 &
        = Index1D_Lin( LogX4(k), LogCoordsX4, SIZE( LogCoordsX4 ) )
      dX(4) &
        = ( LogX4(k) - LogCoordsX4(iX4) ) &
            / ( LogCoordsX4(iX4+1) - LogCoordsX4(iX4) )

      iX3 &
        = Index1D_Lin( LogX3(k), LogCoordsX3, SIZE( LogCoordsX3 ) )
      dX(3) &
        = ( LogX3(k) - LogCoordsX3(iX3) ) &
            / ( LogCoordsX3(iX3+1) - LogCoordsX3(iX3) )

      DO j = 1, SIZE( LogX2 )

        iX2 &
          = Index1D_Lin( LogX2(j), LogCoordsX2, SIZE( LogCoordsX2 ) )
        dX(2) &
          = ( LogX2(j) - LogCoordsX2(iX2) ) &
              / ( LogCoordsX2(iX2+1) - LogCoordsX2(iX2) )

        DO i = 1, SIZE( LogX1 )

          iX1 &
            = Index1D_Lin( LogX1(i), LogCoordsX1, SIZE( LogCoordsX1 ) )
          dX(1) &
            = ( LogX1(i) - LogCoordsX1(iX1) ) &
                / ( LogCoordsX1(iX1+1) - LogCoordsX1(iX1) )

          DO p4 = 0, 1
            DO p3 = 0, 1
              DO p2 = 0, 1
                DO p1 = 0, 1

                  p(p1,p2,p3,p4) &
                    = TABLE(iX1+p1,iX2+p2,iX3+p3,iX4+p4)

                END DO
              END DO
            END DO
          END DO

          Interpolant(i,j,k) &
            = TetraLinear &
                ( p(0,0,0,0), p(1,0,0,0), p(0,1,0,0), p(1,1,0,0), &
                  p(0,0,1,0), p(1,0,1,0), p(0,1,1,0), p(1,1,1,0), &
                  p(0,0,0,1), p(1,0,0,1), p(0,1,0,1), p(1,1,0,1), &
                  p(0,0,1,1), p(1,0,1,1), p(0,1,1,1), p(1,1,1,1), &
                  dX(1), dX(2), dX(3), dX(4) )

        END DO ! i
      END DO ! j
    END DO ! k

    Interpolant(:,:,:) &
      = Ten**( Interpolant(:,:,:) ) - Offset

  END SUBROUTINE LogInterpolateSingleVariable_2D2D


  SUBROUTINE ComputeTempFromIntEnergy_Lookup &
               ( LogD, LogE, LinY, LogDs, LogTs, LinYs, LogEs, OS, T, &
                 Mask_Option )

    REAL(dp), DIMENSION(:),     INTENT(in)  :: LogD
    REAL(dp), DIMENSION(:),     INTENT(in)  :: LogE
    REAL(dp), DIMENSION(:),     INTENT(in)  :: LinY
    REAL(dp), DIMENSION(:),     INTENT(in)  :: LogDs
    REAL(dp), DIMENSION(:),     INTENT(in)  :: LogTs
    REAL(dp), DIMENSION(:),     INTENT(in)  :: LinYs
    REAL(dp), DIMENSION(:,:,:), INTENT(in)  :: LogEs
    REAL(dp),                   INTENT(in)  :: OS
    REAL(dp), DIMENSION(:),     INTENT(out) :: T
    LOGICAL,  DIMENSION(:),     INTENT(in), OPTIONAL :: Mask_Option

    LOGICAL, DIMENSION(SIZE(LogD)) :: Mask
    INTEGER  :: iP
    INTEGER  :: sizeLogDs, sizeLogTs, sizeLinYs
    INTEGER  :: iD, iY, iE, iE1, iE2, iE3, iE4, iEa, iEb
    INTEGER  :: nPtsE
    REAL(dp) :: tmpE
    REAL(dp), DIMENSION(:), ALLOCATABLE :: ptsLogD, ptsLogT, ptsLinY, ptsLogE

    Mask = .TRUE.
    IF( PRESENT( Mask_Option ) ) &
      Mask = Mask_Option

    sizeLogDs = SIZE( LogDs )
    sizeLogTs = SIZE( LogTs )
    sizeLinYs = SIZE( LinYs )

    ALLOCATE &
      ( ptsLogD(sizeLogTs), ptsLogT(sizeLogTs), &
        ptsLinY(sizeLogTs), ptsLogE(sizeLogTs) )

    DO iP = 1, SIZE( LogD )

      IF( .NOT. Mask(iP) ) CYCLE

      iD = Index1D_Lin( LogD(iP), LogDs, sizeLogDs )
      iY = Index1D_Lin( LinY(iP), LinYs, sizeLinYs )

      iE1 = Index1D( logE(iP), LogEs(iD,  :,iY  ), sizeLogTs )
      iE2 = Index1D( logE(iP), LogEs(iD+1,:,iY  ), sizeLogTs )
      iE3 = Index1D( logE(iP), LogEs(iD,  :,iY+1), sizeLogTs )
      iE4 = Index1D( logE(iP), LogEs(iD+1,:,iY+1), sizeLogTs )

      iEa = MAX( MINVAL( [iE1,iE2,iE3,iE4] ) - 1, 1 )
      iEb = MIN( MAXVAL( [iE1,iE2,iE3,iE4] ) + 2, sizeLogTs )

      nPtsE = SIZE( LogTs(iEa:iEb) )

      ptsLogD(1:nPtsE) = LogD (iP)
      ptsLogT(1:nPtsE) = LogTs(iEa:iEb)
      ptsLinY(1:nPtsE) = LinY (iP)

      CALL LogInterpolateSingleVariable_3D &
             ( ptsLogD(1:nPtsE), ptsLogT(1:nPtsE), ptsLinY(1:nPtsE), &
               LogDs(iD:iD+1), LogTs(iEa:iEb), LinYs(iY:iY+1), OS,   &
               LogEs(iD:iD+1,iEa:iEb,iY:iY+1), ptsLogE(1:nPtsE), &
               ReturnLog_Option = .TRUE. )

      IF( (ptsLogE(1)-LogE(iP)) * (ptsLogE(nPtsE)-LogE(iP)) > 0.0_DP )THEN

        WRITE(*,*)
        WRITE(*,'(A4,A)') &
          '', 'Warning: ComputeTempFromIntEnergy_Lookup'
        WRITE(*,'(A6,A20,ES10.4E2,A9,ES10.4E2)') &
          '', 'No Root Between T = ', Ten**ptsLogT(1), &
          ' and T = ', Ten**ptsLogT(nPtsE)
        WRITE(*,*)
        WRITE(*,*) '  nPtsE = ', nPtsE
        WRITE(*,*) '   ptsE = ', Ten**ptsLogE(1:nPtsE) - OS
        WRITE(*,*) '    ia  = ', iEa
        WRITE(*,*) '    ib  = ', iEb
        WRITE(*,*) '    Ta  = ', Ten**LogTs(iEa)
        WRITE(*,*) '    Tb  = ', Ten**LogTs(iEb)
        WRITE(*,*) '    Ea  = ', Ten**ptsLogE(1) - OS
        WRITE(*,*) '    Eb  = ', Ten**ptsLogE(nPtsE) - OS
        WRITE(*,*) '    iP  = ', iP
        WRITE(*,*) '     E  = ', Ten**LogE(iP) - OS
        WRITE(*,*) '     D  = ', Ten**LogD(iP)
        WRITE(*,*) '     Y  = ', LinY(iP)
        WRITE(*,*)
        STOP

      END IF

      iE = Index1D( LogE(iP), ptsLogE, nPtsE )

      T(iP) = Ten**( ptsLogT(iE) &
                     + ( ptsLogT(iE+1)-ptsLogT(iE) ) &
                       * ( LogE(iP)-ptsLogE(iE) ) &
                       / ( ptsLogE(iE+1)-ptsLogE(iE) ) )

    END DO

    DEALLOCATE( ptsLogD, ptsLogT, ptsLinY, ptsLogE )

  END SUBROUTINE ComputeTempFromIntEnergy_Lookup


END MODULE wlInterpolationModule
