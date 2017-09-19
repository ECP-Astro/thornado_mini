MODULE EulerEquationsUtilitiesModule

  USE KindModule, ONLY: &
    DP
  USE FluidFieldsModule, ONLY: &
    nCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne, &
    nPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Conserved
  PUBLIC :: ComputeConserved
  PUBLIC :: Primitive
  PUBLIC :: ComputePrimitive

CONTAINS


  PURE FUNCTION Conserved( Primitive )

    REAL(DP), DIMENSION(1:nPF), INTENT(in) :: Primitive
    REAL(DP), DIMENSION(1:nCF)             :: Conserved

    Conserved(iCF_D)  &
      = Primitive(iPF_D)
    Conserved(iCF_S1) &
      = Primitive(iPF_D) * Primitive(iPF_V1)
    Conserved(iCF_S2) &
      = Primitive(iPF_D) * Primitive(iPF_V2)
    Conserved(iCF_S3) &
      = Primitive(iPF_D) * Primitive(iPF_V3)
    Conserved(iCF_E)  &
      = Primitive(iPF_E) &
          + 0.5_DP * Primitive(iPF_D) &
              * ( Primitive(iPF_V1)**2 + Primitive(iPF_V2)**2 &
                    + Primitive(iPF_V3)**2 )
    Conserved(iCF_Ne) &
      = Primitive(iPF_Ne)

    RETURN
  END FUNCTION Conserved


  PURE FUNCTION Primitive( Conserved )

    REAL(DP), DIMENSION(1:nCF), INTENT(in) :: Conserved
    REAL(DP), DIMENSION(1:nPF)             :: Primitive

    REAL(DP) :: CF_D

    CF_D = MAX( Conserved(iCF_D), TINY( 1.0_DP ) )

    Primitive(iPF_D)  &
      = CF_D
    Primitive(iPF_V1) &
      = Conserved(iCF_S1) / CF_D
    Primitive(iPF_V2) &
      = Conserved(iCF_S2) / CF_D
    Primitive(iPF_V3) &
      = Conserved(iCF_S3) / CF_D
    Primitive(iPF_E)  &
      = Conserved(iCF_E) &
          - 0.5_DP * ( Conserved(iCF_S1)**2 + Conserved(iCF_S2)**2 &
                       + Conserved(iCF_S3)**2 ) / CF_D
    Primitive(iPF_Ne)  &
      = Conserved(iCF_Ne)

    RETURN
  END FUNCTION Primitive


  SUBROUTINE ComputeConserved( uPF, uCF )

    REAL(DP), DIMENSION(:,:), INTENT(in)  :: uPF
    REAL(DP), DIMENSION(:,:), INTENT(out) :: uCF

    uCF(:,iCF_D) &
      = uPF(:,iPF_D)

    uCF(:,iCF_S1) &
      = uCF(:,iCF_D) * uPF(:,iPF_V1)

    uCF(:,iCF_S2) &
      = uCF(:,iCF_D) * uPF(:,iPF_V2)

    uCF(:,iCF_S3) &
      = uCF(:,iCF_D) * uPF(:,iPF_V3)

    uCF(:,iCF_E) &
      = uPF(:,iPF_E) &
        + 0.5_DP * uPF(:,iPF_D) &
            * ( uPF(:,iPF_V1)**2 + uPF(:,iPF_V2)**2 + uPF(:,iPF_V3)**2 )

    uCF(:,iCF_Ne) &
      = uPF(:,iPF_Ne)

  END SUBROUTINE ComputeConserved


  SUBROUTINE ComputePrimitive( uCF, uPF )

    REAL(DP), DIMENSION(:,:), INTENT(in)  :: uCF
    REAL(DP), DIMENSION(:,:), INTENT(out) :: uPF

    uPF(:,iPF_D) &
      = uCF(:,iCF_D)

    uPF(:,iPF_V1) &
      = uCF(:,iCF_S1) / uCF(:,iCF_D)

    uPF(:,iPF_V2) &
      = uCF(:,iCF_S2) / uCF(:,iCF_D)

    uPF(:,iPF_V3) &
      = uCF(:,iCF_S3) / uCF(:,iCF_D)

    uPF(:,iPF_E) &
      = uCF(:,iCF_E) &
        - 0.5_DP * ( uCF(:,iCF_S1)**2 + uCF(:,iCF_S2)**2 &
                     + uCF(:,iCF_S3)**2 ) / uCF(:,iCF_D)

    uPF(:,iPF_Ne) &
      = uCF(:,iCF_Ne)

  END SUBROUTINE ComputePrimitive


END MODULE EulerEquationsUtilitiesModule
