MODULE BoundaryConditionsModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nX, nNodesX, bcX, nE
  USE UtilitiesModule, ONLY: &
    NodeNumberX
  USE RadiationFieldsModule, ONLY: &
    uCR, nCR, nSpecies
  USE ApplicationBoundaryConditionsModule, ONLY: &
    ApplyApplicationBoundaryConditions_Radiation_X1

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ApplyBoundaryConditions_Radiation

CONTAINS


  SUBROUTINE ApplyBoundaryConditions_Radiation( Time, LimiterBC_Option )

    REAL(DP), INTENT(in) :: Time
    LOGICAL,  INTENT(in), OPTIONAL :: LimiterBC_Option

    IF( bcX(1) == 10 )THEN

      CALL ApplyApplicationBoundaryConditions_Radiation_X1 &
             ( Time, LimiterBC_Option )

    ELSE

      CALL ApplyBoundaryConditions_Radiation_X1

    END IF

  END SUBROUTINE ApplyBoundaryConditions_Radiation


  SUBROUTINE ApplyBoundaryConditions_Radiation_X1

    INTEGER :: iS, iX2, iX3, iE, iCR

    SELECT CASE ( bcX(1) )

      CASE ( 0 ) ! No Boundary Condition

      CASE ( 1 ) ! Periodic

        DO iS = 1, nSpecies
          DO iX3 = 1, nX(3)
            DO iX2 = 1, nX(2)
              DO iE = 1, nE

                DO iCR = 1, nCR

                  uCR(:,iE,0,iX2,iX3,iCR,iS) &
                    = uCR(:,iE,nX(1),iX2,iX3,iCR,iS)

                  uCR(:,iE,nX(1)+1,iX2,iX3,iCR,iS) &
                    = uCR(:,iE,1,iX2,iX3,iCR,iS)

                END DO

              END DO
            END DO
          END DO
        END DO

      CASE ( 2 ) ! Homogeneous

        DO iS = 1, nSpecies
          DO iX3 = 1, nX(3)
            DO iX2 = 1, nX(2)
              DO iE = 1, nE

                DO iCR = 1, nCR

                  uCR(:,iE,0,iX2,iX3,iCR,iS) &
                    = uCR(:,iE,1,iX2,iX3,iCR,iS)

                  uCR(:,iE,nX(1)+1,iX2,iX3,iCR,iS) &
                    = uCR(:,iE,nX(1),iX2,iX3,iCR,iS)

                END DO

              END DO
            END DO
          END DO
        END DO

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A5,A45,I2.2)') &
          '', 'Invalid Boundary Condition for Radiation X1: ', bcX(1)
        STOP

    END SELECT

  END SUBROUTINE ApplyBoundaryConditions_Radiation_X1


END MODULE BoundaryConditionsModule
