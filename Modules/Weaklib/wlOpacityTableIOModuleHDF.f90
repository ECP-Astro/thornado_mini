MODULE wlOpacityTableIOModuleHDF
!-----------------------------------------------------------------------
!
!    File:         wlOpacityTableIOModuleHDF.f90
!    Module:       wlOpacityTableIOModuleHDF
!    Type:         Module w/ Subroutines
!    Author:       R. Chu, Dept. Phys. & Astronomy
!                  U. Tennesee, Knoxville
!
!    Created:      3/22/16
!    WeakLib ver:  
!
!    Purpose:
!      Subroutines needed for reading, printing OpacityTable
!
!    CONTAINS:
!       DescribeOpacityTable
!
!    Modules used:
!       wlOpacityTableModule, ONLY: OpacityTableType, OpacityTypeA
!       wlEOSIOModuleHDF, ONLY: DescribeEquationOfStateTable
!       wlGridModule, ONLY: EnergyGridType
!       wlKindModule, ONLY:dp
!       HDF5
!       wlIOModuleHDF
!       wlEquationOfStateTableModule
!
!-----------------------------------------------------------------------
!  NOTE: Only Type A interaction applied. Type B and Type C interaction 
!        needs to be added for future use.
!-----------------------------------------------------------------------

  USE wlKindModule, ONLY:         &
    dp
  USE wlGridModule, ONLY:         &
    GridType
  USE wlOpacityTableModule, ONLY: &
    OpacityTableType,             &
    AllocateOpacityTable
  USE wlOpacityFieldsModule, ONLY:&
    OpacityTypeA,                 &
    OpacityTypeB,                 &
    OpacityTypeC
  USE wlIOModuleHDF, ONLY:        &
    ReadHDF,                      &
    WriteHDF,                     &
    OpenFileHDF,                  &
    CloseFileHDF,                 &
    OpenGroupHDF,                 &
    CloseGroupHDF,                &
    WriteThermoStateHDF,          &
    ReadThermoStateHDF
  USE wlEquationOfStateTableModule
  USE HDF5

  IMPLICIT NONE
  PRIVATE

  INTEGER :: hdferr

  PUBLIC WriteOpacityTableHDF
  PUBLIC ReadOpacityTableHDF

CONTAINS

  SUBROUTINE WriteOpacityTableHDF( OpacityTable, FileName )
 
    TYPE(OpacityTableType), INTENT(inout)       :: OpacityTable
    CHARACTER(len=*), INTENT(in)                :: FileName

    INTEGER(HID_T)                              :: file_id
    INTEGER(HID_T)                              :: group_id

    CHARACTER(LEN=32), DIMENSION(1)             :: tempString
    INTEGER, DIMENSION(1)                       :: tempInteger
    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d
   
    CALL OpenFileHDF( FileName, .true., file_id )

    datasize1d(1) = 1

    tempInteger(1) = OpacityTable % nOpacitiesA
    CALL WriteHDF&
         ( "nOpacitiesA", tempInteger, file_id, datasize1d )

    tempInteger(1) = OpacityTable % nOpacitiesB 
    CALL WriteHDF&
         ( "nOpacitiesB", tempInteger, file_id, datasize1d )
  
    tempInteger(1) = OpacityTable % nMomentsB     
    CALL WriteHDF&
         ( "nMomentsB", tempInteger, file_id, datasize1d )

    tempInteger(1) = OpacityTable % nOpacitiesB_NES
    CALL WriteHDF&
         ( "nOpacitiesB_NES", tempInteger, file_id, datasize1d )

    tempInteger(1) = OpacityTable % nMomentsB_NES
    CALL WriteHDF&
         ( "nMomentsB_NES", tempInteger, file_id, datasize1d )

    tempInteger(1) = OpacityTable % nOpacitiesC     
    CALL WriteHDF&
         ( "nOpacitiesC", tempInteger, file_id, datasize1d )

    tempInteger(1) = OpacityTable % nMomentsC   
    CALL WriteHDF&
         ( "nMomentsC", tempInteger, file_id, datasize1d )

    tempInteger(1) = OpacityTable % nPointsE  
    CALL WriteHDF&
         ( "nPointsE", tempInteger, file_id, datasize1d )

    tempInteger(1) = OpacityTable % nPointsEta
    CALL WriteHDF&
         ( "nPointsEta", tempInteger, file_id, datasize1d )

    datasize1d = 3
    CALL WriteHDF&
         ( "nPointsTS", OpacityTable % nPointsTS, file_id, datasize1d )

    CALL OpenGroupHDF( "EnergyGrid", .true., file_id, group_id )
    CALL WriteGridHDF( OpacityTable % EnergyGrid, group_id )
    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( "EtaGrid", .true., file_id, group_id )
    CALL WriteGridHDF( OpacityTable % EtaGrid, group_id )
    CALL CloseGroupHDF( group_id )
  
    CALL OpenGroupHDF( "ThermoState", .true., file_id, group_id )
    CALL WriteThermoStateHDF( OpacityTable % TS, group_id )
    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( "thermEmAb", .true., file_id, group_id )
    CALL WriteOpacityTableTypeAHDF( OpacityTable % thermEmAb, group_id )
    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( "scatt_Iso", .true., file_id, group_id )
    CALL WriteOpacityTableTypeBHDF( OpacityTable % scatt_Iso, group_id )
    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( "scatt_NES", .true., file_id, group_id )
    CALL WriteOpacityTableTypeBHDF( OpacityTable % scatt_NES, group_id )
    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( "scatt_nIso", .true., file_id, group_id )
    CALL WriteOpacityTableTypeCHDF( OpacityTable % scatt_nIso, group_id )
    CALL CloseGroupHDF( group_id )   

    CALL CloseFileHDF( file_id )

  END SUBROUTINE WriteOpacityTableHDF

  SUBROUTINE WriteGridHDF( Grid, group_id )

    TYPE(GridType), INTENT(in)           :: Grid
    INTEGER(HID_T), INTENT(in)                 :: group_id

    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d
    INTEGER                                     :: i
    INTEGER, DIMENSION(1)                       :: buffer

    CHARACTER(LEN=32), DIMENSION(1)             :: tempString
    INTEGER, DIMENSION(1)                       :: tempInteger          

    datasize1d(1) = 1

    tempString(1) = Grid % Name
    CALL WriteHDF( "Name", tempString, &
                             group_id, datasize1d )
    
    tempString(1) = Grid % Unit
    CALL WriteHDF( "Unit", tempString, &
                            group_id, datasize1d )

    tempInteger(1) = Grid % nPoints  
    CALL WriteHDF( "nPoints", tempInteger, &
                            group_id, datasize1d )
   
    tempInteger(1) = Grid % LogInterp 
    CALL WriteHDF( "LogInterp", tempInteger, &
                             group_id, datasize1d )
   
    datasize1d(1) = Grid % nPoints
    CALL WriteHDF( "Values", Grid % Values(:), &
                              group_id, datasize1d )

  END SUBROUTINE WriteGridHDF

  SUBROUTINE WriteOpacityTableTypeAHDF( thermEmAb, group_id )

    TYPE(OpacityTypeA), INTENT(in)              :: thermEmAb
    INTEGER(HID_T), INTENT(in)                  :: group_id

    INTEGER(HSIZE_T)                            :: datasize1d
    INTEGER(HSIZE_T), DIMENSION(3)              :: datasize3d   
    INTEGER(HSIZE_T), DIMENSION(4)              :: datasize4d
    INTEGER                                     :: i
    INTEGER, DIMENSION(1)                       :: buffer

    CHARACTER(LEN=32), DIMENSION(1)             :: tempString
    INTEGER, DIMENSION(1)                       :: tempInteger
    REAL(dp), DIMENSION(1)                      :: tempReal
    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1dtemp
    INTEGER(HID_T)                              :: subgroup_id
    
    datasize1dtemp(1) = 1
    tempInteger(1) = thermEmAb % nOpacities
    CALL WriteHDF&
         ( "nOpacities", tempInteger, group_id, datasize1dtemp )

    datasize1dtemp(1) = 4
    CALL WriteHDF&
         ( "nPoints", thermEmAb % nPoints, group_id, datasize1dtemp )

    datasize1dtemp(1) = thermEmAb % nOpacities
    CALL WriteHDF&
         ( "Names", thermEmAb % Names, group_id, datasize1dtemp ) 

    CALL WriteHDF&
         ( "Species", thermEmAb % Species, group_id, datasize1dtemp ) 

    CALL WriteHDF&
         ( "Units", thermEmAb % Units, group_id, datasize1dtemp ) 

    CALL WriteHDF&
         ( "Offsets", thermEmAb % Offsets, group_id, datasize1dtemp )

    datasize1d = thermEmAb % nOpacities 
    datasize4d = thermEmAb % nPoints

    CALL OpenGroupHDF( "Absorptivity", .true., group_id, subgroup_id )
    DO i = 1, datasize1d
      CALL Write4dHDF_double&
         ( thermEmAb % Names(i), thermEmAb % Absorptivity(i) % Values(:,:,:,:),&
                              subgroup_id, datasize4d )
    END DO
    CALL CloseGroupHDF( subgroup_id )

    datasize3d = thermEmAb % nPoints(2:4)
    CALL OpenGroupHDF( "GreyOpacity_Number_FD", .true., group_id, subgroup_id )
    DO i = 1, datasize1d
    CALL WriteHDF&
         ( thermEmAb % Names(i), thermEmAb % GreyOpacity_Number_FD(i) % Values,&
           subgroup_id, datasize3d )
    END DO
    CALL CloseGroupHDF( subgroup_id )

    CALL OpenGroupHDF( "GreyOpacity_Energy_FD", .true., group_id, subgroup_id )
    DO i = 1, datasize1d
    CALL WriteHDF&
         ( thermEmAb % Names(i), thermEmAb % GreyOpacity_Energy_FD(i) % Values,&
           subgroup_id, datasize3d )
    END DO
    CALL CloseGroupHDF( subgroup_id )

    CALL OpenGroupHDF( "GreyMoment_Energy_FD", .true., group_id, subgroup_id )
    DO i = 1, datasize1d
    CALL WriteHDF&
         ( thermEmAb % Names(i), thermEmAb % GreyMoment_Energy_FD(i) % Values,&
           subgroup_id, datasize3d )
    END DO
    CALL CloseGroupHDF( subgroup_id )

    CALL OpenGroupHDF( "GreyMoment_Number_FD", .true., group_id, subgroup_id )
    DO i = 1, datasize1d
    CALL WriteHDF&
         ( thermEmAb % Names(i), thermEmAb % GreyMoment_Number_FD(i) % Values,&
           subgroup_id, datasize3d )
    END DO
    CALL CloseGroupHDF( subgroup_id )

  END SUBROUTINE WriteOpacityTableTypeAHDF

  SUBROUTINE WriteOpacityTableTypeBHDF( scatt_Iso , group_id )

    TYPE(OpacityTypeB), INTENT(in)              :: scatt_Iso
    INTEGER(HID_T), INTENT(in)                  :: group_id

    INTEGER(HSIZE_T)                            :: datasize1d
    INTEGER(HSIZE_T), DIMENSION(2)              :: datasize2d
    INTEGER(HSIZE_T), DIMENSION(4)              :: datasize4d
    INTEGER(HSIZE_T), DIMENSION(5)              :: datasize5d
    INTEGER                                     :: i
    INTEGER, DIMENSION(1)                       :: buffer

    CHARACTER(LEN=32), DIMENSION(1)             :: tempString
    INTEGER, DIMENSION(1)                       :: tempInteger
    REAL(dp), DIMENSION(1)                      :: tempReal
    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1dtemp
    INTEGER(HID_T)                              :: subgroup_id

    datasize1dtemp(1) = 1
    tempInteger(1) = scatt_Iso % nOpacities
    CALL WriteHDF&
         ( "nOpacities", tempInteger, group_id, datasize1dtemp )

    tempInteger(1) = scatt_Iso % nMoments
    CALL WriteHDF&
         ( "nMoments", tempInteger, group_id, datasize1dtemp )

    datasize1dtemp(1) = 4
    CALL WriteHDF&
         ( "nPoints", scatt_Iso % nPoints, group_id, datasize1dtemp )

    datasize1dtemp(1) = scatt_Iso % nOpacities
    CALL WriteHDF&
         ( "Names", scatt_Iso % Names, group_id, datasize1dtemp )

    CALL WriteHDF&
         ( "Species", scatt_Iso % Species, group_id, datasize1dtemp )

    CALL WriteHDF&
         ( "Units", scatt_Iso % Units, group_id, datasize1dtemp )

    datasize2d = (/scatt_Iso % nOpacities, scatt_Iso % nMoments/)
    CALL WriteHDF&
         ( "Offsets", scatt_Iso % Offsets, group_id, datasize2d )

    datasize1d = scatt_Iso % nOpacities
    datasize5d(1:4) = scatt_Iso % nPoints
    datasize5d(5) = scatt_Iso % nMoments

    CALL OpenGroupHDF( "Kernel", .true., group_id, subgroup_id )
    DO i = 1, datasize1d
      CALL Write5dHDF_double&
         ( scatt_Iso % Names(i), scatt_Iso % Kernel(i) % Values(:,:,:,:,:),&
                              subgroup_id, datasize5d )
    END DO
    CALL CloseGroupHDF( subgroup_id )

    datasize4d(1:3) = scatt_Iso % nPoints(2:4)
    datasize4d(4) = scatt_Iso % nMoments
    CALL OpenGroupHDF( "GreyOpacity_Number_FD", .true., group_id, subgroup_id )
    DO i = 1, datasize1d
    CALL Write4dHDF_double&
         ( scatt_Iso % Names(i), scatt_Iso % GreyOpacity_Number_FD(i) % Values,&
           subgroup_id, datasize4d )
    END DO
    CALL CloseGroupHDF( subgroup_id )

    CALL OpenGroupHDF( "GreyOpacity_Energy_FD", .true., group_id, subgroup_id )
    DO i = 1, datasize1d
    CALL Write4dHDF_double&
         ( scatt_Iso % Names(i), scatt_Iso % GreyOpacity_Energy_FD(i) % Values,&
           subgroup_id, datasize4d )
    END DO
    CALL CloseGroupHDF( subgroup_id )

  END SUBROUTINE WriteOpacityTableTypeBHDF

  SUBROUTINE WriteOpacityTableTypeCHDF( scattn, group_id )

    TYPE(OpacityTypeC), INTENT(in)              :: scattn
    INTEGER(HID_T), INTENT(in)                  :: group_id

  END SUBROUTINE WriteOpacityTableTypeCHDF

  SUBROUTINE Write4dHDF_double &
              ( name, values, group_id, datasize, desc_option, unit_option )

    CHARACTER(*), INTENT(in)                    :: name
    CHARACTER(*), INTENT(in), OPTIONAL          :: unit_option
    CHARACTER(*), INTENT(in), OPTIONAL          :: desc_option
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), DIMENSION(4), INTENT(in)  :: datasize
    REAL(dp), DIMENSION(:,:,:,:), INTENT(in)      :: values
   
    INTEGER(HID_T)                              :: dataset_id
    INTEGER(HID_T)                              :: dataspace_id
    INTEGER(HID_T)                              :: atype_id
    INTEGER(HID_T)                              :: attr_id
    INTEGER(SIZE_T)                             :: attr_len
    INTEGER(HSIZE_T), DIMENSION(1)              :: adims = (/1/)
   
    CALL h5screate_simple_f( 4, datasize, dataspace_id, hdferr )

    CALL h5dcreate_f( group_id, name, H5T_NATIVE_DOUBLE, &
           dataspace_id, dataset_id, hdferr )

    CALL h5dwrite_f( dataset_id, H5T_NATIVE_DOUBLE, &
           values, datasize, hdferr )

    CALL h5sclose_f( dataspace_id, hdferr ) 

    CALL h5dclose_f( dataset_id, hdferr )

  END SUBROUTINE Write4dHDF_double

  SUBROUTINE Write5dHDF_double &
              ( name, values, group_id, datasize, desc_option, unit_option )

    CHARACTER(*), INTENT(in)                    :: name
    CHARACTER(*), INTENT(in), OPTIONAL          :: unit_option
    CHARACTER(*), INTENT(in), OPTIONAL          :: desc_option
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), DIMENSION(5), INTENT(in)  :: datasize
    REAL(dp), DIMENSION(:,:,:,:,:), INTENT(in)  :: values

    INTEGER(HID_T)                              :: dataset_id
    INTEGER(HID_T)                              :: dataspace_id
    INTEGER(HID_T)                              :: atype_id
    INTEGER(HID_T)                              :: attr_id
    INTEGER(SIZE_T)                             :: attr_len
    INTEGER(HSIZE_T), DIMENSION(1)              :: adims = (/1/)

    CALL h5screate_simple_f( 5, datasize, dataspace_id, hdferr )

    CALL h5dcreate_f( group_id, name, H5T_NATIVE_DOUBLE, &
           dataspace_id, dataset_id, hdferr )

    CALL h5dwrite_f( dataset_id, H5T_NATIVE_DOUBLE, &
           values, datasize, hdferr )

    CALL h5sclose_f( dataspace_id, hdferr )

    CALL h5dclose_f( dataset_id, hdferr )

  END SUBROUTINE Write5dHDF_double


  SUBROUTINE ReadOpacityTableHDF( OpacityTable, FileName )
 
    TYPE(OpacityTableType), INTENT(inout)   :: OpacityTable
    CHARACTER(len=*), INTENT(in)            :: FileName

    INTEGER, DIMENSION(3)                         :: nPointsTS
    INTEGER                                       :: nPointsE
    INTEGER                                       :: nPointsEta
    INTEGER                                       :: nOpacA
    INTEGER                                       :: nOpacB, nMomB
    INTEGER                                       :: nOpacB_NES, nMomB_NES
    INTEGER                                       :: nOpacC, nMomC
    INTEGER(HID_T)                                :: file_id
    INTEGER(HID_T)                                :: group_id
    INTEGER(HID_T)                                :: subgroup_id
    INTEGER(HSIZE_T), DIMENSION(1)                :: datasize1d
    INTEGER, DIMENSION(1)                         :: buffer
    CHARACTER(LEN=32), DIMENSION(1)               :: buffer_string

    CALL OpenFileHDF( FileName, .false., file_id )

    datasize1d(1) = 1
    CALL ReadHDF( "nOpacitiesA", buffer, file_id, datasize1d )
    nOpacA = buffer(1)   

    CALL ReadHDF( "nOpacitiesB", buffer, file_id, datasize1d )
    nOpacB = buffer(1)

    CALL ReadHDF( "nMomentsB", buffer, file_id, datasize1d )
    nMomB = buffer(1)

    CALL ReadHDF( "nOpacitiesB_NES", buffer, file_id, datasize1d )
    nOpacB_NES = buffer(1)

    CALL ReadHDF( "nMomentsB_NES", buffer, file_id, datasize1d )
    nMomB_NES = buffer(1)

    CALL ReadHDF( "nOpacitiesC", buffer, file_id, datasize1d )
    nOpacC = buffer(1)

    CALL ReadHDF( "nMomentsC", buffer, file_id, datasize1d )
    nMomC = buffer(1)

    CALL ReadHDF( "nPointsE", buffer, file_id, datasize1d )
    nPointsE = buffer(1)

    CALL ReadHDF( "nPointsEta", buffer, file_id, datasize1d )
    nPointsEta = buffer(1)

    CALL ReadHDF( "nPointsTS", nPointsTS, file_id, datasize1d )

    CALL AllocateOpacityTable &
               ( OpacityTable, nOpacA, nOpacB, nMomB, nOpacB_NES, nMomB_NES, nOpacC, nMomC, nPointsE, nPointsEta )  

    IF( ( OpacityTable % EOSTable % TS % nPoints(1) .EQ. nPointsTS(1) ) .AND. &
        ( OpacityTable % EOSTable % TS % nPoints(2) .EQ. nPointsTS(2) ) .AND. &
        ( OpacityTable % EOSTable % TS % nPoints(3) .EQ. nPointsTS(3) ) ) THEN

    CALL OpenGroupHDF( "EnergyGrid", .false., file_id, group_id )
    CALL ReadGridHDF( OpacityTable % EnergyGrid, group_id )
    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( "EtaGrid", .false., file_id, group_id )
    CALL ReadGridHDF( OpacityTable % EtaGrid, group_id )
    CALL CloseGroupHDF( group_id )
 
    CALL ReadThermoStateHDF( OpacityTable % TS, file_id )

    CALL OpenGroupHDF( "thermEmAb", .false., file_id, group_id )
    CALL ReadOpacityTypeAHDF( OpacityTable % thermEmAb, group_id )
    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( "scatt_Iso", .false., file_id, group_id )
    CALL ReadOpacityTypeBHDF( OpacityTable % scatt_Iso, group_id )
    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( "scatt_NES", .false., file_id, group_id )
    CALL ReadOpacityTypeBHDF( OpacityTable % scatt_NES, group_id )
    CALL CloseGroupHDF( group_id )

    CALL OpenGroupHDF( "scatt_nIso", .false., file_id, group_id )
    CALL ReadOpacityTypeCHDF( OpacityTable % scatt_nIso , group_id )
    CALL CloseGroupHDF( group_id )

    CALL CloseFileHDF( file_id )
    
    ELSE 
      WRITE(*,*) "ERROR!"
      WRITE(*,*) "EquationOfStateTable is not consistent with OpacityTable!"
      STOP
    END IF

  END SUBROUTINE ReadOpacityTableHDF


  SUBROUTINE ReadOpacityTypeAHDF( thermEmAb, group_id )

    TYPE(OpacityTypeA),INTENT(inout)                 :: thermEmAb
    INTEGER(HID_T), INTENT(in)                       :: group_id

    INTEGER(HSIZE_T), DIMENSION(1)                   :: datasize1d
    INTEGER(HSIZE_T), DIMENSION(4)                   :: datasize4d
    INTEGER                                          :: i
    INTEGER, DIMENSION(1)                            :: buffer
    REAL(dp), DIMENSION(1)                           :: bufferReal
    INTEGER(HID_T)                                   :: subgroup_id

    datasize1d(1) = 1
    CALL ReadHDF( "nOpacities", buffer, group_id, datasize1d )
    thermEmAb % nOpacities = buffer(1)

    datasize1d = buffer(1)
    CALL ReadHDF( "Offsets", thermEmAb % Offsets, group_id, datasize1d )

    Call ReadHDF( "Names", thermEmAb % Names, group_id, datasize1d )

    Call ReadHDF( "Units", thermEmAb % Units, group_id, datasize1d )

    CALL ReadHDF( "Species", thermEmAb % Species, group_id, datasize1d )

    datasize1d(1) = 4
    CALL ReadHDF( "nPoints", thermEmAb % nPoints, group_id, datasize1d )

    datasize4d = thermEmAb % nPoints

    CALL OpenGroupHDF( "Absorptivity", .false., group_id, subgroup_id )
    DO i = 1, thermEmAb % nOpacities
      CALL Read4dHDF_double&
             ( thermEmAb % Names(i), thermEmAb % Absorptivity(i) % Values,&
               subgroup_id, datasize4d )
    END DO ! nOpacities
    CALL CloseGroupHDF( subgroup_id )

    CALL OpenGroupHDF( "GreyOpacity_Number_FD", .false., group_id, subgroup_id )
    DO i = 1, thermEmAb % nOpacities
      CALL ReadHDF&
             ( thermEmAb % Names(i), &
               thermEmAb % GreyOpacity_Number_FD(i) % Values, &
               subgroup_id, datasize4d )
    END DO ! nOpacities
    CALL CloseGroupHDF( subgroup_id )

    CALL OpenGroupHDF( "GreyOpacity_Energy_FD", .false., group_id, subgroup_id )
    DO i = 1, thermEmAb % nOpacities
      CALL ReadHDF&
             ( thermEmAb % Names(i), &
               thermEmAb % GreyOpacity_Energy_FD(i) % Values, &
               subgroup_id, datasize4d )
    END DO ! nOpacities
    CALL CloseGroupHDF( subgroup_id )

    CALL OpenGroupHDF( "GreyMoment_Number_FD", .false., group_id, subgroup_id )
    DO i = 1, thermEmAb % nOpacities
      CALL ReadHDF&
             ( thermEmAb % Names(i), &
               thermEmAb % GreyMoment_Number_FD(i) % Values, &
               subgroup_id, datasize4d )
    END DO ! nOpacities
    CALL CloseGroupHDF( subgroup_id )

    CALL OpenGroupHDF( "GreyMoment_Energy_FD", .false., group_id, subgroup_id )
    DO i = 1, thermEmAb % nOpacities
      CALL ReadHDF&
             ( thermEmAb % Names(i), &
               thermEmAb % GreyMoment_Energy_FD(i) % Values, &
               subgroup_id, datasize4d )
    END DO ! nOpacities
    CALL CloseGroupHDF( subgroup_id )

  END SUBROUTINE ReadOpacityTypeAHDF


  SUBROUTINE ReadOpacityTypeBHDF( scatt_Iso, group_id )

    TYPE(OpacityTypeB),INTENT(inout)                 :: scatt_Iso
    INTEGER(HID_T), INTENT(in)                       :: group_id

    INTEGER(HSIZE_T), DIMENSION(1)                   :: datasize1d
    INTEGER(HSIZE_T), DIMENSION(2)                   :: datasize2d
    INTEGER(HSIZE_T), DIMENSION(4)                   :: datasize4d
    INTEGER(HSIZE_T), DIMENSION(5)                   :: datasize5d
    INTEGER                                          :: i
    INTEGER, DIMENSION(1)                            :: buffer
    REAL(dp), DIMENSION(1)                           :: bufferReal
    INTEGER(HID_T)                                   :: subgroup_id

    datasize1d(1) = 1
    CALL ReadHDF( "nOpacities", buffer, group_id, datasize1d )
    scatt_Iso % nOpacities = buffer(1)

    CALL ReadHDF( "nMoments", buffer, group_id, datasize1d )
    scatt_Iso % nMoments   = buffer(1)

    datasize1d = buffer(1)
    Call ReadHDF( "Names", scatt_Iso % Names, group_id, datasize1d )

    Call ReadHDF( "Units", scatt_Iso % Units, group_id, datasize1d )

    CALL ReadHDF( "Species", scatt_Iso % Species, group_id, datasize1d )

    datasize1d(1) = 4
    CALL ReadHDF( "nPoints", scatt_Iso % nPoints, group_id, datasize1d )

    datasize2d = (/scatt_Iso % nOpacities, scatt_Iso % nMoments/)
    CALL ReadHDF( "Offsets", scatt_Iso % Offsets, group_id, datasize2d )

    datasize5d(1:4) = scatt_Iso % nPoints
    datasize5d(5) = scatt_Iso % nMoments

    CALL OpenGroupHDF( "Kernel", .false., group_id, subgroup_id )
    DO i = 1, scatt_Iso % nOpacities
      CALL Read5dHDF_double &
             ( scatt_Iso % Names(i), &
               scatt_Iso % Kernel(i) % Values, &
               subgroup_id, datasize5d )
    END DO ! nOpacities
    CALL CloseGroupHDF( subgroup_id )

    datasize4d(1:3) = scatt_Iso % nPoints(2:4)
    datasize4d(4)   = scatt_Iso % nMoments
    CALL OpenGroupHDF( "GreyOpacity_Energy_FD", .false., group_id, subgroup_id )
    DO i = 1, scatt_Iso % nOpacities
      CALL Read4dHDF_double&
             ( scatt_Iso % Names(i), &
               scatt_Iso % GreyOpacity_Energy_FD(i) % Values, &
               subgroup_id, datasize4d )
    END DO ! nOpacities
    CALL CloseGroupHDF( subgroup_id )

    CALL OpenGroupHDF( "GreyOpacity_Number_FD", .false., group_id, subgroup_id )
    DO i = 1, scatt_Iso % nOpacities
      CALL Read4dHDF_double&
             ( scatt_Iso % Names(i), &
               scatt_Iso % GreyOpacity_Number_FD(i) % Values,&
               subgroup_id, datasize4d )
    END DO ! nOpacities
    CALL CloseGroupHDF( subgroup_id )

  END SUBROUTINE ReadOpacityTypeBHDF

  SUBROUTINE ReadOpacityTypeCHDF( thermEmAb, group_id )

    TYPE(OpacityTypeC),INTENT(inout)                 :: thermEmAb
    INTEGER(HID_T), INTENT(in)                       :: group_id

  END SUBROUTINE ReadOpacityTypeCHDF

  SUBROUTINE Read4dHDF_double( name, values, group_id, datasize )

    CHARACTER(*), INTENT(in)                    :: name
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), DIMENSION(4), INTENT(in)  :: datasize
    REAL(dp), DIMENSION(:,:,:,:), INTENT(out)   :: values
   
    INTEGER(HID_T)                               :: dataset_id
 
    CALL h5dopen_f( group_id, name, dataset_id, hdferr )
    CALL h5dread_f( dataset_id, H5T_NATIVE_DOUBLE, values, datasize, hdferr )
    CALL h5dclose_f( dataset_id, hdferr )

  END SUBROUTINE Read4dHDF_double


  SUBROUTINE Read5dHDF_double( name, values, group_id, datasize )

    CHARACTER(*), INTENT(in)                    :: name
    INTEGER(HID_T)                              :: group_id
    INTEGER(HSIZE_T), DIMENSION(5), INTENT(in)  :: datasize
    REAL(dp), DIMENSION(:,:,:,:,:), INTENT(out) :: values

    INTEGER(HID_T)                               :: dataset_id

    CALL h5dopen_f( group_id, name, dataset_id, hdferr )
    CALL h5dread_f( dataset_id, H5T_NATIVE_DOUBLE, values, datasize, hdferr )
    CALL h5dclose_f( dataset_id, hdferr )

  END SUBROUTINE Read5dHDF_double


  SUBROUTINE ReadGridHDF( Grid, group_id )

    TYPE(GridType), INTENT(inout)               :: Grid
    INTEGER(HID_T), INTENT(in)                  :: group_id

    INTEGER(HSIZE_T), DIMENSION(1)              :: datasize1d
    INTEGER, DIMENSION(1)                       :: buffer
    CHARACTER(LEN=32), DIMENSION(1)             :: buffer_string

    datasize1d(1) = 1
    Call ReadHDF( "Name", buffer_string, group_id, datasize1d )
    Grid % Name = buffer_string(1)

    Call ReadHDF( "Unit", buffer_string, group_id, datasize1d )
    Grid % Unit = buffer_string(1)

    CALL ReadHDF( "nPoints", buffer, group_id, datasize1d )
    Grid % nPoints = buffer(1)

    CALL ReadHDF( "LogInterp", buffer, group_id, datasize1d )
    Grid % LogInterp = buffer(1)
 
    datasize1d = Grid % nPoints
    CALL ReadHDF( "Values", Grid % Values, &
                              group_id, datasize1d )

    Grid % minValue = MINVAL( Grid % Values )
    
    Grid % maxValue = MAXVAL( Grid % Values )

  END SUBROUTINE ReadGridHDF


END MODULE wlOpacityTableIOModuleHDF
