module mts_utils

   use kinds, only: real_8
   use parac, only: paral, parai
   use system, only: cntl, cnti, cntr
   use cp_xc_utils, only: cp_xc, cp_xc_functional, cp_xc_mts_low_func, cp_xc_functional_env, cp_xc_mts_low_func_env
   use error_handling, only: stopgm
   use inscan_utils, only: inscan
   use readsr_utils, only: keyword_contains, input_string_len, xstring
   use mp_interface, only: mp_bcast
  USE func,                            ONLY: func1,&
                                             func2,&
                                             func3

#include "sizeof.h"

   implicit none

   type :: mts_t
      integer :: timestep_factor
      character(len=10) :: low_level ! string defining the model for high level forces 
      character(len=10) :: high_level ! string defining the model for low level forces  
      character(len=50) :: low_dft_func
      character(len=50) :: high_dft_func
      logical :: print_forces
   end type mts_t

   type (mts_t),public ::  mts

   public :: read_mts_input
   public :: set_mts_functional

   integer, parameter :: output_unit = 6

contains


   ! Purpose: Reads the section &MTS ... &END from file with unit IUNIT
   ! Author: Pablo Baudin (based on dftin routine)
   ! Date: May 2018
   subroutine read_mts_input
      ! format:
      !
      ! &MTS
      !    TIMESTEP_FACTOR
      !       mts%timestep_factor
      !    PRINT_FORCES [OFF]
      !    LOW_LEVEL_FORCES {DFT, EXTERNAL}
      !    HIGH_LEVEL_FORCES {DFT, EXTERNAL}
      ! &END
      !
      CHARACTER(*), PARAMETER                  :: procedureN = 'read_mts_input'
      INTEGER, PARAMETER                       :: max_unknown_lines = 30 

      CHARACTER(len=input_string_len)          :: line, previous_line
      CHARACTER(len=input_string_len)          :: error_message, unknown(max_unknown_lines) 
      INTEGER                                  :: i, first, last, ierr, iunit, nbr_unknown_lines
      LOGICAL                                  :: something_went_wrong, go_on_reading

      !
      ! The read loop is only accessed by the io_parent, therefore, the
      ! io_parent checks within the loop have been removed (redundancy).
      !
      IF (paral%io_parent) THEN
         iunit = 5
         !
         ! Variables for reading
         !
         nbr_unknown_lines     = 0
         line                  = ' '
         previous_line         = ' '
         error_message         = ' '
         go_on_reading         = .true.
         something_went_wrong  = .false.
         !
         ! Defaults
         !
         mts%timestep_factor = 1
         mts%low_level = 'DFT'
         mts%high_level = 'DFT'
         mts%print_forces = .false.
         mts%low_dft_func = ' '
         mts%high_dft_func = ' '
         !
         ! If the MTS section is not there, we simply move on.
         !
         ierr = inscan(iunit,'&MTS')
         !
         IF (ierr == 0) THEN
            !
            ! Check that MTS is actually switched on
            if (.not. cntl%use_mts) then
               go_on_reading = .false.
               write(output_unit,*) ' WARNING: MTS section will be ignored:', &
                  ' USE_MTS keyword missing in CPMD section!'
            end if
            !
            ! Main loop
            !
            DO WHILE(go_on_reading)
               !
               ! Read a line and store the old one
               !
               previous_line = line
               READ(iunit,'(A80)',iostat=ierr) line
               IF (ierr /= 0) THEN
                  something_went_wrong = .TRUE.
                  go_on_reading        = .FALSE.

               ELSE IF (keyword_contains(line,'&END')) THEN
                  go_on_reading = .FALSE.

               ELSE IF ( keyword_contains(line,'TIMESTEP_FACTOR') ) THEN
                  READ(iunit,*,iostat=ierr) mts%timestep_factor
                  if (ierr /= 0) then
                     something_went_wrong = .TRUE.
                     go_on_reading = .FALSE.
                     error_message = 'Could not read TIMESTEP_FACTOR'
                  end if

               ELSE IF ( keyword_contains(line,'PRINT_FORCES') ) THEN
                  IF ( keyword_contains(line,'OFF') ) THEN
                     mts%print_forces = .false.
                  ELSE
                     mts%print_forces = .true.
                  ENDIF

               ELSE IF ( keyword_contains(line,'LOW_LEVEL_FORCES') ) THEN

                  mts%low_level = read_mts_input_level(line)

                  if (mts%low_level == ' ') then
                     something_went_wrong = .TRUE.
                     go_on_reading = .FALSE.
                     error_message = 'Unknown MTS model for low level forces'
                  end if

               ELSE IF ( keyword_contains(line,'HIGH_LEVEL_FORCES') ) THEN

                  mts%high_level = read_mts_input_level(line)

                  if (mts%high_level == ' ') then
                     something_went_wrong = .TRUE.
                     go_on_reading = .FALSE.
                     error_message = 'Unknown MTS model for high level forces'
                  end if

               ELSE
                  ! Unknown Keyword. store and continue
                  IF (' ' /= line) THEN
                     IF (nbr_unknown_lines < max_unknown_lines) THEN
                        nbr_unknown_lines = nbr_unknown_lines+1
                        unknown(nbr_unknown_lines) = line
                     ELSE
                        DO i=1,max_unknown_lines-1
                           unknown(i) = unknown(i+1)
                        ENDDO
                        unknown(nbr_unknown_lines) = line
                     ENDIF
                  ENDIF

               ENDIF
            END DO ! reading

         ENDIF


         ! print error message
         IF (something_went_wrong) THEN
            WRITE(output_unit,'(/,1X,64("!"))')
            WRITE(output_unit,'(1X, A, 1X, A)') 'ERROR:', 'PROBLEM WHILE READING &MTS SECTION:' 
            WRITE(output_unit,'(8X, A)') TRIM(ADJUSTL(error_message))
            IF (line /= ' ' .or. previous_line /= ' ') THEN
               WRITE(output_unit,'(8X, A)') 'THE LAST TWO LINES READ WITHIN THE SECTION WERE:'
               WRITE(output_unit,'(/,1X, A)') TRIM(ADJUSTL(previous_line))
               WRITE(output_unit,'(1X, A)') TRIM(ADJUSTL(line))
            END IF
            WRITE(output_unit,'(1X,64("!"))')
            CALL stopgm(procedureN,'Error while reading &MTS section, cf. output file',& 
               __LINE__,__FILE__)
         ENDIF


         ! print list of unknown keywords
         IF (nbr_unknown_lines /= 0) THEN
            WRITE(output_unit,'(/,1X,64("="))')
            WRITE(output_unit,'(1X,A,14X,A,15X,A)') '= ','UNKNOWN KEYWORDS IN SECTION &MTS','='
            DO i=1,nbr_unknown_lines
               previous_line = unknown(i)
               CALL xstring(previous_line,first,last)
               WRITE(output_unit,'(A,A)') ' = ',previous_line(first:last)
            ENDDO
            WRITE(output_unit,'(1X,64("="),/)')
         ENDIF


         ! check for consitency problems
         IF(cntl%tdipd.AND.(mod(cnti%npdip,mts%timestep_factor).NE.0)) THEN
            WRITE(6,'(1x,A,i5)') 'WARNING MTS: CHANGING DIPOLE SAMPLE VALUE:',cnti%npdip
            cnti%npdip=mts%timestep_factor*(cnti%npdip/mts%timestep_factor+1)
            WRITE(6,'(5x,A,T50,I10,A6)') 'NOW STORE DIPOLE MOMENTS EVERY ',cnti%npdip,' STEPS'
         ENDIF


         ! print MTS input information to CPMD output
         if (cntl%use_mts) then
            write(output_unit,'(/,1x,a)') 'MULTIPLE TIME STEP SCHEME:'
            write(output_unit,'(4x,a,t40,i4)') 'TIMESTEP FACTOR:', mts%timestep_factor
            write(output_unit,'(4x,a,t40,a)')  'LOW LEVEL FROM :', mts%low_level
            write(output_unit,'(4x,a,t40,a)')  'HIGH LEVEL FROM:', mts%high_level
            if (mts%print_forces) then
               write(output_unit,'(4x,a)')  'HIGH AND LOW LEVEL FORCES WILL BE PRINTED TO FILE'
            end if
         end if


      END IF ! select only master node

      CALL mp_bcast(cnti%npdip,parai%io_source,parai%cp_grp)
      CALL mp_bcast_byte(mts, size_in_bytes_of(mts),parai%io_source,parai%cp_grp)

   end subroutine read_mts_input


   function read_mts_input_level(line) result(level)
      implicit none
      character(len=input_string_len), intent(in) :: line
      character(len=10) :: level 

      ! This is the place to add new options for the calculation 
      ! of the forces in the MTS scheme
      if (keyword_contains(line,'DFT')) then
         level = 'DFT'
      else if (keyword_contains(line,'EXTERNAL')) then
         level = 'EXTERNAL'
      else
         level = ' '
      end if
      
   end function read_mts_input_level

   ! Purpose: set functional to be high or low level depending on input
   ! Author: Pablo Baudin
   ! Date: May 2018
   subroutine set_mts_functional(level)
      implicit none
      character(len=*), intent(in) :: level
      character(len=*), parameter :: proceduren = 'set_mts_functional'

      if (cntl%use_xc_driver) then

         select case(level) 
         case('HIGH')
            ! set current functional to the high level one
            cp_xc => cp_xc_functional

         case('LOW')
            ! set current functional to the low level one
            cp_xc => cp_xc_mts_low_func

         case default
            call stopgm(proceduren,'Invalid level request options are HIGH or LOW',&
               __LINE__,__FILE__)

         end select

         CALL cp_xc%get( tgc=cntl%tgc, tgcx=cntl%tgcx, tgcc=cntl%tgcc, ttau=cntl%ttau, &
            thybrid=cntl%thybrid, mhfx=func1%mhfx, phfx=func3%phfx, &
            msrx=func1%msrx, srxa=func2%srxa, cam_alpha=func2%cam_alpha, cam_beta=func2%cam_beta )
      else

         call stopgm(procedureN,' MTS: functionals have to be set up with XC_DRIVER',&
            __LINE__,__FILE__)

      end if

   end subroutine set_mts_functional


end module mts_utils
