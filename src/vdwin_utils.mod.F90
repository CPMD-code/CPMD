#ifdef __SR8000
!option MP(P(0)), LANGLVL(SAVE(0))
#endif

MODULE vdwin_utils
  USE adat,                            ONLY: covrad,&
                                             elem
  USE cnst,                            ONLY: fbohr
  USE dcacp_utils,                     ONLY: dcacp_env,&
                                             dcacp_custom_tag,&
                                             dcacp_init
  USE dftd3_api,                       ONLY: dftd3_init,&
                                             dftd3_set_functional
  USE error_handling,                  ONLY: stopgm
  USE inscan_utils,                    ONLY: inscan
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: parai,&
                                             paral
  USE readsr_utils,                    ONLY: index_of_delimiter,&
                                             input_string_len,&
                                             keyword_contains,&
                                             readsi,&
                                             readsr,&
                                             xstring
  USE system,                          ONLY: cntl,&
                                             maxsys
  USE timer,                           ONLY: tiset,&
                                             tihalt
  USE vdwcmod,                         ONLY: &
       boadwf, empvdwc, empvdwi, empvdwr, idvdw, ifragdata, ivdw, jvdw, &
       nfragx, radfrag, vdwbe, vdwi, vdwl, vdwrm, vdwst, vdwwfi, vdwwfl, &
       vdwwfr, tdftd3, tdftd3noc, dftd3f, dftd3i, dftd3c
  USE wann,                            ONLY: wannl
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: output_unit       =  6
  INTEGER, PARAMETER :: max_unknown_lines = 30

  PUBLIC :: vdwin

CONTAINS

  ! ==================================================================
  SUBROUTINE vdwin
    ! ==--------------------------------------------------------------==
    ! ==  THIS ROUTINE READS THE SECTION &VDW &END ON UNIT IUNIT      ==
    ! ==--------------------------------------------------------------==
    ! ==  THE INPUT HAS TO BE IN THE FOLLOWING FORMAT                 ==
    ! ==     &VDW                                                     ==
    ! ==        Options                                               ==
    ! ==     &END                                                     ==
    ! ==--------------------------------------------------------------==
    ! ==  LIST OF OPTIONS                                             ==
    ! ==                                                              ==
    ! ==    EMPIRICAL CORRECTION                                      ==
    ! ==    WANNIER CORRECTION                                        ==
    ! ==    DCACP Z=z          ! Set custom parameter for DCACP       ==
    ! ==     sigma_1, sigma_2                                         ==
    ! ==    INCLUDE_METALS     ! Apply DCACP to metal centres, too    ==
    ! ==    NO_CONTRIBUTION    ! Do not apply DCACP to species listed ==
    ! ==     species_1_w/o_contrib species_2_w/o_contrib ...          ==
    ! ==                                                              ==
    ! ==  Added support for DCACP.                                    ==
    ! ==                          12.06.2017 M.P.Bircher @ LCBC/EPFL  ==
    ! ==  Present release:                                            ==
    ! ==             Strasbourg/Tokyo/Zurich/Hyogo, 6 May 2019        ==
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'vdwin'

    CHARACTER(len=input_string_len)          :: line, previous_line, &
                                                error_message, &
                                                unknown(max_unknown_lines)
    INTEGER                                  :: i, ierr, isub, iunit, isp, first, last, &
                                                nbr_unknown_lines, z
    LOGICAL                                  :: erread, go_on_reading, &
                                                something_went_wrong

    !
    ! Set routine tracking for proper call-trees in DCACP-init
    !
    CALL tiset(procedureN,isub)

    IF (vdwl%vdwc .OR. vdwl%vdwd .OR. vdwl%dcacp) THEN
       IF ( (vdwl%vdwc.AND.vdwl%vdwd) .OR. &
            (vdwl%vdwc.AND.vdwl%dcacp) .OR. &
            (vdwl%vdwd.AND.vdwl%dcacp) ) THEN
          WRITE(output_unit,'(/,A,/,A,/,A,/)')&
               ' THE FOLLOWING OPTIONS ARE MUTUALLY  EXCLUSIVE ',&
               '                                         DCACP ',&
               '                          EMPIRICAL CORRECTION ',&
               '                            WANNIER CORRECTION '
          CALL stopgm(procedureN,'Conflicting run options.',& 
               __LINE__,__FILE__)
       ENDIF
       IF (paral%io_parent) THEN
          iunit = 5
          !
          nbr_unknown_lines = 0
          line              = ' ' 
          previous_line     = ' '
          error_message     = ' '
          !
          ! Search for &VDW section
          !
          ierr  = inscan(iunit,'&VDW')
          !
          IF (ierr == 0) THEN
             !
             ! Main loop
             !
             go_on_reading        = .true.
             something_went_wrong = .false.
             !
             DO WHILE(go_on_reading)
                previous_line = line
                READ(iunit,'(A80)',iostat=ierr) line
                IF (ierr /= 0) THEN
                   something_went_wrong = .TRUE.
                   error_message        = 'ERROR WHILE READING LINE'
                   go_on_reading        = .FALSE.
                ELSEIF ( keyword_contains(line,'&END') ) THEN
                   go_on_reading        = .FALSE.
                ELSEIF (vdwl%vdwc .AND. &
                        keyword_contains(line,'EMPIRICAL',and='CORRECTION') ) THEN
                   ! Empirical van der Waals corrections
                   CALL empvdwin(iunit)
                ELSEIF (vdwl%vdwd .AND. &
                        keyword_contains(line,'WANNIER',and='CORRECTION') ) THEN
                   ! vdW corrections using Wannier functions
                   CALL wannvdwin(iunit)
                !
                ! DCACP
                !
                ELSEIF (vdwl%dcacp .AND. &
                        keyword_contains(line,'NO_CONTRIBUTION') ) THEN
                   READ(iunit,'(A80)',iostat=ierr) line
                   IF (ierr /= 0) THEN
                      error_message        = 'ERROR WHILE READING NO_CONTRIBUTION LIST'
                      something_went_wrong = .TRUE.
                      go_on_reading        = .FALSE.
                   ENDIF
                   first=1
                   contrib_loop: DO 
                      erread=.FALSE.
                      CALL readsi(line,first,last,isp,erread)
                      IF (erread) EXIT contrib_loop
                      dcacp_env%no_contribution(isp) = .TRUE.
                      first = last
                   ENDDO contrib_loop
                ELSE IF(vdwl%dcacp .AND. &
                        keyword_contains(line,'INCLUDE_METALS') ) THEN
                   dcacp_env%include_metals = .TRUE.
                ELSE IF(vdwl%dcacp .AND. &
                        keyword_contains(line,'DCACP',and='Z',cut_at='=') ) THEN
                   dcacp_env%has_custom_sigma = .TRUE.
                   first = index_of_delimiter(line,'Z','=')
                   CALL readsi(line,first,last,z,erread)
                   IF (erread) THEN
                      error_message        = 'ERROR WHILE READING CORE CHARGE VALUE Z'
                      something_went_wrong = .TRUE.
                      go_on_reading        = .FALSE.
                   ENDIF
                   previous_line  = line
                   READ(iunit,'(A80)',iostat=ierr) line
                   first = 1
                   CALL readsr(line,first,last,dcacp_env%data_from_input(z)%sigma_1,erread)
                   first = last
                   IF (.NOT. erread) CALL readsr(line,first,last,dcacp_env%data_from_input(z)%sigma_2,erread)
                   IF (erread) THEN
                      error_message        = 'ERROR WHILE READING DCACP PARAMETERS'
                      something_went_wrong = .TRUE.
                      go_on_reading        = .FALSE.
                   ENDIF
                   dcacp_env%data_from_input(z)%citation = dcacp_custom_tag
                ELSE
                   ! Unknown keyword
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
             ENDDO
          ELSE
             !
             ! If the section &VDW is not there (this is legit)  
             !
             something_went_wrong = .false.
          ENDIF
          !
          IF (something_went_wrong) THEN
             WRITE(output_unit,'(/,1X,64("!"))')
             WRITE(output_unit,'(1X, A, 1X, A)') 'ERROR:', 'PROBLEM WHILE READING &VDW SECTION:' 
             WRITE(output_unit,'(8X, A)') TRIM(ADJUSTL(error_message))
             IF (line /= ' ' .or. previous_line /= ' ') THEN
                WRITE(output_unit,'(8X, A)') 'THE LAST TWO LINES READ WITHIN THE SECTION WERE:'
                WRITE(output_unit,'(/,1X, A)') TRIM(ADJUSTL(previous_line))
                WRITE(output_unit,'(1X, A)') TRIM(ADJUSTL(line))
             END IF
             WRITE(output_unit,'(1X,64("!"))')
             CALL stopgm(procedureN,'Error while reading &VDW section, cf. output file',& 
                  __LINE__,__FILE__)
          ENDIF
          !
          IF (nbr_unknown_lines /= 0) THEN
             WRITE(output_unit,'(/,1X,64("="))')
             WRITE(output_unit,'(1X,A,14X,A,14X,A)') '= ','UNKNOWN KEYWORDS IN SECTION &VDW','='
             DO i=1,nbr_unknown_lines
                previous_line = unknown(i)
                CALL xstring(previous_line,first,last)
                WRITE(output_unit,'(A,A)') ' = ',previous_line(first:last)
             ENDDO
             WRITE(output_unit,'(1X,64("="),/)')
          ENDIF
          !
       ENDIF
       IF (vdwl%dcacp) THEN
          CALL dcacp_init()
       ELSEIF (vdwl%vdwc) THEN
          CALL empvdw_init()
       ELSEIF (vdwl%vdwd) THEN
          CALL vdwwf_init()
       ENDIF
    ENDIF

    CALL tihalt(procedureN,isub)

    CONTAINS

    ! ==--------------------------------------------------------------==
    SUBROUTINE empvdw_init()
 
       ! First broadcast logical variable tdftd3
       CALL mp_bcast(tdftd3,parai%io_source,parai%cp_grp)
       IF (tdftd3) THEN
          ! Broadcast dftd3 stuff
          CALL mp_bcast(tdftd3noc,SIZE(tdftd3noc),parai%io_source,parai%cp_grp)
          CALL mp_bcast(dftd3i%threebody,parai%io_source,parai%cp_grp)
          CALL mp_bcast(dftd3i%numgrad,parai%io_source,parai%cp_grp)
          CALL mp_bcast(dftd3i%cutoff,parai%io_source,parai%cp_grp)
          CALL mp_bcast(dftd3i%cutoff_cn,parai%io_source,parai%cp_grp)
          CALL mp_bcast_byte(dftd3f,size_in_bytes_of(dftd3f),parai%io_source,parai%cp_grp)
          ! vdwi is also needed for dftd3
          CALL mp_bcast_byte(vdwi,size_in_bytes_of(vdwi),parai%io_source,parai%cp_grp)
          ! Initialize dftd3c
          CALL dftd3_init(dftd3c,dftd3i)
          ! Choose functional
          CALL dftd3_set_functional(dftd3c,func=TRIM(ADJUSTL(dftd3f%func)),version=dftd3f%version,tz=dftd3f%tz)
          ! Report selected scaling & damping parameters
          IF (paral%io_parent) THEN
             WRITE(output_unit,'(A,T39,2(A,F8.4))')&
                  '    GLOBAL SCALING FACTORS:',&
                  '  s6=',dftd3c%s6,'  s18=',dftd3c%s18
             SELECT CASE(dftd3f%version)
             CASE (2)
                WRITE(output_unit,'(A,T39,2(A,F8.4))')&
                  '    FERMI-TYPE DAMPING:',&
                  ' rs6=',dftd3c%rs6,' alp6=',dftd3c%alp
             CASE (3,5)
                WRITE(output_unit,'(A,T39,2(A,F8.4)/T38,2(A,F8.4))')&
                  '    ZERO DAMPING:',&
                  ' rs6=',dftd3c%rs6,' alp6=',dftd3c%alp,&
                  ' rs18=',dftd3c%rs18,' alp8=',dftd3c%alp+2.0_real_8
             CASE (4,6)
                WRITE(output_unit,'(A,T39,2(A,F8.4)/T38,2(A,F8.4))')&
                  '    BECKE-JOHNSON DAMPING:',&
                  ' rs6=',dftd3c%rs6,' alp6=',dftd3c%alp,&
                  ' rs18=',dftd3c%rs18,' alp8=',dftd3c%alp+2.0_real_8
             CASE DEFAULT
                CALL stopgm(procedureN,'unknown version number of dft-d3',&
                   __LINE__,__FILE__)
             END SELECT
          ENDIF
       ELSE
          ! Broadcast empirical vdW stuff
          CALL mp_bcast_byte(vdwi,size_in_bytes_of(vdwi),parai%io_source,parai%cp_grp)
          CALL mp_bcast_byte(empvdwi,size_in_bytes_of(empvdwi),parai%io_source,parai%cp_grp)
          CALL mp_bcast_byte(empvdwr,size_in_bytes_of(empvdwr),parai%io_source,parai%cp_grp)
          IF (empvdwi%nvdw>0) THEN
             IF (.NOT.paral%io_parent) THEN
                ALLOCATE(idvdw(maxsys%ncorx),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                ALLOCATE(ivdw(maxsys%ncorx),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                ALLOCATE(jvdw(maxsys%ncorx),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                ALLOCATE(vdwst(maxsys%ncorx),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                ALLOCATE(vdwrm(maxsys%ncorx),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
                ALLOCATE(vdwbe(maxsys%ncorx),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ENDIF
             CALL mp_bcast(idvdw,SIZE(idvdw),parai%io_source,parai%cp_grp)
             CALL mp_bcast(ivdw,SIZE(ivdw),parai%io_source,parai%cp_grp)
             CALL mp_bcast(jvdw,SIZE(jvdw),parai%io_source,parai%cp_grp)
             CALL mp_bcast(vdwst,SIZE(vdwst),parai%io_source,parai%cp_grp)
             CALL mp_bcast(vdwrm,SIZE(vdwrm),parai%io_source,parai%cp_grp)
             CALL mp_bcast(vdwbe,SIZE(vdwbe),parai%io_source,parai%cp_grp)
          ENDIF
       ENDIF

    END SUBROUTINE empvdw_init
    ! ==--------------------------------------------------------------==
    SUBROUTINE vdwwf_init()

       ! Broadcast vdW-WF stuff
       CALL mp_bcast_byte(vdwi,size_in_bytes_of(vdwi),parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(vdwwfl,size_in_bytes_of(vdwwfl),parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(vdwwfi,size_in_bytes_of(vdwwfi),parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(vdwwfr,size_in_bytes_of(vdwwfr),parai%io_source,parai%cp_grp)
       CALL mp_bcast(wannl%twann,parai%io_source,parai%cp_grp)
       IF (vdwwfi%icriteri==0) THEN
          ! No additional variables to be broadcasted
       ELSEIF (vdwwfi%icriteri==1) THEN
          IF (.NOT.paral%io_parent) THEN
             ALLOCATE(ifragdata(vdwwfi%multifrag),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             ALLOCATE(radfrag(vdwwfi%multifrag),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
          ENDIF
          CALL mp_bcast(ifragdata,SIZE(ifragdata),parai%io_source,parai%cp_grp)
          CALL mp_bcast(radfrag,SIZE(radfrag),parai%io_source,parai%cp_grp)
       ELSEIF (vdwwfi%icriteri==2) THEN
          CALL mp_bcast(covrad,SIZE(covrad),parai%io_source,parai%cp_grp)
          IF (vdwwfi%nboadwf>0) THEN
             IF (.NOT.paral%io_parent) THEN
                ALLOCATE(boadwf(3,vdwwfi%nboadwf),STAT=ierr)
                IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                     __LINE__,__FILE__)
             ENDIF
             CALL mp_bcast(boadwf,SIZE(boadwf),parai%io_source,parai%cp_grp)
          ENDIF
       ENDIF
       ! Allocation of dynamical arrays is moved in vdw_wf_alloc routine

    END SUBROUTINE vdwwf_init
    ! ==--------------------------------------------------------------==
  END SUBROUTINE vdwin
  ! ==================================================================
  SUBROUTINE empvdwin(iunit)
    ! ==--------------------------------------------------------------==
    ! ==  THIS ROUTINE READS THE EMPIRICAL CORRECTION SUBSECTION      ==
    ! ==--------------------------------------------------------------==
    ! ==  THE INPUT HAS TO BE IN THE FOLLOWING FORMAT                 ==
    ! ==    EMPIRICAL CORRECTION                                      ==
    ! ==      Options                                                 ==
    ! ==    END EMPIRICAL CORRECTION                                  ==
    ! ==--------------------------------------------------------------==
    ! ==    For empirical vdW corrections either of keywords ALL or   ==
    ! ==    CUSTUM must be present in EMPIRICAL CORRECTION subsection ==
    ! ==                                                              ==
    ! ==  LIST OF OPTIONS                                             ==
    ! ==                                                              ==
    ! ==    ALL {DFT-D2,DFT-D3}                                       ==
    ! ==    CUSTOM                                                    ==
    ! ==      nvdw                                                    ==
    ! ==      vdwtyp ivdw jvdw vdwst vdwrm vdwbe                      ==
    ! ==      vdwtyp ivdw jvdw vdwst vdwrm vdwbe s6grim               ==
    ! ==      {C6,TANH} typ ivdw jvdw vdwst vdwrm vdwbe               ==
    ! ==      {DFT-D2} typ ivdw jvdw vdwst vdwrm vdwbe s6grim         ==
    ! ==      {DFT-D2} ivdw jvdw                                      ==
    ! ==    CUTOFF                                                    ==
    ! ==      vdweps                                                  ==
    ! ==    CELL                                                      ==
    ! ==      nxvdw nyvdw nzvdw                                       ==
    ! ==    S6GRIMME                                                  ==
    ! ==      s6_str                                                  ==
    ! ==    DFT-D3 PARAMETERS [THREEBODY,NUMGRAD]                     ==
    ! ==      cutoff cutoff_cn                                        ==
    ! ==    DFT-D3 FUNCTIONAL VERSION=version [TZ]                    ==
    ! ==      func                                                    ==
    ! ==    DFT-D3 NO_CONTRIBUTION                                    ==
    ! ==      species_1_w/o_contrib species_2_w/o_contrib ...         ==
    ! ==--------------------------------------------------------------==
    ! VDWC6 and VDWR0 according to Grimme not in a.u. !
    INTEGER                                  :: iunit

    CHARACTER(*), PARAMETER                  :: procedureN = 'empvdwin'
    CHARACTER(len=6), DIMENSION(4), PARAMETER :: &
      vdwtyp = (/'C6    ','TANH  ','DFT-D2','DFT-D3'/)
    REAL(real_8), DIMENSION(99), PARAMETER :: vdwc6 = (/0.14_real_8, &
      0.08_real_8, 1.61_real_8, 1.61_real_8,3.31_real_8, 1.75_real_8, &
      1.23_real_8, 0.70_real_8,0.75_real_8, 0.63_real_8, 5.71_real_8, &
      5.71_real_8, 10.79_real_8,9.23_real_8, 7.84_real_8, 5.57_real_8, &
      5.07_real_8,4.61_real_8, 10.80_real_8, 10.80_real_8, 10.80_real_8,&
      10.80_real_8, 10.80_real_8, 10.80_real_8, 10.80_real_8,10.80_real_8, &
      10.80_real_8, 10.80_real_8, 10.80_real_8,10.80_real_8, 16.99_real_8, &
      17.10_real_8, 16.37_real_8,12.64_real_8, 12.47_real_8, 12.01_real_8, &
      24.67_real_8,24.67_real_8, 24.67_real_8, 24.67_real_8, 24.67_real_8,&
      24.67_real_8, 24.67_real_8, 24.67_real_8, 24.67_real_8,24.67_real_8, &
      24.67_real_8, 24.67_real_8, 37.32_real_8,38.71_real_8, 38.44_real_8, &
      31.74_real_8, 31.50_real_8, 29.99_real_8,315.275_real_8, 226.994_real_8,&
      176.252_real_8,140.68_real_8, 140.68_real_8, 140.68_real_8, &
      140.68_real_8,140.68_real_8, 140.68_real_8, 140.68_real_8, 140.68_real_8&
      ,140.68_real_8, 140.68_real_8, 140.68_real_8, 140.68_real_8,&
      140.68_real_8, 140.68_real_8, 105.112_real_8, 81.24_real_8,81.24_real_8,&
      81.24_real_8, 81.24_real_8, 81.24_real_8, 81.24_real_8,81.24_real_8, &
      57.364_real_8, 57.254_real_8, 63.162_real_8, 63.540_real_8,55.283_real_8&
      , 57.171_real_8, 56.64_real_8, 0._real_8,0._real_8,0._real_8,0._real_8,&
      0._real_8,0._real_8,0._real_8,0._real_8,0._real_8,0._real_8,0._real_8,&
      0._real_8,0._real_8/)
    REAL(real_8), DIMENSION(99), PARAMETER :: vdwr0 = (/1.001_real_8, &
      1.012_real_8, 0.825_real_8, 1.408_real_8,1.485_real_8, 1.452_real_8, &
      1.397_real_8,1.342_real_8, 1.287_real_8, 1.243_real_8, 1.144_real_8,&
      1.364_real_8, 1.639_real_8, 1.716_real_8, 1.705_real_8,1.683_real_8, &
      1.639_real_8, 1.595_real_8, 1.485_real_8,1.474_real_8, 1.562_real_8, &
      1.562_real_8, 1.562_real_8,1.562_real_8, 1.562_real_8, 1.562_real_8, &
      1.562_real_8,1.562_real_8, 1.562_real_8, 1.562_real_8,1.650_real_8, &
      1.727_real_8, 1.760_real_8, 1.771_real_8,1.749_real_8, 1.727_real_8, &
      1.628_real_8, 1.606_real_8,1.639_real_8, 1.639_real_8, 1.639_real_8, &
      1.639_real_8,1.639_real_8, 1.639_real_8, 1.639_real_8, 1.639_real_8,&
      1.639_real_8, 1.639_real_8, 1.672_real_8, 1.804_real_8,1.881_real_8, &
      1.892_real_8, 1.892_real_8, 1.881_real_8,1.802_real_8, 1.762_real_8, &
      1.720_real_8, 1.753_real_8,1.753_real_8, 1.753_real_8, 1.753_real_8, &
      1.753_real_8,1.753_real_8, 1.753_real_8, 1.753_real_8, 1.753_real_8,&
      1.753_real_8, 1.753_real_8, 1.753_real_8, 1.753_real_8,1.753_real_8, &
      1.788_real_8, 1.772_real_8, 1.772_real_8,1.772_real_8, 1.772_real_8, &
      1.772_real_8, 1.772_real_8,1.772_real_8, 1.758_real_8, 1.986_real_8, &
      1.944_real_8,1.898_real_8, 2.005_real_8, 1.991_real_8, 1.924_real_8,&
      1.0e10_real_8,1.0e10_real_8,1.0e10_real_8,1.0e10_real_8,1.0e10_real_8,&
      1.0e10_real_8,1.0e10_real_8,1.0e10_real_8,1.0e10_real_8,1.0e10_real_8,&
      1.0e10_real_8,1.0e10_real_8,1.0e10_real_8/)
    REAL(real_8), PARAMETER                  :: convfac = 1.733979e1_real_8

    CHARACTER(len=80)                        :: line, previous_line, error_message,&
                                                index_list
    INTEGER                                  :: i, i1, ierr, iout, j, k, first, last
    LOGICAL                                  :: erread, tdftd2, telst,&
                                                go_on_reading, something_went_wrong
    REAL(real_8)                             :: a6grim

    ! Safety check
    !
    IF (.not. paral%io_parent) CALL stopgm(procedureN,'Must not be called by '//&
                                           'anything but the io_parent.', &
                                           __LINE__,__FILE__)
    !
    ! Defaults
    !
    vdwi%nxvdw=0
    vdwi%nyvdw=0
    vdwi%nzvdw=0
    tdftd2=.FALSE.
    tdftd3=.FALSE.
    telst=.FALSE.
    empvdwi%nvdw=0
    empvdwr%vdweps=1.0e-2_real_8
    a6grim=20.0_real_8
    empvdwr%s6grim=0.0_real_8   ! No default for s6 in Grimme's vdW
    dftd3i%threebody=.FALSE.
    dftd3i%numgrad=.FALSE.
    dftd3i%cutoff=SQRT(9000.0_real_8)
    dftd3i%cutoff_cn=SQRT(1600.0_real_8)
    dftd3f%func=' '
    dftd3f%version=4
    dftd3f%tz=.FALSE.
    !
    ! Allocate dynamical arrays
    ALLOCATE(idvdw(maxsys%ncorx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ivdw(maxsys%ncorx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(jvdw(maxsys%ncorx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vdwst(maxsys%ncorx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vdwrm(maxsys%ncorx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vdwbe(maxsys%ncorx),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    something_went_wrong = .FALSE.
    line          = ' ' 
    previous_line = ' ' 
    error_message = ' '
    go_on_reading = .TRUE.
    !
    DO WHILE(go_on_reading)
       previous_line = line
       READ(iunit,'(A80)',iostat=ierr) line
       IF (ierr /= 0) THEN
          error_message        = 'ERROR WHILE READING LINE'
          something_went_wrong = .TRUE.
          go_on_reading        = .FALSE.
       ELSEIF (keyword_contains(line,'END')) THEN
          go_on_reading        = .FALSE.
       ! Use hard-coded parameters for all
       ELSEIF (keyword_contains(line,'ALL')) THEN
          IF (keyword_contains(line,'DFT-D2')) THEN
             tdftd2=.TRUE.
             empvdwi%nvdw=ions1%nsp*(ions1%nsp+1)/2
             ! Simulate full input for all types.
             k=0
             DO i=1,ions1%nsp
                ! Reality Check. Parametrization exists only for "H" up to "RN"
                IF ((ions0%iatyp(i) > 86).AND.(ions0%iatyp(i) <= 99)) THEN
                   WRITE(output_unit,'(A,11X,A)') " WARNING : "//&
                        "DFT-D2 VDW ONLY SUPPORTS ELEMENTS FROM H TO RN",&
                        "ELEMENTS BEYOND RN GIVE A ZERO CONTRIBUTION !!"
                ELSE IF (ions0%iatyp(i) > 99) THEN
                   CALL stopgm('VDWIN','DFT-D2 VDW ONLY SUPPORTS'//&
                        'ELEMENTS UP TO A.N. 99 !!',& 
                        __LINE__,__FILE__)
                ENDIF
                DO j=i,ions1%nsp
                   k=k+1
                   ivdw(k)=i
                   jvdw(k)=j
                ENDDO
             ENDDO
             DO i=1,empvdwi%nvdw
                vdwst(i)=SQRT(vdwc6(ions0%iatyp(ivdw(i)))&
                     *vdwc6(ions0%iatyp(jvdw(i))))*convfac
                vdwrm(i)=(vdwr0(ions0%iatyp(ivdw(i)))&
                     +vdwr0(ions0%iatyp(jvdw(i))))*fbohr
                vdwbe(i)=a6grim
                idvdw(i)=3
             ENDDO
          ELSE IF(keyword_contains(line,'DFT-D3')) THEN
             tdftd3=.TRUE.
          ELSE
             CALL stopgm(procedureN,&
                  '"ALL" OPTION ONLY AVAILABLE FOR TYPE "DFT-D2" OR "DFT-D3"',& 
                  __LINE__,__FILE__)
          ENDIF
       ELSEIF (keyword_contains(line,'CUSTOM')) THEN
          ! Read individual pairs
          !bug   i1=1
          !bug   CALL readsi(line,i1,iout,empvdwi%nvdw,erread)
          READ(iunit,*,iostat=ierr) empvdwi%nvdw
          IF (ierr /= 0) THEN
             error_message        = 'ERROR WHILE READING LINE'
             something_went_wrong = .TRUE.
             go_on_reading        = .FALSE.
          ENDIF
          DO i=1,empvdwi%nvdw
             READ(iunit,'(A)',iostat=ierr) line
             IF (ierr /= 0) THEN
                error_message        = 'ERROR WHILE READING LINE'
                something_went_wrong = .TRUE.
                go_on_reading        = .FALSE.
             ENDIF
             IF (INDEX(line,'C6') /= 0) THEN
                telst=.TRUE.
                idvdw(i)=1
                i1=INDEX(line,'C6')+2
             ELSE IF (INDEX(line,'TANH') /= 0) THEN
                telst=.TRUE.
                idvdw(i)=2
                i1=INDEX(line,'TANH')+4
             ELSE IF (INDEX(line,'DFT-D2') /= 0) THEN
                tdftd2=.TRUE.
                idvdw(i)=3
                i1=INDEX(line,'DFT-D2')+6
             ELSE
                WRITE(output_unit,*) ' ERROR WHILE READING VDW CORRECTIONS'
                CALL stopgm(procedureN,'UNKNOWN VDW TYPE',& 
                     __LINE__,__FILE__)
             ENDIF
             CALL readsi(line,i1,iout,ivdw(i),erread)
             i1=iout
             CALL readsi(line,i1,iout,jvdw(i),erread)
             i1=iout
             ! Grimme parameters are hardcoded, but we allow overrides.
             IF (idvdw(i) == 3) THEN
                CALL readsr(line,i1,iout,vdwst(i),erread)
                IF (erread) THEN! ok use defaults
                   IF (MAX(ions0%iatyp(ivdw(i)),ions0%iatyp(jvdw(i))) > 86) THEN
                      CALL stopgm(procedureN,'DFT-D2 VDW ONLY SUPPORTS'//&
                           'ELEMENTS "H" TO "RN"',& 
                           __LINE__,__FILE__)
                   ENDIF
                   vdwst(i)=SQRT(vdwc6(ions0%iatyp(ivdw(i)))&
                        *vdwc6(ions0%iatyp(jvdw(i))))*convfac    ! bug-fix
                   vdwrm(i)=(vdwr0(ions0%iatyp(ivdw(i)))&
                        +vdwr0(ions0%iatyp(jvdw(i))))*fbohr      ! bug-fix
                   vdwbe(i)=a6grim
                ELSE
                   i1=iout
                   CALL readsr(line,i1,iout,vdwrm(i),erread)
                   i1=iout
                   CALL readsr(line,i1,iout,vdwbe(i),erread)
                ENDIF
             ELSE
                CALL readsr(line,i1,iout,vdwst(i),erread)
                i1=iout
                CALL readsr(line,i1,iout,vdwrm(i),erread)
                i1=iout
                CALL readsr(line,i1,iout,vdwbe(i),erread)
             ENDIF
          ENDDO
       ELSEIF (keyword_contains(line,'CUTOFF',alias='VDW-CUTOFF')) THEN
          READ(iunit,*,iostat=ierr) empvdwr%vdweps
          IF (ierr /= 0) THEN
             error_message        = 'ERROR WHILE READING LINE'
             something_went_wrong = .TRUE.
             go_on_reading        = .FALSE.
          ENDIF
       ELSEIF (keyword_contains(line,'CELL',alias='VDW-CELL')) THEN
          READ(iunit,*,iostat=ierr) vdwi%nxvdw,vdwi%nyvdw,vdwi%nzvdw
          IF (ierr /= 0) THEN
             error_message        = 'ERROR WHILE READING LINE'
             something_went_wrong = .TRUE.
             go_on_reading        = .FALSE.
          ENDIF
          ! mb-s6 >>> s6 as input
       ELSEIF (keyword_contains(line,'S6GRIMME')) THEN
          READ(iunit,*,iostat=ierr) empvdwc%s6_str
          IF (ierr /= 0) THEN
             error_message        = 'ERROR WHILE READING LINE'
             something_went_wrong = .TRUE.
             go_on_reading        = .FALSE.
          ENDIF
          empvdwr%s6grim=set_s6(empvdwc%s6_str)
          ! SI: S6_STR should not in general contradict with the functional used for DFT,
          ! SI: but one should be able to do so
          IF (empvdwc%dft_func /= empvdwc%s6_str.AND.empvdwc%dft_func /= " ") THEN
             WRITE(output_unit,*)&
                  '!!! WARNING !!! GRIMME CORRECTION USES', empvdwc%s6_str
             WRITE(output_unit,*)&
                  '!!! WARNING !!! WHEREAS ', empvdwc%dft_func, 'IS USED'
          ENDIF
       ELSEIF (keyword_contains(line,'DFT-D3',and='PARAMETERS')) THEN
          IF (keyword_contains(line,'THREEBODY')) dftd3i%threebody=.TRUE.
          IF (keyword_contains(line,'NUMGRAD')) dftd3i%numgrad=.TRUE.
          READ(iunit,*,iostat=ierr) dftd3i%cutoff,dftd3i%cutoff_cn
          IF (ierr /= 0) THEN
             error_message        = 'ERROR WHILE READING LINE'
             something_went_wrong = .TRUE.
             go_on_reading        = .FALSE.
          ENDIF
          IF (.NOT.cntl%bohr) THEN
             dftd3i%cutoff=dftd3i%cutoff*fbohr
             dftd3i%cutoff_cn=dftd3i%cutoff_cn*fbohr
          ENDIF
       ELSEIF (keyword_contains(line,'DFT-D3',and='FUNCTIONAL')) THEN
          IF (keyword_contains(line,'TZ')) dftd3f%tz=.TRUE.
          IF (keyword_contains(line,'VERSION',cut_at='=')) THEN
             first = index_of_delimiter(line,'VERSION','=')
             CALL readsi(line,first,last,dftd3f%version,erread)
             IF (erread) THEN
                error_message        = "ERROR WHILE READING VALUE"
                something_went_wrong = .true.
                go_on_reading        = .false.
             ENDIF
          ENDIF
          READ(iunit,'(A80)',iostat=ierr) line
          IF (ierr /= 0) THEN
             error_message        = 'ERROR WHILE READING LINE'
             something_went_wrong = .TRUE.
             go_on_reading        = .FALSE.
          ENDIF
          CALL xstring(line,first,last)
          dftd3f%func=line(first:last)
       ELSEIF (keyword_contains(line,'DFT-D3',and='NO_CONTRIBUTION')) THEN
          READ(iunit,'(A80)',iostat=ierr) line
          IF (ierr /= 0) THEN
             error_message        = 'ERROR WHILE READING NO_CONTRIBUTION LIST'
             something_went_wrong = .TRUE.
             go_on_reading        = .FALSE.
          ENDIF
          first=1
          contrib_loop: DO
             erread=.FALSE.
             CALL readsi(line,first,last,i1,erread)
             IF (erread) EXIT contrib_loop
             tdftd3noc(i1) = .TRUE.
             first = last
          ENDDO contrib_loop
       ENDIF
    ENDDO
    !
    IF (something_went_wrong) THEN
       WRITE(output_unit,'(/,1X,64("!"))')
       WRITE(output_unit,'(1X, A, 1X, A)') 'ERROR:', 'PROBLEM WHILE READING EMPIRICAL CORRECTION SECTION:' 
       WRITE(output_unit,'(8X, A)') TRIM(ADJUSTL(error_message))
       IF (line /= ' ' .or. previous_line /= ' ') THEN
          WRITE(output_unit,'(8X, A)') 'THE LAST TWO LINES READ WITHIN THE SECTION WERE:'
          WRITE(output_unit,'(/,1X, A)') TRIM(ADJUSTL(previous_line))
          WRITE(output_unit,'(1X, A)') TRIM(ADJUSTL(line))
       END IF
       WRITE(output_unit,'(1X,64("!"))')
       CALL stopgm(procedureN,'Error while reading EMPIRICAL CORRECTION section, cf. output file',& 
            __LINE__,__FILE__)
    ENDIF
    !
    CALL empvdw_report()

    CONTAINS

    ! ==--------------------------------------------------------------==
    SUBROUTINE empvdw_report()

      WRITE(output_unit,*)
      WRITE(output_unit,'(A)')&
           ' EMPIRICAL VAN DER WAALS CORRECTION AFTER:'
      IF (tdftd3) THEN
         ! Set functional name if not given
         IF (dftd3f%func==' ') dftd3f%func=dftd3func(empvdwc%dft_func)
         WRITE(output_unit,'(A)')&
              '    S. GRIMME ET AL., J. CHEM. PHYS. 132, 154104 (2010)'
         WRITE(output_unit,'(A)')&
              '    INPUT VARIABLES USED FOR SCALING AND DAMPING PARAMETERS:'
         WRITE(output_unit,'(A,1X,A,1X,A,I3,1X,A,L3)')&
              '       FUNCTIONAL=',TRIM(ADJUSTL(dftd3f%func)),' VERSION=',dftd3f%version,' TZ=',dftd3f%tz
         WRITE(output_unit,'(A,T46,2F10.3)')&
              '    DISTANCE CUTOFF FOR VDW AND CN:',dftd3i%cutoff,dftd3i%cutoff_cn
         IF (dftd3i%threebody) WRITE(output_unit,'(A)')&
              '    INCLUDE THREE BODY TERMS' 
         IF (dftd3i%numgrad) WRITE(output_unit,'(A)')&
              '    USE NUMERICAL GRADIENTS'
         IF (ANY(tdftd3noc(:))) THEN
            index_list = list_of_indices()
            WRITE(output_unit,'(A,1X,A)')&
              '    NO CONTRIBUTION FROM SPECIES',TRIM(ADJUSTL(index_list))
         ENDIF
      ENDIF
      IF (telst) THEN
         WRITE(output_unit,'(A)')&
              '    R. LE SAR, J. PHYS. CHEM. 88, 4272 (1984)'
         WRITE(output_unit,'(A)')&
              '    M. ELSTNER ET AL., J. CHEM. PHYS. 114, 5149 (2001)'
      ENDIF
      IF (tdftd2) THEN
         WRITE(output_unit,'(A)')&
              '    S. GRIMME, J. COMP. CHEM. 27, 1787 (2006)'
         ! mb-s6 >>> s6 as input
         ! SI: If 0 then try to set it according to the functional in DFT section
         IF (empvdwr%s6grim == 0._real_8) THEN
            empvdwr%s6grim=set_s6(empvdwc%dft_func)
            empvdwc%s6_str=empvdwc%dft_func
         ENDIF
         IF (empvdwr%s6grim == 0._real_8) THEN
            CALL stopgm('VDWIN',&
                 'ONE MUST PROVIDE A FUNCTIONAL FOR GRIMME VDW',& 
                 __LINE__,__FILE__)
         ELSE
            WRITE(output_unit,'(A,1F5.2,A,T47,A,T55,A)')&
                 '    GLOBAL SCALING FACTOR S6 =', empvdwr%s6grim,':',&
                 ADJUSTR(empvdwc%s6_str),' FUNCTIONAL'
         ENDIF
         ! mb-s6 <<< s6 as input
      ENDIF
      IF (.NOT.tdftd3) THEN
         WRITE(output_unit,'(A/,T8,A,T21,A,T30,A,T44,A,T62,A)')&
              '    PARAMETERS FOR PAIRS:','ELEMENTS','C6(A.U.)',&
              'RADIUS(A.U.)','ALPHA(A.U.)','TYPE'
         DO i=1,empvdwi%nvdw
            WRITE(output_unit,'(T9,A,T13,A,F14.7,2(1X,F12.7),T60,A)')&
                 elem%el(ions0%iatyp(ivdw(i))),elem%el(ions0%iatyp(jvdw(i))),vdwst(i),vdwrm(i),&
                 vdwbe(i),ADJUSTR(vdwtyp(idvdw(i)))
         ENDDO
      ENDIF
      WRITE(output_unit,'(A,A,T52,I2,A,I2,A,I2,A)') '    VAN DER WAALS SUM',&
           ' IN REAL SPACE OVER ',vdwi%nxvdw+1,'*',vdwi%nyvdw+1,'*',vdwi%nzvdw+1,' CELLS'
      IF(.NOT.tdftd3) WRITE(output_unit,*)

    END SUBROUTINE empvdw_report
    ! ==--------------------------------------------------------------==
    ELEMENTAL FUNCTION list_of_indices() RESULT (indices)

       CHARACTER(len=80)         :: indices
       !
       ! (I3) + (,) + (space) = space for 16 indices
       !
       CHARACTER(len=3)          :: atom_index
       LOGICAL                   :: first_in_list
       INTEGER                   :: isp

       indices       = ' '
       first_in_list = .true.

       DO isp=1,ions1%nsp
          atom_index = ' '
          IF (tdftd3noc(isp)) THEN
             WRITE(atom_index,'(I3)') isp
             IF (first_in_list) THEN
                indices       = atom_index
                first_in_list = .false.
             ELSE
                indices = TRIM(ADJUSTL(indices))//&
                          ', '//TRIM(ADJUSTL(atom_index))
             ENDIF
          ENDIF
       ENDDO

    END FUNCTION list_of_indices
    ! ==--------------------------------------------------------------==
  END SUBROUTINE empvdwin
  ! ==================================================================
  SUBROUTINE wannvdwin(iunit)
    ! ==--------------------------------------------------------------==
    ! ==  THIS ROUTINE READS THE WANNIER CORRECTION SUBSECTION        ==
    ! ==--------------------------------------------------------------==
    ! ==  THE INPUT HAS TO BE IN THE FOLLOWING FORMAT                 ==
    ! ==    WANNIER CORRECTION                                        ==
    ! ==       Options                                                ==
    ! ==    END WANNIER CORRECTION                                    ==
    ! ==--------------------------------------------------------------==
    ! ==    Currently, three criteria for fragmenting the system are  ==
    ! ==    implemented:                                              ==
    ! ==    (1) the system is subdivided into 2 fragments separated   ==
    ! ==        by 'zlevel', namely a z coordinate that separates     ==
    ! ==        the fragments                                         ==
    ! ==          FRAGMENT ZLEVEL                                     ==
    ! ==            zlevel                                            ==
    ! ==    (2) the system is subdivided into N fragments selected    ==
    ! ==        by choosing a given ion and a given radius            ==
    ! ==          FRAGMENT RADIUS                                     ==
    ! ==            multifrag                                         ==
    ! ==            i radius                                          ==
    ! ==            ....                                              ==
    ! ==    (3) the system is subdivided into fragments automatically ==
    ! ==        detected by using predefined covalent bond radii      ==
    ! ==          FRAGMENT BOND (default)                             ==
    ! ==            xmfacwf                                           ==
    ! ==                                                              ==
    ! ==  LIST OF OPTIONS                                             ==
    ! ==                                                              ==
    ! ==    VERSION                                                   ==
    ! ==      iswitchvdw                                              ==
    ! ==    FRAGMENT ZLEVEL                                           ==
    ! ==      zlevel                                                  ==
    ! ==    FRAGMENT RADIUS                                           ==
    ! ==      multifrag                                               ==
    ! ==      i  radius                                               ==
    ! ==      ....                                                    ==
    ! ==    FRAGMENT BOND [COVRAD]                                    ==
    ! ==      xmfacwf                                                 ==
    ! ==      covrad2 ...                                             ==
    ! ==    DAMPING [DIPOLE]                                          ==
    ! ==      a6                                                      ==
    ! ==    RESTART WANNIER                                           ==
    ! ==    ENERGY MONOMER                                            ==
    ! ==      enmonomer                                               ==
    ! ==    TOLERANCE WANNIER                                         ==
    ! ==      tolwann                                                 ==
    ! ==    TOLERANCE REFERENCE                                       ==
    ! ==      tolref                                                  ==
    ! ==    CHANGE BONDS                                              ==
    ! ==      nboadwf                                                 ==
    ! ==      i  j   +/- 1                                            ==
    ! ==      ....                                                    ==
    ! ==    CELL                                                      ==
    ! ==      nxvdw nyvdw nzvdw                                       ==
    ! ==    PRINT [INFO,FRAGMENT,C6,FORCES]                           ==
    ! ==    REFERENCE [FITTED CENTERS]                                ==
    ! ==                                                              ==
    ! ==  Note: the unit of zlevel, radfrag, covrad2, and tollwan     ==
    ! ==        follows a logical variable cntl%bohr                  ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: iunit

    CHARACTER(*), PARAMETER                  :: procedureN = 'wannvdwin'

    CHARACTER(len=1)                         :: csign
    CHARACTER(len=80)                        :: line, previous_line, &
                                                error_message
    INTEGER                                  :: i, ierr, l
    LOGICAL                                  :: resetcb, tflag, &
                                                go_on_reading,  &
                                                something_went_wrong
    REAL(real_8), ALLOCATABLE                :: covrad2(:)

    ! Safety check
    !
    IF (.not. paral%io_parent) CALL stopgm(procedureN,'Must not be called by '//&
                                           'anything but the io_parent.', &
                                           __LINE__,__FILE__)
    !
    ! Defaults
    !
    vdwi%nxvdw=0
    vdwi%nyvdw=0
    vdwi%nzvdw=0
    IF (cntl%tlsd) THEN
       vdwwfi%nelpwf=1
    ELSE
       vdwwfi%nelpwf=2
    ENDIF
    vdwwfi%iswitchvdw=1     ! vdW-WF
    vdwwfi%icriteri=2       ! Bond
    vdwwfr%xmfacwf=1.35_real_8
    vdwwfr%a6=20.0_real_8
    vdwwfl%trwannc=.FALSE.
    vdwwfr%enmonomer=0.0_real_8
    vdwwfl%twannup=.TRUE.
    vdwwfr%tolwann=0.0_real_8
    vdwwfr%tolref=1.0_real_8
    vdwwfi%nboadwf=0
    vdwwfi%multifrag=SUM(ions0%na(1:ions1%nsp))
    vdwwfl%tpinfo=.TRUE.
    vdwwfl%tpfrag=.FALSE.
    vdwwfl%tpc6=.FALSE.
    vdwwfl%tpforce=.FALSE.
    vdwwfl%treffit=.FALSE.
    vdwwfl%tdampda=.FALSE.
    !
    !
    !
    something_went_wrong = .FALSE.
    line          = ' ' 
    previous_line = ' ' 
    error_message = ' '
    go_on_reading = .TRUE.
    !
    DO WHILE(go_on_reading)
       previous_line = line
       READ(iunit,'(A80)',iostat=ierr) line
       IF (ierr /= 0) THEN
          error_message        = 'ERROR WHILE READING LINE'
          something_went_wrong = .TRUE.
          go_on_reading        = .FALSE.
       ELSEIF (keyword_contains(line,'END')) THEN
          go_on_reading        = .FALSE.
       ELSEIF (keyword_contains(line,'VERSION')) THEN
          READ(iunit,*,iostat=ierr) vdwwfi%iswitchvdw
          IF (ierr /= 0) THEN
             error_message        = 'ERROR WHILE READING LINE'
             something_went_wrong = .TRUE.
             go_on_reading        = .FALSE.
          ENDIF
          IF ((vdwwfi%iswitchvdw < 1).OR.(vdwwfi%iswitchvdw > 2)) THEN
             CALL stopgm(procedureN,'WRONG ISWITCHVDW! ',& 
                  __LINE__,__FILE__)
          ENDIF
       ELSEIF(keyword_contains(line,'FRAGMENT',alias='FRAGMENTS') .AND. &
            keyword_contains(line,'ZLEVEL')) THEN
          vdwwfi%icriteri=0
          vdwwfi%multifrag=2
          READ(iunit,*,iostat=ierr) vdwwfr%zlevel
          IF (ierr /= 0) THEN
             error_message        = 'ERROR WHILE READING LINE'
             something_went_wrong = .TRUE.
             go_on_reading        = .FALSE.
          ENDIF
          IF(.NOT.cntl%bohr) vdwwfr%zlevel=vdwwfr%zlevel*fbohr
       ELSEIF(keyword_contains(line,'FRAGMENT',and='RADIUS')) THEN
          vdwwfi%icriteri=1
          READ(iunit,*,iostat=ierr) vdwwfi%multifrag
          IF (ierr /= 0) THEN
             error_message        = 'ERROR WHILE READING LINE'
             something_went_wrong = .TRUE.
             go_on_reading        = .FALSE.
          ENDIF
          ALLOCATE(ifragdata(vdwwfi%multifrag),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(radfrag(vdwwfi%multifrag),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          DO i=1,vdwwfi%multifrag
             READ(iunit,*,iostat=ierr) ifragdata(i),radfrag(i)
             IF (ierr /= 0) THEN
                error_message        = 'ERROR WHILE READING LINE'
                something_went_wrong = .TRUE.
                go_on_reading        = .FALSE.
             ENDIF
          ENDDO
          IF(.NOT.cntl%bohr) CALL dscal(vdwwfi%multifrag,fbohr,radfrag,1)
       ELSEIF(keyword_contains(line,'FRAGMENT',and='BOND') .OR. &
              keyword_contains(line,'FRAGMENT',and='BONDS')) THEN
          vdwwfi%icriteri=2
          resetcb=.FALSE.
          READ(iunit,*,iostat=ierr) vdwwfr%xmfacwf
          IF (ierr /= 0) THEN
             error_message        = 'ERROR WHILE READING LINE'
             something_went_wrong = .TRUE.
             go_on_reading        = .FALSE.
          ENDIF
          IF(keyword_contains(line,'COVRAD')) THEN
             ALLOCATE(covrad2(ions1%nsp),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             READ(iunit,*,iostat=ierr) (covrad2(i),i=1,ions1%nsp)
             IF (ierr /= 0) THEN
                error_message        = 'ERROR WHILE READING LINE'
                something_went_wrong = .TRUE.
                go_on_reading        = .FALSE.
             ENDIF
             IF(cntl%bohr) CALL dscal(ions1%nsp,1.0_real_8/fbohr,covrad2,1)
             resetcb=.TRUE.
          ENDIF
       ELSEIF (keyword_contains(line,'DAMPING')) THEN
          IF (keyword_contains(line,'DIPOLE')) THEN
             vdwwfl%tdampda=.TRUE.
          ELSE
             READ(iunit,*,iostat=ierr) vdwwfr%a6
             IF (ierr /= 0) THEN
                error_message        = 'ERROR WHILE READING LINE'
                something_went_wrong = .TRUE.
                go_on_reading        = .FALSE.
             ENDIF
             IF (vdwwfr%a6 < 1.0e-05_real_8) THEN
                WRITE(output_unit,'(/,1X,A,1X,A,/)')&
                     'DAMPING FACTOR TOO SMALL. SET TO 20.0',&
                     'ACCORDING TO: MOL. PHYS. 107, 999 (2009)'
                vdwwfr%a6=20.0_real_8
             ENDIF
          ENDIF
       ELSEIF(keyword_contains(line,'RESTART') .AND.&
            keyword_contains(line,'WANN',alias='WANNIER')) THEN
          vdwwfl%trwannc=.TRUE.
       ELSEIF(keyword_contains(line,'ENERGY') .AND.&
            keyword_contains(line,'MONOMER')) THEN
          READ(iunit,*,iostat=ierr) vdwwfr%enmonomer
          IF (ierr /= 0) THEN
             error_message        = 'ERROR WHILE READING LINE'
             something_went_wrong = .TRUE.
             go_on_reading        = .FALSE.
          ENDIF
       ELSEIF(keyword_contains(line,'TOLERANCE',and='WANNIER')) THEN
          READ(iunit,*,iostat=ierr) vdwwfr%tolwann
          IF (ierr /= 0) THEN
             error_message        = 'ERROR WHILE READING LINE'
             something_went_wrong = .TRUE.
             go_on_reading        = .FALSE.
          ENDIF
          IF(.NOT.cntl%bohr) vdwwfr%tolwann=vdwwfr%tolwann*fbohr
       ELSEIF(keyword_contains(line,'TOLERANCE',and='REFERENCE')) THEN
          READ(iunit,*,iostat=ierr) vdwwfr%tolref
          IF (ierr /= 0) THEN
             error_message        = 'ERROR WHILE READING LINE'
             something_went_wrong = .TRUE.
             go_on_reading        = .FALSE.
          ENDIF
       ELSEIF(keyword_contains(line,'CHANGE',and='BONDS')) THEN
          READ(iunit,*,iostat=ierr) vdwwfi%nboadwf
          IF (ierr /= 0) THEN
             error_message        = 'ERROR WHILE READING LINE'
             something_went_wrong = .TRUE.
             go_on_reading        = .FALSE.
          ENDIF
          IF (vdwwfi%nboadwf > 0) THEN
             ALLOCATE(boadwf(3,vdwwfi%nboadwf),STAT=ierr)
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
                  __LINE__,__FILE__)
             DO i=1,vdwwfi%nboadwf
                READ(iunit,*,iostat=ierr) (boadwf(l,i),l=1,3)
                IF (ierr /= 0) THEN
                   error_message        = 'ERROR WHILE READING LINE'
                   something_went_wrong = .TRUE.
                   go_on_reading        = .FALSE.
                ENDIF
             ENDDO
          ENDIF
       ELSEIF (keyword_contains(line,'CELL')) THEN
          READ(iunit,*,iostat=ierr) vdwi%nxvdw,vdwi%nyvdw,vdwi%nzvdw
          IF (ierr /= 0) THEN
             error_message        = 'ERROR WHILE READING LINE'
             something_went_wrong = .TRUE.
             go_on_reading        = .FALSE.
          ENDIF
       ELSEIF (keyword_contains(line,'PRINT')) THEN
          IF (keyword_contains(line,'OFF')) THEN
             tflag=.FALSE.
          ELSE
             tflag=.TRUE.
          ENDIF
          IF (keyword_contains(line,'INFO')) vdwwfl%tpinfo=tflaG
          IF (keyword_contains(line,'FRAGMENT',alias='FRAGMENTS')) vdwwfl%tpfrag=tflaG
          IF (keyword_contains(line,'C6')) vdwwfl%tpc6=tflag
          IF (keyword_contains(line,'FORCES',alias='FORCE')) vdwwfl%tpforce=tflaG
       ELSEIF (keyword_contains(line,'REFERENCE')) THEN
          IF (keyword_contains(line,'FIT')) vdwwfl%treffit=.TRUE.
          IF (keyword_contains(line,'CENTER')) vdwwfl%treffit=.FALSE.
       ENDIF
    ENDDO
    !
    IF (something_went_wrong) THEN
       WRITE(output_unit,'(/,1X,64("!"))')
       WRITE(output_unit,'(1X, A, 1X, A)') 'ERROR:', 'PROBLEM WHILE READING WANNIER CORRECTION SECTION:' 
       WRITE(output_unit,'(8X, A)') TRIM(ADJUSTL(error_message))
       IF (line /= ' ' .or. previous_line /= ' ') THEN
          WRITE(output_unit,'(8X, A)') 'THE LAST TWO LINES READ WITHIN THE SECTION WERE:'
          WRITE(output_unit,'(/,1X, A)') TRIM(ADJUSTL(previous_line))
          WRITE(output_unit,'(1X, A)') TRIM(ADJUSTL(line))
       END IF
       WRITE(output_unit,'(1X,64("!"))')
       CALL stopgm(procedureN,'Error while reading WANNIER CORRECTION section, cf. output file',& 
            __LINE__,__FILE__)
    ENDIF
    !
    !
!!!cmb-This check can be eliminated. HUGE(0) allocation does already the
!!!cmb-job and the suitable size is set later in vdw_wf_alloc. 
!!!cmb-Having this here can cause the code hanging/crashing for large systems.
!!!    IF (vdwwfi%multifrag > nfragx)&
!!!         CALL stopgm(procedureN,'NFRAGX TOO SMALL ! ',& 
!!!         __LINE__,__FILE__)
    IF (.NOT.wannl%twann) wannl%twann=.TRUE.
    !
    CALL wannvdw_report()

    CONTAINS

    ! ==--------------------------------------------------------------==
    SUBROUTINE wannvdw_report()
    
      WRITE(output_unit,*)
      WRITE(output_unit,'(A)')&
           ' VAN DER WAALS CORRECTION USING WANNIER CENTERS AFTER:'
      IF (vdwwfi%iswitchvdw == 1) THEN
         WRITE(output_unit,'(A)')&
              '    P.L. SILVESTRELLI, PHYS. REV. LETT. 100, 053002 (2008)'
         WRITE(output_unit,'(A)')&
              '                         J. PHYS. CHEM. A 113, 5224 (2009)'
      ENDIF
      IF (vdwwfi%iswitchvdw == 2) THEN
         WRITE(output_unit,'(A)')&
              '    A. AMBROSETTI & P.L. SILVESTRELLI, PRB 85, 073101 (2012)'
      ENDIF
      IF (vdwwfi%icriteri == 0) THEN
         WRITE(output_unit,'(A)')&
              '    ADOPT Z-LEVEL CRITERION FOR FRAGMENTATION'
         WRITE(output_unit,'(A,T56,F10.3)')&
              '       ZLEVEL SEPARATING TWO FRAGMENTS(A.U.): ',vdwwfr%zlevel
      ENDIF
      IF (vdwwfi%icriteri == 1) THEN
         WRITE(output_unit,'(A)')&
              '    ADOPT RADIUS CRITERION FOR FRAGMENTATION'
         WRITE(output_unit,'(A,T61,I5)')&
              '       NUMBER OF FRAGMENTS:',vdwwfi%multifrag
         WRITE(output_unit,'(A,1X,A,3X,A)')&
              '       IFRAG','ATOM','RADIUS(A.U.)'
         DO i=1,vdwwfi%multifrag
            WRITE(output_unit,'(7X,2I5,3X,F12.5)') i,ifragdata(i),radfrag(i)
         ENDDO
      ENDIF
      IF (vdwwfi%icriteri == 2) THEN
         WRITE(output_unit,'(A)')&
              '    ADOPT BOND CRITERION FOR FRAGMENTATION'
         WRITE(output_unit,'(A,T56,F10.3)')&
              '       TOLERANCE RATIO FOR COVALENT BONDS:',vdwwfr%xmfacwf
         IF (resetcb) THEN
            DO i=1,ions1%nsp
               IF (covrad2(i)>0.0_real_8) THEN
                  covrad(ions0%iatyp(i))=covrad2(i)
               ENDIF
            ENDDO
            DEALLOCATE(covrad2,STAT=ierr)
            IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
                 __LINE__,__FILE__)
         ENDIF
         WRITE(output_unit,'(T11,A,T24,A)') 'ELEMENT','RADIUS(ANGSTROM)'
         DO i=1,ions1%nsp
            WRITE(output_unit,'(T13,A,T24,F11.4)') elem%el(ions0%iatyp(i)),covrad(ions0%iatyp(i))
         ENDDO
         IF (vdwwfi%nboadwf > 0) THEN
            WRITE(output_unit,'(A,T56,I10)')&
                 '     NUMBER OF CHANGING BOND ASSIGNMENT:',vdwwfi%nboadwf
            WRITE(output_unit,'(A,1X,A,1X,A,5X,A)')&
                 '        IBOND','ATOM','ATOM','+/-'
            DO i=1,vdwwfi%nboadwf
               IF (boadwf(3,i) > 0) THEN
                  csign='+'
               ELSE IF (boadwf(3,i) < 0) THEN
                  csign='-'
               ENDIF
               WRITE(output_unit,'(8X,3I5,6X,A)')&
                    I,BOADWF(1,I),BOADWF(2,I),CSIGN
            ENDDO
         ENDIF
      ENDIF
      IF (vdwwfl%trwannc) THEN
         WRITE(output_unit,'(A)')&
              '    RESTART WITH OLD WF CENTERS AND SPREAD'
      ENDIF
      IF (vdwwfl%treffit) THEN
         WRITE(output_unit,'(A)')&
              '    USE REFERENCES FITTED TO WF CENTERS'
      ELSE
         WRITE(output_unit,'(A)')&
              '    USE REFERENCES OF ATOM/BOND CENTERS'
      ENDIF
      IF (vdwwfl%tdampda) THEN
         WRITE(output_unit,'(A)')&
              '    USE DIPOLE APPROXIMATION PAULI EXCHANGE OF DAMPING FUNCTION'
         WRITE(output_unit,'(A)')&
              '    SILVESTRELLI-AMBROSETTI J. CHEM. PHYS. 150, 164109 (2019)'
      ELSE
         WRITE(output_unit,'(A,T56,F10.3)')&
              '    USE SEMIEMPIRICAL DAMPING FUNCTION WITH A6 FACTOR:',vdwwfr%a6
      ENDIF
      WRITE(output_unit,'(A,T56,F10.3)')&
           '    TOLERANCE FOR EXTRAPOLATION OF WF CENTERS(A.U.):',vdwwfr%tolwann
      WRITE(output_unit,'(A,T56,F10.3)')&
           '    TOLERANCE RATIO FOR REFERENCE IONS OF WF CENTERS:',vdwwfr%tolref
      WRITE(output_unit,'(A,A,T52,I2,A,I2,A,I2,A)') '    VAN DER WAALS SUM',&
           ' IN REAL SPACE OVER ',vdwi%nxvdw+1,'*',vdwi%nyvdw+1,'*',vdwi%nzvdw+1,' CELLS'
      WRITE(output_unit,*)

    END SUBROUTINE wannvdw_report
    ! ==--------------------------------------------------------------==
  END SUBROUTINE wannvdwin
  ! ==================================================================
  FUNCTION set_s6(str)
    ! ==--------------------------------------------------------------==
    ! == This function returns proper S6 value w.r.t. the functional  ==
    ! ==--------------------------------------------------------------==
    CHARACTER(len=*), PARAMETER              :: procedureN='set_s6'
    CHARACTER(len=8)                         :: str
    REAL(real_8)                             :: set_s6

    IF (str == " ") THEN
       set_s6=0.0_real_8
       RETURN
    ENDIF
    ! Functional-dependent factor for DFT-D2
    IF (str == "BLYP" ) THEN
       set_s6 = 1.20_real_8
    ELSE IF (str == "B3LYP") THEN
       set_s6 = 1.05_real_8
    ELSE IF (str == "BP86") THEN
       set_s6 = 1.05_real_8
    ELSE IF (str == "TPSS") THEN
       set_s6 = 1.00_real_8
    ELSE IF (str == "PBE") THEN
       set_s6 = 0.75_real_8
    ELSE IF (str == "REVPBE") THEN
       set_s6 = 1.25_real_8
    ELSE IF (str == "PBE0") THEN
       set_s6 = 0.60_real_8
    ELSE
       set_s6=0._real_8
       IF (paral%io_parent)&
            WRITE(output_unit,'(A,A)') ' SELECTED VDW FUNCTIONAL = ', str
       CALL stopgm(procedureN,&
            'Grimme vdw does not support this functional: '//trim(adjustl(str)),& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION set_s6
  ! ==================================================================
  FUNCTION dftd3func(dftfunc)
    ! ==--------------------------------------------------------------==
    ! == This function converts functional to that used in dftd3 code ==
    ! ==--------------------------------------------------------------==
    CHARACTER(len=*), INTENT(in)               :: dftfunc
    CHARACTER(len=input_string_len)            :: dftd3func
    CHARACTER(len=*), PARAMETER                :: procedureN='dftd3func'

    SELECT CASE(dftfunc)
    CASE('BLYP')
       dftd3func='b-lyp'
    CASE('B3LYP')
       dftd3func='b3-lyp'
    CASE('BP86')
       dftd3func='b-p'
    CASE('TPSS')
       dftd3func='tpss'
    CASE('PBE')
       dftd3func='pbe'
    CASE('REVPBE')
       dftd3func='revpbe'
    CASE('PBE0')
       dftd3func='pbe0'
    CASE('HSE06')
       dftd3func='hse06'
    CASE('REVPBE0')
       dftd3func='revpbe0'
    CASE('CAM-B3LYP')
       dftd3func='cam-b3lyp'
    CASE('HCTH')
       dftd3func='hcth120'
    CASE('PW91')
       dftd3func='pwgga'
    CASE('OLYP')
       dftd3func='o-lyp'
    CASE('PBES')
       dftd3func='pbesol'
    CASE('OPBE')
       dftd3func='opbe'
    CASE DEFAULT
       CALL stopgm(procedureN,&
            'unknown functional ('//TRIM(ADJUSTL(dftfunc))//') to convert for dftd3',&
            __LINE__,__FILE__)
    END SELECT

  END FUNCTION dftd3func
  ! ==================================================================

END MODULE vdwin_utils
