MODULE fileopenmod
  IMPLICIT NONE

  ! FLAGS AND STATUS INFORMATION FOR FILEOPEN AND RELATED CALLS.
  ! DEFAULT SETTINGS: OPEN WITH ACCESS='UNKNOWN', FORM='FORMATTED' AND QUIET.
  ! 
  ! FILE ACCESS MODE FLAGS (EXCLUSIVE)
  ! - FO_DEF:   USE DEFAULT
  ! - FO_NEW:   NEW FILE, OLD FILE MUST NOT EXIST. FERROR=.TRUE. IF IT DOES
  ! - FO_OLD:   OLD FILE, OLD FILE MUST EXISTS. FERROR=.TRUE. IF IT DOES NOT
  ! - FO_APP:   CREATE OR APPEND. OLD FILE WILL BE APPENDED, ELSE CREATE NEW
  ! 
  ! MODIFIERS (ADDITIVE)
  ! - FO_UFO:   WRITE UNFORMATTED OUTPUT
  ! - FO_VERB:  BE VERBOSE
  ! - FO_MARK:  ADD << NEW DATA >> MARKER
  ! - FO_NOCHK: DO NOT CHECK WHETHER THIS CHANNEL IS ALREADY OPEN.
  ! - FO_NOFP:  DO NOT PREPEND FILEPATH (= FPATH(IAFPATH:IEFPATH))
  ! - FO_DEBUG: DEBUG FLAG, TRACE ALL ACTIONS VERY VERBOSELY
  ! 
  ! - FO_VMARK: = FO_VERB+FO_MARK
  ! 
  INTEGER, PARAMETER :: fo_def=0,fo_new=1,fo_old=2,fo_app=3,&
       fo_ufo=4,fo_verb=8,fo_mark=16,fo_nochk=32,&
       fo_nofp=64,fo_debug=128,fo_scratch=256,fo_vmark=(fo_verb+fo_mark)

  ! STORE THE OPEN STATUS OF UNITS.
  INTEGER, PARAMETER :: fo_stat_max=100,fo_fpath_max=1024 

  ! ==================================================================
  TYPE :: fo_info_t
     LOGICAL :: fo_tdebug
     LOGICAL :: fo_stat(fo_stat_max)
     INTEGER :: iapath
     INTEGER :: iepath
     CHARACTER(len=fo_fpath_max) :: fpath
  END TYPE fo_info_t
  TYPE(fo_info_t) :: fo_info

END MODULE fileopenmod
