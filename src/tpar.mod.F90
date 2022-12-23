MODULE tpar
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: maxsp

  IMPLICIT NONE

  ! ==================================================================
  ! == STEPS FOR GEOMETRY OPTIMIZATION OR MOLECULAR DYNAMICS        ==
  ! ==--------------------------------------------------------------==
  ! == DELT: Value given by the option TIMESTEP                     ==
  ! ==       used by the steepest descent                           ==
  ! ==--------------------------------------------------------------==
  ! == All these quantities are initialized in dynit.F              ==
  ! == DT_ELEC = DELT_ELEC                                          ==
  ! == DT_IONS = DELT_IONS                                          ==
  ! == DTB2ME = DELT_ELEC/(2.0_real_8*EMASS)                             ==
  ! == DT2_ELEC  = DELT_ELEC*DELT_ELEC                              ==
  ! == DT2_IONS  = DELT_IONS*DELT_IONS                              ==
  ! == DT2BYE = DT2_ELEC/EMASS                                      ==
  ! == DT2HBE = DT2BY2_ELEC/EMASS                                   ==
  ! == DT2BYM(1:NSP) step for steepest descent                      ==
  ! ==        = DT2BY2_IONS/PMA(IS)                                 ==
  ! == DTB2MI(1:NSP)                                                ==
  ! ==        = DELT_IONS/(2.0_real_8*PMA(IS))                           ==
  ! ==--------------------------------------------------------------==
  REAL(real_8), SAVE :: dt_elec,&
       dt2_elec,&
       dtb2me,&
       dt2hbe
  REAL(real_8), SAVE :: dt2bye=0.0_real_8  !vw should that set to zero as default?

  ! ==--------------------------------------------------------------==
  REAL(real_8), SAVE :: dt_ions,&
       dt2_ions,&
       dt2bym(maxsp)=0.0_real_8,&
       dtb2mi(maxsp)
  ! ==================================================================

END MODULE tpar
