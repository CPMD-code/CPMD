MODULE control_bcast_utils
  USE andr,                            ONLY: andr2,&
                                             andr3
  USE benc,                            ONLY: ibench
  USE broy,                            ONLY: broy1
  USE cdftmod,                         ONLY: cdftci,&
                                             cdftcom,&
                                             cdftlog,&
                                             cm_dir,&
                                             cm_dr,&
                                             czones,&
                                             sccomm
  USE comvelmod,                       ONLY: comvl
  USE cotr,                            ONLY: sdpl
  USE cp_cuda_types,                   ONLY: cp_cuda_env
  USE error_handling,                  ONLY: stopgm
  USE fileopenmod,                     ONLY: fo_info
  USE g_loc,                           ONLY: gloc_list,&
                                             glocal,&
                                             gloci,&
                                             glocr
  USE glemod,                          ONLY: glepar
  USE ions,                            ONLY: coord_fdiff,&
                                             tref_fdiff
  USE linres,                          ONLY: tshl,&
                                             xfmqc
  USE mergemod,                        ONLY: merge01,&
                                             merge02
  USE mm_input,                        ONLY: clc
  USE mp_interface,                    ONLY: mp_bcast
  USE mw,                              ONLY: mwi,&
                                             tmw
  USE nabdy_types,                     ONLY: nabdyfric,&
                                             nabdyvar
  USE nort,                            ONLY: nort_com
  USE nose,                            ONLY: cafesini,&
                                             cafesinr,&
                                             lctrng,&
                                             loct,&
                                             loctpin,&
                                             ncafesgrp,&
                                             nosl,&
                                             tcafes,&
                                             tnosepc
  USE para_global,                     ONLY: para_buff_size,&
                                             para_stack_buff_size,&
                                             para_use_mpi_in_place
  USE parac,                           ONLY: parai,&
                                             paral
  USE prden,                           ONLY: elfcb,&
                                             mwfn,&
                                             numpr
  USE qspl,                            ONLY: qspl1,&
                                             qsrang
  USE shop,                            ONLY: tshopold
  USE shop_rest,                       ONLY: sh03
  USE spin,                            ONLY: clsd,&
                                             lspin1,&
                                             lspin2,&
                                             lspin3
  USE store_types,                     ONLY: cprint,&
                                             iface1,&
                                             restart1,&
                                             rout1,&
                                             store1,&
                                             trajsmall,&
                                             trajsmalln
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr,&
                                             cp_trace,&
                                             group,&
                                             locpot2,&
                                             restf
  USE time,                            ONLY: tname
  USE vdwcmod,                         ONLY: vdwl
  USE wann,                            ONLY: sw_list,&
                                             wan05,&
                                             wannc,&
                                             wanni,&
                                             wannl,&
                                             wannr
  USE xinr,                            ONLY: inr_logical,&
                                             rmixsd
#include "sizeof.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: control_bcast

CONTAINS

  ! ==================================================================
  SUBROUTINE control_bcast
    ! ==--------------------------------------------------------------==
    ! ==  BROADCAST THE VARIABLES SET IN CONTROL                      ==
    ! ==--------------------------------------------------------------==
    ! Variables
    CHARACTER(*), PARAMETER                  :: procedureN = 'control_bcast'

    INTEGER                                  :: ierr

    CALL mp_bcast_byte(cntl, size_in_bytes_of(cntl), parai%io_source, parai%cp_grp)
    ! CNTI
    CALL mp_bcast_byte(cnti, size_in_bytes_of(cnti), parai%io_source, parai%cp_grp)
    ! CNTR
    CALL mp_bcast_byte(cntr, size_in_bytes_of(cntr), parai%io_source, parai%cp_grp)
    ! RESTART
    CALL mp_bcast_byte(restart1, size_in_bytes_of(restart1), parai%io_source, parai%cp_grp)
    ! ROUT
    CALL mp_bcast_byte(rout1, size_in_bytes_of(rout1), parai%io_source, parai%cp_grp)
    ! TRAJSMALL, TRAJSMALLN
    CALL mp_bcast(trajsmall,parai%io_source,parai%cp_grp)
    CALL mp_bcast(trajsmalln,parai%io_source,parai%cp_grp)
    ! ISTORE
    CALL mp_bcast_byte(store1, size_in_bytes_of(store1), parai%io_source, parai%cp_grp)
    ! QSPL
    CALL mp_bcast_byte(qspl1, size_in_bytes_of(qspl1), parai%io_source, parai%cp_grp)
    ! BENC
    CALL mp_bcast(ibench,SIZE(ibench),parai%io_source,parai%cp_grp)
    ! PRDEN
    CALL mp_bcast(numpr,parai%io_source,parai%cp_grp)
    IF (numpr.NE.0)  THEN
       IF (.NOT.paral%io_parent) THEN
          ALLOCATE(mwfn(numpr),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       CALL mp_bcast(mwfn,SIZE(mwfn),parai%io_source,parai%cp_grp)
    ENDIF
    ! ANDR
    CALL mp_bcast_byte(andr2, size_in_bytes_of(andr2), parai%io_source, parai%cp_grp)
    ! DENSNRDIIS
    CALL mp_bcast_byte(andr3, size_in_bytes_of(andr3), parai%io_source, parai%cp_grp)
    ! BROY
    CALL mp_bcast_byte(broy1, size_in_bytes_of(broy1), parai%io_source, parai%cp_grp)
    ! SPIN (NLSD NLSX)
    CALL mp_bcast_byte(clsd, size_in_bytes_of(clsd), parai%io_source, parai%cp_grp)
    ! ROKS
    CALL mp_bcast_byte(lspin2, size_in_bytes_of(lspin2),parai%io_source,parai%cp_grp)
    CALL mp_bcast(lspin3%mgab,SIZE(lspin3%mgab),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(lspin1, size_in_bytes_of(lspin1), parai%io_source, parai%cp_grp)
    ! NOGRP
    CALL mp_bcast(group%nogrp,parai%io_source,parai%cp_grp)
    CALL mp_bcast(qsrang,parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(elfcb, size_in_bytes_of(elfcb),parai%io_source,parai%cp_grp)
    ! CPRINT
    CALL mp_bcast_byte(cprint, size_in_bytes_of(cprint),parai%io_source,parai%cp_grp)
    ! WANNIER FUNCTION INFORMATION
    CALL mp_bcast_byte(wannc, size_in_bytes_of(wannc),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(wannl, size_in_bytes_of(wannl),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(wanni, size_in_bytes_of(wanni),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(wannr, size_in_bytes_of(wannr),parai%io_source,parai%cp_grp)
    IF (wanni%sw_orb.GT.0) THEN
       IF (.NOT.paral%parent) THEN
          ALLOCATE(sw_list(wanni%sw_orb),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       CALL mp_bcast(sw_list,SIZE(sw_list),parai%io_source,parai%cp_grp)
    ENDIF
    !NABDY[
    CALL mp_bcast(nabdyvar%nabdy_zmax,parai%io_source, parai%cp_grp)
    CALL mp_bcast(nabdyvar%tnarepart,parai%io_source, parai%cp_grp)
    CALL mp_bcast(nabdyvar%scalep,parai%io_source, parai%cp_grp)
    CALL mp_bcast(nabdyvar%naforce_screened,parai%io_source, parai%cp_grp)
    CALL mp_bcast(nabdyvar%nasoftening,parai%io_source, parai%cp_grp)
    CALL mp_bcast(nabdyvar%tnafric,parai%io_source, parai%cp_grp)
    CALL mp_bcast(nabdyvar%nafritc,parai%io_source, parai%cp_grp)
    CALL mp_bcast(nabdyfric%natempcm,parai%io_source, parai%cp_grp)
    CALL mp_bcast(nabdyfric%natempbp,parai%io_source, parai%cp_grp)
    CALL mp_bcast(nabdyvar%ntrajbd,parai%io_source, parai%cp_grp)
    !NABDY]
    ! WAN05
    CALL mp_bcast_byte(wan05,size_in_bytes_of(wan05),parai%io_source,parai%cp_grp)
    ! PRNG SEED
    CALL mp_bcast(cnti%iprng,parai%io_source,parai%cp_grp)
    ! GLE
    CALL mp_bcast_byte(glepar, size_in_bytes_of(glepar),parai%io_source,parai%cp_grp)


    ! RESTF
    CALL mp_bcast_byte(restf%nstepwr,size_in_bytes_of(restf),parai%io_source,parai%cp_grp)

    ! SDPL
    CALL mp_bcast_byte(sdpl, size_in_bytes_of(sdpl),parai%io_source,parai%cp_grp)
    ! CAFES
    CALL mp_bcast(tcafes,parai%io_source,parai%cp_grp)
    IF (tcafes) THEN
       CALL mp_bcast(ncafesgrp,parai%io_source,parai%cp_grp)
       IF (.NOT. paral%parent) THEN
          ALLOCATE(cafesini(2,ncafesgrp),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(cafesinr(2,ncafesgrp),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       CALL mp_bcast(cafesini,SIZE(cafesini),parai%io_source,parai%cp_grp)
       CALL mp_bcast(cafesinr,SIZE(cafesinr),parai%io_source,parai%cp_grp)
    ENDIF
    ! LOCAL TEMPERATURE
    CALL mp_bcast(loct%tloct,parai%io_source,parai%cp_grp)
    IF (loct%tloct) THEN
       CALL mp_bcast(loct%nloct,parai%io_source,parai%cp_grp)
       CALL mp_bcast(loct%nlocrng,parai%io_source,parai%cp_grp)
       IF (.NOT. paral%parent) THEN
          ALLOCATE(loctpin(2,loct%nloct),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          ALLOCATE(lctrng(3,loct%nlocrng),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       CALL mp_bcast(loctpin,SIZE(loctpin),parai%io_source,parai%cp_grp)
       CALL mp_bcast(lctrng,SIZE(lctrng),parai%io_source,parai%cp_grp)
    ENDIF
    ! SURFACE HOPPING  
    IF (cntl%tshop) THEN
       CALL mp_bcast(tshopold,parai%io_source,parai%cp_grp)
       CALL mp_bcast(sh03%tshopres,parai%io_source,parai%cp_grp)
    ENDIF
    ! CLASSICAL
    CALL mp_bcast_byte(clc, size_in_bytes_of(clc),parai%io_source,parai%cp_grp)
    ! TURBOCICCO GEOMETRY OPTIMIZATION
    CALL mp_bcast_byte(inr_logical, size_in_bytes_of(inr_logical),parai%io_source,parai%cp_grp)
    CALL mp_bcast(rmixsd,parai%io_source,parai%cp_grp)
    ! FDIFF
    CALL mp_bcast(tref_fdiff,parai%io_source,parai%cp_grp)
    CALL mp_bcast(coord_fdiff,SIZE(coord_fdiff),parai%io_source,parai%cp_grp)
    ! G_LOCALIZATION OF FUNCTION INFORMATION
    CALL mp_bcast_byte(glocal, size_in_bytes_of(glocal),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(gloci, size_in_bytes_of(gloci),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(glocr, size_in_bytes_of(glocr),parai%io_source,parai%cp_grp)
    IF (gloci%gloc_orb.GT.0) THEN
       IF (.NOT.paral%parent)  THEN
          ALLOCATE(gloc_list(gloci%gloc_orb),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       CALL mp_bcast_byte(gloc_list,SIZE(gloc_list),parai%io_source,parai%cp_grp)
    ENDIF
    ! QMMM (FPATH AND RELATED STUFF)
    CALL mp_bcast_byte(fo_info,size_in_bytes_of(fo_info),parai%io_source,parai%cp_grp)
    ! INTERFACE TO OTHER PW CODES
    CALL mp_bcast_byte(iface1, size_in_bytes_of(iface1),parai%io_source,parai%cp_grp)
    ! COM/ROT SUBTRACTION
    CALL mp_bcast_byte(comvl, size_in_bytes_of(comvl),parai%io_source,parai%cp_grp)
    ! SLIMIT
    CALL mp_bcast_byte(nort_com, size_in_bytes_of(nort_com),parai%io_source,parai%cp_grp)
    ! cntl%cdft
    CALL mp_bcast_byte(cdftlog, size_in_bytes_of(cdftlog),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(cdftci, size_in_bytes_of(cdftci),parai%io_source,parai%cp_grp)
    CALL mp_bcast(czones,SIZE(czones),parai%io_source,parai%cp_grp)
    CALL mp_bcast(cdftcom%cdft_nc,parai%io_source,parai%cp_grp)
    CALL mp_bcast(cdftcom%cdft_v,SIZE(cdftcom%cdft_v),parai%io_source,parai%cp_grp)
    CALL mp_bcast(cdftcom%nother,parai%io_source,parai%cp_grp)
    CALL mp_bcast(cdftcom%cdft_ns,parai%io_source,parai%cp_grp)
    CALL mp_bcast(cdftcom%nsother,parai%io_source,parai%cp_grp)
    ! SYSCOMB
    CALL mp_bcast_byte(sccomm, size_in_bytes_of(sccomm),parai%io_source,parai%cp_grp)
    IF (cntl%tscombm)THEN
       CALL mp_bcast(cm_dr,parai%io_source,parai%cp_grp)
       CALL mp_bcast(cm_dir,parai%io_source,parai%cp_grp)
    ENDIF
    ! VDW CORRECTIONS
    CALL mp_bcast_byte(vdwl, size_in_bytes_of(vdwl),parai%io_source,parai%cp_grp)
    ! FILE MERGE
    CALL mp_bcast_byte(merge01, size_in_bytes_of(merge01),parai%io_source,parai%cp_grp)
    CALL mp_bcast_byte(merge02, size_in_bytes_of(merge02),parai%io_source,parai%cp_grp)
    ! mb-kk -Printing local potential
    CALL mp_bcast_byte(locpot2, size_in_bytes_of(locpot2),parai%io_source,parai%cp_grp)
    ! tracing
    CALL mp_bcast_byte(cp_trace,size_in_bytes_of(cp_trace),parai%io_source,parai%cp_grp)
    CALL mp_bcast(tname%trace_nbr_procedure,parai%io_source,parai%cp_grp)
    CALL mp_bcast(tname%trace_max_depth,parai%io_source,parai%cp_grp)
    CALL mp_bcast(tname%trace_max_calls,parai%io_source,parai%cp_grp)
    CALL mp_bcast(tname%trace_procedure,tname%trace_nbr_procedure,&
         parai%io_source,parai%cp_grp)
    ! nose thermostat
    CALL mp_bcast(nosl%tultra,parai%io_source,parai%cp_grp)
    CALL mp_bcast(nosl%tmnose,parai%io_source,parai%cp_grp)
    CALL mp_bcast(tnosepc,parai%io_source,parai%cp_grp)
    ! PARA_*
    CALL mp_bcast(para_use_mpi_in_place,parai%io_source,parai%cp_grp)
    CALL mp_bcast(para_stack_buff_size,parai%io_source,parai%cp_grp)
    CALL mp_bcast(para_buff_size,parai%io_source,parai%cp_grp)
    ![EXACT FACTORIZATION
    CALL mp_bcast(tshl%txfmqc,parai%io_source,parai%cp_grp)
    IF (tshl%txfmqc) THEN
       CALL mp_bcast(xfmqc%n_xftraj,parai%io_source,parai%cp_grp)
       CALL mp_bcast(tmw,parai%io_source,parai%cp_grp)
       CALL mp_bcast_byte(mwi, size_in_bytes_of(mwi),parai%io_source,parai%cp_grp)
    ENDIF
    !]EXACT FACTORIZATION
    CALL mp_bcast_byte(cp_cuda_env,size_in_bytes_of(cp_cuda_env),parai%io_source,parai%cp_grp)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE control_bcast
  ! ==================================================================

END MODULE control_bcast_utils
