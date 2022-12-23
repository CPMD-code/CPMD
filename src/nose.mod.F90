MODULE nose
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: maxsp

  IMPLICIT NONE
  ! --------------------------------------------------------------------
  ! stuff for the CAFES algorithm, different temperatures for every atom
  ! --------------------------------------------------------------------
  REAL(real_8), ALLOCATABLE :: tempwr(:,:,:)
  REAL(real_8), ALLOCATABLE :: wnosepr(:,:,:)
  REAL(real_8), ALLOCATABLE :: wnosep2r(:,:,:)

  REAL(real_8), ALLOCATABLE :: chainer(:,:,:)

  INTEGER, ALLOCATABLE :: cafescl(:,:)

  INTEGER, ALLOCATABLE :: cafesini(:,:)
  INTEGER :: ncafesgrp

  REAL(real_8), ALLOCATABLE :: cafestmp(:)
  REAL(real_8), ALLOCATABLE :: cafesnse(:)
  REAL(real_8), ALLOCATABLE :: cafesinr(:,:)

  LOGICAL :: tcafes
  ! ==================================================================
  ! ==                  LOCAL TEMPERATUE STUFF                      ==
  ! ==================================================================
  INTEGER, ALLOCATABLE :: lctrng(:,:)
  INTEGER, ALLOCATABLE :: lcttab(:,:)
  INTEGER, ALLOCATABLE :: lctmemb(:,:,:)

  INTEGER, ALLOCATABLE :: nlctmbm(:)

  REAL(real_8), ALLOCATABLE :: loctpin(:,:)
  REAL(real_8), ALLOCATABLE :: dofsl(:)
  REAL(real_8), ALLOCATABLE :: loctt0(:)

  TYPE :: loct_t
     LOGICAL :: tloct
     INTEGER :: nloct
     INTEGER :: nlocrng
  END TYPE loct_t
  TYPE(loct_t) :: loct
  ! ==================================================================
  ! ==                  NOSE STUFF                                  ==
  ! ==================================================================
  INTEGER, PARAMETER :: nchx=20 
  INTEGER, PARAMETER :: msuzx=5 
  INTEGER, PARAMETER :: ncallsx=5**(msuzx-1) 

  TYPE :: nosl_t
     LOGICAL :: tultra
     LOGICAL :: tmnose
     INTEGER :: nose_nrepl
  END TYPE nosl_t
  TYPE(nosl_t) :: nosl
  ! ==--------------------------------------------------------------==
  REAL(real_8) :: wnosee,dofst,wnosep,wnosec,gckt
  REAL(real_8) :: qnosee(nchx),dtsuz(ncallsx),dofsu(maxsp),qnoscc(nchx)
  ! ==--------------------------------------------------------------==
  ! ==--------------------------------------------------------------==
  REAL(real_8), ALLOCATABLE :: eta(:,:)
  REAL(real_8), ALLOCATABLE :: etadot(:,:)
  REAL(real_8), ALLOCATABLE :: etap1(:,:)
  REAL(real_8), ALLOCATABLE :: etap1dot(:,:)
  REAL(real_8), ALLOCATABLE :: gkt(:,:)
  REAL(real_8), ALLOCATABLE :: qnospm(:,:,:)
  REAL(real_8), ALLOCATABLE :: qnosp(:,:,:)
  REAL(real_8), ALLOCATABLE :: qnospc(:,:)
  REAL(real_8), ALLOCATABLE :: etap(:,:,:)
  REAL(real_8), ALLOCATABLE :: etapdot(:,:,:)
  REAL(real_8), ALLOCATABLE :: etapm(:,:,:)
  REAL(real_8), ALLOCATABLE :: etapmdot(:,:,:)
  REAL(real_8), ALLOCATABLE :: etc(:,:)
  REAL(real_8), ALLOCATABLE :: etcdot(:,:)

  LOGICAL :: tnosepc                     ! for centroid of 
  REAL(real_8), ALLOCATABLE :: gkt1(:,:) ! constant NPT CMD

  INTEGER, ALLOCATABLE :: ndfnt(:,:)
  INTEGER, ALLOCATABLE :: mapdof(:,:,:)


  INTEGER :: nit,msuz,ncalls,nche,nedof,nchain,nchc,ncdof
  INTEGER :: ntherm(2)
  REAL(real_8) :: glib
  ! ==================================================================


END MODULE nose
