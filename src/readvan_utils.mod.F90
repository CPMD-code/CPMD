MODULE readvan_utils
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE pslo,                            ONLY: pslo_com
  USE readsr_utils,                    ONLY: xstring
  USE system,                          ONLY: maxsys,&
                                             nbrx
  USE vdbp,                            ONLY: &
       betar, dion, ncpr1, qfunc, qqq, qrl, r, rab, rsatom, rscore, ru, &
       rucore, vdb_pawf, vdb_r
  USE vdbt,                            ONLY: itmax,&
                                             vdbti
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: readvan

CONTAINS

  ! ==================================================================
  SUBROUTINE readvan(is,filename)
    ! ==--------------------------------------------------------------==
    ! READ ULTRASOFT VANDERBILT PSEUDOPOTENTIALS
    ! ==--------------------------------------------------------------==
    ! 
    ! ==--------------------------------------------------------------==
    ! explanatory comments
    ! ==--------------------------------------------------------------==
    ! 
    ! title    name of chemical element (character*20)
    ! z        atomic number of element
    ! zv       pseudo atomic number (net charge of bare ion core)
    ! exfact   encodes type of exchange-correlation to use
    ! if (exfact.eq.-1.) xctit='              wigner'
    ! if (exfact.eq.-2.) xctit='     hedin-lundqvist'
    ! if (exfact.eq.-3.) xctit=' gunnarson-lundqvist'
    ! if (exfact.gt.0.)  xctit='      slater x-alpha'
    ! nvales   number of l states to include
    ! nvales=1  -->  s
    ! nvales=2  -->  s,p
    ! nvales=3  -->  s,p,d
    ! mesh     number of radial mesh points
    ! etot     total energy of pseudo atom in reference config
    ! nnlz       100 s place --> n quantum number
    ! 10 s place --> l quantum number
    ! 1 s place --> (m quantum number) (really always 0)
    ! wwnl     occupation
    ! ee       eigenvalue
    ! 
    ! Note: nnlz,wwnl,ee give info on pseudo eigenstates; these are
    ! always listed in the order s, p, d in current version
    ! 
    ! mesh     number of radial mesh points
    ! keyps    encoding of type of pseudopotential:
    ! 0 --> standard hsc pseudopotential with exponent 4.0
    ! 1 --> standard hsc pseudopotential with exponent 3.5
    ! 2 --> vanderbilt modifications using defaults
    ! 3 --> new generalized eigenvalue pseudopotentials
    ! 4 --> frozen core all-electron case
    ! ifpcor   1 if "nonlinear core correction" of louie, froyen,
    ! and cohen to be used; 0 otherwise
    ! rinner   radius at which to cut off partial core or q_ij
    ! 
    ! For true frozen core case, use keyps=4, ifpcor=1, rinner=0.
    ! For keyps=3:
    ! 
    ! rc       cutoff radii for s,p,d respectively
    ! nbeta    number of beta functions (sum over all l)
    ! kkbeta   last radial mesh point used to describe functions
    ! which vanish outside core
    ! lll      lll(j) is l quantum number of j th beta function
    ! eee      energy at which construction was done
    ! betar    beta function
    ! dion     bare pseudopotential matrix (ionic and screening
    ! parts subtracted out)
    ! ddd      screened pseudopotential matrix (reference config)
    ! qqq      Q_ij matrix
    ! qfunc    Q_ij(r) function
    ! qfcoef   coefficients to pseudize qfunc for different total
    ! angular momentum (see below)
    ! rcloc    cutoff radius used to construct local potential
    ! rucore   bare local potential (see note below)
    ! rscore   partial core charge
    ! ru       screened local potential (see note below)
    ! rsatom   charge density of pseudo atom (reference config)
    ! 
    ! Note: For consistency with HSC pseudopotentials etc., the
    ! program carries a bare local potential rucore for
    ! each l value.  For keyps=3 they are all the same.
    ! In general, ru is the screened s potential (again
    ! for keyps=3 it doesn t matter).
    ! 
    ! ------------------------------------------------------
    ! Important:
    ! ------------------------------------------------------
    ! potentials, e.g. rucore, are really r*v(r)
    ! wave funcs, e.g. snl, are really proportional to r*psi(r)
    ! and are normalized so int dr (snl**2) = 1
    ! thus psi(r-vec)=(1/r)*snl(r)*y_lm(theta,phi)
    ! conventions carry over to beta, etc
    ! charge dens, e.g. rscore, really 4*pi*r**2*rho
    ! ------------------------------------------------------
    ! 
    ! ------------------------------------------------------
    ! Notes on qfunc and qfcoef:
    ! ------------------------------------------------------
    ! Since Q_ij(r) is the product of two orbitals like
    ! psi_{l1,m1}^star * psi_{l2,m2}, it can be decomposed by
    ! total angular momentum L, where L runs over | l1-l2 | ,
    ! | l1-l2 | +2 , ... , l1+l2.  (L=0 is the only component
    ! needed by the atomic program, which assumes spherical
    ! charge symmetry.)
    ! 
    ! Recall  qfunc(r) = y1(r) * y2(r)  where y1 and y2 are the
    ! radial parts of the wave functions defined according to
    ! 
    ! psi(r-vec) = (1/r) * y(r) * Y_lm(r-hat)  .
    ! 
    ! For each total angular momentum L, we pseudize qfunc(r)
    ! inside rc as:
    ! 
    ! qfunc(r) = r**(L+2) * [ a_1 + a_2*r**2 + a_3*r**4 ]
    ! 
    ! in such a way as to match qfunc and its 1 st derivative at
    ! rc, and to preserve
    ! 
    ! integral dr r**L * qfunc(r)   ,
    ! 
    ! i.e., to preserve the L th moment of the charge.  The array
    ! qfunc has been set inside rc to correspond to this pseudized
    ! version using the minimal L, namely L = | l1-l2 | (e.g., L=0
    ! for diagonal elements).  The coefficients a_i (i=1,2,3)
    ! are stored in the array qfcoef(i,L+1,j,k) for each L so that
    ! the correctly pseudized versions of qfunc can be reconstructed
    ! for each L.  (Note that for given l1 and l2, only the values
    ! L = | l1-l2 | , | l1-l2 | +2 , ... , l1+l2 are ever used.)
    ! ------------------------------------------------------
    ! 
    ! Note that some of the variables included in the pseudo
    ! file (e.g. rc, wwnl) will not be used by the solid state
    ! program, but are included to help identify the pseudopo-
    ! tential (i.e. they are printed out in subroutine psrprt).
    ! 
    ! Also arrays like ru, ddd, and rsatom are non-essential,
    ! but are provided so that they can be used for generating
    ! a starting guess at the potential.
    ! 
    ! ==--------------------------------------------------------------==
    ! 
    INTEGER                                  :: is
    CHARACTER(len=*)                         :: filename

    CHARACTER(len=20)                        :: title, xctype
    CHARACTER(len=80)                        :: testerr
    INTEGER :: i, i2, ib, idmy(3), ifip2, ifn, ifqopt, iptype(nbrx), iqf, &
      irel, it, iv, iver(3), j, jv, k, key(3), keyps, lll(nbrx), lloc, lp, &
      nang, natwf, nblock, nlc, nnlz(5), npf, nqf, testbeta(0:4), testlmax, &
      xc_flag
    LOGICAL                                  :: exists
    REAL(real_8) :: aa, aaa(3), ddd(nbrx,nbrx), deltax, ee(5), eee(nbrx), &
      eloc, etot, exfact, ptryc, qf, qfcoef(13,5,nbrx,nbrx), qtryc, rcloc, &
      rinner(5), rpcor, rr, wwnl(5), z

! Variables
! 
! ==--------------------------------------------------------------==

    DO i=1,3
       ions0%rc(is,i) = 0._real_8
    ENDDO
    ! ==--------------------------------------------------------------==
    ifn=22
    IF (paral%io_parent)&
         INQUIRE(file=filename,exist=exists)
    IF (.NOT.exists) THEN
       IF (paral%io_parent)&
            WRITE(6,*) ' READVAN: FILE NOT FOUND ',filename
       CALL stopgm(' READVAN',' ',& 
            __LINE__,__FILE__)
    ENDIF
    ! 
    IF (.NOT.pslo_com%tbin(is)) GOTO 100
    ! 
    ! ==--------------------------------------------------------------==
    ! READ UNFORMATTED INPUT PSEUDOPOTENTIAL FILE
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         OPEN(unit=ifn,file=filename,status='UNKNOWN',form='UNFORMATTED')
    IF (paral%io_parent)&
         REWIND(ifn)
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         READ(ifn) (iver(i),i=1,3),(idmy(i),i=1,3)
    pslo_com%tlog(is) = (iver(1).GT.2)
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         READ(ifn) title,z,ions0%zv(is),exfact,ncpr1%nvales(is),ncpr1%meshva(is),etot
    ions0%iatyp(is) = NINT(z)
    xc_flag = NINT(exfact)
    IF (ncpr1%meshva(is).GT.maxsys%mmaxx)&
         CALL stopgm(' READVAN: ','mesh.gt.maxsys%mmaxx ',& 
         __LINE__,__FILE__)
    IF (ncpr1%nvales(is).GT.5)&
         CALL stopgm(' READVAN: ','nvales.gt.5 ',& 
         __LINE__,__FILE__)
    IF (paral%io_parent)&
         READ(ifn) (nnlz(i),wwnl(i),ee(i),i=1,ncpr1%nvales(is))
    IF (paral%io_parent)&
         READ(ifn) keyps,ncpr1%ifpcor(is),rinner(1)
    IF (iver(1).GE.3) THEN
       IF (paral%io_parent)&
            READ(ifn) nang,lloc,eloc,ifqopt,nqf,qtryc
       IF (nqf.GT.13)&
            CALL stopgm(' READVAN: ','nqf.gt.13 ',& 
            __LINE__,__FILE__)
       IF ( (ifqopt.GE.2) .AND. (iver(1).LT.5) )&
            CALL stopgm(' READVAN: ','unsupported qfcoef',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (10*iver(1)+iver(2).GE.51) THEN
       IF (paral%io_parent)&
            READ(ifn) (rinner(i),i=1,2*nang-1)
    ELSE
       IF (nang.GT.1) THEN
          DO i=2,2*nang-1
             rinner(i) = rinner(1)
          ENDDO
       ENDIF
    ENDIF
    IF (iver(1).GE.4) READ(ifn) irel
    ! ==--------------------------------------------------------------==
    ! set the number of angular momentum terms in q_ij to read in
    IF (iver(1).EQ.1) THEN
       ! no distinction between nang and nvalps
       nang = ncpr1%nvales(is)
       ! no optimisation of q_ij so 3 term taylor series
       nqf = 3
       nlc = 5
    ELSEIF (iver(1).EQ.2) THEN
       ! no distinction between nang and nvalps
       nang = ncpr1%nvales(is)
       ! no optimisation of q_ij so 3 term taylor series
       nqf = 3
       nlc = 2*ncpr1%nvales(is)-1
    ELSE
       nlc = 2*nang-1
    ENDIF
    IF (nlc.GT.5)&
         CALL stopgm(' READVAN: ','NLC.gt.5 ',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF (keyps.LE.2) THEN
       IF (paral%io_parent)&
            READ(ifn) aa,key,aaa,(ions0%rc(is,i),i=1,ncpr1%nvales(is))
       IF (paral%io_parent)&
            READ(ifn) ((rucore(i,lp,is),i=1,ncpr1%meshva(is)),lp=1,ncpr1%nvales(is))
    ELSEIF (keyps.EQ.3) THEN
       IF (paral%io_parent)&
            READ(ifn) (ions0%rc(is,i),i=1,nang)
       IF (paral%io_parent)&
            READ(ifn) ncpr1%nbeta(is),ncpr1%kkbeta(is)
       IF (ncpr1%nbeta(is).GT.nbrx)&
            CALL stopgm(' READVAN: ','nbeta.gt.NBRX',& 
            __LINE__,__FILE__)
       DO j=1,ncpr1%nbeta(is)
          IF (paral%io_parent)&
               READ(ifn) lll(j),eee(j),(betar(i,j,is),i=1,ncpr1%kkbeta(is))
          DO k=j,ncpr1%nbeta(is)
             IF (paral%io_parent)&
                  READ(ifn) dion(j,k,is),ddd(j,k),qqq(j,k,is),&
                  (qfunc(i,j,k,is),i=1,ncpr1%kkbeta(is)),&
                  ((qfcoef(i,lp,j,k),i=1,nqf),lp=1,nlc)
             dion(k,j,is) = dion(j,k,is)
             ddd(k,j) = ddd(j,k)
             qqq(k,j,is) = qqq(j,k,is)
             DO i=1,nqf
                DO lp=1,nlc
                   qfcoef(i,lp,k,j) = qfcoef(i,lp,j,k)
                ENDDO
             ENDDO
             DO i=1,ncpr1%kkbeta(is)
                qfunc(i,k,j,is) = qfunc(i,j,k,is)
             ENDDO
          ENDDO
       ENDDO
       IF (10*iver(1)+iver(2).GE.72) THEN
          IF (paral%io_parent)&
               READ(ifn) (iptype(j),j=1,ncpr1%nbeta(is)),npf,ptryc
       ENDIF
       IF (paral%io_parent)&
            READ(ifn) rcloc,(rucore(i,1,is),i=1,ncpr1%meshva(is))
    ELSEIF (keyps.EQ.4) THEN
       IF (paral%io_parent)&
            READ(ifn) ((rucore(i,lp,is),i=1,ncpr1%meshva(is)),lp=1,ncpr1%nvales(is))
    ELSE
       CALL stopgm(' READVAN','KEYPS',& 
            __LINE__,__FILE__)
    ENDIF
    ! old
    ! IF(IFPCOR(IS).GT.0) READ(IFN) (RSCORE(I,IS),I=1,MESHVA(IS))
    ! old
    IF (ncpr1%ifpcor(is).GT.0) THEN
       IF (paral%io_parent)&
            READ(ifn) rpcor
       IF (paral%io_parent)&
            READ(ifn) (rscore(i,is),i=1,ncpr1%meshva(is))
    ENDIF
    IF (paral%io_parent)&
         READ(ifn) (ru(i,is),i=1,ncpr1%meshva(is))
    IF (paral%io_parent)&
         READ(ifn) (rsatom(i,is),i=1,ncpr1%meshva(is))
    ! ==--------------------------------------------------------------==
    ! possible further read in of information in log meshva(is) case
    IF ( pslo_com%tlog(is) ) THEN
       IF (paral%io_parent)&
            READ(ifn) (r(i,is),i=1,ncpr1%meshva(is))
       ! 
       DO i=1,ncpr1%meshva(is)
          vdb_r(i,is)=r(i,is)
       ENDDO
       ! 
       IF (paral%io_parent)&
            READ(ifn) (rab(i,is),i=1,ncpr1%meshva(is))
    ENDIF
    ! Read in pseudo atomic wavefunction
    IF (iver(1) .GE. 6) THEN
       natwf = ncpr1%nvales(is)
       IF (iver(1) .GE. 7) READ (ifn) natwf
       IF (paral%io_parent)&
            READ(ifn) ((vdb_pawf(is,i,j), i=1,ncpr1%meshva(is)),j=1,natwf)
    ENDIF
    ! 
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         CLOSE(ifn)
    GOTO 200
    ! 
100 CONTINUE
    ! ==--------------------------------------------------------------==
    ! READ PSEUDOPOTENTIAL FILE (FREE ASCII FORMAT)
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         OPEN(unit=ifn,file=filename,status='UNKNOWN',form='FORMATTED')
    IF (paral%io_parent)&
         REWIND(ifn)
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         READ(ifn,*) (iver(i),i=1,3),(idmy(i),i=1,3)
    pslo_com%tlog(is) = (iver(1).GT.2)
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         READ(ifn,'(A20,3F15.9/2I5,1PE19.11)') title,z,ions0%zv(is),exfact,&
         ncpr1%nvales(is),ncpr1%meshva(is),etot
    ions0%iatyp(is) = NINT(z)
    xc_flag = NINT(exfact)
    IF (ncpr1%meshva(is).GT.maxsys%mmaxx)&
         CALL stopgm(' READVAN: ','mesh.gt.maxsys%mmaxx ',& 
         __LINE__,__FILE__)
    IF (ncpr1%nvales(is).GT.5)&
         CALL stopgm(' READVAN: ','nvales.gt.5 ',& 
         __LINE__,__FILE__)
    IF (paral%io_parent)&
         READ(ifn,*) (nnlz(i),wwnl(i),ee(i),i=1,ncpr1%nvales(is))
    IF (paral%io_parent)&
         READ(ifn,*) keyps,ncpr1%ifpcor(is),rinner(1)
    IF (iver(1).GE.3) THEN
       IF (paral%io_parent)&
            READ(ifn,*) nang,lloc,eloc,ifqopt,nqf,qtryc
       IF (nqf.GT.13)&
            CALL stopgm(' READVAN: ','nqf.gt.13 ',& 
            __LINE__,__FILE__)
       IF ( (ifqopt.GE.2) .AND. (iver(1).LT.5) )&
            CALL stopgm(' READVAN: ','unsupported qfcoef',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (10*iver(1)+iver(2).GE.51) THEN
       IF (paral%io_parent)&
            READ(ifn,*) (rinner(i),i=1,2*nang-1)
    ELSE
       IF (nang.GT.1) THEN
          DO i=2,2*nang-1
             rinner(i) = rinner(1)
          ENDDO
       ENDIF
    ENDIF
    IF (iver(1).GE.4) READ(ifn,*) irel
    ! ==--------------------------------------------------------------==
    ! set the number of angular momentum terms in q_ij to read in
    IF (iver(1).EQ.1) THEN
       ! no distinction between nang and nvalps
       nang = ncpr1%nvales(is)
       ! no optimisation of q_ij so 3 term taylor series
       nqf = 3
       nlc = 5
    ELSEIF (iver(1).EQ.2) THEN
       ! no distinction between nang and nvalps
       nang = ncpr1%nvales(is)
       ! no optimisation of q_ij so 3 term taylor series
       nqf = 3
       nlc = 2*ncpr1%nvales(is)-1
    ELSE
       nlc = 2*nang-1
    ENDIF
    IF (nlc.GT.5)&
         CALL stopgm(' READVAN: ','NLC.gt.5 ',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF (keyps.LE.2) THEN
       IF (paral%io_parent)&
            READ(ifn,*) aa,key,aaa,(ions0%rc(is,i),i=1,ncpr1%nvales(is))
       IF (paral%io_parent)&
            READ(ifn,*) ((rucore(i,lp,is),i=1,ncpr1%meshva(is)),lp=1,ncpr1%nvales(is))
    ELSEIF (keyps.EQ.3) THEN
       IF (paral%io_parent)&
            READ(ifn,*) (ions0%rc(is,i),i=1,nang)
       IF (paral%io_parent)&
            READ(ifn,*) ncpr1%nbeta(is),ncpr1%kkbeta(is)
       IF (ncpr1%nbeta(is).GT.nbrx)&
            CALL stopgm(' READVAN: ','nbeta.gt.NBRX',& 
            __LINE__,__FILE__)
       DO j=1,ncpr1%nbeta(is)
          IF (paral%io_parent)&
               READ(ifn,*) lll(j),eee(j),(betar(i,j,is),i=1,ncpr1%kkbeta(is))
          DO k=j,ncpr1%nbeta(is)
             IF (paral%io_parent)&
                  READ(ifn,*) dion(j,k,is),ddd(j,k),qqq(j,k,is),&
                  (qfunc(i,j,k,is),i=1,ncpr1%kkbeta(is)),&
                  ((qfcoef(i,lp,j,k),i=1,nqf),lp=1,nlc)
             dion(k,j,is) = dion(j,k,is)
             ddd(k,j) = ddd(j,k)
             qqq(k,j,is) = qqq(j,k,is)
             DO i=1,nqf
                DO lp=1,nlc
                   qfcoef(i,lp,k,j) = qfcoef(i,lp,j,k)
                ENDDO
             ENDDO
             DO i=1,ncpr1%kkbeta(is)
                qfunc(i,k,j,is) = qfunc(i,j,k,is)
             ENDDO
          ENDDO
       ENDDO
       IF (10*iver(1)+iver(2).GE.72) THEN
          IF (paral%io_parent)&
               READ(ifn,*) (iptype(j),j=1,ncpr1%nbeta(is)),npf,ptryc
       ENDIF
       IF (paral%io_parent)&
            READ(ifn,*) rcloc,(rucore(i,1,is),i=1,ncpr1%meshva(is))
    ELSEIF (keyps.EQ.4) THEN
       IF (paral%io_parent)&
            READ(ifn,*) ((rucore(i,lp,is),i=1,ncpr1%meshva(is)),lp=1,ncpr1%nvales(is))
    ELSE
       CALL stopgm(' READVAN','KEYPS',& 
            __LINE__,__FILE__)
    ENDIF
    ! old
    ! IF(IFPCOR(IS).GT.0) READ(IFN,*) (RSCORE(I,IS),I=1,MESHVA(IS))
    ! old
    IF (ncpr1%ifpcor(is).GT.0) THEN
       IF (paral%io_parent)&
            READ(ifn,*) rpcor
       IF (paral%io_parent)&
            READ(ifn,*) (rscore(i,is),i=1,ncpr1%meshva(is))
    ENDIF
    IF (paral%io_parent)&
         READ(ifn,*) (ru(i,is),i=1,ncpr1%meshva(is))
    IF (paral%io_parent)&
         READ(ifn,*) (rsatom(i,is),i=1,ncpr1%meshva(is))
    ! ==--------------------------------------------------------------==
    ! possible further read in of information in log meshva(is) case
    IF ( pslo_com%tlog(is) ) THEN
       IF (paral%io_parent)&
            READ(ifn,*) (r(i,is),i=1,ncpr1%meshva(is))
       DO i=1,ncpr1%meshva(is)
          vdb_r(i,is)=r(i,is)
       ENDDO
       IF (paral%io_parent)&
            READ(ifn,*) (rab(i,is),i=1,ncpr1%meshva(is))
    ENDIF
    ! Read in pseudo atomic wavefunction
    IF (iver(1) .GE. 6) THEN
       natwf = ncpr1%nvales(is)
       IF (iver(1) .GE. 7) READ (ifn,*) natwf
       IF (paral%io_parent)&
            READ(ifn,*) ((vdb_pawf(is,i,j), i=1,ncpr1%meshva(is)),j=1,natwf)
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent)&
         CLOSE(ifn)
    ! 
200 CONTINUE
    ! 
    ! ==--------------------------------------------------------------==
    ! ==  G E N E R A T E  H E R M A N  S K I L L M A N  M E S H      ==
    ! ==--------------------------------------------------------------==
    IF ( .NOT. pslo_com%tlog(is) ) THEN
       nblock = (ncpr1%meshva(is))/40
       i = 1
       r(i,is) = 0.0_real_8
       ncpr1%cmesh(is) = 0.88534138_real_8/z**(1._real_8/3._real_8)
       deltax = 0.0025_real_8*ncpr1%cmesh(is)
       DO j=1,nblock
          DO k=1,40
             i = i+1
             r(i,is) = r(i-1,is) + deltax
          ENDDO
          deltax = deltax + deltax
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == SMOOTHING
    ! ==--------------------------------------------------------------==
    pslo_com%tpseu(is) = .FALSE.
    IF ( (iver(1).GT.5) .AND. (ifqopt.GE.2) ) pslo_com%tpseu(is) = .TRUE.
    ! 
    IF (pslo_com%tpseu(is)) THEN
       ! 
       DO lp=1,nlc
          ! 
          DO i2 = ncpr1%kkbeta(is),1,-1
             IF (r(i2,is) .LT. rinner(lp)) GOTO 300
          ENDDO
          CALL stopgm(' READVAN','rinner too small',& 
               __LINE__,__FILE__)
300       CONTINUE
          ! 
          DO iv=1,ncpr1%nbeta(is)
             DO jv=iv,ncpr1%nbeta(is)
                ! 
                DO i=1,ncpr1%kkbeta(is)
                   qrl(i,iv,jv,lp,is) = qfunc(i,iv,jv,is)
                ENDDO
                DO i=1,i2
                   rr = r(i,is)*r(i,is)
                   qf = qfcoef(1,lp,iv,jv)
                   DO iqf =2,nqf
                      qf = qf + qfcoef(iqf,lp,iv,jv)*rr**(iqf-1)
                   ENDDO
                   qrl(i,iv,jv,lp,is) = qf*r(i,is)**(lp+1)
                ENDDO
                DO i=1,ncpr1%kkbeta(is)
                   qrl(i,jv,iv,lp,is) = qrl(i,iv,jv,lp,is)
                ENDDO
                ! 
             ENDDO
          ENDDO
          ! 
       ENDDO
       ! 
    ENDIF
    ! ==--------------------------------------------------------------==
    ! ==       P S E U D O P O T E N T I A L  R E P O R T             ==
    ! ==--------------------------------------------------------------==
    ! 
    xctype = ' '
    IF (xc_flag.EQ. 0) xctype = '      ceperley-alder'
    IF (xc_flag.EQ.-1) xctype = '              wigner'
    IF (xc_flag.EQ.-2) xctype = '     hedin-lundqvist'
    IF (xc_flag.EQ.-3) xctype = ' gunnarson-lundqvist'
    IF (xc_flag.EQ. 1) xctype = ' C-A + B88gx + LYPgc'
    IF (xc_flag.EQ. 2) xctype = ' C-A + B88gx        '
    IF (xc_flag.EQ. 3) xctype = ' C-A + B88gx + P86gc'
    IF (xc_flag.EQ. 4) xctype = ' Perdew Wang 1991   '
    IF (xc_flag.EQ. 5) xctype = ' PBE - GGA          '
    IF (xc_flag.EQ. 6) xctype = ' revPBE - GGA       '
    IF (xc_flag.GT.99) THEN
       IF (paral%io_parent)&
            WRITE(xctype,'(5x,I6)') xc_flag
    ENDIF
    ! ==--------------------------------------------------------------==
    it=1
    IF (paral%io_parent)&
         WRITE(vdbti(it,is),'(3X,60("="))')
    it=it+1
    IF (paral%io_parent)&
         WRITE(vdbti(it,is),1000) (iver(i),i=1,3),idmy(2),idmy(1),idmy(3)
1000 FORMAT(3x,'|  pseudopotential report: version',i3,&
         '.',i1,'.',i1,' date',i3,'-',i2,'-',i4,2x,'|')
    it=it+1
    IF (paral%io_parent)&
         WRITE(vdbti(it,is),'(3X,60("-"))')
    it=it+1
    IF (paral%io_parent)&
         WRITE(vdbti(it,is),1010) title,xctype
1010 FORMAT(3x,'|  ',2a20,' exchange-corr  |')
    it=it+1
    IF (paral%io_parent)&
         WRITE(vdbti(it,is),1020) z,ions0%zv(is),exfact
1020 FORMAT(3x,'|  z =',f7.2,2x,'zv =',f7.2,2x,'exfact =',f10.5,&
         13x,'|')
    it=it+1
    IF (paral%io_parent)&
         WRITE(vdbti(it,is),1030) etot
1030 FORMAT(3x,'|     ',9x,'    ',9x,' etot  =',f10.5,13x,'|')
    it=it+1
    IF (paral%io_parent)&
         WRITE(vdbti(it,is),1040)
1040 FORMAT(3x,'|  index    orbital      occupation    energy',14x,'|')
    DO i=1,ncpr1%nvales(is)
       it=it+1
       IF (paral%io_parent)&
            WRITE(vdbti(it,is),1050) i,nnlz(i),wwnl(i),ee(i)
    ENDDO
1050 FORMAT(3x,'|',i5,i11,5x,f10.2,f12.2,15x,'|')
    ! 
    it=it+1
    IF (10*iver(1)+iver(2).GE.51) THEN
       IF (paral%io_parent)&
            WRITE(vdbti(it,is),1061) keyps,ncpr1%ifpcor(is)
1061   FORMAT(3x,'|  keyps =',i2,5x,'ifpcor =',i2,32x,'|')
       DO i = 1,2*nang-1
          it=it+1
          IF (paral%io_parent)&
               WRITE(vdbti(it,is),1065) rinner(i),i
       ENDDO
1065   FORMAT(3x,'|  rinner =',f10.2,5x,'for L=',i5,22x,'|')
    ELSE
       IF (paral%io_parent)&
            WRITE(vdbti(it,is),1060) keyps,ncpr1%ifpcor(is),rinner(1)
1060   FORMAT(3x,'|  keyps =',i2,5x,'ifpcor =',i2,5x,'rinner =',&
            f10.4,9x,'|')
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (keyps.LE.2) THEN
       ! ==--------------------------------------------------------------==
    ELSEIF (keyps.EQ.3) THEN
       ! new scheme
       it=it+1
       IF (paral%io_parent)&
            WRITE(vdbti(it,is),1120)
1120   FORMAT(3x,'|    new generation scheme:',32x,'|')
       it=it+1
       IF (paral%io_parent)&
            WRITE(vdbti(it,is),1130) ncpr1%nbeta(is),ncpr1%kkbeta(is),rcloc
1130   FORMAT(3x,'|    nbeta = ',i2,5x,'kkbeta =',i5,5x,&
            'rcloc =',f10.4,4x,'|')
       it=it+1
       IF (paral%io_parent)&
            WRITE(vdbti(it,is),1135)
1135   FORMAT(3x,'|    ibeta     l     epsilon   rcut iptype',17x,'|')
       CALL zeroing(testbeta)!,5)
       testlmax=-1
       DO ib=1,ncpr1%nbeta(is)
          lp=lll(ib)+1
          it=it+1
          IF (lll(ib).GT.testlmax) testlmax=lll(ib)
          testbeta(lll(ib))=testbeta(lll(ib))+1
          IF (10*iver(1)+iver(2).LT.72) iptype(ib)=-1
          ! -1 means iptype was not recorded and is unknown
          IF (paral%io_parent)&
               WRITE(vdbti(it,is),1150) ib,lll(ib),eee(ib),ions0%rc(is,lp),&
               iptype(ib)
       ENDDO
1150   FORMAT(3x,'|',6x,i2,5x,i2,4x,2f7.2,6x,i1,18x,'|')

       ! AK: the current (28.05.2004) vanderbilt implementation requires
       ! the same number of projectors for all angular momentums and
       ! does not allow to skip one. check this.
       DO ib=1,testlmax
          IF (testbeta(ib).NE.testbeta(0)) THEN
             IF (paral%parent) THEN
                CALL xstring(filename,i,i2)
                testerr = filename(i:i2) // ' INCOMPATIBLE WITH CPMD'
                CALL stopgm('READVAN', testerr,& 
                     __LINE__,__FILE__)
             ENDIF
          ENDIF
       ENDDO
       ! 
       IF (10*iver(1)+iver(2).GE.72) THEN
          ifip2=0
          DO ib=1,ncpr1%nbeta(is)
             IF (iptype(ib).EQ.2) ifip2=1
          ENDDO
          IF (ifip2.EQ.1) THEN
             it=it+1
             IF (paral%io_parent)&
                  WRITE(vdbti(it,is),1155) npf,ptryc
1155         FORMAT(3x,'|  npf    =',i2,'  ptryc =',f8.3,29x,'|')
          ENDIF
       ENDIF
       ! 
       IF (iver(1).GT.2) THEN
          it=it+1
          IF (paral%io_parent)&
               WRITE(vdbti(it,is),1160) lloc,eloc
1160      FORMAT(3x,'|  lloc   =',i2,'  eloc   =',f8.3,28x,'|')
          it=it+1
          IF (paral%io_parent)&
               WRITE(vdbti(it,is),1170) ifqopt,nqf,qtryc
1170      FORMAT(3x,'|  ifqopt =',i2,'  nqf    =',i2,'  qtryc =',f8.3,&
               17x,'|')
       ENDIF
       ! 
       it=it+1
       IF ( (iver(1).GT.3) .AND. (irel.EQ.2) ) THEN
          IF (paral%io_parent)&
               WRITE(vdbti(it,is),1180) 'koelling-harmon equation'
       ELSE
          IF (paral%io_parent)&
               WRITE(vdbti(it,is),1180) 'schroedinger equation   '
       ENDIF
1180   FORMAT(3x,'|',2x,'all electron calculation used ',a24,2x,'|')
       ! 
       it=it+1
       IF ( .NOT. pslo_com%tlog(is) ) THEN
          IF (paral%io_parent)&
               WRITE(vdbti(it,is),2000) 'herman skillman mesh'
       ELSE
          IF (paral%io_parent)&
               WRITE(vdbti(it,is),2000) '**logarithmic mesh**'
       ENDIF
2000   FORMAT(3x,'|',9x,10('*'),a20,10('*'),9x,'|')
       ! ==--------------------------------------------------------------==
    ELSEIF (keyps.EQ.4) THEN
       ! frozen core
       ! write (6,2500)
       ! 2500   format(3x,'|    frozen core all-electron case',25x,'|')
    ENDIF
    ! ==--------------------------------------------------------------==
    it=it+1
    IF (paral%io_parent)&
         WRITE(vdbti(it,is),3000)
3000 FORMAT(3x,60("="))
    itmax(is) = it
    ncpr1%nvales(is) = nang
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE readvan

END MODULE readvan_utils
