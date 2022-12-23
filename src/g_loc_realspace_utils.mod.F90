MODULE g_loc_realspace_utils
  USE coor,                            ONLY: tau0
  USE cppt,                            ONLY: indzs,&
                                             nzhs
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: invfftn
  USE fftutil_utils,                   ONLY: phase
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def
  USE g_loc,                           ONLY: gloc_list,&
                                             glocal,&
                                             gloci,&
                                             ind_st,&
                                             lostate
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_recv,&
                                             mp_send,&
                                             mp_sync
  USE nmr_position_p_utils,            ONLY: putongrid,&
                                             rgrid
  USE parac,                           ONLY: parai,&
                                             paral
  USE readsr_utils,                    ONLY: xstring
  USE system,                          ONLY: fpar,&
                                             maxsys,&
                                             ncpw,&
                                             parap,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: p_real_space

CONTAINS

  ! ==================================================================
  SUBROUTINE p_real_space(infi,ikind,c0,nstate,psi,nstate_loc)
    ! ==--------------------------------------------------------------==

    ! ==--------------------------------------------------------------==
    INTEGER                                  :: infi, ikind, nstate
    COMPLEX(real_8)                          :: c0(ncpw%ngw*2,nstate), psi(:)
    INTEGER                                  :: nstate_loc

    CHARACTER(*), PARAMETER                  :: procedureN = 'p_real_space'

    CHARACTER(len=15)                        :: filbod, filen, filtr
    CHARACTER(len=8)                         :: numb
    INTEGER :: i, i1, i2, i3, ia, iat, icont, ie, ierr, ig, ii1, ii2, ii3, &
      ii3p, ilowerleft(3), ir, isp, istate, isub, length, msgid, no, source_pe
    INTEGER, ALLOCATABLE, SAVE               :: orbpri(:)
    INTEGER, SAVE                            :: ifirst = 0
    LOGICAL                                  :: ferror
    REAL(real_8)                             :: cubecenter(3), lowerleft(3), &
                                                rrr, rrx, rry, rrz, xl, xmax, &
                                                xmin, xs, yl, ymax, ymin, ys, &
                                                zl, zmax, zmin, zs
    REAL(real_8), ALLOCATABLE                :: array(:), dest(:,:), r0(:,:,:)

! ==--------------------------------------------------------------==

    CALL tiset('  P_REAL_SPACE',isub)

    length = fpar%kr1*fpar%kr2s*fpar%kr3s + fpar%kr2s*fpar%kr3s

    ALLOCATE(array(fpar%kr1*fpar%kr2s*fpar%kr3s),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(dest(fpar%kr2s, fpar%kr3s),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)



    IF (ifirst.EQ.0) THEN
       no=nstate/2+1
       ALLOCATE(orbpri(no),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(orbpri)!,nstate)
       IF (gloci%gloc_all.NE.0) THEN
          DO i=1,nstate
             orbpri(i)=1
          ENDDO
       ELSEIF (gloci%gloc_first.GT.0) THEN
          DO i=gloci%gloc_first,gloci%gloc_last
             orbpri(i)=1
          ENDDO
       ELSEIF (gloci%gloc_orb.GT.0) THEN
          DO i=1,gloci%gloc_orb
             orbpri(gloc_list(i))=1
          ENDDO
       ENDIF
       ifirst=1
    ENDIF

    ALLOCATE(r0(3,maxsys%nax,ions1%nsp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! SHIFT THE CELL WHEN NECESSARY
    xmin = tau0(1,1,1)
    xmax = tau0(1,1,1)
    ymin = tau0(2,1,1)
    ymax = tau0(1,1,1)
    zmin = tau0(3,1,1)
    zmax = tau0(1,1,1)
    DO isp=1,ions1%nsp
       DO iat=1,ions0%na(isp)
          xmin = MIN(xmin,tau0(1,iat,isp))
          xmax = MAX(xmax,tau0(3,iat,isp))! ZMAX -> XMAX cmb-bugfix
          ymin = MIN(ymin,tau0(2,iat,isp))
          ymax = MAX(ymax,tau0(3,iat,isp))
          zmin = MIN(zmin,tau0(3,iat,isp))
          zmax = MAX(zmax,tau0(3,iat,isp))
       ENDDO
    ENDDO
    DO isp=1,ions1%nsp
       DO iat=1,ions0%na(isp)
          r0(1,iat,isp) = tau0(1,iat,isp) - xmin
          r0(2,iat,isp) = tau0(2,iat,isp) - ymin
          r0(3,iat,isp) = tau0(3,iat,isp) - zmin
       ENDDO
    ENDDO

    xl = parm%a1(1)+parm%a2(1)+parm%a3(1)
    yl = parm%a1(2)+parm%a2(2)+parm%a3(2)
    zl = parm%a1(3)+parm%a2(3)+parm%a3(3)

    xs = (xl - (xmax - xmin))/(2._real_8*xl) ! (4._real_8*XL)
    ys = (yl - (ymax - ymin))/(2._real_8*yl) ! (4._real_8*YL)
    zs = (zl - (zmax - zmin))/(2._real_8*zl) ! (4._real_8*ZL)

    ! write(6,*) XMAX - XMIN,YMAX - YMIN,ZMAX - ZMIN
    ! write(6,*) XL,YL,ZL,XS,YS,ZS


    cubecenter(1) =  (0.5_real_8-xs)*(parm%a1(1)+parm%a2(1)+parm%a3(1))
    cubecenter(2) =  (0.5_real_8-ys)*(parm%a1(2)+parm%a2(2)+parm%a3(2))
    cubecenter(3) =  (0.5_real_8-zs)*(parm%a1(3)+parm%a2(3)+parm%a3(3))
    ! write(6,*) CUBECENTER(1),CUBECENTER(2),CUBECENTER(3) 
    ! do ir=1,3
    ! LOWERLEFT(ir) = CUBECENTER(ir) 
    ! &        - 0.5_real_8*(a1(ir)+a2(ir)+a3(ir))
    ! enddo
    lowerleft(1) = cubecenter(1)&
         - (0.5_real_8-xs)*(parm%a1(1)+parm%a2(1)+parm%a3(1))
    lowerleft(2) = cubecenter(2)&
         - (0.5_real_8-ys)*(parm%a1(2)+parm%a2(2)+parm%a3(2))
    lowerleft(3) = cubecenter(3)&
         - (0.5_real_8-zs)*(parm%a1(3)+parm%a2(3)+parm%a3(3))

    DO ir = 1,3
       ilowerleft(ir) = rgrid(lowerleft,ir)
    ENDDO
    ! write(6,*) ILOWERLEFT(1),ILOWERLEFT(2),ILOWERLEFT(3)
    CALL putongrid(lowerleft)
    ! write(6,*) LOWERLEFT(1),LOWERLEFT(2),LOWERLEFT(3)
    ! stop
    IF (paral%parent) THEN
       filtr='RWFN_'
       IF (paral%io_parent)&
            WRITE(numb,'(I8)') infi
       CALL xstring(numb,ia,ie)
       filbod=filtr(1:5)//numb(ia:ie)
       IF (glocal%tg_complex) THEN
          CALL xstring(filbod,i1,i2)
          filtr = filbod(i1:i2)
          IF (paral%io_parent)&
               WRITE(numb,'(I8)') ikind
          CALL xstring(numb,ia,ie)
          CALL xstring(filtr,i1,i2)
          filbod=filtr(i1:i2)//'_k'//numb(ia:ie)
       ENDIF
    ENDIF

    DO istate=1,nstate
       IF (orbpri(istate).NE.0) THEN
          IF (paral%parent) THEN

             IF (paral%io_parent)&
                  WRITE(numb,'(I8)') istate
             CALL xstring(numb,ia,ie)
             CALL xstring(filbod,i1,i2)
             filen=filbod(i1:i2)//'.'//numb(ia:ie)

             ! AK 2005/04/08. FIXME: this may be offloaded to the cubefile routine
             ! in util_p.F
             IF (paral%io_parent)&
                  CALL fileopen(2604,filen,fo_def,ferror)
             IF (paral%io_parent)&
                  WRITE (2604,*) 'CPMD CUBE FILE.'
             IF (paral%io_parent)&
                  WRITE (2604,*) 'OUTER : X, MIDDLE : Y, INNER : Z'
             IF (paral%io_parent)&
                  WRITE (2604,200) ions1%nat, (lowerleft(ir),ir=1,3)
             IF (paral%io_parent)&
                  WRITE (2604,200) spar%nr1s,&
                  (parm%a1(ir)/REAL(spar%nr1s,kind=real_8),ir=1,3)
             IF (paral%io_parent)&
                  WRITE (2604,200) spar%nr2s,&
                  (parm%a2(ir)/REAL(spar%nr2s,kind=real_8),ir=1,3)
             IF (paral%io_parent)&
                  WRITE (2604,200) spar%nr3s,&
                  (parm%a3(ir)/REAL(spar%nr3s,kind=real_8),ir=1,3)

             DO isp=1,ions1%nsp
                DO iat=1,ions0%na(isp)
                   IF (paral%io_parent)&
                        WRITE (2604,201) ions0%iatyp(isp),0.0_real_8,&
                        (r0(ir,iat,isp),ir=1,3)
                ENDDO
             ENDDO

          ENDIF

          CALL zeroing(psi)!,maxfft)

          DO ig=1,ncpw%ngw
             psi(nzhs(ig))  = c0(ig,istate)
             psi(indzs(ig)) = c0(ig+ncpw%ngw,istate)
          ENDDO
          CALL  invfftn(psi,.TRUE.,parai%allgrp)
          CALL phase(psi)

          DO ir=1,fpar%nnr1
             array(ir)  = REAL(psi(ir))*REAL(psi(ir))+&
                  AIMAG(psi(ir))*AIMAG(psi(ir))
             array(ir)  = SQRT(array(ir))
          ENDDO

          DO i1=ilowerleft(1),spar%nr1s+ilowerleft(1)-1
             CALL mp_sync(parai%allgrp)! to avoid long transfer queues
             ii1 = i1
             IF (ii1 .GT. spar%nr1s) ii1 = ii1 - spar%nr1s
             ! find the PE who has the plane i1:
             source_pe = 0
             DO WHILE (.NOT.&
                  (ii1.GE.parap%nrxpl(source_pe,1) .AND.&
                  ii1.LE.parap%nrxpl(source_pe,2)))
                source_pe = source_pe + 1
             ENDDO
             source_pe = parap%pgroup(source_pe+1)

             IF (parai%mepos .EQ. source_pe) THEN
                ! &              ARRAY(ii1-NRXPL(MEPOS,1)+1,1,1),KR1,
                CALL dcopy(fpar%kr2s*fpar%kr3s,&
                     array(ii1-parap%nrxpl(parai%mepos,1)+1),fpar%kr1,&
                     dest,1)
                IF (.NOT. paral%parent) THEN
                   !msglen = 8*kr2s*kr3s! one X-plane.
                   msgid = 2! MPI message tag.
                   CALL mp_send(dest,fpar%kr2s*fpar%kr3s,parap%pgroup(1),msgid,parai%allgrp)
                ENDIF
             ELSEIF (paral%parent) THEN
                !msglen = 8*kr2s*kr3s! one X-plane.
                msgid = 2! MPI message tag.
                CALL mp_recv(dest,fpar%kr2s*fpar%kr3s,source_pe,msgid,parai%allgrp)
             ENDIF


             ! Now, the parent should have the plane in PSI.
             IF (paral%parent) THEN
                DO i2=ilowerleft(2),spar%nr2s+ilowerleft(2) - 1
                   ii2 = i2
                   IF (ii2 .GT. spar%nr2s) ii2 = ii2 -spar%nr2s
                   DO i3=ilowerleft(3),spar%nr3s+ilowerleft(3) - 1 -1,2
                      ii3  = i3
                      ii3p = i3 + 1
                      IF (ii3 .GT. spar%nr3s) ii3 = ii3 -spar%nr3s
                      IF (ii3p .GT. spar%nr3s) ii3p = ii3p -spar%nr3s
                      IF (paral%io_parent)&
                           WRITE (2604,202) dest(ii2,ii3),dest(ii2,ii3p)
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
          IF ((paral%parent).AND.paral%io_parent)&
               CALL fileclose(2604)
       ENDIF
    ENDDO


    ! ************************************************************************

    ! CUBECENTER(1) =  (0.5_real_8-XS)*(a1(1)+a2(1)+a3(1))
    ! CUBECENTER(2) =  (0.5_real_8-YS)*(a1(2)+a2(2)+a3(2))
    ! CUBECENTER(3) =  (0.5_real_8-ZS)*(a1(3)+a2(3)+a3(3))


    cubecenter(1) = 0.0_real_8 ! (0.5_real_8-XS)*(a1(1)+a2(1)+a3(1))
    cubecenter(2) = 0.0_real_8 ! (0.5_real_8-YS)*(a1(2)+a2(2)+a3(2))
    cubecenter(3) = 0.0_real_8 ! 

    lowerleft(1) = cubecenter(1)&
         - (0.5_real_8-xs)*(parm%a1(1)+parm%a2(1)+parm%a3(1))
    lowerleft(2) = cubecenter(2)&
         - (0.5_real_8-ys)*(parm%a1(2)+parm%a2(2)+parm%a3(2))
    lowerleft(3) = cubecenter(3)&
         - (0.5_real_8-zs)*(parm%a1(3)+parm%a2(3)+parm%a3(3))

    DO ir = 1,3
       ilowerleft(ir) = rgrid(lowerleft,ir)
    ENDDO

    CALL putongrid(lowerleft)

    IF (infi .NE. 1) THEN
       IF (paral%parent) THEN

          IF (paral%io_parent)&
               WRITE(numb,'(A3)') 'TOT'
          CALL xstring(numb,ia,ie)
          CALL xstring(filbod,i1,i2)
          filen=filbod(i1:i2)//'.'//numb(ia:ie)

          IF (paral%io_parent)&
               CALL fileopen(2604,filen,fo_def,ferror)
          IF (paral%io_parent)&
               WRITE (2604,*) 'CPMD CUBE FILE.'
          IF (paral%io_parent)&
               WRITE (2604,*) 'OUTER : X, MIDDLE : Y, INNER : Z'
          IF (paral%io_parent)&
               WRITE (2604,200) ions1%nat, (lowerleft(ir),ir=1,3)
          IF (paral%io_parent)&
               WRITE (2604,200) spar%nr1s,&
               (parm%a1(ir)/REAL(spar%nr1s,kind=real_8),ir=1,3)
          IF (paral%io_parent)&
               WRITE (2604,200) spar%nr2s,&
               (parm%a2(ir)/REAL(spar%nr2s,kind=real_8),ir=1,3)
          IF (paral%io_parent)&
               WRITE (2604,200) spar%nr3s,&
               (parm%a3(ir)/REAL(spar%nr3s,kind=real_8),ir=1,3)

          DO isp=1,ions1%nsp
             DO iat=1,ions0%na(isp)
                IF (paral%io_parent)&
                     WRITE (2604,201) ions0%iatyp(isp),0.0_real_8,&
                     (r0(ir,iat,isp),ir=1,3)
             ENDDO
          ENDDO

       ENDIF

       CALL zeroing(psi)!,maxfft)

       DO icont = 1,nstate_loc

          IF (lostate%state_all) THEN
             istate = icont
          ELSE
             istate = ind_st(icont)
          ENDIF

          DO ig=1,ncpw%ngw
             ! psi(nzhs(ig))  = psi(nzhs(ig))  + C0(ig,ISTATE)
             ! psi(indzs(ig)) = psi(indzs(ig)) + C0(ig+NGW,ISTATE)
             psi(nzhs(ig))  = psi(nzhs(ig))  +&
                  SQRT(REAL(c0(ig,istate))*REAL(c0(ig,istate))+&
                  AIMAG(c0(ig,istate))*AIMAG(c0(ig,istate)))
             psi(indzs(ig)) = psi(indzs(ig)) +&
                  SQRT(REAL(c0(ig+ncpw%ngw,istate))*REAL(c0(ig+ncpw%ngw,istate))+&
                  AIMAG(c0(ig+ncpw%ngw,istate))*AIMAG(c0(ig+ncpw%ngw,istate)))
          ENDDO
       ENDDO
       CALL  invfftn(psi,.TRUE.,parai%allgrp)
       CALL phase(psi)

       DO ir=1,fpar%nnr1
          array(ir)  = REAL(psi(ir))*REAL(psi(ir))+&
               AIMAG(psi(ir))*AIMAG(psi(ir))
          array(ir)  = SQRT(array(ir))
       ENDDO

       ! PLOT DECAY FROM (000) --->  (1/21/2 1/2)

       IF (paral%io_parent)&
            CALL fileopen(1215,'tot_wan_decay_111.dat',fo_def,ferror)
       DO i1 = 1,spar%nr1s/2+1
          i2 = i1
          i3 = i1
          rrx = parm%a1(1)*i1/REAL(spar%nr1s,kind=real_8)+ parm%a2(1)*i2/REAL(spar%nr2s,kind=real_8) +&
               parm%a3(1)*i3/REAL(spar%nr3s,kind=real_8)
          rry = parm%a1(2)*i1/REAL(spar%nr1s,kind=real_8)+ parm%a2(2)*i2/REAL(spar%nr2s,kind=real_8) +&
               parm%a3(2)*i3/REAL(spar%nr3s,kind=real_8)
          rrz = parm%a1(3)*i1/REAL(spar%nr1s,kind=real_8)+ parm%a2(3)*i2/REAL(spar%nr2s,kind=real_8) +&
               parm%a3(3)*i3/REAL(spar%nr3s,kind=real_8)

          rrr = SQRT(rrx*rrx+rry*rry+rrz*rrz)

          IF (paral%io_parent)&
               WRITE(1215,'(2(1PE12.6))')&
               rrr,array(i1+(i2-1)*spar%nr1s+(i3-1)*spar%nr1s*spar%nr2s)
       ENDDO
       IF (paral%io_parent)&
            CALL fileclose(1215)


       DO i1=ilowerleft(1),spar%nr1s+ilowerleft(1)-1
          CALL mp_sync(parai%allgrp)! to avoid long transfer queues
          ii1 = i1
          IF (ii1 .GT. spar%nr1s) ii1 = ii1 - spar%nr1s
          ! find the PE who has the plane i1:
          source_pe = 0
          DO WHILE (.NOT.&
               (ii1.GE.parap%nrxpl(source_pe,1) .AND.&
               ii1.LE.parap%nrxpl(source_pe,2)))
             source_pe = source_pe + 1
          ENDDO
          source_pe = parap%pgroup(source_pe+1)

          IF (parai%mepos .EQ. source_pe) THEN
             ! &              ARRAY(ii1-NRXPL(MEPOS,1)+1,1,1),KR1,
             CALL dcopy(fpar%kr2s*fpar%kr3s,&
                  array(ii1-parap%nrxpl(parai%mepos,1)+1),fpar%kr1,&
                  dest,1)
             IF (.NOT. paral%parent) THEN
                !msglen = 8*kr2s*kr3s! one X-plane.
                msgid = 2! MPI message tag.
                CALL mp_send(dest,fpar%kr2s*fpar%kr3s,parap%pgroup(1),msgid,parai%allgrp)
             ENDIF
          ELSEIF (paral%parent) THEN
             !msglen = 8*kr2s*kr3s! one X-plane.
             msgid = 2    ! MPI message tag.
             CALL mp_recv(dest,fpar%kr2s*fpar%kr3s,source_pe,msgid,parai%allgrp)
          ENDIF


          ! Now, the parent should have the plane in PSI.
          IF (paral%parent) THEN
             DO i2=ilowerleft(2),spar%nr2s+ilowerleft(2) - 1
                ii2 = i2
                IF (ii2 .GT. spar%nr2s) ii2 = ii2 -spar%nr2s
                DO i3=ilowerleft(3),spar%nr3s+ilowerleft(3) - 1 -1,2
                   ii3  = i3
                   ii3p = i3 + 1
                   IF (ii3 .GT. spar%nr3s) ii3 = ii3 -spar%nr3s
                   IF (ii3p .GT. spar%nr3s) ii3p = ii3p -spar%nr3s
                   IF (paral%io_parent)&
                        WRITE (2604,202) dest(ii2,ii3),dest(ii2,ii3p)
                ENDDO
             ENDDO
          ENDIF
       ENDDO
       IF ((paral%parent).AND.paral%io_parent)&
            CALL fileclose(2604)
    ENDIF


    ! ==--------------------------------------------------------------==
    DEALLOCATE(r0,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    DEALLOCATE(array,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(dest,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
200 FORMAT (5x,i4,6x,3f13.5)
201 FORMAT (5x,i4,f6.2,3f13.8)
202 FORMAT (2x,2f9.5)

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE p_real_space

  ! ==================================================================


END MODULE g_loc_realspace_utils
