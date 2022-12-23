MODULE write_pp_utils
  USE dpot,                            ONLY: dpot_mod
  USE fileopen_utils,                  ONLY: fileclose,&
                                             fileopen
  USE fileopenmod,                     ONLY: fo_def,&
                                             fo_verb
  USE func,                            ONLY: func1,&
                                             func2
  USE ions,                            ONLY: ions0
  USE kinds,                           ONLY: real_8
  USE parac,                           ONLY: paral
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: write_pp

CONTAINS

  SUBROUTINE write_pp(isp,ecpname)
    INTEGER                                  :: isp
    CHARACTER(len=*)                         :: ecpname

    CHARACTER(len=1), DIMENSION(5), PARAMETER :: &
      lblsp = (/'s','p','d','f','g'/)
    CHARACTER(len=40), DIMENSION(0:4), PARAMETER :: cfmtstr = (/&
      '(1X,I3,"   #C" )                        ',&
      '(1X,I3,F16.9,"   #C  C1")               ',&
      '(1X,I3,2(F16.9),"   #C  C1 C2")         ',&
      '(1X,I3,3(F16.9),"   #C  C1 C2 C3")      ',&
      '(1X,I3,4(F16.9),"   #C  C1 C2 C3 C4")   '/)
    CHARACTER(len=52), DIMENSION(0:3), PARAMETER :: hfmtstr = (/&
      '(1X,F16.9,I3," H(",A,") no projector")           ',&
      '(1X,F16.9,I3,1(F14.9)," H(",A,") 11")            ',&
      '(1X,F16.9,I3,3(F14.9)," H(",A,") 1112 22")       ',&
      '(1X,F16.9,I3,6(F14.9)," H(",A,") 1112 1322 2333")'/)

    INTEGER                                  :: excf, i, j, k, l
    LOGICAL                                  :: foerr

! Generate four digit code for functional

    excf=func1%mfxcx*1000+func1%mfxcc*100+func1%mgcx*10+func1%mgcc

    IF (paral%io_parent)&
         CALL fileopen(94,ecpname,fo_def+fo_verb,foerr)
    IF (paral%io_parent)&
         WRITE(94,*) '&ATOM'
    IF (paral%io_parent)&
         WRITE(94,*) ' Z =',REAL(ions0%iatyp(isp),kind=real_8)
    IF (paral%io_parent)&
         WRITE(94,*) ' ZV =',ions0%zv(isp)
    IF (paral%io_parent)&
         WRITE(94,'(A,I4,F15.10)') '  XC = ',excf,func2%salpha
    IF (paral%io_parent)&
         WRITE(94,*) ' TYPE = NORMCONSERVING GOEDECKER'
    IF (paral%io_parent)&
         WRITE(94,*) '&END'
    IF (paral%io_parent)&
         WRITE(94,*) '&INFO'
    IF (paral%io_parent)&
         WRITE(94,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
    IF (paral%io_parent)&
         WRITE(94,*) '>    OACP pseudopotential (Goedecker type)     >'
    IF (paral%io_parent)&
         WRITE(94,'(A,4X,A,F6.2,4X,A,F6.2,21X,A)')&
         ' >','Z=',REAL(ions0%iatyp(isp),kind=real_8),'ZV=',ions0%zv(isp),'>'
    IF (paral%io_parent)&
         WRITE(94,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
    IF (paral%io_parent)&
         WRITE(94,*) '&END'
    IF (paral%io_parent)&
         WRITE(94,*) '&POTENTIAL'
    IF (paral%io_parent)&
         WRITE(94,*) '    GOEDECKER'
    IF (paral%io_parent)&
         WRITE(94,'(1X,I3,46X,A)') dpot_mod%lmax(isp)-1, 'LMAX'
    IF (paral%io_parent)&
         WRITE(94,'(1X,F16.9,33X,A)') sgpp2%rcsg(isp),'RC'

    l=sgpp1%nclsg(isp)
    IF (paral%io_parent)&
         WRITE(94,cfmtstr(l)) l,(sgpp2%clsg(i,isp),i=1,l)

    DO l=1,dpot_mod%lmax(isp)-1
       k=sgpp2%npro(l,isp)
       IF (paral%io_parent)&
            WRITE(94,hfmtstr(k)) sgpp2%rcnl(l,isp),k,&
            ((sgpp2%hlsg(i,j,l,isp),j=i,k),i=1,k),lblsp(l)
    ENDDO

    IF (paral%io_parent)&
         CALL fileclose(94)
    RETURN
  END SUBROUTINE write_pp

END MODULE write_pp_utils
