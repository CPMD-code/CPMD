MODULE ffsum_utils
  USE cppt,                            ONLY: hg,&
                                             inyh
  USE fitpack_utils,                   ONLY: curvd
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE qspl,                            ONLY: ggnh,&
                                             nsplpo,&
                                             voo
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb
  USE system,                          ONLY: cntl,&
                                             ncpw,&
                                             parm

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ffsum

CONTAINS

  ! ==================================================================
  SUBROUTINE ffsum(dsumps)
    ! ==--------------------------------------------------------------==
    ! == Calculate DSUMPS(1:NHG)                                      ==
    ! ==         = derivative of local pseudopotential part           ==
    ! ==           * structure factor                                 ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: dsumps(ncpw%nhg)

    COMPLEX(real_8)                          :: sfacs
    INTEGER                                  :: ia, ig, is, isa, isa0
    REAL(real_8)                             :: dsuml, vol

    vol=1._real_8/parm%omega
    IF (cntl%bigmem) THEN
       !$omp parallel do private(IG,IS,IA,ISA0,ISA,SFACS,DSUML)
#ifdef __SR8000
       !poption parallel, tlocal(IG,IS,IA,ISA0,ISA,SFACS,DSUML)
#endif 
       DO ig=1,ncpw%nhg
          isa0=0
          DO is=1,ions1%nsp
             sfacs=(0.0_real_8,0.0_real_8)
             dsuml=-vol*curvd(hg(ig),nsplpo,ggnh(1),voo(1,1,is),&
                  voo(1,2,is),0.0_real_8)/parm%tpiba2*2._real_8
             DO ia=1,ions0%na(is)
                isa=isa0+ia
                sfacs=sfacs+eigrb(ig,isa)
             ENDDO
             dsumps(ig) = dsumps(ig) + sfacs*dsuml
             isa0=isa0+ions0%na(is)
          ENDDO
       ENDDO
    ELSE
       !$omp parallel do private(IG,IS,IA,ISA0,ISA,SFACS,DSUML)
#ifdef __SR8000
       !poption parallel, tlocal(IG,IS,IA,ISA0,ISA,SFACS,DSUML)
#endif 
       DO ig=1,ncpw%nhg
          isa0=0
          DO is=1,ions1%nsp
             sfacs=(0.0_real_8,0.0_real_8)
             dsuml=-vol*curvd(hg(ig),nsplpo,ggnh(1),voo(1,1,is),&
                  voo(1,2,is),0.0_real_8)/parm%tpiba2*2._real_8
             DO ia=1,ions0%na(is)
                isa=isa0+ia
                sfacs=sfacs+ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                     ei3(isa,inyh(3,ig))
             ENDDO
             dsumps(ig) = dsumps(ig) + sfacs*dsuml
             isa0=isa0+ions0%na(is)
          ENDDO
       ENDDO
    ENDIF
    IF (geq0) dsumps(1)=0._real_8
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ffsum
  ! ==================================================================

END MODULE ffsum_utils
