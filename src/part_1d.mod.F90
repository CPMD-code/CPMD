MODULE part_1d
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_dims_create

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: part_1d_nbr_el_in_blk
  PUBLIC :: part_1d_get_el_in_blk
  PUBLIC :: part_1d_nbr_elems
  PUBLIC :: part_1d_get_elem
  PUBLIC :: part_1d_elem_to_proc
  PUBLIC :: part_1d_symm_holds_pair
  PUBLIC :: part_1d_get_blk_bounds
  PUBLIC :: proc_to_grid2d

CONTAINS

  ! ==================================================================
  INTEGER FUNCTION part_1d_nbr_el_in_blk(n_elem,proc,nproc)
    ! ==================================================================
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
    IMPLICIT NONE
    ! Arguments
    INTEGER :: n_elem,proc,nproc
    ! Variables
    INTEGER :: nbr,res
    ! ==--------------------------------------------------------------==
    res = MOD(n_elem,nproc)
    nbr = (n_elem-res)/nproc
    IF (proc.LT.res) nbr = nbr + 1
    part_1d_nbr_el_in_blk = nbr
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION part_1d_nbr_el_in_blk

  ! ==================================================================
  INTEGER FUNCTION part_1d_get_el_in_blk(i_elem,n_elem,proc,&
       nproc)
    ! ==================================================================
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
    IMPLICIT NONE
    ! Arguments
    INTEGER :: i_elem,n_elem,proc,nproc
    ! Variables
    INTEGER :: nbr,res
    ! ==--------------------------------------------------------------==
    res = MOD(n_elem,nproc)
    nbr = (n_elem-res)/nproc
    part_1d_get_el_in_blk = i_elem + nbr * proc + MIN(proc,res)
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION part_1d_get_el_in_blk

  ! ==================================================================
  INTEGER FUNCTION part_1d_nbr_elems(n_elem,proc,nproc)
    ! ==================================================================
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
    IMPLICIT NONE
    ! Arguments
    INTEGER :: n_elem,proc,nproc
    ! Variables
    INTEGER :: nbr
    ! ==--------------------------------------------------------------==
    nbr = (n_elem-proc)/nproc
    IF (MOD(n_elem-proc,nproc).GT.0) nbr = nbr + 1
    part_1d_nbr_elems = nbr
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION part_1d_nbr_elems

  ! ==================================================================
  INTEGER FUNCTION part_1d_get_elem(i_elem,proc,nproc)
    ! ==================================================================
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
    IMPLICIT NONE
    ! Arguments
    INTEGER :: i_elem,proc,nproc
    ! Variables
    ! ==--------------------------------------------------------------==
    part_1d_get_elem = (i_elem - 1) * nproc + proc + 1
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION part_1d_get_elem

  ! ==================================================================
  SUBROUTINE part_1d_get_blk_bounds(n_elem,proc,nproc,first,last)
    ! ==================================================================
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n_elem, proc, nproc, first, &
                                                last

    INTEGER                                  :: iproc, m, old

! Variables
! ==--------------------------------------------------------------==
! can do better !

    m = n_elem
    old = 0
    DO iproc = 0,proc
       first = old + 1
       last = old + CEILING( REAL( m ,kind=real_8) / REAL( nproc - iproc ,kind=real_8) )
       m = m - ( last - first + 1 )
       old = last
    ENDDO
    IF (first.GT.n_elem) THEN
       first = 1
       last = 0
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE part_1d_get_blk_bounds

  ! ==================================================================
  INTEGER FUNCTION part_1d_elem_to_proc(i_elem,nproc)
    ! ==================================================================
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
    IMPLICIT NONE
    ! Arguments
    INTEGER :: i_elem, nproc
    ! ==--------------------------------------------------------------==
    part_1d_elem_to_proc = MOD(i_elem-1, nproc)
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION part_1d_elem_to_proc

  ! ==================================================================
  LOGICAL FUNCTION part_1d_symm_holds_pair(i_elem,j_elem)
    ! ==================================================================
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
    IMPLICIT NONE
    ! Arguments
    INTEGER :: i_elem,j_elem
    ! ==--------------------------------------------------------------==
    part_1d_symm_holds_pair = BTEST(i_elem+j_elem,0).EQV.&
         i_elem.GE.j_elem
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION part_1d_symm_holds_pair

  ! ==================================================================
  SUBROUTINE proc_to_grid2d(nproc,me,nprows,npcols,prow,pcol)
    ! ==================================================================
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nproc, me, nprows, npcols, &
                                                prow, pcol

    INTEGER                                  :: dims(2)

    dims = 0
    CALL mp_dims_create(nproc,2,dims)
    nprows = dims(1)
    npcols = dims(2)
    pcol = me / nprows
    prow = me - nprows * pcol
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE proc_to_grid2d

END MODULE part_1d
