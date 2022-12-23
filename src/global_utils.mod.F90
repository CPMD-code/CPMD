MODULE global_utils
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: paral
  USE pimd,                            ONLY: num_np,&
                                             parentgroup,&
                                             pimd3

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: global

CONTAINS

  ! ==================================================================
  SUBROUTINE global(array,lda)
    ! ==--------------------------------------------------------------==
    ! ==  Globalizes data in Array to all parents                     ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lda
    REAL(real_8)                             :: array(lda,pimd3%np_total)

    INTEGER                                  :: ip, iparent, isrc

    IF (paral%io_parent) THEN
       DO ip=1,pimd3%np_total
          iparent=num_np(ip)
          ! DM1          ISRC=IP-1
          isrc=iparent-1
          ! DM2
          CALL mp_bcast(array(:,ip),lda,isrc,parentgroup)
       ENDDO
    ENDIF

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE global
  ! ==================================================================

END MODULE global_utils
