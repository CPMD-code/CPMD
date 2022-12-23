MODULE ldosmod
  IMPLICIT NONE




  TYPE :: cldos_t
     INTEGER :: nlayer
     INTEGER :: ldosn
     LOGICAL :: tldos
  END TYPE cldos_t
  TYPE(cldos_t) :: cldos

END MODULE ldosmod
