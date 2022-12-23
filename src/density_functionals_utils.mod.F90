MODULE density_functionals_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: pade_lda
  PUBLIC :: slater_lda

CONTAINS

  ! ==================================================================
  FUNCTION pade_lda(x,order)
    ! ==--------------------------------------------------------------==
    ! Pade approximation to LDA functional
    REAL(real_8)                             :: x
    INTEGER                                  :: order
    REAL(real_8)                             :: pade_lda

    REAL(real_8), PARAMETER :: a0 = 0.4581652932831429_real_8, &
      a1 = 2.217058676663745_real_8, a2 = 0.7405551735357053_real_8, &
      a3 = 0.01968227878617998_real_8 , b1 = 1.0_real_8, &
      b2 = 4.504130959426697_real_8, b3 = 1.110667363742916_real_8, &
      b4 = 0.02359291751427506_real_8 , c = 0.6203504908994001_real_8 

    REAL(real_8) :: t1, t10, t103, t11, t110, t111, t12, t14, t15, t16, t17, &
      t18, t19, t2, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t3, &
      t30, t31, t32, t33, t34, t35, t36, t37, t39, t4, t40, t41, t42, t43, &
      t46, t48, t49, t5, t50, t52, t58, t59, t6, t61, t62, t65, t67, t69, t7, &
      t72, t75, t8, t80, t83, t89, t9, t92, t93, t999

! ==--------------------------------------------------------------==

    IF (order.EQ.0) THEN

       t2 = x**(1._real_8/3._real_8)
       t3 = 1/t2
       t5 = c**2
       t7 = t2**2
       t8 = 1/t7
       t10 = t5*c
       t12 = 1/x
       t22 = t5**2
       t999 = -x*(a0+a1*c*t3+a2*t5*t8+a3*t10*t12)/&
            (b1*c*t3+b2*t5*t8+b3*t10*t12+b4*t22/t2/x)

    ELSEIF (order.EQ.1) THEN

       t1 = a1*c
       t2 = x**(1._real_8/3._real_8)
       t3 = 1/t2
       t5 = c**2
       t6 = a2*t5
       t7 = t2**2
       t8 = 1/t7
       t10 = t5*c
       t11 = a3*t10
       t12 = 1/x
       t14 = a0+t1*t3+t6*t8+t11*t12
       t15 = b1*c
       t17 = b2*t5
       t19 = b3*t10
       t21 = t5**2
       t22 = b4*t21
       t24 = 1/t2/x
       t26 = t15*t3+t17*t8+t19*t12+t22*t24
       t27 = 1/t26
       t32 = 1/t7/x
       t35 = x**2
       t36 = 1/t35
       t42 = t26**2
       t999 = -t14*t27-x*(-t1*t24/3-2._real_8/3._real_8*t6*t32-t11*t36)*t27+&
            x*t14/t42*(-t15*t24/3-2._real_8/3._real_8*t17*t32-t19*t36-4._real_8/&
            3._real_8*t22/t2/t35)

    ELSEIF (order.EQ.2) THEN

       t1 = a1*c
       t2 = x**(1._real_8/3._real_8)
       t4 = 1/t2/x
       t7 = c**2
       t8 = a2*t7
       t9 = t2**2
       t11 = 1/t9/x
       t14 = t7*c
       t15 = a3*t14
       t16 = x**2
       t17 = 1/t16
       t19 = -t1*t4/3-2._real_8/3._real_8*t8*t11-t15*t17
       t20 = b1*c
       t21 = 1/t2
       t23 = b2*t7
       t24 = 1/t9
       t26 = b3*t14
       t27 = 1/x
       t29 = t7**2
       t30 = b4*t29
       t32 = t20*t21+t23*t24+t26*t27+t30*t4
       t33 = 1/t32
       t39 = a0+t1*t21+t8*t24+t15*t27
       t40 = t32**2
       t41 = 1/t40
       t49 = 1/t2/t16
       t52 = -t20*t4/3-2._real_8/3._real_8*t23*t11-t26*t17-4._real_8/3._real_8*t30*t49
       t58 = 1/t9/t16
       t61 = t16*x
       t62 = 1/t61
       t72 = x*t39
       t75 = t52**2
       t999 = -2*t19*t33+2*t39*t41*t52-x*(4._real_8/9._real_8*t1*t49+10._real_8/&
            9._real_8*t8*t58+2*t15*t62)*t33+2*x*t19*t41*t52-2*t72/t40/&
            t32*t75+t72*t41*(4._real_8/9._real_8*t20*t49+10._real_8/9._real_8*t23*t58+&
            2*t26*t62+28._real_8/9._real_8*t30/t2/t61)

    ELSEIF (order.EQ.3) THEN

       t1 = a1*c
       t2 = x**2
       t3 = x**(1._real_8/3._real_8)
       t5 = 1/t3/t2
       t8 = c**2
       t9 = a2*t8
       t10 = t3**2
       t12 = 1/t10/t2
       t15 = t8*c
       t16 = a3*t15
       t17 = t2*x
       t18 = 1/t17
       t21 = 4._real_8/9._real_8*t1*t5+10._real_8/9._real_8*t9*t12+2*t16*t18
       t22 = b1*c
       t23 = 1/t3
       t25 = b2*t8
       t26 = 1/t10
       t28 = b3*t15
       t29 = 1/x
       t31 = t8**2
       t32 = b4*t31
       t34 = 1/t3/x
       t36 = t22*t23+t25*t26+t28*t29+t32*t34
       t37 = 1/t36
       t43 = 1/t10/x
       t46 = 1/t2
       t48 = -t1*t34/3-2._real_8/3._real_8*t9*t43-t16*t46
       t49 = t36**2
       t50 = 1/t49
       t59 = -t22*t34/3-2._real_8/3._real_8*t25*t43-t28*t46-4._real_8/3._real_8*t32*t5
       t65 = a0+t1*t23+t9*t26+t16*t29
       t67 = 1/t49/t36
       t69 = t59**2
       t80 = 1/t3/t17
       t83 = 4._real_8/9._real_8*t22*t5+10._real_8/9._real_8*t25*t12+2*t28*t18+28._real_8/&
            9._real_8*t32*t80
       t89 = 1/t10/t17
       t92 = t2**2
       t93 = 1/t92
       t103 = x*t48
       t110 = x*t65
       t111 = t49**2
       t999 = -3*t21*t37+6*t48*t50*t59-6*t65*t67*t69+3*t65*t50*t83-x*&
            (-28._real_8/27._real_8*t1*t80-80._real_8/27._real_8*t9*t89-6*t16*t93)*&
            t37+3*x*t21*t50*t59-6*t103*t67*t69+3*t103*t50*t83+6*&
            t110/t111*t69*t59-6*t110*t67*t59*t83+t110*t50*(-28._real_8/&
            27._real_8*t22*t80-80._real_8/27._real_8*t25*t89-6*t28*t93-&
            280._real_8/27._real_8*t32/t3/t92)

    ELSE

       CALL stopgm('PADE_LDA','ORDER .GT. 3 NOT PROGRAMMED',& 
            __LINE__,__FILE__)

    ENDIF

    pade_lda = t999

    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION pade_lda
  ! ==================================================================
  FUNCTION slater_lda(x,order)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: x
    INTEGER                                  :: order
    REAL(real_8)                             :: slater_lda

    REAL(real_8), PARAMETER :: alpha = 2._real_8/3._real_8 , &
      f1 = -1.10783814957303361_real_8 

    REAL(real_8)                             :: t2, t3, t999

! ==--------------------------------------------------------------==

    IF (order.EQ.0) THEN

       t2 = x**(1._real_8/3._real_8)
       t999 = f1*alpha*t2*x

    ELSEIF (order.EQ.1) THEN

       t2 = x**(1._real_8/3._real_8)
       t999 = 4._real_8/3._real_8*f1*alpha*t2

    ELSEIF (order.EQ.2) THEN

       t2 = x**(1._real_8/3._real_8)
       t3 = t2**2
       t999 = 4._real_8/9._real_8*f1*alpha/t3

    ELSEIF (order.EQ.3) THEN

       t2 = x**(1._real_8/3._real_8)
       t3 = t2**2
       t999 = -8._real_8/27._real_8*f1*alpha/t3/x

    ELSE

       CALL stopgm('SLATER_LDA','ORDER .GT. 3 NOT PROGRAMMED',& 
            __LINE__,__FILE__)

    ENDIF

    slater_lda = t999

    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION slater_lda
  ! ==================================================================

END MODULE density_functionals_utils
