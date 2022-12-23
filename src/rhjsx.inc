      Ax = 0.738558766382022336_real_8
      kf = 3.09366772628013553_real_8*rho**(1.0_real_8/3.0_real_8)
      y0 = kf*r0
      s = 0.5_real_8*SQRT(grho)/(kf*rho)
      zeta = s**2*(a2*s**2 + a3*s**3 + a4*s**4 + a5*s**5 + a6*s**6 + a7*&
     & s**7)/(b1*s + b2*s**2 + b3*s**3 + b4*s**4 + b5*s**5 + b6*s**6 +  &
     & b7*s**7 + b8*s**8 + b9*s**9 + 1.0_real_8)
      lambada = d + zeta
      eta = a + zeta
      F = 1.0_real_8 - 0.037037037037037037_real_8*s**2/(c*(s**2/s0**2 + 1.0_real_8)) -&
     & 0.5_real_8*zeta/c
      EG = -0.4_real_8*F*c*lambada - 1.2_real_8*a*lambada**3 -                    &
     & 0.266666666666666667_real_8*b*lambada**2 - 2.4_real_8*lambada**(7.0_real_8/2.0_real_8&
     & )*(-SQRT(eta) + SQRT(zeta)) - 1.41796308072441279_real_8*lambada**(   &
     & 7.0_real_8/2.0_real_8)
      Fx = a*DEI(-lambada*y0**2) - eta*DEI(-eta*y0**2) + zeta*DEI(-y0**2&
     & *zeta) + (-0.444444444444444444_real_8*EG*(y0**4/lambada + 2.0_real_8*y0**2&
     & /lambada**2 + 2.0_real_8/lambada**3) - 0.444444444444444444_real_8*F*c*(   &
     & lambada*y0**2 + 1.0_real_8)/lambada**2 - 0.444444444444444444_real_8*b/    &
     & lambada)*EXP(-lambada*y0**2) + (EXP(-y0**2*zeta) - EXP(-eta*y0**2&
     & ))/y0**2
      Ener = -Ax*Fx*rho**(4.0_real_8/3.0_real_8)
      dkfdr = 1.0312225754267117_real_8/rho**(2.0_real_8/3.0_real_8)
      dy0dr = dkfdr*r0
      dsdr = -0.5_real_8*dkfdr*SQRT(grho)/(kf**2*rho) - 0.5_real_8*SQRT(grho)/(kf*&
     & rho**2)
      dzetadr = 2.0_real_8*dsdr*s*(a2*s**2 + a3*s**3 + a4*s**4 + a5*s**5 + a6&
     & *s**6 + a7*s**7)/(b1*s + b2*s**2 + b3*s**3 + b4*s**4 + b5*s**5 + &
     & b6*s**6 + b7*s**7 + b8*s**8 + b9*s**9 + 1.0_real_8) + s**2*(a2*s**2 + &
     & a3*s**3 + a4*s**4 + a5*s**5 + a6*s**6 + a7*s**7)*(-b1*dsdr -     &
     & 2.0_real_8*b2*dsdr*s - 3.0_real_8*b3*dsdr*s**2 - 4.0_real_8*b4*dsdr*s**3 - 5.0_real_8&
     & *b5*dsdr*s**4 - 6.0_real_8*b6*dsdr*s**5 - 7.0_real_8*b7*dsdr*s**6 - 8.0_real_8* &
     & b8*dsdr*s**7 - 9.0_real_8*b9*dsdr*s**8)/(b1*s + b2*s**2 + b3*s**3 + b4&
     & *s**4 + b5*s**5 + b6*s**6 + b7*s**7 + b8*s**8 + b9*s**9 + 1.0_real_8) &
     & **2 + s**2*(2.0_real_8*a2*dsdr*s + 3.0_real_8*a3*dsdr*s**2 + 4.0_real_8*a4*dsdr*&
     & s**3 + 5.0_real_8*a5*dsdr*s**4 + 6.0_real_8*a6*dsdr*s**5 + 7.0_real_8*a7*dsdr*s &
     & **6)/(b1*s + b2*s**2 + b3*s**3 + b4*s**4 + b5*s**5 + b6*s**6 + b7&
     & *s**7 + b8*s**8 + b9*s**9 + 1.0_real_8)
      dlambadadr = dzetadr
      detadr = dzetadr
      dFdr = 0.0740740740740740741_real_8*dsdr*s**3/(c*s0**2*(s**2/s0**2 +   &
     & 1.0_real_8)**2) - 0.0740740740740740741_real_8*dsdr*s/(c*(s**2/s0**2 +     &
     & 1.0_real_8)) - 0.5_real_8*dzetadr/c
      dEGdr = -0.4_real_8*F*c*dlambadadr - 3.6_real_8*a*dlambadadr*lambada**2 -   &
     & 0.533333333333333333_real_8*b*dlambadadr*lambada - 0.4_real_8*c*dFdr*      &
     & lambada - 8.4_real_8*dlambadadr*lambada**(5.0_real_8/2.0_real_8)*(-SQRT(eta) +  &
     & SQRT(zeta)) - 4.962870782535445_real_8*dlambadadr*lambada**(5.0_real_8/    &
     & 2.0_real_8) - 2.4_real_8*lambada**(7.0_real_8/2.0_real_8)*(-0.5_real_8*detadr/SQRT(eta) + &
     & 0.5_real_8*dzetadr/SQRT(zeta))
      dFxdr = -a*(-dlambadadr*y0**2 - 2.0_real_8*dy0dr*lambada*y0)*EXP(      &
     & -lambada*y0**2)/(lambada*y0**2) - detadr*DEI(-eta*y0**2) - 2.0_real_8*&
     & dy0dr*(EXP(-y0**2*zeta) - EXP(-eta*y0**2))/y0**3 + dzetadr*DEI(  &
     & -y0**2*zeta) + (-dlambadadr*y0**2 - 2.0_real_8*dy0dr*lambada*y0)*(    &
     & -0.444444444444444444_real_8*EG*(y0**4/lambada + 2.0_real_8*y0**2/lambada**&
     & 2 + 2.0_real_8/lambada**3) - 0.444444444444444444_real_8*F*c*(lambada*y0**2&
     & + 1.0_real_8)/lambada**2 - 0.444444444444444444_real_8*b/lambada)*EXP(     &
     & -lambada*y0**2) + (-0.444444444444444444_real_8*EG*(-dlambadadr*y0**4/&
     & lambada**2 - 4.0_real_8*dlambadadr*y0**2/lambada**3 - 6.0_real_8*dlambadadr&
     & /lambada**4 + 4.0_real_8*dy0dr*y0**3/lambada + 4.0_real_8*dy0dr*y0/lambada &
     & **2) + 0.888888888888888889_real_8*F*c*dlambadadr*(lambada*y0**2 +    &
     & 1.0_real_8)/lambada**3 - 0.444444444444444444_real_8*F*c*(dlambadadr*y0**2 &
     & + 2.0_real_8*dy0dr*lambada*y0)/lambada**2 + 0.444444444444444444_real_8*b* &
     & dlambadadr/lambada**2 - 0.444444444444444444_real_8*c*dFdr*(lambada*y0&
     & **2 + 1.0_real_8)/lambada**2 - 0.444444444444444444_real_8*dEGdr*(y0**4/   &
     & lambada + 2.0_real_8*y0**2/lambada**2 + 2.0_real_8/lambada**3))*EXP(       &
     & -lambada*y0**2) + (-detadr*y0**2 - 2.0_real_8*dy0dr*eta*y0)*EXP(-eta* &
     & y0**2)/y0**2 + (-(-detadr*y0**2 - 2.0_real_8*dy0dr*eta*y0)*EXP(-eta*y0&
     & **2) + (-2.0_real_8*dy0dr*y0*zeta - dzetadr*y0**2)*EXP(-y0**2*zeta))/ &
     & y0**2 - (-2.0_real_8*dy0dr*y0*zeta - dzetadr*y0**2)*EXP(-y0**2*zeta)/ &
     & y0**2
      dEnerdr = -1.33333333333333333_real_8*Ax*Fx*rho**(1.0_real_8/3.0_real_8) - Ax*   &
     & dFxdr*rho**(4.0_real_8/3.0_real_8)
      dsdn = 0.25_real_8/(SQRT(grho)*kf*rho)
      dzetadn = 2.0_real_8*dsdn*s*(a2*s**2 + a3*s**3 + a4*s**4 + a5*s**5 + a6&
     & *s**6 + a7*s**7)/(b1*s + b2*s**2 + b3*s**3 + b4*s**4 + b5*s**5 + &
     & b6*s**6 + b7*s**7 + b8*s**8 + b9*s**9 + 1.0_real_8) + s**2*(a2*s**2 + &
     & a3*s**3 + a4*s**4 + a5*s**5 + a6*s**6 + a7*s**7)*(-b1*dsdn -     &
     & 2.0_real_8*b2*dsdn*s - 3.0_real_8*b3*dsdn*s**2 - 4.0_real_8*b4*dsdn*s**3 - 5.0_real_8&
     & *b5*dsdn*s**4 - 6.0_real_8*b6*dsdn*s**5 - 7.0_real_8*b7*dsdn*s**6 - 8.0_real_8* &
     & b8*dsdn*s**7 - 9.0_real_8*b9*dsdn*s**8)/(b1*s + b2*s**2 + b3*s**3 + b4&
     & *s**4 + b5*s**5 + b6*s**6 + b7*s**7 + b8*s**8 + b9*s**9 + 1.0_real_8) &
     & **2 + s**2*(2.0_real_8*a2*dsdn*s + 3.0_real_8*a3*dsdn*s**2 + 4.0_real_8*a4*dsdn*&
     & s**3 + 5.0_real_8*a5*dsdn*s**4 + 6.0_real_8*a6*dsdn*s**5 + 7.0_real_8*a7*dsdn*s &
     & **6)/(b1*s + b2*s**2 + b3*s**3 + b4*s**4 + b5*s**5 + b6*s**6 + b7&
     & *s**7 + b8*s**8 + b9*s**9 + 1.0_real_8)
      dlambadadn = dzetadn
      detadn = dzetadn
      dFdn = 0.0740740740740740741_real_8*dsdn*s**3/(c*s0**2*(s**2/s0**2 +   &
     & 1.0_real_8)**2) - 0.0740740740740740741_real_8*dsdn*s/(c*(s**2/s0**2 +     &
     & 1.0_real_8)) - 0.5_real_8*dzetadn/c
      dEGdn = -0.4_real_8*F*c*dlambadadn - 3.6_real_8*a*dlambadadn*lambada**2 -   &
     & 0.533333333333333333_real_8*b*dlambadadn*lambada - 0.4_real_8*c*dFdn*      &
     & lambada - 8.4_real_8*dlambadadn*lambada**(5.0_real_8/2.0_real_8)*(-SQRT(eta) +  &
     & SQRT(zeta)) - 4.962870782535445_real_8*dlambadadn*lambada**(5.0_real_8/    &
     & 2.0_real_8) - 2.4_real_8*lambada**(7.0_real_8/2.0_real_8)*(-0.5_real_8*detadn/SQRT(eta) + &
     & 0.5_real_8*dzetadn/SQRT(zeta))
      dFxdn = a*dlambadadn*EXP(-lambada*y0**2)/lambada - detadn*DEI(-eta&
     & *y0**2) - detadn*EXP(-eta*y0**2) - dlambadadn*y0**2*(            &
     & -0.444444444444444444_real_8*EG*(y0**4/lambada + 2.0_real_8*y0**2/lambada**&
     & 2 + 2.0_real_8/lambada**3) - 0.444444444444444444_real_8*F*c*(lambada*y0**2&
     & + 1.0_real_8)/lambada**2 - 0.444444444444444444_real_8*b/lambada)*EXP(     &
     & -lambada*y0**2) + dzetadn*DEI(-y0**2*zeta) + dzetadn*EXP(-y0**2* &
     & zeta) + (-0.444444444444444444_real_8*EG*(-dlambadadn*y0**4/lambada**2&
     & - 4.0_real_8*dlambadadn*y0**2/lambada**3 - 6.0_real_8*dlambadadn/lambada**4&
     & ) - 0.444444444444444444_real_8*F*c*dlambadadn*y0**2/lambada**2 +     &
     & 0.888888888888888889_real_8*F*c*dlambadadn*(lambada*y0**2 + 1.0_real_8)/   &
     & lambada**3 + 0.444444444444444444_real_8*b*dlambadadn/lambada**2 -    &
     & 0.444444444444444444_real_8*c*dFdn*(lambada*y0**2 + 1.0_real_8)/lambada**2 &
     & - 0.444444444444444444_real_8*dEGdn*(y0**4/lambada + 2.0_real_8*y0**2/     &
     & lambada**2 + 2.0_real_8/lambada**3))*EXP(-lambada*y0**2) + (detadn*y0 &
     & **2*EXP(-eta*y0**2) - dzetadn*y0**2*EXP(-y0**2*zeta))/y0**2
      dEnerdn = -Ax*dFxdn*rho**(4.0_real_8/3.0_real_8)
