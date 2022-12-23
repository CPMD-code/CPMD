      Ener = 0.0_real_8
      dEnerdra = 0.0_real_8
      dEnerdrb = 0.0_real_8
      dEnerdnaa = 0.0_real_8
      dEnerdnbb = 0.0_real_8
      IF(RHOA.GT.TOL) THEN
      IF(RHOB.GT.TOL) THEN
      Ax = 0.738558766382022336_real_8
      kfa = 3.89777708972075354_real_8*RHOA**(1.0_real_8/3.0_real_8)
      kfb = 3.89777708972075354_real_8*RHOB**(1.0_real_8/3.0_real_8)
      y0a = R0*kfa
      y0b = R0*kfb
      sa = 0.5_real_8*SQRT(GRHOAA)/(RHOA*kfa)
      sb = 0.5_real_8*SQRT(GRHOBB)/(RHOB*kfb)
      zetaa = sa**2*(a2*sa**2 + a3*sa**3 + a4*sa**4 + a5*sa**5 + a6*sa**&
     & 6 + a7*sa**7)/(b1*sa + b2*sa**2 + b3*sa**3 + b4*sa**4 + b5*sa**5 &
     & + b6*sa**6 + b7*sa**7 + b8*sa**8 + b9*sa**9 + 1.0_real_8)
      zetab = sb**2*(a2*sb**2 + a3*sb**3 + a4*sb**4 + a5*sb**5 + a6*sb**&
     & 6 + a7*sb**7)/(b1*sb + b2*sb**2 + b3*sb**3 + b4*sb**4 + b5*sb**5 &
     & + b6*sb**6 + b7*sb**7 + b8*sb**8 + b9*sb**9 + 1.0_real_8)
      lambadaa = D + zetaa
      lambadab = D + zetab
      etaa = A + zetaa
      etab = A + zetab
      Fa = 1.0_real_8 - 0.037037037037037037_real_8*sa**2/(C*(1.0_real_8 + sa**2/s0**2 &
     & )) - 0.5_real_8*zetaa/C
      Fb = 1.0_real_8 - 0.037037037037037037_real_8*sb**2/(C*(1.0_real_8 + sb**2/s0**2 &
     & )) - 0.5_real_8*zetab/C
      EGa = -1.2_real_8*A*lambadaa**3 - 0.266666666666666667_real_8*B*lambadaa**2 &
     & - 0.4_real_8*C*Fa*lambadaa - 2.4_real_8*lambadaa**(7.0_real_8/2.0_real_8)*(-SQRT(etaa&
     & ) + SQRT(zetaa)) - 1.41796308072441279_real_8*lambadaa**(7.0_real_8/2.0_real_8)
      EGb = -1.2_real_8*A*lambadab**3 - 0.266666666666666667_real_8*B*lambadab**2 &
     & - 0.4_real_8*C*Fb*lambadab - 2.4_real_8*lambadab**(7.0_real_8/2.0_real_8)*(-SQRT(etab&
     & ) + SQRT(zetab)) - 1.41796308072441279_real_8*lambadab**(7.0_real_8/2.0_real_8)
      Fxa = A*DEI(-lambadaa*y0a**2) - etaa*DEI(-etaa*y0a**2) + zetaa*DEI&
     & (-y0a**2*zetaa) + (-0.444444444444444444_real_8*B/lambadaa -          &
     & 0.444444444444444444_real_8*C*Fa*(lambadaa*y0a**2 + 1.0_real_8)/lambadaa**2&
     & - 0.444444444444444444_real_8*EGa*(y0a**4/lambadaa + 2.0_real_8*y0a**2/    &
     & lambadaa**2 + 2.0_real_8/lambadaa**3))*EXP(-lambadaa*y0a**2) + (EXP(  &
     & -y0a**2*zetaa) - EXP(-etaa*y0a**2))/y0a**2
      Fxb = A*DEI(-lambadab*y0b**2) - etab*DEI(-etab*y0b**2) + zetab*DEI&
     & (-y0b**2*zetab) + (-0.444444444444444444_real_8*B/lambadab -          &
     & 0.444444444444444444_real_8*C*Fb*(lambadab*y0b**2 + 1.0_real_8)/lambadab**2&
     & - 0.444444444444444444_real_8*EGb*(y0b**4/lambadab + 2.0_real_8*y0b**2/    &
     & lambadab**2 + 2.0_real_8/lambadab**3))*EXP(-lambadab*y0b**2) + (EXP(  &
     & -y0b**2*zetab) - EXP(-etab*y0b**2))/y0b**2
      Ener = -1.25992104989487316_real_8*Ax*(Fxa*RHOA**(4.0_real_8/3.0_real_8) + Fxb*  &
     & RHOB**(4.0_real_8/3.0_real_8))
      dkfadra = 1.29925902990691777_real_8/RHOA**(2.0_real_8/3.0_real_8)
      dy0adra = R0*dkfadra
      dsadra = -0.5_real_8*SQRT(GRHOAA)*dkfadra/(RHOA*kfa**2) - 0.5_real_8*SQRT(  &
     & GRHOAA)/(RHOA**2*kfa)
      dzetaadra = 2.0_real_8*dsadra*sa*(a2*sa**2 + a3*sa**3 + a4*sa**4 + a5* &
     & sa**5 + a6*sa**6 + a7*sa**7)/(b1*sa + b2*sa**2 + b3*sa**3 + b4*sa&
     & **4 + b5*sa**5 + b6*sa**6 + b7*sa**7 + b8*sa**8 + b9*sa**9 +     &
     & 1.0_real_8) + sa**2*(a2*sa**2 + a3*sa**3 + a4*sa**4 + a5*sa**5 + a6*sa&
     & **6 + a7*sa**7)*(-b1*dsadra - 2.0_real_8*b2*dsadra*sa - 3.0_real_8*b3*     &
     & dsadra*sa**2 - 4.0_real_8*b4*dsadra*sa**3 - 5.0_real_8*b5*dsadra*sa**4 -   &
     & 6.0_real_8*b6*dsadra*sa**5 - 7.0_real_8*b7*dsadra*sa**6 - 8.0_real_8*b8*dsadra* &
     & sa**7 - 9.0_real_8*b9*dsadra*sa**8)/(b1*sa + b2*sa**2 + b3*sa**3 + b4*&
     & sa**4 + b5*sa**5 + b6*sa**6 + b7*sa**7 + b8*sa**8 + b9*sa**9 +   &
     & 1.0_real_8)**2 + sa**2*(2.0_real_8*a2*dsadra*sa + 3.0_real_8*a3*dsadra*sa**2 +  &
     & 4.0_real_8*a4*dsadra*sa**3 + 5.0_real_8*a5*dsadra*sa**4 + 6.0_real_8*a6*dsadra* &
     & sa**5 + 7.0_real_8*a7*dsadra*sa**6)/(b1*sa + b2*sa**2 + b3*sa**3 + b4*&
     & sa**4 + b5*sa**5 + b6*sa**6 + b7*sa**7 + b8*sa**8 + b9*sa**9 +   &
     & 1.0_real_8)
      dlambadaadra = dzetaadra
      detaadra = dzetaadra
      dFadra = -0.0740740740740740741_real_8*dsadra*sa/(C*(1.0_real_8 + sa**2/s0**&
     & 2)) + 0.0740740740740740741_real_8*dsadra*sa**3/(C*s0**2*(1.0_real_8 + sa**&
     & 2/s0**2)**2) - 0.5_real_8*dzetaadra/C
      dEGadra = -3.6_real_8*A*dlambadaadra*lambadaa**2 -                     &
     & 0.533333333333333333_real_8*B*dlambadaadra*lambadaa - 0.4_real_8*C*Fa*     &
     & dlambadaadra - 0.4_real_8*C*dFadra*lambadaa - 8.4_real_8*dlambadaadra*     &
     & lambadaa**(5.0_real_8/2.0_real_8)*(-SQRT(etaa) + SQRT(zetaa)) -            &
     & 4.962870782535445_real_8*dlambadaadra*lambadaa**(5.0_real_8/2.0_real_8) - 2.4_real_8*&
     & lambadaa**(7.0_real_8/2.0_real_8)*(-0.5_real_8*detaadra/SQRT(etaa) + 0.5_real_8*     &
     & dzetaadra/SQRT(zetaa))
      dFxadra = -A*(-dlambadaadra*y0a**2 - 2.0_real_8*dy0adra*lambadaa*y0a)* &
     & EXP(-lambadaa*y0a**2)/(lambadaa*y0a**2) - detaadra*DEI(-etaa*y0a &
     & **2) - 2.0_real_8*dy0adra*(EXP(-y0a**2*zetaa) - EXP(-etaa*y0a**2))/y0a&
     & **3 + dzetaadra*DEI(-y0a**2*zetaa) + (-dlambadaadra*y0a**2 -     &
     & 2.0_real_8*dy0adra*lambadaa*y0a)*(-0.444444444444444444_real_8*B/lambadaa -&
     & 0.444444444444444444_real_8*C*Fa*(lambadaa*y0a**2 + 1.0_real_8)/lambadaa**2&
     & - 0.444444444444444444_real_8*EGa*(y0a**4/lambadaa + 2.0_real_8*y0a**2/    &
     & lambadaa**2 + 2.0_real_8/lambadaa**3))*EXP(-lambadaa*y0a**2) + (      &
     & 0.444444444444444444_real_8*B*dlambadaadra/lambadaa**2 +              &
     & 0.888888888888888889_real_8*C*Fa*dlambadaadra*(lambadaa*y0a**2 + 1.0_real_8&
     & )/lambadaa**3 - 0.444444444444444444_real_8*C*Fa*(dlambadaadra*y0a**2 &
     & + 2.0_real_8*dy0adra*lambadaa*y0a)/lambadaa**2 -                      &
     & 0.444444444444444444_real_8*C*dFadra*(lambadaa*y0a**2 + 1.0_real_8)/       &
     & lambadaa**2 - 0.444444444444444444_real_8*EGa*(-dlambadaadra*y0a**4/  &
     & lambadaa**2 - 4.0_real_8*dlambadaadra*y0a**2/lambadaa**3 - 6.0_real_8*     &
     & dlambadaadra/lambadaa**4 + 4.0_real_8*dy0adra*y0a**3/lambadaa + 4.0_real_8*&
     & dy0adra*y0a/lambadaa**2) - 0.444444444444444444_real_8*dEGadra*(y0a**4&
     & /lambadaa + 2.0_real_8*y0a**2/lambadaa**2 + 2.0_real_8/lambadaa**3))*EXP(  &
     & -lambadaa*y0a**2) + (-detaadra*y0a**2 - 2.0_real_8*dy0adra*etaa*y0a)* &
     & EXP(-etaa*y0a**2)/y0a**2 + (-(-detaadra*y0a**2 - 2.0_real_8*dy0adra*  &
     & etaa*y0a)*EXP(-etaa*y0a**2) + (-2.0_real_8*dy0adra*y0a*zetaa -        &
     & dzetaadra*y0a**2)*EXP(-y0a**2*zetaa))/y0a**2 - (-2.0_real_8*dy0adra*  &
     & y0a*zetaa - dzetaadra*y0a**2)*EXP(-y0a**2*zetaa)/y0a**2
      dEnerdra = -1.25992104989487316_real_8*Ax*(1.33333333333333333_real_8*Fxa*  &
     & RHOA**(1.0_real_8/3.0_real_8) + RHOA**(4.0_real_8/3.0_real_8)*dFxadra)
      dkfbdrb = 1.29925902990691777_real_8/RHOB**(2.0_real_8/3.0_real_8)
      dy0bdrb = R0*dkfbdrb
      dsbdrb = -0.5_real_8*SQRT(GRHOBB)*dkfbdrb/(RHOB*kfb**2) - 0.5_real_8*SQRT(  &
     & GRHOBB)/(RHOB**2*kfb)
      dzetabdrb = 2.0_real_8*dsbdrb*sb*(a2*sb**2 + a3*sb**3 + a4*sb**4 + a5* &
     & sb**5 + a6*sb**6 + a7*sb**7)/(b1*sb + b2*sb**2 + b3*sb**3 + b4*sb&
     & **4 + b5*sb**5 + b6*sb**6 + b7*sb**7 + b8*sb**8 + b9*sb**9 +     &
     & 1.0_real_8) + sb**2*(a2*sb**2 + a3*sb**3 + a4*sb**4 + a5*sb**5 + a6*sb&
     & **6 + a7*sb**7)*(-b1*dsbdrb - 2.0_real_8*b2*dsbdrb*sb - 3.0_real_8*b3*     &
     & dsbdrb*sb**2 - 4.0_real_8*b4*dsbdrb*sb**3 - 5.0_real_8*b5*dsbdrb*sb**4 -   &
     & 6.0_real_8*b6*dsbdrb*sb**5 - 7.0_real_8*b7*dsbdrb*sb**6 - 8.0_real_8*b8*dsbdrb* &
     & sb**7 - 9.0_real_8*b9*dsbdrb*sb**8)/(b1*sb + b2*sb**2 + b3*sb**3 + b4*&
     & sb**4 + b5*sb**5 + b6*sb**6 + b7*sb**7 + b8*sb**8 + b9*sb**9 +   &
     & 1.0_real_8)**2 + sb**2*(2.0_real_8*a2*dsbdrb*sb + 3.0_real_8*a3*dsbdrb*sb**2 +  &
     & 4.0_real_8*a4*dsbdrb*sb**3 + 5.0_real_8*a5*dsbdrb*sb**4 + 6.0_real_8*a6*dsbdrb* &
     & sb**5 + 7.0_real_8*a7*dsbdrb*sb**6)/(b1*sb + b2*sb**2 + b3*sb**3 + b4*&
     & sb**4 + b5*sb**5 + b6*sb**6 + b7*sb**7 + b8*sb**8 + b9*sb**9 +   &
     & 1.0_real_8)
      dlambadabdrb = dzetabdrb
      detabdrb = dzetabdrb
      dFbdrb = -0.0740740740740740741_real_8*dsbdrb*sb/(C*(1.0_real_8 + sb**2/s0**&
     & 2)) + 0.0740740740740740741_real_8*dsbdrb*sb**3/(C*s0**2*(1.0_real_8 + sb**&
     & 2/s0**2)**2) - 0.5_real_8*dzetabdrb/C
      dEGbdrb = -3.6_real_8*A*dlambadabdrb*lambadab**2 -                     &
     & 0.533333333333333333_real_8*B*dlambadabdrb*lambadab - 0.4_real_8*C*Fb*     &
     & dlambadabdrb - 0.4_real_8*C*dFbdrb*lambadab - 8.4_real_8*dlambadabdrb*     &
     & lambadab**(5.0_real_8/2.0_real_8)*(-SQRT(etab) + SQRT(zetab)) -            &
     & 4.962870782535445_real_8*dlambadabdrb*lambadab**(5.0_real_8/2.0_real_8) - 2.4_real_8*&
     & lambadab**(7.0_real_8/2.0_real_8)*(-0.5_real_8*detabdrb/SQRT(etab) + 0.5_real_8*     &
     & dzetabdrb/SQRT(zetab))
      dFxbdrb = -A*(-dlambadabdrb*y0b**2 - 2.0_real_8*dy0bdrb*lambadab*y0b)* &
     & EXP(-lambadab*y0b**2)/(lambadab*y0b**2) - detabdrb*DEI(-etab*y0b &
     & **2) - 2.0_real_8*dy0bdrb*(EXP(-y0b**2*zetab) - EXP(-etab*y0b**2))/y0b&
     & **3 + dzetabdrb*DEI(-y0b**2*zetab) + (-dlambadabdrb*y0b**2 -     &
     & 2.0_real_8*dy0bdrb*lambadab*y0b)*(-0.444444444444444444_real_8*B/lambadab -&
     & 0.444444444444444444_real_8*C*Fb*(lambadab*y0b**2 + 1.0_real_8)/lambadab**2&
     & - 0.444444444444444444_real_8*EGb*(y0b**4/lambadab + 2.0_real_8*y0b**2/    &
     & lambadab**2 + 2.0_real_8/lambadab**3))*EXP(-lambadab*y0b**2) + (      &
     & 0.444444444444444444_real_8*B*dlambadabdrb/lambadab**2 +              &
     & 0.888888888888888889_real_8*C*Fb*dlambadabdrb*(lambadab*y0b**2 + 1.0_real_8&
     & )/lambadab**3 - 0.444444444444444444_real_8*C*Fb*(dlambadabdrb*y0b**2 &
     & + 2.0_real_8*dy0bdrb*lambadab*y0b)/lambadab**2 -                      &
     & 0.444444444444444444_real_8*C*dFbdrb*(lambadab*y0b**2 + 1.0_real_8)/       &
     & lambadab**2 - 0.444444444444444444_real_8*EGb*(-dlambadabdrb*y0b**4/  &
     & lambadab**2 - 4.0_real_8*dlambadabdrb*y0b**2/lambadab**3 - 6.0_real_8*     &
     & dlambadabdrb/lambadab**4 + 4.0_real_8*dy0bdrb*y0b**3/lambadab + 4.0_real_8*&
     & dy0bdrb*y0b/lambadab**2) - 0.444444444444444444_real_8*dEGbdrb*(y0b**4&
     & /lambadab + 2.0_real_8*y0b**2/lambadab**2 + 2.0_real_8/lambadab**3))*EXP(  &
     & -lambadab*y0b**2) + (-detabdrb*y0b**2 - 2.0_real_8*dy0bdrb*etab*y0b)* &
     & EXP(-etab*y0b**2)/y0b**2 + (-(-detabdrb*y0b**2 - 2.0_real_8*dy0bdrb*  &
     & etab*y0b)*EXP(-etab*y0b**2) + (-2.0_real_8*dy0bdrb*y0b*zetab -        &
     & dzetabdrb*y0b**2)*EXP(-y0b**2*zetab))/y0b**2 - (-2.0_real_8*dy0bdrb*  &
     & y0b*zetab - dzetabdrb*y0b**2)*EXP(-y0b**2*zetab)/y0b**2
      dEnerdrb = -1.25992104989487316_real_8*Ax*(1.33333333333333333_real_8*Fxb*  &
     & RHOB**(1.0_real_8/3.0_real_8) + RHOB**(4.0_real_8/3.0_real_8)*dFxbdrb)
      dsadnaa = 0.25_real_8/(SQRT(GRHOAA)*RHOA*kfa)
      dzetaadnaa = 2.0_real_8*dsadnaa*sa*(a2*sa**2 + a3*sa**3 + a4*sa**4 + a5&
     & *sa**5 + a6*sa**6 + a7*sa**7)/(b1*sa + b2*sa**2 + b3*sa**3 + b4* &
     & sa**4 + b5*sa**5 + b6*sa**6 + b7*sa**7 + b8*sa**8 + b9*sa**9 +   &
     & 1.0_real_8) + sa**2*(a2*sa**2 + a3*sa**3 + a4*sa**4 + a5*sa**5 + a6*sa&
     & **6 + a7*sa**7)*(-b1*dsadnaa - 2.0_real_8*b2*dsadnaa*sa - 3.0_real_8*b3*   &
     & dsadnaa*sa**2 - 4.0_real_8*b4*dsadnaa*sa**3 - 5.0_real_8*b5*dsadnaa*sa**4 -&
     & 6.0_real_8*b6*dsadnaa*sa**5 - 7.0_real_8*b7*dsadnaa*sa**6 - 8.0_real_8*b8*      &
     & dsadnaa*sa**7 - 9.0_real_8*b9*dsadnaa*sa**8)/(b1*sa + b2*sa**2 + b3*sa&
     & **3 + b4*sa**4 + b5*sa**5 + b6*sa**6 + b7*sa**7 + b8*sa**8 + b9* &
     & sa**9 + 1.0_real_8)**2 + sa**2*(2.0_real_8*a2*dsadnaa*sa + 3.0_real_8*a3*dsadnaa&
     & *sa**2 + 4.0_real_8*a4*dsadnaa*sa**3 + 5.0_real_8*a5*dsadnaa*sa**4 + 6.0_real_8*&
     & a6*dsadnaa*sa**5 + 7.0_real_8*a7*dsadnaa*sa**6)/(b1*sa + b2*sa**2 + b3&
     & *sa**3 + b4*sa**4 + b5*sa**5 + b6*sa**6 + b7*sa**7 + b8*sa**8 +  &
     & b9*sa**9 + 1.0_real_8)
      dlambadaadnaa = dzetaadnaa
      detaadnaa = dzetaadnaa
      dFadnaa = -0.0740740740740740741_real_8*dsadnaa*sa/(C*(1.0_real_8 + sa**2/s0&
     & **2)) + 0.0740740740740740741_real_8*dsadnaa*sa**3/(C*s0**2*(1.0_real_8 +  &
     & sa**2/s0**2)**2) - 0.5_real_8*dzetaadnaa/C
      dEGadnaa = -3.6_real_8*A*dlambadaadnaa*lambadaa**2 -                   &
     & 0.533333333333333333_real_8*B*dlambadaadnaa*lambadaa - 0.4_real_8*C*Fa*    &
     & dlambadaadnaa - 0.4_real_8*C*dFadnaa*lambadaa - 8.4_real_8*dlambadaadnaa*  &
     & lambadaa**(5.0_real_8/2.0_real_8)*(-SQRT(etaa) + SQRT(zetaa)) -            &
     & 4.962870782535445_real_8*dlambadaadnaa*lambadaa**(5.0_real_8/2.0_real_8) - 2.4_real_8&
     & *lambadaa**(7.0_real_8/2.0_real_8)*(-0.5_real_8*detaadnaa/SQRT(etaa) + 0.5_real_8*   &
     & dzetaadnaa/SQRT(zetaa))
      dFxadnaa = A*dlambadaadnaa*EXP(-lambadaa*y0a**2)/lambadaa -       &
     & detaadnaa*DEI(-etaa*y0a**2) - detaadnaa*EXP(-etaa*y0a**2) -      &
     & dlambadaadnaa*y0a**2*(-0.444444444444444444_real_8*B/lambadaa -       &
     & 0.444444444444444444_real_8*C*Fa*(lambadaa*y0a**2 + 1.0_real_8)/lambadaa**2&
     & - 0.444444444444444444_real_8*EGa*(y0a**4/lambadaa + 2.0_real_8*y0a**2/    &
     & lambadaa**2 + 2.0_real_8/lambadaa**3))*EXP(-lambadaa*y0a**2) +        &
     & dzetaadnaa*DEI(-y0a**2*zetaa) + dzetaadnaa*EXP(-y0a**2*zetaa) + (&
     & 0.444444444444444444_real_8*B*dlambadaadnaa/lambadaa**2 -             &
     & 0.444444444444444444_real_8*C*Fa*dlambadaadnaa*y0a**2/lambadaa**2 +   &
     & 0.888888888888888889_real_8*C*Fa*dlambadaadnaa*(lambadaa*y0a**2 +     &
     & 1.0_real_8)/lambadaa**3 - 0.444444444444444444_real_8*C*dFadnaa*(lambadaa* &
     & y0a**2 + 1.0_real_8)/lambadaa**2 - 0.444444444444444444_real_8*EGa*(       &
     & -dlambadaadnaa*y0a**4/lambadaa**2 - 4.0_real_8*dlambadaadnaa*y0a**2/  &
     & lambadaa**3 - 6.0_real_8*dlambadaadnaa/lambadaa**4) -                 &
     & 0.444444444444444444_real_8*dEGadnaa*(y0a**4/lambadaa + 2.0_real_8*y0a**2/ &
     & lambadaa**2 + 2.0_real_8/lambadaa**3))*EXP(-lambadaa*y0a**2) + (      &
     & detaadnaa*y0a**2*EXP(-etaa*y0a**2) - dzetaadnaa*y0a**2*EXP(-y0a**&
     & 2*zetaa))/y0a**2
      dEnerdnaa = -1.25992104989487316_real_8*Ax*RHOA**(4.0_real_8/3.0_real_8)*dFxadnaa
      dsbdnbb = 0.25_real_8/(SQRT(GRHOBB)*RHOB*kfb)
      dzetabdnbb = 2.0_real_8*dsbdnbb*sb*(a2*sb**2 + a3*sb**3 + a4*sb**4 + a5&
     & *sb**5 + a6*sb**6 + a7*sb**7)/(b1*sb + b2*sb**2 + b3*sb**3 + b4* &
     & sb**4 + b5*sb**5 + b6*sb**6 + b7*sb**7 + b8*sb**8 + b9*sb**9 +   &
     & 1.0_real_8) + sb**2*(a2*sb**2 + a3*sb**3 + a4*sb**4 + a5*sb**5 + a6*sb&
     & **6 + a7*sb**7)*(-b1*dsbdnbb - 2.0_real_8*b2*dsbdnbb*sb - 3.0_real_8*b3*   &
     & dsbdnbb*sb**2 - 4.0_real_8*b4*dsbdnbb*sb**3 - 5.0_real_8*b5*dsbdnbb*sb**4 -&
     & 6.0_real_8*b6*dsbdnbb*sb**5 - 7.0_real_8*b7*dsbdnbb*sb**6 - 8.0_real_8*b8*      &
     & dsbdnbb*sb**7 - 9.0_real_8*b9*dsbdnbb*sb**8)/(b1*sb + b2*sb**2 + b3*sb&
     & **3 + b4*sb**4 + b5*sb**5 + b6*sb**6 + b7*sb**7 + b8*sb**8 + b9* &
     & sb**9 + 1.0_real_8)**2 + sb**2*(2.0_real_8*a2*dsbdnbb*sb + 3.0_real_8*a3*dsbdnbb&
     & *sb**2 + 4.0_real_8*a4*dsbdnbb*sb**3 + 5.0_real_8*a5*dsbdnbb*sb**4 + 6.0_real_8*&
     & a6*dsbdnbb*sb**5 + 7.0_real_8*a7*dsbdnbb*sb**6)/(b1*sb + b2*sb**2 + b3&
     & *sb**3 + b4*sb**4 + b5*sb**5 + b6*sb**6 + b7*sb**7 + b8*sb**8 +  &
     & b9*sb**9 + 1.0_real_8)
      dlambadabdnbb = dzetabdnbb
      detabdnbb = dzetabdnbb
      dFbdnbb = -0.0740740740740740741_real_8*dsbdnbb*sb/(C*(1.0_real_8 + sb**2/s0&
     & **2)) + 0.0740740740740740741_real_8*dsbdnbb*sb**3/(C*s0**2*(1.0_real_8 +  &
     & sb**2/s0**2)**2) - 0.5_real_8*dzetabdnbb/C
      dEGbdnbb = -3.6_real_8*A*dlambadabdnbb*lambadab**2 -                   &
     & 0.533333333333333333_real_8*B*dlambadabdnbb*lambadab - 0.4_real_8*C*Fb*    &
     & dlambadabdnbb - 0.4_real_8*C*dFbdnbb*lambadab - 8.4_real_8*dlambadabdnbb*  &
     & lambadab**(5.0_real_8/2.0_real_8)*(-SQRT(etab) + SQRT(zetab)) -            &
     & 4.962870782535445_real_8*dlambadabdnbb*lambadab**(5.0_real_8/2.0_real_8) - 2.4_real_8&
     & *lambadab**(7.0_real_8/2.0_real_8)*(-0.5_real_8*detabdnbb/SQRT(etab) + 0.5_real_8*   &
     & dzetabdnbb/SQRT(zetab))
      dFxbdnbb = A*dlambadabdnbb*EXP(-lambadab*y0b**2)/lambadab -       &
     & detabdnbb*DEI(-etab*y0b**2) - detabdnbb*EXP(-etab*y0b**2) -      &
     & dlambadabdnbb*y0b**2*(-0.444444444444444444_real_8*B/lambadab -       &
     & 0.444444444444444444_real_8*C*Fb*(lambadab*y0b**2 + 1.0_real_8)/lambadab**2&
     & - 0.444444444444444444_real_8*EGb*(y0b**4/lambadab + 2.0_real_8*y0b**2/    &
     & lambadab**2 + 2.0_real_8/lambadab**3))*EXP(-lambadab*y0b**2) +        &
     & dzetabdnbb*DEI(-y0b**2*zetab) + dzetabdnbb*EXP(-y0b**2*zetab) + (&
     & 0.444444444444444444_real_8*B*dlambadabdnbb/lambadab**2 -             &
     & 0.444444444444444444_real_8*C*Fb*dlambadabdnbb*y0b**2/lambadab**2 +   &
     & 0.888888888888888889_real_8*C*Fb*dlambadabdnbb*(lambadab*y0b**2 +     &
     & 1.0_real_8)/lambadab**3 - 0.444444444444444444_real_8*C*dFbdnbb*(lambadab* &
     & y0b**2 + 1.0_real_8)/lambadab**2 - 0.444444444444444444_real_8*EGb*(       &
     & -dlambadabdnbb*y0b**4/lambadab**2 - 4.0_real_8*dlambadabdnbb*y0b**2/  &
     & lambadab**3 - 6.0_real_8*dlambadabdnbb/lambadab**4) -                 &
     & 0.444444444444444444_real_8*dEGbdnbb*(y0b**4/lambadab + 2.0_real_8*y0b**2/ &
     & lambadab**2 + 2.0_real_8/lambadab**3))*EXP(-lambadab*y0b**2) + (      &
     & detabdnbb*y0b**2*EXP(-etab*y0b**2) - dzetabdnbb*y0b**2*EXP(-y0b**&
     & 2*zetab))/y0b**2
      dEnerdnbb = -1.25992104989487316_real_8*Ax*RHOB**(4.0_real_8/3.0_real_8)*dFxbdnbb
      ELSE ! RHOB.GT.TOL
      Ax = 0.738558766382022336_real_8
      kfa = 3.89777708972075354_real_8*RHOA**(1.0_real_8/3.0_real_8)
      y0a = R0*kfa
      sa = 0.5_real_8*SQRT(GRHOAA)/(RHOA*kfa)
      zetaa = sa**2*(a2*sa**2 + a3*sa**3 + a4*sa**4 + a5*sa**5 + a6*sa**&
     & 6 + a7*sa**7)/(b1*sa + b2*sa**2 + b3*sa**3 + b4*sa**4 + b5*sa**5 &
     & + b6*sa**6 + b7*sa**7 + b8*sa**8 + b9*sa**9 + 1.0_real_8)
      lambadaa = D + zetaa
      lambadab = D
      etaa = A + zetaa
      etab = A
      Fa = 1.0_real_8 - 0.037037037037037037_real_8*sa**2/(C*(1.0_real_8 + sa**2/s0**2 &
     & )) - 0.5_real_8*zetaa/C
      Fb = 1.00000000000000000_real_8
      EGa = -1.2_real_8*A*lambadaa**3 - 0.266666666666666667_real_8*B*lambadaa**2 &
     & - 0.4_real_8*C*Fa*lambadaa - 2.4_real_8*lambadaa**(7.0_real_8/2.0_real_8)*(-SQRT(etaa&
     & ) + SQRT(zetaa)) - 1.41796308072441279_real_8*lambadaa**(7.0_real_8/2.0_real_8)
      Fxa = A*DEI(-lambadaa*y0a**2) - etaa*DEI(-etaa*y0a**2) + zetaa*DEI&
     & (-y0a**2*zetaa) + (-0.444444444444444444_real_8*B/lambadaa -          &
     & 0.444444444444444444_real_8*C*Fa*(lambadaa*y0a**2 + 1.0_real_8)/lambadaa**2&
     & - 0.444444444444444444_real_8*EGa*(y0a**4/lambadaa + 2.0_real_8*y0a**2/    &
     & lambadaa**2 + 2.0_real_8/lambadaa**3))*EXP(-lambadaa*y0a**2) + (EXP(  &
     & -y0a**2*zetaa) - EXP(-etaa*y0a**2))/y0a**2
      Ener = -1.25992104989487316_real_8*Ax*Fxa*RHOA**(4.0_real_8/3.0_real_8)
      dkfadra = 1.29925902990691777_real_8/RHOA**(2.0_real_8/3.0_real_8)
      dy0adra = R0*dkfadra
      dsadra = -0.5_real_8*SQRT(GRHOAA)*dkfadra/(RHOA*kfa**2) - 0.5_real_8*SQRT(  &
     & GRHOAA)/(RHOA**2*kfa)
      dzetaadra = 2.0_real_8*dsadra*sa*(a2*sa**2 + a3*sa**3 + a4*sa**4 + a5* &
     & sa**5 + a6*sa**6 + a7*sa**7)/(b1*sa + b2*sa**2 + b3*sa**3 + b4*sa&
     & **4 + b5*sa**5 + b6*sa**6 + b7*sa**7 + b8*sa**8 + b9*sa**9 +     &
     & 1.0_real_8) + sa**2*(a2*sa**2 + a3*sa**3 + a4*sa**4 + a5*sa**5 + a6*sa&
     & **6 + a7*sa**7)*(-b1*dsadra - 2.0_real_8*b2*dsadra*sa - 3.0_real_8*b3*     &
     & dsadra*sa**2 - 4.0_real_8*b4*dsadra*sa**3 - 5.0_real_8*b5*dsadra*sa**4 -   &
     & 6.0_real_8*b6*dsadra*sa**5 - 7.0_real_8*b7*dsadra*sa**6 - 8.0_real_8*b8*dsadra* &
     & sa**7 - 9.0_real_8*b9*dsadra*sa**8)/(b1*sa + b2*sa**2 + b3*sa**3 + b4*&
     & sa**4 + b5*sa**5 + b6*sa**6 + b7*sa**7 + b8*sa**8 + b9*sa**9 +   &
     & 1.0_real_8)**2 + sa**2*(2.0_real_8*a2*dsadra*sa + 3.0_real_8*a3*dsadra*sa**2 +  &
     & 4.0_real_8*a4*dsadra*sa**3 + 5.0_real_8*a5*dsadra*sa**4 + 6.0_real_8*a6*dsadra* &
     & sa**5 + 7.0_real_8*a7*dsadra*sa**6)/(b1*sa + b2*sa**2 + b3*sa**3 + b4*&
     & sa**4 + b5*sa**5 + b6*sa**6 + b7*sa**7 + b8*sa**8 + b9*sa**9 +   &
     & 1.0_real_8)
      dlambadaadra = dzetaadra
      detaadra = dzetaadra
      dFadra = -0.0740740740740740741_real_8*dsadra*sa/(C*(1.0_real_8 + sa**2/s0**&
     & 2)) + 0.0740740740740740741_real_8*dsadra*sa**3/(C*s0**2*(1.0_real_8 + sa**&
     & 2/s0**2)**2) - 0.5_real_8*dzetaadra/C
      dEGadra = -3.6_real_8*A*dlambadaadra*lambadaa**2 -                     &
     & 0.533333333333333333_real_8*B*dlambadaadra*lambadaa - 0.4_real_8*C*Fa*     &
     & dlambadaadra - 0.4_real_8*C*dFadra*lambadaa - 8.4_real_8*dlambadaadra*     &
     & lambadaa**(5.0_real_8/2.0_real_8)*(-SQRT(etaa) + SQRT(zetaa)) -            &
     & 4.962870782535445_real_8*dlambadaadra*lambadaa**(5.0_real_8/2.0_real_8) - 2.4_real_8*&
     & lambadaa**(7.0_real_8/2.0_real_8)*(-0.5_real_8*detaadra/SQRT(etaa) + 0.5_real_8*     &
     & dzetaadra/SQRT(zetaa))
      dFxadra = -A*(-dlambadaadra*y0a**2 - 2.0_real_8*dy0adra*lambadaa*y0a)* &
     & EXP(-lambadaa*y0a**2)/(lambadaa*y0a**2) - detaadra*DEI(-etaa*y0a &
     & **2) - 2.0_real_8*dy0adra*(EXP(-y0a**2*zetaa) - EXP(-etaa*y0a**2))/y0a&
     & **3 + dzetaadra*DEI(-y0a**2*zetaa) + (-dlambadaadra*y0a**2 -     &
     & 2.0_real_8*dy0adra*lambadaa*y0a)*(-0.444444444444444444_real_8*B/lambadaa -&
     & 0.444444444444444444_real_8*C*Fa*(lambadaa*y0a**2 + 1.0_real_8)/lambadaa**2&
     & - 0.444444444444444444_real_8*EGa*(y0a**4/lambadaa + 2.0_real_8*y0a**2/    &
     & lambadaa**2 + 2.0_real_8/lambadaa**3))*EXP(-lambadaa*y0a**2) + (      &
     & 0.444444444444444444_real_8*B*dlambadaadra/lambadaa**2 +              &
     & 0.888888888888888889_real_8*C*Fa*dlambadaadra*(lambadaa*y0a**2 + 1.0_real_8&
     & )/lambadaa**3 - 0.444444444444444444_real_8*C*Fa*(dlambadaadra*y0a**2 &
     & + 2.0_real_8*dy0adra*lambadaa*y0a)/lambadaa**2 -                      &
     & 0.444444444444444444_real_8*C*dFadra*(lambadaa*y0a**2 + 1.0_real_8)/       &
     & lambadaa**2 - 0.444444444444444444_real_8*EGa*(-dlambadaadra*y0a**4/  &
     & lambadaa**2 - 4.0_real_8*dlambadaadra*y0a**2/lambadaa**3 - 6.0_real_8*     &
     & dlambadaadra/lambadaa**4 + 4.0_real_8*dy0adra*y0a**3/lambadaa + 4.0_real_8*&
     & dy0adra*y0a/lambadaa**2) - 0.444444444444444444_real_8*dEGadra*(y0a**4&
     & /lambadaa + 2.0_real_8*y0a**2/lambadaa**2 + 2.0_real_8/lambadaa**3))*EXP(  &
     & -lambadaa*y0a**2) + (-detaadra*y0a**2 - 2.0_real_8*dy0adra*etaa*y0a)* &
     & EXP(-etaa*y0a**2)/y0a**2 + (-(-detaadra*y0a**2 - 2.0_real_8*dy0adra*  &
     & etaa*y0a)*EXP(-etaa*y0a**2) + (-2.0_real_8*dy0adra*y0a*zetaa -        &
     & dzetaadra*y0a**2)*EXP(-y0a**2*zetaa))/y0a**2 - (-2.0_real_8*dy0adra*  &
     & y0a*zetaa - dzetaadra*y0a**2)*EXP(-y0a**2*zetaa)/y0a**2
      dEnerdra = -1.25992104989487316_real_8*Ax*(1.33333333333333333_real_8*Fxa*  &
     & RHOA**(1.0_real_8/3.0_real_8) + RHOA**(4.0_real_8/3.0_real_8)*dFxadra)
      dsadnaa = 0.25_real_8/(SQRT(GRHOAA)*RHOA*kfa)
      dzetaadnaa = 2.0_real_8*dsadnaa*sa*(a2*sa**2 + a3*sa**3 + a4*sa**4 + a5&
     & *sa**5 + a6*sa**6 + a7*sa**7)/(b1*sa + b2*sa**2 + b3*sa**3 + b4* &
     & sa**4 + b5*sa**5 + b6*sa**6 + b7*sa**7 + b8*sa**8 + b9*sa**9 +   &
     & 1.0_real_8) + sa**2*(a2*sa**2 + a3*sa**3 + a4*sa**4 + a5*sa**5 + a6*sa&
     & **6 + a7*sa**7)*(-b1*dsadnaa - 2.0_real_8*b2*dsadnaa*sa - 3.0_real_8*b3*   &
     & dsadnaa*sa**2 - 4.0_real_8*b4*dsadnaa*sa**3 - 5.0_real_8*b5*dsadnaa*sa**4 -&
     & 6.0_real_8*b6*dsadnaa*sa**5 - 7.0_real_8*b7*dsadnaa*sa**6 - 8.0_real_8*b8*      &
     & dsadnaa*sa**7 - 9.0_real_8*b9*dsadnaa*sa**8)/(b1*sa + b2*sa**2 + b3*sa&
     & **3 + b4*sa**4 + b5*sa**5 + b6*sa**6 + b7*sa**7 + b8*sa**8 + b9* &
     & sa**9 + 1.0_real_8)**2 + sa**2*(2.0_real_8*a2*dsadnaa*sa + 3.0_real_8*a3*dsadnaa&
     & *sa**2 + 4.0_real_8*a4*dsadnaa*sa**3 + 5.0_real_8*a5*dsadnaa*sa**4 + 6.0_real_8*&
     & a6*dsadnaa*sa**5 + 7.0_real_8*a7*dsadnaa*sa**6)/(b1*sa + b2*sa**2 + b3&
     & *sa**3 + b4*sa**4 + b5*sa**5 + b6*sa**6 + b7*sa**7 + b8*sa**8 +  &
     & b9*sa**9 + 1.0_real_8)
      dlambadaadnaa = dzetaadnaa
      detaadnaa = dzetaadnaa
      dFadnaa = -0.0740740740740740741_real_8*dsadnaa*sa/(C*(1.0_real_8 + sa**2/s0&
     & **2)) + 0.0740740740740740741_real_8*dsadnaa*sa**3/(C*s0**2*(1.0_real_8 +  &
     & sa**2/s0**2)**2) - 0.5_real_8*dzetaadnaa/C
      dEGadnaa = -3.6_real_8*A*dlambadaadnaa*lambadaa**2 -                   &
     & 0.533333333333333333_real_8*B*dlambadaadnaa*lambadaa - 0.4_real_8*C*Fa*    &
     & dlambadaadnaa - 0.4_real_8*C*dFadnaa*lambadaa - 8.4_real_8*dlambadaadnaa*  &
     & lambadaa**(5.0_real_8/2.0_real_8)*(-SQRT(etaa) + SQRT(zetaa)) -            &
     & 4.962870782535445_real_8*dlambadaadnaa*lambadaa**(5.0_real_8/2.0_real_8) - 2.4_real_8&
     & *lambadaa**(7.0_real_8/2.0_real_8)*(-0.5_real_8*detaadnaa/SQRT(etaa) + 0.5_real_8*   &
     & dzetaadnaa/SQRT(zetaa))
      dFxadnaa = A*dlambadaadnaa*EXP(-lambadaa*y0a**2)/lambadaa -       &
     & detaadnaa*DEI(-etaa*y0a**2) - detaadnaa*EXP(-etaa*y0a**2) -      &
     & dlambadaadnaa*y0a**2*(-0.444444444444444444_real_8*B/lambadaa -       &
     & 0.444444444444444444_real_8*C*Fa*(lambadaa*y0a**2 + 1.0_real_8)/lambadaa**2&
     & - 0.444444444444444444_real_8*EGa*(y0a**4/lambadaa + 2.0_real_8*y0a**2/    &
     & lambadaa**2 + 2.0_real_8/lambadaa**3))*EXP(-lambadaa*y0a**2) +        &
     & dzetaadnaa*DEI(-y0a**2*zetaa) + dzetaadnaa*EXP(-y0a**2*zetaa) + (&
     & 0.444444444444444444_real_8*B*dlambadaadnaa/lambadaa**2 -             &
     & 0.444444444444444444_real_8*C*Fa*dlambadaadnaa*y0a**2/lambadaa**2 +   &
     & 0.888888888888888889_real_8*C*Fa*dlambadaadnaa*(lambadaa*y0a**2 +     &
     & 1.0_real_8)/lambadaa**3 - 0.444444444444444444_real_8*C*dFadnaa*(lambadaa* &
     & y0a**2 + 1.0_real_8)/lambadaa**2 - 0.444444444444444444_real_8*EGa*(       &
     & -dlambadaadnaa*y0a**4/lambadaa**2 - 4.0_real_8*dlambadaadnaa*y0a**2/  &
     & lambadaa**3 - 6.0_real_8*dlambadaadnaa/lambadaa**4) -                 &
     & 0.444444444444444444_real_8*dEGadnaa*(y0a**4/lambadaa + 2.0_real_8*y0a**2/ &
     & lambadaa**2 + 2.0_real_8/lambadaa**3))*EXP(-lambadaa*y0a**2) + (      &
     & detaadnaa*y0a**2*EXP(-etaa*y0a**2) - dzetaadnaa*y0a**2*EXP(-y0a**&
     & 2*zetaa))/y0a**2
      dEnerdnaa = -1.25992104989487316_real_8*Ax*RHOA**(4.0_real_8/3.0_real_8)*dFxadnaa
      ENDIF ! RHOB.GT.TOL
      ELSE ! RHOA.GT.TOL
      IF(RHOB.GT.TOL) THEN
      Ax = 0.738558766382022336_real_8
      kfb = 3.89777708972075354_real_8*RHOB**(1.0_real_8/3.0_real_8)
      y0b = R0*kfb
      sb = 0.5_real_8*SQRT(GRHOBB)/(RHOB*kfb)
      zetab = sb**2*(a2*sb**2 + a3*sb**3 + a4*sb**4 + a5*sb**5 + a6*sb**&
     & 6 + a7*sb**7)/(b1*sb + b2*sb**2 + b3*sb**3 + b4*sb**4 + b5*sb**5 &
     & + b6*sb**6 + b7*sb**7 + b8*sb**8 + b9*sb**9 + 1.0_real_8)
      lambadaa = D
      lambadab = D + zetab
      etaa = A
      etab = A + zetab
      Fa = 1.00000000000000000_real_8
      Fb = 1.0_real_8 - 0.037037037037037037_real_8*sb**2/(C*(1.0_real_8 + sb**2/s0**2 &
     & )) - 0.5_real_8*zetab/C
      EGb = -1.2_real_8*A*lambadab**3 - 0.266666666666666667_real_8*B*lambadab**2 &
     & - 0.4_real_8*C*Fb*lambadab - 2.4_real_8*lambadab**(7.0_real_8/2.0_real_8)*(-SQRT(etab&
     & ) + SQRT(zetab)) - 1.41796308072441279_real_8*lambadab**(7.0_real_8/2.0_real_8)
      Fxb = A*DEI(-lambadab*y0b**2) - etab*DEI(-etab*y0b**2) + zetab*DEI&
     & (-y0b**2*zetab) + (-0.444444444444444444_real_8*B/lambadab -          &
     & 0.444444444444444444_real_8*C*Fb*(lambadab*y0b**2 + 1.0_real_8)/lambadab**2&
     & - 0.444444444444444444_real_8*EGb*(y0b**4/lambadab + 2.0_real_8*y0b**2/    &
     & lambadab**2 + 2.0_real_8/lambadab**3))*EXP(-lambadab*y0b**2) + (EXP(  &
     & -y0b**2*zetab) - EXP(-etab*y0b**2))/y0b**2
      Ener = -1.25992104989487316_real_8*Ax*Fxb*RHOB**(4.0_real_8/3.0_real_8)
      dkfbdrb = 1.29925902990691777_real_8/RHOB**(2.0_real_8/3.0_real_8)
      dy0bdrb = R0*dkfbdrb
      dsbdrb = -0.5_real_8*SQRT(GRHOBB)*dkfbdrb/(RHOB*kfb**2) - 0.5_real_8*SQRT(  &
     & GRHOBB)/(RHOB**2*kfb)
      dzetabdrb = 2.0_real_8*dsbdrb*sb*(a2*sb**2 + a3*sb**3 + a4*sb**4 + a5* &
     & sb**5 + a6*sb**6 + a7*sb**7)/(b1*sb + b2*sb**2 + b3*sb**3 + b4*sb&
     & **4 + b5*sb**5 + b6*sb**6 + b7*sb**7 + b8*sb**8 + b9*sb**9 +     &
     & 1.0_real_8) + sb**2*(a2*sb**2 + a3*sb**3 + a4*sb**4 + a5*sb**5 + a6*sb&
     & **6 + a7*sb**7)*(-b1*dsbdrb - 2.0_real_8*b2*dsbdrb*sb - 3.0_real_8*b3*     &
     & dsbdrb*sb**2 - 4.0_real_8*b4*dsbdrb*sb**3 - 5.0_real_8*b5*dsbdrb*sb**4 -   &
     & 6.0_real_8*b6*dsbdrb*sb**5 - 7.0_real_8*b7*dsbdrb*sb**6 - 8.0_real_8*b8*dsbdrb* &
     & sb**7 - 9.0_real_8*b9*dsbdrb*sb**8)/(b1*sb + b2*sb**2 + b3*sb**3 + b4*&
     & sb**4 + b5*sb**5 + b6*sb**6 + b7*sb**7 + b8*sb**8 + b9*sb**9 +   &
     & 1.0_real_8)**2 + sb**2*(2.0_real_8*a2*dsbdrb*sb + 3.0_real_8*a3*dsbdrb*sb**2 +  &
     & 4.0_real_8*a4*dsbdrb*sb**3 + 5.0_real_8*a5*dsbdrb*sb**4 + 6.0_real_8*a6*dsbdrb* &
     & sb**5 + 7.0_real_8*a7*dsbdrb*sb**6)/(b1*sb + b2*sb**2 + b3*sb**3 + b4*&
     & sb**4 + b5*sb**5 + b6*sb**6 + b7*sb**7 + b8*sb**8 + b9*sb**9 +   &
     & 1.0_real_8)
      dlambadabdrb = dzetabdrb
      detabdrb = dzetabdrb
      dFbdrb = -0.0740740740740740741_real_8*dsbdrb*sb/(C*(1.0_real_8 + sb**2/s0**&
     & 2)) + 0.0740740740740740741_real_8*dsbdrb*sb**3/(C*s0**2*(1.0_real_8 + sb**&
     & 2/s0**2)**2) - 0.5_real_8*dzetabdrb/C
      dEGbdrb = -3.6_real_8*A*dlambadabdrb*lambadab**2 -                     &
     & 0.533333333333333333_real_8*B*dlambadabdrb*lambadab - 0.4_real_8*C*Fb*     &
     & dlambadabdrb - 0.4_real_8*C*dFbdrb*lambadab - 8.4_real_8*dlambadabdrb*     &
     & lambadab**(5.0_real_8/2.0_real_8)*(-SQRT(etab) + SQRT(zetab)) -            &
     & 4.962870782535445_real_8*dlambadabdrb*lambadab**(5.0_real_8/2.0_real_8) - 2.4_real_8*&
     & lambadab**(7.0_real_8/2.0_real_8)*(-0.5_real_8*detabdrb/SQRT(etab) + 0.5_real_8*     &
     & dzetabdrb/SQRT(zetab))
      dFxbdrb = -A*(-dlambadabdrb*y0b**2 - 2.0_real_8*dy0bdrb*lambadab*y0b)* &
     & EXP(-lambadab*y0b**2)/(lambadab*y0b**2) - detabdrb*DEI(-etab*y0b &
     & **2) - 2.0_real_8*dy0bdrb*(EXP(-y0b**2*zetab) - EXP(-etab*y0b**2))/y0b&
     & **3 + dzetabdrb*DEI(-y0b**2*zetab) + (-dlambadabdrb*y0b**2 -     &
     & 2.0_real_8*dy0bdrb*lambadab*y0b)*(-0.444444444444444444_real_8*B/lambadab -&
     & 0.444444444444444444_real_8*C*Fb*(lambadab*y0b**2 + 1.0_real_8)/lambadab**2&
     & - 0.444444444444444444_real_8*EGb*(y0b**4/lambadab + 2.0_real_8*y0b**2/    &
     & lambadab**2 + 2.0_real_8/lambadab**3))*EXP(-lambadab*y0b**2) + (      &
     & 0.444444444444444444_real_8*B*dlambadabdrb/lambadab**2 +              &
     & 0.888888888888888889_real_8*C*Fb*dlambadabdrb*(lambadab*y0b**2 + 1.0_real_8&
     & )/lambadab**3 - 0.444444444444444444_real_8*C*Fb*(dlambadabdrb*y0b**2 &
     & + 2.0_real_8*dy0bdrb*lambadab*y0b)/lambadab**2 -                      &
     & 0.444444444444444444_real_8*C*dFbdrb*(lambadab*y0b**2 + 1.0_real_8)/       &
     & lambadab**2 - 0.444444444444444444_real_8*EGb*(-dlambadabdrb*y0b**4/  &
     & lambadab**2 - 4.0_real_8*dlambadabdrb*y0b**2/lambadab**3 - 6.0_real_8*     &
     & dlambadabdrb/lambadab**4 + 4.0_real_8*dy0bdrb*y0b**3/lambadab + 4.0_real_8*&
     & dy0bdrb*y0b/lambadab**2) - 0.444444444444444444_real_8*dEGbdrb*(y0b**4&
     & /lambadab + 2.0_real_8*y0b**2/lambadab**2 + 2.0_real_8/lambadab**3))*EXP(  &
     & -lambadab*y0b**2) + (-detabdrb*y0b**2 - 2.0_real_8*dy0bdrb*etab*y0b)* &
     & EXP(-etab*y0b**2)/y0b**2 + (-(-detabdrb*y0b**2 - 2.0_real_8*dy0bdrb*  &
     & etab*y0b)*EXP(-etab*y0b**2) + (-2.0_real_8*dy0bdrb*y0b*zetab -        &
     & dzetabdrb*y0b**2)*EXP(-y0b**2*zetab))/y0b**2 - (-2.0_real_8*dy0bdrb*  &
     & y0b*zetab - dzetabdrb*y0b**2)*EXP(-y0b**2*zetab)/y0b**2
      dEnerdrb = -1.25992104989487316_real_8*Ax*(1.33333333333333333_real_8*Fxb*  &
     & RHOB**(1.0_real_8/3.0_real_8) + RHOB**(4.0_real_8/3.0_real_8)*dFxbdrb)
      dsbdnbb = 0.25_real_8/(SQRT(GRHOBB)*RHOB*kfb)
      dzetabdnbb = 2.0_real_8*dsbdnbb*sb*(a2*sb**2 + a3*sb**3 + a4*sb**4 + a5&
     & *sb**5 + a6*sb**6 + a7*sb**7)/(b1*sb + b2*sb**2 + b3*sb**3 + b4* &
     & sb**4 + b5*sb**5 + b6*sb**6 + b7*sb**7 + b8*sb**8 + b9*sb**9 +   &
     & 1.0_real_8) + sb**2*(a2*sb**2 + a3*sb**3 + a4*sb**4 + a5*sb**5 + a6*sb&
     & **6 + a7*sb**7)*(-b1*dsbdnbb - 2.0_real_8*b2*dsbdnbb*sb - 3.0_real_8*b3*   &
     & dsbdnbb*sb**2 - 4.0_real_8*b4*dsbdnbb*sb**3 - 5.0_real_8*b5*dsbdnbb*sb**4 -&
     & 6.0_real_8*b6*dsbdnbb*sb**5 - 7.0_real_8*b7*dsbdnbb*sb**6 - 8.0_real_8*b8*      &
     & dsbdnbb*sb**7 - 9.0_real_8*b9*dsbdnbb*sb**8)/(b1*sb + b2*sb**2 + b3*sb&
     & **3 + b4*sb**4 + b5*sb**5 + b6*sb**6 + b7*sb**7 + b8*sb**8 + b9* &
     & sb**9 + 1.0_real_8)**2 + sb**2*(2.0_real_8*a2*dsbdnbb*sb + 3.0_real_8*a3*dsbdnbb&
     & *sb**2 + 4.0_real_8*a4*dsbdnbb*sb**3 + 5.0_real_8*a5*dsbdnbb*sb**4 + 6.0_real_8*&
     & a6*dsbdnbb*sb**5 + 7.0_real_8*a7*dsbdnbb*sb**6)/(b1*sb + b2*sb**2 + b3&
     & *sb**3 + b4*sb**4 + b5*sb**5 + b6*sb**6 + b7*sb**7 + b8*sb**8 +  &
     & b9*sb**9 + 1.0_real_8)
      dlambadabdnbb = dzetabdnbb
      detabdnbb = dzetabdnbb
      dFbdnbb = -0.0740740740740740741_real_8*dsbdnbb*sb/(C*(1.0_real_8 + sb**2/s0&
     & **2)) + 0.0740740740740740741_real_8*dsbdnbb*sb**3/(C*s0**2*(1.0_real_8 +  &
     & sb**2/s0**2)**2) - 0.5_real_8*dzetabdnbb/C
      dEGbdnbb = -3.6_real_8*A*dlambadabdnbb*lambadab**2 -                   &
     & 0.533333333333333333_real_8*B*dlambadabdnbb*lambadab - 0.4_real_8*C*Fb*    &
     & dlambadabdnbb - 0.4_real_8*C*dFbdnbb*lambadab - 8.4_real_8*dlambadabdnbb*  &
     & lambadab**(5.0_real_8/2.0_real_8)*(-SQRT(etab) + SQRT(zetab)) -            &
     & 4.962870782535445_real_8*dlambadabdnbb*lambadab**(5.0_real_8/2.0_real_8) - 2.4_real_8&
     & *lambadab**(7.0_real_8/2.0_real_8)*(-0.5_real_8*detabdnbb/SQRT(etab) + 0.5_real_8*   &
     & dzetabdnbb/SQRT(zetab))
      dFxbdnbb = A*dlambadabdnbb*EXP(-lambadab*y0b**2)/lambadab -       &
     & detabdnbb*DEI(-etab*y0b**2) - detabdnbb*EXP(-etab*y0b**2) -      &
     & dlambadabdnbb*y0b**2*(-0.444444444444444444_real_8*B/lambadab -       &
     & 0.444444444444444444_real_8*C*Fb*(lambadab*y0b**2 + 1.0_real_8)/lambadab**2&
     & - 0.444444444444444444_real_8*EGb*(y0b**4/lambadab + 2.0_real_8*y0b**2/    &
     & lambadab**2 + 2.0_real_8/lambadab**3))*EXP(-lambadab*y0b**2) +        &
     & dzetabdnbb*DEI(-y0b**2*zetab) + dzetabdnbb*EXP(-y0b**2*zetab) + (&
     & 0.444444444444444444_real_8*B*dlambadabdnbb/lambadab**2 -             &
     & 0.444444444444444444_real_8*C*Fb*dlambadabdnbb*y0b**2/lambadab**2 +   &
     & 0.888888888888888889_real_8*C*Fb*dlambadabdnbb*(lambadab*y0b**2 +     &
     & 1.0_real_8)/lambadab**3 - 0.444444444444444444_real_8*C*dFbdnbb*(lambadab* &
     & y0b**2 + 1.0_real_8)/lambadab**2 - 0.444444444444444444_real_8*EGb*(       &
     & -dlambadabdnbb*y0b**4/lambadab**2 - 4.0_real_8*dlambadabdnbb*y0b**2/  &
     & lambadab**3 - 6.0_real_8*dlambadabdnbb/lambadab**4) -                 &
     & 0.444444444444444444_real_8*dEGbdnbb*(y0b**4/lambadab + 2.0_real_8*y0b**2/ &
     & lambadab**2 + 2.0_real_8/lambadab**3))*EXP(-lambadab*y0b**2) + (      &
     & detabdnbb*y0b**2*EXP(-etab*y0b**2) - dzetabdnbb*y0b**2*EXP(-y0b**&
     & 2*zetab))/y0b**2
      dEnerdnbb = -1.25992104989487316_real_8*Ax*RHOB**(4.0_real_8/3.0_real_8)*dFxbdnbb
      ENDIF ! RHOB.GT.TOL
      ENDIF ! RHOA.GT.TOL
