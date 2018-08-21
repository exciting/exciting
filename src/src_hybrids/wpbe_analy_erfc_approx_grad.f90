      SUBROUTINE wpbe_analy_erfc_approx_grad(rho,s,omega,Fx_wpbe, &
                      d1rfx,d1sfx)
!--------------------------------------------------------------------
!
!     wPBE Enhancement Factor (erfc approx.,analytical, gradients)
!
!--------------------------------------------------------------------

      Implicit None

      Real*8 rho,s,omega,Fx_wpbe,d1sfx,d1rfx

      Real*8 f12,f13,f14,f18,f23,f43,f32,f72,f34,f94,f1516,f98
      Real*8 pi,pi2,pi_23,srpi
      Real*8 Three_13

      Real*8 ea1,ea2,ea3,ea4,ea5,ea6,ea7,ea8
      Real*8 eb1
      Real*8 A,B,C,D,E
      Real*8 Ha1,Ha2,Ha3,Ha4,Ha5
      Real*8 Fc1,Fc2
      Real*8 EGa1,EGa2,EGa3
      Real*8 EGscut,wcutoff,expfcutoff

      Real*8 xkf, xkfrho
      Real*8 w,w2,w3,w4,w5,w6,w7,w8
      Real*8 d1rw
      Real*8 A2,A3,A4,A12,A32,A52,A72
      Real*8 X
      Real*8 s2,s3,s4,s5,s6

      Real*8 H,F
      Real*8 Hnum,Hden,d1sHnum,d1sHden
      Real*8 d1sH,d1sF
      Real*8 G_a,G_b,EG
      Real*8 d1sG_a,d1sG_b,d1sEG

      Real*8 Hsbw,Hsbw2,Hsbw3,Hsbw4,Hsbw12,Hsbw32,Hsbw52,Hsbw72
      Real*8 DHsbw,DHsbw2,DHsbw3,DHsbw4,DHsbw5
      Real*8 DHsbw12,DHsbw32,DHsbw52,DHsbw72,DHsbw92
      Real*8 d1sHsbw,d1rHsbw
      Real*8 d1sDHsbw,d1rDHsbw
      Real*8 HsbwA94,HsbwA9412
      Real*8 HsbwA942,HsbwA943,HsbwA945
      Real*8 piexperf,expei
      Real*8 piexperfd1,expeid1
      Real*8 d1spiexperf,d1sexpei
      Real*8 d1rpiexperf,d1rexpei
      Real*8 expei1,expei2,expei3,expei4

      Real*8 DHs,DHs2,DHs3,DHs4,DHs72,DHs92,DHsw,DHsw2,DHsw52,DHsw72
      Real*8 d1sDHs,d1rDHsw

      Real*8 np1,np2
      Real*8 d1rnp1,d1rnp2
      Real*8 t1,t2t9,t10,t10d1
      Real*8 f2,f3,f4,f5,f6,f7,f8,f9
      Real*8 f2d1,f3d1,f4d1,f5d1,f6d1,f8d1,f9d1
      Real*8 d1sf2,d1sf3,d1sf4,d1sf5,d1sf6,d1sf7,d1sf8,d1sf9
      Real*8 d1rf2,d1rf3,d1rf4,d1rf5,d1rf6,d1rf7,d1rf8,d1rf9
      Real*8 d1st1,d1rt1
      Real*8 d1st2t9,d1rt2t9
      Real*8 d1st10,d1rt10
      Real*8 d1sterm1,d1rterm1,term1d1
      Real*8 d1sterm2
      Real*8 d1sterm3,d1rterm3
      Real*8 d1sterm4,d1rterm4
      Real*8 d1sterm5,d1rterm5

      Real*8 term1,term2,term3,term4,term5

      Real*8 ax,um,uk,ul
      Real*8 gc1,gc2

      Real*8, external :: qe_erf, qe_erfc
!      Real*8 ei
      Real*8, external :: expint

      Real*8 Zero,One,Two,Three,Four,Five,Six,Seven,Eight,Nine,Ten
      Real*8 Fifteen,Sixteen
      Real*8 r12,r64,r36,r81,r256,r384,r864,r1944,r4374
      Real*8 r20,r25,r27,r48,r120,r128,r144,r288,r324,r512,r729
      Real*8 r30,r32,r75,r243,r2187,r6561,r40,r105,r54,r135
      Real*8 r1215,r15309

      Save Zero,One,Two,Three,Four,Five,Six,Seven,Eight,Nine,Ten
      Data Zero,One,Two,Three,Four,Five,Six,Seven,Eight,Nine,Ten &
        / 0D0,1D0,2D0,3D0,4D0,5D0,6D0,7D0,8D0,9D0,10D0 /
      Save Fifteen,Sixteen
      Data Fifteen,Sixteen / 1.5D1, 1.6D1 /
      Save r36,r64,r81,r256,r384,r864,r1944,r4374
      Data r36,r64,r81,r256,r384,r864,r1944,r4374 &
        / 3.6D1,6.4D1,8.1D1,2.56D2,3.84D2,8.64D2,1.944D3,4.374D3 /
      Save r27,r48,r120,r128,r144,r288,r324,r512,r729
      Data r27,r48,r120,r128,r144,r288,r324,r512,r729 &
        / 2.7D1,4.8D1,1.2D2,1.28D2,1.44D2,2.88D2,3.24D2,5.12D2,7.29D2 /
      Save r20,r32,r243,r2187,r6561,r40
      Data r20,r32,r243,r2187,r6561,r40 &
        / 2.0d1,3.2D1,2.43D2,2.187D3,6.561D3,4.0d1 /
      Save r12,r25,r30,r54,r75,r105,r135,r1215,r15309
      Data r12,r25,r30,r54,r75,r105,r135,r1215,r15309 &
        / 1.2D1,2.5d1,3.0d1,5.4D1,7.5d1,1.05D2,1.35D2,1.215D3,1.5309D4 /


!     General constants

      f12    = 0.5d0
      f13    = One/Three
      f14    = 0.25d0
      f18    = 0.125d0

      f23    = Two * f13
      f43    = Two * f23

      f32    = 1.5d0
      f72    = 3.5d0
      f34    = 0.75d0
      f94    = 2.25d0
      f98    = 1.125d0
      f1516  = Fifteen / Sixteen

      pi     = ACos(-One)
      pi2    = pi*pi
      pi_23  = pi2**f13
      srpi   = sqrt(pi)

      Three_13 = Three**f13

!     Constants from fit

      ea1 = -1.128223946706117d0
      ea2 = 1.452736265762971d0
      ea3 = -1.243162299390327d0
      ea4 = 0.971824836115601d0
      ea5 = -0.568861079687373d0
      ea6 = 0.246880514820192d0
      ea7 = -0.065032363850763d0
      ea8 = 0.008401793031216d0

      eb1 = 1.455915450052607d0

!     Constants for PBE hole

      A      =  1.0161144d0
      B      = -3.7170836d-1
      C      = -7.7215461d-2
      D      =  5.7786348d-1
      E      = -5.1955731d-2
      X      = - Eight/Nine

!     Constants for fit of H(s) (PBE)

      Ha1    = 9.79681d-3
      Ha2    = 4.10834d-2
      Ha3    = 1.87440d-1
      Ha4    = 1.20824d-3
      Ha5    = 3.47188d-2

!     Constants for F(H) (PBE)

      Fc1    = 6.4753871d0
      Fc2    = 4.7965830d-1

!     Constants for polynomial expansion for EG for small s

      EGa1   = -2.628417880d-2
      EGa2   = -7.117647788d-2
      EGa3   =  8.534541323d-2

!     Constants for large x expansion of exp(x)*ei(-x)

      expei1 = 4.03640D0
      expei2 = 1.15198D0
      expei3 = 5.03627D0
      expei4 = 4.19160D0

!     Cutoff criterion below which to use polynomial expansion

      EGscut     = 8.0d-2
      wcutoff    = 1.4D1
      expfcutoff = 7.0D2

!     Calculate prelim variables

      xkf    = (Three*pi2*rho) ** f13
      xkfrho = xkf * rho

      A2 = A*A
      A3 = A2*A
      A4 = A3*A
      A12 = Sqrt(A)
      A32 = A12*A
      A52 = A32*A
      A72 = A52*A
      w      = omega / xkf
      w2    = w * w
      w3    = w2 * w
      w4    = w2 * w2
      w5    = w3 * w2
      w6    = w5 * w
      w7    = w6 * w
      w8    = w7 * w

      d1rw  = -(One/(Three*rho))*w

      X      = - Eight/Nine

      s2     = s*s
      s3     = s2*s
      s4     = s2*s2
      s5     = s4*s
      s6     = s5*s

!     Calculate wPBE enhancement factor

      Hnum    = Ha1*s2 + Ha2*s4
      Hden    = One + Ha3*s4 + Ha4*s5 + Ha5*s6

      H       = Hnum/Hden

      d1sHnum = Two*Ha1*s + Four*Ha2*s3
      d1sHden = Four*Ha3*s3 + Five*Ha4*s4 + Six*Ha5*s5

      d1sH    = (Hden*d1sHnum - Hnum*d1sHden) / (Hden*Hden)

      F      = Fc1*H + Fc2
      d1sF   = Fc1*d1sH

!     Change exponent of Gaussian if we're using the simple approx.

      if(w .gt. wcutoff) then

        eb1 = 2.0d0

      endif

!     Calculate helper variables (should be moved later on...)

      Hsbw = s2*H + eb1*w2
      Hsbw2 = Hsbw*Hsbw
      Hsbw3 = Hsbw2*Hsbw
      Hsbw4 = Hsbw3*Hsbw
      Hsbw12 = Sqrt(Hsbw)
      Hsbw32 = Hsbw12*Hsbw
      Hsbw52 = Hsbw32*Hsbw
      Hsbw72 = Hsbw52*Hsbw

      d1sHsbw  = d1sH*s2 + Two*s*H
      d1rHsbw  = Two*eb1*d1rw*w

      DHsbw = D + s2*H + eb1*w2
      DHsbw2 = DHsbw*DHsbw
      DHsbw3 = DHsbw2*DHsbw
      DHsbw4 = DHsbw3*DHsbw
      DHsbw5 = DHsbw4*DHsbw
      DHsbw12 = Sqrt(DHsbw)
      DHsbw32 = DHsbw12*DHsbw
      DHsbw52 = DHsbw32*DHsbw
      DHsbw72 = DHsbw52*DHsbw
      DHsbw92 = DHsbw72*DHsbw

      HsbwA94   = f94 * Hsbw / A
      HsbwA942  = HsbwA94*HsbwA94
      HsbwA943  = HsbwA942*HsbwA94
      HsbwA945  = HsbwA943*HsbwA942
      HsbwA9412 = Sqrt(HsbwA94)

      DHs    = D + s2*H
      DHs2   = DHs*DHs
      DHs3   = DHs2*DHs
      DHs4   = DHs3*DHs
      DHs72  = DHs3*sqrt(DHs)
      DHs92  = DHs72*DHs

      d1sDHs = Two*s*H + s2*d1sH

      DHsw   = DHs + w2
      DHsw2  = DHsw*DHsw
      DHsw52 = sqrt(DHsw)*DHsw2
      DHsw72 = DHsw52*DHsw

      d1rDHsw = Two*d1rw*w

      if(s .gt. EGscut) then

        G_a    = srpi * (Fifteen*E + Six*C*(One+F*s2)*DHs + &
                         Four*B*(DHs2) + Eight*A*(DHs3))    &
                      * (One / (Sixteen * DHs72))           &
                       - f34*pi*sqrt(A) * exp(f94*H*s2/A) * &
                         (One - qe_erf(f32*s*sqrt(H/A)))

        d1sG_a = (One/r32)*srpi *                           &
                 ((r36*(Two*H + d1sH*s) / (A12*sqrt(H/A)))  &
                  + (One/DHs92) *                           &
                     (-Eight*A*d1sDHs*DHs3 - r105*d1sDHs*E  &
                      -r30*C*d1sDHs*DHs*(One+s2*F)          &
                      +r12*DHs2*(-B*d1sDHs + C*s*(d1sF*s + Two*F)))  &
                  - ((r54*exp(f94*H*s2/A)*srpi*s*(Two*H+d1sH*s)*     &
                     qe_erfc(f32*sqrt(H/A)*s))                         &
                     / A12))

        G_b    = (f1516 * srpi * s2) / DHs72

        d1sG_b = (Fifteen*srpi*s*(Four*DHs - Seven*d1sDHs*s)) &
                 / (r32*DHs92)

        EG     = - (f34*pi + G_a) / G_b

        d1sEG  = (-Four*d1sG_a*G_b + d1sG_b*(Four*G_a + Three*pi)) &
                 / (Four*G_b*G_b)

      else

        EG    = EGa1 + EGa2*s2 + EGa3*s4
        d1sEG = Two*EGa2*s + Four*EGa3*s3

      endif

!    Calculate the terms needed in any case


      term2 =       (DHs2*B + DHs*C + Two*E + DHs*s2*C*F + Two*s2*EG) / &
                    (Two*DHs3)

      d1sterm2 = (-Six*d1sDHs*(EG*s2 + E)                     &
                  + DHs2 * (-d1sDHs*B + s*C*(d1sF*s + Two*F)) &
                  + Two*DHs * (Two*EG*s - d1sDHs*C            &
                  + s2 * (d1sEG - d1sDHs*C*F)))               &
                 / (Two*DHs4)

      term3 = - w  * (Four*DHsw2*B + Six*DHsw*C + Fifteen*E &
                      + Six*DHsw*s2*C*F + Fifteen*s2*EG) /  &
                     (Eight*DHs*DHsw52)

      d1sterm3 = w * (Two*d1sDHs*DHsw * (Four*DHsw2*B         &
                         + Six*DHsw*C + Fifteen*E             &
                         + Three*s2*(Five*EG + Two*DHsw*C*F)) &
                      + DHs * (r75*d1sDHs*(EG*s2 + E)         &
                         + Four*DHsw2*(d1sDHs*B               &
                              - Three*s*C*(d1sF*s + Two*F))   &
                         - Six*DHsw*(-Three*d1sDHs*C          &
                              + s*(Ten*EG + Five*d1sEG*s      &
                                  - Three*d1sDHs*s*C*F))))    &
                 / (Sixteen*DHs2*DHsw72)

      d1rterm3 = (-Two*d1rw*DHsw * (Four*DHsw2*B              &
                         + Six*DHsw*C + Fifteen*E             &
                         + Three*s2*(Five*EG + Two*DHsw*C*F)) &
                      + w * d1rDHsw * (r75*(EG*s2 + E)        &
                         + Two*DHsw*(Two*DHsw*B + Nine*C      &
                                     + Nine*s2*C*F)))         &
                 / (Sixteen*DHs*DHsw72)

      term4 = - w3 * (DHsw*C + Five*E + DHsw*s2*C*F + Five*s2*EG) /  &
                     (Two*DHs2*DHsw52)

      d1sterm4 = (w3 * (Four*d1sDHs*DHsw * (DHsw*C + Five*E   &
                             + s2 * (Five*EG + DHsw*C*F))     &
                        + DHs * (r25*d1sDHs*(EG*s2 + E)       &
                             - Two*DHsw2*s*C*(d1sF*s + Two*F) &
                             + DHsw * (Three*d1sDHs*C + s*(-r20*EG  &
                                   - Ten*d1sEG*s              &
                                   + Three*d1sDHs*s*C*F)))))  &
                 / (Four*DHs3*DHsw72)

      d1rterm4 = (w2 * (-Six*d1rw*DHsw * (DHsw*C + Five*E   &
                             + s2 * (Five*EG + DHsw*C*F))   &
                        + w * d1rDHsw * (r25*(EG*s2 + E) +  &
                             Three*DHsw*C*(One + s2*F))))  &
                 / (Four*DHs2*DHsw72)

      term5 = - w5 * (E + s2*EG) / &
                     (DHs3*DHsw52)

      d1sterm5 = (w5 * (Six*d1sDHs*DHsw*(EG*s2 + E)               &
                        + DHs * (-Two*DHsw*s * (Two*EG + d1sEG*s) &
                             + Five*d1sDHs * (EG*s2 + E))))       &
                 / (Two*DHs4*DHsw72)

      d1rterm5 = (w4 * Five*(EG*s2 + E) * (-Two*d1rw*DHsw   &
                                           + d1rDHsw * w))  &
                 / (Two*DHs3*DHsw72)


      if((s.gt.0.0d0).or.(w.gt.0.0d0)) then

        t10    = (f12)*A*Log(Hsbw / DHsbw)
        t10d1  = f12*A*(One/Hsbw - One/DHsbw)
        d1st10 = d1sHsbw*t10d1
        d1rt10 = d1rHsbw*t10d1

      endif

!     Calculate exp(x)*f(x) depending on size of x

      if(HsbwA94 .lt. expfcutoff) then

        piexperf = pi*Exp(HsbwA94)*qe_erfc(HsbwA9412)
!        expei    = Exp(HsbwA94)*Ei(-HsbwA94)
        expei    = Exp(HsbwA94)*(-expint(1,HsbwA94))

      else

!       print *,rho,s," LARGE HsbwA94"

        piexperf = pi*(One/(srpi*HsbwA9412)          &
                   - One/(Two*Sqrt(pi*HsbwA943))     &
                   + Three/(Four*Sqrt(pi*HsbwA945)))

        expei  = - (One/HsbwA94) *                         &
                   (HsbwA942 + expei1*HsbwA94 + expei2) /  &
                   (HsbwA942 + expei3*HsbwA94 + expei4)

      endif

!     Calculate the derivatives (based on the orig. expression)
!     --> Is this ok? ==> seems to be ok...

      piexperfd1  = - (Three*srpi*sqrt(Hsbw/A))/(Two*Hsbw)  &
                    + (Nine*piexperf)/(Four*A)
      d1spiexperf = d1sHsbw*piexperfd1
      d1rpiexperf = d1rHsbw*piexperfd1

      expeid1  = f14*(Four/Hsbw + (Nine*expei)/A)
      d1sexpei = d1sHsbw*expeid1
      d1rexpei = d1rHsbw*expeid1

      if (w .eq. Zero) then

!       Fall back to original expression for the PBE hole

        t1 = -f12*A*expei
        d1st1 = -f12*A*d1sexpei
        d1rt1 = -f12*A*d1rexpei

!        write(*,*) s, t1, t10, d1st1,d1rt1,d1rt10

        if(s .gt. 0.0D0) then

          term1    = t1 + t10
          d1sterm1 = d1st1 + d1st10
          d1rterm1 = d1rt1 + d1rt10

          Fx_wpbe = X * (term1 + term2)

          d1sfx = X * (d1sterm1 + d1sterm2)
          d1rfx = X * d1rterm1

        else

          Fx_wpbe = 1.0d0

! TODO    This is checked to be true for term1
!         How about the other terms???

          d1sfx   = 0.0d0
          d1rfx   = 0.0d0

        endif


      elseif(w .gt. wcutoff) then

!       Use simple Gaussian approximation for large w

!       print *,rho,s," LARGE w"

        term1 = -f12*A*(expei+log(DHsbw)-log(Hsbw))

        term1d1  = - A/(Two*DHsbw) - f98*expei
        d1sterm1 = d1sHsbw*term1d1
        d1rterm1 = d1rHsbw*term1d1

        Fx_wpbe = X * (term1 + term2 + term3 + term4 + term5)

        d1sfx = X * (d1sterm1 + d1sterm2 + d1sterm3  &
                              + d1sterm4 + d1sterm5)

        d1rfx = X * (d1rterm1 + d1rterm3 + d1rterm4 + d1rterm5)

      else

!       For everything else, use the full blown expression

!       First, we calculate the polynomials for the first term

         np1    = -f32*ea1*A12*w + r27*ea3*w3/(Eight*A12)     &
                  - r243*ea5*w5/(r32*A32) + r2187*ea7*w7/(r128*A52)

        d1rnp1 = - f32*ea1*d1rw*A12 + (r81*ea3*d1rw*w2)/(Eight*A12) &
                 - (r1215*ea5*d1rw*w4)/(r32*A32)                    &
                 + (r15309*ea7*d1rw*w6)/(r128*A52)

        np2 = -A + f94*ea2*w2 - r81*ea4*w4/(Sixteen*A)        &
              + r729*ea6*w6/(r64*A2) - r6561*ea8*w8/(r256*A3)


        d1rnp2 =   f12*(Nine*ea2*d1rw*w)         &
                 - (r81*ea4*d1rw*w3)/(Four*A)    &
                 + (r2187*ea6*d1rw*w5)/(r32*A2)  &
                 - (r6561*ea8*d1rw*w7)/(r32*A3)

!       The first term is

        t1    = f12*(np1*piexperf + np2*expei)
        d1st1 = f12*(d1spiexperf*np1 + d1sexpei*np2)
        d1rt1 = f12*(d1rnp2*expei + d1rpiexperf*np1 +  &
                     d1rexpei*np2 + d1rnp1*piexperf)

!       The factors for the main polynomoal in w and their derivatives

        f2    = (f12)*ea1*srpi*A / DHsbw12
        f2d1  = - ea1*srpi*A / (Four*DHsbw32)
        d1sf2 = d1sHsbw*f2d1
        d1rf2 = d1rHsbw*f2d1

        f3    = (f12)*ea2*A / DHsbw
        f3d1  = - ea2*A / (Two*DHsbw2)
        d1sf3 = d1sHsbw*f3d1
        d1rf3 = d1rHsbw*f3d1

        f4    =  ea3*srpi*(-f98 / Hsbw12     &
                 + f14*A / DHsbw32)
        f4d1  = ea3*srpi*((Nine/(Sixteen*Hsbw32))-   &
                          (Three*A/(Eight*DHsbw52)))
        d1sf4 = d1sHsbw*f4d1
        d1rf4 = d1rHsbw*f4d1

        f5    = ea4*(One/r128) * (-r144*(One/Hsbw)   &
                 + r64*(One/DHsbw2)*A)
        f5d1  = ea4*((f98/Hsbw2)-(A/DHsbw3))
        d1sf5 = d1sHsbw*f5d1
        d1rf5 = d1rHsbw*f5d1

        f6    = ea5*(Three*srpi*(Three*DHsbw52*(Nine*Hsbw-Two*A) &
                 + Four*Hsbw32*A2))                              &
                 / (r32*DHsbw52*Hsbw32*A)
        f6d1  = ea5*srpi*((r27/(r32*Hsbw52))-        &
                    (r81/(r64*Hsbw32*A))-            &
                    ((Fifteen*A)/(Sixteen*DHsbw72)))
        d1sf6 = d1sHsbw*f6d1
        d1rf6 = d1rHsbw*f6d1

        f7    = ea6*(((r32*A)/DHsbw3                 &
                 + (-r36 + (r81*s2*H)/A)/Hsbw2)) / r32
        d1sf7 = ea6*(Three*(r27*d1sH*DHsbw4*Hsbw*s2 +           &
                Eight*d1sHsbw*A*(Three*DHsbw4 - Four*Hsbw3*A) + &
                r54*DHsbw4*s*(Hsbw - d1sHsbw*s)*H))/            &
                (r32*DHsbw4*Hsbw3*A)
        d1rf7 = ea6*d1rHsbw*((f94/Hsbw3)-((Three*A)/DHsbw4)     &
                           -((r81*s2*H)/(Sixteen*Hsbw3*A)))

        f8    = ea7*(-Three*srpi*(-r40*Hsbw52*A3                &
                 +Nine*DHsbw72*(r27*Hsbw2-Six*Hsbw*A+Four*A2))) &
                 / (r128 * DHsbw72*Hsbw52*A2)
        f8d1  = ea7*srpi*((r135/(r64*Hsbw72)) + (r729/(r256*Hsbw32*A2))  &
                         -(r243/(r128*Hsbw52*A))                         &
                         -((r105*A)/(r32*DHsbw92)))
        d1sf8 = d1sHsbw*f8d1
        d1rf8 = d1rHsbw*f8d1

        f9    = (r324*ea6*eb1*DHsbw4*Hsbw*A                      &
                + ea8*(r384*Hsbw3*A3 + DHsbw4*(-r729*Hsbw2       &
                + r324*Hsbw*A - r288*A2))) / (r128*DHsbw4*Hsbw3*A2)
        f9d1  = -((r81*ea6*eb1)/(Sixteen*Hsbw3*A))               &
                + ea8*((r27/(Four*Hsbw4))+(r729/(r128*Hsbw2*A2)) &
                      -(r81/(Sixteen*Hsbw3*A))                   &
                      -((r12*A/DHsbw5)))
        d1sf9 = d1sHsbw*f9d1
        d1rf9 = d1rHsbw*f9d1

        t2t9    = f2*w  + f3*w2 + f4*w3 + f5*w4 + f6*w5          &
                        + f7*w6 + f8*w7 + f9*w8
        d1st2t9 = d1sf2*w + d1sf3*w2 + d1sf4*w3 + d1sf5*w4       &
                          + d1sf6*w5 + d1sf7*w6 + d1sf8*w7       &
                          + d1sf9*w8
        d1rt2t9 = d1rw*f2 + d1rf2*w + Two*d1rw*f3*w   &
                  + d1rf3*w2 + Three*d1rw*f4*w2       &
                  + d1rf4*w3 + Four*d1rw*f5*w3        &
                  + d1rf5*w4 + Five*d1rw*f6*w4        &
                  + d1rf6*w5 + Six*d1rw*f7*w5         &
                  + d1rf7*w6 + Seven*d1rw*f8*w6       &
                  + d1rf8*w7 + Eight*d1rw*f9*w7 + d1rf9*w8

!       The final value of term1 for 0 < omega < wcutoff is:

        term1 = t1 + t2t9 + t10

        d1sterm1 = d1st1 + d1st2t9 + d1st10
        d1rterm1 = d1rt1 + d1rt2t9 + d1rt10

!       The final value for the enhancement factor and its
!       derivatives is:

        Fx_wpbe = X * (term1 + term2 + term3 + term4 + term5)

        d1sfx = X * (d1sterm1 + d1sterm2 + d1sterm3    &
                              + d1sterm4 + d1sterm5)

        d1rfx = X * (d1rterm1 + d1rterm3 + d1rterm4 + d1rterm5)

      endif

      END SUBROUTINE wpbe_analy_erfc_approx_grad


