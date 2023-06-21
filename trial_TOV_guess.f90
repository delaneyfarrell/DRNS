!--------------------------------------------------------------------------
! solve TOV equation for spherically symmetric stars
!--------------------------------------------------------------------------

      subroutine TOV_guess

          use vars

          implicit none
          integer:: i, s, m, s_temp
          double precision :: r_is_s, lambda_s, nu_s, gamma_eq, rho_eq
          double precision :: interp 

          call integrate(1) 
          call integrate(2) 
          call integrate(3) 

          if (s_max .eq. 1.) then

              s_temp = s_div-1

          else 

              s_temp = s_div

          end if

          n_nearest = r_div/2

          do s = 1, s_temp 

            r_is_s = r_is_final*(s_gp(s)/(1.-s_gp(s)))

            ! transform spherical solutions to the s radial coordinate 
            if (r_is_s .lt. r_is_final) then

                lambda_s = interp(r_is_gp, lambda_gp, r_div, r_is_s)
                nu_s = interp(r_is_gp, nu_gp, r_div, r_is_s) 

            else 

                lambda_s = 2.*log(1.+m_final/(2.*r_is_s))
                nu_s = log((1.-m_final/(2.*r_is_s))/(1.+m_final/ & 
                    (2*r_is_s)))

            end if

            gam(s,1) = nu_s+lambda_s
            rho(s,1) = nu_s-lambda_s

            ! tranform to metrics gamma & rho
            do m = 1, m_div

                gamma_guess(s,m) = gam(s,1)
                rho_guess(s,m) = rho(s,1)
                alpha_guess(s,m) = (gam(s,1)-rho(s,1))/2.
                omega_guess(s,m) = 0.d0

            end do

            gamma_mu_0(s) = gam(s,1)
            rho_mu_0(s) = rho(s,1) 

          end do

          if (s_max .eq. 1.) then 

              do m = 1, m_div

                gamma_guess(s_div,m) = 0.d0
                rho_guess(s_div,m) = 0.d0
                alpha_guess(s_div,m) = 0.d0
                omega_guess(s_div,m) = 0.d0

              end do

          end if 

          n_nearest = s_div/2
          gamma_eq = interp(s_gp, gamma_mu_0, s_div, s_e) 
          rho_eq = interp(s_gp, rho_mu_0, s_div, s_e) 
          r_e_guess = r_final*exp(0.5*(rho_eq-gamma_eq)) 
          r_ec = r_final*sqrt(kappa)
          mass = m_final*sqrt(kappa)*(c*c/G) 
          z_p = exp(-0.5*(gamma_eq+rho_eq))-1.

          print *, 'TOV calc, r_e_guess=', r_e_guess

      end subroutine 

!--------------------------------------------------------------------------
! integration of TOV subroutine
!--------------------------------------------------------------------------

      subroutine integrate(n) 

          use vars

          implicit none
          integer :: i, n
          double precision :: r, r_is, m, p, e_d, r_is_est, e_at_p
          double precision :: dr_is_save, r_is_check, nu_s, h, hh
          double precision :: a1, a2, a3, a4, b1, b2, b3, b4
          double precision :: c1, c2, c3, c4, rho0, dr_dr_is
          double precision :: dp_dr_is, dm_dr_is, rtsec_c
          double precision :: p_at_e, h_at_p 

          if (n .eq. 1) then 
              
              r_is_est = 1.5e6/sqrt(kappa)
              h = r_is_est/100.d0

          else

              r_is_est = r_is_final
              h = r_is_est/10000.d0

          end if

          dr_is_save = r_is_final/r_div
          r_is_check = dr_is_save

          ! initialize (is = isotropic)
          r_is = 0.
          r = 0.
          m = 0.
          p = p_center

          r_is_gp(1) = 0.d0
          r_gp(1) = 0.d0
          m_gp(1) = 0.d0
          lambda_gp(1) = 0.d0
          e_d_gp(1) = e_center 

          i = 2
          do while (p .ge. p_surface) 

            e_d = e_at_p(p) 

            if (n .eq. 3 .and. (r_is .gt. r_is_check) .and. &
                (i .lt. r_div)) then 

                r_is_gp(i) = r_is
                r_gp(i) = r
                m_gp(i) = m
                e_d_gp(i) = e_d
                r_is_check = r_is_check+dr_is_save

                i = i+1

            end if

            r_is_final = r_is
            r_final = r
            m_final = m

            a1 = dr_dr_is(r_is,r,m)
            b1 = dm_dr_is(r_is,r,m,p)
            c1 = dp_dr_is(r_is,r,m,p) 

            a2 = dr_dr_is(r_is+h/2., r+h*a1/2., m+h*b1/2.)
            b2 = dm_dr_is(r_is+h/2., r+h*a1/2., m+h*b1/2., p+h*c1/2.)
            c2 = dp_dr_is(r_is+h/2., r+h*a1/2., m+h*b1/2., p+h*c1/2.)

            a3 = dr_dr_is(r_is+h/2., r+h*a2/2., m+h*b2/2.)
            b3 = dm_dr_is(r_is+h/2., r+h*a2/2., m+h*b2/2., p+h*c2/2.)
            c3 = dp_dr_is(r_is+h/2., r+h*a2/2., m+h*b2/2., p+h*c2/2.) 

            a4 = dr_dr_is(r_is+h, r+h*a3, m+h*b3)
            b4 = dm_dr_is(r_is+h, r+h*a3, m+h*b3, p+h*c3)
            c4 = dp_dr_is(r_is+h, r+h*a3, m+h*b3, p+h*c3) 

            r = r + (h/6.)*(a1+2*a2+2*a3+a4)
            m = m + (h/6.)*(b1+2*b2+2*b3+b4)
            p = p + (h/6.)*(c1+2*c2+2*c3+c4) 

            r_is = r_is+h

          end do

          r_is_gp(r_div) = r_is_final
          r_gp(r_div) = r_final
          m_gp(r_div) = m_final 

          ! rescale r_is & compute lambda 
          if (n .eq. 3) then 

              k_rescale = 0.5*(r_final/r_is_final)*(1.-m_final/r_final &
                  +sqrt(1.-2.*m_final/r_final)) 

              r_is_final = r_is_final*k_rescale

              nu_s = log((1.-m_final/(2.*r_is_final))/(1.+m_final/(2.* &
                  r_is_final)))

              do i = 1, r_div

                r_is_gp(i) = r_is_gp(i)*k_rescale

                if (i .eq. 1) then

                    lambda_gp(i) = log(1./k_rescale) 

                else 

                    lambda_gp(i) = log(r_gp(i)/r_is_gp(i)) 
                    
                end if

                if (e_d_gp(i) .lt. e_surface) then

                    hh = 0.d0 

                else 

                    p = p_at_e(e_d_gp(i)) 
                    hh = h_at_p(p) 

                end if

                nu_gp(i) = nu_s-hh

              end do

              nu_gp(r_div) = nu_s 

          end if

      end subroutine

!--------------------------------------------------------------------------
! dr_dr_is function
!--------------------------------------------------------------------------

      double precision function dr_dr_is(r_is, r, m)

          use vars

          implicit none
          double precision :: r_is, r, m

          if (r_is .lt. r_min) then

              dr_dr_is = 1.

          else 

              dr_dr_is = (r/r_is)*sqrt(1.-2.*m/r) 

          end if

      end function 

!--------------------------------------------------------------------------
! dm_dr_is function
!--------------------------------------------------------------------------

      double precision function dm_dr_is(r_is, r, m, p) 

          use vars

          implicit none
          double precision :: r_is, r, m, p, e_d, e_at_p

          if (p .lt. p_surface) then 

              e_d = 0.d0 

          else 

              e_d = e_at_p(p) 

          end if

          if (r_is .lt. r_min) then

              dm_dr_is = 4.*pi*e_center*r*r*(1.+4*pi*e_center*r*r/3.) 

          else 

              dm_dr_is = 4.*pi*e_d*r*r*r*sqrt(1.-2.*m/r)/r_is

          end if

      end function 

!--------------------------------------------------------------------------
! dp_dr_is function 
!--------------------------------------------------------------------------

      double precision function dp_dr_is(r_is, r, m, p)

          use vars

          implicit none
          double precision :: r_is, r, m, p, e_d, e_at_p

          if (p .lt. p_surface) then 

              e_d = 0.d0 

          else 

              e_d = e_at_p(p) 

          end if

          if (r_is .lt. r_min) then 

              dp_dr_is = -4.*pi*(e_center+p)*(e_center+3.*p)*r* & 
                  (1.+4.*e_center*r*r/3.)/3.

          else 

              dp_dr_is = -(e_d+p)*(m+4.*pi*r*r*r*p)/(r*r_is* &
                  sqrt(1.-2.*m/r))

          end if

      end function


