!--------------------------------------------------------------------------
! computation of physical quantities
!--------------------------------------------------------------------------

      subroutine comp_phys

          use vars 

          implicit none 

          integer :: s_surf_1, s_temp, s, m, i 

          double precision :: velocity(0:s_div,0:m_div), gamma_pole
          double precision :: gamma_equator, rho_pole, rho_equator 
          double precision :: omega_equator, d_m(0:s_div)
          double precision :: d_m_0(0:s_div), d_m_p(0:s_div)
          double precision :: d_j(0:s_div), rho_0(0:s_div+1,0:m_div) 
          double precision :: d_o_e(0:s_div), d_g_e(0:s_div)
          double precision :: d_r_e(0:s_div), d_v_e(0:s_div) 
          double precision :: doe, dge, dre, dve, vek, d_gamma_s
          double precision :: d_rho_s, d_omega_s, d_gamma_m
          double precision :: d_rho_m, d_omega_m, d_alpha_s
          double precision :: d_alpha_m, sqrt_v, s1, s_1, s_plus
          double precision :: s_minus, r_plus, gamma_plus, rho_plus
          double precision :: gamma_minus, rho_minus, b_st_p(0:s_div)
          double precision :: b_st_m(0:s_div), d_v_plus_s, d_v_minus_s
          double precision :: gamma_m(0:s_div), rho_m(0:s_div)
          double precision :: alpha_m(0:s_div), enthaply_m(0:s_div) 
          double precision :: r_surface(0:m_div), s_surface(0:m_div) 
          double precision :: dpi_m, d_r_s_dm, d_gamma_s_dm
          double precision :: d_rho_s_dm, srt_1_m2, f_mu(0:m_div)
          double precision :: virial1, virial2, virial3
          double precision :: s_virial1(0:s_div,0:m_div) 
          double precision :: s_virial2(0:s_div,0:m_div)
          double precision :: s_virial3(0:s_div,0:m_div)
          double precision :: u_t, u_t_plus, omega_plus, vel_p
          double precision :: d_virial1(0:s_div), d_virial2(0:s_div) 
          double precision :: m1, temp1, temp2, temp3, temp4, temp5
          double precision :: temp6, temp_x(5), temp_y1(5)
          double precision :: temp_y2(5), s_ad1(0:s_div,4)
          double precision :: s_ad2(0:s_div,4), mu_ch(0:m_div)
          double precision :: t_rho(0:m_div), t_alpha(0:m_div) 
          double precision :: t_rho_s(0:m_div), t_gamma_s(0:m_div)
          double precision :: t_rho_m(0:m_div), t_gamma_m(0:m_div) 
          double precision :: t_omega_s(0:m_div), t_omega_m(0:m_div) 
          double precision :: t_pressure(0:m_div), t_energy(0:m_div) 
          double precision :: t_v2(0:m_div), rho_ch, alpha_ch
          double precision :: rho_s_ch, gamma_s_ch, gamma_m_ch
          double precision :: rho_m_ch, omega_s_ch, omega_m_ch
          double precision :: pressure_ch, energy_ch, v2_ch
          double precision :: sv1_ch(0:s_div,0:m_div)
          double precision :: sv2_ch(0:s_div,0:m_div)
          double precision :: dv1_ch(0:s_div), dv2_ch(0:s_div) 
          double precision :: gamma_surface(0:m_div)
          double precision :: rho_surface(0:m_div) 
          double precision :: alpha_surface(0:m_div), pi_bar(0:m_div)
          double precision :: z_emb(0:m_div), grv2_spherical
          double precision :: b_st_p_surface, b_st_m_surface
          double precision :: b_st_p_out(0:s_div-(s_div/2)+2) 
          double precision :: b_st_m_out(0:s_div-(s_div/2)+2)
          double precision :: s_gp_out(0:s_div-(s_div/2)+2)
          double precision :: interp, deriv_s, n0_at_e, sq

          ! polar calulations (radius & s-coordinate) 
          r_p = r_ratio*r_e
          s_p = r_p/(r_p+r_e)
          s_e = 0.5

          do s = 1, s_div
            do m = 1, m_div
                velocity(s,m) = sqrt(velocity_sq(s,m))
            end do 
          end do


          do s = 1, s_div
            gamma_mu_1(s) = gam(s,m_div) 
            gamma_mu_0(s) = gam(s,1)
            rho_mu_1(s) = rho(s,m_div) 
            rho_mu_0(s) = rho(s,1)
            omega_mu_0(s) = omega(s,1) 
          end do

          n_nearest = s_div/2
          gamma_pole = interp(s_gp,gamma_mu_1,s_div,s_p)
          gamma_equator = interp(s_gp,gamma_mu_0,s_div,s_e) 
          rho_pole = interp(s_gp, rho_mu_1,s_div,s_p) 
          rho_equator = interp(s_gp, rho_mu_0,s_div,s_e) 

          if (r_ratio .eq. 1.) then 
              velocity_equator = 0.d0
              omega_equator = 0.d0
          else 
              velocity_equator = sqrt(1.-exp(gamma_pole+rho_pole- & 
                  gamma_equator-rho_equator+SQ(Ac)*SQ(Omega_rot_eq &
                  -Omega_rot_c)))
              omega_equator = interp(s_gp,omega_mu_0,s_div,s_e) 

          end if

          ! circumferential radius
          r_ec = sqrt(kappa)*r_e*exp((gamma_equator-rho_equator)/2.)

          mass = 0.
          mass_0 = 0.
          mass_p = 0.
          jam = 0.
          n_nearest = n_tab/2

          do s = 1, s_div

            do m = 1, m_div

                if (energy(s,m) .gt. e_surface) then

                    rho_0(s,m) = n0_at_e(energy(s,m))*mb*c*c*k_scale

                else 

                    rho_0(s,m) = 0.d0

                end if

            end do

          end do

          if (s_max .eq. 1.) then

              s_temp = s_div-1

          else 

              s_temp = s_div

          end if

          do s = 1, s_temp

            d_m(s) = 0.
            d_m_0(s) = 0.
            d_m_p(s) = 0.
            d_j(s) = 0. 

            do m = 1, m_div-2, 2

            d_m(s) = d_m(s) + (1./(3.*(m_div-1)))*(exp(2.*alpha(s,m)+ &
                gam(s,m))*(((energy(s,m)+pressure(s,m))/(1.- &
                velocity_sq(s,m)))*(1.+velocity_sq(s,m)+(2.*s_gp(s)* &
                sqrt(velocity_sq(s,m))/(1.-s_gp(s)))*sqrt(1.-mu(m)* &
                mu(m))*r_e*omega(s,m)*exp(-rho(s,m)))+2.*pressure(s,m))&
                +4.*exp(2.*alpha(s,m+1)+gam(s,m+1))* &
                (((energy(s,m+1)+pressure(s,m+1))/(1.- &
                velocity_sq(s,m+1)))*(1.+velocity_sq(s,m+1)+(2.*s_gp(s)&
                *sqrt(velocity_sq(s,m+1))/(1.-s_gp(s)))*sqrt(1.-mu(m+1)&
                *mu(m+1))*r_e*omega(s,m+1)*exp(-rho(s,m+1)))+2.* &
                pressure(s,m+1))+exp(2.*alpha(s,m+2)+gam(s,m+1))* &
                (((energy(s,m+2)+pressure(s,m+2))/(1.- &
                velocity_sq(s,m+2)))*(1+velocity_sq(s,m+2)+(2.* &
                s_gp(s)*sqrt(velocity_sq(s,m+2))/(1-s_gp(s)))*sqrt(1.- &
                mu(m+2)*mu(m+2))*r_e*omega(s,m+2)*exp(-rho(s,m+2)))+ &
                2.*pressure(s,m+2)))

            d_m_0(s) = d_m_0(s) + (1./(3.*(m_div-1)))*(exp(2.* & 
                alpha(s,m)+(gam(s,m)-rho(s,m))/2.)*rho_0(s,m)/ &
                sqrt(1.-velocity_sq(s,m))+4.* & 
                exp(2.*alpha(s,m+1)+(gam(s,m+1) & 
                -rho(s,m+1))/2.)*rho_0(s,m+1)/sqrt(1.- &
                velocity_sq(s,m+1))+exp(2.*alpha(s,m+2)+(gam(s,m+2) &
                -rho(s,m+2))/2.)*rho_0(s,m+2)/sqrt(1.- &
                velocity_sq(s,m+2)))

            d_m_p(s) = d_m_p(s) + (1./(3.*(m_div-1)))*(exp(2.* &
                alpha(s,m)+(gam(s,m)-rho(s,m))/2.)*energy(s,m)/ &
                sqrt(1.-velocity_sq(s,m))+4.* &
                exp(2.*alpha(s,m+1)+(gam(s,m+1) &
                -rho(s,m+1))/2.)*energy(s,m+1)/sqrt(1.- &
                velocity_sq(s,m+1))+exp(2.*alpha(s,m+2)+(gam(s,m+2) &
                -rho(s,m+2))/2.)*energy(s,m+2)/sqrt(1.- &
                velocity_sq(s,m+2)))

            d_j(s) = d_j(s) + (1./(3.*(m_div-1)))*(sqrt(1.-mu(m)* &
                mu(m))*exp(2.*alpha(s,m)+gam(s,m)-rho(s,m))* &
                (energy(s,m)+pressure(s,m))*sqrt(velocity_sq(s,m))/ &
                (1.-velocity_sq(s,m))+4.*sqrt(1.-mu(m+1)*mu(m+1))* &
                exp(2.*alpha(s,m+1)+gam(s,m+1)-rho(s,m+1))* &
                (energy(s,m+1)+pressure(s,m+1))* &
                sqrt(velocity_sq(s,m+1))/(1.-velocity_sq(s,m+1))+ &
                sqrt(1.-mu(m+2)*mu(m+2))*exp(2.*alpha(s,m+2)+ &
                gam(s,m+2)-rho(s,m+2))*(energy(s,m+2)+ &
                pressure(s,m+2))*sqrt(velocity_sq(s,m+2))/(1.- &
                velocity_sq(s,m+2)))

            end do

          end do

          if (s_max .eq. 1.) then

              d_m(s_div) = 0.
              d_m_0(s_div) = 0.
              d_m_p(s_div) = 0.
              d_j(s_div) = 0.

          end if

          do s = 1, s_div-4, 2

            mass = mass + (s_max/(3.*(s_div-1)))*(((sqrt(s_gp(s))/ &
                (1.-s_gp(s)))**4.)*d_m(s)+4.*((sqrt(s_gp(s+1))/ &
                (1.-s_gp(s+1)))**4.)*d_m(s+1)+((sqrt(s_gp(s+2))/ & 
                (1.-s_gp(s+2)))**4.)*d_m(s+2))

            mass_0 = mass_0 + (s_max/(3.*(s_div-1)))*(((sqrt(s_gp(s))/ &
                (1.-s_gp(s)))**4.)*d_m_0(s)+4.*((sqrt(s_gp(s+1))/ &
                (1.-s_gp(s+1)))**4.)*d_m_0(s+1)+((sqrt(s_gp(s+2))/ &
                (1.-s_gp(s+2)))**4.)*d_m_0(s+2))

            mass_p = mass_p + (s_max/(3.*(s_div-1)))*(((sqrt(s_gp(s))/ &
                (1.-s_gp(s)))**4.)*d_m_p(s)+4.*((sqrt(s_gp(s+1))/ &
                (1.-s_gp(s+1)))**4.)*d_m_p(s+1)+((sqrt(s_gp(s+2))/ &
                (1.-s_gp(s+2)))**4.)*d_m_p(s+2))

            jam = jam + (s_max/(3.*(s_div-1)))*((((s_gp(s))**3.)/ &
                ((1.-s_gp(s))**5.))*d_j(s)+4.*(((s_gp(s+1))**3.)/ &
                ((1.-s_gp(s+1))**5.))*d_j(s+1)+(((s_gp(s+2))**3.)/ &
                ((1.-s_gp(s+2))**5.))*d_j(s+2))

          end do

          mass = sqrt(kappa)*c*c*mass*4.*pi*((r_e)**3.)/G
          mass_0 = sqrt(kappa)*c*c*mass_0*4.*pi*((r_e)**3.)/G
          mass_p = sqrt(kappa)*c*c*mass_p*4.*pi*((r_e)**3.)/G

          if (r_ratio .eq. 1.) then 

              jam = 0.d0
              omega_f = 0.d0

          else 

              jam = (jam*4.*pi*(r_e**4.))*(C**3.)*kappa/G

              do s = 1, s_div
              
                do m = 1, m_div

                    omega_f(s,m) = (C/sqrt(kappa))*(omega_h(s,m)/r_e)

                end do

              end do

          end if

          do s = 1, s_div

            do m = 1, m_div 

                T(s,m) = 0.5*jam*omega_f(s,m)

                    if (r_ratio .eq. 1.) then 

                        im(s,m) = 0.d0

                    else 

                        im(s,m) = jam/omega_f(s,m)

                    end if 

                W(s,m) = mass_p-mass+T(s,m)

            end do

          end do

          ! redshifts 
          z_p = exp(-0.5*(gamma_pole))-1.
          z_b = sqrt((1.+velocity_equator)/(1.-velocity_equator))* &
              (exp(-0.5*(gamma_equator+rho_equator))/(1.- & 
              omega_equator*r_e*exp(-rho_equator)))-1.
          z_f = sqrt((1.-velocity_equator)/(1.+velocity_equator))* &
              (exp(-0.5*(gamma_equator+rho_equator))/(1.+ &
              omega_equator*r_e*exp(-rho_equator)))-1.

          ! kepler frequency 
          do s = 1, s_div

            d_o_e(s) = deriv_s(omega,s,1)
            d_g_e(s) = deriv_s(gam,s,1)
            d_r_e(s) = deriv_s(rho,s,1)
            d_v_e(s) = deriv_s(velocity,s,1) 

            doe = interp(s_gp,d_o_e,s_div,s_e) 
            dge = interp(s_gp,d_g_e,s_div,s_e)
            dre = interp(s_gp,d_r_e,s_div,s_e)
            dve = interp(s_gp,d_v_e,s_div,s_e)

            vek = (doe/(8.+dge-dre))*r_e*exp(-rho_equator)+ &
                sqrt(((dge+dre)/(8.+dge-dre))+((doe/(8.+dge-dre))* &
                r_e*exp(-rho_equator))**2.)

          end do

          omega_k = (c/sqrt(kappa))*(omega_equator+vek* &
              exp(rho_equator)/r_e)
          
          !USE THIS VALUE AS INITIAL GUESS
         
                   
          
                   
                   
              
              

      end subroutine





