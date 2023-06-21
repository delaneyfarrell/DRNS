!--------------------------------------------------------------------------
! kepler models subroutine
!--------------------------------------------------------------------------

       subroutine kepler_models

           use vars 

           implicit none
           integer :: s, m
           double precision :: dr, sgn
           double precision :: d_omega(0:s_div,0:m_div) 
           double precision :: diff_omega(0:s_div,0:m_div)

           ! first rotating model
           call center(e_center)
           call TOV_guess
           d_omega = 1.0
           sgn = 1.0
           dr = 0.1
           r_ratio = 1.0
           diff_omega = 1.0

           ! routine: decrease r_p. whenever d_omega changes sign or
           ! iteration doesn't converge, cut stepsize in half and
           ! reverse direction.

           do s = 1, s_div

            do m = 1, m_div

           do while ((diff_omega(s,m) .gt. omega_error .or. d_omega(s,m)
     $     .lt. 0.d0) .and. r_ratio .le. 1.0) 

            if (d_omega(s,m)*sgn .lt. 0.0) then

                sgn = d_omega(s,m)
                dr = -dr/2.0

            end if

            r_ratio = r_ratio-dr
            a_check = 0
            call iterate
       
            if (a_check .eq. 200) then 
                
                d_omega(s,m) = -1.0

            else 

                call comp_phys
                d_omega(s,m) = omega_k-omega_f(s,m)
                diff_omega(s,m) = dabs(omega_f(s,m)-omega_k)/omega_k

            end if

            if (r_ratio .gt. 1.0) then

                print*, 'r_ratio >= 1.0'

            end if

          end do

          end do

          end do

          call comp_phys
          call print_str_kepler

      end subroutine
      
!--------------------------------------------------------------------------
! kepler models (single star) subroutine
!--------------------------------------------------------------------------

      subroutine kepler_models_single 

          use vars

          implicit none
          integer :: s,m 
          double precision :: dr, sgn
          double precision :: d_omega(0:s_div,0:m_div)
           double precision :: diff_omega(0:s_div,0:m_div)

          ! first rotating model
          call center(e_center) 
          call TOV_guess
          d_omega = 1.0
          sgn = 1.0
          dr = 0.1
          r_ratio = 1.0
          diff_omega = 1.0

          ! routine: decrease r_p. whenever d_omega changes sign or
          ! iteration doesn't converge, cut stepsize in half and
          ! reverse direction.

          do s = 1, s_div

            do m = 1, m_div

          do while ((diff_omega(s,m) .gt. omega_error .or. d_omega(s,m)
     $    .lt. 0.d0) .and. r_ratio .le. 1.0) 

            if (d_omega(s,m)*sgn .lt. 0.0) then

                sgn = d_omega(s,m)
                dr = -dr/2.0

            end if

            r_ratio = r_ratio-dr
            a_check = 0
            call iterate
         
            if (a_check .eq. 200) then 

                d_omega(s,m) = -1.0

            else

                call comp_phys
                d_omega(s,m) = omega_k-omega_f(s,m)
                diff_omega(s,m) = dabs(omega_f(s,m)-omega_k)/omega_k

            end if

            if (r_ratio .gt. 1.0) then 
                
                print *,'r_ratio >= 1.0.'

            end if

          end do

          end do

          end do 
       
          call comp_phys
          call print_str_kepler_single

      end subroutine
       
!--------------------------------------------------------------------------
! kepler model printing subroutine
!--------------------------------------------------------------------------

       subroutine print_str_kepler

           use vars

           implicit none
           integer :: i, j, ji, s, m
           double precision :: r_dim, e_dim, mass0_print, r_eprint
           double precision :: e_center_print, mass_print
           double precision :: omegak_print, mass_0_print
           double precision :: omegaf_print(0:s_div,0:m_div)
         
           mass0_print = mass_0/m_sun
           r_eprint = r_ec/1.d5
           e_center_print = (e_center/(c*c*k_scale))/(1.7827d12)
           mass_print = mass/m_sun
           omegak_print = omega_k/(2.*pi)

           do s = 1, s_div
            
            do m = 1, m_div

                omegaf_print(s,m) = omega_f(s,m)/(2.*pi)

            end do

           end do
      
       write(*,1010) e_center_print, mass_print, mass0_print, r_eprint,
     $ r_ratio*r_eprint, omegak_print
              
1010  format(1F15.5,xx,1E15.5,xx,1F15.5,xx,1F15.5,xx,1F15.5,
     $xx,1F15.5)
       
      end subroutine

!--------------------------------------------------------------------------
! kepler model printing subroutine - single star
!--------------------------------------------------------------------------

      subroutine print_str_kepler_single

          use vars

          implicit none
          integer :: i, j, ji, s, m
          double precision :: r_dim, e_dim, mass0_print, r_eprint
          double precision :: e_center_print, mass_print
          double precision :: omegak_print, mass_0_print, p_dim
          double precision :: omegaf_print(0:s_div,0:m_div) 

          open(unit=10, file='output_struct.dat')

          mass0_print = mass_0/m_sun
          r_eprint = r_ec/1.d5 
          e_center_print = (e_center/(c*c*k_scale))/(1.7827d12)
          mass_print = mass/m_sun
          omegak_print = omega_k/(2.*pi)

          do s = 1, s_div
            
            do m = 1, m_div

                omegaf_print(s,m) = omega_f(s,m)/(2.*pi)

            end do

          end do
     
          print *, 'Central Density (MeV/fm^3) =', e_center_print
          print *, 'Stellar Mass (Msun) =', mass_print
          print *, 'Baryonic Mass (Msun) =', mass0_print
          print *, 'Equatorial Radius (km) =', r_eprint
          print *, 'Polar Radius (km) =', r_ratio*r_eprint
          !print *, 'Frequency (Hz) =', omegaf_print
          print *, 'Kepler Frequency (Hz) =', omegak_print
       
          ! detailed print output
          if (print_option .eq. 1) then

              do j = 1, m_div, 1

                do i = 1, s_div/2 + 2
                
                    r_dim = (((r_ec*(s_gp(i)/(1.d0 -s_gp(i))))))/1.d5
                    e_dim = ((energy(i,j))/(c*c*k_scale))/(1.7827d12)
                    p_dim = ((pressure(i,j))/(c*c*k_scale))/(1.7827d12)


           write(10,1012) theta(j), r_dim, (p_dim), alpha(i,j), 
     $ gam(i,j), rho(i,j), omega(i,j), velocity_sq(i,j)
           
                end do
                
                write(10,*) " "

              end do

          end if
       
1012  format(1F15.5,xxx,1E15.5,xxx,1F15.5,xxx,1F15.5,xxx,1F15.5,
     $xxx,1F15.5,xxx,1F15.5,xxx,1F15.5)
               
      end subroutine 
         
