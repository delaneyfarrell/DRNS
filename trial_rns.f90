!-------------------------------------------------------------------------
! delaney farrell 9/1/21
! differential rotation code for neutron stars - version 1.5
! needs to be compiled with: iteration.f90, TOV_guess.f90, comp.f90,
! comp_f_p.f90, kepler.F, and spin_up_freq.f90 (see makefile)  
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! global variables 
!--------------------------------------------------------------------------

      module vars

          implicit none 

          integer :: s_div = 329, m_div = 125, r_div = 900
          integer :: l_max = 10 
          integer :: n_nearest, a_check
          integer :: print_dif = 0, n_tab = 0, print_option 

          double precision :: s_max = 0.9999d0, gamma_p, cf
          double precision :: e_center, p_center, ac 
          double precision :: r_final = 0.d0, m_final = 0.d0
          double precision :: pi = dacos(-1.d0) 
          double precision :: r_min = 1.d-15
          double precision :: mb = 1.66d-24
          double precision :: r_is_final = 0.d0
          double precision :: e_surface, p_surface, k_rescale 
          double precision :: c = 2.9979d10
          double precision :: G = 6.6832d-8
          double precision :: kappa = 1.346790806509621d+13
          double precision :: k_scale = 1.112668301525780d-36
          double precision :: s_e = 0.5d0
          double precision :: m_sun = 1.987e33
          !double precision :: ac = 0.01d0
          double precision :: mass, r_e, r_e_guess, r_ec, z_p 
          double precision :: gamma_center_old, rho_center_old
          double precision :: gamma_pole_h, gamma_center_h
          double precision :: gamma_equator_h
          double precision :: rho_pole_h, rho_center_h, rho_equator_h 
          double precision :: r_p, s_p, h_center
          double precision :: omega_h_c, omega_h_equator
          double precision :: Omega_rot_eq, Omega_rot_c
          double precision :: omega_equator_h
          double precision :: enthalpy_min
          double precision :: accu, ds, dm 
          double precision :: r_ratio, velocity_equator 
          double precision :: mass_p, mass_0, jam, n_p
          double precision :: z_b, z_f, omega_k, kappa_pole
          double precision :: omega_error, a, m0_error, m0, j_const

          double precision, allocatable :: s_gp(:), s1_s(:), one_s(:) 
          double precision, allocatable :: mu(:), one_m2(:), theta(:) 
          double precision, allocatable :: sin_theta(:) 
          double precision, allocatable :: f_rho(:,:,:) 
          double precision, allocatable :: f_gamma(:,:,:)
          double precision, allocatable :: f_omega(:,:,:)
          double precision, allocatable :: p_2n(:,:), p1_2n_1(:,:) 
          double precision, allocatable :: r_gp(:), r_is_gp(:), m_gp(:) 
          double precision, allocatable :: lambda_gp(:), e_d_gp(:)
          double precision, allocatable :: nu_gp(:) 
          double precision, allocatable :: gamma_guess(:,:)
          double precision, allocatable :: rho_guess(:,:)
          double precision, allocatable :: omega_guess(:,:) 
          double precision, allocatable :: alpha_guess(:,:) 
          double precision, allocatable :: alpha(:,:), gam(:,:) 
          double precision, allocatable :: rho(:,:), omega(:,:)
          double precision, allocatable :: gamma_mu_0(:), rho_mu_0(:) 
          double precision, allocatable :: omega_mu_0(:)
          double precision, allocatable :: omega_h_mu_0(:) 
          double precision, allocatable :: gamma_mu_1(:), rho_mu_1(:) 
          double precision, allocatable :: velocity_sq(:,:)
          double precision, allocatable :: enthalpy(:,:)
          double precision, allocatable :: energy(:,:) 
          double precision, allocatable :: pressure(:,:) 
          double precision, allocatable :: s_rho(:,:), s_gamma(:,:) 
          double precision, allocatable :: s_omega(:,:) 
          double precision, allocatable :: d1_rho(:,:), d1_gamma(:,:) 
          double precision, allocatable :: d1_omega(:,:) 
          double precision, allocatable :: d2_rho(:,:), d2_gamma(:,:)
          double precision, allocatable :: d2_omega(:,:)
          double precision, allocatable :: dgds(:,:), dgdm(:,:)
          double precision, allocatable :: dadm(:,:)
          double precision, allocatable :: log_e_tab(:), log_p_tab(:) 
          double precision, allocatable :: log_h_tab(:), log_n0_tab(:)
          double precision, allocatable :: omega_h(:,:), omega_f(:,:) 
          double precision, allocatable :: W(:,:), T(:,:), im(:,:) 
          double precision, allocatable :: Omg_h_eq(:)
          double precision              :: Omg_EQU
          
          character :: dummy 

      end module 

!--------------------------------------------------------------------------
! main program 
!--------------------------------------------------------------------------

      program dr_rns

          use vars

          implicit none 
          integer :: ialloc, n_of_models, option, j, I 
          double precision :: fix_error
          double precision :: err_tol
          double precision :: e_min, e_max, m_fix
          double precision :: m0_const, omega_const, jconst 
          double precision :: ratio_max, ratio_min, ratio_loop
          double precision :: dratio 
          
          double precision :: mass0_print, r_eprint, e_center_print
          double precision :: mass_print, Omega_rot_c_print
          double precision :: ra_max, ra_min
          double precision :: Omega_rot_eq_print, interp
          
  
          n_nearest = 1 

          ! allocate space for variables 
          allocate(dadm(0:s_div, 0:m_div), stat=ialloc) 
          allocate(dgds(0:s_div, 0:m_div), stat=ialloc)
          allocate(dgdm(0:s_div, 0:m_div), stat=ialloc)
          allocate(d2_rho(0:s_div, 0:l_max), stat=ialloc)
          allocate(d2_gamma(0:s_div, 0:l_max), stat=ialloc)
          allocate(d2_omega(0:s_div, 0:l_max), stat=ialloc)
          allocate(d1_rho(0:l_max, 0:s_div), stat=ialloc)
          allocate(d1_gamma(0:l_max, 0:s_div), stat=ialloc)
          allocate(d1_omega(0:l_max, 0:s_div), stat=ialloc)
          allocate(s_rho(0:s_div, 0:m_div), stat=ialloc)
          allocate(s_omega(0:s_div, 0:m_div), stat=ialloc)
          allocate(s_gamma(0:s_div, 0:m_div), stat=ialloc)
          allocate(pressure(0:s_div, 0:m_div), stat=ialloc)
          allocate(energy(0:s_div, 0:m_div), stat=ialloc)
          allocate(enthalpy(0:s_div, 0:m_div), stat=ialloc)
          allocate(velocity_sq(0:s_div, 0:m_div), stat=ialloc)
          allocate(omega_guess(0:s_div, 0:m_div), stat=ialloc)
          allocate(alpha_guess(0:s_div, 0:m_div), stat=ialloc)
          allocate(gamma_guess(0:s_div, 0:m_div), stat=ialloc)
          allocate(rho_guess(0:s_div, 0:m_div), stat=ialloc)
          allocate(gam(0:s_div, 0:m_div), stat=ialloc)
          allocate(rho(0:s_div, 0:m_div), stat=ialloc)
          allocate(alpha(0:s_div, 0:m_div), stat=ialloc)
          allocate(omega(0:s_div, 0:m_div), stat=ialloc)
          allocate(rho_mu_0(0:s_div), stat=ialloc)
          allocate(rho_mu_1(0:s_div), stat=ialloc)
          allocate(gamma_mu_0(0:s_div), stat=ialloc)
          allocate(gamma_mu_1(0:s_div), stat=ialloc)
          allocate(omega_mu_0(0:s_div), stat=ialloc)
          allocate(omega_h_mu_0(0:s_div), stat=ialloc) 
          allocate(r_is_gp(0:r_div), stat=ialloc)
          allocate(r_gp(0:r_div), stat=ialloc)
          allocate(m_gp(0:r_div), stat=ialloc)
          allocate(lambda_gp(0:r_div), stat=ialloc)
          allocate(e_d_gp(0:r_div), stat=ialloc)
          allocate(nu_gp(0:r_div), stat=ialloc)
          allocate(s_gp(0:s_div), stat=ialloc)
          allocate(s1_s(0:s_div), stat=ialloc)
          allocate(one_s(0:s_div), stat=ialloc)
          allocate(mu(0:m_div), stat=ialloc)
          allocate(one_m2(0:m_div), stat=ialloc)
          allocate(theta(0:m_div), stat=ialloc)
          allocate(sin_theta(0:m_div), stat=ialloc)
          allocate(f_rho(0:s_div, 0:l_max, 0:s_div), stat=ialloc)
          allocate(f_gamma(0:s_div, 0:l_max, 0:s_div), stat=ialloc)
          allocate(f_omega(0:s_div, 0:l_max, 0:s_div), stat=ialloc)
          allocate(p_2n(0:m_div, 0:l_max), stat=ialloc)
          allocate(p1_2n_1(0:m_div, 0:l_max), stat=ialloc)
          allocate(omega_h(0:s_div, 0:m_div), stat=ialloc)
          allocate(omega_f(0:s_div, 0:m_div), stat=ialloc)
          allocate(W(0:s_div, 0:m_div), stat=ialloc)
          allocate(T(0:s_div, 0:m_div), stat=ialloc)
          allocate(im(0:s_div, 0:m_div), stat=ialloc)
          allocate(Omg_h_eq(0:s_div), stat=ialloc)

          ds = s_max/(s_div-1.d0)
          dm = 1.d0/(m_div-1.d0) 

          ! open parameter file 
          open(unit=21, file='parameters.dat', status='unknown') 

          ! parameters for tabulated EOS
          e_surface = 7.8*c*c*k_scale 
          p_surface = 1.01e8*k_scale 
          enthalpy_min = 1.0/(c*c) 

          ! call functions
          call make_grid
          call load_eos
          call comp_f_p

          ! polytropic stars default 
          cf = 0.5d0                     ! relaxation constant
          accu = 1.d-4                  ! accuracy 
          err_tol = 1.d-4               ! error tolerance 
          n_of_models = 1               ! # of models to compute
          e_min = 2.66e15*c*c*k_scale   ! min energy density 
          e_max = 1.e16*c*c*k_scale     ! max energy density 
          r_ratio = 0.75                ! ratio r_p / r_e 
          m_fix = 1.4*m_sun             ! gravitational mass
          m0_const = 1.4*m_sun          ! baryonic mass
          omega_const = 1.e-5           ! omega
          j_const = 0.                  ! angular momentum 

          ! reading options (dummy = skipping a line) 
          read(21,*) dummy      
          read(21,*) dummy      
          read(21,*) dummy      
          read(21,*) dummy      
          read(21,*) dummy      
          read(21,*) dummy     
          read(21,*) dummy      
          read(21,*) dummy      
          read(21,*) dummy      
          read(21,*) dummy      
          read(21,*) dummy      
          read(21,*) dummy      
          read(21,*) option

          ! OPTION 1: single star)
          if (option .eq. 1) then 

              ! read suboptions 
              read(21,*) dummy
              read(21,*) e_center, r_ratio, ac, print_option 

              ! screen output 
              print*, '---------------------------------------------'
              print*, 'begin program for tabulated EOS'
              print*, 'option: single star'
              print*, 'central density = ', e_center, ' MeV/fm^3'
              print*, 'ratio r_p/r_e = ', r_ratio 
              print *, '1/A =', ac
              print*, '---------------------------------------------'

              ! convert units 
              e_center = e_center*1.7828d12       ! MeV/fm^3 to g/cm^3
              e_center = e_center*c*c*k_scale     ! g/cm^3 to dimensionless

              ac = 1./ac
              ! call routines for single star
              call center(e_center)
              call TOV_guess
              call iterate
              call comp_phys
              call print_str

          ! option 2: Constant ratio (r) - Range of densities
          else if (option .EQ. 2) then

              ! reading suboptions
              read(21,*) dummy      !SKIPPING LINE
              read(21,*) e_min, e_max, r_ratio, ac, n_of_models

              ! screen output
              print*, '---------------------------------------------'
              print *, 'solving a sequence of constant e_c, r_ratio =', &
                  r_ratio
              print *, 'central density range = ',e_min, ' to ', e_max
              print *, 'number of stars =', n_of_models 
              print *, '1/A =', ac 
              print*, '---------------------------------------------'

              ! converting units 
              e_min = e_min*1.7827d12       ! MeV/fm3 to g/cm3
              e_min = e_min*C*C*K_SCALE     ! g/cm3 to dimensionless
              e_max = e_max*1.7827d12       ! MeV/fm3 to g/cm3
              e_max = e_max*C*C*K_SCALE     ! g/cm3 to dimensionless

              ! output file: OPT_2_OUTPUT.DAT
              OPEN(110, FILE = 'OPT_2_OUTPUT.DAT')          
              write(110,*) 'Central Density (MeV/Fm3), Stellar & 
                  Mass(M_sun), Baryonic Mass(M_sun), Circmf. Radius & 
                  (km), Omega_c (hz), Omega_eq (hz)'

              ! calculations use A, not 1/A
              ac = 1./ac
       
              ! density steps
              a = (e_max/e_min)**(1.0/(n_of_models -1.0))  

              ! main routine
              do j= 1, n_of_models

                ! energy density at center
                e_center = (a**(1.0*j -1.0))*e_min 
                
                ! call routines
                call center(e_center)
                call TOV_guess
                call iterate
                call comp_phys
             
                ! print results
                mass0_print = mass_0/m_sun
                r_eprint = r_ec/1.d5
                e_center_print = (e_center/(c*c*k_scale))/(1.7827d12) 
                mass_print = mass/m_sun
                Omega_rot_c_print = Omega_rot_c*(c/sqrt(kappa))
              
                do i = 1, s_div
                    Omg_h_eq(i) = Omega_h(i,1)
                end do
                Omg_EQU = interp(s_gp, Omg_h_eq, s_div, s_e)
               
              
                 Omega_rot_eq_print = Omg_EQU*(c/sqrt(kappa))

                 write(110,1010) e_center_print, mass_print, & 
                     mass0_print, r_eprint, Omega_rot_c_print, & 
                     Omega_rot_eq_print

              end do              
              close(110) 

      ! option 3: Constant Baryonic Mass - Range of densities
      else if (option .EQ. 3) then 

          ! reading suboptions
          read(21,*) dummy      !SKIPPING LINE
          read(21,*) M0, e_min, e_max, ac, n_of_models
          
          ! frequency tolerance
          M0_error = fix_error
          
          ! converting mass
          M0 = M0*M_sun
           
          ! screen output
          print*, '---------------------------------------------'
          print *, 'solving sequenc of constant baryon mass stars' 
          print *, 'central density range = ',e_min, ' to ', e_max
          print *, 'number of stars = ', n_of_models 
          print *, '1/A =', ac 
          print*, '---------------------------------------------'
            
          ! converting units 
          e_min = e_min*1.7827d12  !MeV/fm3 to g/cm3
          e_min = e_min*C*C*K_SCALE !g/cm3 to dimensionless
          e_max = e_max*1.7827d12  !MeV/fm3 to g/cm3
          e_max = e_max*C*C*K_SCALE !g/cm3 to dimensionless

          ! output file: OPT_3_OUTPUT.DAT
          OPEN(110, FILE = 'OPT_3_OUTPUT.DAT')          
          write(110,*) 'Central Density (MeV/Fm3), Stellar Mass(M_sun), &
              Baryonic Mass(M_sun), Circmf. Radius (km), Omega_c (hz),&
              Omega_eq (hz)'          
          
          ! inverting A
          ac = 1./ac

          ! density steps
          a = (e_max/e_min)**(1.0/(n_of_models -1.0))
          
          ! main routine
          do j = 1, n_of_models

            e_center =   (a**(1.0*j -1.0))*e_min
          
            ! call routine
            call m0_models

            ! print results
            mass0_print = mass_0/m_sun
            r_eprint = r_ec/1.d5
            e_center_print = (e_center/(c*c*k_scale))/(1.7827d12) 
            mass_print = mass/m_sun
            Omega_rot_c_print = Omega_rot_c*(c/sqrt(kappa))
            Omega_rot_eq_print = Omega_rot_eq*(c/sqrt(kappa))
      
            write(110,1010) e_center_print, mass_print, mass0_print, & 
                r_eprint, Omega_rot_c_print, Omega_rot_eq_print
            
          end do
          close(110)

      ! option 4: Constant central dems (r) - Range of ratios 
      else if (option .EQ. 4) then

          ! reading suboptions
          read(21,*) dummy      !SKIPPING LINE
          read(21,*) ra_min, ra_max, e_center, ac, n_of_models

          ! screen output
          print*, '---------------------------------------------'
          print *, 'solving sequence of constant central den.'
          print *, 'e_c = ', e_center
          print *, 'r_ratio range = ',ra_min, ' to ', ra_max
          print *, 'number of stars =', n_of_models 
          print *, '1/A =', ac 

          ! inverting A
          ac = 1./ac

          ! convert units 
          e_center = e_center*1.7828d12       ! MeV/fm^3 to g/cm^3
          e_center = e_center*c*c*k_scale     ! g/cm^3 to dimensionless

          ! output file: OPT_4_OUTPUT.DAT
          OPEN(110, FILE = 'OPT_4_OUTPUT.DAT')          
          write(110,*) '#Central Density (MeV/Fm3), Stellar Mass(M_sun),&
              Baryonic Mass(M_sun), Circmf. Radius (km), Omega_c (hz),&
              Omega_eq (hz)'
       
          ! r_ratio steps
          a = (ra_max/ra_min)**(1.0/(n_of_models -1.0))  
          
          ! main program
          do j = 1, n_of_models
          
             ! energy density at center
             r_ratio = (a**(1.0*j -1.0))*ra_min 
             
             ! call routines
             call center(e_center)
             call TOV_guess
             call iterate
             call comp_phys
             
             ! printing results
             mass0_print = mass_0/m_sun
             r_eprint = r_ec/1.d5
             e_center_print = (e_center/(c*c*k_scale))/(1.7827d12) 
             mass_print = mass/m_sun
             Omega_rot_c_print = Omega_rot_c*(c/sqrt(kappa))
              
             do i = 1, s_div
                Omg_h_eq(i) = Omega_h(i,1)
             end do
             Omg_EQU = interp(s_gp, Omg_h_eq, s_div, s_e)  
              
              Omega_rot_eq_print = Omg_EQU*(c/sqrt(kappa))
      
              write(110,1011) e_center_print, mass_print, mass0_print, &
                  r_eprint, Omega_rot_c_print, Omega_rot_eq_print, r_ratio
             
          end do
          close(110)

1010  format(1F15.5,xx,1F15.5,xx,1F15.5,xx,1F15.5,xx,1F15.5,&
      xx,1F15.5)
1011  format(1F15.5,xx,1F15.5,xx,1F15.5,xx,1F15.5,xx,1F15.5,&
      xx,1F15.5,xx,1F15.5)  

      end if 

      end program

!--------------------------------------------------------------------------
! grid subroutine
!--------------------------------------------------------------------------

      subroutine make_grid

          use vars

          implicit none
          integer :: i

          do i = 1, s_div
            s_gp(i) = s_max*(1.*real(i)-1.d0)/(s_div-1.d0)
            s1_s(i) = s_gp(i)*(1.d0-s_gp(i))
            one_s(i) = 1.d0-s_gp(i)
          end do

          do i = 1, m_div
            mu(i) = (1.*real(i)-1.d0)/(m_div-1.d0) 
            one_m2(i) = 1.d0-mu(i)**2.d0
            theta(i) = acos(mu(i)) 
            sin_theta(i) = sqrt(one_m2(i)) 
          end do

      end subroutine

!--------------------------------------------------------------------------
! center subroutine
!--------------------------------------------------------------------------

      subroutine center(e_c) 

          use vars 

          implicit none
          double precision :: p_at_e, h_at_p, e_c 

          n_nearest = n_tab/2
          p_center = p_at_e(e_c) 
          h_center = h_at_p(p_center)

      end subroutine 

!--------------------------------------------------------------------------
! functions for center (p_at_e, h_at_p) 
!--------------------------------------------------------------------------

      double precision function p_at_e(ec)

          use vars

          implicit none
          double precision :: ec, argtemp, interp

          argtemp = interp(log_e_tab, log_p_tab, n_tab, log10(ec))
          p_at_e = 10.d0**(argtemp) 

      end function 

      double precision function h_at_p(pc) 

          use vars 

          implicit none
          double precision :: pc, argtemp, interp

          argtemp = interp(log_p_tab, log_h_tab, n_tab, log10(pc))
          h_at_p = 10.d0**(argtemp) 

      end function 

!--------------------------------------------------------------------------
! function rtsec_c
!--------------------------------------------------------------------------

      double precision function rtsec_c(x1, x2, xacc, e_c) 

          implicit none
          double precision :: x1, x2, xacc, e_c
          double precision :: fl, f, dx, swap, xl, rts
          double precision :: e_of_rho0 
          integer :: j, max_iter 

          max_iter = 100
          fl = e_of_rho0(x1)-e_c
          f = e_of_rho0(x2)-e_c 

          if (dabs(fl) .lt. dabs(f)) then 
              rts = x1
              xl = x2
              swap = fl
              fl = f
              f = swap 
          else 
              xl = x1
              rts = x2
          end if 

          do j = 1, max_iter

            dx = (xl-rts)*f/(f-fl) 
            xl = rts
            fl = f
            rts = rts+dx
            f = e_of_rho0(rts)-e_c

            if (dabs(dx) .lt. xacc .or. f .eq. 0.d0) then 
                rtsec_c = rts
                go to 10 
            end if 

          end do 

          print*, 'max # of iterations exceeded in rtsec_c'
          rtsec_c = 0.d0

10        continue 

      end function

!--------------------------------------------------------------------------
! function e_of_rho0 
!--------------------------------------------------------------------------

      double precision function e_of_rho0(x) 

          use vars 

          implicit none
          double precision :: x

          e_of_rho0 = (x**gamma_p)/(gamma_p-1.d0)+x

      end function

!--------------------------------------------------------------------------
! function e_at_p
!--------------------------------------------------------------------------

      double precision function e_at_p(pc) 

          use vars

          implicit none
          double precision :: pc, argtemp, interp

          argtemp = interp(log_p_tab, log_e_tab, n_tab, log10(pc))
          e_at_p = 10.d0**argtemp 

      end function

!--------------------------------------------------------------------------
! function n0_at_e
!--------------------------------------------------------------------------

      double precision function n0_at_e(ec)

          use vars

          implicit none
          double precision :: ec, argtemp, interp

          argtemp = interp(log_e_tab, log_n0_tab, n_tab,log10(ec))
          n0_at_e = 10.d0**argtemp

      end function

!--------------------------------------------------------------------------
! function p_at_h
!--------------------------------------------------------------------------

      double precision function p_at_h(h_c)

          use vars

          implicit none
          double precision :: h_c, argtemp, interp

          argtemp = interp(log_h_tab, log_p_tab, n_tab, log10(h_c))
          p_at_h = 10.d0**argtemp

      end function

!--------------------------------------------------------------------------
! function for interpolation (interp) 
!--------------------------------------------------------------------------

      double precision function interp(xp, yp, np, xb) 

          use vars

          implicit none
          integer :: k, m, np
          double precision :: y, xp(0:np), yp(0:np), xb

          m = 4
          call hunt(xp, np, xb, n_nearest) 

          k = min(max(n_nearest-(m-1)/2, 1), np+1-m) 

          if (xb .eq. xp(k) .or. xb .eq. xp(k+1) .or. xb .eq. xp(k+2) & 
              .or. xb .eq. xp(k+3)) then 
              xb = xb+1.d-12 
          end if

          y = (xb-xp(k+1))*(xb-xp(k+2))*(xb-xp(k+3))*yp(k)/ &
              ((xp(k)-xp(k+1))*(xp(k)-xp(k+2))*(xp(k)-xp(k+3))) & 
              +(xb-xp(k))*(xb-xp(k+2))*(xb-xp(k+3))*yp(k+1)/ &
              ((xp(k+1)-xp(k))*(xp(k+1)-xp(k+2))*(xp(k+1)-xp(k+3))) &
              +(xb-xp(k))*(xb-xp(k+1))*(xb-xp(k+3))*yp(k+2)/ & 
              ((xp(k+2)-xp(k))*(xp(k+2)-xp(k+1))*(xp(k+2)-xp(k+3))) &
              +(xb-xp(k))*(xb-xp(k+1))*(xb-xp(k+2))*yp(k+3)/ & 
              ((xp(k+3)-xp(k))*(xp(k+3)-xp(k+1))*(xp(k+3)-xp(k+2)))

          interp = y

      end function

!--------------------------------------------------------------------------
! hunt (nearest gridpoint) subroutine
!--------------------------------------------------------------------------

      subroutine hunt(xx, n, x, j_lo) 

          use vars

          implicit none 
          integer :: n, j_lo, j_hi, j_m, iter
          double precision :: x, xx(0:n) 
          logical :: ascnd 

          ascnd = xx(n) .ge. xx(1)

          if (j_lo .le. 0 .or. j_lo .gt. n) then 
              j_lo = 0
              j_hi = n+1
              go to 3
          end if

          iter = 1
          if (x .ge. xx(j_lo) .eqv. ascnd) then 

              1 continue 
              j_hi = j_lo+iter
              
              if (j_hi .gt. n) then 
                  j_hi = n+1
              else if (x .ge. xx(j_hi) .eqv. ascnd) then 
                  j_lo = j_hi
                  iter = iter + iter
                  go to 1
              end if 

          else 

              j_hi = j_lo
              2 continue 
              j_lo = j_hi-iter

              if (j_lo .lt. 1) then 
                  j_lo = 0
              elseif (x .lt. xx(j_lo) .eqv. ascnd) then 
                  j_hi = j_lo
                  iter = iter + iter 
                  go to 2
              end if 

          end if 

          3 continue
          if (j_hi-j_lo .eq. 1) then 

              if (x .eq. xx(n)) then
                  j_lo = n-1
              end if

              if (x .eq. xx(1)) then
                  j_lo = 1
              end if
              return 

          end if

          j_m = (j_hi+j_lo)/2

          if (x .ge. xx(j_m) .eqv. ascnd) then 
              j_lo = j_m
          else 
              j_hi = j_m 
          end if 
          goto 3

      end subroutine 
      
!--------------------------------------------------------------------------
! print subroutine for spherically sym.
!--------------------------------------------------------------------------

      subroutine print_str_sph 

          use vars

          implicit none 
          integer :: i, j, ji, m, s, ialloc 
          double precision :: r_dim, e_dim, mass0_print 
          double precision :: r_eprint, e_center_print, mass_print
          double precision :: omegak_print 
          double precision, allocatable :: omegaf_print(:,:)

          allocate(omegaf_print(0:s_div, 0:m_div), stat=ialloc)

          mass0_print = mass_0/m_sun
          r_eprint = r_ec/1.d5
          e_center_print = (e_center/(c*c*k_scale))/(1.7828d12)
          mass_print = mass/m_sun
          omegak_print = omega_k/(2.*pi) 

          do s = 1, s_div
            do m = 1, m_div
                omegaf_print(s,m) = (c/sqrt(kappa))*omega_f(s,m)/(2.*pi)
            end do
          end do

          write(*,1010) e_center_print, mass_print, mass0_print, &
              r_eprint, r_ratio*r_eprint, omegak_print 

1010      format(1F15.5, xx, 1F15.5, xx, 1F15.5, xx, 1F15.5, xx, & 
              1F15.5, xx, 1F15.5)

      end subroutine

!--------------------------------------------------------------------------
! load EOS subroutine 
!--------------------------------------------------------------------------

      subroutine load_eos 

          use vars

          implicit none
          integer :: i, ialloc
          double precision :: pr, rhor, hr, n0r

          ! open file
          open(unit=91, file='eos.dat', status='unknown') 

          ! read total number of points
          read(91,*) n_tab

          ! allocate variables
          allocate(log_e_tab(1:n_tab+1), stat=ialloc)
          allocate(log_p_tab(1:n_tab+1), stat=ialloc)
          allocate(log_h_tab(1:n_tab+1), stat=ialloc)
          allocate(log_n0_tab(1:n_tab+1), stat=ialloc)

          ! read EOS - CGS system
          do i = 1, n_tab
            read(91,*) rhor, pr, hr, n0r
            log_e_tab(i) = log10(rhor*c*c*k_scale)    ! dimensionless
            log_p_tab(i) = log10(pr*k_scale)          ! dimensionless
            log_h_tab(i) = log10(hr/(c*c))            ! dimensionless
            log_n0_tab(i) = log10(n0r)                ! still w/ dims
          end do 

      end subroutine

!--------------------------------------------------------------------------
! print subroutine
!--------------------------------------------------------------------------

      subroutine print_str

          use vars

          implicit none
          integer :: i, j, ji, m, s, ialloc 
          double precision :: r_dim, e_dim, mass0_print, r_eprint
          double precision :: e_center_print, mass_print
          double precision :: omegak_print, mass_0_print, p_dim
          double precision :: omega_h_dim
          double precision, allocatable :: omegaf_print(:,:)
          DOUBLE PRECISION :: Xg_dim(1000,1000)
          DOUBLE PRECISION :: Yg_dim(1000,1000)
          DOUBLE PRECISION :: eprint(1000,1000)
          DOUBLE PRECISION :: RadH(1000),RadSURF(1000)
          DOUBLE PRECISION :: omegaprint(1000,1000)
          DOUBLE PRECISION :: RMAT(1000,1000)
          DOUBLE PRECISION :: interp, r_temp,entRa(1000)
          DOUBLE PRECISION :: pprint(1000,1000)
          

          open(unit = 11, file = 'output_struct.dat') 

          allocate(omegaf_print(0:s_div, 0:m_div), stat=ialloc)

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

          print*, 'central density (MeV/fm^3) = ', e_center_print
          print*, 'stellar mass (m_sun) = ', mass_print
          print*, 'baryonic mass (m_sun) = ', mass0_print
          print*, 'equatorial radius (km) = ', r_eprint
          print*, 'polar radius (km) = ', r_ratio*r_eprint
          !print*, 'frequency (hz) = ', omegaf_print
          print*, 'kepler frequency (km) = ', omegak_print

          if (print_option .eq. 1) then 
              do j = 1, m_div
                do i = 1, s_div
                    r_dim = (((sqrt(kappa)*r_e*(s_gp(i)/(1.d0-s_gp(i))))))/1.d5
                    e_dim = ((energy(i,j))/(c*c*k_scale))/(1.7827d12)
                    eprint(i,j) = e_dim
                    p_dim = ((pressure(i,j))/(c*c*k_scale))/(1.7827d12)
                    
                    pprint(i,j) = enthalpy(i,j)
                    omega_h_dim = omega_h(i,j)*(c/sqrt(kappa))
                    omegaprint(i,j) = omega_h_dim

                    ! coordinates
                    Xg_dim(i,j) = r_dim*sin(theta(j))
                    Yg_dim(i,j) = r_dim*cos(theta(j))

                    write(11,1010) theta(j), r_dim, p_dim, alpha(i,j), &
                        gam(i,j), rho(i,j), omega(i,j), velocity_sq(i,j)&
                        , omega_h_dim
                end do
                write(11,*) ' '
              end do
          end if

          ! routines for printing stellar maps!
          open(23, FILE = 'X.DAT')
          open(24, FILE = 'Y.DAT')
          open(28, FILE = 'ED.DAT')
          open(29, FILE = 'Freq.DAT')
          open(30, FILE = 'vsq.DAT')
          open(31, FILE = 'Pdim.DAT')
          
          do s = 1, S_DIV/2 +20!+50  
            ! x coordinate
            write(23,*) Xg_dim(s,1:M_DIV), Xg_dim(s,M_DIV:1:-1),&
                -Xg_dim(s,1:M_DIV), -Xg_dim(s,M_DIV:1:-1)

            ! y coordinate
            write(24,*) Yg_dim(s,1:M_DIV), -Yg_dim(s,M_DIV:1:-1),&
                -Yg_dim(s,1:M_DIV), Yg_dim(s,M_DIV:1:-1)

            ! energy density 
            write(28,*) eprint(s,1:M_DIV),eprint(s,M_DIV:1:-1), &
                eprint(s,1:M_DIV),eprint(s,M_DIV:1:-1)

            ! frequency
            write(29,*) omegaprint(s,1:M_DIV),omegaprint(s,M_DIV:1:-1)&
                ,omegaprint(s,1:M_DIV),omegaprint(s,M_DIV:1:-1)

            ! velocity squared
             WRITE(30,*) velocity_sq(s,1:M_DIV), & 
                 velocity_sq(s,M_DIV:1:-1),velocity_sq(s,1:M_DIV), & 
                 velocity_sq(s,M_DIV:1:-1)

             ! pdim
             write(31,*) Pprint(s,1:M_DIV),Pprint(s,M_DIV:1:-1), & 
                 Pprint(s,1:M_DIV),Pprint(s,M_DIV:1:-1)
          end do
          
1010      format(1F15.8,xxx,1E15.8,xxx,1F15.8,xxx,1F15.8,xxx,1F15.8, &
              xxx,1F15.8,xxx,1F15.8,xxx,1F15.8,xxx,1F15.8)

      end subroutine

!--------------------------------------------------------------------------
! print for constant ec
!--------------------------------------------------------------------------

      subroutine print_str_const_ec

          use vars

          implicit none
          integer :: i, j, ji, m, s, ialloc
          double precision :: r_dim, e_dim, mass0_print, r_eprint
          double precision :: e_center_print, mass_print
          double precision :: omegak_print, mass_0_print, p_dim
          double precision, allocatable :: omegaf_print(:,:)

          open(unit=10, file='output_struct.dat') 

          allocate(omegaf_print(0:s_div, 0:m_div), stat=ialloc)

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

          print*, e_center_print, mass_print, mass0_print, r_eprint, &
              r_ratio*r_eprint, omegak_print 

      end subroutine

!--------------------------------------------------------------------------
! print for ratio
!--------------------------------------------------------------------------

      subroutine print_str_rat

          use vars

          implicit none
          integer :: i, j, ji, m, s, ialloc
          double precision :: r_dim, e_dim, mass0_print, r_eprint
          double precision :: e_center_print, mass_print
          double precision :: omegak_print, mass_0_print, p_dim
          double precision, allocatable :: omegaf_print(:,:)
       

          allocate(omegaf_print(0:s_div, 0:m_div), stat=ialloc)

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

          write(*,1013) e_center_print, mass_print, r_eprint, & 
              Im/(1.d45), r_p/r_e

1013      format(1F25.15,xx,1F25.15,xx,1F25.15,xx,1F25.15,xx,1F25.15) 

      end subroutine

