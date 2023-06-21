!--------------------------------------------------------------------------
! spin up subroutine
!--------------------------------------------------------------------------

      subroutine spin_up(n_of_models)

          use vars 

          implicit none 

          integer :: n_ok, n_sp, i, n_acc, ialloc, n_of_models
          integer :: n_of_models2, j, k, s, m 

          double precision :: ec_ok(100), m0_ok(100), ec_sph(100)
          double precision :: m0_sph(100), y2_ok(100), y2_sph(100)
          double precision :: y1_ok, yN_ok, y1_sph, yN_sph, m0_d
          double precision :: eck_int, ecsp_int, y1_facc, yN_facc
          double precision :: y1_tacc, yN_tacc, time_int, freq_int
          double precision :: m0_const_seq, ec_const_seq
          double precision :: freq_k_const_seq, y1_oe, yN_oe, y1_oer
          double precision :: yN_oer, m0_const_seq_vec(30)
          double precision :: ec_const_seq_vec(30), y2_oe(30)
          double precision :: freq_const_seq_vec(30), y2_oer(30)
          double precision :: freq_k_const_seq_vec(30), ec_fin_int
          double precision :: r_eprint_const_seq_vec(30)
          double precision :: r_eprint_const_seq, r_fin_int
          double precision :: m0_max, m0_min, del_m, e_min, e_max
          double precision :: freq_const_seq(0:s_div,0:m_div)

          double precision, allocatable :: time(:), m_acc(:)
          double precision, allocatable :: freq_acc(:), y2_facc(:)
          double precision, allocatable :: y2_tacc(:), m0_vec(:)
          double precision, allocatable :: freq_vec(:), time_vec(:)

          allocate(m0_vec(1:n_of_models), stat=ialloc)
          allocate(time_vec(1:n_of_models), stat=ialloc)
          allocate(freq_vec(1:n_of_models), stat=ialloc) 

          n_of_models2 = 10

          ! open omega_k sequence file
          open(unit=28, file='omega_k.dat') 
          read(28,*) n_ok

          do i = 1, n_ok

            read(28,*) ec_ok(i), m0_ok(i) 

          end do

          ! install and call interpolation routine
          y1_ok = 1.d30
          yN_ok = 1.d30
          call spline(m0_ok,ec_ok,n_ok,y1_ok,yN_ok,y2_ok) 

          ! open spherical sequence file
          open(unit=29, file='spherical.dat') 
          read(29,*) n_sp

          do i = 1, n_sp

            read(29,*) ec_sph(i), m0_sph(i) 

          end do

          ! install and call interpolation routine
          y1_sph = 1.d30
          yN_sph = 1.d30
          call spline(m0_sph,ec_sph,n_sp,y1_sph,yN_sph,y2_sph)

          ! open re-sampled accreting sequence file
          open(unit=30, file='accr_out_sampled.dat') 
          read(30,*) n_acc

          allocate(m_acc(1:n_acc), stat=ialloc)
          allocate(time(1:n_acc), stat=ialloc)
          allocate(freq_acc(1:n_acc), stat=ialloc)
          allocate(y2_facc(1:n_acc), stat=ialloc)
          allocate(y2_tacc(1:n_acc), stat=ialloc) 

          do i = 1, n_acc
            
            read(30,*) time(i), freq_acc(i), m_acc(i) 

          end do

          ! install and call interpolation routine
          y1_facc = 1.d30
          yN_facc = 1.d30
          y1_tacc = 1.d30
          yN_tacc = 1.d30
          call spline(m_acc,freq_acc,n_acc,y1_facc,yN_facc,y2_facc)
          call spline(m_acc,time,n_acc,y1_tacc,yN_tacc,y2_tacc)

          ! define range of mass
          m0_max = m_acc(n_acc)
          m0_min = m_acc(1) 
          del_m = (m0_max-m0_min)/(n_of_models-1)
          m0_vec(1) = m0_min 

          do i = 2, n_of_models

            m0_vec(i) = m0_vec(i-1)+del_m

          end do
          
          ! starting sequence loop 
          do i = 1, n_of_models

            ! find frequency & time of stellar mass
            call splint(m_acc,freq_acc,y2_facc,n_acc,m0_vec(i), & 
                freq_vec(i))
            call splint(m_acc,time,y2_tacc,n_acc,m0_vec(i), & 
                time_vec(i)) 

          end do

          ! screen output 
          print*, '---------------------------------------------'
          print*, 'begin program for tabulated EOS'
          print*, 'option: spin up sequence'
          print*, 'baryon mass range = ', m0_vec(1), 'to ', &
              m0_vec(n_of_models) 
          print*, 'number of stars = ', n_of_models
          print*, '---------------------------------------------'
          print*, 'm0/m_sun, freq. (hz), central density (MeV/fm^3)'

          do i = 1, n_of_models

            ! find kepler & spherical central densities
            call splint(m0_ok,ec_ok,y2_ok,n_ok,m0_vec(i),eck_int)
            call splint(m0_sph,ec_sph,y2_ok,n_sp,m0_vec(i),ecsp_int)

            ! converting mass 
            m0 = m0_vec(i)*m_sun

            ! convert units
            e_min = (eck_int+10.)*1.7827d12   ! MeV/fm^3 to g/cm^3
            e_min = e_min*c*c*k_scale         ! g/cm^3 to dimensionless
            e_max = (ecsp_int-10.)*1.7827d12  ! MeV/fm^3 to g/cm^3
            e_max = e_max*c*c*k_scale         ! g/cm^3 to dimensionless

            ! density increments 
            a = (e_max/e_min)**(1./n_of_models2-1.)

            k = n_of_models2

            ! main routine 
            do j = 1, n_of_models2-1

                e_center = (a**(1.*j-1.))*e_min

                call m0_models_spin_up(m0_const_seq,ec_const_seq, &
                    freq_const_seq,freq_k_const_seq,r_eprint_const_seq) 

                m0_const_seq_vec(k) = m0_const_seq
                ec_const_seq_vec(k) = ec_const_seq
                freq_const_seq_vec(k) = freq_const_seq
                freq_k_const_seq_vec(k) = freq_k_const_seq
                r_eprint_const_seq_vec(k) = r_eprint_const_seq

                k = k-1

            end do

            ! solving spherical sequence 
            e_center = e_max
            r_ratio = 1.e0

            call center(e_center)     ! defining variables at center
            call TOV_guess            ! solve for sph. sym. star
            call iterate              ! main iteration
            call comp_phys            ! physical computations

            m0_const_seq = mass_0/m_sun
            ec_const_seq = (e_center/(c*c*k_scale))/(1.7827d12)
            r_eprint_const_seq = r_ec/1.d5
            freq_k_const_seq = omega_k/(2.*pi)
            m0_const_seq_vec(1) = m0_const_seq
            ec_const_seq_vec(1) = ec_const_seq
            freq_const_seq_vec(1) = freq_const_seq
            freq_k_const_seq_vec(1) = freq_k_const_seq
            r_eprint_const_seq_vec(1) = r_eprint_const_seq


            do s = 1, s_div

                do m = 1, m_div

                    freq_const_seq(s,m) = omega_f(s,m)/(2.*pi)

                end do

            end do

            ! install interpolation routine
            y1_oe = 1.d30
            yN_oe = 1.d30
            y1_oer = 1.d30
            yN_oer = 1.d30

            call spline(freq_const_seq_vec,ec_const_seq_vec, &
                n_of_models2,y1_oe,yN_oe,y2_oe) 

            call spline(freq_const_seq_vec,r_eprint_const_seq_vec, &
                n_of_models2,y1_oer,yN_oer,y2_oer) 

            if (freq_vec(i) .le. freq_const_seq_vec(n_of_models2)) then

                ! call interpolation routine
                call splint(freq_const_seq_vec,ec_const_seq_vec,y2_oe, &
                    n_of_models2,freq_vec(i),ec_fin_int) 

                call splint(freq_const_seq_vec,r_eprint_const_seq_vec, &
                    y2_oer,n_of_models2,freq_vec(i),r_fin_int) 

                print*, m0_vec(i), freq_vec(i), ec_fin_int, r_fin_int

            else 

                print*, m0_vec(i), 'above kepler', ec_fin_int

            end if

          end do

       end subroutine

!--------------------------------------------------------------------------
! spline subroutine
!--------------------------------------------------------------------------

      subroutine spline(x,y,n,yp1,ypN,y2)

           implicit none
           integer :: n, i, k, n_max
           double precision :: x, y, yp1, ypN, y2, u, sig, p, qn, un
           parameter (n_max=11000) 
           dimension :: x(n), y(n), y2(n), u(n_max)

           if (yp1 .gt. .99e30) then

               y2(1) = 0
               u(1) = 0

           else 

               y2(1) = -0.5
               u(1) = (3./(x(2)-x(2)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
               
           end if

           do 11, i = 2, n-1

            sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
            p = sig*y2(i-1)+2.
            y2(i) = (sig-1.)/p
            u(i) = (6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/ &
                (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p

11         continue

            if (ypN .gt. .99e30) then

                qn = 0.
                un = 0.

            else 

                qn = 0.5
                un = (3./(x(n)-x(n-1)))*(ypN-(y(n)-y(n-1))/ & 
                    (x(n)-x(n-1)))

            end if

            y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.)

            do 12, k = n-1, 1, -1

                y2(k) = y2(k)*y2(k+1)+u(k)

12          continue

            return

      end subroutine

!--------------------------------------------------------------------------
! spline interpolation subroutine
!--------------------------------------------------------------------------

      subroutine splint(xa, ya, y2a, n, x, y)

          implicit none
          integer :: k_lo, k_hi, k, n
          double precision :: xa, ya, y2a, h, a, b, x, y
          dimension :: xa(n), ya(n), y2a(n) 

          k_lo = 1
          k_hi = n

1         if (k_hi-k_lo .gt. 2) then

              k = (k_hi+k_lo)/2

              if (xa(k) .gt. x) then 

                  k_hi = k

              else 

                  k_lo = k

              end if 

              goto 1

          end if

          h = xa(k_hi)-xa(k_lo) 

          if (h .eq. 0) print*, 'bad xa input'

          a = (xa(k_hi)-x)/h
          b = (x-xa(k_lo))/h
          y = a*ya(k_lo)+b*ya(k_hi)+((a**3-a)*y2a(k_lo)+ &
              (b**3-b)*y2a(k_hi))*(h**2)/6.

          return 

      end subroutine

!--------------------------------------------------------------------------
! m0 models spin up
!--------------------------------------------------------------------------

      subroutine m0_models_spin_up(m0_const_seq, ec_const_seq, & 
              freq_const_seq, freq_k_const_seq, r_eprint_const_seq) 

          use vars

          implicit none
          integer :: s, m
          double precision :: dr, diff_m0, d_ratio_m0, sgn
          double precision :: m0_const_seq, ec_const_seq
          double precision :: freq_k_const_seq, r_eprint_const_seq
          double precision :: freq_const_seq(0:s_div,0:m_div)

          dr = 0.1
          r_ratio = 1.-dr

          ! first rotating model
          call center(e_center) 
          call TOV_guess
          a_check = 0.d0
          call iterate 

          if (a_check .eq. 200) then 

              diff_m0 = -1.
              sgn = -1.

          else 

              call comp_phys
              diff_m0 = m0-mass_0
              sgn = diff_m0
              d_ratio_m0 = dabs(diff_m0)/m0

          end if

          ! if rest mass is greater than desired, reverse direction &
          ! cut stepsize in half

          if (diff_m0 .lt. 0.d0) then 

              dr = -dr/2.d0

          end if

          do while(d_ratio_m0 .gt. m0_error .and. r_ratio .le. 1.) 

            if (diff_m0*sgn .lt. 0.) then 

                sgn = diff_m0
                dr = -dr/2.

            end if

            r_ratio = r_ratio-dr
            a_check = 0
            call iterate

            if (a_check .eq. 200) then 

                diff_m0 = -1.

            else 

                call comp_phys

                do s = 1, s_div

                    do m = 1, m_div

                        if (omega_f(s,m) .gt. omega_k) then 

                            diff_m0 = -1.

                        else

                            diff_m0 = m0-mass_0
                            d_ratio_m0 = dabs(diff_m0)/m0

                        end if

                    end do

                end do

            end if

          end do

          ! storing relevant variables
          m0_const_seq = mass_0/m_sun
          ec_const_seq = (e_center/(c*c*k_scale))/(1.7827d12)
          freq_k_const_seq = omega_k/(2.*pi)
          r_eprint_const_seq = r_ec/1.d5

          do s = 1, s_div

            do m = 1, m_div

                freq_const_seq(s,m) = omega_f(s,m)/(2.*pi)

            end do

          end do

      end subroutine

