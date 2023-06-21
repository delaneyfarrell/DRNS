!--------------------------------------------------------------------------
! m0_sequence subroutine
!--------------------------------------------------------------------------

          subroutine m0_models

              use vars

              implicit none
              integer :: s, m
              double precision :: dr, dif_m0, d_ratio_m0, sgn

              dr = 0.1
              r_ratio = 1.-dr
              call center(e_center)
              call TOV_guess
              a_check = 0.d0
              call iterate

              if (a_check .eq. 200) then

                  dif_m0 = -1.
                  sgn = -1.

              else

                  call comp_phys
                  dif_m0 = m0-mass_0
                  d_ratio_m0 = dabs(dif_m0)/m0

              end if

              ! if rest mass is greater than desired, reverse direction
              ! and cut stepsize in half

              if (dif_m0 .lt. 0.d0) then

                  dr = -dr/2.d0

              end if

              do while (d_ratio_m0 .gt. m0_error .and. r_ratio .le. 1.)

                if (dif_m0*sgn .lt. 0.d0) then

                    sgn = dif_m0
                    dr = -dr/2.

                end if

                r_ratio = r_ratio-dr
                a_check = 0
                call iterate

                if (a_check .eq. 200) then

                    dif_m0 = -1.

                else

                    call comp_phys

                    do s = 1, s_div
                        
                        do m = 1, m_div

                            if (omega_f(s,m) .gt. omega_k) then

                                dif_m0 = -1.

                            else

                                dif_m0 = m0-mass_0
                                d_ratio_m0 = dabs(dif_m0)/m0

                            end if

                        end do 

                    end do

                end if

              end do

!               call print_str_m0

      end subroutine

