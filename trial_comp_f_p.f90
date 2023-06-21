!--------------------------------------------------------------------------
! legendre polynomial subroutine
!--------------------------------------------------------------------------

      subroutine comp_f_p

          use vars

          implicit none 
          integer :: n, k, j, i, stemp
          double precision :: sk, sk1, sj, sj1, legendre, plgndr

          if (s_max .eq. 1) then

              stemp = s_div-1

          else 

              stemp = s_div

          end if

          j = 1
          n = 0

          do k = 2, s_div

            sk = s_gp(k)
            f_rho(j,n,k) = 1./(sk*(1.-sk))

          end do

          n = 1

          do k = 2, s_div

            sk = s_gp(k)
            sk1 = 1.-s_gp(k)
            f_rho(j,n,k) = 0.d0
            f_gamma(j,n,k) = 1./(sk*sk1)
            f_omega(j,n,k) = 1./(sk*sk1) 

          end do

          do n = 2, l_max

            do k = 1, s_div

                f_rho(j,n,k) = 0.d0
                f_gamma(j,n,k) = 0.d0
                f_omega(j,n,k) = 0.d0 

            end do  

          end do

          k = 1
          n = 0

          do j = 1, s_div

            f_rho(j,n,k) = 0.d0 

          end do

          do j = 1, s_div

            do n = 1, l_max

                f_rho(j,n,k) = 0.d0
                f_gamma(j,n,k) = 0.d0
                f_omega(j,n,k) = 0.d0

            end do

          end do

          n = 0

          do j = 2, s_div

            do k = 2, s_div

                if (s_max .eq. 1. .and. (k .eq. s_div .or. j .eq. &
                    s_div)) then 

                    f_rho(j,n,k) = 0.d0
                    print*, f_rho(j,n,k) 

                else 

                    sk = s_gp(k) 
                    sj = s_gp(j) 
                    sk1 = 1.-sk
                    sj1 = 1.-sj

                    if (k .lt. j) then 

                        f_rho(j,n,k) = ((sj1/sj)**(2.d0*n+1.))* &
                            (sk**(2.*n))/((sk1**(2.*n+2)))

                    else 

                        f_rho(j,n,k) = ((sj/sj1)**(2.d0*n))* &
                            ((sk1**(2.*n-1.)))/((sk**(2.*n+1))) 

                    end if

                end if
                
            end do

          end do

          do j = 2, s_div

            do n = 1, l_max 

                do k = 2, s_div

                    if (s_max .eq. 1. .and. (k .eq. s_div .or. j .eq. &
                        s_div)) then 

                        f_rho(j,n,k) = 0.d0
                        f_gamma(j,n,k) = 0.d0
                        f_omega(j,n,k) = 0.d0 

                    else 

                        sk = s_gp(k) 
                        sj = s_gp(j) 
                        sk1 = 1.-sk
                        sj1 = 1.-sj

                        if (k .lt. j) then

                            f_rho(j,n,k) = ((sj1/sj)**(2.*n+1.))* &
                                (sk**(2.*n))/(sk1**(2.*n+2))

                            f_gamma(j,n,k) = ((sj1/sj)**(2.*n))* &
                                (sk**(2.*n-1.))/(sk1**(2.*n+1))

                            f_omega(j,n,k) = ((sj1/sj)**(2.*n+1.))* &
                                (sk**(2.*n))/(sk1**(2.*n+2.))

                        else 

                            f_rho(j,n,k) = ((sj/sj1)**(2.*n))*(sk1** &
                                (2.*n-1.))/(sk**(2.*n+1)) 

                            f_gamma(j,n,k) = ((sj/sj1)**(2.*n-2.))* &
                                (sk1**(2.*n-3.))/(sk**(2.*n-1.)) 

                            f_omega(j,n,k) = ((sj/sj1)**(2.*n-2.))* &
                                (sk1**(2.*n-3.))/(sk**(2.*n-1.)) 

                        end if

                    end if

                end do

            end do

          end do

          n = 0 

          do i = 1, m_div

            p_2n(i,n) = legendre(2*n,mu(i)) 

          end do

          do i = 1, m_div

            do n = 1, l_max

                p_2n(i,n) = legendre(2*n,mu(i)) 
                p1_2n_1(i,n) = plgndr(2*n-1,1,mu(i))

            end do

          end do

      end subroutine 

!--------------------------------------------------------------------------
! legendre polynomial of degree n
!--------------------------------------------------------------------------

      double precision function legendre(n,x) 

          implicit none
          integer :: n, i
          double precision :: x, p_2, p_1, p

          p_2 = 1.
          p_1 = x

          if (n .ge. 2) then 

              do i = 2, n

                p = (x*(2.*i-1.)*p_1-(i-1.)*p_2)/(1.*i)
                p_2 = p_1
                p_1 = p 

              end do

              legendre = p

          else if (n .eq. 1) then 

              legendre = p_1

          else 

              legendre = p_2 

          end if

      end function

!--------------------------------------------------------------------------
! associated legendre polynomial 
!--------------------------------------------------------------------------

      double precision function plgndr(l,m,x) 

          implicit none
          integer :: l, m, ll, i
          double precision :: x, somx2, fact, pmm, pmmp1, pll

          if (m .lt. 0 .or. m .gt. l .or. dabs(x) .gt. 1.) then

              print*, 'bad argument in routine plgndr' 

          end if

          pmm = 1. 

          if (m .gt. 0) then 

              somx2 = sqrt((1.-x)*(1.+x)) 
              fact = 1. 

              do i = 1, m 

                pmm = pmm*(-fact*somx2) 
                fact = fact+2. 

              end do

          end if

          if (l .eq. m) then 

              plgndr = pmm 

          else 

              pmmp1 = x*(2*m+1)*pmm 

              if (l .eq. (m+1)) then 

                  plgndr = pmmp1

                  do ll = (m+2), l 

                    pll = (x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(1.*ll-1.*m) 
                    pmm = pmmp1 
                    pmmp1 = pll

                  end do

                  plgndr = pll

              end if

          end if

      end function

