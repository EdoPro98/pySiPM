!---------------------------------------------------------------------------
!> @brief Generation of SiPM signal.
!> Generation of SiPM signal given starting time, signal shape and noise values.
!> Computes: @f$e^{-frac{t}{t_f}}-e^{-{t}{t_r}}@f$
!
!
!> @param t Time position of the signal
!> @param tr Rising time of the signal shape
!> @param tf Falling time of the signal shape
!> @param h  Relative signal height
!> @param sigpts Number of samples
!> @param gvar Multiplicative factor accounting for gain variation
!
!> @param s: Signal
!---------------------------------------------------------------------------
subroutine signalgenfortran(s, t, h, tf, tr, sigpts, gvar)
  implicit none
  integer(4)  :: t, sigpts, i
  real(4)     :: h, gvar, tf, tr, s(sigpts)
!f2py intent(in) t, h, gvar, tf, tr, sigpts
!f2py intent(out) s

  s = 0.
  forall (i = 1 : sigpts - t)
    s(i + t) = exp(-i / tf) - exp(-i / tr)
  end forall
  s = gvar * h * s
end subroutine signalgenfortran


subroutine rollfortran(out, vect, t, h, npt)
  implicit none
  integer(4)     :: t, npt
  real(4)        :: gvar, h, vect(npt), out(npt)
!f2py intent(in) t, gvar, h, vect
!f2py intent(hide), depend(vect) npt = len(vect)
!f2py intent(out) out

  out = cshift(vect, -t, dim=1)
  out(1 : t) = 0
  out = h * out
end subroutine rollfortran


!----------------------------------------------------------------------
!> @author
!> Edoardo Proserpio
!> @brief F2PY extension for Python. Random arrays generation in Fortran
!-----------------------------------------------------------------------
module frandom
  implicit none
  real(4), parameter, private  :: pi2 = 6.28318530718E+00
contains

  !-----------------------------------------------------------
  !> @brief Generation of random integers in range [0 - sup].
  !
  !> @param sup Superior limit for random generation
  !> @param n Number of elements to generate
  !
  !> @param out
  !----------------------------------------------------------
subroutine randint(out, sup, n)
  implicit none
  integer(4)     :: sup, n, out(n)
  real(4)        :: u(n)
!f2py intent(in) ncells, nsig
!f2py intent(out) out

  call random_number(u)
  out = floor(sup * u) + 1
end subroutine randint


!-------------------------------------------------------------------------
!> @brief Generation of gaussian random values using Box-Muller transform.
!
!> @param mu Mean value of gaussian distribution
!> @param sigma Standard deviation of gaussian distribution
!> @param n Number of elements to generate
!
!> @param out
!-------------------------------------------------------------------------
subroutine randn(out, mu, sigma, n)
  implicit none
  integer(4)          :: n
  real(4)             :: mu, sigma, out(n), temp(n)
!f2py intent(in) mu, sigma, n
!f2py intent(out) out

  call random_number(out)
  call random_number(temp)

  where (out .lt. tiny(out)) out = tiny(out)

  out = sqrt( -2.0E+00 * log(out)) * cos(pi2 * temp) * sigma + mu
end subroutine randn


!------------------------------------------------------
!> @brief Generation of poissonian random values.
!
!> @param mu Mean value of poissonian distribution
!> @param n Number of elements to generate
!
!> @param out
!-------------------------------------------------------
subroutine randpoiss(out, mu, n)
  implicit none
  integer(4)     :: n, out(n),i
  real(4)        :: mu, L, p(n), u
!f2py intent(in) mu, n
!f2py intent(out) out

  L = exp(-mu)
  out = 0
  p = 1
  i = 1

  do i=1,n
    do while(p(i) .gt. L)
      call random_number(u)
      out(i) = out(i) + 1
      p(i) = p(i) * u
    end do
  end do


  out = out - 1
end subroutine randpoiss


!-------------------------------------------------------
!> @brief Generation of exponential random values.
!
!> @param mu Mean value of exponential distribution
!> @param n Number of elements to generate
!
!> @param out
!--------------------------------------------------------
subroutine randexp(out, mu, n)
  implicit none
  integer(4)   :: n
  real(4)      :: mu, out(n)
!f2py intent(in) mu, n
!f2py intent(out) out

  call random_number(out)
  out = -log(out) * mu
end subroutine randexp

end module frandom


!-------------------------------------------------------
!> @brief Array sorting done in Fortran
!
!> @param array Array to be sorted
!> @param last Number of elements in the array
!--------------------------------------------------------
subroutine sortfortran(array,last)
  implicit none
  integer(4)       :: i, j, left, right, last
  real(4)          :: array(last)
  real(4)          :: temp, p, next
!f2py intent(inout) array
!f2py intent(hide), depend(array) last = len(array)

  p = 0.5 * (array(1) + array(last))
  if (array(1) .gt. array(last)) then
     temp = array(last)
     array(last) = array(1)
     array(1) = temp
  endif

  left = 1
  right = last
  temp = array(2)

  do i = 2,last - 1
     if (temp .lt. p) then
        do j = left, 1, -1
           if (array(j) .le. temp) exit
           array(j + 1) = array(j)
        end do
        array(j + 1) = temp
        temp = array(left + 2)
        left = left + 1
     else
        next = array(right - 1)
        do j = right, last
           if (array(j) .ge. temp) exit
           array(j - 1) = array(j)
        end do
        array(j - 1) = temp
        temp = next
        right = right - 1
     endif
  end do
end subroutine sortfortran

subroutine signalanalysisfortran(integral, peak, toa, tot, top, signalingate, sampling)

  real(4)     :: signalingate(:), sampling
  real(4)     :: integral, peak, toa, tot, top
!f2py intent(in) signalingate(:), sampling
!f2py intent(out) integral, peak, toa, tot, top

  integral = sum(signalingate)
  peak = maxval(signalingate)
  toa = findloc(signalingate > 1.5, .true., dim = 1)
  tot = count(signalingate > 1.5, dim = 1)
  top = maxloc(signalingate, dim = 1)
  integral = integral * sampling
  toa = toa * sampling
  tot = tot * sampling
  top = top * sampling
end subroutine signalanalysisfortran
