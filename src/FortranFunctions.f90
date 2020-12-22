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
pure subroutine fsignal(s, t, tf, tr, h, sigpts)
  implicit none
  integer(8),intent(in)  :: t, sigpts
  real(8),intent(in)     :: tf, tr, h
  real(8),intent(out)    :: s(sigpts)
  integer(8)             :: i
!f2py intent(in) t, h, gvar, tf, tr, sigpts
!f2py intent(out) s

  s = 0.
  forall (i = 1 : sigpts - t)
    s(i + t) = exp(-i / tf) - exp(-i / tr)
  end forall
  s = s * h
end subroutine fsignal


pure subroutine froll(out, vect, t, h, npt)
  implicit none
  integer(8),intent(in)     :: t, npt
  real(8),intent(in)        :: vect(npt), h
  real(8),intent(out)       :: out(npt)
!f2py intent(in) t, vect
!f2py intent(hide), depend(vect) npt = vect.size
!f2py intent(out) out

  out = cshift(vect, -t, dim=1)
  out(1 : t) = 0
  out = out * h
end subroutine froll


!----------------------------------------------------------------------
!> @author
!> Edoardo Proserpio
!> @brief F2PY extension for Python. Random arrays generation in Fortran
!-----------------------------------------------------------------------
module frandom
  implicit none
  real(8), parameter, private  :: pi2 = 6.28318530718E+00
contains

  !-----------------------------------------------------------
  !> @brief Generation of random integers in range [0 - sup].
  !
  !> @param sup Superior limit for random generation
  !> @param n Number of elements to generate
  !
  !> @param out
  !----------------------------------------------------------
subroutine integer(out, sup, n)
  implicit none
  integer(8)     :: sup, n, out(n)
  real(8)        :: u(n)
!f2py intent(in) sup,n
!f2py intent(out) out
  sup = sup + 1
  call random_number(u)
  out = floor(sup * u)
end subroutine integer


!-------------------------------------------------------------------------
!> @brief Generation of gaussian random values using Box-Muller transform.
!
!> @param mu Mean value of gaussian distribution
!> @param sigma Standard deviation of gaussian distribution
!> @param n Number of elements to generate
!
!> @param out
!-------------------------------------------------------------------------
subroutine normal(out, mu, sigma, n)
  implicit none
  integer(8)          :: n
  real(8)             :: mu, sigma
  real(8)             :: out(n), v(n), u(n), logs(n)
!f2py intent(in) mu, sigma, n
!f2py intent(out) out

  call random_number(v)
  call random_number(u)

  where (v .lt. tiny(v)) v = tiny(v)
  u = u * pi2

  logs = log(v)
  logs = sqrt(-2.00E0 * logs)

  out = logs * cos(u)
  out = out * sigma + mu
end subroutine normal


!------------------------------------------------------
!> @brief Generation of poissonian random values.
!
!> @param mu Mean value of poissonian distribution
!> @param n Number of elements to generate
!
!> @param out
!-------------------------------------------------------
subroutine poisson(out, mu, n)
  implicit none
  integer(8)     :: n, out(n)
  real(8)        :: mu, L, p(n), u
!f2py intent(in) mu, n
!f2py intent(out) out

  L = exp(-mu)
  out = -1
  p = 1

  do while(any(p .gt. L))
    where (p .gt. L) out = out + 1
    call random_number(u)
    p = p * u
  end do
end subroutine poisson


!-------------------------------------------------------
!> @brief Generation of exponential random values.
!
!> @param mu Mean value of exponential distribution
!> @param n Number of elements to generate
!
!> @param out
!--------------------------------------------------------
subroutine exponential(out, mu, n)
  implicit none
  integer(8)   :: n
  real(8)      :: mu, out(n)
!f2py intent(in) mu, n
!f2py intent(out) out

  call random_number(out)
  out = -log(out) * mu
end subroutine exponential

end module frandom


!-------------------------------------------------------
!> @brief Array sorting done in Fortran
!
!> @param array Array to be sorted
!> @param n Number of elements in the array
!--------------------------------------------------------
pure subroutine fsort(array,n)
  implicit none
  integer(8)            :: i, j, left, right
  integer(8),intent(in) :: n
  real(8),intent(inout) :: array(n)
  real(8)               :: temp, p, next
!f2py intent(inout) array
!f2py intent(hide), depend(array) n = len(array)

  p = 0.5 * (array(1) + array(n))
  if (array(1) .gt. array(n)) then
     temp = array(n)
     array(n) = array(1)
     array(1) = temp
  endif

  left = 1
  right = n
  temp = array(2)

  do i = 2,n - 1
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
        do j = right, n
           if (array(j) .ge. temp) exit
           array(j - 1) = array(j)
        end do
        array(j - 1) = temp
        temp = next
        right = right - 1
     endif
  end do
end subroutine fsort


subroutine signalanalysisfortran(integral, peak, toa, tot, top, signalingate, sampling)

  real(8)     :: signalingate(:), sampling
  real(8)     :: integral, peak, toa, tot, top
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
