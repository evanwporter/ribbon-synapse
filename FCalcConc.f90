#include "fintrf.h"

module iter_mod
  use omp_lib
  implicit none
  
  private
  public :: iter_func

contains

  subroutine iter_func(r, nr, D, dt, current, it, C)
    real(kind=8), intent(in) :: r(nr), D, dt, current(it)
    integer(kind=8), intent(in) :: nr, it
    real(kind=8), intent(out) :: C(nr)
    
    real(kind=8) :: four_pi_D, t, r_i, r_i_squared, t_prime, diff_time, G
    integer(kind=8) :: i, j
    
    four_pi_D = 4.0d0 * acos(-1.0d0) * D
    t = it * dt
    
    !$omp parallel do private(r_i, r_i_squared, t_prime, diff_time, G) shared(C, r, current)
    do i = 1, nr
      r_i = r(i)
      r_i_squared = r_i * r_i
      C(i) = 0.0d0
      do j = 1, it
        t_prime = (j-1) * dt
        diff_time = t - t_prime
        G = (four_pi_D * diff_time)**(-1.5d0) * exp(-r_i_squared / (four_pi_D * diff_time))
        C(i) = C(i) + current(j) * G * dt
      end do
    end do
    !$omp end parallel do
    
  end subroutine iter_func

end module iter_mod

! Gateway function for MATLAB
subroutine mexFunction(nlhs, plhs, nrhs, prhs)
  use iter_mod
  implicit none

  integer*4, intent(in) :: nlhs, nrhs
  mwPointer, intent(inout) :: plhs(*), prhs(*)

  mwPointer :: r_ptr, current_ptr, C_ptr
  real(kind=8), pointer :: r(:), current(:), C(:)
  real(kind=8) :: D, dt
  integer(kind=8) :: nr, it

  if (nrhs .ne. 5) then
    call mexErrMsgIdAndTxt("MATLAB:iter_func:nrhs", "Five inputs required.")
  endif

  if (nlhs .ne. 1) then
    call mexErrMsgIdAndTxt("MATLAB:iter_func:nlhs", "One output required.")
  endif

  r_ptr = mxGetPr(prhs(1))
  nr = mxGetM(prhs(1)) * mxGetN(prhs(1))
  D = mxGetScalar(prhs(2))
  dt = mxGetScalar(prhs(3))
  current_ptr = mxGetPr(prhs(4))
  it = mxGetScalar(prhs(5))

  call mxCopyPtrToReal8(r_ptr, r, nr)
  call mxCopyPtrToReal8(current_ptr, current, it)
  
  C_ptr = mxCreateDoubleMatrix(nr, 1, 0)
  call mxCopyPtrToReal8(mxGetPr(C_ptr), C, nr)

  call iter_func(r, nr, D, dt, current, it, C)
  
  call mxCopyReal8ToPtr(C, mxGetPr(C_ptr), nr)
  plhs(1) = C_ptr
end subroutine mexFunction
