module param
    integer, parameter :: rk= selected_real_kind(18)
    real(kind=rk), parameter :: pi = 4*atan(1.0)
    real(kind=rk), parameter :: hbar = 1
    real(kind=rk), parameter :: beta = 1
    real(kind=rk), parameter :: G = 1
    real(kind=rk), parameter :: gam = 1
    real(kind=rk), parameter :: R0 = 1
    real(kind=rk), parameter :: m = 1
    real(kind=rk), parameter :: w = 1
    real(kind=rk), parameter :: qzpf = sqrt(hbar/(2.0*m*w))
    real(kind=rk), parameter :: pzpf = sqrt(hbar*m*w/2.0)
end module param

module matrixOp
    use param
    implicit none
    contains

    ! Expansion of determinants using Laplace formula
    recursive function determinant(matrix) result(laplace_det)
    real(kind=rk), dimension(:,:) :: matrix
    integer :: msize(2), i, n
    real(kind=rk) :: laplace_det, det
    real(kind=rk), dimension(:,:), allocatable :: cf

    msize = shape(matrix)
    n = msize(1)

    if (n .eq. 1) then
      det = matrix(1,1)
    else
      det = 0
      do i=1, n
        allocate(cf(n-1, n-1))
        cf = cofactor(matrix, i, 1)
        det = det + ((-1)**(i+1))* matrix(i,1) * determinant(cf)
        deallocate(cf)
      end do
    end if
    laplace_det = det
    end function determinant

    function cofactor(matrix, mI, mJ)
    use param
      real(kind=rk), dimension(:,:) :: matrix
      integer :: mI, mJ
      integer :: msize(2), i, j, k, l, n
      real(kind=rk), dimension(:,:), allocatable :: cofactor
      msize = shape(matrix)
      n = msize(1)

      allocate(cofactor(n-1, n-1))
      l=0
      k = 1
      do i=1, n
       if (i .ne. mI) then
         l = 1
         do j=1, n
           if (j .ne. mJ) then
             cofactor(k,l) = matrix(i,j)
             l = l+ 1
           end if
         end do
         k = k+ 1
       end if
      end do
    return
    end function cofactor

    subroutine inverse(a,c,n)
    use param
    !============================================================
    ! Inverse matrix
    ! Method: Based on Doolittle LU factorization for Ax=b
    ! Alex G. December 2009
    !-----------------------------------------------------------
    ! input ...
    ! a(n,n) - array of coefficients for matrix A
    ! n      - dimension
    ! output ...
    ! c(n,n) - inverse matrix of A
    ! comments ...
    ! the original matrix a(n,n) will be destroyed
    ! during the calculation
    !===========================================================
    implicit none
    integer n
    real(kind=rk) a(n,n), c(n,n)
    real(kind=rk) L(n,n), U(n,n), b(n), d(n), x(n)
    real(kind=rk) coeff
    integer i, j, k

    ! step 0: initialization for matrices L and U and b
    ! Fortran 90/95 allows such operations on matrices
    L=0.0
    U=0.0
    b=0.0

    ! step 1: forward elimination
    do k=1, n-1
       do i=k+1,n
          coeff=a(i,k)/a(k,k)
          L(i,k) = coeff
          do j=k+1,n
             a(i,j) = a(i,j)-coeff*a(k,j)
          end do
       end do
    end do

    ! Step 2: prepare L and U matrices
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i=1,n
      L(i,i) = 1.0
    end do
    ! U matrix is the upper triangular part of A
    do j=1,n
      do i=1,j
        U(i,j) = a(i,j)
      end do
    end do

    ! Step 3: compute columns of the inverse matrix C
    do k=1,n
      b(k)=1.0
      d(1) = b(1)
    ! Step 3a: Solve Ld=b using the forward substitution
      do i=2,n
        d(i)=b(i)
        do j=1,i-1
          d(i) = d(i) - L(i,j)*d(j)
        end do
      end do
    ! Step 3b: Solve Ux=d using the back substitution
      x(n)=d(n)/U(n,n)
      do i = n-1,1,-1
        x(i) = d(i)
        do j=n,i+1,-1
          x(i)=x(i)-U(i,j)*x(j)
        end do
        x(i) = x(i)/u(i,i)
      end do
    ! Step 3c: fill the solutions x(n) into column k of C
      do i=1,n
        c(i,k) = x(i)
      end do
      b(k)=0.0
    end do
    end subroutine inverse

end module matrixOp

module initialization
    use param
    use matrixOp
    implicit none
    contains

    function gaussianWF(r,sigma) result(res)
        real(kind=rk), dimension(6,1), intent(in):: r
        real(kind=rk), dimension(6,1) :: rd
        real(kind=rk), dimension(1,6) :: rdT
        real(kind=rk), dimension(1,1) :: expo_array
        real(kind=rk) :: expo
        real(kind=rk), dimension(6,6), intent(in) :: sigma !INVERSE of the covariance matrix
        real(kind=rk) :: res
        integer :: i,j


        !switching to quadratures
        do i = 1,3
            rd(i,1) = r(i,1)/ qzpf
        end do

        do j = 4,6
            rd(j,1) = r(j,1)/ pzpf
        end do

        rdT = transpose(rd)

        expo_array = matmul(matmul(rdT,sigma),rd)
        expo = expo_array(1,1)

        res = sqrt(determinant(sigma))/sqrt(2*pi)**3 * 1.0 * exp(-expo)

    end function gaussianWF

    subroutine const(model, f,D1,D2,D3,D4) ! model = 1 -> CSL, model = 2 -> DP
        integer :: model
        real(kind=rk) :: I, f, D1, D2, D3, D4

        if (model .eq. 1) then
            I = 3*hbar**2*gam/sqrt(pi**(3))/16.0/R0**5
            else if (model .eq. 2) then
                I = hbar*G/2.0/sqrt(pi)/R0**3
                else
                    print*, "No valid model index was selected"
        end if

        D1 = 2*m**2.0/3.0 * I + hbar**4*beta**2*35.0/96.0 * I - hbar**2*beta*m*5.0/6.0 * I
        D2 = hbar**2*beta**2/48.0 * I
        D3 = 2.0*D2
        D4 = D3
        f = m*beta/3.0 * I + hbar**2*beta**2/12.0*5.0/2.0*I

    end subroutine const

end module initialization

module heatEq !test
    implicit none
    contains

    !gaussian initial datum
    function initialDatum(x, mu, var) result(res)
        use param
        real(kind=rk), intent(in) :: x,mu,var
        real(kind=rk) :: res
        res = 1.0/sqrt(2.0*4*atan(1.0)*var)*exp(-(x-mu)**2/2.0/var)
    end function initialDatum

    function pdeHE(NumSteps, GrdPts, dt, dx) result(f_of_t)
        !Solving the pde of the heat + dissipation equation
        use param
        integer, intent(in) :: NumSteps, GrdPts
        real(kind=rk), intent(in) :: dt, dx
        real(kind=rk), dimension(NumSteps,GrdPts):: f_of_t
        real(kind=rk) :: ave, s2, r, f
        real(kind=rk), dimension(GrdPts):: Grd
        real(kind=rk), dimension(NumSteps) :: tAxis
        integer :: i,j,k,l,n

        !creating the space grid
        do i=1,GrdPts
            Grd(i) = (i-GrdPts/2)*dx
        end do

        !creating the time axis
        do j=1,NumSteps
            tAxis(j) = j*dt
        end do

        ave = 0
        s2 = 0.1

        !initialization of the initial datum
        do k=1,GrdPts
            f_of_t(1,k) = initialDatum(Grd(k),ave,s2)
        end do

        !heat equation with the finite difference method

        r = 0.01*dt/dx**2
        f = 0.3

        do l=2,NumSteps
            do n=2,GrdPts-1
                f_of_t(l,n) = (1.0-2*r + f*dt)*f_of_t(l-1,n) + &
                (r+f*dt*Grd(n)/(2*dx))*f_of_t(l-1,n+1) + &
                (r-f*dt*Grd(n)/(2*dx))*f_of_t(l-1,n-1)
            end do
        end do

    end function pdeHE

end module heatEq

module dissipativeQFP
    use param
    use initialization
    implicit none
    contains

    function dyn(NumSteps, GrdPts, dt, dx) result(f_of_t)
        integer, intent(in) :: NumSteps, GrdPts
        real(kind=rk), intent(in) :: dt, dx
        real(kind=rk), dimension(NumSteps,GrdPts,GrdPts,GrdPts,GrdPts,GrdPts,GrdPts):: f_of_t
        real(kind=rk), dimension(6):: r
        real(kind=rk), dimension(6,6) :: sigma !inverse of the covariance matrix of the initial datum
        real(kind=rk) :: f, D1, D2, D3, D4
        real(kind=rk), dimension(GrdPts):: Grd
        real(kind=rk), dimension(NumSteps) :: tAxis

        integer :: i,j,t
        integer :: i1, i2, i3, i4, i5, i6
        integer :: j1, j2, j3, j4, j5, j6

        call const(1, f, D1, D2, D3, D4)

        !creating the grid of one of the dimensions
        do i=1,GrdPts
            Grd(i) = (i-GrdPts/2)*dx
        end do

        !creating the time axis
        do j=1,NumSteps
            tAxis(j) = j*dt
        end do

        !INVERSE OF THE COVARIANCE MATRIX OF THE INITIAL DATUM
        do i=1,6
            do j=1,6
                sigma(i,j) = 10.0
            end do
        end do

        !initialization of the initial datum
        do i1=1,GrdPts
            do i2=1,GrdPts
                do i3=1,GrdPts
                    do i4=1,GrdPts
                        do i5=1,GrdPts
                            do i6=1,GrdPts
                                r = [Grd(i1), Grd(i2), Grd(i3), Grd(i4), Grd(i5), Grd(i6)]
                                f_of_t(1,i1,i2,i3,i4,i5,i6) = gaussianWF(r,sigma)
                            end do
                        end do
                    end do
                end do
            end do
        end do

        do t=2,NumSteps
            do j1=2,GrdPts-1
                do j2=2,GrdPts-1
                    do j3=2,GrdPts-1
                        do j4=2,GrdPts-1
                            do j5=2,GrdPts-1
                                do j6=2,GrdPts-1

                                    f_of_t(t,j1,j2,j3,j4,j5,j6) = f_of_t(t-1,j1,j2,j3,j4,j5,j6) + dt*()!ROBE

                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do

    end function dyn

end module dissipativeQFP



program main
    use param
    use dissipativeQFP
    implicit none

    character(len=100) :: filepath

    real(kind=rk) :: deltat, d !the spacing between the points of the phase-space grid is isotropic
    integer :: i, n
    integer, parameter :: GPts = 10
    integer, parameter :: Nsteps = 500

    real(kind=rk), dimension(GPts):: Grid
    real(kind=rk), dimension(Nsteps) :: times
    real(kind=rk), dimension(Nsteps, Gpts, GPts) :: sol !the solution projected on the (q1,p1) surface

    deltat = 0.01
    d = 0.02
    !the idea is to use the less number of 6-dim arrays as possible
    !creating the space grid (one-dimensional)
    do i=1,GPts
        Grid(i) = (i-Gpts/2)*d
    end do

    !creating the time axis
    do n=1,Nsteps
        times(n) = n*deltat
    end do

    ! Specify directory and file name
    !filepath = 'C:\Users\simon\OneDrive\Desktop\fortran programming\solt.txt'

    !open(unit=10, file=trim(filepath), status='replace')

    !do i = 1,GPts
    !    write(10,*) Grid(i), sol(1,i)
    !end do
    !close(unit=10)

end program main
