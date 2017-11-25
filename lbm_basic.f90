! Basic D2Q9 Lattice Boltzmann in fortran
! although some of this is generalized, a lot of it is hard coded for D2Q9
program lbm

    use io_routines, only : io_write
    implicit none

    integer :: n_iterations = 20000

    real    :: Re     = 200.0  ! Reynolds number.
    integer, parameter :: nx     = 520
    integer, parameter :: ny     = 180
    ! ly=ny-1.0;
    integer, parameter :: q      = 9      ! lattice type/size
    real    :: u_lattice    = 0.04   ! lattice fluid velocity
    real    :: radius = 20
    real    :: tau


    real :: f(q,nx,ny)
    real :: f_eq(q,nx,ny)
    real :: c_u(q,nx,ny)
    real :: u(2,nx,ny)

    real :: velocity(2,nx,ny)
    real :: rho(nx,ny) = 1
    logical :: obstacle(nx,ny) = .False.

    ! Lattice Constants
    integer :: c(2,q)        ! Lattice particle velocities.
    real :: w(q) = 1./36. ! Lattice weights

    ! points to the no slip lattice element for each index
    integer :: noslip(q) = [1,3,2,7,9,8,4,6,5]


    integer :: c_vals(3)    = [0,-1,1]
    integer :: i

    obstacle(:,1) = .True.

    tau    = 1.0 / (3. * (u_lattice * radius/Re) + 0.5) ! time relaxation
    print*, 'Tau',tau
    ! Not the most typical c configuration, will have to think about this some
    !
    ! c =  1 [[ 0,  0],       6  3  9
    !      2  [ 0, -1],        \ | /
    !      3  [ 0,  1],         \|/
    !      4  [-1,  0],       4--1--7
    !      5  [-1, -1],         /|\
    !      6  [-1,  1],        / | \
    !      7  [ 1,  0],       5  2  8
    !      8  [ 1, -1],
    !      9  [ 1,  1]]
    !
    do i=1,q
        c(1,i) = c_vals((i-1)/3+1)
        c(2,i) = c_vals(mod(i-1,3)+1)

        if (sum(abs(c(:,i))) < 1.1) w(i) = 1/9.0
    end do
    w(1) = 4.0/9.0
    ! w = [ 4/9,  1/9,  1/9,  1/9,  1/36,
    !      1/36,  1/9, 1/36, 1/36]

    velocity(1,:,:) = u_lattice
    velocity(2,:,:) = 0

    call equilibrium(rho, velocity, f_eq, c_u, w, c)

    f = f_eq
    u = velocity

    ! Time loop
    do i=1,n_iterations
        if (mod(i,100)==0) print*, i, "  of", n_iterations
        call update(f, u, c_u, velocity, f_eq, noslip, obstacle, w, c)

    end do

    call io_write("output_U.nc","u",  &
                  reshape(u, [ny,nx,2], order=[3,2,1]))

contains

    subroutine equilibrium(rho, vel, feq, cu, w, c)
        implicit none
        real,    intent(in)    :: rho(:,:)
        real,    intent(in)    :: vel(:,:,:)
        real,    intent(inout) :: feq(:,:,:)
        real,    intent(inout) :: cu(:,:,:)
        real,    intent(in)    :: w(:)
        integer, intent(in)    :: c(:,:)

        integer :: i,j,k, nc,nx,ny

        nc = size(feq,1)
        nx = size(feq,2)
        ny = size(feq,3)

        do j=1,ny
            do i=1,nx
                cu(:,i,j) = 3 * (c(1,:) * vel(1,i,j) + c(2,:) * vel(2,i,j))
            enddo
        enddo

        do k=1,ny
            do j=1,nx
                feq(:,j,k) = rho(j,k) * w * (1. + cu(:,j,k) + 0.5 * cu(:,j,k)**2  &
                              - 3./2. * (vel(1,j,k)**2 + vel(2,j,k)**2))
            enddo
        enddo

        ! print*, 'feq',feq(:,3,3)

    end subroutine equilibrium


    subroutine update(f, u, cu, boundary_velocity, feq, noslip, obstacle, w, c)
        implicit none
        real,    intent(inout) :: f(:,:,:)
        real,    intent(inout) :: u(:,:,:)
        real,    intent(inout) :: cu(:,:,:)
        real,    intent(inout) :: boundary_velocity(:,:,:)
        real,    intent(inout) :: feq(:,:,:)
        integer, intent(in)    :: noslip(:)
        logical, intent(in)    :: obstacle(:,:)
        real,    intent(in)    :: w(:)
        integer, intent(in)    :: c(:,:)

        integer :: i, j, k, n, nq, nx, ny
        real, allocatable :: f_temp(:,:,:), rho(:,:)

        integer :: on_right(3)  = [4,5,6]
        integer :: in_middle(3) = [1,2,3]
        integer :: on_left(3)   = [7,8,9]

        nq = size(f,1)
        nx = size(f,2)
        ny = size(f,3)
        allocate(f_temp(nq,nx,ny))
        allocate(rho(nx,ny))

        ! Set open outflow boundary condition on right wall
        do i=1,3
            f(on_right(i),nx,:) = f(on_right(i),nx-1,:)
        enddo

        ! compute density
        rho = sum(f,1)

        do k = 1, ny
            do j = 1, nx
                do i = 1, 2
                    u(i,j,k) = 0

                    do n = 1, nq
                        u(i,j,k) = u(i,j,k) + (c(i,n) * f(n,j,k))
                    enddo
                    u(i,j,k) = u(i,j,k) / rho(j,k)
                enddo
            enddo
        enddo

        ! Left wall: compute density from known populations.
        u(:,1,:) = boundary_velocity(:,1,:)
        rho(1,:) = 0
        do i = 1, 3
            rho(1,:) = rho(1,:) + 1./(1.-u(1,1,:)) * (f(in_middle(i),1,:) + 2. * f(on_right(i),1,:))
        enddo

        call equilibrium(rho, u, feq, cu, w, c)

        ! Left wall: Zou/He boundary condition.
        do i = 1, 3
            f(on_left(i),1,:) = feq(on_left(i),1,:)
            ! f(on_left(i),1,:) = f(on_right(i),1,:) + feq(on_left(i),1,:) - f(on_right(i),1,:)
        enddo

        f_temp = f - tau * (f - feq) ! Collision

        ! bounce back no slip wall boundaries
        do k = 1, ny
            do j = 1, nx
                do i = 1, nq
                    if (obstacle(j,k)) then
                        f_temp(i,j,k) = f(noslip(i),j,k)
                    endif
                enddo
            enddo
        enddo

        do k = 1, ny
            do j = 1, nx
                do i = 1, nq
                    f(i,j,k) = f_temp(i,mod((j-1)-c(1,i)+nx,nx)+1, mod((k-1)-c(2,i)+ny,ny)+1)
                enddo
            enddo
        enddo
        ! print*, 'f',f(:,1,1)
        ! print*, ""
    end subroutine update

end program lbm
