! Basic D2Q9 Lattice Boltzmann in fortran
! although some of this is generalized, a lot of it is hard coded for D2Q9
program lbm

    use io_routines, only : io_write, io_read
    implicit none

    integer :: n_iterations = 100000

    real,    parameter :: Re        = 300.0  ! Reynolds number. = u L / v (wind speed * length / viscosity)
    real,    parameter :: radius    = 20     ! Effective object radius in grid-cells for use with Reynolds number

    real,    parameter :: viscosity = 1.81E-5! Viscosity of Air at 15 C           [kg / (m s)]
    real,    parameter :: kviscosity= 1.48E-5! Kinematic Viscosity of Air at 15 C [m^2 / s]

    real,    parameter :: dx        = 1.0    ! grid cell width          [m]
    real,    parameter :: u_sound   = 343.0  ! Speed of sound in air    [m/s]
    real,    parameter :: dt        = dx / u_sound ! technically the lattice speed of sound is ~1 / sqrt(3) (?)

    real,    parameter :: u_lattice = 0.04   ! lattice fluid velocity = wind speed / dx * dt
    integer, parameter :: q         = 9      ! lattice type/size
    integer, parameter :: nx        = 1024   ! lattice length
    integer, parameter :: ny        = 384    ! lattice height

    integer, parameter :: output_dt = 1      ! output interval [s]
    integer, parameter :: output_steps = nint(output_dt / dt)    ! output interval time steps
    character(len=1024):: output_filename

    ! Lattice Constants
    integer :: c(2,q)        ! Lattice particle velocities.
    real    :: w(q) = 1./36. ! Lattice weights
    ! points to the no slip lattice element for each index
    integer :: noslip(q) = [1,3,2,7,9,8,4,6,5]

    integer :: c_vals(3)    = [0,-1,1]

    ! lattice relaxation time used in BGK collision operator
    real    :: tau

    ! convenience variables
    integer :: start_x, ymax, i

    ! lattice distribution function
    real :: f(q,nx,ny)
    ! equilibrium distribution function
    real :: f_eq(q,nx,ny)
    ! lattice speed in each direction (q) at all grid points
    real :: c_u(q,nx,ny)
    ! real velocity at each grid cell in x and y directions
    real :: u(2,nx,ny)
    real :: output_u(2,nx,ny)

    ! initial velocity (used for boundary conditions)
    real :: velocity(2,nx,ny)
    ! density
    real :: rho(nx,ny) = 1
    ! lower no-slip boundary location
    real, allocatable :: topography(:)
    ! defines the no-slip surface
    logical :: obstacle(nx,ny) = .False.

    print*, (dt * kviscosity) / (0.577**2) + 0.5,  &
            dt * kviscosity / (0.577**2)

    ! set up the noslip surface
    call io_read("profile.nc","data",topography)
    obstacle(:,1) = .True.

    topography = topography - minval(topography)
    start_x = 250 !size(topography,1) - nx
    do i=1,nx
        ymax = nint(topography(min(i + start_x,size(topography,1))))
        if (ymax > 0) obstacle(i,:ymax)=.True.
    enddo

    tau    = 1.0 / (3. * (u_lattice * radius/Re) + 0.5) ! 1 / time relaxation (collision operator)
    ! tau = 3 * kviscosity * (dt / dx**2) + 0.5
    print*, 'Tau',1/tau

    ! Not the most typical c configuration, will have to think about this some, maybe it doesn't matter...
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

        if (sum(abs(c(:,i))) > 1.1)     w(i) = 1 /36.0
        if ((sum(abs(c(:,i))) < 1.1)    &
        .and.(sum(abs(c(:,i))) > 0.1))  w(i) = 1 / 9.0
        if (sum(abs(c(:,i))) < 0.1)     w(i) = 4 / 9.0
    end do
    ! w = [ 4/9,  1/9,  1/9,  1/9,  1/36,
    !      1/36,  1/9, 1/36, 1/36]

    velocity(1,:,:) = u_lattice
    velocity(2,:,:) = 0

    call equilibrium(rho, velocity, f_eq, c_u, w, c)

    f = f_eq
    u = velocity
    output_u = 0

    print*, "Output every ", output_steps
    ! Time loop
    do i=1,n_iterations

        call update(f, u, c_u, velocity, f_eq, noslip, obstacle, w, c)

        output_u = output_u + u

        if (mod(i,output_steps)==0) then

            if (i/output_steps > 250) then
                write(output_filename,"(A,I4.4,A)") "output/output_",i/output_steps,".nc"
                print*, "Writing outputfile:",trim(output_filename)

                call io_write(trim(output_filename),"u",  &
                              reshape(reshape((output_u / output_steps) * (dx/dt), [nx,ny,2], order=[3,1,2]), [nx,ny,2,1]))
            endif
            output_u = 0
        endif

    end do

contains

    subroutine equilibrium(rho, vel, feq, cu, w, c)
        implicit none
        real,    intent(in)    :: rho(:,:)
        real,    intent(in)    :: vel(:,:,:)
        real,    intent(inout) :: feq(:,:,:)
        real,    intent(inout) :: cu(:,:,:)
        real,    intent(in)    :: w(:)
        integer, intent(in)    :: c(:,:)

        integer :: i,j, nx,ny

        nx = size(feq,2)
        ny = size(feq,3)

        do j=1,ny
            do i=1,nx
                cu(:,i,j) = 3 * (c(1,:) * vel(1,i,j) + c(2,:) * vel(2,i,j))

                feq(:,i,j) = rho(i,j) * w * (1. + cu(:,i,j) + 0.5 * cu(:,i,j)**2  &
                              - 3./2. * (vel(1,i,j)**2 + vel(2,i,j)**2))
            enddo
        enddo

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

        !$omp parallel default(shared) private(i,j,k)
        !$omp do
        do k = 1, ny
            ! Left wall: compute density from known populations.
            u(:,1,k) = boundary_velocity(:,1,k)
            rho(1,k) = 0

            do i = 1,3
                rho(1,k) = rho(1,k) + 1./(1.-u(1,1,k)) * (f(in_middle(i),1,k) + 2. * f(on_right(i),1,k))
            enddo
            do j = 1,nx
                if (j>1) then
                    if (.not.obstacle(j,k)) then
                        do i = 1,2
                            u(i,j,k) = 0

                            do n = 1,nq
                                u(i,j,k) = u(i,j,k) + (c(i,n) * f(n,j,k))
                            enddo
                            u(i,j,k) = u(i,j,k) / rho(j,k)
                        enddo
                    else
                        u(:,j,k)=0
                    endif
                endif

                ! call equilibrium(rho, u, feq, cu, w, c) is effectively this
                cu(:,j,k) = 3 * (c(1,:) * u(1,j,k) + c(2,:) * u(2,j,k))

                feq(:,j,k) = rho(j,k) * w * (1. + cu(:,j,k) + 0.5 * cu(:,j,k)**2  &
                              - 3./2. * (u(1,j,k)**2 + u(2,j,k)**2))


                if (j==1) then
                    ! Left wall: Zou/He boundary condition.
                    do i = 1,3
                      f(on_left(i),1,k) = feq(on_left(i),1,k)
                      ! f(on_left(i),1,:) = f(on_right(i),1,:) + feq(on_left(i),1,:) - f(on_right(i),1,:)
                    enddo
                endif

                ! BGK Collision
                f_temp(:,j,k) = f(:,j,k) - tau * (f(:,j,k) - feq(:,j,k))

                ! Bounce back no slip wall boundaries
                if (obstacle(j,k)) then
                    do i = 1,nq
                        f_temp(i,j,k) = f(noslip(i),j,k)
                    enddo
                endif

            enddo
        enddo
        !$omp end do
        !$omp barrier
        ! separate loop because we have to know that all collisions / bounceback processes are accounted for before streaming
        !$omp do
        do k = 1, ny
            do j = 1,nx
                do i = 1,nq
                    f(i,j,k) = f_temp(i,mod((j-1)-c(1,i)+nx,nx)+1, mod((k-1)-c(2,i)+ny,ny)+1)
                enddo
            enddo
        enddo
        !$omp end do
        !$omp end parallel
        ! print*, 'f',f(:,1,1)
        ! print*, ""
    end subroutine update

end program lbm
