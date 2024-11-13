MODULE global_parameters
IMPLICIT NONE
DOUBLE PRECISION, PARAMETER :: pi=3.1415926535897932384D0
DOUBLE PRECISION, PARAMETER :: TOL=5.D-4 
type node
	DOUBLE PRECISION :: x
	DOUBLE PRECISION :: y
end type
type(node),allocatable :: particle(:)     !particle center 
type(node),allocatable :: force_p(:)     ! inter particle force
type(node),allocatable :: v_p(:)     ! velocity

Integer :: NMCs, Npre, N_pre, Npart, N_total, na, nb, nc, nt, r_3, nmax, Ntype, layout_type
Integer :: nx, ny 
Integer :: seed
Integer,allocatable :: i_x(:), i_y(:)     ! cell number of each particle
Integer,allocatable :: charge(:)     !charge
Integer,allocatable :: cell(:,:,:)     !cell list
Integer,allocatable :: pointor(:,:)     !pointor
integer, DIMENSION(:), ALLOCATABLE :: nn

DOUBLE PRECISION :: Lx, Ly, r_b, r_c 
DOUBLE PRECISION, allocatable :: shift_x(:,:), shift_y(:,:)
DOUBLE PRECISION, allocatable :: phiyy(:), phi_y(:)
DOUBLE PRECISION :: axis(3)
DOUBLE PRECISION :: ff, f_b  ! external field
DOUBLE PRECISION :: dd, dd_2, kappa, kWCA
DOUBLE PRECISION :: mobility, vy, vxx, vyy
DOUBLE PRECISION :: vvx, vvy
DOUBLE PRECISION :: deltat, eta, dt_f, dt_r, sigma 

END MODULE global_parameters