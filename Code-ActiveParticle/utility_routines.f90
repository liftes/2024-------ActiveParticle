MODULE utility_routines
IMPLICIT NONE
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE initialize_C_particles(a_plane_2, layout_type)
    USE global_parameters, ONLY: particle, na, nb, N_total, Lx, Ly
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: layout_type
    DOUBLE PRECISION :: a_plane_2, r_c
    INTEGER :: i
    REAL :: dx, dy, rand_val

    r_c = a_plane_2
    SELECT CASE(layout_type)
    CASE(1)  ! Square grid
        dx = a_plane_2
        dy = a_plane_2
        DO i = na+nb+1, N_total
            particle(i)%x = MOD(i - na - nb - 1, CEILING((Lx - 2*r_c)/dx)) * dx - Lx/2.0 + r_c
            particle(i)%y = (i - na - nb - 1) / CEILING(Lx/dx) * dy
            IF (particle(i)%x < 0.0) THEN
                particle(i)%x = particle(i)%x + Lx
            END IF
            IF (particle(i)%x > Lx) THEN
                particle(i)%x = particle(i)%x - Lx
            END IF
            particle(i)%y = MOD(particle(i)%y, Ly - Lx)
        END DO

    CASE(2)  ! Hexagonal grid
        dx = a_plane_2
        dy = a_plane_2 * SQRT(3.0) / 2.0
        DO i = na+nb+1, N_total
            particle(i)%x = MOD(i - na - nb - 1, CEILING((Lx - 2*r_c)/dx)) * dx + &
            MOD((i - na - nb - 1) / CEILING(Lx/dx), 2) * dx / 2.0 + r_c - Lx/2.0 
            particle(i)%y = (i - na - nb - 1) / CEILING(Lx/dx) * dy
            IF (particle(i)%x < 0.0) THEN
                particle(i)%x = particle(i)%x + Lx
            END IF
            particle(i)%y = MOD(particle(i)%y, Ly - Lx)
        END DO

    CASE(3)  ! Rectangular grid
        dx = a_plane_2
        dy = a_plane_2 * 1.5
        DO i = na+nb+1, N_total
            particle(i)%x = MOD(i - na - nb - 1, CEILING((Lx - 2*r_c)/dx)) * dx - Lx/2.0 + r_c
            particle(i)%y = (i - na - nb - 1) / CEILING(Lx/dx) * dy
            IF (particle(i)%x < 0.0) THEN
                particle(i)%x = particle(i)%x + Lx
            END IF
            particle(i)%y = MOD(particle(i)%y, Ly - Lx)
        END DO

    CASE(4)  ! Centered rectangular grid
        dx = a_plane_2
        dy = a_plane_2 * 1.5
        DO i = na+nb+1, N_total
            ! Correct calculation for particles per row considering boundary adjustment
            particle(i)%x = MOD(i - na - nb - 1, CEILING((Lx - 2*r_c)/dx)) * dx + r_c - Lx/2.0
            particle(i)%y = (i - na - nb - 1) / CEILING((Lx - 2*r_c)/dx) * dy
            ! Staggering rows
            IF (MOD((i - na - nb - 1) / CEILING((Lx - 2*r_c)/dx), 2) == 1) THEN
                particle(i)%x = particle(i)%x + dx / 2.0
            END IF
            ! Boundary conditions
            IF (particle(i)%x < 0.0) THEN
                particle(i)%x = particle(i)%x + Lx
            ELSEIF (particle(i)%x > Lx) THEN
                particle(i)%x = particle(i)%x - Lx
            END IF
            ! Modulo for y to handle wrapping
            particle(i)%y = MOD(particle(i)%y, Ly - Lx)
        END DO

    CASE(5)  ! Skewed grid
        dx = a_plane_2
        dy = a_plane_2
        DO i = na+nb+1, N_total
            particle(i)%x = MOD(i - na - nb - 1, CEILING((Lx - 2*r_c)/dx)) * dx + &
            MOD((i - na - nb - 1) / CEILING(Lx/dx), 2) * dx / 4.0 + r_c - Lx/2.0
            particle(i)%y = (i - na - nb - 1) / CEILING(Lx/dx) * dy
            IF (particle(i)%x < 0.0) THEN
                particle(i)%x = particle(i)%x + Lx
            END IF
            particle(i)%y = MOD(particle(i)%y, Ly - Lx)
        END DO

    CASE(6)  ! Random distribution
        DO i = na+nb+1, N_total
            CALL RANDOM_NUMBER(rand_val)
            particle(i)%x = rand_val * (Lx - 2*r_c) + r_c - Lx/2.0
            CALL RANDOM_NUMBER(rand_val)
            particle(i)%y = rand_val * (Ly - Lx)
            IF (particle(i)%x < 0.0) THEN
                particle(i)%x = particle(i)%x + Lx
            END IF
        END DO
    END SELECT
END SUBROUTINE initialize_C_particles


SUBROUTINE pcf ()
USE global_parameters
IMPLICIT NONE
DOUBLE PRECISION :: dis,vbar
INTEGER :: i, j 
vbar = 0 
vy = 0
do i=1, nn(1)
    dis = vbar + v_p( i )%x
    vy = vy + v_p( i )%y                 
end do ! all particle

vbar = 1.0d0*vbar/nn(1)
mobility = vbar/ff/eta
vy = 1.0d0*vy/nn(1) 

vxx = 0 
vyy = 0
do i=1, nn(1)
    vxx = vxx + (v_p( i )%x-vbar)**2
    vyy = vyy + (v_p( i )%y-vy)**2                 
end do ! all particle
vxx = 1.0d0*vxx/nn(1)
vyy = 1.0d0*vyy/nn(1)

vvx = 0 
vvy = 0
do i=1, nn(1)
    do j = 1, nn(1)
        vvx = vvx + v_p( i )%x*v_p( i )%x
        vvy = vvy + v_p( i )%y*v_p( i )%y
    end do                 
end do ! all particle
vvx = 1.0d0*vvx/nn(1)/nn(1) - vbar*vbar
vvy = 1.0d0*vvy/nn(1)/nn(1) - vy*vy

End SUBROUTINE pcf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE v_correlate ()
USE global_parameters
IMPLICIT NONE
DOUBLE PRECISION :: vbar
INTEGER :: i, j 
vbar = 0 
vy = 0
do i=1, nn(1)
    vbar = vbar + v_p( i )%x
    vy = vy + v_p( i )%y                 
end do ! all particle

vbar = 1.0d0*vbar/nn(1)
mobility = vbar/ff/eta
vy = 1.0d0*vy/nn(1) 

vxx = 0 
vyy = 0
do i=1, nn(1)
    vxx = vxx + (v_p( i )%x-vbar)**2
    vyy = vyy + (v_p( i )%y-vy)**2                 
end do ! all particle
vxx = 1.0d0*vxx/nn(1)
vyy = 1.0d0*vyy/nn(1)

vvx = 0 
vvy = 0
do i=1, nn(1)
    do j = 1, nn(1)
        vvx = vvx + v_p( i )%x*v_p( i )%x
        vvy = vvy + v_p( i )%y*v_p( i )%y
    end do                 
end do ! all particle
vvx = 1.0d0*vvx/nn(1)/nn(1) - vbar*vbar
vvy = 1.0d0*vvy/nn(1)/nn(1) - vy*vy

End SUBROUTINE v_correlate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE order_w (w, vbar)
USE global_parameters
IMPLICIT NONE
DOUBLE PRECISION :: v_x(1:N_total), v_y(1:N_total)
DOUBLE PRECISION :: rr, w, vbar
INTEGER :: i 

vbar = 0 
do i=1, nn(1)
    vbar = vbar + v_p( i )%x                  
end do ! all particle
vbar = 1.0d0*vbar/nn(1)
mobility = vbar/ff/eta

do i = 1, N_total-nt
    rr = dsqrt( v_p( i )%x*v_p( i )%x + v_p( i )%y*v_p( i )%y )
    v_x(i) = v_p( i )%x/rr
    v_y(i) = v_p( i )%y/rr                    
end do ! all particle

w = 0
do i = 1, nn(1)
    w = w + v_x(i)*charge(i)/( abs(charge(i)) )                    
end do ! all particle
w = w/nn(1)

End SUBROUTINE order_w


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE statis (fun,n_x,funbar,funvar)
IMPLICIT NONE
INTEGER :: n_x 
DOUBLE PRECISION :: fun(1:n_x)
DOUBLE PRECISION :: funbar,funvar
INTEGER :: i 
funbar = 0
do i=1, n_x
    funbar = funbar + fun(i)
end do
funbar = funbar/n_x

funvar = 0
do i=1, n_x
    funvar = funvar + (fun(i) - funbar)**2
end do
funvar = funvar/n_x

End SUBROUTINE statis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_position ()
USE global_parameters
IMPLICIT NONE

INTEGER :: i, j, ix, iy
 
DOUBLE PRECISION :: rr, deltax, deltay

do i=1, na
    particle( i )%x = particle( i )%x + v_p( i )%x 
    particle( i )%y = particle( i )%y + v_p( i )%y
    if (particle(i)%x>Lx ) then
        particle(i)%x = particle(i)%x - Lx
    else if (particle(i)%x<0 ) then
        particle(i)%x = particle(i)%x + Lx
    end if 
    if (particle(i)%y>Ly ) then
        particle(i)%y = particle(i)%y - Ly
    else if (particle(i)%y<0 ) then
        particle(i)%y = particle(i)%y + Ly
    end if 
!     print*,"p",i, particle( i )%x,particle( i )%y                  
end do ! all particle

End SUBROUTINE update_position

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! should be used right after calling update 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE position_corelation ()
USE global_parameters
IMPLICIT NONE

INTEGER :: i, j, m, n
DOUBLE PRECISION :: phibar

phi_y = 0
do i=1, nn(1)
    phi_y( i_y(i) ) = phi_y( i_y(i) ) + 1                  
end do ! all particle

phibar = 0
do i = 1, Ny
    phibar = phibar + phi_y(i)
end do 
phibar = phibar/ny

phiyy = 0
do m = 1, ny
    do n = 1, ny
        phiyy( abs(m-n) ) = phiyy( abs(m-n) ) + phi_y(m)*phi_y(n)
    end do     
end do
phiyy = phiyy/ny/ny - phibar*phibar

End SUBROUTINE position_corelation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE velocity ()
USE global_parameters
IMPLICIT NONE

INTEGER :: i 
 
do i=1, na
    v_p( i )%x = ( force_p(i)%x ) * dt_f + gasdev(seed)*dt_r 
    v_p( i )%y = ( force_p(i)%y + ff*charge(i) ) * dt_f + gasdev(seed)*dt_r
!     if (force_p(i)%x /= 0 ) then
!         print*,"v1",force_p(i)%x, force_p(i)%y, gasdev(seed)*dt_r 
!     end if 
!     print*,"v",i, v_p( i )%x,v_p( i )%y                 
end do ! all particle

End SUBROUTINE velocity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE F_poten ()
USE global_parameters
IMPLICIT NONE

INTEGER :: i, j, ix, iy
INTEGER :: ii, jj,iprime
DOUBLE PRECISION :: rr, deltax, deltay
DOUBLE PRECISION :: force_temp
 
do i=1, na
    force_p(i)%x = 0
    force_p(i)%y = 0
    do ii = -1, 1
        do jj = -1,1
            ix = i_x(i) + ii
            iy = i_y(i) + jj 
            if ( ii == 0 .and. jj==0 ) then
                if ( pointor(ix,iy)>1) then
                    do j = 1, pointor(ix,iy)                          
                        iprime = cell( j, ix, iy )
                        if( i /= iprime )then                                                     
                            deltax = particle( i )%x - particle( iprime )%x
                            deltay = particle( i )%y - particle( iprime )%y
                            rr = dsqrt( deltax*deltax + deltay*deltay ) 
                            if (iprime<=na) then
                                if( rr<=sigma ) then
! wca interaction
                                    force_temp = kWCA*(2.0d0/rr**13-1.0d0/rr**7)/rr
                                    if ( force_temp.ne.force_temp ) then
                                       force_temp = 0.0d0 
                                    end if                             
                                    if ( isnan(force_temp) .or. force_temp>100000 ) then 
                                        print*, force_temp, rr
                                        print*, particle( i )%x, particle( iprime )%x
                                        print*, particle( i )%y, particle( iprime )%y
                                        print*, deltax, deltay
                                        print*, i, iprime
                                        print*, ix,iy
                                        print*, i_x(i),i_y(i)                                        
                                      stop "force_temp1"
                                    end if
                                    force_p(i)%x = force_p(i)%x + force_temp*deltax
                                    force_p(i)%y = force_p(i)%y + force_temp*deltay
!                                     print*,"f1",i, iprime,force_p(i)%x, force_p(i)%y
                                end if        
                            else if (iprime>=na+nb+1) then
                                if( rr<=sigma*0.5d0+r_c ) then
! wca interaction
                                    force_temp = kWCA*(2.0d0/rr**13-1.0d0/rr**7)/rr
                                    if ( force_temp.ne.force_temp ) then
                                       force_temp = 0.0d0 
                                    end if                             
                                    if ( isnan(force_temp) .or. force_temp>100000 ) then
                                        print*, force_temp, rr
                                        print*, particle( i )%x, particle( iprime )%x
                                        print*, particle( i )%y, particle( iprime )%y
                                        print*, deltax, deltay
                                        print*, i, iprime
                                        print*, ix,iy
                                        print*, i_x(i),i_y(i)                                     
                                        stop "force_temp2"
                                    end if
                                    force_p(i)%x = force_p(i)%x + force_temp*deltax
                                    force_p(i)%y = force_p(i)%y + force_temp*deltay
!                                     print*,"f2",i, iprime,force_p(i)%x, force_p(i)%y
                                end if
                            else 
                                if( rr<=sigma*0.5d0+r_b ) then 
                                    if(particle( i )%x <= Lx*0.5d0) then
                                        force_p(i)%x = -f_b
                                        force_p(i)%y = 0
                                    else
                                        force_p(i)%x = f_b
                                        force_p(i)%y = 0
                                    end if
                                end if
                            end if
                        end if
                    end do                        
                end if                    
            else !if ( ii /= 0 .and. jj /=0 )

                if ( pointor(ix,iy)>0 ) then
!                    print*, ix,iy,"ix,iy"
                    do j = 1, pointor(ix,iy)  
                        iprime = cell( j, ix, iy )                                                    
                        deltax = particle( i )%x - ( particle( iprime )%x + shift_x(ix,iy) )
                        deltay = particle( i )%y - ( particle( iprime )%y + shift_y(ix,iy) )
                        rr = dsqrt( deltax*deltax + deltay*deltay )  
                        if (iprime<=na) then
                            if( rr<=sigma ) then
! wca interaction
                                force_temp = kWCA*(2.0d0/rr**13-1.0d0/rr**7)/rr
                                if ( force_temp.ne.force_temp ) then
                                   force_temp = 0.0d0 
                                end if                             
                                if ( isnan(force_temp) .or. force_temp>100000 ) then
                                    print*, force_temp, rr
                                    print*, particle( i )%x, particle( iprime )%x
                                    print*, particle( i )%y, particle( iprime )%y
                                    print*, deltax, deltay
                                    print*, i, iprime
                                    print*, ix,iy
                                    print*, i_x(i),i_y(i)                                      
                                    stop "force_temp3"
                                end if
                                force_p(i)%x = force_p(i)%x + force_temp*deltax
                                force_p(i)%y = force_p(i)%y + force_temp*deltay
!                                 print*,"f3",i, iprime,force_p(i)%x, force_p(i)%y
                            end if        
                        else if (iprime>=na+nb+1) then
                            if( rr<=sigma*0.5d0+r_c ) then
! wca interaction
                                force_temp = kWCA*(2.0d0/rr**13-1.0d0/rr**7)/rr
                                if ( force_temp.ne.force_temp ) then
                                   force_temp = 0.0d0 
                                end if                             
                                if ( isnan(force_temp) .or. force_temp>100000 ) then
                                    print*, force_temp, rr
                                    print*, particle( i )%x, particle( iprime )%x
                                    print*, particle( i )%y, particle( iprime )%y
                                    print*, deltax, deltay
                                    print*, i, iprime
                                    print*, ix,iy
                                    print*, i_x(i),i_y(i)                                     
                                    stop "force_temp4"
                                end if
                                force_p(i)%x = force_p(i)%x + force_temp*deltax
                                force_p(i)%y = force_p(i)%y + force_temp*deltay
!                                 print*,"f4",i, iprime,force_p(i)%x, force_p(i)%y
                            end if
                        else 
                            if( rr<=sigma*0.5d0+r_b ) then 
                                if(particle( i )%x <= Lx*0.5d0) then
                                    force_p(i)%x = -f_b
                                    force_p(i)%y = 0
                                else
                                    force_p(i)%x = f_b
                                    force_p(i)%y = 0
                                end if
                            end if
                        end if
                    end do    
                end if ! no particle in cell( ix,iy ) 
                           
            end if ! ( ii /= 0 .and. jj /=0 ) 
                           
        end do  ! jj
    end do  ! ii

end do ! all particle

End SUBROUTINE F_poten

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE update_cell ()
USE global_parameters
IMPLICIT NONE

INTEGER :: i, j
!open(24,file='ixiy.dat')
cell = 0
pointor = 0 
i_x = 0
i_y = 0

do i=1, N_total
        i_x(i) = floor( particle(i)%x / dd ) + 1
        i_y(i) = floor( particle(i)%y / dd ) + 1
        pointor (i_x(i),i_y(i)) = pointor (i_x(i),i_y(i)) + 1        
        cell(pointor(i_x(i),i_y(i)),i_x(i),i_y(i)) = i
        
!        write(24,"(2E20.8,2I5)") particle(i)%x,particle(i)%y,i_x(i),i_y(i)
end do
!close(24)

!!initialize the ghoust cell
do i = 1, nx
    do j = 1, Nmax
        cell(j, i, 0) = cell(j, i, ny)
        cell(j, i, ny+1) = cell(j, i, 1)   
    end do
    pointor(i, 0) = pointor(i, ny)
    pointor(i, ny+1) = pointor(i, 1)
end do 
do i = 1, ny
    do j = 1, Nmax
        cell(j, 0, i) = cell(j, nx, i)
        cell(j, nx+1, i) = cell(j, 1, i)   
    end do
    pointor(0, i) = pointor(nx, i)
    pointor(nx+1, i) = pointor(1, i)
end do 

do j= 1, Nmax
    cell(j, 0, 0) = cell(j, nx, ny)
    cell(j, 0, ny+1) = cell(j, nx, 1)
    cell(j, nx+1, 0) = cell(j, 1, ny)
    cell(j, nx+1, ny+1) = cell(j, 1, 1)    
end do
pointor(0, 0) = pointor(nx, ny)
pointor(0, ny+1) = pointor(nx, 1)
pointor(nx+1, 0) = pointor(1, ny)
pointor(nx+1, ny+1) = pointor(1, 1) 

End SUBROUTINE update_cell
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE dump_force( tt )
USE global_parameters
IMPLICIT NONE

INTEGER :: i, tt 
character*9 res
character res0,res1,res2,res3
character(len=50) :: filename

res0 = achar(48+mod(tt,10))
res1 = achar(48+mod(int(tt/10),10))
res2 = achar(48+mod(int(tt/100),10))
res3 = achar(48+int(tt/1000))

write(filename, '(A,I1,A)') 'Result-', layout_type, '/ProcessLocation/'
res = trim(filename)//'f'// res3 // res2 // res1 // res0 // '.dat'
! res = 'f'// res3 // res2 // res1 // res0 // '.dat'
open(68,file=res)
do i=1, N_total
    write(68,"(4E20.8)") particle( i )%x, particle( i )%y, particle( i )%x+force_p( i )%x, particle( i )%y+force_p(i)%y                   
end do ! all particle
close(68)

End SUBROUTINE dump_force

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE dump_position( tt )
USE global_parameters
IMPLICIT NONE

INTEGER :: i, tt 
character*40 resl,resr,resi,rest
character res0,res1,res2,res3
character(len=50) :: filename

res0 = achar(48+mod(tt,10))
res1 = achar(48+mod(int(tt/10),10))
res2 = achar(48+mod(int(tt/100),10))
res3 = achar(48+int(tt/1000))

! resl = 'p1'// res3 // res2 // res1 // res0 // '.dat'
! resr = 'p2'// res3 // res2 // res1 // res0 // '.dat'
! resi = 'p3'// res3 // res2 // res1 // res0 // '.dat'
write(filename, '(A,I1,A)') 'Result-', layout_type, '/ProcessLocation/'
rest = trim(filename)//'p'// res3 // res2 // res1 // res0 // '.dat'
! rest = 'Result/ProcessLocation/p'// res3 // res2 // res1 // res0 // '.dat'
! open(68,file=resl)
! open(69,file=resr)
! open(70,file=resi)
open(71,file=rest)
! do i=1,na
!     write(68,*) particle(i)%x, particle(i)%y
! end do
! do i=na+1,na+nb
!     write(69,*) particle(i)%x, particle(i)%y
! end do
! do i=na+nb+1,N_total
!     write(70,*) particle(i)%x, particle(i)%y
! end do
do i=1,N_total
    write(71,*) particle(i)%x, particle(i)%y
end do
! close(69)
! close(68)
! close(70)
close(71)


End SUBROUTINE dump_position

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE dump_velocity( tt )
USE global_parameters
IMPLICIT NONE

INTEGER :: i, tt 
character*40 resl,resr,resi
character res0,res1,res2,res3
character(len=50) :: filename

res0 = achar(48+mod(tt,10))
res1 = achar(48+mod(int(tt/10),10))
res2 = achar(48+mod(int(tt/100),10))
res3 = achar(48+int(tt/1000))

write(filename, '(A,I1,A)') 'Result-', layout_type, '/ProcessLocation/'
resl = trim(filename)//'v1'// res3 // res2 // res1 // res0 // '.dat'

! resr = 'v2'// res3 // res2 // res1 // res0 // '.dat'
! resi = 'v3'// res3 // res2 // res1 // res0 // '.dat'
open(68,file=resl)
! open(69,file=resr)
! open(70,file=resi)
do i=1,na
    write(68,*) v_p( i )%x , v_p( i )%y
end do
! do i=na+1,na+nb
!     write(69,*) v_p( i )%x , v_p( i )%y
! end do
! do i=na+nb+1,N_total
!     write(70,*) v_p( i )%x , v_p( i )%y
! end do


close(68)
! close(69)
! close(70)


End SUBROUTINE dump_velocity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION gasdev(idum)
INTEGER idum
REAL gasdev

!Returns a normally distributed deviate with zero mean and unit variance, using ran1(idum)
!as the source of uniform deviates.
INTEGER iset
REAL fac,gset,rsq,v1,v2,ran1
SAVE iset,gset
DATA iset/0/
if (idum.lt.0) iset=0 !Reinitialize.
if (iset.eq.0) then !We don��t have an extra deviate handy, so
    rsq = 10.
    do while(rsq.ge.1..or.rsq.eq.0.)
        v1=2.*ran2(idum)-1. !pick two uniform numbers in the square extending from -1 to +1 in each direction, 
        v2=2.*ran2(idum)-1.
        rsq=v1**2+v2**2 !see if they are in the unit circle,
       
    end do
    fac=sqrt(-2.*log(rsq)/rsq) !Now make the Box-Muller transformation to get two normal deviates. Return one and save the other for next time.
    gset=v1*fac
    gasdev=v2*fac
    iset=1 !Set flag.
    
else ! We have an extra deviate handy,
    gasdev=gset !so return it,
    iset=0 !and unset the flag.
endif
return
END function

!!!!!!!!!!This is the function that creat a random number of [0-1]!!!!!!!!!!!!!!
      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
        IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,	 &
        NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
        end do
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=AM*iy
      if(ran2.gt.RNMX) ran2=RNMX
      return
      END function  



END MODULE utility_routines
