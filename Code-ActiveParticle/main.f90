program MC
USE global_parameters         ! 使用全局参数模块
USE utility_routines          ! 使用实用工具模块
IMPLICIT none                 ! 禁用隐式类型

Integer :: i, j, m, ic, ii, MCS, tt   ! 定义整数变量
Integer :: Nmove, move, num,n_vo      ! 定义移动相关的整数变量
integer :: ddt(8)                      ! 定义整数数组，存储日期和时间信息
Integer, allocatable :: nstart(:),nend(:)    ! 定义可分配整数数组，用于存储起始和结束索引
Integer, allocatable :: delta_t(:)     ! 定义时间步的可分配整数数组
DOUBLE PRECISION :: rho, charge3, distance2   ! 定义双精度变量
DOUBLE PRECISION :: dx, number, ftemp         ! 定义双精度变量
DOUBLE PRECISION :: a_plane,a_plane_2,a_plane_3  ! 定义平面相关的双精度变量
DOUBLE PRECISION :: rx, ry, rr                 ! 定义轨迹相关的双精度变量
DOUBLE PRECISION :: wbar, wvar, vbar, vvar     ! 定义统计变量
DOUBLE PRECISION, allocatable :: pyy(:), ww(:),vv(:), distance(:)  ! 定义可分配的双精度数组

type(node),allocatable :: r_end(:)   ! 定义节点类型的可分配数组
character*9 resmsd                   ! 定义字符变量
character*7 res                      ! 定义字符变量
character res0,res1,res2,res3        ! 定义字符变量
character(len=50) :: filename

open(unit=10,file='input.txt')       ! 打开输入文件
read(10,*) Nmax                      ! 读取最大粒子数
read(10,*) Ntype                     ! 读取粒子类型数
read(10,*) na                        ! 读取第一种粒子的数量
read(10,*) nb                        ! 读取第二种粒子的数量
read(10,*) nc                        ! 读取第三种粒子的数量
read(10,*) a_plane_2                 ! 读取底部粒子B的参数
read(10,*) a_plane_3                 ! 读取墙壁粒子C的参数
read(10,*) r_b                       ! 读取第二种类型粒子的半径
read(10,*) r_c                       ! 读取第三种类型粒子的半径
read(10,*) deltat                    ! 读取时间步长
read(10,*) dd                        ! 读取单元格大小
read(10,*) eta                       ! 读取粘度
read(10,*) kappa                     ! 读取粒子的刚度
read(10,*) Lx                        ! 读取x方向长度
read(10,*) Ly                        ! 读取y方向长度
read(10,*) Npre                      ! 读取预旋转数
read(10,*) ftemp                     ! 读取场强
read(10,*) f_b                       ! 读取触碰墙时的力
read(10,*) Nmove                     ! 读取一步中的移动数
read(10,*) ii                        ! 读取NMCs的幂
read(10,*) kWCA                      ! 读取WCA参数
read(10,*) charge3                   ! 读取第三种类型粒子的电荷
read(10,*) layout_type               ! 读取晶格初始化模式
close(10)                            ! 关闭文件

! 初始化参数
open(unit=30,file='delta.txt')       ! 打开delta文件
read(30,*) ic                        ! 读取ic值
close(30)                            ! 关闭文件
call date_and_time(values=ddt)       ! 调用日期和时间函数
seed = (ddt(8)+10*ic)                ! 设置随机种子

ftemp = ftemp * 2**ic                ! 调整场强
kWCA = kWCA*24.0d0                   ! 根据f调整WCA参数
dt_f = deltat/eta                    ! 计算流体时间步长
dt_r = dsqrt(2.0d0*deltat/eta)       ! 计算随机时间步长

print*, dt_f,"dt_f"                  ! 输出流体时间步长
print*, dt_r,"dt_r"                  ! 输出随机时间步长

sigma = 2.0d0**(1.0d0/6)             ! 计算sigma值
dd = 3.0d0*sigma                     ! 计算单元格尺寸
print*, sigma,"sigma"                ! 输出sigma
print*, dd,"dd"                      ! 输出单元格尺寸

allocate( nn(1:ntype),nstart(1:ntype),delta_t(1:n_vo) )   ! 分配数组
allocate( distance(1:ntype),Nend(1:ntype) )              ! 分配距离和结束数组

delta_t(1) = 10                                           ! 设置时间增量
delta_t(2) = 100
delta_t(3) = 1000
delta_t(4) = 10000

nn(1) = na                                                ! 设置第一种粒子的数量
nstart(1) = 1                                             ! 设置第一种粒子的起始位置
nend(1) = nn(1)                                           ! 设置第一种粒子的结束位置

nn(2) = nb+nc                                             ! 设置第二和第三种粒子的总数
nstart(2) = nn(1)+1                                       ! 设置第二种粒子的起始位置
nend(2) = nn(1)+nn(2)                                     ! 设置第二种粒子的结束位置

if(ntype==3)then
    nn(3) = na                                            ! 如果粒子类型为3，则设置第三种粒子的数量
    nstart(3) = nn(2)+1                                   ! 设置第三种粒子的起始位置
    nend(3) = nn(2) + nn(3)                               ! 设置第三种粒子的结束位置
    print*,nstart(3),"nstart3"                            ! 输出第三种粒子的起始位置
end if

N_total = na+nb+nc                                        ! 计算总粒子数
rho = N_total/Lx/Ly                                       ! 计算密度

nx = floor(0.999999*Lx/dd)+1                              ! 计算x方向的单元格数
ny = floor(0.999999*Ly/dd)+1                              ! 计算y方向的单元格数
print*,"nxny",nx,ny                                       ! 输出单元格数

NMCs = 2*10**ii                                           ! 计算总的MC步数

allocate( particle(1:N_total), i_x(1:N_total), i_y(1:N_total) )  ! 分配粒子和坐标数组
allocate( r_end(1:N_total),v_p(1:N_total) )               ! 分配末端向量和速度数组
allocate( charge(1:N_total), force_p(1:N_total) )         ! 分配电荷和力数组
allocate( pointor(0:nx+1,0:ny+1), cell(1:Nmax,0:nx+1,0:ny+1) )   ! 分配指针和单元格数组
allocate( shift_x(0:nx+1,0:ny+1), shift_y(0:nx+1,0:ny+1)) ! 分配位移数组
allocate( phiyy(0:ny+1),pyy(0:ny+1), phi_y(0:ny+1))       ! 分配势能数组
allocate( ww(1:NMCs),vv(1:NMCs) )                         ! 分配统计数组

print*,"allocated"                                        ! 输出已分配提示

dx = dsqrt( lx*ly/N_total )                                ! 计算x方向位移

shift_x = 0                                               ! 初始化x方向位移
shift_y = 0                                               ! 初始化y方向位移

!! 初始化边界单元格
do i = 1, nx
    shift_x(i, 0) = 0
    shift_y(i, 0) = -ly
    shift_x(i, ny+1) = 0
    shift_y(i, ny+1) = ly
end do 
do i = 1, ny
    shift_x(0, i) = -lx
    shift_y(0, i) = 0
    shift_x(nx+1, i) = lx
    shift_y(nx+1, i) = 0
end do 

shift_x(0, 0) = -lx
shift_y(0, 0) = -ly

shift_x(0, ny+1) = -lx
shift_y(0, ny+1) = ly

shift_x(nx+1, 0) = lx
shift_y(nx+1, 0) = -ly

shift_x(nx+1, ny+1) = lx
shift_y(nx+1, ny+1) = ly

do i=1, nn(1)
    charge(i) = 1                                             ! 设置第一种粒子的电荷
end do
do i=1, nn(2)
    charge(nn(1)+i) = 0                                       ! 设置第二种粒子的电荷
end do
if (ntype==3)then
    do i=1, nn(3)
        charge(nn(2)+i) = charge3                             ! 设置第三种粒子的电荷
    end do
end if
print*,charge(nstart(1)),"charge1"                            ! 输出第一种粒子的电荷
print*,charge(nstart(2)),"charge2"                            ! 输出第二种粒子的电荷
print*,charge(N_total),"charge3"                              ! 输出最后一种粒子的电荷

print*,"initializing"                                         ! 输出初始化提示

! 初始化粒子位置
a_plane = 1.01d0*(1.0d0*Lx*Ly/(N_total+1))**(1.0d0/2.0d0)     ! 计算平面参数
print*, a_plane,"1"                                           ! 输出平面参数
particle(1)%x = a_plane*0.5d0                                 ! 初始化第一个粒子的x坐标
particle(1)%y = Ly-a_plane                                    ! 初始化第一个粒子的y坐标
do i=2, na
    particle(i)%x = particle(i-1)%x + a_plane                 ! 初始化第一种粒子的x坐标
    particle(i)%y = particle(i-1)%y                           ! 初始化第一种粒子的y坐标
    if (particle(i)%x >= Lx*0.5d0 - a_plane_3 - a_plane .and. particle(i)%x <= Lx*0.5d0 + a_plane_3 ) then
        particle(i)%x = Lx*0.5d0 + a_plane_3*2.0d0 +a_plane   ! 调整第一种粒子的x坐标
    end if
    if (particle(i)%x >= (Lx-sigma)) then
        particle(i)%x = a_plane*0.5d0                         ! 重置第一种粒子的x坐标
        particle(i)%y = particle(i)%y - a_plane               ! 调整第一种粒子的y坐标
    end if
end do

particle(na+1)%x = Lx*0.5d0 - a_plane_3*0.01d0                       ! 初始化第二种粒子的x坐标
particle(na+1)%y = 0                                          ! 初始化第二种粒子的y坐标
do i=na+2,na+nb
    particle(i)%x = particle(i-1)%x                           ! 初始化第二种粒子的x坐标
    particle(i)%y = particle(i-1)%y + a_plane_3               ! 初始化第二种粒子的y坐标
    if (particle(i)%y >= Ly + a_plane_3) then
        particle(i)%x = Lx*0.5d0 + a_plane_3*0.01d0                  ! 调整第二种粒子的x坐标
        particle(i)%y = 0.0d0                                 ! 重置第二种粒子的y坐标
    end if 
end do  


! particle(na+nb+1)%x = a_plane_2*0.5d0                         ! 初始化第三种粒子的x坐标
! particle(na+nb+1)%y = a_plane_2                               ! 初始化第三种粒子的y坐标
! do i=na+nb+2, N_total
!     particle(i)%x = particle(i-1)%x + a_plane_2               ! 初始化第三种粒子的x坐标
!     particle(i)%y = particle(i-1)%y                           ! 初始化第三种粒子的y坐标
!     if (particle(i)%x >= Lx*0.5d0 - a_plane_3-a_plane_2 .and. particle(i)%x <= Lx*0.5d0 + a_plane_3 ) then
!         particle(i)%x = Lx*0.5d0 + a_plane_3*2.0d0 +a_plane_2 ! 调整第三种粒子的x坐标
!     end if
!     if (particle(i)%x >= (Lx-sigma)) then
!         particle(i)%x = a_plane_2*0.5d0                       ! 重置第三种粒子的x坐标
!         particle(i)%y = particle(i)%y + a_plane_2             ! 调整第三种粒子的y坐标
!     end if
! end do 


! layout_type = 1
CALL initialize_C_particles(a_plane_2, layout_type)

write(filename, '(A,I1,A)') 'Result-', layout_type, '/InitLocation/A.txt'
open(20,file=filename) 
do i=1, na
    write(20,*)  particle(i)%x, particle(i)%y                ! 将第一种粒子的位置写入A.txt
end do
close (20) 
write(filename, '(A,I1,A)') 'Result-', layout_type, '/InitLocation/B.txt'
open(21,file=filename) 
do i=na+1,na+nb
    write(21,*)  particle(i)%x, particle(i)%y                ! 将第二种粒子的位置写入B.txt
end do
close (21)
write(filename, '(A,I1,A)') 'Result-', layout_type, '/InitLocation/C.txt'
open(22,file=filename) 
do i=na+nb+1, N_total
    write(22,*)  particle(i)%x, particle(i)%y                ! 将第三种粒子的位置写入C.txt
end do
close (22)


call update_cell ()                                           ! 更新单元格信息
call dump_position( 0 )                                       ! 输出位置信息
! 预跑
ff = ftemp                                                    ! 设置场强
MCS = 0                                                       ! 初始化MC步数
do while(MCS < Npre)

    call F_poten ()                                           ! 计算势能
    call velocity ()                                          ! 计算速度
    call update_position ()                                   ! 更新位置
    call update_cell ()                                       ! 更新单元格  
    Mcs = Mcs + 1                                             ! 更新MC步数
    
end do

print*, "initial is finished"                                 ! 输出初始化完成提示

! MD模拟
MCS = 0                                                       ! 重置MC步数
ff = ftemp                                                    ! 重新设置场强
tt = 0                                                        ! 初始化时间
do while(MCS < NMCs)
    Mcs = Mcs + 1    
    move = 0                                                  ! 初始化移动次数
! 
    do while ( move<Nmove ) 
        move = move + 1                                       ! 更新移动次数
        call F_poten ()                                       ! 计算势能
        call velocity ()                                      ! 计算速度
        call update_position ()                               ! 更新位置
        call update_cell ()                                   ! 更新单元格
    end do

    if ( mod(MCS,2000)==0 ) then
        tt = tt + 1                                           ! 更新时间
        print*,tt,"tt"
        call dump_position(tt)                                ! 输出位置
        call dump_velocity(tt)                                ! 输出速度      
    end if     
end do

print*, "steady state is reached"                             ! 输出达到稳态提示

write(filename, '(A,I1,A)') 'Result-', layout_type, '/FinalLocation/A.txt'
open(20,file=filename) 
do i=1, na
    write(20,*)  particle(i)%x, particle(i)%y                ! 将最终的第一种粒子位置写入Afinal.txt
end do
close (20) 
write(filename, '(A,I1,A)') 'Result-', layout_type, '/FinalLocation/B.txt'
open(21,file=filename) 
do i=na+1,na+nb
    write(21,*)  particle(i)%x, particle(i)%y                ! 将最终的第二种粒子位置写入Bfinal.txt
end do
close (21)
write(filename, '(A,I1,A)') 'Result-', layout_type, '/FinalLocation/C.txt'
open(22,file=filename) 
do i=na+nb+1, N_total
    write(22,*)  particle(i)%x, particle(i)%y                ! 将最终的第三种粒子位置写入Cfinal.txt
end do
close (22)
! 模拟结束

print*, "sim finished"                                       ! 输出模拟完成提示

End program MC                                               ! 结束程序
