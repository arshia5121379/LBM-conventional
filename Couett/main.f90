
program couett
    implicit none
    integer , parameter :: n=101,m=41
    doubleprecision :: f(0:8,1:n,1:m),feq(0:8,1:n,1:m)
    doubleprecision :: rho(1:n,1:m),rho_old(1:n,1:m)
    doubleprecision :: u(1:n,1:m),v(1:n,1:m),u_old(1:n,1:m),v_old(1:n,1:m)
    real :: X(1:n,1:m),Y(1:n,1:m),tminv(0:8,0:8),sm(0:8),tm(0:8,0:8),stmiv(0:8,0:8)
    real :: w(0:8),cx(0:8),cy(0:8),ev(0:8,0:8)
    real :: nu , dx , dy , dt , uo , Re , omega , rhoo , eps = 1e-7 &
            & , err_rho , err_u , err_v , a1 , sumcc,tau,PI ,Cs
    integer :: i,j,k, kk , mstep

    PI = 4*atan(1.d0)

    err_rho = 0.0
    err_u   = 0.0
    err_v   = 0.0
    mstep=10000

    !   6 2 5
    !    \|/
    !   3-0-1
    !    /|\
    !   7 4 8

    w(0)=4./9.
    do i=1,4
        w(i)=1./9.
    end do
    do i=5,8
        w(i)=1./36.
    end do

    ! compute the direction of lattice velocities
    cx(:) = (/0.0,1.0,0.0,-1.0,0.0,1.0,-1.0,-1.0,1.0/)
    cy(:) = (/0.0,0.0,1.0,0.0,-1.0,1.0,1.0,-1.0,-1.0/)

    !cx(0) = 0.0
    !cy(0) = 0.0
    !do i=1,4
     !   cx(i) = cos(PI*(i-1.0)/(2.0))
     !   cy(i) = sin(PI*(i-1.0)/(2.0))
    !end do

     !do i=5,8
     !   cx(i) = sqrt(2.0)*cos(PI*(2.0*i-9.0)/(4.0))
     !   cy(i) = sqrt(2.0)*sin(PI*(2.0*i-9.0)/(4.0))
    !end do


    CALL init(n,m,u,v,rho,omega,tau,dx,dy,uo,nu)
    CALL mesh(n,m,X,Y,dx,dy)
    CALL MRT(n,m,tminv,sm,tm,stmiv,tau,nu)

        open(unit=2,file="unsteady.dat")
        write(2, '(A)') 'VARIABLES = "U" , "Y"'

! main loop
    do kk=1,mstep

        rho_old = rho
        u_old   = u
        v_old   = v

        call collesion(u,v,f,feq,rho,omega,w,cx,cy,n,m,tm,tminv,stmiv)
        call streaming(f,n,m)
        ! ——————————–
        call boundary(f,n,m,uo)
        call macro(f,rho,u,v,cx,cy,n,m,uo)

        ! convergence checking

        err_rho = maxval(abs(rho - rho_old))
        err_u   = maxval(abs(u   - u_old))
        err_v   = maxval(abs(v   - v_old))


        if (mod(kk,10)==0) then
            print * , "========================================================="
            print *, "the rho residual : " , err_rho
            print *, "the u residual   : " , err_u
            print *, "the v residual   : " , err_v
            print *, "time step        : " , kk
            print * , "========================================================="
            CALL unsteady(m,u,Y,uo,kk)
        end if

        if (err_rho .LE. eps .AND. err_u .LE. eps .AND. err_v .LE. eps) then
            print * , "the solution has been converged at iteration : " , kk
            exit
        end if
    END DO

! end of the main loop

    CALL output(n,m,dx,dy,rho,u,v,X,Y)

    CALL section(n,m,u,X,Y,uo)

    contains

! end of the main program
    subroutine init(n,m,u,v,rho,omega,tau,dx,dy,uo,nu)
        integer,INTENT(IN)::n,m
        doubleprecision,INTENT(INOUT)::u(1:n,1:m),v(1:n,1:m),rho(1:n,1:m)
        real,INTENT(INOUT):: omega,tau,dx,dy,uo,nu
        real:: rhoo,dt,Re,H,L
        integer::i,j

        H = 0.04
        L = 0.08
        uo=0.01
        rhoo=1.00
        dx=1.0
        dy=dx
        dt=dx
        Re = 5000.0
        nu =uo*m/Re


        print *, "Re = ", Re , "nu = ", nu
        omega=1.0/(3.*nu+0.5)
        tau  = 3.*nu/dt+0.5

        do j=1,m
            do i=1,n
                rho(i,j)=rhoo
                u(i,j)  =0.0
                v(i,j)  =0.0
            end do
        end do

        do i=1,n
            u(i,m)=uo
            v(i,m)=0.0
        end do

    end subroutine


subroutine MRT(n,m,tminv,sm,tm,stmiv,tau,nu)
        integer,INTENT(IN) ::n,m
        real,INTENT(INOUT) :: tminv(0:8,0:8),sm(0:8),tm(0:8,0:8),stmiv(0:8,0:8),tau,nu
        real::s3,s5,s7,s8,dt=1.0
        integer:: i,j,k,l


            tm(0,:) = (/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/)
            tm(1,:) = (/-4.0,-1.0,-1.0,-1.0,-1.0,2.0,2.0,2.0,2.0/)
            tm(2,:) = (/4.0,-2.0,-2.0,-2.0,-2.0,1.0,1.0,1.0,1.0/)
            tm(3,:) = (/0.0,1.0,0.0,-1.0,0.0,1.0,-1.0,-1.0,1.0/)
            tm(4,:) = (/0.0,-2.0,0.0,2.0,0.0,1.0,-1.0,-1.0,1.0/)
            tm(5,:) = (/0.0,0.0,1.0,0.0,-1.0,1.0,1.0,-1.0,-1.0/)
            tm(6,:) = (/0.0,0.0,-2.0,0.0,2.0,1.0,1.0,-1.0,-1.0/)
            tm(7,:) = (/0.0,1.0,-1.0,1.0,-1.0,0.0,0.0,0.0,0.0/)
            tm(8,:) = (/0.0,0.0,0.0,0.0,0.0,1.0,-1.0,1.0,-1.0/)

            a1 = 1.0/36.0

            tminv(0,:) = (/4.0*a1,-4.0*a1,4.0*a1,0.0,0.0,0.0,0.0,0.0,0.0/)
            tminv(1,:) = (/4.0*a1,-a1,-2.0*a1,6.0*a1,-6.0*a1,0.0,0.0,9.0*a1,0.0/)
            tminv(2,:) = (/4.0*a1,-a1,-2.0*a1,0.0,0.0,6.0*a1,-6.0*a1,-9.0*a1,0.0/)
            tminv(3,:) = (/4.0*a1,-a1,-2.0*a1,-6.0*a1,6.0*a1,0.0,0.0,9.0*a1,0.0/)
            tminv(4,:) = (/4.0*a1,-a1,-2.0*a1,0.0,0.0,-6.0*a1,6.0*a1,-9.0*a1,0.0/)
            tminv(5,:) = (/4.0*a1,2.0*a1,a1,6.0*a1,3.0*a1,6.0*a1,3.0*a1,0.0,9.0*a1/)
            tminv(6,:) = (/4.0*a1,2.0*a1,a1,-6.0*a1,-3.0*a1,6.0*a1,3.0*a1,0.0,-9.0*a1/)
            tminv(7,:) = (/4.0*a1,2.0*a1,a1,-6.0*a1,-3.0*a1,-6.0*a1,-3.0*a1,0.0,9.0*a1/)
            tminv(8,:) = (/4.0*a1,2.0*a1,a1,6.0*a1,3.0*a1,-6.0*a1,-3.0*a1,0.0,-9.0*a1/)

            Cs = 1.0/sqrt(3.0)

            s3 = 1.0
            s5 = s3
            s7 = 1.0/tau
            s8 = s7

            sm(:) = (/0.5,0.8,0.8,s3,1.2,s5,1.2,s7,s8/)

            do i=0,8
                do j=0,8
                    stmiv(i,j) = tminv(i,j)*sm(j)
                end do
            end do


            do i=0,8
                do j=0,8
                    sumcc = 0.0
                    do k=0,8
                        sumcc = sumcc+tminv(i,k)*tm(k,j)
                    end do
                    ev(i,j) = sumcc
                end do
            end do

    end subroutine

    subroutine collesion(u,v,f,feq,rho,omega,w,cx,cy,n,m,tm,tminv,stmiv)
            integer , INTENT(IN) :: n , m
            doubleprecision , INTENT(INOUT) :: f(0:8,1:n,1:m),feq(0:8,1:n,1:m),rho(1:n,1:m),u(1:n,1:m),v(1:n,1:m)
            doubleprecision :: fmom(0:8,1:n,1:m),fmeq(0:8,1:n,1:m)
            real , INTENT(IN) :: w(0:8),cx(0:8),cy(0:8),omega,tminv(0:8,0:8),tm(0:8,0:8),stmiv(0:8,0:8)
            real :: suma,sumb
            integer :: i,j,k,l

            ! calculate the feq moments:
            do i=1,n
                do j=1,m
                    fmeq(0,i,j) = rho(i,j)
                    fmeq(1,i,j) = rho(i,j)*(-2.0+3.0*rho(i,j)*(u(i,j)*u(i,j)+v(i,j)*v(i,j)))
                    fmeq(2,i,j) = rho(i,j)*(1.0-3.0*rho(i,j)*(u(i,j)*u(i,j)+v(i,j)*v(i,j)))
                    fmeq(3,i,j) = rho(i,j)*u(i,j)
                    fmeq(4,i,j) = -rho(i,j)*u(i,j)
                    fmeq(5,i,j) = rho(i,j)*v(i,j)
                    fmeq(6,i,j) = -rho(i,j)*v(i,j)
                    fmeq(7,i,j) = rho(i,j)*(u(i,j)*u(i,j)-v(i,j)*v(i,j))
                    fmeq(8,i,j) = rho(i,j)*u(i,j)*v(i,j)
                end do
            end do

            ! calculate the moment:
            do i=1,n
                do j=1,m
                    do k=0,8
                        suma  = 0.0
                        do l=0,8
                            suma = suma + tm(k,l)*f(l,i,j)
                        end do
                        fmom(k,i,j) = suma
                    end do
                end do
            end do

            ! collision

            do i=1,n
                do j=1,m
                    do k=0,8
                        sumb = 0.0
                        do l=0,8
                            sumb = sumb +stmiv(k,l)*(fmom(l,i,j)-fmeq(l,i,j))
                        end do
                        f(k,i,j) = f(k,i,j)-sumb
                    end do
                end do
            end do

    end subroutine collesion

    subroutine streaming(f,n,m)
        integer , INTENT(IN) :: n , m
        doubleprecision , INTENT(INOUT) :: f(0:8,1:n,1:m)
        doubleprecision :: ff(0:8,1:n,1:m)
        integer :: i , j , ip, jp , in, jn

        do i=1,n
            do j=1,m

                ip = mod(i,n)+1
                in = n - mod(n+1-i,n)
                jp = mod(j,m)+1
                jn = m - mod(m+1-j,m)

                ff(0,i,j) = f(0,i,j)
                ff(1,ip,j) = f(1,i,j)
                ff(2,i,jp) = f(2,i,j)
                ff(3,in,j) = f(3,i,j)
                ff(4,i,jn) = f(4,i,j)
                ff(5,ip,jp) = f(5,i,j)
                ff(6,in,jp) = f(6,i,j)
                ff(7,in,jn) = f(7,i,j)
                ff(8,ip,jn) = f(8,i,j)

            end do
        end do

        do i=1,n
            do j=1,m
                do k=0,8
                    f(k,i,j) = ff(k,i,j)
                end do
            end do
        end do

    end subroutine streaming

    subroutine boundary(f,n,m,uo)
        doubleprecision , INTENT(INOUT) :: f(0:8,1:n,1:m)
        integer , INTENT(IN) :: n , m
        real , INTENT (IN)   :: uo
        integer :: i, j , k
        real :: rhon , t1 , t2

        ! bounce back on south boundary
        do i=1,n
            f(2,i,1)=f(4,i,1)
            f(5,i,1)=f(7,i,1)
            f(6,i,1)=f(8,i,1)
        end do
        ! moving lid, north boundary
        do i=1,n
            rhon=f(0,i,m)+f(1,i,m)+f(3,i,m)+2.*(f(2,i,m)+f(6,i,m)+f(5,i,m))
            f(4,i,m)=f(2,i,m)
            f(8,i,m)=f(6,i,m)+rhon*uo/6.0
            f(7,i,m)=f(5,i,m)-rhon*uo/6.0
        end do

    end subroutine boundary


    subroutine macro(f,rho,u,v,cx,cy,n,m,uo)
        integer , INTENT(IN) :: m , n
        doubleprecision , INTENT(IN) :: f(0:8,1:n,1:m)
        doubleprecision , INTENT(INOUT) :: rho(1:n,1:m),u(1:n,1:m),v(1:n,1:m)
        real , INTENT(IN) :: cx(0:8),cy(0:8),uo
        integer :: i , j , k
        real :: ssum , usum , vsum

        do j=1,m
            do i=1,n
                ssum=0.0
            do k=0,8
                ssum=ssum+f(k,i,j)
            end do
                rho(i,j)=ssum
            end do
        end do

        DO i=1,n
            DO j=1,m
                usum=0.0
                vsum=0.0
                DO k=0,8
                    usum=usum+f(k,i,j)*cx(k)
                    vsum=vsum+f(k,i,j)*cy(k)
                END DO
                u(i,j)=usum/rho(i,j)
                v(i,j)=vsum/rho(i,j)
            END DO
        END DO

        do j=1,m
            v(n,j) = 0.0
        end do

    end subroutine macro

    subroutine mesh(n,m,X,Y,dx,dy)
        integer , INTENT(IN) :: m , n
        real , INTENT(IN) :: dx , dy
        real ,INTENT(OUT) :: X(1:n,1:m) , Y(1:n,1:m)
        integer :: i,j

        do i=1,n
            do j=1,m
                X(i,j) = (real(i)-1.0)*dx
                Y(i,j) = (real(j)-1.0)*dy
            end do
        end do

    end subroutine mesh

    subroutine output(n,m,dx,dy,rho,u,v,X,Y)
        integer , INTENT(IN) :: m , n
        real , INTENT(IN) :: dx , dy
        doubleprecision , INTENT(IN) :: rho(1:n,1:m),u(1:n,1:m),v(1:n,1:m)
        real , INTENT(IN):: X(1:n,1:m) , Y(1:n,1:m)
        integer :: i,j

        open(unit=1,file="output.dat")

        write(1, '(A)') 'TITLE = "output"'
        write(1, '(A)') 'VARIABLES = "X", "Y" , "U" , "V"'

        write(1, '(A, I0, A, I0, A)') 'ZONE I=', n, ', J=', m, ', F=POINT'

        do j=1,m
            do i=1,n
                write(1,*) X(i,j) , Y(i,j) , u(i,j) , v(i,j)
            end do
        end do

    end subroutine output

    subroutine unsteady(m,u,Y,uo,kk)
        integer , INTENT(IN) :: m,kk
        doubleprecision , INTENT(IN) :: u(1:n,1:m)
        real , INTENT(IN):: Y(1:n,1:m) , uo
        integer :: j,j1

        write(2, *) 'ZONE t="T=',kk,'"'

        do j=0,m-1
            j1 = m-j
            write(2,*) Y(int(n/2),j1)/real(m-1) , u(int(n/2),j+1)/uo
        end do

        print* , "write iteration has been done!"

    end subroutine

    subroutine section(n,m,u,X,Y,uo)
        integer,INTENT(IN) :: n,m
        doubleprecision,INTENT(IN)::u(1:n,1:m)
        real,INTENT(IN):: X(1:n,1:m) , Y(1:n,1:m),uo
        real::t,H,nu_vis,eta,eta1,Re_vis,U_plate,deltaT,deltaY,U_unsteady(1:m,0:5000),Yn(1:m),&
            & suma
        integer ::i,j,k,Tmax,Nmax

        Nmax = 1000
        Tmax = 1000
        deltaT = 0.01
        H = 0.04
        deltaY = H/m
        Re_vis = 10.0
        U_plate = 0.01

        do j=1,m
            Yn(j) = deltaY*(j-1.0)
        end do

        nu_vis = U_plate*H/Re_vis

        open(unit=1,file="varticalSEC.dat") ! x=0.5
        write(1,*) 'VARIABLES = "U/Uo" , "Y/H" '

        t = 0.01

        do k=0,Tmax
            eta1 = 0.5*(H/sqrt(nu_vis*t))
            do j=1,m
                suma = 0.0
                eta = 0.5*(Yn(j)/sqrt(nu_vis*t))
                do i=0,Nmax
                    suma = suma + erfc(2*i*eta1+eta)-erfc(2*(i+1)*eta1-eta)
                end do
            U_unsteady(j,k) = suma
            end do
            t = t + deltaT
        end do

        t=0.0
        do k=0,Tmax
            write(1,*) 'ZONE t="T = ' , t ,' " '
            do j=1,m
                write(1,*) Yn(j)/H , U_unsteady(j,k)
            end do
            t = t + deltaT
        end do

    end subroutine

end program couett

