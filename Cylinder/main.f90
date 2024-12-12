!	 ========================================================
!	 The main program, implementing a flow past a cylinder
!	 ========================================================

PROGRAM unsteady
    USE simParam, ONLY: xDim, yDim, tMax
    implicit none

    double precision:: omega, time1, time2, timeTot
    double precision, dimension(:,:,:), allocatable:: f, fEq, u,u_old
    double precision, dimension(:,:), allocatable:: rho,rho_old,uSqr
    integer, dimension(:,:), allocatable:: image
    integer:: tStep

    allocate(f(yDim,xDim,0:8))
    allocate(fEq(yDim,xDim,0:8))
    allocate(u_old(yDim,xDim,0:1))
    allocate(u(yDim,xDim,0:1))
    allocate(uSqr(yDim,xDim))
    allocate(rho(yDim,xDim))
    allocate(rho_old(yDim,xDim))
    allocate(image(yDim,xDim))

    CALL constructImage(image)
    CALL computeOmega(omega)
    CALL writeInput(omega)
    CALL initMacro(rho,u,uSqr)
    CALL computeFeq(fEq,rho,u,uSqr)

    f = fEq

    open(14,file='output.dat')
    write(14,*) "VARIABLE = X , Y , U , V , image"

    timeTot = 0.0d0
    do tStep = 1, tMax
        CALL CPU_TIME(time1)

        rho_old = rho
        u_old(:,:,0) = u(:,:,0)
        u_old(:,:,1) = u(:,:,1)

        CALL inletOutlet(f,rho,u,image)

        CALL boundaries(f,image)

        CALL computeMacros(f,rho,u,uSqr)

        CALL computeFeq(fEq,rho,u,uSqr)

        CALL collide(f,fEq,omega,image)

        CALL stream(f)

        IF(mod(tStep , 20) .EQ. 0 ) THEN
            write(*,*) "====================================================="
            write(*,*) "rho residual : " , maxval(abs(rho-rho_old))
            write(*,*) "U   residual : " , maxval(abs(u(:,:,0)-u_old(:,:,0)))
            write(*,*) "V   residual : " , maxval(abs(u(:,:,1)-u_old(:,:,1)))
            write(*,*) "time step    : " , tStep
            write(*,*) "====================================================="
            CALL writeOutput(u,image,tStep)
        ENDIF
        CALL CPU_TIME(time2)
        timeTot = timeTot + (time2-time1)
    end do

    write(*,*) dble(tMax) * (dble(yDim * xDim)) / timeTot ,'cells per second'

    deallocate(f)
    deallocate(fEq)
    deallocate(u)
    deallocate(uSqr)
    deallocate(rho)
    deallocate(image)

END PROGRAM unsteady


!	 ========================================================
!	 Compute the relaxation parameter from the Reynolds number
!	 ========================================================
SUBROUTINE computeOmega(omega)
    USE simParam, ONLY: Re,uMax,obstR

    implicit none

    double precision, INTENT(INOUT):: omega
    double precision:: nu

    nu    =  uMax * 2.0d0 * dble(obstR) / Re
    omega = 1.0d0 / (3.0d0*nu+0.5d0)
END SUBROUTINE computeOmega


!	 ========================================================
!	 Construct an array the defines the flow geometry
!	 ========================================================
SUBROUTINE constructImage(image)
    USE cellConst
    USE simParam, ONLY: xDim, yDim, obstX, obstY, obstR
    USE D2Q9Const, ONLY: v

    implicit none

    integer, INTENT(INOUT):: image(yDim,xDim)
    integer:: x,y

    v(0:8,0) = (/0,1,0,-1,0,1,-1,-1,1/)
    v(0:8,1) = (/0,0,1,0,-1,1,1,-1,-1/)

    image          = fluid
    image(:,1)     = inlet
    image(:,xDim)  = outlet
    image(1,:)     = wall
    image(yDim,:)  = wall
    do x = 1, xDim
        do y = 1, yDim
            if (((x-obstX)**2 + (y-obstY)**2) <= (obstR**2) ) image(y,x) = wall
        end do
    end do

END SUBROUTINE constructImage


!	 ========================================================
!	 Initialize the simulation to Poiseuille profile at
!        an equilibrium distribution
!	 ========================================================
SUBROUTINE initMacro(rho,u,uSqr)
    USE simParam, ONLY: xDim, yDim

    implicit none

    double precision, INTENT(INOUT):: rho(yDim,xDim), u(yDim,xDim,0:1), uSqr(yDim,xDim)
    double precision:: uProf
    integer:: y

    do y = 1, yDim
        u(y,:,0) = uProf(y)
        u(y,:,1) = 0.0d0
    end do
    rho  = 1.0d0
    uSqr = u(:,:,0) * u(:,:,0) + u(:,:,1) * u(:,:,1)
END SUBROUTINE initMacro


!	 ========================================================
!	 Compute equilibrium distribution
!	 ========================================================
SUBROUTINE computeFeq(fEq,rho,u,uSqr)
    USE D2Q9COnst, ONLY: t, v
    USE simParam, ONLY: xDim, yDim
    implicit none

    double precision, INTENT(IN):: rho(yDim,xDim), uSqr(yDim,xDim), u(yDim,xDim,0:1)
    double precision, INTENT(INOUT):: fEq(yDim,xDim,0:8)
    integer:: i, x, y
    double precision:: uxy

    do i = 0, 8
        do x = 1, xDim
            do y = 1, yDim
                uxy = u(y,x,0) * v(i,0) + u(y,x,1) * v(i,1)
                fEq(y,x,i) = t(i) * rho(y,x) * (1.0d0 + 3.0d0 * uxy + 4.5d0 * uxy * uxy - 1.5d0 * uSqr(y,x))
            end do
        end do
    end do
END SUBROUTINE computeFeq


!	 ========================================================
!	 Compute density and velocity from distribution functions
!	 ========================================================
SUBROUTINE computeMacros(f,rho,u,uSqr)
    USE simParam, ONLY: xDIm, yDim
    implicit none

    double precision, INTENT(IN):: f(yDim,xDim,0:8)
    double precision, INTENT(INOUT):: u(yDim,xDim,0:1), rho(yDim, xDim), uSqr(yDim, xDim)
    integer:: x,y

    do x = 1, xDim
        do y = 1, yDim
            rho(y,x)  = f(y,x,0) + f(y,x,1) + f(y,x,2) + f(y,x,3) + f(y,x,4) + f(y,x,5) + f(y,x,6) + f(y,x,7) + f(y,x,8)
            u(y,x,0)  = (f(y,x,1) - f(y,x,3) + f(y,x,5) - f(y,x,6) - f(y,x,7) + f(y,x,8)) / rho(y,x)
            u(y,x,1)  = (f(y,x,2) - f(y,x,4) + f(y,x,5) + f(y,x,6) - f(y,x,7) - f(y,x,8)) / rho(y,x)
            uSqr(y,x) = u(y,x,0) * u(y,x,0) + u(y,x,1) * u(y,x,1)
        end do
    end do
END SUBROUTINE computeMacros


!	 ========================================================
!	 Implement Bounce-back on upper/lower boundaries
!	 ========================================================
SUBROUTINE boundaries(f,image)
    USE D2Q9Const, ONLY: opposite
    USE cellConst, ONLY: wall
    USE simParam, ONLY: xDim, yDim
    implicit none

    integer, INTENT(IN):: image(yDim,xDim)
    double precision, INTENT(INOUT):: f(yDim,xDim,0:8)
    double precision:: fTmp(0:8)
    integer:: i, x, y

    do x = 1, xDim
        do y = 1, yDim
            if (image(y,x) == wall) then
                do i = 0, 8
                    fTmp(i) = f(y,x,opposite(i))
                end do
                do i = 0, 8
                    f(y,x,i) = fTmp(i)
                end do
            end if
        end do
    end do
END SUBROUTINE boundaries


!	 ========================================================
!	 Use Zou/He boundary condition to implement Dirichlet
!        boundaries on inlet/outlet
!	 ========================================================
SUBROUTINE inletOutlet(f,rho,u,image)
    USE cellConst, ONLY: inlet, outlet
    USE simParam

    implicit none

    double precision, INTENT(INOUT):: f(yDim,xDim,0:8), u(yDim,xDim,0:1), rho(yDim,xDim)
    integer, INTENT(IN):: image(yDim,xDim)

    double precision:: uProf
    integer:: x, y

    do x = 1, xDim
        do y = 1, yDim
            if (image(y,x) == inlet) then
                u(y,x,0) = uProf(y)
                u(y,x,1) = 0.0d0
                CALL inletZou(f(y,x,:),u(y,x,:),rho(y,x))
            else if (image(y,x) == outlet) then
                u(y,x,0) = uProf(y)
                u(y,x,1) = 0.0d0
                CALL outletZou(f(y,x,:),u(y,x,:),rho(y,x))
            end if
        end do
    end do

CONTAINS


    !	 ========================================================
    !	 Zou/He boundary on inlet
    !	 ========================================================
    SUBROUTINE inletZou(f,u,rho)
        implicit none

        double precision, INTENT(INOUT):: f(0:8),rho
        double precision, INTENT(IN):: u(0:1)
        double precision:: fInt, fInt2

        fInt   = f(0) + f(2) + f(4)
        fInt2  = f(3) + f(6) + f(7)
        rho    = (fInt + 2.0d0 * fInt2) / (1.0d0 - u(0))
        CALL zouWestWall(f,rho,u)
    END SUBROUTINE inletZou

    SUBROUTINE zouWestWall(f,rho,u)
        implicit none

        double precision, INTENT(INOUT):: f(0:8)
        double precision, INTENT(IN):: rho, u(0:1)
        double precision:: fDiff, rhoUx, rhoUy

        fDiff = 0.5d0 * (f(2) - f(4))
        rhoUx = rho * u(0) / 6.0d0
        rhoUy = 0.5d0 * rho * u(1)

        f(1) = f(3) + 4.0d0 * rhoUx
        f(5) = f(7) - fDiff + rhoUx + rhoUy
        f(8) = f(6) + fDiff + rhoUx - rhoUy
    END SUBROUTINE zouWestWall


    !	 ========================================================
    !	 Zou/He boundary on outlet
    !	 ========================================================
    SUBROUTINE outletZou(f,u,rho)
        implicit none

        double precision, INTENT(INOUT):: f(0:8),rho,u(0:1)
        double precision:: fInt, fInt2

        fInt  = f(0) + f(2) + f(4)
        fInt2 = f(1) + f(8) + f(5)
        rho   = (fInt + 2.0d0 * fInt2) / (1.0d0 + u(0))
        CALL zouEastWall(f,rho,u)
    END SUBROUTINE outletZou

    SUBROUTINE zouEastWall(f,rho,u)
        implicit none

        double precision, INTENT(INOUT):: f(0:8)
        double precision, INTENT(IN):: rho, u(0:1)
        double precision:: fDiff, rhoUx, rhoUy

        fDiff = 0.5d0 * (f(2) - f(4))
        rhoUx = rho * u(0) / 6.0d0
        rhoUy = 0.5d0 * rho * u(1)

        f(3) = f(1) - 4.0d0 * rhoUx
        f(7) = f(5) + fDiff - rhoUx - rhoUy
        f(6) = f(8) - fDiff - rhoUx + rhoUy
    END SUBROUTINE zouEastWall

END SUBROUTINE inletOutlet


!	 ========================================================
!	 Computation of Poiseuille profile for the inlet/outlet
!	 ========================================================
FUNCTION uProf(y)
    USE simParam, ONLY: yDIm, uMax
    implicit none

    integer, INTENT(IN):: y
    double precision:: radius, uProf

    radius = dble(yDim-1) * 0.5d0
    uProf  = -uMax * ((abs(1 - dble(y-1) / radius))**2 - 1.0d0)
END FUNCTION uProf


!	 ========================================================
!	 Streaming step: the population functions are shifted
!        one site along their corresponding lattice direction
!        (no temporary memory is needed)
!	 ========================================================
SUBROUTINE stream(f)
    USE simParam
    implicit none

    double precision, INTENT(INOUT):: f(yDim,xDim,0:8)
    double precision:: periodicHor(yDim), periodicVert(xDim)

!	 -------------------------------------
!	 right direction
    periodicHor   = f(:,xDim,1)
    f(:,2:xDim,1) = f(:,1:xDim-1,1)
    f(:,1,1)      = periodicHor
!	 -------------------------------------
!	 up direction
    periodicVert    = f(1,:,2)
    f(1:yDim-1,:,2) = f(2:yDim,:,2)
    f(yDim,:,2)     = periodicVert
!	 -------------------------------------
!	 left direction
    periodicHor     = f(:,1,3)
    f(:,1:xDim-1,3) = f(:,2:xDim,3)
    f(:,xDim,3)     = periodicHor
!	 -------------------------------------
!	 down direction
    periodicVert  = f(yDim,:,4)
    f(2:yDim,:,4) = f(1:yDim-1,:,4)
    f(1,:,4)      = periodicVert
!	 -------------------------------------
!	 up-right direction
    periodicVert = f(1,:,5)
    periodicHor  = f(:,xDim,5)
    f(1:yDim-1,2:xDim,5) = f(2:yDim,1:xDim-1,5)
    f(yDim,2:xDim,5)     = periodicVert(1:xDim-1)
    f(yDim,1,5)          = periodicVert(xDim)
    f(1:yDim-1,1,5)      = periodicHor(2:yDim)
!	 -------------------------------------
!	 up-left direction
    periodicVert = f(1,:,6)
    periodicHor  = f(:,1,6)
    f(1:yDim-1,1:xDim-1,6) = f(2:yDim,2:xDim,6)
    f(yDim,1:xDim-1,6)     = periodicVert(2:xDim)
    f(yDim,xDim,6)         = periodicVert(1)
    f(1:yDim-1,xDim,6)     = periodicHor(2:yDim)
!	 -------------------------------------
!	 down-left direction
    periodicVert = f(yDim,:,7)
    periodicHor  = f(:,1,7)
    f(2:yDim,1:xDim-1,7) = f(1:yDim-1,2:xDim,7)
    f(1,1:xDim-1,7)      = periodicVert(2:xDim)
    f(1,xDim,7)          = periodicVert(1)
    f(2:yDim,xDim,7)     = periodicHor(1:yDim-1)
!	 -------------------------------------
!	 down-right direction
    periodicVert = f(yDim,:,8)
    periodicHor  = f(:,xDim,8)
    f(2:yDim,2:xDim,8) = f(1:yDim-1,1:xDim-1,8)
    f(1,2:xDim,8)      = periodicVert(1:xDim-1)
    f(1,1,8)           = periodicVert(xDim)
    f(2:yDim,1,8)      = periodicHor(1:yDim-1)
END SUBROUTINE stream


!	 ========================================================
!	 LBGK collision step
!	 ========================================================
SUBROUTINE collide(f,fEq,omega,image)
    USE simParam, ONLY: xDim, yDim
    USE cellConst, ONLY: wall
    implicit none

    integer, INTENT(IN):: image(yDim,xDim)
    double precision, INTENT(IN):: fEq(yDim,xDim,0:8), omega
    double precision, INTENT(INOUT):: f(yDim,xDim,0:8)

    integer:: x,y,i

    do i = 0, 8
        do x = 1, xDim
            do y = 1, yDim
                if (image(y,x) /= wall) f(y,x,i) = (1.0d0 - omega) * f(y,x,i) + omega * feq(y,x,i)
            end do
        end do
    end do
END SUBROUTINE collide


!	 ========================================================
!	 Write the components of the velocity to a text file,
!        with indices (x,y)
!	 ========================================================
SUBROUTINE writeOutput(u,image,tStep)
    USE simParam, ONLY: xDim, yDim
    implicit none

    integer, INTENT(IN):: tStep
    double precision, INTENT(IN):: u(yDim,xDim,0:1)
    integer, INTENT(IN):: image(yDim,xDim)
    integer:: x,y

    write(14,*) 'ZONE ' , ',I = ' , xDim , ', J = ' , yDim

    do y=1, yDim
        do x=1, xDim
            write(14,*) x,y,u(y,x,0),u(y,x,1),image(y,x)
        end do
    end do

END SUBROUTINE writeOutput

!	 ========================================================
!	 Print out simulation parameters to screen
!	 ========================================================
SUBROUTINE writeInput(omega)
    USE simParam
    implicit none

    double precision, INTENT(IN):: omega

    write(*,*) 'xDim                 = ', xDim
    write(*,*) 'yDim                 = ', yDim
    write(*,*) 'Obstacle X           = ', obstX
    write(*,*) 'Obstacle Y           = ', obstY
    write(*,*) 'Obstacle Radius      = ', obstR
    write(*,*) 'tMax                 = ', tMax
    write(*,*) 'uMax                 = ', uMax
    write(*,*) 'Re                   = ', Re
    write(*,*) 'omega                = ', omega
END SUBROUTINE writeInput
