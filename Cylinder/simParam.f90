
!	 ========================================================
!	 Constants for simulation setup
!	 ========================================================
MODULE simParam
    integer, parameter:: xDim = 250
    integer, parameter:: yDim = 60
    integer, parameter:: obstX = xDim/5
    integer, parameter:: obstY = yDim/2
    integer, parameter:: obstR = 5

    integer, parameter:: tMax = 5000

    double precision, parameter:: uMax = 0.05d0
    double precision, parameter:: Re = 50.0d0
END MODULE simParam
