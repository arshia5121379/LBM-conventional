!	 ========================================================
!	 Lattice constants for the D2Q9 lattice
!	 ========================================================
MODULE D2Q9Const
!	 D2Q9 Weights
    double precision,parameter:: t(0:8) = (/4.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0&
                                           &,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0/)
!	D2Q9 Directions
    integer:: v(0:8,0:1)
!       = (/(/0,1,0,-1,0,1,-1,-1,1/),(/0,0,1,0,-1,1,1,-1,-1/)/)

    integer, parameter:: opposite(0:8) = (/0,3,4,1,2,7,8,5,6/)
END MODULE D2Q9Const

