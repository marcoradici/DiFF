module evo_grid   ! calculate grid in 1< Q< 100 GeV, 0.3< Mh< 1.3 GeV, and 0.2< z< 1 for evolution

	implicit none
	
	integer, parameter:: nQ=50, nM=50, nz=30
	double precision, parameter:: qmin=1.d0, qmax=1.d2
	double precision, parameter:: zmin=0.2d0, Mmin=0.3d0, Mmax=1.3d0
	double precision:: Qvec(nQ), Mvec(nM), zvec(nz)

	save
	contains
	
	subroutine compute_grids
	
		integer:: i
	    
	    do i=1,nQ
		  Qvec(i) = 0.d0
	    end do
	    do i=1,nM
		  Mvec(i) = 0.d0
	    end do
	    do i=1,nz
		  zvec(i) = 0.d0
	    end do
	
	    do i=1,nQ
		  Qvec(i) = qmin * (qmax/qmin)**(dble(i-1)/dble(nQ))
	    end do
	
	    do i=1,nM      ! 50 equidistant points
		  Mvec(i) = Mmin + (Mmax-Mmin) * dble(i-1)/dble(nM)
	    end do
	
	    do i=1,nz      ! 30 equidistant points 
		  zvec(i) = zmin + (1.d0-zmin) * dble(i-1)/dble(nz)
	    end do
	
	end subroutine compute_grids

end module evo_grid
	
