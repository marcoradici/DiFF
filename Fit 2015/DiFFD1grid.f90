! Code to produce a grid for unpolarized DiFF D1
!
! Reference: Phys. Rev. D85 (2012) 114023, arXiv:1202.0323 [hep-ph]
!            JHEP 1505 (2015) 123, arXiv:1503.03495v2 [hep-ph] 
! 
!  output format:   zD1_PV15_XX_YY.dat    where XX = LO / NLO using MSTW08 parametrization and 
!                                         YY = 1    D1g(Q0) = 0
!                                              2    D1g(Q0) = D1u(Q0) / 4
!                                              3    D1g(Q0) = D1u(Q0)
! 
!  for each file the output is: grid in Mh , grid in Q , grid in z , then zD1q  with q = t-bar,b-bar,c-bar,s-bar,u-bar,d-bar,g,d,u,s,c,b,t 
!  where  1< Q< 100 GeV , 0.3< Mh< 1.3 GeV , 0.2< z< 1
!
!  At Q0=1 GeV, u=d=u-bar=d-bar & s=s-bar & c=c-bar 		  
!
!  compilation: see Makefile 

program DiFFD1grid
	
	use evo_grid

	implicit none
	
	integer:: replica, i, j, k, l
	double precision:: p(1), zD1(-6:6), zd1u, zd1d, zd1s, zd1c
	character(len=27):: fname   ! LO
!	character(len=29):: fname   ! NLO

	

!    prepare grids in M , Q , z
	call compute_grids

! 	prepare output in zD1_PV15_XX_YY.dat with  XX=LO/NLO and YY= 1 (D1g=0) / 2 (D1g=D1u/4) / 3 (D1g=D1u) 
    write(fname,'(a,i1,a)')  'MSTW08/LO/zD1_PV15_LO_',1,'.dat'  !  D1g(Q0)=0   
!    write(fname,'(a,i1,a)')  'MSTW08/NLO/zD1_PV15_NLO_',1,'.dat'  !  D1g(Q0)=0   
!    write(fname,'(a,i1,a)')  'MSTW08/LO/zD1_PV15_LO_',2,'.dat'  !  D1g(Q0)=D1u(Q0)/4   
!    write(fname,'(a,i1,a)')  'MSTW08/NLO/zD1_PV15_NLO_',2,'.dat'  !  D1g(Q0)=D1u(Q0)/4  
!    write(fname,'(a,i1,a)')  'MSTW08/LO/zD1_PV15_LO_',3,'.dat'  !  D1g(Q0)=D1u(Q0) 
!    write(fname,'(a,i1,a)')  'MSTW08/NLO/zD1_PV15_NLO_',3,'.dat'  !  D1g(Q0)=D1u(Q0)  

	open(unit=1,file=fname,status='unknown') 

    write(1,*) ' Mh '
    write(1,'(5(e11.5,1x))') (Mvec(i),i=1,nM)
    write(1,*)
    write(1,*) ' Q ' 
    write(1,'(5(e11.5,1x))') (Qvec(i),i=1,nQ)
    write(1,*)
    write(1,*) ' z '
    write(1,'(5(e11.5,1x))') (zvec(i),i=1,nz)
    write(1,*)
    write(1,*)
    write(1,*) ' anti-t       anti-b       anti-c       anti-s      ', &
   &           ' anti-u       anti-d       g            d           ', &
   &           ' u            s            c            b            t'
	
	
!   loop on invariant mass, treated as input parameter (not affected by evolution)
	do i=1,nM

	  p(1) = Mvec(i)
	  
! 	Initialize APFEL++
	  call InitialiseEvolution(p)

!   evolve zD1 at predefined grid in (Q,z) and write output 
	
	  do j=1,nQ
        do k=1,nz
 
		  do l=-6,6
	  	    zD1(l) = 0.d0
		  end do
         
		  if(j == 1) then
            zD1(1) = zd1d(zvec(k),p(1))
			zD1(-1) = zD1(1)
            zD1(2) = zd1u(zvec(k),p(1))
			zD1(-2) = zD1(2)
            zD1(3) = zd1s(zvec(k),p(1))
			zD1(-3) = zD1(3)
            zD1(4) = zd1c(zvec(k),p(1))
			zD1(-4) = zD1(4)
			zD1(0) = 0.d0
!			zD1(0) = zD1(2) / 4.d0
!			zD1(0) = zD1(2)	
          else	 
            call EvolveTransversities(zvec(k), Qvec(j), zD1)
          end if

		  write(1,'(13(1x,e12.5))') (zD1(l),l=-6,6)

    	end do
	  end do
	end do

	close(unit=1)

	stop

end program DiFFD1grid