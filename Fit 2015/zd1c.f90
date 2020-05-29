! =================================================
!  parametrization for chiral-even DiFF D1c(z,Mh;Q0) with Q0 = 1 GeV
!  from  Phys. Rev. D85 (2012) 114023, arXiv:1202.0323 [hep-ph]
!  input z , parameter Mh 
!  output  z D1c(z,Mh;Q0)
! =================================================

double precision function zd1c(z,m)

	implicit none

    double precision, intent(in):: z, m
	
	double precision:: bkg, rho, omega, kaon

	double precision, parameter:: pi = 4.d0*atan(1.d0)
	double precision, parameter:: mpi = 0.13957d0 
	double precision, parameter:: mrho = 0.776d0, grho = 0.150d0
	double precision, parameter:: mome = 0.783d0, gome = 0.008d0  

!  parameters for background channel
	double precision, parameter:: N_BG3=1.4373d0, B_BG2=0.94002d0, A_BG4=-2.3096d0, A_BG5=1.7020d0, &
   &                              G_BG5=0.63359d0, G_BG6=0.81632d0, G_BG7=-0.64550d0
   
!  parameters for rho resonant channel
	double precision, parameter:: N_RH3=0.44995d0, B_RH2=0.69703d0, A_RH4=1.8497d0, A_RH5=2.4739d0, &
   &                              G_RH4=3.9579d0, E_RH2=2.2231d0, D_RH3=-1.2198d0, D_RH4=3.7210d0
   
!  parameters for the omega channel
	double precision, parameter:: N_OM3=1.7577d13, B_OM3=11.326d0, B_OM4=0.38219d0, A_OM3=1.8371d0, &
   &                              G_OM4=33.268d0, E_OM2=-0.027674d0, D_OM3=0.3381d0, D_OM4=7.7996d0
   
!  parameters for the kaon resonant channel
	double precision, parameter:: N_KA4=0.59628d0, E_KA2=0.10853d0, &
   &                              G_KA5=0.43489d0, G_KA6=1.9868d0, G_KA7=3.6238d0, G_KA8=-11.641d0, &
   &                              D_KA4=2.7231d0, D_KA5=-5.1218d0, D_KA6=-0.18046d0

   
    bkg = N_BG3 * z**A_BG4 * (1.d0-z)**(A_BG5**2) *                   &
   &      dsqrt(m**2-4.d0*mpi**2)**(B_BG2**2) *                       &
   &      dexp(-(G_BG5/z+G_BG6*z+G_BG7/m)**2 * (m**2-4.d0*mpi**2))

    rho = N_RH3**2 * z**A_RH4 * (1.d0-z)**(A_RH5**2) *               &
   &      dsqrt(m**2-4.d0*mpi**2)**(B_RH2**2) *                      &
   &      (  dexp(-(G_RH4-G_RH4*z**3) * m**2) *                      &
   &         dexp(-(D_RH3/z+D_RH4*z)) +                              &
   &         E_RH2**2 / ((m**2-mrho**2)**2 + mrho**2 * grho**2)  )

    omega = (1.d0-z)**(A_OM3**2) *                                          &
   &        dsqrt(m**2-4.d0*mpi**2)**B_OM3 / (1.d0+dexp(5.d0*(m-1.2d0))) *  &
   &        (  N_OM3 * dexp(-(D_OM3/z+D_OM4*z)) *                           &
   &           dexp(-(G_OM4+0.d0*B_OM4*m)*(m**2-4.d0*mpi**2)**B_OM4) +      &
   &           E_OM2**2 / ((m**2-mome**2)**2 + mome**2 * gome**2)  )

    if( m >= 0.49d0   .and.   m < (0.49d0+0.01d0) ) then
      kaon = dexp(G_KA5/z+G_KA6+G_KA7*z+G_KA8*z**2) * N_KA4**2 * 2.d0
    else
      kaon = dsqrt(m**2-4.d0*mpi**2) *                           &
   &         dexp(G_KA5/z+G_KA6+G_KA7*z+G_KA8*z**2) *            &
   &         E_KA2**2 * dexp(1.d0+D_KA4*m+D_KA5*m**2+D_KA6*m*z)
    end if 
	
	zd1c = z * ( bkg + rho + omega + kaon )
	
end function zd1c
