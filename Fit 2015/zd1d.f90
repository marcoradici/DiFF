! =================================================
!  parametrization for chiral-even DiFF D1d(z,Mh;Q0) with Q0 = 1 GeV
!  from  Phys. Rev. D85 (2012) 114023, arXiv:1202.0323 [hep-ph]
!  input z , parameter Mh 
!  output  z D1d(z,Mh;Q0)
! =================================================

double precision function zd1d(z,m)

	implicit none

    double precision, intent(in):: z, m
	
	double precision:: bkg, rho, omega, kaon

	double precision, parameter:: pi = 4.d0*atan(1.d0)
	double precision, parameter:: mpi = 0.13957d0 
	double precision, parameter:: mrho = 0.776d0, grho = 0.150d0
	double precision, parameter:: mome = 0.783d0, gome = 0.008d0  

!  parameters for background channel
	double precision, parameter:: N_BG1=0.60126d0, B_BG1=0.84461d0, A_BG1=-2.2824d0, A_BG2=1.0012d0, &
   &                              G_BG1=0.71333d0, G_BG2=-0.15473d0, G_BG3=1.1798d0, G_BG4=-1.0506d0
   
!  parameters for rho resonant channel
	double precision, parameter:: N_RH1=0.20851d0, B_RH1=0.99914d0, A_RH1=0.10423d0, A_RH2=-1.2095d0, &
   &                              G_RH1=4.0446d0, G_RH2=-15.679d0, G_RH3=20.582d0, E_RH1=1.1035d0, &
   &                              D_RH1=-1.0668d0, D_RH2=-1.3569d0
   
!  parameters for the omega channel
	double precision, parameter:: N_OM1=3.234d14, B_OM1=12.539d0, B_OM2=0.28995d0, A_OM1=1.2197d0, &
   &                              G_OM1=1.970d0, G_OM2=31.032d0, G_OM3=10.228d0, E_OM1=0.038805d0, &
   &                              D_OM1=-0.86237d0, D_OM2=-0.27902d0
   
!  parameters for the kaon resonant channel
	double precision, parameter:: N_KA1=0.19076d0, E_KA1=0.063383d0, &
   &                              G_KA1=0.20997d0, G_KA2=5.2428d0, G_KA3=-2.795d0, G_KA4=-5.2705d0, &
   &                              D_KA1=2.3843d0, D_KA2=-5.0427d0, D_KA3=0.63341d0, &
   &                              N_KA2=1.373d0, A_KA1=0.42636d0

   
    bkg = N_BG1 * z**A_BG1 * (1.d0-z)**(A_BG2**2) *   &
   &      dexp(-(G_BG1/z+G_BG2+G_BG3*z+G_BG4/m)**2 *  &
   &            (m**2-4.d0*mpi**2)) *                 &
   &      dsqrt(m**2-4.d0*mpi**2)**(B_BG1**2)

    rho = N_RH1**2 * z**A_RH1 * (1.d0-z)**(A_RH2**2) *                       &
   &      dsqrt(m**2-4.d0*mpi**2)**(B_RH1**2) *                              &
   &      (  dexp(-(G_RH1/z+G_RH2+G_RH3*z-(G_RH1+G_RH2+G_RH3)*z**3)*m**2) *  &
   &         dexp(-(D_RH1/z+D_RH2*z)) +                                      &
   &         E_RH1**2 / ((m**2-mrho**2)**2 + mrho**2 * grho**2)  ) 

    omega = (1.d0-z)**(A_OM1**2) *                                        &
   &        dsqrt(m**2-4.d0*mpi**2)**B_OM1/(1.d0+dexp(5.d0*(m-1.2d0))) *  &
   &        (  N_OM1 * dexp(-(D_OM1/z+D_OM2*z)) *                         &
   &           dexp(-(G_OM1/z+G_OM2+G_OM3*z+0.d0*B_OM2*m) *               &
   &                 (m**2-4.d0*mpi**2)**B_OM2) +                         &
   &           E_OM1**2 / ((m**2-mome**2)**2 + mome**2 * gome**2)  )

    if( m >= 0.49d0  .and.  m < (0.49d0+0.01d0) ) then
      kaon = N_KA1**2 *4.d0* dexp(G_KA1/z+G_KA2+G_KA3*z+G_KA4*z**2)
    else
      kaon = E_KA1**2 * dexp(1.d0+D_KA1*m+D_KA2*m**2+D_KA3*m*z) *       &
   &         dexp(G_KA1/z+G_KA2+G_KA3*z+G_KA4*z**2) * dsqrt(m**2-4.d0*mpi**2)
    end if 
	
	zd1d = z * ( bkg + rho + omega + N_KA2**2 * z**A_KA1 * kaon )
	
end function zd1d
