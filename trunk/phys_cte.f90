

module phys_cte

  implicit none


  integer,parameter :: sp = kind(1.e0) !single precision
  integer,parameter :: dp = kind(1.d0) !double precision


  integer,parameter :: stdout = 6

  real(dp),parameter :: s2 = 1.414213562d0 
  real(dp),parameter :: s3 = 1.732050808d0
	real(dp),parameter :: s8=2.8284271d0
	real(dp),parameter :: s11=3.31662479d0 
  real(dp),parameter :: pi    = 3.141592654d0 

  real(dp),parameter :: ecb = 1.60217662d-19 ! electron charge (C)
  real(dp),parameter :: hbar = 1.0545718d-34 ! reduced Planck cte (J.s)
  real(dp),parameter :: kb = 8.617330342d-5  ! Boltzmann constant (eV.K-1)
	real(dp),parameter ::eps0=8.854188d-12     ! vaccum permittivity (F.m-1)
	real(dp),parameter ::me0=9.10938356d-31    ! electron mass (kg)
	real(dp),parameter ::clight=2.99792458d8   ! speed of light in vaccum (m.s-1)
	real(dp),parameter ::rydcst=13.605693d0    ! Rydberg constant (eV)	

	real(dp),parameter ::hbarnorm=1.0545716d0  
	real(dp),parameter ::me0norm=9.10938356d0
	real(dp),parameter ::ecbnorm=1.602176487d0
	real(dp),parameter ::Atomet=1.d30          ! to convert between A-3 and m-3



end module phys_cte
