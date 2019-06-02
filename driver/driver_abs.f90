program driver_abs

	use mpi
	use phys_cte
	use setup_epm
	use abs_routine
	
	implicit none
	
	!init params
	real(dp) ::muSO_in,zetaSO_in,Vloc_in(11),ryd_unit
	real(dp) ::a0_in,rbasis_in(3,3),pos1_in(3),pos2_in(3)
	real(dp) ::nr_in,broad_in,TK_in,ndop_in,pdop_in
	integer :: nbvb_in,nbcb_in
	integer,allocatable ::bndvb_in(:),bndcb_in(:)
	
	!params for absorption
	real(dp) ::gcutratio,pol(3),Elow,Eup,lbdsF(2),ubdsF(2)
	integer ::flagabs(3)
	integer ::ii
	
	!offset to highest VB at G
	real(dp) ::kG(3)
	real(dp),allocatable :: EGarr(:)
	real(dp) ::bndoff_in
	complex(dp),allocatable ::EvGarr(:,:)
	!MPI
	integer ::rank,root,nproc,ierr
	
	!Fermi_scan variables
	integer ::nbFc,nbFv
	real(dp),allocatable ::Fcls(:),Fvls(:)
	real(dp) ::Fcbds(2),Fvbds(2)
	
	!gainninj scan
	integer :: nbninj,nEls
	real(dp) ::ninjdown,ninjup
	real(dp),allocatable :: ninjls(:),Els(:)
	
	!gainEph scan
	integer ::nEmap,ndensgEph
	real(dp) ::Emapdown,Emapup
	real(dp),allocatable ::Emap(:),densgEph(:)
	
	!NNNscan
	integer ::nscan
	integer,allocatable ::N1ls(:),N2ls(:),N3ls(:)
	
	!prefix
	character(len=120) ::suffix
	
	!strain params
	real(dp) ::epsval,c11,c12,epsmat(3,3),Iden(3,3),strmode
	real(dp) ::rbasis_raw(3,3)
	

	!open MPI comm world (call in main program)
	call MPI_INIT(ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
	

	!********************************************
	!********************************************
	!DECLARATION OF INPUT PARAMETERS

	!value used in driver_fit during the parameters fitting process
	epscib=0.0125d0 

	!root of MPI
	root=0

	
	!params for EPM model
	!example:Ge
	ryd_unit=13.605693d0 !Rydberg constant, in eV
	muSO_in=ryd_unit*0.00043866d0 !SO parameters, in eV
	zetaSO_in=10.09 !SO parameters, in A0-1
	
	!local empirical pseudopotential, in eV
	Vloc_in=ryd_unit*(/0.d0,-0.26033304d0,-0.25887354d0,-0.25740279d0,&
									&0.03636844d0,0.03991221d0,0.04318571d0,&
									&0.02953685d0,0.03224141d0,0.03521343d0,0.d0/) 
	
	!----------------------------------
	!params for unit cell configuration
	
	! type of applied strain (strmode=1 for biaxe [001],=2 for uniaxe [001]) 
	strmode=1
	!strain value applied
	epsval=0.d0
	!C11 and C12 (GPa)
	c11=129.d0
	c12=48.d0
	
	!Identical matrix
	Iden(1,:)=(/1.d0,0.d0,0.d0/)
	Iden(2,:)=(/0.d0,1.d0,0.d0/)
	Iden(3,:)=(/0.d0,0.d0,1.d0/)


	if (strmode==1) then
		epsmat(:,1)=(/1.d0,0.d0,0.d0/)
		epsmat(:,2)=(/0.d0,1.d0,0.d0/)
		epsmat(:,3)=(/0.d0,0.d0,-2.d0*(c12/c11)/)
	else if (strmode==2) then
		epsmat(:,1)=(/1.d0,0.d0,0.d0/)
		epsmat(:,2)=(/0.d0,-c12/(c12+c11),0.d0/)
		epsmat(:,3)=(/0.d0,0.d0,-c12/(c12+c11)/)
	end if
	
	!cell parameters
	a0_in=5.658d0 !A

	rbasis_raw(:,1)=a0_in*(/0.d0,0.5d0,0.5d0/)
	rbasis_raw(:,2)=a0_in*(/0.5d0,0.d0,0.5d0/)
	rbasis_raw(:,3)=a0_in*(/0.5d0,0.5d0,0.d0/)
	
	rbasis_in=matmul(Iden+epsval*epsmat,rbasis_raw)
	
	!2 positions of atoms in the diamond unit cell
	pos1_in=(0.d0)*(rbasis_in(:,1)+rbasis_in(:,2)+rbasis_in(:,3))
	pos2_in=(1.d0/4.d0)*(rbasis_in(:,1)+rbasis_in(:,2)+rbasis_in(:,3))
	
	!----------------------------------
	!params for calculation of absorption

	nr_in=4.02 !refractive index
	broad_in=0.04d0 !broadening of the Gaussian function (in eV)
	nbvb_in=6 ! Number of valence band taken account
	nbcb_in=8 ! Number of conduction band taken account

	allocate(bndvb_in(nbvb_in),bndcb_in(nbcb_in))
	bndvb_in=(/3,4,5,6,7,8/) !list of valence band
	bndcb_in=(/9,10,11,12,13,14,15,16/) !list of conduction band


	TK_in=25.d0 !Temperature (K)
	ndop_in=0.d0! input n-doping concentration as cm-3 
	pdop_in=0.d0 ! input p-doping concentration as cm-3 
	
	!params for the calculation of absorption

	!Cutoff for wavevector=gcutratio*|g1|
	gcutratio=3.8d0

	!Polarization of incoming photon
	pol=(/1.d0,0.d0,0.d0/)

	flagabs=(/1,1,1/) !to test only Fermi function, set flagabs to (/0,0,0/)
	
	
	!Fermi_scan
	nbFc=200 !number of scan for quasi-Fermi level Fc
	nbFv=200 !number of scan for quasi-Fermi level Fv
	
	
	nbninj=100
	ninjdown=15.d0
	ninjup=20.d0
	
	!parameters for gain distribution as function of energy:

	nEmap=140 !Number of energy to calculate gain
	Emapdown=0.1d0 !Lower bound of energy
	Emapup=1.5d0 !Upper bound of energy
	ndensgEph=2 !Number of case for injected density carrier

	allocate (densgEph(ndensgEph))
	densgEph=(/17.d0,18.d0/) !list injected density carrier (in exponent unit)
	
	!Scan of Monkhorst-Pack grid resolution
	nscan=1
	allocate (N1ls(nscan),N2ls(nscan),N3ls(nscan))
	N1ls=(/200/)
	N2ls=(/200/)
	N3ls=(/200/)
	
	!suffix
	suffix="scan_absorption.dat"
	
	allocate (ninjls(nbninj))
	do ii=1,nbninj
		ninjls(ii)=ninjdown+(ninjup-ninjdown)*(ii-1)/(1.d0*(nbninj-1))
	end do
	
	allocate (Emap(nEmap))
	do ii=1,nEmap
		Emap(ii)=Emapdown+(Emapup-Emapdown)*(ii-1)/(1.d0*(nEmap-1))
	end do
	

	call init_params_EPM(muSO_in,zetaSO_in,Vloc_in)
	
	call init_geo_EPM(a0_in,rbasis_in,pos1_in,pos2_in)
	

			
	!******************************************************
	
	call init_sep(nr_in,broad_in,nbvb_in,nbcb_in,bndvb_in,bndcb_in,TK_in,ndop_in,pdop_in)

	!offset the band to highest VB at G
	kG=(/0.d0,0.d0,0.d0/)
	!if (rank==root) then
	call pw_generation(kG,gcutratio)
	allocate (EGarr(2))
	call Hfull_solver(2,(/8,9/),EGarr,EvGarr,.false.)
	bndoff_in=EGarr(1)
	
	!181107: Ggap is now module variable
	Ggap=EGarr(2)-EGarr(1)
	
	call set_bandoffset(bndoff_in)
	!write (stdout, '("bndoff: ",f10.5)')bndoff
	
	deallocate(EGarr,outkbas)

	!gainninj scan (energy position of which we scan gain as function of injected density carrier)
	nEls=1
	allocate (Els(nEls))
	Els=(/Ggap+0.06d0,Ggap+0.12d0 /)
	
	!Fermi_scan (lower bounds and upper bounds for Fc list and Fv list)
	Fcbds=(/0.d0,Ggap+0.6d0/)
	Fvbds=(/-0.6d0,Ggap/)
	
	allocate(Fcls(nbFc),Fvls(nbFv))
	
	do ii=1,nbFc
		Fcls(ii)=Fcbds(1)+(Fcbds(2)-Fcbds(1))*(ii-1)/(1.d0*(nbFc-1))
	end do

	do ii=1,nbFv
		Fvls(ii)=Fvbds(1)+(Fvbds(2)-Fvbds(1))*(ii-1)/(1.d0*(nbFv-1))
	end do
	
	!only Fermi scan and gain distribution as function of energy
	call convergence_scan(nscan,N1ls,N2ls,N3ls,rank,nproc,root,gcutratio,pol,&
													&nEmap,Emap,nEls,Els,flagabs,&
													&nbFc,nbFv,Fcls,Fvls,nbninj,ninjls,ndensgEph,densgEph,&
													&suffix,1,0) 
	
	!Activate this call of convergence_scan and deactivate the previous if we want to scan gain as function of injected density carrier

	!call convergence_scan(nscan,N1ls,N2ls,N3ls,rank,nproc,root,gcutratio,pol,&
	!												&nEmap,Emap,nEls,Els,flagabs,&
	!												&nbFc,nbFv,Fcls,Fvls,nbninj,ninjls,ndensgEph,densgEph,&
	!												&suffix,1,1) 
	
	if (rank==root) then
		write(stdout,'("----------------Scan done-----------------")')
	end if	

	call MPI_FINALIZE(ierr)

end program driver_abs
