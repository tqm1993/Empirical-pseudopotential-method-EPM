program driver_fit
	
	use phys_cte
	use setup_epm
	use fit_routine

	
	implicit none
	
	!list of parameters
	!geometrical parameters
	real(dp) ::a0_i,rbasis_i(3,3),pos1_i(3),pos2_i(3)
	!EPM parameters
	real(dp) ::muSO_i,zetaSO_i,Vloc_i(11),ryd_unit,gcutratio
	real(dp) ::xval(10),lbdsfit(10),ubdsfit(10),c11_i,c12_i
	real(dp) :: aG_i,aL_i,Ecib_relax(3),Ewei_i(3,3),epstar_i

	
	!DECLARE THE PARAMETERS
	!geometrical (pos1 and pos2 follow EPM convention)
	
	a0_i=5.658d0 !Ge pure

	rbasis_i(:,1)=a0_i*(/0.d0,0.5d0,0.5d0/)
	rbasis_i(:,2)=a0_i*(/0.5d0,0.d0,0.5d0/)
	rbasis_i(:,3)=a0_i*(/0.5d0,0.5d0,0.d0/)
	pos1_i=(0.d0)*(rbasis_i(:,1)+rbasis_i(:,2)+rbasis_i(:,3))
	pos2_i=(1.d0/4.d0)*(rbasis_i(:,1)+rbasis_i(:,2)+rbasis_i(:,3))
	ryd_unit=13.605693d0 !eV
	
	!Cutoff for wavevector=gcutratio*|g1|
	gcutratio=3.8d0

	!muSO_i,zetaSO_i,Vloc_i will be used to initialize init_params_EPM (it will not affected the fit routine, you can keep these values)
	muSO_i=ryd_unit*0.00045d0 !eV
	zetaSO_i=10.09
	Vloc_i=ryd_unit*(/0.d0,-0.28d0,-0.27d0,-0.26d0,0.04d0,0.05d0,0.06d0,0.02d0,0.0163038d0,0.005d0,0.d0/) 
	

	Ecib_relax=(/0.742d0,0.89d0,0.29d0/) !Fit target for Ge (EL,EG,ESO)
	!Deformation potential for direct gap and indirect gap
	aG_i=-9.48d0 !eV
	aL_i=-2.78d0 !eV
	
	!C11 and C12
	c11_i=129.d0 !GPa
	c12_i=48.d0 !GPa
	
	Ewei_i(1,:)=(/2.d0,2.d0,0.8d0/) !weight for relax case (EL,EG,ESO)
	Ewei_i(2,:)=(/1.d0,1.d0,0.d0/)  !weight for hydrocompression case (EL,EG,ESO)
	Ewei_i(3,:)=(/1.d0,12.d0,0.d0/) !weight for hydrotension case (EL,EG,ESO)
	
	
	!Lower bounds and upper bounds of fit parameters
	dvarfit=ryd_unit*1.d-4
	lbdsfit=ryd_unit*(/-0.3d0,-0.3d0,-0.3d0,-0.03d0,-0.03d0,-0.03d0,-0.03d0,-0.03d0,-0.03d0,0.0001d0/)
	ubdsfit=ryd_unit*(/-0.22d0,-0.22d0,-0.22d0,0.08d0,0.08d0,0.08d0,0.08d0,0.08d0,0.08d0,0.001d0/)
	

	epscib=0.0125d0  !REMEMBER TO UPDATE EPSCIB IF FIT VALUE IS DIFFERENT THAN DEFAULT VALUE
	!initialize the parameters
	
	call init_geo_EPM(a0_i,rbasis_i,pos1_i,pos2_i)
	call init_params_EPM(muSO_i,zetaSO_i,Vloc_i)
	
	call load_cible(aG_i,aL_i,Ecib_relax,Ewei_i)
	
	call bfgs_fit_strain(xval,gcutratio,lbdsfit,ubdsfit)
	

end program driver_fit

















