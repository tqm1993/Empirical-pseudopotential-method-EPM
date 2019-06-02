
module bandplot


	use phys_cte 
	use setup_epm
	
	contains
	
	subroutine band_visualize (nksym,ksymls1,ksymls2,flagsym,nkstep,nband,outfile,gcut)
		
		implicit none
		integer,intent(in) :: nksym,nkstep,nband
		integer,intent(in) :: flagsym(nksym-1)
		real(dp),intent(in) ::ksymls1(nksym,3),ksymls2(nksym-1,3),gcut
		character(len=120),intent(in) :: outfile
		
		!local variables
		integer ::ii,jj,bndls(nband),aa
		real(dp) :: tmpkold(3),kpath,kzero(3),tmpkraw(3),tmpk(3)
		real(dp) ::kstep(nksym-1,3),Earr(nband),kspec(3)
		logical :: exist
		complex(dp),allocatable ::Evarr(:,:)
		
		kpath=0.d0
		tmpkold=ksymls1(1,1)*gbasis(:,1)+ksymls1(1,2)*gbasis(:,2)+ksymls1(1,3)*gbasis(:,3)

		inquire(file=outfile,exist=exist)
		if (exist) then
			open(1, file=outfile, status="replace", action="write",form="formatted")
		else
			open(1, file=outfile, status="new", action="write",form="formatted")

		end if
		do ii=1,nband
			bndls(ii)=ii
		end do
		
		do ii=1,nksym-1
			kstep(ii,:)=(1.d0/(1.d0*nkstep))*(ksymls2(ii,:)-ksymls1(ii,:))
		end do
		
		do ii=1,nksym-1
			kzero=ksymls1(ii,:)
			print *,"step: ",ii
			do jj=0,nkstep-1
				tmpkraw=kzero+1.d0*jj*kstep(ii,:)
				tmpk=tmpkraw(1)*gbasis(:,1)+tmpkraw(2)*gbasis(:,2)+tmpkraw(3)*gbasis(:,3)
				if ((jj==0) .and.(flagsym(ii)==1).and. (ii .ne. 1)) then
					kspec=ksymls2(ii-1,1)*gbasis(:,1)+ksymls2(ii-1,2)*gbasis(:,2)+ksymls2(ii-1,3)*gbasis(:,3)
					kpath=kpath+sqrt(dot_product(kspec-tmpkold,kspec-tmpkold))
				else
					kpath=kpath+sqrt(dot_product(tmpk-tmpkold,tmpk-tmpkold))
				end if	
				!print *,kpath,tmpkraw
				call pw_generation(tmpk,gcut)
				call Hfull_solver(nband,bndls,Earr,Evarr,.false.)
				call dealloc_kbasis()
				write (1,'(100e16.8)')kpath,(Earr(aa),aa=1,nband)
				tmpkold=tmpk
			end do
		end do
		
		close(1)

	end subroutine band_visualize
	
end module bandplot

program driver_scan_strain
	
	use phys_cte
	use setup_epm
	use bandplot
	!use num_fnc

	
	implicit none
	
	!list of parameters
	!geometrical parameters
	real(dp) ::a0_i,rbasis_raw(3,3),pos1_i(3),pos2_i(3),rbasis_i(3,3)
	!EPM parameters
	real(dp) ::muSO_i,zetaSO_i,Vloc_i(11),ryd_unit
	!other parameters
	integer ::ii,jj,uu
	integer ::nksym,nkstep,nband,ierr
	integer,allocatable ::bndls(:),flagsym(:)
	real(dp),allocatable ::ksymls1(:,:),ksymls2(:,:)
	real(dp),allocatable ::Earr(:) 
	complex(dp),allocatable :: Evarr(:,:)
	real(dp) :: gcut,kzero(3),tmpk(3),tmpx,Vloc2(5),tmpkred(3),VBlev
	logical ::flagEv,exist
	real(dp)::kpath,tmpkold(3),deltakred(3),deltak(3),bb(5),cc(5),dd(5),coeff,fff
	complex(dp),allocatable ::hsomat(:,:)
	
	!deformation parameters
	real(dp) ::epsstr,c11,c12,epsmat(3,3),Iden(3,3),epstmp
	integer ::flag,dr_mode,nbscan
	real(dp) ::kGpts(3),kLpts(3),epsinit,epsfinal,EGarr(12),ELarr(12)
	character(len=120) ::outfile,outband
	!logical ::exist
	
	!select mode: 1 for scan mode, 2 for band visualize
	dr_mode=1
	
	!DECLARE THE PARAMETERS
	!flag=1 biaxe [001], flag=2 uniaxe [001]
	!elastic coeff c11,c12
	flag=1
	c11=129.d0
	c12=48.d0
	
	!strain matrix
	Iden(:,1)=(/1.d0,0.d0,0.d0/)
	Iden(:,2)=(/0.d0,1.d0,0.d0/)
	Iden(:,3)=(/0.d0,0.d0,1.d0/)
	
	if (flag==1) then
		epsmat(:,1)=(/1.d0,0.d0,0.d0/)
		epsmat(:,2)=(/0.d0,1.d0,0.d0/)
		epsmat(:,3)=(/0.d0,0.d0,-2.d0*(c12/c11)/)
	else if (flag==2) then
		epsmat(:,1)=(/1.d0,0.d0,0.d0/)
		epsmat(:,2)=(/0.d0,-c12/(c12+c11),0.d0/)
		epsmat(:,3)=(/0.d0,0.d0,-c12/(c12+c11)/)
	end if
	
	
	!geometrical (pos1 and pos2 follow EPM convention)
	
	a0_i=5.658d0! cell parameter, here for example Ge
	
	!real space basis vector in relaxed case
	rbasis_raw(:,1)=a0_i*(/0.d0,0.5d0,0.5d0/)
	rbasis_raw(:,2)=a0_i*(/0.5d0,0.d0,0.5d0/)
	rbasis_raw(:,3)=a0_i*(/0.5d0,0.5d0,0.d0/)
	
	!EPM params
	!example Ge
	ryd_unit=13.605693d0 !eV
	
	!EPM parameters
	muSO_i=ryd_unit*0.00043866d0 !eV
	zetaSO_i=10.09 !A0-1
	Vloc_i=ryd_unit*(/0.d0,-0.26033304d0,-0.25887354d0,-0.25740279d0,&
									&0.03636844d0,0.03991221d0,0.04318571d0,&
									&0.02953685d0,0.03224141d0,0.03521343d0,0.d0/) !example_Ge
	!print *,Vloc_i
	!cutoff
	gcut=5.d0
	flagEv=.false. !.false. means no storage of wavevector
	

	!lower and upper bound for the scan of applied strain
	epsinit=-0.014
	epsfinal=0.014
	
	!Number of scan points
	nbscan=28
	!Output filename
	outfile="scan_strain.dat"
		
	
	call init_params_EPM(muSO_i,zetaSO_i,Vloc_i)	
	!initialize the parameters
	
	if (dr_mode==1) then
		
		inquire(file=outfile,exist=exist)
		if (exist) then
			open(1, file=outfile, status="replace", action="write",form="formatted")
		else
			open(1, file=outfile, status="new", action="write",form="formatted")

		end if
		
		allocate (bndls(12))
		do ii=1,12
			bndls(ii)=ii
		end do
		
		do ii=0,nbscan
			
			epstmp=epsinit+(epsfinal-epsinit)*ii/(1.d0*nbscan)
			rbasis_i=matmul(Iden+epstmp*epsmat,rbasis_raw)
			pos1_i=(0.d0)*(rbasis_i(:,1)+rbasis_i(:,2)+rbasis_i(:,3))
			pos2_i=(1.d0/4.d0)*(rbasis_i(:,1)+rbasis_i(:,2)+rbasis_i(:,3))
			call init_geo_EPM(a0_i,rbasis_i,pos1_i,pos2_i)
			
			kGpts= (/0.d0,0.d0,0.d0/)
			call pw_generation(kGpts,gcut)
			call Hfull_solver(12,bndls,EGarr,Evarr,flagEv)
			call dealloc_kbasis()
			
			kLpts=0.5d0*gbasis(:,1)+0.5d0*gbasis(:,2)+0.5d0*gbasis(:,3)
			call pw_generation(kLpts,gcut)
			call Hfull_solver(12,bndls,ELarr,Evarr,flagEv)
			call dealloc_kbasis()
			
			!epstmp,LH,HH,CB-G,CB-L,G-LH,G-HH,L-LH,L-HH
			write (1,'(9f16.7)')epstmp,EGarr(6),EGarr(8),EGarr(9),ELarr(9),&
													&EGarr(9)-EGarr(6),EGarr(9)-EGarr(8),ELarr(9)-EGarr(6),ELarr(9)-EGarr(8)
		
		end do
		
		close(1)
	else if (dr_mode==2) then
		nksym=5
		nkstep=49
		nband=30
		outband="band_structure.dat"
		epsstr=-0.015
		
		rbasis_i=matmul(Iden+epsstr*epsmat,rbasis_raw)
		pos1_i=(0.d0)*(rbasis_i(:,1)+rbasis_i(:,2)+rbasis_i(:,3))
		pos2_i=(1.d0/4.d0)*(rbasis_i(:,1)+rbasis_i(:,2)+rbasis_i(:,3))
		call init_geo_EPM(a0_i,rbasis_i,pos1_i,pos2_i)	
		allocate (ksymls1(nksym,3),ksymls2(nksym-1,3),flagsym(nksym-1))
		
		!ksymls1
		ksymls1(1,:)=(/0.5d0,0.5d0,0.5d0/)
		ksymls1(2,:)=(/0.d0,0.d0,0.d0/)
		ksymls1(3,:)=(/0.5d0,0.5d0,0.d0/)
		ksymls1(4,:)=(/0.d0,0.d0,0.d0/)
		ksymls1(5,:)=(/0.5d0,0.5d0,0.d0/)
		
		!ksymls2
		ksymls2(1,:)=(/0.d0,0.d0,0.d0/)
		ksymls2(2,:)=(/0.d0,0.5d0,0.5d0/)
		ksymls2(3,:)=(/0.d0,0.d0,0.d0/)
		ksymls2(4,:)=(/0.5d0,0.5d0,0.d0/)
		
		!flagsym
		flagsym=(/0,0,1,0/)
		
		call band_visualize (nksym,ksymls1,ksymls2,flagsym,nkstep,nband,outband,gcut)
	end if
		
		
	


	
	
end program driver_scan_strain

















