! This module provides setup for atom positions in the crystal and 
! calculate band structure energy using local potential Empirical Pseudopotential Method (EPM)

module setup_epm

	use phys_cte
	use spline

	
	implicit none

	real(dp) ::gcuttol=1.0001d0
	real(dp) ::tolrelax=0.0001d0
	real(dp),allocatable :: outkbas(:,:)
	integer ::npwave

	real(dp) ::epscib=0.01  
	!computational params
	complex(dp) ::pauli_x(2,2),pauli_y(2,2),pauli_z(2,2)
	real(dp) ::rbasis(3,3),gbasis(3,3),Vcell_r,Vcell_g
	real(dp) ::pos1(3),pos2(3) !position of atom, in cartesian coordinate
	
	!EPM parameters
	real(dp) ::muSO,zetaSO
	real(dp) ::Vloc(11),x_fit(11) !x_fit for special norms of kpts
	real(dp) ::a0param
	
	!spline parameters 
	integer ::btype(2)
	real(dp) ::dy2j(11)
	
	!BFGS variables 
	real(dp) :: pgtol_set=1.d-5,factr_set=1.d+7
	integer :: m_BFGS_set=5
	
	contains

	!cross product vv= vv1 x vv2 
	subroutine cross_prod(vv1,vv2,vv)
		
		implicit none
		real(dp),intent(in) ::vv1(3),vv2(3)
		real(dp),intent(inout) :: vv(3)
		
		vv(1)=vv1(2)*vv2(3)-vv1(3)*vv2(2)
		vv(2)=vv1(3)*vv2(1)-vv1(1)*vv2(3)
		vv(3)=vv1(1)*vv2(2)-vv1(2)*vv2(1)

	end subroutine cross_prod
	
	!initializing parameters in EPM and Pauli matrices 
	!For a definition of each parameters, see PRB 48, 14276 (1993)	
 
	subroutine init_params_EPM(muSO_i,zetaSO_i,Vloc_i)

		implicit none
		real(dp),intent(in) ::muSO_i,zetaSO_i
		real(dp),intent(in) ::Vloc_i(11)
		
		
		muSO=muSO_i
		zetaSO=zetaSO_i
		Vloc=Vloc_i

		pauli_x(1,:)=(/dcmplx(0.d0,0.d0),dcmplx(1.d0,0.d0)/)
		pauli_x(2,:)=(/dcmplx(1.d0,0.d0),dcmplx(0.d0,0.d0)/)
		pauli_y(1,:)=(/dcmplx(0.d0,0.d0),dcmplx(0.d0,-1.d0)/)
		pauli_y(2,:)=(/dcmplx(0.d0,1.d0),dcmplx(0.d0,0.d0)/)
		pauli_z(1,:)=(/dcmplx(1.d0,0.d0),dcmplx(0.d0,0.d0)/)
		pauli_z(2,:)=(/dcmplx(0.d0,0.d0),dcmplx(-1.d0,0.d0)/)
		
	end subroutine init_params_EPM


	subroutine init_geo_EPM(a0_i,rbasis_i,pos1_i,pos2_i)
		
		implicit none
		real(dp),intent(in) ::rbasis_i(3,3),pos1_i(3),pos2_i(3),a0_i
		
		!local variable
		real(dp) ::r1cr2(3),g1cg2(3)
		real(dp) ::coefffit
		

		
		rbasis=rbasis_i
		call cross_prod(rbasis_i(:,1),rbasis_i(:,2),r1cr2)
		Vcell_r=abs(dot_product(r1cr2,rbasis_i(:,3)))
		call cross_prod(rbasis_i(:,2),rbasis_i(:,3),gbasis(:,1))
		call cross_prod(rbasis_i(:,3),rbasis_i(:,1),gbasis(:,2))
		call cross_prod(rbasis_i(:,1),rbasis_i(:,2),gbasis(:,3))

		gbasis=gbasis*(2.d0*pi/Vcell_r)
		
		call cross_prod(gbasis(:,1),gbasis(:,2),g1cg2)
		Vcell_g=abs(dot_product(g1cg2,gbasis(:,3)))


		pos1=pos1_i
		pos2=pos2_i
		a0param=a0_i
		
		coefffit=2.d0*pi/a0param
		x_fit=coefffit*(/0.d0,s3/(1.d0+epscib),s3,s3/(1.d0-epscib),s8/(1.d0+epscib),s8,s8/(1.d0-epscib),&
									&s11/(1.d0+epscib),s11,s11/(1.d0-epscib),3.8d0/) 
		btype=(/1,1/)
		!extract array of second derivatives dy2j
		call second_deriv(x_fit,Vloc,11,dy2j,btype)
		
	end subroutine init_geo_EPM
	
	!generate plane wave basis (without spin component)
	!input k-vector k in cartesian coordinates
	!gcutratio define the wavevector cutoff

	subroutine pw_generation(k,gcutratio)
		
		implicit none
		real(dp),intent(in)::k(3),gcutratio

		!local variable
		integer ::ii,jj,uu
		integer ::lim1,lim2,lim3,limmax
		real(dp) ::norm1,norm2,norm3,tmpk(3),kmodule,gcut
		
		gcut=gcutratio*2.d0*pi/a0param

		!check carefully gbasis definition
		norm1=sqrt(dot_product(gbasis(:,1),gbasis(:,1)))
		norm2=sqrt(dot_product(gbasis(:,2),gbasis(:,2)))
		norm3=sqrt(dot_product(gbasis(:,3),gbasis(:,3)))
		lim1=nint(gcut/norm1)+2
		lim2=nint(gcut/norm2)+2
		lim3=nint(gcut/norm3)+2
		limmax=max(lim1,lim2,lim3)
		allocate(outkbas(8*limmax*limmax,7))
		npwave=0

		do ii=-limmax,limmax
			do jj=-limmax,limmax
				do uu=-limmax,limmax

					!tmpk in cartesian coordinate
					tmpk=k+1.d0*ii*gbasis(:,1)+1.d0*jj*gbasis(:,2)+1.d0*uu*gbasis(:,3)
					kmodule=sqrt(dot_product(tmpk,tmpk))
					if (kmodule < gcut*gcuttol) then
						npwave=npwave+1
						outkbas(npwave,:)=(/1.d0*ii,1.d0*jj,1.d0*uu,tmpk(1),tmpk(2),tmpk(3),kmodule/)
					end if
					
				end do
			end do
		end do
		
	end subroutine pw_generation
	
	!deallocating the plane wave basis
	subroutine dealloc_kbasis()
		deallocate(outkbas)
	end subroutine dealloc_kbasis
	
	!generating the spin-orbit hamiltonian in EPM
	!basis vector with spin up appear first, then with spin down

	subroutine Hso_gen(hsomat)
		implicit none
		complex(dp),intent(inout) ::hsomat(2*npwave,2*npwave)
		
		!local variable
		integer ::ii,jj
		real(dp) ::kirkj(3),tmpki(3),tmpkj(3),normi,normj
		integer ::spin_i,spin_j
		real(dp) ::tmp2
		complex(dp) ::tmp1,tmp3,coeff,tmpval

		coeff=dcmplx(0.d0,-2.d0*muSO*((Vcell_r/(2.d0))**(2.d0/3.d0))/(pi*pi))
		hsomat=dcmplx(0.d0,0.d0)

		do ii=1,2*npwave
		
			if (ii> npwave) then
				tmpki=outkbas(ii-npwave,4:6)
				normi=outkbas(ii-npwave,7)
				spin_i=2
			else
				tmpki=outkbas(ii,4:6)
				normi=outkbas(ii,7)
				spin_i=1
			end if
			
			do jj=ii,2*npwave

				if (jj> npwave) then
					tmpkj=outkbas(jj-npwave,4:6)
					normj=outkbas(jj-npwave,7)
					spin_j=2
				else
					tmpkj=outkbas(jj,4:6)
					normj=outkbas(jj,7)
					spin_j=1
				end if
				
				call cross_prod(tmpki,tmpkj,kirkj)
				tmp3=kirkj(1)*pauli_x(spin_i,spin_j)+kirkj(2)*pauli_y(spin_i,spin_j)+&
						&kirkj(3)*pauli_z(spin_i,spin_j)
				
				tmp2=(5.d0-(normi/zetaSO)**2.d0)*(5.d0-(normj/zetaSO)**2.d0)/&
							&(25.d0*((1.d0+(normi/zetaSO)**2.d0)**4.d0)*&
							&((1.d0+(normj/zetaSO)**2.d0)**4.d0))
				
				tmp1=0.5d0*(exp(dcmplx(0.d0,-dot_product(tmpki-tmpkj,pos1)))+&
										&exp(dcmplx(0.d0,-dot_product(tmpki-tmpkj,pos2))))
				tmpval=coeff*tmp1*tmp2*tmp3
				hsomat(ii,jj)=tmpval
				
				if (ii .ne. jj) then
					hsomat(jj,ii)=dconjg(tmpval)
				end if

			end do
		end do
		
	end subroutine Hso_gen
	
	!generating the hamiltonian which contains both the contribution from the kinetic part and the
	!local potential part

	subroutine Hloc_gen(hlocmat)
		
		implicit none
		complex(dp),intent(inout) ::hlocmat(2*npwave,2*npwave)
		
		!local variables
		integer ::ii,jj,ierr
		complex(dp) ::sgfact
		real(dp) ::tmpki(3),tmpkj(3),normij,vlocval,coeffkin,normi
		real(dp) ::coeff
		
		hlocmat=dcmplx(0.d0,0.d0)
		coeffkin=(hbar**2.d0)*(1.d20)/(2.d0*me0*ecb) !eV*A02

		do ii=1,npwave
			tmpki=outkbas(ii,4:6)
			normi=outkbas(ii,7)
			do jj=ii,npwave
				tmpkj=outkbas(jj,4:6)
				sgfact=0.5d0*(exp(dcmplx(0.d0,-dot_product(tmpki-tmpkj,pos1)))+&
										&exp(dcmplx(0.d0,-dot_product(tmpki-tmpkj,pos2))))
					
				normij=sqrt(dot_product(tmpki-tmpkj,tmpki-tmpkj))
				vlocval=0.d0

					
				if (normij <x_fit(11)) then
					call spline_interp(x_fit,Vloc,dy2j,11,normij,vlocval)
				end if
				

				
				hlocmat(ii,jj)=vlocval*sgfact

				
				if (ii .ne. jj) then
					hlocmat(jj,ii)=dconjg(hlocmat(ii,jj))
				else 
					hlocmat(ii,jj)=hlocmat(ii,jj)+coeffkin*(normi**2.d0)
				end if
				
			end do
		end do

		hlocmat((npwave+1):2*npwave,(npwave+1):2*npwave)=hlocmat(1:npwave,1:npwave)
	
	end subroutine Hloc_gen
	
	
	!Solving the full hamiltonian H=H SO + H loc + H kinetic
	!nband: number of bands for which we extract the energy
	!bndls(nband): list of bands
	!Earr(nband): list of energy for each band
	!Evarr(nband,2*npwave): list of eigenvector for each band (optional)
	!flagEv: .true. means eigenvectors will be stored, .false. otherwise

	subroutine Hfull_solver(nband,bndls,Earr,Evarr,flagEv)
		
		implicit none
		integer,intent(in) :: nband
		integer,intent(in) ::bndls(nband)
		real(dp),intent(out) ::Earr(nband)
		complex(dp),allocatable,intent(out) ::Evarr(:,:)
		logical,intent(in) ::flagEv
		
		!local variable
		complex(dp) ::hlocmat(2*npwave,2*npwave),hsomat(2*npwave,2*npwave)
		complex(dp),allocatable :: hbulk(:,:)
		integer ::lwork,info,hfulldim,ii
		complex(dp),allocatable ::work(:)
		real(dp),allocatable :: rwork(:)
		real(dp),allocatable ::tmpEres(:)
		
		allocate(hbulk(2*npwave,2*npwave))
		if (flagEv .eqv. .true.) then
			allocate(Evarr(nband,2*npwave))
		end if
		
		call Hloc_gen(hlocmat)
		call Hso_gen(hsomat)
		hbulk=hlocmat+hsomat

		hfulldim=2*npwave
		lwork=2*hfulldim-1
		allocate(tmpEres(hfulldim),work(lwork),rwork(3*hfulldim-2))
		
		call zheev('V','L',hfulldim,hbulk,hfulldim,tmpEres,work,lwork,rwork,info)
		do ii=1,nband
			Earr(ii)=tmpEres(bndls(ii))
			if (flagEv .eqv. .true.) then
				Evarr(ii,:)=hbulk(:,bndls(ii)) !column convention
			end if
		end do
		
		deallocate(hbulk,tmpEres,work,rwork)
	
	end subroutine Hfull_solver

	
	
end module setup_epm











