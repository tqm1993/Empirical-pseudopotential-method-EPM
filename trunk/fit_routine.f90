!This module contains fitting routine for EPM process, with input of direct, indirect gap and SO gap
!More fit option will be added in the future

module fit_routine

	use phys_cte
	use setup_epm
	use spline
	
	implicit none
	
	real(dp) ::aG_cib,aL_cib,c11fit,c12fit
	real(dp) ::Ecibstr(3,3),Ewei(3,3)
	real(dp) ::Imat(3,3),epsmat(3,3)
	real(dp) ::rbas_copy(3,3),epstar

	!fitting
	real(dp) ::dvarfit(11),dbndfit(5)
	
	contains
	
	!Loading reference data of energy gap
	!aG_i,aL_i(eV): Deformation potential of the direct gap at Gamma and the indirect gap at L
	!Ecib_relax(3) (eV): Indirect gap at L, Direct gap at Gamma and SO gap for the relaxed material
	!Ewei_i(3,3):Weight of the gap in the fit function, in the relaxed case, hydrocompression and hydrotension	

	subroutine load_cible(aG_i,aL_i,Ecib_relax,Ewei_i)
		
		implicit none
		real(dp),intent(in) ::Ecib_relax(3),aG_i,aL_i
		real(dp),intent(in) ::Ewei_i(3,3)
		
		Ecibstr(1,:)=Ecib_relax
		aG_cib=aG_i
		aL_cib=aL_i
		Ewei=Ewei_i
		!hydrostatic compression

		Ecibstr(2,:)=Ecib_relax-epscib*3.d0*(/aL_cib,aG_cib,0.d0/)
		!hydrostatic strain
		Ecibstr(3,:)=Ecib_relax+epscib*3.d0*(/aL_cib,aG_cib,0.d0/)
		
				
		Imat(:,1)=(/1.d0,0.d0,0.d0/)
		Imat(:,2)=(/0.d0,1.d0,0.d0/)
		Imat(:,3)=(/0.d0,0.d0,1.d0/)
		
		epsmat=Imat
		
		!epsmat=epscib*Imat
		rbas_copy=rbasis
		
		
	end subroutine load_cible
	
	

	!function used in the fit process
	!x(10): List of trial parameters. The first 9 parameters belong to the local potential (see definition in setup_epm module),
	!the last parameter in the SO coefficient muSO
	!gcutratio: Defining the wavevector cutoff
	!f: Function value
	!df(10): First derivativeof the function, calculated at x
	!Egap(3,3): Gap energy (Gamma, L, SO) calculated at x

	subroutine ffit_strain(x,gcutratio,f,df,Egap)
		
		implicit none
		real(dp),intent(in) ::x(10),gcutratio
		real(dp),intent(out) ::f,df(10),Egap(3,3)
		
		!local variable
		real(dp) ::fbf(2,10),xx(10)
		real(dp) ::EcL,EcG,EvG,EsoG,rbastmp(3,3)
		real(dp) ::tmpk(3),EarrL(1),EarrG(4),pos1tmp(3),pos2tmp(3),ryd_to_eV
		integer  ::bndlsL(1),bndlsG(4)

		complex(dp),allocatable ::Evarr(:,:)
		integer ::ii,jj,aa
		ryd_to_eV=13.605693d0
		bndlsL=(/9/)
		bndlsG=(/3,5,7,9/)
		f=0.d0
		fbf=0.d0
		
		do aa=1,3
		
			if (aa==1) then
				rbastmp=rbas_copy
			else if (aa==2) then
				rbastmp=matmul(Imat-epscib*epsmat,rbas_copy) 
			else if (aa==3) then
				rbastmp=matmul(Imat+epscib*epsmat,rbas_copy)
			end if
			pos1tmp=(0.d0)*(rbastmp(:,1)+rbastmp(:,2)+rbastmp(:,3))
			pos2tmp=(1.d0/4.d0)*(rbastmp(:,1)+rbastmp(:,2)+rbastmp(:,3))
			call init_geo_EPM(a0param,rbastmp,pos1tmp,pos2tmp)
			
			
			
			
			muSO=x(10)
			Vloc=(/0.d0,x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),0.d0/)
			call second_deriv(x_fit,Vloc,11,dy2j,(/1,1/))
			!L pts
			tmpk=0.5d0*gbasis(:,1)+0.5d0*gbasis(:,2)+0.5d0*gbasis(:,3)
			call pw_generation(tmpk,gcutratio)
			call Hfull_solver(1,bndlsL,EarrL,Evarr,.false.)
			call dealloc_kbasis()
			
			!Gpts
			tmpk=0.d0
			call pw_generation(tmpk,gcutratio)
			call Hfull_solver(4,bndlsG,EarrG,Evarr,.false.)
			call dealloc_kbasis()
			
			f=f+Ewei(aa,1)*((EarrL(1)-EarrG(3)-Ecibstr(aa,1))**2.d0)+&
					&Ewei(aa,2)*((EarrG(4)-EarrG(3)-Ecibstr(aa,2))**2.d0)+&
					&Ewei(aa,3)*((EarrG(3)-EarrG(1)-Ecibstr(aa,3))**2.d0)			
					
			Egap(aa,:)=(/EarrL(1)-EarrG(3),EarrG(4)-EarrG(3),EarrG(3)-EarrG(1)/)
	
			do ii=1,10
				
				do jj=1,2
					xx=x
					if (jj==1) then
						xx(ii)=x(ii)-dvarfit(ii)
					else if (jj==2) then
						xx(ii)=x(ii)+dvarfit(ii)
					end if
					muSO=xx(10)
					Vloc=(/0.d0,xx(1),xx(2),xx(3),xx(4),xx(5),xx(6),xx(7),xx(8),xx(9),0.d0/)
					call second_deriv(x_fit,Vloc,11,dy2j,(/1,1/))
					!L pts
					tmpk=0.5d0*gbasis(:,1)+0.5d0*gbasis(:,2)+0.5d0*gbasis(:,3)
					call pw_generation(tmpk,gcutratio)
					call Hfull_solver(1,bndlsL,EarrL,Evarr,.false.)
					call dealloc_kbasis()
					
					!Gpts
					tmpk=0.d0
					call pw_generation(tmpk,gcutratio)
					call Hfull_solver(4,bndlsG,EarrG,Evarr,.false.)
					call dealloc_kbasis()
					
					fbf(jj,ii)=fbf(jj,ii)+Ewei(aa,1)*((EarrL(1)-EarrG(3)-Ecibstr(aa,1))**2.d0)+&
										&Ewei(aa,2)*((EarrG(4)-EarrG(3)-Ecibstr(aa,2))**2.d0)+&
										&Ewei(aa,3)*((EarrG(3)-EarrG(1)-Ecibstr(aa,3))**2.d0)
					
					
					
				end do
				
				df(ii)=(fbf(2,ii)-fbf(1,ii))/(2.d0*dvarfit(ii))
				!print *,"ii,df(ii)",ii,df(ii)
			end do
		
		end do

	end subroutine ffit_strain
	

	!Fit process, using BFGS algorithm
	!x(10): the first 9 parameters belong to the local potential (see definition in setup_epm module),
	!the last parameter in the SO coefficient muSO
	!gcutratio: Defining the wavevector cutoff
	!lbdsfit(10),ubdsfit(10): Lower and upper bounds for each parameter

	subroutine bfgs_fit_strain (x,gcutratio,lbdsfit,ubdsfit)
	
		implicit none
		real(dp),intent(in) ::gcutratio,lbdsfit(10),ubdsfit(10)
		real(dp),intent(inout) ::x(10)
		
		!local variable
		integer :: iprint=3,n=10
		character(len=60) :: task,csave
		integer ::isave(44)
		logical ::lsave(4)
		real(dp) ::dsave(29),f,Egap(3,3)
		integer,  allocatable  :: nbd(:), iwa(:)
		real(dp), allocatable  :: l(:), u(:), df(:), wa(:)
		real(dp) ::ryd_to_eV
		allocate ( nbd(n), l(n), u(n), df(n))
		allocate ( iwa(3*n) )
		allocate ( wa(2*m_BFGS_set*n + 5*n + 11*m_BFGS_set*m_BFGS_set + 8*m_BFGS_set) )
		
		ryd_to_eV=13.605693d0
		nbd=2
		l=lbdsfit
		u=ubdsfit
		x=0.5d0*(lbdsfit+ubdsfit)
		
		task='START'
		!L-BFGS-B call.
		do while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or.task.eq.'START') 
			
			call setulb (n,m_BFGS_set,x,l,u,nbd,f,df,factr_set,pgtol_set,wa,iwa,task,iprint,csave,lsave,isave,dsave)
			
			if (task(1:2) .eq. 'FG') then
					!print *,"toto"
					call ffit_strain(x,gcutratio,f,df,Egap) !x is the exponent of the trial density
					
					!verbose
					write(stdout,'("----------------------------------------------------")')
					write(stdout,'("----------------------------------------------------")')
					write(stdout,'("Vloc(s3): ",3f16.8)')x(1)/ryd_to_eV,x(2)/ryd_to_eV,x(3)/ryd_to_eV
					write(stdout,'("Vloc(s8): ",3f16.8)')x(4)/ryd_to_eV,x(5)/ryd_to_eV,x(6)/ryd_to_eV
					write(stdout,'("Vloc(s11): ",3f16.8)')x(7)/ryd_to_eV,x(8)/ryd_to_eV,x(9)/ryd_to_eV
					write(stdout,'("muSO: ",f16.8)')x(10)/ryd_to_eV
					
					write(stdout,'("*********")')
					write(stdout,'("Case relaxed: ")')
					write(stdout,'("EcL: ",f14.7," | Target: ",f14.7)')Egap(1,1),Ecibstr(1,1)
					write(stdout,'("EcG: ",f14.7," | Target: ",f14.7)')Egap(1,2),Ecibstr(1,2)
					write(stdout,'("EsoG: ",f14.7," | Target: ",f14.7)')Egap(1,3),Ecibstr(1,3)
					
					write(stdout,'("*********")')
					write(stdout,'("Case hydro compression: ")')
					write(stdout,'("EcL: ",f14.7," | Target: ",f14.7)')Egap(2,1),Ecibstr(2,1)
					write(stdout,'("EcG: ",f14.7," | Target: ",f14.7)')Egap(2,2),Ecibstr(2,2)
					write(stdout,'("EsoG: ",f14.7," | Target: ",f14.7)')Egap(2,3),Ecibstr(2,3)
					
					write(stdout,'("*********")')
					write(stdout,'("Case hydro strain: ")')
					write(stdout,'("EcL: ",f14.7," | Target: ",f14.7)')Egap(3,1),Ecibstr(3,1)
					write(stdout,'("EcG: ",f14.7," | Target: ",f14.7)')Egap(3,2),Ecibstr(3,2)
					write(stdout,'("EsoG: ",f14.7," | Target: ",f14.7)')Egap(3,3),Ecibstr(3,3)
					
					write(stdout,'("f residu: ",f14.7)')f

			end if
			
		end do
		
		deallocate (nbd,l,u,df,iwa,wa)
		
	end subroutine bfgs_fit_strain
	

end module fit_routine
