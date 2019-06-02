!This module calculates the absorption using formula from Fermi golden rule
!Parallelization with MPI is included

module abs_routine
	
	use mpi 
	use phys_cte
	use setup_epm

	
	implicit none

	real(dp) ::nr,dVbz,c_alpha,broad,TK   
	integer ::N1,N2,N3,Ntotal ! BZ sampling with Monkhorst-Pack scheme
	integer ::nbvb,nbcb
	integer,allocatable ::bndvb(:),bndcb(:)
	real(dp),allocatable ::Evbls(:,:),Ecbls(:,:)
	real(dp) ::ninj,ndop,pdop
	real(dp) ::bndoff,sigmag
	real(dp),allocatable ::NNNls(:,:),NNNraw(:,:)
	real(dp) ::Ggap
	real(dp) ::offsetCB=0.d0
	
	real(dp) ::dEheav=1.d-3
	real(dp) ::dEfermi=2.d-3,Teps=0.01d0
	real(dp) :: pgtol=1.d-5,factr=1.d+1
	integer :: m_BFGS=5
	logical ::have_spin=.false.
	

	contains
	
	!linear interpolation for gain value
	!Nbpts: number of points in reference data
	!inx(Nbpts): Reference array of x in increasing order
	!inval(Nbpts): Reference array of y
	!x: input value of x to find interpolated value of y

	real(dp) function inter1p(Nbpts,inx,inval,x)
		
		implicit none
		integer,intent(in) ::Nbpts
		real(dp),intent(in) ::inx(Nbpts),inval(Nbpts)
		real(dp),intent(in) ::x
		
		!local variables
		integer ::ii
		
		if ((x < inx(1)) .or. (x > inx(Nbpts))) then
			write (stdout,'("Warning: Density value is not found in input range.&
					&inter1p value will be set to 0")')

			inter1p=0.d0
		
		!check exit statement
		else
			do ii=1,Nbpts
				if ((x >=inx(ii)).and.(x <=inx(ii+1))) then
					inter1p=inval(ii)+(inval(ii+1)-inval(ii))*(x-inx(ii))/&
							&(inx(ii+1)-inx(ii))
					exit
				end if
			end do
			
		end if
	
	end function inter1p
	
	!Read Fermi levels and carrier densities from file
	!infile: Filename
	!fmtout: Line format to read data from file
	!Nbpts: Number of value will be read from data for each data column
	!densls(Nbpts) : Output array of carrier densities
	!Efls(Nbpts) :Output array of Fermi levels
	
	subroutine read_fermi_inter1p(infile,fmtout,Nbpts,densls,Efls)
		
		implicit none
		integer,intent(in) ::Nbpts
		character(len=120),intent(in) ::infile,fmtout
		real(dp),intent(out) ::densls(Nbpts),Efls(Nbpts)
		
		!local variables
		integer ::ii,ierr
		real(dp) ::dummy
		open(1,file=infile,action="read")
		
		do ii=1,Nbpts
			read(1,fmtout,iostat=ierr)Efls(ii),dummy,densls(ii)
			!print *,dummy
		end do
		
		close(1)
		
	end subroutine read_fermi_inter1p
	
	
	!Parameters initialization
	!nr_in: refractive index of materials
	!broad_in(eV): Broadening of the Gaussian function, used to replace the Dirac distribution in the absorption formula
	!nbvb_in: Number of valence bands
	!nbcb_in: Number of conduction bands
	!bndvb_in(nbvb_in): List of valence bands
	!bndcb_in(nbcb_in): List of conduction bands	
	!TK_in(K): Temperature
	!ndop_in(cm-3): Donor doping density
	!pdop_in(cm-3): Acceptor doping density

	subroutine init_sep(nr_in,broad_in,nbvb_in,nbcb_in,bndvb_in,bndcb_in,TK_in,ndop_in,pdop_in)
		
		implicit none
		real(dp),intent(in) :: nr_in,broad_in
		integer,intent(in) :: nbvb_in,nbcb_in
		integer,intent(in) :: bndvb_in(nbvb_in),bndcb_in(nbcb_in)
		real(dp),intent(in) ::TK_in,ndop_in,pdop_in
		
		nr=nr_in
		broad=broad_in
		ndop=ndop_in
		pdop=pdop_in
		nbcb=nbcb_in
		nbvb=nbvb_in
		allocate(bndvb(nbvb),bndcb(nbcb))
		bndvb=bndvb_in
		bndcb=bndcb_in
		TK=TK_in
		sigmag=broad/(2.d0*sqrt(2.d0*log(2.d0)))
	
	end subroutine init_sep
	
	!Mesh initialization (Monkhorst-Pack mesh)
	!rank: MPI processor rank
	!root: MPI processor root
	!N1_in,N2_in,N3_in: number of division in each basis vector g1,g2,g3
	subroutine init_sepN1N2N3(rank,root,N1_in,N2_in,N3_in)
		
		implicit none
		integer,intent(in) ::rank,root
		integer,intent(in) :: N1_in,N2_in,N3_in
		
		!local variable
		integer ::cnt,tt,uu,vv
		real(dp) ::expo,tmpk(3),tmpkred(3)
		
		N1=N1_in
		N2=N2_in
		N3=N3_in
		Ntotal=N1*N2*N3
		dVbz=Vcell_g/(1.d0*N1*N2*N3*1.d0)
		expo=1.d0*(50-34*3+31*2-2) 
		c_alpha=pi*(hbarnorm**3.d0)*dVbz*(10.d0**expo)/(nr*clight*eps0*(me0norm**2.d0)*8.d0*(pi**3.d0))
		
		allocate (NNNls(Ntotal,3),NNNraw(Ntotal,3))
		
		if (rank==root) then
			allocate (Evbls(nbvb,Ntotal),Ecbls(nbcb,Ntotal))
			!nr=2.d0
			write (stdout,'("hello",i5,e12.5)')rank
			write (stdout,'("sigmag",f14.6)')sigmag
		end if
		
		!form list of k points for mpi process
		cnt=1
		do tt=1,N1
			do uu=1,N2
				do vv=1,N3

					tmpkred=(/1.d0*(2*tt-N1-1)/(2.d0*N1),1.d0*(2*uu-N2-1)/(2.d0*N2),1.d0*(2*vv-N3-1)/(2.d0*N3)/)
					tmpk=tmpkred(1)*gbasis(:,1)+tmpkred(2)*gbasis(:,2)+tmpkred(3)*gbasis(:,3)
					NNNls(cnt,:)=tmpk
					NNNraw(cnt,:)=tmpkred
					cnt=cnt+1		
				end do
			end do
		end do
	
	end subroutine init_sepN1N2N3
	
	!dealloc all array related to N1,N2,N3
	subroutine dealloc_N1N2N3(rank,root)
		
		implicit none
		integer,intent(in) ::rank,root
		
		deallocate (NNNls,NNNraw)
		if (rank==root) then
			deallocate (Evbls,Ecbls)

		end if
		
	end subroutine dealloc_N1N2N3

	!set band offset (used most of the time to set the highest valence band to 0 eV)
	subroutine set_bandoffset(bndoff_in)
		
		implicit none
		real(dp),intent(in) ::bndoff_in
		bndoff=bndoff_in

	end subroutine set_bandoffset
	
	!deallocating band data
	subroutine dealloc_list()
		
		deallocate(bndvb,bndcb,Evbls,Ecbls)
	end subroutine dealloc_list
	
	!********************************************************
	!calculating electron density  for defined value of conduction band quasi-Fermi level 
	!Fc(eV): conduction band quasi-Fermi level

	real(dp) function cb_density(Fc)
		
		implicit none
		real(dp),intent(in) ::Fc
		
		!local var
		real(dp) :: val
		real(dp) ::tmpE
		integer :: i,j,cnt

		cnt=0
		val=0.d0
		do j=1,nbcb
		
			do i=1,Ntotal
				tmpE=Ecbls(j,i)
				if (TK>Teps) then
					val=val+Atomet*dVbz/((8.d0*(pi**3.d0))*(1.d0+exp((tmpE-Fc)/(kb*TK))))
					
				else
					if (tmpE>Fc)then
						val=val+Atomet*dVbz/(8.d0*(pi**3.d0))
						print *,"toto"
					end if
				end if
			end do
			
		end do
		
		if (have_spin) then
			cb_density=2.d0*val
		else
			cb_density=val
		end if
		
	end function cb_density
	
	!calculating hole density  for defined value of valence band quasi-Fermi level 
	!Fc(eV): valence band quasi-Fermi level

	real(dp) function vb_density(Fv)
		
		implicit none
		real(dp),intent(in) ::Fv
		
		!local var
		real(dp) :: val
		real(dp) ::tmpE
		integer :: i,j
		
		val=0.d0
		do j=1,nbvb
		
			do i=1,Ntotal
				tmpE=Evbls(j,i)
				if (TK>Teps) then
					val=val+Atomet*dVbz/((8.d0*(pi**3.d0))*(1.d0+exp((Fv-tmpE)/(kb*TK))))
					!print *,"toto"
				else
					if (tmpE<Fv)then
						val=val+Atomet*dVbz/(8.d0*(pi**3.d0))
						!print *,"toto"
					end if
				end if
				
			end do
		end do
		
		if (have_spin) then
			vb_density=2.d0*val
		else
			vb_density=val
		end if
	end function vb_density
	

	
	!********************************************************
	!Test subroutine for the calculation of optical matrix
	!Calculating optical matrix for a list of k-points
	!outfile: Output file name
	!gcutratio: Defined the cutoff for wave vector, see setup_epm module
	!pol(3): Polarization of incoming photon
	!nbkpts: Number of k-points for which the optical matrix will be calculated
	!kls(nbkpts,3): List of k-points
	!nstate: Number of electronic states. Optical matrix will be calculated for transitions between these states
	!ab_state(nstate): Output results

	subroutine test_optmat_k(outfile,gcutratio,pol,nbkpts,kls,nstate,ab_state)
		
		implicit none
		character(len=120),intent(in) ::outfile
		real(dp),intent(in) ::gcutratio,pol(3)
		integer,intent(in) ::nbkpts,nstate
		integer,intent(in) ::ab_state(nstate)
		real(dp),intent(in) ::kls(nbkpts,3)
		
		!local variable
		integer ::ii,aa,bb,jj,ndim,nblap
		complex(dp),allocatable::diag_op(:),ab_lap(:),Evarr(:,:)
		complex(dp) ::zdotc
		logical ::exist
		real(dp) ::tmpk(3),tmpdist,dk(3)
		real(dp),allocatable ::tmpval(:),Earr(:)
		
		inquire(file=outfile,exist=exist)
		if (exist) then
			open(1, file=outfile, status="replace", action="write",form="formatted")
		else
			open(1, file=outfile, status="new", action="write",form="formatted")

		end if	
		
		nblap=(nstate/2)-1
		!******************************************

		!main code
		allocate (tmpval(nblap),Earr(nstate))
		tmpdist=0.d0
		do ii=1,nbkpts
			tmpval=0.d0
			tmpk=kls(ii,:)
			call pw_generation(tmpk,gcutratio)
			call Hfull_solver(nstate,ab_state,Earr,Evarr,.true.)
			ndim=2*npwave
			allocate(diag_op(ndim),ab_lap(ndim))
			do aa=1,npwave
				diag_op(aa)=dcmplx(dot_product(pol,outkbas(aa,4:6)),0.d0)
			end do
				diag_op((npwave+1):2*npwave)=diag_op(1:npwave)

			
			do aa=1,2
				do bb=1,nblap
					
					ab_lap=dconjg(Evarr(aa,:))*Evarr(bb*2+1,:)
					tmpval(bb)=tmpval(bb)+abs(zdotc(ndim,diag_op,1,ab_lap,1))**2.d0
					ab_lap=dconjg(Evarr(aa,:))*Evarr(bb*2+2,:)
					tmpval(bb)=tmpval(bb)+abs(zdotc(ndim,diag_op,1,ab_lap,1))**2.d0
					
				end do
			end do			
			
			if (ii==1) then
				tmpdist=0.d0
			else
				dk=tmpk-kls(ii-1,:)
				tmpdist=tmpdist+sqrt(dot_product(dk,dk))
			end if
			write(1,'(100f20.6)') tmpdist,(tmpval(bb),bb=1,nblap),(Earr(aa),aa=1,nstate)
			
			deallocate(outkbas,Evarr,diag_op,ab_lap)
		end do

		close(1)
		
	end subroutine test_optmat_k

	!Gaussian function

	real(dp) function gaussian(deltaE)
		
		implicit none
		real(dp),intent(in):: deltaE
		
		gaussian=(1.d0/(sigmag*sqrt(2.d0*pi)))*exp(-(deltaE**2.d0)/(2.d0*(sigmag**2.d0)))
	end function gaussian
	
	
	!Calculation of absorption, with scan over parameters
	subroutine convergence_scan(nscan,N1ls,N2ls,N3ls,rank,nproc,root,gcutratio,pol,&
													&nEmap,Emap,nEls,Els,flagabs,&
													&nbFc,nbFv,Fcls,Fvls,nbninj,ninjls,ndensgEph,densgEph,&
													&suffix,flagsep,flag_ninj)
		
		implicit none
		
		integer,intent(in) ::rank,nproc,root,nEmap,nEls,nscan,flagabs(3)
		integer,intent(in) ::nbFc,nbFv,nbninj,flagsep,flag_ninj
		integer,intent(in) ::N1ls(nscan),N2ls(nscan),N3ls(nscan)
		real(dp),intent(in) ::gcutratio,pol(3)
		real(dp),intent(in) ::Emap(nEmap),Els(nEls),Fcls(nbFc),Fvls(nbFv)
		real(dp),intent(in) ::ninjls(nbninj) !in exponent unit
		integer,intent(in) ::ndensgEph
		real(dp),intent(in) ::densgEph(ndensgEph) !in exponent unit
		character(len=120),intent(in) ::suffix
		
		!local variables:
		character(len=120) :: outFermi_Fc,outFermi_Fv,outgain_ninj,outgain_Eph
		character(len=120) ::dummy
		logical ::exist
		integer ::ii,tt
		real(dp) ::Fvlsinv(nbFv)
		
		!band structure variable
		integer ::offcb,bndls(nbvb+nbcb)
		real(dp),allocatable ::Earr(:)
		real(dp),allocatable ::VBVBproc(:,:),VBCBproc(:,:),CBCBproc(:,:)
		integer ::Npart,Nres,Ntotpr
		real(dp),allocatable ::Evbproc(:,:),Ecbproc(:,:)
		integer ::Nelvb,Nelcb,NVBVB,NVBCB,NCBCB,ndim,jj
		integer ::aa,bb,scnt,cnt
		real(dp) ::tmpk(3)
		complex(dp),allocatable ::diag_op(:),ab_lap(:),Evarr(:,:)
		
		!Fermi_level scan
		real(dp) ::tmpdensFc,tmpdensFv
		real(dp) ::cbdenslog(nbFc),vbdenslog(nbFv)
		
		!Gain-ninj scan
		real(dp) ::tmpninj
		real(dp) ::tmpE,Fc(nbninj),Fv(nbninj)
		real(dp) ::coeff,cbdval,vbdval
		real(dp),allocatable ::gain_table(:,:),rgain_table(:,:) !at root
		real(dp) ::valVBVB,valVBCB,valCBCB,valLHCB,valHHCB
		real(dp) ::tmpEa,tmpEb,fer
		
		!Gain-Eph scan
		!190131: add LH,HH-CB gain
		real(dp),allocatable ::gainVBVB(:,:),gainVBCB(:,:),&
							&gainCBCB(:,:),fullgain(:,:),gainLHCB(:,:),gainHHCB(:,:)
		real(dp),allocatable ::rgainVBVB(:,:),rgainVBCB(:,:),&
							&rgainCBCB(:,:),rfullgain(:,:),rgainLHCB(:,:),rgainHHCB(:,:)
		real(dp) ::Fc2(ndensgEph),Fv2(ndensgEph)
		
		!MPI and other
		integer ::ierr
		complex(dp) ::zdotc
		
		!check scan variable
		scanloop:		do tt=1,nscan
			
			call init_sepN1N2N3(rank,root,N1ls(tt),N2ls(tt),N3ls(tt))
			if (rank==root) then
				write(stdout,'("************************************************************")')
				write(stdout,'("************************************************************")')
				write(stdout,'("************************************************************")')
				write(stdout,'("Starting scan ................................")')
				write(stdout,'(A70)')suffix
				write(stdout,'("---------------------------------------")')
				write(stdout,'("scan ",i2," | Ntotal: ",i12)')tt,Ntotal
				write(stdout,'("offsetCB: ",f12.4)')offsetCB
				!initialize all params related to N1,N2,N3 (check carefully)
				
				
				write(dummy,'(A14,I1,A1)')"scan_Fermi_Fc_",tt,"_"
				outFermi_Fc= trim(dummy) //suffix

				write(dummy,'(A14,I1,A1)')"scan_Fermi_Fv_",tt,"_"
				outFermi_Fv= trim(dummy) //suffix
				
				if (flag_ninj==1) then
					write(dummy,'(A15,I1,A1)')"scan_gain_ninj_",tt,"_"
					outgain_ninj= trim(dummy) //suffix
				end if
				
				write(dummy,'(A14,I1,A1)')"scan_gain_Eph_",tt,"_"
				outgain_Eph=trim(dummy) //suffix
				
			end if
			
			!*******************************************
			!Band structure part
			offcb=nbvb
			bndls(1:nbvb)=bndvb
			bndls(offcb+1:offcb+nbcb)=bndcb
			
			if (rank==root) then
				!nbvb and nbcb are always even number
				if (flagabs(1)==1) then
					write (stdout,'("Calculating VB-VB transition")')
				end if
				if (flagabs(2)==1) then
					write (stdout,'("Calculating VB-CB transition")')
				end if
				if (flagabs(3)==1) then
					write (stdout,'("Calculating CB-CB transition")')
				end if
			end if
			
			!simple version (with overhead on last processor)
			!divide work area for each processor
			Npart=Ntotal/nproc
			Nres=mod(Ntotal,nproc)
			
			if (rank .ne. (nproc-1)) then
				Ntotpr=Npart
			else
				Ntotpr=Npart+Nres
			end if
			
			allocate (Evbproc(nbvb,Ntotpr),Ecbproc(nbcb,Ntotpr))
			if (flagabs(1)==1) then
				allocate(VBVBproc((nbvb/2)*(nbvb-1),Ntotpr))
			end if
			if (flagabs(2)==1) then
				allocate(VBCBproc(nbvb*nbcb,Ntotpr))
			end if
			if (flagabs(3)==1) then
				allocate(CBCBproc((nbcb/2)*(nbcb-1),Ntotpr))
			end if
			!print *,"shape ", shape(CBCBproc),rank
			!write (stdout,'("rank: ",i5 " Ntotpr: ",i8)')rank,Ntotpr
			!181027: changed all 5
			Nelvb=Ntotpr*nbvb
			Nelcb=Ntotpr*nbcb
			NVBVB=Ntotpr*(nbvb/2)*(nbvb-1)
			NVBCB=Ntotpr*nbvb*nbcb
			NCBCB=Ntotpr*(nbcb/2)*(nbcb-1)
			!Loop of band structure solver for each processor
			cnt=1
			do ii=1,Ntotpr
				if ((rank==root) .and. (mod(ii,10000)==0)) then
					write (stdout,'("Band str calculation, step: ",i10)'),ii
				end if
				tmpk=NNNls(Npart*rank+ii,:)
				call pw_generation(tmpk,gcutratio)
				allocate(Earr(nbcb+nbvb))
				call Hfull_solver(nbcb+nbvb,bndls,Earr,Evarr,.true.)
				
				Evbproc(:,cnt)=Earr(1:nbvb)
				Ecbproc(:,cnt)=Earr(offcb+1:offcb+nbcb)
				ndim=2*npwave
				
				allocate (diag_op(ndim),ab_lap(ndim))
				do jj=1,npwave
					diag_op(jj)=dcmplx(dot_product(pol,outkbas(jj,4:6)),0.d0)
				end do
				diag_op((npwave+1):2*npwave)=diag_op(1:npwave)
				
				
				!VB-VB transition
				if (flagabs(1)==1) then
					scnt=1
					do aa=1,(nbvb-1)
						do bb=(aa+1),nbvb
							ab_lap=dconjg(Evarr(bb,:))*Evarr(aa,:)
							VBVBproc(scnt,cnt)=abs(zdotc(ndim,diag_op,1,ab_lap,1))**2.d0
							scnt=scnt+1
						end do
					end do

				end if
				
				!VB-CB transition
				if (flagabs(2)==1) then
					scnt=1
					do aa=1,nbvb
						do bb=1,nbcb
							ab_lap=dconjg(Evarr(offcb+bb,:))*Evarr(aa,:)
							VBCBproc(scnt,cnt)=abs(zdotc(ndim,diag_op,1,ab_lap,1))**2.d0
							scnt=scnt+1
						end do
					end do
				end if
					
				!CB-CB transition
				if (flagabs(3)==1) then
					scnt=1
					do aa=1,(nbcb-1)
						do bb=(aa+1),nbcb
							ab_lap=dconjg(Evarr(offcb+bb,:))*Evarr(offcb+aa,:)
							CBCBproc(scnt,cnt)=abs(zdotc(ndim,diag_op,1,ab_lap,1))**2.d0
							scnt=scnt+1
						end do
					end do
				end if
				cnt=cnt+1
				deallocate (outkbas,diag_op,ab_lap,Earr,Evarr)
			end do
			!End of band structure solver for each tmpk
			
			!set HH band as 0
			Evbproc=Evbproc-bndoff
			Ecbproc=Ecbproc-bndoff-offsetCB
			
			!check MPI again to spot if there some arrangement error (double checking with single proc if possible)
			!only true if the number of step divide the number of processor
			!Gather Evbproc,Ecbproc
			call MPI_GATHER(Evbproc,Nelvb,MPI_DOUBLE_PRECISION,Evbls,Nelvb,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
			call MPI_GATHER(Ecbproc,Nelcb,MPI_DOUBLE_PRECISION,Ecbls,Nelcb,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
			!deallocate(Evbproc,Ecbproc)
			
			call MPI_Barrier(MPI_COMM_WORLD,ierr) !block until all processor reaches this routine
			!------------------------------------------------
			!scan Fermi on root
			if (rank==root) then
				inquire(file=outFermi_Fc,exist=exist)
				if (exist) then
					open(1, file=outFermi_Fc, status="replace", action="write",form="formatted")
				else
					open(1, file=outFermi_Fc, status="new", action="write",form="formatted")
				end if
					!print *,Evbls(:,6)
				do ii=1,nbFc
					tmpdensFc=cb_density(Fcls(ii))*1.d-6 !convert to cm-3
				
					write (1,'(f12.4,e14.5,f12.4)'),Fcls(ii),tmpdensFc,log10(tmpdensFc)
					cbdenslog(ii)=log10(tmpdensFc)
				end do
				close(1)
				
				inquire(file=outFermi_Fv,exist=exist)
				if (exist) then
					open(1, file=outFermi_Fv, status="replace", action="write",form="formatted")
				else
					open(1, file=outFermi_Fv, status="new", action="write",form="formatted")
				end if
				
				do ii=1,nbFv
					tmpdensFv=vb_density(Fvls(ii))*1.d-6
					vbdenslog(nbFv-(ii-1))=log10(tmpdensFv) !inverse vbdens array order
					write (1,'(f12.4,e14.5,f12.4)'),Fvls(ii),tmpdensFv,log10(tmpdensFv)
					Fvlsinv(nbFv-(ii-1))=Fvls(ii)
				end do
				close(1)
				
			end if
			call MPI_Bcast(cbdenslog,nbFc,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
			call MPI_Bcast(vbdenslog,nbFv,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
			call MPI_Bcast(Fvlsinv,nbFv,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,ierr)
			call MPI_Barrier(MPI_COMM_WORLD,ierr) !block until all processor reaches this routine
			
			!**************************----------------------
			!Gain-ninj scan
			
			coeff=(clight/nr)*1.d2 !convert to Gp (s-1)
			if (flag_ninj==1) then
				!interpolation
				!BEWARE: ninjls in log10 exponent unit
				do ii=1,nbninj
					tmpninj=10.d0**ninjls(ii)
					Fc(ii)=inter1p(nbFc,cbdenslog,Fcls,log10(tmpninj+ndop))
					Fv(ii)=inter1p(nbFv,vbdenslog,Fvlsinv,log10(tmpninj+pdop))
				end do
				! End of calculation for Fermi level

				allocate(gain_table(nbninj,nEls)) !to be deallocated later
				if (rank==root) then
					allocate (rgain_table(nbninj,nEls))
				end if

				!assembling
				!TODO: check carefully
				!here the gain coefficient with unit of s-1 will be calculated, hence the introduction of coeff
				gainninj_&
				&Els_loop : do jj=1,nEls
					tmpE=Els(jj)
					
					gainninj_&
					&ninj_loop : do ii=1,nbninj
					
						if ((rank==root) .and. (mod(ii,10)==0)) then
							write (stdout,'("Gain-ninj calculation, step: ",i5,i5)'),jj,ii
						end if

						valVBVB=0.d0
						valVBCB=0.d0
						valCBCB=0.d0
						do cnt=1,Ntotpr

							!print *,cnt
							!check fer for each case
							!VB-VB transition
							if (flagabs(1)==1) then
								scnt=1
								do aa=1,(nbvb-1)
									tmpEa=Evbproc(aa,cnt)
									do bb=(aa+1),nbvb
										tmpEb=Evbproc(bb,cnt)
										fer=-(1.d0/(1.d0+exp((Fv(ii)-tmpEa)/(kb*TK))))+(1.d0/(1.d0+exp((Fv(ii)-tmpEb)/(kb*TK))))
										valVBVB=valVBVB+(coeff*c_alpha/tmpE)*VBVBproc(scnt,cnt)*fer*gaussian(tmpEb-tmpEa-tmpE)
										
										scnt=scnt+1
									end do
								end do
							end if
							
							!VB-CB transition
							if (flagabs(2)==1) then
								scnt=1
								do aa=1,nbvb
									tmpEa=Evbproc(aa,cnt)
									do bb=1,nbcb
										tmpEb=Ecbproc(bb,cnt)
										fer=1.d0-(1.d0/(1.d0+exp((Fv(ii)-tmpEa)/(kb*TK))))-(1.d0/(1.d0+exp((tmpEb-Fc(ii))/(kb*TK))))
										valVBCB=valVBCB+(coeff*c_alpha/tmpE)*VBCBproc(scnt,cnt)*fer*gaussian(tmpEb-tmpEa-tmpE)

										scnt=scnt+1
									end do
								end do
							end if

							!CB-CB transition
							if (flagabs(3)==1) then
								scnt=1
								do aa=1,(nbcb-1)
									tmpEa=Ecbproc(aa,cnt)
									do bb=(aa+1),nbcb
										tmpEb=Ecbproc(bb,cnt)
										fer=(1.d0/(1.d0+exp((tmpEa-Fc(ii))/(kb*TK))))-(1.d0/(1.d0+exp((tmpEb-Fc(ii))/(kb*TK))))
										valCBCB=valCBCB+(coeff*c_alpha/tmpE)*CBCBproc(scnt,cnt)*fer*gaussian(tmpEb-tmpEa-tmpE)
														
										scnt=scnt+1
									end do
								end do
							end if

						end do
					!write on table
					gain_table(ii,jj)=valVBVB+valVBCB+valCBCB
					
					end do 	gainninj_ninj_loop
					
				end do gainninj_Els_loop
			
				call MPI_Barrier(MPI_COMM_WORLD,ierr) !wait all processes
				!rgain_table contains the final result
				call MPI_Reduce(gain_table,rgain_table,nbninj*nEls,MPI_DOUBLE_PRECISION,MPI_SUM,root,&
												&MPI_COMM_WORLD,ierr)
				call MPI_Barrier(MPI_COMM_WORLD,ierr)
				deallocate (gain_table)
			
			end if
			
			!**************************----------------------
			!Gain-Eph scan
			!BEWARE: ninjls in log10 exponent unit
			!190131: add separation between LH and HH to CB
			do ii=1,ndensgEph
				
				!print *,densgEph(ii),vbdenslog(1),vbdenslog(nbFv)
				tmpninj=10.d0**densgEph(ii)
				Fc2(ii)=inter1p(nbFc,cbdenslog,Fcls,log10(tmpninj+ndop))
				Fv2(ii)=inter1p(nbFv,vbdenslog,Fvlsinv,log10(tmpninj+pdop))
			end do
			
			if (rank==root) then
				write (stdout,'("Fc list: ",12f10.3)')(Fc2(ii),ii=1,ndensgEph)
				write (stdout,'("Fv list: ",12f10.3)')(Fv2(ii),ii=1,ndensgEph)
			end if
			! End of calculation for Fermi level

			allocate(gainVBVB(ndensgEph,nEmap),gainVBCB(ndensgEph,nEmap),&
							&gainCBCB(ndensgEph,nEmap),fullgain(ndensgEph,nEmap)) !to be deallocated later
							
			if (flagsep==1) then
				allocate (gainLHCB(ndensgEph,nEmap),gainHHCB(ndensgEph,nEmap))
			end if
			
			if (rank==root) then
				allocate(rgainVBVB(ndensgEph,nEmap),rgainVBCB(ndensgEph,nEmap),&
							&rgainCBCB(ndensgEph,nEmap),rfullgain(ndensgEph,nEmap)) !to be deallocated later
				if (flagsep==1) then
					allocate (rgainLHCB(ndensgEph,nEmap),rgainHHCB(ndensgEph,nEmap))
				end if
			end if

			!assembling
			!TODO: check carefully
			!here the gain coefficient with unit of s-1 will be calculated, hence the introduction of coeff
			gainEph_&
			&Eph_loop : do jj=1,nEmap
				tmpE=Emap(jj)
				if ((rank==root) .and. (mod(jj,10)==0)) then
					write (stdout,'("Gain-Eph calculation, step: ",i5)'),jj
				end if
				gainEph_&
				&ninj_loop : do ii=1,ndensgEph

					valVBVB=0.d0
					valVBCB=0.d0
					valCBCB=0.d0
					valLHCB=0.d0
					valHHCB=0.d0
					do cnt=1,Ntotpr
					
						!VB-VB transition
						if (flagabs(1)==1) then
							scnt=1
							do aa=1,(nbvb-1)
								tmpEa=Evbproc(aa,cnt)
								do bb=(aa+1),nbvb
									tmpEb=Evbproc(bb,cnt)
									fer=-(1.d0/(1.d0+exp((Fv2(ii)-tmpEa)/(kb*TK))))+(1.d0/(1.d0+exp((Fv2(ii)-tmpEb)/(kb*TK))))
									valVBVB=valVBVB+(c_alpha/tmpE)*VBVBproc(scnt,cnt)*fer*gaussian(tmpEb-tmpEa-tmpE)
									
									scnt=scnt+1
								end do
							end do
						end if
						
						!VB-CB transition
						if (flagabs(2)==1) then
							scnt=1
							do aa=1,nbvb
								tmpEa=Evbproc(aa,cnt)
								do bb=1,nbcb
									tmpEb=Ecbproc(bb,cnt)
									fer=1.d0-(1.d0/(1.d0+exp((Fv2(ii)-tmpEa)/(kb*TK))))-(1.d0/(1.d0+exp((tmpEb-Fc2(ii))/(kb*TK))))
									valVBCB=valVBCB+(c_alpha/tmpE)*VBCBproc(scnt,cnt)*fer*gaussian(tmpEb-tmpEa-tmpE)
									
									if (flagsep==1) then
										if ((aa==3).or.(aa==4)) then
											valLHCB=valLHCB+(c_alpha/tmpE)*VBCBproc(scnt,cnt)*fer*gaussian(tmpEb-tmpEa-tmpE)
										else if ((aa==5).or.(aa==6)) then
											valHHCB=valHHCB+(c_alpha/tmpE)*VBCBproc(scnt,cnt)*fer*gaussian(tmpEb-tmpEa-tmpE)
										end if	
									end if
									
									scnt=scnt+1
								end do
							end do
						end if

						!CB-CB transition
						if (flagabs(3)==1) then
							scnt=1
							do aa=1,(nbcb-1)
								tmpEa=Ecbproc(aa,cnt)
								do bb=(aa+1),nbcb
									tmpEb=Ecbproc(bb,cnt)
									fer=(1.d0/(1.d0+exp((tmpEa-Fc2(ii))/(kb*TK))))-(1.d0/(1.d0+exp((tmpEb-Fc2(ii))/(kb*TK))))
									valCBCB=valCBCB+(c_alpha/tmpE)*CBCBproc(scnt,cnt)*fer*gaussian(tmpEb-tmpEa-tmpE)
													
									scnt=scnt+1
								end do
							end do
						end if

					end do
					gainVBVB(ii,jj)=valVBVB
					gainVBCB(ii,jj)=valVBCB
					gainCBCB(ii,jj)=valCBCB
					fullgain(ii,jj)=valVBVB+valVBCB+valCBCB
					
					if (flagsep==1) then
						gainLHCB(ii,jj)=valLHCB
						gainHHCB(ii,jj)=valHHCB
					end if
				
				end do 	gainEph_ninj_loop
				
			end do gainEph_Eph_loop
			
			call MPI_Barrier(MPI_COMM_WORLD,ierr) !wait all processes
			call MPI_Reduce(gainVBVB,rgainVBVB,ndensgEph*nEmap,MPI_DOUBLE_PRECISION,MPI_SUM,root,&
											&MPI_COMM_WORLD,ierr)
			call MPI_Reduce(gainVBCB,rgainVBCB,ndensgEph*nEmap,MPI_DOUBLE_PRECISION,MPI_SUM,root,&
											&MPI_COMM_WORLD,ierr)
			call MPI_Reduce(gainCBCB,rgainCBCB,ndensgEph*nEmap,MPI_DOUBLE_PRECISION,MPI_SUM,root,&
											&MPI_COMM_WORLD,ierr)
			call MPI_Reduce(fullgain,rfullgain,ndensgEph*nEmap,MPI_DOUBLE_PRECISION,MPI_SUM,root,&
											&MPI_COMM_WORLD,ierr)
			
			if (flagsep==1) then
				call MPI_Reduce(gainLHCB,rgainLHCB,ndensgEph*nEmap,MPI_DOUBLE_PRECISION,MPI_SUM,root,&
								&MPI_COMM_WORLD,ierr)
				call MPI_Reduce(gainHHCB,rgainHHCB,ndensgEph*nEmap,MPI_DOUBLE_PRECISION,MPI_SUM,root,&
								&MPI_COMM_WORLD,ierr)
			end if
			
			call MPI_Barrier(MPI_COMM_WORLD,ierr)
			
			deallocate (gainVBVB,gainVBCB,gainCBCB,fullgain)
			
			if (flagsep==1) then
				deallocate (gainLHCB,gainHHCB)
			end if
			!write in file
			if (rank==root) then
			
				!write gain-ninj
				
				if (flag_ninj==1) then
					inquire(file=outgain_ninj,exist=exist)
					if (exist) then
						open(1, file=outgain_ninj, status="replace", action="write",form="formatted")
					else
						open(1, file=outgain_ninj, status="new", action="write",form="formatted")
					end if
					
					do ii=1,nbninj
						tmpninj=10.d0**ninjls(ii)
						write(1,'(e14.6,e14.6,e14.6,f12.5,f12.5,80e16.6)')ninjls(ii),&
																								&log10(tmpninj+ndop),log10(tmpninj+pdop),&
																								&Fc(ii),Fv(ii),(rgain_table(ii,jj),jj=1,nEls),&
																								&(rgain_table(ii,jj)/coeff,jj=1,nEls)
						
					end do
					close(1)
					
					deallocate (rgain_table)
				
				end if
				
				!write gain-Eph
				inquire(file=outgain_Eph,exist=exist)
				if (exist) then
					open(1, file=outgain_Eph, status="replace", action="write",form="formatted")
				else
					open(1, file=outgain_Eph, status="new", action="write",form="formatted")
				end if
				
				if (flagsep .ne. 1) then 
					do ii=1,nEmap
						write (1,'(f12.4,200e16.6)')Emap(ii),(rgainVBVB(aa,ii),aa=1,ndensgEph),&
																			&(rgainVBCB(aa,ii),aa=1,ndensgEph),&
																			&(rgainCBCB(aa,ii),aa=1,ndensgEph),&
																			&(rfullgain(aa,ii),aa=1,ndensgEph)
					end do
					
				else
					do ii=1,nEmap
						write (1,'(f12.4,200e16.6)')Emap(ii),(rgainVBVB(aa,ii),aa=1,ndensgEph),&
																			&(rgainVBCB(aa,ii),aa=1,ndensgEph),&
																			&(rgainCBCB(aa,ii),aa=1,ndensgEph),&
																			&(rfullgain(aa,ii),aa=1,ndensgEph),&
																			&(rgainLHCB(aa,ii),aa=1,ndensgEph),&
																			&(rgainHHCB(aa,ii),aa=1,ndensgEph)
					end do
				end if
				
				deallocate (rgainCBCB,rgainVBVB,rgainVBCB,rfullgain)
				
				if (flagsep==1) then
					deallocate (rgainLHCB,rgainHHCB)
				end if
				
			end if
			deallocate (Evbproc,Ecbproc,VBVBproc,CBCBproc,VBCBproc)
			call dealloc_N1N2N3(rank,root)
			call MPI_Barrier(MPI_COMM_WORLD,ierr)
		end do scanloop 
		
		
	end subroutine convergence_scan
	
end module abs_routine













