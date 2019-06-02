!This module calculate the coefficient for the spline interpolation from a dataset
!and perform the spline interpolation on a given point

module spline

	use phys_cte

	contains
	
	!extract the coefficient for spline interpolation
	!nbdata: Number of data in the reference array
	!xx(nbdata),yy(nbdata): Arrays of x and y points used to calculated the coefficient for the spline interpolation
	!dy2j(nbdata): Output spline coefficient
	!btype(2): Boundary condition at the initial and end points: A value of 0 indicate a zero value for the second derivative and
	!a value of 1 indicate a zero value for the first derivative.	

	subroutine second_deriv(xx,yy,nbdata,dy2j,btype)

		implicit none
		integer,intent(in) ::nbdata,btype(2)
		real(dp),intent(in) ::xx(nbdata),yy(nbdata)
		real(dp),intent(out) ::dy2j(nbdata)
		!local variables
		integer ::ii,info
		real(dp) ::diag(nbdata-2),udiag(nbdata-3),ldiag(nbdata-3),bb(nbdata-2,1)
		
		dy2j(1)=0.d0
		dy2j(nbdata)=0.d0
		
		if (btype(1)==1) then
			dy2j(1)=3.d0*(yy(2)-yy(1))/((xx(2)-xx(1))**2.d0)
		end if
		
		if (btype(2)==1) then
			dy2j(nbdata)=-3.d0*(yy(nbdata)-yy(nbdata-1))/((xx(nbdata)-xx(nbdata-1))**2.d0)
		end if		
		
		do ii=2,nbdata-1
			
			diag(ii-1)=(1.d0/3.d0)*(xx(ii+1)-xx(ii-1))
			bb(ii-1,1)=((yy(ii+1)-yy(ii))/(xx(ii+1)-xx(ii)))-((yy(ii)-yy(ii-1))/(xx(ii)-xx(ii-1)))
			
			if (ii ==2) then
				bb(ii-1,1)=bb(ii-1,1)-(1.d0/6.d0)*(xx(ii)-xx(ii-1))*(dy2j(ii-1))
			end if
			
			if (ii ==(nbdata-1)) then
				bb(ii-1,1)=bb(ii-1,1)-(1.d0/6.d0)*(xx(ii+1)-xx(ii))*(dy2j(ii+1))
			end if
			
			if (ii .ne. 2) then
				ldiag(ii-2)=(1.d0/6.d0)*(xx(ii)-xx(ii-1))
			end if
			
			if (ii .ne. (nbdata-1)) then
				udiag(ii-1)=(1.d0/6.d0)*(xx(ii+1)-xx(ii))
			end if
			
		end do

		call dgtsv(nbdata-2,1,ldiag,diag,udiag,bb,nbdata-2,info)
		
		if (info .ne. 0) then
			write (stdout,'("Error in tridiagonal subroutine. Process halts")')
		end if
		!print *,"bb",bb
		dy2j(2:nbdata-1)= bb(:,1)

	end subroutine second_deriv
	
	!Calculating the spline interpolation at a given point, knowing the spline coefficients dy2j(nbdata)
	!nbdata: Number of points in the array of reference xx and yy
	!xx(nbdata),yy(nbdata):Arrays of x and y points used to calculated the coefficient for the spline interpolation
	!dy2j(nbdata): Array of spline coefficients
	!x: Point to calculate the spline interpolation
	!y: Output spline interpolation 	

	subroutine spline_interp(xx,yy,dy2j,nbdata,x,y)
		
		implicit none
		integer,intent(in) ::nbdata
		real(dp),intent(in) :: xx(nbdata),yy(nbdata),dy2j(nbdata),x
		real(dp),intent(out) ::y
		!local variables
		integer::j,jp1,ii,info
		real(dp) ::AA,BB,CC,DD
		
		info=0
		do ii=1,nbdata-1
			if ((x >= xx(ii)) .and. (x <= xx(ii+1))) then
				info=1
				j=ii
				jp1=ii+1
				exit
			end if
		end do
		
		if (info==0) then
			write (stdout,'("x-value is out of range. Process halts")')
		end if
		
		AA=(xx(jp1)-x)/(xx(jp1)-xx(j))
		BB=(x-xx(j))/(xx(jp1)-xx(j))
		CC=(1.d0/6.d0)*((AA**3.d0)-AA)*((xx(jp1)-xx(j))**2.d0)
		DD=(1.d0/6.d0)*((BB**3.d0)-BB)*((xx(jp1)-xx(j))**2.d0)
		
		y=AA*yy(j)+BB*yy(jp1)+CC*dy2j(j)+DD*dy2j(jp1)
	
	
	end subroutine spline_interp


end module spline
