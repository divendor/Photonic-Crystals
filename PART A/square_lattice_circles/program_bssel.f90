program test_my_bessel1



contains

subroutine my_bessel1(n, x, output)

	real(kind = sp), intent(in)		:: x
	integer							:: n
	real(kind = sp), intent(out)	:: output
	
	real(kind = sp)					:: aux
	integer							:: i, k

	
	output = 0.0
	aux = 1
	do i = 1, n
		
		k = 2*(i-1) + 1
		output = output + (x**k)*aux
	
		aux = 1/(real(aux**2.0)*real(aux + 1))
		
		aux
		
	end do



end subroutine my_bessel1




end program my_bessel1
