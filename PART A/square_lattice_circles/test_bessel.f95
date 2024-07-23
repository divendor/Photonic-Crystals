program test_bessel

! gfortran -o test_bessel.o ../modulos/tipos.f90 test_bessel.f95


use mcf_tipos

integer(kind = long)				:: i
real(kind = sp)						:: dx


dx = 0.1

open(unit = 1, file = "bessel1.dat", action = "write", status = "replace")

do i = 1, 1510	
	write(unit = 1, fmt = *) dx*(real(i)-1.0), bessel_j1(dx*(i-1.0))
	!print *, bessel_j1(dx*(i-1.0))
end do

close(unit = 1)


end program test_bessel

