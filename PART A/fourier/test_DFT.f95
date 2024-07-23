program test_DFT

! F -o out ../modulos/tipos.f90 test_DFT.f95

use mcf_tipos

real(kind = sp), dimension(:), allocatable			:: y_original, & 
													   x_original
real(kind = sp)										:: dx													
integer												:: N, m
real(kind = sp), dimension(:,:), allocatable		:: DFT


dx = 0.1
N = 100
allocate(y_original(N))
allocate(x_original(N))
allocate(DFT(N,N))

open(unit = 1, file = "datos.dat", action = "write", status = "replace")
do m = 1, N/2
	
	x_original(m) = -(N/2-(m-1))*dx
	x_original(N-(m-1)) = (N/2-(m-1))*dx
	y_original(m) = -0.5*x_original(m) 
	y_original(N-(m-1)) = 0.5*x_original(N-(m-1))
	write(unit = 1, fmt = *) x_original(m), y_original(m)
	write(unit = 1, fmt = *) x_original(N-(m-1)), y_original(N-(m-1))
end do
close(unit = 1)

do m = 1, N
	print *, x_original(m)
end do






end program test_DFT










