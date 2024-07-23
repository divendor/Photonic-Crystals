module module_utilities
! F ../modulos/tipos.f90 module_utilities.f95

use mcf_tipos

! Physical constants:
real(kind = sp), parameter, public :: q_e  = 1.60217733*(10**(-19.0)), &
									  hbar = 1.05457266*(10**(-34.0)), &
								      c0   = 2.99792458*(10**(8.0)),   &
								      PI   = 4.0*atan(1.0) 
			   
public		:: make_ordered_Gpoints, fourier_coef, 			           &
			   kvector, vectormodule, kunitaryvector, linspace,        &
			   get_kvalues, file_to_integer_matrix,                     & 
			   print_integer_matrix, print_real_matrix
			 
private		:: count_valid_Gpoints, make_Gpoints, quick_sort

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

subroutine file_to_integer_matrix(file_name, Nlines, matrix)

	
	character(len = *), intent(in)			:: file_name
	integer, intent(in)						:: Nlines
	integer, dimension(:,:), intent(out)	:: matrix
	integer									:: i
	
	open(unit = 1, file = file_name, status = "old", action = "read")
	do i = 1, Nlines
		read(unit = 1, fmt = *) matrix(i, 1), matrix(i, 2)
	end do
	close(unit = 1)

end subroutine file_to_integer_matrix

subroutine print_integer_matrix(matrix)
	
	integer, dimension(:,:), intent(in)	:: matrix
	integer								:: i, Npoints
	integer, dimension(2)				:: shape_matrix
	
	shape_matrix = shape(matrix)
	Npoints = shape_matrix(1)
	do i = 1, Npoints
		print *, "[", matrix(i, :), "]"
	end do 

end subroutine print_integer_matrix

subroutine print_real_matrix(matrix)
	
	real(kind = sp), dimension(:,:), intent(in)	:: matrix
	integer										:: i, Npoints
	integer, dimension(2)						:: shape_matrix
	
	shape_matrix = shape(matrix)
	Npoints = shape_matrix(1)
	do i = 1, Npoints
		print *, "[", matrix(i, :), "]"
	end do 

end subroutine print_real_matrix
		
subroutine make_ordered_Gpoints(m_max, file_name, valid_points)
! Outputs a file with the indexes (m,n) and the number of valid_points
! that are within a radius given by m_max

	integer, intent(in)						:: m_max
	character(len = *), intent(in)			:: file_name
	integer, intent(out)					:: valid_points
	integer									:: m, n
	real									:: radius
	
	open(unit = 1, file = "Gpoints.dat", status = "replace", &
		 action = "write")
	valid_points = 0
	do m = -m_max, m_max
		do n = -m_max, m_max
			radius = sqrt(m**2.0 + n**2.0)
			if (radius <= m_max) then
				write(unit = 1, fmt= *) m, n 
				valid_points = valid_points + 1
			end if	
		end do
	end do
	close(unit = 1)

end subroutine make_ordered_Gpoints
	
subroutine get_kvalues(G_start, G_end, Npoints, step, k_values)
! returns a matrix containing the equally spaced points that want to be 
! used along the
! desired path.
! G_path is the vector representing the desired path
! G_start is a vector that points the start of the path: 
!					G_path = G_end - G_start

	real(kind = sp), dimension(:), intent(in)		:: G_start, G_end
	integer, intent(in)								:: Npoints
	real(kind = sp), dimension(:,:), intent(out)	:: k_values
	real(kind = sp)	,intent(out)					:: step
	real(kind = sp), dimension(2)					:: unitary, G_path
	real(kind = sp), dimension(Npoints) 			:: space, kx, ky
	real(kind = sp)									:: G_path_module
	integer											:: i
			
	G_path = G_end - G_start											
	call kunitaryvector(G_path, unitary)
	call vectormodule(G_path, G_path_module)
	call linspace(0.0, G_path_module, Npoints, space, step)
	
	kx = G_start(1) + unitary(1) * space
	ky = G_start(2) + unitary(2) * space
	
	!print *, unitary
	!print *, G_path_module
	!print *, space
	
	!print *, "kx", kx
	!print *, "ky", ky
	
	do i = 1, Npoints
		k_values(i, 1) = kx(i)
		k_values(i, 2) = ky(i)
	end do

end subroutine get_kvalues
	
!subroutine get_kvalues(G_start, G_path, Npoints, step, k_values)
!! returns a matrix containing the equally spaced points that want to be 
!! used along the
!! desired path.
!! G_path is the vector representing the desired path
!! G_start is a vector that points the start of the path: 
!!					G_path = G_end - G_start

!	real(kind = sp), dimension(:), intent(in)		:: G_start, G_path
!	integer, intent(in)								:: Npoints
!	real(kind = sp), dimension(:,:), intent(out)	:: k_values
!	real(kind = sp)	,intent(out)					:: step
!	real(kind = sp), dimension(2)					:: unitary
!	real(kind = sp), dimension(Npoints) 			:: space, kx, ky
!	real(kind = sp)									:: G_path_module
!	integer											:: i
														
!	call kunitaryvector(G_path, unitary)
!	call vectormodule(G_path, G_path_module)
!	call linspace(0.0, G_path_module, Npoints, space, step)
	
!	kx = G_start(1) + unitary(1) * space
!	ky = G_start(2) + unitary(2) * space
	
!	print *, unitary
!	print *, G_path_module
!	print *, space
	
!	print *, "kx", kx
!	print *, "ky", ky
	
!	do i = 1, Npoints
!		k_values(i, 1) = kx(i)
!		k_values(i, 2) = ky(i)
!	end do

!end subroutine get_kvalues

subroutine kvector(coord_x, coord_y, vector)
!retunrs a 2D k-space vector, given the x,y coordinates on each axis

	real(kind = sp), intent(in)					:: coord_x, coord_y
	real(kind = sp), dimension(:), intent(out)	:: vector

	vector(1) = coord_x
	vector(2) = coord_y

end subroutine kvector
	
subroutine vectormodule(vector, vector_module)
!returns the module of a 2D vector

	real(kind = sp), dimension(:), intent(in)	:: vector
	real(kind = sp), intent(out)				:: vector_module
	
	vector_module = sqrt(vector(1)**2 + vector(2)**2)
	
end subroutine vectormodule

subroutine kunitaryvector(vector, unitary)
! returns the unitary vector, given the vector and its module

	real(kind = sp), dimension(:), intent(in)	:: vector
	real(kind = sp), dimension(:), intent(out)	:: unitary
	real(kind = sp)								:: vector_module
	
	call vectormodule(vector, vector_module)
	unitary = vector / vector_module

end subroutine kunitaryvector

subroutine linspace(first_point, final_point, Npoints, space, step)
! divides the space between first_point and final_point
! in Npoints points and returns a vector containing each segment

	real(kind = sp), intent(in)						:: first_point, &
													   final_point
	integer, intent(in)								:: Npoints
	real(kind = sp), dimension(:), intent(out)		:: space
	real(kind = sp)	,intent(out)					:: step
	integer											:: i
	
	step = (final_point - first_point)/real((Npoints - 1))
	
	space(1) = first_point
	space(Npoints) = final_point

	do i = 2, Npoints -1
		space(i) = (i-1) * step
	end do

end subroutine linspace

subroutine fourier_coef(m, n, gam_1, gam_2, a, b, coef)
! Coefficient (m,n) of 2d-Fourier transform of gamma = 1/eps
! a is the square lattice's constant
! b is the thickness of the vein
	
	integer, intent(in)				:: m, n
	real(kind = sp), intent(in)		:: gam_1, gam_2, a, b
	real(kind = sp), intent(out)	:: coef
	
	if (m == 0) then
		if (n == 0) then
			coef = gam_2 + (gam_1 - gam_2)*(b/a)**2
		else
			coef = ((gam_1-gam_2)*b*((-1)**abs(n))*sin(n*PI*b/a)) & 
			       /(a*n*PI)
		end if
	else 
		if (n == 0) then
			coef = ((gam_1-gam_2)*b*((-1)**abs(m))*sin(m*PI*b/a)) & 
			       /(a*m*PI)
		else
			coef = ((gam_1-gam_2) & 
			       *((-1)**abs(m+n))*sin(m*PI*b/a)*sin(n*PI*b/a)) & 
			       /(m*n*PI**2.0)
		end if
	end if

end subroutine fourier_coef

end module module_utilities
