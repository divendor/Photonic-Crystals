module module_solver

! F ../modulos/tipos.f90 ../modulos/mcf_matrices.f95 module_utilities.f95 module_solver.f95

use mcf_tipos
use module_utilities
use mcf_diagonalizacion


public		:: make_matrixes, get_eigenvalues


contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine make_matrixes(G_points_matrix, k, TM, TE, b1, b2, &
						 gam_1, gam_2, a, r)
						 
	integer, dimension(:,:), intent(in)			:: G_points_matrix
	real(kind = sp), dimension(:), intent(in)	:: k, b1, b2
	real(kind = sp), intent(in)					:: gam_1, gam_2, a, r
	real(kind = sp), dimension(:,:), intent(out):: TE, TM
	integer										:: NGpoints, i, j, &
												   o, p, m, n
	real(kind = sp), dimension(2)				:: k_G_op, k_G_mn, &
												   shape_G_points_matrix
	real										:: module_k_G_op, &
												   module_k_G_mn, &
												   coef_fourier
	
	shape_G_points_matrix = shape(G_points_matrix)
	NGpoints = shape_G_points_matrix(1)
	
	do i = 1, NGpoints
		o = G_points_matrix(i,1)
		p = G_points_matrix(i,2)
		do j = i, NGpoints
			m = G_points_matrix(j,1)
			n = G_points_matrix(j,2)	
			
			k_G_op = k + o*b1 + p*b2
			k_G_mn = k + m*b1 + n*b2
			
			!print *, "k_G_op: ", "[", k_G_op,"]", "[",k_G_mn,"]"

			
			call vectormodule(k_G_op, module_k_G_op)
			call vectormodule(k_G_mn, module_k_G_mn)
			
			call fourier_coef(o - m, p - n, & 
							  gam_1, gam_2, a, r, coef_fourier)
			!print *, coef_fourier
			
			! TM
			TM(i,j) = coef_fourier * module_k_G_mn * module_k_G_op
			TM(j,i) = TM(i,j) 
			
			! TM
			TE(i,j) = coef_fourier * dot_product(k_G_mn, k_G_op)
			TE(j,i) = TE(i,j)
		end do
	end do 

end subroutine make_matrixes

subroutine get_eigenvalues(matrix,real_eigenvalues,complex_eigenvalues)
	
	real(kind = sp), dimension(:,:), intent(in):: matrix
	real(kind = sp), dimension(:), intent(out)	  :: real_eigenvalues, &
													 complex_eigenvalues
	real(kind = sp), dimension(:,:), allocatable  :: matrix_copy
	integer, dimension(2)						  :: matrix_shape
	integer									      :: matrix_size
	
	matrix_shape = shape(matrix)
	matrix_size = matrix_shape(1)
	allocate(matrix_copy(matrix_size, matrix_size))
	matrix_copy = matrix
	
	call eival(matrix_copy, real_eigenvalues, complex_eigenvalues)
	
	real_eigenvalues = c0*sqrt(real_eigenvalues)

end subroutine get_eigenvalues


end module module_solver
