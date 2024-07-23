program square_lattice_circles
! F -o out ../modulos/tipos.f90 ../modulos/mcf_matrices.f95 module_utilities.f95 module_solver.f95 square_lattice_circles.f95
! gfortran -o out ../modulos/tipos.f90 ../modulos/mcf_matrices.f95 module_utilities.f95 module_solver.f95 square_lattice_circles.f95



use mcf_tipos
use mcf_diagonalizacion
use module_utilities
use module_solver


real(kind = sp)									:: a, r_a, d, r,       &
												   epsr_1, gam_1,      & 
												   epsr_2, gam_2,      &
												   dk, Gmax_factor,    &  
												   Gmax, coef, 		   &
												   last_point, aux
integer											:: Nkpoints, Nbands,   & 
												   m_max, 				&		 	    
												   considered_points, i
real(kind = sp), dimension(2)					:: b1, b2, aux_vec
												   
integer, dimension(:,:), allocatable			:: G_points 
												 
integer, dimension(:), allocatable				:: G_x_array, 		   & 
												   G_y_array
												   
												   
real(kind = sp), dimension(:,:), allocatable	:: TM, TE,			   &
												   path_GAMMA_X,       &
												   path_X_M,		   &
												   path_M_GAMMA
real(kind = sp), dimension(:), allocatable		:: TM_eigenval_real,   & 
												   TM_eigenval_complex,&
												   TE_eigenval_real,   & 
											       TE_eigenval_complex
											   									       
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! MATERIAL: A Square Lattice of Dielectric Veins

! Grid settings
a	= 10 				! lattice constant 			[nm]
r_a	= 0.2				! d/a 						[dimless]
r	= a * r_a 			! thickness of grid veins	[nm]

! Material 1 == "Hole" material:

! mu_r1 = 1.0  			  relative magnetic   permeability
epsr_1 = 8.9			! relative dielectric permittivity
gam_1   = 1.0/(epsr_1)

! Material 2 == "Grid" material:

! mu_r2 = 1.0             relative magnetic   permeability
epsr_2 = 1.0 			! relative dielectric permittivity
gam_2 = 1.0/(epsr_2)  

! Vectores primitivos

b1 = (/2*pi/a, 0.0/)
b2 = (/0.0, 2*pi/a/)

!print *, "b1: ", b1
!print *, "b2: ", b2
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! Calculation settings

Nkpoints = 30		! Number of k-points on path:
Nbands = 10				! Unused yet
dk = PI/(a*(Nkpoints - 1))	! distance between k-states.  

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! Plane wave settings

Gmax_factor = 15.0 ! it controls the number of considered plane waves
Gmax = Gmax_factor*2*PI/a
	m_max = floor(Gmax_factor) ! closest integer below
!m_max = 17 ! funciona
m_max = 10

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! 1. WE MAKE THE RECIPROCAL LATTICE WITH THE CONSIDERED
!    RECIPROCAL LATTICE VECTORS 

call make_ordered_Gpoints(m_max, "Gpoints.dat", considered_points)
allocate(G_points(considered_points, 2))
call file_to_integer_matrix("Gpoints.dat", considered_points,& 
							G_points)
							
! 2. WE MAKE THE DESIRED k-PATHS ALONG HIGH SYMMETRY DIRECTIONS and
!    solve the eigenvalue problem

allocate(TM(considered_points, considered_points))
allocate(TM_eigenval_real(considered_points))
allocate(TM_eigenval_complex(considered_points))

allocate(TE(considered_points, considered_points))
allocate(TE_eigenval_real(considered_points))
allocate(TE_eigenval_complex(considered_points))


last_point = 0.0

open(unit = 30, file = "TM.dat", status = "replace", &
		 action = "write")
open(unit = 40, file = "TE.dat", status = "replace", &
		 action = "write")

! Path GAMMA --> X......................................................

allocate(path_GAMMA_X(Nkpoints, 2))

call get_kvalues((/0.0, 0.0/), (/PI/a, 0.0/), Nkpoints, dk, &
                 path_GAMMA_X)

print *, last_point

open(unit = 10, file = "TM_GAMMA_X.dat", status = "replace", &
		 action = "write")
open(unit = 20, file = "TE_GAMMA_X.dat", status = "replace", &
		 action = "write")
print *, "Computing path: ", 1
do i = 1, Nkpoints
	call make_matrixes(G_points, path_GAMMA_X(i,:), TM, TE, b1, b2, &
					   gam_1, gam_2, a, r)
	
	call get_eigenvalues(TM, TM_eigenval_real, TM_eigenval_complex)
	call get_eigenvalues(TE, TE_eigenval_real, TE_eigenval_complex)
	
	write(unit = 10, fmt = *) (i-1)*dk, & 
	                          TM_eigenval_real(:)*a/(2.0*PI*c0)
	write(unit = 20, fmt = *) (i-1)*dk, & 
	                          TE_eigenval_real(:)*a/(2.0*PI*c0)
	
	write(unit = 30, fmt = *) (i-1)*dk, & 
	                          TM_eigenval_real(:)*a/(2.0*PI*c0)
	write(unit = 40, fmt = *) (i-1)*dk, & 
	                          TE_eigenval_real(:)*a/(2.0*PI*c0)
	
	print *, "points computed: ", i

end do 
close(unit = 10)
close(unit = 20)

! Path X --> M..........................................................

allocate(path_X_M(Nkpoints, 2))

call get_kvalues((/PI/a, 0.0/), (/PI/a, PI/a/), Nkpoints, dk, &
                 path_X_M)
                
 
call vectormodule(path_GAMMA_X(Nkpoints,:), aux)
last_point = last_point + aux  
print *, last_point      
        
print *, "k_values obtained"

open(unit = 11, file = "TM_X_M.dat", status = "replace", &
		 action = "write")
open(unit = 21, file = "TE_X_M.dat", status = "replace", &
	 action = "write")

print *, "Computing path: ", 2
do i = 1, Nkpoints
	call make_matrixes(G_points, path_X_M(i,:), TM, TE, b1, b2, &
					   gam_1, gam_2, a, r)
	
	call get_eigenvalues(TM, TM_eigenval_real, TM_eigenval_complex)
	call get_eigenvalues(TE, TE_eigenval_real, TE_eigenval_complex)
	
	write(unit = 11, fmt = *) (i-1)*dk, TM_eigenval_real(:)*a/(2*PI*c0)
	write(unit = 21, fmt = *) (i-1)*dk, & 
	                          TE_eigenval_real(:)*a/(2.0*PI*c0)	  
	                          
	write(unit = 30, fmt = *) (i-1)*dk + last_point, & 
	                          TM_eigenval_real(:)*a/(2.0*PI*c0)
	write(unit = 40, fmt = *) (i-1)*dk + last_point, & 
	                          TE_eigenval_real(:)*a/(2.0*PI*c0)                       
	                        
	print *, "points computed: ", i

end do 
close(unit = 11)
close(unit = 21)


! Path M --> GAMMA......................................................

allocate(path_M_GAMMA(Nkpoints, 2))

call get_kvalues((/PI/a, PI/a/), (/0.0, 0.0/), Nkpoints, dk, &
                 path_M_GAMMA)

aux_vec = path_X_M(Nkpoints,:) - path_X_M(1,:)

call vectormodule(aux_vec, aux)

print *, "ultimo punto: ", path_X_M(Nkpoints,:)
print *, aux
last_point = last_point + aux

print *, last_point

open(unit = 12, file = "TM_M_GAMMA.dat", status = "replace", &
		 action = "write")
open(unit = 22, file = "TE_M_GAMMA.dat", status = "replace", &
		 action = "write")
print *, "Computing path: ", 3
do i = 1, Nkpoints
	call make_matrixes(G_points, path_M_GAMMA(i,:), TM, TE, b1, b2, &
					   gam_1, gam_2, a, r)
	
	call get_eigenvalues(TM, TM_eigenval_real, TM_eigenval_complex)
	call get_eigenvalues(TE, TE_eigenval_real, TE_eigenval_complex)
	
	write(unit = 12, fmt = *) (i-1)*dk, & 
	                          TM_eigenval_real(:)*a/(2.0*PI*c0)
	write(unit = 22, fmt = *) (i-1)*dk, & 
	                          TE_eigenval_real(:)*a/(2.0*PI*c0)
	
	write(unit = 30, fmt = *) (i-1)*dk + last_point, & 
	                          TM_eigenval_real(:)*a/(2.0*PI*c0)
	write(unit = 40, fmt = *) (i-1)*dk + last_point, & 
	                          TE_eigenval_real(:)*a/(2.0*PI*c0)
	
	print *, "points computed: ", i

end do 
close(unit = 12)
close(unit = 22)

close(unit = 30)
close(unit = 40)


end program square_lattice_circles
