module gauss_legendre_sp
!
!Cuadratura de Gauss-Legendre unidimensional. Trasladada
!de Numerical Recipes Fortran77 a F. Simple precision.
!
use mcf_tipos

public :: qgauss
interface qgauss
	module procedure qgauss_sp
end interface
private :: qgauss_sp, gauleg_sp

contains

subroutine gauleg_sp(x1,x2,x,w,n)
integer, intent(in)                      :: n
real(kind=sp), intent(in)                :: x1,x2
real(kind=sp), intent(out), dimension(:) :: x,w
integer                                  :: i,j,m
!precision relativa
real(kind=dp), parameter                 :: eps=3.0e-14   
real(kind=dp)                            :: p1,p2,p3,pp,xl,xm,z,z1
!
!Dados los limites de integracion x1 (inferior) y x2 (superior), y dado n,
!esta subrutina devuelve las matrices x(1:n) y w(1:n) que contienen las
!abscisas y pesos necesarios para la cuadratura de Gauss-Legendre con n puntos.
!
!Los calculos son realizados en doble precision aunque el resultado se devuelva
!en simple precision.
!

!Las raices son simétricas en el intervalo, por lo que solo es preciso calcular 
!la mitad
m=(n+1)/2         

xm=0.5_dp*(x2+x1)
xl=0.5_dp*(x2-x1)

do i=1,m      !bucle sobre todas las raices
!
!Primera aproximacion de la raiz i-esima.
!
        z=cos(3.141592654_dp*(i-0.25_dp)/(n+0.5_dp))

!
!Refinamiento de la raiz (Metodo de Newton)
!
newton: do
            p1=1.0_dp
            p2=0.0_dp
!
!Recurrencia de los polinomios de Legendre evaluados en z.
!
            do j=1,n
                 p3=p2
                 p2=p1
                 p1=((2.0_dp*j-1.0_dp)*z*p2-(j-1.0_dp)*p3)/j
            end do
!
!p1 es el polinomio de Legendre. A continuacion calculamos su derivada, pp,
!mediante una relacion que require p2, el polinomio un orden mas bajo.
!
            pp=n*(z*p1-p2)/(z*z-1.0_dp)
            z1=z
!
!Metodo de Newton        
!
            z=z1-p1/pp
            
            if (abs(z-z1) <= eps) then
                exit newton
            end if
            
        end do newton
!
!Escalado de las raices al intervalo deseado
!
             x(i)=real(xm-xl*z,sp)
!
!Raiz simetrica
!             
             x(n+1-i)=real(xm+xl*z,sp)
!
!Calculo de los pesos
!
             w(i)=real(2.0_dp*xl/((1.0_dp-z*z)*pp*pp),sp)
!
!y su contrapartida simetrica
!
             w(n+1-i)=w(i)
end do
    
end subroutine gauleg_sp

subroutine qgauss_sp (func,a,b,ss,n_points)
integer, intent(in), optional            :: n_points
real(kind=sp), intent(in)                :: a,b
real(kind=sp), intent(out)               :: ss
integer, parameter                       :: n_default=5
integer                                  :: j,n
real(kind=sp), dimension(:), allocatable :: w,x
!
!Realiza la integral usando las coordenadas y pesos calculados 
!por gauleg_sp. Por defecto utiliza 5 puntos de integracion.
!
interface
function func(x) result (func_result)
use mcf_tipos
real(kind=sp), intent(in) :: x
real(kind=sp)             :: func_result
end function func
end interface

if (present(n_points)) then
     n=n_points
else
     n=n_default
end if

allocate (x(n))
allocate (w(n))

call gauleg_sp(a,b,x,w,n)

ss=0.0

do j=1,n
     ss=ss+w(j)*func(x(j))
end do

end subroutine qgauss_sp

end module gauss_legendre_sp
!
!---------------------------------------------------------------
!
module gauss_legendre_dp
!
!Cuadratura de Gauss-Legendre unidimensional. Trasladada
!de Numerical Recipes Fortran77 a F. Doble precision.
!
use mcf_tipos

public :: qgauss
interface qgauss
	module procedure qgauss_dp
end interface
private :: qgauss_dp, gauleg_dp

contains

subroutine gauleg_dp(x1,x2,x,w,n)
integer, intent(in)                      :: n
real(kind=dp), intent(in)                :: x1,x2
real(kind=dp), intent(out), dimension(:) :: x,w
integer                                  :: i,j,m
!precision relativa
real(kind=dp), parameter                 :: eps=3.0e-14   
real(kind=dp)                            :: p1,p2,p3,pp,xl,xm,z,z1
!
!Dados los limites de integracion x1 (inferior) y x2 (superior), y dado n,
!esta subrutina devuelve las matrices x(1:n) y w(1:n) que contienen las
!abscisas y pesos necesarios para la cuadratura de Gauss-Legendre con n puntos.
!
!Los calculos son realizados en doble precision y el resultado es devuelto
!en doble precision.
!

!Las raices son simétricas en el intervalo, por lo que solo es preciso calcular 
!la mitad
m=(n+1)/2         

xm=0.5_dp*(x2+x1)
xl=0.5_dp*(x2-x1)

do i=1,m      !bucle sobre todas las raices
!
!Primera aproximacion de la raiz i-esima
!
        z=cos(3.141592654_dp*(i-0.25_dp)/(n+0.5_dp))

!
!Refinamiento de la raiz (Metodo de Newton)
!
newton: do
            p1=1.0_dp
            p2=0.0_dp
!
!Recurrencia de os polinomios de Legendre evaluados en z.
!
            do j=1,n
                 p3=p2
                 p2=p1
                 p1=((2.0_dp*j-1.0_dp)*z*p2-(j-1.0_dp)*p3)/j
            end do
!
!p1 es el polinomio de Legendre. A continuacion calculamos su derivada, pp,
!mediante una relacion que require p2, el polinomio un orden mas bajo.
!
            pp=n*(z*p1-p2)/(z*z-1.0_dp)
            z1=z
!
!Metodo de Newton        
!
            z=z1-p1/pp
            
            if (abs(z-z1) <= eps) then
                exit newton
            end if
            
        end do newton
!
!Escalado de las raices al intervalo deseado
!
             x(i)=xm-xl*z
!
!Raiz simetrica
!             
             x(n+1-i)=xm+xl*z
!
!Calculo de los pesos
!
             w(i)=2.0_dp*xl/((1.0_dp-z*z)*pp*pp)
!
!y su contrapartida simetrica
!
             w(n+1-i)=w(i)
end do
    
end subroutine gauleg_dp

subroutine qgauss_dp (func,a,b,ss,n_points)
integer, intent(in), optional            :: n_points
real(kind=dp), intent(in)                :: a,b
real(kind=dp), intent(out)               :: ss
integer, parameter                       :: n_default=5
integer                                  :: j,n
real(kind=dp), dimension(:), allocatable :: w,x
!
!Realiza la integral usando las coordenadas y pesos calculados 
!por gauleg_dp. Por defecto utiliza 5 puntos de integracion.
!
interface
function func(x) result (func_result)
use mcf_tipos
real(kind=dp), intent(in) :: x
real(kind=dp)             :: func_result
end function func
end interface

if (present(n_points)) then
     n=n_points
else
	 n=n_default
end if

allocate (x(n))
allocate (w(n))

call gauleg_dp(a,b,x,w,n)

ss=0.0

do j=1,n
print "(i3,2f15.10)",j,x(j),w(j)
     ss=ss+w(j)*func(x(j))
end do

end subroutine qgauss_dp

end module gauss_legendre_dp
!
!---------------------------------------------------------------
!
module gauss_legendre
!
!Cuadratura de Gauss-Legendre
!
use gauss_legendre_sp
use gauss_legendre_dp

public

end module gauss_legendre
!
!---------------------------------------------------------------
!
