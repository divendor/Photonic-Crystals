module bessel

use mcf_tipos



public :: besselj0,besselj1,besseljn,sinc

private

real(kind=sp), parameter  :: eps = epsilon(1.0_sp), pi = acos(-1.0_sp), &
                             pi2 = pi/2.0_sp, pi4 = pi/4.0_sp 

real(kind=sp), parameter  :: asymptotic_j0 = 40.0_sp, &
                             asymptotic_j1 = 30.0_sp

integer, public           :: n

contains

function besselj0 (x) result (j0_result)
!***************************************************************
! Esta funcion  calcula la funcion de Bessel J de orden 0.     *
! Empiricamente esta dividido el intervalo entre [0,40]        *
! donde se utiliza una relacion de recurrencia, mientras       *
! que para valores de x superiores a 40 se utiliza una         *
! ecuacion asintotica. No es el mejor algoritmo dado que,      *
! en general se utilizan aproximaciones polinomicas en         *
! diferente intervalos. Jn puede calcularse iterativamente a   *
! partir de J0 y J1.                                           *
!***************************************************************
real(kind=sp), intent(in) :: x
real(kind=sp)             :: j0_result

real(kind=sp)             :: contador, correccion, factor, xlocal

j0_result = 1.0_sp

if (x == 0.0_sp) then

    return

else

    xlocal = abs(x)

    if (xlocal <= asymptotic_j0) then
!
! Procedimiento iterativo basado en la serie de J0.
!
        xlocal     = xlocal/2.0_sp
        contador   = 0.0_sp
        correccion = 1.0_sp

        do

            contador   = contador + 1.0_sp
            factor     = xlocal/contador
            correccion = -correccion*factor*factor

            if (abs(correccion) < eps) then

                exit

            else

                j0_result = j0_result + correccion

            end if

        end do

    else

        j0_result = sqrt(1.0_sp/(pi2*xlocal))*cos(xlocal-pi4)

    end if

end if

end function besselj0

function besselj1 (x) result (j1_result)
!***************************************************************
! Esta funcion  calcula la funcion de Bessel J de orden 1.     *
! Empiricamente esta dividido el intervalo entre [0,30]        *
! donde se utiliza una relacion de recurrencia, mientras       *
! que para valores de x superiores a 30 se utiliza una         *
! ecuacion asintotica. No es el mejor algoritmo dado que,      *
! en general se utilizan aproximaciones polinomicas en         *
! diferente intervalos. Jn puede calcularse iterativamente a   *
! partir de J0 y J1.                                           *
!***************************************************************
real(kind=sp), intent(in) :: x
real(kind=sp)             :: j1_result

real(kind=sp)             :: contador, correccion, factor, xlocal


if (x == 0.0_sp) then

        j1_result = 0.0_sp

    return

else

    xlocal = abs(x)

    if (xlocal <= asymptotic_j1) then
!
! Procedimiento iterativo basado en la serie de J1.
!
        xlocal     = xlocal/2.0_sp

        j1_result  = xlocal

        correccion = xlocal
        factor     = xlocal*xlocal
        contador   = 0.0_sp

        do

            contador   = contador + 1.0_sp
            correccion = -correccion*factor/contador/(contador+1.0_sp)

            if (abs(correccion) < eps) then

                exit

            else

                j1_result = j1_result + correccion

            end if

        end do

    else

        j1_result = sqrt(1.0_sp/(pi2*xlocal))*cos(xlocal-pi2-pi4)

    end if

end if

j1_result = sign(1.0_sp,x)*j1_result

end function besselj1

function besseljn(x) result (jn_result)
!***************************************************************
! Esta funcion  calcula la funcion de Bessel J de orden n.     *
! El calculo se realiza iterativamente a partir de J0 y J1.    *
!                                                              *
!            J(n,x) = [2(n-1)/x] J(n-1,x) - J(n-2,x)           *
!                                                              *
! Aunque simple este algoritmo es inestable porque la propia   *
! recurrencia lo es.                                           *
!                                                              *
! Hasta n = 24 funciona razonablemente bien                    *
!***************************************************************
real(kind=sp), intent(in) :: x
real(kind=sp)             :: jn_result

real(kind=sp)             :: jn_1, jn_2, xlocal, xlimit
integer                   :: i, nlocal
integer, parameter        :: nlimit = 24


xlocal = x/2.0_sp
nlocal = abs(n)

!Si |n| > nlimit el calculo no es bueno. No se debe seguir

if (nlocal > nlimit) then

    stop "besseljn: Si |n|>24 el algoritmo es divergente en zonas de valor nulo. " // &
         "Es inadmisible."

end if

!Limite para el rango asintotico para x pequenas. De manera empirica este
!control es necesario para |n| > 8.

xlimit = 1.0_sp*sqrt(real(nlocal,kind=sp)+1.0_sp)

!Punto de partida

jn_2 = besselj0(x)
jn_1 = besselj1(x)

if (n == 0) then

    jn_result = jn_2

else if (n == 1) then

    jn_result = jn_1

else if (abs(x) <= xlimit .and. nlocal > 8) then

    jn_result = 1.0_sp

    do i=1,nlocal

        jn_result = jn_result*(xlocal/real(i,kind=sp))

    end do

else

    do i=2,nlocal-1

        jn_result = (real(i,kind=sp)-1.0_sp)*jn_1/xlocal - jn_2
        jn_2      = jn_1
        jn_1      = jn_result

    end do

    jn_result = (nlocal-1)*jn_1/xlocal - jn_2

end if

jn_result = jn_result*sign(1,n)**nlocal

end function besseljn

function sinc (x) result (sinc_result)

real(kind=sp), intent(in) :: x
real(kind=sp)             :: sinc_result

if (x == 0.0_sp) then

    sinc_result = 1.0_sp

else

    sinc_result = sin(x)/x

end if

end function sinc

end module bessel
