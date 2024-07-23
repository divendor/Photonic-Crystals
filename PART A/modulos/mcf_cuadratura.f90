
module m_cuadratura_sp

use mcf_tipos
private
  !
  ! Routines to integrate functions
  !
public  :: trapecio, romberg
interface trapecio
   module procedure trapecio_sp
end interface
interface romberg
   module procedure romberg_sp
end interface
private :: trapecio_sp
private :: romberg_sp
!
CONTAINS !=====================================

SUBROUTINE TRAPECIO_sp(FUNC,A,B,S,N)
!
! Calcula una aproximacion S a la integral de func(x) entre a y b. 
! Usa el metodo del trapecio con 2**N puntos
!
! Ha de ser llamada sucesivas veces, sin modificar S,
! con N empezando en 1 y 
! llegando hasta el valor necesario, segun la precision deseada.
!
real(kind=sp), intent(in)    :: A,B
real(kind=sp), intent(inout) :: S
integer, intent(in)              :: N 

interface
  function func(x) result(res_func)
  use mcf_tipos
  real(kind=sp), intent(in)   :: x
  real(kind=sp)               :: res_func
  end function func
end interface

real(kind=sp)         ::   DEL,SUM,X
INTEGER                   ::   J, TNM
integer, save             ::  IT

IF (N == 1) THEN  
    ! Primera entrada. h=(b-a)
    !
    S = 0.5_sp* (B-A)* (FUNC(A)+FUNC(B))
    IT = 1
ELSE
    !
    ! Tomamos mas puntos intermedios
    !
    TNM = IT
    DEL = (B-A)/TNM
    X = A + 0.5_sp*DEL
    SUM = 0.0_sp
    DO J = 1,IT
        SUM = SUM + FUNC(X)
        X = X + DEL
    ENDDO
    !
    ! Sumamos la contribucion de estos puntos 
    !
    S = 0.50_sp* (S+ (B-A)*SUM/TNM)
    !
    ! Preparamos IT para la siguiente iteracion
    !
    IT = 2*IT

END IF

END SUBROUTINE TRAPECIO_sp

SUBROUTINE romberg_sp(FUNC,A,B,SS,eps,order,debug,converged)
use mcf_interpoli

! Calcula una aproximacion SS a la integral
! de func(x) entre a y b. Usa el metodo
! del trapecio con extrapolacion Romberg en h**2 para alcanzar
! la precision eps.
!
! El argumento opcional order determina el numero de puntos en
! la extrapolacion (por defecto es cinco)
! Si el argumento opcional debug es .true., la rutina da informacion
! detallada acerca de su funcionamiento.
! El argumento opcional converged sera .true. a la salida 
! si la rutina tiene exito y .false. en caso contrario.

real(kind=sp), intent(in)    :: A,B
real(kind=sp), intent(out)   :: SS
real(kind=sp), intent(in)    :: eps
integer, intent(in), optional    :: order
logical, intent(in), optional    :: debug
logical, intent(out), optional   :: converged

interface
  function func(x) result(res_func)
  use mcf_tipos
  real(kind=sp), intent(in)   :: x
  real(kind=sp)               :: res_func
  end function func
end interface

INTEGER, parameter   ::  JMAX = 20 ,JMAXP = jmax+1
integer, parameter   ::  K_default=5  ! Numero de puntos en la extrapolacion

real(kind=sp)         ::  DSS
INTEGER                   ::  K, J, jmin
logical                   ::  verbose
logical                   ::  converged_flag

real(kind=sp), dimension(jmaxp)     ::    H1, H, S

verbose = (present(debug) .and. debug)

if (present(order)) then
   k = order
else
   k = k_default
endif

H1(1) = 1.0_sp
H(1) = 1.0_sp
DO J = 1,JMAX
    call trapecio_sp(FUNC,A,B,S(J),J)
    if (verbose) then
        print *, "j, S: ", j, S(j)
    end if
    IF (J >= K) THEN
        ! Usamos una extrapolacion a h=0 con k puntos
        jmin = j - k + 1
        CALL POLINT(H(jmin:j),S(jmin:j),K,0.0_sp,SS,DSS)
        if (verbose) then
           print "(a14,5f11.7)", "h usados:",H(jmin:j)
           print "(a14,5f11.7)", "S usados:",S(jmin:j)
           print *, "Extrapolacion hacia h=0: ", SS
        endif
        IF (ABS(DSS) < eps*ABS(SS)) then
           exit
        end if
    END IF

    S(J+1) = S(J)
    H(J+1) = 0.25_sp*H(J)    ! La extrapolacion es en h**2...
ENDDO

if (j>jmax) then
   if (verbose) then
     print *,"ROMBERG did not converge in ", Jmax, " iterations."
   end if
   converged_flag = .false.
else
   converged_flag = .true.
endif
if (present(converged)) then
        converged = converged_flag
end if
END subroutine romberg_sp


end module m_cuadratura_sp




module m_cuadratura_dp

use mcf_tipos
private
  !
  ! Routines to integrate functions
  !
public  :: trapecio, romberg
interface trapecio
   module procedure trapecio_dp
end interface
interface romberg
   module procedure romberg_dp
end interface
private :: trapecio_dp
private :: romberg_dp
!
CONTAINS !=====================================

SUBROUTINE TRAPECIO_dp(FUNC,A,B,S,N)
!
! Calcula una aproximacion S a la integral de func(x) entre a y b. 
! Usa el metodo del trapecio con 2**N puntos
!
! Ha de ser llamada sucesivas veces, sin modificar S,
! con N empezando en 1 y 
! llegando hasta el valor necesario, segun la precision deseada.
!
real(kind=dp), intent(in)    :: A,B
real(kind=dp), intent(inout) :: S
integer, intent(in)              :: N 

interface
  function func(x) result(res_func)
  use mcf_tipos
  real(kind=dp), intent(in)   :: x
  real(kind=dp)               :: res_func
  end function func
end interface

real(kind=dp)         ::   DEL,SUM,X
INTEGER                   ::   J, TNM
integer, save             ::  IT

IF (N == 1) THEN  
    ! Primera entrada. h=(b-a)
    !
    S = 0.5_dp* (B-A)* (FUNC(A)+FUNC(B))
    IT = 1
ELSE
    !
    ! Tomamos mas puntos intermedios
    !
    TNM = IT
    DEL = (B-A)/TNM
    X = A + 0.5_dp*DEL
    SUM = 0.0_dp
    DO J = 1,IT
        SUM = SUM + FUNC(X)
        X = X + DEL
    ENDDO
    !
    ! Sumamos la contribucion de estos puntos 
    !
    S = 0.50_dp* (S+ (B-A)*SUM/TNM)
    !
    ! Preparamos IT para la siguiente iteracion
    !
    IT = 2*IT

END IF

END SUBROUTINE TRAPECIO_dp

SUBROUTINE romberg_dp(FUNC,A,B,SS,eps,order,debug,converged)
use mcf_interpoli

! Calcula una aproximacion SS a la integral
! de func(x) entre a y b. Usa el metodo
! del trapecio con extrapolacion Romberg en h**2 para alcanzar
! la precision eps.
!
! El argumento opcional order determina el numero de puntos en
! la extrapolacion (por defecto es cinco)
! Si el argumento opcional debug es .true., la rutina da informacion
! detallada acerca de su funcionamiento.
! El argumento opcional converged sera .true. a la salida 
! si la rutina tiene exito y .false. en caso contrario.

real(kind=dp), intent(in)    :: A,B
real(kind=dp), intent(out)   :: SS
real(kind=dp), intent(in)    :: eps
integer, intent(in), optional    :: order
logical, intent(in), optional    :: debug
logical, intent(out), optional   :: converged

interface
  function func(x) result(res_func)
  use mcf_tipos
  real(kind=dp), intent(in)   :: x
  real(kind=dp)               :: res_func
  end function func
end interface

INTEGER, parameter   ::  JMAX = 20 ,JMAXP = jmax+1
integer, parameter   ::  K_default=5  ! Numero de puntos en la extrapolacion

real(kind=dp)         ::  DSS
INTEGER                   ::  K, J, jmin
logical                   ::  verbose
logical                   ::  converged_flag

real(kind=dp), dimension(jmaxp)     ::    H1, H, S

verbose = (present(debug) .and. debug)

if (present(order)) then
   k = order
else
   k = k_default
endif

H1(1) = 1.0_dp
H(1) = 1.0_dp
DO J = 1,JMAX
    call trapecio_dp(FUNC,A,B,S(J),J)
    if (verbose) then
        print *, "j, S: ", j, S(j)
    end if
    IF (J >= K) THEN
        ! Usamos una extrapolacion a h=0 con k puntos
        jmin = j - k + 1
        CALL POLINT(H(jmin:j),S(jmin:j),K,0.0_dp,SS,DSS)
        if (verbose) then
           print "(a14,5f11.7)", "h usados:",H(jmin:j)
           print "(a14,5f11.7)", "S usados:",S(jmin:j)
           print *, "Extrapolacion hacia h=0: ", SS
        endif
        IF (ABS(DSS) < eps*ABS(SS)) then
           exit
        end if
    END IF

    S(J+1) = S(J)
    H(J+1) = 0.25_dp*H(J)    ! La extrapolacion es en h**2...
ENDDO

if (j>jmax) then
   if (verbose) then
        print *,"ROMBERG did not converge in ", Jmax, " iterations."
   end if
   converged_flag = .false.
else
   converged_flag = .true.
endif
if (present(converged)) then
        converged = converged_flag
end if
END subroutine romberg_dp


end module m_cuadratura_dp


module mcf_cuadratura
  !
  ! Integracion de funciones
  !
use m_cuadratura_sp
use m_cuadratura_dp

public

end module mcf_cuadratura





