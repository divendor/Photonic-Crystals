module m_slineales_sp

use mcf_tipos
private
!
! Routines to solve systems of linear equations
!
public :: tridag, ctridag, lu_descomposicion, lu_resolucion, gaussj
public :: ludcmp, lubksb     ! Aliases
public :: print_matrix

interface tridag
   module procedure tridag_sp
end interface
private :: tridag_sp

interface ctridag
   module procedure ctridag_sp
end interface
private :: ctridag_sp

interface gaussj
   module procedure gaussj_sp
end interface
private :: gaussj_sp

interface ludcmp
   module procedure lu_descomposicion_sp
end interface

interface lubksb
   module procedure lu_resolucion_sp
end interface
!
interface lu_descomposicion
   module procedure lu_descomposicion_sp
end interface
private :: lu_descomposicion_sp

interface lu_resolucion
   module procedure lu_resolucion_sp
end interface
private :: lu_resolucion_sp

interface print_matrix
   module procedure print_matrix_sp
end interface
private :: print_matrix_sp

CONTAINS  !==================================================

SUBROUTINE TRIDAG_sp(A,D,C,B,X,n)
!
! Resuelve un sistema tridiagonal de ecuaciones lineales.
! D(1:n) contiene los elementos diagonales
! a(1:n) contiene los elementos infradiagonales (a(1) no se usa)
! c(1:n) contiene los elementos supradiagonales (c(n) no se usa)
! B(1:n) contiene los terminos independientes del sistema
!
! A la salida, x(1:n) contiene la solucion del sistema
!
!
integer, intent(in)              :: n
real(kind=sp), intent(in), dimension(:)   :: a, d, c, b
real(kind=sp), intent(out), dimension(:)  :: x

real(kind=sp)    ::  BET
INTEGER :: J

real(kind=sp), allocatable, dimension(:)  ::  GAM     ! Array de trabajo

IF (D(1) == 0.0_sp) THEN
    STOP "D(1)=0 in tridag"
END IF

allocate(gam(n))

BET = D(1)
X(1) = B(1)/BET
DO J = 2,N
    GAM(J) = C(J-1)/BET
    BET = D(J) - A(J)*GAM(J)
    IF (BET == 0.0_sp) THEN
        deallocate(gam)
        STOP "BET=0.0, Tridag failed."
    END IF

    X(J) = (B(J)-A(J)*X(J-1))/BET
ENDDO
DO J = N - 1,1,-1
    X(J) = X(J) - GAM(J+1)*X(J+1)
ENDDO

deallocate(gam)

END subroutine tridag_sp

SUBROUTINE ctridag_sp(A,D,C,B,X,n)
!
! Resuelve un sistema tridiagonal de ecuaciones lineales.
! D(1:n) contiene los elementos diagonales
! a(1:n) contiene los elementos infradiagonales (a(1) no se usa)
! c(1:n) contiene los elementos supradiagonales (c(n) no se usa)
! B(1:n) contiene los terminos independientes del sistema
!
! A la salida, x(1:n) contiene la solucion del sistema
!
!
integer, intent(in)                         :: n
complex(kind=sp), intent(in), dimension(:)  :: a, d, c, b
complex(kind=sp), intent(out), dimension(:) :: x

complex(kind=sp)                            :: BET
INTEGER                                     :: J

complex(kind=sp), allocatable, dimension(:) :: GAM     ! Array de trabajo

IF (D(1) == 0.0_sp) THEN
    STOP "D(1)=0 in Ctridag"
END IF

allocate(gam(n))

BET = D(1)
X(1) = B(1)/BET
DO J = 2,N
    GAM(J) = C(J-1)/BET
    BET = D(J) - A(J)*GAM(J)
    IF (BET == 0.0_sp) THEN
        deallocate(gam)
        STOP "BET=0.0, Ctridag failed."
    END IF

    X(J) = (B(J)-A(J)*X(J-1))/BET
ENDDO
DO J = N - 1,1,-1
    X(J) = X(J) - GAM(J+1)*X(J+1)
ENDDO

deallocate(gam)

END subroutine ctridag_sp
!===================================================================

SUBROUTINE LU_DESCOMPOSICION_sp(A,INDX,D)
!
! Performs the LU decomposition of a square matrix A
! indx is an array holding permutation indexes. It must
! be kept unchanged before calling lu_resolucion.
! d is the signature of the permutation.
! On exit, A contains the L and U matrices.
!
real(kind=sp), dimension(:,:) , intent(inout)   :: A
INTEGER, dimension(:), intent(out)                  :: INDX
real(kind=sp), intent(out)                      :: D

real(kind=sp)                     ::  AAMAX,DUM,SUM
INTEGER                               ::  n,I,IMAX,J,K

real(kind=sp), PARAMETER      :: TINY=1.0e-20_sp
                                  ! To avoid division by zero

real(kind=sp), dimension(:), allocatable ::  VV

n = size(a,dim=1)
if (size(a,dim=2) /= n) then
        STOP "Matrix is not square in lu_descomposicion"
end if
if (size(indx) < n) then
        STOP "Size of Indx is too small in lu_descomposicion"
end if
allocate(vv(n))

D = 1.0_sp
DO I = 1,N
    !! aamax = maxval(abs(a(i,1:n)))
    AAMAX = 0.0_sp
    DO J = 1,N
        IF (ABS(A(I,J)) > AAMAX) then
                AAMAX = ABS(A(I,J))
        end if
    ENDDO
    IF (AAMAX == 0.0_sp) THEN
        deallocate(vv)
        STOP "Singular matrix in lu_descomposicion"
    END IF
    VV(I) = 1.0_sp/AAMAX
ENDDO

DO J = 1,N
    IF (J > 1) THEN
        DO I = 1,J - 1
            SUM = A(I,J)
            IF (I > 1) THEN
                DO K = 1,I - 1
                    SUM = SUM - A(I,K)*A(K,J)
                ENDDO
                A(I,J) = SUM
            END IF
        ENDDO
    END IF

    AAMAX = 0.0_sp
    DO I = J,N
        SUM = A(I,J)
        IF (J > 1) THEN
            DO K = 1,J - 1
                SUM = SUM - A(I,K)*A(K,J)
            ENDDO
            A(I,J) = SUM
        END IF

        DUM = VV(I)*ABS(SUM)
        IF (DUM >= AAMAX) THEN
            IMAX = I
            AAMAX = DUM
        END IF
    ENDDO
    IF (J /= IMAX) THEN
        DO K = 1,N
            DUM = A(IMAX,K)
            A(IMAX,K) = A(J,K)
            A(J,K) = DUM
        ENDDO
        D = -D
        VV(IMAX) = VV(J)
    END IF

    INDX(J) = IMAX
    IF (J /= N) THEN
        IF (A(J,J) == 0.0_sp) then
                A(J,J) = TINY
        end if
        DUM = 1.0_sp/A(J,J)
        DO I = J + 1,N
            A(I,J) = A(I,J)*DUM
        ENDDO
    END IF

ENDDO
IF (A(N,N) == 0.0_sp) then
        A(N,N) = TINY
end if
deallocate(vv)

END subroutine lu_descomposicion_sp

SUBROUTINE LU_RESOLUCION_sp(A,INDX,B)
!
! Performs backsubstitution to solve a linear system. 
! A is a matrix holding the LU decomposition in the standard manner.
! Indx is an array holding permutation indexes
! B is the right-hand-side vector, which is replaced by the solution
!
real(kind=sp), dimension(:,:) , intent(inout)   :: A
INTEGER, dimension(:), intent(in)                   :: INDX
real(kind=sp), dimension(:), intent(inout)      :: B

real(kind=sp)          :: SUM
INTEGER                    :: n, I,II,J,LL

n = size(a,dim=1)
if (size(a,dim=2) /= n) then
        STOP "Matrix is not square in lu_resolucion"
end if
if (size(indx) < n) then
        STOP "Size of Indx is too small in lu_resolucion"
end if
if (size(b) < n) then
        STOP "Size of B is too small in lu_resolucion"
end if

II = 0
DO I = 1,N
    LL = INDX(I)
    SUM = B(LL)
    B(LL) = B(I)
    IF (II /= 0) THEN
        DO J = II,I - 1
            SUM = SUM - A(I,J)*B(J)
        ENDDO
    ELSE IF (SUM /= 0.0_sp) THEN
        II = I
    END IF
    B(I) = SUM
ENDDO

DO I = N,1,-1
    SUM = B(I)
    IF (I < N) THEN
        DO J = I + 1,N
            SUM = SUM - A(I,J)*B(J)
        ENDDO
    END IF
    B(I) = SUM/A(I,I)
ENDDO

END subroutine lu_resolucion_sp

SUBROUTINE GAUSSJ_sp(A,B)
!
! Performs Gaussian elimination in matrix A to solve the
! multiple linear system Ax=B
!
! A is replaced by its inverse, and B by the solution vectors
!
real(kind=sp), intent(inout), dimension(:,:) ::  A
real(kind=sp), intent(inout), dimension(:,:) ::  B

real(kind=sp)           :: BIG,DUM,PIVINV
INTEGER        :: I,ICOL,IROW,J,K,L,LL
integer        :: n, m

INTEGER, dimension(:), allocatable ::  INDXC,INDXR,IPIV

n = size(a,dim=2)
if (size(a,dim=1) /= n) then
        STOP "A is not square in gaussj"
end if
if (size(b,dim=1) /= n) then
        STOP "Array mismatch in gaussj"
end if
m = size(b,dim=2)

allocate(ipiv(n),indxc(n),indxr(n))

IPIV(1:N) = 0

DO I = 1,N
    BIG = 0.0_sp
    DO J = 1,N
        IF (IPIV(J) /= 1) THEN
            DO K = 1,N
                IF (IPIV(K) == 0) THEN
                    IF (ABS(A(J,K)) >= BIG) THEN
                        BIG = ABS(A(J,K))
                        IROW = J
                        ICOL = K
                    END IF
                ELSE IF (IPIV(K) > 1) THEN
                    deallocate(ipiv,indxc,indxr)
                    STOP "Singular matrix in gaussj"
                END IF
            ENDDO
        END IF
    ENDDO

    IPIV(ICOL) = IPIV(ICOL) + 1
    IF (IROW /= ICOL) THEN
        DO L = 1,N
            DUM = A(IROW,L)
            A(IROW,L) = A(ICOL,L)
            A(ICOL,L) = DUM
        ENDDO
        DO L = 1,M
            DUM = B(IROW,L)
            B(IROW,L) = B(ICOL,L)
            B(ICOL,L) = DUM
        ENDDO
    END IF

    INDXR(I) = IROW
    INDXC(I) = ICOL
    IF (A(ICOL,ICOL) == 0.0_sp) THEN
       deallocate(ipiv,indxc,indxr)
       STOP "Singular matrix in gaussj"
    END IF

    PIVINV = 1.0_sp/A(ICOL,ICOL)
    A(ICOL,ICOL) = 1.0_sp
    A(ICOL,1:N) = A(ICOL,1:N) * PIVINV
    B(ICOL,1:M) = B(ICOL,1:M) * PIVINV

    DO LL = 1,N
        IF (LL /= ICOL) THEN
            DUM = A(LL,ICOL)
            A(LL,ICOL) = 0.0_sp
            A(LL,1:N) = A(LL,1:N) - A(ICOL,1:N)*DUM
            B(LL,1:M) = B(LL,1:M) - B(ICOL,1:M)*DUM
        END IF
    ENDDO
ENDDO

DO L = N,1,-1
    IF (INDXR(L) /= INDXC(L)) THEN
        DO K = 1,N
            DUM = A(K,INDXR(L))
            A(K,INDXR(L)) = A(K,INDXC(L))
            A(K,INDXC(L)) = DUM
        ENDDO
    END IF
ENDDO
deallocate(ipiv,indxc,indxr)

END SUBROUTINE GAUSSJ_sp

subroutine print_matrix_sp(a,unit)
real(kind=sp), intent(in), dimension(:,:)  :: a
integer, intent(in)               :: unit

integer  :: n, m, i, j
n = size(a,dim=1)
m = size(a,dim=2)
do i = 1, n
   do j = 1, m-1
      write(unit=unit,fmt="(f10.4)", advance="no") a(i,j)
   enddo
   write(unit=unit,fmt="(f10.4)") a(i,m)
enddo
write(unit=unit,fmt=*) "--------"
end subroutine print_matrix_sp

end module m_slineales_sp






module m_slineales_dp

use mcf_tipos
private
!
! Routines to solve systems of linear equations
!
public :: tridag, ctridag, lu_descomposicion, lu_resolucion, gaussj
public :: ludcmp, lubksb     ! Aliases
public :: print_matrix

interface tridag
   module procedure tridag_dp
end interface
private :: tridag_dp

interface ctridag
   module procedure ctridag_dp
end interface
private :: ctridag_dp

interface gaussj
   module procedure gaussj_dp
end interface
private :: gaussj_dp

interface ludcmp
   module procedure lu_descomposicion_dp
end interface

interface lubksb
   module procedure lu_resolucion_dp
end interface
!
interface lu_descomposicion
   module procedure lu_descomposicion_dp
end interface
private :: lu_descomposicion_dp

interface lu_resolucion
   module procedure lu_resolucion_dp
end interface
private :: lu_resolucion_dp

interface print_matrix
   module procedure print_matrix_dp
end interface
private :: print_matrix_dp

CONTAINS  !==================================================

SUBROUTINE TRIDAG_dp(A,D,C,B,X,n)
!
! Resuelve un sistema tridiagonal de ecuaciones lineales.
! D(1:n) contiene los elementos diagonales
! a(1:n) contiene los elementos infradiagonales (a(1) no se usa)
! c(1:n) contiene los elementos supradiagonales (c(n) no se usa)
! B(1:n) contiene los terminos independientes del sistema
!
! A la salida, x(1:n) contiene la solucion del sistema
!
!
integer, intent(in)              :: n
real(kind=dp), intent(in), dimension(:)   :: a, d, c, b
real(kind=dp), intent(out), dimension(:)  :: x

real(kind=dp)    ::  BET
INTEGER :: J

real(kind=dp), allocatable, dimension(:)  ::  GAM     ! Array de trabajo

IF (D(1) == 0.0_dp) THEN
    STOP "D(1)=0 in tridag"
END IF

allocate(gam(n))

BET = D(1)
X(1) = B(1)/BET
DO J = 2,N
    GAM(J) = C(J-1)/BET
    BET = D(J) - A(J)*GAM(J)
    IF (BET == 0.0_dp) THEN
        deallocate(gam)
        STOP "BET=0.0, Tridag failed."
    END IF

    X(J) = (B(J)-A(J)*X(J-1))/BET
ENDDO
DO J = N - 1,1,-1
    X(J) = X(J) - GAM(J+1)*X(J+1)
ENDDO

deallocate(gam)

END subroutine tridag_dp

SUBROUTINE ctridag_dp(A,D,C,B,X,n)
!
! Resuelve un sistema tridiagonal de ecuaciones lineales.
! D(1:n) contiene los elementos diagonales
! a(1:n) contiene los elementos infradiagonales (a(1) no se usa)
! c(1:n) contiene los elementos supradiagonales (c(n) no se usa)
! B(1:n) contiene los terminos independientes del sistema
!
! A la salida, x(1:n) contiene la solucion del sistema
!
!
integer, intent(in)                         :: n
complex(kind=dp), intent(in), dimension(:)  :: a, d, c, b
complex(kind=dp), intent(out), dimension(:) :: x

complex(kind=dp)                            :: BET
INTEGER                                     :: J

complex(kind=dp), allocatable, dimension(:) :: GAM     ! Array de trabajo

IF (D(1) == 0.0_dp) THEN
    STOP "D(1)=0 in Ctridag"
END IF

allocate(gam(n))

BET = D(1)
X(1) = B(1)/BET
DO J = 2,N
    GAM(J) = C(J-1)/BET
    BET = D(J) - A(J)*GAM(J)
    IF (BET == 0.0_dp) THEN
        deallocate(gam)
        STOP "BET=0.0, Ctridag failed."
    END IF

    X(J) = (B(J)-A(J)*X(J-1))/BET
ENDDO
DO J = N - 1,1,-1
    X(J) = X(J) - GAM(J+1)*X(J+1)
ENDDO

deallocate(gam)

END subroutine ctridag_dp
!===================================================================

SUBROUTINE LU_DESCOMPOSICION_dp(A,INDX,D)
!
! Performs the LU decomposition of a square matrix A
! indx is an array holding permutation indexes. It must
! be kept unchanged before calling lu_resolucion.
! d is the signature of the permutation.
! On exit, A contains the L and U matrices.
!
real(kind=dp), dimension(:,:) , intent(inout)   :: A
INTEGER, dimension(:), intent(out)                  :: INDX
real(kind=dp), intent(out)                      :: D

real(kind=dp)                     ::  AAMAX,DUM,SUM
INTEGER                               ::  n,I,IMAX,J,K

real(kind=dp), PARAMETER      :: TINY=1.0e-20_dp
                                  ! To avoid division by zero

real(kind=dp), dimension(:), allocatable ::  VV

n = size(a,dim=1)
if (size(a,dim=2) /= n) then
        STOP "Matrix is not square in lu_descomposicion"
end if
if (size(indx) < n) then
        STOP "Size of Indx is too small in lu_descomposicion"
end if

allocate(vv(n))

D = 1.0_dp
DO I = 1,N
    !! aamax = maxval(abs(a(i,1:n)))
    AAMAX = 0.0_dp
    DO J = 1,N
        IF (ABS(A(I,J)) > AAMAX) then
                AAMAX = ABS(A(I,J))
        end if
    ENDDO
    IF (AAMAX == 0.0_dp) THEN
        deallocate(vv)
        STOP "Singular matrix in lu_descomposicion"
    END IF
    VV(I) = 1.0_dp/AAMAX
ENDDO

DO J = 1,N
    IF (J > 1) THEN
        DO I = 1,J - 1
            SUM = A(I,J)
            IF (I > 1) THEN
                DO K = 1,I - 1
                    SUM = SUM - A(I,K)*A(K,J)
                ENDDO
                A(I,J) = SUM
            END IF
        ENDDO
    END IF

    AAMAX = 0.0_dp
    DO I = J,N
        SUM = A(I,J)
        IF (J > 1) THEN
            DO K = 1,J - 1
                SUM = SUM - A(I,K)*A(K,J)
            ENDDO
            A(I,J) = SUM
        END IF

        DUM = VV(I)*ABS(SUM)
        IF (DUM >= AAMAX) THEN
            IMAX = I
            AAMAX = DUM
        END IF
    ENDDO
    IF (J /= IMAX) THEN
        DO K = 1,N
            DUM = A(IMAX,K)
            A(IMAX,K) = A(J,K)
            A(J,K) = DUM
        ENDDO
        D = -D
        VV(IMAX) = VV(J)
    END IF

    INDX(J) = IMAX
    IF (J /= N) THEN
        IF (A(J,J) == 0.0_dp) then
                A(J,J) = TINY
        end if
        DUM = 1.0_dp/A(J,J)
        DO I = J + 1,N
            A(I,J) = A(I,J)*DUM
        ENDDO
    END IF

ENDDO
IF (A(N,N) == 0.0_dp) then
        A(N,N) = TINY
end if
deallocate(vv)

END subroutine lu_descomposicion_dp

SUBROUTINE LU_RESOLUCION_dp(A,INDX,B)
!
! Performs backsubstitution to solve a linear system. 
! A is a matrix holding the LU decomposition in the standard manner.
! Indx is an array holding permutation indexes
! B is the right-hand-side vector, which is replaced by the solution
!
real(kind=dp), dimension(:,:) , intent(inout)   :: A
INTEGER, dimension(:), intent(in)                   :: INDX
real(kind=dp), dimension(:), intent(inout)      :: B

real(kind=dp)          :: SUM
INTEGER                    :: n, I,II,J,LL

n = size(a,dim=1)
if (size(a,dim=2) /= n) then
        STOP "Matrix is not square in lu_resolucion"
end if
if (size(indx) < n) then
        STOP "Size of Indx is too small in lu_resolucion"
end if
if (size(b) < n) then
        STOP "Size of B is too small in lu_resolucion"
end if

II = 0
DO I = 1,N
    LL = INDX(I)
    SUM = B(LL)
    B(LL) = B(I)
    IF (II /= 0) THEN
        DO J = II,I - 1
            SUM = SUM - A(I,J)*B(J)
        ENDDO
    ELSE IF (SUM /= 0.0_dp) THEN
        II = I
    END IF
    B(I) = SUM
ENDDO

DO I = N,1,-1
    SUM = B(I)
    IF (I < N) THEN
        DO J = I + 1,N
            SUM = SUM - A(I,J)*B(J)
        ENDDO
    END IF
    B(I) = SUM/A(I,I)
ENDDO

END subroutine lu_resolucion_dp

SUBROUTINE GAUSSJ_dp(A,B)
!
! Performs Gaussian elimination in matrix A to solve the
! multiple linear system Ax=B
!
! A is replaced by its inverse, and B by the solution vectors
!
real(kind=dp), intent(inout), dimension(:,:) ::  A
real(kind=dp), intent(inout), dimension(:,:) ::  B

real(kind=dp)           :: BIG,DUM,PIVINV
INTEGER        :: I,ICOL,IROW,J,K,L,LL
integer        :: n, m

INTEGER, dimension(:), allocatable ::  INDXC,INDXR,IPIV

n = size(a,dim=2)
if (size(a,dim=1) /= n) then
        STOP "A is not square in gaussj"
end if
if (size(b,dim=1) /= n) then
        STOP "Array mismatch in gaussj"
end if

m = size(b,dim=2)

allocate(ipiv(n),indxc(n),indxr(n))

IPIV(1:N) = 0

DO I = 1,N
    BIG = 0.0_dp
    DO J = 1,N
        IF (IPIV(J) /= 1) THEN
            DO K = 1,N
                IF (IPIV(K) == 0) THEN
                    IF (ABS(A(J,K)) >= BIG) THEN
                        BIG = ABS(A(J,K))
                        IROW = J
                        ICOL = K
                    END IF
                ELSE IF (IPIV(K) > 1) THEN
                    deallocate(ipiv,indxc,indxr)
                    STOP "Singular matrix in gaussj"
                END IF
            ENDDO
        END IF
    ENDDO

    IPIV(ICOL) = IPIV(ICOL) + 1
    IF (IROW /= ICOL) THEN
        DO L = 1,N
            DUM = A(IROW,L)
            A(IROW,L) = A(ICOL,L)
            A(ICOL,L) = DUM
        ENDDO
        DO L = 1,M
            DUM = B(IROW,L)
            B(IROW,L) = B(ICOL,L)
            B(ICOL,L) = DUM
        ENDDO
    END IF

    INDXR(I) = IROW
    INDXC(I) = ICOL
    IF (A(ICOL,ICOL) == 0.0_dp) THEN
       deallocate(ipiv,indxc,indxr)
       STOP "Singular matrix in gaussj"
    END IF

    PIVINV = 1.0_dp/A(ICOL,ICOL)
    A(ICOL,ICOL) = 1.0_dp
    A(ICOL,1:N) = A(ICOL,1:N) * PIVINV
    B(ICOL,1:M) = B(ICOL,1:M) * PIVINV

    DO LL = 1,N
        IF (LL /= ICOL) THEN
            DUM = A(LL,ICOL)
            A(LL,ICOL) = 0.0_dp
            A(LL,1:N) = A(LL,1:N) - A(ICOL,1:N)*DUM
            B(LL,1:M) = B(LL,1:M) - B(ICOL,1:M)*DUM
        END IF
    ENDDO
ENDDO

DO L = N,1,-1
    IF (INDXR(L) /= INDXC(L)) THEN
        DO K = 1,N
            DUM = A(K,INDXR(L))
            A(K,INDXR(L)) = A(K,INDXC(L))
            A(K,INDXC(L)) = DUM
        ENDDO
    END IF
ENDDO
deallocate(ipiv,indxc,indxr)

END SUBROUTINE GAUSSJ_dp

subroutine print_matrix_dp(a,unit)
real(kind=dp), intent(in), dimension(:,:)  :: a
integer, intent(in)               :: unit

integer  :: n, m, i, j
n = size(a,dim=1)
m = size(a,dim=2)
do i = 1, n
   do j = 1, m-1
      write(unit=unit,fmt="(f10.4)", advance="no") a(i,j)
   enddo
   write(unit=unit,fmt="(f10.4)") a(i,m)
enddo
write(unit=unit,fmt=*) "--------"
end subroutine print_matrix_dp

end module m_slineales_dp






module m_diagonalizacion_sp
  !
  ! Este modulo contiene rutinas para diagonalizar matrices simetricas y no simetricas
  ! (En versiones 'single' y 'double', pero con el mismo nombre abreviado)
  ! En el caso de matrices no simetrícas sólo se calculan los valores propios.
  !
use mcf_tipos
private

public :: tred2   ! Convierte una matriz simetrica a tridiagonal
interface tred2
  module procedure tred2_sp
end interface

public :: tqli    ! Diagonaliza una matriz tridiagonal
interface tqli
  module procedure tqli_sp
end interface

public :: eigsrt  ! Ordena autovalores 
interface eigsrt
  module procedure eigsrt_sp
end interface

public :: eival
interface eival
  module procedure eival_sp
end interface
!--------------------------------------------------------------------
!--------------------------------------------------------------------
private :: pythag
interface pythag
  module procedure pythag_sp
end interface
private :: balanc
interface balanc
  module procedure balanc_sp
end interface
private :: elmhes
interface elmhes
  module procedure elmhes_sp
end interface
private :: hqr
interface hqr
  module procedure hqr_sp
end interface
private :: tred2_sp
private :: tqli_sp
private :: pythag_sp
private :: eigsrt_sp
private :: balanc_sp
private :: elmhes_sp
private :: hqr_sp
private :: eival_sp

CONTAINS !===========================

SUBROUTINE TQLI_sp(D,E,Z)
  ! Calcula los autovalores de una matriz simetrica tridiagonal
  ! D(1:n) es la diagonal
  ! E(1:n) (e(1) no se usa) son los valores infra- y supradiagonales
  ! Si Z aparece como argumento, se calculan tambien los autovectores.
  !  (En este caso Z ha de ser la matriz unidad a la entrada)
  !
real(kind=sp), dimension(:), intent(inout) :: D,E
real(kind=sp), dimension(:,:), intent(inout), OPTIONAL :: Z

real(kind=sp)       ::   B,C,DD,F,G,P,R,S
INTEGER    ::   I,ITER,K,L,M
integer    ::   n
logical    ::   premature_exit

integer, parameter :: maxit = 30

n = size(d)
if (size(e) < n) then
        STOP "E too short in tqli"
end if
if (present(z)) then
   if (size(z,dim=1) < n) then
        STOP "1st dim of z too short in tqli"
   end if
   if (size(z,dim=2) < n) then
        STOP "2nd dim of z too short in tqli"
   end if
endif

IF (N <= 1) then
        RETURN
end if

DO I = 2,N
   E(I-1) = E(I)
ENDDO
E(N) = 0.0_sp

DO L = 1,N
    ITER = 0

 loop:    DO      ! iterations
         
        premature_exit = .false.    ! Kludge to work around horrible goto
        DO  M = L,N - 1
            DD = ABS(D(M)) + ABS(D(M+1))
            IF (ABS(E(M))+DD == DD) then
               premature_exit = .true.
               exit         ! of  "do m" loop
            endif
        ENDDO
        if (.not. premature_exit) then
                M = N
        end if

        IF (M == L) then
                EXIT
        end if

        IF (ITER == MAXIT) then
                STOP "too many iterations in TQLI"
        end if
        ITER = ITER + 1
        G = (D(L+1)-D(L))/ (2.0_sp*E(L))
        R = pythag(G,1.0_sp)
        G = D(M) - D(L) + E(L)/ (G+SIGN(R,G))
        S = 1.0_sp
        C = 1.0_sp
        P = 0.0_sp
        DO I = M - 1,L,-1
                F = S*E(I)
                B = C*E(I)
                r = pythag(f,g)
                e(i+1) = r
                ! Recover from underflow
                if (r == 0.0_sp) then
                   d(i+1) = d(i+1) - p
                   e(m) = 0.0_sp
                   cycle loop
                endif
                s = f/r
                c = g/r
                G = D(I+1) - P
                R = (D(I)-G)*S + 2.0*C*B
                P = S*R
                D(I+1) = G + P
                G = C*R - B
                if (present(z)) then  !------------ Autovectores
                   DO K = 1,N
                      F = Z(K,I+1)
                      Z(K,I+1) = S*Z(K,I) + C*F
                      Z(K,I) = C*Z(K,I) - S*F
                   ENDDO
                endif                 !------------
        ENDDO
        D(L) = D(L) - P
        E(L) = G
        E(M) = 0.0_sp
     ENDDO  loop !------------ outer loop

   ENDDO      ! loop over L


END subroutine tqli_sp


!
function pythag_sp(a,b) result(res)
real(kind=sp), intent(in)  :: a, b
real(kind=sp)              :: res
!
! Computes sqrt(a*a+b*b) without overflow or underflow
!
real(kind=sp)    :: absa, absb

absa=abs(a)
absb=abs(b)
if (absa  > absb) then
   res = absa*sqrt(1.0_sp + (absb/absa)**2)
else
   if (absb == 0.0_sp) then
      res = 0.0_sp
   else
      res = absb*sqrt(1.0_sp + (absa/absb)**2)
   endif
endif
end function pythag_sp


!SUBROUTINE EIGSRT_sp(D,V,N,NP)
SUBROUTINE EIGSRT_sp(D,V)

!INTEGER, intent(in) ::  N,NP

real(kind=sp), intent(inout), dimension(:)             ::  D
real(kind=sp), intent(inout), dimension(:,:), optional ::  V

real(kind=sp)        :: P
INTEGER              :: I,J,K,n

n = size(d)

DO I = 1,N - 1
    K = I
    P = D(I)
    DO J = I + 1,N
        IF (D(J) <= P) THEN
            K = J
            P = D(J)
        END IF
    enddo
    IF (K /= I) THEN
        D(K) = D(I)
        D(I) = P
        if (present(v)) then
            DO J = 1,N
                P = V(J,I)
                V(J,I) = V(J,K)
                V(J,K) = P
            enddo
        end if
    END IF
ENDDO

END subroutine eigsrt_sp

SUBROUTINE TRED2_sp(A,D,E)
! Reduce una matriz simetrica a tridiagonal

real(kind=sp), intent(inout), dimension(:,:)  :: a
real(kind=sp), intent(out), dimension(:)      :: d, e



real(kind=sp)      ::  F,G,H,HH,SCALE
INTEGER            ::  I,J,K,L, N


n = size(a,dim=2)
if (size(a,dim=1) /= n) then
        stop "A is not square"
end if
if (size(d) < n) then
        stop "d is not big enough in tred2"
end if
if (size(e) < n) then
        stop "e is not big enough in tred2"
end if

IF (N > 1) THEN
    DO I = N,2,-1
        L = I - 1
        H = 0.0_sp
        SCALE = 0.0_sp
        IF (L > 1) THEN
            DO K = 1,L
                SCALE = SCALE + ABS(A(I,K))
            ENDDO
            IF (SCALE == 0.0_sp) THEN
                E(I) = A(I,L)

            ELSE
                DO K = 1,L
                    A(I,K) = A(I,K)/SCALE
                    H = H + A(I,K)**2
                ENDDO
                F = A(I,L)
                G = -SIGN(SQRT(H),F)
                E(I) = SCALE*G
                H = H - F*G
                A(I,L) = F - G
                F = 0.0_sp
                DO J = 1,L
                    A(J,I) = A(I,J)/H
                    G = 0.0_sp
                    DO K = 1,J
                        G = G + A(J,K)*A(I,K)
                    ENDDO
                    IF (L > J) THEN
                        DO K = J + 1,L
                            G = G + A(K,J)*A(I,K)
                        ENDDO
                    END IF
                    E(J) = G/H
                    F = F + E(J)*A(I,J)
                ENDDO
                HH = F/ (H+H)
                DO J = 1,L
                    F = A(I,J)
                    G = E(J) - HH*F
                    E(J) = G
                    DO K = 1,J
                        A(J,K) = A(J,K) - F*E(K) - G*A(I,K)
                    ENDDO
                ENDDO
            END IF

        ELSE
            E(I) = A(I,L)
        END IF

        D(I) = H
  ENDDO
END IF

D(1) = 0.0_sp
E(1) = 0.0_sp
DO I = 1,N
    L = I - 1
    IF (D(I) /= 0.0_sp) THEN
        DO J = 1,L
            G = 0.0_sp
            DO K = 1,L
                G = G + A(I,K)*A(K,J)
            ENDDO
            DO K = 1,L
                A(K,J) = A(K,J) - G*A(K,I)
            ENDDO
        ENDDO
    END IF

    D(I) = A(I,I)
    A(I,I) = 1.0_sp
    IF (L >= 1) THEN
        DO J = 1,L
            A(I,J) = 0.0_sp
            A(J,I) = 0.0_sp
        ENDDO
    END IF

ENDDO

END subroutine tred2_sp

subroutine eival_sp (a,wr,wi)

real(kind=sp), intent(inout), dimension(:,:) :: a
real(kind=sp), intent(out), dimension(:)     :: wr,wi
integer :: n,i

n = size(a,dim=2)

if (size(a,dim=1) /= n) then

    stop "A is not square"

end if

if (size(wi) < n) then

        stop "wi is not big enough in eival"

end if

if (size(wr) < n) then

        stop "wr is not big enough in eival"

end if


! Matrix balancing

call balanc(a)

! Hessenberg form of the matrix

call elmhes(a)

! Eigenvalues

call hqr(a,wr,wi)

! Now we can order the eigenvalues using eigsrt(D,V). 
! D are the eigenvalue moduli. In V the first and the second
! column are the real and imaginary parts of the eigen values.
! Dirty. V is A. Also wr but it and wi are finally recovered in 
! the desired order. The matrix A is completely destroyed.

a      = 0.0_sp
a(1,:) = wr
a(2,:) = wi

do i=1,n

    wr(i) = sqrt(wr(i)*wr(i) + wi(i)*wi(i))

end do

call eigsrt(wr,a)

wr = a(1,:)
wi = a(2,:)

end subroutine eival_sp

subroutine balanc_sp(a)

real(kind=sp), intent(inout), dimension(:,:) :: a
real(kind=sp), parameter                     :: RADIX=2.0_sp,SQRDX=RADIX**2.0_sp
integer                                      :: i,j,last,n
real(kind=sp)                                :: c,f,g,r,s

n = size(a,dim=2)
if (size(a,dim=1) /= n) then
    stop "A is not square"
end if

do
    last=1
    do i=1,n

        c=0.0_sp
        r=0.0_sp

        do j=1,n

            if(j /= i)then

                c=c+abs(a(j,i))
                r=r+abs(a(i,j))

            end if

        end do

        if(c /= 0.0_sp .and. r /= 0.0_sp)then

            g=r/RADIX
            f=1.0_sp
            s=c+r

            do

                if(c < g)then

                    f=f*RADIX
                    c=c*SQRDX

                else

                    exit

                end if

            end do

            g=r*RADIX

            do

                if(c > g)then

                    f=f/RADIX
                    c=c/SQRDX
          
                else
              
                    exit

                end if

            end do

                if((c+r)/f < 0.95_sp*s)then

                    last=0
                    g=1.0_sp/f

                do j=1,n

                    a(i,j)=a(i,j)*g
              
                end do
              
                do j=1,n
                  
                    a(j,i)=a(j,i)*f

                end do

            end if

        end if

    end do
      
    if(last /= 0) then
      
        exit

    end if

end do

end subroutine balanc_sp
!  (C) Copr. 1986-92 Numerical Recipes Software Y].

subroutine elmhes_sp(a)

real(kind=sp), intent(inout), dimension(:,:) :: a
integer                                      :: i,j,m,n
real                                         :: x,y

n = size(a,dim=2)

if (size(a,dim=1) /= n) then

    stop "A is not square"

end if

do m=2,n-1
        
    x=0.0_sp
    i=m

    do j=m,n

        if(abs(a(j,m-1)) > abs(x))then

            x=a(j,m-1)
            i=j

        end if

    end do

    if(i /= m)then

        do j=m-1,n

            y=a(i,j)
            a(i,j)=a(m,j)
            a(m,j)=y

        end do

        do j=1,n

            y=a(j,i)
            a(j,i)=a(j,m)
            a(j,m)=y

        end do

    end if

    if(x /= 0.0_sp)then

        do i=m+1,n

            y=a(i,m-1)

            if(y /= 0.0_sp)then

                y=y/x
                a(i,m-1)=y

                do j=m,n

                    a(i,j)=a(i,j)-y*a(m,j)

                end do

                do j=1,n

                    a(j,m)=a(j,m)+y*a(j,i)

                end do

            end if

        end do

    end if

end do

end subroutine elmhes_sp      
!  (C) Copr. 1986-92 Numerical Recipes Software Y].

subroutine hqr_sp(a,wr,wi)


real(kind=sp), intent(inout), dimension(:,:) :: a
real(kind=sp), intent(out), dimension(:)     :: wi,wr
integer                                      :: i,its,j,k,l,m,nn,n
real(kind=sp)                                :: anorm,p,q,r,s,t,u,v,w,x,y,z
logical                                      :: premature_exit

n = size(a,dim=2)

if (size(a,dim=1) /= n) then

    stop "A is not square"

end if

if (size(wi) < n) then

        stop "wi is not big enough in hqr"

end if

if (size(wr) < n) then

        stop "wr is not big enough in hqr"

end if

anorm=abs(a(1,1))

do i=2,n

    do j=i-1,n

        anorm=anorm+abs(a(i,j))

    end do

end do

nn=n
t=0.0_sp

do_1: do

    if(nn >= 1)then

        its=0

do_2:   do

            premature_exit = .false. 
   
do_l:       do l=nn,2,-1

                s=abs(a(l-1,l-1))+abs(a(l,l))

                if(s == 0.0_sp) then

                    s=anorm

                end if

                if(abs(a(l,l-1))+s == s)then

                    premature_exit = .true.
                    exit do_l

                end if

            end do do_l

            if (.not. premature_exit) then

                l=1

            end if

            x=a(nn,nn)

            if(l == nn)then

                wr(nn)=x+t
                wi(nn)=0.0_sp
                nn=nn-1

            else
   
                y=a(nn-1,nn-1)
                w=a(nn,nn-1)*a(nn-1,nn)

                if(l == nn-1)then
           
                    p=0.5_sp*(y-x)
                    q=p**2.0_sp+w
                    z=sqrt(abs(q))
                    x=x+t
            
                    if(q >= 0.0_sp)then
                
                        z=p+sign(z,p)
                        wr(nn)=x+z
                        wr(nn-1)=wr(nn)
 
                        if(z /= 0.0_sp)then
     
                            wr(nn)=x-w/z

                        end if
        
                        wi(nn)=0.0_sp
                        wi(nn-1)=0.0_sp
 
                    else
               
                        wr(nn)=x+p
                        wr(nn-1)=wr(nn)
                        wi(nn)=z
                        wi(nn-1)=-z
            
                    end if
            
                    nn=nn-2
          
                else
        
                    if(its == 30)then
  
                        stop "too many iterations in hqr"
         
                    end if

                    if(its == 10 .or. its == 20)then
      
                        t=t+x
              
                        do i=1,nn
                
                            a(i,i)=a(i,i)-x
              
                        end do

                        s=abs(a(nn,nn-1))+abs(a(nn-1,nn-2))
                        x=0.75_sp*s
                        y=x
                        w=-0.4375_sp*s**2.0_sp

                    end if
      
                    its=its+1
            
do_m:               do m=nn-2,l,-1
              
                        z=a(m,m)
                        r=x-z
                        s=y-z
                        p=(r*s-w)/a(m+1,m)+a(m,m+1)
                        q=a(m+1,m+1)-z-r-s
                        r=a(m+2,m+1)
                        s=abs(p)+abs(q)+abs(r)
                        p=p/s
                        q=q/s
                        r=r/s
             
                        if(m == l) then
                  
                            exit do_m

                        end if
          
                        u=abs(a(m,m-1))*(abs(q)+abs(r))
                        v=abs(p)*(abs(a(m-1,m-1))+abs(z)+abs(a(m+1,m+1)))
              
                        if(u+v == v)then
                  
                            exit do_m
             
                        end if

                    end do do_m

                    do i=m+2,nn
              
                        a(i,i-2)=0.0_sp
                    
                        if (i /= m+2) then
                  
                            a(i,i-3)=0.0_sp
              
                        end if
            
                    end do
            
                    do k=m,nn-1
              
                        if(k /= m)then
                
                            p=a(k,k-1)
                            q=a(k+1,k-1)
                            r=0.0_sp
              
                            if(k /= nn-1)then
                    
                                r=a(k+2,k-1)
                
                            end if
               
                            x=abs(p)+abs(q)+abs(r)
               
                            if(x /= 0.0_sp)then
                  
                                p=p/x
                                q=q/x
                                r=r/x
                
                            end if
              
                        end if

                        s=sign(sqrt(p**2.0_sp+q**2.0_sp+r**2.0_sp),p)
        
                        if(s /= 0.0_sp)then
                
                            if(k == m)then
                  
                                if(l /= m)then
                      
                                    a(k,k-1)=-a(k,k-1)
                  
                                end if
               
                            else
                  
                                a(k,k-1)=-s*x
                
                            end if
                
                            p=p+s
                            x=p/s
                            y=q/s
                            z=r/s
                            q=q/p
                            r=r/p
                
                            do j=k,nn
                  
                                p=a(k,j)+q*a(k+1,j)
                  
                                if(k /= nn-1)then
                    
                                    p=p+r*a(k+2,j)
                                    a(k+2,j)=a(k+2,j)-p*z
                  
                                end if
                
                                a(k+1,j)=a(k+1,j)-p*y
                                a(k,j)=a(k,j)-p*x
                
                            end do
                
                            do i=l,min(nn,k+3)
                   
                                p=x*a(i,k)+y*a(i,k+1)
                  
                                if(k /= nn-1)then
                 
                                    p=p+z*a(i,k+2)
                                    a(i,k+2)=a(i,k+2)-p*r
                  
                                end if
                  
                                a(i,k+1)=a(i,k+1)-p*q
                                a(i,k)=a(i,k)-p
                
                            end do
                        
                        end if

                    end do

                    cycle do_2

                end if

            end if

            exit do_2

            end do do_2
      
    else

            exit do_1
   
    end if

end do do_1

end subroutine hqr_sp
!  (C) Copr. 1986-92 Numerical Recipes Software Y].

end module m_diagonalizacion_sp






module m_diagonalizacion_dp
  !
  ! Este modulo contiene rutinas para diagonalizar matrices simetricas
  ! (En versiones 'single' y 'double', pero con el mismo nombre abreviado)
  !
use mcf_tipos
private

public :: tred2   ! Convierte una matriz simetrica a tridiagonal
interface tred2
  module procedure tred2_dp
end interface

public :: tqli    ! Diagonaliza una matriz tridiagonal
interface tqli
  module procedure tqli_dp
end interface

public :: eigsrt  ! Ordena autovalores 
interface eigsrt
  module procedure eigsrt_dp
end interface

public :: eival
interface eival
  module procedure eival_dp
end interface

!--------------------------------------------------------------------
!--------------------------------------------------------------------
private :: pythag
interface pythag
  module procedure pythag_dp
end interface
private :: balanc
interface balanc
  module procedure balanc_dp
end interface
private :: elmhes
interface elmhes
  module procedure elmhes_dp
end interface
private :: hqr
interface hqr
  module procedure hqr_dp
end interface

private :: tred2_dp
private :: tqli_dp
private :: pythag_dp
private :: eigsrt_dp
private :: balanc_dp
private :: elmhes_dp
private :: hqr_dp
private :: eival_dp
!--------------------------------------------------------------------
!--------------------------------------------------------------------

CONTAINS !===========================

SUBROUTINE TQLI_dp(D,E,Z)
  ! Calcula los autovalores de una matriz simetrica tridiagonal
  ! D(1:n) es la diagonal
  ! E(1:n) (e(1) no se usa) son los valores infra- y supradiagonales
  ! Si Z aparece como argumento, se calculan tambien los autovectores.
  !  (En este caso Z ha de ser la matriz unidad a la entrada)
  !
real(kind=dp), dimension(:), intent(inout) :: D,E
real(kind=dp), dimension(:,:), intent(inout), OPTIONAL :: Z

real(kind=dp)       ::   B,C,DD,F,G,P,R,S
INTEGER    ::   I,ITER,K,L,M
integer    ::   n
logical    ::   premature_exit

integer, parameter :: maxit = 30

n = size(d)
if (size(e) < n) then
        STOP "E too short in tqli"
end if
if (present(z)) then
   if (size(z,dim=1) < n) then
        STOP "1st dim of z too short in tqli"
   end if
   if (size(z,dim=2) < n) then
        STOP "2nd dim of z too short in tqli"
   end if
endif

IF (N <= 1) then
        RETURN
end if

DO I = 2,N
   E(I-1) = E(I)
ENDDO
E(N) = 0.0_dp

DO L = 1,N
    ITER = 0

 loop:    DO      ! iterations
         
        premature_exit = .false.    ! Kludge to work around horrible goto
        DO  M = L,N - 1
            DD = ABS(D(M)) + ABS(D(M+1))
            IF (ABS(E(M))+DD == DD) then
               premature_exit = .true.
               exit         ! of  "do m" loop
            endif
        ENDDO
        if (.not. premature_exit) then
                M = N
        end if

        IF (M == L) then
                EXIT
        end if

        IF (ITER == MAXIT) then
                STOP "too many iterations in TQLI"
        end if
        ITER = ITER + 1
        G = (D(L+1)-D(L))/ (2.0_dp*E(L))
        R = pythag(G,1.0_dp)
        G = D(M) - D(L) + E(L)/ (G+SIGN(R,G))
        S = 1.0_dp
        C = 1.0_dp
        P = 0.0_dp
        DO I = M - 1,L,-1
                F = S*E(I)
                B = C*E(I)
                r = pythag(f,g)
                e(i+1) = r
                ! Recover from underflow
                if (r == 0.0_dp) then
                   d(i+1) = d(i+1) - p
                   e(m) = 0.0_dp
                   cycle loop
                endif
                s = f/r
                c = g/r
                G = D(I+1) - P
                R = (D(I)-G)*S + 2.0*C*B
                P = S*R
                D(I+1) = G + P
                G = C*R - B
                if (present(z)) then  !------------ Autovectores
                   DO K = 1,N
                      F = Z(K,I+1)
                      Z(K,I+1) = S*Z(K,I) + C*F
                      Z(K,I) = C*Z(K,I) - S*F
                   ENDDO
                endif                 !------------
        ENDDO
        D(L) = D(L) - P
        E(L) = G
        E(M) = 0.0_dp
     ENDDO  loop !------------ outer loop

   ENDDO      ! loop over L


END subroutine tqli_dp


!
function pythag_dp(a,b) result(res)
real(kind=dp), intent(in)  :: a, b
real(kind=dp)              :: res
!
! Computes sqrt(a*a+b*b) without overflow or underflow
!
real(kind=dp)    :: absa, absb

absa=abs(a)
absb=abs(b)
if (absa  > absb) then
   res = absa*sqrt(1.0_dp + (absb/absa)**2)
else
   if (absb == 0.0_dp) then
      res = 0.0_dp
   else
      res = absb*sqrt(1.0_dp + (absa/absb)**2)
   endif
endif
end function pythag_dp


!SUBROUTINE EIGSRT_dp(D,V,N,NP)
SUBROUTINE EIGSRT_dp(D,V)

!INTEGER, intent(in) ::  N,NP

real(kind=dp), intent(inout), dimension(:)             ::  D
real(kind=dp), intent(inout), dimension(:,:), optional ::  V

real(kind=dp)    :: P
INTEGER              :: I,J,K,n

n = size(d)

DO I = 1,N - 1
    K = I
    P = D(I)
    DO J = I + 1,N
        IF (D(J) <= P) THEN
            K = J
            P = D(J)
        END IF
    enddo
    IF (K /= I) THEN
        D(K) = D(I)
        D(I) = P
        if (present(v)) then
            DO J = 1,N
                P = V(J,I)
                V(J,I) = V(J,K)
                V(J,K) = P
            enddo
        end if
    END IF
ENDDO

END subroutine eigsrt_dp

SUBROUTINE TRED2_dp(A,D,E)
! Reduce una matriz simetrica a tridiagonal

real(kind=dp), intent(inout), dimension(:,:)  :: a
real(kind=dp), intent(out), dimension(:)      :: d, e



real(kind=dp)      ::  F,G,H,HH,SCALE
INTEGER            ::  I,J,K,L, N


n = size(a,dim=2)
if (size(a,dim=1) /= n) then
        stop "A is not square"
end if
if (size(d) < n) then
        stop "d is not big enough in tred2"
end if
if (size(e) < n) then
        stop "e is not big enough in tred2"
end if

IF (N > 1) THEN
    DO I = N,2,-1
        L = I - 1
        H = 0.0_dp
        SCALE = 0.0_dp
        IF (L > 1) THEN
            DO K = 1,L
                SCALE = SCALE + ABS(A(I,K))
            ENDDO
            IF (SCALE == 0.0_dp) THEN
                E(I) = A(I,L)

            ELSE
                DO K = 1,L
                    A(I,K) = A(I,K)/SCALE
                    H = H + A(I,K)**2
                ENDDO
                F = A(I,L)
                G = -SIGN(SQRT(H),F)
                E(I) = SCALE*G
                H = H - F*G
                A(I,L) = F - G
                F = 0.0_dp
                DO J = 1,L
                    A(J,I) = A(I,J)/H
                    G = 0.0_dp
                    DO K = 1,J
                        G = G + A(J,K)*A(I,K)
                    ENDDO
                    IF (L > J) THEN
                        DO K = J + 1,L
                            G = G + A(K,J)*A(I,K)
                        ENDDO
                    END IF
                    E(J) = G/H
                    F = F + E(J)*A(I,J)
                ENDDO
                HH = F/ (H+H)
                DO J = 1,L
                    F = A(I,J)
                    G = E(J) - HH*F
                    E(J) = G
                    DO K = 1,J
                        A(J,K) = A(J,K) - F*E(K) - G*A(I,K)
                    ENDDO
                ENDDO
            END IF

        ELSE
            E(I) = A(I,L)
        END IF

        D(I) = H
  ENDDO
END IF

D(1) = 0.0_dp
E(1) = 0.0_dp
DO I = 1,N
    L = I - 1
    IF (D(I) /= 0.0_dp) THEN
        DO J = 1,L
            G = 0.0_dp
            DO K = 1,L
                G = G + A(I,K)*A(K,J)
            ENDDO
            DO K = 1,L
                A(K,J) = A(K,J) - G*A(K,I)
            ENDDO
        ENDDO
    END IF

    D(I) = A(I,I)
    A(I,I) = 1.0_dp
    IF (L >= 1) THEN
        DO J = 1,L
            A(I,J) = 0.0_dp
            A(J,I) = 0.0_dp
        ENDDO
    END IF

ENDDO

END subroutine tred2_dp

subroutine eival_dp (a,wr,wi)

real(kind=dp), intent(inout), dimension(:,:) :: a
real(kind=dp), intent(out), dimension(:)     :: wr,wi
integer :: n,i

n = size(a,dim=2)

if (size(a,dim=1) /= n) then

    stop "A is not square"

end if

if (size(wi) < n) then

        stop "wi is not big enough in eival"

end if

if (size(wr) < n) then

        stop "wr is not big enough in eival"

end if


! Matrix balancing

call balanc(a)

! Hessenberg form of the matrix

call elmhes(a)

! Eigenvalues

call hqr(a,wr,wi)

! Now we can order the eigenvalues using eigsrt(D,V). 
! D are the eigenvalue moduli. In V the first and the second
! column are the real and imaginary parts of the eigen values.
! Dirty. V is A. Also wr but it and wi are finally recovered in 
! the desired order. The matrix A is completely destroyed.

a      = 0.0_dp
a(1,:) = wr
a(2,:) = wi

do i=1,n

    wr(i) = sqrt(wr(i)*wr(i) + wi(i)*wi(i))

end do

call eigsrt(wr,a)

wr = a(1,:)
wi = a(2,:)

end subroutine eival_dp

subroutine balanc_dp(a)

real(kind=dp), intent(inout), dimension(:,:) :: a
real(kind=dp), parameter                     :: RADIX=2.0_dp,SQRDX=RADIX**2.0_dp
integer                                      :: i,j,last,n
real(kind=dp)                                :: c,f,g,r,s

n = size(a,dim=2)
if (size(a,dim=1) /= n) then
    stop "A is not square"
end if

do
    last=1
    do i=1,n

        c=0.0_dp
        r=0.0_dp

        do j=1,n

            if(j /= i)then

                c=c+abs(a(j,i))
                r=r+abs(a(i,j))

            end if

        end do

        if(c /= 0.0_dp .and. r /= 0.0_dp)then

            g=r/RADIX
            f=1.0_dp
            s=c+r

            do

                if(c < g)then

                    f=f*RADIX
                    c=c*SQRDX

                else

                    exit

                end if

            end do

            g=r*RADIX

            do

                if(c > g)then

                    f=f/RADIX
                    c=c/SQRDX
          
                else
              
                    exit

                end if

            end do

                if((c+r)/f < 0.95_dp*s)then

                    last=0
                    g=1.0_dp/f

                do j=1,n

                    a(i,j)=a(i,j)*g
              
                end do
              
                do j=1,n
                  
                    a(j,i)=a(j,i)*f

                end do

            end if

        end if

    end do
      
    if(last /= 0) then
      
        exit

    end if

end do

end subroutine balanc_dp
!  (C) Copr. 1986-92 Numerical Recipes Software Y].

subroutine elmhes_dp(a)

real(kind=dp), intent(inout), dimension(:,:) :: a
integer                                      :: i,j,m,n
real                                         :: x,y

n = size(a,dim=2)

if (size(a,dim=1) /= n) then

    stop "A is not square"

end if

do m=2,n-1
        
    x=0.0_dp
    i=m

    do j=m,n

        if(abs(a(j,m-1)) > abs(x))then

            x=a(j,m-1)
            i=j

        end if

    end do

    if(i /= m)then

        do j=m-1,n

            y=a(i,j)
            a(i,j)=a(m,j)
            a(m,j)=y

        end do

        do j=1,n

            y=a(j,i)
            a(j,i)=a(j,m)
            a(j,m)=y

        end do

    end if

    if(x /= 0.0_dp)then

        do i=m+1,n

            y=a(i,m-1)

            if(y /= 0.0_dp)then

                y=y/x
                a(i,m-1)=y

                do j=m,n

                    a(i,j)=a(i,j)-y*a(m,j)

                end do

                do j=1,n

                    a(j,m)=a(j,m)+y*a(j,i)

                end do

            end if

        end do

    end if

end do

end subroutine elmhes_dp      
!  (C) Copr. 1986-92 Numerical Recipes Software Y].

subroutine hqr_dp(a,wr,wi)


real(kind=dp), intent(inout), dimension(:,:) :: a
real(kind=dp), intent(out), dimension(:)     :: wi,wr
integer                                      :: i,its,j,k,l,m,nn,n
real(kind=dp)                                :: anorm,p,q,r,s,t,u,v,w,x,y,z
logical                                      :: premature_exit

n = size(a,dim=2)

if (size(a,dim=1) /= n) then

    stop "A is not square"

end if

if (size(wi) < n) then

        stop "wi is not big enough in hqr"

end if

if (size(wr) < n) then

        stop "wr is not big enough in hqr"

end if

anorm=abs(a(1,1))

do i=2,n

    do j=i-1,n

        anorm=anorm+abs(a(i,j))

    end do

end do

nn=n
t=0.0_dp

do_1: do

    if(nn >= 1)then

        its=0

do_2:   do

            premature_exit = .false. 

do_l:       do l=nn,2,-1

                s=abs(a(l-1,l-1))+abs(a(l,l))

                if(s == 0.0_dp) then

                    s=anorm

                end if

                if(abs(a(l,l-1))+s == s)then

                    premature_exit = .true.
                    exit do_l

                end if

            end do do_l

            if (.not. premature_exit) then

                l=1

            end if

            x=a(nn,nn)

            if(l == nn)then

                wr(nn)=x+t
                wi(nn)=0.0_dp
                nn=nn-1

            else
   
                y=a(nn-1,nn-1)
                w=a(nn,nn-1)*a(nn-1,nn)

                if(l == nn-1)then
           
                    p=0.5_dp*(y-x)
                    q=p**2.0_dp+w
                    z=sqrt(abs(q))
                    x=x+t
            
                    if(q >= 0.0_dp)then
                
                        z=p+sign(z,p)
                        wr(nn)=x+z
                        wr(nn-1)=wr(nn)
 
                        if(z /= 0.0_dp)then
     
                            wr(nn)=x-w/z

                        end if
        
                        wi(nn)=0.0_dp
                        wi(nn-1)=0.0_dp
 
                    else
               
                        wr(nn)=x+p
                        wr(nn-1)=wr(nn)
                        wi(nn)=z
                        wi(nn-1)=-z
            
                    end if
            
                    nn=nn-2
          
                else
        
                    if(its == 30)then
  
                        stop "too many iterations in hqr"
         
                    end if

                    if(its == 10 .or. its == 20)then
      
                        t=t+x
              
                        do i=1,nn
                
                            a(i,i)=a(i,i)-x
              
                        end do

                        s=abs(a(nn,nn-1))+abs(a(nn-1,nn-2))
                        x=0.75_dp*s
                        y=x
                        w=-0.4375_dp*s**2.0_dp

                    end if
      
                    its=its+1
            
do_m:               do m=nn-2,l,-1
              
                        z=a(m,m)
                        r=x-z
                        s=y-z
                        p=(r*s-w)/a(m+1,m)+a(m,m+1)
                        q=a(m+1,m+1)-z-r-s
                        r=a(m+2,m+1)
                        s=abs(p)+abs(q)+abs(r)
                        p=p/s
                        q=q/s
                        r=r/s
             
                        if(m == l) then
                  
                            exit do_m

                        end if
          
                        u=abs(a(m,m-1))*(abs(q)+abs(r))
                        v=abs(p)*(abs(a(m-1,m-1))+abs(z)+abs(a(m+1,m+1)))
              
                        if(u+v == v)then
                  
                            exit do_m
             
                        end if

                    end do do_m

                    do i=m+2,nn
              
                        a(i,i-2)=0.0_dp
                    
                        if (i /= m+2) then
                  
                            a(i,i-3)=0.0_dp
              
                        end if
            
                    end do
            
                    do k=m,nn-1
              
                        if(k /= m)then
                
                            p=a(k,k-1)
                            q=a(k+1,k-1)
                            r=0.0_dp
              
                            if(k /= nn-1)then
                    
                                r=a(k+2,k-1)
                
                            end if
               
                            x=abs(p)+abs(q)+abs(r)
               
                            if(x /= 0.0_dp)then
                  
                                p=p/x
                                q=q/x
                                r=r/x
                
                            end if
              
                        end if

                        s=sign(sqrt(p**2.0_dp+q**2.0_dp+r**2.0_dp),p)
        
                        if(s /= 0.0_dp)then
                
                            if(k == m)then
                  
                                if(l /= m)then
                      
                                    a(k,k-1)=-a(k,k-1)
                  
                                end if
               
                            else
                  
                                a(k,k-1)=-s*x
                
                            end if
                
                            p=p+s
                            x=p/s
                            y=q/s
                            z=r/s
                            q=q/p
                            r=r/p
                
                            do j=k,nn
                  
                                p=a(k,j)+q*a(k+1,j)
                  
                                if(k /= nn-1)then
                    
                                    p=p+r*a(k+2,j)
                                    a(k+2,j)=a(k+2,j)-p*z
                  
                                end if
                
                                a(k+1,j)=a(k+1,j)-p*y
                                a(k,j)=a(k,j)-p*x
                
                            end do
                
                            do i=l,min(nn,k+3)
                   
                                p=x*a(i,k)+y*a(i,k+1)
                  
                                if(k /= nn-1)then
                 
                                    p=p+z*a(i,k+2)
                                    a(i,k+2)=a(i,k+2)-p*r
                  
                                end if
                  
                                a(i,k+1)=a(i,k+1)-p*q
                                a(i,k)=a(i,k)-p
                
                            end do
                        
                        end if

                    end do

                    cycle do_2

                end if

            end if

            exit do_2

            end do do_2
      
    else

            exit do_1
   
    end if

end do do_1

end subroutine hqr_dp
!  (C) Copr. 1986-92 Numerical Recipes Software Y].

end module m_diagonalizacion_dp






module mcf_slineales
 
use m_slineales_sp
use m_slineales_dp

public

end module mcf_slineales







module mcf_diagonalizacion

use m_diagonalizacion_sp
use m_diagonalizacion_dp

public

end module mcf_diagonalizacion






