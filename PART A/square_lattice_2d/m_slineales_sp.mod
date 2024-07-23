G95 module created on Fri Nov 11 16:12:13 2022 from ../modulos/mcf_matrices.f95
If you edit this, you'll get what you deserve.
module-version 9
(() () () () () () () () () () () () () () () () () () () () ())

()

(('ctridag' 2) ('gaussj' 3) ('lu_descomposicion' 4) ('lu_resolucion' 5)
('lubksb' 5) ('ludcmp' 4) ('print_matrix' 6) ('tridag' 7))

()

()

(8 'ctridag' '(global)' 1 ((PROCEDURE UNKNOWN UNKNOWN UNKNOWN NONE NONE)
(UNKNOWN) 0 0 () () () '' () ())
9 'gaussj' '(global)' 1 ((PROCEDURE UNKNOWN UNKNOWN UNKNOWN NONE NONE) (
UNKNOWN) 0 0 () () () '' () ())
10 'lu_descomposicion' '(global)' 1 ((PROCEDURE UNKNOWN UNKNOWN UNKNOWN
NONE NONE) (UNKNOWN) 0 0 () () () '' () ())
11 'lu_resolucion' '(global)' 1 ((PROCEDURE UNKNOWN UNKNOWN UNKNOWN NONE
NONE) (UNKNOWN) 0 0 () () () '' () ())
12 'lubksb' '(global)' 1 ((PROCEDURE UNKNOWN UNKNOWN UNKNOWN NONE NONE)
(UNKNOWN) 0 0 () () () '' () ())
13 'ludcmp' '(global)' 1 ((PROCEDURE UNKNOWN UNKNOWN UNKNOWN NONE NONE)
(UNKNOWN) 0 0 () () () '' () ())
14 'print_matrix' '(global)' 1 ((PROCEDURE UNKNOWN UNKNOWN UNKNOWN NONE
NONE) (UNKNOWN) 0 0 () () () '' () ())
15 'tridag' '(global)' 1 ((PROCEDURE UNKNOWN UNKNOWN UNKNOWN NONE NONE)
(UNKNOWN) 0 0 () () () '' () ())
7 'tridag_sp' 'm_slineales_sp' 1 ((PROCEDURE UNKNOWN MODULE-PROC DECL
NONE NONE SUBROUTINE INVOKED) (PROCEDURE 0) 0 0 (16 NONE 17 NONE 18 NONE
19 NONE 20 NONE 21 NONE) () () '' () ())
21 'n' '' 22 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (INTEGER 8)
0 0 () () () '' () ())
20 'x' '' 22 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
REAL 4) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
19 'b' '' 22 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
REAL 4) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
18 'c' '' 22 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
REAL 4) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
17 'd' '' 22 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
REAL 4) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
16 'a' '' 22 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
REAL 4) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
6 'print_matrix_sp' 'm_slineales_sp' 1 ((PROCEDURE UNKNOWN MODULE-PROC
DECL NONE NONE SUBROUTINE INVOKED) (PROCEDURE 0) 0 0 (23 NONE 24 NONE) ()
() '' () ())
24 'unit' '' 25 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (INTEGER
8) 0 0 () () () '' () ())
23 'a' '' 25 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
REAL 4) 0 0 () (2 ASSUMED_SHAPE () () () ()) () '' () ())
5 'lu_resolucion_sp' 'm_slineales_sp' 1 ((PROCEDURE UNKNOWN MODULE-PROC
DECL NONE NONE SUBROUTINE INVOKED) (PROCEDURE 0) 0 0 (26 NONE 27 NONE 28
NONE) () () '' () ())
28 'b' '' 29 ((VARIABLE INOUT UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY)
(REAL 4) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
27 'indx' '' 29 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY)
(INTEGER 8) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
26 'a' '' 29 ((VARIABLE INOUT UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY)
(REAL 4) 0 0 () (2 ASSUMED_SHAPE () () () ()) () '' () ())
4 'lu_descomposicion_sp' 'm_slineales_sp' 1 ((PROCEDURE UNKNOWN
MODULE-PROC DECL NONE NONE SUBROUTINE INVOKED) (PROCEDURE 0) 0 0 (30
NONE 31 NONE 32 NONE) () () '' () ())
32 'd' '' 33 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DUMMY) (REAL 4) 0
0 () () () '' () ())
31 'indx' '' 33 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY)
(INTEGER 8) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
30 'a' '' 33 ((VARIABLE INOUT UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY)
(REAL 4) 0 0 () (2 ASSUMED_SHAPE () () () ()) () '' () ())
3 'gaussj_sp' 'm_slineales_sp' 1 ((PROCEDURE UNKNOWN MODULE-PROC DECL
NONE NONE SUBROUTINE INVOKED) (PROCEDURE 0) 0 0 (34 NONE 35 NONE) () () ''
() ())
35 'b' '' 36 ((VARIABLE INOUT UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY)
(REAL 4) 0 0 () (2 ASSUMED_SHAPE () () () ()) () '' () ())
34 'a' '' 36 ((VARIABLE INOUT UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY)
(REAL 4) 0 0 () (2 ASSUMED_SHAPE () () () ()) () '' () ())
2 'ctridag_sp' 'm_slineales_sp' 1 ((PROCEDURE UNKNOWN MODULE-PROC DECL
NONE NONE SUBROUTINE INVOKED) (PROCEDURE 0) 0 0 (37 NONE 38 NONE 39 NONE
40 NONE 41 NONE 42 NONE) () () '' () ())
42 'n' '' 43 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DUMMY) (INTEGER 8)
0 0 () () () '' () ())
41 'x' '' 43 ((VARIABLE OUT UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
COMPLEX 4) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
40 'b' '' 43 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
COMPLEX 4) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
39 'c' '' 43 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
COMPLEX 4) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
38 'd' '' 43 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
COMPLEX 4) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
37 'a' '' 43 ((VARIABLE IN UNKNOWN UNKNOWN NONE NONE DIMENSION DUMMY) (
COMPLEX 4) 0 0 () (1 ASSUMED_SHAPE () ()) () '' () ())
)

('ctridag' 0 8 'gaussj' 0 9 'lu_descomposicion' 0 10 'lu_resolucion' 0
11 'lubksb' 0 12 'ludcmp' 0 13 'print_matrix' 0 14 'tridag' 0 15)
