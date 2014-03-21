        SUBROUTINE aa (a)
	       IMPLICIT NONE
           DOUBLE PRECISION :: a
           CALL DBLEPR("a+1", -1, a+1.0, 1)
           CALL DBLEPR("a+1_int", -1, a+1, 1)
        END SUBROUTINE aa	   