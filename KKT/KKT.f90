      SUBROUTINE KKT(q, b0, beta, nobs, nvars, &
      & lam1, lam2, l, x, y, res, thr, vio_count, quiet)
         IMPLICIT NONE
!--------arg types---------------------------------------------
         INTEGER :: nobs, nvars, q, l, vio_count
         DOUBLE PRECISION :: b0(l), beta(nvars, l), thr
         DOUBLE PRECISION :: dl(nobs, l), y(nobs), lam1(l), lam2
         DOUBLE PRECISION :: x(nobs,nvars), res(nvars, l)    
         LOGICAL :: quiet 
		  
!--------local declarations-------------------------------------
         INTEGER :: i, j
		 DOUBLE PRECISION :: decib, fdr, r(nobs,l), dly(nobs, l)
		 DOUBLE PRECISION :: vio1, vio2
		 
!--------some initial setup-------------------------------------
         res = -5777.0D0
		 dly = -5777.0D0
	     vio_count = 0
!--------compute r: u_i-----------------------------------------		
! r
! In R: r = y * (x %*% beta + matrix(rep(b0,n), n, p, byrow=T))			 
         loop_r: DO i = 1, l
           r(:,i) = y * (Matmul(x, beta(:,i)) + b0(i))
         END DO loop_r		
		 
!--------compute dl: V'(u_i)------------------------------------
! dl
! In R: 
!        dlfun = function(r){
!                ifelse(r > decib, r ^ (-q-1) * fdr, -1)
!                }
!        dl = apply(r, c(1, 2), dlfun)
	
         decib = DBLE(q)/ (DBLE(q) + 1.0D0)    
         fdr =  decib ** (DBLE(q) + 1.0D0)                               
 
		 loop_dl_col: DO j = 1, l
            loop_dl_row: DO i = 1, nobs
               IF (r(i, j) > decib) THEN
                  dl (i, j) = -1.0 * r(i, j) ** (- q - 1) * fdr
                  ELSE
                     dl (i, j) = -1.0D0
               END IF
            END DO loop_dl_row
         END DO loop_dl_col
		 
!--------compute res-------------------------------------------
! For every column 		 
! res: \frac{1}{n} \sum_{i=1}^n V'(u_i)y_i x_{ij}
! res
! In R:  res = t(x) %*% (dl * y) / n
         loop_dly: DO i = 1, l
            dly(:, i) = y * dl(:, i)
         END DO loop_dly	
         res = Matmul(Transpose(x), dly) / nobs

!--------KKT check---------------------------------------------
        loop_KKT_col: DO j = 1, l
           loop_KKT_row: DO i = 1, nvars
		      IF (beta(i, j) == 0) THEN
			     vio1 = Abs(res(i, j)) - lam1(j)
			     IF (vio1 > thr) THEN
				    vio_count = vio_count + 1
					IF (quiet .EQV. .FALSE.) CALL DBLEPR("b=0", -1, vio1, 1)
				 END IF
			     ELSE  !beta_j \neq 0
				    vio2 = Abs(res(i, j) + Sign(lam1(j), beta(i, j)) + lam2 * beta(i, j))
                    IF (vio2 > thr) THEN
				      vio_count = vio_count + 1
					IF (quiet .EQV. .FALSE.) CALL DBLEPR("b!=0", -1, vio2, 1)
				    END IF
		      END IF
		   END DO loop_KKT_row
		END DO loop_KKT_col
		IF (quiet .EQV. .FALSE.) CALL INTPR("vio_count", -1, vio_count, 1)
        RETURN
      END SUBROUTINE KKT

	  
	  