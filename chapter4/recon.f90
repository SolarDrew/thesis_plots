SUBROUTINE reconstruct(t, resp, logt, x, y, emiss)

INTEGER, INTENT(IN) :: x, y
REAL, DIMENSION(x, y), INTENT(IN) :: t
REAL, DIMENSION(301), INTENT(IN) :: resp, logt
REAL, DIMENSION(x, y), INTENT(OUT) :: emiss
REAL :: mean, width, height
REAL, DIMENSION(301) :: dem, f

  PRINT*, MINVAL(resp), MAXVAL(resp)
  PRINT*, MINVAL(t), MAXVAL(t)
  PRINT*, logt

  DO j = 1, y
    DO i = 1, x
      mean = t(i, j)
      height = 1.0
      width = 0.1

      dem = height * EXP(-((logt - mean) ** 2.0) / (2.0 * (width ** 2.0)))
      f = resp * dem
      emiss(i, j) = SUM(f)
    END DO
  END DO

END SUBROUTINE
