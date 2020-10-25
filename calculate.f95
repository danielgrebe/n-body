MODULE calculate
    USE constants, ONLY: nl_name, G
    IMPLICIT NONE
    PRIVATE
    PUBLIC initial_guess, construct_matrix, construct_b, solve, rel_err, &
        calc_b_prime
    
CONTAINS

    SUBROUTINE initial_guess(n_body, n, m, x0, xf, u)
        INTEGER :: n_body, n, m, i
        DOUBLE PRECISION :: x0(3,n_body), xf(3, n_body), u(m), dx(3,n_body) 
        dx = (xf - x0)/n
        u = 0
        DO i=0,n-1
            u(i*3*n_body+1:(i+1)*3*n_body) = RESHAPE(x0+i*dx, (/3*n_body/)) 
        END DO
!       DO i = n*n_body*3+1,m
!           u(i) = 1
!       END DO
    END SUBROUTINE

    SUBROUTINE construct_matrix(masses, dt, n_body, n, m, A, x0, xf)
        INTEGER :: n_body, n, m, i, offset, body_idx
        DOUBLE PRECISION :: masses(n_body), dt, A(m,m), mass, x0(3*n_body), &
            xf(3*n_body)
        A = 0
        offset = 3*n_body*n
        !initialize boundary conditions
        DO i=1,n_body*3
            A(i,i) = 1/x0(i)
            A(i+offset,i) = 1/dt
            A(i+offset,i+3*n_body) = -1/dt
        END DO
        DO i=n_body*(n-1)*3+1,n*n_body*3
            A(i,i) = 1/xf(i-(n-1)*3*n_body)
            A(i+offset,i) = -1/dt
            A(i+offset,i-3*n_body) = 1/dt
        END DO

        DO i = n_body*3+1, offset - n_body*3
        ! initialize inertia
            body_idx = MOD(i, n_body)/3
            IF (body_idx .EQ. 0) THEN
                body_idx = n_body
            END IF
            mass = masses(body_idx)
            A(i, offset+i-n_body*3) = -mass/(2*dt)
            A(i, offset+i+n_body*3) = mass/(2*dt)
        ! initialize  derivative
            A(offset+i, i-n_body*3) = 1/(2*dt)
            A(offset+i, i+n_body*3) = -1/(2*dt)
        END DO

        DO i = offset+1, m
            A(i,i) = 1
        END DO
!       DO i = 1, m
!           PRINT *, A(i,:)
!       END DO


    END SUBROUTINE


    SUBROUTINE construct_b(u_bar, masses, n_body, m, n, b)
        INTEGER :: n_body, m, n, body_idx, t_idx, i, k, coord_body_idx
        DOUBLE PRECISION :: u_bar(m), masses(n_body), b(m), r(3), norm_r
        DOUBLE PRECISION, EXTERNAL :: DNRM2
        b = 0
        ! boundary condition
!       b(1:3*n_body) = u_bar(1:3*n_body)
        b(1:3*n_body) = 1
!       b(3*n_body*(n-1)+1:3*n_body*n) = u_bar(3*n_body*(n-1)+1:3*n_body*n)
        b(3*n_body*(n-1)+1:3*n_body*n) = 1

        ! gravitational attraction
        DO i = n_body*3+1, 3*n_body*(n-1), 3
            body_idx = MOD(i, n_body)/3 + 1
            t_idx = i/(n_body*3)
            IF (body_idx .EQ. 0) THEN
                body_idx = n_body
            END IF
            coord_body_idx = MOD(i,3)-1
            IF (coord_body_idx .EQ. -1) THEN
                coord_body_idx = 2
            END IF
            DO k = 1,n_body
                IF (body_idx .NE. k) THEN
                    r = u_bar(3*n_body*t_idx+k*3:3*n_body*t_idx+k*3+2) &
                        -u_bar(i:i+2)
                    norm_r = DNRM2(3, r, 1)
                    b(i:i+2) = b(i:i+2) + r/norm_r*G*masses(body_idx)/norm_r* &
                        masses(k)/norm_r
                END IF
            END DO
        END DO
    END SUBROUTINE


    SUBROUTINE calc_b_prime(u_bar, masses, n_body, m, n, b, h_u, b_prime)
        INTEGER :: n_body, m, n, body_idx, t_idx, i, k, coord_body_idx
        DOUBLE PRECISION :: u_bar(m), masses(n_body), b(m), &
            h_u, b_prime(m,m), b_star(m), u_prime(m)
        DOUBLE PRECISION, EXTERNAL :: DNRM2
        !$OMP PARALLEL
        !$OMP DO
        DO i = 1,m
            u_prime = u_bar
            u_prime(i) = u_prime(i) + h_u
            CALL construct_b(u_prime, masses, n_body, m, n, b_star)
            b_prime(:,i) = (b_star - b)/h_u
        END DO
        !$OMP END DO
        !$OMP END PARALLEL

    END SUBROUTINE

    SUBROUTINE solve(A, b, m, u_star, FERR)
        INTEGER :: m, i
        INTEGER :: pivot(m), IWORK(m), INFO, RANK
        DOUBLE PRECISION :: A(m,m), AF(m,m), b(m), u_star(m), R(m), C(m), &
            RCOND, FERR, BERR(1), WORK(4*m), S(m)
        CHARACTER :: EQUED
        CALL DGESVX('E', 'N', m, 1, A, m, AF, m, pivot, EQUED, R, C, b, m, &
            u_star, m, RCOND, FERR, BERR, WORK, IWORK, INFO)
!       CALL DGELSD(m, m, 1, A, m, b, m, S, RCOND, RANK, WORK, m*m, IWORK, INFO)
        PRINT *, INFO
!       PRINT *, SUM(A(INFO,:))
!       PRINT *, SUM(AF(220,:))
!       DO i = 1, m
!           IF (MAXVAL(A(i,:)) .EQ. 0) THEN
!               PRINT *, i
!           END IF
!       END DO
        PRINT *, FERR
        PRINT *, RCOND
    END SUBROUTINE

    SUBROUTINE rel_err(u_star, u_bar, m, error)
        INTEGER :: m
        DOUBLE PRECISION :: u_star(m), u_bar(m), error
        error = NORM2(u_star-u_bar)/NORM2(u_bar) - 1
    END SUBROUTINE

END MODULE
