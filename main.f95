PROGRAM main
    USE io, ONLY: read_config, read_data, write_u
    USE calculate, ONLY: initial_guess, construct_matrix, construct_b, solve, &
        rel_err, calc_b_prime
    IMPLICIT NONE
    ! Initialize
    INTEGER :: n_body, n, m
    DOUBLE PRECISION :: dt, t0, tf, error, eps, relax, h_u, gmma
    DOUBLE PRECISION, ALLOCATABLE :: masses(:), x0(:,:), xf(:,:), u_bar(:), &
        A(:,:), d_u(:), u(:), b(:), b_prime(:,:)
    CALL read_config(n_body, dt, t0, tf, eps, relax, h_u, gmma)
    n = (tf-t0)/dt + 1
    m = 2*3*n_body*n
    ALLOCATE(masses(n_body), x0(3, n_body), xf(3, n_body), u_bar(m), A(m,m), &
        d_u(m), u(m), b(m), b_prime(m,m))
    CALL read_data(n_body, masses, x0, xf)
    CALL initial_guess(n_body, n, m, x0, xf, u_bar)
    CALL construct_matrix(masses, dt, n_body, n, m, A, &
        RESHAPE(x0, (/3*n_body/)), RESHAPE(xf, (/3*n_body/)))
    
    ! Solve
    error = 1
    DO WHILE (error > eps)
        CALL construct_b(u_bar, masses, n_body, m, n, b)
        CALL calc_b_prime(u_bar, masses, n_body, m, n, b, h_u, b_prime)
        CALL solve(A - gmma*b_prime, MATMUL(A, u_bar) + b, m, d_u, error)
        u = relax*d_u + u_bar
        error = NORM2(d_u)
!       CALL rel_err(d_u, u_bar, m, error)
        CALL write_u(u, n, n_body, dt, t0,  m, error)
        u_bar = u
!       error = eps
    END DO

END PROGRAM
