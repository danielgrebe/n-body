PROGRAM main
    USE io, ONLY: read_config, read_data, write_u
    USE calculate, ONLY: initial_guess, construct_matrix, construct_b, solve, &
        rel_err
    IMPLICIT NONE
    ! Initialize
    INTEGER :: n_body, n, m
    DOUBLE PRECISION :: dt, t0, tf, error, eps, relax
    DOUBLE PRECISION, ALLOCATABLE :: masses(:), x0(:,:), xf(:,:), u_bar(:), &
        A(:,:), u_star(:), u(:), b(:)
    CALL read_config(n_body, dt, t0, tf, eps, relax)
    n = (tf-t0)/dt + 1
    m = 2*3*n_body*n
    ALLOCATE(masses(n_body), x0(3, n_body), xf(3, n_body), u_bar(m), A(m,m), &
        u_star(m), u(m), b(m))
    CALL read_data(n_body, masses, x0, xf)
    CALL initial_guess(n_body, n, m, x0, xf, u_bar)
    CALL construct_matrix(masses, dt, n_body, n, m, A)
    
    ! Solve
    error = 1
    DO WHILE (error > eps)
        CALL construct_b(u_bar, masses, n_body, m, n, b)
        CALL solve(A, b, m, u_star, error)
        u = relax*u_star + (1-relax)*u_bar
!       CALL rel_err(u_star, u_bar, m, error)
        CALL write_u(u, n, n_body, dt, t0,  m, error)
        u_bar = u
!       error = eps
    END DO

END PROGRAM
