MODULE io
    USE constants, ONLY: nl_name, u_name, x_name
    IMPLICIT NONE
    PRIVATE
    PUBLIC read_config, read_data, write_u
    REAL, PARAMETER :: G = 6.6743e-10
    
CONTAINS

    SUBROUTINE read_config(n_body, dt, t0, tf, eps, relax, d_u, gmma)
        INTEGER :: n_body
        DOUBLE PRECISION :: dt, t0, tf, eps, relax, d_u, gmma
        NAMELIST /config/ n_body, dt, t0, tf, eps, relax, d_u, gmma
        OPEN(10, FILE=nl_name)
        READ(10, NML=config)
        WRITE(*, NML=config)
        CLOSE(10)
    END SUBROUTINE

    SUBROUTINE read_data(n_body, masses, x0, xf)
        INTEGER :: n_body, n
        DOUBLE PRECISION :: masses(n_body), x0(3, n_body), xf(3, n_body)
        NAMELIST /init/ masses, x0, xf
        OPEN(10, FILE=nl_name)
        READ(10, NML=init)
        WRITE(*, NML=init)
        CLOSE(10)
    END SUBROUTINE

    SUBROUTINE write_u(u, n, n_body, dt, t0, m, error)
        INTEGER :: m, n, n_body, i
        DOUBLE PRECISION :: u(m), error, dt, t0
        PRINT *, 'The error is ', error
        OPEN(10, FILE=u_name)
        WRITE(10,*) u
        close(10)
        OPEN(10, FILE=x_name)
        WRITE(10,"(2A)", ADVANCE="no") 't '
        DO i=1,n_body
            WRITE(10,"(2A)", ADVANCE='no') " X"
            WRITE(10,"(I2.2)", ADVANCE='no') i
            WRITE(10,"(2A)", ADVANCE='no') " Y"
            WRITE(10,"(I2.2)", ADVANCE='no') i
            WRITE(10,"(2A)", ADVANCE='no') " Z"
            WRITE(10,"(I2.2)", ADVANCE='no') i
        END DO
        WRITE(10, *) ' '
        DO i=1,n
            WRITE(10,"(F15.2)", ADVANCE='no') dt*(i-1)
            WRITE(10,*) u((i-1)*3*n_body+1:i*3*n_body)
        END DO
    END SUBROUTINE


END MODULE
