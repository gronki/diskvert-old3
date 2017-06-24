module lapack

    implicit none

    interface GESV
        SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
            INTEGER            INFO, LDA, LDB, N, NRHS
            INTEGER            IPIV( * )
            DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
        END SUBROUTINE
    end interface


end module
