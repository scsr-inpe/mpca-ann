MODULE normalR8
    contains
    subroutine r8vec_normal_01(n, seed, x)
        !*****************************************************************************80
        !
        !! R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
        !
        !  Discussion:
        !
        !    An R8VEC is an array of double precision real values.
        !
        !    The standard normal probability distribution function (PDF) has
        !    mean 0 and standard deviation 1.
        !
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    18 May 2014
        !
        !  Author:
        !
        !    John Burkardt
        !
        !  Parameters:
        !
        !    Input, integer ( kind = 4 ) N, the number of values desired.
        !
        !    Input/output, integer ( kind = 4 ) SEED, a seed for the random
        !    number generator.
        !
        !    Output, real ( kind = 8 ) X(N), a sample of the standard normal PDF.
        !
        !  Local parameters:
        !
        !    Local, real ( kind = 8 ) R(N+1), is used to store some uniform
        !    random values.  Its dimension is N+1, but really it is only needed
        !    to be the smallest even number greater than or equal to N.
        !
        !    Local, integer ( kind = 4 ) X_LO_INDEX, X_HI_INDEX, records the range
        !    of entries of X that we need to compute
        !
        use uniformR8

        implicit none
        integer ( kind = 4) n

        integer ( kind = 4) m
        real ( kind = 8) r(n + 1)
        real ( kind = 8), parameter :: r8_pi = 3.141592653589793D+00
        integer ( kind = 4) seed
        real ( kind = 8) x(n)
        integer ( kind = 4) x_hi_index
        integer ( kind = 4) x_lo_index
        !
        !  Record the range of X we need to fill in.
        !
        x_lo_index = 1
        x_hi_index = n
        !
        !  If we need just one new value, do that here to avoid null arrays.
        !
        if (x_hi_index - x_lo_index + 1 == 1) then

            r(1) = r8_uniform_01(seed)

            if (r(1) == 0.0D+00) then
                write ( *, '(a)') ' '
                write ( *, '(a)') 'R8VEC_NORMAL_01 - Fatal error!'
                write ( *, '(a)') '  R8_UNIFORM_01 returned a value of 0.'
                stop 1
            end if

            r(2) = r8_uniform_01(seed)

            x(x_hi_index) = &
            sqrt(-2.0D+00 * log(r(1))) * cos(2.0D+00 * r8_pi * r(2))
            !
            !  If we require an even number of values, that's easy.
            !
        elseif (mod(x_hi_index - x_lo_index, 2) == 1) then

            m = (x_hi_index - x_lo_index + 1) / 2

            call r8vec_uniform_01(2 * m, seed, r)

            x(x_lo_index:x_hi_index - 1:2) = &
            sqrt(-2.0D+00 * log(r(1:2 * m - 1:2))) &
            *cos(2.0D+00 * r8_pi * r(2:2 * m:2))

            x(x_lo_index + 1:x_hi_index:2) = &
            sqrt(-2.0D+00 * log(r(1:2 * m - 1:2))) &
            *sin(2.0D+00 * r8_pi * r(2:2 * m:2))
            !
            !  If we require an odd number of values, we generate an even number,
            !  and handle the last pair specially, storing one in X(N), and
            !  saving the other for later.
            !
        else

            x_hi_index = x_hi_index - 1

            m = (x_hi_index - x_lo_index + 1) / 2 + 1

            call r8vec_uniform_01(2 * m, seed, r)

            x(x_lo_index:x_hi_index - 1:2) = &
            sqrt(-2.0D+00 * log(r(1:2 * m - 3:2))) &
            *cos(2.0D+00 * r8_pi * r(2:2 * m - 2:2))

            x(x_lo_index + 1:x_hi_index:2) = &
            sqrt(-2.0D+00 * log(r(1:2 * m - 3:2))) &
            *sin(2.0D+00 * r8_pi * r(2:2 * m - 2:2))

            x(n) = sqrt(-2.0D+00 * log(r(2 * m - 1))) &
            *cos(2.0D+00 * r8_pi * r(2 * m))

        end if

        return
  !  end
END SUBROUTINE
 end MODULE normalR8
