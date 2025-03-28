
module Rand
        use, intrinsic :: iso_fortran_env, only: int32, dp => real64 
        type random 
                integer :: seed, ix, iy  
        end type random 
        public :: random, makeRandom, algor_uniform_random, Normal      
       !ix=ieor(777755555_int32,iseed)                   !Marsaglia generator        
        contains
        type(random) function makeRandom(seed)
                integer(kind=int32), intent(in) :: seed
                type(random) :: output
                output%seed = seed 
                output%ix = ieor(777755555_int32,seed)
                output%iy = ior(ieor(888889999_int32,seed),1_int32)     !Parks-Miller generator
                makeRandom = output
        end function makeRandom 
 function algor_uniform_random(inRand)
    !=========================================================================!
    ! Return a single random deviate ~ uniform [0,1].                         !
    ! Based on Park-Miller "minimal standard" generator with Schrage's method !
    !  to do 32-bit multiplication without requiring higher precision, plus   !
    !  additional Marsaglia shift to suppress any weaknesses & correlations.  !
    ! Using two independent methods greatly increases the period of the       !
    !  generator, s.t. resulting period ~2*10^18                              !
    ! NB Routine is only set to work with 32 bit integers!                    !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   none                                                                  !
    !-------------------------------------------------------------------------!
    ! References:                                                             !
    !   S.K. Park and K.W. Miller, Commun. ACM, 31, p1192-1201 (1988)         !
    !   L. Schrage, ACM Trans. Math. Soft., 5, p132-138 (1979)                !
    !   G. Marsaglia, Linear Algebra and its Applications, 67, p147-156 (1985)!
    !-------------------------------------------------------------------------!
    ! Return value:                                                           !
    !   algor_uniform_random => required random deviate                       !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   ix as next 32-bit integer in Marsaglia generator (updated)            !
    !   iy as next 32-bit integer in Park-Miller generator (updated)          !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   none                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   none                                                                  !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   none                                                                  !
    !-------------------------------------------------------------------------!
    ! Written by Matt Probert, v0.1, 01/07/2000                               !
    !=========================================================================!

    implicit none
    real(kind=dp)                 :: algor_uniform_random

    !NB We use 3 logical XOR operations to do Marsaglia shift
    !=> hard-wire 3 shift values (not all triplets any good)
    !=> entire routine preset to only work with 32 bit integers.

    !local variables ...
    type(random), intent(inout) :: inRand 
    integer(kind=int32)            :: iy_tmp       !working value to force integer division
    integer(kind=int32), parameter :: iy_max=2147483647_int32 !2^31-1
    real(kind=dp), parameter        :: inv_iy_max=1.0_dp/2147483647.0_dp
    integer                         :: seed

    !Catch uninitialised random number generator and set to random seed
    integer                         :: status
    integer                         :: iy,ix ! SAM: putting this here as a stop gap to manage state


    seed=inRand%seed 
    iy = inRand%iy 
    !do Marsaglia shift sequence, period 2^32-1, to get new ix
    ix = inRand%ix 
    ix=ieor(ix,ishft(ix, 13_int32))
    ix=ieor(ix,ishft(ix,-17_int32))
    ix=ieor(ix,ishft(ix,  5_int32))


    !Park-Miller sequence, period iy_max-1, to get new iy
    iy_tmp=iy/127773_int32                         !NB integer division
    iy=16807_int32*(iy-iy_tmp*127773_int32)-iy_tmp*2836_int32  !next value of iy
    if (iy < 0_int32) iy=iy+iy_max                 !integer modulo iy_max

    !Combine ix and iy to get new random number, rescale onto range [0,1]
    algor_uniform_random=inv_iy_max*ior(iand(iy_max,ieor(ix,iy)),1_int32) !with masking to ensure non-zero
    inRand%ix = ix 
    inRand%iy = iy 

    return
  end function algor_uniform_random

  function Normal(rand,mean,stdev) result(z1)
        type(random), intent(inout) :: rand 
        real(kind=8), intent(in) :: mean, stdev
        real(kind=8) :: u1,u2, z1, z2
        real(kind=8), parameter :: pi = 3.14159265358979323846_8

        u1 = algor_uniform_random(rand)
        u2 = algor_uniform_random(rand)

        z1 = sqrt(-2*log(u1)) * cos(2*pi*u2)
        z2 = sqrt(-2*log(u1)) * sin(2*pi*u2)
         
        z1 = mean + stdev*z1 
        z2 = mean + stdev*z2
  end function Normal
end module Rand 

