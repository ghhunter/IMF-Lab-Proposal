program IMF
	implicit none
	real*8, external :: INTEGRAL2, probdist			!Calls funcitons
	real*8 :: integral, error

	integral = INTEGRAL2(probdist,0d0,10d0,10000)	!Integrates to get probability
	
	error = abs(integral-1d0)						!Calculats the error in the value
	
	print *, "The value of the integral is ", integral
	print *, "This gives an error of ", error, " from the value of 1."


end program IMF

!----------------------------------------------------------------------------
real*8 function probdist(x)
!----------------------------------------------------------------------------
! *** Wave function of a particle in an infinite deep square well ***
!----------------------------------------------------------------------------	
	implicit none
	real*8, intent(in) :: x			!Input values of x
	real*8 :: L,n,pi,root,y,wave	!Constants
	
	L = 10d0						!Length of the 
	n = 1d0							!Energy state where n=1 for the ground state
	pi = ACOS(-1d0)					!Just pi
	
	root = SQRT(2d0/L) 				
	y = n*pi*x/L
	
	wave = root * SIN(y)			!Calculates the wave function

	probdist = wave**2				!Calculates the probability distribution
end function probdist


!----------------------------------------------------------------------------
real*8 function INTEGRAL2(FUNC,a,b,N)
!----------------------------------------------------------------------------
! ***  numerical integration (Simpson-rule) with equidistant spacing      ***
!----------------------------------------------------------------------------
	implicit none
	real*8,external :: FUNC                    ! the function to be integrated
	real*8,intent(in) :: a,b                   ! boundary values
	integer,intent(in) :: N                    ! number of sub-intervals
	real*8 :: dx,x1,x2,xm,f1,f2,fm,int         ! local variables
	integer :: i
	dx  = (b-a)/DBLE(N)                        ! x subinterval
	x1  = a                                    ! left   
	f1  = FUNC(a)
	int = 0.d0
	do i=1,N
		x2  = a+DBLE(i)*dx                       ! right
		xm  = 0.5d0*(x1+x2)                      ! midpoint
		f2  = FUNC(x2)
		fm  = FUNC(xm)
		int = int + (f1+4.d0*fm+f2)/6.d0*(x2-x1) ! Simpson rule
		x1  = x2
		f1  = f2                                 ! save for next subinterval
	enddo
	INTEGRAL2 = int
end 
