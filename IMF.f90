!-------------------------------------------------------------------
module constants
!-------------------------------------------------------------------
! ***  Defines the constants needed in the IMF fucntions ***
!-------------------------------------------------------------------
	implicit none
	real*8 :: t
	real*8,parameter :: mc = 0.079d0
	real*8,parameter :: sigma = 0.69d0
	real*8,parameter :: x = 2.3d0
	real*8,parameter :: C1 = 2.3419119626187666d0
	real*8,parameter :: C2 = 0.65359841526902296d0
	real*8 :: C3
end module constants

!------------------------------------------------------------------
program initial_mass
!------------------------------------------------------------------
	implicit none
	
	call stars_alive
	
end program initial_mass

!-------------------------------------------------------------------
real*8 function IMF(m)
!-------------------------------------------------------------------
! *** Calulcates the IMF depending on the mass **8
!-------------------------------------------------------------------
	use constants, only : C1,C2,mc,sigma,x		!Calls constants
	implicit none
	real*8,intent (in) :: m				!Input mass
	real*8 :: y

	if (m<1.d0) then
		y = -((LOG10(m/mc))**2)/(2d0*sigma*sigma)
		IMF = C1 * EXP (y)/m			!IMF for m<1Msol
	else
		IMF = C2 * m**(-x)		!IMF for m>=1Msol
	end if
end function IMF

!-------------------------------------------------------------------
real*8 function tms(m)
!-------------------------------------------------------------------
! *** Caulcates the main sequence liftime of a star *** 
!-------------------------------------------------------------------
	use constants, only : C3		!Call C3
	implicit none
	real*8, intent(in) :: m
	
	C3 = 7.d9 * 10.d0**(-0.5d0)		!Exact values of C3
	
	if (m < 10.d0) then
		tms = 7.d9 * m**(-3)		!tms for m<10Msol
	else 
		tms = C3 * m**(-2.5d0)		!tms for m>10Msol
	end if

end function tms

!-------------------------------------------------------------------
real*8 function IMFt(m)
!-------------------------------------------------------------------
! *** Function to be called in calculating the number of  stars per Msol
! *** that are currently alive ***
!-------------------------------------------------------------------
	use constants, only : t	
	implicit none
	real*8, intent(in) :: m
	real*8, external :: IMF, tms
	
	t = tms(m)
	
	if (t > 13.d9) then
		IMFt = IMF(m) * 13.d09		!No. of stars/Msol alive if tms > 13Gyr
	else 
		IMFt = IMF(m) * tms(m)		!No. of stars/Msol alive if tms < 13Gyr
	end if

end function IMFt

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
end function INTEGRAL2

!-------------------------------------------------------------------
subroutine stars_alive
!-------------------------------------------------------------------
! *** Computes the current stars alive where ***
! *** dn/dm = 13Gyr if tms > 13Gy or
! ***       = tms if tms < 13Gyr
!-------------------------------------------------------------------
	implicit none
	real*8,external :: IMFt,INTEGRAL2
	real*8 :: total
	real*8 :: O,B,A,F,G,K,Mstar,BD
	
	
	total = INTEGRAL2(IMFt,0.01d0,150.0d0,10**5)		!Number of stars alive
	
	print '(A29,1pE10.4,A6)', "Total number of stars alive= ", total, " stars"
	
	BD = INTEGRAL2(IMFt,0.01d0,0.07d0,10**5)		!Number of brown dwarfs alive
	Mstar = INTEGRAL2(IMFt,0.07d0,0.45d0,10**5)		!Number of M stars alive
	K = INTEGRAL2(IMFt,0.45d0,0.8d0,10**5)			!Number of K stars alive
	G = INTEGRAL2(IMFt,0.8d0,1.04d0,10**5) 			!Number of G stars alive
	F = INTEGRAL2(IMFt,1.04d0,1.4d0,10**5) 			!Number of F stars alive
	A = INTEGRAL2(IMFt,1.4d0,2.1d0,10**5)			!Number of A stars alive
	B = INTEGRAL2(IMFt,2.1d0,16.d0,10**5)			!Number of B stars alive
	O = INTEGRAL2(IMFt,16.d0,150.d0,10**5)			!Number of O stars alive
	
	print '(A19,1pE10.4,A6)', "Number of O stars= ", O," stars"
	print '(A19,1pE10.4,A6)', "Number of B stars= ", B," stars"
	print '(A19,1pE10.4,A6)', "Number of A stars= ", A," stars"
	print '(A19,1pE10.4,A6)', "Number of F stars= ", F," stars"
	print '(A19,1pE10.4,A6)', "Number of G stars= ", G," stars"
	print '(A19,1pE10.4,A6)', "Number of K stars= ", K," stars"
	print '(A19,1pE10.4,A6)', "Number of M stars= ", Mstar," stars"
	print '(A24,1pE10.4,A13)', "Number of brown dwarfs= ", BD," brown dwarfs"

end subroutine stars_alive

