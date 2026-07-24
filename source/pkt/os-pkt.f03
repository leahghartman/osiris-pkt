!-----------------------------------------------------------------------------------------
! pkt calculation functions
!-----------------------------------------------------------------------------------------


module m_pkt

  use m_system
  use m_parameters
  use m_math
  
  private
  
  real( p_double ), parameter :: M_PI    = 3.14159265358979323846264338327950288d0


  interface syncpkt
    module procedure syncpkt
  end interface
  
  interface find_photon_chi
    module procedure find_photon_chi
  end interface

  interface specpair_fast
    module procedure specpair_fast
  end interface

  interface interp_spec_pkt
    module procedure interp_spec_pkt
  end interface

  interface interp_spec_pairpkt
    module procedure interp_spec_pairpkt
  end interface

  
  public :: syncpkt, find_photon_chi, specpair_fast, interp_spec_pkt, interp_spec_pairpkt
  
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  
contains

!-------------------------------------------------------------------------------
! generates a gaussian random number 
!-------------------------------------------------------------------------------
	function  gauss()
	
    implicit none
    
	real(p_double) :: gauss
	real(p_double) :: x, P, varand(2), vmin, vmax, v
	
	vmin = -6.0d0
	vmax =  6.0d0	
	x = 1.0d0
	P = 0.0d0
		
	do while (x > P)	
    	call random_number( varand )
    	v = vmin + (vmax - vmin)*varand(1)
    	x = varand(2)
    	P = exp(-v**2/2.0d0)    
    enddo
    
    gauss = v

	end function gauss
	
	
!-------------------------------------------------------------------------------
! generates a random number distributed according to the pkt spectrum
!-------------------------------------------------------------------------------
	function  syncpkt(eta)
	
    implicit none
    
	real(p_double) :: syncpkt
    real(p_double) :: x, P, varand(2), chimax, chi0, chi, eta, logchimin, logchimax, logchi
	
	chi0 = 0.214d0 * eta**2 / (1 + eta**4) + (0.5d0*eta - 1.d-6) / (1 / (eta**3) + 1)	
	chimax = 0.5d0*eta - 1.d-6	
	logchimin = log10(chi0)-12.d0
	logchimax = log10(chimax)	
	x = 1.0d0
	P = 0.0d0
		
	do while (x > P)
    	call random_number( varand )
    	logchi = logchimin + (logchimax - logchimin)*varand(1)

   	 	chi = exp(log(10.d0)*logchi)  
   	 	
   	 	! Test if this is faster (probably implemented the same way)
   	 	! chi = 10.0d0 ** logchi  

    	x = varand(2)*peak_pkt(eta)  
    	P = F_pkt(eta,chi)            
    enddo
    
    syncpkt = chi

	end function syncpkt
	
	
!-------------------------------------------------------------------------------
! generates a random number distributed according to the pair emission spectrum
!-------------------------------------------------------------------------------
	function  specpair(chi_g)
	
    implicit none
    
	real(p_double) :: specpair
	real(p_double) :: r, P, varand(2), etamin, etamax, chi_g, eta	
	
	etamin = 1.0d-12
	etamax = chi_g-1e-9	
	r = 1.0d0
	P = 0.0d0
	
	do while (r > P)
    	call random_number( varand )
    	eta = etamin + (etamax - etamin)*varand(1)  
    	r = varand(2)*peak_pairpkt(chi_g) 
    	P = F_pair(eta,chi_g) 	
    enddo
        
    specpair = eta

	end function specpair	
	
!-----------------------------------------------------------------------------------------
! generates a random eta distributed according to the pair emission spectrum
! This version calculates the samee as specpair, but using already interpolated values 
! of the Bessel functions -> one random is enough to emit a pair for 0.1<chi_g<100.0
!-----------------------------------------------------------------------------------------
	function  specpair_fast(chi_g)
	
	real(p_double) :: specpair_fast
	real(p_double) :: ratio, ratiomin, ratiomax, chi_g, eta, lt_c, ltr	
	real(p_double) :: pc1, pc2, pc3, pc4, pc5, pc
	integer :: n
	logical :: mirror

	! coefficients for 0.1 < chi_g < 10     
	real(p_double), dimension(25), parameter :: p1_array = (/ -0.000005225401859,     -0.000080299602413,     0.000140584960763,     0.009000197490255,     0.035872232367153, &
	     -0.000512836744468,     -0.005836167050800,     -0.024756793004610,     -0.039904597410351,     -0.001077077386309, &
	     -0.004136487977769,     -0.030318650457853,     -0.087607137642237,     -0.112531074265945,     -0.032921973247876,  &
	      0.057693258982845,     0.179283783937442,     0.180705696041231,     0.041353652234700,     0.000475472998592,  &
	     -0.145328538070807,     -0.233082359479153,     -0.124513341694971,     -0.055982996943532,     -0.010737979156167 /)
	real(p_double), dimension(25), parameter :: p2_array = (/ -0.000017186457941,     -0.000501308272203,     -0.005835292277119,     -0.031158652183816,     -0.014223734979921,  &
          0.000017025291305,     -0.000183506798540,     -0.004884676164642,     -0.030445605558507,     -0.014858163874475,   &
          0.002544961647828,     0.013669710411411,     0.023407691578925,     -0.004987492811994,     -0.006377414300204,  &
          0.020533649186576,     0.088984339299603,     0.140104063967745,     0.074713915215919,     0.013935964256991,  &
         -0.793006417498808,     -1.555472861127672,     -1.115360068815503,     -0.354708337189465,     -0.041628012743737 /)
	real(p_double), dimension(25), parameter :: p3_array = (/ 0.000053676346116,     0.001342395538715,     0.011878135658543,     0.033049144144598,     -0.129859763407402,  &
	      0.001864793520613,     0.022926907676793,     0.111050470233787,     0.241997090137700,     0.040713187936894,  &
	      0.010915660757696,     0.081889940287695,     0.257794761863463,     0.407301837742367,     0.111748829824497,  &
	     -0.094704410042950,     -0.259225069042258,     -0.155491759584177,     0.184605599105285,     0.066710295541102,  &
	     -1.321702458506052,     -2.669062636880125,     -1.968481008170236,     -0.435526991491079,     -0.014677433853649 /)
	real(p_double), dimension(25), parameter :: p4_array = (/ 0.000146287990639,     0.004170117816959,     0.047793530425434,     0.281487950503571,     0.891407062392472,  &
	      0.002355625428871,     0.031357451125372,     0.177383183278795,     0.565709748126881,     1.133279049852064,   &
	      0.008894669086513,     0.075669333596765,     0.291976904248376,     0.699564295108904,     1.192750654018492,   &
	     -0.052884210441859,     -0.128784161738488,     0.037015124061894,     0.557408165404765,     1.162815353169648,  &
	     -0.626653322990238,     -1.231597890005989,     -0.776953017271845,     0.283251642562608,     1.127225781030429 /)
	real(p_double), dimension(25), parameter :: p5_array = (/ 0.000101035895204,     0.003034961791079,     0.038517000340175,     0.296395206350618,     -0.340255742643951,  &
	      0.001367534920433,     0.018706596307303,     0.113713739198754,     0.462593023080581,     -0.197612683567629,  &
	      0.007675199160302,     0.058921747887693,     0.211874594432208,     0.571362244498481,     -0.151477476660753,  &
	      0.036908328292682,     0.162109824654076,     0.350672986287757,     0.655784776438876,     -0.131875713605264,  &
	      0.421923148153023,     0.869248252496938,     0.845970255776770,     0.813262278773464,     -0.112629271074978 /)
	! coefficients for 10 < chi_g < 100      
	real(p_double), dimension(25), parameter :: p6_array = (/ 0.000006082040645,     0.000564737585995,     0.008242128861585,     0.046025770183968,     0.084647285368009,  &
	     -0.002485658336182,     -0.021312822865797,     -0.052708137449931,     0.005176422736026,     0.122282126102417,  &
	     -0.029376830591882,     -0.213113658843963,     -0.568537193360927,     -0.614609623468511,     -0.158288550851653,  &
	      0.111656465232389,     0.257195253365006,     0.026946304622848,     -0.274468967945381,     -0.084172472851054,   &
	     -0.810417539310961,     -0.880137768795069,     -0.281075095233155,     -0.188412694151709,     -0.048622693553012 /)
	real(p_double), dimension(25), parameter :: p7_array = (/ 0.000190332705603,     0.004331676342178,     0.041328744292588,     0.203626212858701,     0.533900410274208,  &
	      0.010500649205440,     0.111418057985674,     0.453780672610153,     0.897182903247248,     0.958322396366310,  &
	      0.067311632302158,     0.525103861335212,     1.586743785301362,     2.280266329245493,     1.593141572200228,  &
	     -0.109511244628452,     -0.038419575312760,     0.903676736272272,     1.904879603203579,     1.513766659274223,  &
	     -5.524238464421281,     -11.990343682335268,     -8.993570148855301,     -1.742244248896068,     1.008817566301562 /)
	real(p_double), dimension(25), parameter :: p8_array = (/ 0.000081962355663,     0.003080647795245,     0.043137842803378,     0.340288131821678,     -0.174328302877703,  &
	     -0.003003229443531,     -0.023844138917779,     -0.030457433806767,     0.297061696048246,     -0.117978287081051,  &
	     -0.012253448661188,     -0.114401777093548,     -0.328095829070981,     -0.113331976962520,     -0.322994083205383,  &
	     -0.037353987824001,     -0.187562068001494,     -0.396651342055515,     -0.131745222357037,     -0.320882921260851,   &
	      3.978736012514611,     8.264364245327915,     6.286822482443694,     2.223401647058546,     -0.008653167964833 /)

	! Boolean to mirror the value obtained
	mirror = .false.
	
	ratiomin = 1.e-8
	ratiomax = 1.0-1.e-8
	if ( ( chi_g < 0.1 ) .or. ( chi_g > 100.0 ) ) then
		eta = specpair(chi_g)
	else
		call random_number( ratio )
		lt_c = log10(chi_g)

		! Check if we use the symmetry
        if (ratio > 0.5) then
            mirror = .true.
            ratio = 1.0 - ratio
		endif

		if ( ratio < ratiomin ) then
			ratio = ratiomin
		else if ( ratio > ratiomax ) then
			ratio = ratiomax
		endif
				
	    if ( ratio < 0.007 ) then 
    		n=1
    	else if ( ratio < 0.025 ) then
    		n = 6
    	else if ( ratio < 0.12 ) then
    		n = 11
    	else if (ratio < 0.3) then  
    		n = 16
    	else
			n = 21

    	endif
		
		if 	(lt_c < 1.0) then 
			pc1 = p1_array(n) * ( lt_c**4.0 ) + p2_array(n) * ( lt_c**3.0 ) + p3_array(n) * ( lt_c**2.0 ) + p4_array(n) * lt_c + p5_array(n) 
    		pc2 = p1_array(n+1) * ( lt_c**4.0 ) + p2_array(n+1) * ( lt_c**3.0 ) + p3_array(n+1) * ( lt_c**2.0 ) + p4_array(n+1) * lt_c + p5_array(n+1) 
    		pc3 = p1_array(n+2) * ( lt_c**4.0 ) + p2_array(n+2) * ( lt_c**3.0 ) + p3_array(n+2) * ( lt_c**2.0 ) + p4_array(n+2) * lt_c + p5_array(n+2) 
    		pc4 = p1_array(n+3) * ( lt_c**4.0 ) + p2_array(n+3) * ( lt_c**3.0 ) + p3_array(n+3) * ( lt_c**2.0 ) + p4_array(n+3) * lt_c + p5_array(n+3) 
    		pc5 = p1_array(n+4) * ( lt_c**4.0 ) + p2_array(n+4) * ( lt_c**3.0 ) + p3_array(n+4) * ( lt_c**2.0 ) + p4_array(n+4) * lt_c + p5_array(n+4) 
			
			ltr = log10(ratio)
			pc = pc1 * ( ltr**4.0 ) + pc2 * ( ltr**3.0 ) + pc3 * ( ltr**2.0 ) + pc4 * ltr + pc5
				
		else
    		pc1 = p6_array(n) * ( lt_c**2.0 ) + p7_array(n) * lt_c + p8_array(n) 
    		pc2 = p6_array(n+1) * ( lt_c**2.0 ) + p7_array(n+1) * lt_c + p8_array(n+1) 
    		pc3 = p6_array(n+2) * ( lt_c**2.0 ) + p7_array(n+2) * lt_c + p8_array(n+2) 
    		pc4 = p6_array(n+3) * ( lt_c**2.0 ) + p7_array(n+3) * lt_c + p8_array(n+3) 
    		pc5 = p6_array(n+4) * ( lt_c**2.0 ) + p7_array(n+4) * lt_c + p8_array(n+4) 
			
			ltr = log10(ratio)
			pc = pc1 * ( ltr**4.0 ) + pc2 * ( ltr**3.0 ) + pc3 * ( ltr**2.0 ) + pc4 * ltr + pc5
		

		endif 
		
		eta = 10**pc

    endif

	! Use the symmetry axis if needed
    if (mirror) then
        specpair_fast = chi_g - eta
	else
	    specpair_fast = eta
	endif

	end function specpair_fast	

	
	
!-------------------------------------------------------------------------------
! Bessel function of fractional order
!-------------------------------------------------------------------------------	
	
	
subroutine bessik(x,xnu,ri,rk,rip,rkp)

    implicit none

    integer :: MAXIT
    real(p_double) :: ri,rip,rk,rkp,x,xnu,XMIN
    real(p_double) :: EPS,FPMIN,PI
parameter (EPS=1.d-10,FPMIN=1.d-30,MAXIT=10000,XMIN=2.d0,PI=3.141592653589793d0)
	
    integer :: i,l,nl
    real(p_double) :: a,a1,b,c,d,del,del1,delh,dels,e,f,fact, &
fact2,ff,gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2, &
qnew,ril,ril1,rimu,rip1,ripl,ritemp,rk1,rkmu,rkmup, &
rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2

!print *, "x = ", x
!print *, "xnu = ", xnu

if (x < 100.0d0) then 

    ! if(x <= 0. .or. xnu .lt. 0.) pause "bad arguments in bessik"

nl=int(xnu+.5d0) 				!nl is the number of downward recurrences
								!of the I’s and upward recurrences
								!of K’s. xmu lies between −1/2 and
								!1/2.
xmu=xnu-nl
xmu2=xmu*xmu
xi=1.d0/x
xi2=2.d0*xi	
h=xnu*xi 						!Evaluate CF1 by modiﬁed Lentz’s method

if(h.lt.FPMIN)h=FPMIN 			!(§5.2).
b=xi2*xnu
d=0.d0
c=h
do i=1,MAXIT
	b=b+xi2
	d=1.d0/(b+d) 					!Denominators cannot be zero here, so no
	c=b+1.d0/c 						!need for special precautions.
	del=c*d
	h=del*h
	if(abs(del-1.d0).lt.EPS)goto 1
enddo 
!pause "x too large in bessik; try asymptotic expansion"
print *,"x too large in bessik; try asymptotic expansion"
1 continue
	ril=FPMIN
	ripl=h*ril						!Initialize Iν and Iν for downward recurrence.
	ril1=ril 						!Store values for later rescaling.
	rip1=ripl
	fact=xnu*xi
	do l=nl,1,-1
		ritemp=fact*ril+ripl
		fact=fact-xi
		ripl=fact*ritemp+ril
		ril=ritemp
	enddo 
	f=ripl/ril							!Now have unnormalized Iµ and Iµ.
    if(x.lt.XMIN) then					!Use series.
		x2=.5d0*x
		pimu=PI*xmu
		if(abs(pimu).lt.EPS)then
			fact=1.d0
		else
			fact=pimu/sin(pimu)
		endif
		d=-log(x2)
		e=xmu*d
		if(abs(e).lt.EPS)then
			fact2=1.d0
		else
			fact2=sinh(e)/e
		endif
		call beschb(xmu,gam1,gam2,gampl,gammi) 		!Chebyshev evaluation of Γ1 and Γ2.
		ff=fact*(gam1*cosh(e)+gam2*fact2*d) 		!f0.
		sum=ff
		e=exp(e)
		p=0.5d0*e/gampl 							!p0.
		q=0.5d0/(e*gammi) 							!q0.
		c=1.d0
		d=x2*x2
		sum1=p
		do i=1,MAXIT
			ff=(i*ff+p+q)/(i*i-xmu2)
			c=c*d/i
			p=p/(i-xmu)
			q=q/(i+xmu)
			del=c*ff
			sum=sum+del
			del1=c*(p-i*ff)  
			sum1=sum1+del1
			if(abs(del).lt.abs(sum)*EPS) goto 2
		enddo 
!		pause "bessk series failed to converge"
		print *,"bessk series failed to converge"
2 continue
		rkmu=sum
		rk1=sum1*xi2  
	else 							!Evaluate CF2 by Steed’s algorithm (§5.2),
									!which is OK because there can be no
									!zero denominators.
		b=2.d0*(1.d0+x)
		d=1.d0/b
		delh=d
		h=delh
		q1=0.d0 					!Initializations for recurrence (6.7.35).
		q2=1.d0
		a1=.25d0-xmu2
		c=a1
		q=c 						!First term in equation (6.7.34).
		a=-a1
		s=1.d0+q*delh
		do i=2,MAXIT
			a=a-2*(i-1)
			c=-a*c/i
			qnew=(q1-b*q2)/a
			q1=q2
			q2=qnew
			q=q+c*qnew
			b=b+2.d0
			d=1.d0/(b+a*d)
			delh=(b*d-1.d0)*delh
			h=h+delh
			dels=q*delh
			s=s+dels
			if(abs(dels/s).lt.EPS)goto 3 	!Need only test convergence of sum since
		enddo    							!CF2 itself converges more quickly.
!		pause "pause 3"
		print *,"bessik: failure to converge in cf2"
3 continue
		h=a1*h
		rkmu=sqrt(PI/(2.d0*x))*exp(-x)/s 	!Omit the factor exp(−x) to scale all the
											!returned functions by exp(x) for x ≥XMIN.
		rk1=rkmu*(xmu+x+.5d0-h)*xi
	endif
	rkmup=xmu*xi*rkmu-rk1
	rimu=xi/(f*rkmu-rkmup) 					!Get Iµ from Wronskian.
	ri=(rimu*ril1)/ril 						!Scale original Iν and Iν

	rip=(rimu*rip1)/ril
	do  i=1,nl 							!Upward recurrence of Kν.
		rktemp=(xmu+i)*xi2*rk1+rkmu
		rkmu=rk1
		rk1=rktemp
enddo 
rk=rkmu
rkp=xnu*xi*rkmu-rk1
else
	ri = 0.0d0
	rk = 0.0d0
	rip = 0.0d0
	rkp = 0.0d0
endif
    
end  subroutine bessik

!-------------------------------------------------------------------------------
! Chebyshev expansion
!-------------------------------------------------------------------------------	
	
subroutine beschb(x,gam1,gam2,gampl,gammi)

    implicit none

    integer :: NUSE1,NUSE2
    real(p_double) :: gam1,gam2,gammi,gampl,x
parameter (NUSE1=5,NUSE2=5)
! USES chebev
!Evaluates \Gamma1 and \Gamma2 by Chebyshev expansion for |x| ≤ 1/2. Also returns 1/\Gamma(1 + x) and
!1/\Gamma(1 − x). If converting to double precision, set NUSE1 = 7, NUSE2 = 8.
    real(p_double) :: xx,c1(7),c2(8)
SAVE c1,c2
DATA c1/-1.142022680371168d0,6.5165112670737d-3, &
3.087090173086d-4,-3.4706269649d-6,6.9437664d-9, &
3.67795d-11,-1.356d-13/
DATA c2/1.843740587300905d0,-7.68528408447867d-2, &
1.2719271366546d-3,-4.9717367042d-6,-3.31261198d-8, &
2.423096d-10,-1.702d-13,-1.49d-15/
xx=8.d0*x*x-1.d0 						!Multiply x by 2 to make range be −1 to 1, and then
										!apply transformation for evaluating even Chebyshev series.
gam1=chebev(-1.d0,1.d0,c1,NUSE1,xx)
gam2=chebev(-1.d0,1.d0,c2,NUSE2,xx)
gampl=gam2-x*gam1
gammi=gam2+x*gam1

end	subroutine beschb


!-------------------------------------------------------------------------------
! Chebyshev evaluation
!-------------------------------------------------------------------------------
  

function chebev(a,b,c,m,x)

    implicit none

    integer :: m
    real(p_double) :: chebev,a,b,x,c(m)
!Chebyshev evaluation: All arguments are input. c(1:m) is an array of Chebyshev coeﬃcients, 
!the ﬁrst m elements of c output from chebft (which must have been called with the same a and b). 
!The Chebyshev polynomial is evaluated at a point and the result is returned as the funccao value.
    integer :: j
    real(p_double) :: d,dd,sv,y,y2

! if ((x-a)*(x-b).gt.0.) pause " not in range in chebev"

d=0.
dd=0.
y=(2.*x-a-b)/(b-a) 				!Change of variable.
y2=2.*y
do j=m,2,-1 					!Clenshaw’s recurrence.
	sv=d
	d=y2*d-dd+c(j)
	dd=sv
enddo 
chebev=y*d-dd+0.5*c(1) 			!Last step is diﬀerent.

end function chebev
  

function arth_d(first,increment,n)      

    implicit none

integer, parameter :: NPAR_ARTH=16,NPAR2_ARTH=8

real(p_double), intent(in) :: first,increment
integer, intent(in) :: n
real(p_double), dimension(n) :: arth_d

integer :: k,k2
real(p_double) :: temp

if (n > 0) arth_d(1)=first

if (n <= NPAR_ARTH) then
  do k=2,n
    arth_d(k)=arth_d(k-1)+increment
  end do
else
   do k=2,NPAR2_ARTH
	 arth_d(k)=arth_d(k-1)+increment
   end do
   temp=increment*NPAR2_ARTH
   k=NPAR2_ARTH
   do
	  if (k >= n) exit
	  k2=k+k
	  arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
	  temp=temp+temp
	  k=k2
   end do
end if

end function arth_d


!-------------------------------------------------------------------------------
! Numerical integration using an extended trapezoidal rule for a function f(x,y)
!-------------------------------------------------------------------------------
  
  
subroutine trapzd1(func,a,b,s,n)

implicit none

real(p_double), intent(in) :: a,b
real(p_double), intent(inout) :: s
integer, intent(in) :: n
interface
	function func(x,y)
	use m_system
	real(p_double), dimension(:), intent(in) :: x
	real(p_double), intent(in) :: y
	real(p_double), dimension(size(x)) :: func
	end function func
end interface
!This routine computes the nth stage of reﬁnement of an extended trapezoidal rule. func is
!input as the name of the function to be integrated between limits a and b, also input. When
!called with n=1, the routine returns as s the crudest estimate of int(a,b)f(x)dx. Subsequent
!calls with n=2,3,... (in that sequential order) will improve the accuracy of s by adding 2
!n-2 additional interior points. s should not be modiﬁed between sequential calls.
real(p_double) :: del,fsum
integer :: it

s=0.0d0

if (n == 1) then
	s=0.5d0*(b-a)*sum(func( (/ a,b /),a ))
else
	it=2**(n-2)
	del=(b-a)/it 							!This is the spacing of the points to be added.
	fsum=sum(func(arth_d(a+0.5d0*del,del,it),a))
	s=0.5d0*(s+del*fsum)*2.0d0 					!This replaces s by its reﬁned value.
	
end if

end subroutine trapzd1 

!-------------------------------------------------------------------------------
! Numerical integration using an extended trapezoidal rule for a function f(x)
!-------------------------------------------------------------------------------
  
  
subroutine trapzd2(func,a,b,s,n)

implicit none

real(p_double), intent(in) :: a,b
real(p_double), intent(inout) :: s
integer, intent(in) :: n
interface
	function func(x)
	use m_system
	real(p_double), dimension(:), intent(in) :: x
	real(p_double), dimension(size(x)) :: func
	end function func
end interface
!This routine computes the nth stage of reﬁnement of an extended trapezoidal rule. func is
!input as the name of the function to be integrated between limits a and b, also input. When
!called with n=1, the routine returns as s the crudest estimate of int(a,b)f(x)dx. Subsequent
!calls with n=2,3,... (in that sequential order) will improve the accuracy of s by adding 2
!n-2 additional interior points. s should not be modiﬁed between sequential calls.
real(p_double) :: del,fsum
integer :: it

s=0.0d0

if (n == 1) then
	s=0.5d0*(b-a)*sum(func( (/ a,b /) ))
else
	it=2**(n-2)
	del=(b-a)/it 							!This is the spacing of the points to be added.
	fsum=sum(func(arth_d(a+0.5d0*del,del,it)))
	s=0.5d0*(s+del*fsum)*2.0d0 					!This replaces s by its reﬁned value.
	
end if

end subroutine trapzd2  



!-------------------------------------------------------------------------------
!  intermediate function  : J_1(y)
!-------------------------------------------------------------------------------
 
function J_1(t,y)

    implicit none

real(p_double), dimension(:), intent(in) :: t
real(p_double), intent(in) :: y
real(p_double), dimension(size(t)) :: J_1
real(p_double) :: ri1,rk1,rip1,rkp1
integer :: k

do k=1, size(t)
  call bessik(t(k),2.d0/3.d0,ri1,rk1,rip1,rkp1)
  J_1(k) = (rk1**2)*t(k)/sqrt( (t(k)/y)**(2.d0/3.d0) -1.0d0 )/(3*y**2)
enddo    
   

end function J_1


!-------------------------------------------------------------------------------
!  intermediate function  : J_2(y)
!-------------------------------------------------------------------------------
 
function J_2(t,y)

    implicit none

real(p_double), dimension(:), intent(in) :: t
real(p_double), intent(in) :: y
real(p_double), dimension(size(t)) :: J_2
real(p_double) :: ri1,rk1,rip1,rkp1
integer :: k

do k=1, size(t)
  call bessik(t(k),1.d0/3.d0,ri1,rk1,rip1,rkp1)
  J_2(k) = (rk1**2)*((t(k)/y)**(1.d0/3.0d0))*sqrt( (t(k)/y)**(2.d0/3.d0) -1.0d0 )/(3*y)
enddo    
   

end function J_2


!-------------------------------------------------------------------------------
!   intermediate function  : J_3(y)
!-------------------------------------------------------------------------------
 
function J_3(t,y)

    implicit none

real(p_double), dimension(:), intent(in) :: t
real(p_double), intent(in) :: y
real(p_double), dimension(size(t)) :: J_3
real(p_double) :: ri1,rk1,rip1,rkp1
integer :: k

do k=1, size(t)
  call bessik(t(k),1.d0/3.d0,ri1,rk1,rip1,rkp1)
  J_3(k) = (rk1**2)*t(k)/sqrt( (t(k)/y)**(2.d0/3.d0) -1.0d0 )/(3*y**2)
enddo    
   

end function J_3


!-------------------------------------------------------------------------------
!   intermediate function  : K_13(y)
!-------------------------------------------------------------------------------
 
function K_13(t)

    implicit none

real(p_double), dimension(:), intent(in) :: t
real(p_double), dimension(size(t)) :: K_13
real(p_double) :: ri1,rk1,rip1,rkp1
integer :: k

do k=1, size(t)
  call bessik(t(k),1.d0/3.d0,ri1,rk1,rip1,rkp1)
  K_13(k) = rk1
enddo  
   

end function K_13


!-------------------------------------------------------------------------------
!  function synchrotron emissivity  : F(eta,chi)
!-------------------------------------------------------------------------------

function F_pkt(eta,chi)

    implicit none

real(p_double), intent(in) :: eta,chi
real(p_double) :: F_pkt, y, x, J1eval, J2eval, J3eval, sum_bess
    integer :: n

n = 10

y=2.d0*chi/(3.d0*eta*(eta-2.d0*chi))
x=2.0d0*chi/eta

call trapzd1(J_1,y,y+10.0d0,J1eval,n)
call trapzd1(J_2,y,y+10.0d0,J2eval,n)
call trapzd1(J_3,y,y+10.0d0,J3eval,n)
sum_bess=(1.d0+1.d0/(1.d0-x)**2)*J1eval+2.d0*J2eval/(1.d0-x)+J3eval*(x**2)/((1.d0-x)**2)

    F_pkt=sum_bess*8.d0*(chi**2)/(3.d0*sqrt(3.d0)*pi*(eta**4))

end function F_pkt


!-------------------------------------------------------------------------------
!  function pair emission probability  : F(eta,chi)
!-------------------------------------------------------------------------------

function F_pair(eta,chi)

    implicit none

real(p_double), intent(in) :: eta,chi
real(p_double) :: F_pair, nu, delta, K_13eval, ri2, rk2, rip2, rkp2

    integer :: n

n = 10

nu = eta / chi
delta = 2.d0/(3.d0*chi*nu*(1-nu))

call bessik(delta,2.d0/3.d0,ri2,rk2,rip2,rkp2)

call trapzd2(K_13,delta,delta+10.0d0,K_13eval,n)

F_pair=((1-nu)/nu + nu/(1-nu))*rk2+K_13eval

end function F_pair
  


!-------------------------------------------------------------------------------
!  function Peak of pkt synchrotron spectrum
!-------------------------------------------------------------------------------

function peak_pkt(eta)

    implicit none

integer :: i, loc
real(p_double), intent(in) :: eta
real(p_double) :: peak_pkt, coef

real(p_double), dimension(19), parameter :: eta_array = (/1.0d-4,1.0d-2, 1.0d-1, 0.2d0, 0.5d0, 1.0d0, 2.0d0, 3.0d0, 5.0d0, 10.d0, 15.0d0, 20.0d0, &
			30.0d0, 40.0d0, 60.0d0, 80.0d0, 100.0d0, 200.0d0, 300.0d0/)
real(p_double), dimension(19), parameter :: peak_array = (/0.918d0,0.914d0, 0.883d0, 0.855d0, 0.795d0, 0.732d0, 0.656d0, 0.610d0, 0.568d0, 0.574d0, 0.582d0, 0.587d0, &
			0.591d0, 0.596d0, 0.599d0, 0.6d0, 0.6d0, 0.6d0, 0.6d0/)

    if (eta <= 1.0d-4) then
	peak_pkt = 0.918d0
else 
	if ((eta.gt.1.0d-4).and.(eta.lt.300.0d0)) then
           loc = 15
	   do i=1,18
              if (eta <= eta_array(i+1) .and. eta >= eta_array(i)) then
			 loc = i
			exit
		  endif
	   enddo   
	   coef = (peak_array(loc+1)-peak_array(loc))/(eta_array(loc+1)-eta_array(loc))
	   peak_pkt = peak_array(loc) + coef*(eta-eta_array(loc))
	else   
		peak_pkt = 0.605d0
	endif	
endif	

end function peak_pkt

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  function interpolation of integration of the pkt synchrotron spectrum  
!-------------------------------------------------------------------------------

function interp_spec_pkt(eta)

    implicit none

integer :: i, loc
real(p_double), intent(in) :: eta
real(p_double) :: interp_spec_pkt, coef

real(p_double), dimension(18), parameter :: eta_array = (/1.0d-4, 1.0d-3, 1.0d-2, 1.0d-1, 0.5d0, 1.0d0, 2.0d0, 5.0d0, 1.0d1, 15.0d0, 20.0d0, &
			30.0d0, 40.0d0, 60.0d0, 80.0d0, 100.0d0, 200.0d0, 300.0d0/)
real(p_double), dimension(18), parameter :: int_eval_array = (/5.236d0, 5.231d0, 5.19d0, 4.87d0, 4.177d0, 3.75d0, 3.284d0, 2.659d0, 2.217d0, 1.98d0, 1.82d0, &
			1.61d0, 1.48d0, 1.3d0, 1.19d0, 1.11d0, 0.89d0, 0.78d0/)

    if (eta <= 1.0d-4) then
	interp_spec_pkt = 5.236d0
else
        loc = 17
	do i=1,17
            if (eta <= eta_array(i+1) .and. eta >= eta_array(i)) then
			loc = i
			exit
		endif
	enddo
	coef = (int_eval_array(loc+1)-int_eval_array(loc))/(eta_array(loc+1)-eta_array(loc))
	interp_spec_pkt = int_eval_array(loc) + coef*(eta-eta_array(loc))
endif

end function interp_spec_pkt

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  function polynomial expansion for synchrontron cumulative probability for y<<1
!  			where y = 2*chi/(3*eta(eta-2*chi))           ; Marija
!-------------------------------------------------------------------------------

function poly_sync_pkt(eta, chi)

real(p_double), intent(in) :: eta, chi
real(p_double) :: poly_sync_pkt
real(p_double) :: a, b, d, I12, I34

    a = 0.921
    b = 0.307
    d = 2.0*chi/eta
    I12 = 6.0 - d + (1.0/3.0)* (d**2.0) + 0.069023569023569 * (d**3.0)
    I34 =  6.0 - d + (1.0/3.0)* (d**2.0) + 0.072016460905350 * (d**3.0)
	poly_sync_pkt = 0.764540211031963 * (eta**(-1.0/3.0)) * (d**(1.0/3.0)) *( a * I12 + b * I34 ) 
	! this is (2/pi)*3^(1/6) * eta^(-1/3) d^(1/3) * [ a(I1+I2) + b (I3+I4) ]

end function poly_sync_pkt;

!-------------------------------------------------------------------------------
!  function estimates the integral based on improved polynomial expansion; Marija
!             The output is normalized to the integral from 0 to eta/2
!-------------------------------------------------------------------------------

function integral_alpha_pkt(eta, chi)

real(p_double), intent(in) :: eta, chi
real(p_double) :: integral_alpha_pkt
real(p_double) :: chi_y1, delta, total, int

	total = interp_spec_pkt(eta)
	
	! this is the chi for which y = 1
	chi_y1 = 3.0 * eta * eta / ( 2.0*(1 + 3.0 * eta) )
	
	if ( chi > chi_y1 ) then
			int = 1.0
	else  
			
			delta = poly_sync_pkt(eta, chi_y1) - total
			int = ( poly_sync_pkt(eta, chi) - (chi * delta / chi_y1) ) / total
						
	endif
	
	integral_alpha_pkt = int

end function integral_alpha_pkt

!-------------------------------------------------------------------------------
!  set of functions to get new chi from eta and one random value; Marija;
!  There are 3 functions: eta> 3.3 (high), 0.01< eta < 3.3 (midrange) and eta < 0.01 (low);   
!-------------------------------------------------------------------------------
function interpolate_photon_high(eta, ratio)	

real(p_double), intent(in) :: eta, ratio
real(p_double) :: interpolate_photon_high
real(p_double) :: chi, lt_eta, ltr, pc1, pc2, pc3, pc4, pc5, pc
integer :: n

real(p_double), dimension(20), parameter :: p1_array = (/ 0.000394589503051,     0.001758682038381,     0.002893874817675,     0.002104730355267,     -0.003634049745981,  &
		   -0.055937387233642,     -0.065886320279835,     -0.026147110013602,     -0.002940621306168,     -0.003894502232776,  &
		   -1.048747279807481,     -0.872789273537110,     -0.267400248879477,     -0.034452832734123,     -0.005414601503678,   &
		  1875.071802019952656,     300.495149184593060,     16.762405490950311,     0.373442696138872,     -0.001652051103320 /)
real(p_double), dimension(20), parameter :: p2_array = (/ -0.003971150489006,     -0.018055330405033,     -0.030513120343084,     -0.023030679161369,     0.042273956423612, &
           0.329951816637183,     0.370401763115319,     0.126156506804666,     0.000407195758592,     0.042908607010348,   &
          -0.354179154844088,     1.286949411802133,     0.739449801699535,     0.114599218464971,     0.049666212297222,   &
          -16743.456812191390782,     -2690.227130010071960,     -151.629508314305980,     -3.554314465956319,     0.015338323224990 /)
real(p_double), dimension(20), parameter :: p3_array = (/ 0.015624480455891,     0.072569603108912,     0.126108887187262,     0.098804713961978,     -0.196613435156459,  &
	      -0.092237867707433,     0.054258855096189,     0.207495688575473,     0.146153041062157,     -0.189146562438188,  &
	      51.069523904044814,     29.535660711749539,     6.237888502630087,     0.654622916464811,     -0.174825363262565,   &
	      56887.183787550107809,     9184.043012321235437,     526.266972325473603,     13.271756641941396,     -0.054733202346443 /)
real(p_double), dimension(20), parameter :: p4_array = (/ -0.029199070329492,     -0.138548698396030,     -0.247451431218712,     -0.200886031023574,     1.438475028852042,   &
         -2.638987961872001,     -3.620798753363849,     -2.020650203097647,     -0.611879714396155,     1.401780723462706,  &
         -195.122632154715916,     -125.818836023634518,     -30.690761393156599,     -3.569844737674069,     1.287758082995097,  &
         -89508.217916543173487,     -14572.859228520872421,     -856.519716321191368,     -23.819455788619784,     1.091395618969734 /)
real(p_double), dimension(20), parameter :: p5_array = (/ 0.048890836069625,     0.244727685140817,     0.454442708483244,     3.374002199572533,     -0.719404122634595,   &
		  6.615523470324536,     9.858398511739829,     6.018148533542798,     4.902198980643337,     -0.551056458816563,   &
	      117.323516878384979,     83.233549894578488,     24.205231838998493,     6.911554027078656,     -0.466984430164604,   &
	      55498.371969807310961,     9008.148552239785204,     526.826181903363249,     18.719241080848093,     -0.363030323053338 /)
    
    lt_eta = log10(eta)
    ltr = log10(ratio)
    
   ! Here we look for the right section of coefficients to use; 

    if ( ratio < 0.4 ) then 
    	n=1
    else if ( ratio < 0.7 ) then
    		n = 6
    	  else if ( ratio < 0.9 ) then
    				n = 11
    	  		else  
    				n = 16
!    	  		endif
  ! these two are commented because the compiler complained. Marija  		
 !   	  endif
    endif
    
    pc1 = p1_array(n) * ( lt_eta**4.0 ) + p2_array(n) * ( lt_eta**3.0 ) + p3_array(n) * ( lt_eta**2.0 ) + p4_array(n) * lt_eta + p5_array(n) 
    pc2 = p1_array(n+1) * ( lt_eta**4.0 ) + p2_array(n+1) * ( lt_eta**3.0 ) + p3_array(n+1) * ( lt_eta**2.0 ) + p4_array(n+1) * lt_eta + p5_array(n+1) 
    pc3 = p1_array(n+2) * ( lt_eta**4.0 ) + p2_array(n+2) * ( lt_eta**3.0 ) + p3_array(n+2) * ( lt_eta**2.0 ) + p4_array(n+2) * lt_eta + p5_array(n+2) 
    pc4 = p1_array(n+3) * ( lt_eta**4.0 ) + p2_array(n+3) * ( lt_eta**3.0 ) + p3_array(n+3) * ( lt_eta**2.0 ) + p4_array(n+3) * lt_eta + p5_array(n+3) 
    pc5 = p1_array(n+4) * ( lt_eta**4.0 ) + p2_array(n+4) * ( lt_eta**3.0 ) + p3_array(n+4) * ( lt_eta**2.0 ) + p4_array(n+4) * lt_eta + p5_array(n+4) 

	pc = pc1 * ( ltr**4.0 ) + pc2 * ( ltr**3.0 ) + pc3 * ( ltr**2.0 ) + pc4 * ltr + pc5

    chi = 10**pc

	interpolate_photon_high = chi
	
end function interpolate_photon_high

!-------------------------------------------------------------------------------
!  Find chi for midrange photons 0.01< eta < 3.3  ; Marija
!-------------------------------------------------------------------------------
function interpolate_photon_midrange(eta, ratio)	

real(p_double), intent(in) :: eta, ratio
real(p_double) :: interpolate_photon_midrange
real(p_double) :: chi, lt_eta, ltr, pc1, pc2, pc3, pc4, pc5, pc
integer :: n

real(p_double), dimension(30), parameter :: p1_array = (/ 0.001773836232853,     0.004159773105190,     0.005909583673540,     0.006022338374199,     0.012028184722626,  &
	      -0.004317070831837,     0.009330893070144,     0.022133372122301,     0.015757086764085,     0.013817061729165,   &
	      -0.980486796007142,     -0.682697620064631,     -0.161298903150042,     -0.005809661704548,     0.012866929690746,  &
	      -0.000482475937002,     -0.000985616750807,     0.006209647636546,     0.013027947138266,     0.019074621074484,  &
	       0.001265577345342,     0.013336749438294,     0.048791963835746,     0.067308775582128,     0.044043687269626,   &
	      -0.000055358248880,     -0.001348721803424,     -0.012275052087300,     -0.045221323100608,     -0.033430575699799 /)
real(p_double), dimension(30), parameter :: p2_array = (/ -0.000170639654639,     0.007170182865574,     0.029058662140918,     0.034982833305507,     0.016674896115385,  &
	       0.268338782993809,     0.469351068095490,     0.351483214319252,     0.141973073260982,     0.030569503524505,   &
	      -4.909508871257551,     -2.996785964714983,     -0.507820495166353,     0.048413754117628,     0.026789393408900,   &
	      -0.004329396889427,     -0.024295168496621,     -0.033052751294962,     -0.017905509970837,     0.018581438999043,   &
	       0.004381022575082,     0.048108792452574,     0.186206593872071,     0.268302058587707,     0.154364638546128,  &
	      -0.000148670396644,     -0.003846576613533,     -0.038205552382314,     -0.164785712598556,     -0.161041306668578 /)
real(p_double), dimension(30), parameter :: p3_array = (/ 0.018175348373514,     0.072438926956503,     0.105531563563381,     0.067040495583811,     -0.170000259566593,   &
		  0.223262330890935,     0.471684923546082,     0.388922426066648,     0.153796985221281,     -0.160391408919875,   &
		  -31.680084424698883,     -20.084722542030761,  	     -4.510229318398666,     -0.359559263573851,     -0.180449508408224,   &
		  -0.012347312978628,     -0.090172836468116,     -0.213250060611825,     -0.192351240589369,    -0.229109462035024,   &
		   0.003111666251559,     0.039690828311627,     0.185564686028165,     0.338027042839089,     0.028759127464761,    &
		   0.000023481513433,     0.000022801554947,     -0.007165145535375,     -0.081707990054295,     -0.316879491793383 /)
real(p_double), dimension(30), parameter :: p4_array = (/ -0.031911342972680,     -0.144318629353999,     -0.245920478028274,     -0.190427282615784,     1.421796540765927,   &
		  -2.545220550191603,     -3.553366289404980,     -2.017380580443948,     -0.610503032302649,     1.383423724377177,  &
	      -109.316616098549332,     -72.579453825704149,     -18.597839190749774,     -2.373170424438470,     1.312928625714792,   &
	      -0.015073626703808,     -0.125181013600805,     -0.346874427727219,     -0.241763740356303,     1.383511050857337,    &
	      -0.002791369594754,     -0.022584555300707,     -0.032865372915183,     0.175959378659173,     1.587766611263539,   &
	       0.000269829602639,     0.006650538450129,     0.064272326168450,     0.299244401499319,     1.624536255589060 /)
real(p_double), dimension(30), parameter :: p5_array = (/ 0.049421028623077,     0.245823831180409,     0.453711481196624,     3.371207968157663,     -0.715551819439772,  &
           6.515531213872445,     9.740563653619228,     5.966701464021353,     4.890315707496041,     -0.547768232060367,   &
           93.612358800233736,     68.274849458436137,     20.721252596432443,     6.552637908179546,     -0.476703814260746,  &
          -0.038414492310977,     -0.361223808578362,     -1.312294486900120,     -2.265078646565482,     -2.015389446419847,   &
          -0.011376879355673,     -0.133755830526832,     -0.611451162601370,     -1.328331100768032,     -1.556871006962085,   &
          -0.000243937957341,     -0.005785523688170,     -0.052781127864644,     -0.226485066799127,     -0.725872671153799 /)
    
    lt_eta = log10(eta)
    
    ! Here we look for the right section of interpolation coefficients to use; 
    ! Please note that second half of the coefficients (above n=16) is done for log10(1-ratio)
    ! due to the numerical precision when ratio => 1.0
    
    if ( ratio < 0.4 ) then 
    	n=1
    	ltr = log10(ratio)
    else if ( ratio < 0.7 ) then
    		n = 6
    		ltr = log10(ratio)
    	  else if ( ratio < 0.88 ) then
    				n = 11
    				ltr = log10(ratio)
    	  		else  if (ratio < 0.975) then
    					   n = 16
    					   ltr = log10(1.0-ratio)
    				  else if (ratio < 0.9993) then
    					   		n = 21
    					   		ltr = log10(1.0-ratio)
    				  		else 
    				  			n = 26
    					   		ltr = log10(1.0-ratio)
    				  		!endif
    				  !endif
    	  		!endif
    		! all this commented for the compiler not to complain
    	  !endif
    endif
    
    pc1 = p1_array(n) * ( lt_eta**4.0 ) + p2_array(n) * ( lt_eta**3.0 ) + p3_array(n) * ( lt_eta**2.0 ) + p4_array(n) * lt_eta + p5_array(n) 
    pc2 = p1_array(n+1) * ( lt_eta**4.0 ) + p2_array(n+1) * ( lt_eta**3.0 ) + p3_array(n+1) * ( lt_eta**2.0 ) + p4_array(n+1) * lt_eta + p5_array(n+1) 
    pc3 = p1_array(n+2) * ( lt_eta**4.0 ) + p2_array(n+2) * ( lt_eta**3.0 ) + p3_array(n+2) * ( lt_eta**2.0 ) + p4_array(n+2) * lt_eta + p5_array(n+2) 
    pc4 = p1_array(n+3) * ( lt_eta**4.0 ) + p2_array(n+3) * ( lt_eta**3.0 ) + p3_array(n+3) * ( lt_eta**2.0 ) + p4_array(n+3) * lt_eta + p5_array(n+3) 
    pc5 = p1_array(n+4) * ( lt_eta**4.0 ) + p2_array(n+4) * ( lt_eta**3.0 ) + p3_array(n+4) * ( lt_eta**2.0 ) + p4_array(n+4) * lt_eta + p5_array(n+4) 

	pc = pc1 * ( ltr**4.0 ) + pc2 * ( ltr**3.0 ) + pc3 * ( ltr**2.0 ) + pc4 * ltr + pc5

    chi = 10**pc

	interpolate_photon_midrange = chi
	
end function interpolate_photon_midrange

!-------------------------------------------------------------------------------
!  Find chi for midrange photons  eta < 0.01 ; Marija
!-------------------------------------------------------------------------------
function interpolate_photon_low(eta, ratio)	

real(p_double), intent(in) :: eta, ratio
real(p_double) :: interpolate_photon_low
real(p_double) :: chi, lt_eta, ltr, pc1, pc2, pc3, pc4, pc5, pc
integer :: n

real(p_double), dimension(30), parameter :: p1_array = (/ 0.000263811201959,     0.000553145340450,     0.000095369097674,     -0.000411875299010,     -0.001853132487289,  &
		   -0.034968036509529,     -0.048613728699486,     -0.027108238057496,     -0.007689397564446,     -0.002662591941810,  &
		   -0.471389071640141,     -0.366882206365127,     -0.112218758479580,     -0.017612011812319,     -0.003089636664014,  &
		   -0.000026474752394,     -0.000220282443738,     -0.000682523450780,     0.000433200917317,     -0.001573728870609,   &
		   -0.000005530223021,     -0.000071572557034,     -0.000306087538604,     0.000833871524539,     -0.001422824476591,   &
		   -0.000000140633974,     -0.000004012965775,     0.000014133707467,     0.001516930650954,     -0.000868627294210 /)
real(p_double), dimension(30), parameter :: p2_array = (/ 0.001759803729180,     0.002219230270537,     -0.003710995832019,     -0.007443369492250,     -0.025954464715172,   &
	      -0.481290048732986,     -0.671043559591993,     -0.375418889566618,     -0.106850616220667,     -0.037061410137816,   &
	      -7.275620737513042,     -5.506669738559794,     -1.643778940803329,     -0.252565594681990,     -0.043269949679557,   &
	      -0.000356276154489,     -0.003002703096092,     -0.009471182536737,     0.006046013708690,     -0.021917477241251,    &
	      -0.000082444971341,     -0.001046577943293,     -0.004492265877667,     0.011365874881429,     -0.019913641794140,   &
	      -0.000002547037527,     -0.000068162145320,     0.000053060484331,     0.020896093063842,     -0.012292213053433 /)
real(p_double), dimension(30), parameter :: p3_array = (/ 0.006986723640187,     0.009997910101316,     -0.014612961168528,     -0.034208021629582,     -0.135324497875304,   &
          -2.503213435436670,     -3.500587251076293,     -1.964980868889572,     -0.561163702345697,     -0.194881296072327,  &
          -42.093636131306923,     -31.107872464612527,     -9.085101831836093,     -1.368415148610610,     -0.228959888084581,  &
          -0.001824697740166,     -0.015564332548708,     -0.049964958271390,     0.031574604130121,     -0.115395098457585,    &
          -0.000456844928734,     -0.005715449721112,     -0.024696987365157,     0.058773043081745,     -0.105091818810049,    &
          -0.000016740343327,     -0.000428850551023,     -0.000563308599435,     0.108577421492149,     -0.065824907236507 /)
real(p_double), dimension(30), parameter :: p4_array = (/ 0.005335131360358,     -0.022501648387618,     -0.110725979737021,     -0.136216820923786,     1.666869897886332,   &
          -5.855109866397301,     -8.212834724838951,     -4.625499915752781,     -1.325329840321259,     1.539440714165051,   &
          -108.585158882810205,     -78.638367024033329,     -22.526356253671199,     -3.331216319051826,     1.455465999810217,  &
          -0.004237355165717,     -0.036525658976146,     -0.119237916752667,     0.073381346509617,     1.726770340767212,   &
          -0.001122548215616,     -0.013891027040214,     -0.060598648484234,     0.137146528360352,     1.751177760085899,  &
          -0.000048009457794,     -0.001192614380272,     -0.003528934286669,     0.253157886207936,     1.841329935045384 /)
real(p_double), dimension(30), parameter :: p5_array = (/ 0.012440512264447,     0.058208450503953,     0.135433713807426,     3.159377819076660,     -0.576999863431748,  &
	      5.328981984837045,     8.165820883577529,     5.169630202193487,     4.713848412633423,     -0.375348777377691,  &
	     109.908794069604866,     75.340162347559087,     21.372576915740748,     6.468554084565493,     -0.302616117826916,  &
         -0.034560041731897,     -0.326548751981733,     -1.219198942872906,     -2.148249186408203,     -1.782748307370904,  &
         -0.009402387169714,     -0.115966993497007,     -0.573762433354587,     -1.290381599174277,     -1.365467788910251,  &
         -0.000418029920662,     -0.010630206487499,     -0.103007290874155,     -0.336591670303331,     -0.624414723674595 /)
    
    lt_eta = log10(eta)
    
    ! Here we look for the right section of interpolation coefficients to use; 
    ! Please note that second half of the coefficients (above n=16) is done for log10(1-ratio)
    ! due to the numerical precision when ratio => 1.0
    
    if ( ratio < 0.4 ) then 
    	n=1
    	ltr = log10(ratio)
    else if ( ratio < 0.7 ) then
    		n = 6
    		ltr = log10(ratio)
    	  else if ( ratio < 0.88 ) then
    				n = 11
    				ltr = log10(ratio)
    	  		else  if (ratio < 0.975) then
    					   n = 16
    					   ltr = log10(1.0-ratio)
    				  else if (ratio < 0.9993) then
    					   		n = 21
    					   		ltr = log10(1.0-ratio)
    				  		else 
    				  			n = 26
    					   		ltr = log10(1.0-ratio)
    				  		!endif
    				  !endif
    	  		!endif
    		! commented for compiler not to complain
    	  !endif
    endif
    
    pc1 = p1_array(n) * ( lt_eta**4.0 ) + p2_array(n) * ( lt_eta**3.0 ) + p3_array(n) * ( lt_eta**2.0 ) + p4_array(n) * lt_eta + p5_array(n) 
    pc2 = p1_array(n+1) * ( lt_eta**4.0 ) + p2_array(n+1) * ( lt_eta**3.0 ) + p3_array(n+1) * ( lt_eta**2.0 ) + p4_array(n+1) * lt_eta + p5_array(n+1) 
    pc3 = p1_array(n+2) * ( lt_eta**4.0 ) + p2_array(n+2) * ( lt_eta**3.0 ) + p3_array(n+2) * ( lt_eta**2.0 ) + p4_array(n+2) * lt_eta + p5_array(n+2) 
    pc4 = p1_array(n+3) * ( lt_eta**4.0 ) + p2_array(n+3) * ( lt_eta**3.0 ) + p3_array(n+3) * ( lt_eta**2.0 ) + p4_array(n+3) * lt_eta + p5_array(n+3) 
    pc5 = p1_array(n+4) * ( lt_eta**4.0 ) + p2_array(n+4) * ( lt_eta**3.0 ) + p3_array(n+4) * ( lt_eta**2.0 ) + p4_array(n+4) * lt_eta + p5_array(n+4) 

	pc = pc1 * ( ltr**4.0 ) + pc2 * ( ltr**3.0 ) + pc3 * ( ltr**2.0 ) + pc4 * ltr + pc5

    chi = 10**pc

	interpolate_photon_low = chi
	
end function interpolate_photon_low


!-------------------------------------------------------------------------------
!  function find photon chi with eta and random value between 0 and 1
!-------------------------------------------------------------------------------

! function find_photon_chi(eta)
! 
! real(p_double), intent(in) :: eta
! real(p_double) :: find_photon_chi
! real(p_double) :: chi, varand, first_guess, a, b, d, total, delta_y, derivative, int
! 
! 	a = 0.921
!     b = 0.307
!     total = interp_spec_pkt(eta)
! 	call random_number( varand )	
! 	
! 	! inverse from the non-improved polynomial expansion
! 	first_guess =0.5*(eta**2 )* ( varand * total /( 0.764540211031963 * (a + b) *6.0 )) **3
!     int = integral_alpha_pkt(eta, first_guess)
!     delta_y = varand - int
!     d = 2.0 * first_guess / eta
!     derivative = 0.764540211031963*(a+b)* (eta**(-2.0/3.0))* 2.0*(2.0* (d**(-2.0/3.0))-(4.0/3.0)*(d**(1.0/3.0))) /total
!     ! this derivative is calculated by = (2/pi) * 3^(1/6) * (a+b) * 2 * eta^(-2/3)* (2 d^(-2/3)- (4/3)d^(1/3) )   
!     chi = first_guess + delta_y / derivative
!     
! 	find_photon_chi = chi
! 
! end function find_photon_chi


!-------------------------------------------------------------------------------
!  function find photon chi with eta and random value between 0 and 1
!-------------------------------------------------------------------------------

function find_photon_chi(eta)

real(p_double), intent(in) :: eta
real(p_double) :: find_photon_chi
real(p_double) :: chi, varand, first_guess, a, b, d, total, delta_y, derivative, int, chi_y1, chi_new

	a = 0.921
    b = 0.307
    total = interp_spec_pkt(eta)
	call random_number( varand )	
	
	! inverse from the non-improved polynomial expansion
	first_guess =0.5*(eta**2 )* ( varand * total /( 0.764540211031963 * (a + b) *6.0 )) **3
    int = integral_alpha_pkt(eta, first_guess)    
    delta_y = varand - int
    
    ! this is the chi for which y = 1
	chi_y1 = 3.0 * eta * eta / ( 2.0*(1 + 3.0 * eta) )
	
	! if chi << chi_y1, then we can use this expansion; otherwise, we use interpolation
	if (first_guess < 0.001 * chi_y1) then
		d = 2.0 * first_guess / eta
    	derivative = 0.764540211031963*(a+b)* (eta**(-4.0/3.0))* 2.0*(2.0* (d**(-2.0/3.0))-(4.0/3.0)*(d**(1.0/3.0))) /total
    	! this derivative is calculated by = (2/pi) * 3^(1/6) * (a+b) * 2 * eta^(-4/3)* (2 d^(-2/3)- (4/3)d^(1/3) )   
    	chi = first_guess + delta_y / derivative
    else
    	if (eta > 3.3) then
    			 chi_new = interpolate_photon_high(eta, varand)	
    	else if (eta > 0.01) then
    		           chi_new = interpolate_photon_midrange(eta, varand)
    		 else 
    		 		   chi_new = interpolate_photon_low(eta, varand)
    		 !endif
        endif
        chi = chi_new
	endif
    
	find_photon_chi = chi

end function find_photon_chi


!-------------------------------------------------------------------------------
!  function Peak of pkt pair emission spectrum
!-------------------------------------------------------------------------------

function peak_pairpkt(chi_g)

    implicit none

integer :: i, loc
real(p_double), intent(in) :: chi_g
real(p_double) :: peak_pairpkt, coef

real(p_double), dimension(46), parameter :: chi_g_array = (/ 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, &
			0.7d0, 0.8d0, 0.9d0, 1.0d0, 1.2d0, 1.5d0, 2.0d0, 2.5d0, 3.0d0, 3.5d0, 4.0d0, 4.5d0, 5.0d0, 5.5d0, 6.0d0, 6.5d0, &
			7.0d0, 7.5d0, 8.0d0, 8.5d0, 9.0d0, 9.5d0, 10.0d0, 15.0d0, 20.0d0, &
			 30.0d0, 40.0d0, 60.0d0, 80.0d0, 100.0d0, 125.0d0, 150.0d0, 175.0d0, 200.0d0, 300.0d0, 500.0d0, 1.0d3, 5.0d3, 1.0d4, 5.0d4/)
real(p_double), dimension(46), parameter :: peak_array = (/ 1.9d-12, 1.65d-6, 1.7d-4, 1.8d-3, 7.7d-3, 2.0d-2, &
			4.2d-2, 7.2d-2, 1.1d-1, 1.55d-1, 0.26d0, 0.46d0, 0.82d0, 1.2d0, 1.57d0, 1.96d0, 2.36d0, 2.76d0, 3.17d0, 3.59d0, 4.01d0, 4.43d0, &
			4.86d0, 5.29d0, 5.72d0, 6.15d0, 6.59d0, 7.03d0, 7.46d0, 11.89d0, 16.37d0, &
			25.39d0, 34.45d0, 52.6d0, 70.8d0, 88.9d0, 111.7d0, 134.4d0, 157.2d0, 180.0d0, 271.0d0, 453.0d0, 908.0d0, 4547.0d0, 9100.0d0, 45530.0d0/)


    loc = 45
	do i=1,45
        if (chi_g <= chi_g_array(i+1) .and. chi_g >= chi_g_array(i)) then
			loc = i
			exit
		endif
	enddo	
	if (loc < 6) then
		coef = (log(peak_array(loc+1)) - log(peak_array(loc)))/(chi_g_array(loc+1) - chi_g_array(loc))
		peak_pairpkt = exp(log(peak_array(loc)) + coef*(chi_g - chi_g_array(loc)))
	else	
		coef = (peak_array(loc+1) - peak_array(loc))/(chi_g_array(loc+1) - chi_g_array(loc))
		peak_pairpkt = peak_array(loc) + coef*(chi_g - chi_g_array(loc))
	endif


end function peak_pairpkt

!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!  function interpolation of integration of the pkt pair prob rate  
!-------------------------------------------------------------------------------

function interp_spec_pairpkt(chi_g)

	implicit none

integer :: n
real(p_double), intent(in) :: chi_g
real(p_double) :: interp_spec_pairpkt
! integer :: i, loc
! real(p_double) :: coef
real(p_double) :: ri1,rk1,rip1,rkp1, p, ltc

real(p_double), dimension(15), parameter :: p_array = (/   -4.683083764792118 ,  -1.091211988944788 , -4.309558261657857 ,  &
		4.491779676559683,   -1.107281703238826,   -0.461518946888616 ,  1.754341491184969 , -2.973628079492484 ,   4.567415354988209 , &
		 -1.114863165044028,   -0.019692135577373 ,  0.199292591618420 , -0.774418835297167  , 3.062688187056846 , -0.696902875956990 /)



!real(p_double), dimension(40), parameter :: chi_g_array = (/ 0.1d0, 0.2d0, 0.3d0, 0.4d0, 0.5d0, 0.6d0, &
!				0.7d0, 0.8d0, 0.9d0, 1.0d0, 1.2d0, 1.5d0, 2.0d0, 2.5d0, 3.0d0, 3.5d0, 4.0d0, 4.5d0, 5.0d0, &
!				5.5d0, 6.0d0, 6.5d0, 7.0d0, 7.5d0, 8.0d0, 8.5d0, 9.0d0, 9.5d0, 10.0d0, 15.0d0, 20.0d0, &
!				30.0d0, 60.0d0, 80.0d0, 100.0d0, 125.0d0, 150.0d0, 175.0d0, 200.0d0, 300.0d0/)
!real(p_double), dimension(40), parameter :: int_eval_array = (/ 8.6d12, 6.7d6, 5.0d4, 4.0d3, 821.0d0, 275.0d0, & 
!				125d0, 66d0, 40.9d0, 26.4d0, 13.8d0, 6.94d0, 3.25d0, 1.96d0, 1.35d0, 1.02d0, 0.8d0, 0.66d0, 0.56d0, &
!				0.49d0, 0.43d0, 0.38d0, 0.35d0, 0.32d0, 0.293d0, 0.272d0, 0.253d0, 0.237d0, 0.22d0, 0.14d0, 0.10d0, &
!				0.072d0, 0.04d0, 0.0314d0, 0.0264d0, 0.022d0, 0.019d0, 0.0173d0, 0.0156d0, 0.011d0/)
!
!	do i=1,39
!		if (chi_g.le.chi_g_array(i+1) .and. chi_g.ge.chi_g_array(i)) then
!			loc = i
!			exit
!		endif
!	enddo
!	if (loc < 6) then
!		coef = (log(int_eval_array(loc+1))-log(int_eval_array(loc)))/(chi_g_array(loc+1)-chi_g_array(loc))
!		interp_spec_pairpkt = exp( log(int_eval_array(loc)) + coef*(chi_g-chi_g_array(loc)))
!	else	
!		coef = (int_eval_array(loc+1)-int_eval_array(loc))/(chi_g_array(loc+1)-chi_g_array(loc))
!		interp_spec_pairpkt = int_eval_array(loc) + coef*(chi_g-chi_g_array(loc))
!	endif
	
 if ( chi_g > 1000.0 ) 	then
 	call bessik(4.0d0/(3.0d0*chi_g),1.d0/3.d0,ri1,rk1,rip1,rkp1)
	! function given by Erber, Review of Modern physics 1966
	interp_spec_pairpkt = sqrt(3.0d0)*M_PI*0.16d0*(rk1**2); 
	
	else
		ltc = log10(chi_g)
		if ( chi_g < 1.0 ) then
			n = 1 
			else if (chi_g < 10.0) then
					n = 6
				 else 
					n = 11
						
		endif
		p = p_array(n)*(ltc**4) + p_array(n+1)*(ltc**3) + p_array(n+2)*(ltc**2) + p_array(n+3)*ltc + p_array(n+4) 
		! polynomial expansion of the total emission rate integral - more precise than Erber's approximate exp.	
		interp_spec_pairpkt = ( 10**p ) / chi_g
    endif
    		


end function interp_spec_pairpkt

!-------------------------------------------------------------------------------




end module m_pkt


