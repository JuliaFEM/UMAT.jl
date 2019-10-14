! This file is a part of JuliaFEM.
! License is MIT: see https://github.com/JuliaFEM/UMAT.jl/blob/master/LICENSE
SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,&
	RPL,DDSDDT,DRPLDE,DRPLDT,&
	STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
	NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,&
	CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

implicit none
! Rate-independent Gurson-Tveergaad-Needleman porous plasticity model
!INCLUDE'ABA_PARAM.INC'

CHARACTER*80 CMNAME
REAL*8,INTENT(INOUT) :: STRESS(NTENS),STATEV(NSTATV),&
DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),&
STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),&
PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

REAL*8,INTENT(INOUT) :: SSE,SPD,SCD,RPL,DRPLDT,DTIME,TEMP,DTEMP,&
	PNEWDT,CELENT
integer,intent(inout) :: NDI,NSHR,NOEL,NPT,LAYER,NTENS,NSTATV,&
	KSPT,KSTEP,KINC,NPROPS

INTEGER :: I,J,MAXITERS,DISP
REAL*8 :: tol,res
REAL*8 :: f0,fc,ff,fn,sn,en,qinf,b,C,ga,EE,nu,syield
REAL*8 :: K,G,la
REAL*8 :: q1,q2,q3,fu
REAL*8 :: pi,A
REAL*8 :: DEp,DEq,ep,f
REAL*8 :: ft,dftdf,epn,fnn
REAL*8 :: p,q,sflow,phi,cosht,sinht,dphidq,dphidp
REAL*8,DIMENSION(6) :: STRIAL,X,EELAS,EPLAS,BETA,DEin,N,onevec,DSTRESS,BETAiter
REAL*8,DIMENSION(6,6) :: CELAS
REAL*8 :: eratio,qiter,piter,phiter,dpderatio,dqderatio
! Constants
pi=3.14159265358979
tol = 1e-6
MAXITERS = 20
IF (KSTEP==2) THEN
	DISP=0
ELSE
	DISP=0
END IF
IF (DISP>=1) THEN
	PRINT *, "***** NEW INCREMENT IN MAIN ROUTINE! STEP:", KSTEP, " *****"
	PRINT *, "TIME:", TIME(1), "DTIME:", DTIME
END IF
!Load properties
call LoadPROPS(f0,fc,ff,fn,sn,en,qinf,b,C,ga,EE,nu,syield,PROPS,NPROPS,DISP)

!Bulk and shear moduli
K=EE/(3.0*(1.0-nu))
G=EE/(2.0*(1.0+nu))
la=K-2.0/3.0*G

!Tvergaard values for periodic void distributions
q1=1.5
q2=1.0
q3=q1**2

!Effective porosity at failure
fu=(q1-SQRT(q1**2-q3))/q3

! Load state variables
call LoadSTATEV(f,ep,DEp,DEq,EELAS,EPLAS,X,STATEV,NSTATV,&
	NTENS,NDI,NSHR,DROT,DISP)

IF (STATEV(NSTATV)==0) THEN
	! Initialize porosity
	f = f0
	STATEV(NSTATV)=1
END IF
!Effective linear elasticity tensor CELAS
CELAS(1:NTENS,1:NTENS)=0.0
DO I=1,NDI
	DO J=1,NDI
		CELAS(I,J)=la
	END DO 
	CELAS(I,I)=2.0*G+la
END DO 
DO I=NDI+1,NTENS
	CELAS(I,I)=G
END DO 

! Trial stress
STRIAL=MATMUL(CELAS,EELAS+DSTRAN)
! Relative stress
BETA=STRIAL-X

! Calculate pressure and von mises
	call pressure(BETA,p)
	call vonmises(BETA,q)

! Flow stress
sflow = syield+qinf*(1.0-exp(-b*ep))

! Update effective porosity and its derivative
	call Updateft(f,ft,dftdf,fu,fc,ff)

! Yield function
	call mcosh(1.5*q2*p/sflow,cosht)
	call msinh(1.5*q2*p/sflow,sinht)
phi = (q/sflow)**2 + 2.0*ft*q1*cosht-(1.0+q3*ft**2)
IF (DISP>=2) THEN
	print *, "Phi:", phi
	print *, "Pressure:", p, "Von Mises:", q
END IF
IF (phi<0.0) THEN
	! Elastic step
	! Update state variables
	IF (DISP>=2) THEN
		print *, "EELAS(n):"
		print *, EELAS
		print *, "DSTRAN:"
		print *, DSTRAN
		print *, "STRESS(n):"
		print *, STRESS
		print *, "STRIAL:"
		print *, STRIAL
	END IF
	EELAS = EELAS+DSTRAN
	EPLAS = EPLAS
	X = X
	call UpdateSTATEV(f,ep,DEp,DEq,EELAS,EPLAS,X,STATEV,NSTATV,NTENS,DISP)
	! Update STRESS
	STRESS=STRIAL
	! Update DDSDDE
	DDSDDE = CELAS
ELSE
	IF (DISP>=2) THEN
		print *, "***** Yielding *****"
	END IF
	fnn = f
	! Plastic step
	! 22.11.2013 Added estimation for ratio of elasticity
	DSTRESS = MATMUL(CELAS,DSTRAN)
	onevec=(/1.0,1.0,1.0,0.0,0.0,0.0/)
	N = 1.5*(BETA+p*onevec)/q
	dpderatio = -1.0/3.0*(DSTRESS(1)+DSTRESS(2)+DSTRESS(3))
	dqderatio = DOT_PRODUCT(N(1:3),DSTRESS(1:3))+2*DOT_PRODUCT(N(4:6),DSTRESS(4:6))
	
	res = 1000.0
	! Initial guess: ratio of elasticity = 0.1
	eratio=0.1
	DO I=1,MAXITERS
		BETAiter=STRESS+(1.0-eratio)*DSTRESS-X
			call pressure(BETAiter,piter)
			call vonmises(BETAiter,qiter)
		! Yield function
		call mcosh(1.5*q2*piter/sflow,cosht)
		call msinh(1.5*q2*piter/sflow,sinht)
		phiter = (qiter/sflow)**2 + 2.0*ft*q1*cosht-(1.0+q3*ft**2)
		dphidq=2.0*qiter/sflow**2
		dphidp=3.0*q1*q2*ft/sflow*sinht
		eratio = eratio + phiter/(dphidp*dpderatio + dphidq*dqderatio)
		res = ABS(phiter)
		IF (res<tol) EXIT
	END DO
	IF ((res>tol) .OR. (eratio<0) .OR. (eratio>1) .OR. (eratio/=eratio)) THEN
		IF (DISP>=2) THEN
			print *, "Failed elasticity ratio:", eratio
		END IF
		eratio = 0.0
	END IF
	IF (DISP>=2) THEN
		print *, "Elasticity ratio:", eratio
	END IF
	
	res = 1000.0
	! 12.11.2013 Changed plastic starting-guess
	! Close to ideal-plastic starting guess
	! Least-squares method, N:N = 3/2
	DEp = (1.0-eratio)*(DSTRAN(1)+DSTRAN(2)+DSTRAN(3))
	DEq = 2.0/3.0*DOT_PRODUCT(N,(1.0-eratio)*DSTRAN-DEp/3.0*onevec)
	!DEq = 2.0/3.0*(DOT_PRODUCT(N(1:3),(1.0-eratio)*DSTRAN(1:3)-DEp/3.0)+2*DOT_PRODUCT(N(4:6),(1.0-eratio)*DSTRAN(4:6)))
	IF ((DEq==0) .and. (DEp==0)) THEN
		DEq = tol
	END IF
	
	! 19.11.2013 Extended starting guess to f and ep as well
	epn = ep
	fnn = f
	ep = ep + (-p*DEp + q*DEq)/((1.0-f)*sflow)
	IF (p<0) THEN
		A=fn*exp(-0.5*(-en + ep)**2/sn**2)/(sqrt(2.0*pi)*sn)
		f = f + (1-f)*DEp + A*(ep-epn)
	ELSE
		f = f + (1-f)*DEp
	END IF
	
	! Returns iteration variables: f, ep, DEp, DEq and residual res
	IF (DISP>=2) THEN
		print *, "EELAS(n):"
		print *, EELAS
		print *, "DSTRAN:"
		print *, DSTRAN
		print *, "STRESS(n):"
		print *, STRESS
		print *, "STRIAL:"
		print *, STRIAL
		print *, "X(n):"
		print *, X
		print *, "f(n), ep(n)"
		print *, fnn, epn
		print *, "DEp(n), DEq(n)"
		print *, DEp, DEq
	END IF
	! Analytical tangent added 19.11.2013, Calculated inside PlasticIteration
	call PlasticIteration(f,ep,DEp,DEq,DDSDDE,STRIAL,X,epn,fnn,MAXITERS,TOL,PROPS,&
		NPROPS,NTENS,NDI,res,DISP)
	IF (res>tol) THEN
		PRINT *, "Plastic iteration failed at",TIME(1),"in element", NOEL, " in point", NPT
		PNEWDT=0.5
	END IF
	IF (DISP>=2) THEN
		print *, "f(n+1), ep(n+1)"
		print *, f, ep
		print *, "DEp(n+1), DEq(n+1)"
		print *, DEp, DEq
	END IF
		
	DEin = EPLAS
	! Update STRESS,X,DEin
	call UpdateStress(f,ep,DEp,DEq,STRIAL,X,PROPS,NPROPS,STRESS,DEin,fnn,NTENS,NDI,DISP)
	IF (DISP>=2) THEN
		print *, "STRESS(n+1):"
		print *, STRESS
		print *, "X(n+1):"
		print *, X
	END IF
	
	! Update internal strains using return mapping algorithms
	EELAS = EELAS+DSTRAN-DEin
	EPLAS = EPLAS+DEin
	IF (DISP>=2) THEN
		print *, "EELAS(n+1):"
		print *, EELAS
	END IF
	
	! Update state variables
	call UpdateSTATEV(f,ep,DEp,DEq,EELAS,EPLAS,X,STATEV,NSTATV,NTENS,DISP)
	
	! Disabled 19.11.2013 for the Analytical consistent tangent to be enabled
	! Calculate Numerical consistent tangent DDSDDE
	! using Forward-difference method
	!DISP=5
	!call NumericalConsistentTangent(STRESS,STATEV,DSTRAN,DDSDDE,NDI,&
	!	NSHR,NTENS,NSTATV,PROPS,NPROPS,DROT,DISP)
	IF (DISP>=2) THEN
		print *, "DDSDDE:"
		print *, DDSDDE
	END IF
	DISP=0
END IF
END SUBROUTINE
subroutine MVector(S)
	implicit none
	real(8),dimension(6), intent(inout) :: S
	S(4:6)=S(4:6)*sqrt(2.0)
end subroutine
subroutine MiVector(S)
	implicit none
	real(8),dimension(6), intent(inout) :: S
	S(4:6)=S(4:6)/sqrt(2.0)
end subroutine
subroutine MMatrix(S)
	implicit none
	real(8),dimension(6,6), intent(inout) :: S
	S(4:6,1:3) = S(4:6,1:3)*sqrt(2.0)
	S(1:3,4:6) = S(1:3,4:6)*sqrt(2.0)
	S(4:6,4:6) = S(4:6,4:6)*2.0
end subroutine
subroutine MiMatrix(S)
	implicit none
	real(8),dimension(6,6), intent(inout) :: S
	S(4:6,1:3) = S(4:6,1:3)/sqrt(2.0)
	S(1:3,4:6) = S(1:3,4:6)/sqrt(2.0)
	S(4:6,4:6) = S(4:6,4:6)/2.0
end subroutine
subroutine UpdateStress(f,ep,DEp,DEq,STRIAL,X,PROPS,&
	NPROPS,STRESS,DEin,fnn,NTENS,NDI,DISP)
	implicit none
	
	real(8),intent(in) :: DEp,DEq,ep,f
	real(8),dimension(6),intent(in) :: STRIAL
	real(8),dimension(6),intent(inout) :: STRESS,X,DEin
	integer,intent(in) :: NPROPS
	real(8),dimension(NPROPS),intent(in) :: PROPS
	real(8),intent(in) :: fnn
	integer,intent(in) :: NTENS,NDI,DISP
	real(8),dimension(6) :: onevec,N,BETAplus
	real(8),dimension(6,6) :: CELAS
	integer :: I, J
	real(8) :: f0,fc,ff,fn,sn,en,qinf,b,C,ga,EE,nu,syield,K,G,la
	real(8) :: q1,q2,q3,fu
	real(8) :: T,s,Ds,Deinscl,q,p,qplus,pplus
	real(8) :: ft,dftdf,ftn,dftndf
	IF (DISP>=1) THEN
		print *, "***** Updating STRESS *****"
	END IF
	call LoadPROPS(f0,fc,ff,fn,sn,en,qinf,b,C,ga,EE,nu,syield,PROPS,NPROPS,DISP)
	!Bulk and shear moduli
	K=EE/(3.0*(1.0-nu))
	G=EE/(2.0*(1.0+nu))
	!la=K-2.0/3.0*G
	
	!Tvergaard values for periodic void distributions
	q1=1.5
	q2=1.0
	q3=q1**2

	! Elasticity tensor CELAS
	!CELAS(1:NTENS,1:NTENS)=0.0
	!DO I=1,NDI
	!	DO J=1,NDI
	!		CELAS(I,J)=la
	!	END DO 
	!	CELAS(I,I)=2.0*G+la
	!END DO 
	!DO I=NDI+1,NTENS
	!	CELAS(I,I)=G
	!END DO
	
	!Effective porosity at failure
	fu=(q1-SQRT(q1**2-q3))/q3
	
	call Updateft(fnn,ftn,dftndf,fu,fc,ff)
	call Updateft(f,ft,dftdf,fu,fc,ff)
	
	s=1.0-ft/fu
	Ds=(ftn-ft)/fu
	Deinscl=SQRT(2.0/9.0*DEp**2+DEq**2)
	T=1.0-Ds/s+s*ga*Deinscl
	
	BETAplus = STRIAL-X/T
	! Calculate p and q
	call pressure(BETAplus,pplus)
	call vonmises(BETAplus,qplus)
	
	p = pplus+(K+2.0/9.0*C/T*s)*DEp
	q = qplus-(G+1.0/3.0*C/T*s)*3.0*DEq
	
	onevec = (/1.0,1.0,1.0,0.0,0.0,0.0/)
	
	!N=(BETAplus+pplus*onevec)/(2.0/3.0*q+(G+1.0/3.0*C/T*s)*2.0*DEq)
	! Changed 12.11.2013 - more compact form
	N = 1.5*(BETAplus+pplus*onevec)/(q+(3.0*G+C/T*s)*DEq)
	DEin=1.0/3.0*DEp*onevec+DEq*N
	DO I=NDI+1,NTENS
		DEin(I) = DEin(I)*2.0
	END DO
	X = (2.0/3.0*C*DEin*s+X)/T
	STRESS = STRIAL-K*DEp*onevec-2.0*G*DEq*N
	!STRESS = STRIAL-MATMUL(CELAS,DEin)
	IF (DISP>=1) THEN
		print *, "***** STRESS updated *****"
	END IF
end subroutine
subroutine PlasticIteration(f,ep,DEp,DEq,DDSDDE,STRIAL,X,epn,fnn,MAXITERS,TOL,&
	PROPS,NPROPS,NTENS,NDI,res,DISP)
	implicit none
	real(8),intent(inout) :: DEp,DEq,ep,f
	real(8),dimension(6,6),intent(inout) :: DDSDDE
	real(8),dimension(6),intent(in) :: STRIAL,X
	integer,intent(in) :: NPROPS,MAXITERS,DISP,NTENS,NDI
	real(8),dimension(NPROPS),intent(in) :: PROPS
	real(8),intent(in) :: TOL,epn,fnn
	real(8),intent(inout) :: res
	real(8), dimension(4,4) :: ALIN
	real(8), dimension(4) :: BLIN, dXLIN
	real(8), dimension(6) :: BETAplus
	real(8) :: f0,fc,ff,fn,sn,en,qinf,b,C,ga,EE,nu,syield,K,G,la
	real(8) :: q1,q2,q3,fu,pi
	real(8) :: T,s,Ds,Deinscl,q,p,Depscl
	real(8) :: ft,dftdf,ftn,dftndf,Df
	real(8) :: sinht,cosht
	real(8) :: pplus,qplus
	real(8) :: sflow, dsflowdep, A, dAdep
	real(8) :: dTdDEp, dTdDEq, dTdep, dTdf
	real(8) :: dpdDEp, dpdDEq, dpdep,dpdf
	real(8) :: dqdDEp, dqdDEq, dqdep,dqdf
	real(8) :: mult1,mult2,mult3,mult4
	real(8) :: dphidq, dphidq2, dphidp, dphidpdf, dphidp2
	real(8) :: dphidqdsflow,dphidpdsflow,dphidsflow
	real(8) :: F1, dF1dDEp, dF1dDEq, dF1dep, dF1df
	real(8) :: F2, dF2dDEp, dF2dDEq, dF2dep, dF2df
	real(8) :: G1, dG1dDEp, dG1dDEq, dG1dep, dG1df
	real(8) :: G2, dG2dDEp, dG2dDEq, dG2dep, dG2df
	integer :: I,J
	real(8), dimension(2,2) :: Amat2,Bmat2,Cmat2
	real(8), dimension(2) :: bvec2
	real(8), dimension(6) :: N, onevec, dG1dstr, dG2dstr, dF1dstr, dF2dstr
	real(8), dimension(6) :: dpdstr, dqdstr, dH1dstr, dH2dstr
	real(8), dimension(6) :: dNdT
	real(8), dimension(6) :: B1,B2,Bp,Bq,Bep,Bf
	real(8), dimension(6,6) :: dNdstr,PP,tempmat6,QQ,tempmat7,CELAS
	real(8), dimension(3,3,3,3) :: Itensor, Ptensor,Qtensor,Ctensor,Dtensor,tensor1,tensor2
	real(8) :: qe,dH1dDEp,dH1dDEq,dH2dDEp,dH2dDEq
	real(8) :: resn
	
	IF (DISP>=1) THEN
		print *, "***** Inside PlasticIteration *****"
	END IF
	call LoadPROPS(f0,fc,ff,fn,sn,en,qinf,b,C,ga,EE,nu,syield,PROPS,NPROPS,DISP)
	pi=3.14159265358979
	
	!Bulk and shear moduli
	K=EE/(3.0*(1.0-nu))
	G=EE/(2.0*(1.0+nu))
	la=K-2.0/3.0*G
	!Effective linear elasticity tensor CELAS
	CELAS(1:NTENS,1:NTENS)=0.0
	DO I=1,NDI
		DO J=1,NDI
			CELAS(I,J)=la
		END DO 
		CELAS(I,I)=2.0*G+la
	END DO 
	DO I=NDI+1,NTENS
		CELAS(I,I)=G
	END DO 
		
	!Tvergaard values for periodic void distributions
	q1=1.5
	q2=1.0
	q3=q1**2

	!Effective porosity at failure
	fu=(q1-SQRT(q1**2-q3))/q3
	
	call Updateft(f,ftn,dftndf,fu,fc,ff)
	if (res<tol) THEN
		res=1.5*tol
	END IF
	IF (DISP>=1) THEN
		print *, "***** Starting PlasticIteration *****"
	END IF
	IF (DISP>=3) THEN
		print *, "f, ep"
		print *, f, ep
		print *, "DEp, DEq"
		print *, DEp, DEq
		print *, "STRIAL:"
		print *, STRIAL
		print *, "X(n):"
		print *, X
	END IF
	! Plastic iterations
	do I=1,MAXITERS
		IF (DISP>=3) THEN
			print *, "***** Iteration", I, " *****"
			print *, "f, ep"
			print *, f, ep
			print *, "DEp, DEq"
			print *, DEp, DEq
		END IF
		! 22.11.2013 -  Added divergence control
		resn=res
		! Update ft and dftdf
		call Updateft(f,ft,dftdf,fu,fc,ff)
		s=1.0-ft/fu
		Ds=(ftn-ft)/fu
		Deinscl=SQRT(2.0/9.0*DEp**2+DEq**2)
		!print *, "s:", s, "Ds:", Ds
		T=1.0-Ds/s+s*ga*Deinscl
		!print *, "Deinscl:", Deinscl, "dftdf:", dftdf
		BETAplus = STRIAL-X/T
		sflow=syield+qinf*(1.0-exp(-b*ep))
		dsflowdep=b*qinf*exp(-b*ep)
		
		!print *, "Sflow:", sflow, "T:", T
		
		A=fn*exp(-0.5*(-en + ep)**2/sn**2)/(sqrt(2.0*pi)*sn)
		dAdep=-(ep-en)/sn**2*A
		
		!print *, "A:", A
		
		! Calculate p and q
		call pressure(BETAplus,pplus)
		call vonmises(BETAplus,qplus)
		
		p = pplus+(K+2.0/9.0*C/T*s)*DEp
		q = qplus-(G+1.0/3.0*C/T*s)*3.0*DEq
		IF (DISP>=3) THEN
			print *, "Pressure:", p, "Von Mises:", q
		END IF
		call msinh(1.5*q2*p/sflow,sinht)
		call mcosh(1.5*q2*p/sflow,cosht)
		!print *, "sinh*:", sinht, "cosh*:", cosht
		
		! Derivatives
		dTdDEp = 2.0*ga*s*DEp/(9.0*Deinscl)
		dTdDEq = ga*s*DEq/Deinscl
		dTdep = 0.0
		dTdf = dftdf*(1.0/(fu-ft)-(ftn-ft)/(fu-ft)**2-ga*Deinscl/fu) ! Different in Metzger 2009
		
		dpdDEp = K-(X(1)+X(2)+X(3))/(3.0*T**2)*dTdDEp-2.0*C*DEp*s/(9.0*T**2)*dTdDEp+2.0*C*s/(9.0*T)
		dpdDEq = -(X(1)+X(2)+X(3))/(3.0*T**2)*dTdDEq-2.0*C*DEp*s/(9.0*T**2)*dTdDEq
		dpdep = -(X(1)+X(2)+X(3))/(3.0*T**2)*dTdep-2.0*C*DEp*s/(9.0*T**2)*dTdep
		dpdf = -(X(1)+X(2)+X(3))/(3.0*T**2)*dTdf-2.0*C*DEp*s/(9.0*T**2)*dTdf-2.0*C*DEp/(9.0*T*fu)*dftdf
		
		!mult1 = STRIAL(1)**2+STRIAL(2)**2+STRIAL(3)**2+2.0*(STRIAL(4)**2+STRIAL(5)**2+STRIAL(6)**2)
		!mult2 = 2.0*(STRIAL(1)*X(1)+STRIAL(2)*X(2)+STRIAL(3)*X(3)+2.0*(STRIAL(4)*X(4)+STRIAL(5)*X(5)+STRIAL(6)*X(6)))
		!mult3 = X(1)**2+X(2)**2+X(3)**2+2.0*(X(4)**2+X(5)**2+X(6)**2)
		!mult4 = (1.5*mult2/T**2-3.0*mult3/T**3)/SQRT(6.0*(mult1-mult2/T+mult3/T**2))

		mult4 = (0.25*(2*X(1)/T**2 - 2*X(2)/T**2)*(STRIAL(1) - STRIAL(2) - X(1)/T + X(2)/T)&
			+ 0.25*(2*X(1)/T**2 - 2*X(3)/T**2)*(STRIAL(1) - STRIAL(3) - X(1)/T + X(3)/T)&
			+ 0.25*(2*X(2)/T**2 - 2*X(3)/T**2)*(STRIAL(2) - STRIAL(3) - X(2)/T + X(3)/T)&
			+ 3.0*X(4)*(STRIAL(4) - X(4)/T)/T**2 + 3.0*X(5)*(STRIAL(5) - X(5)/T)/T**2&
			+ 3.0*X(6)*(STRIAL(6) - X(6)/T)/T**2)/sqrt(3.0*(STRIAL(4) - X(4)/T)**2&
			+ 3.0*(STRIAL(5) - X(5)/T)**2 + 3.0*(STRIAL(6) - X(6)/T)**2&
			+ 0.5*(STRIAL(1) - STRIAL(2) - X(1)/T + X(2)/T)**2&
			+ 0.5*(STRIAL(1) - STRIAL(3) - X(1)/T + X(3)/T)**2&
			+ 0.5*(STRIAL(2) - STRIAL(3) - X(2)/T + X(3)/T)**2)
		
		dqdDEp = mult4*dTdDEp&
			+C*s*DEq/T**2*dTdDEp
		dqdDEq = mult4*dTdDEq&
			+C*s*DEq/T**2*dTdDEq-3.0*G-C*s/T
		dqdep = mult4*dTdep&
			+C*s*DEq/T**2*dTdep
		dqdf = mult4*dTdf&
			+C*DEq*(s/T**2*dTdf+dftdf/(T*fu))
		
		! Equation 1
		dphidq=2.0*q/sflow**2
		dphidq2=2.0/sflow**2
		dphidp=3.0*q1*q2*ft/sflow*sinht
		dphidpdf=3.0*q1*q2*dftdf/sflow*sinht
		dphidp2=4.5*q1*q2**2*ft/sflow**2*cosht
		dphidqdsflow = -4.0*q/sflow**3
		dphidpdsflow = -4.5*ft*p*q1*q2**2*cosht/sflow**3&
			-3.0*ft*q1*q2*sinht/sflow**2
		F1 = DEp*dphidq+DEq*dphidp
		
		dF1dDEp = dphidq+DEp*dphidq2*dqdDEp+DEq*dphidp2*dpdDEp
		dF1dDEq = DEp*dphidq2*dqdDEq+dphidp+DEq*dphidp2*dpdDEq
		dF1dep = DEp*(dphidq2*dqdep+dphidqdsflow*dsflowdep)&
			+ DEq*(dphidp2*dpdep+dphidpdsflow*dsflowdep)
		dF1df = DEp*(dphidq2*dqdf) + DEq*(dphidp2*dpdf+dphidpdf)
		
		! Equation 2
		F2 = (q/sflow)**2+2.0*ft*q1*cosht-(1.0+q3*ft**2)
		dphidsflow = -3.0*ft*p*q1*q2*sinht/sflow**2 - 2.0*q**2/sflow**3
		dF2dDEp = dphidq*dqdDEp + dphidp*dpdDEp
		dF2dDEq = dphidq*dqdDEq + dphidp*dpdDEq
		dF2dep = dphidq*dqdep + dphidp*dpdep + dphidsflow*dsflowdep
		dF2df = dphidq*dqdf + 2.0*dftdf*q1*cosht + dphidp*dpdf-2.0*q3*dftdf*ft
		
		! Equation 3
		Depscl = (-p*DEp+q*DEq)/((1.0-f)*sflow)
		G1 = (ep-epn)-Depscl
		dG1dDEp = -(-dpdDEp*DEp-p+dqdDEp*DEq)/((1.0-f)*sflow)
		dG1dDEq = -(-dpdDEq*DEp+dqdDEq*DEq+q)/((1.0-f)*sflow)
		dG1dep = 1.0 - (-dpdep*DEp+dqdep*DEq)/((1.0-f)*sflow) + (-p*DEp+q*DEq)/((1.0-f)*sflow**2)*dsflowdep
		dG1df = -(-dpdf*DEp+dqdf*DEq)/((1.0-f)*sflow) - (-p*DEp+q*DEq)/((1.0-f)**2*sflow)
		
		! Equation 4
		IF (p<0) THEN
			! Nucleation term is only active under hydrostatic tension p<0
			! Condition added 11.11.2013
			Df = (1.0-f)*DEp+A*Depscl
			!Df = (1.0-f)*DEp+A*(ep-epn)
			G2 = (f-fnn)-Df
			!dG2dDEp = f-1.0
			!dG2dDEq = 0.0
			!dG2dep = -dAdep*(ep-epn)-A
			!dG2df = 1.0+DEp
			dG2dDEp = -(1.0-f)+A*dG1dDEp
			dG2dDEq = A*dG1dDEq
			dG2dep = -dAdep*Depscl + A*(dG1dep-1.0)
			dG2df = 1.0 + DEp + A*dG1df
		ELSE
			! Only growth is active under hydrostatic compression p>0
			! Condition added 11.11.2013
			Df = (1.0-f)*DEp
			G2 = (f-fnn)-Df
			dG2dDEp = f-1.0
			dG2dDEq = 0.0
			dG2dep = 0.0
			dG2df = 1.0+DEp
		END IF
		
		!Linearized system for Newton's iterations
		ALIN(1,1)=dF1dDEp
		ALIN(2,1)=dF1dDEq
		ALIN(3,1)=dF1dep
		ALIN(4,1)=dF1df
		ALIN(1,2)=dF2dDEp
		ALIN(2,2)=dF2dDEq
		ALIN(3,2)=dF2dep
		ALIN(4,2)=dF2df
		ALIN(1,3)=dG1dDEp
		ALIN(2,3)=dG1dDEq
		ALIN(3,3)=dG1dep
		ALIN(4,3)=dG1df
		ALIN(1,4)=dG2dDEp
		ALIN(2,4)=dG2dDEq
		ALIN(3,4)=dG2dep
		ALIN(4,4)=dG2df
		
		BLIN=-1.0*(/F1,F2,G1,G2/)
		
		!Convergence check
		res=0.0
		do J=1,4
			if (abs(BLIN(J))>res) then
				res=abs(BLIN(J))
			end if
		end do
		IF (DISP>=3) THEN
			print *, "Residual:", res
		END IF
		
		
		! Solve linear system
		IF (DISP>=3) THEN
			print *, "ALIN"
			print *, ALIN
			print *, "BLIN"
			print *, BLIN
		END IF
		call LinearSolve(TRANSPOSE(ALIN),4,BLIN,dXLIN)
		
		IF (DISP>=3) THEN
			print *, "dXLIN"
			print *, dXLIN
		END IF
		
		!Update variables
		DEp=DEp+dXLIN(1)
		DEq=DEq+dXLIN(2)
		ep=ep+dXLIN(3)
		f=f+dXLIN(4)
		! 22.11.2013 - Divergence control
		IF ((res<=TOL) .OR. (res>resn)) EXIT
	end do
	! 18.11.2013 Added Analytical Consistent Tangent
	BETAplus = STRIAL-X/T
	! Calculate p and q
	call pressure(BETAplus,pplus)
	call vonmises(BETAplus,qplus)
	
	p = pplus+(K+2.0/9.0*C/T*s)*DEp
	q = qplus-(G+1.0/3.0*C/T*s)*3.0*DEq
	qe = q + 3.0*G*DEq
	
	onevec = (/1.0,1.0,1.0,0.0,0.0,0.0/)
	tempmat6(1:6,1:6) = 0.0
	tempmat7(1:6,1:6) = 0.0
	N = 1.5*(BETAplus+pplus*onevec)/qplus
	
	! Determine derivatives d/dstr
	dpdstr = -1.0/3.0*onevec
	dqdstr = N
	dF1dstr = dphidp2*dpdstr*DEq + dphidq2*dqdstr*DEp
	dF2dstr = dphidp*dpdstr + dphidq*dqdstr
	dG1dstr = -(-dpdstr*DEp+dqdstr*DEq)/((1.0-f)*sflow)
	IF (p<0) THEN
		dG2dstr = A*dG1dstr
	ELSE
		dG2dstr = 0.0*onevec
	END IF
	! Total differentiation of equations 3 and 4
	! Solve variations varep, varf using static condensation from 2x2 system
	Amat2 = ALIN(3:4,3:4)
	Bmat2 = -1.0*ALIN(1:2,3:4)
	call MatrixInverse2(Amat2, Cmat2)
	Amat2 = MATMUL(Cmat2,Bmat2)
	dH1dDEp = Amat2(1,1)
	dH1dDEq = Amat2(2,1)
	dH1dstr = -(Cmat2(1,1)*dG1dstr+Cmat2(2,1)*dG2dstr)
	dH2dDEp = Amat2(1,2)
	dH2dDEq = Amat2(2,2)
	dH2dstr = -(Cmat2(1,2)*dG1dstr+Cmat2(2,2)*dG2dstr)
	! Variations of internal variables are expressed using varDEp, varDEq and varstr:
	! varep = dH1dDEp*varDEp + dH1dDEq*varDEq + DOT_PRODUCT(dH1dstr,varstr)
	! varf = dH2dDEp*varDEp + dH2dDEq*varDEq + DOT_PRODUCT(dH2dstr,varstr)
	
	! Total differentiation of equations 1 and 2 with respect to DEp, DEq, ep, f, STRIAL
	! Expressing varep and varf using the above relations => varDEp, varDEq, varstr
	Amat2 = ALIN(1:2,1:2)
	Amat2(1,1) = Amat2(1,1) + dF1dep*dH1dDEp + dF1df*dH2dDEp
	Amat2(2,1) = Amat2(2,1) + dF1dep*dH1dDEq + dF1df*dH2dDEq
	Amat2(1,2) = Amat2(1,2) + dF2dep*dH1dDEp + dF2df*dH2dDEp
	Amat2(2,2) = Amat2(2,2) + dF2dep*dH1dDEq + dF2df*dH2dDEq
	B1 = -1.0*(dF1dstr + dF1dep*dH1dstr + dF1df*dH2dstr)
	B2 = -1.0*(dF2dstr + dF2dep*dH1dstr + dF2df*dH2dstr)
	call MatrixInverse2(Amat2,Cmat2)
	! Variation mapping operators
	! dDEp = Bp:dSTRIAL
	! dDEq = Bq:dSTRIAL
	! dep = Bep:dSTRIAL
	! df = Bf:dSTRIAL
	Bp = Cmat2(1,1)*B1 + Cmat2(2,1)*B2
	Bq = Cmat2(1,2)*B1 + Cmat2(2,2)*B2
	Bep = dH1dDEp*Bp + dH1dDEq*Bq + dH1dstr
	Bf = dH2dDEp*Bp + dH2dDEq*Bq + dH2dstr
	
	dNdT = -1.5*X/(qplus*T**2)-N/qplus*mult4
	
	!call OUTER_PRODUCT(N,N,tempmat6,6)
	!PP(1:6,1:6) = 0.0
	!PP(1:3,1) = (/2.0/3.0, -1.0/3.0, -1.0/3.0/)
	!PP(1:3,2) = (/-1.0/3.0, 2.0/3.0, -1.0/3.0/)
	!PP(1:3,3) = (/-1.0/3.0, -1.0/3.0, 2.0/3.0/)
	!PP(4,4) = 1.0
	!PP(5,5) = 1.0
	!PP(6,6) = 1.0
	!dNdstr = (1.5*PP-tempmat6)/qe
	!call OUTER_PRODUCT(onevec,Bp,tempmat6,6)
	!QQ = 1.0/3.0*tempmat6
	!call OUTER_PRODUCT(N,Bq,tempmat6,6)
	!QQ = QQ + tempmat6
	!! Build Dn/Dstr
	!! dNdT*dTdDEp*varDEp
	!call OUTER_PRODUCT(dNdT,Bp,tempmat7,6)
	!tempmat6 = dTdDEp*tempmat7
	!! dNdT*dTdDEq*varDEq
	!call OUTER_PRODUCT(dNdT,Bq,tempmat7,6)
	!tempmat6 = tempmat6 + dTdDEq*tempmat7
	!! dNdT*dTdep*varep
	!call OUTER_PRODUCT(dNdT,Bep,tempmat7,6)
	!tempmat6 = tempmat6 + dTdep*tempmat7
	!! dNdT*dTdf*varf
	!call OUTER_PRODUCT(dNdT,Bf,tempmat7,6)
	!tempmat6 = tempmat6 + dTdf*tempmat7
	!! Add dNdstr
	!tempmat6 = tempmat6 + dNdstr
	!! Combine to Q matrix
	!QQ = QQ + tempmat6*DEq
	!DDSDDE = CELAS - MATMUL(TRANSPOSE(CELAS),MATMUL(QQ,CELAS))
	
	! 20.11.2013 Added tensor calculus here
	
	! Fourth order identity tensor "Itensor"
	call IdentityTensor(Itensor)
	! Projection tensor "Ptensor"
	call TensorDyad(onevec,onevec,tensor1)
	Ptensor = Itensor-1.0/3.0*tensor1
	! Elasticity tensor "Ctensor"
	CTensor = K*tensor1+2.0*G*Ptensor
	
	! Calculate variation of N
	! dN/dstr
	call TensorDyad(N,N,tensor1)
	tensor2 = (1.5*Ptensor-tensor1)/qe
	call TensorDyad(dNdT,dTdDEp*Bp+dTdDEq*Bq+dTdep*Bep+dTdf*Bf,tensor1)
	! variation of var(N)/var(str)
	tensor2 = tensor2 + tensor1
	
	! Build Q tensor
	Qtensor = tensor2*DEq
	call TensorDyad(N,Bq,tensor1)
	Qtensor = Qtensor + tensor1
	call TensorDyad(onevec,Bp,tensor1)
	Qtensor = Qtensor + 1.0/3.0*tensor1
	
	! Calculate C_ijst Q_stpq C_pqkl
	call TensorQuad(Ctensor,Qtensor,Ctensor,tensor1)
	DTensor = CTensor - tensor1
	
	! Return 4th order tensor to 6x6 matrix
	call reduce4tensor(DTensor, DDSDDE)
	
	IF (DISP>=1) THEN
		print *, '***** Finished PlasticIteration *****'
	END IF
end subroutine
subroutine LoadPROPS(f0,fc,ff,fn,sn,en,qinf,b,C,ga,EE,nu,syield,PROPS,NPROPS,DISP)
	implicit none
	real(8),intent(inout) :: f0,fc,ff,fn,sn,en,qinf,b,C,ga,EE,nu,syield
	integer, intent(in) :: NPROPS,DISP
	real(8),dimension(NPROPS), intent(in) :: PROPS
	IF (DISP>=1) THEN
		print *, '***** Loading PROPS *****'
	END IF
	!Properties for porosity
	f0=PROPS(1)!f0,initial porosity[GJS-700:0.02,GJV-450:0.05,GJL-250:0.125]
	fc=PROPS(2)!fc,critical porosity[0.0325,0.08,1.0]
	ff=PROPS(3)!ff,failure porosity[0.25,0.1,1.0]
	fn=PROPS(4)!fn,volume fraction of graphite inclusions[0.04,0.1,0.25]
	sn=PROPS(5)!sn,standard deviation of graphite inclusions[0.001,0.001,-]
	en=PROPS(6)!en,mean value of graphite inclusions[0.0,0.0,-]
	qinf=PROPS(7)!Isotropic hardening amplitude
	b=PROPS(8)!Isotropic hardening exponent sflow=syield+qinf*(1-EXP(-b*ep))
	C=PROPS(9)!Cyclic hardening parameter
	ga=PROPS(10)!Cyclic saturation parameter. In von Mises plasticity C/ga=qinf,ga=b
	EE=PROPS(11)!Young's modulus
	nu=PROPS(12)!Poisson's ratio
	syield=PROPS(13)!Yield stress
	IF (DISP>=1) THEN
		print *, '***** PROPS Loaded *****'
	END IF
end subroutine
subroutine LoadSTATEV(f,ep,DEp,DEq,EELAS,EPLAS,X,STATEV,NSTATV,NTENS,NDI,NSHR,DROT,DISP)
	implicit none
	real(8),intent(inout) :: f,ep,DEp,DEq
	integer,intent(in) :: NSTATV,NTENS,NDI,NSHR,DISP
	real(8),dimension(6),intent(inout) :: EELAS,EPLAS,X
	real(8),dimension(NSTATV),intent(in) :: STATEV
	real(8),dimension(3,3),intent(in) :: DROT
	IF (DISP>=1) THEN
		print *, '***** Loading STATEV *****'
	END IF
	!STATEV(1)=f:Porosity
	f = STATEV(1)
	!STATEV(2)=ep:equivalent plastic strain
	ep=STATEV(2)
	!STATEV(3)=DEp: Hydrostatic inelastic strain increment
	DEp=STATEV(3)
	!STATEV(4)=DEp: Deviatoric inelastic strain increment
	DEq=STATEV(4)
	!ROTATE TENSORS
	!STATEV(5)...STATEV(10)=EELAS
	CALL ROTSIG(STATEV(5),DROT,EELAS,2,NDI,NSHR)
	!STATEV(11)...STATEV(16)=EPLAS
	CALL ROTSIG(STATEV(5+NTENS),DROT,EPLAS,2,NDI,NSHR)
	!STATEV(17)...STATEV(22)=X
	CALL ROTSIG(STATEV(5+2*NTENS),DROT,X,1,NDI,NSHR)
	IF (DISP>=1) THEN
		print *, '***** STATEV Loaded *****'
	END IF
end subroutine
subroutine UpdateSTATEV(f,ep,DEp,DEq,EELAS,EPLAS,X,STATEV,NSTATV,NTENS,DISP)
	implicit none
	real(8),intent(in) :: f,ep,DEp,DEq
	integer,intent(in) :: NSTATV,NTENS,DISP
	real(8),dimension(6),intent(in) :: EELAS,EPLAS,X
	real(8),dimension(NSTATV),intent(inout) :: STATEV
	integer :: I
	IF (DISP>=1) THEN
		print *, '***** Updating STATEV *****'
	END IF
	STATEV(1)=f
	STATEV(2)=ep
	STATEV(3)=DEp
	STATEV(4)=DEq
	DO I=1,NTENS
		STATEV(4+I)=EELAS(I)
		STATEV(4+I+NTENS)=EPLAS(I)
		STATEV(4+I+2*NTENS)=X(I)
	END DO
	IF (DISP>=1) THEN
		print *, '***** STATEV Updated *****'
	END IF
end subroutine
subroutine Updateft(f,ft,dftdf,fu,fc,ff)
	implicit none
	real(8),intent(in) :: f,fu,fc,ff
	real(8),intent(inout) :: ft,dftdf
	real(8) :: df,k,a,b,c
	k = (fu-fc)/(ff-fc)
	df = 0.1*fc
	if (f<=(fc-df)) then
		ft = f
		dftdf = 1.0
	elseif ((f>(fc-df)) .and. (f<=(fc+df))) then
		a = (k-1.0)/(4.0*df)
		b = fc+(1.0+k)/(1.0-k)*df
		c = fc+k/(1-k)*df
		ft = a*(f-b)**2+c
		dftdf = 2*a*(f-b)
	else
		ft = fc + k*(f-fc)
		dftdf = k
	end if
end subroutine
subroutine IdentityTensor(II)
	implicit none
	real(8), dimension(3,3,3,3), intent(inout) :: II
	integer :: i,j,k,l
	DO i=1,3
		DO j=1,3
			DO k=1,3
				DO l=1,3
					II(i,j,k,l) = 0.5*(Delta(i,k)*Delta(j,l) + Delta(i,l)*Delta(j,k))
				END DO
			END DO
		END DO
	END DO
CONTAINS
	REAL function Delta(i,j)
		implicit none
		integer, intent(in) :: i,j
		IF (i==j) THEN
			Delta = 1.0
		ELSE
			Delta = 0.0
		END IF
	END FUNCTION Delta
end subroutine
subroutine TensorDyad(a,b,AAt)
	implicit none 
	real(8), dimension(6), intent(in) :: a,b
	real(8), dimension(3,3,3,3), intent(inout) :: AAt
	real(8), dimension(3,3) :: at,bt
	integer :: i,j,k,l
	
	! Convert vectors to tensor
	call vector2tensor(a,at)
	call vector2tensor(b,bt)
	
	! Calculate A_ijkl = a_ij b_kl
	DO i=1,3
		DO j=1,3
			DO k=1,3
				DO l=1,3
					AAt(i,j,k,l) = at(i,j)*bt(k,l)
				END DO
			END DO
		END DO
	END DO
end subroutine
subroutine TensorProd(A,B,RES)
	implicit none
	real(8), dimension(3,3,3,3), intent(inout) :: A,B,RES
	integer :: i,j,k,l,s,t
	! Calculates C_ijkl = A_ijst B_stkl
	RES(1:3,1:3,1:3,1:3) = 0.0
	DO i=1,3
		DO j=1,3
			DO k=1,3
				DO l=1,3
					DO s=1,3
						DO t=1,3
							RES(i,j,k,l) = RES(i,j,k,l) + A(i,j,s,t)*B(s,t,k,l)
						END DO
					END DO
				END DO
			END DO
		END DO
	END DO	
end subroutine
subroutine TensorQuad(A,B,C,RES)
	implicit none
	real(8), dimension(3,3,3,3), intent(inout) :: A,B,C,RES
	real(8), dimension(3,3,3,3) :: temp
	! Calculates D_ijkl = A_ijst B_stpq C_pqkl
	call TensorProd(B,C,temp)
	call TensorProd(A,temp,RES)
end subroutine
subroutine vector2tensor(a,AA)
	implicit none
	real(8), dimension(6), intent(in) :: a
	real(8), dimension(3,3), intent(inout) :: AA
	AA(1,1) = a(1)
	AA(2,2) = a(2)
	AA(3,3) = a(3)
	AA(3,2) = a(4)
	AA(2,3) = a(4)
	AA(1,3) = a(5)
	AA(3,1) = a(5)
	AA(1,2) = a(6)
	AA(2,1) = a(6)
end subroutine
subroutine reduce4tensor(AAt, AA)
	implicit none
	real(8), dimension(3,3,3,3), intent(in) :: AAt
	real(8), dimension(6,6), intent(inout) :: AA
	integer, dimension(2,6) :: idx
	integer :: row, col,i,j,k,l
	! Mapping for the indices
	idx(1,:) = (/1,2,3,2,1,1/)
	idx(2,:) = (/1,2,3,3,3,2/)
	DO row=1,6
		i = idx(1,row)
		j = idx(2,row)
		DO col=1,6
			k = idx(1,col)
			l = idx(2,col)
			AA(row,col) = AAt(i,j,k,l)
		END DO
	END DO
end subroutine
subroutine MatrixInverse2(A,Ainv)
	implicit none
	real(8), dimension(2,2), intent(in) :: A
	real(8), dimension(2,2), intent(inout) :: Ainv
	real(8) :: detA
	detA = A(1,1)*A(2,2)-A(1,2)*A(2,1)
	Ainv(1,1) = A(2,2)/detA
	Ainv(2,1) = -A(2,1)/detA
	Ainv(1,2) = -A(1,2)/detA
	Ainv(2,2) = A(1,1)/detA
end subroutine
subroutine OUTER_PRODUCT(a,b,C,n)
	implicit none
	real(8), intent(in) :: a(n),b(n)
	integer, intent(in) :: n
	real(8), intent(inout) :: C(n,n)
	integer :: i,j
	do i=1,n
		do j=1,n
			C(i,j) = a(j)*b(i)
		end do
	end do
end subroutine
subroutine LinearSolve(A,n,b,x)
	implicit none

	real(8),intent(in) :: A(n,n),b(n)
	integer,intent(in) :: n
	real(8),intent(inout) :: x(n)
	real(8) :: c(1,n)

	integer :: ipiv(n),info,nrhs,lda,ldb

	nrhs=1
	lda=n
	ldb=n
	c(1,1:n)=b(1:n)
	call dgesv(n,nrhs,A,lda,ipiv,c,ldb,info)
	x(1:n)=c(1,1:n)
end subroutine
subroutine pressure(S,p)
	implicit none
	real(8), dimension(6), intent(in) :: S
	real(8), intent(inout) :: p
	p = -(S(1)+S(2)+S(3))/3.0
end subroutine
subroutine vonmises(S,q)
	implicit none
	real(8), dimension(6), intent(in) :: S
	real(8), intent(inout) :: q
	!real(8) :: p
	!real(8),dimension(6)::onevec
	q = sqrt(0.5*((S(1)-S(2))**2+(S(2)-S(3))**2+(S(1)-S(3))**2&
		+6.0*(S(4)**2+S(5)**2+S(6)**2)))
	!call pressure(S,p)
	!onevec=(/1.0,1.0,1.0,0.0,0.0,0.0/)
	!q = sqrt(1.5*DOT_PRODUCT(S+p*onevec,S+p*onevec))
end subroutine
subroutine mcosh(x,c)
	implicit none
	real(8),intent(in) :: x
	real(8),intent(inout) :: c
	real(8) :: xc
	xc = 30.0
	if (abs(x)<=xc) then
		c = cosh(x)
	else
		c = cosh(xc)+sinh(xc)*(abs(x)-xc)+0.5*cosh(xc)*(abs(x)-xc)**2
	end if
end subroutine
subroutine msinh(x,s)
	implicit none
	real(8), intent(in) :: x
	real(8), intent(inout) :: s
	real(8) :: xc
	xc = 30.0
	if (abs(x)<=xc) then
		s = sinh(x)
	else
		s = (cosh(xc)+sinh(xc)*(abs(x)-xc)+0.5*cosh(xc)*(abs(x)-xc)**2)*x/abs(x)
	end if
end subroutine
subroutine NumericalConsistentTangent(STRESS,STATEV,DSTRAN,DDSDDE,NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,DROT,DISP)
	implicit none
	REAL*8,INTENT(IN) :: STRESS(NTENS),STATEV(NSTATV),&
		DSTRAN(NTENS),PROPS(NPROPS),DROT(3,3)
	REAL*8,INTENT(INOUT) :: DDSDDE(NTENS,NTENS)
	integer,intent(in) :: NDI,NSHR,NTENS,NSTATV,NPROPS,DISP
	INTEGER :: I,J,MAXITERS
	REAL*8 :: tol,res
	REAL*8 :: f0,fc,ff,fn,sn,en,qinf,b,C,ga,EE,nu,syield
	REAL*8 :: K,G,la
	REAL*8 :: q1,q2,q3,fu
	REAL*8 :: pi,cosht,sinht,dphidq,dphidp
	REAL*8 :: DEp,DEq,ep,f,fnn
	REAL*8 :: ft, dftdf,ftn
	REAL*8 :: p,q,sflow,phi,qplus,pplus
	REAL*8 :: eps,eta,typx
	REAL*8,DIMENSION(6) :: STRIAL,X,EELAS,EPLAS,BETA,DEin,DSTRAIN,STRESSI,BETAplus,STRAIN
	REAL*8,DIMENSION(6) :: N,onevec
	REAL*8,DIMENSION(6,6) :: CELAS
	REAL*8,DIMENSION(3,3) :: IdentityMatrix
	! Constants
	pi=3.14159265358979
	tol = 1e-9
	MAXITERS = 20
	
	IF (DISP>=1) THEN
		print *, '***** Inside NumericalConsistentTangent *****'
	END IF
	! Metzger values
	!typx=1.0e-9
	typx=1.0
	eta=1.0e-5
	IdentityMatrix(1:3,1) = (/1.0, 0.0, 0.0/)
	IdentityMatrix(1:3,2) = (/0.0, 1.0, 0.0/)
	IdentityMatrix(1:3,3) = (/0.0, 0.0, 1.0/)
	
	!Load properties
	call LoadPROPS(f0,fc,ff,fn,sn,en,qinf,b,C,ga,EE,nu,syield,PROPS,NPROPS,DISP)
	
	!Bulk and shear moduli
	K=EE/(3.0*(1.0-nu))
	G=EE/(2.0*(1.0+nu))
	la=K-2.0/3.0*G

	!Tvergaard values for periodic void distributions
	q1=1.5
	q2=1.0
	q3=q1**2

	!Effective porosity at failure
	fu=(q1-SQRT(q1**2-q3))/q3

	! Load state variables
	call LoadSTATEV(f,ep,DEp,DEq,EELAS,EPLAS,X,STATEV,NSTATV,&
		NTENS,NDI,NSHR,IdentityMatrix,DISP)
	
	fnn = f
	!Effective linear elasticity tensor CELAS
	CELAS(1:NTENS,1:NTENS)=0.0
	DO I=1,NDI
		DO J=1,NDI
			CELAS(I,J)=la
		END DO 
		CELAS(I,I)=2.0*G+la
	END DO 
	DO I=NDI+1,NTENS
		CELAS(I,I)=G
	END DO
	
	!DEp = DEp/10.0
	!DEq = DEq/10.0

	
	DO I=1,NTENS
		IF (ABS(EELAS(I)+EPLAS(I))>typx) THEN
			eps=eta*ABS(EELAS(I)+EPLAS(I))*SIGN(1.0,DSTRAN(I))
		ELSE
			eps=eta*typx*SIGN(1.0,DSTRAN(I))
		END IF
		! Metzger values for perturbation
		!IF (ABS(DSTRAN(I))>typx) THEN
		!	eps = sqrt(eta)*DSTRAN(I)
		!ELSE
		!	eps = sqrt(eta)*typx*SIGN(1.0,DSTRAN(I))
		!END IF
		DSTRAIN(1:NTENS)=0.0
		DSTRAIN(I)=eps
		! Trial stress
		STRIAL=MATMUL(CELAS,EELAS+DSTRAIN)
		IF (DISP>=3) THEN
			print *, "DSTRAIN"
			print *, DSTRAIN
			print *, "STRESS(n):"
			print *, STRESS
			print *, "STRIAL:"
			print *, STRIAL
			print *, "X(n):"
			print *, X
		END IF
		! Relative stress
		BETA=STRIAL-X
		
		call pressure(BETA,p)
		call vonmises(BETA,q)
		IF (DISP>=3) THEN
			print *, "Pressure:", p, "Von Mises:", q
		END IF
		! Flow stress
		sflow = syield+qinf*(1.0-exp(-b*ep))

		! Update effective porosity and its derivative
		call Updateft(f,ft,dftdf,fu,fc,ff)

		! Yield function
		call mcosh(1.5*q2*p/sflow,cosht)
		call msinh(1.5*q2*p/sflow,sinht)
		phi = (q/sflow)**2 + 2.0*ft*q1*cosht-(1.0+q3*ft**2)
		IF (phi<0.0) THEN
			! Elastic step
			! Update DDSDDE
			DDSDDE(1:NTENS,I) = CELAS(1:NTENS,I)
			IF (DISP>=2) THEN
				PRINT *, "***** Elastic step *****"
			END IF
		ELSE
			! Plastic step
			res = 1.0
			! 12.11.2013 Changed plastic starting-guess
			! Close to ideal-plastic starting guess
			! Least-squares method, N:N = 3/2
			DEp = DSTRAIN(1)+DSTRAIN(2)+DSTRAIN(3)
			!IF (ABS(DEp)<TOL) THEN
			!	DEp = tol*sign(1.0,DEp)
			!END IF
			onevec=(/1.0,1.0,1.0,0.0,0.0,0.0/)
			N = 1.5*(BETA+p*onevec)/q
			DEq = 2.0/3.0*DOT_PRODUCT(N,DSTRAIN-DEp/3.0*onevec)
			! Returns iteration variables: f, ep, DEp, DEq and residual res
			call PlasticIteration(f,ep,DEp,DEq,DDSDDE,STRIAL,X,MAXITERS,TOL,&
				PROPS,NPROPS,NTENS,NDI,res,DISP-1)
			! Check for convergence
			IF (res>tol) THEN
				PRINT *, "***** Plastic iteration failed at consistent tangent iteration direction", I, "*****"
			!	PRINT *, "***** Using elastic tangent *****"
			!	DDSDDE(1:NTENS,I) = CELAS(1:NTENS,I)
			ELSE
				DEin = EPLAS
				! Update STRESS,X,DEin
				call UpdateStress(f,ep,DEp,DEq,STRIAL,X,PROPS,NPROPS,&
					STRESSI,DEin,fnn,NTENS,NDI,DISP)
				DDSDDE(1:NTENS,I) = (STRESSI-STRESS)/eps
			END IF
		END IF
	END DO
	IF (DISP>=1) THEN
		print *, '***** Finished NumericalConsistentTangent *****'
	END IF
end subroutine
