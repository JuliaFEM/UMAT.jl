! This file is a part of JuliaFEM.
! License is MIT: see https://github.com/JuliaFEM/UMAT.jl/blob/master/LICENSE
!***********************************************************************
! File drucker_prager_plasticity_deriv.f90
!
!   Subroutine DF_Derivatives
!   Subroutine DFDSIGMA
!   Subroutine DEVSTRESS
!   Function   YIELDF
!   Function   J2INVARIANT
!   Function   HARDENING_K
!   Function   HARDENING_K_DERIV
!   Function   DPHIDF
!   subroutine invert
!   Subroutine WriteOut
!
!***********************************************************************
!=======================================================================
! Subroutine DF_Derivatives
!
! to calculate Derivates of Yield function F
!-----------------------------------------------------------------------
!
! Input parameters
! ----------------
! M           = number of stress components
! STRESS      = stress vector
! DEV_STRESS  = deviatoric stress vector
! PARAM       = some material parameters
! HIDDEN      = hidden variables
! IFLAG       = flag for yield function type
!             = 'von_Mises'
!             = 'Drucker-Prager_Overstress_A'
!             = 'Drucker-Prager_Consistency_A'
! HFLAG       = flag for Hardening function type
!     1       = Linear
!     2       = Exponential
! LambdaDot   = Lambda/DT
!
! Output parameters
! ------------------
!   DF_DLambda
!   DF_DLambdaDot
!   DDF_DDSigma
!   DDF_DSigmaDLambda
!   DDF_DSigmaDLambdaDot
!   ALPHA
!-----------------------------------------------------------------------
! Originally coded 19/03/2003 KKO
! last modified
!=======================================================================
!
Subroutine DF_Derivatives(M, IFLAG,HFLAG,STRESS,DEV_STRESS,PARAM,      &
                          HIDDEN,LambdaDot,DF_DLambda,DF_DLambdaDot,   &
                          DDF_DDSigma,DDF_DSigmaDLambda,               &
                          DDF_DSigmaDLambdaDot,ALPHA)
  Implicit None
  Integer, Parameter      :: dp = Selected_real_kind(12,50)
  Integer, Intent(in)     :: M,HFLAG
  Character(*), Intent(in):: IFLAG
  Real(dp), Intent(in)    :: STRESS(M),DEV_STRESS(M),PARAM(*),HIDDEN(*),&
                             LambdaDot
  Real(dp), Intent(out)   :: DF_DLambda,DF_DLambdaDot,DDF_DDSigma(M,M),  &
                             DDF_DSigmaDLambda(M),DDF_DSigmaDLambdaDot(M),&
                             ALPHA
!- local variables -----------------------------------------------------
  Integer :: II,JJ
  Real(dp) :: nvec(M),dum,SLEN,J2INVARIANT,DHardeningDeffStrain,&
              Sigma_t,Sigma_C,BEETA,N_Power, A_Const, B_Const,ETA

  IF (M /= 6) Then
    Write (*,*) ' *** Error **** in Subroutine DF_Derivatives, 6 Stress'&
                ' components needed'
    Call XIT
  Endif

  If(INDEX(IFLAG,'von_Mises') /= 0) Then

    SLEN    = SQRT(2.0_dp * J2INVARIANT(M, DEV_STRESS) )
    DF_DLambda = DHardeningDeffStrain(PARAM, HIDDEN, HFLAG)
    DF_DLambdaDot        = 0.0_dp
    DDF_DSigmaDLambda    = 0.0_dp
    DDF_DSigmaDLambdaDot = 0.0_dp

  Elseif (INDEX(IFLAG,'Drucker-Prager_Overstress_A') /= 0 ) Then

    SIGMA_T = PARAM(1)
    SIGMA_C = PARAM(2)
    ALPHA   = (SIGMA_C - SIGMA_T) / (SIGMA_C + SIGMA_T)
    BEETA   = (1.0_dp + ALPHA) * SIGMA_T
    SLEN    = SQRT(2.0_dp * J2INVARIANT(M, DEV_STRESS) )

    DF_DLambda = DHardeningDeffStrain(PARAM, HIDDEN, HFLAG)*(1.0_dp+Alpha)
    DF_DLambdaDot        = 0.0_dp
    DDF_DSigmaDLambda    = 0.0_dp
    DDF_DSigmaDLambdaDot = 0.0_dp

  Elseif (INDEX(IFLAG,'Drucker-Prager_Consistency_A') /= 0    ) Then

    SIGMA_T = PARAM(1)
    SIGMA_C = PARAM(2)
    ETA     = PARAM(6)
    N_Power = PARAM(5)
    A_Const = PARAM(7)
    B_Const = PARAM(8)
    ALPHA   = (SIGMA_C - SIGMA_T) / (SIGMA_C + SIGMA_T)
    BEETA   = (1.0_dp + ALPHA) * SIGMA_T
    SLEN    = SQRT(2.0_dp * J2INVARIANT(M, DEV_STRESS) )

    DF_DLambda = DHardeningDeffStrain(PARAM, HIDDEN, HFLAG)*(1.0_dp+Alpha)
    DF_DLambdaDot = -A_Const/ETA/N_Power*(LambdaDot/ETA)**(1.0_dp/N_Power-1.0_dp)
    DDF_DSigmaDLambda    = 0.0_dp
    DDF_DSigmaDLambdaDot = 0.0_dp

  Endif

  If(INDEX(IFLAG,'von_Mises') /= 0                   .or. &
     INDEX(IFLAG,'Drucker-Prager_Overstress_A') /= 0 .or. &
     INDEX(IFLAG,'Drucker-Prager_Consistency_A') /= 0    ) Then

    nvec(1) = DEV_STRESS(1)/SLEN
    nvec(2) = DEV_STRESS(2)/SLEN
    nvec(3) = DEV_STRESS(3)/SLEN
    nvec(4) = 2.0_dp*DEV_STRESS(4)/SLEN
    nvec(5) = 2.0_dp*DEV_STRESS(5)/SLEN
    nvec(6) = 2.0_dp*DEV_STRESS(6)/SLEN
    dum = SQRT(1.5)/SLEN
    DO JJ=1,M
      DO II = 1,M
        DDF_DDSigma(II,JJ) = -dum*nvec(II)*nvec(JJ)
      Enddo
    Enddo
    dum = SQRT(1.5)/SLEN/3.0_dp
    DDF_DDSigma(1,1) = DDF_DDSigma(1,1) + dum * 2.0_dp
    DDF_DDSigma(2,2) = DDF_DDSigma(2,2) + dum * 2.0_dp
    DDF_DDSigma(3,3) = DDF_DDSigma(3,3) + dum * 2.0_dp
    DDF_DDSigma(4,4) = DDF_DDSigma(4,4) + dum * 6.0_dp
    DDF_DDSigma(5,5) = DDF_DDSigma(5,5) + dum * 6.0_dp
    DDF_DDSigma(6,6) = DDF_DDSigma(6,6) + dum * 6.0_dp

    dum = SQRT(1.5)/SLEN/3.0_dp
    DDF_DDSigma(1,2) = DDF_DDSigma(1,2) - dum
    DDF_DDSigma(1,3) = DDF_DDSigma(1,3) - dum
    DDF_DDSigma(2,1) = DDF_DDSigma(2,1) - dum
    DDF_DDSigma(2,3) = DDF_DDSigma(2,3) - dum
    DDF_DDSigma(3,1) = DDF_DDSigma(3,1) - dum
    DDF_DDSigma(3,2) = DDF_DDSigma(3,2) - dum

  Endif
!
End Subroutine DF_Derivatives
!
!***********************************************************************
!=======================================================================
! Subroutine DFDSIGMA
!
! to evaluate plastic potential derivate over stress
!-----------------------------------------------------------------------
! input parameters
! ----------------
! M           = number of stress components
! STRESS      = stress vector
! DEV_STRESS  = deviatoric stress vector
! PARAM       = some material parameters
! HIDDEN      = hidden variables
! IFLAG       = flag for yield function type
!             = 'von_Mises'
!             = 'Drucker-Prager_Overstress_A'
!             = 'Drucker-Prager_Consistency_A'
! HFLAG  = flag for Hardening function type
!     1  = Linear
!     2  = Exponential
! Output parameters
!
! ------------------
! DF_DSIGMA(M)
!-----------------------------------------------------------------------
! Originally coded 19/03/2003 KKO
! last modified
!=======================================================================
!
Subroutine DFDSIGMA(M, STRESS, DEV_STRESS, PARAM, HIDDEN, &
                  IFLAG,HFLAG,DF_DSIGMA)
  Implicit None
  Integer, Parameter        :: dp = Selected_real_kind(12,50)
  Integer, Intent(in)       :: M,HFLAG
  Character(*), Intent(in)  :: IFLAG
  Real(dp), Intent(in)      :: STRESS(M), DEV_STRESS(M), PARAM(*), &
                               HIDDEN(*)
  Real(dp), Intent(out)     :: DF_DSIGMA(M)
!- local variables -----------------------------------------------------
  Real(dp) :: SIGMA_T,SIGMA_C,ALPHA,BEETA,dum,HARDENING_K,&
              SM,K,SLEN,J2INVARIANT
!
  IF (M /= 6) Then
    Write (*,*) ' *** Error **** in Function DFDSIGMA, 6 Stress '&
                'components needed'
    Call XIT
  Endif

  SM    = (STRESS(1) + STRESS(2) + STRESS(3)) / 3.0_dp

  If(INDEX(IFLAG,'von_Mises') /= 0) Then

    K = HARDENING_K(PARAM, HIDDEN, HFLAG)
    SLEN    = SQRT(2.0_dp * J2INVARIANT(M, DEV_STRESS) )
    dum     = SQRT(1.5_dp)/SLEN
    DF_DSIGMA(1) = dum*DEV_STRESS(1)
    DF_DSIGMA(2) = dum*DEV_STRESS(2)
    DF_DSIGMA(3) = dum*DEV_STRESS(3)
    DF_DSIGMA(4) = dum*DEV_STRESS(4)  * 2.0_dp
    DF_DSIGMA(5) = dum*DEV_STRESS(5)  * 2.0_dp
    DF_DSIGMA(6) = dum*DEV_STRESS(6)  * 2.0_dp

  Elseif (INDEX(IFLAG,'Drucker-Prager_Overstress_A') /= 0 .or. &
          INDEX(IFLAG,'Drucker-Prager_Consistency_A') /= 0    ) Then

    SIGMA_T = PARAM(1)
    SIGMA_C = PARAM(2)
    !C_const = PARAM(4)
    !K_Zero  = PARAM(3)
    !N_Power = PARAM(5)
    !A_Const = PARAM(7)
    !B_Const = PARAM(8)
    K = HARDENING_K(PARAM, HIDDEN, HFLAG)
    ALPHA   = (SIGMA_C - SIGMA_T) / (SIGMA_C + SIGMA_T)
    BEETA   = (1.0_dp + ALPHA) * SIGMA_T
    SLEN    = SQRT(2.0_dp * J2INVARIANT(M, DEV_STRESS) )
    dum     = SQRT(1.5_dp)/SLEN
    DF_DSIGMA(1) = dum*DEV_STRESS(1)  + ALPHA
    DF_DSIGMA(2) = dum*DEV_STRESS(2)  + ALPHA
    DF_DSIGMA(3) = dum*DEV_STRESS(3)  + ALPHA
    DF_DSIGMA(4) = dum*DEV_STRESS(4)  * 2.0_dp
    DF_DSIGMA(5) = dum*DEV_STRESS(5)  * 2.0_dp
    DF_DSIGMA(6) = dum*DEV_STRESS(6)  * 2.0_dp

  End If

End Subroutine DFDSIGMA

!***********************************************************************
!=======================================================================
! Subroutine DEVSTRESS
!
! to calculate deviatoric stress vector
!-----------------------------------------------------------------------
!
! Input parameters
! ----------------
!   M      = number of stress components
!   STRESS = stress vector
!
! Output parameters
! ------------------
!   DEV_STRESS
!   SM          = Mean Stress
!-----------------------------------------------------------------------
! Originally coded 19/03/2003 KKO
! last modified
!=======================================================================
!
Subroutine DEVSTRESS(M, STRESS,SM,DEV_STRESS)
  Implicit None
  Integer, Parameter        :: dp = Selected_real_kind(12,50)
  Integer, Intent(in)       :: M
  Real(dp), Intent(in)      :: STRESS(M)
  Real(dp), Intent(out)     :: SM,DEV_STRESS(M)
!- local variables -----------------------------------------------------

  IF (M /= 6) Then
    Write (*,*) ' *** Error **** in Subroutine DEVSTRESS, 6 Stress '&
                'components needed'
    Call XIT
  Endif

  SM    = (STRESS(1) + STRESS(2) + STRESS(3)) / 3.0_dp
  DEV_STRESS(1) = STRESS(1) - SM
  DEV_STRESS(2) = STRESS(2) - SM
  DEV_STRESS(3) = STRESS(3) - SM
  DEV_STRESS(4) = STRESS(4)
  DEV_STRESS(5) = STRESS(5)
  DEV_STRESS(6) = STRESS(6)
!
End Subroutine DEVSTRESS
!
!***********************************************************************
!=======================================================================
! function YIELDF
!
! to evaluate the value of yield function
!-----------------------------------------------------------------------
! input parameters
! ----------------
! M      = number of stress components
! STRESS = stress vector
! PARAM  = some material parameters
! HIDDEN = hidden variables
! HFLAG       = flag for Hardening function type
!     1       = Linear
!     2       = Exponential
! IFLAG  = flag for yield function type
!        = 'von_Mises'
!        = 'Drucker-Prager_Overstress_A'
!        = 'Drucker-Prager_Consistency_A'
!-----------------------------------------------------------------------
! Originally coded RK (2001)
! last modified 19/03/2003 KKO
!=======================================================================
!
Function YIELDF(M, STRESS, PARAM, HIDDEN, IFLAG,HFLAG) Result(YIELD)
  Implicit None
  Integer, Parameter        :: dp = Selected_real_kind(12,50)
  Integer, Intent(in)       :: M,HFLAG
  Character(*), Intent(in)  :: IFLAG
  Real(dp), Intent(in)      :: STRESS(M), PARAM(*), HIDDEN(*)
!- local variables -----------------------------------------------------
  Real(dp) :: SM,SX,SY,SZ,SXY,SYZ,SXZ,YLD2,YLD,YIELD,SIGMA_T,&
              SIGMA_C,ALPHA,BEETA,RJ2,K,K_Zero,C_const,HARDENING_K
!
  IF (M /= 6) Then
    Write (*,*) ' *** Error **** in Function YIELDF, 6 Stress '&
                'components needed'
    Call XIT
  Endif

  SM    = (STRESS(1) + STRESS(2) + STRESS(3)) / 3.0_dp
  SX    = STRESS(1) - SM
  SY    = STRESS(2) - SM
  SZ    = STRESS(3) - SM
  SXY   = STRESS(4)
  SYZ   = STRESS(5)
  SXZ   = STRESS(6)

  If(INDEX(IFLAG,'von_Mises') /= 0) Then
    YLD   = PARAM(1)
    YLD2  = YLD*YLD
    YIELD = 1.5_dp*(SX*SX+SY*SY+SZ*SZ+2.0_dp*(SXY*SXY+SYZ*SYZ+SXZ*SXZ))-YLD2

  Elseif (INDEX(IFLAG,'Drucker-Prager_Overstress_A') /= 0 .or. &
          INDEX(IFLAG,'Drucker-Prager_Consistency_A') /= 0) Then
    SIGMA_T = PARAM(1)
    SIGMA_C = PARAM(2)
    K_Zero  = PARAM(3)
    C_const = PARAM(4)
    !N_Power = PARAM(5)
    !A_Const = PARAM(7)
    !B_Const = PARAM(8)
    K = HARDENING_K(PARAM, HIDDEN, HFLAG)
    ALPHA   = (SIGMA_C - SIGMA_T) / (SIGMA_C + SIGMA_T)
    BEETA    = (1.0_dp + ALPHA) * SIGMA_T
    RJ2     = 0.5*(SX*SX+SY*SY+SZ*SZ)+SXY*SXY+SXZ*SXZ+SYZ*SYZ
    YIELD  = Sqrt(3.0_dp*RJ2) + ALPHA*3.0_dp*SM - BEETA - K
  End If
End Function YIELDF
!
!***********************************************************************
!=======================================================================
! function J2INVARIANT
!
! to evaluate Deviatoric stress 2nd invariant
!-----------------------------------------------------------------------
! input parameters
! ----------------
! M           = number of stress components
! DEV_STRESS  = deviatoric stress vector
!-----------------------------------------------------------------------
! Originally coded 19/03/2003 KKO
! last modified
!=======================================================================
!
Function J2INVARIANT(M, DEV_STRESS) Result(J2Inv)
  Implicit None
  Integer, Parameter        :: dp = Selected_real_kind(12,50)
  Integer, Intent(in)       :: M
  Real(dp), Intent(in)      :: DEV_STRESS(M)
!- local variables -----------------------------------------------------
  Real(dp) :: J2Inv
!
  IF (M /= 6) Then
    Write (*,*) ' *** Error **** in Function J_INVARIANTS, 6 Stress '&
                'components needed'
    Call XIT
  Endif

    J2Inv = 0.5_dp*(Dev_Stress(1)**2+Dev_Stress(2)**2+Dev_Stress(3)**2)&
                 +  Dev_Stress(4)**2+Dev_Stress(5)**2+Dev_Stress(6)**2

End Function J2INVARIANT
!
!***********************************************************************
!=======================================================================
! function HARDENING_K
!
! to evaluate Hardening parameter K
!-----------------------------------------------------------------------
! input parameters
! ----------------
! PARAM  = some material parameters
! HIDDEN = hidden variables
! HFLAG  = flag for Hardening function type
!     1  = Linear
!     2  = Exponential
!-----------------------------------------------------------------------
! Originally coded 19/03/2003 KKO
! last modified
!=======================================================================
!
Function HARDENING_K(PARAM, HIDDEN, HFLAG)  Result(K)
  Implicit None
  Integer, Parameter        :: dp = Selected_real_kind(12,50)
  Integer, Intent(in)       :: HFLAG
  Real(dp), Intent(in)      :: PARAM(*),HIDDEN(*)
!- local variables -----------------------------------------------------
  Real(dp) :: K
!

  If(HFLAG == 1) Then
    K = -PARAM(3) *HIDDEN(1)
  Elseif(HFLAG == 2) Then
    K = -PARAM(3) *(1.0_dp-EXP(-PARAM(4)*HIDDEN(1)))
  Endif

End Function HARDENING_K
!
!***********************************************************************
!=======================================================================
! function DHardeningDeffStrain
!
! to evaluate derivatives of parameter K
!-----------------------------------------------------------------------
! input parameters
! ----------------
! PARAM  = some material parameters
! HIDDEN = hidden variables
! HFLAG  = flag for Hardening function type
!     1  = Linear
!     2  = Exponential
!-----------------------------------------------------------------------
! Originally coded 19/03/2003 KKO
! last modified
!=======================================================================
!
Function DHardeningDeffStrain(PARAM, HIDDEN, HFLAG)  Result(DK)
  Implicit None
  Integer, Parameter        :: dp = Selected_real_kind(12,50)
  Integer, Intent(in)       :: HFLAG
  Real(dp), Intent(in)      :: PARAM(*),HIDDEN(*)
!- local variables -----------------------------------------------------
  Real(dp) :: DK
!

  If(HFLAG == 1) Then
    DK = -PARAM(3)
  Elseif(HFLAG == 2) Then
    DK = -PARAM(4)*PARAM(3) *EXP(-PARAM(4)*HIDDEN(1))
  Endif

End Function DHardeningDeffStrain
!
!***********************************************************************
!=======================================================================
! function DPHIDF
!
! to evaluate derivative of the function Phi inside the Heaviside step
! function in Perzyna's viscoplastic modell
!-----------------------------------------------------------------------
! input parameters
! ----------------
!
! PARAM  = some material parameters
! STRESS = stress vector
! M      = number of stress components
! HIDDEN = hidden variables
! HFLAG  = flag for Hardening function type
!     1  = Linear
!     2  = Exponential
! IFLAG  = flag for yield function type
!        = 'von_Mises'
!        = 'Drucker-Prager_Overstress_A'
!        = 'Drucker-Prager_Consistency_A'
!-----------------------------------------------------------------------
! Originally coded 19/03/2003 KKO
! last modified
!=======================================================================
!
Function DPHIDF(F,PARAM,STRESS,M, HIDDEN, IFLAG,HFLAG)  Result(DPhi_Df)
  Implicit None
  Integer, Parameter        :: dp = Selected_real_kind(12,50)
  Character(*), Intent(in)  :: IFLAG
  Integer , Intent(in)      :: M,HFLAG
  Real(dp), Intent(in)      :: PARAM(*),HIDDEN(*),STRESS(M)
!- local variables -----------------------------------------------------
  Real(dp) :: DPhi_Df,N_Power, A_Const, B_Const,F,YIELDF
!

  If(INDEX(IFLAG,'von_Mises') /= 0) Then

    DPhi_Df = 0.0_dp ! ei m��ritelty

  Elseif (INDEX(IFLAG,'Drucker-Prager_Overstress_A') /= 0 ) Then

    N_Power = PARAM(5)
    A_Const = PARAM(7)
    B_Const = PARAM(8)
    F = YIELDF(M, STRESS, PARAM, HIDDEN, IFLAG,HFLAG)
    DPhi_Df = N_Power / A_Const *(F/A_Const-B_Const)**(N_Power-1.0_dp)

  Elseif (INDEX(IFLAG,'Drucker-Prager_Consistency_A') /= 0    ) Then
  Endif
End Function DPHIDF
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c    routine to invert a matrix by Gaussian elimination
!c     Ainv=inverse(A)    A(n,n)   Ainv(n,n)
!c           D matrix is extended.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine invert( A, Ainv, n )

       Implicit none
       Integer, Parameter :: dp = Selected_real_kind(12,50)
       Integer :: n
       Real(dp), Intent(in) :: A(n,n)
       Real(dp), Intent(out) :: Ainv(n,n)

       Integer :: n2,i,j,k
       Real(dp) :: D(n,2*n),alpha,beta

!  initialize the reduction matrix
       n2 = 2*n
       do 1 i = 1,n
        do 2 j = 1,n
         D(i,j) = A(i,j)
         d(i,n+j) = 0.
 2      continue
        D(i,n+i) = 1.
 1     continue

!  do the reduction
       do 3 i = 1,n
        alpha = D(i,i)
        if(alpha .eq. 0.) go to 300
        do 4 j = 1,n2
         D(i,j) = D(i,j)/alpha
 4      continue
        do 5 k = 1,n
         if((k-i).eq.0) go to 5
         beta = D(k,i)
         do 6 j = 1,n2
          D(k,j) = D(k,j) - beta*D(i,j)
 6       continue
 5      continue
 3     continue

!  copy result into output matrix
       do 7 i = 1,n
        do 8 j = 1,n
         Ainv(i,j) = D(i,j+n)
 8      continue
 7     continue

       return
 300   print *,'*** ERROR: Singular matrix ***'
       return
       end
!  ============================================================================================
!
!   Write to dat,log and msg-files
!
!  ============================================================================================
!
  SUBROUTINE WriteOut(string,form,IND,KSTEP,KINC,NOEL,NPT)

  character*(*), Intent(in) :: string,form
  Integer, Intent(in) :: IND,KSTEP,KINC,NOEL,NPT
  Integer, Parameter :: DATF=6,LOGF=5,MSGF=7

  !Form = '(/,'''//Trim(Task)//LCHR(1:JL)//''','//XCHR(1:XL)//'X'
  !Form = '('''//Trim(Task)//LCHR(1:JL)//''','//XCHR(1:XL)//'X'

  If(KSTEP > 0 ) Then

    If(IND == 1 .or. IND == 4) &
    Write (DATF,'(/,A,I2,A,I5,A,I6,A,I2)') &
      ' Step =',KSTEP,' Increment =',KINC,' Element =',NOEL,' Integration point =',NPT
!    If(IND == 2 .or. IND == 4) &
!    Write (LOGF,'(/,A,I2,A,I5,A,I6,A,I2)') &
!      ' Step =',KSTEP,' Increment =',KINC,' Element =',NOEL,' Integration point =',NPT
    If(IND == 2 .or. IND == 4) &
    Write (MSGF,'(/,A,I2,A,I5,A,I6,A,I2)') &
      ' Step =',KSTEP,' Increment =',KINC,' Element =',NOEL,' Integration point =',NPT
  Endif

  If(IND == 1 .or. IND == 4) &
    Write(DATF,form),string
! If(IND == 2 .or. IND == 4) &
!   Write(LOGF,form),string
  If(IND == 3 .or. IND == 4) &
    Write(MSGF,form),string

  RETURN
  END SUBROUTINE WriteOut
