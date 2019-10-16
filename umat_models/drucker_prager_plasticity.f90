! This file is a part of JuliaFEM.
! License is MIT: see https://github.com/JuliaFEM/UMAT.jl/blob/master/LICENSE
!***********************************************************************
! File drucker_prager_plasticity.f90
! Including:
!
!   subroutine UMAT
!   subroutine UMAT_Overstress
!   subroutine UMAT_Consistency
!

!***********************************************************************
!=======================================================================
! Subroutine UMAT
!
! ABAQUS User Material Subroutine, Calls real material subroutines
! depending on the material type
!-----------------------------------------------------------------------
!
!     Variables to be defined in all situtions
!     - - - - - - - - - - - - - - - - - - - - -
!      DDSDDE(NTENS,NTENS) Jacobian matrix of the constitutive model
!      STRESS(NTENS) This array is passed in as the stress tensor at the
!           beginning of the
!           increment and must be updated in this routine to be the stress
!           tensor at the end of the increment. If initial stresses are
!           specified using the *INITIAL CONDITIONS, TYPE=STRESS option,
!           this array will contain the initial stresses at the start of
!           the analysis. The size of this array depends on the value of
!           NTENS as defined below. In finite-strain problems the stress
!           tensor has already been rotated to account or rigid body motion
!           in the increment before UMAT is called, so that only the
!           corotational part of the stress integration should be done
!           in UMAT. The measure of stress used is "true" (Cauchy) stress.
!      STATEV(NSTATV) An array containing the solution-dependent state
!           variables. These are passed in as the values at the beginning of
!           the increment unless they are updated in user subroutines
!           USDFLD (``USDFLD,'' Section 24.2.37) or UEXPAN (``UEXPAN,''
!           Section 24.2.20), in which case the updated values are passed
!           in. In all cases STATEV must be returned as the values at
!           the end of the increment. The size of the array is that
!           defined on the *DEPVAR option associated with this material.
!           In finite-strain problems any vector-valued or tensor-valued
!           state variables must be rotated to account for rigid body
!           motion of the material, in addition to any update in the
!           values associated with constitutive behavior. The rotation
!           increment matrix, drott, is provided for this purpose.
!      SSE, SPD, SCD
!           Specific elastic strain energy, plastic dissipation, and
!           "creep" dissipation, respectively. These are passed in as
!           the values at the start of the increment and should be
!           updated to the corresponding specific energy values at the
!           end of the increment. They have no effect on the solution,
!           except that they are used for the energy output options.
!
!  The following only in a fully coupled temperature-displacement analysis:
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      RPL  Volumetric heat generation per unit time at the end of the
!           increment caused by mechanical working of the material.
!      DDSDDT(NTENS) Variation of the stress increments with respect to
!           the temperature.
!      DRPLDE(NTENS) Variation of RPL with respect to the strain increments.
!      DRPLDT Variation of RPL with respect to the temperature.
!
!      Variable that can be updated
!      - - - - - - - - - - - - - - -
!       PNEWDT  Ratio of suggested new time increment to the time increment
!           being used (DTIME, see discussion later in this section). This
!           variable allows the user to provide input to the automatic
!           time incrementation algorithms in ABAQUS (if automatic time
!           incrementation is chosen)....
!
!  Variables passed in for information
!  - - - - - - - - - - - - - - - - - -
!        STRAN(NTENS)
!           An array containing the total strains at the beginning of
!           the increment. If the *EXPANSION option is used in the same
!           material definition, the strains passed into UMAT are the
!           mechanical strains only (that is, the thermal expansion strains
!           defined by the *EXPANSION option have been subtracted from the
!           total strains). These strains are available for output as the
!           "elastic" strains. In finite-strain problems the strain
!           components have been rotated to account for rigid body motion
!           in the increment before UMAT is called and are approximations
!           to logarithmic strain.
!        DSTRAN(NTENS)
!           Array of strain increments. If the *EXPANSION option is used
!           in the same material definition, these are the mechanical
!           strain increments (the total strain increments minus the thermal
!           strain increments).
!        TIME(1)
!           Value of step time at the beginning of the current increment.
!        TIME(2)
!           Value of total time at the beginning of the current increment.
!        DTIME
!           Time increment.
!        TEMP
!           Temperature at the start of the increment.
!        DTEMP
!           Increment of temperature.
!        PREDEF
!           Array of interpolated values of predefined field variables at
!           this point at the start of the increment, based on the values
!           read in at the nodes in the *FIELD option.
!        DPRED
!           Array of increments of predefined field variables.
!        CMNAME
!           Name given on the *MATERIAL option, left justified.
!        NDI
!           Number of direct stress components at this point.
!        NSHR
!           Number of engineering shear stress components at this point.
!        NTENS
!           Size of the stress or strain component array (NDI + NSHR).
!        NSTATV
!           Number of solution-dependent state variables that are associated
!           with this material type (as defined in the *DEPVAR option).
!        PROPS(NPROPS)
!           Array of material constants entered in the *USER MATERIAL option
!           for this material.
!        NPROPS
!           Number of material constants (the value given to the CONSTANTS
!           parameter on the *USER MATERIAL option).
!        COORDS
!           An array containing the coordinates of this point. These are
!           the current coordinates if geometric nonlinearity is accounted
!           for during the step (see ``Procedures: overview,'' Section 6.1.1);
!           otherwise, the array contains
!           the original coordinates of the point.
!        drott(3,3)
!           Rotation increment matrix. This matrix represents the increment
!           of rigid body rotation of the basis system in which the components
!           of stress (STRESS) and strain (STRAN) are stored. It is provided so
!           that vector- or tensor-valued state variables can be rotated
!           appropriately in this subroutine: stress and strain components are
!           already rotated by this amount before UMAT is called. This matrix
!           is passed in as a unit matrix for small-displacement analysis and
!           for large-displacement analysis if the basis system for the material
!           point rotates with the material (as in a shell element or when
!           the *ORIENTATION option is used).
!        CELENT
!           Characteristic element length, which is a typical length of a
!           line across an element for a first-order element; it is half of
!           the same typical length for a second-order element. For beams and
!           trusses it is a characteristic length along the element axis.
!           For membranes and shells it is a characteristic length in the
!           reference surface. For axisymmetric elements it is a characteristic
!           length in the  plane only.
!        DFGRD0(3,3)
!           Array containing the deformation gradient at the beginning of
!           the increment. See the discussion regarding the availability of
!           the deformation gradient
!           for various element types.
!        DFGRD1(3,3)
!           Array containing the deformation gradient at the end of the
!           increment. The components of this array are set to zero if
!           the NLGEOM parameter is not used on the *STEP option associated
!           with this increment. See the discussion regarding the availability
!           of the deformation gradient for various element types.
!        NOEL
!           Element number.
!        NPT
!           Integration point number.
!        LAYER
!           Layer number (for composite shells and layered solids).
!        KSPT
!           Section point number within the current layer.
!        KSTEP
!           Step number.
!        KINC
!           Increment number.
!
!-----------------------------------------------------------------------
! Originally coded 19/03/2003 KKO
! last modified
!=======================================================================
  subroutine umat(stress, statev, ddsdde, sse, spd, scd, rpl,           &
          ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp,     &
          dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatv, props,&
          nprops, coords, drott, pnewdt, celent, dfgrd0, dfgrd1, noel,  &
          npt, layer, kspt, kstep, kinc)
!
	implicit none
!  include 'aba_param.inc'
! implicit real*8(a-h,o-z) !aba_param.inc
  INTEGER nprecd           !aba_param.inc
  parameter (nprecd=2)
!
  character*8 cmname
  real*8 :: stress, statev, ddsdde, ddsddt, drplde, stran, dstran,     &
      predef, dpred, props, coords, drott, dfgrd0, dfgrd1,sse,spd,scd, &
      rpl,drpldt,time,dtime,temp,dtemp,pnewdt,celent
	Integer ntens,nstatv,nprops,noel,npt,layer,kspt,kstep,kinc,ndi,      &
          nshr
  dimension stress(ntens), statev(nstatv), ddsdde(ntens, ntens),       &
      ddsddt(ntens), drplde(ntens), stran(ntens), dstran(ntens),       &
      predef(1), dpred(1), props(nprops), coords(3), drott(3, 3),      &
      dfgrd0(3, 3), dfgrd1(3, 3)
!- local variables -----------------------------------------------------
  Integer, Parameter :: DATF=6,MSGF=7,TESTFIL1=105,TESTFIL2=106
  Integer :: II,LENGTH,ISCAN,ELN,IPN,MTYP
  character*80 :: OPEN_File_Name1, OPEN_File_Name2
  Logical :: LOPEN,LNAME,From_Back
  real*8 :: stress_guess(ntens)

! Open files for printing results

  IF (KINC <= 1 .and. KSTEP == 1 .and. npt == 1) Then
    INQUIRE (UNIT=TESTFIL1, OPENED = LOPEN, NAMED=LNAME,NAME=OPEN_File_Name1)
    IF(.NOT. LOPEN)   Then
      INQUIRE (UNIT=6, OPENED = LOPEN, NAMED=LNAME,NAME=OPEN_File_Name1)
      From_Back = .True.
      ISCAN = SCAN(OPEN_File_Name1,'.',From_Back)
      If(ISCAN /= 0) Then
        OPEN_File_Name1(ISCAN:ISCAN+5) ='-1.txt'
        OPEN_File_Name2= OPEN_File_Name1
        OPEN_File_Name1(ISCAN:ISCAN+5) ='-2.txt'
      Else
        LENGTH = LEN_TRIM(OPEN_File_Name1)
        OPEN_File_Name1(LENGTH+1:LENGTH+6) = '-1.txt'
        OPEN_File_Name2= OPEN_File_Name1
        OPEN_File_Name2(LENGTH+1: LENGTH+6) = '-1.txt'
      Endif
      Open(unit=TESTFIL1,status='UNKNOWN',file=OPEN_File_Name1)
      Open(unit=TESTFIL2,status='UNKNOWN',file=OPEN_File_Name2)
    Endif

   ! DO II=1,100
   !   INQUIRE (UNIT=II, OPENED = LOPEN, NAMED=LNAME,NAME=OPEN_File_Name1)
   !   IF(LOPEN) THEN
   !     Write(MSGF,*) 'II,LOPEN',II,LOPEN
   !     Write(MSGF,*) 'LNAME',LNAME
   !     Write(MSGF,*) 'FNAME',OPEN_File_Name1
   !   ENDIF
   ! ENDDO


  Endif


  ELN  =  INT(PROPS(15))
  IPN  =  INT(PROPS(16))
  MTYP =  INT(PROPS(17))

! Update stresses for printing:
  IF (KINC >= 1 .and. NOEL == ELN .and. NPT == IPN ) &
    stress_guess = stress

  If(MTYP == 1) Then
    CALL UMAT_Overstress(stress, statev, ddsdde, sse, spd, scd, stran,  &
          dstran, time, dtime, cmname, ndi, nshr, ntens, nstatv, props, &
          nprops, coords, drott, pnewdt, celent, noel,                  &
          npt, layer, kspt, kstep, kinc)

  ElseIf(MTYP == 2 .or. MTYP == 4)Then
    CALL UMAT_Consistency(stress, statev, ddsdde, sse, spd, scd, stran, &
          dstran, time, dtime, cmname, ndi, nshr, ntens, nstatv, props, &
          nprops, coords, drott, pnewdt, celent, noel,                  &
          npt, layer, kspt, kstep, kinc)
  Endif

! Updete stresses for printing:
  IF (KINC >= 1 .and. NOEL == ELN .and. NPT == IPN ) &
    stress_guess = stress_guess + matmul(ddsdde,dstran)


  IF (KINC >= 1 .and. NOEL == ELN .and. NPT == IPN ) &
    WRITE(TESTFIL1,'(1P,21E13.4)') (stran(II)+dstran(II),II=1,6),&
                (stress(II),II=1,6),statev(1),(stress_guess(II),II=1,6)

  RETURN
	end subroutine umat
!
!***********************************************************************
!=======================================================================
! Subroutine UMAT_Overstress
!
! ABAQUS User Material Subroutine
!-----------------------------------------------------------------------
!
! input parameters
! ----------------
!  PROPS(1)  - Sigma_t
!  PROPS(2)  - Sigma_c
!  PROPS(3)  - K_Zero
!  PROPS(4)  - C_const
!  PROPS(5)  - N_Power
!  PROPS(6)  - Eta
!  PROPS(7)  - A_Const
!  PROPS(8)  - B_Const
!  PROPS(9)  - E
!  PROPS(10) - NU
!  PROPS(11) - DT (If real DTIME = 0, gien as negative)
!  PROPS(12) - 1.0       ! Print All to DAT
!  PROPS(13) - FTOLER (tolerance for local iteration
!  PROPS(14) - MAXITER
!  PROPS(15) - Element number for printing
!  PROPS(16) - Integration point number for printing
!  PROPS(17) - Plasticity type
!         1  - Perzyna Drucker-Prager (overstress)
!         2  - Drucker-Prager Consistency model for Perzyna's overstress
!              modell
!  PROPS(18) - Hardening type
!         1  - Linear hardening
!         2  - Exponential hardening
!
! input/output parameters
! ----------------
!
!  STATEV(1) - Equivalent plastic strain
!  stress    - Stress vector, to be updated
!  ddsdde    - Consistent tangent, to be updated
!
!-----------------------------------------------------------------------
! Originally coded 19/03/2003 KKO
! last modified
!=======================================================================
!-----------------------------------------------------------------------
subroutine UMAT_Overstress(stress, statev, ddsdde, sse, spd, scd, stran,    &
          dstran, time, dtime, cmname, ndi, nshr, ntens, nstatv, props, &
          nprops, coords, drott, pnewdt, celent, noel,                  &
          npt, layer, kspt, kstep, kinc)
  implicit none
!  include 'aba_param.inc'
! implicit real*8(a-h,o-z) !aba_param.inc
  integer nprecd           !aba_param.inc
  parameter (nprecd=2)
  Integer, Parameter :: dp = Selected_real_kind(12,50)
  character*8 cmname
  Real(dp) :: stress, statev, ddsdde, ddsddt, drplde, stran, dstran,   &
      predef, dpred, props, coords, drott, dfgrd0, dfgrd1,sse,spd,scd, &
      rpl,drpldt,time,dtime,temp,dtemp,pnewdt,celent
	Integer ntens,nstatv,nprops,noel,npt,layer,kspt,kstep,kinc,ndi,      &
	    nshr
  dimension stress(ntens), statev(nstatv), ddsdde(ntens, ntens),       &
      ddsddt(ntens), drplde(ntens), stran(ntens), dstran(ntens),       &
      predef(1), dpred(1), props(nprops), coords(3), drott(3, 3),      &
      dfgrd0(3, 3), dfgrd1(3, 3)
!- local variables -----------------------------------------------------
  Real(dp) :: EMOD,ENU,EBULK3,EG2,EG,ELAM,ONE,TWO,THREE
  Real(dp) :: Stress_Old(ntens),YIELDF,F_Tilde,Dev_Stress(ntens),SM, &
              ALPHA,A_const,B_const, DF_DSigma(ntens),               &
              DDF_DDSigma(ntens,ntens),Lambda,Eta,H(ntens,ntens),    &
              N_Power,CINV(ntens,ntens),Resid,HN(ntens)
  Real(dp) :: TOLER,DT,DLambda_DF,DPhi_DF,DKappa_DLambda,DF_DKappa,  &
              LambdaDot,DF_DLambda,DDF_DSigmaDLambda(ntens),         &
              DF_DLambdaDot,DDF_DSigmaDLambdaDot(ntens),             &
              Statev_Old(nstatv),DPHIDF,AA
  Character*80 :: string,form,CFLAG,DumCHR
	Integer K1,K2,II,JJ,ITER,MAXITER,ICONVE, Hardening_type,PRNTYP,ELN,IPN
  Integer, Parameter :: DATF=6,LOGF=5,MSGF=7,TESTFIL1=105,TESTFIL2=106
  Parameter(ONE=1.0_dp,TWO=2.0_dp,THREE=3.0_dp)

  Real(dp) :: SLEN,BEETA,SIGMA_T,SIGMA_C,J2INVARIANT,HARDENING_K,DUM


  If (NDI.NE.3) then
    Write(6,1)
 1  Format(//,30X,'***ERROR - THIS UMAT MAY ONLY BE USED FOR ',&
               'ELEMENTS WITH THREE DIRECT STRESS COMPONENTS')
  Endif

  N_Power = PROPS(5)
  Eta     = PROPS(6)
  A_Const = PROPS(7)
  B_Const = PROPS(8)
  PRNTYP  = INT(PROPS(12))
  TOLER   = PROPS(13)
  ELN     = INT(PROPS(15))
  IPN     = INT(PROPS(16))
  !MTYP    = INT(PROPS(17))
  MAXITER = INT(PROPS(14))
  IF(MAXITER < 2) MAXITER = 2
  DT      = DTIME
  IF(PROPS(11) < 0.0_dp ) DT = -PROPS(11)
  CFLAG = 'Drucker-Prager_Overstress_A'

! Check Hardening type
  If(INT(PROPS(18)) == 1) Then
    Hardening_type = 1
  ElseIf(INT(PROPS(18)) == 2) Then
    Hardening_type = 2
  Else
    Write(DumCHR,*) PROPS(18)
    string = ' Unknown Hardening type: No = '//Trim(Adjustl(DumCHR))
    form = '(/,'//''' *** ERROR in UMAT_Vers01 *** '''//',/,A)'
    Call WriteOut(string,form,4,KSTEP,KINC,NOEL,NPT)
    Call xit
  Endif

!     ELASTIC PROPERTIES

  EMOD=PROPS(9)
  ENU=PROPS(10)
  IF(ENU.GT.0.4999.AND.ENU.LT.0.5001) ENU=0.499
  EBULK3=EMOD/(ONE-TWO*ENU)
  EG2=EMOD/(ONE+ENU)
  EG=EG2/TWO
  ELAM=(EBULK3-EG2)/THREE

!     ELASTIC STIFFNESS

  Do K1=1,NTENS
    Do K2=1,NTENS
       DDSDDE(K2,K1)=0.0
    Enddo
  Enddo

  Do K1=1,NDI
    Do K2=1,NDI
       DDSDDE(K2,K1)=ELAM
    Enddo
    DDSDDE(K1,K1)=EG2+ELAM
  Enddo

  Do K1=NDI+1,NTENS
    DDSDDE(K1,K1)=EG
  Enddo

  Stress_Old = stress
  stress = stress + MATMUL(DDSDDE,dstran)

! Check if Yielding or not?
  F_Tilde = YIELDF(NTENS, STRESS, PROPS, STATEV, CFLAG,Hardening_type)

  If(F_Tilde <= 0.0_dp) Then
!   tee elastic temput
    RETURN
  Endif


! Print if PRNTYP == 1

  IF(PRNTYP == 1 .and. noel == ELN) Then
    string = ' Subroutine UMAT_Overstress, material type: '//Trim(cmname)
    form = '(A)'
    Call WriteOut(string,form,4,KSTEP,KINC,NOEL,NPT)
    Write(DumCHR,'(1P,E13.4)') F_Tilde
    string = ' '
    form = '('''//' F_Tilde = '//Trim(Adjustl(DumCHR))//''','//' A)'
    Call WriteOut(string,form,4,0,KINC,NOEL,NPT)
  Endif
  IF(PRNTYP /= 0  .and. NPT == IPN .and. noel == ELN) &
     Write(MSGF,'(A,1P,4E12.3)') '** DT, DTIME',DT, DTIME
  IF(PRNTYP /=0  .and. NPT == IPN .and. noel == ELN) &
    WRITE(MSGF,'(2A)') &
      'ITER,F           Lambda      Kappa,      K,',&
             '          DLambda_DF, slen        R_norm'

! Local Iteration ( Wang et al )

  Lambda  = 0.0_dp
  ICONVE  = 0
  Statev_Old = Statev
  Call invert( DDSDDE, CINV, ntens )

! > > > > > > > > > > > > > > > > > > > >

  Do ITER = 1, MAXITER

    F_Tilde = YIELDF(NTENS,STRESS,PROPS,STATEV,CFLAG,Hardening_type)
    Resid = (F_Tilde / A_Const+B_Const)**N_Power -Lambda/Eta/DT
    If(Abs(Resid) < TOLER) Then
      ICONVE = 1
      EXIT
    Endif

    Call DEVSTRESS(NTENS, STRESS,SM,DEV_STRESS)
    Call DF_Derivatives(NTENS, CFLAG,Hardening_type,STRESS,DEV_STRESS,&
                          PROPS,                                      &
                          STATEV,LambdaDot,DF_DLambda,DF_DLambdaDot,  &
                          DDF_DDSigma,DDF_DSigmaDLambda,              &
                          DDF_DSigmaDLambdaDot,Alpha)
    H = CINV + Lambda*DDF_DDSigma
    Call invert( H, H, ntens )

    Call DFDSIGMA(NTENS,STRESS,DEV_STRESS,PROPS,STATEV,CFLAG,         &
                Hardening_type,DF_DSIGMA)
    HN = MATMUL(H, Lambda*DDF_DSigmaDLambda + DF_DSigma)
    DPhi_DF= DPHIDF(F_Tilde,PROPS,STRESS,NTENS, STATEV, CFLAG,        &
                    Hardening_type)
    AA = DPhi_DF*Dot_Product(DF_DSigma,HN) +1.0_dp/Eta/DT             &
        -DPhi_DF*DF_DLambda
    Lambda = Lambda + Resid/AA
    Stress = Stress_Old + Matmul(DDSDDE,dstran-Lambda*DF_DSigma)
    Statev(1) = Statev_Old(1) + (1+Alpha)*Lambda

!   Printing:

    IF(PRNTYP /= 0  .and. NPT == IPN .and. noel == ELN ) Then
      SLEN    = SQRT(2.0_dp * J2INVARIANT(NTENS, DEV_STRESS) )
      WRITE(MSGF,'(1P,I2,1x,9E12.4)') ITER,F_Tilde,Lambda,Statev(1),&
             HARDENING_K(PROPS, STATEV, Hardening_type),1.0/DF_DLambda,&
             SLEN,Resid
    EndIf

  Enddo
! < < < < < < < < < < < < < < < < < < < <

!   Printing:

    IF(PRNTYP /= 0  .and. NPT == IPN .and. noel == ELN) Then
      SLEN    = SQRT(2.0_dp * J2INVARIANT(NTENS, DEV_STRESS) )
      WRITE(MSGF,'(1P,I2,1x,9E12.4)') ITER,F_Tilde,Lambda,Statev(1),&
             HARDENING_K(PROPS, STATEV, Hardening_type),1.0/DF_DLambda,&
             SLEN,Resid
    EndIf
    IF(PRNTYP == 2  .and. NPT == IPN .and. noel == ELN) Then
      WRITE(MSGF,'(/,A)') 'ITER,F,Lambda,Kappa,K,DLambda_DF,0.0, &
                        Dot_Product(Resid,Resid),DPhi_DF,DKappa_DLambda'
      WRITE(MSGF,'(1P,I2,1x,13E12.4)') ITER,F_Tilde,Lambda,Statev(1), &
                           HARDENING_K(PROPS, STATEV, Hardening_type),&
                           1.0/DF_DLambda,slen,Resid,DPhi_DF,1.0+alpha
      WRITE(MSGF,'(/,A,2I3)') ' STRESS, KINC,ITER',KINC,ITER
      WRITE(MSGF,'(1P,6E13.4)') (STRESS(ii),ii=1,6)
      WRITE(MSGF,'(/,A,2I3)') ' Dev_STRESS, KINC,ITER',KINC,ITER
      WRITE(MSGF,'(1P,6E13.4)') (Dev_STRESS(ii),ii=1,6)
      WRITE(MSGF,'(/,A,2I3)') ' DF_DSigma, KINC,ITER',KINC,ITER
      WRITE(MSGF,'(1P,6E13.4)') (DF_DSigma(ii),ii=1,6)
      WRITE(MSGF,'(/,A,2I3)') ' DDF_DDSigma, KINC,ITER',KINC,ITER
      WRITE(MSGF,'(1P,6E13.4)') ((DDF_DDSigma(ii,jj),ii=1,6),jj=1,6)
    ENDIF

    IF(PRNTYP /= 0.and. NPT==IPN .and. noel==ELN .and. ICONVE == 1) Then
      SIGMA_T = PROPS(1)
      SIGMA_C = PROPS(2)
      Call DEVSTRESS(NTENS, STRESS,SM,DEV_STRESS)
      ALPHA   = (SIGMA_C - SIGMA_T) / (SIGMA_C + SIGMA_T)
      BEETA   = (1.0_dp + ALPHA) * SIGMA_T
      SLEN    = SQRT(2.0_dp * J2INVARIANT(NTENS, DEV_STRESS) )

      WRITE(TESTFIL2,'(1P,6E13.4)') &
                    3.0*SM,SQRT(1.5_dp)*SLEN,BEETA-ALPHA*3*SM,&
                    HARDENING_K(PROPS, STATEV, Hardening_type),Statev(1)
    ENDIF


  IF (ICONVE == 1) Then

!   "Algorithmic Tangent"

    If(ITER > 1 ) Then
      Stress_Old = MATMUL(H,DF_DSigma) ! Stress_Old used as dummy vector
      AA = DPhi_DF/AA                  ! Now: HN = Stress_Old
      DO II=1,ntens
        DO JJ=1,ntens
          DDSDDE(JJ,II) = H(JJ,II) - AA*HN(JJ)*Stress_Old(II)
        ENDDO
      ENDDO
    Endif

!   Suggestion for new increment length
    DUM = ITER*1.0/(1.0*MAXITER)
    IF(DUM>0.8) PNEWDT = 0.8D0

  Else
    write(DumCHR,'(A,1P,E13.4,A,E13.4)') ' ABS(F)',ABS(F_Tilde),' > ',TOLER
    string =' Material type: '//Trim(cmname)//' '//Trim(Adjustl(DumCHR))
    form = '(/,'//''' *** ERROR in UMAT_Vers01 *** '''//',/,A)'
    Call WriteOut(string,form,4,KSTEP,KINC,NOEL,NPT)
    !Call xit
    PNEWDT = 0.75D0
  Endif

End subroutine UMAT_Overstress
!***********************************************************************
!=======================================================================
! Subroutine UMAT_Consistency
!
! ABAQUS User Material Subroutine For Drucker Prager (Perzyna Consistency
! model)
!-----------------------------------------------------------------------
!
! input parameters
! ----------------
!  PROPS(1)  - Sigma_t
!  PROPS(2)  - Sigma_c
!  PROPS(3)  - K_Zero
!  PROPS(4)  - C_const
!  PROPS(5)  - N_Power
!  PROPS(6)  - Eta
!  PROPS(7)  - A_Const
!  PROPS(8)  - B_Const
!  PROPS(9)  - E
!  PROPS(10) - NU
!  PROPS(11) - DT (If real DTIME = 0, gien as negative)
!  PROPS(12) - 1.0       ! Print All to DAT
!  PROPS(13) - FTOLER
!  PROPS(14) - MAXITER
!  PROPS(15) - Element number for printing
!  PROPS(16) - Integration point number for printing
!  PROPS(17) - Plasticity type
!         1  - Perzyna Drucker-Prager (overstress)
!         2  - Drucker-Prager Consistency model for Perzyna's overstress
!              modell
!  PROPS(18) - Hardening type
!         1  - Linear hardening
!         2  - Exponential hardening
!
! input/output parameters
! ----------------
!
!  STATEV(1) - Equivalent plastic strain
!  stress    - Stress vector, to be updated
!  ddsdde    - Consistent tangent, to be updated
!
!-----------------------------------------------------------------------
! Originally coded 19/03/2003 KKO
! last modified
!=======================================================================
!-----------------------------------------------------------------------
subroutine UMAT_Consistency(stress, statev, ddsdde, sse, spd, scd, stran,&
          dstran, time, dtime, cmname, ndi, nshr, ntens, nstatv, props,  &
          nprops, coords, drott, pnewdt, celent, noel,                   &
          npt, layer, kspt, kstep, kinc)
  implicit none
!  include 'aba_param.inc'
! implicit real*8(a-h,o-z) !aba_param.inc
  integer nprecd           !aba_param.inc
  parameter (nprecd=2)
  Integer, Parameter :: dp = Selected_real_kind(12,50)
  character*8 cmname
  Real(dp) :: stress, statev, ddsdde, ddsddt, drplde, stran, dstran,   &
      predef, dpred, props, coords, drott, dfgrd0, dfgrd1,sse,spd,scd, &
      rpl,drpldt,time,dtime,temp,dtemp,pnewdt,celent
	Integer ntens,nstatv,nprops,noel,npt,layer,kspt,kstep,kinc,ndi,      &
	    nshr
  dimension stress(ntens), statev(nstatv), ddsdde(ntens, ntens),       &
      ddsddt(ntens), drplde(ntens), stran(ntens), dstran(ntens),       &
      predef(1), dpred(1), props(nprops), coords(3), drott(3, 3),      &
      dfgrd0(3, 3), dfgrd1(3, 3)
!- local variables -----------------------------------------------------
  Real(dp) :: EMOD,ENU,EBULK3,EG2,EG,ELAM,ONE,TWO,THREE
  Real(dp) :: Stress_Old(ntens),YIELDF,F_Tilde,Dev_Stress(ntens),SM, &
              ALPHA,A_const,B_const, DF_DSigma(ntens),               &
              DDF_DDSigma(ntens,ntens),Lambda,Eta,H(ntens,ntens),    &
              N_Power,CINV(ntens,ntens),Resid,HN(ntens)
  Real(dp) :: TOLER,DT,DLambda_DF,DPhi_DF,DKappa_DLambda,DF_DKappa,  &
              LambdaDot,DF_DLambda,DDF_DSigmaDLambda(ntens),         &
              DF_DLambdaDot,DDF_DSigmaDLambdaDot(ntens),             &
              Statev_Old(nstatv),DPHIDF,AA,Lambda_Dot
  Character*80 :: string,form,CFLAG,DumCHR
	Integer K1,K2,II,JJ,ITER,MAXITER,ICONVE, Hardening_type,PRNTYP,ELN,IPN
  Integer, Parameter :: DATF=6,LOGF=5,MSGF=7,TESTFIL1=105,TESTFIL2=106
  Parameter(ONE=1.0_dp,TWO=2.0_dp,THREE=3.0_dp)

  Real(dp) :: SLEN,BEETA,SIGMA_T,SIGMA_C,J2INVARIANT,HARDENING_K,DUM


  If (NDI.NE.3) then
    Write(6,1)
 1  Format(//,30X,'***ERROR - THIS UMAT MAY ONLY BE USED FOR ',&
               'ELEMENTS WITH THREE DIRECT STRESS COMPONENTS')
  Endif


  N_Power = PROPS(5)
  Eta     = PROPS(6)
  A_Const = PROPS(7)
  B_Const = PROPS(8)
  PRNTYP  = INT(PROPS(12))
  TOLER   = PROPS(13)
  MAXITER = INT(PROPS(14))
  ELN     = INT(PROPS(15))
  IPN     = INT(PROPS(16))
  IF(MAXITER < 2) MAXITER = 2
  DT      = DTIME
  IF(PROPS(11) < 0.0_dp ) DT = -PROPS(11)
  CFLAG = 'Drucker-Prager_Consistency_A'

! Check Hardening type

  If(INT(PROPS(18)) == 1) Then
    Hardening_type = 1
  ElseIf(INT(PROPS(18)) == 2) Then
    Hardening_type = 2
  Else
    Write(DumCHR,*) PROPS(18)
    string = ' Unknown Hardening type: No = '//Trim(Adjustl(DumCHR))
    form = '(/,'//''' *** ERROR in UMAT_Vers01 *** '''//',/,A)'
    Call WriteOut(string,form,4,KSTEP,KINC,NOEL,NPT)
    Call xit
  Endif

!     ELASTIC PROPERTIES

  EMOD=PROPS(9)
  ENU=PROPS(10)
  IF(ENU.GT.0.4999.AND.ENU.LT.0.5001) ENU=0.499
  EBULK3=EMOD/(ONE-TWO*ENU)
  EG2=EMOD/(ONE+ENU)
  EG=EG2/TWO
  ELAM=(EBULK3-EG2)/THREE

!     ELASTIC STIFFNESS

  Do K1=1,NTENS
    Do K2=1,NTENS
       DDSDDE(K2,K1)=0.0
    Enddo
  Enddo

  Do K1=1,NDI
    Do K2=1,NDI
       DDSDDE(K2,K1)=ELAM
    Enddo
    DDSDDE(K1,K1)=EG2+ELAM
  Enddo

  Do K1=NDI+1,NTENS
    DDSDDE(K1,K1)=EG
  Enddo

  Stress_Old = stress
  stress = stress + MATMUL(DDSDDE,dstran)

! Check if Yielding or not?
  F_Tilde = YIELDF(NTENS, STRESS, PROPS, STATEV, CFLAG,Hardening_type)

  If(F_Tilde <= 0.0_dp) Then
!   tee elastic temput
    RETURN
  Endif

  Lambda    = 0.0_dp
  Lambda_Dot= 0.0_dp
  Statev_Old = Statev
  ICONVE = 0
  Call invert( DDSDDE, CINV, ntens )

! Print if PRNTYP == 1

  IF(PRNTYP == 1 .and. noel == ELN) Then
    string = ' Subroutine UMAT_Consistency, material type: '//Trim(cmname)
    form = '(A)'
    Call WriteOut(string,form,4,KSTEP,KINC,NOEL,NPT)
    Write(DumCHR,'(1P,E13.4)') F_Tilde
    string = ' '
    form = '('''//' F_Tilde = '//Trim(Adjustl(DumCHR))//''','//' A)'
    Call WriteOut(string,form,4,0,KINC,NOEL,NPT)
  Endif
  IF(PRNTYP /= 0  .and. NPT == IPN .and. noel == ELN) &
     Write(MSGF,'(A,1P,4E12.3)') '** DT, DTIME',DT, DTIME
  IF(PRNTYP /=0  .and. NPT == IPN .and. noel == ELN) &
    WRITE(MSGF,'(2A)') &
      'ITER,F           Lambda      Kappa,      K,',&
             '          DLambda_DF, slen        R_norm'

! > > > > > > > > > > > > > > > > > > > >

! Local Iteration ( Wang et al )

  Do ITER = 1, MAXITER

    Lambda_Dot = Lambda/DT
    F_Tilde = YIELDF(NTENS,STRESS,PROPS,STATEV,CFLAG,Hardening_type) &
              -  A_Const*(B_Const+(Lambda_Dot/Eta)**(1.0_dp/N_Power))
    Resid = F_Tilde
    If(Abs(Resid) < TOLER) Then
      ICONVE = 1
      EXIT
    Endif

    Call DEVSTRESS(NTENS, STRESS,SM,DEV_STRESS)
    Call DF_Derivatives(NTENS, CFLAG,Hardening_type,STRESS,DEV_STRESS,&
                          PROPS,                                      &
                          STATEV,LambdaDot,DF_DLambda,DF_DLambdaDot,  &
                          DDF_DDSigma,DDF_DSigmaDLambda,              &
                          DDF_DSigmaDLambdaDot,Alpha)
    H = CINV + Lambda*DDF_DDSigma
    Call invert( H, H, ntens )

    Call DFDSIGMA(NTENS,STRESS,DEV_STRESS,PROPS,STATEV,CFLAG,         &
                Hardening_type,DF_DSIGMA)
    HN = MATMUL(H, Lambda*DDF_DSigmaDLambda + DF_DSigma               &
                + Lambda/DT*DDF_DSigmaDLambdaDot)
    AA = Dot_Product(DF_DSigma,HN) - DF_DLambda - DF_DLambdaDot/DT
    Lambda = Lambda + Resid/AA
    Stress = Stress_Old + Matmul(DDSDDE,dstran-Lambda*DF_DSigma)
    Statev(1) = Statev_Old(1) + (1+Alpha)*Lambda

!   Printing:

    IF(PRNTYP /= 0  .and. NPT == IPN .and. noel == ELN ) Then
      SLEN    = SQRT(2.0_dp * J2INVARIANT(NTENS, DEV_STRESS) )
      WRITE(MSGF,'(1P,I2,1x,9E12.4)') ITER,F_Tilde,Lambda,Statev(1),&
             HARDENING_K(PROPS, STATEV, Hardening_type),1.0/DF_DLambda,&
             SLEN,Resid
    EndIf

  Enddo

! < < < < < < < < < < < < < < < < < < < <

!   Printing:

    IF(PRNTYP /= 0  .and. NPT == IPN .and. noel == ELN) Then
      SLEN    = SQRT(2.0_dp * J2INVARIANT(NTENS, DEV_STRESS) )
      WRITE(MSGF,'(1P,I2,1x,9E12.4)') ITER,F_Tilde,Lambda,Statev(1),&
             HARDENING_K(PROPS, STATEV, Hardening_type),1.0/DF_DLambda,&
             SLEN,Resid
    EndIf
    IF(PRNTYP == 2  .and. NPT == IPN .and. noel == ELN) Then
      WRITE(MSGF,'(/,A)') 'ITER,F,Lambda,Kappa,K,DLambda_DF,0.0, &
                        Dot_Product(Resid,Resid),DF_DLambdaDot,DKappa_DLambda'
      WRITE(MSGF,'(1P,I2,1x,13E12.4)') ITER,F_Tilde,Lambda,Statev(1), &
                           HARDENING_K(PROPS, STATEV, Hardening_type),&
                           1.0/DF_DLambda,slen,Resid,DF_DLambdaDot,1.0+alpha
      WRITE(MSGF,'(/,A,2I3)') ' STRESS, KINC,ITER',KINC,ITER
      WRITE(MSGF,'(1P,6E13.4)') (STRESS(ii),ii=1,6)
      WRITE(MSGF,'(/,A,2I3)') ' Dev_STRESS, KINC,ITER',KINC,ITER
      WRITE(MSGF,'(1P,6E13.4)') (Dev_STRESS(ii),ii=1,6)
      WRITE(MSGF,'(/,A,2I3)') ' DF_DSigma, KINC,ITER',KINC,ITER
      WRITE(MSGF,'(1P,6E13.4)') (DF_DSigma(ii),ii=1,6)
      WRITE(MSGF,'(/,A,2I3)') ' DDF_DDSigma, KINC,ITER',KINC,ITER
      WRITE(MSGF,'(1P,6E13.4)') ((DDF_DDSigma(ii,jj),ii=1,6),jj=1,6)
    ENDIF

    IF(PRNTYP /= 0.and. NPT==IPN .and. noel==ELN .and. ICONVE == 1) Then
      SIGMA_T = PROPS(1)
      SIGMA_C = PROPS(2)
      Call DEVSTRESS(NTENS, STRESS,SM,DEV_STRESS)
      ALPHA   = (SIGMA_C - SIGMA_T) / (SIGMA_C + SIGMA_T)
      BEETA   = (1.0_dp + ALPHA) * SIGMA_T
      SLEN    = SQRT(2.0_dp * J2INVARIANT(NTENS, DEV_STRESS) )

      WRITE(TESTFIL2,'(1P,6E13.4)') &
                    3.0*SM,SQRT(1.5_dp)*SLEN,BEETA-ALPHA*3*SM,&
                    HARDENING_K(PROPS, STATEV, Hardening_type),Statev(1)
    ENDIF

  IF (ICONVE == 1) Then

!   "Algorithmic Tangent"

    If(ITER > 1 ) Then
      Stress_Old = MATMUL(H,DF_DSigma) ! Stress_Old used as dummy vector
      AA = 1.0_dp/AA                   ! Now: HN = Stress_Old
      DO II=1,ntens
        DO JJ=1,ntens
          DDSDDE(JJ,II) = H(JJ,II) - AA*HN(JJ)*Stress_Old(II)
        ENDDO
      ENDDO
    Endif

!   Suggestion for new increment length
    DUM = ITER*1.0/(1.0*MAXITER)
    IF(DUM>0.8) PNEWDT = 0.8D0

  Else
    write(DumCHR,'(A,1P,E13.4,A,E13.4)') ' ABS(F)',ABS(F_Tilde),' > ',TOLER
    string =' Material type: '//Trim(cmname)//' '//Trim(Adjustl(DumCHR))
    form = '(/,'//''' *** ERROR in UMAT_Vers01 *** '''//',/,A)'
    Call WriteOut(string,form,4,KSTEP,KINC,NOEL,NPT)
    !Call xit
    PNEWDT = 0.75D0
  Endif

End subroutine UMAT_Consistency
