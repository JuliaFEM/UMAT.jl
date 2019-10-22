# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/UMAT.jl/blob/master/LICENSE

"""
!cmname    - DPC_LH1  = Drucker-Prager Linear hardening
!cmname    - DPC_EXP1 = Drucker-Prager Exponential hardening

! ------------------- ABAQUS:
!
  PROPS(1)  = 200.0d0   ! Sigma_t
  PROPS(2)  = 200.0d0   ! Sigma_c
  PROPS(3)  =  100.00d2   ! K_Zero
  PROPS(4)  =   100.0d0   ! C_const
  PROPS(5)  =   1.0     ! N_Power
  PROPS(6)  =   5.0d0    ! Eta
  PROPS(7)  =   200.0d0   ! A_Const
  PROPS(8)  =   0.0d0   ! B_Const
  PROPS(9)  =   2.0D5   ! E
  PROPS(10) =   0.3d0   ! NU
  PROPS(11) =   0.0d0   ! DT (jos todellinen DT=0)
  PROPS(12) = 1.0       ! Print All to DAT
  PROPS(13) = 1.e-13      ! FTOLER
  PROPS(14) = 20.0      ! MAXITER
  PROPS(15) = 1.0       ! ELNUM
  PROPS(16) = 1.0       ! IPNUM
  PROPS(17) = 2.0       ! MTYP
  PROPS(18) = 1.0       ! HARDTYP
  DTIME = .00011      ! Aika-askel DT
  !DTIME = 1.0        ! Aika-askel DT

  cmname  = 'DPC_LH1'
  !cmname  = 'DPO_EXP1'
  ndi = 3
  nshr = 3
  Stress = 0.0d0
  ddsdde = 0.0d0
  dstran = 0.0d0
  stran  = 0.0d0
  Stran_Old = 0.0D0
  dstran(3) =  5.0d-5                 ! Venymainkrementti (vakio)
  dstran(2) = -PROPS(10) * dstran(3)  ! Venymainkrementti (vakio)
  dstran(1) = -PROPS(10) * dstran(3)  ! Venymainkrementti (vakio)
"""
function DruckerPragerMaterial(;Sigma_t = 200.,
                               Sigma_c = 200.,
                               K_Zero = 100.,
                               C_const = 100.,
                               N_Power = 1.,
                               Eta = 5.,
                               A_Const = 200.,
                               B_Const = 0.,
                               E = 2.e5,
                               NU = 0.3,
                               DT = 0.,
                               PRNTYP = false, # remember to convert float
                               FTOLER = 1.e-13,
                               MAXITER = 20,   # remember to convert float
                               ELNUM = 1,      # remember to convert float
                               IPNUM = 1,      # remember to convert float
                               MTYP = 2,       # remember to convert float
                               HARDTYP = 1,    # remember to convert float
                               DPC="DPC_LH1")


    umat_other = UmatOtherState(NTENS=6, NSTATV=1, NPROPS=18, CMNAME=DPC)
    umat_params = [Sigma_t, Sigma_c, K_Zero, C_const, N_Power, Eta, A_Const,
                   B_Const, E, NU, DT, float(PRNTYP), FTOLER, float(MAXITER),
                   float(ELNUM), float(IPNUM), float(MTYP), float(HARDTYP)]

    NTENS=6
    NSTATV=1
    NPROPS=18
    return UmatMaterial(NTENS=NTENS,
                        NSTATV=NSTATV,
                        NPROPS=NPROPS,
                        #drivers=UmatDriverState(NTENS=NTENS),
                        #ddrivers=UmatDriverState(NTENS=NTENS),
                        parameters=UmatParameterState(NPROPS=NPROPS, PROPS=umat_params),
                        #variables=UmatVariableState(NTENS=NTENS, NSTATV=NSTATV),
                        #variables_new=UmatVariableState(NTENS=NTENS, NSTATV=NSTATV),
                        #dparameters=UmatParameterState(NPROPS=NPROPS),
                        umat_other=umat_other,
                        lib_path=joinpath(lib_dir,"libdrucker_prager_plasticity." * dlext))
end
