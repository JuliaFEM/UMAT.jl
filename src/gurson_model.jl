# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/UMAT.jl/blob/master/LICENSE


"""
!STATEV(1)=f:Porosity
!STATEV(2)=ep:equivalent plastic strain
!STATEV(3)=DEp: Hydrostatic inelastic strain increment
!STATEV(4)=DEp: Deviatoric inelastic strain increment
!STATEV(5)...STATEV(10)=EELAS
!STATEV(11)...STATEV(16)=EPLAS
!STATEV(17)...STATEV(22)=X

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
"""
function GursonMaterial(initial_porosity=0.05, critical_porosity=0.08, failure_porosity=0.1,
                        volume_fraction_of_graphite=0.1, sdev_of_graphite=0.001,
                        mean_of_graphite=0.0, isotropic_hardening_amplitude=100.0,
                        isotropic_hardening_exponent=1.0, cyclic_hardening_parameter=1.0,
                        cyclic_saturation_parameter=1.0, Youngs_modulus=169000.0,
                        Poissons_ratio=0.3, yield_stress=150.0)

    umat_params= [initial_porosity,critical_porosity,failure_porosity,volume_fraction_of_graphite,
                  sdev_of_graphite, mean_of_graphite, isotropic_hardening_amplitude,
                  isotropic_hardening_exponent, cyclic_hardening_parameter, cyclic_saturation_parameter,
                  Youngs_modulus, Poissons_ratio, yield_stress]

    return UmatMaterial(NTENS=6, NSTATV=22, NPROPS=13,
                        parameters=UmatParameterState(NPROPS=13, PROPS=umat_params),
                        lib_path=joinpath(lib_dir,"libgurson_porous_plasticity." * dlext))
end
