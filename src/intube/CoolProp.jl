using CoolProp
using Interpolations

export nondi_PtoT,nondi_TtoP,nondi_PtoD,nondi_DtoP

fluid_type = "Butane"

Tcrit = CoolProp.PropsSI("Tcrit",fluid_type);
Tmin = CoolProp.PropsSI("Tmin",fluid_type);

Trange = LinRange(Tmin, Tcrit, 1000)
Prange = CoolProp.PropsSI.("P","T",Trange,"Q",1.0,fluid_type);
Drange = CoolProp.PropsSI.("D","T",Trange,"Q",1.0,fluid_type);

T0=295.0;
D0=CoolProp.PropsSI("D","T",T0,"Q",1.0,fluid_type);
P0=CoolProp.PropsSI("P","T",T0,"Q",1.0,fluid_type);

nondi_Trange = Trange./T0
nondi_Prange = Prange./P0
nondi_Drange = Drange./D0;

nondi_PtoT = LinearInterpolation(nondi_Prange, nondi_Trange);
nondi_TtoP = LinearInterpolation(nondi_Trange, nondi_Prange);
nondi_PtoD = LinearInterpolation(nondi_Prange, nondi_Drange);
nondi_DtoP = LinearInterpolation(nondi_Drange, nondi_Prange);
# function nondi_PtoT(nondi_P;P0=103225,T0=273.15,quality=1.0,fluid_type="Butane")
#     CoolProp.PropsSI("T","P",nondi_P*P0,"Q",quality,fluid_type)/T0
# end
#
# function nondi_TtoP(nondi_T;T0=273.15,P0=103225,quality=1.0,fluid_type="Butane")
#     CoolProp.PropsSI("P","T",nondi_T*T0,"Q",quality,fluid_type)/P0
# end
#
# function nondi_PtoD(nondi_P;P0=103225,D0=2.75673,quality=1.0,fluid_type="Butane")
#     CoolProp.PropsSI("D","P",nondi_P*P0,"Q",quality,fluid_type)/D0
# end
#
# function nondi_DtoP(nondi_D;D0=2.75673,P0=103225,quality=1.0,fluid_type="Butane")
#     CoolProp.PropsSI("P","D",nondi_D*D0,"Q",quality,fluid_type)/P0
# end
