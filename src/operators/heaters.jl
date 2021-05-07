"""
    PrescribedHeatFluxRegion(q::Real,b::Body)

Set up a prescribed heat flux (per unit area) `q` in a region described by the
interior of body `b`.
"""
struct PrescribedHeatFluxRegion{BT}
    q :: Float64
    body :: BT
end


struct PrescribedHeatFlux{GT,MT,BT}
    q :: GT
    mask :: MT
    body :: BT
    cache :: GT
end

PrescribedHeatFlux(p::PrescribedHeatFluxRegion,u::GridData,g::PhysicalGrid) =
            (o = similar(u); fill!(o,1); PrescribedHeatFlux(p.q*o,Mask(p.body,g,u),p.body,o))

(e::PrescribedHeatFlux)() = e.mask(e.cache,e.q)

"""
    PrescribedHeatModelRegion(hc::Real,Tc::Real,b::Body)

Set up a prescribed heat flux model `hc*(Tc - T)` when acting on data `T` in a region
described by the interior of body `b`.
"""
struct PrescribedHeatModelRegion{BT}
    hc :: Float64
    Tc :: Float64
    body :: BT
end

struct PrescribedHeatModel{GT,MT,BT}
    hc :: Float64
    Tc :: GT
    mask :: MT
    body :: BT
    cache1 :: GT
    cache2 :: GT
end

PrescribedHeatModel(p::PrescribedHeatModelRegion,u::GridData,g::PhysicalGrid) =
            (o = similar(u); fill!(o,1); PrescribedHeatModel(p.hc,p.Tc*o,Mask(p.body,g,u),p.body,o,similar(o)))

(c::PrescribedHeatModel{GT})(T::GT) where {GT} = (c.cache1 .= c.hc.*(c.Tc.-T); c.mask(c.cache2,c.cache1))
