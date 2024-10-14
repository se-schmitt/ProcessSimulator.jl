using ProcessSimulator
using ModelingToolkit, DifferentialEquations
using ModelingToolkit: t_nounits as t

const PS = ProcessSimulator

# Create material sources
Mw = 0.004          # kg/mol
R = 2.1e3*Mw        # J/(mol K)
cₚ = 5.2e3*Mw       # J/(mol K)
cᵥ = cₚ - R

matsource = PS.MaterialSource(
    "helium",
    ["helium"],[0.004],
    (ϱ,T,n) -> ϱ*R*T, 
    (p,T,n;kwargs...) -> p/(R*Mw*T),
    (ϱ,T,n) -> cᵥ*T,
    (ϱ,T,n) -> cₚ*T,
    (ϱ,T,n) -> cᵥ*log(T) + R*log(1/ϱ),
    (p,T,n) -> NaN
)

# Create flowsbheet
# Create flowsbheet
@mtkmodel Flowsheet begin
    @structural_parameters begin
        ms = missing
    end

    @components begin
        comp_12 = PS.SimpleAdiabaticCompressor(ms=ms)
        heat_22⁺ = PS.SimpleIsobaricHeatExchanger(ms=ms)  
        heat_2⁺3 = PS.SimpleIsobaricHeatExchanger(ms=ms)
        turb_34 = PS.SimpleAdiabaticCompressor(ms=ms)
        heat_44⁺ = PS.SimpleIsobaricHeatExchanger(ms=ms)
        heat_4⁺1 = PS.SimpleIsobaricHeatExchanger(ms=ms)
    end

    @equations begin
        connect(comp_12.cv.f2.c,  heat_22⁺.cv.f1.c)
        connect(heat_22⁺.cv.f2.c, heat_2⁺3.cv.f1.c)
        connect(heat_2⁺3.cv.f2.c, turb_34.cv.f1.c)
        connect(turb_34.cv.f2.c,  heat_44⁺.cv.f1.c)
        connect(heat_44⁺.cv.f2.c, heat_4⁺1.cv.f1.c)
        connect(heat_4⁺1.cv.f2.c, comp_12.cv.f1.c)
        heat_22⁺.cv.q1.Q ~ -heat_4⁺1.cv.q1.Q
    end
end

@named flowsheet_ = Flowsheet(ms=matsource)

inp_str = [
    "(comp_12₊cv₊f1₊s₊nᵢ(t))[1]",
    # "comp_12₊cv₊f1₊c₊n(t)",
    "comp_12₊cv₊f1₊s₊T(t)",
    "comp_12₊cv₊f1₊c₊p(t)",
    "comp_12₊cv₊f2₊c₊p(t)",
    "heat_22⁺₊cv₊f2₊s₊T(t)",
    "turb_34₊cv₊f1₊s₊T(t)",
    "turb_34₊cv₊f2₊c₊p(t)",
    "heat_44⁺₊cv₊f2₊s₊T(t)",
]
out_str = [
    "heat_44⁺₊cv₊q1₊Q(t)",
    "comp_12₊cv₊w1₊W(t)",
    "turb_34₊cv₊w1₊W(t)"
]

unk = unknowns(flowsheet_)
inp = unk[[findfirst(str .== string.(unk)) for str in inp_str]]
out = unk[[findfirst(str .== string.(unk)) for str in out_str]]

flowsheet = structural_simplify(flowsheet_,(inp,out))

inp = [
    (comp_12.cv.f1.c.xᵢ)[1] => 1.0,     # T1
    comp_12.cv.f1.c.n => 1.0,           # T1
    comp_12.cv.f1.s.T => 298.15,        # T1
    comp_12.cv.f1.c.p => 1e5,           # p1
    comp_12.cv.f2.c.p => 1e6,           # p2
    heat_22⁺.cv.f2.s.T => 298.15,       # T2⁺ 
    turb_34.cv.f1.s.T => 253.15,        # T3
    turb_34.cv.f2.c.p => 1e5,           # p4
    heat_44⁺.cv.f2.s.T => 253.15,       # T4⁺
    comp_12.ηᴱ => 1.0,                  # ηᴱ
    turb_34.ηᴱ => 1.0                   # ηᴱ
]


# @mtkbuild flowsheet = Flowsheet() (first.(giv),unk)        # or: 'io=(first.(giv),unk)'?

@named flowsheet = Flowsheet()
structural_simplify!(flowsheet,(first.(giv),unk))

u0 = [
    comp_12.cv.s2.T => 1.0,
    comp_12.cv_s.s2.T => 1.0,
    comp_12.cv.w1.W => -1.0,
    turb_34.cv.s2.T => 1.0,
    turb_34.cv_s.s2.T => 1.0,
    turb_34.cv.w1.W => 1.0,
    heat_4⁺1.cv.s2.T => 1.0,
]

prob = SteadyStateProblem(flowsheet,giv,u0)
sol = solve(prob)

ε_KM = abs(sol[heat_44⁺.cv.Qs[1]])/abs(sol[turb_34.cv.Ws[1]] + sol[comp_12.cv.Ws[1]])

(T₁,T₂,T₂₊,T₃,T₄,T₄₊) = (
    sol.ps[comp_12.cv.s1.T],
    sol[comp_12.cv.s2.T],
    sol.ps[heat_22⁺.cv.s2.T],
    sol.ps[turb_34.cv.s1.T],
    sol[turb_34.cv.s2.T],
    sol.ps[heat_44⁺.cv.s2.T]
)

@test round(T₂,digits=2) ≈ 755.58
@test round(T₄,digits=2) ≈ 99.89
@test round(ε_KM,digits=3) ≈ 0.504