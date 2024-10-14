# Base components
@mtkmodel MaterialFlow begin
    @structural_parameters begin
        ms = missing
        N_c::Int = 1
        phase::String = "unknown"
    end
    
    @components begin
        s = TϱState(N_c=N_c,phase=phase)
        c = PHxConnector(N_c=N_c,phase=phase)
    end
    
    @equations begin
        c.p ~ ms.pressure(s.ϱ,s.T,c.xᵢ)
        c.h ~ ms.VT_enthalpy(s.ϱ,s.T,c.xᵢ)
        1.0 ~ sum([xᵢ for xᵢ in c.xᵢ])
        scalarize(c.xᵢ .~ s.nᵢ ./ c.n)
    end
end

@mtkmodel TϱState begin
    @structural_parameters begin
        N_c::Int = 1
        phase::String = "unknown"
    end
    
    @variables begin
        T(t),               [description="temperature", output=true]    #, unit=u"K"]
        ϱ(t),               [description="density", output=true]        #, unit=u"mol m^-3"]
        nᵢ(t)[1:N_c]=1/N_c, [description="mole fractions", output=true] #, unit=u"mol s^-1"]
    end
end

@connector PHxConnector begin 
    @structural_parameters begin 
        N_c::Int = 1
        phase::String = "unknown"
    end
    
    @variables begin
        p(t),               [description="pressure"]                        #, unit=u"Pa"]
        h(t),               [description="molar enthalpy", connect=Stream]  #, unit=u"J mol^-1"]
        xᵢ(t)[1:N_c]=1/N_c, [description="mole fractions", connect=Stream]  #, unit=u"mol mol^-1"]
        n(t),               [description="total molar flow", connect=Flow]  #, unit=u"mol s^-1"]
    end
end

@connector HeatConnector begin
    @variables begin
        Q(t)#,              [description="heat flux", connect=Flow]  #, unit=u"J s^-1"]
    end
end

@connector WorkConnector begin
    @variables begin
        W(t)#,              [description="power", connect=Flow]      #, unit=u"J s^1"]
    end
end

@mtkmodel SimpleControlVolume begin

    @structural_parameters begin
        ms = missing
        N_c::Int = 1
        N_flows::Int = 0
        N_heats::Int = 0
        N_works::Int = 0
        N_phases::Int = 1
        phases = ["unknown"]
        stationary::Bool = true
    end

    @components begin
        flows = [MaterialFlow(ms=ms,name=Symbol("f$i")) for i in 1:N_flows]
        works = [WorkConnector(name=Symbol("w$i")) for i in 1:N_works]
        heats = [HeatConnector(name=Symbol("q$i")) for i in 1:N_heats]
    end

    @variables begin
        ΔH(t), [description="enthalpy difference inlets/outlets", output=true]  #, unit=u"J s^-1"]
        ΔE(t), [description="added/removed heat or work", output=true]          #, unit=u"J s^-1"]
    end

    @equations begin
        ΔH ~ sum([f.c.h*f.c.n for f in flows])
        ΔE ~ (isempty(heats) ? 0.0 : sum([q.Q for q in heats])) + (isempty(works) ? 0.0 : sum([w.W for w in works]))
        # Energy balance
        (stationary ? 0.0 : D(U)) ~ ΔH + ΔE
        # Mole balance
        [(stationary ? 0.0 : D(sum(nᵢ[:,j]))) ~ sum([f.c.n*f.c.xᵢ[j] for f in flows]) for j in 1:N_c][:]
    end
end

# TODO: extending models with array components does not work
# @mtkmodel TPControlVolume begin

#     @extend SimpleControlVolume(;ms, N_c, N_flows, N_heats, N_works, N_phases, phases)

#     @variables begin
#         T(t),                       [description="temperature", output=true]          #, unit=u"K"]
#         p(t),                       [description="pressure", output=true]             #, unit=u"Pa"]
#         ϱ(t)[1:N_phases],           [description="density", output=true]              #, unit=u"mol m^-3"]
#         (nᵢ(t))[1:N_phases,1:N_c],  [description="molar holdup", output=true]         #, unit=u"mol"]
#         n(t),                       [description="total molar holdup", output=true]   #, unit=u"mol"]
#         U(t),                       [description="internal energy", output=true]      #, unit=u"J"]
#     end

#     # reactions
#     #TODO: concept for modeling reactions (model reation enthalpy implicitly in MateiralSource function for enthalpy, here only the reaction rates and stochiometry)
#     # @parameters/@component begin     
#     # end

#     @equations begin
#         # Mole balance
#         n ~ sum(nᵢ)
#         # Thermodynamic system properties
#         [ϱ[i] ~ ms.molar_density(p,T,nᵢ[i,:];phase=phases[i]) for i in 1:N_ph]
#         U ~ sum([ms.VT_internal_energy(ϱ[i],T,nᵢ[i,:]) for i in 1:N_ph])
#     end
# end