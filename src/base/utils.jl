@mtkmodel Inlet begin
    @structural_parameters begin
        ms = missing
    end

    @components begin
        f = MaterialFlow(ms=ms)
    end

    @variables begin
        p(t),               [description="inlet pressure"]      #, unit=u"Pa"]
        T(t),               [description="inlet temperature"]   #, unit=u"K"]
        n(t),               [description="inlet molar flow"]    #, unit=u"mol s^-1"]
        m(t),               [description="inlet mass flow"]     #, unit=u"kg s^-1"]
        xᵢ(t)[1:ms.N_c],    [description="inlet mole fractions"]#, unit=u"mol mol^-1"] 
    end

    @equations begin
        f.s.T ~ T
        f.s.ϱ ~ ms.molar_density(p,T,xᵢ)
        scalarize(f.s.nᵢ .~ n .* xᵢ)
        f.c.p ~ p
        f.c.h ~ ms.VT_enthalpy(f.s.ϱ,T,xᵢ)
        f.c.xᵢ ~ xᵢ
        n ~ m / (ms.Mw' * xᵢ)
        f.c.n ~ n
    end
end