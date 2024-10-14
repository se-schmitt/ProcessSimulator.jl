@mtkmodel SimpleAdiabaticCompressor begin
    @structural_parameters begin
        ms = missing
    end

    @components begin
        cv = SimpleControlVolume(ms=ms,N_flows=2,N_works=1)
    end

    @variables begin
        # W(t),               [description="power", output=true] #, unit=u"J s^1"]    
        # T1(t),              [description="inlet temperature", output=true] #, unit=u"K"]
        # p1(t),              [description="inlet pressure", output=true] #, unit=u"Pa"]
        # T2(t),              [description="outlet temperature", output=true] #, unit=u"K"]
        # p2(t),              [description="outlet pressure", output=true] #, unit=u"Pa"]
        T2_s(t),            [description="isentropic temperature", output=true] #, unit=u"K"]
        ϱ2_s(t),            [description="isentropic density", output=true]     #, unit=u"mol m^-3"]
    end

    @parameters begin
        ηᴱ(t),              [description="efficiency", output=true] #, unit=u"1"]        
    end

    # Equations
    @equations begin
        # W ~ cv.w1.W
        # T1 ~ cv.f1.s.T 
        # T2 ~ cv.f2.s.T
        # p1 ~ cv.f1.c.p 
        # p2 ~ cv.f2.c.p
        ms.VT_entropy(cv.f1.s.ϱ,cv.f1.s.T,cv.f1.c.xᵢ) ~ ms.VT_entropy(cv.f2.s.ϱ,cv.f2.s.T,cv.f2.c.xᵢ)
        ϱ2_s ~ ms.molar_density(cv.f2.c.p,T2_s,cv.f2.c.xᵢ)
        ηᴱ ~ cv.w1.W / (ms.VT_enthalpy(ϱ2_s,T2_s,cv.f2.c.xᵢ) - ms.VT_enthalpy(cv.f1.s.ϱ,cv.f1.s.T,cv.f1.c.xᵢ))
    end
end