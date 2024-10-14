@mtkmodel SimpleIsobaricHeatExchanger begin
    @structural_parameters begin
        ms = missing
    end

    @components begin
        cv = SimpleControlVolume(ms=ms,N_flows=2,N_heats=1)
    end

    @equations begin
        cv.f1.c.p ~ cv.f2.c.p
    end
end