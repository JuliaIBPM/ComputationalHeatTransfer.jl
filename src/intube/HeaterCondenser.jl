export OHPConfiguration

function OHPConfiguration(configure_type::String,power::Real,Tc::Real,hc::Real,Δx::Real;hc2ratio=1/20)

    if configure_type == "ASETS-II OHP 1 LARGE HEATER"
        total_heater_area = 2.0inch*2.0inch;
        qe = power/total_heater_area;

        eb1 = Rectangle(0.5inch,1.0inch,1.5*Δx)
        Tfe = RigidTransform((0.7inch,-0.0),0.0)
        Tfe(eb1)

        eb2 = Rectangle(0.5inch,1.0inch,1.5*Δx)
        Tfe = RigidTransform((-0.7inch,-0.0),0.0)
        Tfe(eb2)

        cb1 = Rectangle(0.55inch,0.0648/2 ,1.5*Δx) # 0.02916 = 0.0648/2 
        Tfc = RigidTransform((-2.5inch,-0.0),0.0)
        Tfc(cb1)

        cb2 = Rectangle(0.55inch,0.0648/2 ,1.5*Δx)
        Tfc = RigidTransform((2.5inch,-0.0),0.0)
        Tfc(cb2)

        eparams1 = PrescribedHeatFluxRegion(qe,eb1);
        eparams2 = PrescribedHeatFluxRegion(qe,eb2);
        cparams1 = PrescribedHeatModelRegion(hc,Tc,cb1);
        cparams2 = PrescribedHeatModelRegion(hc,Tc,cb2);

    return [eparams1,eparams2], [cparams1,cparams2]
    end

    if (configure_type == "ASETS-II OHP 2 LARGE HEATER") || (configure_type == "ASETS-II OHP 3 LARGE HEATER")
        total_heater_area = 2.0inch*2.0inch;
        qe = power/total_heater_area;

        eb1 = Rectangle(0.5inch,1.0inch,1.5*Δx)
        Tfe = RigidTransform((0.7inch,-0.0),0.0)
        Tfe(eb1)

        eb2 = Rectangle(0.5inch,1.0inch,1.5*Δx)
        Tfe = RigidTransform((-0.7inch,-0.0),0.0)
        Tfe(eb2)

        cb1 = Rectangle(0.55inch,0.0648/2 ,1.5*Δx) # 0.02916 = 0.0648/2 
        Tfc = RigidTransform((-2.5inch,-0.0),0.0)
        Tfc(cb1)

        cb2 = Rectangle(0.55inch,0.0648/2,1.5*Δx)
        Tfc = RigidTransform((2.5inch,-0.0),0.0)
        Tfc(cb2)

        eparams1 = PrescribedHeatFluxRegion(qe,eb1);
        eparams2 = PrescribedHeatFluxRegion(qe,eb2);
        cparams1 = PrescribedHeatModelRegion(hc*hc2ratio,Tc,cb1);
        cparams2 = PrescribedHeatModelRegion(hc,Tc,cb2);

    return [eparams1,eparams2], [cparams1,cparams2]
    end

    if configure_type == "ASETS-II OHP 1 SMALL HEATER"
        total_heater_area = 0.5inch*0.5inch;
        qe = power/total_heater_area;

        eb1 = Rectangle(0.25inch,0.25inch,1.5*Δx)
        Tfe = RigidTransform((0.0inch,-0.0),0.0)
        Tfe(eb1)

        cb1 = Rectangle(0.55inch,0.0648/2,1.5*Δx)
        Tfc = RigidTransform((-2.5inch,-0.0),0.0)
        Tfc(cb1)

        cb2 = Rectangle(0.55inch,0.0648/2,1.5*Δx)
        Tfc = RigidTransform((2.5inch,-0.0),0.0)
        Tfc(cb2)

        eparams1 = PrescribedHeatFluxRegion(qe,eb1);
        cparams1 = PrescribedHeatModelRegion(hc,Tc,cb1);
        cparams2 = PrescribedHeatModelRegion(hc,Tc,cb2);

    return [eparams1], [cparams1,cparams2]
    end

    if (configure_type == "ASETS-II OHP 2 SMALL HEATER") || (configure_type == "ASETS-II OHP 3 SMALL HEATER")
        total_heater_area = 0.5inch*0.5inch;
        qe = power/total_heater_area;

        eb1 = Rectangle(0.25inch,0.25inch,1.5*Δx)
        Tfe = RigidTransform((0.0inch,-0.0),0.0)
        Tfe(eb1)

        cb1 = Rectangle(0.55inch,0.0648/2,1.5*Δx)
        Tfc = RigidTransform((-2.5inch,-0.0),0.0)
        Tfc(cb1)

        cb2 = Rectangle(0.55inch,0.0648/2,1.5*Δx)
        Tfc = RigidTransform((2.5inch,-0.0),0.0)
        Tfc(cb2)

        eparams1 = PrescribedHeatFluxRegion(qe,eb1);
        cparams1 = PrescribedHeatModelRegion(hc*hc2ratio,Tc,cb1);
        cparams2 = PrescribedHeatModelRegion(hc,Tc,cb2);

    return [eparams1], [cparams1,cparams2]
    end


    return "configuration not recognized"
end

