using Test, GasChromatographySimulator, GasChromatographyTools

"""
    test_parameters(nonuniform="temperature", control_mode="Pressure")

Function to construct the parameters for a test system with two options:
    
    * nonuniform ... which parameter should have a non-uninform profile. Selection "temperature" (non-uniform temperature along column), "diameter" (non-uniform diameter of the colum, not yet implemented), "film_thickness" (non-uniform film ticknes of the column, not yet implemented)
    * control_mode ... "Pressure" control of the inlet pressure directly, "Flow" control of the flow, inlet pressure is calculated
"""
function test_parameters(nonuniform="temperature", control_mode="Pressure")
	L = 10.0
    if nonuniform == "diamater"
    else
	    d = 0.25e-3
    end
    if nonuniform == "film_thickness"
    else
        df = 0.25e-6
    end
	sp = "SLB5ms"
	gas = "He"

    time_steps = [0.0, 60.0, 600.0, 300.0]
    temp_steps = [40.0, 40.0, 300.0, 300.0]
    if control_mode == "Pressure"
        Fpin_steps = [200.0, 200.0, 300.0, 300.0].*1000.0 .+ 101300
    else
        Fpin_steps = 1.0.*ones(length(time_steps))./(60e6)
    end
    pout_steps = [0.0, 0.0, 0.0, 0.0].*1000.0
    if nonuniform == "temperature"
        ΔT_steps = [20.0, 30.0, 50.0, 40.0]
    else
        ΔT_steps = zeros(length(time_steps))
    end
    x₀_steps = zeros(length(time_steps))
    L₀_steps = L.*ones(length(time_steps))
    α_steps = zeros(length(time_steps))
	
    db_path = string(@__DIR__, "/data")
    db_file = "Database_test.csv"

	#db_path = "/Users/janleppert/Documents/GitHub/GasChromatographyTools/test/data/" 
	#db_file = "Database_test.csv"
	first_alkane = "C8"
	last_alkane = "C15"


    col = GasChromatographySimulator.constructor_System(L, d, df, sp, gas)
    sub = GasChromatographySimulator.load_solute_database(db_path, db_file, sp, gas, [first_alkane, last_alkane], zeros(2), zeros(2))
    opt = GasChromatographySimulator.Options(control=control_mode)
	prog = GasChromatographySimulator.constructor_Program(time_steps,temp_steps, Fpin_steps, pout_steps, ΔT_steps, x₀_steps, L₀_steps, α_steps, opt.Tcontrol, L)
	
	par = GasChromatographySimulator.Parameters(col, prog, sub, opt)

    return par
end

testpar = test_parameters() 
Tst = 273.15
Tref = 150.0 + Tst

tMref = GasChromatographySimulator.holdup_time(Tref, testpar.prog.Fpin_steps[1], testpar.prog.pout_steps[1], testpar.col.L, testpar.col.a_d[1], testpar.col.gas; vis="Blumberg", control="Pressure")
tR_lock = 12*tMref
tR_tol = 1e-3

n = GasChromatographyTools.RT_locking(testpar, tR_lock, tR_tol, testpar.sub[end].name; opt_itp="linear")
testpar1 = GasChromatographyTools.stretched_program(n , testpar)

pl, sol = GasChromatographySimulator.simulate(testpar1)

@test abs(pl.tR[end] - tR_lock)/tR_lock<tR_tol