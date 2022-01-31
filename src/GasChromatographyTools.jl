module GasChromatographyTools

#using Reexport
using GasChromatographySimulator
using Interpolations
using Intervals
using Dierckx
using HypertextLiteral
using PlutoUI
using Plots

#----RT-locking----------------------------------------------------
# helper function for 'RT-locking' and 'stretching' of GC-program
function stretched_program(n::Float64, par::GasChromatographySimulator.Parameters)
	# stretch the temperature program in 'par' by a factor 'n'
	if isa(par, Array)==true
		error("Select an element of the array of GC-system parameters.")
	else
		new_tsteps = n.*par.prog.time_steps
		new_T_itp = GasChromatographySimulator.temperature_interpolation(new_tsteps, par.prog.temp_steps, par.prog.gf, par.sys.L)
		new_pin_itp = GasChromatographySimulator.pressure_interpolation(new_tsteps, par.prog.pin_steps)
		new_pout_itp = GasChromatographySimulator.pressure_interpolation(new_tsteps, par.prog.pout_steps)
		new_prog = GasChromatographySimulator.Program(new_tsteps, par.prog.temp_steps, par.prog.pin_steps, par.prog.pout_steps, par.prog.gf, par.prog.a_gf, new_T_itp, new_pin_itp, new_pout_itp)
		new_par = GasChromatographySimulator.Parameters(par.sys, new_prog, par.sub, par.opt)
		return new_par
	end	
end

function initial_n(n::Float64, tR_lock::Float64, ii::Int, par::GasChromatographySimulator.Parameters)
    par_n = stretched_program(n, par)
    sol_n = GasChromatographySimulator.solve_system_multithreads(par_n)
    tR_n = sol_n[ii].u[end][1]
    if n>1
        while tR_n-tR_lock<0 && n<130.0
            n = n*2.0
            par_n = stretched_program(n, par)
            sol_n = GasChromatographySimulator.solve_system_multithreads(par_n)
            tR_n = sol_n[ii].u[end][1]
        end
    elseif n<1
        while tR_n-tR_lock>0 && n>0.01
            n = n*0.5
            par_n = stretched_program(n, par)
            sol_n = GasChromatographySimulator.solve_system_multithreads(par_n)
            tR_n = sol_n[ii].u[end][1]
        end
    end
    if n>130.0
        error("The choosen retention time for locking is to big.")
    elseif n<0.01
        error("The choosen retention time for locking is to small.")
    else
        return n
    end
end

function RT_locking(par::GasChromatographySimulator.Parameters, tR_lock::Float64, tR_reltol::Float64, solute_RT::String; opt_itp="linear")
	# estimate the factor 'n' for the temperature program to achieve the retention time 'tR_lock' for 'solute_RT' with the GC-system defined by 'par' 
	if isa(par, Array)==true
		error("Select an element of the array of GC-systems.")
	elseif tR_reltol<par.opt.reltol
		error("The tolerance for retention time locking tR_reltol is to small. Use tR_reltol >= par.opt.reltol ($(par.opt.reltol)).")
	else
		# find 'solute_RT' in the substances of 'par'
		name = Array{String}(undef, length(par.sub))
		for i=1:length(par.sub)
			name[i] = par.sub[i].name
		end
		ii = findfirst(name.==solute_RT)
		# calculate the retention time for the original (un-stretched) program 
		sol₀ = GasChromatographySimulator.solving_odesystem_r(par.sys, par.prog, par.sub[ii], par.opt)
		tR₀ = sol₀.u[end][1]
		# start value for the factor 'n'
		if tR₀-tR_lock<0
			n₁ = initial_n(2.0, tR_lock, ii, par)
		else
			n₁ = initial_n(0.5, tR_lock, ii, par)
		end
		# using a recursive function to estimate 'n'
		n = recur_RT_locking(n₁, [1.0], [tR₀], par, tR_lock, tR_reltol, ii; opt_itp=opt_itp)
	end
	return n		
end

function recur_RT_locking(n::Float64, n_vec::Array{Float64,1}, tR_vec::Array{Float64,1}, par::GasChromatographySimulator.Parameters, tR_lock::Float64, tR_reltol::Float64, ii::Int64; opt_itp="linear")
	# recursive function to find the factor 'n' for the temperature program to achieve the retention time 'tR_lock' for solute index 'ii' with the GC-system defined by 'par'
	# calculate the retention time with the input guess 'n'
	par₁ = stretched_program(n, par)
	if par.opt.odesys==true
		sol₁ = GasChromatographySimulator.solving_odesystem_r(par₁.sys, par₁.prog, par₁.sub[ii], par₁.opt)
		tR₁ = sol₁.u[end][1]
	else
		sol₁ = GasChromatographySimulator.solving_migration(par₁.sys, par₁.prog, par₁.sub[ii], par₁.opt)
		tR₁ = sol₁.u[end]
	end
	if abs(tR₁-tR_lock)/tR_lock<tR_reltol
		# if retention time is less than 'tR_tol' from 'tR_lock' we found the factor 'n'
		return n
	else
		# estimate a new factor 'new_n' by linear interpolation of the factors 'n_vec' + 'n' over the corresponding retention times
		new_n_vec = sort([n_vec; n])
		new_tR_vec = sort([tR_vec; tR₁])
		if opt_itp=="spline"
			# Dierckx.jl
			if length(new_tR_vec)<4
				k = length(new_tR_vec)-1
			else
				k = 3
			end
			itp = Spline1D(new_tR_vec, new_n_vec, k=k)
		else # opt_itp=="linear"
			# Interpolations.jl
			itp = LinearInterpolation(sort([tR_vec; tR₁]), sort([n_vec; n]))
		end
		new_n = itp(tR_lock)
		#println("new_n=$(new_n), tR₁=$(tR₁)")
		# use the new factor 'new_n' and call the recursive function again
		return recur_RT_locking(new_n, new_n_vec, new_tR_vec, par, tR_lock, tR_reltol, ii; opt_itp=opt_itp)
	end
end
#---end-RT-locking-----------------------------------------------------------------------------------

#---misc-functions-----------------------------------------------------------------------------------
function change_initial(par::GasChromatographySimulator.Parameters, init_τ, init_t)
	# copys the parameters `par` and changes the values of par.sub[i].τ₀ and par.sub[i].t₀ to init_τ[i] resp. init_t[i]
	newsub = Array{GasChromatographySimulator.Substance}(undef, length(par.sub))
	for i=1:length(par.sub)
		newsub[i] = GasChromatographySimulator.Substance(par.sub[i].name, par.sub[i].CAS, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀, par.sub[i].ann, par.sub[i].Dag, init_τ[i], init_t[i])
	end
	newpar = GasChromatographySimulator.Parameters(par.sys, par.prog, newsub, par.opt)
	return newpar
end

function general_step(x::Float64, L::Array{Float64,1}, a::Array{<:Any,1})
    # a generalized step function
    #
    # x ... variable (x-positon)
    # L ... Array with the different lengths of the column segments
    # a ... Array with the values for the different column segments
    #       values can also be functions
    if length(L)!=length(a)
        error("Parameters `L` and `a` must have the same length.")
        return
    end
	intervals = Array{Intervals.Interval}(undef, length(L))
	cumL = round.(cumsum([0; L]), digits=4)
	for i=1:length(intervals)
		if i==length(intervals)
			intervals[i] = Interval{Closed, Closed}(cumL[i],cumL[i+1])
		else
			intervals[i]=Interval{Closed, Open}(cumL[i],cumL[i+1])
		end
		if x in intervals[i]
			return a[i]
		end
	end
end

#----end-misc-functions-----------------------------------------------------------------------------

#----begin-notebooks-functions----------------------------------------------------------------------
##---begin-UI-functions-----------------------------------------------------------------------------
"""
    sys_set_UI(sp)

Construct a combined PlutoUI widget for the settings of the GC system with then selectable stationary phases `sp`. 
	
# UI fields
* ``L``: column length in m.
* ``d``: column diameter in mm.
* ``d_f``: film thickness in µm.
* stat. phase: stationary phase of the column
* Gas: mobile phase
"""
function sys_set_UI(sp)
		PlutoUI.combine() do Child
			@htl("""
			<h3>System settings</h3>
			``L`` [m]: $(
				Child(NumberField(0.1:0.1:100.0; default=10.0))
			) ``d`` [mm]: $(
				Child(NumberField(0.01:0.01:1.00; default=0.25))
			) ``d_f`` [µm]: $(
				Child(NumberField(0.01:0.01:1.00; default=0.25))
			) stat. phase: $(
				Child(Select(sp; default="SLB5ms"))
			) Gas: $(
				Child(Select(["He", "H2", "N2"]; default="He"))
			) 
			
			""")
	end
end

"""
    prog_set_UI()

Construct a combined PlutoUI widget for the settings of the program of a GC
system with a thermal gradient.

# UI fields
`time steps`: the time steps after which duration the values of temperature,
inlet pressure, ΔT and α are achieved by linear interpolation (in s).
`temperature steps`: the temperature steps (in °C). 
`ΔT steps`: the steps of the temperature difference (in °C) between column inlet
and outlet.
`α steps`: the steps of the gradient profile (α = 0 ... linear change of
temperature along column, α < 0 ... concave exponential profile, α > 0 ...
convexe exponential profile).
``p_{in}`` steps: the steps of the inlet pressure (in kPa(g))
``p_{out}`` steps: the steps of the outlet pressure (in kPa(a))
"""
function prog_set_UI()
	PlutoUI.combine() do Child
		@htl("""
		<h3>Program settings</h3> 
		_Note: Same number of entrys for every text field._
		
		$(
			Child(TextField((50,1); default="0 60 300 300 120"))
		) time steps [s] 
		
		$(
			Child(TextField((50,1); default="40 40 170 300 300"))
		) temperature steps [°C]
		
		$(
			Child(TextField((50,1); default="0 0 40 60 0"))
		) ``ΔT`` steps [°C]
		
		$(
			Child(TextField((50,1); default="-3 -3 -3 -3 -3"))
		) ``α`` steps

		$(
			Child(TextField((50,1); default="18 18 58 98 98"))
		) ``p_{in}`` steps [kPa(g)]

		$(
			Child(TextField((50,1); default="0 0 0 0 0"))
		)``p_{out}`` steps [kPa(a)]
			
		""")
	end
end

"""
    sub_set_UI(sol)

Construct a combined PlutoUI widget for the settings of the substances separated
in the simulated GC system with the selectable substances `subs`. 
	
# UI fields
* Select Substances: Selection of the substances, which will be simulated.
* Injection time: Start time (in s) of the simulation. The same for all selected
  substances.
* Injection width: Peak width (in s) of all selected substances at the time of injection.
"""
function sub_set_UI(sol)
	if length(sol)>10
		select_size = 10
	else
		select_size = length(sol)
	end
	PlutoUI.combine() do Child
		@htl("""
		<h3>Substance settings</h3> 
		
		Select Substances: $(
			Child(MultiSelect(sol; default=sol[1:4], size=select_size))
		) 
		
		Injection time [s]: $(
			Child(NumberField(0.0:0.1:100.0; default=0.0))
		) and Injection width [s]: $(
			Child(NumberField(0.00:0.01:10.0; default=0.0))
		) 
		""")
	end
end

"""
    opt_set_UI()

Construct a combined PlutoUI widget for the settings of the options for the simulation.    
"""
function opt_set_UI()
	PlutoUI.combine() do Child
		@htl("""
		<h3>Option settings</h3>
		
		abstol: 1e $(
			Child(NumberField(-10:1:-3; default=-8))
		) reltol: 1e $(
			Child(NumberField(-8:1:-2; default=-5))
		) Tcontrol: $(
			Child(Select(["inlet", "outlet"]; default="inlet"))
		)
		""")
	end
end
##---end-UI-functions-------------------------------------------------------------------------------

##---begin-Plot-functions---------------------------------------------------------------------------
"""
	local_plots(xx, yy, sol, par)

Show additional 'local' plots of selected `yy` quantities over selected `xx`
quantities.

# Arguments
* `xx`: Selected quantity shown on the x-axis. Possible values: "z", "t", "T",
  "τ", "σ" and "u".
* `yy`: Selected quantity shown on the y-axis. Possible values: "z", "t", "T",
  "τ", "σ" and "u".
* `sol`: The solution of the simulation.
* `par`: The parameters of the simulated GC-system.
"""   
function local_plots(xx, yy, sol, par)
	n = size(sol)[1]

	df_sol = GasChromatographySimulator.sol_extraction(sol, par)
	xvalues = Array{Array{Float64,1}}(undef, n)
	yvalues = Array{Array{Float64,1}}(undef, n)
	
	p_add = plot(legend=false)
	for i=1:n
		if xx=="z"
			xvalues[i] = df_sol.z[i]
			xlabel = "position z in m"
		elseif xx=="t"
			xvalues[i] = df_sol.t[i]
			xlabel = "time t in s"
		elseif xx=="T"
			xvalues[i] = par.prog.T_itp.(df_sol.z[i], df_sol.t[i]).-273.15
			xlabel = "temperature T in °C"
		elseif xx=="τ"
			xvalues[i] = sqrt.(df_sol.τ²[i])
			xlabel = "peak width τ in s"
		elseif xx=="σ"
			xvalues[i] = velocity(df_sol, i, par).*sqrt.(df_sol.τ²[i])
			xlabel = "band width in m"
		elseif xx=="u"
			xvalues[i] = velocity(df_sol, i, par)
			xlabel = "solute velocity in m/s"
		end
		if yy=="z"
			yvalues[i] = df_sol.z[i]
			ylabel = "position z in m"
		elseif yy=="t"
			yvalues[i] = df_sol.t[i]
			ylabel = "time t in s"
		elseif yy=="T"
			yvalues[i] = par.prog.T_itp.(df_sol.z[i], df_sol.t[i]).-273.15
			ylabel = "temperature T in °C"
		elseif yy=="τ"
			yvalues[i] = sqrt.(df_sol.τ²[i])
			ylabel = "peak width in s"
		elseif yy=="σ"
			yvalues[i] = velocity(df_sol, i, par).*sqrt.(df_sol.τ²[i])
			ylabel = "band width in m"
		elseif yy=="u"
			yvalues[i] = velocity(df_sol, i, par)
			ylabel = "solute velocity in m/s"
		end
		plot!(p_add, xvalues[i], yvalues[i], xlabel=xlabel, ylabel=ylabel, label=par.sub[i].name, m=:o)
	end
	return p_add
end

"""
	velocity(df_sol, i, par)

Calculate the velocity (in m/s) coressponding to solution of the `i-th` sunstance of a
GC-system defined by `par`.
""" 
function velocity(df_sol, i, par)
	x = df_sol.z[i]
	t = df_sol.t[i]
	T_itp = par.prog.T_itp
	pin_itp = par.prog.pin_itp
	pout_itp = par.prog.pout_itp
	L = par.sys.L
	d = par.sys.d
	df = par.sys.df
	gas = par.sys.gas
	ΔCp = par.sub[i].ΔCp
	Tchar = par.sub[i].Tchar
	θchar = par.sub[i].θchar
	φ₀ = par.sub[i].φ₀
	u = Array{Float64}(undef, length(x))
	for j=1:length(x)
		u[j] = 1/GasChromatographySimulator.residency(x[j], t[j], T_itp, pin_itp, pout_itp, L, d, df, gas, ΔCp, Tchar, θchar, φ₀)
	end
	return u
end
#----notebooks-functions----------------------------------------------------------------------------
end # module
