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
include("./RTlocking.jl")
#---end-RT-locking-------------------------------------------------



#----begin-notebooks-functions----------------------------------------------------------------------
## UI-functions
include("./UI.jl")

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

"""
	common(s_1, s_1)

Compare two arrays and return the common elements.
"""
function common(s_1, s_2)
	common_s = String[]
	if length(s_1) >= length(s_2)
		for i=1:length(s_1)
			if s_1[i] in s_2
				push!(common_s, s_1[i])
			end
		end
	else
		for i=1:length(s_2)
			if s_2[i] in s_1
				push!(common_s, s_2[i])
			end
		end
	end
	return common_s
end

"""
	compare_peaklist(pl_1, pl_2)

Compare two peaklists (results of GasChromatographySimulator.jl) and calculate
absolute and relative differences of retention times and peak widths.
"""
function compare_peaklist(pl_1, pl_2)
	name = pl_1.Name
	tR1 = pl_1.tR
	τR1 = pl_1.τR
	tR2 = Array{Float64}(undef, size(pl_1)[1])
	τR2 = Array{Float64}(undef, size(pl_1)[1])
	for i=1:size(pl_1)[1]
		i2 = findfirst(name[i].==pl_2.Name)
		tR2[i] = pl_2.tR[i2]
		τR2[i] = pl_2.τR[i2]
	end
	ΔtR = tR1 .- tR2
	ΔτR = τR1 .- τR2
	rel_tR = ΔtR.*100.0./tR1
	rel_τR = ΔτR.*100.0./τR1
	compare_pl = DataFrame(Name=name, tR1=tR1, tR2=tR2, ΔtR=ΔtR, rel_tR=rel_tR, τR1=τR1, τR2=τR2, ΔτR=ΔτR, rel_τR=rel_τR)
	return compare_pl
end

"""
	compare_measurement_simulation(meas, peaklist)

Compare the retention times of measured and simulated substances.

# Arguments
* `meas`: DataFrame consisting at least of the columns `:Name` and `:RT`
  (measured retention time in s)
* `peaklist`: DataFrame as result from GasChromatographySimulator.jl with the
  columns `:Name` and `:tR` (simulated retention time in s).
  
The comparison is done by searching the same `Name` of the substance in both
DataFrames and calculating the absolute difference of the retention times (in s)
and the relative difference (in %).
"""
function compare_measurement_simulation(meas, peaklist)
	name = meas.Name
	tRm = meas.RT
	tRs = Array{Float64}(undef, size(meas)[1])
	for i=1:size(meas)[1]
		i2 = findfirst(name[i].==peaklist.Name)
		if typeof(i2) == Nothing
			tRs[i] = NaN
		else
			tRs[i] = peaklist.tR[i2]
		end
	end
	compare_df = DataFrame(Name=name, measured_tR=tRm, simulated_tR=tRs, ΔtR=tRm.-tRs, rel_tR=(tRm.-tRs)./tRm.*100.0)
	return compare_df
end

#----notebooks-functions----------------------------------------------------------------------------


#---misc-functions-----------------------------------------------------------------------------------
"""
	change_initial(par, init_τ, init_t)

Change the inital time and peak widths for the substances in a defined GC-system
`par` to the values `init_τ` and `init_t`.
""" 
function change_initial(par::GasChromatographySimulator.Parameters, init_t, init_τ)
	# copys the parameters `par` and changes the values of par.sub[i].τ₀ and par.sub[i].t₀ to init_τ[i] resp. init_t[i]
	newsub = Array{GasChromatographySimulator.Substance}(undef, length(par.sub))
	for i=1:length(par.sub)
		newsub[i] = GasChromatographySimulator.Substance(par.sub[i].name, par.sub[i].CAS, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀, par.sub[i].ann, par.sub[i].Dag, init_t[i], init_τ[i])
	end
	newpar = GasChromatographySimulator.Parameters(par.sys, par.prog, newsub, par.opt)
	return newpar
end

"""
	general_step(x, L, a)

A generalized step function

# Arguments
* `x::Float64`: variable (e.g. x-position)
* `L::Array{Float64,1}`: Array with different length of the constant segments
* `a::Array{<:Any,1}`: Array with the values for the segments, can also be
  functions.
"""
function general_step(x::Float64, L::Array{Float64,1}, a::Array{<:Any,1})
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
end # module
