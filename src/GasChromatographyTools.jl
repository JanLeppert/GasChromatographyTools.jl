module GasChromatographyTools

using GasChromatographySimulator
using Interpolations
using Dierckx

# helper function for 'RT-locking' and 'stretching' of GC-program
function stretched_program(n::Float64, par::GasChromatographySimulator.Parameters)
	# stretch the temperature program in 'par' by a factor 'n'
	if isa(par, Array)==true
		error("Select an element of the array of GC-system parameters.")
	else
		new_tsteps = n.*par.prog.time_steps
		new_T_itp = GasChromatographySimulator.temperature_interpolation(new_tsteps, par.prog.temp_steps, par.prog.gf, par.sys.L)
		new_pin_itp = LinearInterpolation((cumsum(new_tsteps), ), par.prog.pin_steps, extrapolation_bc=Flat())
		new_pout_itp = LinearInterpolation((cumsum(new_tsteps), ), par.prog.pout_steps, extrapolation_bc=Flat())
		new_prog = GasChromatographySimulator.Program(new_tsteps, par.prog.temp_steps, par.prog.pin_steps, par.prog.pout_steps, par.prog.gf, par.prog.a_gf, new_T_itp, new_pin_itp, new_pout_itp)
		new_par = GasChromatographySimulator.Parameters(par.sys, new_prog, par.sub, par.opt)
		return new_par
	end	
end

function RT_locking(par::GasChromatographySimulator.Parameters, tR_lock::Float64, tR_tol::Float64, solute_RT::String)
	# estimate the factor 'n' for the temperature program to achieve the retention time 'tR_lock' for 'solute_RT' with the GC-system defined by 'par' 
	if isa(par, Array)==true
		error("Select an element of the array of GC-systems.")
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
			n₁ = 1.5
		else
			n₁ = 0.5
		end
		# using a recursive function to estimate 'n'
		n = recur_RT_locking(n₁, [1.0], [tR₀], par, tR_lock, tR_tol, ii)
	end
	return n		
end

function recur_RT_locking(n::Float64, n_vec::Array{Float64,1}, tR_vec::Array{Float64,1}, par::GasChromatographySimulator.Parameters, tR_lock::Float64, tR_tol::Float64, ii::Int64)
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
	if abs(tR₁-tR_lock)<tR_tol
		# if retention time is less than 'tR_tol' from 'tR_lock' we found the factor 'n'
		return n
	else
		# estimate a new factor 'new_n' by linear interpolation of the factors 'n_vec' + 'n' over the corresponding retention times
		new_n_vec = sort([n_vec; n])
		new_tR_vec = sort([tR_vec; tR₁])
		# Dierckx.jl
		if length(new_tR_vec)<4
			k = length(new_tR_vec)-1
		else
			k = 3
		end
		itp = Spline1D(new_tR_vec, new_n_vec, k=k)
		# Interpolations.jl
		#itp = LinearInterpolation(sort([tR_vec; tR₁]), sort([n_vec; n]))
		new_n = itp(tR_lock)
		#println("new_n=$(new_n), tR₁=$(tR₁)")
		# use the new factor 'new_n' and call the recursive function again
		return recur_RT_locking(new_n, new_n_vec, new_tR_vec, par, tR_lock, tR_tol, ii)
	end
end
#----------------------------------

function change_initial(par::GasChromatographySimulator.Parameters, init_τ, init_t)
	# copys the parameters `par` and changes the values of par.sub[i].τ₀ and par.sub[i].t₀ to init_τ[i] resp. init_t[i]
	newsub = Array{GasChromatographySimulator.Substance}(undef, length(par.sub))
	for i=1:length(par.sub)
		newsub[i] = GasChromatographySimulator.Substance(par.sub[i].name, par.sub[i].CAS, par.sub[i].Tchar, par.sub[i].θchar, par.sub[i].ΔCp, par.sub[i].φ₀, par.sub[i].ann, par.sub[i].Dag, init_τ[i], init_t[i])
	end
	newpar = GasChromatographySimulator.Parameters(par.sys, par.prog, newsub, par.opt)
	return newpar
end

end # module
