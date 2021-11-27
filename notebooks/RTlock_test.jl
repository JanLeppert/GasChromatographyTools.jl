### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ 6f9a94bd-dd28-4110-a7ca-0ed84e9c7c3f
begin
	import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
	using Interpolations
	using Plots
	using Dierckx
	using GasChromatographySimulator
	using GasChromatographyTools
	using Plots
	using PlutoUI
	TableOfContents()
end

# ╔═╡ 0b51a792-d9b8-4c29-928f-57aaf41f1c20
plotly()

# ╔═╡ 93ba6afc-4e9a-11ec-08d9-81c0a9bc502e
md"""
# Test of the RT lock functions
The test is necessary, because of failing (infinit loop in case of opt_ipt="linear", NaN abort in case of opt_itp="spline").
"""

# ╔═╡ 51ec0223-199a-4a25-8768-cd0b7e9f864d
begin
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

	function RT_locking(par::GasChromatographySimulator.Parameters, tR_lock::Float64, tR_tol::Float64, solute_RT::String; opt_itp="linear")
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
				n₁ = initial_n(2.0, tR_lock, ii, par)
			else
				n₁ = initial_n(0.5, tR_lock, ii, par)
			end
			# using a recursive function to estimate 'n'
			n = recur_RT_locking(n₁, [1.0], [tR₀], par, tR_lock, tR_tol, ii; opt_itp=opt_itp)
		end
		return n		
	end

	function recur_RT_locking(n::Float64, n_vec::Array{Float64,1}, tR_vec::Array{Float64,1}, par::GasChromatographySimulator.Parameters, tR_lock::Float64, tR_tol::Float64, ii::Int64; opt_itp="linear")
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
			return recur_RT_locking(new_n, new_n_vec, new_tR_vec, par, tR_lock, tR_tol, ii; opt_itp=opt_itp)
		end
	end
	
	md"""
	## Copies of the RT lock functions
	In the package GasChromatographyTools.jl the functions where adapted to solve this problem.
	"""
end

# ╔═╡ 5e642f09-9a8f-4cca-b61f-b27c8433a2e5
begin
	opt = GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-3, "inlet", true)
	L = 4.0
	d = 0.1e-3
	df = 0.1e-6
	sp = "SLB5ms" # ["Rxi17SilMS" -> ok, "SLB5ms" -> ERR, "SPB50" -> ok, "Wax" -> ok, "DB5ms" -> ok, "Rxi5MS" -> ok, "genericLB", "genericJL"]
	gas = "He"
	sys = GasChromatographySimulator.constructor_System(L, d, df, sp, gas)

	db_path = "/Users/janleppert/Documents/GitHub/Publication_GCsim/data/Databases/" 
	db_file = "Database_append.csv"
	first_alkane = "C8"
	last_alkane = "C15"
	sub = GasChromatographySimulator.load_solute_database(db_path, db_file, sp, gas, [first_alkane, last_alkane], zeros(2), zeros(2))
	
	Tst = 273.15
	Tref = 150.0 + Tst
	pn = 101300.0
	pin = 300000.0 + pn
	pout = 0.0
	dimless_rate = 0.4
	Theat = 1000.0
	Tshift = 40.0
	tMref = GasChromatographySimulator.holdup_time(Tref, pin, pout, L, d, gas)
	rate = dimless_rate*30/tMref
	theat = Theat/rate
	Tstart = sub[1].Tchar-Tst-Tshift
	
	ΔT = 30.0
	α = -3.0
	prog0 = GasChromatographySimulator.constructor_Program([0.0, theat],[Tstart, Tstart+Theat], [pin, pin],[pout, pout],[ΔT, ΔT], [0.0, 0.0], [L, L], [α, α], opt.Tcontrol, L)
	
	par0 = GasChromatographySimulator.Parameters(sys, prog0, sub, opt)

	tR_lock = 12*tMref
	tR_tol = 1e-3
	md"""
	## Settings
	The following settings produce the problem:
	"""
end

# ╔═╡ 3de48bb6-8eb7-4a33-9a98-d5fe3a19f6c6
md"""
## The ERROR
"""

# ╔═╡ 19d1abb6-32b6-414f-bc95-55635cbaa73a
n = RT_locking(par0, tR_lock, 1e-3, last_alkane; opt_itp="spline")

# ╔═╡ 2cad8185-4fc5-4009-98df-07a2be2133c6
md"""
## Reproducing the error by stepwise calculations
"""

# ╔═╡ 30207849-cafe-4ec6-af8d-ee7b2e2e6de0
begin
	# initial simulation
	sol₀ = GasChromatographySimulator.solving_odesystem_r(par0.sys, par0.prog, par0.sub[2], par0.opt)
	tR₀ = sol₀.u[end][1]
end

# ╔═╡ 84e7e869-0fbf-4590-a8de-28855856661f
tR₀-tR_lock

# ╔═╡ c5ace7ce-cfa3-4c15-bcc1-8b66ad1c16e9
begin
	# first stretch of the program
	if tR₀-tR_lock<0
	    n₁ = initial_n(2.0, tR_lock, 2, par0)
	else
	    n₁ = initial_n(0.5, tR_lock, 2, par0)
	end
	# recur function
	par₁ = stretched_program(n₁, par0)
	sol₁ = GasChromatographySimulator.solving_odesystem_r(par₁.sys, par₁.prog, par₁.sub[2], par₁.opt)
	tR₁ = sol₁.u[end][1]
	n_vec₁ = sort([1.0; n₁])
	tR_vec₁ = sort([tR₀; tR₁])
	# first interpolation
	itp₁ = Spline1D(tR_vec₁, n_vec₁, k=1)
end

# ╔═╡ c4d53222-03b1-4ce4-a49a-690945347432
tR₁-tR_lock

# ╔═╡ 6acfbf9b-f8ac-4483-8c93-17920d0d9f0e
begin
	# second stretch 
	n₂ = itp₁(tR_lock)
	par₂ = stretched_program(n₂, par0)
	sol₂ = GasChromatographySimulator.solving_odesystem_r(par₂.sys, par₂.prog, par₂.sub[2], par₂.opt)
	tR₂ = sol₂.u[end][1]
	n_vec₂ = sort([n_vec₁; n₂])
	tR_vec₂ = sort([tR_vec₁; tR₂])
	itp₂ = Spline1D(tR_vec₂, n_vec₂, k=2)

	p1 = plot([tR_lock, tR_lock], [0.95, 0.98], label="tR_lock")
	scatter!(p1, [tR₂, tR₂], [n₂, n₂], label="2")
	plot!(p1, 51.0:0.001:52.0, itp₂.(51.0:0.001:52.0), xlims=(51.0,52.0), ylims=(0.95, 0.98), label="itp₂")
	p1
end

# ╔═╡ 3b31f580-efe2-4a7a-89dc-228b38f2a71e
tR₂-tR_lock

# ╔═╡ b6c5126a-3d41-4809-ab9b-d554262e0668
begin
	n₃ = itp₂(tR_lock)
	par₃ = stretched_program(n₃, par0)
	sol₃ = GasChromatographySimulator.solving_odesystem_r(par₃.sys, par₃.prog, par₃.sub[2], par₃.opt)
	tR₃ = sol₃.u[end][1]
	n_vec₃ = sort([n_vec₂; n₃])
	tR_vec₃ = sort([tR_vec₂; tR₃])
	k = length(tR_vec₃)-1
	itp₃ = Spline1D(tR_vec₃, n_vec₃, k=k)
	
	scatter!(p1, [tR₃, tR₃], [n₃, n₃], label="3")
	plot!(p1, 51.0:0.001:52.0, itp₃.(51.0:0.001:52.0), label="itp₃")
end

# ╔═╡ 6c5443bc-3dab-4317-99c9-1bb32b184fbc
tR₃-tR_lock

# ╔═╡ 3b1bb187-66f9-40bd-a0a5-a3c4e1a23819
begin
	n₄ = itp₃(tR_lock)
	par₄ = stretched_program(n₄, par0)
	sol₄ = GasChromatographySimulator.solving_odesystem_r(par₄.sys, par₄.prog, par₄.sub[2], par₄.opt)
	tR₄ = sol₄.u[end][1]
	n_vec₄ = sort([n_vec₃; n₄])
	tR_vec₄ = sort([tR_vec₃; tR₄])
	itp₄ = Spline1D(tR_vec₄, n_vec₄, k=3)
	
	scatter!(p1, [tR₄, tR₄], [n₄, n₄], label="4")
	plot!(p1, 51.0:0.001:52.0, itp₄.(51.0:0.001:52.0), label="itp₄")
end

# ╔═╡ b81b2905-6756-445b-b1be-c54298df8b3f
tR₄-tR_lock

# ╔═╡ eb2d9f6a-8384-427a-a036-b4f0019d8251
md"""
The proposed settings seam to lead to jumping estimations around the searched retention time.
"""

# ╔═╡ 734bae50-12fa-4717-8446-addb963b8673
begin 
	n₅ = itp₄(tR_lock)
	par₅ = GasChromatographyTools.stretched_program(n₅, par0)
	sol₅ = GasChromatographySimulator.solving_odesystem_r(par₅.sys, par₅.prog, par₅.sub[2], par₅.opt)
	tR₅ = sol₅.u[end][1]
	n_vec₅ = sort([n_vec₄; n₅])
	tR_vec₅ = sort([tR_vec₄; tR₅])
	itp₅ = Spline1D(tR_vec₅, n_vec₅, k=3)

	scatter!(p1, [tR₅, tR₅], [n₅, n₅], label="5")
	plot!(p1, 51.0:0.001:52.0, itp₅.(51.0:0.001:52.0), xlims=(51.45,51.55), ylims=(0.969, 0.972), label="itp₅")
end

# ╔═╡ 8de716b9-ca56-42f5-aebf-39ad467f4613
tR₅-tR_lock

# ╔═╡ dadbefbf-a107-4c89-a4ef-0eb756517c1e
begin
	n₆ = itp₅(tR_lock)
	par₆ = GasChromatographyTools.stretched_program(n₆, par0)
	sol₆ = GasChromatographySimulator.solving_odesystem_r(par₆.sys, par₆.prog, par₆.sub[2], par₆.opt)
	tR₆ = sol₆.u[end][1]
	n_vec₆ = sort([n_vec₅; n₆])
	tR_vec₆ = sort([tR_vec₅; tR₆])
	tR_vec₆.-tR_lock
	itp₆ = Spline1D(tR_vec₆, n_vec₆, k=3)

	scatter!(p1, [tR₆, tR₆], [n₆, n₆], label="6")
	plot!(p1, 51.0:0.001:52.0, itp₆.(51.0:0.001:52.0), xlims=(51.45,51.55), ylims=(0.969, 0.972), label="itp₆")
end

# ╔═╡ 313d1431-ec39-4d1e-91fd-e025efc1f5c3
tR₆-tR_lock

# ╔═╡ fab59b1e-ba32-49d3-b0f1-564db037400c
begin
	n₇ = itp₆(tR_lock)
	par₇ = GasChromatographyTools.stretched_program(n₇, par0)
	sol₇ = GasChromatographySimulator.solving_odesystem_r(par₇.sys, par₇.prog, par₇.sub[2], par₇.opt)
	tR₇ = sol₇.u[end][1]
	n_vec₇ = sort([n_vec₆; n₇])
	tR_vec₇ = sort([tR_vec₆; tR₇])
	tR_vec₇.-tR_lock
	itp₇ = Spline1D(tR_vec₇, n_vec₇, k=3)

	scatter!(p1, [tR₇, tR₇], [n₇, n₇], label="7")
	plot!(p1, 51.0:0.001:52.0, itp₇.(51.0:0.001:52.0), xlims=(51.47,51.49), ylims=(0.9699, 0.9701), label="itp₇")
end

# ╔═╡ 25725568-ae9a-4ab0-8f98-87fec12c867a
begin
	n₈ = itp₇(tR_lock)
	par₈ = GasChromatographyTools.stretched_program(n₈, par0)
	sol₈ = GasChromatographySimulator.solving_odesystem_r(par₈.sys, par₈.prog, par₈.sub[2], par₈.opt)
	tR₈ = sol₈.u[end][1]
	n_vec₈ = sort([n_vec₇; n₈])
	tR_vec₈ = sort([tR_vec₇; tR₈])
	tR_vec₈.-tR_lock
	itp₈ = Spline1D(tR_vec₈, n_vec₈, k=3)

	scatter!(p1, [tR₈, tR₈], [n₈, n₈], label="8")
	plot!(p1, 51.0:0.001:52.0, itp₈.(51.0:0.001:52.0), xlims=(51.47,51.49), ylims=(0.9699, 0.9701), label="itp₈")
end

# ╔═╡ 404c1d12-5217-485e-b3a6-6f024d29b544
md"""
The estimates continue to jump around tR_lock. The last estimate is further away than the one before.

Also, there seems to be some form of discontiuity of the simulation result around tR_lock, which produces the problem in the first place. 
"""

# ╔═╡ 57634f00-45d5-4cf9-94f5-76319d4e5436
begin
	n₉ = itp₈(tR_lock)
	par₉ = GasChromatographyTools.stretched_program(n₉, par0)
	sol₉ = GasChromatographySimulator.solving_odesystem_r(par₉.sys, par₉.prog, par₉.sub[2], par₉.opt)
	tR₉ = sol₉.u[end][1]
	n_vec₉ = sort([n_vec₈; n₉])
	tR_vec₉ = sort([tR_vec₈; tR₉])
	tR_vec₉.-tR_lock
	itp₉ = Spline1D(tR_vec₉, n_vec₉, k=3)

	scatter!(p1, [tR₉, tR₉], [n₉, n₉], label="9")
	plot!(p1, 51.0:0.001:52.0, itp₉.(51.0:0.001:52.0), xlims=(51.47,51.49), ylims=(0.9699, 0.9701), label="itp₉")
end

# ╔═╡ ae803a26-c436-497a-b245-a0afec92e46f
md"""
Make simulations for a range of stretch factors ``n``.
"""

# ╔═╡ 1983258e-e84b-4f39-9ce8-0e20e78a0893
begin
	nn = [0.969970, 0.969971, 0.969972, 0.969973, 0.969974, 0.969975, 0.969976, 0.969977]
	ttR = Array{Float64}(undef, length(nn))
	for i=1:length(nn)
	    pars = stretched_program(nn[i], par0)
	    sols = GasChromatographySimulator.solving_odesystem_r(pars.sys, pars.prog, pars.sub[2], pars.opt)
	    ttR[i] = sols.u[end][1]
	end
	p2 = plot([tR_lock, tR_lock], [0.95, 0.98], label="tR_lock")
	plot!(p2, ttR, nn, line=(2,:solid), markers=:square, xlims=(51.479, 51.484), ylims=(0.969969,0.96998))	
end

# ╔═╡ 89cbd5c0-0df6-43d3-bd5c-d3d55d669f33
begin
	nnn = 0.969974.+collect(0.0000001:0.0000001:0.0000009)
	tttR = Array{Float64}(undef, length(nnn))
	for i=1:length(nnn)
	    pars = stretched_program(nnn[i], par0)
	    sols = GasChromatographySimulator.solving_odesystem_r(pars.sys, pars.prog, pars.sub[2], pars.opt)
	    tttR[i] = sols.u[end][1]
	end
	plot!(p2, tttR, nnn, ylims=(0.9699739,0.9699751), markers=:circle)
end

# ╔═╡ f541b210-3ef8-4003-b358-065e5c5949ad
par8 = stretched_program(nnn[8], par0)

# ╔═╡ 524b256b-0c1e-48bc-ac0b-2ee7fc1eab2b
par9 = stretched_program(nnn[9], par0)

# ╔═╡ 6d15c2a4-6d1e-4b89-b9e8-ecff004c4730
sol8 = GasChromatographySimulator.solving_odesystem_r(par8.sys, par8.prog, par8.sub[2], par8.opt)

# ╔═╡ 665f302d-204b-47ac-90df-c5979350707c
sol9 = GasChromatographySimulator.solving_odesystem_r(par9.sys, par9.prog, par9.sub[2], par9.opt)

# ╔═╡ 56095d71-6169-44b3-89d1-7ea7f1b6ddfb
sol8.destats

# ╔═╡ a51cd1dc-50bb-4444-a0f7-c0e4229b1257
sol9.destats

# ╔═╡ 225434d5-c753-4476-8f68-f5761a454852
sol8.retcode

# ╔═╡ 9a6f436e-e50e-4697-a1d1-3d2d43fc62fc
sol9.retcode

# ╔═╡ 8a6eccf0-8740-4d69-b96b-1248557d4c4d
begin 
	nnnn = sort!(rand(0.9699745:0.000000001:0.9699755, 100))
	ttttR = Array{Float64}(undef, length(nnnn))
	for i=1:length(nnnn)
	    pars = stretched_program(nnnn[i], par0)
	    sols = GasChromatographySimulator.solving_odesystem_r(pars.sys, pars.prog, pars.sub[2], pars.opt)
	    ttttR[i] = sols.u[end][1]
	end
	plot!(p2, ttttR, nnnn, ylims=(0.9699739,0.9699756), markers=:circle)
	md"""
	In the range around tR_lock the solution seems to be unstable. A small change in the temperature program leads to a bigger change in the retention time. There seem to be two branches of n(tR), which does not meet.

	$(embed_display(p2))
	"""
end

# ╔═╡ 08a0972a-55ff-42c9-838b-1a649afe9a46
md"""
## Decreased relative tolerance

The problem is not originated in the RT_lock algorithm but in the simulation itself. For certain settings a small change (like small strectch in the program) can result in bigger changes in retention time.
"""

# ╔═╡ 2d7da55e-2711-44a4-b8b9-bda6321a4c48
begin
	# decrease the relative tolerance
	opt_1 = GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-4, "inlet", true)
	par_1 = GasChromatographySimulator.Parameters(sys, prog0, sub, opt_1)
	# repeat the simulation from above
	nnnn_1 = sort!(rand(0.9699745:0.000000001:0.9699755, 100))
	ttttR_1 = Array{Float64}(undef, length(nnnn_1))
	for i=1:length(nnnn_1)
	    pars = stretched_program(nnnn_1[i], par_1)
	    sols = GasChromatographySimulator.solving_odesystem_r(pars.sys, pars.prog, pars.sub[2], pars.opt)
	    ttttR_1[i] = sols.u[end][1]
	end
	plot!(p2, ttttR_1, nnnn_1, ylims=(0.9699739,0.9699756), markers=:v, label="reltol=1e-4")
end

# ╔═╡ 89fa3139-1d10-4dd3-90b1-d9c0f65745c6
begin
	n_1 = RT_locking(par_1, tR_lock, 1e-3, last_alkane; opt_itp="spline")
	par_1_s = stretched_program(n_1, par_1)
	sol_1_s = GasChromatographySimulator.solving_odesystem_r(par_1_s.sys, par_1_s.prog, par_1_s.sub[2], par_1_s.opt)
	tR_1_s = sol_1_s.u[end][1]
	scatter!(p2, [tR_1_s, tR_1_s], [n_1, n_1], ylims=(0.9699739, n_1*1.000001))
end

# ╔═╡ e9accaac-587f-498d-ab0c-a8064f573678
begin
	# simulation arround n_1
	n_range = (n_1-0.0001):0.000001:(n_1+0.0001)
	tR_range = Array{Float64}(undef, length(n_range))
	for i=1:length(n_range)
	    pars = stretched_program(n_range[i], par_1)
	    sols = GasChromatographySimulator.solving_odesystem_r(pars.sys, pars.prog, pars.sub[2], pars.opt)
	    tR_range[i] = sols.u[end][1]
	end
	plot!(p2, tR_range, n_range, ylims=(n_range[1], n_range[end]), xlims=(tR_range[1], tR_range[end]))
end

# ╔═╡ 4cd795b1-efeb-430b-a025-77d0ee9d2b1b
md"""
## Modified Functions with relative difference
"""

# ╔═╡ 0bf4903c-074d-41b8-a4e8-5834ce8659d5
function recur_RT_locking_rel(n::Float64, n_vec::Array{Float64,1}, tR_vec::Array{Float64,1}, par::GasChromatographySimulator.Parameters, tR_lock::Float64, tR_tol::Float64, ii::Int64; opt_itp="linear")
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
	if abs(tR₁-tR_lock)/tR_lock<tR_tol
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
		return recur_RT_locking_rel(new_n, new_n_vec, new_tR_vec, par, tR_lock, tR_tol, ii; opt_itp=opt_itp)
	end
end

# ╔═╡ be7d3bf5-d682-40e7-8161-4243354373f8
function RT_locking_rel(par::GasChromatographySimulator.Parameters, tR_lock::Float64, tR_tol::Float64, solute_RT::String; opt_itp="linear")
	# estimate the factor 'n' for the temperature program to achieve the retention time 'tR_lock' for 'solute_RT' with the GC-system defined by 'par' 
	if isa(par, Array)==true
		error("Select an element of the array of GC-systems.")
	elseif tR_tol<par.opt.reltol
		error("The relative tolerance for retention time locking tR_tol is to small. Use tR_tol > par.opt.reltol ($(par.opt.reltol)).")
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
		n = recur_RT_locking_rel(n₁, [1.0], [tR₀], par, tR_lock, tR_tol, ii; opt_itp=opt_itp)
	end
	return n		
end

# ╔═╡ 8a7a24e3-3896-4012-8e50-473a70ec3a44
n_rel = RT_locking_rel(par0, tR_lock, 1e-3, last_alkane; opt_itp="spline")

# ╔═╡ 4cbd8054-d1c7-4132-ae61-ef09ce6ba9e7
par0_rel = stretched_program(n_rel, par0)

# ╔═╡ b842960f-ef61-4a70-8a17-5405ec3e5ed3
sol_rel = GasChromatographySimulator.solving_odesystem_r(par0_rel.sys, par0_rel.prog, par0_rel.sub[2], par0_rel.opt)

# ╔═╡ 24269e34-fc1d-4ee2-87b6-afb6a3081596
tR_rel = sol_rel.u[end][1]

# ╔═╡ 0f7e36cb-c4a7-43b5-928f-d364b82d218b
tR_rel-tR_lock

# ╔═╡ 3bdf621e-969b-42b8-85f8-5ecd261646d3
(tR_rel-tR_lock)/tR_lock

# ╔═╡ 42a1fccd-e782-42fb-91de-18c8b33d507d
md"""
# End
"""

# ╔═╡ Cell order:
# ╠═6f9a94bd-dd28-4110-a7ca-0ed84e9c7c3f
# ╠═0b51a792-d9b8-4c29-928f-57aaf41f1c20
# ╠═93ba6afc-4e9a-11ec-08d9-81c0a9bc502e
# ╠═51ec0223-199a-4a25-8768-cd0b7e9f864d
# ╠═5e642f09-9a8f-4cca-b61f-b27c8433a2e5
# ╠═3de48bb6-8eb7-4a33-9a98-d5fe3a19f6c6
# ╠═19d1abb6-32b6-414f-bc95-55635cbaa73a
# ╠═2cad8185-4fc5-4009-98df-07a2be2133c6
# ╠═30207849-cafe-4ec6-af8d-ee7b2e2e6de0
# ╠═84e7e869-0fbf-4590-a8de-28855856661f
# ╠═c5ace7ce-cfa3-4c15-bcc1-8b66ad1c16e9
# ╠═c4d53222-03b1-4ce4-a49a-690945347432
# ╠═6acfbf9b-f8ac-4483-8c93-17920d0d9f0e
# ╠═3b31f580-efe2-4a7a-89dc-228b38f2a71e
# ╠═b6c5126a-3d41-4809-ab9b-d554262e0668
# ╠═6c5443bc-3dab-4317-99c9-1bb32b184fbc
# ╠═3b1bb187-66f9-40bd-a0a5-a3c4e1a23819
# ╠═b81b2905-6756-445b-b1be-c54298df8b3f
# ╟─eb2d9f6a-8384-427a-a036-b4f0019d8251
# ╠═734bae50-12fa-4717-8446-addb963b8673
# ╠═8de716b9-ca56-42f5-aebf-39ad467f4613
# ╠═dadbefbf-a107-4c89-a4ef-0eb756517c1e
# ╠═313d1431-ec39-4d1e-91fd-e025efc1f5c3
# ╠═fab59b1e-ba32-49d3-b0f1-564db037400c
# ╠═25725568-ae9a-4ab0-8f98-87fec12c867a
# ╠═404c1d12-5217-485e-b3a6-6f024d29b544
# ╠═57634f00-45d5-4cf9-94f5-76319d4e5436
# ╠═ae803a26-c436-497a-b245-a0afec92e46f
# ╠═1983258e-e84b-4f39-9ce8-0e20e78a0893
# ╠═89cbd5c0-0df6-43d3-bd5c-d3d55d669f33
# ╠═f541b210-3ef8-4003-b358-065e5c5949ad
# ╠═524b256b-0c1e-48bc-ac0b-2ee7fc1eab2b
# ╠═6d15c2a4-6d1e-4b89-b9e8-ecff004c4730
# ╠═665f302d-204b-47ac-90df-c5979350707c
# ╠═56095d71-6169-44b3-89d1-7ea7f1b6ddfb
# ╠═a51cd1dc-50bb-4444-a0f7-c0e4229b1257
# ╠═225434d5-c753-4476-8f68-f5761a454852
# ╠═9a6f436e-e50e-4697-a1d1-3d2d43fc62fc
# ╠═8a6eccf0-8740-4d69-b96b-1248557d4c4d
# ╟─08a0972a-55ff-42c9-838b-1a649afe9a46
# ╠═2d7da55e-2711-44a4-b8b9-bda6321a4c48
# ╠═89fa3139-1d10-4dd3-90b1-d9c0f65745c6
# ╠═e9accaac-587f-498d-ab0c-a8064f573678
# ╠═4cd795b1-efeb-430b-a025-77d0ee9d2b1b
# ╠═be7d3bf5-d682-40e7-8161-4243354373f8
# ╠═0bf4903c-074d-41b8-a4e8-5834ce8659d5
# ╠═8a7a24e3-3896-4012-8e50-473a70ec3a44
# ╠═4cbd8054-d1c7-4132-ae61-ef09ce6ba9e7
# ╠═b842960f-ef61-4a70-8a17-5405ec3e5ed3
# ╠═24269e34-fc1d-4ee2-87b6-afb6a3081596
# ╠═0f7e36cb-c4a7-43b5-928f-d364b82d218b
# ╠═3bdf621e-969b-42b8-85f8-5ecd261646d3
# ╠═42a1fccd-e782-42fb-91de-18c8b33d507d
