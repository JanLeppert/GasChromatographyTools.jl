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

# ╔═╡ d3c2c2f4-1100-438d-88d6-d7502bcb9310
using Statistics

# ╔═╡ 0b51a792-d9b8-4c29-928f-57aaf41f1c20
plotly()

# ╔═╡ 93ba6afc-4e9a-11ec-08d9-81c0a9bc502e
md"""
# Test of the RT lock functions
The test is necessary, because of failing (infinit loop in case of opt\_ipt="linear", NaN abort in case of opt\_itp="spline").
"""

# ╔═╡ 51ec0223-199a-4a25-8768-cd0b7e9f864d
begin
	function stretched_program(n::Float64, par::GasChromatographySimulator.Parameters)
		# stretch the temperature program in 'par' by a factor 'n'
		if isa(par, Array)==true
			error("Select an element of the array of GC-system parameters.")
		else
			new_tsteps = n.*par.prog.time_steps
			new_T_itp = GasChromatographySimulator.temperature_interpolation(new_tsteps, par.prog.temp_steps, par.prog.gf, par.col.L)
			new_pin_itp = GasChromatographySimulator.pressure_interpolation(new_tsteps, par.prog.pin_steps)
			new_pout_itp = GasChromatographySimulator.pressure_interpolation(new_tsteps, par.prog.pout_steps)
			new_prog = GasChromatographySimulator.Program(new_tsteps, par.prog.temp_steps, par.prog.pin_steps, par.prog.pout_steps, par.prog.gf, par.prog.a_gf, new_T_itp, new_pin_itp, new_pout_itp)
			new_par = GasChromatographySimulator.Parameters(par.col, new_prog, par.sub, par.opt)
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
			sol₀ = GasChromatographySimulator.solving_odesystem_r(par.col, par.prog, par.sub[ii], par.opt)
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
			sol₁ = GasChromatographySimulator.solving_odesystem_r(par₁.col, par₁.prog, par₁.sub[ii], par₁.opt)
			tR₁ = sol₁.u[end][1]
		else
			sol₁ = GasChromatographySimulator.solving_migration(par₁.col, par₁.prog, par₁.sub[ii], par₁.opt)
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
	opt = GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-4, "inlet", true)
	L = 4.0
	d = 0.1e-3
	df = 0.1e-6
	sp = "Rxi17SilMS" # ["Rxi17SilMS" -> ERR, "SLB5ms", "SPB50", "Wax", "DB5ms", "Rxi5MS", "genericLB", "genericJL"]
	gas = "He"
	col = GasChromatographySimulator.constructor_System(L, d, df, sp, gas)

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
	
	ΔT = 135.0
	α = 9.0
	prog0 = GasChromatographySimulator.constructor_Program([0.0, theat],[Tstart, Tstart+Theat], [pin, pin],[pout, pout],[ΔT, ΔT], [0.0, 0.0], [L, L], [α, α], opt.Tcontrol, L)
	
	par0 = GasChromatographySimulator.Parameters(col, prog0, sub, opt)

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
	sol₀ = GasChromatographySimulator.solving_odesystem_r(par0.col, par0.prog, par0.sub[2], par0.opt)
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
	sol₁ = GasChromatographySimulator.solving_odesystem_r(par₁.col, par₁.prog, par₁.sub[2], par₁.opt)
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
	sol₂ = GasChromatographySimulator.solving_odesystem_r(par₂.col, par₂.prog, par₂.sub[2], par₂.opt)
	tR₂ = sol₂.u[end][1]
	n_vec₂ = sort([n_vec₁; n₂])
	tR_vec₂ = sort([tR_vec₁; tR₂])
	itp₂ = Spline1D(tR_vec₂, n_vec₂, k=2)

	p1 = plot([tR_lock, tR_lock], [0.4, 0.6], label="tR_lock")
	scatter!(p1, [tR₂, tR₂], [n₂, n₂], label="2")
	plot!(p1, 40.0:0.001:60.0, itp₂.(40.0:0.001:60.0), xlims=(50.0,55.0), ylims=(0.40, 0.60), label="itp₂")
	p1
end

# ╔═╡ 3b31f580-efe2-4a7a-89dc-228b38f2a71e
tR₂-tR_lock

# ╔═╡ b6c5126a-3d41-4809-ab9b-d554262e0668
begin
	n₃ = itp₂(tR_lock)
	par₃ = stretched_program(n₃, par0)
	sol₃ = GasChromatographySimulator.solving_odesystem_r(par₃.col, par₃.prog, par₃.sub[2], par₃.opt)
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
	sol₄ = GasChromatographySimulator.solving_odesystem_r(par₄.col, par₄.prog, par₄.sub[2], par₄.opt)
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
	sol₅ = GasChromatographySimulator.solving_odesystem_r(par₅.col, par₅.prog, par₅.sub[2], par₅.opt)
	tR₅ = sol₅.u[end][1]
	n_vec₅ = sort([n_vec₄; n₅])
	tR_vec₅ = sort([tR_vec₄; tR₅])
	itp₅ = Spline1D(tR_vec₅, n_vec₅, k=3)

	scatter!(p1, [tR₅, tR₅], [n₅, n₅], label="5")
	plot!(p1, 51.0:0.001:52.0, itp₅.(51.0:0.001:52.0), xlims=(51.47,51.485), ylims=(0.4995, 0.49975), label="itp₅")
end

# ╔═╡ 8de716b9-ca56-42f5-aebf-39ad467f4613
tR₅-tR_lock

# ╔═╡ dadbefbf-a107-4c89-a4ef-0eb756517c1e
begin
	n₆ = itp₅(tR_lock)
	par₆ = GasChromatographyTools.stretched_program(n₆, par0)
	sol₆ = GasChromatographySimulator.solving_odesystem_r(par₆.col, par₆.prog, par₆.sub[2], par₆.opt)
	tR₆ = sol₆.u[end][1]
	n_vec₆ = sort([n_vec₅; n₆])
	tR_vec₆ = sort([tR_vec₅; tR₆])
	tR_vec₆.-tR_lock
	itp₆ = Spline1D(tR_vec₆, n_vec₆, k=3)

	scatter!(p1, [tR₆, tR₆], [n₆, n₆], label="6")
	plot!(p1, 51.0:0.001:52.0, itp₆.(51.0:0.001:52.0), xlims=(51.47,51.485), ylims=(0.4995, 0.49975), label="itp₆")
end

# ╔═╡ 313d1431-ec39-4d1e-91fd-e025efc1f5c3
tR₆-tR_lock

# ╔═╡ fab59b1e-ba32-49d3-b0f1-564db037400c
begin
	n₇ = itp₆(tR_lock)
	par₇ = GasChromatographyTools.stretched_program(n₇, par0)
	sol₇ = GasChromatographySimulator.solving_odesystem_r(par₇.col, par₇.prog, par₇.sub[2], par₇.opt)
	tR₇ = sol₇.u[end][1]
	n_vec₇ = sort([n_vec₆; n₇])
	tR_vec₇ = sort([tR_vec₆; tR₇])
	tR_vec₇.-tR_lock
	itp₇ = Spline1D(tR_vec₇, n_vec₇, k=3)

	scatter!(p1, [tR₇, tR₇], [n₇, n₇], label="7")
	plot!(p1, 51.0:0.001:52.0, itp₇.(51.0:0.001:52.0), xlims=(51.47,51.485), ylims=(0.4995, 0.49975), label="itp₇")
end

# ╔═╡ 25725568-ae9a-4ab0-8f98-87fec12c867a
begin
	n₈ = itp₇(tR_lock)
	par₈ = GasChromatographyTools.stretched_program(n₈, par0)
	sol₈ = GasChromatographySimulator.solving_odesystem_r(par₈.col, par₈.prog, par₈.sub[2], par₈.opt)
	tR₈ = sol₈.u[end][1]
	n_vec₈ = sort([n_vec₇; n₈])
	tR_vec₈ = sort([tR_vec₇; tR₈])
	tR_vec₈.-tR_lock
	itp₈ = Spline1D(tR_vec₈, n_vec₈, k=3)

	scatter!(p1, [tR₈, tR₈], [n₈, n₈], label="8")
	plot!(p1, 51.0:0.001:52.0, itp₈.(51.0:0.001:52.0), xlims=(51.47,51.485), ylims=(0.4995, 0.49975), label="itp₈")
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
	sol₉ = GasChromatographySimulator.solving_odesystem_r(par₉.col, par₉.prog, par₉.sub[2], par₉.opt)
	tR₉ = sol₉.u[end][1]
	n_vec₉ = sort([n_vec₈; n₉])
	tR_vec₉ = sort([tR_vec₈; tR₉])
	tR_vec₉.-tR_lock
	itp₉ = Spline1D(tR_vec₉, n_vec₉, k=3)

	scatter!(p1, [tR₉, tR₉], [n₉, n₉], label="9")
	plot!(p1, 51.0:0.001:52.0, itp₉.(51.0:0.001:52.0), xlims=(51.47,51.485), ylims=(0.4995, 0.49975), label="itp₉")
end

# ╔═╡ eba919a1-4c2c-4eaf-ba12-2c02a8d0d6fb
n_vec₉

# ╔═╡ ae803a26-c436-497a-b245-a0afec92e46f
md"""
Make simulations for a range of stretch factors ``n``.
"""

# ╔═╡ 1983258e-e84b-4f39-9ce8-0e20e78a0893
begin
	nn = [0.49954, 0.49956, 0.49958, 0.49960, 0.49962, 0.49964, 0.49966]
	ttR = Array{Float64}(undef, length(nn))
	for i=1:length(nn)
	    pars = stretched_program(nn[i], par0)
	    sols = GasChromatographySimulator.solving_odesystem_r(pars.col, pars.prog, pars.sub[2], pars.opt)
	    ttR[i] = sols.u[end][1]
	end
	p2 = plot([tR_lock, tR_lock], [0.4, 0.6], label="tR_lock")
	plot!(p2, ttR, nn, line=(2,:solid), markers=:square, xlims=(51.47,51.49), ylims=(0.4995, 0.49975))	
end

# ╔═╡ 89cbd5c0-0df6-43d3-bd5c-d3d55d669f33
begin
	nnn = 0.499625.+collect(0.0000001:0.0000001:0.0000019)
	tttR = Array{Float64}(undef, length(nnn))
	for i=1:length(nnn)
	    pars = stretched_program(nnn[i], par0)
	    sols = GasChromatographySimulator.solving_odesystem_r(pars.col, pars.prog, pars.sub[2], pars.opt)
	    tttR[i] = sols.u[end][1]
	end
	plot!(p2, tttR, nnn, markers=:circle, ylims=(0.499625, 0.499627), xlims=(51.482, 51.485))
end

# ╔═╡ f541b210-3ef8-4003-b358-065e5c5949ad
par8 = stretched_program(nnn[8], par0)

# ╔═╡ 524b256b-0c1e-48bc-ac0b-2ee7fc1eab2b
par9 = stretched_program(nnn[9], par0)

# ╔═╡ 6d15c2a4-6d1e-4b89-b9e8-ecff004c4730
sol8 = GasChromatographySimulator.solving_odesystem_r(par8.col, par8.prog, par8.sub[2], par8.opt)

# ╔═╡ 665f302d-204b-47ac-90df-c5979350707c
sol9 = GasChromatographySimulator.solving_odesystem_r(par9.col, par9.prog, par9.sub[2], par9.opt)

# ╔═╡ 56095d71-6169-44b3-89d1-7ea7f1b6ddfb
sol8.destats

# ╔═╡ a51cd1dc-50bb-4444-a0f7-c0e4229b1257
sol9.destats

# ╔═╡ 225434d5-c753-4476-8f68-f5761a454852
sol8.retcode

# ╔═╡ 9a6f436e-e50e-4697-a1d1-3d2d43fc62fc
sol9.retcode

# ╔═╡ d165d269-58a9-4ad1-b5b2-6704ade60936
tR8 = sol8.u[end][1]

# ╔═╡ 61171239-d6f3-4e1b-9334-349633ee61a6
tR9 = sol9.u[end][1]

# ╔═╡ 3a759d13-d67a-4e9b-a166-dbd4eeb7b61b
tR8-tR9

# ╔═╡ 6040a096-38bd-4789-98d3-a02166e73f49
(tR8-tR9)/tR8

# ╔═╡ 08a0972a-55ff-42c9-838b-1a649afe9a46
md"""
## Decreased relative tolerance

The problem is not originated in the RT_lock algorithm but in the simulation itself. For certain settings a small change (like small strectch in the program) can result in bigger changes in retention time.
"""

# ╔═╡ 2d7da55e-2711-44a4-b8b9-bda6321a4c48
begin
	# decrease the relative tolerance
	opt_1 = GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-4, "inlet", true)
	par_1 = GasChromatographySimulator.Parameters(col, prog0, sub, opt_1)
	# repeat the simulation from above
	nnnn_1 = sort!(rand(0.499625:0.000000001:0.499627, 100))
	ttttR_1 = Array{Float64}(undef, length(nnnn_1))
	for i=1:length(nnnn_1)
	    pars = stretched_program(nnnn_1[i], par_1)
	    sols = GasChromatographySimulator.solving_odesystem_r(pars.col, pars.prog, pars.sub[2], pars.opt)
	    ttttR_1[i] = sols.u[end][1]
	end
	plot!(p2, ttttR_1, nnnn_1, ylims=(0.499625, 0.499627), xlims=(51.482, 51.485), markers=:v, label="reltol=1e-4")
end

# ╔═╡ 544b9b37-e047-49c1-b0d4-c58b4fde741c
Statistics.mean(ttttR_1)

# ╔═╡ b4a880cb-6138-432a-be88-a556908ab207
Statistics.std(ttttR_1)

# ╔═╡ ccd65c3a-2e6a-4540-9d82-e2bc23ef7230
begin
	# decrease the relative tolerance
	opt_2 = GasChromatographySimulator.Options(OwrenZen5(), 1e-7, 1e-4, "inlet", true)
	par_2 = GasChromatographySimulator.Parameters(col, prog0, sub, opt_2)
	# repeat the simulation from above
	nnnn_2 = sort!(rand(0.499625:0.000000001:0.499627, 100))
	ttttR_2 = Array{Float64}(undef, length(nnnn_2))
	for i=1:length(nnnn_2)
	    pars = stretched_program(nnnn_2[i], par_2)
	    sols = GasChromatographySimulator.solving_odesystem_r(pars.col, pars.prog, pars.sub[2], pars.opt)
	    ttttR_2[i] = sols.u[end][1]
	end
	plot!(p2, ttttR_2, nnnn_2, ylims=(0.499625, 0.499627), xlims=(51.482, 51.485), markers=:v, label="abstol=1e-7, reltol=1e-4")
end

# ╔═╡ 70d47e99-fcee-42bf-b502-b311bd347904
Statistics.mean(ttttR_2)

# ╔═╡ ef410f01-1bdb-4e1f-9e09-e1d309d1c546
Statistics.std(ttttR_2)

# ╔═╡ f012e638-ecdb-41f6-8cfa-16a43b503970
begin
	# decrease the relative tolerance
	opt_3 = GasChromatographySimulator.Options(OwrenZen5(), 1e-7, 1e-5, "inlet", true)
	par_3 = GasChromatographySimulator.Parameters(col, prog0, sub, opt_3)
	# repeat the simulation from above
	nnnn_3 = sort!(rand(0.499625:0.000000001:0.499627, 100))
	ttttR_3 = Array{Float64}(undef, length(nnnn_3))
	for i=1:length(nnnn_3)
	    pars = stretched_program(nnnn_3[i], par_3)
	    sols = GasChromatographySimulator.solving_odesystem_r(pars.col, pars.prog, pars.sub[2], pars.opt)
	    ttttR_3[i] = sols.u[end][1]
	end
	plot!(p2, ttttR_3, nnnn_3, ylims=(0.499625, 0.499627), xlims=(51.482, 51.485), markers=:v, label="abstol=1e-7, reltol=1e-5")
end

# ╔═╡ 012230e2-33f1-4658-8f5c-067b3df23a28
Statistics.mean(ttttR_3)

# ╔═╡ 834bfb61-b3cf-4ff9-b041-78b397b254b0
Statistics.std(ttttR_3)

# ╔═╡ 0aeb3ece-93a4-4638-bd90-1fdf13231e2c
md"""
**This is a good test inviroment for variations of abstol and reltol.**

It seems, that, if the result is 'noisy' (which is not always the case), than we have two branches of results for retention times for variing program stretch factors n. 

-> new notebook in GasChromatographySimulator.jl
"""

# ╔═╡ 89fa3139-1d10-4dd3-90b1-d9c0f65745c6
begin
	n_3 = RT_locking(par_3, tR_lock, 1e-3, last_alkane; opt_itp="spline")
	par_3_s = stretched_program(n_3, par_3)
	sol_3_s = GasChromatographySimulator.solving_odesystem_r(par_3_s.col, par_3_s.prog, par_3_s.sub[2], par_3_s.opt)
	tR_3_s = sol_3_s.u[end][1]
	#scatter!(p2, [tR_1_s, tR_1_s], [n_1, n_1], ylims=(0.9699739, n_1*1.000001))
	tR_3_s-tR_lock
end

# ╔═╡ e9accaac-587f-498d-ab0c-a8064f573678
begin
	# simulation arround n_1
	n_range = (n_3-0.0001):0.000001:(n_3+0.0001)
	tR_range = Array{Float64}(undef, length(n_range))
	for i=1:length(n_range)
	    pars = stretched_program(n_range[i], par_3)
	    sols = GasChromatographySimulator.solving_odesystem_r(pars.col, pars.prog, pars.sub[2], pars.opt)
	    tR_range[i] = sols.u[end][1]
	end
	#plot!(p2, tR_range, n_range, ylims=(n_range[1], n_range[end]), xlims=(tR_range[1], tR_range[end]))
end

# ╔═╡ a0e85c4b-7332-42a6-8fa0-77a9bb0c89a8
md"""
# Conclusion
The restriction has to be more strict, than only tR\_tol<reltol. Properly tR\_tol<1000*abstol (tR\_tol=1e-3 with abstol=1e-7).
"""

# ╔═╡ f38f4866-0f1b-4a41-a370-42dcb0c9e907
function RT_locking_mod(par::GasChromatographySimulator.Parameters, tR_lock::Float64, tR_tol::Float64, solute_RT::String; opt_itp="linear")
		# estimate the factor 'n' for the temperature program to achieve the retention time 'tR_lock' for 'solute_RT' with the GC-system defined by 'par' 
		if isa(par, Array)==true
			error("Select an element of the array of GC-systems.")
			elseif tR_tol<=par.opt.reltol || tR_tol<=1000*par.opt.abstol
		error("The tolerance for retention time locking tR_tol is to small. Use tR_tol > par.opt.reltol ($(par.opt.reltol)), resp. tR_tol > 1000*par.opt.abstol ($(1000*par.opt.abstol)).")
		else
			# find 'solute_RT' in the substances of 'par'
			name = Array{String}(undef, length(par.sub))
			for i=1:length(par.sub)
				name[i] = par.sub[i].name
			end
			ii = findfirst(name.==solute_RT)
			# calculate the retention time for the original (un-stretched) program 
			sol₀ = GasChromatographySimulator.solving_odesystem_r(par.col, par.prog, par.sub[ii], par.opt)
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

# ╔═╡ dfd1adb5-4681-446f-9c2e-25f89e56f6fc
RT_locking_mod(par_1, tR_lock, tR_tol, last_alkane; opt_itp="spline")

# ╔═╡ 420cff3b-f302-4751-8c6c-c11a6a45e2a5
RT_locking_mod(par_2, tR_lock, tR_tol, last_alkane; opt_itp="spline")

# ╔═╡ e37c9975-6c51-46a8-8695-c633566a841b
RT_locking_mod(par_3, tR_lock, tR_tol, last_alkane; opt_itp="spline")

# ╔═╡ abade442-bb18-406f-b9f7-6c6f6f04c515
RT_locking_mod(par_1, tR_lock, tR_tol*10, last_alkane; opt_itp="spline")

# ╔═╡ cfe9d6ec-5b22-46b5-b30f-ceed51d49644
function recur_RT_locking_rel(n::Float64, n_vec::Array{Float64,1}, tR_vec::Array{Float64,1}, par::GasChromatographySimulator.Parameters, tR_lock::Float64, tR_tol::Float64, ii::Int64; opt_itp="linear")
		# recursive function to find the factor 'n' for the temperature program to achieve the retention time 'tR_lock' for solute index 'ii' with the GC-system defined by 'par'
		# calculate the retention time with the input guess 'n'
		par₁ = stretched_program(n, par)
		if par.opt.odesys==true
			sol₁ = GasChromatographySimulator.solving_odesystem_r(par₁.col, par₁.prog, par₁.sub[ii], par₁.opt)
			tR₁ = sol₁.u[end][1]
		else
			sol₁ = GasChromatographySimulator.solving_migration(par₁.col, par₁.prog, par₁.sub[ii], par₁.opt)
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

# ╔═╡ 56a2001a-b09b-4a52-a889-c6cb599be242
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
			sol₀ = GasChromatographySimulator.solving_odesystem_r(par.col, par.prog, par.sub[ii], par.opt)
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

# ╔═╡ bb863e5c-5830-41a0-929a-001b15df7606
n_rel = RT_locking_rel(par0, tR_lock, 1e-4, last_alkane; opt_itp="spline")

# ╔═╡ 49c332ae-ae5f-4625-a693-b698ea678569
par0_rel = stretched_program(n_rel, par0)

# ╔═╡ 9356daef-2045-4391-afe8-2f67bb646dea
sol_rel = GasChromatographySimulator.solving_odesystem_r(par0_rel.col, par0_rel.prog, par0_rel.sub[2], par0_rel.opt)

# ╔═╡ 1422fb36-1933-4074-aa2c-556b9814fc42
	tR_rel = sol_rel.u[end][1]

# ╔═╡ da60b623-4c87-4798-899c-ec6e0b991cf4
tR_rel-tR_lock

# ╔═╡ 6dc84c22-042e-42ac-aab6-d048f879e228
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
# ╠═eba919a1-4c2c-4eaf-ba12-2c02a8d0d6fb
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
# ╠═d165d269-58a9-4ad1-b5b2-6704ade60936
# ╠═61171239-d6f3-4e1b-9334-349633ee61a6
# ╠═3a759d13-d67a-4e9b-a166-dbd4eeb7b61b
# ╠═6040a096-38bd-4789-98d3-a02166e73f49
# ╟─08a0972a-55ff-42c9-838b-1a649afe9a46
# ╠═2d7da55e-2711-44a4-b8b9-bda6321a4c48
# ╠═d3c2c2f4-1100-438d-88d6-d7502bcb9310
# ╠═544b9b37-e047-49c1-b0d4-c58b4fde741c
# ╠═b4a880cb-6138-432a-be88-a556908ab207
# ╠═ccd65c3a-2e6a-4540-9d82-e2bc23ef7230
# ╠═70d47e99-fcee-42bf-b502-b311bd347904
# ╠═ef410f01-1bdb-4e1f-9e09-e1d309d1c546
# ╠═f012e638-ecdb-41f6-8cfa-16a43b503970
# ╠═012230e2-33f1-4658-8f5c-067b3df23a28
# ╠═834bfb61-b3cf-4ff9-b041-78b397b254b0
# ╠═0aeb3ece-93a4-4638-bd90-1fdf13231e2c
# ╠═89fa3139-1d10-4dd3-90b1-d9c0f65745c6
# ╠═e9accaac-587f-498d-ab0c-a8064f573678
# ╠═a0e85c4b-7332-42a6-8fa0-77a9bb0c89a8
# ╠═f38f4866-0f1b-4a41-a370-42dcb0c9e907
# ╠═dfd1adb5-4681-446f-9c2e-25f89e56f6fc
# ╠═420cff3b-f302-4751-8c6c-c11a6a45e2a5
# ╠═e37c9975-6c51-46a8-8695-c633566a841b
# ╠═abade442-bb18-406f-b9f7-6c6f6f04c515
# ╠═56a2001a-b09b-4a52-a889-c6cb599be242
# ╠═cfe9d6ec-5b22-46b5-b30f-ceed51d49644
# ╠═bb863e5c-5830-41a0-929a-001b15df7606
# ╠═49c332ae-ae5f-4625-a693-b698ea678569
# ╠═9356daef-2045-4391-afe8-2f67bb646dea
# ╠═1422fb36-1933-4074-aa2c-556b9814fc42
# ╠═da60b623-4c87-4798-899c-ec6e0b991cf4
# ╠═6dc84c22-042e-42ac-aab6-d048f879e228
# ╠═42a1fccd-e782-42fb-91de-18c8b33d507d
