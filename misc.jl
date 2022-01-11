using GasChromatographyTools
using GasChromatographySimulator

# Options
opt = GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-3, "inlet", true)
GasChromatographySimulator.Options()
# System
L = 10.0
a_d = [0.25e-3]
a_df = [0.25e-6]
d(x) = GasChromatographySimulator.gf_const(x, a_d)
df(x) = GasChromatographySimulator.gf_const(x, a_df)
sp = "SLB5ms"
gas = "He"
sys = GasChromatographySimulator.constructor_System(L, a_d[1], a_df[1], sp, gas)

# Program
time_steps = [0.0, 60.0, 600.0, 300.0]
temp_steps = [40.0, 40.0, 300.0, 300.0]
pin_steps = [200.0, 200.0, 300.0, 300.0].*1000.0 .+ 101300
pout_steps = [101.3, 101.3, 101.3, 101.3].*1000.0
ΔT_steps = zeros(length(time_steps))
x₀_steps = zeros(length(time_steps))
L₀_steps = L.*ones(length(time_steps))
α_steps = zeros(length(time_steps))
a_gf = [ΔT_steps x₀_steps L₀_steps α_steps]
gf(x) = GasChromatographySimulator.gradient(x, a_gf)
T_itp = GasChromatographySimulator.temperature_interpolation(time_steps, temp_steps, gf, sys.L)
pin_itp = GasChromatographySimulator.pressure_interpolation(time_steps, pin_steps)
pout_itp = GasChromatographySimulator.pressure_interpolation(time_steps, pout_steps)
prog =  GasChromatographySimulator.constructor_Program(time_steps, temp_steps, pin_steps, pout_steps, ΔT_steps, x₀_steps, L₀_steps, α_steps, opt.Tcontrol, L)

# Substance
sub = Array{GasChromatographySimulator.Substance}(undef, 2)
sub[1] = GasChromatographySimulator.Substance("solute 1",
                                                "0-0-0",
                                                400.0,
                                                30.0,
                                                100.0,
                                                1e-3,
                                                "test",
                                                1e-4,
                                                0.0,
                                                0.0)
sub[2] = GasChromatographySimulator.Substance("solute 2",
                                                "1-0-0",
                                                420.0,
                                                31.0,
                                                100.0,
                                                1e-3,
                                                "test",
                                                1e-4,
                                                0.0,
                                                0.0)

# Parameters
par = GasChromatographySimulator.Parameters(sys, prog, sub, opt)

#----test RT-locking functions--------------------------------
sol_n1 = GasChromatographySimulator.solve_system_multithreads(par)
tR_lock_n1 = sol_n1[2].u[end][1]

stretched_prog = GasChromatographyTools.stretched_program(2.0, par)
sol_n2 = GasChromatographySimulator.solve_system_multithreads(stretched_prog)
tR_lock_n2 = sol_n2[2].u[end][1]

n_factor = GasChromatographyTools.RT_locking(par, 300.0, 1e-3, "solute 2")

#---test nb_solutes_RTlock.jl
opt_C30 = GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-3, "outlet", true)
sys_C30 = GasChromatographySimulator.constructor_System(4.0, 0.1e-3, 0.1e-6, "Rxi5MS", "He")
prog_C30 =  GasChromatographySimulator.constructor_Program([0.0, 10.0, 60.0, 120.0], [50.1518, 50.1518, 310.152, 310.152], [401300.0, 401300.0, 576300.0, 576300.0], [0.0, 0.0, 0.0, 0.0], [20.0, 20.0, 20.0, 20.0], [0.0, 0.0, 0.0, 0.0], [4.0, 4.0, 4.0, 4.0], [0.0, 0.0, 0.0, 0.0], opt_C30.Tcontrol, sys_C30.L)
sub_C30 = [GasChromatographySimulator.Substance("C30",
                                                "638-68-6",
                                                600.477,
                                                37.5137,
                                                181.004,
                                                1e-3,
                                                "Gaida.2021",
                                                5.17619e-5,
                                                0.0,
                                                0.0)]
par_C30 = GasChromatographySimulator.Parameters(sys_C30, prog_C30, sub_C30, opt_C30)

tR_lock_C30 = 60.0
tR_tol_C30 = 1e-3
n = GasChromatographyTools.RT_locking(par_C30, tR_lock_C30, tR_tol_C30, "C30"; opt_itp="linear")
n2 = GasChromatographyTools.RT_locking(par_C30, tR_lock_C30, tR_tol_C30, "C30"; opt_itp="spline")