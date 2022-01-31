##---begin-UI-functions-----------------------------------------------------------------------------
## functions defining PlutoUI widgets for Pluto notebooks
"""
    UI_System(sp)

Construct a combined PlutoUI widget for the settings of the GC system with then selectable stationary phases `sp`. 
	
# UI fields
* ``L``: column length in m.
* ``d``: column diameter in mm.
* ``d_f``: film thickness in µm.
* stat. phase: stationary phase of the column
* Gas: mobile phase
"""
function UI_System(sp)
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
    UI_Program()

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
function UI_Program()
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
    UI_Substance(sol)

Construct a combined PlutoUI widget for the settings of the substances separated
in the simulated GC system with the selectable substances `subs`. 
	
# UI fields
* Select Substances: Selection of the substances, which will be simulated.
* Injection time: Start time (in s) of the simulation. The same for all selected
  substances.
* Injection width: Peak width (in s) of all selected substances at the time of injection.
"""
function UI_Substance(sol)
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
    UI_Options()

Construct a combined PlutoUI widget for the settings of the options for the simulation.    
"""
function UI_Options()
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