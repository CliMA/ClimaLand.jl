# example of parameter value, from https://gmd.copernicus.org/articles/11/1909/2018/#&gid=1&pid=1

Pz = FT(101.0) # 1 [kPa] pressure just above the soil surface at time t
Cᵣ = FT(10.0) # 2 [mg C cm-3] Cᵣ is root biomass carbon, see figure S5
Csom = FT(5.0) # 3 [mg C cm-3] soil organic C content at depth z
Cmic = FT(1.0) # 4 [mg C cm-3] Microbial C pool, ~ 1 % of Csom at DETECT site

Rᵦ = FT(6e-5) # 5 [mg C cm-2] total root biomass C in a 1 m deep by 1 cm2 soil column
α₁ᵣ = FT(11.65) # 6 [-]
α₂ᵣ = FT(20.7) # 7 [-]
α₃ᵣ = FT(-164.2) # 8 [-]

Sᶜ = FT(711.6) # 9 [mg C cm-2] total soil organic C in a 1 m deep by 1 cm2 soil column
Mᶜ = FT(12.3) # 10 [mg C cm-2] total microbial biomass C in a 1 m deep by 1 cm2 soil column
Vᵦ = FT(0.0015) # 11 [mg C cm-3 h-1] value of Vmax at 10 °C and mean environmental conditions
α₁ₘ = FT(14.05) # 12 [-]
α₂ₘ = FT(11.05) # 13 [-]
α₃ₘ = FT(-87.6) # 14 [-]
Kₘ = FT(10e-5) # 15 [mg C cm-3 h-1] Michaelis-Menten half saturation constant
CUE = FT(0.8) # 16 [mg C mg-1 C-1] microbial carbon use efficiency
pf = FT(0.004) # 17 [-] fraction of soil organic C that is soluble
Dₗᵢ = FT(3.17) # 18 [-] Diffusivity of soil C substrate in liquid

E₀ₛ = FT(324.6) # 19 [Kelvin] temperature sensitivity parameter
T₀ = FT(227.5) # 20 [Kelvin] temperature sensitivity-related parameter
α₄ = FT(4.7) # 21 [-]
Tᵣₑ = FT(283.0) # 22 [Kelvin] ref temperature for other param e.g., Rᵦ

α₅ = FT(4.547) # 23 [-]
BD = FT(1.12) # 24 [g cm-3] soil bulk density
ϕ₁₀₀ = FT(18.46) # 25 [%] air filled porosity at soil water potential of -100 cm H2O (~ 10 kPa)
PD = FT(2.52) # 26 [g cm-3] soil particle density
Dstp = FT(1.39e-5) # 27 [m2 s-1] diffusion coefficient for CO2 in air at standard T and P
P₀ = FT(101.325) # 28 [kPa] standard pressure
b = FT(4.547) # 29 [-] parameter related to pore size distribution

