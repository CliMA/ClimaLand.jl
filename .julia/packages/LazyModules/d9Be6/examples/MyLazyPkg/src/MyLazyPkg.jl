module MyLazyPkg

export generate_data, draw_figure
using LazyModules
@lazy import Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80"

generate_data(n) = sin.(range(start=0, stop=5, length=n) .+ 0.1.*rand(n))
draw_figure(data) = Plots.plot(data, title="MyPkg Plot")

end
