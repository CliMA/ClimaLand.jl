module MyPkg

export generate_data, draw_figure
import Plots

generate_data(n) = sin.(range(start=0, stop=5, length=n) .+ 0.1.*rand(n))
draw_figure(data) = Plots.plot(data, title="MyPkg Plot")

end
