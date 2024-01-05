
  #= GLMakie layout
  ax3D = Axis3(fig[5,2], xlabel = inputs.drivers.names[1], ylabel = inputs.drivers.names[2])
  ax_d1 = Axis(fig[6,1], xlabel = inputs.drivers.names[1])
  ax_d2 = Axis(fig[6,2], xlabel = inputs.drivers.names[2])
  =#


  #= GLMakie slider  
  sliders = [Slider(fig[i, 1], range = -5:1:5, startvalue = 0) for i in 1:4]
  s_vals = [sliders[i].value for i in 1:4]   
  s1_d = Slider(fig[1, 1], range = -5:1:5, startvalue = 0) # driver
  s2_d = Slider(fig[2, 1], range = -5:1:5, startvalue = 0)
  s1_p = Slider(fig[3, 1], range = -5:1:5, startvalue = 0) # parameter
  s2_p = Slider(fig[4, 1], range = -5:1:5, startvalue = 0)
  s1_d_v = s1_d.value
  s2_d_v = s2_d.value
  s1_p_v = s1_p.value
  s2_p_v = s2_p.value
  =#


  # Update parameters and rescale x and y limits  
  #=
  on(s1_p_v) do val 
    autolimits!(ax3D)
    autolimits!(ax_d1)
    autolimits!(ax_d2)
  end

  on(s2_p_v) do val 
    autolimits!(ax3D)
    autolimits!(ax_d1)
    autolimits!(ax_d2)
  end

  on(s1_d_v) do val
    autolimits!(ax_d1)
    autolimits!(ax_d2)
  end

  on(s2_d_v) do val
    autolimits!(ax_d1)
    autolimits!(ax_d2)
  end
  =#

