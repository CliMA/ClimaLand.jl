# Graphics.jl

Graphics.jl is an [abstraction layer](https://en.wikipedia.org/wiki/Abstraction_layer)
for graphical operations in Julia.
Its goal is to allow developers to write graphical programs in a manner independent
of the particular graphical backend.
One needs to load a backend package that implements the operations in its API;
currently, [Cairo.jl](https://github.com/JuliaGraphics/Cairo.jl)
is the only such backend.

To get an organized overview of the API, try typing `?Graphics` at the Julia REPL.
You can see the same information in greater detail on the next page.
