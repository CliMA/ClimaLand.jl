# Getting Started with Matrix-Free Operators in Julia

SciMLOperators.jl is a package for defining operators for use in solvers.
One of the major use cases is to define matrix-free operators in cases where
using a matrix would be too memory expensive. In this tutorial we will walk
through the main features of SciMLOperators and get you going with matrix-free
and updating operators.

## Simplest Operator: MatrixOperator

Before we get into the deeper operators, let's show the simplest SciMLOperator:
`MatrixOperator`. `MatrixOperator` just turns a matrix into an `AbstractSciMLOperator`,
so it's not really a matrix-free operator but it's a starting point that is good for
understanding the interface and testing. To create a `MatrixOperator`, simply call the
constructor on a matrix:

```@example getting_started
using SciMLOperators, LinearAlgebra
A = [-2.0 1 0 0 0
     1 -2 1 0 0
     0 1 -2 1 0
     0 0 1 -2 1
     0 0 0 1 -2]

opA = MatrixOperator(A)
```

The operators can do [operations as defined in the operator interface](@ref operator_interface), for example,
matrix multiplication as the core action:

```@example getting_started
v = [3.0, 2.0, 1.0, 2.0, 3.0]
opA*v
```

```@example getting_started
opA(v, nothing, nothing, nothing) # Call = opA*v
```

```@example getting_started
w = zeros(5)
mul!(w, opA, v)
```

```@example getting_started
α = 1.0;
β = 1.0
mul!(w, opA, v, α, β) # α*opA*v + β*w
```

and the inverse operation:

```@example getting_started
opA \ v
```

```@example getting_started
ldiv!(w, lu(opA), v)
```

## State, Parameter, and Time-Dependent Operators

Now let's define a `MatrixOperator` the is dependent on state, parameters, and time.
For example, let's make the operator `A .* u + dt*I` where `dt` is a parameter
and `u` is a state vector:

```@example getting_started
A = [-2.0 1 0 0 0
     1 -2 1 0 0
     0 1 -2 1 0
     0 0 1 -2 1
     0 0 0 1 -2]

function update_function!(B, u, p, t)
    dt = p
    B .= A .* u + dt*I
end

u = Array(1:1.0:5);
p = 0.1;
t = 0.0
opB = MatrixOperator(copy(A); update_func! = update_function!)
```

To update the operator, you would use `update_coefficients!(opB, u, p, t)`:

```@example getting_started
update_coefficients!(opB, u, p, t)
```

We can use the interface to see what the current matrix is by converting to a standard matrix:

```@example getting_started
convert(AbstractMatrix, opB)
```

And now applying the operator applies the updated one:

```@example getting_started
opB*v
```

Or if you use the operator application, it will update and apply in one step:

```@example getting_started
opB(v, Array(2:1.0:6), 0.5, nothing) # opB(u,p,t)*v
```

This is how for example, when an ODE solver asks for an operator `L(u,p,t)*u`, this is how
such an operator can be defined. Notice that the interface can be queried to understand
the traits of the operator, such as for example whether an operator is constant (does not
change w.r.t. `(u,p,t)`):

```@example getting_started
isconstant(opA)
```

```@example getting_started
isconstant(opB)
```

## Matrix-Free Operators via FunctionOperator

Now let's define the operators from above in a matrix-free way using `FunctionOperator`.
With `FunctionOperator`, we directly define the operator application function `opA(w,v,u,p,t)`
which means `w = opA(u,p,t)*v`. For example we can do the following:

```@example getting_started
function Afunc!(w, v, u, p, t)
    w[1] = -2v[1] + v[2]
    for i in 2:4
        w[i] = v[i - 1] - 2v[i] + v[i + 1]
    end
    w[5] = v[4] - 2v[5]
    nothing
end

function Afunc!(v, u, p, t)
    w = zeros(5)
    Afunc!(w, v, u, p, t)
    w
end

mfopA = FunctionOperator(Afunc!, zeros(5), zeros(5))
```

Now `mfopA` acts just like `A*v` and thus `opA`:

```@example getting_started
mfopA*v - opA*v
```

```@example getting_started
mfopA(v, u, p, t) - opA(v, u, p, t)
```

We can also create the state-dependent operator as well:

```@example getting_started
function Bfunc!(w, v, u, p, t)
    dt = p
    w[1] = -(2*u[1]-dt)*v[1] + v[2]*u[1]
    for i in 2:4
        w[i] = v[i - 1]*u[i] - (2*u[i]-dt)*v[i] + v[i + 1]*u[i]
    end
    w[5] = v[4]*u[5] - (2*u[5]-dt)*v[5]
    nothing
end

function Bfunc!(v, u, p, t)
    w = zeros(5)
    Bfunc!(w, v, u, p, t)
    w
end

mfopB = FunctionOperator(Bfunc!, zeros(5), zeros(5); u, p, t, isconstant = false)
```

```@example getting_started
opB(v, Array(2:1.0:6), 0.5, nothing) - mfopB(v, Array(2:1.0:6), 0.5, nothing)
```

## Operator Algebras

While the operators are lazy operations and thus are not full matrices, you can still
do algebra on operators. This will construct a new lazy operator that will be able to
compute the same action as the composed function. For example, let's create `mfopB`
using `mfopA`. Recall that we defined this via `A .* u + dt*I`. Let's first create an
operator for `A .* u` (since right now there is not a built in operator for vector scaling,
but that would be a fantastic thing to add!):

```@example getting_started
function Cfunc!(w, v, u, p, t)
    w[1] = -2v[1] + v[2]
    for i in 2:4
        w[i] = v[i - 1] - 2v[i] + v[i + 1]
    end
    w[5] = v[4] - 2v[5]
    w .= w .* u
    nothing
end

function Cfunc!(v, u, p, t)
    w = zeros(5)
    Cfunc!(w, v, u, p, t)
    w
end

mfopC = FunctionOperator(Cfunc!, zeros(5), zeros(5))
```

And now let's create the operator `mfopC + dt*I`. We can just directly build it:

```@example getting_started
mfopD = mfopC + 0.5*I
```

SciMLOperators.jl uses an `IdentityOperator` and `ScalarOperator` instead of the Base
utilities, but the final composed operator acts just like the operator that was built:

```@example getting_started
mfopB(v, Array(2:1.0:6), 0.5, nothing) - mfopD(v, Array(2:1.0:6), 0.5, nothing)
```

There are many cool things you can do with operator algebras, such as `kron` (Kronecker products),
adjoints, inverses, and more. For more information, see the [operator algebras tutorial](@ref operator_algebras).

## Where to go next?

Great! You now know how to be state/parameter/time-dependent operators and make them matrix-free, along with
doing algebras on operators. What's next?

  - Interested in more examples of building operators? See the example of [making a fast fourier transform linear operator](@ref fft)
  - Interested in more operators ready to go? See the [Premade Operators page](@ref premade_operators) for all of the operators included with SciMLOperators. Note that there are also downstream packages that make new operators.
  - Want to make your own SciMLOperator? See the [AbstractSciMLOperator interface page](@ref operator_interface) which describes the full interface.

How do you use SciMLOperators? Check out the following downstream pages:

  - [Using SciMLOperators in LinearSolve.jl for matrix-free Krylov methods](https://docs.sciml.ai/LinearSolve/stable/tutorials/linear/)
  - [Using SciMLOperators in OrdinaryDiffEq.jl for semi-linear ODE solvers](https://docs.sciml.ai/DiffEqDocs/stable/solvers/nonautonomous_linear_ode/)
