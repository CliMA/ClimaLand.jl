#####
##### 1st order
#####
import SciMLBase

# TODO: is it better to use `first_order_tableau = Union{ARS111,ARS121}`? to
#       reduce the number of methods?
first_order_tableau() = [ARS111, ARS121]

#####
##### 2nd order
#####

second_order_tableau() = [
    ARS122,
    ARS232,
    ARS222,
    IMKG232a,
    IMKG232b,
    IMKG242a,
    IMKG242b,
    IMKG243a,
    IMKG252a,
    IMKG252b,
    IMKG253a,
    IMKG253b,
    IMKG254a,
    IMKG254b,
    IMKG254c,
    HOMMEM1,
    SSP222,
    SSP322,
    SSP332,
    SSPKnoth,
]

#####
##### 3rd order
#####
third_order_tableau() = [ARS233, ARS343, ARS443, IMKG342a, IMKG343a, DBM453, SSP333, SSP433]

for m in first_order_tableau()
    @eval SciMLBase.alg_order(::$m) = 1
end
for m in second_order_tableau()
    @eval SciMLBase.alg_order(::$m) = 2
end
for m in third_order_tableau()
    @eval SciMLBase.alg_order(::$m) = 3
end
