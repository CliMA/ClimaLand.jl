using IfElse
using Test

x = 2
@test IfElse.ifelse(x>0,1,-1) == Core.ifelse(x>0,1,-1)
