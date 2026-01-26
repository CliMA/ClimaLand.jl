module IfElse

@static if Base.ifelse === Core.ifelse
ifelse(args...) = Core.ifelse(args...)
else
const ifelse = Base.ifelse
end
end
