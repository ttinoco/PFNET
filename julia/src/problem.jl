
mutable struct Problem

    alloc::Bool
    ptr::Ptr{Void}
    functions::Array{Function,1}
    constraints::Array{Constraint,1}
    net::Network

    function Problem(a, ptr, f, c, n)
        this = new(a, ptr, f, c, n)
        finalizer(this, dealloc)
        this
    end
end

Problem(net::Network) = Problem(true,
                                ccall((:PROB_new, libpfnet), Ptr{Void}, (Ptr{Void},), net.ptr),
                                Function[],
                                Constraint[],
                                net)

function add_constraint(prob::Problem, constr::Constraint)
    constr.alloc = false
    push!(prob.constraints, constr)
    ccall((:PROB_add_constr, libpfnet), Void, (Ptr{Void}, Ptr{Void},), prob.ptr, constr.ptr)
end

function add_function(prob::Problem, func::Function)
    func.alloc = false
    push!(prob.functions, func)
    ccall((:PROB_add_func, libpfnet), Void, (Ptr{Void}, Ptr{Void},), prob.ptr, func.ptr)
end
    
analyze(prob::Problem) = ccall((:PROB_analyze, libpfnet), Void, (Ptr{Void},), prob.ptr)

function dealloc(prob::Problem)
    if prob.alloc
        ccall((:PROB_del, libpfnet), Void, (Ptr{Void},), prob.ptr)
    end
    prob.alloc = false
    prob.ptr = C_NULL
end

show_problem(prob::Problem) = ccall((:PROB_show, libpfnet), Void, (Ptr{Void},), prob.ptr)
