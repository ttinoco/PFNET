
mutable struct Function

    alloc::Bool
    ptr::Ptr{Void}
    net::Network

    function Function(a, ptr, n)
        this = new(a, ptr, n)
        finalizer(this, dealloc)
        this
    end
end

function Function(name::String, weight::Float64, net::Network)
    
    if name == "generation cost"
        Function(true, ccall((:FUNC_GEN_COST_new, libpfnet), Ptr{Void}, (Float64, Ptr{Void},), weight, net.ptr), net)

    elseif name == "consumption utility"
        Function(true, ccall((:FUNC_LOAD_util_new, libpfnet), Ptr{Void}, (Float64, Ptr{Void},), weight, net.ptr), net)

    elseif name == "net consumption cost"
        Function(true, ccall((:FUNC_NETCON_COST_new, libpfnet), Ptr{Void}, (Float64, Ptr{Void},), weight, net.ptr), net)

    elseif name == "phase shift regularization"
        Function(true, ccall((:FUNC_REG_PHASE_new, libpfnet), Ptr{Void}, (Float64, Ptr{Void},), weight, net.ptr), net)

    elseif name == "generator powers regularization"
        Function(true, ccall((:FUNC_REG_PQ_new, libpfnet), Ptr{Void}, (Float64, Ptr{Void},), weight, net.ptr), net)

    elseif name == "tap ratio regularization"
        Function(true, ccall((:FUNC_REG_RATIO_new, libpfnet), Ptr{Void}, (Float64, Ptr{Void},), weight, net.ptr), net)

    elseif name == "susceptance regularization"
        Function(true, ccall((:FUNC_REG_SUSC_new, libpfnet), Ptr{Void}, (Float64, Ptr{Void},), weight, net.ptr), net)

    elseif name == "voltage angle regularization"
        Function(true, ccall((:FUNC_REG_VANG_new, libpfnet), Ptr{Void}, (Float64, Ptr{Void},), weight, net.ptr), net)

    elseif name == "voltage magnitude regularization"
        Function(true, ccall((:FUNC_REG_VMAG_new, libpfnet), Ptr{Void}, (Float64, Ptr{Void},), weight, net.ptr), net)

    elseif name == "variable regularization"
        Function(true, ccall((:FUNC_REG_VAR_new, libpfnet), Ptr{Void}, (Float64, Ptr{Void},), weight, net.ptr), net)

    elseif name == "soft voltage magnitude limts"
        Function(true, ccall((:FUNC_SLIM_VMAG_new, libpfnet), Ptr{Void}, (Float64, Ptr{Void},), weight, net.ptr), net)

    elseif name == "sparse control penalty"
        Function(true, ccall((:FUNC_SP_CONTROLs_new, libpfnet), Ptr{Void}, (Float64, Ptr{Void},), weight, net.ptr), net)

    else
        throw(ArgumentError("invalid function name"))
    end
end 

function dealloc(func::Function)
    if func.alloc
        ccall((:FUNC_del, libpfnet), Void, (Ptr{Void},), func.ptr)
    end
    func.alloc = false
    func.ptr = C_NULL
end
