
mutable struct Constraint

    alloc::Bool
    ptr::Ptr{Void}
    net::Network

    function Constraint(a, ptr, n)
        this = new(a, ptr, n)
        finalizer(this, dealloc)
        this
    end
end

function Constraint(name::String, net::Network)
    
    if name == "AC power balance"
        Constraint(true, ccall((:CONSTR_ACPF_new, libpfnet), Ptr{Void}, (Ptr{Void},), net.ptr), net)

    elseif name == "DC power balance"
        Constraint(true, ccall((:CONSTR_DCPF_new, libpfnet), Ptr{Void}, (Ptr{Void},), net.ptr), net)

    elseif name == "linearized AC power balance"
        Constraint(true, ccall((:CONSTR_LINPF_new, libpfnet), Ptr{Void}, (Ptr{Void},), net.ptr), net)

    elseif name == "variable fixing"
        Constraint(true, ccall((:CONSTR_FIX_new, libpfnet), Ptr{Void}, (Ptr{Void},), net.ptr), net)

    elseif name == "variable bounds"
        Constraint(true, ccall((:CONSTR_LBOUND_new, libpfnet), Ptr{Void}, (Ptr{Void},), net.ptr), net)

    elseif name == "generator active power participation"
        Constraint(true, ccall((:CONSTR_PAR_GEN_P_new, libpfnet), Ptr{Void}, (Ptr{Void},), net.ptr), net)

    elseif name == "PVPQ switching"
        Constraint(true, ccall((:CONSTR_PVPQ_SWITCHING_new, libpfnet), Ptr{Void}, (Ptr{Void},), net.ptr), net)

    elseif name == "generator ramp limits"
        Constraint(true, ccall((:CONSTR_GEN_RAMP_new, libpfnet), Ptr{Void}, (Ptr{Void},), net.ptr), net)

    elseif name == "voltage regulation by generators"
        Constraint(true, ccall((:CONSTR_REG_GEN_new, libpfnet), Ptr{Void}, (Ptr{Void},), net.ptr), net)

    elseif name == "voltage regulation by transformers"
        Constraint(true, ccall((:CONSTR_REG_TRAN_new, libpfnet), Ptr{Void}, (Ptr{Void},), net.ptr), net)

    elseif name == "voltage regulation by shunts"
        Constraint(true, ccall((:CONSTR_REG_SHUNT_new, libpfnet), Ptr{Void}, (Ptr{Void},), net.ptr), net)

    elseif name == "DC branch flow limits"
        Constraint(true, ccall((:CONSTR_DC_FLOW_LIM_new, libpfnet), Ptr{Void}, (Ptr{Void},), net.ptr), net)

    elseif name == "AC branch flow limits"
        Constraint(true, ccall((:CONSTR_AC_FLOW_LIM_new, libpfnet), Ptr{Void}, (Ptr{Void},), net.ptr), net)

    elseif name == "linearized AC branch flow limits"
        Constraint(true, ccall((:CONSTR_AC_LIN_FLOW_LIM_new, libpfnet), Ptr{Void}, (Ptr{Void},), net.ptr), net)

    elseif name == "battery dynamics"
        Constraint(true, ccall((:CONSTR_BAT_DYN_new, libpfnet), Ptr{Void}, (Ptr{Void},), net.ptr), net)

    elseif name == "load constant power factor"
        Constraint(true, ccall((:CONSTR_LOAD_PF_new, libpfnet), Ptr{Void}, (Ptr{Void},), net.ptr), net)
        
    else
        throw(ArgumentError("invalid constraint name"))
    end
end 

function dealloc(constr::Constraint)
    if constr.alloc
        ccall((:CONSTR_del, libpfnet), Void, (Ptr{Void},), constr.ptr)
    end
    constr.alloc = false
    constr.ptr = C_NULL
end
