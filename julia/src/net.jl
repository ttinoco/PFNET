
mutable struct Network

    alloc::Bool
    ptr::Ptr{Void}

    function Network(a, ptr)
        this = new(a, ptr)
        finalizer(this, dealloc)
        this
    end
end

Network() = Network(true, ccall((:NET_new, libpfnet), Ptr{Void}, (Int,), 1))
Network(num_periods::Int) = Network(true, ccall((:NET_new, libpfnet), Ptr{Void}, (Int,), num_periods))
Network(ptr::Ptr{Void}) = Network(false, ptr)

clear_flags(net::Network) = ccall((:NET_clear_flags, libpfnet), Void, (Ptr{Void},), net.ptr)

function dealloc(net::Network)
    if net.alloc
        ccall((:NET_del, libpfnet), Void, (Ptr{Void},), net.ptr)
    end
    net.alloc = false
    net.ptr = C_NULL
end

num_bounded(net::Network) = ccall((:NET_get_num_bounded, libpfnet), Int, (Ptr{Void},), net.ptr)
num_buses(net::Network) = ccall((:NET_get_num_buses, libpfnet), Int, (Ptr{Void},), net.ptr)
num_generators_not_on_outage(net::Network) = ccall((:NET_get_num_gens_not_on_outage, libpfnet), Int, (Ptr{Void},), net.ptr)
num_slack_buses(net::Network) = ccall((:NET_get_num_slack_buses, libpfnet), Int, (Ptr{Void},), net.ptr)
num_periods(net::Network) = ccall((:NET_get_num_periods, libpfnet), Int, (Ptr{Void},), net.ptr)
num_vars(net::Network) = ccall((:NET_get_num_vars, libpfnet), Int, (Ptr{Void},), net.ptr)

#function set_flags(net::Network, 

show_components(net::Network) = ccall((:NET_show_components, libpfnet), Void, (Ptr{Void},), net.ptr)



