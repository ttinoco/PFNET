
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

num_periods(net::Network) = ccall((:NET_get_num_periods, libpfnet), Int, (Ptr{Void},), net.ptr)
show_components(net::Network) = ccall((:NET_show_components, libpfnet), Void, (Ptr{Void},), net.ptr)

function dealloc(net::Network)
    if net.alloc
        ccall((:NET_del, libpfnet), Void, (Ptr{Void},), net.ptr)
    end
    net.alloc = false
    net.ptr = C_NULL
end

