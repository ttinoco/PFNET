
mutable struct Network

    alloc::Bool
    _c_net::Ptr{Void}

    function Network(a, ptr)
        this = new(a, ptr)
        finalizer(this, dealloc)
        this
    end
end

Network() = Network(true, ccall((:NET_new, libpfnet), Ptr{Void}, (Int,), 1))
Network(num_periods::Int) = Network(true, ccall((:NET_new, libpfnet), Ptr{Void}, (Int,), num_periods))
Network(ptr::Ptr{Void}) = Network(false, ptr)

show_components(net::Network) = ccall((:NET_show_components, libpfnet), Cstring, (Ptr{Void},), net._c_net)

function dealloc(net::Network)
    print("dealloc net")
    if net.alloc
        print("alloc true")
        ccall((:NET_del, libpfnet), Void, (Ptr{Void},), net._c_net)
    else
        print("alloc false")
    end
end

