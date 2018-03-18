
mutable struct Parser

    alloc::Bool
    ptr::Ptr{Void}

    function Parser(a, ptr)
        this = new(a, ptr)
        finalizer(this, dealloc)
        this
    end
end

Parser(ext::String) = Parser(true, ccall((:PARSER_new_for_file, libpfnet), Ptr{Void}, (Cstring,), ext))

parse_case(parser::Parser, filename::String; num_periods::Int=1) = Network(ccall((:PARSER_parse, libpfnet),
                                                                                 Ptr{Void},
                                                                                 (Ptr{Void}, Cstring, Int,),
                                                                                 parser.ptr,
                                                                                 filename,
                                                                                 num_periods))

function dealloc(parser::Parser)
    if parser.alloc
        ccall((:PARSER_del, libpfnet), Void, (Ptr{Void},), parser.ptr)
    end
    parser.alloc = false
    parser.ptr = C_NULL
end
