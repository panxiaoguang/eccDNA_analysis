function check_eofs(fs::String)
    open(fs, "r") do IO
        seekend(IO)
        eof_positions = position(IO) - 28
        seek(IO, eof_positions)
        eof_markers = UInt8[0x1f, 0x8b, 0x08, 0x04, 0x00,
            0x00, 0x00, 0x00, 0x00, 0xff,
            0x06, 0x00, 0x42, 0x43, 0x02,
            0x00, 0x1b, 0x00, 0x03, 0x00,
            0x00, 0x00, 0x00, 0x00, 0x00,
            0x00, 0x00, 0x00]
        actual_markers = Vector{UInt8}(undef, 28)
        read!(IO, actual_markers)
        if isequal(eof_markers, actual_markers)
            @info "file have BGZF EOF Markers and is corrected!" filename = fs
        else
            @warn "file maybe truncted,please check your file!" filename = fs
        end
    end
end

function main(ARGS)
    if length(ARGS) > 1
        check_eofs.(ARGS)
    else
        check_eofs(ARGS[1])
    end
end

main(ARGS)