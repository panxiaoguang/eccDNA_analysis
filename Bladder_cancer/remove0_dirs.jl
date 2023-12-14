

function remove_no_file(filePath::String)
    fs = readdir(filePath, join=true)
    for i in 1:lastindex(fs)
        fsSize = filesize(fs[i])
        if fsSize == 0
            rm(fs[i], recursive=true, force=true)
        end
    end
end

remove_no_file(ARGS[1])