function reformat(fs::String)
    for line in eachline(fs)
        hao = split(line, "\t")
        if hao[9] == "genome-circle"
            newline = join([hao[1], hao[4], hao[5], hao[2], hao[3], hao[6], string(hao[8], hao[7]), "circle-genome"], "\t")
            println(newline)
        else
            newline = join([hao[1], hao[2], hao[3], hao[4], hao[5], hao[6], string(hao[7], hao[8]), hao[9]], "\t")
            println(newline)
        end
    end
end

reformat("test9.2.txt")