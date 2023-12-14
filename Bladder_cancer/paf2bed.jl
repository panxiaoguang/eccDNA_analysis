function paf2bed()
    for line in eachline("alignData/5T.align.paf")
        infors = split(line, "\t")
        println(infors[6], "\t", infors[8], "\t", infors[9], "\t", infors[1])
    end
end

paf2bed()