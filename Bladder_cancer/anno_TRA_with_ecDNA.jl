using GenomicFeatures

function get_intervals(file::String)
    bkp1 = IntervalCollection{Vector{SubString{String}}}()
    bkp2 = IntervalCollection{Vector{SubString{String}}}()
    for line in eachline(file)
        if startswith(line, "chrom1")
            continue
        else
            shuju = split(line, "\t")
            if shuju[11] == "TRA"
                push!(bkp1, Interval(shuju[1], parse(Int64, shuju[2]), parse(Int64, shuju[2]) + 2, '.', shuju[1:end-1]))
                push!(bkp2, Interval(shuju[4], parse(Int64, shuju[5]), parse(Int64, shuju[5]) + 2, '.', shuju[1:end-1]))
            end
        end
    end
    (bkp1, bkp2)
end

function get_database(file::String)
    col = IntervalCollection{Vector{SubString{String}}}()
    for line in eachline(file)
        shuju = split(line, "\t")
        push!(col, Interval(shuju[1], parse(Int64, shuju[2]), parse(Int64, shuju[3]), '.', shuju))
    end
    col
end


function get_overlap(case::IntervalCollection{Vector{SubString{String}}}, db::IntervalCollection{Vector{SubString{String}}})
    marks = Vector{String}()
    for (a, b) in eachoverlap(case, db)
        if leftposition(b) <= leftposition(a) <= rightposition(b)
            push!(marks, string(GenomicFeatures.seqname(a), ":", string(GenomicFeatures.leftposition(a))))
        end
    end
    marks
end


function Anno()
    open("test9.2.txt", "w+") do IO
        (bkp1, bkp2) = get_intervals("SV/9.merged3.sv.txt")
        dbs = get_database("tmp_shuju/9.ecDNA.bed")
        bkp1_s = get_overlap(bkp1, dbs)
        bkp2_s = get_overlap(bkp2, dbs)
        for line in eachline("SV/9.merged3.sv.txt")
            if startswith(line, "chrom1")
                continue
            else
                hao = split(line, "\t")
                if hao[11] == "TRA"
                    need1 = string(hao[2], ":", hao[3])
                    need2 = string(hao[4], ":", hao[5])
                    if need1 in bkp1_s && need2 in bkp2_s
                        println(IO, join(["9", hao[1], hao[2], hao[4], hao[5], hao[7], hao[9], hao[10], "circle-circle"], "\t"))
                    elseif need1 in bkp1_s && need2 ∉ bkp2_s
                        println(IO, join(["9", hao[1], hao[2], hao[4], hao[5], hao[7], hao[9], hao[10], "circle-genome"], "\t"))
                    elseif need1 ∉ bkp1_s && need2 in bkp2_s
                        println(IO, join(["9", hao[1], hao[2], hao[4], hao[5], hao[7], hao[9], hao[10], "genome-circle"], "\t"))
                    elseif need1 ∉ bkp1_s && need2 ∉ bkp2_s
                        println(IO, join(["9", hao[1], hao[2], hao[4], hao[5], hao[7], hao[9], hao[10], "genome-genome"], "\t"))
                    end
                end
            end
        end
    end
end

Anno()


