using CSV
using DataFrames

### merge part 
eachCt(x::Int64, co::Vector{Int64})::Int64 = count(haha -> haha == x, co)

function getDUP(x::Vector{Int64})::Vector{Int64}
    dups = Int64[]
    for ni in unique(x)
        if eachCt(ni, x) > 1
            push!(dups, ni)
        end
    end
    dups
end

function mg!(haha::Vector{Vector{Int64}})
    oneline = vcat(haha...)
    rongyu = getDUP(oneline)
    for kk in rongyu
        ind = findall(x -> (kk in x), haha)
        if length(ind) > 0
            push!(haha, union(haha[ind]...))
            deleteat!(haha, ind)
        end
    end
    haha
end
################################

function mergeSV(sample::String, IO)
    delly = DataFrame(CSV.File("format_delly/format.$(sample).sv.delly.txt", delim="\t"))
    novo = DataFrame(CSV.File("format_novo/format.$(sample).sv.novobreak.txt", delim="\t"))
    manta = DataFrame(CSV.File("format_manta/format.$(sample).sv.manta.txt", delim="\t"))
    svaba = DataFrame(CSV.File("format_svaba/format.$(sample).allsvs.svaba.txt", delim="\t"))
    alls = foldl((x, y) -> vcat(x, y), [delly, novo, manta, svaba])
    sort!(alls, [:svclass, :chrom1, :start1])
    alls.label = collect(1:nrow(alls))
    #alls |> CSV.write("$(sample).merge.original.txt", delim="\t")
    shunxu = Dict("Manta" => 1, "Svaba" => 2, "Delly" => 3, "NovoBreak" => 4)
    HEADER = "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsv_id\tpe_support\tstrand1\tstrand2\tsvclass\tsvmethod\torignal"
    println(IO, HEADER)
    waiting_merge = Vector{Vector{Int64}}()
    for row in eachrow(alls)
        panxuan = (alls.chrom1 .== row.chrom1) .& (alls.chrom2 .== row.chrom2) .& (row.start1 - 100 .< alls.start1 .< row.start1 + 100) .& (row.start2 - 100 .< alls.start2 .< row.start2 + 100) .& (alls.svclass .== row.svclass) .& (alls.strand1 .== row.strand1) .& (alls.strand2 .== row.strand2)
        need = alls[panxuan, :]
        if length(need.svmethod) > 1
            push!(waiting_merge, need.label)
        else
            println(IO, row.chrom1, "\t", row.start1, "\t", row.end1, "\t", row.chrom2, "\t", row.start2, "\t", row.end2, "\t", row.sv_id, "\t", row.pe_support, "\t", row.strand1, "\t", row.strand2, "\t", row.svclass, "\t", row.svmethod, "\t", row.label)
        end
    end
    mg!(waiting_merge)
    for rst in waiting_merge
        houhou = alls[rst, :]
        houhou = sort(houhou, order(:svmethod, by=(x) -> shunxu[x]))
        houhou.mix = fill(join(unique(houhou.svmethod), ","), nrow(houhou))
        houhou.labels = fill(join(unique(string.(houhou.label)), ","), nrow(houhou))
        woyao = houhou[1, :]
        println(IO, woyao.chrom1, "\t", woyao.start1, "\t", woyao.end1, "\t", woyao.chrom2, "\t", woyao.start2, "\t", woyao.end2, "\t", woyao.sv_id, "\t", woyao.pe_support, "\t", woyao.strand1, "\t", woyao.strand2, "\t", woyao.svclass, "\t", woyao.mix, "\t", woyao.labels)
    end
end

function main()
    for sm in eachline("samples")
        open("merged/$(sm).merged.allsv.txt", "w+") do io
            mergeSV(sm, io)
        end
    end
end

main()
