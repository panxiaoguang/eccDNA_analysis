using GenomicFeatures
using Glob
function get_intervals(file::String)
    col = IntervalCollection{Vector{SubString{String}}}()
    for line in eachline(file)
        shuju = split(line, "\t")
        push!(col, Interval(shuju[1], parse(Int64, shuju[2]), parse(Int64, shuju[3]), '.', shuju))
    end
    col
end

function get_overlap(case::IntervalCollection{Vector{SubString{String}}}, db::IntervalCollection{Vector{SubString{String}}}, IO::IOStream)
    for (a, b) in eachoverlap(case, db)
        if leftposition(b) <= leftposition(a) <= rightposition(b)
            println(IO, join(GenomicFeatures.metadata(a), "\t"), "\t", join(GenomicFeatures.metadata(b), "\t"))
        end
    end
end



function Anno(sample::String)
    case = get_intervals("filter_bed/$(sample).ecc.bed")
    dbs = get_intervals("/home/panxiaoguang/Project/Bladder/dbs/hg38.genetic_elements.merge.bed")
    w = open("element_anno/$(sample).startAnno.bed", "w+")
    get_overlap(case, dbs, w)
    close(w)
end

##for sample in Glob.glob("filter_bed/*.ecc.bed")
#    sample = replace(sample, "filter_bed/" => "", ".ecc.bed" => "")
#    Anno(sample)
#end

for sample in ["filter_bed/cKca_141T.ecc.bed","filter_bed/cKca_141N.ecc.bed","filter_bed/cKca_142T.ecc.bed","filter_bed/cKca_142N.ecc.bed"]
    sample = replace(sample, "filter_bed/" => "", ".ecc.bed" => "")
    Anno(sample)
end
