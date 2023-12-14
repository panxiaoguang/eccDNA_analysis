using GenomicFeatures


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


function Anno(input::String,dbFile::String,sample::String,outdir::String)
    case = get_intervals(input)
    dbs = get_intervals(dbFile)
    w = open("$(outdir)/$(sample).startAnno.bed", "w+")
    get_overlap(case, dbs, w)
    close(w)
end


Anno(ARGS[1],ARGS[2],ARGS[3],ARGS[4])

