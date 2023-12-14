const miRNA_db_gff = "/home/panxiaoguang/下载/hsa.gff3"
function parse_primary(miRNAdb, outfile)
    f = open(outfile, "w")
    for line in eachline(miRNAdb)
        if startswith(line, "#")
            continue
        else
            chrom, source, type, start, stop, score, strand, phase, attributes = split(line, '\t')
            if type == "miRNA_primary_transcript"
                id = split(attributes, ';')[1]
                id = split(id, '=')[2]
                name = split(attributes, ';')[3]
                name = split(name, '=')[2]
                println(f, join([chrom, start, stop, strand, id, name], '\t'))
            end
        end
    end
    close(f)
end

parse_primary(miRNA_db_gff, "miRNA_primary_transcript.bed")