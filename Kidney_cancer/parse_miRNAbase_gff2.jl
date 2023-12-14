const miRNA_db_gff = "/home/panxiaoguang/下载/hsa.gff3"
function parse_miRNA(miRNAdb, outfile)
    f = open(outfile, "w")
    for line in eachline(miRNAdb)
        if startswith(line, "#")
            continue
        else
            chrom, source, type, start, stop, score, strand, phase, attributes = split(line, '\t')
            if type == "miRNA"
                name = split(attributes, ';')[3]
                name = split(name, '=')[2]
                derive = split(attributes, ';')[4]
                derive = split(derive, '=')[2]
                println(f, "$derive\t$name")
            end
        end
    end
    close(f)
end

parse_miRNA(miRNA_db_gff, "miRNA_duiying.txt")