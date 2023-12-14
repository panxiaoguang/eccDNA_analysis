using Glob
function parse_svaba(fs::String, IO::IOStream, cutoff::Int64)
    svtype = ""
    POS2 = 0
    pesupport = 0
    strand = ""
    HEADER = "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsv_id\tpe_support\tstrand1\tstrand2\tsvclass\tsvmethod"
    println(IO, HEADER)
    for line in eachline(fs)
        if startswith(line, "#")
            continue
        else
            CHROM, POS, ID, ref, alt, _, _, INFO, _, _, _ = split(line, "\t")
            POS = parse(Int64, POS)
            sv_length = abs(length(ref) - length(alt))
            if sv_length > cutoff
                if length(ref) > length(alt)
                    svtype = "DEL"
                    POS2 = POS + sv_length
                    strand = "+\t-"
                elseif length(ref) < length(alt)
                    svtype = "INS/DUP"
                    POS2 = POS
                    strand = "-\t+"
                end
                println(IO, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", strand, "\t", svtype, "\t", "Svaba")
            end
        end
    end
end

function main()
    cutoff = 50
    for fs in Glob.glob("SVABA/*.svaba.smallIndel.vcf")
        prx = replace(fs, ".svaba.smallIndel.vcf" => "", "SVABA/" => "")
        open("format_svaba/format.$(prx).smallIndel.svaba.txt", "w+") do IO
            parse_svaba(fs, IO, cutoff)
        end
    end
end

main()