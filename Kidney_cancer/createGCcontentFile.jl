using BioSequences
using FASTX

function cal(sequence::LongDNASeq)
    nGC = count(x -> ((x == DNA_G) || (x == DNA_C)), sequence)
    nAT = count(x -> ((x == DNA_A) || (x == DNA_T)), sequence)
    (nGC, nAT)
end

function getPosGCs(IO::IOStream, CHRseq::BioSequence, totalLen::Int64, Chr::SubString, POS::Int64, WINDOW::Vector{Int64}, THRESH::Int64, idx::SubString)
    print(IO, idx, "\t", Chr, "\t", POS, "\t")
    for window in WINDOW
        window = window % 2 == 0 ? window + 1 : window
        startPOS = POS - Int64(floor(window / 2))
        tailPOS = POS + Int64(floor(window / 2))
        tailPOS = tailPOS > totalLen ? totalLen : tailPOS
        startPOS = startPOS <= 0 ? 1 : startPOS
        gc, at = cal(CHRseq[startPOS:tailPOS])
        if gc + at > THRESH
            print(IO, round(gc / (gc + at); digits=6), "\t")
        else
            print(IO, "NA", "\t")
        end
    end
    print(IO, "\n")
end

function main(snpLoci::String, fasta::String)
    open("GC_correct.file.txt", "w") do io
        println("Loding reference genome!")
        Genome = open(FASTA.Reader, fasta, index=string(fasta, ".fai"))
        WINDOWS = Int64[25, 50, 100, 200, 500, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6]
        println(io, "\tChr\tPosition\t25bp\t50bp\t100bp\t200bp\t500bp\t1kb\t2kb\t5kb\t10kb\t20kb\t50kb\t100kb\t200kb\t500kb\t1Mb")
        chrMarker = ""
        seq = dna"NNNNN"
        totalLen = 0
        for line in eachline(snpLoci)
            if startswith(line, "\t")
                continue
            else
                idx, chr, pos = split(line, "\t")
                pos = parse(Int64, pos)
                if chr != chrMarker
                    seq = FASTA.sequence(Genome[chr])
                    totalLen = length(seq)
                    println("processing chromosome ", chr, "...")
                end
                chrMarker = chr
                getPosGCs(io, seq, totalLen, chr, pos, WINDOWS, 20, idx)
            end
        end
    end
end

main(ARGS[1], ARGS[2])