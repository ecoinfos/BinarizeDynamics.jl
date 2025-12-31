module TestIO

using Test
using BinarizeDynamics
using DataFrames
using CSV

@testset "IO Adapters" begin
    # Create temp files
    fasta_content = """
>seq1
AAAA
TTTT
>seq2
CCCC
GGGG
"""
    fasta_path = "test_seqs.fasta"
    write(fasta_path, fasta_content)
    
    csv_content = """
id,seq
1,AAAA
2,TTTT
3,CCCC
"""
    csv_path = "test_seqs.csv"
    write(csv_path, csv_content)
    
    try
        @testset "FASTA Parsing" begin
            headers, seqs = read_fasta(fasta_path)
            @test length(headers) == 2
            @test length(seqs) == 2
            @test headers[1] == "seq1"
            @test seqs[1] == "AAAATTTT"
            @test headers[2] == "seq2"
            @test seqs[2] == "CCCCGGGG"
        end
        
        @testset "CSV Parsing" begin
            # Test by column name
            seqs_csv = read_csv_sequences(csv_path, "seq")
            @test length(seqs_csv) == 3
            @test seqs_csv[1] == "AAAA"
            
            # Test by column index
            seqs_idx = read_csv_sequences(csv_path, 2)
            @test seqs_idx == seqs_csv
            
            # Test missing column
            @test_throws  ArgumentError read_csv_sequences(csv_path, "nonexistent")
            @test_throws  ArgumentError read_csv_sequences(csv_path, 5)
        end
        
    finally
        rm(fasta_path, force=true)
        rm(csv_path, force=true)
    end
end

end # module
