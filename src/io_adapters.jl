module IOAdapters

using CSV
using DataFrames
using ..BinarizeDynamics: validate_sequences

export read_fasta, read_csv_sequences

"""
    read_fasta(path::String)

Reads a FASTA file and returns a tuple `(headers, sequences)`.
- `headers`: Vector of strings containing the sequence descriptions (without the `>`).
- `sequences`: Vector of strings containing the sequences.

Handles multi-line sequences by concatenating lines until the next header or end of file.
"""
function read_fasta(path::String)
    headers = String[]
    sequences = String[]
    
    current_header = nothing
    current_seq_buffer = IOBuffer()
    
    open(path, "r") do io
        for line in eachline(io)
            line = strip(line)
            if isempty(line)
                continue
            end
            
            if startswith(line, ">")
                # If we were processing a sequence, save it
                if current_header !== nothing
                    push!(headers, current_header)
                    push!(sequences, String(take!(current_seq_buffer)))
                end
                
                # Start new record
                current_header = String(strip(line[2:end]))
            else
                # Append to current sequence
                if current_header === nothing
                    # This might happen if file doesn't start with >
                    # For now, we can skip or error. Let's error strictly.
                    # Or be lenient and treat as raw sequences? 
                    # Standard FASTA requires header first.
                    throw(ArgumentError("FASTA file must start with a header line (>)"))
                end
                print(current_seq_buffer, line)
            end
        end
        
        # Save the last record
        if current_header !== nothing
            push!(headers, current_header)
            push!(sequences, String(take!(current_seq_buffer)))
        end
    end
    
    return headers, sequences
end

"""
    read_csv_sequences(path::String, col::Union{String, Int, Symbol})

Reads a CSV or Excel-exported CSV file and extracts sequences from the specified column.
Returns `sequences` as a Vector{String}.

Arguments:
- `path`: Path to the CSV file.
- `col`: Column name (String or Symbol) or index (Int) containing the sequences.
"""
function read_csv_sequences(path::String, col::Union{String, Int, Symbol})
    df = CSV.read(path, DataFrame)
    
    if col isa Int
        if col < 1 || col > ncol(df)
            throw(ArgumentError("Column index $col out of bounds for DataFrame with $(ncol(df)) columns"))
        end
        # Get column by index
        raw_seqs = df[!, col]
    else
        # validation handled by DataFrame indexing
        if !hasproperty(df, Symbol(col))
             throw(ArgumentError("Column '$col' not found in CSV headers: $(names(df))"))
        end
        raw_seqs = df[!, Symbol(col)]
    end
    
    # Ensure they are strings
    return String.(string.(raw_seqs))
end

end
