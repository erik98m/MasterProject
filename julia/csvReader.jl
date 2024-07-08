using CSV
using DataFrames
using Allocations

# Function to read the CSV and process it into additive matrices
function process_csv_to_additive_matrices(filename)
    iterations = []
    current_iteration = []

    # Open the file and read it line by line
    open(filename, "r") do file
        for line in eachline(file)
            if occursin("Iteration", line)
                if !isempty(current_iteration)
                    push!(iterations, current_iteration)
                    current_iteration = []
                end
            elseif !isempty(strip(line))
                push!(current_iteration, parse.(Int, split(line, ",")))
            end
        end
    end

    if !isempty(current_iteration)
        push!(iterations, current_iteration)
    end

    additive_matrices = []

    for iter in iterations
        matrix = reduce(vcat, (row' for row in iter))  # Create a matrix by vertically concatenating rows
        push!(additive_matrices, matrix)
    end

    return additive_matrices
end

# Function to process additive matrices and store results to a CSV file
function process_and_store_additive_matrices(matrices, output_filename)
    results = []

    for (i, matrix) in enumerate(matrices)
        println("Additive Matrix for Iteration ", i - 1)
        println(matrix)
        
        # Create the valuation profile using the Additive function from the Allocations library
        V = Additive(matrix)
        
        # Run alloc_mm(V).mm
        result = alloc_mm(V).mm
        
        # Store the result in the results array
        push!(results, result)
        
        # Print the result
        println("Result of alloc_mm(V).mm: ", result)
        println()
    end

    # Write results to CSV file
    results_df = DataFrame(Result = results)
    CSV.write(output_filename, results_df)
end

# Specify the CSV file path
input_filename = "/home/erik98m/masteroppgave/TS/eqmaxi/valprofile/10agents_valuations.csv"
output_filename = "pythonPlots/TS/eqmaxi/n10.csv"

# Process the CSV and get the additive matrices
#additive_matrices = process_csv_to_additive_matrices(input_filename)

# Process the additive matrices and store results to a CSV file
#process_and_store_additive_matrices(additive_matrices, output_filename)
