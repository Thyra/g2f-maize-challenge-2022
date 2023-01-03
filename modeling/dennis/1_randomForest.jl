using DataFrames
using CSV
using DecisionTree
using Statistics

println("Reading CSV...")
d = DataFrame(CSV.File("../../processed_data/combined_mat.csv"))
println("  done!")
# Only take wt and sl columns for now (no pheno and no EC)
x = d[!, 796:ncol(d)]

# The following two are just for safety, they should actually not throw out any columns
x = x[:, sum.(ismissing, eachcol(x)) .== 0] # Remove cols with missing values
x = x[:, Not(names(x, AbstractString))] # Remove string variables (cannot work with them so far)

# Transform Strings to Categorical Variables
# transform!(x, names(x, AbstractString) .=> categorical, renamecols=false)

# Adjust yield to 15.5% moisture
y = d.Yield_Mg_ha - d.Yield_Mg_ha .* (d.Grain_Moisture .- 15.5)/100

import JSON

println("Reading JSON...")
split = JSON.parsefile("../../processed_data/train_test_split_v2.json")

rmses = DataFrame(Set_name = String[], RMSE = Float64[])

for(key, val) in split
    print(key)
    print(" ... ")
    x_train = x[convert(Vector{Int}, val["train"] .+ 1), :]
    y_train = y[convert(Vector{Int}, val["train"] .+ 1)]
    x_test  = x[convert(Vector{Int}, val["test"] .+ 1), :]
    y_test  = y[convert(Vector{Int}, val["test"] .+ 1)]
    
    # We're using R defaults here
    model = build_forest(y_train, Matrix(x_train), ceil(ncol(x)/3), 500, 0.632)
    preds = apply_forest(model, Matrix(x_test))

    rmse = sqrt(mean((y_test - preds) .^ 2))
    println(rmse)
    push!(rmses, [key, rsme])
    CSV.write("./1_rmses.csv", rmses)
end
println("Done.")
