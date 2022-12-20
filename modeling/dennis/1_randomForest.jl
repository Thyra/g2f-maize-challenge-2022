using DataFrames
using CSV
using DecisionTree
using Statistics

println("Reading CSV...")
d = DataFrame(CSV.File("../../processed_data/combined_mat.csv"))
println("  done!")
x = d[:, 12:ncol(d)]
x = x[:, sum.(ismissing, eachcol(x)) .== 0]

# Remove string variables (cannot work with them so far)
x = x[!, Not(names(x, AbstractString))]

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
end
CSV.write("./1_rmses.csv", rmses)
rmses
