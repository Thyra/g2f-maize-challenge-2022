using DataFrames
using CSV
using Evolutionary
using DecisionTree
using Statistics
import StatsBase

println("Reading pheno CSV...")
d = DataFrame(CSV.File("../../processed_data/combined_mat.csv"))
println("  done!")
println("Reading geno CSV...")
g = DataFrame(CSV.File("../../processed_data/geno_processed.miss.1.mac.1.biallelic.txt"))
println("  done!")
println("Joining them...")
d = leftjoin(d, g, on = :Hybrid)
println("  done!")
# Only take wt and sl columns for now (no pheno and no EC) + genotypic data
x = d[!, 796:ncol(d)]

# The following two are just for safety, they should actually not throw out any columns
x = x[:, sum.(ismissing, eachcol(x)) .== 0] # Remove cols with missing values
x = x[:, Not(names(x, AbstractString))] # Remove string variables (cannot work with them so far)

# Transform Strings to Categorical Variables
# transform!(x, names(x, AbstractString) .=> categorical, renamecols=false)

# Adjust yield to 15.5% moisture
y = d.Yield_Mg_ha - d.Yield_Mg_ha .* (d.Grain_Moisture .- 15.5)/100

# Free some space
d = g = nothing

test_rows = StatsBase.sample(1:nrow(x), ceil(Int, nrow(x)*0.3), replace=false)
train_rows = setdiff(1:nrow(x), test_rows)

x_train = x[train_rows, :]
x_test  = x[test_rows, :]
y_train = y[train_rows]
y_test  = y[test_rows]

# Now for the feature selection thingy

# features_active = BitVector with the length = number of columns in x; 0 for this column is not used, 1 for this column is used.
# We want to optimize the columns that are used (i.e. feature selection), that means the resulting rmse should be as low as possible
function evaluate_feature_subset(features_active) 
    x_train_active = x_train[:, features_active]
    
    model = build_forest(y_train, Matrix(x_train_active), ceil(ncol(x_train_active)/3), 200, 0.632)
    preds = apply_forest(model, Matrix(x_test))
    rmse = sqrt(mean((y_test - preds) .^ 2))
    return(rmse)
end

# Starting individual
fa = BitVector(zeros(ncol(x_train)))
for k in StatsBase.sample(1:ncol(x_train), 30)
    fa[k] = 1
end

res = Evolutionary.optimize(x->evaluate_feature_subset(x),
                            fa,
                            GA(selection=susinv, populationSize=120,
                            mutation=flip, crossover=BSX(1000), ɛ=10),
                            Evolutionary.Options(time_limit=86400.0, show_trace=true, parallelization = :thread))

# @TODO to save features after every iteration, do something with the trace https://github.com/wildart/Evolutionary.jl/issues/104

open("selected_features.txt", "w") do f
    for n in names(x)[res.minimizer]
        write(f, n * "\n")
    end
end

