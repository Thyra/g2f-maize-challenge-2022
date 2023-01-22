using DataFrames
using CSV
using DecisionTree
using Statistics

println("Reading pheno and soil CSV...")
d = DataFrame(CSV.File("../../processed_data/combined_mat_w_o_BLUEs_final_v2.csv"))
s = DataFrame(CSV.File("soils_combined_fixed.csv"))
d = leftjoin(d, s, on = :Env)
println("  done!")
println("Reading geno CSV...")
g = DataFrame(CSV.File("../../processed_data/geno_processed.miss.1.mac.1.biallelic.txt"))
# @TODO for now: only use selected genetic features
 g = g[:, readlines("top_markers2.txt")]
println("  done!")
println("Joining them...")
d = leftjoin(d, g, on = :Hybrid)
println("  done!")

x = d[:, 5:ncol(d)]
x_train = x[x.type .== "train", : ]
x_sub = x[x.type .== "submission", : ]

# The following two are just for safety, they should actually not throw out any columns
x_train = x_train[:, sum.(ismissing, eachcol(x_train)) .== 0] # Remove cols with missing values
x_train = x_train[:, Not(names(x_train, AbstractString))] # Remove string variables (cannot work with them so far)

x_sub = x_sub[:, Not(names(x_sub, AbstractString))] # Remove string variables (cannot work with them so far)

labels_sub = d[d.type .== "submission", ["Env", "Hybrid"]]

y = d[d.type .== "train", :Yield_Mg_ha]
y = convert(Vector{Float64}, y)

model = build_forest(y, Matrix(x_train), ceil(ncol(x_train)/3), 500, 0.632)
preds = apply_forest(model, Matrix(x_sub))

labels_sub.Yield_Mg_ha = preds
CSV.write("./SUBMISSION.csv", labels_sub)
