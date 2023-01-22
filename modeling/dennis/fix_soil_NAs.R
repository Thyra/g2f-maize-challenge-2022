library(stringr)
# To combine i just removed columns and then did rbind
co = read.csv("soils_combined.csv")

co$loc = str_split_fixed(co$Env, "_", 2)[,1]
for(i in co$Env) {
  row = co[co$Env == i,]
  loc.envs = co[co$loc == row$loc,]
  for(j in names(co)) {
    if(is.na(co[co$Env == i, j])) {
      if(nrow(loc.envs) > 1) {
        # For locations where data is available from another year, just use that!
        co[co$Env == i, j] = median(loc.envs[, j], na.rm=T)
      } else {
        # Otherwise use the overall median to have something
        co[co$Env == i, j] = median(co[, j], na.rm=T)
      }
      print(paste(i, j, co[co$Env == i, j]))
    }
  }
}
# Remove the loc column again
co = co[1:(ncol(co)-1)]

write.csv(co, "soils_combined_fixed.csv", row.names = F)
