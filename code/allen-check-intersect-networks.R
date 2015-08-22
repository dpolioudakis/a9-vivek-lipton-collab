# Interesect gene modules in networks made with different parameters

load("../processed_data/allen_BW_modules.rda")

# softPower 7 minModSize 30 deepSplit 2 MEmergeCutHeight 0.25 maxBlockSize 12000
MinMod30 <- split(colnames(exprData), bwModulesLL[[11]]$colors)

#softPower 9 minModSize 100 deepSplit 2 MEmergeCutHeight 0.25 maxBlockSize 12000
MinMod100 <- split(colnames(exprData), bwModulesLL[[18]]$colors)

toMatch <- MinMod30$plum1
toMatch <- MinMod30$grey60
toMatch <- MinMod30$brown

toMatch <- MinMod30$cyan
toMatch <- MinMod30$green
toMatch <- MinMod30$sienna3
toMatch <- MinMod30$royalblue

sapply(MinMod100, function (x) toMatch[toMatch %in% x])

