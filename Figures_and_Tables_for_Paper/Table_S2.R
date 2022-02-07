
### Numbers for Table S2


D1 <- read.table("data/D1/D1_quant/preprocessed/FC_NA_mind2.txt", sep = "\t", header = TRUE)

## quantified peptides
nrow(D1) # 22939

## mean number of valid ratios
mean(colSums(!is.na(D1[, -c(1:2)]))) # 18092.8




D2 <- read.table("data/D2_without_isoforms/D2_quant/preprocessed/FC_NA_mind2.txt", sep = "\t", header = TRUE)

## quantified peptides
nrow(D2)  # 8101

## mean number of valid ratios
mean(colSums(!is.na(D2[, -c(1:2)])))  # 5618.639




D3 <- read.table("data/D3_without_isoforms/D3_quant/preprocessed/FC_NA_mind2.txt", sep = "\t", header = TRUE)

## quantified peptides
nrow(D3)  # 46544

## mean number of valid ratios
sum(!is.na(D3[, -c(1:2)]))  # 30369


D3_iso <- read.table("data/D3/D3_quant/preprocessed/FC_NA_mind2.txt", sep = "\t", header = TRUE)

## quantified peptides
nrow(D3_iso)  # 46544

## mean number of valid ratios
sum(!is.na(D3_iso[, -c(1:2)]))  # 29490

