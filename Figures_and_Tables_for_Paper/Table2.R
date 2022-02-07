################################################################################
### durch Isoformen neu entstehende Peptide und Uniqueness

library(limma)

M_iso <- readRDS("data/D3/D3_fasta/preprocessed/01_matrix_D3_fasta_min7AA.rds")
M <- readRDS("data/D3_without_isoforms/D3_fasta/preprocessed/01_matrix_D3_without_isoforms_fasta_min7AA.rds")

pep_iso <- rownames(M_iso)
pep <- rownames(M)


pep_new <- setdiff(pep_iso, pep)
ind_pep_new <- which(pep_iso %in% pep_new)
x <- rowSums(M_iso[ind_pep_new,])  ## wie viele der neuen Peptide sind unique bzw shared?
#     1     2     3     4     5     6     7     8     9    10    12    13    19
# 93856  7953  1455   515    97    80    42    29     5     5     1     5     3

ind_pep_new_shared



rs_M <- rowSums(M)
ind_unique <- pep[which(rs_M == 1)]
ind_shared <- pep[which(rs_M > 1)]


rs_M_iso <- rowSums(M_iso)
ind_unique_iso <- pep_iso[which(rs_M_iso == 1)]
ind_shared_iso <- pep_iso[which(rs_M_iso > 1)]


### isoform-unique Peptide: d.h. Peptide, die zwar geshart sind zwischen verschiedenen
### accessions, die aber jeweils nur Isoformen voneinander sind

unique_iso <- logical(length(ind_shared_iso))
for (i in seq_along(ind_shared_iso)) {
  print(i)
  j <- ind_shared_iso[i]
  ind_tmp <- which(M_iso[j,] == 1)
  x <- names(ind_tmp)
  unique_iso[i] <- (length(unique(limma::strsplit2(x, "-")[,1])) == 1)
}


ind_isounique_iso <- ind_shared_iso[unique_iso]
ind_shared2_iso <- ind_shared_iso[!unique_iso]



# M_iso2 <- M_iso
# colnames(M_iso2) <- limma::strsplit2(colnames(M_iso2), "-")[,1]
# table(duplicated(colnames(M_iso2)))
# rs_M_iso2 <- rowSums(M_iso2[, !duplicated(colnames(M_iso2))])
# ind_shared_iso2 <- pep_iso[which(rs_M_iso2 > 1)]
# ind_unique_iso2 <- pep_iso[which(rs_M_iso2 == 1)]



### Wie viele Peptide bleiben geshared bzw unique, wenn Isoformen dazugenommen werden?
(unique_unique <- length(intersect(ind_unique, ind_unique_iso))) # 1218900
(unique_shared <- length(intersect(ind_unique, ind_shared2_iso))) # 26027
(unique_isounique <- length(intersect(ind_unique, ind_isounique_iso))) # 485723

(shared_unique <- length(intersect(ind_shared, ind_unique_iso))) # 0
(shared_shared <- length(intersect(ind_shared, ind_shared2_iso))) # 1319690
(shared_isounique <- length(intersect(ind_shared, ind_isounique_iso))) # 0

(new_unique <- length(intersect(pep_new, ind_unique_iso))) # 93856
(new_unique <- length(intersect(pep_new, ind_shared2_iso))) #  103
(new_unique <- length(intersect(pep_new, ind_isounique_iso))) # 10087


