library('textreuse')
library('NameNeedle')
library('readr')

musMSH2 <- read_file('MSH2Mus.txt')
humMSH2 <- read_file('MSH2Homo.txt')


# myParam <- defaultNeedleParams
# myParam$MATCH <- 2
# myParam$MISMATCH <- -2
# 
# scores <- needleScores(musMSH2, humMSH2, myParam)
# ## This gene is too large to use
# 
# GamaDMus <- read_file('GamaDMus.txt')
# GamaDHomo <- read_file('GamaDHomo.txt')
# 
# 
# myParam <- defaultNeedleParams
# myParam$MATCH <- 2
# myParam$MISMATCH <- -2
# 
# scores <- needleScores(GamaDHomo, GamaDMus, myParam)
# ## This gene is too large to use
# 
# 
# scores
# 
# 
# result <- align_local(GamaDMus, GamaDHomo)
# 
# str(result)

library(Biostrings)
s1 <- DNAString(GamaDMus)
s2 <- DNAString(GamaDHomo)
mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
globalAlign <-
  pairwiseAlignment(s1, s2, substitutionMatrix = mat,
                    gapOpening = 5, gapExtension = 2)
localAlign <-
  pairwiseAlignment(s1, s2, type = "local", substitutionMatrix = mat,
                    gapOpening = 5, gapExtension = 2)


GamaDMus <- read_file('GamaDMus.mRNA.txt')
GamaDHomo <- read_file('GamaDHomo.mRNA.txt')

s3 <- DNAString(GamaDMus)
s4 <- DNAString(GamaDHomo)

mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
globalAlign <-
  pairwiseAlignment(s3, s4, substitutionMatrix = mat,
                    gapOpening = 5, gapExtension = 2)
localAlign <-
  pairwiseAlignment(s3, s4, type = "local", substitutionMatrix = mat,
                    gapOpening = 5, gapExtension = 2)


GamaDMus <- read_file('GamaDMus.Pro.txt')
GamaDHomo <- read_file('GamaDHomo.Pro.txt')

s5 <- AAString(GamaDMus)
s6 <- AAString(GamaDHomo)

globalAlign <-
  pairwiseAlignment(s5, s6, substitutionMatrix = 'BLOSUM62',
                    gapOpening = 5, gapExtension = 2)
localAlign <-
  pairwiseAlignment(s5, s6, type = "local", substitutionMatrix = 'BLOSUM62',
                    gapOpening = 5, gapExtension = 2)




