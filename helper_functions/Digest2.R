#### modified version of OrgMassSpecR::Digest
#### - delete functionality to calculate peptide masses & enzymes other than trypsin
#### - interpret "missed" agrument as maximum number of allowed missed cleavages and 
####   only warn if nr or missed cleavages is not possible

Digest2 <- function (sequence, enzyme = "trypsin", missed = 0, IAA = TRUE, 
          N15 = FALSE, custom = list()) {
  seq_vector <- strsplit(sequence, split = "")[[1]]
  end_position <- length(seq_vector)
  if (enzyme == "trypsin") {
    if (seq_vector[end_position] == "K" | seq_vector[end_position] == 
        "R") {
      seq_vector[end_position] <- "!"
      seq_string <- paste(seq_vector, collapse = "")
    }
    else seq_string <- sequence
    seq_string <- gsub("KP", "!P", seq_string)
    seq_string <- gsub("RP", "!P", seq_string)
    seq_vector <- strsplit(seq_string, split = "")[[1]]
    stop <- grep("K|R", seq_vector)
    start <- stop + 1
  }
  if (enzyme == "trypsin.strict") {
    if (seq_vector[end_position] == "K" | seq_vector[end_position] ==
        "R") {
      seq_vector[end_position] <- "!"
      seq_string <- paste(seq_vector, collapse = "")
    }
    else seq_string <- sequence
    seq_vector <- strsplit(seq_string, split = "")[[1]]
    stop <- grep("K|R", seq_vector)
    start <- stop + 1
  }
  if (enzyme != "trypsin" & enzyme != "trypsin.strict") 
    stop("undefined enzyme, defined enzymes are trypsin, trypsin.strict")
  if (length(stop) == 0) {
    warning("sequence does not contain cleavage sites")
    return(data.frame(sequence = sequence, start = 1, stop = nchar(sequence), mc = 0))
  }
  
  if (missed > length(stop)) {
    warning("number of specified missed cleavages is greater than the maximum possible")
  }
  
  cleave <- function(sequence, start, stop, misses) {
    peptide <- substring(sequence, start, stop)
    mc <- rep(misses, times = length(peptide))
    result <- data.frame(sequence = peptide, start, stop, mc, stringsAsFactors = FALSE)
    return(result)
  }
  stop_ <- stop
  start <- c(1, start)
  stop <- c(stop, end_position)
  results <- cleave(sequence, start, stop, 0)
  if (missed > 0) {
    for (i in 1:min(missed, length(stop_))) {
      start_tmp <- start[1:(length(start) - i)]
      stop_tmp <- stop[(1 + i):length(stop)]
      peptide <- cleave(sequence, start_tmp, stop_tmp, i)
      results <- rbind(results, peptide)
    }
  }
   return(results)
}


# License: BSD_2_clause + file LICENSE
# Copyright (c) 2011-2017, Nathan Dodder
#   
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
#   Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
# 
#   Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in
# the documentation and/or other materials provided with the
# distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.