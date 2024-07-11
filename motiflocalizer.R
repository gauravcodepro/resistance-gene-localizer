#! usr/bin/env python3
# Author Gaurav 

library(stringr)
suppressPackageStartupMessages(library(Biostrings, pos = "package:base"))
motiflocalizer <- function(fasta_files, motif_pattern){
    # a exact motif localizer based on the string pattern
    # searching and return a exact location of the string
    #matches of those motifs in your sequences. Below are
    #some of the examples of the protein domains of the NB-ARC
    #plant resistance genes. You can give the complete domains
    # of a specific part of it. This is an application of exact 
    # match
    #NB-ARC example
    #>domain_seq
    #DLCAGGAVGVAFNELFVVLKHVIKTIAGFKSTFNRLEATMLAISPIFEDIKR
    #VSKLAAVITQVQNKTRHFKSHLSQIQETLTQINPGYKESEKRNHDLGRGRMQATEVFIGQLKEAEELVRKCEHIS
    #VVEGAALGALFQVLFNSVVNASFRLSQFRSQLDLLKSTLSFIKPAIDDVENLDRFLERPLSETQ
    #AVGNSRIKEGVECVRKCKSVGWYNIKKKFWYIQKLEELDRTLELLIPFAYSKRDIKEIFIRTRKTS
    # > motiflocalizer(fasta_files = "./motif_analyze.fasta", motif_pattern = "ATAT")
    #m64020e_230609_061545/36768059/ccs.start m64020e_230609_061545/36768059/ccs.end
    #                                 143                                    146
    #m64020e_230609_061545/114820129/ccs.start m64020e_230609_061545/114820129/ccs.end
    #                                143                                     146
    #m64020e_230609_061545/114820428/ccs.start m64020e_230609_061545/114820428/ccs.end
    #                                166                                     169.
    fasta = readDNAStringSet(fasta_files)
    pattern = as.character(motif_pattern)
    fasta_names = names(fasta)
    fasta_string = vector(length = length(names(fasta)))
    for (i in seq_along(fasta))
    {
        fasta_string[i] <- as.character(fasta[i])
    }
    motif_locator <- data.frame(length(fasta_string))
    for (i in seq_along(fasta_string))
    {
        motif_locator[i] <- str_locate(fasta_string[i], pattern)
    }
    names(motif_locator) <- fasta_names
    return(motif_locator)
    
}
