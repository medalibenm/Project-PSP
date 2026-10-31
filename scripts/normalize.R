# normalize.R - Normalise les fichiers CEL avec RMA et ecrit une matrice log2.
# Usage : Rscript normalize.R <dossier_cel> <sortie.tsv>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    message("Usage : Rscript normalize.R <dossier_cel> <sortie.tsv>")
    quit(save = "no", status = 2)
}
cel_dir <- args[1]
output <- args[2]

suppressPackageStartupMessages(library(affy))

files <- list.files(cel_dir, pattern = "\\.cel\\.gz$", ignore.case = TRUE,
                    full.names = TRUE)
if (length(files) == 0) {
    message("Aucun fichier CEL dans ", cel_dir)
    quit(save = "no", status = 3)
}

expression <- exprs(rma(ReadAffy(filenames = files)))
colnames(expression) <- sub("^(GSM[0-9]+).*$", "\\1", colnames(expression))
expression <- round(expression, 4)

table <- data.frame(probe = rownames(expression), expression, check.names = FALSE)
write.table(table, output, sep = "\t", quote = FALSE, row.names = FALSE)
message(nrow(expression), " sondes, ", ncol(expression), " echantillons")
