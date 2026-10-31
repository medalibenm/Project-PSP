# rebuild_preprocesscore.R - Recompile preprocessCore sans threads.
# La version conda s'arrete sur "pthread_create() is 22" pendant RMA.
# Le calcul reste identique : il s'execute simplement sur un seul thread.

# L'option est nommee : BiocManager installe aussi BiocVersion, et une option
# sans nom n'atteint pas preprocessCore.
BiocManager::install("preprocessCore",
                     configure.args = list(preprocessCore = "--disable-threading"),
                     type = "source", force = TRUE, update = FALSE, ask = FALSE)

# BiocManager ne signale qu'un avertissement si la compilation echoue :
# on appelle donc la fonction qui plantait pour prouver que la correction marche.
library(preprocessCore)
invisible(rma.background.correct(matrix(c(1, 2, 3, 4, 5, 6) * 100, ncol = 2)))
message("preprocessCore ", packageVersion("preprocessCore"), " recompile sans threads")
