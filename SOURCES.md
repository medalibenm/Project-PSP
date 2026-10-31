# Sources des données

Toutes les données sont publiques, anonymisées et librement téléchargeables. `run.sh` les
télécharge et vérifie leur intégrité. Elles ne sont pas jointes au dépôt.

| Identifiant | Rôle | Adresse | Taille | Téléchargement |
|---|---|---|---|---|
| GSE34248 (`GSE34248_RAW.tar`) | jeu de données principal : 28 fichiers CEL Affymetrix, 14 patients | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE34nnn/GSE34248/suppl/GSE34248_RAW.tar | 158,6 Mo | 2026-10-31 |
| GSE34248 (series matrix) | métadonnées des échantillons | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE34nnn/GSE34248/matrix/GSE34248_series_matrix.txt.gz | 4,5 Mo | 2026-10-31 |
| GPL570 | annotation sonde vers gène, puce HG-U133 Plus 2.0 | https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz | 8,5 Mo | 2026-10-31 |

Page de la série : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE34248

Publication associée : PubMed 23308107. La série est le sous-ensemble « Set 1 » de GSE41664.

Les empreintes SHA-256 attendues sont dans `data/checksums.sha256`.
