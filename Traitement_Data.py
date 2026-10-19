# Pipeline de traitement des données pour l'algo de Mapping:

# 1- importer les séquences des reads (1,5 millions) et du génome de référence (15 chromosomes)

# 2- traitement des références : on met tout en majuscules et on vérifie que tous les caractères sont valides (A, C, G, T)

# 3- traitement des reads : on enlève les morceaux de séquences dont le score de qualité est insuffisante, puis on ne garde
#    que les reads dont la longueur est suffisante (seuils à ajuster)

# 4- encodage : remplace tous les caractères par des entiers en gardant l'ordre lexicographique (utile pour DC3)

# 5- découpage des reads en k-mers (sert pour le mapping approché)


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Génome de référence au format FASTA
record_reference = list(SeqIO.parse("GCF_000002765.5_GCA_000002765_genomic.fna", "fasta"))

# Reads au format FASTQ
record_reads = list(SeqIO.parse("single_Pfal_dat.fq", "fastq"))

#print(record_reference[0].seq)
#a = list(record_reads)
#print(a[0])




def normaliser_sequence(sequence, ALPHABET_ADN = set("ACGT")):
    """
    Normalise une chaine de caractères : met tout en majuscules et vérifie que tous les caractères sont valides
    """
    sequence = str(sequence).upper()

    caracteres_invalides = set(sequence) - ALPHABET_ADN

    if caracteres_invalides:
        raise ValueError(f"Caractères invalides : {caracteres_invalides}")
    return sequence



# Paramètres du trimming
leading_threshold = 30
trailing_threshold = 30
window_size = 2
window_quality = 30
min_length = 30

# Statistiques
nb_initial = len(record_reads)
nb_supprimes = 0
nb_coupes = 0

reads_conserves = []

for record in record_reads:

    sequence = str(record.seq)
    qualites = list(record.letter_annotations["phred_quality"])

    longueur_initiale = len(sequence)

    # --------------------------------------------------
    # 1. LEADING : suppression des bases faibles au début
    # --------------------------------------------------
    debut = 0

    while debut < len(qualites) and qualites[debut] < leading_threshold:
        debut += 1

    sequence = sequence[debut:]
    qualites = qualites[debut:]

    # --------------------------------------------------
    # 2. TRAILING : suppression des bases faibles à la fin
    # --------------------------------------------------
    fin = len(qualites)

    while fin > 0 and qualites[fin - 1] < trailing_threshold:
        fin -= 1

    sequence = sequence[:fin]
    qualites = qualites[:fin]

    # --------------------------------------------------
    # 3. SLIDINGWINDOW depuis le début
    # Coupe au début de la première fenêtre de qualité
    # moyenne inférieure au seuil.
    # --------------------------------------------------
    position_coupe = len(qualites)

    if len(qualites) >= window_size:
        for i in range(len(qualites) - window_size + 1):
            fenetre = qualites[i:i + window_size]

            if sum(fenetre) / window_size < window_quality:
                position_coupe = i
                break

    sequence = sequence[:position_coupe]
    qualites = qualites[:position_coupe]

    # --------------------------------------------------
    # 4. SLIDINGWINDOW depuis la fin
    # Examine les fenêtres en partant de l'extrémité droite.
    # --------------------------------------------------
    position_coupe = 0

    if len(qualites) >= window_size:
        for fin_fenetre in range(
            len(qualites), window_size - 1, -1
        ):
            debut_fenetre = fin_fenetre - window_size
            fenetre = qualites[debut_fenetre:fin_fenetre]

            if sum(fenetre) / window_size < window_quality:
                position_coupe = fin_fenetre
                break

        # Si une fenêtre faible est trouvée, conserver la
        # séquence située à droite de cette fenêtre.
        if position_coupe > 0:
            sequence = sequence[position_coupe:]
            qualites = qualites[position_coupe:]

    # Comptabiliser les reads ayant subi une coupe
    if len(sequence) < longueur_initiale:
        nb_coupes += 1

    # --------------------------------------------------
    # 5. MINLEN : élimination des reads trop courts
    # --------------------------------------------------
    if len(sequence) < min_length:
        nb_supprimes += 1
        continue

    # Mettre à jour le SeqRecord existant
    record.letter_annotations.clear()
    record.seq = Seq(sequence)
    record.letter_annotations["phred_quality"] = qualites

    reads_conserves.append(record)

# Modifier la liste existante, sans remplacer son identité
record_reads[:] = reads_conserves

# Afficher les statistiques
print(f"Reads initiaux       : {nb_initial}")
print(f"Reads coupés         : {nb_coupes}")
print(f"Reads supprimés      : {nb_supprimes}")
print(f"Reads conservés      : {len(record_reads)}")
print(f"Taux de conservation : {len(record_reads) / nb_initial:.2%}"
      if nb_initial else "Aucun read en entrée")



print(f"après : {len(record_reads)}")


#print(record_reads[0].letter_annotations["phred_quality"])