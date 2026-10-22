import random
from array import array
from itertools import accumulate
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Génome de référence au format FASTA
record_reference = list(SeqIO.parse(r"C:\Users\lebou\OneDrive\Documents\COURS\UE2 - GENOMICS\Projet Mapping\GCF_000002765.5_GCA_000002765_genomic.fna", "fasta"))

# Reads au format FASTQ
record_reads = list(SeqIO.parse(r"C:\Users\lebou\OneDrive\Documents\COURS\UE2 - GENOMICS\Projet Mapping\single_Pfal_dat.fq", "fastq"))


class Triplets:


    def __init__(self, T=None, P=None, R=None):
        """
        Prépare les clés à trier, puis calcule leur ordre radix.


        Args:
            T : Texte encodé en entiers et complété par des sentinelles.
            P : Index de départ des triplets à extraire de T.
            R : Liste d'enregistrements (clé, position) à trier directement.
                Les clés peuvent être des entiers ou des tuples.
        """
        if R is not None: #si la liste des triplets est donnée
            self.R = R  
        else: #si la liste des triplets n'est pas donnée on la construit
            self.R = [(tuple(T[i:i + 3]), i) for i in P]


        self.radix_sort()


    def radix_sort(self):
        """Trie les clés par radix et calcule leurs rangs.


        Les enregistrements (R) sont des paires (clé, position). Une clé peut
        être un entier ou un tuple d'entiers. Le tri stable de chaque
        composante produit un ordre lexicographique pour les tuples.


        Returns:
            Rsorted : Les clées triées
        """
        R = list(self.R)


        if not self.R:
            # Initialise aussi les attributs utilisés par les étapes suivantes.
            self.Rsorted = []
            self.index = []
            self.order = []
            self.T_prime = []
            return []




        def composants(cle):
            """Retourne les composantes d'une clé sous forme de tuple.


            Args:
                cle : Clé à normaliser. Peut être un entier ou un tuple d'entiers.


            Returns:
                Un tuple contenant les composantes de la clé. Si 'cle' est déjà un tuple, il est retourné tel quel ;
                sinon, il est enveloppé dans un tuple à un seul élément.
            """
            return cle if isinstance(cle, tuple) else (cle,)


        key_length = len(composants(R[0][0]))


        # Le tri stable de droite à gauche construit l'ordre lexicographique.
        for chiffre in range(key_length-1, -1, -1):
                values = [
                composants(cle)[chiffre]
                for cle, _ in R]
                max_value = max(values)
                counts = [0] * (max_value + 1)
       
                # Compte le nombre d'occurences de chaque valeur
                for value in values:
                    counts[value] += 1
       
                # Convertit les comptes en positions cumulées dans le résultat.
                for i in range(1, len(counts)):
                    counts[i] += counts[i - 1]
       
                sorted_R = [None] * len(R)
       
                # Le parcours inversé préserve la stabilité de cette passe.
                for i in range(len(R)-1, -1, -1):
                    value = values[i]
                    counts[value] -= 1
                    sorted_R[counts[value]] = R[i]
       
                # La passe suivante utilise cet ordre comme critère secondaire.
                R = sorted_R
   
        self.Rsorted = [(cle, position) for cle, position in R]
        self.index = [position for cle, position in R]


        # Attribue un rang à chaque clé distincte selon l'ordre trié.
        rang_par_triplet = {
            triplet: rang
            for rang, triplet in enumerate(
                dict.fromkeys(cle for cle, position in self.Rsorted),
                start=1
            )
        }
        # Attribue un ordre à chaque clé distincte selon l'ordre trié
        self.order = [
            rang_par_triplet[triplet]
            for triplet, position in self.Rsorted
        ]
        # Création d'une liste des rangs dans l'ordre initial
        self.T_prime = [rang_par_triplet[triplet] for triplet, position in self.R]


        return self.Rsorted






class DC3:


    def __init__(self, T):
        """
        Args:
            T : chaîne de caractères ou liste d'entiers


        Construit le tableau des suffixes de T avec l'algorithme DC3.


        DC3 trie d'abord les suffixes dont la position n'est pas multiple
        de trois, puis ceux dont la position est multiple de trois, avant
        de fusionner les deux séquences triées.
        """


        # Conversion de la chaîne de caractères en liste d'entiers
        T = [ord(char) for char in T] if isinstance(T, str) else list(T)
        n = len(T)
        # Trois sentinelles permettent de lire les triplets en fin de texte.
        T += [0,0,0]


        # Répartit les positions selon leur reste modulo trois.
        self.P_0 = [i for i in range(n) if i % 3 == 0]
        self.P_1 = [i for i in range(n) if i % 3 == 1]
        self.P_2 = [i for i in range(n) if i % 3 == 2]
        self.P_12 = self.P_1 + self.P_2


        # Cas particulier DC3 : ajoute une position sentinelle à l'échantillon.
        if n % 3 == 1:
            self.P_12.append(n)


        # Récupère les rangs compressés calculés par Triplets
        self.triplets = Triplets(T, self.P_12)
        T_prime = self.triplets.T_prime


        # Si tous les rangs sont distincts, leur ordre donne directement
        # l'ordre des suffixes échantillonnés; sinon, applique DC3 récursivement.
        if len(set(T_prime)) == len(T_prime):
            T_prime_suffix_array = sorted(
                range(len(T_prime)),
                key=lambda i: T_prime[i]
            )
        else:
            T_prime_suffix_array = DC3(T_prime).suffix_array


        # Convertit l'ordre récursif en positions du texte d'origine et
        # construit les rangs nécessaires au tri et à la fusion.
        R_12 = [self.P_12[index] for index in T_prime_suffix_array]
        rang_12 = {
                position: rang
                for rang, position in enumerate(R_12, start=1)
            }


        cles_R0 = [
            ((T[i], rang_12.get(i + 1, 0)), i)
                for i in self.P_0
            ]


        # Trie les positions multiples de trois à partir de leurs clés
        # (caractère courant, rang du suffixe suivant).
        tri_R0 = Triplets(R=cles_R0)
        ordre_R0 = tri_R0.index


        # Fusionne les deux ordres partiels, puis retire la position sentinelle éventuellement ajoutée dans le cas n % 3 == 1.
        self.suffix_array = self.merge(T, ordre_R0, R_12, rang_12)
        self.suffix_array = [
                position for position in self.suffix_array
            if position < n
        ]






   
    @staticmethod
    def merge(T, R_0_prime_index, R_12_prime, rank_by_position):
        """Fusionne les suffixes triés de classes 0 et 1/2.


        Args:
            T: Texte encodé en entiers et complété par des sentinelles.
            R_0_prime_index: Positions de classe 0 dans l'ordre trié.
            R_12_prime: Positions de classe 1 et 2 dans l'ordre trié.
            rank_by_position: Rang de chaque suffixe échantillonné.


        Returns:
            Les positions de tous les suffixes dans l'ordre lexicographique.
        """
        i = j = 0
        R_merged = []
        while i < len(R_0_prime_index) and j < len(R_12_prime):
            a = R_0_prime_index[i]
            b = R_12_prime[j]


            # La clé de comparaison dépend de la classe modulo trois de b.
            if b % 3 == 1:
                key_0 = (T[a], rank_by_position.get(a + 1, 0))
                key_12 = (T[b], rank_by_position.get(b + 1, 0))
            else:
                key_0 = (
                    T[a],
                    T[a + 1],
                    rank_by_position.get(a + 2, 0)
                )
                key_12 = (
                    T[b],
                    T[b + 1],
                    rank_by_position.get(b + 2, 0)
                )


            if key_0 <= key_12:
                R_merged.append(a)
                i += 1
            else:
                R_merged.append(b)
                j += 1


        # Ajoute les éléments restants de la séquence non épuisée.
        R_merged.extend(R_0_prime_index[i:])
        R_merged.extend(R_12_prime[j:])
        return R_merged






def decouper_en_kmers(sequence, k):
    """Découpe la séquence en k-mers de taille k.
    Entrée : int, la taille des k-mers
    Sortie : list, les k-mers"""
    return [(i, sequence[i:i + k]) for i in range(0, len(sequence) - k + 1, k)]




class Bwt:


    def __init__(self, texte, suffix_array):
        self.texte = texte                  # avec le $ à la fin
        self.suffix_array = suffix_array    # sortie de DC3 : liste d'entiers
        self.n = len(texte)
        self.premiere_colonne = [texte[p] for p in suffix_array]
        self.colonne_bwt = [texte[p - 1] for p in suffix_array]  # Dernière colonne (la BWT) : la lettre qui précède chaque suffixe


        ###Remplace first_last par un dictionnaire pour stocker les positions de chaque lettre dans la première colonne
        self.debut ={} #dictionnaire
        nb = 0
        for lettre in sorted(set(texte)): #on cherche le bloc de la lettre dans la 1ere colonne
            self.debut[lettre] = nb
            nb += texte.count(lettre) #permet d'obtenir le nombre de lettres avant la lettre c dans la 1ere colonne


        alphabet = sorted(set(texte)) #on récupère l'alphabet de la séquence
        self.apparition = {
            a: array('I', accumulate((x == a for x in self.colonne_bwt), initial=0))
            for a in alphabet
        }


    ###Méthodes magiques###
    def __contains__(self, sequence_cherchee): #utile pour les comparaisons ?  
        return len(self.recherche_motif_intervalles(sequence_cherchee)) > 0
    #à voir pour en ajouter notamment avec la gestion de fichiers
    ### ###


    def recherche_motif_intervalles(self, k_mer, max_occ = None):
        """Fonction permettant de trouver un motif ADN dans une séquence ADN, selon la méthode de BWT c'est-à-dire passage par la matrice
        Entrée : str, une sous-chaine de caractères à trouver dans la grande séquence
        Sortie : list, les positions des motifs dans la grande chaîne = le génome"""
        #Recherche rapide exacte
       
        # Initialisation des intervalles pour la recherche du motif
        debut = 0
        fin = self.n #demi ouvert [debut, fin[


        # Parcours du motif de la dernière lettre à la première dans la première colonne
        for lettre in reversed(k_mer): #on commence bien par la gauche


            if lettre not in self.debut : #si la lettre cherchée n'est pas au début on ajoute rien
                return []
           
            ###Définition des bornes reste à définir self.debut et self.apparition
            debut = self.debut[lettre]+ self.apparition[lettre][debut] #lettre prend la valeur str de chaque lettre du motif cherché donc là en complexité on est de l'ordre de la longueur du motif cherché pas du génome entier
            fin = self.debut[lettre]+self.apparition[lettre][fin]


            #à partir de là on a des bornes pour chaque lettre de la sous séquence cherchée


            # Si l'intervalle est vide, le motif n'est pas présent
            if debut >= fin:
                return []


        if max_occ is not None and fin - debut > max_occ:
            return []
        return sorted(self.suffix_array[i] for i in range(debut, fin))


        # Récupération des positions dans le texte original
        #positions = sorted([self.suffix_array[i] for i in range(debut, fin)])
        #return positions


    def recherche_motif_intervalles_avec_mutations(self, sequence_cherchee, max_erreurs=1):
        """Trouve les positions du motif avec au plus max_erreurs substitutions. Si on met max_erreurs = 0, ça revient à la même chose que la fonction d'avant.
        Retourne un dico {position: nombre_de_mutations}.
        Donc, max_erreurs : nombre maximum de substitutions autorisées par read"""
        #Fonction lente mais qui permet de gérer les subsititions -> à utiliser si la version exacte ne suffit pas.
        resultats = {}
        alphabet = [c for c in self.debut if c != "$"]


        def explorer(i, debut, fin, erreurs):
            # i : index de la lettre à traiter dans le motif (on recule vers 0)
            if i < 0:
                for r in range(debut, fin):
                    pos = self.suffix_array[r]
                    # on garde le meilleur score si une position est atteinte par plusieurs branches
                    resultats[pos] = min(erreurs, resultats.get(pos, max_erreurs + 1))
                return


            attendue = sequence_cherchee[i]
            for lettre in alphabet: #exploration sur l'alphabet donc plus long
                cout = erreurs + (lettre != attendue)
                if cout > max_erreurs:
                    continue  # budget d'erreurs dépassé : on élague
                nouveau_debut = self.debut[lettre] + self.apparition[lettre][debut]
                nouvelle_fin = self.debut[lettre] + self.apparition[lettre][fin]
                if nouveau_debut < nouvelle_fin:  # intervalle non vide
                    explorer(i - 1, nouveau_debut, nouvelle_fin, cout)


        explorer(len(sequence_cherchee) - 1, 0, self.n, 0)
        return dict(sorted(resultats.items())) #renvoie un dictionnaire trié par position pour plus de lisibilité


    def aligner_read(self, read, max_erreurs=2, erreurs_kmer=0, k=20, seuil_occ=200):
        """Fonction qui permet d'aligner un read sur le génome avec un nombre d'erreurs donné.
        Entrée : str, un read et un nombre maximal d'erreurs autorisées
        Sortie : list, tuples (positions,erreurs)
        """
        k = min(k, len(read) // (max_erreurs + 1))
        if k < 1:
            return []
       
        candidats = set() #stoke les positions de départ possible sur le génome sans doublons (fonction set())
       
        for i, kmer in decouper_en_kmers(read, k): #on découpe le read en k-mers de taille k
                 
            if erreurs_kmer == 0:
                positions = self.recherche_motif_intervalles(kmer, max_occ=seuil_occ) #on cherche les positions de départ possibles pour chaque k-mer du read
            else :
                positions = self.recherche_motif_intervalles_avec_mutations(kmer, erreurs_kmer) #on cherche les positions de départ possibles pour chaque k-mer du read


            if len(positions) > seuil_occ:
                continue
           
            for pos in positions:
                candidats.add(pos - i) #on ajoute la position de départ possible sur le génome en soustrayant l'index du k-mer dans le read


        resultats = [] #stocke les alignements valides
       
        for debut in candidats: #on parcout les positions de départ possibles sur le génome
            if debut >= 0 and debut+len(read) <= self.n-1: #on vérifie que le read ne dépasse la taille du génome


                segment = self.texte[debut:debut+len(read)] #on récupère le segment du génome correspondant à la position de départ
                erreurs = sum(a != b for a, b in zip(read, segment)) #on compte le nombre d'erreurs entre le segment du génome et le read


                if erreurs <= max_erreurs: #si le nombre d'erreurs est inférieur ou égal au nombre d'erreurs autorisées
                    resultats.append((debut, erreurs)) #on ajoute la position de départ et le nombre d'erreurs à la liste des résultats


       
        return sorted(resultats, key=lambda c : c[1]) #on trie les résultats par nombre d'erreurs croissant


 


def normaliser_sequence(sequence, ALPHABET_ADN = set("ACGT")):
    """
    Normalise une chaine de caractères : met tout en majuscules et vérifie que tous les caractères sont valides
    """
    sequence = str(sequence).upper()

    caracteres_invalides = set(sequence) - ALPHABET_ADN

    if caracteres_invalides:
        raise ValueError(f"Caractères invalides : {caracteres_invalides}")
    return sequence




#####################################################################
# MAIN PROGRAM
#####################################################################

# Construction du génome complet avec son complément inverse
"""
sequence = ""
for chromosome in record_reference:
    sequence += str(chromosome.seq)
sequence = normaliser_sequence(sequence)
sequence += str(Seq(sequence).reverse_complement())
sequence += "$"
"""
sequence = str(record_reference[0].seq)
sequence = normaliser_sequence(sequence)
sequence += "$"


resultat = DC3(sequence).suffix_array
    #print("Obtenu :", resultat)
    #print("Attendu :", attendu)
    #assert resultat == attendu

suffix_array = resultat
objet_bwt = Bwt(sequence, suffix_array)

m = sequence[0:50]
print(objet_bwt.recherche_motif_intervalles(m))


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

reads = list(str(read.seq) for read in record_reads)

 

texte = sequence


"""if __name__ == "__main__":
    resultat = DC3(texte).suffix_array
    print("Obtenu :", resultat)
    print("Attendu :", attendu)
    assert resultat == attendu


    suffix_array = resultat
    objet_bwt = Bwt(texte, suffix_array) #création d'un objet bwt


    for motif in ["CGA", "ACA", "GAT", "TG", "ACATACAGATG"]:
        positions = objet_bwt.recherche_motif_intervalles_avec_mutations(motif)
        if positions: #si le motif est présent
            print(f"{motif} a été trouvée aux positions {positions}.")
        else: #si la séquence est vide
            print(f"La séquence cherchée {motif} n'est pas présente dans le génome.")
"""


if __name__ == "__main__":
    # 1. Vérification du suffix array
    resultat = DC3(texte).suffix_array
    suffix_array = resultat




    # 2. Construction de la BWT
    objet_bwt = Bwt(texte, suffix_array)


    # 3. Alignement des reads
    k = 4
    max_erreurs = 1
    reads = ["ACATACAG", "ACATTCAG", "GATG", "CCCCCCCC"]


    for read in reads:
        alignements = objet_bwt.aligner_read(read, k=k, max_erreurs=max_erreurs)


        if not alignements:
            print(f"Read {read} : NON TROUVÉ dans le génome "
                  f"(avec au plus {max_erreurs} erreur(s)).")
            continue


        print(f"Read {read} : TROUVÉ ({len(alignements)} alignement(s)).")
        for position, erreurs in alignements:
            if erreurs == 0:
                print(f"   - à la position {position}, correspondance exacte.")
            else:
                print(f"   - à la position {position}, avec {erreurs} substitution(s).")
    # 4. Test de cohérence de l'index : un motif pris dans le génome lui-même
    #    doit toujours être retrouvé à sa position d'origine.




    longueur_motif = 4 #20 pour vrai génome
    nb_tests = 5 #20 pour vrai génome
    echecs = 0


    for _ in range(nb_tests):
        # position aléatoire (on exclut le $ final)
        p = random.randint(0, len(texte) - 1 - longueur_motif)
        motif = texte[p:p + longueur_motif]
        positions = objet_bwt.recherche_motif_intervalles(motif)


        if p in positions:
            print(f"OK : motif pris en {p} retrouvé ({len(positions)} occurrence(s)).")
        else:
            echecs += 1
            print(f"ÉCHEC : motif pris en {p} NON retrouvé. Positions renvoyées : {positions[:5]}")


    if echecs == 0:
        print("\nL'index (DC3 + BWT) est cohérent : le problème vient des reads ou du chargement des données.")
    else:
        print(f"\n{echecs} échec(s) sur {nb_tests} : le problème vient de DC3 ou de la BWT.")
