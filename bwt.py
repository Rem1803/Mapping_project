###attributs 
# table des suffixes
#bwt = "ATG$CG"
#1ère colonne de la table 
#tableau des positions des caractères de la bwt

###fonctions 
#first last
#recherche de motif
#utiliser les fonctions 

###gestion des fichiers
#load 
#save

###classes
#bwt = BWT(sequence qui peut etre en .fasta ou .save peu importe)
#pour qu'on puisse ensuite faire bwt.search (motif) et l'appliquer à différents types de fichiers

#j'écris en dessous une fausse séquence qu'on pourra utiliser pour tester le code

#on récupère un suffix array sorted qui peut potentiellement être un objet
array = [
    ["$ACATACAGATG"], #1
    ["ACAGATG$ACAT"], #2
    ["ACATACAGATG$"], #3
    ["AGATG$ACATAC"], #4
    ["ATACAGATG$AC"], #5
    ["ATG$ACATACAG"], #6
    ["CAGATG$ACATA"], #7
    ["CATACAGATG$A"], #8
    ["G$ACATACAGAT"], #9
    ["GATG$ACATACA"], #10
    ["TACAGATG$ACA"], #11
    ["TG$ACATACAGA"] #12
]

class Bwt:

    def __init__(self, suffix_table):
        self.suffix_table = suffix_table #tableau contenant les suffixes ordonnées : c'est la sortie de DC3 et l'entrée de BWT
        self.colonne_bwt = self.extraire_bwt() #on appelle la méthode extraire ce qui évite d'avoir à trouver la dernière colonne en amont pour la donner en paramètres
        #colonne_bwt = liste de la dernière colonne de la suffix array donnant la BWT
        self.liste_indices_bwt = self.first_last(self.colonne_bwt) 
        #liste_indices_bwt = liste de position/indice de chaque caractère de la colonne_BWT dans la 1ère colonne de la suffix array
        #là pareil pour que ce soit plus efficace on regroupe dans un dictionnaire les listes des positions de chaque lettre grâce à une méthode
        self.tableau_pos = self.positions_par_lettre()
        #liste obtenue en parcourant la colonne_BWT et en cherchant la position des A de façon successive ex : colonne_bwt = [ACAGT] A1 en position 1 puis A2 en position 3 ça donne 133 en partant du haut vers le bas

    def extraire_bwt(self): #permet d'avoir accès directement à la dernière colonne qui correspond à la BWT
        return [ligne[-1] for ligne in self.suffix_table]

    def positions_par_lettre(self):
        """Ex : {'A': [1, 3, ...], 'C': [...], 'G': [...], 'T': [...], '$': [...]}
        Position (dans la BWT) de chaque occurrence de chaque lettre, de haut en bas."""
        positions = {}
        for i, lettre in enumerate(self.colonne_bwt):
            positions.setdefault(lettre, []).append(i)
        return positions


#si on a pas les rotations
class Bwt:

    def __init__(self, texte, suffix_array):
        self.texte = texte                  # avec le $ à la fin
        self.suffix_array = suffix_array    # sortie de DC3 : liste d'entiers
        self.n = len(texte)
        self.premiere_colonne = [texte[p] for p in suffix_array]
        self.colonne_bwt = [texte[p - 1] for p in suffix_array]
        self.liste_indices_bwt = self.first_last(self.colonne_bwt)

    def position_dans_texte(self, ligne):
        return self.suffix_array[ligne]

class Bwt :
 
    def __init__(self, suffix_table) :
        self.suffix_table = suffix_table 
        self.colonne_bwt = colonne_bwt
        self.liste_indices_bwt = liste_indices_bwt 
        self.tableau_pos_A = tableau_pos_A 
        self.tableau_pos_T = tableau_pos_T #idem
        self.tableau_pos_C = tableau_pos_C #idem
        self.tableau_pos_G = tableau_pos_G #idem
        self.tableau_pos_dollar = tableau_pos_dollar #idem
        self.colonne_bwt = self.extraire_bwt()
        self.liste_indices_bwt = self.first_last(self.colonne_bwt)
        self.tableau_pos = self.positions_par_lettre()

    ###Méthodes magiques 
    #def __contains__(self, sequence_cherchee): #utile pour les comparaisons ? 



    def first_last(self, colonne_bwt):
        ...  # fonction de ta collègue
   

    ###Méthodes 

    #def first_last(self, colonne_BWT):
    #"""Fonction qui attribue à chaque caractère de la colonne_BWT son indice dans la 1ere colonne de la suffix array
    #Entrée : list, la colonne_BWT
    #Sortie : list, la liste des positions de chaque caractère de l'alphabet A,T,C,G,$
    #"""
     #compter le nombre de chaque caractères


    def recherche_motif(self, sequence_cherchee):

        """Fonction permettant de trouver un motif ADN dans une séquence ADN, selon la méthode de BWT c'est-à-dire passage par la matrice
        Entrée : str, une sous-chaine de caractères à trouver dans la grande séquence 
        Sortie : int, la position du motif dans la grande chaîne
        """
        premiere_lettre = sequence_cherchee[-1]
        print(premiere_lettre) #par exemple là on cherche le A
        table_des_suffixes = self.suffix_table

        pos_lettre_cherchee = self.tableau_pos_[lettre_cherchee] # théoriquement on parcourt la colonne bwt à la recherche du premier A
        #mais pour optimiser on va direct dans tableau_pos_A qui donne la position des A dans la BWT là dans l'exemple c'est position 7
        #donc ici pos_lettre_cherchee est 7 car on commence par le A de haut en bas
        lettre_cherchee = self.colonne_bwt[pos_lettre_cherchee] #là on se positionne sur le A dans la BWT
        
        ###Etape 4
        indice_lettre_cherchee = self.liste_indices_bwt[lettre_cherchee] #là on récup l'indice 2 car c'est à cet indice que se trouve le 1er A dans la première colonne
        #table_des_suffixes[1][indice_lettre_cherchee] correspond à l'endroit ou se trouve le 1er A dans la 1ere colonne 

        lettre_a_gauche = table_des_suffixes[len(table_des_suffixes)][indice_lettre_cherchee] #on récupère la lettre se trouvant sur la meme ligne dans la derniere colonne de la matrice cad la colonne bwt
        #elle correspond à la lettre venant juste avant le A

        #si c'est G on recommence 
        if lettre_a_gauche == sequence_cherchee[-2] :
            #recommencer à partir de l'etape 4

        #si c'est pas G : il faut regarder A2 c'est à dire le 2e A apparaissant dans la colonne BWT

        #si on trouve toujours pas alors c'est que la séquence n'est pas dans la matrice des suffixes
            print(f"La séquence {sequence_cherchee} n'est pas présente dans le génome")


###Main
objet_suffix_table = Bwt(array) #création de l'objet suffix array
sequence_cherchee = 'CGA' #la séquence test que je cherche à retrouver dans la BWT
objet_suffix_table.recherche_motif(sequence_cherchee) #j'applique la méthode recherche motif à la table des suffixes car c'est dedans que l'on trouvera le motif