class Bwt:

    def __init__(self, texte, suffix_array):
        self.texte = texte                  # avec le $ à la fin
        self.suffix_array = suffix_array    # sortie de DC3 : liste d'entiers
        self.n = len(texte)
        self.premiere_colonne = [texte[p] for p in suffix_array]
        self.colonne_bwt = [texte[p - 1] for p in suffix_array]  # Dernière colonne (la BWT) : la lettre qui précède chaque suffixe
        self.liste_indices_bwt = self.first_last(self.colonne_bwt)

        ###Remplace first_last par un dictionnaire pour stocker les positions de chaque lettre dans la première colonne
        self.debut ={} #dictionnaire 
        nb = 0
        for lettre in sorted(set(texte)): #on cherche le bloc de la lettre dans la 1ere colonne 
            self.debut[lettre] = nb
            nb += texte.count(lettre) #permet d'obtenir le nombre de lettres avant la lettre c dans la 1ere colonne

        alphabet = sorted(set(texte)) #on récupère l'alphabet de la séquence
        self.apparition = {lettre: [0] * (self.n + 1) for lettre in set(texte)} #dictionnaire de listes pour compter les occurrences de chaque lettre
        for i, lettre in enumerate(self.colonne_bwt):
            for a in alphabet :
                self.apparition[a][i + 1] = self.apparition[a][i] + (a==lettre) #on incrémente le compteur de la lettre si elle est présente à l'indice i de la bwt

    ###Méthodes magiques###
    def __contains__(self, sequence_cherchee): #utile pour les comparaisons ?  
        return len(self.recherche_motif_intervalles(sequence_cherchee)) > 0
    #à voir pour en ajouter notamment avec la gestion de fichiers
    ### ###

    def first_last(self, colonne_BWT):
        """Fonction qui attribue à chaque caractère de la colonne_BWT son indice dans la 1ere colonne de la suffix array
        Entrée : list, la colonne_BWT
        Sortie : list, la liste des positions de chaque caractère de l'alphabet A,T,C,G,$
        -> renvoie le numéro de la ligne dans la premiere colonne qui correspond à la lettre i de la bwt
        """
        # Compter le nombre de chaque caractère
        L_pos = [] #liste des positions de chaque lettre dans la première colonne
        count_A = 0
        count_T = 0
        count_C = 0
        count_G = 0
        count_doll= 0
        for elem in self.colonne_bwt:
            if elem == 'A' :
                count_A += 1
                L_pos.append(count_A)
            elif elem == 'T' :
                count_T += 1
                L_pos.append(count_T)
            elif elem == 'C' :
                count_C += 1
                L_pos.append(count_C)
            elif elem == 'G' :
                count_G += 1
                L_pos.append(count_G)
            elif elem == '$' :
                count_doll += 1
                L_pos.append(count_doll)
        return self, L_pos
            



    def recherche_motif_intervalles(self, sequence_cherchee):
        """Fonction permettant de trouver un motif ADN dans une séquence ADN, selon la méthode de BWT c'est-à-dire passage par la matrice
        Entrée : str, une sous-chaine de caractères à trouver dans la grande séquence 
        Sortie : list, les positions des motifs dans la grande chaîne = le génome"""
        
        # Initialisation des intervalles pour la recherche du motif
        debut = 0
        fin = self.n #demi ouvert [debut, fin[

        # Parcours du motif de la dernière lettre à la première dans la première colonne
        for lettre in reversed(sequence_cherchee): #on commence bien par la gauche

            if lettre not in self.debut : #si la lettre cherchée n'est pas au début on ajoute rien
                return []
            
            ###Définition des bornes reste à définir self.debut et self.apparition
            debut = self.debut[lettre]+ self.apparition[lettre][debut] #lettre prend la valeur str de chaque lettre du motif cherché donc là en complexité on est de l'ordre de la longueur du motif cherché pas du génome entier
            fin = self.debut[lettre]+self.apparition[lettre][fin]

            #à partir de là on a des bornes pour chaque lettre de la sous séquence cherchée

            # Si l'intervalle est vide, le motif n'est pas présent
            if debut >= fin: 
                return []

        # Récupération des positions dans le texte original
        positions = sorted([self.suffix_array[i] for i in range(debut, fin)])
        return positions


###Main
if __name__ == "__main__":
    texte = "ACATACAGATG$" #exemple de texte pour tester les fonctions
    suffix_array = sorted(range(len(texte)), key=lambda p: texte[p:])  # en attendant DC3
    objet_bwt = Bwt(texte, suffix_array) #création d'un objet bwt

    for motif in ["CGA", "ACA", "GAT", "TG", "ACATACAGATG"]:
        positions = objet_bwt.recherche_motif_intervalles(motif)
        if positions: #si le motif est présent
            print(f"{motif} a été trouvée aux positions {positions}.")
        else: #si la séquence est vide
            print(f"La séquence cherchée {motif} n'est pas présente dans le génome.")

