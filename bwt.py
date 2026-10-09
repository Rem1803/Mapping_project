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
array = [[$ACATACAGATG],
[ACAGATG$ACAT],
[ACATACAGATG$],
[AGATG$ACATAC],
[ATACAGATG$AC],
[ATG$ACATACAG],
[CAGATG$ACATA],
[CATACAGATG$A],
[G$ACATACAGAT],
[GATG$ACATACA],
[TACAGATG$ACA],
[TG$ACATACAGA]]

 class Bwt :
 
    def __init__(self, suffix_table, colonne_bwt, liste_indices_bwt, tableau_pos_A, tableau_pos_T, tableau_pos_C, tableau_pos_G, tableau_pos_dollar)
        self.suffix_table = suffix_table #tableau contenant les suffixes ordonnées : c'est la sortie de DC3 et l'entrée de BWT
        self.colonnne_bwt = colonne_bwt #liste dernière colonne de la suffix array donnant la BWT
        self.liste_indices_bwt = liste_indices_bwt #liste de position/indice de chaque caractère de la colonne_BWT dans la 1ère colonne de la suffix array
        self.tableau_pos_A = tableau_pos_A #liste obtenue en parcourant la colonne_BWT et en cherchant la position des A de façon successive ex : colonne_bwt = [ACAGT] A1 en position 1 puis A2 en position 3 ça donne 133 en partant du haut vers le bas
        self.tableau_pos_T = tableau_pos_T #idem
        self.tableau_pos_C = tableau_pos_C #idem
        self.tableau_pos_G = tableau_pos_G #idem
        self.tableau_pos_dollar = tableau_pos_dollar #idem

    ###Méthodes magiques 
    def __contains__(self, sequence_cherchee):
    

    ###Méthodes 

    def first_last(self, colonne_BWT):
    """Fonction qui attribue à chaque caractère de la colonne_BWT son indice dans la 1ere colonne de la suffix array
    Entrée : list, la colonne_BWT
    Sortie : list, la liste des positions de chaque caractère de l'alphabet A,T,C,G,$
    """

    def recherche_motif(self, sequence_cherchee):
    """Fonction permettant de trouver un motif ADN dans une séquence ADN
    Entrée : str, une sous-chaine de caractères à trouver dans la grande séquence 
    Sortie : int, la position du motif dans la grande chaîne
    """

 

