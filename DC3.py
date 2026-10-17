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
                

 

texte = "ACA$"
attendu = sorted(range(len(texte)), key=lambda i: texte[i:])

if __name__ == "__main__": 
    resultat = DC3(texte).suffix_array
    print("Obtenu :", resultat)
    print("Attendu :", attendu)
    assert resultat == attendu
