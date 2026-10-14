class Triplets:

    def __init__(self, T, P):
        self.R = [(tuple(T[i:i + 3]), i) for i in P]
        self.radix_sort()

    def radix_sort(self):
        R = list(self.R)

        if not self.R :
            return []

        def composants(cle):
            return cle if isinstance(cle, tuple) else (cle,)

        key_length = len(composants(R[0][0]))

        # On trie d'abord le 3e chiffre, puis le 2e, puis le 1er.
        for chiffre in range(key_length-1, -1, -1):
                
                max_value = max(composants(cle)[chiffre] for cle, position in R) #Prend le plus grand des éléments en positions 2 puis 1 puis 0 entre les triplets 
                counts = [0] * (max_value + 1) #Chaque case correspond à une valeur possible du chiffre qu’on est en train de trier, parmi les triplets.
        
                #compte les triplets dont le n-ième élément est égal à max_value 
                for cle, position in R:
                    counts[composants(cle)[chiffre]] += 1
        
                # Transforme les comptes en positions finales dans le tableau trié.
                for i in range(1, len(counts)):
                    counts[i] += counts[i - 1]
        
                sorted_R = [None] * len(R)
        
                # Place les éléments dans l'ordre, en préservant le tri précédent.
                for record in reversed(R):
                    cle, position = record
                    value = composants(cle)[chiffre]
                    counts[value] -= 1
                    sorted_R[counts[value]] = record
        
                # Indispensable : le prochain passage trie le résultat du précédent.
                R = sorted_R
    
        self.Rsorted = [(cle, position) for cle, position in R]
        self.index = [position for cle, position in R]
        rang_par_triplet = {
            triplet: rang
            for rang, triplet in enumerate(
                dict.fromkeys(cle for cle, position in self.Rsorted),
                start=1
            )
        }
        self.order = [
            rang_par_triplet[triplet]
            for triplet, position in self.Rsorted
        ]
        self.T_prime = [rang_par_triplet[triplet] for triplet, position in self.R]

        return self.Rsorted



class DC3:

    def __init__(self, T):
        T = [ord(char) for char in T] if isinstance(T, str) else list(T)
        n = len(T)
        T += [0,0,0]
        self.P_0 = [i for i in range(n) if i % 3 == 0]
        self.P_1 = [i for i in range(n) if i % 3 == 1]
        self.P_2 = [i for i in range(n) if i % 3 == 2]
        self.P_12 = self.P_1 + self.P_2 #concaténation P_12s

        self.triplets = Triplets(T, self.P_12)
        T_prime = self.triplets.T_prime
        R = self.triplets.R

        if len(set(T_prime)) == len(T_prime):
            T_prime_suffix_array = sorted(
        range(len(T_prime)),
        key=lambda i: T_prime[i]
    )
        else:
            T_prime_suffix_array = DC3(T_prime).suffix_array
             
        
            
        R_12 = [self.P_12[index] for index in T_prime_suffix_array]
        rang_12 = {
                position: rang
                for rang, position in enumerate(R_12, start=1)
            }

        cles_R0 = [
            ((T[i], rang_12.get(i + 1, 0)), i)
                for i in self.P_0
            ]

            # Tri lexicographique des clés (caractère, rang du suffixe suivant).
        cles_R0_sorted = sorted(cles_R0, key=lambda element: element[0])
        ordre_R0 = [position for cle, position in cles_R0_sorted]

        self.suffix_array = self.merge(T, ordre_R0, R_12)


    
    @staticmethod
    def merge(T, R_0_prime_index, R_12_prime):
        rank_by_position = {
            position: rank
            for rank, position in enumerate(R_12_prime, start=1)
        }

        i = j = 0
        R_merged = []
        while i < len(R_0_prime_index) and j < len(R_12_prime):
            a = R_0_prime_index[i]
            b = R_12_prime[j]

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
