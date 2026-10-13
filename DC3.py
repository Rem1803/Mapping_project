class Triplets:

    def __init__(self, T, P):
        self.R = [(tuple(T[i:i + 3]), i) for i in P]
        self.radix_sort()

    def radix_sort(self):
        R = self.R
        key_length = len(R[0][0])
    # On trie d'abord le 3e chiffre, puis le 2e, puis le 1er.
        for chiffre in range(key_length):
            if not R:
                break
    
            max_value = max(triplet[chiffre] for triplet, position in R) #Prend le plus grand des éléments en positions 2 puis 1 puis 0 entre les triplets 
            counts = [0] * (max_value + 1) #Chaque case correspond à une valeur possible du chiffre qu’on est en train de trier, parmi les triplets.
    
            #compte les triplets dont le n-ième élément est égal à max_value 
            for triplet, position in R:
                counts[triplet[chiffre]] += 1
    
            # Transforme les comptes en positions finales dans le tableau trié.
            for i in range(1, len(counts)):
                counts[i] += counts[i - 1]
    
            sorted_R = [None] * len(R)

            # Place les éléments dans l'ordre, en préservant le tri précédent.
            for record in reversed(R):
                triplet, position = record
                value = triplet[chiffre]
                counts[value] -= 1
                sorted_R[counts[value]] = record
    
            #le prochain passage trie le résultat du précédent.
            R = sorted_R
    
        self.Rsorted = [triplet for triplet, position in R]
        self.index = [position for triplet, position in R]
        rang_par_triplet = {
        triplet: rang
        for rang, triplet in enumerate(dict.fromkeys(self.Rsorted), start=1)
        }
        self.order = [rang_par_triplet[triplet] for triplet, position in self.Rsorted]
        self.T_prime = [rang_par_triplet[triplet] for triplet, position in self.R]

        return self.Rsorted



class DC3:

    def __init__(self, T):
        T = [ord(char) for char in T]
        n = len(T)
        T += [0,0,0]
        self.P_0 = [i for i in range(n) if i % 3 == 0]
        P_1 = [i for i in range(n) if i % 3 == 1]
        P_2 = [i for i in range(n) if i % 3 == 2]
        self.P_12 = self.P_1 + self.P_2 #concaténation P_12s

        self.triplets = Triplets(T, self.P_12)
        T_prime = self.triplets.T_prime
        R = self.triplets.R

        if len(set(T_prime)) == len(T_prime): #vérifier que les rangs soient tous distincts
            index_0_1_2_prime = sorted(
        range(len(T_prime)),
        key=lambda i: T_prime[i]
    )
            suffix_array_12 = [self.P_12[i] for i in index_0_1_2_prime]
        else :  
            P_0 = DC3(T_prime).P_0
            P_12 = DC3(T_prime).P_12
            R_12_prime = Triplets(T_prime, P_12).Rsorted
            R_0_prime = [(tuple(T_prime[i:i+2]),i) for i in P_0]
            R_0_prime_sorted = Triplets.radix_sort(R_0_prime)
            rang_12 = {
            indice: rang
            for rang, indice in enumerate(R_12_prime)
            }

            i = j = 0
            R_merged = []
            while i < len(R_0_prime_sorted) and j < R_12_prime :
                a = R_0_prime_sorted[i]
                b = R_12_prime[j]
                if T_prime[a] < T_prime[b] : 
                    R_merged.append(a)
                    i += 1
                elif T_prime[b] < T_prime[a]:
                    R_merged.append(b)
                    j += 1
                elif T_prime[a] == T_prime[b] and b%3 == 1:
                    if rang_12[a+1] < rang_12[b+1]:
                        R_merged.append(a)
                        i += 1
                    else : 
                        R_merged.append(b)
                        j += 1
                elif T_prime[a] == T_prime[b] and b%3 == 2:
                    if T_prime[a+1] < T_prime[b+1]:
                        R_merged.append(a)
                        i += 1
                    elif T_prime[b+1] < T_prime[a+1]:
                        R_merged.append(b)
                        j += 1 
                    elif T_prime[a+1] == T_prime[b+1]:
                        if rang_12[a+2] < rang_12[b+2]:
                            R_merged.append(a)
                            i += 1
                        else : 
                            R_merged.append(b)
                            j += 1 
                    

                
            #merge R_12' et R_0'

        return suffix_array 


            