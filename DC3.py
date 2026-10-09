class Triplets:

    def __init__(self, T, P):
        self.R = [(tuple(T[i:i + 3]), i) for i in P]
        self.radix_sort()

    def radix_sort(self):
        R = self.R
        for chiffre in (2, 1, 0): 
            max_value = max(triplet[chiffre] for triplet, _ in R)
            



class DC3:

    def __init__(self, T):
        self.n = len(T)
        self.P_0 = [i for i in range(self.n + 1) if i % 3 == 0]
        self.P_1 = [i+1 for i in self.P_0 if i<self.n] #P_0 + 1
        self.P_2 = [i+1 for i in self.P_1 if i<self.n] #P_1 + 1
        self.P_12 = self.P_1 + self.P_2 #concaténation P_12s