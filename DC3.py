class P :

    def __init__(self, T):
        self.n = len(T)
        self.P_0 = [i for i in range(n+1) if i%3 == 0]
        self.P_1 = [i+1 for i in self.P_0 if i<self.n] #P_0 + 1
        self.P_2 = [i+1 for i in self.P_1 if i<self.n] #P_1 + 1
        self.P_12 = self.P_1 + self.P_2 #concaténation P_12s