import numpy as np


class Genetic:

    num_points = 0
    population_size = 0

    intervals = [[0,0]]
    counts = [0] 

    atomic_intervals = [[0,0]]

    member_matrix=np.ndarray((0,0),dtype=np.int64) #atomic intervals -> intervals. 1 if the interval contains the atomic interval

    population=np.ndarray((0,0),dtype=np.int64)

    def __init__(self, num_points, intervals, counts, population_size):
        self.num_points=num_points
        self.intervals=intervals
        self.counts=counts
        self.population_size=population_size

        endpoints = sorted(set([i for pair in intervals for i in pair]))
        self.atomic_intervals= [[endpoints[i],endpoints[i+1]] for i in range(len(endpoints)-1) ]

        self.member_matrix = np.zeros((len(self.intervals),len(self.atomic_intervals)), dtype=np.int64)

        for i in range(len(self.atomic_intervals)):
            curinterval = self.atomic_intervals[i]
            for j in range(len(intervals)):
                if curinterval[0]>= intervals[j][0] and curinterval[1]<= intervals[j][1]:
                    self.member_matrix[j][i] = 1
        print(self.member_matrix)
        return



    def init_population(self):
        pvals = [abs((tuple[1]-tuple[0])/(self.atomic_intervals[-1][-1]-self.atomic_intervals[0][0])) for tuple in self.atomic_intervals]
        self.population=np.random.multinomial(self.num_points,pvals, self.population_size).astype(np.int64)
        self.update_population(self.population)
        return

    def getbest(self):
        return self.population[-1].tolist()
    def getcounts(self):
        individual = self.population[-1]
        return individual @ self.member_matrix.T

    def update_population(self, children):
        curcounts_m = children @ self.member_matrix.T
        fitness = -np.sum((curcounts_m - self.counts)**2, axis=1)
        pi = np.argsort(fitness)
        self.population = children[pi]

        return

    def step(self):
        #NOTE: Population must be sorted by fitness before calling this function

        children = np.ndarray(self.population.shape, dtype=np.int64)

        # 2 parents, can try more
        for i in range(children.shape[0]):
            parentinds = [0]*2
            for pind in range(len(parentinds)):
                #Using the fast roulette method
                while True:
                    cand = np.random.randint(0,self.population_size)
                    if np.random.rand() < (cand+1)/(self.population_size):
                        break
                parentinds[pind] = cand

            alpha = np.random.rand()
            children[i] = np.round(self.population[parentinds[0]]*alpha + self.population[parentinds[1]]*(1-alpha)).astype(np.int64)
            toadd = self.num_points-sum(children[i])

            assert(abs(toadd) <= children.shape[1])

            sign = (1,-1)[int(toadd<0)]

            pi = np.random.permutation(children.shape[1])
            while toadd != 0:
                for idx in pi:
                    if toadd == 0:
                        break
                    if children[i, idx] + sign >= 0:
                        children[i, idx] += sign
                        toadd -= sign

            #mutation
            if(np.random.randint(0,3)==0):

                if(np.random.randint(0,5) == 0): #swap mutation
                    i1 = np.random.randint(0,children.shape[1])
                    i2 = np.random.randint(0,children.shape[1])
                    children[i, i1], children[i,i2] = children[i,i2], children[i,i1]

                for idx in range(children.shape[1]):
                    if(np.random.randint(0,8) == 0): 
                        offs = abs(round(np.random.normal(0, self.num_points/100)))
                        if children[i,idx] < offs:
                            offs = children[i,idx]
                        children[i,idx] -= offs
                        children[i, np.random.randint(0,children.shape[1])] += offs

        #elitism
        children[-1] = self.population[-1]

        self.update_population(children)

        return
    





NUM_POINTS = 10000
pvals = [np.random.randint(0,1001) for i in range(10000)]
pvals = [i/sum(pvals) for i in pvals]
data = np.random.multinomial(NUM_POINTS,pvals)
intervals = [  [np.random.randint(0,5001), np.random.randint(0,5001)] for i in range(100)]
counts=[sum(data[a:b]) for (a,b) in intervals]

print(data)
print(intervals)
print(counts)

genetic = Genetic(NUM_POINTS,intervals, counts, 200)
genetic.init_population()


for i in range(50000):
    
    genetic.step()
    expected = genetic.counts
    actual = genetic.getcounts()

    # print("OG:", data)
    # print(genetic.getbest())
    # print(expected)
    # print(actual)
    # mse= sum([(a-b)*(a-b) for (a,b) in zip(expected,actual)])
    # print("MSE: ", mse)
    taxidist= sum(abs(expected-actual))

    print("taxidsit: ", taxidist)
    if taxidist == 0:
        break



expected = genetic.counts
actual = genetic.getcounts()

print("OG:", data)
print(genetic.getbest())
print(expected)
print(actual)
mse= sum([(a-b)*(a-b) for (a,b) in zip(expected,actual)])
print("MSE: ", mse)

