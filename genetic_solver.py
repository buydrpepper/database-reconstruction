import numpy as np
import time


class Genetic:

    num_points = 0
    population_size = 0

    intervals = [[0,0]]
    counts = [0] 

    atomic_intervals = [[0,0]]

    member_matrix=np.ndarray((0,0)) #atomic intervals -> intervals. 1 if the interval contains the atomic interval

    population=np.ndarray((0,0))

    def __init__(self, num_points, intervals, counts, population_size):
        self.num_points=num_points
        self.intervals=intervals
        self.counts=counts
        self.population_size=population_size

        endpoints = sorted(set([i for pair in intervals for i in pair]))
        self.atomic_intervals= [[endpoints[i],endpoints[i+1]] for i in range(len(endpoints)-1) ]

        self.member_matrix = np.zeros((len(self.intervals),len(self.atomic_intervals)))

        for i in range(len(self.atomic_intervals)):
            curinterval = self.atomic_intervals[i]
            for j in range(len(intervals)):
                if curinterval[0]>= intervals[j][0] and curinterval[1]<= intervals[j][1]:
                    self.member_matrix[j][i] = 1
        return



    def init_population(self):
        pvals = [abs((tuple[1]-tuple[0])/(self.atomic_intervals[-1][-1]-self.atomic_intervals[0][0])) for tuple in self.atomic_intervals]
        self.population=np.random.multinomial(self.num_points,pvals, self.population_size)
        self.update_population(self.population)
        return

    def getbest(self):
        return self.population[-1].tolist()
    def getcounts(self):
        individual = self.population[-1]
        return individual @ self.member_matrix.T

    def update_population(self, children):
        fitness = np.ndarray((self.population_size))
        curcounts_m = children @ self.member_matrix.T
        fitness = -np.sum((curcounts_m - self.counts)**2, axis=1)
        pi = np.argsort(fitness)
        self.population = children[pi]

        return

    def step(self):
        #NOTE: Population must be sorted by fitness before calling this function

        children = np.ndarray(self.population.shape)

        # 2 parents, can try more
        for i in range(children.shape[0]):
            parentinds = [0]*2
            for pind in range(len(parentinds)):
                #Using the roulette method
                while True:
                    cand = np.random.randint(0,self.population_size)
                    if np.random.rand() < (cand+1)/(self.population_size):
                        break
                parentinds[pind] = cand

            alpha = np.random.rand()
            children[i] = np.round(self.population[parentinds[0]]*alpha + self.population[parentinds[1]]*(1-alpha))
            #delta = self.num_points - sum(children[i])

            #TODO: TOMORROW: FIX MUTATION, AND COMPLETE CROSSOVER

            #mutation
            #TODO: definitely need some kind of swap or permutation.

            if(np.random.randint(0,3)==0):
                for idx in range(children.shape[1]):
                    if(np.random.randint(0,8) == 0):
                        offs = abs(round(np.random.normal(0, self.num_points/100)))
                        if children[i,idx] < offs:
                            offs = children[i,idx]
                        children[i,idx] -= offs
                        children[i, np.random.randint(0,children.shape[1])] += offs


            #fix number of points
            #the choice to add uniformly is completely arbitrary, so this can be improved

            pi = np.argsort(children[i])
            toadd = self.num_points-sum(children[i])

            num_skipped = 0
            for idx in pi:
                curval = children[i, idx]

                if curval + toadd//(children.shape[1] - num_skipped) > 0:
                    break

                children[i, idx] = 0
                toadd += curval
                num_skipped += 1

            for j in range(num_skipped,children.shape[1]):
                children[i, pi[j]] += toadd//(children.shape[1]-num_skipped)

            children[i,np.random.randint(0,children.shape[1])] += toadd%(children.shape[1]-num_skipped)

        self.update_population(children)

        return
    





NUM_POINTS = 10000
pvals = [1,3,10,7,1,2,3,4,0,9]
pvals = [i/sum(pvals) for i in pvals]
#data = np.random.multinomial(NUM_POINTS,pvals)
data = [ 250 , 760, 2593, 1725 , 269 , 496 , 739  ,980,    0, 2188]
intervals = [[0,4],[1,5],[2,6],[3,7],[4,8],[5,9],[6,10],[0,10]]
counts=[sum(data[a:b]) for (a,b) in intervals]

print(data)
print(intervals)
print(counts)

genetic = Genetic(NUM_POINTS,intervals, counts, 100)
genetic.init_population()


for i in range(50000):
    
    genetic.step()
    expected = genetic.counts
    actual = genetic.getcounts()

    # print("OG:", data)
    # print(genetic.getbest())
    # print(expected)
    # print(actual)
    #mse= sum([(a-b)*(a-b) for (a,b) in zip(expected,actual)])
    #print("MSE: ", mse)
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

