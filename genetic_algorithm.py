import numpy as np


class Genetic:

    num_points = 0
    population_size = 0

    intervals = [[0,0]]
    counts = [0] 

    atomic_intervals = [[0,0]]

    member_matrix=np.ndarray((0,0)) #atomic intervals -> intervals. 1 if the interval contains the atomic interval

    population=np.array([[]])

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
        return

    def getbest(self):
        return self.population[-1].tolist()
    def getcounts(self):
        individual = self.population[-1]
        return individual @ self.member_matrix.T

    def step(self):
        #NOTE: Population must be sorted by fitness before calling this function
        
        NUM_SELECT = self.population_size//2

        selected = self.population[self.population_size-NUM_SELECT:]
        

        #have children


        children = np.ndarray(self.population.shape)



        pvals = [i/sum(range(1,NUM_SELECT+1)) for i in range(1,NUM_SELECT+1)]

        # 2 parents, can try more
        parentdist = np.random.multinomial(2,pvals, self.population_size)
        
        for i in range(children.shape[0]):

            #mask is slow
            parentinds = np.nonzero(parentdist[i])[0]
            mask = np.random.randint(0,2,children.shape[1]) == 0
            children[i] = np.where(mask,selected[parentinds[0]],selected[parentinds[-1]])

            # pivot = np.random.randint(0,children.shape[1])
            # children[i, :pivot] = selected[parentinds[0], :pivot]
            # children[i, pivot:] = selected[parentinds[-1], pivot:]


            #mutation
            if(np.random.randint(0,3)==0):
                for idx in range(children.shape[1]):
                    if(np.random.randint(0,5) == 0):
                        offs = round(np.random.normal(0, self.num_points/100))
                        children[i,idx] += offs
                        if children[i,idx] <0:
                            children[i,idx]=0


            #fix number of points

            toadd = self.num_points-sum(children[i])
            pi = np.random.permutation(children.shape[1])
            for j in pi:
                children[i,j] += toadd//children.shape[1]
                if children[i,j] < 0:
                    children[i,(j+1)%children.shape[1]] += children[i,j]
                    children[i,j]=0

            children[i,np.random.randint(0,children.shape[1])] += toadd%children.shape[1]


        #sort children based on fitness

        population_counts = []
        fitness = []

        #compute counts and fitness

        for individual in children:
            population_counts.append(individual @ self.member_matrix.T ) 

        
        for curcounts in population_counts:
            fitness.append(-sum([(a-b)*(a-b) for (a,b) in zip(curcounts,self.counts)]))


        #select
        idx = np.argsort(fitness)
        self.population = children[idx]

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

genetic = Genetic(NUM_POINTS,intervals, counts, 20)
genetic.init_population()


for i in range(50000):

    input()
    genetic.step()
    expected = genetic.counts
    actual = genetic.getcounts()

    print("OG:", data)
    print(genetic.getbest())
    print(expected)
    print(actual)
    mse= sum([(a-b)*(a-b) for (a,b) in zip(expected,actual)])
    print("MSE: ", mse)
    if mse == 0:
        break



expected = genetic.counts
actual = genetic.getcounts()

print("OG:", data)
print(genetic.getbest())
print(expected)
print(actual)
mse= sum([(a-b)*(a-b) for (a,b) in zip(expected,actual)])
print("MSE: ", mse)

