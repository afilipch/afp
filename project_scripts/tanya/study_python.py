import numpy as np


years = list(range(1950, 2100, 1))
pop = np.linspace(2, 15, 151)

#print(pop)

for a, y in zip(pop, years):
    if (a>10):
        print(y)
        break


print([x[1] for x in zip(pop,years) if x[0]>10][0])

#for nomer, a in enumerate(pop):
    #if (a>10):
        #print(year[nomer])
        #break
    
    
mylist = [1, 6, 4, 8]
result = [];
for a in mylist:
    if(a>5):
        result.append(a**2)
print(result)
        
result = [x**2 for x in mylist if x>5]
print(result)


### Central limit theorem -> normal distribution

