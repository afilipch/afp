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

### Dictionaries

countries = ['spain', 'france', 'germany', 'norway']
capitals = ['madrid', 'paris', 'berlin', 'oslo']
# From string in countries and capitals, create dictionary europe
europe = { 'spain':'madrid', 'france':'paris', 'germany':'berlin', 'norway':'oslo' }
europe = dict(zip(countries, capitals))
#print(list(europe.items())[2])
#print(europe.keys())
#print(europe.values())

europe['italy'] = 'rome'
#print('italy' in europe)  =  print('italy' in europe.keys())  - not in values!!

for key, value in list(europe.items()):
    if(value == 'vienna'):
        del(europe[key]);
        
europe = dict([ (key, value) for key, value in europe.items() if (value != 'vienna')])      
europe = dict([ x for x in europe.items() if (x[1] != 'vienna')])

#a = list(range(5,8))   ===  [5, 6, 7]
#b = list(zip(a[::2],a[1::2]))  - 0, 1; 2, 3; 4, 5 etc.
