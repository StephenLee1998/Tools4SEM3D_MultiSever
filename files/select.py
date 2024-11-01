f1 = open('SOURCES_SELECTED','r')

lines = f1.readlines()
data = []

for line in lines:
    data.append(line)

for i,(d) in enumerate(data):
    if i%3 == 0:
        print(d.split()[0]+' '+d.split()[1]+' '+d.split()[2])

f1.close()