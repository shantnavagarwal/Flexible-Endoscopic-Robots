
num_elems=10
len=5
inter=len/num_elems
theta= []
j=1


data=""   
for j in range(1,num_elems+1):
    data += "PLBEAM " + str(j) + " " + str(j*2-1) + " " + str(j*2) + " " + str(j*2+1) + " " + str(j*2+2) + "\n"

j=1
data+='\n'
init_X=0

for j in range(num_elems+1):
    data+="X "+str(j*2+1)+" "+str(init_X+j*inter)+" "+str(0.0)+"\n"
    
data+="\nFIX 2"+"\n"
data+="INPUTX 1\n\n"

for j in range(1,num_elems+1):
    data+="DYNE "+str(j)+" 1 2 3\n"

data+="\nEND\nHALT\n"

for j in range(1,num_elems+1):
    data+="EM "+str(j)+" 0.001532   0.00000000002393\n"
data+="\n"
for j in range(1,num_elems+1):
    data+="ESTIFF "+str(j)+" 69300  0.001614\n"
data+="\n"
for j in range(1,num_elems+1):
    data+="EDAMP "+str(j)+" 0.000001 0.0000001\n"
data+="\n"


data+="INPUTX 1 1 0 2 0\nUSERSIG circ_pot\nTIMESTEP 5 8000\n\n"

for j in range(1,num_elems+1):
    data+="STARTDE "+str(j)+" 1 0 0\n"
    data+="STARTDE "+str(j)+" 2 0 0\n"
    data+="STARTDE "+str(j)+" 3 0 0\n"
data+="\n"

data+="END\nEND"

with open('SingleFF_10E.dat', 'w') as file:
    file.write(data)
    
    