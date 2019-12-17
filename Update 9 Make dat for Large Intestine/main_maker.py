import numpy as np
#import pandas as pd

def func(lon_stiff,ben_stiff):
    
    dat=np.genfromtxt('large_intes_discret.csv',delimiter=',')
    dat=np.nan_to_num(dat)
#    print(dat)
    
    data="""
PLBEAM 1 1 2 3 4
"""
    j=2
    k=5
    n=dat.shape[0]-1
    for i in range(n):
        data+="PLBEAM "+str(j)+" "+str(k)+" "+str(k+1)+" "+str(k+2)+" "+str(k+3)+'\n'
        j+=1
        k+=2
    #PLBEAM 2 5 6 7 8
    #PLBEAM 3 7 8 9 10
    #PLBEAM 4 9 10 11 12
    
    data+='\n'
    data+="X 1 "+"0.5 0\n"
#    data+="X 3 "+str(dist)+" "+str(-0.25-depth)+"\n"
    j=5
    
    for i in range(n+1):
        data+="X "+str(j)+" "+str(dat[i,0])+" "+str(dat[i,1])+'\n'
        j+=2
    
    data+="""

FIX 2
FIX 5
"""

    data+="FIX "+str(j-2)

    data+="""
    
INPUTX 1

DYNE 1 2 3
"""
    for i in range(n-1):
        data+="DYNE "+str(i+2)+" 1 2 3\n"
    
    data+="DYNE "+str(n+1)+" 2 3\n"
    data+="RLSE "+str(n+1)+" 1\n"
#DYNE 2 1 2 3
#DYNE 3 1 2 3
#DYNE 4 2 3 
#RLSE 4 1


    data+="""
END
HALT

"""
    for i in range(1,n+2):
        data+="EM "+str(i)+" 0.001532   0.00000000002393\n"
#EM 1 0.001532   0.00000000002393
#EM 2 0.001532   0.00000000002393
#EM 3 0.001532   0.00000000002393
#EM 4 0.001532   0.00000000002393
    data+='\n'
    data+="ESTIFF 1 6930  0.003614\n"
    for i in range(2,n+2):
        data+="ESTIFF "+str(i)+" "+str(lon_stiff)+" "+str(ben_stiff)+"\n"
    data+="\n"
    for i in range(1,n+2):
        data+="EDAMP "+str(i)+" 0.001 0.001\n"
        
#    """
#    
#    data+="ESTIFF 2 "+str(lon_stiff)+" "+str(ben_stiff)+"\n"
#    data+="ESTIFF 3 "+str(lon_stiff)+" "+str(ben_stiff)+"\n"
#    data+="ESTIFF 4 "+str(lon_stiff)+" "+str(ben_stiff)+"\n"
    
    
    data+="""
    
INPUTX 1 1 0.5 0 0
INPUTX 1 2 0 0 0

USERSIG CIRC_POT
TIMESTEP 1 200

"""
    for i in range(1,n+2):
        data+="STARTDE "+str(i)+" 1 0 0\n"
        data+="STARTDE "+str(i)+" 2 0 0\n"
        data+="STARTDE "+str(i)+" 3 0 0\n"
#STARTDE 1 1 0 0
#STARTDE 1 2 0 0
#STARTDE 1 3 0 0
#STARTDE 2 1 0 0
#STARTDE 2 2 0 0
#STARTDE 2 3 0 0


    data+="""
    
END
END
    """
    
    with open('trial_code.dat', 'w') as file:
        file.write(data)
    
#%%

func(3500, 0.003)
