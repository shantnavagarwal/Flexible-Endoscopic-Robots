
def func(depth,dist,lon_stiff,ben_stiff):
    data="""
PLBEAM 1 1 2 3 4
PLBEAM 2 5 6 7 8
PLBEAM 3 7 8 9 10
PLBEAM 4 9 10 11 12
    
    """
    data+="X 1 "+str(dist)+" 0\n"
    data+="X 3 "+str(dist)+" "+str(-0.25-depth)+"\n"
    
    
    
    data+="""
X 5 0.0 -1.2
X 7 0.33 -1.2
X 9 0.66 -1.2
X 11 1.0 -1.2


FIX 2
FIX 5
FIX 11

INPUTX 1

DYNE 1 2 3
#DYNE 2 2 3

#DYNE 6 1 2 3
DYNE 2 1 2 3
DYNE 3 1 2 3
DYNE 4 2 3 
RLSE 4 1
#DYNE 10 2 3
#RLSE 10 1


END
HALT
EM 1 0.001532   0.00000000002393
EM 2 0.001532   0.00000000002393
EM 3 0.001532   0.00000000002393
EM 4 0.001532   0.00000000002393

ESTIFF 1 6930  0.003614
    """
    
    data+="ESTIFF 2 "+str(lon_stiff)+" "+str(ben_stiff)+"\n"
    data+="ESTIFF 3 "+str(lon_stiff)+" "+str(ben_stiff)+"\n"
    data+="ESTIFF 4 "+str(lon_stiff)+" "+str(ben_stiff)+"\n"
    
    
    data+=
    """
EDAMP 1 0.001 0.001
EDAMP 2 0.001 0.001

EDAMP 3 0.001 0.001
EDAMP 4 0.001 0.001



INPUTX 1 1 0.5 0 0
INPUTX 1 2 0 0 0

USERSIG CIRC_POT
TIMESTEP 1 200

STARTDE 1 1 0 0
STARTDE 1 2 0 0
STARTDE 1 3 0 0
STARTDE 2 1 0 0
STARTDE 2 2 0 0
STARTDE 2 3 0 0



END
END
    """
    
    with open('main_code.dat', 'w') as file:
        file.write(data)