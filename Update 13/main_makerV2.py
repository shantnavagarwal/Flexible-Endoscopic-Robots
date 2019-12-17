import numpy as np


# import pandas as pd

def func(lon_stiff, ben_stiff):
    dat = np.genfromtxt('large_intes_discret.csv', delimiter=',')
    dat = np.nan_to_num(dat)
    #    print(dat)

    data = """# Endoscope \n"""
    # Endoscope
    edsx = 0.0
    edsy = -0.6
    edex = 0.0
    edey = -0.0
    j = 1
    edj = j
    k = 1
    edk = k
    n = 40
    edn = n
    for i in range(n):
        data += "PLBEAM " + str(j) + " " + str(k) + " " + str(k + 1) + " " + str(k + 2) + " " + str(k + 3) + '\n'
        j += 1
        k += 2
    # PLBEAM 1 1 2 3 4
    # PLBEAM 2 5 6 7 8
    # PLBEAM 3 7 8 9 10
    # PLBEAM 4 9 10 11 12
    # Large Intestine
    data += "# Large Intestine \n"
    lij = j
    k += 2
    lik = k
    n = dat.shape[0] - 1
    lin = n
    for i in range(n):
        data += "PLBEAM " + str(j) + " " + str(k) + " " + str(k + 1) + " " + str(k + 2) + " " + str(k + 3) + '\n'
        j += 1
        k += 2
    # PLBEAM 2 5 6 7 8
    # PLBEAM 3 7 8 9 10
    # PLBEAM 4 9 10 11 12

    data += "\n"
    data += "# Fix nodes of Endoscope \n"
    #    data+="X 1 .1 -0.1\n"
    #    data+="X 3 .1 -1\n"
    j = edj
    n = edn
    # edsteps = (edex - edsx) / n
    # eddatx = np.arange(edsx, edex, edsteps)
    eddatx = np.ones(n+1)*edsx
    edstepy = (edey - edsy) / n
    eddaty = np.arange(edsy, edey + edstepy, edstepy)
    for i in range(n + 1):
        data += "X " + str(j) + " " + str(eddatx[i]) + " " + str(eddaty[i]) + '\n'
        j += 2
    data += "\n# Fix nodes of Large Intestine \n"
    j = lij*2+1
    n = lin
    for i in range(n + 1):
        data += "X " + str(j) + " " + str(dat[i, 0]) + " " + str(dat[i, 1]) + '\n'
        j += 2

    data += """

FIX 2
FIX """ + str(lik) + "\n"

    data += "FIX " + str(lik + 28) + "\n"
    data += "FIX " + str(lik + 42)
    data += """
    
INPUTX 1

DYNE 1 2 3
"""
    # Defining DYNE for endoscope nodes
    j = edj+1
    n = edn
    for i in range(j, n+1):
        data += "DYNE " + str(i) + " 1 2 3\n"
    # Defining DYNE for Large Intestine
    j = lij
    n = lin
    for i in range(j, j+n):
        if i==lij or i==lij+14 or i==lij+19:
            data += "DYNE "+str(i)+" 2 3\n"
            data += "RLSE "+str(i)+" 1\n"
        else:
            data += "DYNE " + str(i) + " 1 2 3\n"

    #data += "DYNE " + str(j + n - 1) + " 2 3\n"
    #data += "RLSE " + str(j + n - 1) + " 1\n"
    # DYNE 2 1 2 3
    # DYNE 3 1 2 3
    # DYNE 4 2 3
    # RLSE 4 1

    data += """
END
HALT

"""
    # EM for all elements
    for i in range(1, lin + edn+1):
        data += "EM " + str(i) + " 0.001532   0.00000000002393\n"
    # EM 1 0.001532   0.00000000002393
    # EM 2 0.001532   0.00000000002393
    # EM 3 0.001532   0.00000000002393
    # EM 4 0.001532   0.00000000002393
    data += '\n'
    # Stiffness for endoscope
    for i in range(1, edn+1):
        data += "ESTIFF " + str(i) + " 300  0.0015\n"
    for i in range(lij, edn+lin+1):
        data += "ESTIFF " + str(i) + " " + str(lon_stiff) + " " + str(ben_stiff) + "\n"
    data += "\n"
    #for i in range(1, edn + 1):
     #   data += "EDAMP " + str(i) + " 0.000001 0.00001\n"
    #for i in range(edn+1, edn + lin+1):
     #   data += "EDAMP " + str(i) + " 0.000001 0.000001\n"

    #    """
    #
    #    data+="ESTIFF 2 "+str(lon_stiff)+" "+str(ben_stiff)+"\n"
    #    data+="ESTIFF 3 "+str(lon_stiff)+" "+str(ben_stiff)+"\n"
    #    data+="ESTIFF 4 "+str(lon_stiff)+" "+str(ben_stiff)+"\n"

    data += """
    
INPUTX 1 1 """ +str(edsx) + """ 0 0
INPUTX 1 2 """ +str(edsy) + """ .1 0

USERSIG CIRC_POT
TIMESTEP 5 20000

"""
    for i in range(1, edn + lin + 1):
        data += "STARTDE " + str(i) + " 1 0 0\n"
        data += "STARTDE " + str(i) + " 2 0 0\n"
        data += "STARTDE " + str(i) + " 3 0 0\n"
    # STARTDE 1 1 0 0
    # STARTDE 1 2 0 0
    # STARTDE 1 3 0 0
    # STARTDE 2 1 0 0
    # STARTDE 2 2 0 0
    # STARTDE 2 3 0 0

    data += """
    
END
END
    """

    with open('trial_code.dat', 'w') as file:
        file.write(data)


# %%

func(700, 0.003)
