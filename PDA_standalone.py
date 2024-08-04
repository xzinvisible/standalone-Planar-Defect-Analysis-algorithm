############################################################
#    pda index      0      1     2     3     4     5    6  #
# structure type   other  fcc   ISF   HCP   ESF   TB   BCC #
############################################################


from ovito.io import import_file
from ovito.modifiers import CommonNeighborAnalysisModifier
from ovito.data import *
from time import time
import os
os.environ["OVITO_THREAD_COUNT"] = "1"

# input & output files
fileinput = 'test.xyz'
fileoutput = 'pda.txt'

# number of processes
cores = 24 

pipeline = import_file(fileinput)   
pipeline.modifiers.append(CommonNeighborAnalysisModifier())

def PDAcount(frame):
    data = pipeline.compute(frame)
    finder = NearestNeighborFinder(12, data)
    pda = [0,0,0,0,0,0,0]
    step = 100*frame

    def ESF_TWIN(k, neighbor):
        for i in neighbor:
            index = i
            flag1 = 0
            fcc, hcp = 0, 0
            
            if data.particles['Structure Type'][index] == 1:               
                for neigh in finder.find(index):
                    if data.particles['Structure Type'][neigh.index] == 1:
                        fcc += 1
                    elif data.particles['Structure Type'][neigh.index] == 2:
                        hcp += 1
                if (5<= fcc <= 6) and (5<= hcp <= 6):
                    pda[4] += 1
                    flag1 = 1
                    break

        if flag1 != 1:
            pda[5] += 1            
    
    def HCP(index):
        fcc, hcp =0, 0
        neighbor = []
        flag = 0

        for neigh in finder.find(index):
            neighbor.append(neigh.index)
            if data.particles['Structure Type'][neigh.index] == 1:
                fcc += 1
            elif data.particles['Structure Type'][neigh.index] == 2:
                hcp += 1

        if hcp >= 11:
            pda[3] += 1
            flag = 1
        elif (2<= fcc <= 3) and ( 8<= hcp <= 9):
            if flag != 1:
                pda[2] += 1  
        elif (5<= fcc <= 6) and ( 5<= hcp <= 6):
            ESF_TWIN(index, neighbor)  
    
    for index in range(data.particles.count):
        #yield(index / data.particles.count)
        if data.particles['Structure Type'][index] == 2:
            HCP(index)
        elif data.particles['Structure Type'][index] == 1:
            pda[1] += 1
        elif data.particles['Structure Type'][index] == 3:
            pda[6] += 1
        else:
            pda[0] += 1
    
    pda.append(step)
    print(f'{step:7d} {pda[0]:7d} {pda[1]:7d} {pda[2]:7d} {pda[3]:7d} {pda[4]:7d} {pda[5]:7d} {pda[6]:7d}')
    return pda
    

if __name__ == '__main__':

    import multiprocessing as mp
    mp.set_start_method('spawn')

    t_start = time()

    print('step    other   FCC     ISF     HCP     ESF     TB      BCC\n')
    frame = range(pipeline.source.num_frames)

    with mp.Pool(cores) as pool:
        PDAList = list(pool.map(PDAcount,frame))

    t_end = time()

    with open(fileoutput,'w') as f:
        f.write(f'Planar Defect Analysis for {fileinput}: \n')
        f.write('step    other   FCC     ISF     HCP     ESF     TB      BCC\n')
        for i in range(len(PDAList)):
            f.write(f'{PDAList[i][7]:7d} {PDAList[i][0]:7d} {PDAList[i][1]:7d} {PDAList[i][2]:7d} {PDAList[i][3]:7d} {PDAList[i][4]:7d} {PDAList[i][5]:7d} {PDAList[i][6]:7d}\n')

    print(f"Computation took {t_end - t_start} seconds")


    
    
