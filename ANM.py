from prody import *
from pylab import *
import sys

#construct the ANM network

def readXYZ(frameFile):
    
    # read in all frame coordinates & align each frame to ref
    Nframe = 0  # number of frames
    atom_count = 0
    frame_coord = []  # coordinates for one frame
    frame_list = []   # a list of frames

    frame_fh=open(frameFile)
    
    lines = []
    for line in frame_fh.readlines():
        lines.append(line)
    #print(lines)
    atom_count=int(lines[0])
    print(atom_count)
    Nframe=int(len(lines)/(atom_count+2))
    print(Nframe)
    for i in range(Nframe):
        #print(i)
        for j in range(i*(atom_count+2)+2,(i+1)*(atom_count+2)):
            #print(j)
            temp=lines[j].split()
            frame_coord.append( [ float(temp[1]), float(temp[2]), float(temp[3]) ])
        frame_list.append(frame_coord)
        frame_coord=[]

    frame_fh.close()
    return frame_list


#ion()

def bulidANM():
    

#np.set_printoptions(threshold=10000000)

    #pdbdata = parsePDB('')

    #for all atom
    #ag = parsePDB("1D1D.pdb")

    ag = parsePDB("HNS_bou_min_CG3.pdb")


    #writeArray("mdm2_coords.txt", ag.getCoords())


    #ag2 = AtomGroup()

    #coords2 = parseArray("mdm2_coords.txt")

    #ag2.setCoords(coords2)

    #print(ag2)

    #for CG
    calphas2 = ag.select('all')

    #for all atom
    #calphas2 = ag.select('calpha')

    print(type(calphas2))

    anm = ANM('8E3X')

    anm.buildHessian(calphas2,cutoff=15,gamma=1.0)

#print(anm.getHessian().round(3))

    print(anm.getCutoff())

    print(anm.getGamma())

#print(anm.getModel())
    anm.calcModes(n_modes=20)

    print(anm.getEigvals().round(4))

#print(anm.getEigvecs().round(3))

    anm.getCovariance().round(2)

    slowest_mode = anm[0]

    slowest_mode.getEigval().round(3)

    slowest_mode.getEigvec().round(3)

    writeNMD('./HNS_bou_min_CG3.nmd', anm, calphas2)


if __name__ == "__main__":

    #(junk,frameXYZ)=sys.argv
    
    #frames=readXYZ(frameXYZ)  
    
    bulidANM()
