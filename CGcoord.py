import numpy as np
import sys

#this code is used to get coarse grained coordinates of backbone


def readXYZ(frameFile):

    # read in all frame coordinates & align each frame to ref
    Nframe = 0  # number of frames
    atom_count = 0  
    frame_coord = []  # coordinates for one frame
    frame_list = []   # a list of frames
    
    frame_fh = open(frameFile)
    for line in frame_fh.readlines():
        #print(line)       
        #if 
        if line.strip():  # remove the whitespace at the beginning & end; if not empty
            temp = line.split()  # split the line, and put the words into a list call "temp"
            frame_coord.append( [ float(temp[0]), float(temp[1]), float(temp[2]) ])
            atom_count = atom_count + 1
        else:
            frame_list.append(frame_coord)
            Nframe = Nframe + 1
            frame_coord = []  # reset for a new frame
            atom_count = 0
                
    frame_fh.close()
    #print(frame_list)
    return frame_list

def readBoundari(boundaryFile):
    boundry=np.loadtxt(boundaryFile,delimiter=',',dtype=int)
    return boundry



def CGCoord(coord, boundary):
    #print(boundary)    
    FList=[]
    for frame in coord:
        Fcood=[]
        for i in range(len(boundary)):
            x=0
            y=0
            z=0
            count=0
            if(i==0):
                #for(2,2): 0 1,,,, 2 3 4,,, 5 6 7,,,, (1,4,7,10)
                for j in range(0,boundary[i]+1):
                    x=x+frame[j][0]
                    y=y+frame[j][1]
                    z=z+frame[j][2]
                    count=count+1         

            elif(i>0):
                for j in range(boundary[i-1]+1,boundary[i]+1):
                    x=x+frame[j][0]
                    y=y+frame[j][1]
                    z=z+frame[j][2]
                    count=count+1
            
            x=x/count
            y=y/count
            z=z/count
            Fcood.append([x,y,z])
        
        FList.append(Fcood)

    return FList

def printCoord(Flist):
    for frame in Flist:
        index=1
        print(len(frame))
        print('text')
        for cood in frame:
            print(index,'\t',cood[0],'\t',cood[1],'\t',cood[2])
            index=index+1


if __name__ == "__main__":

    (junk, trjfile, boundary) = sys.argv
    

    gettrj=readXYZ(trjfile)
    getBound=readBoundari(boundary)
    #print('1')
    result=CGCoord(gettrj,getBound)
    #print(result)
    printCoord(result)
