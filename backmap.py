import numpy as np
import sys


#backmap from 6 to 3, 3 to 1 res./site

def readXYZ(frameFile):

    # read in all frame coordinates & align each frame to ref
    Nframe = 0  # number of frames
    atom_count = 0  
    frame_coord = []  # coordinates for one frame
    frame_list = []   # a list of frames
    
    frame_fh = open(frameFile)
    for line in frame_fh.readlines():
        
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
    return frame_list

def readBoundari(boundaryFile):
    boundry=np.loadtxt(boundaryFile,delimiter=',',dtype=int)
    return boundry



def CGCoord(coord, boundary):
    
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


def backmapCG6toCG3(Flist):

    #new coordinate 
     # coordinates for one frame
    frame_list = []   # a list of frames


    for frame in Flist:
        f1=0.3
        f2=0.7
        frame_coord = [] 
        for i in range(len(frame)):
            if i==0:
                tempx=f1*frame[i][0]
                tempy=f1*frame[i][1]
                tempz=f1*frame[i][2]
                frame_coord.append([tempx,tempy,tempz])
                
                tempx=f2*frame[i][0]
                tempy=f2*frame[i][1]
                tempz=f2*frame[i][2]
                frame_coord.append([tempx,tempy,tempz])


            elif i>0:
                #this could be fixed
                tempx=f2*frame[i-1][0]+f1*frame[i][0]
                tempy=f2*frame[i-1][1]+f1*frame[i][1]
                tempz=f2*frame[i-1][2]+f1*frame[i][2]
                frame_coord.append([tempx,tempy,tempz])

                tempx=f1*frame[i-1][0]+f2*frame[i][0]
                tempy=f1*frame[i-1][1]+f2*frame[i][1]
                tempz=f1*frame[i-1][2]+f2*frame[i][2]
                frame_coord.append([tempx,tempy,tempz])

        frame_list.append(frame_coord)
    
    return frame_list


def backmapCG3toCG1(Flist):

    #new coordinate 
    # coordinates for one frame
    frame_list = []   # a list of frames
    f1=0.4
    f2=0.6
    #print(Flist)
    for frame in Flist:
        frame_coord = []
        for i in range(len(frame)):
            if i==0:
                #print(frame[i])
                tempx=0.3*frame[i][0]
                tempy=0.3*frame[i][1]
                tempz=0.3*frame[i][2]
                frame_coord.append([tempx,tempy,tempz])
                
                tempx=0.7*frame[i][0]
                tempy=0.7*frame[i][1]
                tempz=0.7*frame[i][2]
                frame_coord.append([tempx,tempy,tempz])

                tempx=1*frame[i][0]
                tempy=1*frame[i][1]
                tempz=1*frame[i][2]
                frame_coord.append([tempx,tempy,tempz])

            elif i>0:
                #this could be fixed
                tempx=f1*frame[i-1][0]+f2*frame[i][0]
                tempy=f1*frame[i-1][1]+f2*frame[i][1]
                tempz=f1*frame[i-1][2]+f2*frame[i][2]
                frame_coord.append([tempx,tempy,tempz])

                tempx=f2*frame[i-1][0]+f1*frame[i][0]
                tempy=f2*frame[i-1][1]+f1*frame[i][1]
                tempz=f2*frame[i-1][2]+f1*frame[i][2]
                frame_coord.append([tempx,tempy,tempz])

                tempx=1*frame[i][0]
                tempy=1*frame[i][1]
                tempz=1*frame[i][2]
                frame_coord.append([tempx,tempy,tempz])
                
        frame_list.append(frame_coord)
    
    return frame_list

def RMSD(list1,list2):
    Alltemp=[]
    Allcount=0
    #print(list2)
    temp=0
    for i in range(len(list1)):
        #print(len(list1),len(list2))
        #print()
        Allcount=Allcount+1
        temp=0
        count=0
        #print(len(list1[i]),len(list2[i]))
        for j in range(len(list1[i])):
            
            #print(j)
            d2=(list1[i][j][0]-(list2[i][j][0]+list1[i][j][0])*0.5)**2+(list1[i][j][1]-(list2[i][j][1]+list1[i][j][1])*0.5)**2+(list1[i][j][2]-(list2[i][j][2]+list1[i][j][2])*0.5)**2
            #d2=(list1[i][j][0]-list2[i][j][0])**2+(list1[i][j][1]-list2[i][j][1])**2+(list1[i][j][2]-list2[i][j][2])**2
            #print(list1[i][j][0],'\t',list2[i][j][0])
            #print(j,d2)
            count=count+1
            temp=temp+d2
        temp=temp/count
        temp=temp**0.5
        Alltemp.append(temp)
    #print(temp)               
    #temp=temp/count
    #temp=temp**0.5
    
    ave=0
    for i in Alltemp:
        ave=ave+i
    ave=ave/len(Alltemp)    
    return ave


if __name__ == "__main__":

    

    (junk, trjfile, largeBoundary, smallBoundary) = sys.argv
    

    gettrj=readXYZ(trjfile)

    getLargeBound=readBoundari(largeBoundary)
    getSmallBound=readBoundari(smallBoundary)

    largeCG=CGCoord(gettrj,getLargeBound)
    smallCG=CGCoord(gettrj,getSmallBound)

    bpLarge=backmapCG6toCG3(largeCG)

    bpLarge2=backmapCG3toCG1(smallCG)

    printCoord(bpLarge)
    #printCoord(bpLarge2)

    #//compare large and samll
    value=RMSD(smallCG,bpLarge)

    value2=RMSD(gettrj,bpLarge2)

    print(value,value2)

    #printCoord(result)
