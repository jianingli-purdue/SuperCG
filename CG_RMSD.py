import numpy as np
import sys
import math

#deviation of number of res. per site

def readBoundari(boundaryFile):
    boundry=np.loadtxt(boundaryFile,delimiter=',',dtype=int)
    return boundry


def calboundari(boundary):
    blist=[]
    for i in range(len(boundary)):
        if(i==0):
            n=boundary[i]+1
            blist.append(n)
        if(i!=0):
            n=boundary[i]-boundary[i-1]
            blist.append(n)
    return blist


def devCal(list,CGsite):
    sum=0
    for l in list:
        sum=sum+(float(l)-CGsite)*(float(l)-CGsite)
        #print(sum)
    #print(len(list))
    sum=sum/len(list)
    #print(sum)
    sum=math.sqrt(sum)
    return sum


if __name__ == "__main__":

    (junk,boundaryfile,CGsite)=sys.argv
    
    boun=readBoundari(boundaryfile)

    siteList=calboundari(boun)
    #print(siteList)
    CGsite=float(CGsite)
    
    d=devCal(siteList,CGsite)
    print(d)
    


    
    
