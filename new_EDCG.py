import sys
import numpy as np


#change the python2 version into python3 version
#------------------
def calc_cov_mat(reffile, framefile, begin_atom, end_atom, report=False):

    # read in reference coordinates, usually from xtal structure
    ref_coord = []  
    ref_fh = open(reffile)
    for line in ref_fh.readlines():
        temp = line.split()
        if len(temp) == 0:
            continue
        else:
            ref_coord.append( [ float(temp[0]), float(temp[1]), float(temp[2]) ])
    ref_fh.close()

    # find out the number of Ca (C alpha) atoms
    Natom = len(ref_coord) # number of Ca
    if Natom == 0:
        print("Error: Zero atom in the reference structure")
        sys.exit(1)
    elif end_atom >= Natom:
        end_atom = Natom-1 # update the last Ca number

    actual_Natom =  end_atom - begin_atom + 1  # end_atom included (need to +1)
    if actual_Natom == 0:
        print("Error: Actual zero atom to read from the reference structure")
        print("       Please check the begin and end atoms", begin_atom, end_atom)
        sys.exit(1)
    else:
        if actual_Natom < Natom:
            ref_coord = ref_coord[begin_atom: end_atom+1]  # update the the ref coordiates
        ref_mat = np.mat( center(ref_coord, actual_Natom) )

#report data
    if report:
        print("\nNatom =",Natom,) 
        print(", Actual Natom =", actual_Natom, ", (begin & end):", begin_atom, end_atom,"\n")

    # read in all frame coordinates & align each frame to ref
    Nframe = 0  # number of frames
    atom_count = 0  
    frame_coord = []  # coordinates for one frame
    frame_list = []   # a list of frames
    
    frame_fh = open(framefile)
    for line in frame_fh.readlines():
   
        if line.strip():  # remove the whitespace at the beginning & end; if not empty
            temp = line.split()  # split the line, and put the words into a list call "temp"
            if atom_count >= begin_atom and atom_count <= end_atom:  # begin_atom and end_atom are included in the selection
                frame_coord.append( [ float(temp[0]), float(temp[1]), float(temp[2]) ])
            atom_count = atom_count + 1
        else:
            if atom_count == Natom:
                frame_mat = np.mat( center(frame_coord, actual_Natom) ) # convert the list to a matrix

                mat_R = rotate(frame_mat[:], ref_mat[:])  # calculate rotation matrix (note: copy the list with [:])
                f2 = (mat_R*frame_mat.T).T
                frame_list.append(np.ravel(f2))
		#for each frame, rotate and saved into frame list as an array
                Nframe = Nframe + 1

                frame_coord = []  # reset for a new frame
                atom_count = 0
                
    frame_fh.close()

    if report: print("Total frame #: %d\n"%Nframe)

    # construct covariance matrix and solve for eigenvalues and eigenvectors
    newcov = np.cov(np.transpose(np.array(frame_list)), bias=1)   # small sample, use N instead of N-1 for data normalization (bias=1)


    #calculate the eigenvector and eigenvalue
    eig_val_cov, eig_vec_cov = np.linalg.eigh(newcov)
	
    #return the index that sort the array in descending order
    idx = eig_val_cov.argsort()[::-1]   # change eigenvalues from ascending order to descending order
    eig_val_cov = eig_val_cov [idx]
    eig_vec_cov = eig_vec_cov [:, idx]
    
    #from here, Use PCA to reduce the dimension

    #essential modes??
    # find out how many essential modes to include to represent (cutoff)% of dynamics
    cutoff = 0.99
    eig_val_sum = np.sum(eig_val_cov)
    
    for N_ed in range( eig_val_cov.shape[0] ): 
        ratio = np.sum(eig_val_cov[0:N_ed+1])/eig_val_sum
        if ratio > cutoff: break
    N_ed = N_ed + 1
    if report: print("%d PCA modes selected for %d percent of the dynamics\n"%(N_ed, cutoff*100))

    # truncate the eigenvalues and eigenvectors
    vec_ess = eig_vec_cov[:, 0:N_ed]
    val_mat = np.diag(eig_val_cov[0:N_ed]) # val_ess is a vector, needs to change to a n_ed*n_ed diagonal matrix

    cov_ed = np.dot(vec_ess, val_mat) # covariance matrix in essential space
    cov_ed = np.dot(cov_ed, vec_ess.T)

    # reduce cov_ed from 3N*3N to N*N
    red_cov_ed = np.zeros((actual_Natom, actual_Natom))
    for i in range(actual_Natom):
        for j in range(actual_Natom):
            for k in range(3):
                red_cov_ed[i,j] = red_cov_ed[i,j] + cov_ed[3*i+k, 3*j+k]

    return red_cov_ed

#------------------

def MCM_search_min_residual(cov_ed, begin_atom, NCG, niter=100, ba_list=[]):

    import random

    Natom = cov_ed.shape[0]

    Natom_per_CG = 1.0*Natom/NCG  # integer

    if len(ba_list)==0:
        for iCG in range(NCG-1):
            ba_list.append(int( iCG*Natom_per_CG+Natom_per_CG  ))
        ba_list.append(Natom-1)  # the last atom number is always at the end of the list

    #print "debug", ba_list
    resid = calc_residual(cov_ed, ba_list, NCG)
    if NCG == 1: return ba_list, resid
    #print "ba_list", ba_list, ", resid", resid

    best_ba_list = ba_list[:]
    best_resid = resid
    resid_0 = resid
    for one_iter in range(niter):

        negative_beta = -0.05*np.sqrt(one_iter) # a function of one_iter, simulated annealing

        iba_to_change = random.randrange(0, NCG-1)  # the last boundary atom (# NCG) is always the last atom

        if iba_to_change > 0:
            new_boundary = random.randrange(ba_list[iba_to_change-1]+1, ba_list[iba_to_change+1])
        else:
            new_boundary = random.randrange(0, ba_list[iba_to_change+1])
 
        old_ba_list = ba_list[:]  # save ba list first
        ba_list[ iba_to_change ] = new_boundary

        resid = calc_residual(cov_ed, ba_list, NCG)

        if resid < best_resid:
            best_ba_list = ba_list[:]
            best_resid = resid
            #print "%10d  "%one_iter, "ba_list", best_ba_list, ", resid", best_resid

        diff_resid = resid - resid_0
        move_prob = np.exp(diff_resid*negative_beta)
        prob_criteria = random.random()

        if move_prob > prob_criteria:  # accept
            resid_0 = resid
            if False: print(one_iter, negative_beta, diff_resid, move_prob, prob_criteria, ba_list)
        else:  # reject
            ba_list = old_ba_list[:]

#rearrange the offset of atom indexs
    for iba in range(NCG):   # correct boundary atom numbering
        best_ba_list[iba] = best_ba_list[iba] + begin_atom

    return best_ba_list, best_resid

#------------------

def MCM_search_incl_fixba(cov_ed, begin_atom, NCG, niter=100, ba_list=[], fix_ba=[]):

    import random

    Natom = cov_ed.shape[0]
    
    #number of atom per coarse graiend bead
    Natom_per_CG = 1.0*Natom/NCG  # integer
    
    #NCG: number of coarse grained beads
    if len(ba_list)==0:
        for iCG in range(NCG-1):
            #append the last index of atom in each coarse grained bead into the list, for NCG-1
            ba_list.append(int( iCG*Natom_per_CG+Natom_per_CG  ))
	#the same as before, for last CG bead
        ba_list.append(Natom-1)  # the last atom number is always at the end of the list
	

#this section should go through the second or third index
#if the second not in previous ba_list, give it a random index point ba_val list 

    if len(fix_ba)>0:
        for ba_val in fix_ba:
            if ba_val not in ba_list:
                iba_to_change = random.randrange(0, NCG-1)
                ba_list[ iba_to_change ] = ba_val
        ba_list.sort()

#calculate the residual of at current step 
    resid = calc_residual(cov_ed, ba_list, NCG)

#assgin the best value of residule and boudnay arry for later calculation 
    best_ba_list = ba_list[:]
    best_resid = resid
    resid_0 = resid

#inerative until find best solution

#niter: number of iteration
    for one_iter in range(niter):
	#should be some interation constant
        negative_beta = -0.05*np.sqrt(one_iter) # a function of one_iter, simulated annealing

        old_ba_list = ba_list[:]  # save ba list first
	
	#randomly find a iba_to_replace not in fix_ba, if it is in fix_ba, do it again
        iba_to_replace = random.randrange(0, NCG-1)  # the last boundary atom (# NCG) is always the last atom
       	while ba_list[iba_to_replace] in fix_ba:
            iba_to_replace = random.randrange(0, NCG-1) 
	
	#the same, randomly find a new boundary, it this new_boundary is in ba_list, do it again
        new_boundary = random.randrange(0, Natom-1)
        while new_boundary in ba_list:
            new_boundary = random.randrange(0, Natom-1)
   	
	#get a new boudary and replace it with a old one, and sort it 
        ba_list[ iba_to_replace ] = new_boundary
        ba_list.sort()
 	
	#calculate the residual again based on the boundaries 
        resid = calc_residual(cov_ed, ba_list, NCG)
	

	#calculate thebest residual for record, but does not move to new solution
	#the new solution based on current 
        if resid < best_resid:
            best_ba_list = ba_list[:]
            best_resid = resid
    #        print "%10d  "%one_iter, "ba_list", best_ba_list, ", resid", best_resid

	#calculate the difference between previous residual and the current residual
	#and use the difference to calculate the probability function
	#get a random number for criteria, which is similiar to the Monte Carlo
        diff_resid = resid - resid_0
        move_prob = np.exp(diff_resid*negative_beta)
        prob_criteria = random.random()	
	
	
        if move_prob > prob_criteria:  # accept
            resid_0 = resid
            if False: print(one_iter, negative_beta, diff_resid, move_prob, prob_criteria, ba_list)
        else:  # reject
            ba_list = old_ba_list[:]


    for iba in range(NCG):   # correct boundary atom numbering
        best_ba_list[iba] = best_ba_list[iba] + begin_atom

    return best_ba_list, best_resid

#--------------------------------------------------------

#calculate overall residual 
def calc_residual(cov_ed, ba_list, NCG):
    total_resid = 0.0
    
    begin_atom = 0
    for iCG in range(len(ba_list)):
        end_atom = ba_list[iCG]
        total_resid = total_resid + residual_from_oneCG(cov_ed, begin_atom, end_atom)
        begin_atom = end_atom + 1

    total_resid = total_resid/(3*NCG)

    return total_resid

#------------------

#calculate the residual from one coarse grain bead
def residual_from_oneCG(cov_ed, begin_atom, end_atom):

    #print("debug", begin_atom, end_atom)
    #print(cov_ed.size)

    ibegin = begin_atom
    iend = end_atom+1

    partial_residual = 0.0

    for i in range(ibegin, iend):
        for j in range(i, iend):
            partial_residual = partial_residual + cov_ed[i,i]-2*cov_ed[i,j]+cov_ed[j,j]

    return partial_residual

#------------------
#this function is not used
def rigid_transform_3D(A, B, N):
    
    '''
    Input: expects Nx3 matrix of points
    Returns R,t
    R = 3x3 rotation matrix
    t = 3x1 column vector
    fit A to B (ref)
    '''    
    assert len(A) == len(B)

    centroid_A = np.mean(A, axis=0)
    centroid_B = np.mean(B, axis=0)

    #print(centroid_A, centroid_B)

    # centre the points
    AA = A - np.tile(centroid_A, (N, 1))
    BB = B - np.tile(centroid_B, (N, 1))

    # dot is matrix multiplication for array
    H = np.dot(np.transpose(AA), BB)

    U, S, Vt = np.linalg.svd(H)

    R = np.dot(np.transpose(Vt), np.transpose(U))

    # special reflection case
    if np.linalg.det(R) < 0:
       #print("Reflection detected")
       Vt[2,:] *= -1
       R = np.dot(np.transpose(Vt), np.transpose(U))

    t = -1.0*np.dot(R, np.transpose(centroid_A)) + np.transpose(centroid_B)
  
    return R, t

#-------------------------------------------
#center of mass
def center(A, N):

#calculate the center
    centroid_A = np.mean(A, axis=0)
    if False: print("xcm:", centroid_A)
    #repeat the elements
    AA = A - np.tile(centroid_A, (N, 1))
    #substract center
    return AA

#--------------------------------------------

#rotate the transfomation
def rotate(A, B):
    
    '''
    Input: expects Nx3 matrix of points; Returns: R = 3x3 rotation matrix
    fit A to B (ref)
    '''    
    assert len(A) == len(B)

    # dot is matrix multiplication for array
    H = np.dot(np.transpose(A), B)

    U, S, Vt = np.linalg.svd(H)

    R = np.dot(np.transpose(Vt), np.transpose(U))

    # special reflection case
    if np.linalg.det(R) < 0:
       #print("Reflection detected")
       Vt[2,:] *= -1
       R = np.dot(np.transpose(Vt), np.transpose(U))

    return R

##################################################

def new_EDCG(trjfile, reffile, NCG, begin_atom=0, end_atom=10000, niter=100  ):

    import time

    start_time = time.time()

    # construct covariance matrix
    
    (cov_ed) = calc_cov_mat(reffile, trjfile, begin_atom, end_atom)

    # total sampling is (Natom-1)!/[(NCG-1)!(Natom-NCG)!], use Monte Carlo Metropolis smapling
    (balist, chi2) = MCM_search_min_residual(cov_ed, begin_atom, NCG, niter)
    
    print("\nbest:", balist, chi2)

    end_time = time.time()

    print("\nNeed %20.2f s"%(end_time-start_time))

    return True

#####################

if __name__ == "__main__":
    
    ## NCG = 4
    ## reffile = "/Users/jianingli/Documents/Research/CoarseGraining_TMPs_2011-2015/transm_proteins/EDCG/2LOP/new_2lop_reffile"
    ## trjfile = "/Users/jianingli/Documents/Research/CoarseGraining_TMPs_2011-2015/transm_proteins/EDCG/2LOP/new_2lop_trjfile"

    print("\n#############################################\n")

    (junk, trjfile, reffile, NCG_str) = sys.argv
    NCG = int(NCG_str)
    #if True:  print("Number of CG:%8d"%NCG)
    
    new_EDCG(trjfile, reffile, NCG, niter=20000*NCG)

    print("\n#############################################\n")

