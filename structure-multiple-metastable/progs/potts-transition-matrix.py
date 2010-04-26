"""Calculate the transition matrix occurring in the master equation for a small Potts system"""

from numpy import *
#from pylab import *

from scipy import sparse
from scipy.sparse.linalg.eigen.arpack import *

from sys import argv, exit

import cPickle

#import arpack

import matplotlib
#matplotlib.use('ps')
#matplotlib.rc('text', usetex=True)
import matplotlib.pyplot as plt


def neighbour(pos, k):
    """Return value of kth neighbour of site pos=i,j)"""
    
    new_pos = pos+neighbour_vector[k]
    new_pos = new_pos[0]%Lx, new_pos[1]%Ly
    return lattice[new_pos]


def calc_energy():
    interaction_energy = 0
    field_energy = 0.0
    
    # more efficient to use ndenumerate?
    # or ndindex
    
    for pos in pos_list:
        value = lattice[pos]    
        
        for j in range(z):
            neighbour_value = neighbour(pos, j)
            
            interaction_energy -= interactions[value, neighbour_value]
            
        field_energy -= field[value]
            
            
    return 0.5*interaction_energy + field_energy
        


def calc_magnetisation():
    
    magnetisation = zeros(q)
    
    for pos in pos_list:
        value = lattice[pos]
        magnetisation[value] += 1
        
    return magnetisation
    
    
def calc_transition_prob(pos, new_value):
    """Calculate transition probability given current configuration in lattice and the new value of the spin at position pos"""
    
    global sum_off_diag
    
    
    delta_H = 0.0
    
    old_value = lattice[pos]
    
    for j in range(z):
        neighbour_value = neighbour(pos, j)
        
        delta_H += interactions[old_value, neighbour_value] - interactions[new_value, neighbour_value]
        # Can make this more efficient with a table?
        
    delta_H += field[old_value] - field[new_value]
    
    if delta_H <= 0.0:
        transition_prob = 1.0
    else:
        transition_prob = exp(-beta * delta_H)
    
    transition_prob /= N * (q-1)
    
    sum_off_diag += transition_prob
    
    if symmetric:
        transition_prob *= exp(0.5*beta*delta_H)
        
    return transition_prob
    
    
def setup_configuration(config_num):
    """Setup lattice with configuration with number config_num"""
    
    mask = max_power
    
    for pos in pos_list:
        current_value = config_num / mask
        lattice[pos] = current_value
        
        config_num -= current_value*mask
        mask /= q
        
        
def get_configuration_num(configuration):
    
    mask = max_power
    config_num = 0
    
    for pos in pos_list:
        current_value = configuration[pos]
        config_num += current_value*mask
        mask /= q
        
    return config_num
    
    
def find_transitions(config_num):
    """Find transitions starting from config_num"""
    
    global sum_off_diag, magnetisation_list
    
    sum_off_diag = 0.0
    
    setup_configuration(config_num)   # fill lattice with starting configuration
    
    boltzmann.append(exp(-beta*calc_energy()))
    
    magnetisation_list.append(calc_magnetisation())
    
    current_config_num = config_num
    
    mask = max_power
    
    for spin in pos_list:
        current_value = current_config_num / mask	# value of current spin being changed
        current_config_num -= (current_value*mask)
        
        for i in range(1,q):
            new_value = (current_value + i) % q
            #if new_value >= q:
                #new_value -= q
                
            transition_prob = calc_transition_prob(spin, new_value)
            
            new_config = config_num + (new_value - current_value)*mask
            
            #print config_num, new_config, transition_prob
            
            transition_matrix[new_config, config_num] = transition_prob  # if symmetric could reduce calculations in half by exploiting symmetry
            #transition_matrix[config_num, new_config] = transition_prob  # if symmetric could reduce calculations in half by exploiting symmetry
            
        mask /= q
        
    transition_matrix[config_num, config_num] = 1.0 - sum_off_diag
    
        
        
def calc_transition_matrix():
    
    #print "Finding transition matrix..."
    
    for i in xrange(total_num_configs):
        if i%1000 == 0:
            print "i = ", i
        find_transitions(i)
        
    print "Done"
    
def inner_prod(v, w):
    if symmetric:
        return dot(v, w)
    else:
        return sum(v*w / P[0])
    
    
def ensure_positive(v):
    max_pos = argmax(abs(v))
    if v[max_pos] > 0:
        return v
    return -v
    
def calc_C(evecs):	
    
    # Gram--Schmidt orthogonalization on eigenvectors
    # If transition matrix is symmetric, then the eigenvectors are orthogonal w.r.t. the usual norm
    
    global P
    
    P = []
    C = zeros(shape(evecs))
    
    
    if symmetric:
        evecs[0] /= sqrt(inner_prod(evecs[0], evecs[0]))
        P.append(ensure_positive(evecs[0]*evecs[0]))
        
    else:
        evecs[0] /= sum(evecs[0])   # Normalise
        
        P.append(ensure_positive(evecs[0]))
    
    num_vecs = len(evecs)
    
    
    # Gram--Schmidt with inner product adecuado
    
    for i in range(1, num_vecs):
        
        #evecs[i] = ensure_positive(evecs[i])
        
        for j in range(0, i):
            evecs[i] -= inner_prod(evecs[i], evecs[j])*evecs[j]
                
        evecs[i] /= sqrt(inner_prod(evecs[i], evecs[i]))
        evecs[i] = ensure_positive(evecs[i])

    
    if symmetric:
        for i in range(1, num_vecs):
            C[i] = ensure_positive(evecs[i] / evecs[0])      # C[0] is still equal to v[0]
                
            #P.append(ensure_positive(evecs[i]*evecs[0]))
            P.append(ensure_positive(C[i]*P[0]))
        
    else: 
        for i in range(1,num_vecs):
            P.append(evecs[i])
            C[i] = ensure_positive(P[i] / P[0])
                
                
    P = array(P)
        
    return C, P

def read_param_file(q, filename):
    
    global field, interactions	
    
    param_file = open(filename, "r")
    
    for line in param_file:
        if line[0]=="#": continue
        
        if len(line[0])==0: continue   # empty line
        
        data = map(float, line.split())
        
        if len(data)==1:
            if q != int(data[0]):
                print "Mismatch in q!"
                exit()
        
        if len(data)==2:    # only two entries, so must be a field
            field[int(data[0])] = data[1]
        
        elif len(data)==3:
            interactions[int(data[0]), int(data[1])] = data[2]
            interactions[int(data[1]), int(data[0])] = data[2]   # symmetric interaction matrix
                
    param_file.close()


def output_C(i):
    outfile = open("C%d.dat" % i, "w")
    
    if reduced:
        for j, mag in enumerate(magnetisation_list):
            for k in mag:
                outfile.write("%g\t" % k)
            outfile.write("%g\n" % C[i, representative_list[j]] )
            #outfile.write("%g\t%g\t%g\t%g\n" %( mag[0], mag[1], mag[2], C[i, representative_list[j]]) )
    else:
        for mag, j in zip(magnetisation_list, C[i]):
            for k in mag:
                outfile.write("%g\t" % k)
            outfile.write("%g\n" % j)
                
            #outfile.write("%g\t%g\t%g\t%g\n" %( j[0], j[1], j[2], k) )


def output_header(outfile):
    outfile.write("# q: %d\n" % q)
    outfile.write("# size: %d x %d\n" % (Lx, Ly))
    
    outfile.write("#\n# beta: %g\n" % beta)
    
    outfile.write("#\n# fields[i]: \n")
    for i in range(q):
        outfile.write("# %d\t%g\n" % (i, field[i]))
        
    outfile.write("#\n# interactions matrix: \n")
    for i in range(q):
        outfile.write("# ");
        for j in range(q):
            outfile.write("%g\t" %  interactions[i,j])
        outfile.write("\n")
        
    outfile.write("#\n# data: \n#\n")
    

def output_file(filename, data):
    '''Output to filename the array data'''
    
    outfile = open(filename, "w")
    print "# Writing to file ", filename
    
    output_header(outfile)

    
    if reduced:
        for j, mag in enumerate(magnetisation_list):
            for k in mag:
                outfile.write("%g\t" % k)
            outfile.write("%g\n" % data[representative_list[j] ] )
            #outfile.write("%g\t%g\t%g\t%g\n" %( mag[0], mag[1], mag[2], C[i, representative_list[j]]) )
    else:
        for mag, j in zip(magnetisation_list, data):
            for k in mag:
                outfile.write("%g\t" % k)
            outfile.write("%g\n" % j)
    
    #for i in data:
        #outfile.write("%g\n" % i)
    
    
    
    outfile.close()
    

def output(base_filename, num):
    '''Output calculated C and P
    Uses base_filename to construct filename, outputs _Ci and _Pi up to num
    Each file contains the information on field, interactions etc.'''
    
    for i in range(num+1):
        output_file(base_filename+"_P"+str(i)+".dat", P[i])
        
        if i>0:
            output_file(base_filename+"_C"+str(i)+".dat", C[i])
            
        
    
    outfile = open(base_filename+"_evals.dat", "w")
    output_header(outfile)
    
    outfile.write("%s" % str(evals) )
    #for i in evals:
        #outfile.write("%g\t" % i)
    
    #output_file(base_filename+"_evals.dat", evals)
    


# Symmetry operations on configs:


def translate(inc):
    """Translation by inc = (x_inc, y_inc)"""
    global translated_config
    
    #print "Translation by", inc
    
    for pos in pos_list:
        new_pos = (pos[0]+inc[0]) % Lx, (pos[1]+inc[1]) % Ly
        #new_pos = new_pos[0]%Lx, new_pos[1]%Ly
        
        translated_config[new_pos] = lattice[pos]
    #print "translated_config = ", translated_config
        
def reflect_horizontal(config):
    for pos in pos_list:
        new_x = pos[0]
        new_y = Ly - 1 - pos[1]
        
        reflected_config[new_x, new_y] = config[pos]
        
def reflect_vertical(config):
    for pos in pos_list:
        new_x = Lx - 1 - pos[0]
        new_y = pos[1]
        
        reflected_config[new_x, new_y] = config[pos]
        

def rotate():
    # Will rotate several times, so copy the current state:
    
    global rotated_config, new_config
    
    new_config = rotated_config.copy()
    
    
    for pos in pos_list:
        new_x = Lx - 1 - pos[1]
        new_y = pos[0]
        
        rotated_config[new_x, new_y] = new_config[pos]


def do_orbit(config_num):
    """Find the lowest-numbered representative for the orbit of config_num"""
    
    # To do so, apply all symmetry operations allowed by the lattice, and find resulting configurations
    # The one with the minimum number is the representative of all the others
    
    # Assume that everything generated in increasing order, so if find a member of the orbit which is smaller, or more generally which has already
    # been assigned a representative, then can immediately assign that representative to all the things found in the current orbit?
    
    
    # Each configuration has a representative of its orbit
    # The representative is the configuration which is lexicographically least in the group orbit under the symmetry group
    
    # representative_list  is a list of all the representatives, once each
    
    # representatives gives the *position* in representative_list of the representative of each configuration
    
    # For a given configuration, we thus start generating its complete group orbit
    # As soon as we find a configuration which has already been assigned a representative, we know what the representative is for that whole orbit, namely this representative just found?
    
    # [Shouldn't we also check that the new configuration hasn't already been found in the current orbit?]
    
    found = False
    
    global rotated_config, translated_config, reflected_config, representative_list, representatives
    
    if representatives[config_num] >= 0:   # already processed
        return 
    
    
    
    
    orbit = set([])   # orbit = set of images under the symmetry elements
    
    setup_configuration(config_num)
    
    #current_representative = config_num
    
    # Loop over possible increments in the lattice:
    for inc in pos_list:
        translate(inc)   # sets up translated config, starting from non-translated, which we know hasn't been found yet, so orbit has at least one element
        
        new_config_num = get_configuration_num(translated_config)
        #print new_config_num
        #if representatives[new_config_num] >= 0:
            #break
        
    
        if Lx==Ly:   # square lattice has more symmetry
            # for each translated config, do rotations and a single reflection of each rotation
            
            #print "translated_config:"
            #print translated_config
            
            rotated_config = translated_config.copy()
            
            for rotation in range(4):
                # start from unrotated (identity)
                
                #print "rotated_config = \n", rotated_config
                
                new_config_num = get_configuration_num(rotated_config)
                #print new_config_num
                #if representatives[new_config_num] >= 0:
                    #break
                orbit.add(new_config_num)
                
                reflect_horizontal(rotated_config)
                new_config_num = get_configuration_num(reflected_config)
                #print new_config_num
                #if representatives[new_config_num] >= 0:
                    #break
                orbit.add(new_config_num)
                
                rotate()
                
        else:   # rectangular, non-square.  Do translated, HT, VT, HVT, where H and V are horiz and vert
            
            new_config_num = get_configuration_num(translated_config)
            #if representatives[new_config_num] >= 0:
                #break
            orbit.add(new_config_num)
            
            reflect_horizontal(translated_config)
            new_config_num = get_configuration_num(reflected_config)
            #if representatives[new_config_num] >= 0:
                #break
            orbit.add(new_config_num)
            
            reflect_vertical(translated_config)
            new_config_num = get_configuration_num(reflected_config)
            #if representatives[new_config_num] >= 0:
                #break
            orbit.add(new_config_num)
            
            new_config = reflected_config.copy()
            
            reflect_horizontal(new_config)
            new_config_num = get_configuration_num(reflected_config)
            #if representatives[new_config_num] >= 0:
                #break
            orbit.add(new_config_num)
        
        #if representatives[new_config_num] >= 0:
                #break
    
    
    #if representatives[new_config_num] >= 0:
        #print "HOLA"
        # THIS CASE NEVER OCCURS!  SINCE IF IT DID, IT WOULD MEAN THAT THE CURRENT CONFIG HAD ALREADY BEEN FOUND IN A PREVIOUS ORBIT
        # already found a representative
        #representative_position = representatives[new_config_num]
        #class_size_list[representative_position] += len(orbit)
    #else:
    representative = min(orbit)    # Take the minimum configuration found
    # ACTUALLY WE KNOW THIS MUST BE config_num
    
    representative_position = len(representative_list)   # new element
    representative_list.append(representative)
    
    class_size_list.append(len(orbit))
        # This should be the original one
    
    for i in orbit:
        representatives[i] = representative_position
    
    
    
    #print "Config_num = ", config_num
    #print "Representative = ", representative
    #print "Representative position = ", representative_position
    #print "Orbit = ", orbit
    
    
def find_representatives():
    for i in range(total_num_configs):
        if i%100==0:
            print "Finding orbit of config", i
        do_orbit(i)
        
    #print "Representative_list:"
    #print representative_list
    
    print "Number of representatives:", len(representative_list)
    
    print "Ratio of total number of configs to number of representatives: ", float(total_num_configs) / len(representative_list)
    
    print "Max theoretical: ",
    if Lx==Ly: print 8*N
    else: print 4*N
    
    
# The following versions do the symmetry-reduced transition matrix

    
def calc_reduced_transition_prob(pos, new_value):
    """Calculate transition probability given current configuration in lattice and the new value of the spin at position pos
    Goes in the OPPOSITE DIRECTION from calc_transition_prob,
    i.e. calculates the probability of going *from* new_value *to* current value, simply by taking -Delta_H
    """
    
    
    global sum_off_diag
    
    
    delta_H = 0.0
    
    old_value = lattice[pos]
    
    for j in range(z):
        neighbour_value = neighbour(pos, j)
        
        delta_H += interactions[old_value, neighbour_value] - interactions[new_value, neighbour_value]
        # Can make this more efficient with a table?
        
    delta_H += field[old_value] - field[new_value]
    
    # first calculate the contribution the transitions *from* the current config, which will be used to calculate the diagonal part
    if delta_H <= 0.0:
        sum_off_diag_contribution = 1.0
    else:
        sum_off_diag_contribution = exp(-beta * delta_H)
    
    sum_off_diag_contribution /= N * (q-1)
    
    
    # now calculate the transition *to* the current config
    
    delta_H = -delta_H    # transition in *opposite* direction
    
    if delta_H <= 0.0:
        transition_prob = 1.0
    else:
        transition_prob = exp(-beta * delta_H)
    
    transition_prob /= N * (q-1)
    
    #sum_off_diag_contribution = transition_prob
    
    if symmetric:
        transition_prob *= exp(0.5*beta*delta_H)
        
    return transition_prob, sum_off_diag_contribution


def find_reduced_transitions(representative_position):
    """Find transitions TO a given representative
    To do so, start from the given representative, and find all configurations accesible *from* it
    These are the configurations from which transitions are possible!
    Sum over *those* configurations in each equivalence class """
    
    global magnetisation_list
    
    sum_off_diag = 0.0   # sum of off-diagonal elements, i.e. probs of transition to other states
    #new_sum_off_diag = 0.0
    
    config_num = representative_list[representative_position]
    setup_configuration(config_num)   # fill lattice with starting configuration
    
    boltzmann.append(exp(-beta*calc_energy()))
    
    magnetisation_list.append(calc_magnetisation())
    
    current_config_num = config_num
    #current_representative = representatives[current_config_num]
    
    mask = max_power
    
    for spin in pos_list:
        
        current_value = current_config_num / mask	# value of the spin 
        current_config_num -= (current_value*mask)
        
        for i in range(1,q):
            new_value = (current_value + i) % q
            #if new_value >= q:
                #new_value -= q
                
            #transition_prob = class_size_list[representative_position] * calc_transition_prob(spin, new_value)
            transition_prob, sum_off_diag_contribution = calc_reduced_transition_prob(spin, new_value)
            # there are lots of states with the same representative!
            
            new_config_num = config_num + (new_value - current_value)*mask
            
            #print config_num, new_config, transition_prob
            
            new_representative_position = representatives[new_config_num]
            
            #ratio = sqrt(float(class_size_list[new_representative_position]) / class_size_list[representative_position])
            #ratio = float(class_size_list[representative_position]) / class_size_list[new_representative_position]
            # ratio gives ratio of sizes of classes
            
            #ratio = class_size_list[new_representative_position]
            
            #transition_prob *= ratio
            #sum_off_diag_contribution *= ratio
            
            sum_off_diag += sum_off_diag_contribution
            
            reduced_transition_matrix[representative_position, new_representative_position] += transition_prob  
            #if new_representative_position != representative_position:
                #new_sum_off_diag += transition_prob
            
        mask /= q
        
    reduced_transition_matrix[representative_position, representative_position] = 1.0 - sum_off_diag
    #1.0 - class_size_list[representative_position]*sum_off_diag
    # config->config prob is 1.0 - sum_off_diag if symm
    # But could already be other transitions when reduced to representatives?
    # No, I don't think so, since the magnetisation always changes
    
    
    
        
        
def calc_reduced_transition_matrix():
    
    print "Calculating reduced transition matrix..."
    
    for i in range(num_representatives):
        if i%100 == 0:
            print "i = ", i
        find_reduced_transitions(i)
        
    #for i in range(num_representatives):
        #reduced_transition_matrix[i,i] = 1 - (reduced_transition_matrix.getcol(i)).sum()     # works in non-symmetric case?
        
    print "Done"
    
    
def reconstruct_full(reduced_evecs):
    num_evecs = len(reduced_evecs)
    
    full_evecs = zeros( (num_evecs, total_num_configs) )
    
    for j in range(total_num_configs):
        reduced_pos = representatives[j]
        
        #for i in range(num_evecs):
            #full_evecs[i][j] = reduced_evecs[i][reduced_pos]
        full_evecs[:, j] = reduced_evecs[:,reduced_pos]
        
    return full_evecs
            
    
##################
# MAIN


beta = 1.0

# Read parameter file:

#q = 2
#Lx = Ly = 3


try:
    q = int(argv[1])
    param_file = argv[2]
    beta = float(argv[3])
    Lx, Ly = map(int, argv[4:6])
    symmetric = int(argv[6])
    reduced = int(argv[7])
except:
    print "Sintaxis: python potts-transition-matrix.py q parameter_file beta Lx Ly symmetric reduced"
    exit()


field = zeros(q)
interactions = eye(q)  # start with standard Potts model interactions

read_param_file(q, param_file)

N = Lx * Ly
z = 4   # coordination number of lattice 

lattice = zeros((Lx, Ly), dtype='int')
translated_config = zeros((Lx, Ly), dtype='int')
reflected_config = lattice = zeros((Lx, Ly), dtype='int')
rotated_config = lattice = zeros((Lx, Ly), dtype='int')
new_config = lattice = zeros((Lx, Ly), dtype='int')

magnetisation = zeros(q)

magnetisation_list = []

boltzmann = []

sum_off_diag = 0.0

#symmetric = True

neighbour_vector = array([ [-1,0], [0,1], [1,0], [0,-1] ])

max_power = q**(N-1)
total_num_configs = q**N
representatives = -ones(total_num_configs, dtype='int')
representative_list = []
class_size_list = []




m = meshgrid(range(Lx), range(Ly))
pos_list = zip(m[0].flatten(), m[1].flatten())
#pos_list = array(pos_list)



print "Fields:"
print field
print "\nInteractions:"
print interactions
print "\nInverse temp. beta = ", beta
print "\nSystem size: ", Lx, "x", Ly
print "\nTotal number of configurations = ", total_num_configs
if symmetric:
    print "\nUsing symmetric transition matrix"
else:
    print "\nUsing non-symmetric transition matrix"
print

#if not symmetric:
    #transition_matrix = transition_matrix.T

P = []

if reduced:
    
    try:  # check if already found representatives for this case
        
        reps_filename = "reps_q%d_%dx%d.pickle" % (q, Lx, Ly)
        reps_file = open(reps_filename, "r")
        
        print "Already found representatives, so unpickling..."
        (representative_list, representatives) = cPickle.load(reps_file)
        print "Done"
        
    except:  # if not, then find representatives, and save them
        
        print "Haven't found representatives yet"
        
        find_representatives()
        reps_filename = "reps_q%d_%dx%d.pickle" % (q, Lx, Ly)
        reps_file = open(reps_filename, "w")
        
        print "Pickling..."
        cPickle.dump( (representative_list, representatives), reps_file ) 
        print "Done"
        
    num_representatives = len(representative_list)
    
    reduced_transition_matrix = sparse.lil_matrix( (num_representatives, num_representatives) )
    
    calc_reduced_transition_matrix()
    
    class_size_list = array(class_size_list)
    
    boltzmann = reconstruct_full(array([boltzmann]))[0]
    boltzmann /= sum(boltzmann)		# normalise
    
    print "Finding eigen-cosas..."
    
    if num_representatives >= 10:
        evals, evecs = arpack.eigen(reduced_transition_matrix, 10)
    else:
        evals, evecs = arpack.eigen(reduced_transition_matrix, 3)
    evecs = evecs.transpose()    # so eigenvalues are now rows of evecs
    
    evecs = evecs[0:-1]
    evals = evals[0:-1]
    
    evecs = real(evecs)
    evals = real(evals)
    
    
    perm = argsort(-evals)
    
    evals = evals[perm]
    evecs = evecs[perm]
    
    print "Evals: "
    print evals
    print "\nTransition times: ", 1./(1-evals)
    
    print "Reconstructing full evecs..."
    
    full_evecs = reconstruct_full(evecs)
    
    print "Calculating C and P..."
    C, P = calc_C(full_evecs)


else:   # non-reduced (full) case
    
    
    transition_matrix = sparse.lil_matrix( (total_num_configs, total_num_configs) )
    
    print "Calculating transition matrix..."
    calc_transition_matrix()
    
    boltzmann /= sum(boltzmann)
    
    
    
    print "Finding eigen-cosas..."
    #e,v = arpack.eigen(transition_matrix)
    
    if symmetric:
        #evals, evecs = arpack.eigen_symmetric(transition_matrix,10)
        evals, evecs = arpack.eigen(transition_matrix, 10)
    else:
        evals, evecs = arpack.eigen(transition_matrix, 10)
    evecs = evecs.transpose()    # so eigenvalues are now rows of evecs
    
    evecs = evecs[0:-1]
    evals = evals[0:-1]
    
    evecs = real(evecs)
    evals = real(evals)
    
    #d = {}  # dictionary to sort evecs and evals
    
    #for i in range(len(evals)):
        #d[evals[i]] = evecs[i]
        
    ##new_evals = zeros(len(evals))
    ##new_evecs = zeros(len(evecs))
        
    #new_evals = []
    #new_evecs = []
        
    #for i in sorted(d.keys(), reverse=True):
        #new_evals.append(i)
        #new_evecs.append(d[i])
    
    ##z = array(zip(evecs, evals))
    ##sort(z)
    
    #new_evals = array(new_evals)
    #new_evecs = array(new_evecs)
    
    perm = argsort(-evals)
    
    evals = evals[perm]
    evecs = evecs[perm]
    
    
    C, P = calc_C(evecs)
    
    print "Evals: ", evals
    print "\nTransition times: ", 1./(1-evals)

#for i,j in zip(magnetisation_list, C[1]):
    #print i[0], i[1], i[2], j
    
    
def plot_scaled(i):
    t = linspace(0., 1., total_num_configs) 
    sorted_C = sort(C[i])
    
    plt.plot(t, sorted_C/max(sorted_C), 'o-')  