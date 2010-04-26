
# Routines for symmetry-reduced transition matrix

from configurations import *


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
			
	