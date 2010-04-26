
	
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