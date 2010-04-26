
# Routines on configurations

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