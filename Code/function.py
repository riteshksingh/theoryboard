# This is a module of useful functions. See descriptions below.
#==================================================================================

from sympy import zeros, sqrt

#(i)
def makedict(key,value,name):                   	
	if not key in name:
		name[key] = [value]
	else:
		if not value in name[key]:
			name[key].append(value)


#(ii)
def t_add(any_tup1,any_tup2):                 		   
	add=[]
	for i in range(len(list(any_tup1))):
		add.append(any_tup1[i]+any_tup2[i]);
	return tuple(add);


#(iii)
def t_mult(k,any_tup):                        		
	mult = []
	for i in range(len(list(any_tup))):
		mult.append(k*any_tup[i]);
	return tuple(mult)	


#(iv)
def t_subs(any_tup1,any_tup2):                 		   
	subs=[]
	for i in range(len(list(any_tup1))):
		subs.append(any_tup1[i]-any_tup2[i]);
	return tuple(subs);		


#(v)
def makelistoffixedlength(size):					
	listofobjects = list()
	for i in range(size):
		listofobjects.append(list())
	return listofobjects


#(vi)
def ntuple(k,pos):									
	nt = []
	for i in range(k):
		if( i == pos and pos != -1):
			nt.append(1)
		else:
			nt.append(0)
	return nt


#(vii)	spM : scalar product matrix
def gram_schmidt_rotation(spM):	
	dimension = spM.shape[0]
	rotation_matrix = zeros(dimension)
	for i in range(dimension):
		rotation_matrix[i,i] = 1
		for j in range(i):
			scalar_product = sum([rotation_matrix[j,k]*spM[i,k]\
for k in range(j+1)])
			for k in range(j+1):
				rotation_matrix[i,k] += -rotation_matrix[j,k]*scalar_product
		norm_factor = sum([rotation_matrix[i,j]*rotation_matrix[i,k]*\
spM[j,k] for j in range(i+1) for k in range(i+1)])
		for j in range(i+1):
			rotation_matrix[i,j] = rotation_matrix[i,j]/sqrt(norm_factor)
	# resulting matrix is a lower triangular matrix
	return rotation_matrix


#(viii)
def list_product(l1,l2):							
	return sum([a*b for a,b in zip(l1,l2)])

#================================================================================= 			
#                                   Descriptions          
#---------------------------------------------------------------------------------               
#i) makedict is a general dictionary constructor function. Takes 3 arguments in the 
#order - a key of a dict, an 'entry' in the value list for that key and name of the 
#dictionary resp. 
#		makedict first checks whether the key argument exists in the name. If not, 
#then it creates that key with the value list containing the 'entry' . On the other
# hand, if key is in name, then it in turn checks whether the 'entry' argument is in
# the value list for that key. If yes, then it does nothing (we do not want 
# repetitions in value list); if not,then it appends the value list for that key
# with the 'entry' as an item.
#---------------------------------------------------------------------------------
#ii) t_add & iii) t_mult are functions for component-wise addition (for tuples) 
#and scalar multiplication (As we do for elements in vector spaces). iv) t_sub, in a
#similar way, is a function to subtract one tuple from another.
#---------------------------------------------------------------------------------
#v) makelistoffixedlength is self explanatory. Once, a size is chosen, it returns
#a list of objects of that length.
#---------------------------------------------------------------------------------
#vi) ntuple asks for two arguments: size of tuple (say k) and the position (say pos)
#for which we want a non-zero value (value is specifically 1). We can create a zero
#tuple by simply assigning pos==-1 (wlog). CAUTION: misonomer. actually creates list.
#Why list? Because tuple does not allow item assignment
#---------------------------------------------------------------------------------
#vii) Given a matrix of scalar products of linearly independent vectors, returns 
# a rotation matrix that can convert this set of vectors into an orthonormal set.
#---------------------------------------------------------------------------------
#viii) Returns the dot product of two lists/tuples

