from Lagrangian import *
#create the groups
SU5 = Group('SU(5)',False,latex_index_name='a')
SU3 = Group('SU(3)',False,latex_index_name='b')
SU2 = Group('SU(2)',False,latex_index_name='c')
i = Rep('i', SU5, singlet= False, dim=10)
j=Rep('j', SU3, singlet=True)
k=Rep('k', SU2, singlet=False)
#list of non-singlet fields
phi1 = Scalar('phi1',latex_name='\phi',rep=[i, j, k])

phi2 = Scalar('phi2',latex_name='\ph2',rep=[i, j, k])
phi3 = Scalar('phi3',latex_name='\ph3',rep=[i, j, k])

#list of singlet fields
phi4 = Scalar('phi4',latex_name='\phi4',rep=[j])
phi5 = Scalar('phi5',latex_name='\phi5',rep=[j])
phi6 = Scalar('phi6',latex_name='\phi6',rep=[j])
c=MixingTensor('C', latex_name='C', rep=[i])
d=MixingTensor('D', latex_name='D', rep=[i])
print(c)



L = Lagrangian()
L.AddParticles(phi6)
L.AddMass(phi6, "m")
#L.AddParticles( phi1)
#L.AddMass(phi1, 'm')
#L.AddTerm('l', phi1.dag, phi1.field, phi1.dag, phi1.field, contraction_pattern=[[1, 0, 2], [-1, 0, -2], [3, 0, 4], [-3,0, -4]])
for x in L():
    print(x)

print('start replace')
v=Symbol('v', Real=False)
b=Symbol('b')

print(isinstance(phi3, Field))
#how to replace field with set of linear combination of fields
M=L.Replace3(phi6, [['a', 'b' ],[ phi3, phi5]])
#Now we need to implement the replacement with product of fields.
#M=L.Replace3(phi, [[c],[]])
print("output")
for x in M:
    print(x)


