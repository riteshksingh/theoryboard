from Lagrangian import *
SU5 = Group('SU(5)',False,latex_index_name='a')
SU3 = Group('SU(3)',False,latex_index_name='b')
i = Rep('i', SU5, singlet= False, dim=10)
j=Rep('j', SU3, singlet=True)
k=Rep('k', SU5, singlet=False)
phi1 = Scalar('phi1',latex_name='\phi',rep=[i,j])

phi2 = Scalar('phi2',latex_name='\ph2',rep=[i,j])
phi3 = Scalar('phi3',latex_name='\ph3',rep=[i])
c=MixingTensor('C', latex_name='C', rep=[i])
d=MixingTensor('D', latex_name='D', rep=[i])
print(c)
L = Lagrangian()
L.AddParticles( phi1)
L.AddMass(phi1, 'm')
L.AddTerm('l', phi1.dag, phi1.field, phi1.dag, phi1.field, contraction_pattern=[[1,0], [-1,0], [2,0], [-2,0]])
for x in L():
    print(L.contraction_pattern(x))
