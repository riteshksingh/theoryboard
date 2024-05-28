'''
Symbolic computations in QFT Lagrangians
Language: Python



The programme can be divided into 2 parts. the first one corrosponds to the 
creation of all the components required in the programme and the second one 
deals with the creation and modification of the Lagrangian
'''
from sympy import Symbol, sympify, Matrix, zeros, sqrt, I
from Utilities import *
from LieAlgebra import *
from Rep import *

####################Global attributes###########################################
FieldDict={}
################################################################################
class Group():
    '''
    Represents the Lie Group.
    grp_dict stores details of all the instances created in the code.
    this class also deals with structure constants.
    '''
    grp_dict = {}
    #===========================================================================
    def __init__(self, name,abelian=False, **kwargs):
        #major attributes are name, abelian.
        #dim is optional
        #rest of the attributes can be added as kwargs
        self.name = name
        self.abelian = abelian
        

        #if latex name is not given expilicitly we take 'Ad+name' as the 
        #name of adjjoint representation---------------------------------------------> add indices() directly from here
        if 'latex_index_name' in kwargs.keys():
            self.index_name=kwargs['latex_index_name']
        else:
            self.index_name='Ad_'+self.name
        
        #name is the primary key, as no two group instances take same name
        if self.name not in Group.grp_dict:
            Group.grp_dict.update({self.name: self})
        else:
            raise Exception("Already defined name!")

      
        flag=0
        #---------including Lie Algebra Package---------------------------------
        if 'LA' in kwargs.keys():
            flag=1
            arguments=kwargs['LA']
            print(arguments)
            self.LA=LieAlgebra(arguments[0], arguments[1])
            print(self.LA)
            self.adjoint_rep=Rep(self.index_name, self,
                   LARep=[self.LA, self.LA.adj_rep()])
            self.dim=self.adjoint_rep.dim
            self.Nroots=(self.dim-self.LA.rank)/2
        else:
            self.LA=None
        #-------------adjoint rep--------------------------------
            if 'dim' in kwargs.keys():
                print('inside dim')
                if flag==1:
                    pass
                elif kwargs['dim']!=None:
                    self.adjoint_rep = Rep(self.index_name, 
                                       self, 
                                       dim=kwargs['dim'])   
                                 
                                                                                    #why create it?       
            else:
                print('inside else')
                self.adjoint_rep = Rep( self.index_name, self, dim=None)
                print(self.adjoint_rep)
         #-------------adjoint rep-------------------------------
        if 'Nroots' in kwargs.keys():
            if flag==1:
                pass
            self.Nroots=kwargs['Nroots']
        #create the f matrix
        if self.abelian==False and self.name != 'Lorentz':
            self.f= self.f_matrix()   
            FieldDict.update({self.f.name: self.f})
        else:
            self.f=None
                                                                                    #why create it?       
            
    #=========================================================
    def f_matrix(self, idx=[1]): #this need to be written using the constant Tensor 
        #create the structure constants of the adjoint representation
        
        #explicit name----------------------------------------------------------------->define externally
        #explicit_name='f_x_'+self.adjoint_rep.index_name+'_'+\
                   #self.adjoint_rep.index_name+'_'+self.adjoint_rep.index_name

        # First arrange all the properties in a dict.
        
        #--------------This part deals with the Lie algebra Package-----------------------
#        if self.LA:
#            numerical_values=get_numerical_values_mapping( self.adjoint_rep, True)
 #       else:
  #          numerical_values=None
        #---------------------------------------------end of lie algebra---------------------
        properties={'rep':[self.adjoint_rep, self.adjoint_rep, 
                      self.adjoint_rep],
                     'field_type':'constant', 'latex_name': 'f',
                     'symmetry':-1, 
                     #'numerical_values':numerical_values
                     }

        e_name='f'+explicitname(properties['rep'])
        properties.update({'explicit_name':e_name})
        #Create indices based on idx=[0]/[1](covariant/contravariant)       
        indices=Create_indices(properties['rep'],idx)

        #create the field
        f=Field('f_'+ self.name, indices, **properties)
        
        return(f)
       
    #===========================================================================
    def __repr__(self):
        #printing of the group instance in a terminal
        return (self.name)



################################################################################
class Rep():
    '''
    Represents the representation of the Lie Group.
    rep_dict stores all the Rep instances.
    T_name stores all the T matrices created in the code.
    '''
    rep_dict = {}
    T_name = []

    #========================================================
    def __init__(self, name, group,singlet=False,**kwargs):
        #main attributes includes name, group and singlet.
        #dim is optional. rest of the attributes can be added as kwargs
        self.group = group
        self.name = name
        self.singlet = singlet
        

        #modified while doing the MakeExplicit method in the Lagrangian class
        self.pseudo_singlet = False 
        
        #----latex name of the Rep, uses for printing. default value is name of Rep---
        if 'latex_index_name'in kwargs.keys():
            self.index_name=kwargs['latex_index_name']
        else:
            self.index_name=self.name	

        #----name of the Rep is primary as no 2 Rep instances has same name-----
        if self.name not in Rep.rep_dict:
            Rep.rep_dict.update({self.name: self})
        else:
            print(self.name)
            print(Rep.rep_dict)
            raise Exception ("Already defined name!")
        #---Representation from LiePy--------------------------
        if 'LARep' in kwargs.keys():
            arguments=kwargs['LARep']
            #print(arguments)
            self.LARep=Representation(arguments[0], arguments[1])
            self.dim=self.LARep.weyldim
        else:
            self.LARep=None
            
        if 'dim' in kwargs.keys():
            if 'LARep' in kwargs.keys():
                pass
            elif kwargs['dim']!=None: 
                self.dim=kwargs['dim']
            else:
                self.dim=None
        #if self.singlet == False and self.group.name != 'Lorentz':
        #    self.T = self.T_matrix() 
    #===========================================================================
    def T(self, gauge_field, idx=[1]): #This need to be written using Constant Tensor
        # creates T(a,i,j)
        #follows the same structure as f().arrange properties,createindices
        #create the tensor

        #add the name of Tmatrix to the list T_name--------------------------------------> wrong name needs correction
        if 'T'+gauge_field.g_name not in self.T_name:
            self.T_name.append('T'+gauge_field.g_name)

        print(self.group.name)
        #if name group is not ableian : one adjoint rep + 2 rep indices 
        if self.group.abelian == False:
            properties = {'rep': [self.group.adjoint_rep, self, self]}
            #explicit_name='_x_'+self.group.adjoint_rep.index_name+'_'+\      #------check make explicit
                             #self.index_name+'_'+self.index_name
        #if group is abelian: 2 rep indices
        else:
            properties = {'rep': [self, self]}
            #explicit_name='_x_'+self.index_name+'_'+self.index_name         #-------check makeexplicit
        e_name='T'+self.index_name.upper()+\
                       explicitname(properties['rep'])
        if self.LARep:
            numerical_values=get_numerical_values_mapping( self)
        else:
            numerical_values=None
        #combining the properties of the T matrix as a dictionary
        properties = {**properties,'field_type':'constant',
                     'explicit_name':e_name,
                     'latex_name':'(T_{'+self.index_name+'})',
                     'symmetry':0, 
                     'numerical_values': numerical_values}

        #creating indices 
        indices=Create_indices(properties['rep'],idx)
 
        #creating field
        t_matrix = Field('T'+self.index_name.upper(), indices, **properties)
        
        return (t_matrix)
    #==========================================================================
    def delta(self,indices):
        #create the kronecker delta with the indices passed as arguments
        d=Field('KD',indices,field_type='constant')
        return(d)
    #---------------------------------------------------------------------------
    def __repr__(self):
        #printing of the Rep instance
        return (self.name  )

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def get_numerical_values_mapping(rep, if_f=False):  
    #if I pass only matrices  and name as input, I can take this to utilities
    Tvalues={}
    matrices =rep.LARep.matrices()
    roots=rep.LARep.LA.labelled_roots()
    if if_f==True:
        nroots=rep.group.Nroots
        dim=rep.dim
        conversion=dim_to_cartan_weyl(dim, nroots)
        conversion1=conversion[0]
        conversion2=conversion[1]
    for x in roots.keys():
        y=roots[x]
        if isinstance(y, list):
            for z in y:
                w=matrices[z]
                nonzero_values=nonzero_matrix(w)
                for p in nonzero_values:
                    
                    if if_f==True:
                        Tvalues.update({('c'+str(z+1), conversion1[p[0]], conversion1[p[1]]): p[2]})
                    else:
                        Tvalues.update({('c'+str(z+1),str(p[0]+1), str(p[1]+1)): p[2]})
                    
        else:
            z=matrices[y]
            nonzero_values=nonzero_matrix(z)
            if '+' in str(x):
                v=x.replace('+', 'p')
            if '-' in str(x):
                v=x.replace('-', 'm')
            for p in nonzero_values:
                
                if if_f==True:
                    Tvalues.update({(str(v), conversion1[p[0]], conversion1[p[1]]): p[2]})
                else:

                    Tvalues.update({(str(v), str(p[0]+1), str(p[1]+1)):p[2]})
    count=0
    for x in Tvalues:
        if len(x)-len(set(list(x)))!=0:
            count+=1
    return(Tvalues)

################################################################################
class Index():
    '''
    Represents the index of the Representation
    idx_dict stores all the Index instancs--defined to call any instance outside
    the Index class
    '''
    idx_dict = {}
    #===========================================================================
    def __init__(self, name, Rep):
        #attribute is the name string and the Rep instance
        #name string is converted into a symbol
        self.index=Symbol(name, Real=True)
        #type of the index is name of the rep---This way it is easier to get back 
        # to the representation.
        self.index_type=Rep.name
        
        #all the instances created are added to the idx_dict
        #index is stored in the dictionary with keyword 'index'
        if self.index not in Index.idx_dict:
            Index.idx_dict.update({self.index:self})
    #===========================================================================
    def __str__(self):
        #printing of the Index
        return ('%s' % (str(self.index))  )
    #==========================================================================
    def __repr__(self):
        #printing of the Index instnace
        return (str(self.index))

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def Create_indices(rep_list, idx=[0]):
    '''
    To create indices of the fields from the defined representations
    Take the list of reps and the type of the indices(covariant/ contravariant) 
    as arguments. Type of indices can be mixed as well.In that case,
    len(idx)=len(rep_list)

    '''
    N_list=rep_list
    I_list=idx

    #if we give idx=[1]/[0] then multiply [1]/[0] with the len(rep_list)
    #this all the indices will be either covariant or contravarient
    length = len(N_list)
    if len(I_list) == 1:
        if I_list[0] == 1:
            I_list = [1 for i in range(length)]
        elif I_list[0] == 0:
            I_list = [0 for i in range(length)]
        else:
            raise Exception("The input should [0] or [1]")

    #if any of the rep is singlet, then the index should be taken as covariant
    for x in range(len(rep_list)):
        if rep_list[x].singlet==True:
            I_list[x]=0
   
    #stores the index symbols
    index_list = []
    #stores the index position
    index_pos=[] 

    #keep the count of each rep in the rep_list. prevents the creation of same
    #index twice   
    index_names = {}

    for i in range(length):
        try:
            index_names[N_list[i].name] += 1
        except NameError:
            index_names = {N_list[i].name: 1}
        except KeyError:
            index_names.update({N_list[i].name: 1})
        num = index_names[N_list[i].name]
        #if index is lower add -1 to the index_pos. Uses index name attribute
        #from Rep class as the name of index
        if I_list[i] == 0:                                                           
            index_list.append(Index(N_list[i].index_name + str(num), N_list[i]))
            index_pos.append(-1)

        # if index is upper add +1 to the index pos
        else:
            index_list.append(Index(N_list[i].index_name + str(num), N_list[i]))
            index_pos.append(1)
            
    return (index_list,index_pos)

################################################################################
class Field:

    '''
    Represents Field/ Constant tensors
    '''
    #---------------------------------------------------------------------------
    def __init__(self, name,indices,**kwargs):
        #main attributes are name and indices
        self.name=name
        #indices has format of '2 lists in a list'. first list stores the symbol
        #of the indices, second is the sign of the indices(+1 or -1)
        #couldn't find another way to deal with -ve symbols in sympy
        self.indices=indices

        #latex name of the field
        #default value is name itself
        if 'latex_name' in kwargs.keys():
            self.latex_name=kwargs['latex_name']
        else:
            self.latex_name=self.name
       
        #type of the field - scalar, fermion, constant, gauge
        if 'field_type' in kwargs.keys():
            self.field_type=kwargs['field_type'] 

        #explicit name of the field - for MakeExplicit expansion
        if 'explicit_name' in kwargs.keys():
            self.explicit_name=kwargs['explicit_name']

        #self adjoint is the property only for scalar fields
        if 'self_adjoint' in kwargs.keys():
            self.self_adjoint=kwargs['self_adjoint']
        
        #list of representations of the field
        if 'rep' in kwargs.keys():
            self.rep=[x.name for x in kwargs['rep']]

        #For normal field, base field is itsself. It is for fields after explicit
        #expansion.
        if 'base' in kwargs.keys():
            self.base=kwargs['base']
        else:
            self.base=self

        #if the field is symmetric in expansion we, use +1, antisymmetric -1,
        #no symmetry 0.
        if 'symmetry' in kwargs.keys():
            self.symmetry=kwargs['symmetry']
        else:
            self.symmetry=0

        #-----Charge in each representation----------------->not implemented fully
        if 'charge' in kwargs.keys():
            self.charge=kwargs['charge']

        #---------------numerical values for constant tensors
        if 'numerical_values' in kwargs.keys():
            self.numerical_values=kwargs['numerical_values']
        else:
            self.numerical_values=None
    #===========================================================================
    @property
    def get_indices(self):       
        #gives the indices of the field. Not exactly required here. But 
        #introduced to make it similar to the FieldMul class
        return(self.indices)
    #===========================================================================
    @property
    def get_free_indices(self):
        #Gives only the free indices in a field

        #take element wise product of 2 lists in indices
        all_idx=[self.indices[1][x]*self.indices[0][x].index \
                   for x in range(len(self.indices[0]))]
        free_idx_symbol=[]
        free_idx_sign=[]

        #if a index and its -ve appears in the list-- index is not free
        for x in range(len(all_idx)):
            if -1*all_idx[x] in all_idx:
                continue
            else:
                free_idx_symbol.append(self.indices[0][x])
                free_idx_sign.append(self.indices[1][x])
        return(free_idx_symbol, free_idx_sign)
    
    #===========================================================================
    def check_indices(self):
        #check no two indices repeats in the field--------------------------------------> can move to the init 
        all_idx=[self.indices[1][x]*self.indices[0][x].index \
                   for x in range(len(self.indices[0]))]
        for x in all_idx:
            if all_idx.count(x)==1:
                continue
            else:
                raise Exception('Field contain multiple number of same indices')
        #one more command is required                         
            
    #==========================================================================
    def convert(self):
        #covert all the indices to its opposite sign------------------------------------not needed.can be removed
        self.indices[1]=[-1*x for x in self.indices[1]]
    #==========================================================================
    def __str__(self):
        #printing of the field: 'name(indices)'
        index_string='('
        for x in range(len(self.indices[0])):
            index_string+=str(self.indices[1][x]*self.indices[0][x].index)+','
            
        index_string=index_string[:-1]+')'
        return ('%s%s'% (self.name,index_string))
    #==========================================================================
    def __repr__(self):
        #printing of the field
        index_string='('
        for x in range(len(self.indices[0])):
            index_string+=str(self.indices[1][x]*self.indices[0][x].index)+','
        index_string=self.name+index_string[:-1]+')'
        return (index_string)


################################################################################
class FieldMul:  #check the printing of the composite fields
    '''
    Represents composite tensors
    '''
    #==========================================================================
    def __init__(self, *args, **kwargs):
        #main attribute is the list of fields as args
        self.fields=list(args)
        
        #name of the composite field is the combination of names of all fields
        name=''
        for x in args:
            name+=x.name
        self.name=name

        #field typ is composite
        self.field_type='composite'

        #if composite field is symmetric/ antisymmetric in indices 
        if 'symmetry' in kwargs.keys():
            self.symmetry=kwargs['symmetry']
        else:
            self.symmetry=0
   
        
    #==========================================================================
    @property  
    def get_indices(self):
        #combine all the indices of all fields in the composite field.
        #In this way, both field and FieldMul has same attribute which gives
        #the indices
        all_idx_sym=[]
        all_idx_sign=[]
        for x in self.fields:
            all_idx_sym.extend(x.get_indices[0])
            all_idx_sign.extend(x.get_indices[1])
        return (all_idx_sym,all_idx_sign)

    #==========================================================================
    @property  
    def get_free_indices(self):
        #return only the free indices from get_indices()
        all_idx=[self.get_indices[1][y]*self.get_indices[0][y].index \
                    for y in range(len(self.get_indices[0]))]
        all_idx_sym=[]
        all_idx_sign=[]
        for x in self.fields:
            all_idx_sym.extend(x.get_indices[0])
            all_idx_sign.extend(x.get_indices[1])
       
        free_idx_sym=[]
        free_idx_sign=[]
          
        for x in range(len(all_idx)):
            if -1*all_idx[x] in all_idx:
                continue
            else:
                free_idx_sym.append(all_idx_sym[x])
                free_idx_sign.append(all_idx_sign[x])
        return(free_idx_sym,free_idx_sign)
    #==========================================================================
    def __str__(self):
        #printing of the field
        string=''
        for x in self.fields:
            string+=x.name + str(tuple(x.get_indices))
        return ('%s' %(string))
    #==========================================================================
    def __repr__(self):
        #printing of the field
        #structure is name1(indices1)name2(indices2)....
        string=''
        for x in self.fields:
            indices=[x.get_indices[1][y]*x.get_indices[0][y].index for y in \
                       range(len(x.get_indices[0]))]
            string+=x.name + str(tuple(indices))
        return (string) 
#################################################################################

#Global attributes: Lorentz Group and Representation
#_______________________________________________________________________________
Lorentz = Group('Lorentz', Abelian=False,latex_index_name='\\mu')
L = Rep('Lorentz', Lorentz, latex_index_name='\\mu')
###################################################################################gg
class ConstantTensor():   #Is it better to rename it as constant tensor with any kind of indices?
    #============================================================================
    def __init__(self, name, **kwargs):
        self.name=name
        self.field_type='mixing Tensor'
        self.rep=[*kwargs['rep']]
        self.explicit_name=explicitname(self.rep)
        self.properties={'rep': self.rep,'field_type':self.field_type,
                         'explicit_name':self.name+self.explicit_name,
                         'latex_name': kwargs["latex_name"], 
                          }
        self.indices1=Create_indices(self.properties['rep'])
        #indices of the conjugate field 
        self.indices2=Create_indices(self.properties['rep'],idx=[1]) 
        self.field=Field(self.name, self.indices1,**self.properties) 
        self.conj= Field(self.name+'_dag',self.indices2, **self.properties)
            #modifying latex name and explicit name for conjugate field
        self.conj.latex_name=self.field.latex_name+'^{\dagger}'
        self.conj.explicit_name=self.conj.name+self.explicit_name
        FieldDict.update({self.field.name: self.field,
                          self.conj.name: self.field})

    #==========================================================================
    def __str__(self):
        #printing of the scalar instance
        return ('constant_field_name:%s ,field_type: %s'\
                                          %(self.name,self.field_type))
#################Global Attributes################################################
Gamma = ConstantTensor('Gamma', rep=[L], latex_name="\gamma", explicit_name = None)
##################################################################################     
class Scalar():
    '''
    Represents the scalar particles
    '''
    #==========================================================================
    def __init__(self, name, self_adjoint=False,**kwargs):
        #main attributes are name, self_adjoint, list of reps
        self.field_type = 'scalar'
        self.name = name
        self.self_adjoint=self_adjoint

        #all the properties other than name is stored here to pass through the
        #field class        
        self.properties = kwargs

        #charge in each rep----------------------------------------------------------only introduced the abelian charge
        if 'charge' in kwargs.keys():
            self.charge=sympify(kwargs['charge'])
        else:
            self.charge=Symbol('Q_'+self.name)

        #To make sure all the reps belong to different groups----------------------------will move to a different method
        check=check_scalar_fermion_rep(self.properties['rep'])
        if check==1:
            raise Exception ("All representations of the field need to be \
                                different and belong to different groups.")


        #explicit name of the scalar field--------------------------------------------move as global function
        self.explicit_name=explicitname(self.properties['rep'])
        #combining all the remaining properties
        self.properties={**self.properties,'field_type':self.field_type,
                          'self_adjoint':self_adjoint,
                          'explicit_name':self.name+self.explicit_name,
                          'charge':self.charge}
        
        #indicesof the field                    
        self.indices1=Create_indices(self.properties['rep'])
        #indices of the conjugate field 
        self.indices2=Create_indices(self.properties['rep'],idx=[1])      
        #scalar field
        self.field=Field(self.name, self.indices1,**self.properties)
        #conjugate field
        #if self adjoint dag field is field itself--------------------------------check for the indices while printing lagrangian
        if self.properties['self_adjoint']==True:
            self.dag=self.field#------------------------------------------------------check if we should change the sign of indices? 

        else:
            self.dag= Field(self.name+'_dag',self.indices2, **self.properties)
            #modifying latex name and explicit name for conjugate field
            self.dag.latex_name=self.field.latex_name+'^{\dagger}'
            self.dag.explicit_name=self.dag.name+self.explicit_name
        self.conj=self.dag
        self.DelField = Del(self.field)
        self.DelConj= Del(self.conj, idx =[1])
        self.CovDelField =CovDel(self.field)
        self.CovDelConj=CovDel(self.conj, idx=[1])
        FieldDict.update({self.field.name: self.field,
                          self.conj.name : self.conj,
                          self.DelField.name: self.DelField,
                          self.DelConj.name: self.DelConj,
                          self.CovDelField.name: self.CovDelField,
                          self.CovDelConj.name: self.CovDelConj,
                         })

    #==========================================================================
    def __str__(self):
        #printing of the scalar instance
        return ('scalar_field_name:%s ,self_adjoint: %s'\
                                          %(self.name,self.field_type))


################################################################################
class Fermion():
    '''
    Represents the fermion fields
    '''
    #==========================================================================
    def __init__(self, name,**kwargs):
        #main attributes includes name, rep list
        #field type is fermion
        self.field_type = 'fermion'
        self.name = name

        #all the attributes aother than name is stored here to pass to the Field
        #class 
        self.properties = kwargs
    
        #charge in the Field-------------------------------------------------------not implemented fully
        if 'charge' in kwargs.keys():
            self.charge=sympify(kwargs['charge'])
        else:
            self.charge=Symbol('Q_'+self.name)

        #check reps belongs to different groups-------------------------------------different function       
        check=check_scalar_fermion_rep(self.properties['rep'])
        if check==1:
            raise Exception ("All representations of the field need to be \
                                different and belong to different groups.")


        #explicit name of the field
        self.explicit_name=explicitname(self.properties['rep'])

        #fields---same structure as Scalar class
        self.properties={**self.properties,'field_type':self.field_type,
                     'explicit_name':self.name+self.explicit_name,
                     'charge':self.charge}
        self.indices1=Create_indices(self.properties['rep'])      
        self.field=Field(self.name, self.indices1,**self.properties)
        self.indices2=Create_indices(self.properties['rep'],idx=[1])
        self.bar=Field(self.name+'_bar', self.indices2,**self.properties)
        self.bar.latex_name='\overline{'+self.field.latex_name+'}'
        self.bar.explicit_name=self.bar.name+self.explicit_name
        self.conj=self.bar
        self.DelField = Del(self.field)
        self.DelConj= Del(self.conj, idx=[1])
        self.CovDelField =CovDel(self.field)
        self.CovDelConj=CovDel(self.conj, idx=[1])
        FieldDict.update({self.field.name: self.field,
                          self.conj.name : self.conj,
                          self.DelField.name: self.DelField,
                          self.DelConj.name: self.DelConj,
                          self.CovDelField.name: self.CovDelField,
                          self.CovDelConj.name: self.CovDelConj,
                         })
      

    #==========================================================================
    def __str__(self):
        #printing of the Fermion class
        return ('fermionic_field_name:%s ,field_type: %s' \
                         % (self.name, self.field_type))






################################################################################
class VBoson():
    
    #==========================================================================
    def __init__(self, name, **kwargs):
        self.field_type='VBoson'
        self.name= name
        
        self.properties=kwargs

        self.properties['rep'].insert(0, L)
        #To make sure all the reps belong to different groups----------------------------will move to a different method
        check=check_scalar_fermion_rep(self.properties['rep'])
        if check==1:
            raise Exception ("All representations of the field need to be \
                                different and belong to different groups.")


        #explicit name of the scalar field--------------------------------------------move as global function
        self.explicit_name=explicitname(self.properties['rep'])
        #combining all the remaining properties
        self.properties={**self.properties,'field_type':self.field_type,
                          'explicit_name':self.name+self.explicit_name,
                         }
        
        #indicesof the field                    
        self.indices1=Create_indices(self.properties['rep'])
        
        #indices of the conjugate field 
        self.indices2=Create_indices(self.properties['rep'],idx=[1])      
        #scalar field
        self.field=Field(self.name, self.indices1,**self.properties)
        
        #conjugate field
        #if self adjoint dag field is field itself--------------------------------check for the indices while printing lagrangian
        if self.properties['self_adjoint']==True:
            self.dag=self.field#------------------------------------------------------check if we should change the sign of indices? 

        else:
            self.dag= Field(self.name+'_dag',self.indices2, **self.properties)
            #modifying latex name and explicit name for conjugate field
            self.dag.latex_name=self.field.latex_name+'^{\dagger}'
            self.dag.explicit_name=self.dag.name+self.explicit_name
        self.conj=self.dag
        FieldDict.update({self.field.name: self.field,
                          self.conj.name: self.conj})
        

    




################################################################################
class GaugeField():
    '''
    Represents Gauge particles
    Creates Gauge Field and corrosponding Field strength Tensor
    '''

    #stores the fields corrosponding to each group

    

    #==========================================================================
    def __init__(self, name, group,**kwargs):
        #major attributes are name, group

        #name should be single character
        if len(name)>1:
            raise Exception ('Gauge field name should be one character length.')
        self.g_name = name.upper()

        self.group = group
        
        #properties to pass to Field 
        self.properties=kwargs

        #field type is gauge
        self.field_type='gauge'

        #latex name of the Field strength tensor
        if 'latex_name' in kwargs.keys():
            self.fs_latex=kwargs['latex_name']
        #default value is the mathbb(name)
        else:
            self.fs_latex='\mathbb{'+self.g_name+'}'

        self.gauge_field=self.field()
        self.fstr=self.field_strength()  #need to add the derivative fields as well.
        FieldDict.update({self.gauge_field.name: self.gauge_field,
                          self.fstr.name: self.fstr})
    #==========================================================================
    def field(self, idx=[1]):
        # Arrange the properties-> create indices-> Create gauge field

        #if group is non abelian, rep contains-[adjoint rep,  lorentz]
        #if group is abelian - 1 lorentz rep
        if self.group.abelian == False:
            properties= {'rep': [self.group.adjoint_rep, L]}
        else:
            properties = {'rep': [L]}

        ename=self.g_name+explicitname(properties['rep'])
        #all the properties

        properties= {**properties, 'field_type': self.field_type,
                     'latex_name':self.g_name,
                     'explicit_name':ename}
          
        #create indices           
        indices=Create_indices(properties['rep'],idx)
        #gauge field
        g_field = Field(self.g_name, indices, **properties)

        return (g_field)

    #==========================================================================
    def field_strength(self, idx=[1]):
        #Field strength tensor same procedure as field() only change is 
        # addition of one more lorentz rep
        if self.group.abelian == False:
            properties = {'rep': [self.group.adjoint_rep, L, L]}
            
        else:
            properties = {'rep': [L, L]}
            
        ename=self.g_name+explicitname(properties['rep'])
        

        properties = {**properties,'field_type':self.field_type,
                                   'latex_name':self.fs_latex,
                                   'explicit_name':ename,
                                   'symmetry': -1}
                                     
        #create indices
        indices=Create_indices(properties['rep'],idx)
        #field strength tensor
        g_field = Field(self.g_name + self.g_name,indices,**properties)

        return (g_field)
    #==========================================================================
    def __str__(self):
        return ('gauge_field_name:%s ,group: {%s}' % (self.g_name, self.group))

    
################################################################################
#8888888888888888888888888888888888888888888888888888888888888888888888888888888
def Del(field,indices=None,idx=[0]):
    '''
    repressnts the partial derivative of a field symbolically
    takes the field as the main argument.
    We can provide either lorentz index or the type of the index as the other 
    argument.
    If index is given, create Del(index) and combine with Field
    if idx is given create a lorentz index and create del(new index) and combine 
    with the field
    '''
     
    properties={'rep':[L],'field_type':'partial','latex_name':'\partial',
                 'explicit_name':None}
    if indices==None:
        index=Create_indices(properties['rep'],idx)
    else:
        index=indices
    Del=Field('Del', index, **properties)
    Del_field=FieldMul(Del,field,symmetry=0)
    return (Del_field)
     
#8888888888888888888888888888888888888888888888888888888888888888888888888888888
def CovDel(field,idx=[0]):
    '''
    Represents the covariant derivative symmbolically
    Slighty different from Del as it takes the field and index type(idx=0/1) as 
    arguments. The other type not introduced because that requirement didn't 
    occur in the code
    '''
    properties={'rep':[L],'field_type':'partial','latex_name':'D',
                'explicit_name':None}
    index=Create_indices(properties['rep'],idx)
    covdel=Field('CovDel', index, **properties)
    CovDel_field=FieldMul(covdel,field, symmetry=0)
    return (CovDel_field)




    
     

        

