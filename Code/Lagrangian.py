'''
This part represnts the second part of the code which deals with the creation 
and modification of the Lagrangian

It has a total of 35 class methods(23 internal, 12 user) ---- This information is wrong need to be corrected.

'''

#------------------------------------------------------------------------------- check if there any additional modules written here
from sympy import Symbol,Rational,latex, sympify,srepr,simplify, I
from sympy.abc import iota
from Components import*
from IPython.core.display import display, Math
from copy import deepcopy
from itertools import permutations, product
import os,sys
from re import findall
from signPermutation import *
from numpy import prod,array
from Utilities import *
from collections import Counter
import datetime
#--------------------------------------------------------------------------------
################################################################################
class Lagrangian():
    #It takes no arguments.
    #---------------------------------------------------------------------------
    def __init__(self):
        #stores the terms of the Lagrangian
        self.L_exp = []     
        #stores the list of scalar fields
        self.L_scalar = [] 
        #stores the list of fermion fields
        self.L_fermion = []
        #Stores the list of groups which are gauged
        self.GaugedGroups = []  
        
        #we need to create dictionaries which can store details of the fields present in a lagrangian
        # rightnow I am going to create a single field dictioary for all types of fields . 
        #Later we can modify it to separate dictioaries if neccessary
        #the keyword of the argument will be the field name and value corresponds to the key is the field instance
        self.FieldDict={}
        

        #list of Vector Bosons
        self.L_VB=[]
        #list of gaugefield instances corrosponding to each group
        self.GaugeFieldInstances = {}
        #list of field strength tensors corrosponding to each group
        self.fstr={}
        #list of gaugefields corrosponding to each groups
        self.gfield={}
        # Stores representations of the Lagrangian.Only adjoint reps are not 
        #includes
        self.rep = []     
        #If gauge flag is true, it prevents the addition of gauge group
        #to the Lagrangian
        self.gauge_flag = False

        #list of balancing tensors
        self.tensor_names = [] 
        
        self.rename_index_dict = {}#---------------------------------------------fill this part
        self.n_dict = {}

        #for the mpping of field strength tensor and covdel
        self.mapping = {'CovDel': { Lorentz: {'coupling': 1, 'delta': [], 
                                    'fields': []}},
                        'fstr': {}}
        #list of groups of the Lagrangian
        self.groups = []
        
    # This function is completely unneccessary as I included it to avoid printing of certian methods.
    #---------------------------------------------------------------------------
    def blockPrinting(func):
        def func_wrapper(*args, **kwargs):
            # block all printing to the constole
            sys.stdout = open(os.devnull, 'w')
            # call the method in question
            value = func(*args, **kwargs)
            # enable all printing to the console
            sys.stdout = sys.__stdout__
            # pass the return value of the method back
            return (value)

        return (func_wrapper)
    

    #---------------------------------------------------------------------------
    def __call__(self):
        '''
        A fuction to give Lagrangian term on calling the Lagrangian 
        instance
        '''
        return (self.L_exp)

    #next 6 functions are to add terms to the lagrangian.
    #---------------------------------------------------------------------------
    def AddGaugeGroup(self, group, coupling, field_name,**kwargs):
        '''
        A functin that adds Kinetic term of the gauge bosons to the 
        Lagrangian. gauge groups should be added before the addition of other
        particle. It takes group, coupling cinstant and name of the field as 
        main argument. Rest of the arguments such as latex name of the 
        fstr can be added as kwargs.
        '''

        #gauge groups can be added only if gauge flag is False
        if self.gauge_flag != False:
            raise Exception ("Gauge fields can no more be added after \
                                adding matter fields!")
       
        self.kwargs=kwargs

        #add the GaugeField instance to the dictionary
        gauge_field_names=[self.gfield[y].name for y in self.gfield.keys()]
        if field_name in gauge_field_names:
            raise Exception ('Field of given name already exists in the Lagrangian')
        else:
            #Create the gauge field and fstr using GaugeField class
            g = GaugeField(field_name, group,**self.kwargs)
            self.GaugeFieldInstances[group] = g
            self.fstr[group]=g.fstr
            self.gfield[group]=g.gauge_field
            self.FieldDict.update({g.fstr.name: g.fstr, 
                                   g.gauge_field.name: g.gauge_field}) #updates the field dict of the Lagrangian

        #only one gauge field per group is allowed--------------------------------->will change later
        if group in self.GaugedGroups:
            raise Exception ('Field of desired group  already exists')
        else:
            self.GaugedGroups.append(group)

        #create the coupling as a symbol
        coupling = Symbol(coupling, Real=True)
        #Stores the couplimg details as a dict for mapping purposes
        self.mapping['CovDel'][group]={'coupling': -I*coupling, 'delta': [],
            'fields': []}
        #self.mapping['CovDel'][group]={'coupling': Rational(-1)*1j*coupling, 'delta': [],
        #    'fields': []}
        #Stores the field details for apping of fstr
        self.mapping['fstr'][group]={'coupling': 1, 'delta': [],
                    'fields': [g.field_strength(), g.field(), coupling]}

        #Add term to the Lagrangian
        self.L_exp.append({'coupling': -Rational(1,4),'delta': [], 
              'fields': [g.field_strength(idx=[1]), g.field_strength(idx=[0])],
              'dirac_pair':[]})
    

    #---------------------------------------------------------------------------
    def AddParticles(self,*args):
        '''
        This is a wrapper function of AddFermion and AddScalar.
        This method enables the addition of a list of fermion and scalars
        in a single step. It basically check the type of the field and pass it 
        to either fermion or boson.
        '''

        for v in args:
            if v.field_type =='scalar':
                self.AddScalar(v)
            elif v.field_type =='fermion':
                self.AddFermion(v)
            else:
                raise Exception("Particle is not defined")
                  
    
    #---------------------------------------------------------------------------
    def AddScalar(self, scalar):
        '''
        Add a single scalar particle to the field
        Take a single scalar field as argument.
        '''

        #List of all the groups in the Particles.
        #The groups in the first particle is considered as the group of the 
        #Lagrangian
        if not self.groups:
            self.groups = [i.group.name for i in scalar.properties['rep']]
        elif self.groups != [i.group.name for i in scalar.properties['rep']]:
            raise Exception ("Poor input! The theory of the Lagrangian is",
                self.groups)
      
        #scalar field
        sf = scalar.field
        #Conjugate field
        sf_dag=scalar.dag
        
        #if gauge group is added-- take covarint derivative, o.w normal derivative
        if self.GaugeFieldInstances:
            sf_deriv = FieldDict["CovDel"+sf.name]
            sf_deriv_conj = FieldDict["CovDel"+sf_dag.name]
            #sf_deriv = CovDel(sf, idx=[0])
            #sf_deriv_conj = CovDel(sf_dag,idx=[1])
        else:
            sf_deriv = FieldDict["Del"+ sf.name]
            sf_deriv_conj = FieldDict["Del"+sf_dag.name]
            #sf_deriv = Del(sf, idx=[0])
            #sf_deriv_conj = Del(sf_dag,idx=[1])
            print(sf_deriv_conj)

        #add the scalarf field to the L_scalar list
        if len(self.L_scalar) == 0:
            self.L_scalar.append(sf)
        else:
            for key in self.L_scalar:
                if sf == key:
                    raise Exception('The desired field is already present')
            self.L_scalar.append(sf)
        self.FieldDict.update({sf.name: sf,
                               sf_dag.name: sf_dag}) #updates the field_dict

        #put gauge flag to 0, which prevents the further addition of gaugegroups
        self.gauge_flag = True
       
        #create the list of reps in the scalar field
        for i in scalar.properties['rep']:
            if i.name not in self.rep:
                self.rep.append(i.name)
         
        #if self_adjoint--> coupling is 1/2 else 1
        if scalar.properties['self_adjoint']==True:
            coupling=Rational(1,2)
        else:
            coupling=1

        #creating the klein Gordan Lagrangian
        scalar_term={'coupling':coupling, 'delta':[], 
                                'fields':[sf_deriv_conj,sf_deriv],
                                'dirac_pair':[]}
        # add the term to the Lagrangian
        self.L_exp.append(scalar_term)
        
    #---------------------------------------------------------------------------
    def AddFermion(self, fermion):
        '''
        Add a single fermion to the Lagrangian. Similar structure as 
        AddScalar
        '''

        #list the groups in the Lagrangian
        #Groups of th first particle is the groups of lagrangian
        if not self.groups:
            self.groups = [i.group.name for i in fermion.properties['rep']]
        elif self.groups != [i.group.name for i in fermion.properties['rep']]:
            raise Exception ("Poor input! The theory of the Lagrangian is",
                self.groups)


        #fermion field
        ff = fermion.field
        #conjuagte field
        ff_bar=fermion.bar
        self.FieldDict.update({ff.name: ff}) #updates the field dict

        #if gauge group is added--> covdel else--> del
        if self.GaugeFieldInstances:
            #ff_deriv = CovDel(ff, idx=[0]) #need to replace with the default field created
            ff_deriv = FieldDict["CovDel"+ff.name] #need to replace with the default field created
        else:
            #ff_deriv = Del(ff, idx=[0]) #need to replace with the default field created
            ff_deriv = FieldDict["Del"+ff.name] #need to replace with the default field created

        #add the fermion to the list of fermion
        if len(self.L_fermion) == 0:
            self.L_fermion.append(ff)
        else:
            for key in self.L_fermion:
                if ff == key:
                    raise Exception('The desired field is already present')
            self.L_fermion.append(ff)

        #prevents addition of extra gauge groups
        self.gauge_flag = True

        #gamma matrix
        gamma=Gamma.field
        gamma.indices[1][0]=1
        #list of reps is added to the rep list
        for i in fermion.properties['rep']:
            if i.name not in self.rep:
                self.rep.append(i.name)
        #couping constant
        coupling=I
        #creating dirac lagrangian
        fermion_term={'coupling':coupling,'delta':[],
                       'fields':[ff_bar, gamma, ff_deriv],
                       'dirac_pair':[[ff_bar.name,gamma.name, ff_deriv.name]] }

        #add term to the lagrangian
        self.L_exp.append(fermion_term)


    #---------------------------------------------------------------------------
    def AddVBoson(self,vboson): #need to change this part later
        del_vb=Del(vboson.field, idx=[0])
        del_vb_dag=Del(vboson.dag, idx=[0]) #replace with the default field 
        l= len(del_vb.get_indices[0])
        print(l)
        t1_up=[x for x in range(1, l+1)]
        t1_down=[-x for x in t1_up]
        t2_up= [ x for x in t1_up]
        t2_up[0], t2_up[1]= t2_up[1], t2_up[0]
        t2_down=t1_down
        print(t1_up, t2_up, t1_down, t2_down)
        
        if len(self.L_VB)==0:
            self.L_VB.append(vboson.field)
        else: 
            for key in self.L_VB:
                if key==vb:
                    raise Exception('the desired vector boson is already present')
            self.L_VB.append(vb)
        self.AddTerm('-1/2', del_vb_dag, del_vb,
                          contraction_pattern=[t1_up, t1_down])
        self.AddTerm('1/2', del_vb_dag, del_vb, 
                          contraction_pattern=[t2_up, t2_down])

        

    #---------------------------------------------------------------------------
    def AddMass(self, field, mass):
        '''
        depending on the type of field(scalar or fermion), mass 
        term is added.
        for scalar fields we consider whether it is self adjoint or not.
        gauge_flag is set to true, so that no gauge fields can be added 
        after the addition 
        of mass term
        '''
        #mass sting is converted to symbol
        mass = Symbol(mass, Real=True)
        self.gauge_flag = True

        #for scalar field
        if field.properties['field_type']== 'scalar':
            m_field = field.field
            m_field_dag=field.dag
            #coupming constant for scalar field
            if field.properties['self_adjoint'] == True:
                coupling=-Rational(1, 2) * mass ** 2
            else:
                coupling=-mass ** 2
            term={'coupling':coupling,'delta':[],
                  'fields':[m_field_dag,m_field],
                  'dirac_pair':[]}

       
        #for fermion field
        elif field.properties['field_type'] == 'fermion':
            m_field = field.field
            m_field_bar=field.bar
            #coupling constant
            coupling=-mass
            term={'coupling':coupling,'delta':[],
                  'fields':[m_field_bar,m_field],
                  'dirac_pair':[[m_field_bar.name,m_field.name]] }
        self.L_exp.append(term)
    

    #-----------------------------------------------------------
    def count_fermion_pairs(self, list1):
        d_pair=[]
        fermion_count = {"bar":0, "not_bar": 0}
        for x in list1:
            if isinstance(x, FieldMul):
                f = x.fields[1]
            else:
                f = x
            if f.field_type == 'fermion' :
                if '_bar' in f.name:
                    fermion_count['bar'] += 1
                else:
                    fermion_count['not_bar'] += 1
                d_pair.append(x)
        return(fermion_count, d_pair)
    #---------------------------------------------
    def convert_str_to_fields(self, list1):
        #new function to convert fields inputed as strings into their corresponding fields
        field_list = []
        #first deal with single component fields
        for x in list1:
            count=Counter(x)
            print(count)
            bracket_number=count['('] 
            print(bracket_number)
            if bracket_number == 2:
                field1, field2, extra = x.split(')')
                f1_name, f1_indices = field1.split('(')
                f2_name, f2_indices = field2.split('(')
                print(field1, field2)
                print("f1",f1_name, f1_indices)
                print("f2",f2_name, f2_indices)
                #identify the field from FieldDict
                fields= FieldDict[f1_name+f2_name].fields
                print(fields)
                if ',' in f1_indices:
                    print("if")
                    f1_splitted_indices = f1_indices.split(',')
                else:
                    print("else")
                    f1_splitted_indices = [f1_indices]
                if ',' in f2_indices:
                    f2_splitted_indices = f2_indices.split(',')
                else:
                    f2_splitted_indices = [f2_indices]
           
                print("ss",f1_splitted_indices, f2_splitted_indices)
                f1_indices_sym = []
                f1_indices_sign = []
                for y in range(len(f1_splitted_indices)):
                    print(y)
                    if '-' in f1_splitted_indices[y]:
                        z1 = f1_splitted_indices[y].replace('-', '')
                        f1_indices_sign.append(-1)
                    else:
                        z1 = f1_splitted_indices[y]
                        f1_indices_sign.append(1)
                    index1 = Index(z1, \
                             Rep.rep_dict[fields[0]\
                             .get_indices[0][y].index_type])
                    f1_indices_sym.append(index1)
                f1_new_indices = [f1_indices_sym,\
                                  f1_indices_sign]
                f1_new = Field(fields[0].name, f1_new_indices,
                    field_type=fields[0].field_type,
                    latex_name=fields[0].latex_name,
                    explicit_name=fields[0].explicit_name,
                    symmetry=fields[0].symmetry, )

                f2_indices_sym = []
                f2_indices_sign = []
                for y in range(len(f2_splitted_indices)):
                    if '-' in f2_splitted_indices[y]:
                        z2 = f2_splitted_indices[y]\
                             .replace('-', '')
                        f2_indices_sign.append(-1)
                    else:
                        z2 = f2_splitted_indices[y]
                        f2_indices_sign.append(1)
                    index2 = Index(z2, \
                             Rep.rep_dict[fields[1]\
                             .get_indices[0][y].index_type])
                    f2_indices_sym.append(index2)
                f2_new_indices = [f2_indices_sym,\
                                  f2_indices_sign]
                f2_new = Field(fields[1].name, f2_new_indices,
                    field_type=fields[1].field_type,
                    latex_name=fields[1].latex_name,
                    explicit_name=fields[1].explicit_name,
                    symmetry=fields[1].symmetry, )
                new_field = FieldMul(f1_new, f2_new)
                print(new_field)
            else: 
                name, indices = x.split('(')
                #Identify the field from the FieldDict
                field = FieldDict[name]
                field_indices=field.get_indices[0]
                #separate indices and remove bracket
                splited_indices = indices.split(',')
                splited_indices = [s.replace(')', '') for s in\
                               splited_indices]
                #create new indices
                #need to include certain checks which can be added
                # later
                indices_sym = []
                indices_sign = []
                for y in range(len(splited_indices)):
                    if '-' in splited_indices[y]:
                        z = splited_indices[y].replace('-','')
                        indices_sign.append(-1)
                    else:
                        z = splited_indices[y]
                        indices_sign.append(1)
                    index= Index(z, \
                        Rep.rep_dict[field_indices[y].index_type])
                    indices_sym.append(index)
                #need to remove the extra spaces in strings
                new_indices=[indices_sym, indices_sign]
                #need to add other properties to the field
                new_field = Field(field.name, new_indices,
                    field_type=field.field_type,
                    latex_name=field.latex_name,
                    explicit_name=field.explicit_name,
                    symmetry=field.symmetry, )
            field_list.append(new_field)
        return(field_list)      
                
    #----------------------------------------------
    def contraction(self, coupling, field_list):
        term = {'coupling': Symbol(coupling, Real = True),
                'delta' : [],
                'fields': [*field_list],
                'dirac_pair': []
               }
        return(term)
    #---------------------------------------------            
    def add_term_manage(self, coupling, list1):
        #This is the method in which the code generates constant
        #tensors to create contraction pattern.
        #In this method the fields should not be given as string
        term = {'coupling' : Symbol(coupling, Real = True),
                'delta': [],
                'fields': [],
                'dirac_pair': []
               }
        fermion_count = {'bar': 0, 
                         'not_bar': 0
                        }
        index = [[],[]]
        dirac_pair = []
        for x in list1:
            if isinstance(x, Field) or isinstance(x, FieldMul):
                indices = x.get_indices
                index[0].extend(indices[0])
                index[1].extend(indices[1])
            else:
                raise Exception ("only fields are allowed")
        print(index)  
        new_index = [[], []]
        rep_list = []
        for i in range(len(index[0])):
            if Rep.rep_dict[index[0][i].index_type].singlet \
            == False:
                new_index[0].append(index[0][i])
                new_index[1].append(index[0][i])
                rep_list.append(index[0][i].index_type)
        group_list = [Rep.rep_dict[x].group.name \
                      for x in rep_list]
        index_dict = {}
        for i in range(len(group_list)):
            try: 
                index_dict[group_list[i]][rep_list[i]] += 1
            except NameError:
                index_dict = {group_list[i]: {rep_list[i]: 1}}
            except KeyError:
                try:
                    isinstance(index_dict[group_list[i]], dict)
                except KeyError:
                    index_dict.update({group_list[i]: \
                                       {rep_list[i]: 1}\
                                      })
                else: 
                    index_dict[group_list[i]].update(\
                              {rep_list[i]:1})
        if index_dict == {}:
            term['fields'] = [*list1]
            return (term)
        else:
            temp = {}
            for z, y in index_dict.items():
                one_g_rep = []
                bal_indices = [[],[]]
                for i, j in y.items():
                    one_rep = [[Index(Rep.rep_dict[str(i)].\
                                index_name + str (k+1), \
                                Rep.rep_dict[str(i)]) \
                                for k in range(j)], \
                                [1 for k in range(j)]]
                    one_g_rep.extend([i for k in range(j)])
                    temp.update({i: one_rep})
                if len(one_g_rep) == 1:
                    raise Exception ( 'only one index \
                                      for rep present')
                elif len(one_g_rep) == 2:
                    if one_rep[0][0].index_type != \
                    one_rep[0][1].index_type:
                        raise Exception ('Entered term :null')
                    delta = self.create_delta(\
                            [[one_rep[0][0]],\
                             [one_rep[1][0]]],\
                            [[one_rep[0][1]],\
                             [one_rep[1][1]]])
                    term['delta'].extend(delta)
                else:
                    e_name = '_x'
                    for x in one_g_rep:
                        e_name += '_' + x
                    rep_list = [Rep.rep_dict[ str(i)] for i in \
                                             one_g_rep ]
                    if rep_list != []:
                        indices= Create_indices(rep_list, \
                                               idx=[1])
                        if coupling[0] == '-':
                            c= coupling[1:]
                        else:
                            c = coupling
                        f = Field(c + '_'+ str(z),
                                  indices,
                                  field_type = 'constant',
                                  latex_name = '(C_{' +  \
                                                str(z) + '})',
                                  explicit_name = c + str(z) + \
                                                  e_name,
                                  symmetry = 0,
                                  numerical_values = None)
                        
                    term['fields'].append(f)
                    self.tensor_names.append(f)
            for x in list1:
                if isinstance(x, FieldMul):
                    y=x.fields
                    new_field=[]
                else:
                    y=[x]
                    new_field=[]
                for j in range(len(y)):
                    ar = [[],[]]
                    idcs = y[j].get_indices

                    for i in range(len(idcs[0])):
                        if idcs[0][i].index_type in temp.keys():
                            ar[0].append(temp[idcs[0][i].\
                                       index_type][0].pop())
                            ar[1].append(-1)
                        elif idcs[1][i]==-1:
                            ar[0].append(idcs[0][i])
                            ar[1].append(-1)
                        else:
                            ar[0].append(idcs[0][i])
                            ar[1].append(-idcs[1][i]) 

  
                    new_field.append(Field(y[j].name, ar, \
                    field_type=y[j].field_type,
                    latex_name=y[j].latex_name,
                    explicit_name=y[j].explicit_name,
                    symmetry=y[j].symmetry, 
                    numerical_values=y[j].numerical_values))                
                if len(y) == 1:
                    term['fields'].append(new_field[0])
                else:
                    term['fields'].append(FieldMul(*new_field))
        if term['delta']:
            term= self.eliminate_delta(term)
        return(term)
    #-----------------------------------------------------      
    def AddTerm(self, coupling, list1, manage= True):
        #this the method that adds interaction terms to the 
        #Lagragian. It should be done in 2 ways. One the code  
        # will automatically generates the contraction pattern 
        #by creating Constant tensor. Second by user providing
        # the contraction along with the field in  the form of
        # strings.
        #
              
        if manage == False :
            #Since the fields are given as strings the first
            #step is to convert them in to fields 
            field_list=self.convert_str_to_fields( list1)
            #check the number of fermion pairs
            #check the contraction is legit
            term = self.contraction(coupling, field_list)
        else:
            term = self.add_term_manage(coupling, list1)
        #return(term)
        self.L_exp.append(term)

                  
    #-----------------------------------------------------------
    def create_delta(self, nlist1, nlist2):      #no connection with Lagrangian 
        '''
        Create kronecker deltas with set of indices passed.
        Creates only for non singlet indices
        '''
        list1=[[],[]]
        list2=[[],[]]
        for x in range(len(nlist1[0])):
            if Rep.rep_dict[nlist1[0][x].index_type].singlet==False:
                list1[0].append(nlist1[0][x])
                list1[1].append(nlist1[1][x])

        for x in range(len(nlist2[0])):
            if Rep.rep_dict[nlist2[0][x].index_type].singlet==False:
                list2[0].append(nlist2[0][x])
                list2[1].append(nlist2[1][x])

        #length of the arguments should match as
        #KD is created from one index from first and another index from second
        #lists
        if len(list1) != len(list2):
            raise Exception(" kds cannot be created")

        delta_list = []

        for i in range(len(list1[0])):
            index_type=list1[0][i].index_type
            if list1[0][i].index_type != list2[0][i].index_type:
                raise Exception("Index type should match for the Kronecker\
                                 delta creation")   
            indices=[[list1[0][i],list2[0][i]],[list1[1][i],list2[1][i]]]
            #create delta from the method in Rep Class                      
            delta_list.append(Rep.rep_dict[index_type].delta(indices))  
        return (delta_list)
                         
    #---------------------------------------------------------------------------
    def eliminate_delta(self, expression): #--------------------------------------------can be removed and moved to utilities
        '''
        delta and metric tensors are considered equally.At present, the function 
        only considers the 2 indices of the tensors in the delta part of the 
        dictionary(not the name--> whether it is KD or metric tensor)
        For the 2 indices of the delta tensor, code go through all the terms of 
        the lagrangian.
        If it finds a matching index, it will replace that index with the other 
        index of the delta.
        '''

        f = expression['fields']        
        d_holder = []        
        #check each delta against all the fields
        for d in expression['delta']:  
             
            flag1=0
            i1=d.get_indices[0]
            i2=d.get_indices[1]  

            for num2 in range(len(f)):

                #create a field in list---> to consider both Field and FieldMul
                #together
                if isinstance(f[num2],FieldMul):
                    fields=f[num2].fields
                else:
                    fields=[f[num2]]

                for m in range(len(fields)):
                    l1=fields[m].get_indices[0]
                    string_list=[str(x) for x in l1]
                    l2=fields[m].get_indices[1]
                    #compare first index of the KD to indices of the field
                    #if matches, replace with the second KD index
                    if str(i1[0]) in string_list:
                        pos=string_list.index(str(i1[0]))
                        if l2[pos]!=i2[0]:
                            l1[pos]=i1[1]
                            l2[pos]=i2[1]
                            flag1=1
                    #compare second index of the KD to indices of the field
                    #if matches, replace with the first KD index
                    elif str(i1[1]) in string_list:
                        pos=string_list.index(str(i1[1]))
                        if l2[pos]!=i2[1]:
                            l1[pos]=i1[0]
                            l2[pos]=i2[0]
                            flag1=1
                    #create new field with new indices, if indices list is 
                    #modified
                    if flag1==1:
                        new_field=Field(fields[m].name,(l1,l2),
                                        latex_name=fields[m].latex_name,
                                        field_type=fields[m].field_type,
                                        explicit_name=fields[m].explicit_name,
                                        symmetry=fields[m].symmetry,
                                        numerical_values=fields[m].numerical_values)
                        fields[m]=new_field
                        break
                if flag1==1:
                    if len(fields)==1:
                        f[num2]=fields[0]
                    else:
                        f[num2]=FieldMul(*fields)

                    #if KD is removed, that is separated from the original KD list
                    #so that, it can be removed from the delta list of the 
                    #expression
                    d_holder.append(d)
                    break
                
                                   
        temp = []
        for d in expression['delta']:
            if d not in d_holder:
                temp.append(d)
        expression['delta'] = temp
        return (expression)


#-------------------------------------------------------------------------------
    def LatexPrint(self,dict1, **kwargs): 

        '''
        For the latex printing of the output in a jupyter notebook/in a file
        Pass the term the term to be printed as argument.
        If needed user can provide the name of the output file.
        
        '''
        #defualt file name is the combination of groups in the lagrangian
        if 'file_name' in kwargs.keys():
            file_name=kwargs['file_name']
        else:
            
            #if only gauge groups are there, name is 'fstr'------fieldstrength
            if self.groups==[]:
                file_name='fstr.txt'
            else:
                grp=set(self.groups)
                s='_'
                output=s.join(grp)
                file_name=str(output)+'.txt'

        #initial latex commands
        file1=open(file_name,'w')
        file1.write('\\documentclass[9pt]{article}\n')
        file1.write('\\usepackage{amsmath,amssymb,geometry}\n')
        file1.write('\\allowdisplaybreaks\n')
        file1.write('\\begin{document}\n')
        

        file1.write('\\begin{align*}\n')
        #for the printing in jupyter notebook
        final=''
        
        for x in dict1:
            #if str(x['coupling'])==sympify('0'):
            #    continue
            #for coupling-------------------------------------------------------> may need modification if the coupling becomes complicated
            #if ('iota*iota' in str(x['coupling'])) or \
             #   ('iota**2' in str(x['coupling'])):
              #  x['coupling']= x['coupling'].subs(iota*iota,-1) 
            if str(x['coupling'])=='1':
                coupling='+'+' '
            elif srepr(sympify(str(x['coupling'])))[:3]=='Add':
                coupling='+('+latex(sympify(str(x['coupling']))) +')' + ' '
            elif str(sympify(str(x['coupling'])))[0]!='-':
                coupling= '+ '+latex(sympify(str(x['coupling']))) + ' '
            else:
                coupling=latex(sympify(str(x['coupling']))) + ' '
          
           #for the delta part --create delta for each KD in the list
            D=''
            if x['delta']!=[]:
                for n in x['delta']:
                    delta='\delta'
                    d_index=self.up_down_index(n)
                    delta+=d_index[0] +d_index[1]
                    D+=delta
            term=''
           #for the field part:
            for field in x['fields']:      
                if field==None:
                    continue                      
                field_name=' '

                #composite fields are combined together
                if isinstance(field,FieldMul):
                    f=field.fields
                    name1=f[0].latex_name
                    name2=f[1].latex_name
                    index1=self.up_down_index(f[0])
                    index2=self.up_down_index(f[1])
                    up1=index1[0]
                    down1=index1[1]
                    up2=index2[0]
                    down2=index2[1]
                    field_name=name1
                    if up1!='^{}':
                        field_name +=up1
                    if down1!='_{}':
                        field_name+=down1

                    #dagger symbol needed to be put outside the indices 
                    if '\dagger' in name2:
                        dagger='^{\dagger}'
                        name2=name2.replace(dagger,'')
                        field_name+=' ('+name2
                        if up2!='^{}':
                            field_name +=up2
                        if down2!='_{}':
                            field_name+=down2 
                        field_name+=')'+dagger
                        
                      
                    else:
                        field_name+=name2
                        if up2!='^{}':
                            field_name +=up2
                        if down2!='_{}':
                            field_name+=down2
                    if ('D' in field_name )and ('\dagger' in field_name):
                        field_name=field_name.replace('(','')
                        field_name="("+field_name
                    
                #for non composite fields  
                else:
                    name=field.latex_name                                          
                    index_name=self.up_down_index(field)
                    up=index_name[0]
                    down=index_name[1] 
                    #indices=field.get_indices
                   

                    if 'T' in field.latex_name \
                                          and index_name[2]==[]:
                        if '_x_' not in field.name:
                            continue                                             
                    if '\dagger' in name:
                        dagger='^{\dagger}'
                        name=name.replace(dagger,'')
                        field_name+=' ('+name
                        if up!='^{}':
                            field_name +=up
                        if down!='_{}':
                            field_name+=down
                        field_name+=')'+dagger
                    else:
                        field_name+=name
                        if up!='^{}':
                            field_name +=up
                        if down!='_{}':
                           field_name+=down
                  
                term+=field_name
            final+=coupling+D+term
            #in the file code prints each term in a line
            file1.write('&'+coupling+D+term+'\\\\')
            file1.write("\n")
            
            #uses multiple align in latex for printing bigger terms
            if (dict1.index(x)%299==0) and (dict1.index(x)!=0):
                file1.write('\\end{align*}\n')
                file1.write('\\begin{align*}\n')
            
        if final[0]=='+':
            final=final[1:]
        file1.write('\\end{align*}\n')
        file1.write('\\end{document}\n')
        file1.close()
        #output prints only in jupyter if the length terms is less than 300
        if len(dict1)<=300:
           display(Math(final))
           os.remove(file_name)
    #----------------------------------------------------------------------------      
    def up_down_index(self,field):  #can be taken out accessed from latex_print----------------------
        indices=field.get_indices 

        #make a list of non singlet indices
        non_singlet_index=[]
        for i in range(len(indices[0])):
            if Rep.rep_dict[indices[0][i].index_type].singlet==False:
                non_singlet_index.append(indices[0][i].index*indices[1][i])
       

        up='^{'
        down='_{'   
        #if index is covariant add it to down string,  o.w up string               
        for o in non_singlet_index:
            if '-' in str(o):
                down+=str(o)[1:-1]+'_'+str(o)[-1]
            else:
                up+=str(o)[:-1]+'_'+str(o)[-1]
        up+='}'
        down+='}'
        return(up,down, non_singlet_index)#---------------------------------------------------------add one more attribute to combine related parts in latexprint

    #---------------------------------------------------------------------------
    @blockPrinting
    def FuncDeriv(self,lag,*args,constant=False):
        '''
        gets the terms to be differentiated and 
        w.r.t what to differentiate as argument.
        Has the keyword constant=True/False
        True means  return final expression having only coupling 
        False return all terms after differentaition
        
        '''                  
        X,Y=[],[]

        #creates the dictionary of fields passed as arguments
        list1={}
        for v in args:
            #derivative cannot be taken with constant tensors
            if v.field_type=='constant':
                raise Exception("cannot take derivative w.r.t constant field")            
            try:
                list1[v.name]+=1
            except KeyError:
                list1.update({v.name:1})  
        # create the dictionary of fields for each term of lagrangian and
        #compare with the fields in args.
        for i in range(len(lag)):
            list2={}
            for m in lag[i]["fields"]:                               
                try:
                    list2[m.name]+=1
                except KeyError:
                    list2.update({m.name:1})
            #compare the number of fields in each term with the arguments passed
            if  list1.keys()<=list2.keys(): 
                l1=[list1[x] for x in list1.keys()]
                l2=[list2[x] for x in list1.keys()]
                diff=[l2[x]-l1[x] for x in range(len(l1))]
                if min(diff)>=0:
                    X.append(lag[i])
                if list1.keys()==list2.keys():
                    l1=[list1[x] for x in list1.keys()]
                    l2=[list2[x] for x in list2.keys()]
                    if l1==l2:
                        Y.append(lag[i])   
        #if the user demands only constant output
        if constant:
            return self.func_deriv_body(Y,*args)
        else:
            return self.func_deriv_body(X,*args)

    #---------------------------------------------------------------------------
    def func_deriv_body(self,l, *args):
        '''
        Take the fields with respect to which we differentiate the lagrangian 
        expression as arguments.
        For each argument the count increments by one.
        count is used to remame the indices of the fields provided as arguments.
        To rename the  multi component field we first split the field. After the 
        modification is made, we combine them back.
        To take the functional derivative, for each argument we check through 
        the Lagrangian and replace with the delta tensor when we find a match.
        Use eliminate_delta to make the proper contractions at the end.
        '''
        count=0
        #for each field in the args list
        for v in args:
            #increment the count for renaming the indices
            count+=1
            
            if isinstance(v,FieldMul):
                fields=v.fields
            else:
                fields=[v]
            #creating new indices for the field 
            for j in range(len(fields)):
                x=fields[j]
                i1=x.get_indices[0]
                i2=x.get_indices[1]
                new_index_sym=[]
                for y in i1:
                    type1=y.index_type
                    #for lorentz field name sigma is given
                    if type1=='Lorentz':
                        new_name='\sigma'+str(count)
                    else: 
                        new_name=Rep.rep_dict[type1].index_name.upper()\
                                 +str(count)
                    #creating the new index
                    new_index=Index(new_name, Rep.rep_dict[type1])
                    new_index_sym.append(new_index)
                #creating new field with new indices
                new_field=Field(x.name, (new_index_sym,i2),
                                latex_name=x.latex_name,
                                field_type=x.field_type,
                                explicit_name=x.explicit_name,
                                symmetry=x.symmetry,
                                numerical_values=x.numerical_values)
                fields[j]=new_field
            if len(fields)==1:
                v=fields[0]
            else:
                v=FieldMul(*fields)                                          
            f_expression = [] 
            #Taking functional derivative with one field
            for i in range(len(l)):
                for j in range(len(l[i]['fields'])):
                    #comparing the field name with the name in the term
                    if v.name == l[i]['fields'][j].name:                        
                        fermion_flag=0
                        if v.field_type=='fermion':
                            fermion_flag=1
                        if isinstance(v, FieldMul) and  \
                             v.fields[1].field_type=='fermion':
                            fermion_flag=1
                        temp_dict = deepcopy(l[i])
                        if fermion_flag==1:
                            #to get the sign after differentiation
                            a = self.count_fermion(l[i]['fields'],\
                                       l[i]['dirac_pair'], j) 
                       
                            temp_dict=self._func_deriv_dirac_pair(temp_dict,v,j, a[1])  
                            if a[0]%2==0:
                                temp_dict['coupling']*=-1    
                                                                                                      
                        #replace the field with delta function
                        delta =self.create_delta((v.get_indices[0],\
                            [-1*x for x in v.get_indices[1]]),
                            l[i]['fields'][j].get_indices)
                        if delta:
                            temp_dict['delta'].extend(delta)
                        del temp_dict['fields'][j]
                        
                          
                        f_expression.append(temp_dict)
            expr = []
            #removing the delta function
            for fe in f_expression:
                if fe['delta']:
                    expr.append(self.eliminate_delta(fe))
                else:
                    expr.append(fe)
            l = expr
        final_expression = []
        couplings = []
        for exp in l:
            if exp not in final_expression:
                couplings.append(l.count(exp))
                final_expression.append(exp)

        for i in range(len(final_expression)):
            final_expression[i]['coupling'] *= couplings[i]

        return(final_expression)

    #--------------------------------------------------------------------------
    def count_fermion(self, lag, dirac_p, position):
        #To count the number of fermions in the term
        #stores the total umber of fermions
        count_fermion=0
        #Stores the number of free fermions
        count_free_fermion=0
        #position is the position of fermion in the lagrangian term
        for x in range(position + 1):
            #if the required field is fermion-- count imcrements
            if isinstance(lag[x], FieldMul):
                if lag[x].fields[1].field_type=='fermion':
                    count_fermion+=1
            else:
                if lag[x].field_type=='fermion':
                    count_fermion+=1

            #check the dirac pair list if it is free or part of boson pair
            for y in dirac_p:
                #if its part of tuple inthe dirac pair, the free fermion
                #count increments
                if lag[x].name in y:
                    if isinstance(y, list):
                        pos=y.index(lag[x].name)
                        length=len(y)
                        str1=[]
                        if pos==0:
                            for p in range(x, x+length):
                                str1.append(lag[p].name)
                        else:
                             for p in range(x-length+1, x+1):
                                str1.append(lag[p].name)
                        if str1==y:
                             break
                  
                    elif isinstance(y,tuple):
                        z=list(y)
                        pos=z.index(lag[x].name)
                        count_free_fermion+=1
                        break
        return(count_fermion, count_free_fermion)
            
                        
                            
        
    #--------------------------------------------------------------------------
    def _func_deriv_dirac_pair(self, lag,v,j, count_tuple):
        #modification of dirac pairs in a term is required after
        #differentiating with a fermion.
        #dirac pair has 2 type of data inside. a list denotes the a bosonic pair.
        #a tuple list all the free fermions
        for x in lag['dirac_pair']:

            #if the field is inside a list remove it fro there and
            #move rest of the list to the tuple
            if v.name in x:
                if isinstance(x, list):
                    pos=x.index(v.name)
                    length=len(x)
                    str1=[]
                    if pos==0:
                    
                        for p in range(j, j+length):
                            str1.append(lag['fields'][p].name)
                    else:
                                       
                        for p in range(j-length+1, j+1):
                            str1.append(lag['fields'][p].name)
                    y=deepcopy(x)
                    if (y==str1):
                        y.remove(v.name)
                    else: 
                        continue
                    lag['dirac_pair'].remove(x)
                    if y==[]:
                        break
                    else:
                        if lag['dirac_pair']==[]:
                            y=tuple(y)  
                            lag['dirac_pair'].append(y) 
                        elif isinstance(lag['dirac_pair'][-1], list):
                            y=tuple(y)  
                            lag['dirac_pair'].append(y)
                        elif isinstance(lag['dirac_pair'][-1], tuple):
                            z=lag['dirac_pair'][-1]
                            str2=''
                            for a in y:
                                str2+=a + '+'
                            lag['dirac_pair'].remove(z)     
                            z=list(z)
                            z.insert(count_tuple, str2[:-1])
                            z=tuple(z)
                            lag['dirac_pair'].append(z)   
                        break
                    
                #if the field is in a tuple-- remove it   
                else:                                    
                    y=list(x)                                      
                    y.remove(v.name)
                    lag['dirac_pair'].remove(x)
                    if y!=[]:
                        lag['dirac_pair'].append(tuple(y))
                    break
                                                
               
        return(lag) 

    #---------------------------------------------------------------------------
    #@blockPrinting
    def Expand(self, keyword):
        '''
        expand covdel and fstr using mapping function.
        
        '''
        terms_tbd = []
        #to expand the field strength tensor
        if keyword.lower() == 'fstr':
            expression = []
            #list of field strength names
            fstrs = [self.GaugeFieldInstances[y].fstr.name 
                                for y in self.GaugeFieldInstances.keys()]
            #for each  term in the Lagrangian
            for exp in self.L_exp:
                self.rename_index_dict = {}
                temp_list1 = []
                temp_list2 = []
                #for each field in the term
                for field in exp['fields']:
                    temp_flag = 0
                    #if the name of the field strength is in the term
                    if field.name in fstrs:
                        temp_flag=1
                        #add the mapping to the temp_list1
                        temp_list1.append(self.mapping_function('fstr',field))
                    #else add the field  to temp list2
                    if temp_flag == 0:
                        temp_list2.append(field)
                    else:
                        temp_list2.append(None)
                k = 0
                true_list = []
                temp = deepcopy(exp)
                temp['fields'] = []
                #combine both temp list 1and 2 
                for i in range(len(temp_list2)):
                    if temp_list2[i] != None:
                        temp['fields'].append(temp_list2[i])
                        temp['coupling']=1      
                    else:
                        if temp['fields']:
                            t = deepcopy(temp)
                            true_list.append([t])
                            temp['fields'] = []
                        true_list.append(temp_list1[k])
                        k += 1
                if temp['fields']:
                    true_list.append([temp])
                #will write about it later
                if temp_list1:
                    a = []
                    for i in true_list:
                        if a == []:
                            a = i
                        else:
                            a = self.field_multiplication(a, i)
                    terms_tbd.append(exp)
                    for temp_var in a:
                        temp_var['coupling'] = \
                                            exp['coupling']*temp_var['coupling']
                        
                    expression.append(a)
                
            for i in expression:
                for j in i:
                    self.L_exp.append(j)
        #For replacing the covariant derivative
        if keyword.lower()=='covdel':
            expression = []
            for exp in self.L_exp:
                self.rename_index_dict = {}
                temp_list1=[]
                temp_list2=[]
                for field in exp['fields']:
                    temp_flag=0
                    if isinstance(field, FieldMul):
                        y=field.fields
                        if y[0].name=='CovDel':
                            position=exp['fields'].index(field)
                            #for fermions change of dirac pair is required
                            if y[1].field_type=='fermion':
                                for x in exp['dirac_pair']:
                                    if field.name in x:
                                        x.remove(field.name)
                            #follow the same procedure as fstr
                            temp_flag=1
                            temp_list1.append(self.mapping_function(\
                                              'covdel',y[1],y[0]))
                    
                    if temp_flag == 0:
                        temp_list2.append(field)
                    else:
                        temp_list2.append(None)
                
                k = 0
                true_list = []
                temp = deepcopy(exp)
                temp['fields'] = []
                for i in range(len(temp_list2)):
                    if temp_list2[i] != None:
                        temp['fields'].append(temp_list2[i])
                    else:
                        if temp['fields']:
                            t = deepcopy(temp)
                            true_list.append([t])
                            temp['fields'] = []
                        true_list.append(temp_list1[k])
                        k += 1
                if temp['fields']:
                    true_list.append([temp])
                

                if temp_list1:
                    a = []
                    for i in true_list:
                        if a == []:
                            a = i
                        else:
                            a = self.field_multiplication(a, i)
                    terms_tbd.append(exp)
                    expression.append(a)
                
            for i in expression:
                for j in i:
                    self.L_exp.append(j)                                          
        if terms_tbd:
            for i in terms_tbd:
                self.L_exp.remove(i)     

    #---------------------------------------------------------------------------
    #@blockPrinting
    def mapping_function(self,keyword, field, field_lorentz=None): 
        #replacement of field strength and covariant derivative 
        field_list=[]
        if keyword=='fstr':
             for key,values in self.mapping['fstr'].items():
                print(key)
                print(values)
                
                if field.name == values['fields'][0].name:
                    a = None
                    idx_list = field.get_free_indices
                    l_index=[idx_list[0][-2:],idx_list[1][-2:]]
                    if len(idx_list[0]) == 3:
                        a = [[idx_list[0][0]],[idx_list[1][0]]]
                        f=key.f
                        f.symmetry=-1
                        f.indices[1][0]=idx_list[1][0]
                    ori_field=values['fields'][1]
                    

                    if a:
                        #create fields for the non abelian field strength tensor
                        ori_field.indices=(ori_field.indices[0],\
                                        [f.get_indices[1][0],idx_list[1][-1]])
                        field1=Field(ori_field.name,
                                     ([idx_list[0][0],l_index[0][1]],
                                      [idx_list[1][0],l_index[1][1]]),
                                     latex_name=ori_field.latex_name,
                                     explicit_name=ori_field.explicit_name,
                                     field_type=ori_field.field_type,
                                     symmetry=ori_field.symmetry)
                        field2=Field(ori_field.name,ori_field.indices,
                                     latex_name=ori_field.latex_name,
                                     explicit_name=ori_field.explicit_name,
                                     field_type=ori_field.field_type,
                                     symmetry=ori_field.symmetry)
                        field3=Field(ori_field.name,
                                     ([f.get_indices[0][1],
                                       field2.indices[0][1]],\
                                       [-1*f.get_indices[1][1],\
                                       ori_field.get_indices[1][1]]),
                                     latex_name=ori_field.latex_name,
                                     explicit_name=ori_field.explicit_name,
                                     field_type=ori_field.field_type,                                     
                                     symmetry=ori_field.symmetry)
                        field4=Field(ori_field.name,([f.get_indices[0][2],
                                      field1.indices[0][1]],\
                                      [-1*f.get_indices[1][2],\
                                      ori_field.get_indices[1][1]]),
                                     latex_name=ori_field.latex_name,
                                     explicit_name=ori_field.explicit_name,
                                     field_type=ori_field.field_type,
                                     symmetry=ori_field.symmetry)
                        
                        field_del1=Del(field1, indices=([l_index[0][0]],\
                                                        [l_index[1][0]]))
                        field_del2=Del(field2, indices=([l_index[0][1]],\
                                                        [l_index[1][1]]))
                        field_3 = [ f, field3,field4]
                        temp3 = deepcopy(self.mapping['fstr'][key])
                        temp3['coupling'] *= values['fields'][2]
                        temp3['fields'] = field_3
                        temp3.update({'dirac_pair':[]})
                        field_list.append(temp3)
                    else:
                        #create field strength tesnor expansion for abelian
                        ori_field.indices=(ori_field.indices[0],\
                                          [idx_list[1][-1]])
                        field1=Field(ori_field.name,([l_index[0][1]],\
                                                     [l_index[1][1]]),
                                     latex_name=ori_field.latex_name,
                                     explicit_name=ori_field.explicit_name,
                                     field_type=ori_field.field_type,                                     
                                     symmetry=ori_field.symmetry)
                        field2=Field(ori_field.name,ori_field.indices,
                                     latex_name=ori_field.latex_name,
                                     explicit_name=ori_field.explicit_name,
                                     field_type=ori_field.field_type,
                                     symmetry=ori_field.symmetry)
                        
                        field_del1=Del(field1, indices=([l_index[0][0]],\
                                                         [l_index[1][0]]))
                        field_del2=Del(field2, indices=([l_index[0][1]],\
                                                         [l_index[1][1]]))
                    temp1 = deepcopy(self.mapping['fstr'][key])
                    temp2 = deepcopy(self.mapping['fstr'][key])
                    temp1['fields'] = [field_del1]
                    temp2['coupling'] *= -1
                    temp2['fields'] = [field_del2]
                    temp1['dirac_pair']=[]
                    temp2['dirac_pair']=[]
                    field_list.extend([temp1, temp2])
        elif keyword=='covdel':
            for key,values in self.mapping['CovDel'].items():
                #represents the first term in the covdel expansion
                if str(key)=='Lorentz':
                    temp = deepcopy(self.mapping['CovDel'][Lorentz])
                    l_index=field_lorentz.get_indices
                    field_del=Del(field, indices=l_index)
                    temp['fields']=[field_del]
                    #modifiaction of dirac pair
                    if field.field_type=='fermion':
                        temp['dirac_pair']=[[field_del.name]]
                    else:
                        temp['dirac_pair']=[]
                    field_list.append(temp)
                else:
                    #if u1 group has gauge field
                    if key.abelian==True:
                        field_idx=[]
                    #for other non abelian groups
                    else:
                        field_idx=[0]
                    g=self.GaugeFieldInstances[key]
                    gf=deepcopy(g.gauge_field)
                    gf.indices[1][-1]= field_lorentz.indices[1][0]
                    i1=field.get_free_indices[0]                   
                    i2=field.get_free_indices[1]
                    for i in range(len(i1)):
                        #T matrix is created only if non abelian group and 
                        #non singlet rep
                        if Rep.rep_dict[i1[i].index_type].group.abelian==False\
                           and  Rep.rep_dict[i1[i].index_type].singlet==True:
                            continue
                        if g.group==Rep.rep_dict[i1[i].index_type].group:
                            temp = deepcopy(self.mapping['CovDel'][key])
                            if key.abelian == False:
                                idx_list = [0]
                            else:
                                idx_list = []
                            if '_dag' in field.name:
                                temp['coupling']=-1*temp['coupling']
                            if i2[i] == -1:
                                idx_list.extend([1,1])
                            else:
                                idx_list.extend([0,0])
                            #creation of t atrix
                            t_matrix = Rep.rep_dict[i1[i].index_type].T(g, \
                                                                  idx=idx_list)  
                            t_matrix_indices=t_matrix.get_indices[0][-2:]
                            
                            
                            if key.abelian==False:
                                tm_new=[]
                                for x in t_matrix_indices:
                                    name=str(x.index)
                                    type1=Rep.rep_dict[x.index_type]
                                    name=name.replace(name[-1],str(int(name[-1])+2))
                                    tm_new.append(Index(name, type1))
                                    t_matrix.indices=(\
                                    [t_matrix.get_indices[0][0],*tm_new], \
                                    t_matrix.get_indices[1])
                                '''
                                if '_dag' or '_bar' in field.name:
                                    temp['fields'].extend([field,  t_matrix, gf])
                                else:
                                '''
                                temp['fields'].extend([gf, t_matrix, field])
                                delta=self.create_delta(\
                                [[t_matrix.get_indices[0][-2]],\
                                [t_matrix.get_indices[1][-2]*-1]],\
                                                        [[i1[i]],[i2[i]*-1]])
                                temp['delta'].extend(delta)
                                temp=self.eliminate_delta(temp)
                                temp['fields']=self.rename_dummy(temp['fields'])
                                delta2=self.create_delta(\
                                [[t_matrix.get_indices[0][-1]],\
                                [t_matrix.get_indices[1][-1]*-1]],\
                                                        [[i1[i]],[i2[i]]])
                                temp['delta'].extend(delta2)
                                temp=self.eliminate_delta(temp)
                                temp['dirac_pair']=[[]]
                                if field.field_type=='fermion':                                  
                                    temp['dirac_pair'][0].append(field.name)
                                else:
                                    temp['dirac_pair']=[]
                            else:
                                temp['fields'].extend([gf, t_matrix, field])
                                temp['coupling']*=field.charge
                                
                                if field.field_type=='fermion':
                                    temp['dirac_pair']=[[field.name]]
                                else:
                                    temp['dirac_pair']=[]

                            field_list.append(temp) 
        return(field_list)


    #---------------------------------------------------------------------------
    def field_multiplication(self, field1, field2): 

        #mutliplication of 2 terms 

        tlist1_indices,tlist2_indices=[[],[]],[[],[]]
        tlist1 = [x['fields'] for x in field1]
        tlist2 = [x['fields'] for x in field2]
        #before the multiplication remane all the dummy indices.
        #check if there is any same dummy indices
        for y in range(len(tlist1)):
            for x in tlist1[y]:
                tlist1_indices[0].extend(x.get_indices[0])
                tlist1_indices[1].extend(x.get_indices[1])
        for y in range(len(tlist2)):
            for x in tlist2[y]:
                tlist2_indices[0].extend(x.get_indices[0])    
                tlist2_indices[1].extend(x.get_indices[1])
        all_idx1=[tlist1_indices[0][x].index*tlist1_indices[1][x] \
                  for x in range(len(tlist1_indices[0]))]
        all_idx2=[tlist2_indices[0][x].index*tlist2_indices[1][x] \
                  for x in range(len(tlist2_indices[0]))]
        #finding the common indices in the terms
        common_indices=[x for x in all_idx1 if x in all_idx2]
        #if common index is there replace the index with dummy index
        if len(common_indices)!=0:
            tlist2new=[]
            for x in tlist2:
                tlist2new.append(self.rename_dummy(x))
        else:
            tlist2new=tlist2
        #multiply the fields in the term
        a = self.list_multiplication(tlist1,tlist2new)
        field_return = []
        k=0
        for i in field1: 
            for j in field2:
                #modification of the dirac pair
                dirac_pair=deepcopy(i['dirac_pair'] )
                for x in range(len(j['dirac_pair'])):
                    dirac_pair[x].extend(j['dirac_pair'][x])
 

                
                temp_field = {'coupling': i['coupling']*j['coupling'],
                                'delta': [], 'fields': a[k],
                                'dirac_pair': dirac_pair}
                field_return.append(temp_field)
                dirac_pair=[]
                k += 1
        return(field_return)
                       
    #---------------------------------------------------------------------------
    @staticmethod
    def list_multiplication(list1, list2):
        #mutliply 2 set of fields--follows distributive laws
        field_list = []
        for i in list1:
            for j in list2:
                temp_field_list = [*i]
                temp_field_list.extend(j)
                field_list.append(temp_field_list)
        return(field_list) 
    #---------------------------------------------------------------------------
    def rename_dummy(self, fields):
        #Take the product of field to get indices in single step
        Fields=FieldMul(*fields)
        free_idx=Fields.get_free_indices
        all_idx=Fields.get_indices
        temp_list=[]
        
        new_list=[]
        #compare free indices and all indices
        if len(free_idx[0])!=len(all_idx[0]):
            index_list=all_idx[0]
            index_list_sym=[str(x.index) for x in index_list]
            #modify all the dummy indices in the list
            new_index_list=[self.rename_index(x,temp_list)\
                            if index_list_sym.count(str(x.index))==2 else x\
                              for x in index_list]          

            k=0
            l=0
            #create field with the new indices created
            for f in fields:
                l=k+len(f.get_indices[0])
                idx=new_index_list[k:l]
                g=Field(f.name,(idx,f.get_indices[1]),
                    latex_name=f.latex_name,
                    explicit_name=f.explicit_name,
                    field_type=f.field_type,
                    symmetry=f.symmetry,
                    numerical_values=f.numerical_values)
              
                k=l

                new_list.append(g)
        else:
            new_list=[*fields]
        return(new_list)

    #---------------------------------------------------------------------------
    def rename_index(self,index,temp_list):
        #rename the index passed as agrument
        type1=index.index_type
        if index.index not in temp_list:
            temp_list.append(index.index)
            try:
                self.rename_index_dict[type1] += 1
                #name of the index is the capital letter of the index nae of 
                #representations
                name=Rep.rep_dict[type1].index_name.upper()+ \
                                        str(self.rename_index_dict[type1])
                self.n_dict[type1].append(Index(name, Rep.rep_dict[type1]))
            except KeyError:
                self.rename_index_dict[type1] = 1
                name=Rep.rep_dict[type1].index_name.upper()+ \
                                        str(self.rename_index_dict[type1])
                self.n_dict[type1] = [Index(name, Rep.rep_dict[type1])]

            return(self.n_dict[type1][-1])
        else:
            temp_list.remove(index.index)
            return(self.n_dict[type1].pop(0))
    #---------------------------------------------------------------------------
    def MakeExplicit(self, list1, group, mode='cartesian'):
        group_names=[x.name for x in group]
        #print(group_names)
        final_list=[]
        adjoint_reps=[x.adjoint_rep.name for x in group]
        for term in list1:
            explicit_indices_non_adj={}
            explicit_indices_adj={}
            for f in term['fields']:
                indices=f.get_indices
                for idx in indices[0]:
                    if (Rep.rep_dict[idx.index_type].group.name in group_names)\
                    and (Rep.rep_dict[idx.index_type].singlet==False):
                        if Index.idx_dict[idx.index].index_type in adjoint_reps:
                            try:
                                explicit_indices_adj[idx.index].\
                                append(term['fields'].index(f))
                            except KeyError:
                                explicit_indices_adj.update({idx.index:\
                                [term['fields'].index(f)]})
                        else:
                            try:
                                explicit_indices_non_adj[idx.index].\
                                append(term['fields'].index(f))
                            except KeyError:
                                explicit_indices_non_adj.update({idx.index:\
                                [term['fields'].index(f)]})
            #print(explicit_indices_non_adj, explicit_indices_adj)
            #fixing the mode of explicit function
            LA=False
            Numerical=False
            explicit_indices={}
            if mode=='cartesian':
                explicit_indices.update(explicit_indices_adj)
                explicit_indices.update(explicit_indices_non_adj)
                
            elif mode=='adjoint_only':
                #expand adjoint indices in cartan weyl
                explicit_indices.update(explicit_indices_adj)
                LA=True
            elif mode=='cartan_weyl':
                explicit_indices.update(explicit_indices_adj)
                explicit_indices.update(explicit_indices_non_adj)
                LA=True
            elif mode=='cartan_weyl_with_numerical':
                explicit_indices.update(explicit_indices_adj)
                explicit_indices.update(explicit_indices_non_adj)
                LA=True
                Numerical=True

            temp_term_list=[deepcopy(term)]  #I need this line
            for key in explicit_indices.keys():
                #---check if the index is adjoint or not-----
                if Index.idx_dict[key].index_type in adjoint_reps and mode!='cartesian':
                    adjoint=True
                else: 
                    adjoint=False
                new_temp_list=[]
                for t in temp_term_list:
                    pos1, pos2 = explicit_indices[key][0], \
                                 explicit_indices[key][1]
                    new_fields_list=self.make_components(t['fields'][pos1],\
                                    t['fields'][pos2], Index.idx_dict[key], adjoint,LA,Numerical)
                    for x in range(len(new_fields_list[0])):
                        temp_term=deepcopy(t)
                        if isinstance(new_fields_list[0][x], Field) or isinstance(new_fields_list[0][x], FieldMul): 
                            temp_term['fields'][pos1]=new_fields_list[0][x]
                        else:
                            temp_term['coupling']*=new_fields_list[0][x]
                            temp_term['fields'][pos1]=None
                        if isinstance(new_fields_list[1][x], Field) or isinstance(new_fields_list[1][x], FieldMul):
                            temp_term['fields'][pos2]=new_fields_list[1][x]
                        else:
                            temp_term['coupling']*=new_fields_list[1][x]
                            temp_term['fields'][pos2]=None
                        if temp_term['dirac_pair']!=[]:
                            for y in temp_term['dirac_pair']:
                                if term['fields'][pos1].name in y:
                                    f_pos=y.index(term['fields'][pos1].name)
                                    if isinstance(new_fields_list[0][x], Field) or isinstance(new_fields_list[0][x], FieldMul):
                                        y[f_pos]=new_fields_list[0][x].name
                                    else:
                                        y[f_pos]=None
                                if term['fields'][pos2].name in y:
                                    f_pos2=y.index(term['fields'][pos2].name)
                                    if isinstance(new_fields_list[1][x], Field) or isinstance(new_fields_list[1][x], FieldMul):
                                        y[f_pos2]=new_fields_list[1][x].name
                                    else:
                                        y[f_pos2]=None
                        new_temp_list.append(temp_term)
                        #new_temp_list.append(temp_term)
                
                temp_term_list=new_temp_list
            final_list.extend(temp_term_list)
        
        final_list=self.cleaning_up_after_explicit(final_list)
        return(final_list)
                
            

        
        

                
    #==============================================================================
    def cleaning_up_after_explicit(self, list1):
        new_list=[]
        for x in list1:
            dict1={}
            dict1.update({'coupling': x['coupling'], 
                          'delta':x['delta']})
            field_list=[]
            for field in x['fields']:
                if (field !=None):
                    field_list.append(field)

            dict1.update({'fields':field_list})
            dict1.update({'dirac_pair':x['dirac_pair']})
            #dirac_pair=[]    
            new_list.append(dict1)
        return(new_list)                                  
           
    #@blockPrinting
    #==============================================================================
    def make_components(self, field1, field2, index, adjoint=False, useLA=True, numerical=True):
        #print ('start processing make_components')
        #print('mc1', field1, field2, index)
        
        
        #--removing first component if field is multicomponent------------------

        #multicomponent flag
        mc1=0
        mc2=0
        if isinstance(field1, FieldMul):
            f1=field1.fields[1]
            mc1=1
        else:
            f1=field1
        if isinstance(field2, FieldMul):
            f2=field2.fields[1]
            mc2=1
        else:
            f2=field2
        #print('mc2', f1, f2)
    
        #---make pseudo singlet indices for f1 and f2---------------------------
        new_indices1=self.make_pseudo_singlet(f1, index)
        new_indices2=self.make_pseudo_singlet(f2, index)
        #print('indices', new_indices1, new_indices2)

        #-----name of the components--------------------------------------------
        components1=self.create_components_name(f1, index, adjoint, useLA, numerical)
        components2=self.create_components_name(f2, index, adjoint, useLA,  numerical)
        #print('comp_names', components1, components2)
        
        #-------for No LA and antisymmetric tensors---------------------------------
        #print('anti symmetry property', f1.name, f2.name, f1.symmetry, f2.symmetry)
        if adjoint==False:
            #print('yes')
            if f1.symmetry==-1:
                #print('2yes')
                temp_comp1=deepcopy(components1)
                for x in temp_comp1:
                    partname=x.split('_x_')
                    number_string= [i for i in \
                                partname[1].split('_') if i.isalpha()==False] 
                    if len(list(set(number_string)))< len(number_string):
                        pos_temp=temp_comp1.index(x)
                        for y in components1:
                            if y==x:
                                pos_temp2=components1.index(y)
                                components1.remove(y)
 
                                components2.remove(components2[pos_temp2])  
                                           

            if f2.symmetry==-1:
                temp_comp2=deepcopy(components2)
     
                for x in temp_comp2:
                    partname=x.split('_x_')
                    number_string= [i for i in \
                                partname[1].split('_') if i.isalpha()==False] 
                    if len(list(set(number_string)))< len(number_string):
                        pos_temp=temp_comp2.index(x)
                        for y in components2:
                            if y==x:
                                pos_temp2=components2.index(y)
                                components2.remove(y)
 
                                components2.remove(components1[pos_temp2])
        
        #-------flip components for adjoint representation----------------------
        if adjoint==True:
            new_field_list1=self.swap_components(components1, index)
        #-----remove unneccessary components-----------------------------------
        new_components1, new_components2=[], []
        for x in range(len(components1)):
            if components1[x]==None or components2[x]==None:
                pass
            else:
                new_components1.append(components1[x])
                new_components2.append(components2[x])
                
        #print(new_components1, new_components2)
        #-------create latex name of the components----------------------------
        if new_components1!=[]:
            latex1=self.create_latex_name_components(f1, new_components1)
            latex2=self.create_latex_name_components(f2, new_components2)
        
        #----create new field with the components------------------------------
        new_field_list1, new_field_list2=[], []
        for x in range(len(new_components1)):
            
            x1=new_components1[x]            
            x2=new_components2[x]
            if isinstance(x1, str):
                new_field1=Field(x1, new_indices1,
                             field_type=f1.field_type,
                             latex_name=latex1[x],
                             explicit_name=x1,
                             numerical_values=f1.numerical_values,
                             symmetry=f1.symmetry,
                             base=f1.base)
            else:
                new_field1=x1
            if isinstance(x2, str):
                new_field2=Field(x2, new_indices2,
                             field_type=f2.field_type,
                             latex_name=latex2[x],
                             explicit_name=x2,
                             numerical_values=f2.numerical_values,
                             symmetry=f2.symmetry,
                             base=f2.base)
            else:
                new_field2=x2
            #--------adding the untouched field in multicomponent fields------
            if mc1==1:
                new_field_list1.append(FieldMul(field1.fields[0], new_field1))
            else:
                new_field_list1.append(new_field1)
            if mc2==1:
                new_field_list2.append(FieldMul(field2.fields[0], new_field2))
            else:
                new_field_list2.append(new_field2)
        
        #print(new_field_list1, new_field_list2)
        #print('end processing make components')
        return(new_field_list1, new_field_list2)
    #==========================================================================
    def create_latex_name_components(self, field, components_list):
       
        #print('start processing create_latex_name_components')
        #print(field, components_list)
        if isinstance(components_list[0], str)==False:
             new_latex_name_list=[None for x in range(len(components_list))]
             return(new_latex_name_list)
        #print(field, components_list)
        base_latex_name=field.base.latex_name
        new_latex_name_list=[]
        if field.field_type=='constant':            
             latex_name='('+field.base.latex_name+')'
        else:
            latex_name=field.base.latex_name
        for x in components_list:
               
                sp=x.split('_x_')
                sp2=sp[1].replace('_','')
                
                if 'c' or 'p' or 'm' in sp2:
                    #print('yes')
                    if 'p' in sp2:
                        #sp2=sp2.replace('p', 'p_')
                        sp2=sp2.replace('p', '+')
                    if 'm' in sp2:
                        #sp2=sp2.replace('m', 'm_')
                        sp2=sp2.replace('m', '-')
                    if 'c' in sp2:
                        #sp2=sp2.replace('c', 'c_')
                        sp2=sp2.replace('c', ' ')
                
                new_latex_name_list.append('('+latex_name+'_{'+sp2+'})')
        #print('end processing create_latex_name_components')
        return(new_latex_name_list)
    #==========================================================================
    def swap_components(self, list1, index):  #can be moved to the utilities
        '''
        the order of the states and operator in the field is 
        heighest to lowest.
        i.e. p1, p2,..., c1, c2 ,...,...m3, m2, m1
        This function interchanges p_i and m_i 's
        i.e. p1<->m1, p2<->m2 
        This is applicable only for the adjoint representation indices.
        '''
        rep=Rep.rep_dict[index.index_type]
        dim =rep.dim
        Nroots=rep.group.Nroots
        for x in range(Nroots):
            list1[x], list1[-(x+1)]=list1[-(x+1)], list1[x]
        #print('end processing swap_components')
        return(list1)
        
    #==========================================================================
    def create_components_name(self, field, index, adjoint=False,useLA=True, numerical=True):
        #print('start processing create components name')
        #print('ccn1', field, index, adjoint)
        e_name=[]
        numerical_value=[]
        flag=0
        if field.field_type=='constant' and field.numerical_values!=None and useLA==True:
            flag=1
            #print('values' , field.numerical_values)
        #print(flag)
        explicit_name=field.explicit_name
        rep=Rep.rep_dict[index.index_type]
        dim =rep.dim
        names=explicit_name.split('_x')
        #print(names)
        #print('flag', flag)
        
        for y in range(dim):
            if adjoint==False:
                name2=names[1].replace('_' + rep.index_name, '_'+str(y+1), 1)
            else:
                nroots=rep.group.Nroots
                if (y<nroots):
                    n='p'+str(y+1)
                elif dim-nroots<=y<dim:
                    n='m'+str(dim-y)
                else:
                    n='c'+ str(y-nroots+1)
                name2=names[1].replace('_'+rep.index_name, '_'+n, 1)
            #print('name2', name2)
            flag2=0
            if flag==1:
                onlyno=name2.split('_')
                onlyno.remove('')
                #print('onlyno', onlyno)
                possible=check_if_value_available(field.numerical_values, onlyno, numerical)
                if possible[0]==1:
                    if possible[1]!=None:
                        e_name.append(possible[1])
                    else:
                        e_name.append(names[0]+'_x'+name2)
                    flag2=1
                else:
                    e_name.append(None)
                    flag2=1
            if flag2==0:
                e_name.append(names[0]+'_x'+name2)
        #print(e_name)
        #print('end processing create components name')          
        return(e_name)

    #===========================================================================
    
               
      
    #============================================================================
    def make_pseudo_singlet(self, field,index, adjoint=False):  
        '''
        take one field tensor and index as arguments
        inside create the pseudo rep if not already created
        create the components name and latex name
        create new field using the components name and features 
        of the old field
        return new field
        
        '''             
        #print('start processing make_pseudo_singlet')
        index_type=index.index_type        
        indices=field.get_indices
        #print(index_type, indices)
        #list of reps in the field
        theory = [Rep.rep_dict[indices[0][k].index_type] \
                    for k in range(len(indices[0]))]
        #create pseudo rep
        if Rep.rep_dict[index_type].singlet==False:
            if not Rep.rep_dict[index_type].pseudo_singlet:
                Rep.rep_dict[index_type].pseudo_singlet = Rep('Ps_' + \
                                            str(index_type),Rep.rep_dict\
                                            [str(index_type)].\
                                            group, singlet=True)
        #modify the reps with pseudo rep of the index
        theory_new = [theory[i] if (str(indices[0][i]) !=str(index) ) \
                                else theory[i].pseudo_singlet \
                                for i in range(len(theory))] 
        new_index_list=[[],[]]
        #creating the indices of the new field
        for x in range(len(indices[0])):
            if str(indices[0][x]) ==str(index):
                if indices[1][x]==1:
                    new_idx=Index('Ps_' + str(indices[0][x].index),
                                        theory_new[x])
                    new_index_list[0].append(new_idx)
                    new_index_list[1].append(1)
                else:
                    name='Ps_'+str(indices[0][x].index)
                    new_idx=Index(name,theory_new[x])
                    new_index_list[0].append(new_idx)
                    new_index_list[1].append(-1)
            else:
                new_index_list[0].append(indices[0][x])
                new_index_list[1].append(indices[1][x])
        #print('end processing make_pseudo_singlet')
        
        return(new_index_list)
        
    #------------------------------------------------------------------------------------
    def GiveVev( self,lag,replace_field, vev_value, mode='cartan_weyl_with_numerical'):    
        #only scalar field and the VEV as a tuple is allowed as arguments                                
        if replace_field.properties['field_type'] != 'scalar':
            raise Exception ("Only Scalars allowed!")
        tb_replaced = replace_field.field
        tb_replaced_dag = replace_field.dag
               
        #__________exlicit expansion________________________________
        theory = [Rep.rep_dict[tb_replaced.get_indices[0][k].index_type]
                         for k in range(len(tb_replaced.get_indices[0]))]
        temp = [0 if theory[i].singlet else 1 for i in range(len(theory))]
        temp_explicit = [theory[i].group for i in \
                         range(len(temp)) if temp[i] == 1]  
        Lag=self.MakeExplicit(lag,temp_explicit, mode=mode)
        #__________check if tuple is passed__________
        if not isinstance(vev_value,tuple):
            raise Exception("vev should be given as a tuple")
        #__________dimension of the tuple_________________
        dimension=1
        for d in range(len(theory)):
            if theory[d].singlet:
                pass
            else:
                dimension*=theory[d].dim
        if len(vev_value)!=dimension:
            raise Exception("Wrong dimension of the tuple")
        #______finding the non zero components of the tuple______
        non_brok_ind=[]
        broken_value=[]
        broken_indices=[]
        for a in range(len(vev_value)):
            if vev_value[a]!=0:
                broken_indices.append(a)
                broken_value.append(vev_value[a])
            else:
                non_brok_ind.append(a)

        name_list1,name_list2=[],[]
        
        #_____creating the pseudo singlet names___________
        if temp[0]!=0:
            components=['_x_'+str(i+1) for i in range(theory[0].dim)]
        else:
            components=['_x_s']
        if len(theory)>1:
            for i in range(1,len(theory)):                      
                small_list=components
                components=[]
                new_list=[]
                for m in range(len(small_list)):
                    if temp[i]!=0:                                                                                 
                        for j in range(theory[i].dim):
                            new_list.append(small_list[m]+'_'+str(j+1))
                    else:
                        new_list.append(small_list[m]+'_s')
                components.extend(new_list)     
        for i in broken_indices:
            name_list1.append(components[i])
        for i in non_brok_ind:
            name_list2.append(tb_replaced.name+components[i])
            name_list2.append(tb_replaced_dag.name +components[i])
        final ,rest= [],[]
        if temp_explicit==[]:
            name_list1=[tb_replaced.name,tb_replaced_dag.name]
            name_list2=[]
        #________Checking each term and removing terms end up giving zero______
        for i in Lag:
            for j in i['fields']:
                if isinstance(j,FieldMul):
                    f=j.fields[1]
                    if f.name in name_list2:                               
                        rest.append(i)
                        break
                    else:
                        continue                                                       
                elif not isinstance(j,FieldMul) and j.field_type!='constant':
                    if j.name in name_list2:
                        rest.append(i)
                        break 
                    else:
                        continue                            
                else:#for constant
                    continue
        for m in rest:
            Lag.remove(m)
        for n in Lag:
            final.extend(self._give_vev_mapping(n,replace_field,\
                         name_list1,broken_value))
        return(final)
    #---------------------------------------------------------------------------
    #@blockPrinting
    def _give_vev_mapping(self, list1, replace_field, name_list,value_list):
        '''
        take the scalar field and its dagger , components name and vev 
        value as arguments
        for each term of the Lag expr we replace the scalar field as V + H 
        use list multiplication to get final expression
        '''
        field=replace_field.field
        field_dag=replace_field.dag
        if replace_field.properties['self_adjoint']==False:
            norm=Symbol('1/sqrt(2)',positive=True)
        else:
            norm=Symbol('1',positive=True)
        value_list_sym=[]
        for y in value_list:
            Y=sympify(y)
            value_list_sym.append(Y)

        new_list=[]
        if name_list==[field.name,field_dag.name]:
            new_list=name_list
        else:
            for x in name_list:
                new_list.append(field.name+x)
                new_list.append(field_dag.name+x)
        
        final_list=[]
        fields=[]
        return_list=[]
        #__________replacement of the field with its component____________
        for i in list1['fields']:
            if isinstance(i,FieldMul):                 
                f=i.fields[1]               
                if f.name in new_list:
                    pos=new_list.index(f.name)
                    h = self._generate_vev_components(f,new_list[pos],\
                           value_list[int(pos/2)])
                    v=value_list_sym[int(pos/2)]
                    h1=FieldMul(i.fields[0],h)
                    final_list.append([[h1]])
                    final_list.append([[norm]])
                else:
                    final_list.append([[i]])                                             
            elif i.field_type != 'constant':
                if i.name in new_list:
                    pos=new_list.index(i.name)
                    h = self._generate_vev_components(i, new_list[pos],\
                          value_list[int(pos/2)])
                    v=value_list_sym[int(pos/2)]
                    final_list.append([[v],[h]])
                    final_list.append([[norm]])
                else:
                    final_list.append([[i]])
            else:
                    final_list.append([[i]])
        #multiplication of the fields after the replacement
        for m in final_list: 
            if fields: 
                fields=self.list_multiplication(fields,m)
            else:
                fields=m
        for i in fields:            
            temp_dict = deepcopy(list1)           
            t_list = []
            for j in i:                
                if isinstance(j,Symbol):
                    temp_dict['coupling'] *= j
                else:
                    t_list.append(j)
            temp_dict['fields'] = t_list            
            return_list.append(temp_dict)
        return(return_list)    

    #------------------------------------------------------------------------------------
    def _generate_vev_components(self, field1,name1 ,value1): 
        '''
        Create the higgs field and (the vev value as symbol)
        '''
        #v = Symbol(value1, positive=True)
        index_list = field1.indices[0]
        latex_name=field1.latex_name
        
        #creating the name for Higgs field
        if '_dag' in name1:
            name1=name1.replace('_dag','')
        if '_' not in field1.latex_name:
             latex_name='H'
        else:
             latex_index=(field1.latex_name.split('_'))[1]
             latex_name= 'H_'+latex_index[:-1]
        rep_name_list = [Rep.rep_dict[i.index_type] for i in index_list]

        #create indices
        #indices=Create_indices(rep_name_list)
        #creating the field
        h= Field('H' + name1, field1.get_indices, field_type=field1.field_type,
                              explicit_name=field1.explicit_name,
                              latex_name=latex_name,
                              symmetry=0,
                              numerical_values=None)
        
        
        return(h)

    #---------------------------------------------------------------------------
    def Simplify(self,lag):
        #simplification of certain terms of the Lagrangian
        new_lag=[]
        
        name_list=[]
        #name grouping of each term in the Lagrangian
        for term in lag:          
            name_list.append(self.name_grouping(term))

        #sort terms of the lagrangian based on having same set of fields
        similar_term_index=group_same_strings(name_list)

        #each of the sorted list is passed for finding mapping and doing
        #permutations
        for g in similar_term_index:
            similar_terms=[]
            #only one term with particular set of fields--no need of further 
            #simplification                                                      #name_grouping is working fine
            if len(g)==1:
                new_lag.append(lag[g[0]])
                continue
            for h in g: 
                similar_terms.append(lag[h])
            #print('1')
            #for x in similar_terms:
             #   print(x)
            new_lag.extend(self.sort_same_mapping(similar_terms))
           
        return(new_lag)
    

    #---------------------------------------------------------------------------
    def name_grouping(self,term): 
        #count the fields in a term and out it as dictionary with key 'name of 
        # field and value its count
        bal_tensor=[y.name for y in self.tensor_names]
        name={}
        for field in term['fields']:
            #for balancing tensor, consider its latex name.
            if field.name in bal_tensor:
                try:
                    name[field.latex_name]+=1
                except:
                    name.update({field.latex_name:1})
            else:
                try:
                    name[field.name]+=1
                except:
                    name.update({field.name:1})
        #sorting is important as we compare it with other terms
        sorted_keys=sorted(name.keys())
        #convert the dict into a string for managig easily    
        string=''
        for key in sorted_keys:
            string+=str(key)+':'+ str(name[key])+','
        return(string[:-1])
    #---------------------------------------------------------------------------
    def dirac_p_same_mapping(self, lag):
        #print('inside dirac_p_permutation')
        #for doing dirac permuatation for terms having same mapping
        #fix the first term and compare the dirac pairing of the rest of the 
        #terms         
        first=lag[0]
        for x in range(1, len(lag)):
            second=lag[x]
            new_second=self._dirac_permutation(first, second)
            first['coupling']+=new_second['coupling']
            first['coupling']=simplify(first['coupling'])
        if first['coupling']!=simplify('0'):
            return(first)
        else:
            return(None)
    #---------------------------------------------------------------------------
    def indices_mapping(self, term): 
        #create mapping of each term in the Lagrangian
        term_mapping=[]
        fields_to_permute=[]
        index_to_permute=[]

        field_count={}
        index={}
        
        
        for x in range(len(term['fields'])):
            #count the fields in the term
            field_name=term['fields'][x].name
            try:
                field_count[field_name].append(x)
            except:
                field_count.update({field_name:[x]})

        for key in field_count:
            #if the count of one field is greater than 1 name the field with 
            #a number
            if len(field_count[key])>1:
                list1=[]
                for k in range(len(field_count[key])):
                    name1=str(key)+str(k+1)
                    list1.append(name1)
                #list of fields to permute
                fields_to_permute.append(list1)

        #counting the indices of each field in the term
        for field in term['fields']:
            field_name=field.name
            if len(field_count[field.name])>1:
                position_term=term['fields'].index(field)
                position_count=field_count[field_name].index(position_term)
                field_name=field_name+str(position_count +1)

            
            indices=field.get_indices[0]
            field_index_count={}
            
            for x in indices:
                idx_type=x.index_type                   
                try:
                    field_index_count[idx_type].append(x)
                except:
                    field_index_count.update({idx_type:[x]})
           
            temp_name=field_name
            
           
            for key in field_index_count:
                if len(field_index_count[key])>1:
                    #if the field has symmetry and anti symmetry in indices
                    #require index permutations
                    if (field.symmetry==1) or (field.symmetry==-1):
                        list2=[]
                        for k in range(len(field_index_count[key])):
                            name2=str(temp_name)+'_'+str(k+1)
                            list2.append(name2)
                        #indices to permute(field_name+_number)
                        index_to_permute.append(list2)

            #create a dict with key 'index' and values 'fields with that index'
            for x in indices:
                
                idx_type=x.index_type
                
                if len(field_index_count[idx_type])>1:
                    pos=field_index_count[idx_type].index(x)
                    field_name=temp_name+'_'+str(pos +1)
                    
                if Rep.rep_dict[idx_type].singlet==False:    
                    try:
                        index[x.index].append(field_name)
                    except:
                        index.update({x.index:[field_name]})
                
        #creating the mapping
        for key in index:
            index_type=Index.idx_dict[key].index_type
            if len(index[key])==2:
                str1=index[key][0] +'-'+index_type +'-' +index[key][1] 
                str2=index[key][1] +'-'+index_type +'-' +index[key][0]
                    
            elif len(index[key])==1:
                str1=index[key][0] + '-' + str(key)
                str2=str(key)+'-'+index[key][0]
            else:
                raise Exception('incorrect contraction')
            term_mapping.append(str1)
            term_mapping.append(str2) 
        
        term_mapping=sorted(term_mapping)
        str_term_mapping=str(term_mapping)[1:-1]
        return(term_mapping, sorted(fields_to_permute), sorted(index_to_permute))
        

    #---------------------------------------------------------------------------
    def sort_same_mapping(self, lag):
        #terms with same set of fields are passed here
        #print('2')
        #for x in lag:
         #   print(x)
        new_lag=[]
        mapping_list=[]
        mapping_str=[]
        fields_to_permute_list=[]
        index_to_permute_list=[]
        #for each term in the lagrangian
        for term in lag:
            #print('2', lag.index(term))
            #find the mapping
            indices_mapping=self.indices_mapping(term)
            #convert the mapping to string for easy management
            mapping_str.append(str(indices_mapping[0]))
            mapping_list.append(indices_mapping[0])
            fields_to_permute_list.append(indices_mapping[1])
            index_to_permute_list.append(indices_mapping[2])
        index_set=group_same_strings(mapping_str)
        fields_list=fields_to_permute_list[0]
        index_list=index_to_permute_list[0]
        new_mapping_list=[]
        #for terms with same mapping do dirac permutations
        for x in index_set:
            similar_terms=[]
            new_mapping_list.append(mapping_list[x[0]])
            for y in x:
                similar_terms.append(lag[y])
            dirac_term=self.dirac_p_same_mapping(similar_terms)
            if dirac_term!=None:
                new_lag.append(dirac_term)

        lag_after_per=[]
        #perform field and index permutations
        if len(new_lag)>1:
            if fields_list!=[] or index_list!=[]:
                while new_lag!=[]:
                    first=new_lag[0]
                    first_mapping=new_mapping_list[0]
                    temp_new_lag=deepcopy(new_lag)
                    for x in range(1, len(new_lag)):
                        second=new_lag[x]
                        second_mapping=new_mapping_list[x]
                        result=self._permutations(first, second, first_mapping, 
                             second_mapping, index_list, fields_list)
                        if result==1:
                             new_second=self._dirac_permutation(first,second)
                             new_lag.remove(second)
                             first['coupling']+=new_second['coupling']
                             first['coupling']=simplify(first['coupling'])
                             new_mapping_list.remove(second_mapping)
                    new_lag.remove(first)
                    new_mapping_list.remove(first_mapping)
                    if first['coupling']!=sympify('0'):
                        lag_after_per.append(first)
                   
                new_lag=lag_after_per
        return(new_lag)

    #-----------------------------------------------------------------------
    def _dirac_permutation(self, first, second):
        #take 2 terms and permuta the second one w.r.t first one
        f=0
        b=0
        bosonic1,bosonic2=[],[]
        fermionic1=()
        fermionic2=()       
        dirac1=first['dirac_pair']
        dirac2=second['dirac_pair']
        for x in dirac1:
            if isinstance(x, list):
                bosonic1.append(x)
            else:
                fermionic1=x

        for x in dirac2:
            if isinstance(x, list):
                bosonic2.append(x)
            else:
                fermionic2=x
        sign=1
        #checking both terms have same bosonic pairs
        if bosonic1!=[]:
            sort1=sorted(bosonic1)
            sort2=sorted(bosonic2)
            
            if sort1==sort2:
               b=1
        elif bosonic1==[]:
            b=1
        if len(fermionic1)==0:
            f=1
        elif len(fermionic1)!=0:
            #if free fermions' order is same
            if fermionic1==fermionic2:
                f=1
            else:
                number=[]
                #if the orer of free fermion doesn't match do permutations
                for x in range(len(fermionic1)):
                    number.append(x+1)
                perm=list(permutations(fermionic2))
                num=list(permutations(tuple(number)))
                for x in range(len(perm)):
                   
                    if perm[x] == fermionic1:
                        f=1
                        sign=sign_dict[num[x]]
                        break
        if (b ==1) and (f==1):
            new_second=deepcopy(second)
            new_second['coupling']*=sign 
            return (new_second) 
        else:
            return(second)          
       

    #---------------------------------------------------------------------------
    def _permutations(self, term1, term2, mapping1, mapping2, 
                       permuting_indices, permuting_fields):
        #input 2 terms their mapping, fields topermute , indices to permute
        new_lag=[]
        cases1=[]
        sign_list=[]
        for x in permuting_indices:
            sign=[]
            permute=list(permutations(tuple(x)))
            for y in permute:
                sign.append(Sign(y))
            sign_list.append(tuple(sign))
            cases1.append(tuple(permute))
        case=list(product(*cases1))
        sign_tuple=list(product(*sign_list))
        cases_sign=[]
        for x in sign_tuple:
            y=array(x)
            cases_sign.append(prod(y))
        
        cases=[]       
        for x in case:
            one_case=[]
            for  y in x:
                one_case.extend(list(y))
            cases.append(one_case)
        case1_fields=[]
        #creating permutations of all the set of indices needed to be permutated
        for x in permuting_fields:
            fields_permuted=list(permutations(tuple(x)))           
            case1_fields.append(tuple(fields_permuted))
        case_fields=list(product(*case1_fields))
        final_case_fields=[]
        for x in case_fields:
            one_case=[]
            for y in x:
                one_case.extend(list(y))
            final_case_fields.append(one_case)
        flag=0  
        #3 cases needed to be considered 1) only index permutation
        if permuting_fields==[]:
            result=self._index_permutation( mapping1, mapping2, cases)
            flag=result[0]
            if flag==1:
                sign=cases_sign[result[1]]
        else:        
            #2) only field permutation  3) both                            
            for a in range(len(final_case_fields)):
                temp_mapping1=deepcopy(mapping2)
                if a ==0:     
                    #for the first case in field permuatation, perform index
                    #permutation            
                    if permuting_indices!=[]:
                        result=self._index_permutation(mapping1, temp_mapping1, cases)
                        flag=result[0]
                        if flag==1:
                            sign=cases_sign[result[1]]
                            break
                else: 
                    #go to next case and perform index permutation                               
                    q1=final_case_fields[0]
                    q2=final_case_fields[a]
                    dict1f={}
                    dict2f={}
                    for b in range(len(q1)):
                        dict1f.update({q1[b]:q2[b].swapcase()})
                        dict2f.update({q2[b].swapcase(): q2[b]}) 
                    for key in dict1f:   
                        new_mapping=[]            
                        for z in temp_mapping1:                   
                            z=z.replace(key, dict1f[key])
                            new_mapping.append(z)
                        temp_mapping1=new_mapping
                    for key in dict2f:   
                        new_mapping=[]            
                        for z in temp_mapping1:
                            z=z.replace(key, dict2f[key])
                            new_mapping.append(z)
                        temp_mapping1=new_mapping
                    temp_mapping1=sorted(temp_mapping1)
                    if str(mapping1)==str(temp_mapping1):
                        flag=1
                        sign=1
                        break
                    if permuting_indices!=[]:
                        result=self._index_permutation(mapping1, temp_mapping1, cases)
                        flag=result[0]
                        if flag==1:
                            sign=cases_sign[result[1]]
                            break
        #depending on the permutation changing the sign of coupling   
        if flag==1:
            term2['coupling']*=sympify(sign)
            return(1)
        else:
            return(0)

    #----------------------------------------------------------------------------------------------
    def _index_permutation(self, mapping1, mapping2, cases):
        #input mapping 1 and 2 , cases of permutations
        flag=0
        for x in range(1,len(cases)):
            p1=cases[0]
            p2=cases[x]
            temp_mapping=deepcopy(mapping2)
            dict1={}
            dict2={}
            for y in range(len(p1)):
                dict1.update({p1[y]:p2[y].swapcase()})
                dict2.update({p2[y].swapcase(): p2[y]})
            for key in dict1:   
                new_mapping=[]            
                for z in temp_mapping:
                    z=z.replace(key, dict1[key])
                    new_mapping.append(z)
                temp_mapping=new_mapping
            for key in dict2:   
                new_mapping=[]            
                for z in temp_mapping:
                    z=z.replace(key, dict2[key])
                    new_mapping.append(z)
                temp_mapping=new_mapping
                         
            temp_mapping=sorted(temp_mapping)
            if str(mapping1)==str(temp_mapping):
                 flag=1
                 
                 break
        return(flag, x)
    

#------------------------------------------------------------------------------------
    def Replace3(self, field1, field2):
        #field1 is the field to be replaced
        #field2 is the list of fields is using to replace field1
        #we need to make sure that the indices of the bothe field1 and the field 
        #in list field2 matches 
        #first step is to check whether the fields of both field1 and fielsds in field2
        #matches 
        new_lag=[] #new list to store the modified lagrangian
        indices1= field1.field.get_indices[0]
        #print('r1', indices1)
        rep1=[x.index_type for x in indices1]
        singlet1=[Rep.rep_dict[x].singlet for x in rep1]
        #print('r2', rep1, singlet1)
        
        #next step is to verify the that the fields we are replacing with matches 
        #the field we are replacing in terms of indices/ representations
        
        rep2=[]
        for x in field2[1]:
            if x!=1:
                rep2_1=[y.index_type for y in  x.field.get_indices[0]]
                if len(rep2_1) != len(rep1):
                    raise Exception(" replacing is not possible because fields \\                                    are not compatiable")
                if rep2_1 !=rep1:
                    raise Exception("fields are not compatiable as reps \
                                     are different")
                
            
            else:
                if False in singlet1:
                    raise Exception(" no singlet field cannot be replaced with \\                                      a number")
                pass
            
            

        #preparing list of fields to compare
        
        list1=[field1.field.name, field1.conj.name]
        #print('r3', list1)
        field_conj_list=[x.conj if x!=1 else None for x in field2[1]]
        #print('r4', field_conj_list)
        field_list=[x.field if x!=1 else None for x in field2[1]]
        #print('r5', field_list)
        
        #creating symbols for the constants according to its type
        symbol_list=[]
        symbol_conjugate_list=[]
        for x in field2[0]:
            if isinstance(x, str):
                c=Symbol(x, Real=False)
                c_=c.conjugate()
            else:
                c=x
                c_=c.conjugate()
            symbol_list.append(c)
            symbol_conjugate_list.append(c_)
        #print(symbol_list, symbol_conjugate_list)
        for term in self.L_exp:
            temp_term= copy.deepcopy(term)
            contraction_pattern=self.contraction_pattern(term)
            for field in term['fields']:
                multi=0
                if isinstance(field, FieldMul):
                    f=field.fields[1]
                    multi=1
                else:
                    f=field
                if f.name in list1:
                    which_field=list1.index(f.name)
                    flag=1
                    pos=term['fields'].index(field)
                    #if which_field==0:
                    field_list_updated=[]
                    symbol_list_updated=[]
                    for x in range(len(field_list)):
                        if isinstance(field_list[x], Field):
                            if which_field==0:
                                new_f=deepcopy(field_list[x])
                                symbol_list_updated.append(symbol_list[x])
                            else:
                                new_f=deepcopy(field_conj_list[x])
                                symbol_list_updated.append(symbol_conjugate_list[x])
                            if multi==1:
                                new_field=FieldMul(field.fields[0], new_f)
                            else:
                                new_field=new_f
                            field_list_updated.append(new_field)
                            
                        else:
                            continue
                    temp_term['fields'][pos]=[symbol_list_updated, field_list_updated]
            new_term=self.field_associative2(term, temp_term, contraction_pattern)
            new_lag.extend(new_term)
            
        return(new_lag)
#--------------------------------------------------------------------------------
    def field_associative2(self, ori_term, term, contraction):
        pos_list=[]
        
        for x in range(len(term['fields'])):
            if isinstance(term['fields'][x], list):
                pos_list.append(x)
            else:
                continue
        temp_term=deepcopy(ori_term)
        lag=[ori_term]
        for m in pos_list:
            replace=term['fields'][m]
            final_list=[]
            for t in lag:

               for y in range(len(replace[0])):
               
                  temp=deepcopy(t)
                  temp['fields'][m]=replace[1][y]
                  temp['coupling']*=replace[0][y]
                  final_list.append(temp)
            lag=final_list
        new_lag=[]
        print(len(lag))
        print(contraction)
        
        for x in lag:
            new_term=self.AddTerm(str(x['coupling']), *x['fields'],out=True, contraction_pattern=contraction)
            new_lag.append(new_term)
        
        return(new_lag)
        
            
  #--------------------------------------------------------------------------------
    def split_tuple(self, field_tuple, index):
        print("inside split_tuple")
        dummy_field1= field_tuple[0]
        #new_index= Index(dummy_field1.get_indices()[1][1].intex_type
        new_field_constant=Field(dummy_field1.name)
        return(0)        
#------------------------------------------------------------------------------
    def contraction_pattern(self, term):
        explicit_indices={}
        print(term)
        length=len(term['fields'])
        pattern_list=[]
        track_dict={}
        count=1
        for f in term['fields']:
            print(f)
            pattern=[]
            indices=f.get_indices
            #cases=f.get_indices[1]
            #print(indices, cases)
            for idx in range(len(indices[0])):
                #if (Rep.rep_dict[idx.index_type].group.name in group_names)\
                #and (Rep.rep_dict[idx.index_type].singlet==False):
                if Rep.rep_dict[indices[0][idx].index_type].singlet==True:
                    pattern.append(0)
                elif indices[0][idx].index in track_dict.keys():
                    pattern.append(track_dict[indices[0][idx].index]*indices[1][idx])
                else:
                    track_dict.update({indices[0][idx].index:count})  
                    pattern.append(count*indices[1][idx])
                    count+=1
                

            pattern_list.append(pattern)
        print(pattern_list)
      
        #print(pattern_list)
        return(pattern_list)
                          
                        

        

    

       
   
        
             
                            
                    
         
       
             

                                        
