from sympy import sqrt
'''
All the global functions and attributes used in the code
'''
#-------------------------------------------------------------------------------
IndicesCount=0
indices_name_list=['a','b','c','d','e','f','g','h','i','j','k','m','n', 'o',
                    'p','q','r', 's','t','u', 'v', 'w', 'x', 'y','z']
GaugeFieldNames=['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
GaugeFieldCount=0
#-------------------------------------------------------------------------------
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def check_if_value_available(dict1, key, numerical=True):
    #print('start processing check_if_value_available')
    value=None
    length=len(key)
    #print(length)
    #print(dict1, key)
    for x in dict1.keys():
         #   print(x, key)
        count1=0
        count2=0
        for y in range(length):
            if key[y].isalpha():
                pass
            else:
                count1+=1
                if key[y]==x[y]:
                    count2+=1
        if count1==count2:
            if count1==length and numerical==True:
                value=dict1[x]
            result=1
        else:           
            result=0
            #print(count1, count2, result)
        if result==1:
            break

    return(result, value)   
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def dim_to_cartan_weyl(dim, nroots):
    conversion_dict={}
    conversion_dict2={}
    for y in range(dim):

        if (y<nroots):
            n='p'+str(y+1)
            conversion_dict.update({y: n })
            conversion_dict2.update({dim-y-1:n})
        elif dim-nroots<=y<dim:
            n='m'+str(dim-y) 
            conversion_dict.update({y: n })
            conversion_dict2.update({dim-y-1:n})
        else:
            n='c'+ str(y-nroots+1)
            conversion_dict2.update({y:n})
            conversion_dict.update({y: n })
        
        #conversion_dict2.update({dim-y-1:n})
    #print('matrix column',conversion_dict)
    #print('matrix row', conversion_dict2)
    return(conversion_dict, conversion_dict2)

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def nonzero_matrix(mat):
    length=len(mat)
    rowno=sqrt(length)
    non_zero_index=[]
    for i in range(length):
        if mat[i]!=0:
            row=int(i/rowno)
            column=i%rowno
            list1=[row, column, mat[i]]
            non_zero_index.append(list1)
    return(non_zero_index)


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def indices():
    global IndicesCount
    index=indices_name_list[IndicesCount]
    IndicesCount+=1
    return(index)

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def gaugenames():
    global GaugeFieldCount
    name=GaugeFieldNames[GaugeFieldCount]
    GaugeFieldCount+=1
    return(name)


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def group_same_strings( big_list):
        #for making string set index out of string list
        #[aa,bb,aa,cc,bb,dd] gives [[0,2],[1,4],[3],[5]]
    small_set=set(big_list)
    index_set=[]
    for x in small_set:
        idx=[index for index, value in enumerate(big_list) if value == x]
        index_set.append(idx)
    return(index_set)

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def explicitname(rep_list):
    explicit_name='_x'
    for x in rep_list:
        if x.singlet==False:
            if x.name=='Lorentz':
                pass
                #explicit_name+='_l'
            else:
                explicit_name+='_'+x.index_name
        else:
            explicit_name+='_s'
    return(explicit_name)

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def check_scalar_fermion_rep(rep_list):
    G_list = [x.group.name for x in rep_list]
    G_set = set(G_list)
    if len(G_list) != len(G_set):
        return (1)
    else: 
        return(0)



