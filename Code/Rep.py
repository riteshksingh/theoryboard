'''
Module : irreps for SLA
Language : python3
#==================================================
'''

from function import *
import copy
from sympy import eye, together, SparseMatrix,\
Rational, sqrt, Matrix, zeros, Symbol, solve, Eq
from LieAlgebra import LieAlgebra

class Representation(object):
	# HOW? Highest weight theorem governs this module
	# irrep ---> (simple Lie algebra, highest weight)
	def __init__(self,LA,hweight):

		########## LA check 
		if not isinstance(LA,LieAlgebra):
			raise TypeError('The 1st arg of Representation class must be a \
							LieAlgebra instance')
		self.LA=LA

		# WHAT? A positive roots method is run here  (wlog pos2_roots()) 
		# WHY? need info such as klist, positiveroots as raw material  
		#		in class methods like weyldim and multiplicity
		kl={}
		for key,tup in LA.pos2_roots()[0].items():
			makedict(tup[2]+1,key,kl)
		self.klist=kl
		self.positiveroots=LA.pos2_roots()[0]

		########## hweight check 
		if not isinstance(hweight, list):
			raise TypeError("highest weight should be input in list format.")
		if len(hweight) != self.LA.rank:
			raise ValueError('highest weight inconsistent with SLA rank.')
		if not all(isinstance(i, int) for i in hweight):
			raise TypeError("all entries in the highest weight must be integers")
		if not all(i>=0 for i in hweight):
			raise ValueError('all entries must be non negative integer')

		self.hweight=hweight

		########## IDENTIFYING SPINOR REPRESENTATIONS
		# Context : Standard Orthogonal Algebra
		if self.LA.fam=='b':
			if self.hweight==ntuple(self.LA.rank,self.LA.rank-1):
				print('The given representation is a spinor representation')
		elif self.LA.fam=='d':
			if self.hweight==ntuple(self.LA.rank,self.LA.rank-1):
				print('The given representation is a spinor representation')
			elif self.hweight==ntuple(self.LA.rank,self.LA.rank-2):
				print('The given representation is a spinor representation')
		else:
			pass

		########## PRIVATE/INTERNAL ATTRIBUTES (not for user)
		# WHY? Placeholders for internal attr in case explicit printing is sought
		self.NO_path=set()	# WHAT? to store non existent paths of descent
		self.scalar_pdt_list={} # WHAT? to store scalar products b/w states of irrep
		self.sp_list={}
	
		#-------------------methods run Automatically in init----------
		self.multiplicity()

		#--------more------
		#self.showRep_to_User()

	########################## 
	# OPTIONAL METHOD 
	def showRep_to_User(self):
		print('')
		print("Dimension of the representation is {}".format(self.weyldim))
		print('')
		print('Weights of the representation in Dynkin Labels')
		print("{} \t {} \t {}" .format('level','weights', 'connections'))
		print('--------------------------------------------------------')
		for level, weights in self.weights.items():
			for keys, values in weights.items():
				print("{} \t {} \t {}" .format(level, keys , values))
		print('--------------------------------------------------------')
		print('state multiplicity of the weights')
		print("{} \t {}" .format('weight', 'multiplicity'))
		print('--------------------------------------------------------')		
		for keys in self.mlist:
			print("{} \t {}".format(keys, self.mlist[keys]))

	##########################
	# WEIGHT-WEIGHT DOT
	# WHAT? dot product calculator b/w any two weights of irrep.
	# HOW? self.ww_dot(t1,t2) where t1,t2 are weight tuples of irrep in Dynkin format
	def ww_dot(self,t1,t2):
		la=self.LA
		WW=la.WW_Matrix
		k=la.rank
		dot=0
		for i in range(k):
			for j in range(k):
				dot+=Rational(WW[i,j]*t1[i]*t2[j])
		return dot	

	##########################
	# CASIMIR ELEMENT
	# WHY interesting? ---> Casimir element x IdentityMatrix is CASIMIR OPERATOR
	# HOW ? ---> Formula given in Cahn using highest weight and twodelta
	@property
	def Casimir(self):
		wt1=t_add(self.hweight,self.LA.twodelta)
		casimir = self.ww_dot(wt1,self.hweight)
		self.casimir=casimir
		return 'The CASIMIR Element',casimir

	##########################
	# MULTIPLICITY OF WEIGHTS (i.e. dimension of each weight space).
	# WHY?  shorter method w/o explicit states identification.
	# HOW?  by Freudenthal's Recursion. Tools : weights of irrep & SLA positive roots.
	def multiplicity(self): 
		Lam=self.hweight
		twdl=self.LA.twodelta
		#-----------------------
		Wlist=[] # Wlist is the list of all weights of irrep
		for levels in self.weights:
			for key,value in self.weights[levels].items():
				Wlist.append(key)
		self.Wlist=Wlist
		#-----------------------
		self.mlist={} # to store weight vs multiplicity.
		ML=self.mlist
		makedict(tuple(self.hweight),1,ML)		
		dim=1
		mlistprime={}

		# we go level by level in the weight diagram.
		for levl in self.weights:
			k1=levl
			if k1!=0:
				# for loop on the weights in level = k1 (a weights LEVEL)
				for M in self.weights[k1].keys():
					denom=self.ww_dot(t_add(Lam,t_add(twdl,M)),t_subs(Lam,M))
					nume=0	
					# recall: self.klist=> level vs positive roots 
					# We deal with each level one by one. 				
					for k2 in self.klist.keys(): # k2 is a positive roots LEVEL
						if k1>=k2:	
							for v1 in self.klist[k2]:
								for k in range(1,int(Rational(k1,k2))+1):
									keff=k1-k*k2
									if keff>=0:
										M1=t_add(M,t_mult(k,v1))
										if M1 in self.weights[keff].keys() :
											#Recursion Formula
											nume += ML[M1][0]*2*self.ww_dot(v1,M1)											
					# The Multiplicity Formula (using Freudenthal Recursion)
					mlt = Rational(nume,denom)
					makedict(M,mlt,ML)
					# dim by adding successively each multiplicity
					dim = dim+mlt

		# a check for the dimension of irrep 
		self.dim_from_multiplicity=dim

	##########################
	# WEYL's DIMENSION FORMULA (only 2delta & positive roots are used in the formula) 
	# to calculate the dim of irrep w/o explicit weights/states/multiplicities. 
	# WHY? to quickly get an idea of how large an irrep is.
	@property
	def weyldim(self):
		delta=t_mult(Rational(1/2),self.LA.twodelta)
		wt1=t_add(self.hweight,delta)									 
		# involves products of terms. so we initialise as:
		weylD=1

		for e in self.positiveroots.keys():
			weylD=Rational(self.ww_dot(wt1,e),self.ww_dot(delta,e))*weylD
		self.weylD=weylD
		return weylD

	##########################
	# WHAT? path of descent from highest weight state,|mu>. Checks existence (T/F)
	# 		Canoically, paths are lowering chains. For eg:- [3,4,0,7,1]|mu>
	# HOW? self.NO_path book-keeps forbidden paths. If in self.NO_path, then trivially
	#		return FALSE. Else do EXISTENCY CHECK.
	# NOTE: all allowed paths need not be Linearly Independent.	
	def path_exists(self,lowering_chain,weights):		
		if tuple(lowering_chain) in self.NO_path:
			return False
		current_wt=tuple(self.hweight)
		current_level=0
		# EXISTENCY CHECK
		# idea is to trace the whole path respecting connections.
		for SR in reversed(lowering_chain):
			# reversed b/c (eg: In [3,4,0,7,1]|mu>) , rightmost acts first
			try:
				# Is current_wt 'connected' by SR to another weight in the wt diag?
				current_wt=weights[current_level][current_wt]['connections'][SR]
				current_level+=1
			except KeyError:
				# store it in the set NO_path
				self.NO_path.add(tuple(lowering_chain))
				return False
		return True

	##########################
	# WHAT? ---> SCALAR PRODUCT between STATES/PATHS
	# WHY? ---> primarily used to RESOLVE DEGENERACY of STATES in irrep
	def scalar_pdt(self, s1, s2, weights):
		la=self.LA
		SRlsq=la.SR_lensq_list
		CM=la.Cmatrix

		# PATHS s1 and s2 are lowering chains
		if len(s1)==0 and len(s2)==0:							
			return 1  # <mu|mu> = 1, Assumption: hweight state is normalised.
		if self.path_exists(s1,weights)==False:
			return 0
		if self.path_exists(s2,weights)==False:
			return 0

		# self.scalar_pdt_list stores (s1,s2) vs <s1|s2> 
		# Check if given (s1,s2) already exists in self.scalar_pdt_list
		try:
			return self.scalar_pdt_list[(tuple(s1),tuple(s2))]
		except KeyError:
			pass
		try:
			return self.scalar_pdt_list[(tuple(s2),tuple(s1))]
		except KeyError:
			pass		

		# Else, compute <s1|s2> from scratch.Store in self.scalar_pdt_list,thereafter. 
		result = 0
		moving_operator =s1[0]
		for i in range(len(s2)):
			if moving_operator==s2[i]:
				SR_len=SRlsq[moving_operator]
				product = self.hweight[moving_operator]
				for j in range(i+1,len(s2)):
					product -= CM[s2[j],moving_operator]
				product*=Rational(1,2)*SR_len
				new_s1=s1[1:]
				new_s2=s2[:i]+s2[i+1:]
				result+=product*self.scalar_pdt(new_s1,new_s2,weights)
		result=together(result)
		self.scalar_pdt_list[(tuple(s1),tuple(s2))]=result
		return result

	##########################
	# WHAT? calculate WEIGHTS of irrep as a dictionary
	# 		level ~> weight ~> nearest (1 step below) connections
	# HOW? starting from hweight and going down in steps using Master Formula
	@property
	def weights(self):
		la=self.LA
		k=la.rank

		# p(raising) and q(lowering) lists
		ps = [0 for i in range(k)]
		qs = None
		fully_lowered =[False for i in range(k)]
		# level=0 is assigned to highest weight. level vs weight is thus, descending. 
		level=0
		weights = {0 :{tuple(self.hweight):{"ps" : ps, "qs" : qs,\
					"fully_lowered" : fully_lowered,"connections" : {}}} }
		# WL_conn to store p,q info againt each weight
		WL_conn={}	

		while True:
			try:
				wgt_this_lvl=weights[level]
			except KeyError:
				break

			for wgt in wgt_this_lvl:
				ps = wgt_this_lvl[wgt]["ps"]
				fully_lowered = wgt_this_lvl[wgt]["fully_lowered"]
				qs = [w_i+p_i for w_i,p_i in zip(wgt, ps)]
				WL_conn[wgt]=[ps,qs]
				for i in range(k):
					if qs[i] > 0 and not fully_lowered[i]:
						q_i = qs[i]
						current_level=level
						current_weight=wgt
						p_i = ps[i]
						for lower_num in range(q_i):
							new_level= current_level + 1
							new_p_i  = p_i + 1
							this_simple_root = la.simpleroots[i]
							new_weight=t_subs(current_weight,this_simple_root)

							# Baggage for each weight
							first_time={"ps": k*[0], "qs" : None,\
									"fully_lowered" : k*[0],"connections":{}}
							try:
								lvl_check = weights[new_level]
								try:
									lvl_w_check = weights[new_level][new_weight]
								# store the new weight in the new level
								except KeyError:
									weights[new_level][new_weight]=first_time
							# create the new level. Then store the new weight there
							except KeyError:
									weights[new_level]={new_weight :first_time}

							weights[new_level][new_weight]["ps"][i] = new_p_i
							# The True-False key help to avoid redundant calculations
							# since, new_weight lies in the middle of the i chain
							weights[new_level][new_weight]["fully_lowered"][i]=True
							weights[current_level][current_weight]\
													["connections"][i] = new_weight

							current_level = new_level
							current_weight = new_weight
							p_i = new_p_i
						wgt_this_lvl[wgt]["fully_lowered"][i]=True
			level+=1
		for level in weights:
			for wgt in weights[level]:
				weights[level][wgt]={"connections":weights[level][wgt]["connections"]}
		self.WL_conn=WL_conn
		return weights

	##########################
	# EXHAUSTIVE CONNECTIONS
	def Wconn_exhaustive(self):
		self.weights
		M={}
		for wgt in self.WL_conn:
			M[wgt]={}
			for i in range(self.LA.rank):
				if self.WL_conn[wgt][0][i]!=0:
					M[wgt]['+'+str(i)]=t_add(wgt,self.LA.simpleroots[i])
				if self.WL_conn[wgt][1][i]!=0:
					M[wgt]['-'+str(i)]=t_subs(wgt,self.LA.simpleroots[i])
		return(M)
		

	#===========================================================================
	# WHAT? Reality Property of Representations (complex vs real) 
	def reality(self):
		wl=self.Wlist
		nwl=[t_mult(-1,wgt) for wgt in wl]
		if set(wl)==set(nwl):
			print('The given Representation is a real representation')
		else:
			print('The given Representation is a complex representation')

	#===========================================================================
	# WHAT? ---> a Private Function to create 'states' dictionary
	# HOW does the scheme work? ---> 
	# 	- each level has corresponding weights (level denoted by lowering steps);
	#	- each weight has an associated weight space. 
	#	  "states" list records the LI paths for each weight space. In a "states"
	#	  list, each state is denoted by : a lowering chain, a normalizer,
	#	  and a "1-step down" info recorder : "matrix_elt_info"
	#	- in presence of multiplicity, "rotn_to_ob" matrix carries out 
	#					the orthogonalisation of degenerate states.
	def __states(self):

		la=self.LA
		weights= self.weights

		self.scalar_pdt_list={}
		self.NO_path=set()

		# level=0 => hweight
		states= { 0 : {tuple(self.hweight) : { "states" : \
				[ { "lowering_chain" : [], "normalizer" : 1,\
					"matrix_elt_info" : [] } ], "rotn_to_ob" : Matrix([[1]]) }}}

		max_step = max(weights)	# maximum number of lowering steps

		for Level in range(max_step):
			for weight in weights[Level]:
				connections = weights[Level][weight]["connections"]
				for sroot_num in connections:
					new_weight = connections[sroot_num]
					for state_num in range(len(states[Level][weight]["states"])):
						# a state (note: it is a dict) for the current weight
						state = states[Level][weight]["states"][state_num]
						normalizer = state["normalizer"]
						lowering_chain = state["lowering_chain"]

						# new lowering chain by moving in the direction of sroot_num
						new_lowering_chain = [sroot_num] + lowering_chain 
						#------------------
						# SANITY CHECK
						if self.path_exists(new_lowering_chain,weights)==False:
							exit()	
						#------------------
						# other identifiers for a state:
						new_scalar_product = sqrt(self.scalar_pdt(new_lowering_chain,\
												new_lowering_chain, weights))
						if new_scalar_product>0:
							new_normalizer = 1/new_scalar_product
							new_level = Level + 1
							new_state = {"lowering_chain" : new_lowering_chain,\
									"normalizer" : new_normalizer, "matrix_elt_info" :\
								[{"level" : Level, "weight":weight,"state_num" :\
								state_num, "direction" : sroot_num,"matrix_elt": \
								Rational(1,1)*new_normalizer/normalizer}]}
						
							try:
								level_check = states[new_level]

								try:
									weight_check = states[new_level][new_weight]
									states[new_level]\
											[new_weight]["states"].append(new_state)
								except KeyError:
									states[new_level][new_weight] = {"states" :\
											[new_state], "rotn_to_ob" : Matrix([[1]]),}

							except KeyError:
								states[new_level]={new_weight:{"states":\
										[new_state], "rotn_to_ob" : Matrix([[1]]),}}


		# LClist :      a dictionary with allowed paths as key and expression as
		#				linear combination of LI paths as value.						
		LClist={}
		# Keeping only LINEARLY INDEPENDENT states
		LIS=copy.deepcopy(states)

		for new_level in states:	
			# Resolving DEGENERACY of States a Weight at level=Level+1
			for weight in weights[new_level]:
				states_for_this_weight = states[new_level][weight]["states"]
				degeneracy = len(states_for_this_weight)
				if degeneracy == 1:
					# there is nothing to check, hence continue
					LClist[tuple(states_for_this_weight[0]["lowering_chain"])]\
								=[(1,states_for_this_weight[0]["lowering_chain"])]
					continue

				scalar_product_matrix = Matrix(degeneracy, degeneracy,\
										lambda i,j : states_for_this_weight[i]\
										["normalizer"]*states_for_this_weight[j]\
										["normalizer"]*self.scalar_pdt\
										(states_for_this_weight[i]["lowering_chain"],\
										states_for_this_weight[j]["lowering_chain"]\
										,weights) if i > j else 0)
				# scalar product is symmetric. <s1|s2> = <s2|s1>
				scalar_product_matrix += scalar_product_matrix.T
				# The <si|si> part
				scalar_product_matrix += eye(degeneracy)
				#=================Scalar Product Matrix===========

				# converts MATRIX to reduced row-echelon form : rref()
				rref = scalar_product_matrix.rref()
				#=================================================
				# collects all those which are dependent
				dependents = [i for i in range(degeneracy) if i not in rref[1]]
				for i in rref[1]:
					X=states_for_this_weight[i]["lowering_chain"]
					LClist[tuple(X)]=[(1,X)]
				for i in dependents:
					X=states_for_this_weight[i]["lowering_chain"]
					LClist[tuple(X)]=[(rref[0][j,i]*states_for_this_weight[j]\
					["normalizer"]/states_for_this_weight[i]["normalizer"],\
					states_for_this_weight[j]["lowering_chain"]) for j in rref[1]]

				# ADDITIONAL MATRIX elts (if any)

				# The linearly independent column numbers of rref of scalar pdt matrix
				# correspond to the LI states.
				for independent in rref[1]:	
					for state_num in range(degeneracy):
						if state_num != independent:
							# If independent, then no projection
							if state_num in rref[1]:
								continue
							# If dependent, to find proper projection

							# About the state in concern, i.e. 
							# states_for_this_weight[state_num]
							state_num_LC = states_for_this_weight[state_num]\
																	["lowering_chain"]
							state_num_norm = states_for_this_weight[state_num]\
																	["normalizer"]
							state_num_MeI =states_for_this_weight[state_num]\
																["matrix_elt_info"]

							# Parent history
							parent_level = state_num_MeI[0]["level"]
							parent_weight = state_num_MeI[0]["weight"]
							parent_state_num = state_num_MeI[0]["state_num"]
							direction = state_num_MeI[0]["direction"]
							parent_mel= state_num_MeI[0]["matrix_elt"]

							parent_LC = states[parent_level][parent_weight]\
									["states"][parent_state_num]["lowering_chain"]
							parent_norm = states[parent_level][parent_weight]\
										["states"][parent_state_num]["normalizer"]
							
							
							for tup in LClist[tuple(parent_LC)]:
								a1=tup[0]
								LCa1=tup[1]
								nza1=1/sqrt(self.scalar_pdt(LCa1,LCa1,self.weights))
								for snm in range(len(states[parent_level]\
													[parent_weight]["states"])):
									if states[parent_level][parent_weight]\
									["states"][snm]["lowering_chain"]==LCa1:
										st1=snm 
								mel=Rational(1,1)*parent_mel*\
													state_num_norm*a1/parent_norm
								
								if mel == 0:
									continue

								states[new_level][weight]["states"][independent]\
											["matrix_elt_info"].append({"level" : \
											parent_level, "weight":parent_weight,\
											"state_num":st1,"direction":direction,\
												"matrix_elt" : mel })

				# WHAT? ---> "rot_to_ob"
				# WHY? ---> matrix to ORTHOGONALISE LI degenerate states
				if degeneracy == 1:
					# nothing to rotate!
					states[new_level][weight]["rotn_to_ob"] = Matrix([[1]])
				else:
					norm_matrix = scalar_product_matrix.extract(rref[1], rref[1])
					rotation_to_ob = gram_schmidt_rotation(norm_matrix)
					states[new_level][weight]["rotn_to_ob"] = rotation_to_ob
					LIS[new_level][weight]["rotn_to_ob"]=rotation_to_ob
				# Keeping only LINEARLY INDEPENDENT states
				LIS[new_level][weight]["states"] = [states_for_this_weight[i] \
																for i in rref[1]]

		self.LClist=LClist
		# INDEXING each STATE uniquely
		state_index = 0
		states=LIS
		for level in states:
			for weight in states[level]:
				states[level][weight]["start_index"] = state_index
				for state_num in range(len(states[level][weight]["states"])):
					states[level][weight]["states"][state_num]["index"]=(state_index)
					state_index += 1
				states[level][weight]["end_index"] = state_index - 1
		print("total number of states is {}".format(state_index))
		# after del, self will have no attributes called NO_path or scalar_pdt_list
		self.scalar_pdt_list #del self.scalar_pdt_list
		self.NO_path #del self.NO_path


		return states


	#----------------------------------------------------------------------------------
	# Calling STATES [External Function]
	def states(self):
		if not hasattr(self, "stored_states"):
			
			self.stored_states = self.__states()

		# Representation dimension calculation by brute force : i.e. counting states. 
		repdim=0
		for lvl in self.stored_states:
			for wgt in self.stored_states[lvl]:
				degeneracy = len(self.stored_states[lvl][wgt]["states"])
				repdim += degeneracy
		self.repdim=repdim

		return self.stored_states
		

	#================================================================================
	# AUXILLARY FEATURES from STATES calculation

	# Neat Printing of STATES (if the print commands below are uncommented)
	# levind; wind : level vs indices ; weight vs indices, respectively.
	def print_states(self):
		S=self.states()
		levind={}
		wind={}
		
		#print('')
		#print('level','weight','state','index')
		#print('')
		for l in S:
			for w in S[l]:
				for i in range(len(S[l][w]['states'])):
					Z=S[l][w]['states'][i]
					#print(l,w,Z['lowering_chain'],Z['index'])
					makedict(l,Z['index'],levind)
					makedict(w,Z['index'],wind)
		self.levind=levind
		self.wind=wind

	# state index vs lowering chain (LC)
	def indexedLC(self):
		IC={}
		IW={}
		S=self.states()
		for l in S:
			for w in S[l]:
				for i in range(len(S[l][w]['states'])):
					Z=S[l][w]['states'][i]
					IC[Z['index']]=Z['lowering_chain']
					IW[Z['index']]=w
		# IW stores index vs weight
		self.IW=IW
		return IC

	# opposite of indexedLC, i.e., LC vs index
	def opp(self):
		X={}
		Y=self.indexedLC()
		for keys, values in Y.items():
			X[tuple(values)]=keys
		return X
			


	# short cut for states info.
	def fleeting(self):
		IC={}
		S=self.states()
		for l in S:
			for w in S[l]:
				for i in range(len(S[l][w]['states'])):
					Z=S[l][w]['states'][i]
					IC[Z['index']]=\
						(Z['lowering_chain'],Z['normalizer'],Z['matrix_elt_info'])
		# IW stores index vs weight
		return IC 


	#================================================================================
	# RAISING OPERATOR

	def Rai(self,y,Ch):
		# Ch above is a dict {lowering chain: factor, ...}
		# such a chain represents a vector in the irrep.
		l={}
		RL=self.LClist
		print(RL)
		for St in Ch:
			St=list(St)
			for i in range(len(St)):
				if y==St[i]:
					new1=St[0:i]+St[i+1:]
					print('new1',new1)
					if tuple(new1) in RL:
						SR_len=self.LA.SR_lensq_list[y]
						#print(SR_len,'s')
						product = self.hweight[y]
						#print(product)
						print(St[i+1:])
						for j in St[i+1:]:
							product -= self.LA.Cmatrix[St[j],y]
						
						#print('j',product)
						product*=Rational(1,2)*SR_len
						for entry in RL[tuple(new1)]:
							try:
								l[tuple(entry[1])]+=entry[0]*product*Ch[tuple(St)]
							except:
								l[tuple(entry[1])]=entry[0]*product*Ch[tuple(St)]
		return l
	#==============================================================================
	# MATRICES of the Representation of an SLA-------------in an orthonormal basis.
	#------------in dict format. KEY : generator; VALUE : Matrix
	# integer : Cartan Generator, tuple: Ladder Generator

	# NOTE : SparseMatrix creates the final matrix once all non-zero matrix
	#		 elements are set (used to create matrices with a lot of zeroes)

	def matrices(self):
		la=self.LA
		k=la.rank
		CM=la.Cmatrix
		states = self.states()
		Dim=self.weyldim
		
		RL=self.LClist

		representation_matrices = {}

		# -----------------------------------------------------------------------
		# REP MATRICES for Cartan Generators (Basis independent)
		# 					H_{\alpha_i}=representation_matrices[i]
		cartan_matrices_nonzero = [{} for i in range(k)]
		for level in states:
			for weight in states[level]:
				states_for_this_weight = states[level][weight]["states"]
				for state in states_for_this_weight:
					state_index = state["index"]
					for i in range(k):
						cartan_matrices_nonzero[i][(state_index,state_index)]=\
											 self.ww_dot(la.simpleroots[i],weight)

		for i in range(k):
			representation_matrices[i]=\
								SparseMatrix(Dim,Dim,cartan_matrices_nonzero[i])

			# Normalising Factor
			representation_matrices[i]*=Rational(1, la.SR_lensq_list[i])

		# Cartan Gens all have root 0. Hence, degeneracy. 
		# Need to orthogonalise.
		HL=[]
		hl=[{0:1}]
		for n in range(1,k):
			var=[]
			for i in range(0,n+1):
				var.append('a'+str(i))
			var=[Symbol(i) for i in var]
			Var={}
			# Assigning variable co-effs to each SR
			for i in range(0,n+1):
				Var[i]=var[i]

			EQ=[]

			for live_L in hl:
				eq=0
				for key1 in Var:
					for key2 in live_L:
						eq+=Var[key1]*live_L[key2]*\
							self.ww_dot(la.simpleroots[key2],la.simpleroots[key1])
				EQ.append(eq)
			EQ.append(var[0]-1)
			res=solve([Eq(z) for z in EQ],var)

			new_L={}
			for var in Var:
				new_L[var]=res[Var[var]]
			hl.append(new_L)
				
		self.hl=hl
		for L in hl:
			nor=0
			for p in L:
				for q in L:
					nor+=L[p]*L[q]*self.ww_dot(la.simpleroots[p],la.simpleroots[q])
			H_Mat=zeros(self.weyldim)
			for j in L:
				H_Mat+=int(L[j])*representation_matrices[j]
			HL.append(1/sqrt(nor)*H_Mat)

		self.HL=HL			

		# -----------------------------------------------------------------------
		# REP MATRICES for SIMPLE ROOTS GENERATORS				
		for i in range(k):
			simple_root_i_nonzero = {}
			for level in states:
				for weight in states[level]:
					states_for_this_weight = states[level][weight]["states"]
					for state in states_for_this_weight:
						state_nz=state["normalizer"]
						state_LC=state["lowering_chain"]
						state_index1 = state["index"]
						
						new_st=[i]+state_LC
						if tuple(new_st) in RL:
							expansion= RL[tuple(new_st)]
							for entry in expansion:
								state_index2=self.opp()[tuple(entry[1])]
								c1=entry[0]
								n1=self.fleeting()[state_index2][1]

								mel=Rational(1,1)*state_nz*c1/n1
								simple_root_i_nonzero[(state_index2, state_index1)]\
														= Rational(1,1)*mel

			simple_root_i_matrix = SparseMatrix(Dim,Dim,simple_root_i_nonzero)
			
			# negative Simple Roots Generators : E_{-\alpha_i}=simple_root_i_matrix
			representation_matrices[t_mult(-1,tuple(CM.row(i)))] = simple_root_i_matrix

			#representation_matrices[t_mult(-1,tuple(CM.row(i))]*=sqrt(Rational(1, la.SR_lensq_list[i]))

			# -----  ----- ----- ----- ----- ----- ----- ----- ----- ----- 
			# REP MATRICES for Positive Simple Root Generators
			representation_matrices[tuple(CM.row(i))] = (simple_root_i_matrix.T)


		
		# ORTHONORMALISATION 

		rotation_to_ob_nonzero = {}
		for level in states:
			for weight in states[level]:
				rotation_matrix = states[level][weight]["rotn_to_ob"]
				#print(rotation_matrix)
				start_index = states[level][weight]["start_index"]
				end_index = states[level][weight]["end_index"]
				rotation_matrix_nonzero = dict([((i,j),rotation_matrix[i-start_index,\
				j-start_index]) for i in range(start_index, end_index+1) for j in \
					range(start_index, end_index+1)])
				rotation_to_ob_nonzero.update(rotation_matrix_nonzero)

		rotation_to_ob = SparseMatrix(Dim, Dim,rotation_to_ob_nonzero)

		for key in representation_matrices:
			# change of basis via similarity transform
			rotated_matrix = rotation_to_ob.multiply(representation_matrices[key].\
														multiply(rotation_to_ob.T))
			rotated_matrix.simplify()
			representation_matrices[key] = rotated_matrix
		
		# -----------------------------------------------------------------------
		# REP MATRICES for OTHER Generators 
		positive_roots = la.pos2_roots()[0]
		positive_roots_list = [(key,) +  positive_roots[key] \
		for key in positive_roots if isinstance(positive_roots[key][1],list)]
		positive_roots_list = sorted(positive_roots_list,key = lambda item : item[3])		

		for i in range(len(positive_roots_list)):
			matrix1 = representation_matrices[positive_roots_list[i][2][0]]
			matrix2 = representation_matrices[positive_roots_list[i][2][1]]
			# commutator
			positive_root_matrix = matrix1.multiply(matrix2).\
										add(-matrix2.multiply(matrix1))
			factor = positive_roots_list[i][1]
			root = positive_roots_list[i][0]

			# Generator matrix for POSITIVE ROOT
			representation_matrices[root]=positive_root_matrix.\
											applyfunc(lambda i : factor*i)
			# Generator matrix for NEGATIVE ROOT
			representation_matrices[t_mult(-1,root)] = representation_matrices[root].T


		return representation_matrices
		

 							
