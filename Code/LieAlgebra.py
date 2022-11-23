'''
Module : Simple Lie Algebra
Language : python3

#NOTE:
#property decorators, which are used in some methods, let us access those in the 
#fashion of accessing attributes. No particular benefit, just for its own sake.
#======================================================================================
'''

from sympy import Rational, sqrt, Matrix, zeros
from function import *

class LieAlgebra(object):
	def __init__(self,fam,dim):
	#instance variables--------->alg family & dimension. 
	#family ------->standard name/common name <----- both acceptable.
		self.fam = fam
		self.dim = dim

		# DC positiveroots method calulates the positive roots and
		# also expresses each root as commutator of a simple root with 
		# a positive root that is one level below.

		# AB positiveroots method calculates the positive roots and
		# records the information about how many of each simpleroot
		# are present in a given positive root. Exact SR chain not recorded.
		# Only the combination is. 

		# BINARY SWITCHes pmeth1 & pmeth2 for AB method and DC method resp.
		# In the init method, (pmeth1,pmeth2)=(0,0) implies that both are off.
		# positive roots switches are used to make certain calculations method indep.
		self.pmeth1=0
		self.pmeth2=0
		
		

		#LA CONVERSION to standard form ------ family & rank
		fam = self.fam
		dim = self.dim
		#preliminary check
		if fam not in ['u','a','su','c','sp','b','d','so','e','f','g']:
			print("Lie Algebra not Recognized. Choose from 'u','a','su',\
										'c','sp','b','d','so','e','f','g'")
			exit()
		if not isinstance(dim,int):
			print('dimension should be +ve integer valued')
			exit()
		#Unitary 
		if fam == 'u':
			self.fam = 'u'
			if dim != 1:
				print('for u, dimension is invalid')
				exit()
			else:
				print('u(1) has no simple root')
				exit()
		#Standard Unitary
		elif fam == 'su' or fam == 'a':
			if fam=='a':
				self.rank = dim
			else:
				if dim<2:
					print('invalid dimension')
					exit()
				else:
					self.rank = dim -1
					self.fam = 'a'
					
		#Symplectic
		elif fam == 'sp' or fam == 'c':
			if fam=='c':
				self.rank=dim
			else:
				if dim%2 !=0:
					print('invalid dimension')
					exit()
				else:
					self.rank=int(dim/2)
			self.fam='c'
		#Standard Orthogonal
		elif fam=='so' or fam=='b' or fam=='d':
			if fam=='b':
				if self.dim==2:
					self.fam='c'
					self.rank=2
				else:
					self.fam='b'
					self.rank=dim
			elif fam=='d':
				if dim<3:
					print('invalid dimension')
					exit()
				elif dim==3:
					self.fam='a'
					self.rank=3
				else:
					self.fam='d'
					self.rank=dim
			else:
				if dim%2 == 0 and dim>4:
					if dim==6:
						self.rank=3
						self.fam='a'
					else:
						self.rank=int(dim/2)
						self.fam='d'
				elif dim==5:
					self.fam='c'
					self.rank=2
				
				elif dim>=5:
					self.rank=int((dim-1)/2)
					self.fam='b'
				else:
					print('invalid dimension')
					exit()
		#Exceptional lie algebras
		elif fam=='g':
			if dim==2:
				self.fam='g'
				self.rank=dim
			else:
				print('for family g, only dimension 2 is valid')
				exit()
		elif fam=='f':
			if dim==4:
				self.fam='f'
				self.rank=dim
			else:
				print('for family f, only dimension 4 is valid')
				exit()
		else:
			self.fam='e'
			if dim<6 or dim>8:
				print('for family e, only dimension 6,7,8 are valid')
				exit()
			else:
				self.rank=dim 
	
	# Dynkin Rep heavily employed in subsequent methods

	##########################
	# CARTAN MATRIX of SLA
	@property
	def Cmatrix(self):
		family=self.fam
		k=self.rank
		cmatrix=zeros(k)
		if family=='g':
			cmatrix=Matrix([[2,-1],[-3,2]])
		elif family=='f':
			cmatrix=Matrix([[2,-1,0,0],[-1,2,-1,0],[0,-2,2,-1],[0,0,-1,2]])
		elif family=='e':
			if k==6:
				cmatrix=Matrix([[2,-1,0,0,0,0],[-1,2,-1,0,0,0],[0,-1,2,-1,0,-1],\
								[0,0,-1,2,-1,0],[0,0,0,-1,2,0],[0,0,-1,0,0,2]])
			elif k==7:
				cmatrix=Matrix([[2,-1,0,0,0,0,0],[-1,2,-1,0,0,0,0],\
								[0,-1,2,-1,0,0,0],[0,0,-1,2,-1,0,-1],\
								[0,0,0,-1,2,-1,0],[0,0,0,0,-1,2,0],[0,0,0,-1,0,0,2]])
			
			else:
				cmatrix=Matrix([[2,-1,0,0,0,0,0,0],[-1,2,-1,0,0,0,0,0],\
							[0,-1,2,-1,0,0,0,0],[0,0,-1,2,-1,0,0,0],\
								[0,0,0,-1,2,-1,0,-1],[0,0,0,0,-1,2,-1,0],\
									[0,0,0,0,0,-1,2,0],[0,0,0,0,-1,0,0,2]])
		else:
			cmatrix=zeros(k)
			for i in range(k):
				for j in range(k):
					if i==j:
						cmatrix[i,j]=2
					elif abs(i-j)==1:
						if family=='a':
							cmatrix[i,j]=-1
						elif family=='c':
							if i!=k-1:
								cmatrix[i,j]=-1
							else:
								cmatrix[i,j]=-2
						elif family=='b':
							if j!=k-1:
								cmatrix[i,j]=-1
							else:
								cmatrix[i,j]=-2
						else:
							if i==k-1 or j==k-1:
								cmatrix[i,j]=0
							else:
								cmatrix[i,j]=-1
					elif family=='d' and abs(i-j)==2 and (j==k-1 or i==k-1):
						cmatrix[i,j]=-1
		

		return cmatrix

	
	##########################
	#SIMPLE ROOTS of the SLA------> in Dykin Representation
	@property
	def simpleroots(self):
		simpleroots=[]
		for i in range(self.rank):
			#simply the Cartan Matrix rows
			simpleroots.append(tuple(self.Cmatrix.row(i)))
		return simpleroots

	##########################
	# From the Dynkin Diagram info (recall we began with Cmatrix ),
	# we can easily write the norm squares in the language of ratio-proportions: 

	# NORM square--------The length squares of simple roots in a given SLA
	# Can also be thought of as ratios among SR lengths
	@property
	def SR_lensq_list(self):

		k=self.rank
		A=self.Cmatrix
		lenlst=[]
		if self.fam in ['a','d','e']:
			lenlst=k*[1,]
		elif self.fam=='g':
			lenlst=[1,3]
		else:
			slst=[]
			lenlst=k*[1]
			for i in range(k):
				for j in range(i+1,k):
					if A[i,j]!=0:
						slst.append(Rational(A[i,j],A[j,i])) 
			slst.append(1)
			for i in range(len(slst)):
				if slst[i]>1:
					for j in range(i+1):
						lenlst[j] = slst[i]*lenlst[j]
				else:
					if slst[i]<1:
						for j in range(i+1,k):
							lenlst[j] = Rational(lenlst[j],slst[i])
		return lenlst

	##########################
	# FUNDAMENTAL WEIGHTS : Using dynkin idea
	# SS[i]=C[i,j]*FF[j] --------> system of Linear Equations
	@property
	def FW(self):
		SS=self.simpleroots
		k=self.rank
		C=self.Cmatrix
		# inverse of Cartan Matrix
		B=C.inv()
		FF=[]
		for i in range(k):
			FF.append(k*(0,))
			for j in range(k):
				FF[i]=t_add(FF[i],t_mult(B[i,j],SS[j]))
		return FF 

	# Note that when expanded in SR basis, the jth component (co-eff of alpha_j) of
	# fundamental weight mu_i is Cmatrix^(-1)[i,j]. 
	# Used in the following method (decorated)for calculation of wgt-wgt matrix.

	##########################
	#  WW_Matrix. Stores dot products between fundamental weights (number=rank)
	#  as elements of the weight-weight matrix
	@property
	def WW_Matrix(self):	
		k=self.rank
		Mx=zeros(k)
		C=self.Cmatrix
		lenlst=self.SR_lensq_list
		# Worthwhile to note that C is not necessarily symmetric.
		CI=Matrix(C)**(-1)
		# finally the wgt-wgt matrix elements. Symmetric. Hence, we can reduce the
		# running cycles of the for loop.
		for i in range(k):
			for j in range(i,k):
				Mx[i,j]=Rational(1,2)*CI[i,j]*lenlst[j]
				Mx[j,i]=Mx[i,j]
		return Mx

	#===============================================================================
	#           METHODS valuation of POSITIVE ROOTS of a SIMPLE LIE ALGEBRA
	#-------------------------------------------------------------------------------

	# POSITIVE ROOTS ------ meth1
	# To Call: self.pos1_roots() ------ 2 new attr get assigned
	# positiveroots and klist
	#=============================
	#procedure involves raising each root at each level to the maximum possible extent
	#using each simple root. The path is not recorded. Only the combination of SRs is.

	def pos1_roots(self):
		# the moment pos1_roots() is called, pmeth1 switch is turned on
		self.pmeth1=1
		
		self.klist={}
		self.positiveroots={}

		k=self.rank
		KL=self.klist			# dictionary of levels
		PR=self.positiveroots	# dictionary of positiveroots


		newlist={}				# purpose : start from a root. 
		klistprime={}			# prime denotes primary level

		# Starting with Simple Roots
		for i in range(k):
			a=self.simpleroots[i]
			# SRs stored in level 1 of KL
			makedict(1,a,KL)
			makedict(a,ntuple(k,-1),newlist)
			makedict(a,tuple(ntuple(k,i)),PR)

		count =0
		while(True):
			klistprime.update(KL)
			for keys in klistprime:
				# starting keys is 1
				if keys>count:
					# b is the associated value list for key keys. 
					# NOTE: makedict generates value lists for dictionary keys.
					b=klistprime[keys]
					# len(b) is number of proots in the key=keys level in KL
					for i in range(len(b)):
						# MASTER FORMULA : Algorithm
						for j in range(k):	# j gives direction.
							# b[i][j] is the jth dynkin label of the proot b[i]
							if b[i][j]<0 and newlist[b[i]][0][j]==0:
								# abs(b[i][j]) is the number of times proot b[i] 
								# can be raised by SR_j
								# NOTE: newlist is a switch facility to 
								# 		remember paths already traversed.
								newlist[b[i]][0][j]=1
								pr=b[i]

								for c in range(1,abs(b[i][j])+1):
									# new positive root
									fr=t_add(b[i],t_mult(c,self.simpleroots[j]))

									# PR stores only the kind of SRs present in
									# each positive root. Not the actual path.
									# NOTE: adding the entries in the value list 
									#		for each +ve root gives the level.
									makedict(fr,t_add(ntuple(k,j),PR[pr][0]),PR)
									makedict(keys+c,fr,KL)

									# if fr is a key in newlist, turn its j^th switch
									# on to signify that j^th direction is taken care
									# of ; if not, then create newlist for fr and turn
									# j^th switch on.
									try:
										newlist[fr][0][j]=1
									except:
										makedict(fr,ntuple(k,-1),newlist)
										newlist[fr][0][j]=1									
									# rise up a level to pr (cursor is here now !!)
									pr=fr
			if len(KL.keys())-1<count:
				#dimension of adjoint rep
				self.adim=2*len(PR.keys())+k
				break
			else:
				count+=1

	#-----------------------------------------------------------------------

	# POSITIVE ROOTS ------ meth2
	# To Call: self.pos2_roots() ------ returns a 2-tuple
	# 1st elt of the tuple is the +ve roots dict, 2nd elt is hweight of ad Rep

	# each positive root comes with the following 3 info :
	# multiplication factor, expression in terms of commutator, level
	#=============================
	# commutators are included in the root dictionary itself
	def pos2_roots(self):
		self.pmeth2=1	
		k=self.rank
		roots_dict={}
		roots_laststep={}
		roots_thisstep={}
		# RL_conn to store p,q info against each root
		RL_conn={}
		for i in range(k):
			SR_i = self.simpleroots[i]
			q_values=t_mult(2,ntuple(k,i))
			# from MASTER FORMULA
			p_values= t_subs(q_values,SR_i)
			roots_laststep[SR_i]=[p_values,q_values,1,SR_i,0]
			RL_conn[SR_i]=[p_values,q_values]
		while True:
			for (roots,value) in roots_laststep.items():
				# remember that value is a list object
				roots_dict[roots]=(value[2],value[3],value[4])
			for roots in roots_laststep:
				p_values=roots_laststep[roots][0]
				q_values=roots_laststep[roots][1]
				level=roots_laststep[roots][4]
				for i in range(k):
					if p_values[i]>0:
						SR_i=self.simpleroots[i]
						newroot=t_add(roots,SR_i)
						try:
							roots_thisstep[newroot][1][i]=q_values[i]+1
						except KeyError:
							# proportionality factor in commutation
							factor=1/sqrt(Rational((q_values[i] +1)*p_values[i],2)\
															*self.SR_lensq_list[i])

							# The order of writing the commutator is crucial.
							commutator=[SR_i,roots]
							newlevel=level+1
							new_q_values=ntuple(k,-1)
							new_q_values[i]+=q_values[i]+1
							roots_thisstep[newroot]=[None,new_q_values,factor,\
															commutator,newlevel]
	
			for roots in roots_thisstep:
				roots_thisstep[roots][0]=t_subs(roots_thisstep[roots][1],roots)
			if len(roots_thisstep)==0:
				adjoint_rep=list([key for key in roots_laststep][0])
				# RL_conn set as an attr in this final step
				self.RL_conn=RL_conn
				return (roots_dict,adjoint_rep)
			else:
				# RL_conn stores p,q info against each root
				for nr in roots_thisstep:
					RL_conn[nr]=[roots_thisstep[nr][0],roots_thisstep[nr][1]]
				roots_laststep=roots_thisstep
				roots_thisstep={}

	#########################
	def adj_rep(self):
		HW=[]
		hw=self.pos2_roots()[1]
		for i in hw:
			HW.append(int(i))
		return HW

	##########################
	# EXHAUSTIVE CONNECTIONS
	def Rconn_exhaustive(self):
		self.pos2_roots()
		M={}
		for root in self.RL_conn:	# Note : RL_conn only stores +ve roots info
			M[root]={}	# in M, we'd like to store both -ve & +ve roots info
			M[t_mult(-1,root)]={}	# -ve data is simply Mirror image of +ve data
			for i in range(self.rank):
				SRi=self.simpleroots[i]
				if self.RL_conn[root][0][i]!=0:
					M[root]['+'+str(i)]=t_add(root,SRi)
					M[t_mult(-1,root)]['-'+str(i)]=t_mult(-1,t_add(root,SRi))
				if self.RL_conn[root][1][i]!=0:
					M[root]['-'+str(i)]=t_subs(root,SRi)
					M[t_mult(-1,root)]['+'+str(i)]=t_mult(-1,t_subs(root,SRi))
		zero_vec=tuple(ntuple(self.rank,-1))
		M[zero_vec]={}
		for j in range(len(self.simpleroots)): 
			M[zero_vec]['+'+str(j)]=self.simpleroots[j]
			M[zero_vec]['-'+str(j)]=t_mult(-1,self.simpleroots[j])
		return(M)

	##########################
	# INDEX vs ROOTS
	def labelled_roots(self):
		L={}
		for i in range(self.rank):
			makedict(0,i,L)
		self.pos1_roots()
		S=self.klist
		j=1
		for i in range(len(self.klist)):
			kmax=max(S)
			for num in range(len(S[kmax])):
				L['+'+str(j)]=S[kmax][num]
				L['-'+str(j)]=t_mult(-1,S[kmax][num])
				j=j+1
			del S[kmax]
		return L
				
						
	##########################
	# TWODELTA
	# two-delta [Independent of whichever positive roots method is used]
	@property
	def twodelta(self):
		twdl=self.rank*(0,)
		# ab initio CASE
		if self.pmeth1==0 and self.pmeth2==0:
			print(""" Positive Roots are Prerequisites to run twodelta.""")
		elif self.pmeth1 !=0:
			for keys in self.positiveroots:
				twdl=t_add(twdl,keys)
		elif self.pmeth2 !=0:
			for keys in self.pos2_roots()[0]:
				twdl=t_add(twdl,keys)
		return twdl

	##########################
	# ROOT-ROOT DOT
	def rr_dot(self,p1,p2):
		WW=self.WW_Matrix
		k=self.rank
		dot=0
		for i in range(k):
			for j in range(k):
				dot+=Rational(WW[i,j]*p1[i]*p2[j])
		return dot


	##########################
	# Commutators of Generators of the Algebra
	# Separate method needed because, pos1_roots method doesn't encompass this
	def commutator(self):

		k=self.rank
		# commutatotors of positive roots are stored here
		self.comdict={}
		comdict=self.comdict
		
		# simpleroots intactly inserted in comdict. (not expressed as commutators)
		for si in self.simpleroots:
			makedict(si,(1,si),comdict)

		# We call AB's pos1_roots 
		self.pos1_roots()
		
		# we deal level by level
		for i in range(1,max(self.klist.keys())+1):
			b=self.klist[i]
			for root in range(len(b)):
				for l in range(k):
					# length square of SR[l]
					ll=self.SR_lensq_list[l] 
					if b[root][l] <0:
						# c takes care of all roots in this chain.
						for c in range(abs(b[root][l])): 
							flag=0
							TADD=t_add(b[root],t_mult(c+1,self.simpleroots[l]))
							if TADD in comdict:
									#print('whee',\
									#t_add(b[root],t_mult(c+1,self.simpleroots[l])))
									flag = 1
							# flag=1 implies that the chosen root exists
							# in a previously constructed SR[l]-directed chain

							# We always seek a fresh SR[l]-chain starting at root 
							if flag != 1:	
								A=t_add(b[root],t_mult(c+1,self.simpleroots[l]))
								prevroot=t_add(b[root],t_mult(c,self.simpleroots[l]))	
								# recall that for pi systems (root systems),
								# the connections can only be 0,1,2,3
								
								# abs(b[root][l]) denotes the length of a chain
								# in terms of edges. 
								# factor : commutation factors case by case						
								if abs(b[root][l])==1:
									factor = 1/sqrt(ll*Rational(1,2))
									makedict(A,(factor,\
											[self.simpleroots[l],prevroot]),comdict)
																				
								elif abs(b[root][l])==2:
									factor = 1/sqrt(ll)
									makedict(A,(factor,\
											[self.simpleroots[l],prevroot]),comdict)
																			
								elif abs(b[root][l])==3:		 									 
									if (c==0 or c==2):
										factor = 1/sqrt(ll*Rational(3,2))
										makedict(A,(factor,\
											[self.simpleroots[l],prevroot]),comdict)	
										 
									elif (c==1):
										factor = 1/sqrt(2*ll)
										makedict(A,(factor,\
											[self.simpleroots[l],prevroot]),comdict)
		for key in comdict:
			comdict[key]=comdict[key][0]
			
 
#==================================== x =============================================
		 

