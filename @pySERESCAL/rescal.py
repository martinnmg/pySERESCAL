#-*- encoding: ISO-8859-1 -*-
# N. MARTIN 2022/04/29
# (c) Laboratoire Léon Brillouin LLB
# From the 'SERESCAL' MATLAB toolbox by K. Habicht (HZB) 
# References:
# [Habicht2003] K. Habicht et al., J. Appl. Cryst. 36, 1307-1318 (2003)
# [Groitl2018] F. Groitl et al., J. Appl. Cryst. 51, 818-830 (2018)	

from sys import version_info
from numpy import arccos,arcsin,argmin,array,cos,cross,deg2rad,dot,exp, \
	linspace,log,minimum,pi,ones,rad2deg,resize,sin,sign,sqrt,tan,transpose,zeros
from numpy.linalg import det,inv,norm
from scipy.constants import hbar,m_n
from scipy.constants import e as eV_to_J
from math import isnan
from time import sleep

class RESCAL:

	def __init__(self,parlist,method):		
		self.updateParams(parlist,method)
		self.computeCNmatrix()
		self.printTASparams()
	
	def printTASparams(self):
# 		Prints TAS & NRSE parameters in the command line window
		
		A1,A2,A3,A4,A5,A6,phi_S,sigma_S=self.calcAngles()
		ki,kf=self.calcKiKf()
		
		if version_info[:3]<(3,0):
			print '\nTAS parameters:\n---------------'
			if self.fix == 1:
				print 'Q = [',self.Qvec[0],self.Qvec[1],self.Qvec[2],'], hw =',self.En/(1e-3*eV_to_J),'meV (fixed ki)'	
			else:
				print 'Q = [',self.Qvec[0],self.Qvec[1],self.Qvec[2],'], hw =',self.En/(1e-3*eV_to_J),'meV (fixed kf)'	
			print 'ki =',round(ki*1e-10,3),'1/AA, kf = ',round(kf*1e-10,3),'1/A\n' 
			print 'a1 =',round(rad2deg(A1),2),'deg, a2 =',round(rad2deg(A2),2),'deg'		
			print 'a3 =',round(rad2deg(A3),2),'deg, a4 =',round(rad2deg(A4),2),'deg'		
			print 'a5 =',round(rad2deg(A5),2),'deg, a6 =',round(rad2deg(A6),2),'deg\n'	
		else:
			print('\nTAS parameters:\n---------------')
			if self.fix == 1:
				print('Q = [',self.Qvec[0],self.Qvec[1],self.Qvec[2],'], hw =',self.En,'meV (fixed ki)')	
			else:
				print('Q = [',self.Qvec[0],self.Qvec[1],self.Qvec[2],'], hw =',self.En,'meV (fixed kf)')			
			print('a1 =',round(rad2deg(A1),2),'deg, a2 =',round(rad2deg(A2),2),'deg')		
			print('a3 =',round(rad2deg(A3),2),'deg, a4 =',round(rad2deg(A4),2),'deg')		
			print('a5 =',round(rad2deg(A5),2),'deg, a6 =',round(rad2deg(A6),2),'deg\n')	
							
	def updateParams(self,parlist,method):
# 		Gets TAS parameters from the main panel,
#		converts them into S.I. units
#		and stores them into instance attributes
	
		self.method=method
		parlist=[float(i) for i in parlist]
		# Monochromator & analyzer lattice constants (input unit == Å), Å -> m
		self.dM=parlist[0]*1e-10
		self.dA=parlist[1]*1e-10
		# Mosaicities (input unit == arc minutes, FWHM), arcmin -> rad
		FWHMminstoSigdegs=1/(60*2*sqrt(2*log(2))) 
		if parlist[2]<1e-6: # Set mochromator mosaic very small if zero
			self.etaM=1e-6
		else:		
			self.etaM=deg2rad(parlist[2]*FWHMminstoSigdegs)		
		if parlist[3]<1e-6: # Set analyzer mosaic very small if zero
			self.etaA=1e-6
		else:		
			self.etaA=deg2rad(parlist[3]*FWHMminstoSigdegs)			
		if parlist[4]<1e-6: # Set sample mosaic very small if zero
			self.etaS=1e-6
		else:		
			self.etaS=deg2rad(parlist[4]*FWHMminstoSigdegs)
		# Scattering senses
		self.SM=sign(parlist[6])
		self.SS=sign(parlist[7])
		self.SA=sign(parlist[8])
		# Neutron initial/final wavevector (input unit == Å-1), Å-1 -> m-1
		self.kfix=parlist[9]*1e10
		self.fix=parlist[10]
		# Horizontal divergences (input unit == arc minutes, FWHM), arcmin -> rad
		if parlist[11]>1e-6:
			self.alpha0=deg2rad(parlist[11]*FWHMminstoSigdegs)
		else:
			self.alpha0=1e-6
		if parlist[12]>1e-6:
			self.alpha1=deg2rad(parlist[12]*FWHMminstoSigdegs)
		else:
			self.alpha1=1e-6
		if parlist[13]>1e-6:
			self.alpha2=deg2rad(parlist[13]*FWHMminstoSigdegs)
		else:
			self.alpha2=1e-6
		if parlist[14]>1e-6:
			self.alpha3=deg2rad(parlist[14]*FWHMminstoSigdegs)
		else:
			self.alpha3=1e-6			
		# Vertical divergences (input unit == arc minutes, FWHM), arcmin -> rad	
		if parlist[15]>1e-6:
			self.beta0=deg2rad(parlist[15]*FWHMminstoSigdegs)
		else:
			self.beta0=1e-6
		if parlist[16]>1e-6:
			self.beta1=deg2rad(parlist[16]*FWHMminstoSigdegs)
		else:
			self.beta1=1e-6
		if parlist[17]>1e-6:
			self.beta2=deg2rad(parlist[17]*FWHMminstoSigdegs)
		else:
			self.beta2=1e-6
		if parlist[18]>1e-6:
			self.beta3=deg2rad(parlist[18]*FWHMminstoSigdegs)
		else:
			self.beta3=1e-6
		# Sample lattice constants (input unit == Å & degrees), Å -> m & degrees -> rad
		self.a=parlist[19]*1e-10
		self.b=parlist[20]*1e-10
		self.c=parlist[21]*1e-10
		self.alpha=deg2rad(parlist[22])
		self.beta=deg2rad(parlist[23])
		self.gamma=deg2rad(parlist[24])
		# Orienting vectors
		self.Avec=array([parlist[25],parlist[26],parlist[27]])
		self.Bvec=array([parlist[28],parlist[29],parlist[30]])
		# Excitation parameters
		self.Qvec=array([parlist[31],parlist[32],parlist[33]])	# Excitation wavevector
		self.En=parlist[34]*eV_to_J*1e-3						# Energy transfer (input unit == meV), meV -> J
		self.dEdq=parlist[38]*eV_to_J/hbar*1e-13				# Local slope of the dispersion (input unit == meV.Å), meV.Å -> m/s
		self.dE=parlist[85]*eV_to_J*1e-3						# Energy detuning parameter
		# Wavectors & excitation velocity
		self.Gmodvec=array([parlist[35],parlist[36],parlist[37]])	# Direction of the dispersion
		self.qvec=array([parlist[39],parlist[40],parlist[41]])		# Q = G + q
		self.Gvec=self.Qvec-self.qvec 								# Zone center (Bragg)
		# Sample distribution of lattice constants
		if parlist[5]<1e-8: 
			self.deltad=1e-8	# Set delta d/d very small if zero
		else:
			self.deltad=parlist[5]
		# Is there an in-pile neutron guide?	
		self.guide=parlist[58]
		if self.guide == 0:
			self.guidehdiv=deg2rad(parlist[59]*FWHMminstoSigdegs)	# Input unit == arc minutes, FWHM), arcmin -> rad	
			self.guidevdiv=deg2rad(parlist[60]*FWHMminstoSigdegs)	# Input unit == arc minutes, FWHM), arcmin -> rad	
		# Spatial parameters (Popovici)
		if self.method == 'pop':
			# Source shape
			self.sourceshape=parlist[55]
			if self.sourceshape==0:
				self.sourcediameter=parlist[56]*1e-2	# Input unit == cm, cm -> m
			else:
				self.sourcey=parlist[56]*1e-2			# Input unit == cm, cm -> m
				self.sourcez=parlist[57]*1e-2			# Input unit == cm, cm -> m		
			# Sample shape
			self.sampleshape=parlist[61]
			if self.sampleshape == 0:
				self.samplediameter=parlist[62]*1e-2	# Input unit == cm, cm -> m
				self.sampleheight=parlist[63]*1e-2		# Input unit == cm, cm -> m
			else:
				self.samplethickness=parlist[62]*1e-2	# Input unit == cm, cm -> m
				self.samplewidth=parlist[63]*1e-2		# Input unit == cm, cm -> m
				self.sampleheight=parlist[64]*1e-2		# Input unit == cm, cm -> m
			# Detector shape
			self.detshape=parlist[65]
			if self.detshape == 0:
				self.detdiameter=parlist[66]*1e-2		# Input unit == cm, cm -> m
			else:
				self.detwidth=parlist[66]*1e-2			# Input unit == cm, cm -> m
				self.detheight=parlist[67]*1e-2			# Input unit == cm, cm -> m
			# Monochromator geometry
			self.monothickness=parlist[68]*1e-2			# Input unit == cm, cm -> m
			self.monowidth=parlist[69]*1e-2				# Input unit == cm, cm -> m
			self.monoheight=parlist[70]*1e-2			# Input unit == cm, cm -> m
			# Analyzer geometry			
			self.anathickness=parlist[71]*1e-2			# Input unit == cm, cm -> m
			self.anawidth=parlist[72]*1e-2				# Input unit == cm, cm -> m
			self.anaheight=parlist[73]*1e-2				# Input unit == cm, cm -> m
			# Distances [cm]			
			self.L0=parlist[74]*1e-2					# Input unit == cm, cm -> m
			self.L1=parlist[75]*1e-2					# Input unit == cm, cm -> m
			self.L2=parlist[76]*1e-2					# Input unit == cm, cm -> m
			self.L3=parlist[77]*1e-2					# Input unit == cm, cm -> m
			# Monoch/analyzer curvatures
			self.rmh=parlist[78]						# Input unit == 1/m
			self.rmv=parlist[79] 						# Input unit == 1/m
			self.rah=parlist[80] 						# Input unit == 1/m
			self.rav=parlist[81]	 					# Input unit == 1/m
		
	def calcAngles(self):
# 		Calculates TAS parameters (crystal angles, scattering vector basis, etc.)
	
		HKLtoRec=self.getReciprocalBasis()
		Avec=dot(HKLtoRec,self.Avec)
		A=norm(Avec)
		Bvec=dot(HKLtoRec,self.Bvec)
		Cvec=cross(Avec,Bvec)
		Qvec=dot(HKLtoRec,self.Qvec)
		Q=norm(Qvec)
		qvec=dot(HKLtoRec,self.qvec)
		q=norm(qvec)
		Gvec=dot(HKLtoRec,self.Gvec)
		G=norm(Gvec)
		
		ki,kf=self.calcKiKf()
		dM=self.dM
		dA=self.dA
		A1=self.SM*arcsin(pi/(ki*dM))
		A2=2*A1
		if (ki**2+kf**2-Q**2)/(2*ki*kf) > 1 or (ki**2+kf**2-Q**2)/(2*ki*kf) <-1:
			if version_info[:3]<(3,0):
				print '\nError: scattering triangle does not close!!!'
			else:
				print('\nError: scattering triangle does not close!!!')
			return
		else:
			A4=self.SS*arccos((ki**2+kf**2-Q**2)/(2*ki*kf))
		A5=self.SA*arcsin(pi/(kf*dA))
		A6=2*A5
		
		phi_S=-self.SS*arccos((Q**2+ki**2-kf**2)/(2*Q*ki))
		sigma_S=-self.SS*arccos((Q**2+kf**2-ki**2)/(2*Q*kf))
				
		if Q == G or q == 0:
			phi_G=0		
		else:
			phiarg=float((q**2+Q**2-G**2)/(2*q*Q))
			if abs(phiarg) > 1:
				phiarg=sign(phiarg)*1
			phi_G=arccos(phiarg)	
			
		if dot(cross(Gvec,qvec),Cvec)<0:
			phi_G*=-1
			
		gamma=arccos(dot(Qvec,Avec)/(Q*A))		
		if dot(cross(Qvec,Avec),Cvec)<0:
			gamma*=-1
		
		Qvec=array([Q,0,0])*hbar/m_n		
		Avec=A*array([cos(gamma),sin(gamma),0])	
		kivec=ki*array([cos(phi_S),-sin(phi_S),0])*hbar/m_n
		kfvec=kf*array([-cos(sigma_S),-sin(sigma_S),0])*hbar/m_n
		qvec=q*array([cos(phi_G),sin(phi_G),0])*hbar/m_n		
				
		A3=arccos(dot(kivec,Avec)/(norm(kivec)*norm(Avec)))
		if dot(cross(kivec,Avec),array([0,0,1]))<0:
			A3*=-1
		
		return A1,A2,A3,A4,A5,A6,phi_S,sigma_S
	
	def calcKiKf(self):	
# 		Calculates ki and kf (in AA^-1) knowing fixed k and energy transfer values
	
		if self.fix==1:
			ki=self.kfix
			kf=sqrt(ki**2-2*m_n*self.En/hbar**2)
		else:
			kf=self.kfix
			ki=sqrt(2*m_n*self.En/hbar**2+kf**2)
			
		return ki,kf	
			
	def getReciprocalBasis(self):
		
		Gtensor=zeros([3,3]) # Metric tensor 
		Gtensor[0,0]=self.a**2
		Gtensor[0,1]=self.a*self.b*cos(self.gamma)
		Gtensor[0,2]=self.a*self.c*cos(self.beta)
		Gtensor[1,0]=self.a*self.b*cos(self.gamma)
		Gtensor[1,1]=self.b**2
		Gtensor[1,2]=self.b*self.c*cos(self.alpha)
		Gtensor[2,0]=self.a*self.c*cos(self.beta)
		Gtensor[2,1]=self.b*self.c*cos(self.alpha)
		Gtensor[2,2]=self.c**2
		
		V=sqrt(det(Gtensor)) # Unit cell volume
		
		astar=2*pi*self.b*self.c*sin(self.alpha)/V # Reciprocal lattice parameters
		bstar=2*pi*self.a*self.c*sin(self.beta)/V
		cstar=2*pi*self.a*self.b*sin(self.gamma)/V
		alphastar=arccos((cos(self.beta)*cos(self.gamma)-cos(self.alpha)) \
			/(sin(self.beta)*sin(self.gamma)));
		betastar=arccos((cos(self.alpha)*cos(self.gamma)-cos(self.beta)) \
			/(sin(self.alpha)*sin(self.gamma)));
		gammastar=arccos((cos(self.beta)*cos(self.alpha)-cos(self.gamma)) \
			/(sin(self.beta)*sin(self.alpha)));
		
		self.astar=astar
		self.bstar=bstar
		self.cstar=cstar
		self.alphastar=alphastar
		self.betastar=betastar
		self.gammastar=gammastar		
		
		astarvec=array([astar,0,0]) # Tranformation matrix [h,k,l] -> cartesian [a*,b*,c*]
		bstarvec=array([bstar*cos(gammastar),bstar*sin(gammastar),0])
		cstarvec1=cstar*cos(betastar)
		cstarvec2=cstar*(cos(alphastar)-cos(betastar)*cos(gammastar)) \
			/sin(gammastar)
		cstarvec3=sqrt(cstar**2-cstarvec1**2-cstarvec2**2)
		cstarvec=array([cstarvec1,cstarvec2,cstarvec3])
		
		HKLtoRec=array([astarvec,bstarvec,cstarvec]).T
		return HKLtoRec

	def getOrthonormalBasis(self):
	
		HKLtoRec=self.getReciprocalBasis()
		V1=dot(HKLtoRec,self.Avec)
		V2=dot(HKLtoRec,self.Bvec)
		V3=cross(V1,V2)
		V2=cross(V3,V1)
		V3/=norm(V3)
		V2/=norm(V2)
		V1/=norm(V1)
		
		return V1,V2,V3	

	def scalarProd(self,x,y,a,b,c,aa,bb,cc):
# 		Scalar product within crystal system
		s=x[0]*y[0]*a**2+x[1]*y[1]*b**2+x[2]*y[2]*c**2
		s+=(x[0]*y[1]+y[0]*x[1])*a*b*cos(cc)
		s+=(x[0]*y[2]+y[0]*x[2])*a*c*cos(bb)
		s+=(x[2]*y[1]+y[2]*x[1])*c*b*cos(aa)
		return s
	
	def computeCNmatrix(self):	
# 		Determines Cooper-Nathans resolution matrix for host TAS spectrometer
#		using Popovici's compact formalism.
		
		A1,A2,A3,A4,A5,A6,phi_S,sigma_S=self.calcAngles()
		ki,kf=self.calcKiKf()
		ki*=1e-10
		kf*=1e-10
		
		Amat=zeros((6,8))
		Amat[0,0]=0.5*ki/tan(A1)
		Amat[0,1]=-Amat[0,0]
		Amat[1,1]=ki
		Amat[2,3]=ki
		Amat[3,4]=0.5*kf/tan(A5)
		Amat[3,5]=-Amat[3,4]
		Amat[4,4]=kf
		Amat[5,6]=kf

		Bmat=zeros((4,6))
		Bmat[0,0]=cos(phi_S)
		Bmat[0,1]=sin(phi_S)
		Bmat[0,3]=-cos(phi_S-A4)
		Bmat[0,4]=-sin(phi_S-A4)
		Bmat[1,0]=-Bmat[0,1]
		Bmat[1,1]=Bmat[0,0]
		Bmat[1,3]=-Bmat[0,4]
		Bmat[1,4]=Bmat[0,3]
		Bmat[2,2]=1
		Bmat[2,5]=-1
		#Bmat[3,0]=hbar*ki/m_n
		#Bmat[3,3]=-hbar*kf/m_n
		Bmat[3,0]=2*2.07654*ki
		Bmat[3,3]=-2*2.07654*kf
		
		Cmat=zeros((4,8))
		Cmat[0,0]=0.5
		Cmat[0,1]=0.5
		Cmat[1,2]=1/(2*sin(A1))
		Cmat[1,3]=-1/(2*sin(A1))
		Cmat[2,4]=0.5
		Cmat[2,5]=0.5
		Cmat[3,6]=1/(2*sin(A5))
		Cmat[3,7]=-1/(2*sin(A5))

		Fmat=zeros((4,4))
		Fmat[0,0]=1/self.etaM**2        
		Fmat[1,1]=1/self.etaM**2 # should be etaM_v
		Fmat[2,2]=1/self.etaA**2 
		Fmat[3,3]=1/self.etaA**2 # should be etaA_v
		
		Gmat=zeros([8,8])
		if self.guide == 0:
			Gmat[0,0]=1/(self.guidehdiv*2*pi/ki)**2
			Gmat[2,2]=1/(self.guidevdiv*2*pi/ki)**2
		else:
			Gmat[0,0]=1/self.alpha0**2
			Gmat[2,2]=1/self.beta0**2
		Gmat[1,1]=1/self.alpha1**2
		Gmat[3,3]=1/self.beta1**2		
		Gmat[4,4]=1/self.alpha2**2
		Gmat[5,5]=1/self.alpha3**2
		Gmat[6,6]=1/self.beta2**2
		Gmat[7,7]=1/self.beta3**2	

		# Calculate Cooper-Nathans resolution matrix
		# M-1 = B A (G + C'FC)-1 A' B';
		MCNinv=dot(dot(Bmat,dot(dot(Amat,inv(dot(transpose(Cmat),dot(Fmat,Cmat))+Gmat)),transpose(Amat))),transpose(Bmat))
		
		# Sample mosaic contribution (Werner & Pynn)
		HKLtoRec=self.getReciprocalBasis()
		Qvec=dot(HKLtoRec,self.Qvec)
		Q=norm(Qvec)*1e-10
		MCNinv[1,1]=MCNinv[1,1]+(self.etaS*Q)**2
		MCNinv[2,2]=MCNinv[2,2]+(self.etaS*Q)**2
		
		# Rotate towards (Qx,Qy,Qz,hw)-frame, was (X1,X2,X3,X4)
		V1,V2,V3=self.getOrthonormalBasis()
		G=norm(self.Gvec)
		Q=norm(self.Qvec)
		q=norm(self.qvec)
		if Q == G or q == 0:
			phi_G=0		
		else:
			phiarg=float((q**2+Q**2-G**2)/(2*q*Q))
			if abs(phiarg) > 1:
				phiarg=sign(phiarg)*1
			phi_G=arccos(phiarg)		
		if dot(cross(self.Gvec,self.qvec),V3)<0:
			phi_G*=-1
		Rmat=zeros([4,4])
		Rmat[0,0]=cos(phi_G)
		Rmat[0,1]=sin(phi_G)
		Rmat[1,0]=-sin(phi_G)
		Rmat[1,1]=cos(phi_G)
		Rmat[2,2]=1
		Rmat[3,3]=1
		MCNinv=inv(dot(Rmat,dot(inv(MCNinv),Rmat.T)))
						
		# Return CN matrix
		MCN=inv(MCNinv)
		
		print '\nCooper-Nathans matrix :\n-----------------------------------\n'
		print MCN
		return MCN

	def calcTASellipsoid(self,method):
				
		if method=='cn':
			
			MCN=self.computeCNmatrix() # get Cooper-Nathans matrix			
			
			### (Qx,E)-plane
			HKLtoRec=self.getReciprocalBasis()
			Q=norm(dot(HKLtoRec,self.Qvec))*1e-10
			En=self.En/(eV_to_J*1e-3)
			dE=self.dE/(eV_to_J*1e-3)
			th=linspace(0,2*pi,1e4+1)
			# 50 % ellipse coordinates
			r=sqrt(2*log(2)/(MCN[3,3]*(sin(th))**2+MCN[0,0]*(cos(th))**2+2*MCN[0,3]*sin(th)*cos(th))) # 50% ellipsoid		
			xel50=r*cos(th)+Q		
			yel50=r*sin(th)+En		
			# 30 % ellipse coordinates
			r=sqrt(2*log(3)/(MCN[3,3]*(sin(th))**2+MCN[0,0]*(cos(th))**2+2*MCN[0,3]*sin(th)*cos(th))) # 33% ellipsoid		
			xel30=r*cos(th)+Q		
			yel30=r*sin(th)+En
			# 20 % ellipse coordinates (to get comfortable space between ellipse and figure borders)
			r=sqrt(2*log(5)/(MCN[3,3]*(sin(th))**2+MCN[0,0]*(cos(th))**2+2*MCN[0,3]*sin(th)*cos(th))) # 20% ellipsoid
			xel20=r*cos(th)+Q			
			yel20=r*sin(th)+En	
			xmin=min(xel20)
			xmax=max(xel20)
			ymin=min(yel20)
			ymax=max(yel20)				
			# Local dispersion curve
			xvel=linspace(xmin,xmax,101)
			#yvel=abs(xvel)*self.dEdq/rlutoinvA+self.Egap							
			#yvel=sign(xcen)*(xvel-xcen)*self.dEdq/rlutoinvA+En-dE
			yvel=zeros(101)
			# (Qx,E)-widths
			qxwidth=float(sqrt(8*log(2))/sqrt(MCN[0,0]))
			qywidth=float(sqrt(8*log(2))/sqrt(MCN[1,1]))	
			qzwidth=float(sqrt(8*log(2))/sqrt(MCN[2,2]))
			ewidth=float(sqrt(8*log(2))/sqrt(MCN[3,3]))
			#
			xwdthline=[En,En]
			xwdthlinex=[Q-0.5*qxwidth,Q+0.5*qxwidth]
			ywdthlinex=[Q,Q]
			ywdthline=[En-0.5*ewidth,En+0.5*ewidth]			
		
		else:
			print '\nError: unknown calculation method!!!'
			return
		
		return xmin,xmax,ymin,ymax,xvel,yvel,xwdthlinex,xwdthline,ywdthlinex,ywdthline,xel30,yel30,xel50,yel50,qxwidth,qywidth,qzwidth,ewidth