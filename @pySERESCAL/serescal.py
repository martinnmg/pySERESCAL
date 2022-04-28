#-*- encoding: ISO-8859-1 -*-
# N. MARTIN 2022/04/28
# (c) Laboratoire Léon Brillouin LLB
# From the 'SERESCAL' MATLAB toolbox by K. Habicht (HZB) 
# References:
# [Habicht2003] K. Habicht et al., J. Appl. Cryst. 36, 1307-1318 (2003)
# [Groitl2018] F. Groitl et al., J. Appl. Cryst. 51, 818-830 (2018)	

from sys import version_info
from numpy import arccos,arcsin,argmin,array,cos,cross,deg2rad,dot,exp,inf, \
	linspace,log,minimum,ndarray,pi,ones,rad2deg,sin,sign,sqrt,tan,transpose,zeros
from numpy.linalg import det,inv,norm
from scipy.constants import hbar,m_n
from scipy.constants import e as eV_to_J
from scipy.integrate import quad
from math import isnan

class SERESCAL:

	def __init__(self,parlist,method):
		# Creates a 'SERESCAL' object containing information & methods for all NRSE-related calculations
		self.updateParams(parlist,method)
		self.printSEparams()
		
	def updateParams(self,parlist,method):
# 		Gets TAS & NRSE parameters from the main panel,
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
		# Wavectors & excitation velocity
		self.Gmodvec=array([parlist[35],parlist[36],parlist[37]])	# Direction of the dispersion
		self.qvec=array([parlist[39],parlist[40],parlist[41]])		# Q = G + q
		self.Gvec=self.Qvec-self.qvec 								# Zone center (Bragg)
		# Sample distribution of lattice constants
		if parlist[5]<1e-8: 
			self.deltad=1e-8	# Set delta d/d very small if zero
		else:
			self.deltad=parlist[5]
		# Dispersion curvature matrix (input units depend on chosen mode), to be transformed later
		self.H=zeros([3,3])
		self.H[0,0]=parlist[42]
		self.H[0,1]=parlist[43]
		self.H[0,2]=parlist[44]
		self.H[1,0]=parlist[45]
		self.H[1,1]=parlist[46]
		self.H[1,2]=parlist[47]
		self.H[2,0]=parlist[48]
		self.H[2,1]=parlist[49]
		self.H[2,2]=parlist[50]
		self.cumbasis=parlist[86]
		# NRSE parameters
		self.fmin=parlist[51]*1e3	# Minimum RF frequency
		self.fmax=parlist[52]*1e3	# Maximum RF frequency
		self.L4pi=parlist[53]		# Flipper distance (normal)
		self.L8pi=parlist[54]		# Flipper distance (bootstrap)
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
				self.samplethickness=parlist[63]*1e-2	# Input unit == cm, cm -> m
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
		# Magnetic stuff (work in progress...)
		self.extype=parlist[82]							# Magnon/phonon
		self.Myy=parlist[83]							# In-plane moment (Myy + Mzz = 1)
		self.Mzz=parlist[84]							# Out-of-plane moment (Myy + Mzz = 1)

	def calcAngles(self):
# 		Calculates TAS (crystal angles, scattering vector basis, etc.)
#		and NRSE parameters (coil tilt angles, field integral ratio, etc.)
	
		# self.getReciprocalBasis()
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
				
		A3=arccos(dot(kivec,Avec)/(ki*A))
		if dot(cross(kivec,Avec),array([0,0,1]))<0:
			A3*=-1
		
		dEdq=self.dEdq
		if norm(qvec) == 0:
			dEdqvec=dEdq*Qvec/norm(Qvec)
		else:
			dEdqvec=dEdq*qvec/norm(qvec)		
		
		nivec=kivec-dEdqvec
		nfvec=kfvec-dEdqvec
		ni=norm(nivec)
		nf=norm(nfvec)
		self.ni=ni
		self.nf=nf	

		# RF coils tilt angles
		theta1=arccos(dot(kivec,nivec)/(norm(kivec)*norm(nivec)))
		if isnan(theta1) == 1:
			theta1=0
		if dot(cross(kivec,nivec),array([0,0,1]))<0:
			theta1*=-1
		theta2=arccos(dot(kfvec,nfvec)/(norm(kfvec)*norm(nfvec)))
		if isnan(theta2) == 1:
			theta2=0
		if dot(cross(kfvec,nfvec),array([0,0,1]))<0:
			theta2*=-1		

		fratio=(ni*ki**2*cos(theta1))/(nf*kf**2*cos(theta2))
		
		fmin=self.fmin
		fmax=self.fmax
		L4pi=self.L4pi
		L8pi=self.L8pi
		
		if fratio > 1: # f1 > f2
			if norm(qvec) == 0:	
				tmin=m_n/hbar*4*pi*fmin*L4pi/(kf**2)
				tmax=m_n/hbar*8*pi*fmax*L8pi/(ki**2)			
			else:
				tmin=m_n/hbar*4*pi*fmin*L4pi/(nf*kf**2*cos(theta2))
				tmax=m_n/hbar*8*pi*fmax*L8pi/(ni*ki**2*cos(theta1))
		else:
			if norm(qvec) == 0:	
				tmin=m_n/hbar*4*pi*fmin*L4pi/(ki**2)
				tmax=m_n/hbar*8*pi*fmax*L8pi/(kf**2)				
			else:
				tmin=m_n/hbar*4*pi*fmin*L4pi/(ni*ki**2*cos(theta1))
				tmax=m_n/hbar*8*pi*fmax*L8pi/(nf*kf**2*cos(theta2))	
		
		return A1,A2,A3,A4,A5,A6,theta1,theta2,ni,nf,fratio,tmin,tmax,phi_S,sigma_S
	
	def printSEparams(self):
# 		Prints TAS & NRSE parameters in the command line window
		
		A1,A2,A3,A4,A5,A6,theta1,theta2,ni,nf,fratio,tmin,tmax,phi_S,sigma_S=self.calcAngles()
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
			#
			print '\nNRSE parameters:\n----------------'
			print 'Zeta_i =',round(rad2deg(theta1),2),'deg'
			print 'Zeta_f =',round(rad2deg(theta2),2),'deg'
			print 'Frequency ratio f_i/f_f =',round(fratio,3)
			print 'Minimum tau =',round(tmin*1e12,1),'ps'
			print 'Maximum tau =',round(tmax*1e12,1),'ps'
		else:
			print('\nTAS parameters:\n---------------')
			if self.fix == 1:
				print('Q = [',self.Qvec[0],self.Qvec[1],self.Qvec[2],'], hw =',self.En,'meV (fixed ki)')	
			else:
				print('Q = [',self.Qvec[0],self.Qvec[1],self.Qvec[2],'], hw =',self.En,'meV (fixed kf)')			
			print('a1 =',round(rad2deg(A1),2),'deg, a2 =',round(rad2deg(A2),2),'deg')		
			print('a3 =',round(rad2deg(A3),2),'deg, a4 =',round(rad2deg(A4),2),'deg')		
			print('a5 =',round(rad2deg(A5),2),'deg, a6 =',round(rad2deg(A6),2),'deg\n')	
			#
			print('\nNRSE parameters:\n----------------')
			print('Zeta_i =',round(rad2deg(theta1),2),'deg')
			print('Zeta_f =',round(rad2deg(theta2),2),'deg')
			print('Frequency ratio f_i/f_f =',round(fratio,3))
			print('Minimum tau =',round(tmin*1e12,1),'ps')
			print('Maximum tau =',round(tmax*1e12,1),'ps')		
	
	def getReciprocalBasis(self):
		
		# Metric tensor
		Gtensor=zeros([3,3]) 
		Gtensor[0,0]=self.a**2
		Gtensor[0,1]=self.a*self.b*cos(self.gamma)
		Gtensor[0,2]=self.a*self.c*cos(self.beta)
		Gtensor[1,0]=self.a*self.b*cos(self.gamma)
		Gtensor[1,1]=self.b**2
		Gtensor[1,2]=self.b*self.c*cos(self.alpha)
		Gtensor[2,0]=self.a*self.c*cos(self.beta)
		Gtensor[2,1]=self.b*self.c*cos(self.alpha)
		Gtensor[2,2]=self.c**2
		
		# Reciprocal lattice constants
		astar=2*pi*self.b*self.c*sin(self.alpha)/sqrt(det(Gtensor))	
		bstar=2*pi*self.a*self.c*sin(self.beta)/sqrt(det(Gtensor))
		cstar=2*pi*self.a*self.b*sin(self.gamma)/sqrt(det(Gtensor))
		# Reciprocal lattice angles
		alphastar=arccos((cos(self.beta)*cos(self.gamma)-cos(self.alpha)) \
			/(sin(self.beta)*sin(self.gamma)));
		betastar=arccos((cos(self.alpha)*cos(self.gamma)-cos(self.beta)) \
			/(sin(self.alpha)*sin(self.gamma)));
		gammastar=arccos((cos(self.beta)*cos(self.alpha)-cos(self.gamma)) \
			/(sin(self.beta)*sin(self.alpha)));
			
		astarvec=array([astar,0,0])	# Reciprocal lattice vectors
		bstarvec=array([bstar*cos(gammastar),bstar*sin(gammastar),0])
		cstarvec1=cstar*cos(betastar)
		cstarvec2=cstar*(cos(alphastar)-cos(betastar)*cos(gammastar)) \
			/sin(gammastar)
		cstarvec3=sqrt(cstar**2-cstarvec1**2-cstarvec2**2)
		cstarvec=array([cstarvec1,cstarvec2,cstarvec3])
		
		HKLtoRec=array([astarvec,bstarvec,cstarvec]).T	# Tranformation matrix [h,k,l] -> cartesian [a*,b*,c*]
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
						
	def calcKiKf(self):	
# 		Calculates ki and kf (in AA^-1) knowing fixed k and energy transfer values
	
		if self.fix==1:
			ki=self.kfix
			kf=sqrt(ki**2-2*m_n*self.En/hbar**2)
		else:
			kf=self.kfix
			ki=sqrt(2*m_n*self.En/hbar**2+kf**2)
			
		return ki,kf		

	def computeResolutionCurves(self,tau,method):
# 		Computes resolution matrices, taking Cooper-Nathans approximation for the TAS transmission function.
	
		# Instrumental (flat dispersion and perfect sample)
		if method=='cn':
			self.computeInstrumentalResolutionMatrix_CN(tau)
		elif method=='pop':
			self.computeInstrumentalResolutionMatrix_POP(tau)
		yinstr=zeros(tau.size)
		for k in xrange(tau.size):
			yinstr[k]=sqrt(1./abs(det(self.Linstmat[1:,1:,k])))
		yinstr*=sqrt(abs(det(self.Linstmat[1:,1:,0])))
		# Curvature (perfect sample)
		self.computeCurvatureResolutionMatrix(tau)
		ycurv=zeros(tau.size)
		for k in xrange(tau.size):
			ycurv[k]=sqrt(1./abs(det(self.Lcurvmat[1:,1:,k])))
		ycurv*=sqrt(abs(det(self.Lcurvmat[1:,1:,0])))
		# Sample imperfections (flat dispersion)
		self.computeSampleImperfResolutionMatrix(tau)
		yimperf=zeros(tau.size)
		for k in xrange(tau.size):
			yimperf[k]=sqrt(1./abs(det(self.Limperfmat[1:,1:,k])))*self.Limperflinterm[k]
		yimperf*=sqrt(abs(det(self.Limperfmat[1:,1:,0])))
		# Total resolution (curved dispersion and imperfect sample)
		self.computeTotalResolutionMatrix(tau,method)
		ytot=zeros(tau.size)
		for k in xrange(tau.size):
			ytot[k]=sqrt(1./abs(det(self.Ltotmat[1:,1:,k])))*self.Ltotlinterm[k]
		ytot*=sqrt(abs(det(self.Ltotmat[1:,1:,0])))
		
		return yinstr,ycurv,yimperf,ytot

	def computeInstrumentalResolutionMatrix_CN(self,tau): 
# 		Calculates NRSE resolution matrix taking Cooper-Nathans approximation
# 		for the TAS transmission function. Assumes flat dispersion and perfect sample.
		
		A1,A2,A3,A4,A5,A6,theta1,theta2,ni,nf,fratio,tmin,tmax,phi_S,sigma_S=self.calcAngles()
		ki,kf=self.calcKiKf()
		
		## Real part (TAS transmission) 
		
		# Cooper-Nathans "a" constants		
		a1=tan(A1)/self.etaM/ki
		a2=1/self.etaM/ki
		a3=1/self.alpha1/ki
		a4=1/self.alpha2/kf
		a5=tan(A5)/self.etaA/kf
		a6=-1/self.etaA/kf
		if self.guide == 0:
			a7=2*tan(A1)/(self.guidehdiv*2*pi/(ki*1e-10))/ki
			a8=1/(self.guidehdiv*2*pi/(ki*1e-10))/ki		
			a11=sqrt(1/(4*(sin(A1))**2*self.etaM**2*ki**2+(self.guidevdiv*2*pi/(ki*1e-10))**2*ki**2)+1/(self.beta1**2*ki**2))
		else:
			a7=2*tan(A1)/self.alpha0/ki
			a8=1/self.alpha0/ki
			a11=sqrt(1/(4*(sin(A1))**2*self.etaM**2*ki**2+self.beta0**2*ki**2)+1/(self.beta1**2*ki**2))
		a9=2*tan(A5)/self.alpha1/kf
		a10=-1/self.alpha3/kf
		a12=sqrt(1/(4*(sin(A5))**2*self.etaA**2*kf**2+self.beta3**2*kf**2)+1/(self.beta2**2*kf**2))
		
		# Cooper-Nathans "b" constants	
		b0=a1*a2+a7*a8
		b1=a2**2+a3**2+a8**2
		b2=a4**2+a6**2+a10**2
		b3=a5**2+a9**2
		b4=a5*a6+a9*a10
		b5=a1**2+a7**2
		
		# "XI" matrix - 6x6, symmetric (cf. [Habicht2003])
		XImat=zeros([6,6],dtype=complex)
		XImat[0,0]=b3/(cos(theta2)*nf)**2
		XImat[0,1]=-b3*ni/(cos(theta2)*nf)**2
		XImat[0,3]=b3*tan(theta2)/(cos(theta2)*nf)-b4/(cos(theta2)*nf)
		XImat[1,1]=b5/cos(theta1)**2+b3*(ni/(cos(theta2)*nf))**2
		XImat[1,2]=-b5*tan(theta1)/cos(theta1)+b0/cos(theta1)
		XImat[2,1]=XImat[1,2]
		XImat[1,3]=-b3*tan(theta2)*ni/(cos(theta2)*nf)+b4*ni/(cos(theta2)*nf) # Mistake in [Habicht2003], Eq. (56)
		XImat[3,1]=XImat[1,3]
		XImat[2,2]=b5*tan(theta1)**2+b1-2*b0*tan(theta1) # Mistake in [Habicht2003], Eq. (52)
		XImat[3,3]=-2*b4*tan(theta2)+b2+b3*tan(theta2)**2
		XImat[4,4]=a11**2
		XImat[5,5]=a12**2
		
		## Imaginary part (Larmor precession)
		
		# "PSI" matrix - 6x6, symmetric (cf. [Habicht2003])
		PSImat=zeros([6,6],dtype=complex)
		PSImat[0,0]=hbar/m_n/(nf*cos(theta2))**2+2/(nf*kf*cos(theta2))
		PSImat[0,1]=-hbar*ni/m_n/(nf*cos(theta2))**2-2*ni/(nf*kf*cos(theta2))
		PSImat[0,3]=hbar*tan(theta2)/(m_n*cos(theta2)*nf)
		PSImat[1,2]=hbar*tan(theta1)/(m_n*cos(theta1))
		PSImat[2,1]=PSImat[1,2]
		PSImat[1,3]=-hbar*tan(theta2)*ni/(m_n*cos(theta2)*nf)
		PSImat[3,1]=PSImat[1,3]
		PSImat[1,1]=-hbar/m_n/cos(theta1)**2+hbar/m_n*(ni/(nf*cos(theta2)))**2-2*ni/(ki*cos(theta1))+2*ni**2/(nf*kf*cos(theta2))
		PSImat[2,2]=-hbar/m_n/cos(theta1)**2
		PSImat[3,3]=hbar/m_n/cos(theta2)**2
		PSImat[4,4]=-hbar/m_n
		PSImat[5,5]=hbar/m_n
		
		## Instrumental resolution matrix
		
		self.Linstmat=ndarray([6,6,tau.size],dtype=complex)
		for k in xrange(tau.size):
			self.Linstmat[:,:,k]=XImat+1j*tau[k]*PSImat

	def computeInstrumentalResolutionMatrix_POP(self,tau): 
# 		Calculates NRSE resolution matrix taking Popovici approximation
# 		for the TAS transmission function. Assumes flat dispersion and perfect sample.
		
		A1,A2,A3,A4,A5,A6,theta1,theta2,ni,nf,fratio,tmin,tmax,phi_S,sigma_S=self.calcAngles()
			
		HKLtoRec=self.getReciprocalBasis()	
		ki,kf=self.calcKiKf()
			
		## Real part (TAS transmission) 
		
		if self.sourceshape == 0:
			sourced=self.sourcediameter
		else:
			sourcey=self.sourcey
			sourcez=self.sourcez
		
		monoth=self.monothickness
		monow=self.monowidth
		monoh=self.monoheight
		
		if self.sampleshape == 0:
			samth=self.samplethickness
			samd=self.samplediameter
		else:
			samth=self.samplethickness
			samw=self.samplewidth
			samh=self.sampleheight
		
		if self.detshape == 0:
			detd=self.detdiameter
		else:
			dety=self.detwidth
			detz=self.detheight
		
		anath=self.anathickness
		anaw=self.anawidth
		anah=self.anaheight	
				
		AM=zeros([6,8])
		AM[0,0]=0.5*ki/tan(A1)
		AM[0,1]=-0.5*ki/tan(A1)
		AM[1,1]=ki
		AM[2,3]=ki
		AM[3,4]=0.5*kf/tan(A5)
		AM[3,5]=-0.5*kf/tan(A5)
		AM[4,4]=kf
		AM[5,6]=kf
		
		FM=zeros([4,4])
		FM[0,0]=1/self.etaM**2
		FM[1,1]=1/self.etaM**2 # Should be vertical mosaic
		FM[2,2]=1/self.etaA**2
		FM[3,3]=1/self.etaA**2 # Should be vertical mosaic
	
		GM=zeros([8,8])
		if self.guide == 0:
			GM[0,0]=1/(self.guidehdiv*2*pi/(ki*1e-10))**2
			GM[2,2]=1/(self.guidevdiv*2*pi/(ki*1e-10))**2
		else:
			GM[0,0]=1/self.alpha0**2
			GM[2,2]=1/self.beta0**2
		GM[1,1]=1/self.alpha1**2
		GM[3,3]=1/self.beta1**2		
		GM[4,4]=1/self.alpha2**2
		GM[5,5]=1/self.alpha3**2
		GM[6,6]=1/self.beta2**2
		GM[7,7]=1/self.beta3**2		
		
		CM=zeros([4,8])
		CM[0,0]=0.5
		CM[0,1]=0.5
		CM[2,4]=0.5
		CM[2,5]=0.5
		CM[1,2]=0.5/sin(A1)
		CM[1,3]=-0.5/sin(A1)
		CM[3,6]=0.5/sin(A5)
		CM[3,7]=-0.5/sin(A5)

		TM=zeros([4,13])
		TM[0,0]=-0.5/self.L0
		TM[0,2]=0.5*cos(A1)*(1/self.L1-1/self.L0)
		TM[0,3]=0.5*sin(A1)*(1/self.L0+1/self.L1-2*self.rmh*self.SM/abs(sin(A1)))
		TM[0,5]=0.5*sin(0.5*A4)/self.L1
		TM[0,6]=0.5*cos(0.5*A4)/self.L1
		TM[1,1]=-0.5/(self.L0*sin(A1))
		TM[1,4]=0.5*(1/self.L0+1/self.L1-2*self.rmv*sin(A1))/sin(abs(A1))
		TM[1,7]=-0.5/(self.L1*sin(A1))
		TM[2,5]=0.5*sin(0.5*A4)/self.L2
		TM[2,6]=-0.5*cos(0.5*A4)/self.L2
		TM[2,8]=0.5*cos(A5)*(1/self.L3-1/self.L2)
		TM[2,9]=0.5*sin(A5)*(1/self.L2+1/self.L3-2*self.rah*self.SA/abs(sin(A5)))
		TM[2,11]=0.5/self.L3
		TM[3,7]=-0.5/(self.L2*sin(A5))
		TM[3,10]=0.5*(1/self.L2+1/self.L3-2*self.rav*sin(A5))/abs(sin(A5))
		TM[3,12]=-0.5/(self.L3*sin(A5))
		
		DM=zeros([8,13])
		DM[0,0]=-1/self.L0
		DM[0,2]=-cos(A1)/self.L0
		DM[0,3]=sin(A1)/self.L0
		DM[2,1]=-1/self.L0
		DM[2,4]=1/self.L0
		DM[1,2]=cos(A1)/self.L1
		DM[1,3]=sin(A1)/self.L1
		DM[1,5]=sin(0.5*A4)/self.L1
		DM[1,6]=cos(0.5*A4)/self.L1
		DM[3,4]=-1/self.L1
		DM[3,7]=1/self.L1
		DM[4,5]=sin(0.5*A4)/self.L2
		DM[4,6]=-cos(0.5*A4)/self.L2
		DM[4,8]=-cos(A5)/self.L2
		DM[4,9]=sin(A5)/self.L2
		DM[6,7]=-1/self.L2
		DM[6,10]=1/self.L2
		DM[5,8]=cos(A5)/self.L3
		DM[5,9]=sin(A5)/self.L3
		DM[5,11]=1/self.L3
		DM[7,10]=-1/self.L3
		DM[7,12]=1/self.L3
		
		# Spatial covariances matrix
		Smat=zeros([13,13])
		
		# Source covariance matrix
		if self.sourceshape == 0:
			Smat[0,0]=1/16.*sourced**2
			Smat[1,1]=1/16.*sourced**2
		else:
			Smat[0,0]=1/12.*sourcey**2
			Smat[1,1]=1/12.*sourcez**2
		
		# Monochromator covariance matrix
		Smono=zeros([3,3])
		Smono[0,0]=1/12.*monoth**2
		Smono[1,1]=1/12.*monow**2
		Smono[2,2]=1/12.*monoh**2
		Smat[2:5,2:5]=Smono
		
		# Sample covariance matrix
		Ssam=zeros([3,3])
		if self.sampleshape == 0:
			Ssam[0,0]=1/16.*samd**2
			Ssam[1,1]=1/16.*samd**2
			Ssam[2,2]=1/12.*samth**2
		else:
			Ssam[0,0]=1/12.*samth**2
			Ssam[1,1]=1/12.*samw**2
			Ssam[2,2]=1/12.*samh**2		
		Smat[5:8,5:8]=Ssam
				
		# Analyzer covariance matrix
		Sana=zeros([3,3])
		Sana[0,0]=1/12.*anath**2
		Sana[1,1]=1/12.*anaw**2
		Sana[2,2]=1/12.*anah**2
		Smat[8:11,8:11]=Sana

		# Detector covariance matrix
		if self.detshape == 0:
			Smat[11,11]=1/16.*detd**2
			Smat[12,12]=1/16.*detd**2
		else:
			Smat[11,11]=1/12.*dety**2
			Smat[12,12]=1/12.*detz**2		
			
		# Builds full Popovici matrix and transform it to the coordinate space of the resolution matrix (Dw,Dkin,y1,y2,z1,z2)
		
		ICmat=zeros([6,6])
		ICmat[0,1]=1/cos(theta1)
		ICmat[0,2]=-tan(theta1)
		ICmat[3,0]=-1/(nf*cos(theta2))
		ICmat[3,1]=ni/(nf*cos(theta2))
		ICmat[3,3]=-tan(theta2)
		ICmat[1,2]=1
		ICmat[2,4]=1
		ICmat[4,3]=1
		ICmat[5,5]=1
			
		XImat=inv(inv(dot(DM,dot(inv(inv(Smat)+dot(transpose(TM),dot(FM,TM))),transpose(DM))))+GM)
		XImat=dot(AM,dot(XImat,transpose(AM)))
		XImat=inv(dot(inv(ICmat),dot(XImat,transpose(inv(ICmat)))))

		## Imaginary part (Larmor precession)
		
		# "PSI" matrix - 6x6, symmetric (cf. [Habicht2003])
		PSImat=zeros([6,6],dtype=complex)
		PSImat[0,0]=hbar/m_n/(nf*cos(theta2))**2+2/(nf*kf*cos(theta2))
		PSImat[0,1]=-hbar*ni/m_n/(nf*cos(theta2))**2-2*ni/(nf*kf*cos(theta2))
		PSImat[0,3]=hbar*tan(theta2)/(m_n*cos(theta2)*nf)
		PSImat[1,2]=hbar*tan(theta1)/(m_n*cos(theta1))
		PSImat[2,1]=PSImat[1,2]
		PSImat[1,3]=-hbar*tan(theta2)*ni/(m_n*cos(theta2)*nf)
		PSImat[3,1]=PSImat[1,3]
		PSImat[1,1]=-hbar/m_n/cos(theta1)**2+hbar/m_n*(ni/(nf*cos(theta2)))**2-2*ni/(ki*cos(theta1))+2*ni**2/(nf*kf*cos(theta2))
		PSImat[2,2]=-hbar/m_n/cos(theta1)**2
		PSImat[3,3]=hbar/m_n/cos(theta2)**2
		PSImat[4,4]=-hbar/m_n
		PSImat[5,5]=hbar/m_n
		
		## Instrumental resolution matrix
		
		self.Linstmat=ndarray([6,6,tau.size],dtype=complex)
		for k in xrange(tau.size):
			self.Linstmat[:,:,k]=XImat+1j*tau[k]*PSImat
					
	def computeCurvatureResolutionMatrix(self,tau):
# 		Calculates NRSE resolution matrix taking Cooper-Nathans approximation
# 		for the TAS transmission function. Includes dispersion curvature
#		and neglects sample imperfections.

		A1,A2,A3,A4,A5,A6,theta1,theta2,ni,nf,fratio,tmin,tmax,phi_S,sigma_S=self.calcAngles()

		# Input Hessian is by default given in meV.AA^2 and expressed in the reciprocal lattice basis (a*,b*,c*)
		H=self.H*eV_to_J*1e-3/hbar*1e-20	# meV.AA^2 -> m^2/s
		
		# H is transformed to the "orthonormalized" AB basis
		V1,V2,V3=self.getOrthonormalBasis()
		Avec=self.Avec/norm(self.Avec)
		Avec=self.Gvec/norm(self.Gvec)
		Bvec=cross(V3.T,Avec)
		Cvec=cross(Avec,Bvec)
		RectoAB=array([Avec.T,Bvec.T,Cvec.T])
		H=dot(RectoAB.T,dot(H,RectoAB))
		
		# H is transformed to the (Qx,Qy,Qz) basis
		#HKLtoRec=self.getReciprocalBasis()
		#Qvec=dot(HKLtoRec,self.Qvec)
		Qvec=self.Qvec
		phi_AQ=arcsin(norm(cross(Avec,Qvec))/(norm(Avec)*norm(Qvec)))	
		if dot(cross(Avec,Qvec),Cvec)<0:
			phi_AQ*=-1
		ABtoQ=zeros([3,3])
		ABtoQ[0,0]=cos(phi_AQ)
		ABtoQ[0,1]=sin(phi_AQ)
		ABtoQ[1,0]=-sin(phi_AQ)
		ABtoQ[1,1]=cos(phi_AQ)
		ABtoQ[2,2]=1
		H=dot(ABtoQ.T,dot(H,ABtoQ))		
		
		# H is transformed to the coordinate space of the resolution matrix (Dw,Dkin,y1,y2,z1,z2)
		THETACmat=zeros([3,6])
		THETACmat[0,0]=cos(phi_S)
		THETACmat[0,1]=-sin(phi_S)		# Sign correction [Groitl2018], Eq. (66)
		THETACmat[0,3]=-cos(sigma_S)	# Sign correction [Groitl2018], Eq. (66)
		THETACmat[0,4]=sin(sigma_S)		# Sign correction [Groitl2018], Eq. (66)
		THETACmat[1,0]=sin(phi_S)		# Sign correction [Groitl2018], Eq. (66)
		THETACmat[1,1]=cos(phi_S)		# Sign correction [Groitl2018], Eq. (66)
		THETACmat[1,3]=-sin(sigma_S)	# Sign correction [Groitl2018], Eq. (66)
		THETACmat[1,4]=-cos(sigma_S)	# Sign correction [Groitl2018], Eq. (66)
		THETACmat[2,2]=1
		THETACmat[2,5]=-1	
			
		ICmat=zeros([6,6])
		ICmat[0,1]=1/cos(theta1)
		ICmat[0,2]=-tan(theta1)
		ICmat[3,0]=-1/(nf*cos(theta2))
		ICmat[3,1]=ni/(nf*cos(theta2))
		ICmat[3,3]=-tan(theta2)
		ICmat[1,2]=1
		ICmat[2,4]=1
		ICmat[4,3]=1
		ICmat[5,5]=1
				
		H=dot(ICmat.T,dot(THETACmat.T,dot(H,dot(THETACmat,ICmat))))
		
		self.Lcurvmat=ndarray([6,6,tau.size],dtype=complex)
		for k in xrange(tau.size):
			self.Lcurvmat[:,:,k]=self.Linstmat[:,:,k]+1j*tau[k]*H		
	
	def computeSampleImperfResolutionMatrix(self,tau): 
# 		Calculates resolution matrix taking Cooper-Nathans approximation
# 		for the TAS transmission function. Neglects dispersion curvature
#		and includes for sample imperfections (mosaic and lattice spacing spreads)
		
		A1,A2,A3,A4,A5,A6,theta1,theta2,ni,nf,fratio,tmin,tmax,phi_S,sigma_S=self.calcAngles()
		ki,kf=self.calcKiKf()

		HKLtoRec=self.getReciprocalBasis()		
				
		Avec=dot(HKLtoRec,self.Avec)
		A=norm(Avec)
		Bvec=dot(HKLtoRec,self.Bvec)
		Cvec=cross(Avec,Bvec)
		
		Gmodvec=dot(HKLtoRec,self.Gmodvec)
		Cx=dot(Gmodvec/norm(Gmodvec),Avec/norm(Avec))*self.dEdq # m/s
		Cy=dot(Gmodvec/norm(Gmodvec),Bvec/norm(Bvec))*self.dEdq
		Cz=dot(Gmodvec/norm(Gmodvec),Cvec/norm(Cvec))*self.dEdq
		if norm(Cz)>1e-3:
			if version_info[:3]<(3,0):
				print '\nSlope is defined out of the scattering plane!!!'
			else:
				print('\nSlope is defined out of the scattering plane!!!')
		Gvec=dot(HKLtoRec,self.Gvec)
		G=norm(Gvec) # m^-1	
		
		Imat=zeros([6,9])
		Imat[0,0]=1
		Imat[1,1]=1
		Imat[2,2]=1
		Imat[3,3]=1
		Imat[4,4]=1
		Imat[5,5]=1
		Imat[0,6]=-Cx
		Imat[0,7]=-Cy*G
		Imat[0,8]=-Cz*G
		self.Imat=Imat
		
		Nmat=zeros([9,9])
		Nmat[6,6]=1/(self.deltad*G)**2
		Nmat[7,7]=1/self.etaS**2
		Nmat[8,8]=1/self.etaS**2 # Vertical and horizontal mosaic spreads are assumed to be equal
		self.Nmat=Nmat
		
		Wmat=zeros([9,9])
		Wmat[6,7]=Cy
		Wmat[7,6]=Cy
		Wmat[6,8]=Cz
		Wmat[8,6]=Cz
		Wmat[7,7]=-Cx*G # Correction from [Groitl2018], Eq. (55)
		Wmat[8,8]=-Cx*G # Correction from [Groitl2018], Eq. (55)
		self.Wmat=Wmat
		
		Tvec=array([0,0,0,0,0,Cx,Cy*G,Cz*G])
		self.Tvec=Tvec
		
		self.Limperfmat=ndarray([9,9,tau.size],dtype=complex)
		self.Limperflinterm=zeros(tau.size)		
		for k in xrange(tau.size):		
			self.Limperfmat[:,:,k]=dot(Imat.T,dot(self.Linstmat[:,:,k],Imat))+Nmat-1j*tau[k]*Wmat
			self.Limperflinterm[k]=abs(exp(-0.5*tau[k]**2*dot(Tvec.T,dot(inv(self.Limperfmat[1:,1:,k]),Tvec))))

	def computeTotalResolutionMatrix(self,tau,method):  
# 		Calculates resolution matrix taking Cooper-Nathans approximation
# 		for the TAS transmission function

		if method=='cn':
			self.computeInstrumentalResolutionMatrix_CN(tau)
		elif method=='pop':
			self.computeInstrumentalResolutionMatrix_POP(tau)
		self.computeCurvatureResolutionMatrix(tau)
		self.computeSampleImperfResolutionMatrix(tau)
				
		self.Ltotmat=ndarray([9,9,tau.size],dtype=complex)
		self.Ltotlinterm=zeros(tau.size)	
		for k in xrange(tau.size):
			self.Ltotmat[:,:,k]=dot(self.Imat.T,dot(self.Lcurvmat[:,:,k],self.Imat))+self.Nmat-1j*tau[k]*self.Wmat
			self.Ltotlinterm[k]=abs(exp(-0.5*tau[k]**2*dot(self.Tvec.T,dot(inv(self.Limperfmat[1:,1:,k]),self.Tvec))))

	def computeMagneticFactor(self,tau):
#		Calculates magnetic correction factor for a given excitation if defined
#		as a 'magnon'. Needs magnetic terms My (in-precession plane) and Mz (vertical). 
#		Work in progress...
		
		A1,A2,A3,A4,A5,A6,theta1,theta2,ni,nf,fratio,tmin,tmax,phi_S,sigma_S=self.calcAngles()
		ki,kf=self.calcKiKf()
		# ki distribution
		alpha0=self.alpha0*sqrt(8*log(2))
		etaM=self.etaM*sqrt(8*log(2))
		a1=A1*pi/180
		dki=1/sqrt((tan(a1)**2)*(4./alpha0**2+1./etaM**2))*0.0001
		#
		len=tau.size
		ymag_par=ones(len)
		ymag_antipar=ones(len)
		#
		mtot=sqrt(self.Myy**2+self.Mzz**2)
		if mtot ==0: # No magnetic contribution...
			ymag_par=zeros(len)
			ymag_antipar=zeros(len)
		else:
			anisoM=(self.Myy**2-self.Mzz**2)/(self.Myy**2+self.Mzz**2)
			icount=0
			for t in tau:
				cosavg=quad(lambda dk: calcPhiNSE(dk,t,ki,dki),0,inf)
				ymag_par[icount]=0.5*(1+anisoM)+0.5*(1-anisoM)*cosavg[0]
				ymag_antipar[icount]=0.5*(1+anisoM)*cosavg[0]-0.5*(1-anisoM)	
				icount+=1
		return ymag_par,ymag_antipar

def calcPhiNSE(dk,t,ki,dki):	
	return cos(4*6.29595*t*ki*dk)*exp(-4*log(2.)*(dk/dki)**2)
	
if __name__ == "__main__":
	
	fid=open('serescal.par','r')
	params=fid.readlines()
	fid.close()
	parlist=[float(i) for i in params]
	calc=SERESCAL(parlist)	