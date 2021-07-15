#-*- encoding: ISO-8859-1 -*-
# N. MARTIN 19/12/2019
# (c) Laboratoire Léon Brillouin LLB
# From the 'serescal' MATLAB toolbox by K. Habicht (HZB) 
# References:
# K. Habicht et al., J. Appl. Cryst. 36, 1307-1318 (2003)
# F. Groitl et al., J. Appl. Cryst. 51, 818-830 (2018)

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
		A1,A2,A3,A4,A5,A6,theta1,theta2,ni,nf,fratio,tmin,tmax,phi_S,sigma_S=self.calcAngles()		
		print '\nTAS parameters:\n---------------'
		if self.fix == 1:
			print 'Q = [',self.Qvec[0],self.Qvec[1],self.Qvec[2],'], hw =',self.En,'meV (fixed ki)'	
		else:
			print 'Q = [',self.Qvec[0],self.Qvec[1],self.Qvec[2],'], hw =',self.En,'meV (fixed kf)'		
		print 'a1 =',int(1e2*A1)/1.e2,'deg, a2 =',int(1e2*A2)/1.e2,'deg'		
		print 'a3 =',int(1e2*A3)/1.e2,'deg, a4 =',int(1e2*A4)/1.e2,'deg'		
		print 'a5 =',int(1e2*A5)/1.e2,'deg, a6 =',int(1e2*A6)/1.e2,'deg\n'	
							
	def updateParams(self,parlist,method):
		# Gets TAS & NRSE parameters from the main panel
		self.method=method
		parlist=[float(i) for i in parlist]
		# Monoch/analyzer lattice constants [Å]
		self.dM=parlist[0]
		self.dA=parlist[1]
		# Mosaicities [arc minutes, FWHM]
		FWHMminstoSigdegs=1/(60*2*sqrt(2*log(2))) 
		if parlist[2] < 1e-5: # Set mochromator mosaic very small if zero
			self.etaM=1e-5
		else:		
			self.etaM=deg2rad(parlist[2]*FWHMminstoSigdegs)		
		if parlist[3] < 1e-5: # Set analyzer mosaic very small if zero
			self.etaA=1e-5
		else:		
			self.etaA=deg2rad(parlist[3]*FWHMminstoSigdegs)			
		if parlist[4] < 1e-5: # Set sample mosaic very small if zero
			self.etaS=1e-5
		else:		
			self.etaS=deg2rad(parlist[4]*FWHMminstoSigdegs)
		# Scattering senses
		self.SM=sign(parlist[6])
		self.SS=sign(parlist[7])
		self.SA=sign(parlist[8])
		# Neutron initial/final wavevector [Å]
		self.kfix=parlist[9]
		self.fix=parlist[10]
		# Horizontal divergences [arc minutes, FWHM]
		if parlist[11] > 1e-5:
			self.alpha0=deg2rad(parlist[11]*FWHMminstoSigdegs)
		else:
			self.alpha0=1e-5
		if parlist[12] > 1e-5:
			self.alpha1=deg2rad(parlist[12]*FWHMminstoSigdegs)
		else:
			self.alpha1=1e-5
		if parlist[13] > 1e-5:
			self.alpha2=deg2rad(parlist[13]*FWHMminstoSigdegs)
		else:
			self.alpha2=1e-5
		if parlist[14] > 1e-5:
			self.alpha3=deg2rad(parlist[14]*FWHMminstoSigdegs)
		else:
			self.alpha3=1e-5			
		# Vertical divergences [arc minutes, FWHM]		
		if parlist[15] > 1e-5:
			self.beta0=deg2rad(parlist[15]*FWHMminstoSigdegs)
		else:
			self.beta0=1e-5
		if parlist[16] > 1e-5:
			self.beta1=deg2rad(parlist[16]*FWHMminstoSigdegs)
		else:
			self.beta1=1e-5
		if parlist[17] > 1e-5:
			self.beta2=deg2rad(parlist[17]*FWHMminstoSigdegs)
		else:
			self.beta2=1e-5
		if parlist[18] > 1e-5:
			self.beta3=deg2rad(parlist[18]*FWHMminstoSigdegs)
		else:
			self.beta3=1e-5
		# Sample lattice constants [Å & °]
		self.a=parlist[19]
		self.b=parlist[20]
		self.c=parlist[21]
		self.alpha=deg2rad(parlist[22])
		self.beta=deg2rad(parlist[23])
		self.gamma=deg2rad(parlist[24])
		# Orienting vectors
		self.Avec=array([parlist[25],parlist[26],parlist[27]])
		self.Bvec=array([parlist[28],parlist[29],parlist[30]])
		# Excitation parameters
		self.Qvec=array([parlist[31],parlist[32],parlist[33]]) # Excitation wavevector
		self.En=parlist[34]
		self.dEdq=parlist[38] # Local slope of the dispersion [meV.Å]
		# Wavectors & excitation velocity
		self.Gmodvec=array([parlist[35],parlist[36],parlist[37]]) # Direction of the dispersion
		self.qvec=array([parlist[39],parlist[40],parlist[41]]) # Q - G
		self.Gvec=self.Qvec-self.qvec # Zone center (Bragg)
		# Sample distribution of lattice constants
		if parlist[5] < 1e-8: # Set delta d/d very small if zero
			self.deltad	=1e-8
		else:
			self.deltad=parlist[5]
		# Dispersion curvature matrix
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
		# NRSE parameters
		self.fmin=parlist[51] # Minimum RF frequency
		self.fmax=parlist[52] # Maximum RF frequency
		self.L4pi=parlist[53] # Flipper distance (normal)
		self.L8pi=parlist[54] # Flipper distance (bootstrap)
		# Spatial parameters (Popovici)
		if self.method == 'pop':
			# Source shape
			self.sourceshape=parlist[55]
			if self.sourceshape == 0:
				self.sourcediameter=parlist[56] # cm
			else:
				self.sourcey=parlist[56]
				self.sourcez=parlist[57]
			# Is there a neutron guide upstream?	
			self.guide=parlist[58]
			if self.guide == 0:
				self.guidehdiv=deg2rad(parlist[59])*FWHMminstoSigdegs
				self.guidevdiv=deg2rad(parlist[60])*FWHMminstoSigdegs
			# Sample shape
			self.sampleshape=parlist[61]
			if self.sampleshape == 0:
				self.samplediameter=parlist[62] # cm
				self.samplethickness=parlist[63]
			else:
				self.samplethickness=parlist[62]
				self.samplewidth=parlist[63]
				self.sampleheight=parlist[64]
			# Detector shape
			self.detshape=parlist[65]
			if self.detshape == 0:
				self.detdiameter=parlist[66] # cm
			else:
				self.detwidth=parlist[66]
				self.detheight=parlist[67]
			# Monochromator geometry
			self.monothickness=parlist[68] # cm
			self.monowidth=parlist[69]
			self.monoheight=parlist[70]
			# Analyzer geometry			
			self.anathickness=parlist[71] # cm
			self.anawidth=parlist[72]
			self.anaheight=parlist[73]
			# Distances [cm]			
			self.L0=parlist[74] # cm
			self.L1=parlist[75]
			self.L2=parlist[76]
			self.L3=parlist[77]
			# Monoch/analyzer curvatures
			self.rmh=parlist[78]/100. # Conversion form 1/m -> 1/cm
			self.rmv=parlist[79]/100.
			self.rah=parlist[80]/100.
			self.rav=parlist[81]/100.	
		# Zone center gap (for viszualisation of the dispersion curve only)
		self.Egap=parlist[85]
		
	def calcAngles(self):
# 		Calculates TAS (crystal angles, scattering vector basis, etc.)
#		and NRSE parameters (coil tilt angles, field integral ratio, etc.)
	
		C=1e20*hbar**2/(2*m_n*eV_to_J*1e-3)
		
		# self.getReciprocalBasis()
		HKLtoRec=self.getReciprocalBasis()
		Avec=dot(HKLtoRec,self.Avec)
		A=norm(Avec)
		Bvec=dot(HKLtoRec,self.Bvec)
		Cvec=cross(Avec,Bvec)
		Qvec=dot(HKLtoRec,self.Qvec)
		Q=norm(Qvec)
		self.Q=Q
		qvec=dot(HKLtoRec,self.qvec)
		q=norm(qvec)
		Gvec=dot(HKLtoRec,self.Gvec)
		G=norm(Gvec)
		
		ki,kf=self.calcKiKf()
		dM=self.dM
		dA=self.dA
		A1=rad2deg(self.SM*arcsin(pi/(ki*dM)))
		A2=2*A1
		if (ki**2+kf**2-Q**2)/(2*ki*kf) > 1 or (ki**2+kf**2-Q**2)/(2*ki*kf) <-1:
			print '\nError: scattering triangle does not close!!!'
			return
		else:
			A4=rad2deg(self.SS*arccos((ki**2+kf**2-Q**2)/(2*ki*kf)))
		A5=rad2deg(self.SA*arcsin(pi/(kf*dA)))
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
			
		Qvec=array([Q,0,0])		
		Avec=A*array([cos(gamma),sin(gamma),0])	
		kivec=ki*array([cos(phi_S),-sin(phi_S),0])
		kfvec=kf*array([-cos(sigma_S),-sin(sigma_S),0])
		qvec=q*array([cos(phi_G),sin(phi_G),0])		
				
		A3=rad2deg(arccos(dot(kivec,Avec)/(ki*A)))
		if dot(cross(kivec,Avec),array([0,0,1]))<0:
			A3*=-1
		
		dEdq=self.dEdq/(2*C)
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
		
		theta1=arccos(dot(kivec,nivec)/(ki*ni))
		if isnan(theta1) == 1:
			theta1=0
		if dot(cross(kivec,nivec),array([0,0,1]))<0:
			theta1*=-1
		theta2=arccos(dot(kfvec,nfvec)/(kf*nf))
		if isnan(theta2) == 1:
			theta2=0
		if dot(cross(kfvec,nfvec),array([0,0,1]))<0:
			theta2*=-1		

		if norm(qvec) == 0:			
			fratio=(ki**2)/(kf**2)
		else:
			fratio=(ni*ki**2*cos(theta1))/(nf*kf**2*cos(theta2))
		
		fmin=self.fmin
		fmax=self.fmax
		L4pi=self.L4pi
		L8pi=self.L8pi
		
		if fratio > 1: # f1 > f2
			if norm(qvec) == 0:	
				tmin=2.522549e-4*2*pi*2*fmin*1e3*L4pi/(kf**2)
				tmax=2.522549e-4*2*pi*4*fmax*1e3*L8pi/(ki**2)			
			else:
				tmin=2.522549e-4*2*pi*2*fmin*1e3*L4pi/(nf*kf**2*cos(theta2))
				tmax=2.522549e-4*2*pi*4*fmax*1e3*L8pi/(ni*ki**2*cos(theta1))
		else:
			if norm(qvec) == 0:	
				tmin=2.522549e-4*2*pi*2*fmin*1e3*L4pi/(ki**2)
				tmax=2.522549e-4*2*pi*4*fmax*1e3*L8pi/(kf**2)				
			else:
				tmin=2.522549e-4*2*pi*2*fmin*1e3*L4pi/(ni*ki**2*cos(theta1))
				tmax=2.522549e-4*2*pi*4*fmax*1e3*L8pi/(nf*kf**2*cos(theta2))	
		
		return A1,A2,A3,A4,A5,A6,theta1,theta2,ni,nf,fratio,tmin,tmax,phi_S,sigma_S
	
	def calcKiKf(self):	
# 		Calculates ki and kf (in AA^-1) knowing fixed k and energy transfer values
	
		C=1e20*hbar**2/(2*m_n*eV_to_J*1e-3)
		if self.fix==1:
			ki=self.kfix
			kf=sqrt(ki**2-self.En/C)
		else:
			kf=self.kfix
			ki=sqrt(self.En/C+kf**2)
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

	def scalarProd(self,x,y,a,b,c,aa,bb,cc):
	# Scalar product within crystal system
		s=x[0]*y[0]*a**2+x[1]*y[1]*b**2+x[2]*y[2]*c**2
		s+=(x[0]*y[1]+y[0]*x[1])*a*b*cos(cc)
		s+=(x[0]*y[2]+y[0]*x[2])*a*c*cos(bb)
		s+=(x[2]*y[1]+y[2]*x[1])*c*b*cos(aa)
		return s
	
	def computeCNmatrix(self):	
	# Determines Cooper-Nathans resolution matrix for host TAS spectrometer.
		
		A1,A2,A3,A4,A5,A6,theta1,theta2,ni,nf,fratio,tmin,tmax,phi_S,sigma_S=self.calcAngles()
		ki,kf=self.calcKiKf()
	
		AM=zeros((6,8)) # AM is in AA-1
		AM[0,0]=ki/tan(deg2rad(A1))/2
		AM[0,1]=-AM[0,0]
		AM[1,1]=ki
		AM[2,3]=AM[1,1]
		AM[3,4]=kf/tan(deg2rad(A5))/2
		AM[3,5]=-AM[3,4]
		AM[4,4]=kf
		AM[5,6]=AM[4,4]

		BM=zeros((4,6)) # BM is dimensionless
		BM[0,0]=cos(phi_S)
		BM[0,1]=sin(phi_S)
		BM[0,3]=-cos(phi_S-deg2rad(A4))
		BM[0,4]=-sin(phi_S-deg2rad(A4))
		BM[1,0]=-BM[0,1]
		BM[1,1]=BM[0,0]
		BM[1,3]=-BM[0,4]
		BM[1,4]=BM[0,4]
		BM[2,2]=1
		BM[2,5]=-1
		# f_const=151.9266;
		BM[3,0]=hbar/m_n*ki*1e10/151.9266
		BM[3,3]=-hbar/m_n*kf*1e10/151.9266
		
		CM=zeros((4,8)) # CM is dimensionless
		CM[0,0]=0.5
		CM[0,1]=0.5
		CM[1,2]=1/(2*sin(deg2rad(A1)))
		CM[1,3]=-1/(2*sin(deg2rad(A1)))
		CM[2,4]=0.5
		CM[2,5]=0.5
		CM[3,6]=1/(2*sin(deg2rad(A5)))
		CM[3,7]=-1/(2*sin(deg2rad(A5)))

		FM=zeros((4,4)) # FM is dimensionless
		FM[0,0]=1/self.etaM**2        
		FM[1,1]=1/self.etaM**2 # should be etaM_v
		FM[2,2]=1/self.etaA**2 
		FM[3,3]=1/self.etaA**2 # should be etaA_v

		GM=zeros((8,8)) # GM is dimensionless
		GM[0,0]=1/self.alpha0**2     
		GM[1,1]=1/self.alpha1**2
		GM[2,2]=1/self.beta0**2
		GM[3,3]=1/self.beta1**2
		GM[4,4]=1/self.alpha2**2
		GM[5,5]=1/self.alpha3**2
		GM[6,6]=1/self.beta2**2
		GM[7,7]=1/self.beta3**2

		# Calculate resolution matrix
		# MCNINV=BM*(AM*(inv(CM'*FM*CM + GM))*AM')*BM';
		MCNinv=dot(dot(BM,dot(dot(AM,inv(dot(transpose(CM),dot(FM,CM))+GM)),transpose(AM))),transpose(BM))
			
		# Sample mosaic contribution (Werner & Pynn)
		MCNinv[1,1]=MCNinv[1,1]+(self.etaS*self.Q)**2
		MCNinv[3,3]=MCNinv[3,3]+(self.etaS*self.Q)**2
				
		MCN=inv(MCNinv)
		MCNtemp=MCN
		MCN=zeros([4,4])	
		MCN[0,0]=MCNtemp[0,0]
		MCN[1,0]=MCNtemp[1,0]
		MCN[0,1]=MCNtemp[0,1]
		MCN[1,1]=MCNtemp[1,1]	
		MCN[0,2]=MCNtemp[0,3]
		MCN[2,0]=MCNtemp[3,0]
		MCN[2,2]=MCNtemp[3,3]
		MCN[2,1]=MCNtemp[3,1]
		MCN[1,2]=MCNtemp[1,3]
		MCN[0,3]=MCNtemp[0,2]
		MCN[3,0]=MCNtemp[2,0]
		MCN[3,3]=MCNtemp[2,2]
		MCN[3,1]=MCNtemp[2,1]
		MCN[1,3]=MCNtemp[1,2]
		
		# Rotation towards the sample (Q,w) coordinate system
		Uvec=array([self.Qvec[0],self.Qvec[1],self.Qvec[2]])/self.Q
		Avec=self.Avec/norm(self.Avec)	
		proj=self.scalarProd(self.Bvec,Avec,self.astar,self.bstar,self.cstar,self.alphastar,self.betastar,self.gammastar)
		Bvec=self.Bvec-Avec*proj
		Bvec/=sqrt(self.scalarProd(Bvec,Bvec,self.astar,self.bstar,self.cstar,self.alphastar,self.betastar,self.gammastar))
		xq=self.scalarProd(Avec,Uvec,self.astar,self.bstar,self.cstar,self.alphastar,self.betastar,self.gammastar)
		yq=self.scalarProd(Bvec,Uvec,self.astar,self.bstar,self.cstar,self.alphastar,self.betastar,self.gammastar)
		TM=zeros([4,4])
		TM[3,3]=1
		TM[2,2]=1
		TM[0,0]=xq
		TM[0,1]=yq
		TM[1,1]=xq
		TM[1,0]=-yq
		MCN=dot(dot(transpose(TM),MCN),TM)	
		
		# Calculate prefactor
		A1*=pi/180
		A5*=pi/180
		R0=(ki**3/tan(A1))*(kf**3/tan(A5))*(2*pi)**4/(64*pi**2*sin(A1)*sin(A5))*sqrt(det(FM)/det(GM+dot(dot(transpose(CM),FM),CM)))
		R0=R0/(2*pi)**2*sqrt(det(MCN))*kf/ki # Chesser & Axe
		R0=R0/sqrt((1+(self.Q*self.etaS)**2*MCN[3,3])*(1+(self.Q*self.etaS)**2*MCN[1,1])) # Sample mosaic
		
		# Return CN matrix
		print '\nCooper-Nathans matrix:\n----------------------\n'
		print MCN
		print '\nPrefactor = %f\n' % R0
		return R0,MCN