#-*- encoding: ISO-8859-1 -*-
# N. MARTIN 19/12/2019
# (c) Laboratoire Léon Brillouin LLB
# From the 'serescal' MATLAB toolbox by K. Habicht (HZB) 
# References:
# K. Habicht et al., J. Appl. Cryst. 36, 1307-1318 (2003)
# F. Groitl et al., J. Appl. Cryst. 51, 818-830 (2018)

from numpy import arccos,arcsin,argmin,array,cos,cross,deg2rad,dot,exp,inf, \
	linspace,log,minimum,pi,ones,rad2deg,sin,sign,sqrt,tan,transpose,zeros
from numpy.linalg import det,inv,norm
from scipy.constants import hbar,m_n
from scipy.constants import e as eV_to_J
from scipy.integrate import quad
from math import isnan

class SERESCAL:

	def __init__(self,parlist,method):
		# Creates a 'SERESCAL' object containing information & methods for the NRSE-related calculations
		self.updateParams(parlist,method)
		self.printSEparams()
		
	def printSEparams(self):		
		# Prints TAS & NRSE parameters in the command line window
		A1,A2,A3,A4,A5,A6,theta1,theta2,ni,nf,fratio,tmin,tmax,phi_S,sigma_S=self.calcAngles()

		print '\nTAS parameters:\n---------------'
		if self.fix == 1:
			print 'Q = [',self.Qvec[0],self.Qvec[1],self.Qvec[2],'], hw =',self.En,'meV (fixed ki)'	
		else:
			print 'Q = [',self.Qvec[0],self.Qvec[1],self.Qvec[2],'], hw =',self.En,'meV (fixed kf)'			
		print 'a1 =',round(A1,2),'deg, a2 =',round(A2,2),'deg'		
		print 'a3 =',round(A3,2),'deg, a4 =',round(A4,2),'deg'		
		print 'a5 =',round(A5,2),'deg, a6 =',round(A6,2),'deg\n'	
		#
		print '\nNRSE parameters:\n----------------'
		print 'Zeta_i =',round(theta1*180/pi,2),'deg'
		print 'Zeta_f =',round(theta2*180/pi,2),'deg'
		print 'Frequency ratio f_i/f_f =',round(fratio,3)
		print 'Minimum tau =',round(tmin,1),'ps'
		print 'Maximum tau =',round(tmax,1),'ps'
	
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
		# Magnetic stuff
		self.extype=parlist[82] # Magnon/phonon
		self.Myy=parlist[83] # In-plane moment
		self.Mzz=parlist[84] # Out-of-plane moment
			
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
		
	def computeInstrumentalResolutionCN(self,tau): 
# 		Calculates instrumental resolution curve taking Cooper-Nathans approximation
# 		for the TAS transmission function. Assumes flat dispersion and perfect sample.
		
		A1,A2,A3,A4,A5,A6,theta1,theta2,ni,nf,fratio,tmin,tmax,phi_S,sigma_S=self.calcAngles()
		
		Ni=(hbar*1e10/m_n)*ni
		Nf=(hbar*1e10/m_n)*nf		
		
		## Instrumental part (Cooper-Nathans)
		ki,kf=self.calcKiKf()
		thetaM = deg2rad(A1)
		thetaA = deg2rad(A5)
		etaM=self.etaM
		etaA=self.etaA
		alpha0=self.alpha0
		alpha1=self.alpha1
		alpha2=self.alpha2
		alpha3=self.alpha3
		beta0=self.beta0
		beta1=self.beta1
		beta2=self.beta2
		beta3=self.beta3		
		
		# Cooper-Nathans constants		
		a1=tan(thetaM)/etaM/(ki*1e10)
		a2=1/etaM/(ki*1e10)
		a3=1/alpha1/(ki*1e10)
		a4=1/alpha2/(kf*1e10)
		a5=tan(thetaA)/etaA/(kf*1e10)
		a6=-1/etaA/(kf*1e10)
		a7=2*tan(thetaM)/alpha0/(ki*1e10)
		a8=1/alpha0/(ki*1e10)
		a9=2*tan(thetaA)/alpha1/(kf*1e10)
		a10=-1/alpha3/(kf*1e10)
		a11=sqrt(1/(4*(sin(thetaM))**2*etaM**2*ki**2+beta0**2*ki**2)+1/(beta1**2*ki**2))*1e-10
		a12=sqrt(1/(4*(sin(thetaA))**2*etaA**2*kf**2+beta3**2*kf**2)+1/(beta2**2*kf**2))*1e-10
		#
		b0=a1*a2+a7*a8
		b1=a2**2+a3**2+a8**2
		b2=a4**2+a6**2+a10**2
		b3=a5**2+a9**2
		b4=a5*a6+a9*a10
		b5=a1**2+a7**2
		
		## Real part (TAS transmission)		
		M11=-0.5*b5*(tan(theta1))**2-0.5*b1+b0*tan(theta1)
		M12=b5/cos(theta1)*tan(theta1)-b0/cos(theta1)
		J11=0.5*hbar/m_n*1/(cos(theta1))**2
		L12=-0.5*b5/(cos(theta1))**2-0.5*b3/(cos(theta2))**2*Ni**2/Nf**2
		N11=-0.5*b3*(tan(theta2))**2-0.5*b2+b4*tan(theta2)
		N12=b3/cos(theta2)*tan(theta2)*Ni/Nf-b4/cos(theta2)*Ni/Nf
		
		## Imaginary part (Larmor precession)
		H12=hbar/2/m_n*1/cos(theta1)**2-hbar/2/m_n*Ni**2/Nf**2/(cos(theta2))**2+ \
			Ni/(ki*cos(theta1)*1e10)-Ni/(kf*cos(theta2)*1e10)*Ni/Nf
		J12=-hbar/m_n*tan(theta1)/cos(theta1)
		K11=-0.5*hbar/m_n/cos(theta2)**2
		K12=hbar/m_n*tan(theta2)/cos(theta2)*Ni/Nf
		
		U12=-L12-1j*H12*tau*1e-12		
		V11=-M11-1j*J11*tau*1e-12
		V12=-M12-1j*J12*tau*1e-12
		W11=-N11-1j*K11*tau*1e-12
		W12=-N12-1j*K12*tau*1e-12
		L44=a11**2-1j*tau*1e-12*hbar/m_n
		L55=a12**2+1j*tau*1e-12*hbar/m_n
		
		len=tau.size
		Linst=zeros([5,5],dtype=complex)	
		Linst_0=zeros([5,5],dtype=complex)		
		y_tau_instr=zeros(len)
		yinstr=zeros(len)
		
		for k in xrange(len):
			
			Linst[0,0]=2*U12[k]*1e20
			Linst[0,1]=V12[k]*1e20
			Linst[0,2]=W12[k]*1e20
			Linst[1,0]=V12[k]*1e20
			Linst[1,1]=2*V11[k]*1e20
			Linst[2,0]=W12[k]*1e20
			Linst[2,2]=2*W11[k]*1e20
			Linst[3,3]=L44[k]*1e20
			Linst[4,4]=L55[k]*1e20
						
			y_tau_instr[k]=sqrt(abs(1./det(Linst)))
		
			
		# Normalize polarization such that P(tau=0)=1
		tau=0
		U12=-L12-1j*H12*tau*1e-12		
		V11=-M11-1j*J11*tau*1e-12
		V12=-M12-1j*J12*tau*1e-12
		W11=-N11-1j*K11*tau*1e-12
		W12=-N12-1j*K12*tau*1e-12
		L44=a11**2-1j*tau*1e-12*hbar/m_n
		L55=a12**2+1j*tau*1e-12*hbar/m_n
		
		Linst_0[0,0]=2*U12*1e20
		Linst_0[0,1]=V12*1e20
		Linst_0[0,2]=W12*1e20
		Linst_0[1,0]=V12*1e20
		Linst_0[1,1]=2*V11*1e20
		Linst_0[2,0]=W12*1e20
		Linst_0[2,2]=2*W11*1e20
		Linst_0[3,3]=L44*1e20
		Linst_0[4,4]=L55*1e20
				
		yinstr=sqrt(abs(det(Linst_0)))*y_tau_instr
		
		return yinstr
	
	def computeCurvatureResolutionCN(self,tau): 
# 		Calculates resolution curve taking Cooper-Nathans approximation
# 		for the TAS transmission function. Includes dispersion curvature
#		and neglects sample imperfections.
		
		A1,A2,A3,A4,A5,A6,theta1,theta2,ni,nf,fratio,tmin,tmax,phi_S,sigma_S=self.calcAngles()
		
		Ni=(hbar*1e10/m_n)*ni
		Nf=(hbar*1e10/m_n)*nf
				
		HKLtoRec=self.getReciprocalBasis()
		
		# Dispersion curvature matrices
		Avec=dot(HKLtoRec,self.Avec)
		A=norm(Avec)
		Bvec=dot(HKLtoRec,self.Bvec)
		Cvec=cross(Avec,Bvec)
		RectoAB=transpose(array([Avec,Bvec,Cvec]))
		H=dot(dot((RectoAB.T),self.H),RectoAB) # H is now expressed in the AB basis
		Qvec=dot(HKLtoRec,self.Qvec)
		Q=norm(Qvec)
		phi_AQ=arcsin(norm(cross(Avec,Qvec))/(A*Q))	
		if dot(cross(Avec,Qvec),Cvec)<0:
			phi_AQ*=-1
		ABtoQ=zeros([3,3])
		ABtoQ[0,0]=cos(phi_AQ)
		ABtoQ[0,1]=sin(phi_AQ)
		ABtoQ[1,0]=-sin(phi_AQ)
		ABtoQ[1,1]=cos(phi_AQ)
		ABtoQ[2,2]=1
		H=dot(dot((ABtoQ.T),H),ABtoQ) # H is now expressed in the scattering plane basis
		H*=eV_to_J*1e-3/hbar*1e-20 # meV.AA^2 -> m^2.s^-1
		
		Theta=zeros([3,6])
		Theta[0,0]=cos(phi_S)
		Theta[0,1]=-sin(phi_S)		# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[0,3]=-cos(sigma_S)	# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[0,4]=sin(sigma_S)		# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[1,0]=sin(phi_S)		# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[1,1]=cos(phi_S)		# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[1,3]=-sin(sigma_S)	# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[1,4]=-cos(sigma_S)	# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[2,2]=1
		Theta[2,5]=-1	
		
		IC=zeros([6,6])
		IC[0,1]=1/cos(theta1)
		IC[0,2]=-tan(theta1)
		IC[3,0]=-1/(Nf*cos(theta2))
		IC[3,1]=Ni/(Nf*cos(theta2))
		IC[3,3]=-tan(theta2)
		IC[1,2]=1
		IC[2,4]=1
		IC[4,3]=1
		IC[5,5]=1
		
		CUM=dot(dot(dot(dot(transpose(IC),transpose(Theta)),H),Theta),IC)
		
		# Instrumental part (Cooper-Nathans)
		
		ki,kf=self.calcKiKf()
		thetaM = deg2rad(A1)
		thetaA = deg2rad(A5)
		etaM=self.etaM
		etaA=self.etaA
		alpha0=self.alpha0		
		alpha1=self.alpha1
		alpha2=self.alpha2
		alpha3=self.alpha3
		beta0=self.beta0
		beta1=self.beta1
		beta2=self.beta2
		beta3=self.beta3		
		
		# Cooper-Nathans constants
		
		a1=tan(thetaM)/etaM/(ki*1e10)
		a2=1/etaM/(ki*1e10)
		a3=1/alpha1/(ki*1e10)
		a4=1/alpha2/(kf*1e10)
		a5=tan(thetaA)/etaA/(kf*1e10)
		a6=-1/etaA/(kf*1e10)
		a7=2*tan(thetaM)/alpha0/(ki*1e10)
		a8=1/alpha0/(ki*1e10)
		a9=2*tan(thetaA)/alpha1/(kf*1e10)
		a10=-1/alpha3/(kf*1e10)
		a11=sqrt(1/(4*(sin(thetaM))**2*etaM**2*ki**2+beta0**2*ki**2)+1/(beta1**2*ki**2))*1e-10
		a12=sqrt(1/(4*(sin(thetaA))**2*etaA**2*kf**2+beta3**2*kf**2)+1/(beta2**2*kf**2))*1e-10
		
		b0=a1*a2+a7*a8
		b1=a2**2+a3**2+a8**2
		b2=a4**2+a6**2+a10**2
		b3=a5**2+a9**2
		b4=a5*a6+a9*a10
		b5=a1**2+a7**2
		
		# Real part (TAS transmission)		
		
		M11=-0.5*b5*(tan(theta1))**2-0.5*b1+b0*tan(theta1)
		M12=b5/cos(theta1)*tan(theta1)-b0/cos(theta1)
		J11=0.5*hbar/m_n*1/(cos(theta1))**2
		L11=-0.5*b3/(cos(theta2))**2*1/Nf**2
		L12=-0.5*b5/(cos(theta1))**2-0.5*b3/(cos(theta2))**2*Ni**2/Nf**2
		L13=b3/(cos(theta2))**2*Ni/Nf**2
		N11=-0.5*b3*(tan(theta2))**2-0.5*b2+b4*tan(theta2)
		N12=b3/cos(theta2)*tan(theta2)*Ni/Nf-b4/cos(theta2)*Ni/Nf
		N13=-b3/cos(theta2)*tan(theta2)/Nf+b4/cos(theta2)/Nf
		
		# Imaginary part (Larmor precession)
		
		H11=-hbar/2/m_n*1/Nf**2-hbar/2/m_n*(tan(theta2))**2/Nf**2-1/(kf*cos(theta2)*1e10)/Nf
		H12=hbar/2/m_n*1/cos(theta1)**2-hbar/2/m_n*Ni**2/Nf**2/(cos(theta2))**2+ \
			Ni/(ki*cos(theta1)*1e10)-Ni/(kf*cos(theta2)*1e10)*Ni/Nf
		H13=hbar/m_n/Nf**2*Ni/cos(theta2)**2+Ni/(kf*cos(theta2)*1e10)*2/Nf
		J12=-hbar/m_n*tan(theta1)/cos(theta1)
		K11=-hbar/2/m_n/cos(theta2)**2
		K12=hbar/m_n*tan(theta2)/cos(theta2)*Ni/Nf
		K13=-hbar/m_n*tan(theta2)/cos(theta2)/Nf
		
		U11=-L11-1j*H11*tau*1e-12
		U12=-L12-1j*H12*tau*1e-12
		U13=-L13-1j*H13*tau*1e-12		
		V11=-M11-1j*J11*tau*1e-12
		V12=-M12-1j*J12*tau*1e-12
		W11=-N11-1j*K11*tau*1e-12
		W12=-N12-1j*K12*tau*1e-12
		W13=-N13-1j*K13*tau*1e-12
		L44=a11**2-1j*tau*1e-12*hbar/m_n
		L55=a12**2+1j*tau*1e-12*hbar/m_n
		
		len=tau.size
		Linst=zeros([6,6],dtype=complex)	
		Linst_0=zeros([6,6],dtype=complex)		
		Lc=zeros([5,5],dtype=complex)
		y_tau_curv=zeros(len)
		ycurv=zeros(len)
		
		for k in xrange(len):
			
			Linst[0,0]=2*U11[k]*1e24
			Linst[0,1]=U13[k]*1e22
			Linst[0,3]=W13[k]*1e22
			Linst[1,0]=U13[k]*1e22
			Linst[1,1]=2*U12[k]*1e20
			Linst[1,2]=V12[k]*1e20
			Linst[1,3]=W12[k]*1e20
			Linst[2,1]=V12[k]*1e20
			Linst[2,2]=2*V11[k]*1e20
			Linst[3,0]=W13[k]*1e22
			Linst[3,1]=W12[k]*1e20
			Linst[3,3]=2*W11[k]*1e20
			Linst[4,4]=L44[k]*1e20
			Linst[5,5]=L55[k]*1e20
			
			Lc=Linst[1:,1:]+1j*tau[k]*1e8*CUM[1:,1:]
			
			y_tau_curv[k]=sqrt(abs(1./det(Lc)))
		
			
		# Normalize polarization such that P(tau=0)=1
		tau=0
		U11=-L11-1j*H11*tau*1e-12
		U12=-L12-1j*H12*tau*1e-12
		U13=-L13-1j*H13*tau*1e-12		
		V11=-M11-1j*J11*tau*1e-12
		V12=-M12-1j*J12*tau*1e-12
		W11=-N11-1j*K11*tau*1e-12
		W12=-N12-1j*K12*tau*1e-12
		W13=-N13-1j*K13*tau*1e-12
		L44=a11**2-1j*tau*1e-12*hbar/m_n
		L55=a12**2+1j*tau*1e-12*hbar/m_n
		
		Linst_0[0,0]=2*U11*1e24
		Linst_0[0,1]=U13*1e22
		Linst_0[0,3]=W13*1e22
		Linst_0[1,0]=U13*1e22
		Linst_0[1,1]=2*U12*1e20
		Linst_0[1,2]=V12*1e20
		Linst_0[1,3]=W12*1e20
		Linst_0[2,1]=V12*1e20
		Linst_0[2,2]=2*V11*1e20
		Linst_0[3,0]=W13*1e22
		Linst_0[3,1]=W12*1e20
		Linst_0[3,3]=2*W11*1e20
		Linst_0[4,4]=L44*1e20
		Linst_0[5,5]=L55*1e20
		
		ycurv=sqrt(abs(det(Linst_0[1:,1:])))*y_tau_curv
		
		return ycurv
		
	def computeSampleImperfResolutionCN(self,tau): 
# 		Calculates resolution curve taking Cooper-Nathans approximation
# 		for the TAS transmission function. Neglects dispersion curvature
#		and account for sample imperfections (mosaic and lattice spacing spreads)
		
		A1,A2,A3,A4,A5,A6,theta1,theta2,ni,nf,fratio,tmin,tmax,phi_S,sigma_S=self.calcAngles()
		
		Ni=(hbar*1e10/m_n)*ni
		Nf=(hbar*1e10/m_n)*nf
		
		HKLtoRec=self.getReciprocalBasis()		
				
		# Sample imperfections effect matrices
		
		Avec=dot(HKLtoRec,self.Avec)
		A=norm(Avec)
		Bvec=dot(HKLtoRec,self.Bvec)
		Cvec=cross(Avec,Bvec)
		dEdq=self.dEdq*eV_to_J/hbar*1e-13;
		Gmodvec=dot(HKLtoRec,self.Gmodvec)
		Cx=dot(Gmodvec/norm(Gmodvec),Avec/norm(Avec))*dEdq
		Cy=dot(Gmodvec/norm(Gmodvec),Bvec/norm(Bvec))*dEdq
		Cz=dot(Gmodvec/norm(Gmodvec),Cvec/norm(Cvec))*dEdq
		if norm(Cz)>1e-3:
			print '\nSlope is defined out of the scattering plane!!!'
		Gvec=dot(HKLtoRec,self.Gvec)
		G=norm(Gvec)
		
		I=zeros([6,9])
		I[0,0]=1
		I[1,1]=1
		I[2,2]=1
		I[3,3]=1
		I[4,4]=1
		I[5,5]=1
		I[0,6]=-Cx*1e-2
		I[0,7]=-Cy*G*1e-22
		I[0,8]=-Cz*G*1e-22
		
		N=zeros([9,9])
		N[6,6]=1/(self.deltad*G)**2
		N[7,7]=1/self.etaS**2
		N[8,8]=1/self.etaS**2 # Vertical and horizontal mosaic spreads are assumed to be equal
		
		W=zeros([9,9])
		W[6,7]=Cy
		W[7,6]=Cy
		W[6,8]=Cz
		W[8,6]=Cz
		W[7,7]=-Cx*G # Correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (55)
		W[8,8]=-Cx*G # Correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (55)
		W*=1e10
		
		# Instrumental part (Cooper-Nathans)
		
		ki,kf=self.calcKiKf()
		thetaM = deg2rad(A1)
		thetaA = deg2rad(A5)
		etaM=self.etaM
		etaA=self.etaA
		alpha0=self.alpha0
		alpha1=self.alpha1
		alpha2=self.alpha2
		alpha3=self.alpha3
		beta0=self.beta0
		beta1=self.beta1
		beta2=self.beta2
		beta3=self.beta3		
		
		# Cooper-Nathans constants
		
		a1=tan(thetaM)/etaM/(ki*1e10)
		a2=1/etaM/(ki*1e10)
		a3=1/alpha1/(ki*1e10)
		a4=1/alpha2/(kf*1e10)
		a5=tan(thetaA)/etaA/(kf*1e10)
		a6=-1/etaA/(kf*1e10)
		a7=2*tan(thetaM)/alpha0/(ki*1e10)
		a8=1/alpha0/(ki*1e10)
		a9=2*tan(thetaA)/alpha1/(kf*1e10)
		a10=-1/alpha3/(kf*1e10)
		a11=sqrt(1/(4*(sin(thetaM))**2*etaM**2*ki**2+beta0**2*ki**2)+1/(beta1**2*ki**2))*1e-10
		a12=sqrt(1/(4*(sin(thetaA))**2*etaA**2*kf**2+beta3**2*kf**2)+1/(beta2**2*kf**2))*1e-10
		
		b0=a1*a2+a7*a8
		b1=a2**2+a3**2+a8**2
		b2=a4**2+a6**2+a10**2
		b3=a5**2+a9**2
		b4=a5*a6+a9*a10
		b5=a1**2+a7**2
		
		# Real part (TAS transmission)		
		
		M11=-0.5*b5*(tan(theta1))**2-0.5*b1+b0*tan(theta1)
		M12=b5/cos(theta1)*tan(theta1)-b0/cos(theta1)
		J11=0.5*hbar/m_n*1/(cos(theta1))**2
		L11=-0.5*b3/(cos(theta2))**2*1/Nf**2
		L12=-0.5*b5/(cos(theta1))**2-0.5*b3/(cos(theta2))**2*Ni**2/Nf**2
		L13=b3/(cos(theta2))**2*Ni/Nf**2
		N11=-0.5*b3*(tan(theta2))**2-0.5*b2+b4*tan(theta2)
		N12=b3/cos(theta2)*tan(theta2)*Ni/Nf-b4/cos(theta2)*Ni/Nf
		N13=-b3/cos(theta2)*tan(theta2)/Nf+b4/cos(theta2)/Nf
		
		# Imaginary part (Larmor precession)
		
		H11=-hbar/2/m_n*1/Nf**2-hbar/2/m_n*(tan(theta2))**2/Nf**2-1/(kf*cos(theta2)*1e10)/Nf
		H12=hbar/2/m_n*1/cos(theta1)**2-hbar/2/m_n*Ni**2/Nf**2/(cos(theta2))**2+ \
			Ni/(ki*cos(theta1)*1e10)-Ni/(kf*cos(theta2)*1e10)*Ni/Nf
		H13=hbar/m_n/Nf**2*Ni/cos(theta2)**2+Ni/(kf*cos(theta2)*1e10)*2/Nf
		J12=-hbar/m_n*tan(theta1)/cos(theta1)
		K11=-hbar/2/m_n/cos(theta2)**2
		K12=hbar/m_n*tan(theta2)/cos(theta2)*Ni/Nf
		K13=-hbar/m_n*tan(theta2)/cos(theta2)/Nf
		
		U11=-L11-1j*H11*tau*1e-12
		U12=-L12-1j*H12*tau*1e-12
		U13=-L13-1j*H13*tau*1e-12		
		V11=-M11-1j*J11*tau*1e-12
		V12=-M12-1j*J12*tau*1e-12
		W11=-N11-1j*K11*tau*1e-12
		W12=-N12-1j*K12*tau*1e-12
		W13=-N13-1j*K13*tau*1e-12
		L44=a11**2-1j*tau*1e-12*hbar/m_n
		L55=a12**2+1j*tau*1e-12*hbar/m_n
		
		len=tau.size
		Linst=zeros([6,6],dtype=complex)	
		Linst_0=zeros([6,6],dtype=complex)		
		Limperf=zeros([9,9],dtype=complex)
		Limperf_0=zeros([6,6],dtype=complex)
		T=array([0,0,0,0,0,Cx,Cy*G,Cz*G])
		y_tau_imperf=zeros(len)
		yimperf=zeros(len)
		
		for k in xrange(len):
			
			Linst[0,0]=2*U11[k]*1e24
			Linst[0,1]=U13[k]*1e22
			Linst[0,3]=W13[k]*1e22
			Linst[1,0]=U13[k]*1e22
			Linst[1,1]=2*U12[k]*1e20
			Linst[1,2]=V12[k]*1e20
			Linst[1,3]=W12[k]*1e20
			Linst[2,1]=V12[k]*1e20
			Linst[2,2]=2*V11[k]*1e20
			Linst[3,0]=W13[k]*1e22
			Linst[3,1]=W12[k]*1e20
			Linst[3,3]=2*W11[k]*1e20
			Linst[4,4]=L44[k]*1e20
			Linst[5,5]=L55[k]*1e20
			
			Limperf=dot(I.T,dot(Linst,I))+N-1j*tau[k]*1e-12*W
			
			linterm=abs(exp(-0.5*(tau[k]**2)*1e-4*dot(T,dot(inv(Limperf[1:,1:]),T.T))))
			
			y_tau_imperf[k]=sqrt(abs(1./det(Limperf[1:,1:])))*linterm
		
			
		# Normalize polarization such that P(tau=0)=1
		tau=0
		U11=-L11-1j*H11*tau*1e-12
		U12=-L12-1j*H12*tau*1e-12
		U13=-L13-1j*H13*tau*1e-12		
		V11=-M11-1j*J11*tau*1e-12
		V12=-M12-1j*J12*tau*1e-12
		W11=-N11-1j*K11*tau*1e-12
		W12=-N12-1j*K12*tau*1e-12
		W13=-N13-1j*K13*tau*1e-12
		L44=a11**2-1j*tau*1e-12*hbar/m_n
		L55=a12**2+1j*tau*1e-12*hbar/m_n
		
		Linst_0[0,0]=2*U11*1e24
		Linst_0[0,1]=U13*1e22
		Linst_0[0,3]=W13*1e22
		Linst_0[1,0]=U13*1e22
		Linst_0[1,1]=2*U12*1e20
		Linst_0[1,2]=V12*1e20
		Linst_0[1,3]=W12*1e20
		Linst_0[2,1]=V12*1e20
		Linst_0[2,2]=2*V11*1e20
		Linst_0[3,0]=W13*1e22
		Linst_0[3,1]=W12*1e20
		Linst_0[3,3]=2*W11*1e20
		Linst_0[4,4]=L44*1e20
		Linst_0[5,5]=L55*1e20
		
		Limperf_0=dot(I.T,dot(Linst_0,I))+N		
		yimperf=sqrt(abs(det(Limperf_0[1:,1:])))*y_tau_imperf
		
		return yimperf
	
	def computeTotalResolutionCN(self,tau): 
# 		Calculates resolution curve taking Cooper-Nathans approximation
# 		for the TAS transmission function
		
		A1,A2,A3,A4,A5,A6,theta1,theta2,ni,nf,fratio,tmin,tmax,phi_S,sigma_S=self.calcAngles()
		
		Ni=(hbar*1e10/m_n)*ni
		Nf=(hbar*1e10/m_n)*nf
		
		HKLtoRec=self.getReciprocalBasis()
		
		# Dispersion curvature matrices
		
		Avec=dot(HKLtoRec,self.Avec)
		A=norm(Avec)
		Bvec=dot(HKLtoRec,self.Bvec)
		Cvec=cross(Avec,Bvec)
		RectoAB=array([Avec,Bvec,Cvec]).T
		H=dot(dot((RectoAB.T),self.H),RectoAB) # H is now expressed in the AB basis
		Qvec=dot(HKLtoRec,self.Qvec)
		Q=norm(Qvec)
		phi_AQ=arcsin(norm(cross(Avec,Qvec))/(A*Q))	
		if dot(cross(Avec,Qvec),Cvec)<0:
			phi_AQ*=-1
		ABtoQ=zeros([3,3])
		ABtoQ[0,0]=cos(phi_AQ)
		ABtoQ[0,1]=sin(phi_AQ)
		ABtoQ[1,0]=-sin(phi_AQ)
		ABtoQ[1,1]=cos(phi_AQ)
		ABtoQ[2,2]=1
		H=dot(dot((ABtoQ.T),H),ABtoQ) # H is now expressed in the scattering plane basis
		H*=eV_to_J*1e-3/hbar*1e-20 # meV.AA^2 -> m^2.s^-1
		
		Theta=zeros([3,6])
		Theta[0,0]=cos(phi_S)
		Theta[0,1]=-sin(phi_S)		# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[0,3]=-cos(sigma_S)	# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[0,4]=sin(sigma_S)		# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[1,0]=sin(phi_S)		# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[1,1]=cos(phi_S)		# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[1,3]=-sin(sigma_S)	# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[1,4]=-cos(sigma_S)	# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[2,2]=1
		Theta[2,5]=-1	

		IC=zeros([6,6])
		IC[0,1]=1/cos(theta1)
		IC[0,2]=-tan(theta1)
		IC[3,0]=-1/(Nf*cos(theta2))
		IC[3,1]=Ni/(Nf*cos(theta2))
		IC[3,3]=-tan(theta2)
		IC[1,2]=1
		IC[2,4]=1
		IC[4,3]=1
		IC[5,5]=1
		
		CUM=dot(dot(dot(dot(transpose(IC),transpose(Theta)),H),Theta),IC)
		CUM[1:,1:]*=1e20
		CUM[0,0]*=1e24
		CUM[0,1:]*=1e22
		CUM[1:,0]*=1e22
		
		# Sample imperfections effect matrices
		
		dEdq=self.dEdq*eV_to_J/hbar*1e-13;
		Gmodvec=dot(HKLtoRec,self.Gmodvec)
		Cx=dot(Gmodvec/norm(Gmodvec),Avec/norm(Avec))*dEdq
		Cy=dot(Gmodvec/norm(Gmodvec),Bvec/norm(Bvec))*dEdq
		Cz=dot(Gmodvec/norm(Gmodvec),Cvec/norm(Cvec))*dEdq
		if norm(Cz)>1e-3:
			print '\nSlope is defined out of the scattering plane!!!'
		Gvec=dot(HKLtoRec,self.Gvec)
		G=norm(Gvec)
		
		I=zeros([6,9])
		I[0,0]=1
		I[1,1]=1
		I[2,2]=1
		I[3,3]=1
		I[4,4]=1
		I[5,5]=1
		I[0,6]=-Cx*1e-2
		I[0,7]=-Cy*G*1e-22
		I[0,8]=-Cz*G*1e-22
		
		N=zeros([9,9])
		N[6,6]=1/(self.deltad*G)**2
		N[7,7]=1/self.etaS**2
		N[8,8]=1/self.etaS**2 # Vertical and horizontal mosaic spreads are assumed to be equal
		
		W=zeros([9,9])
		W[6,7]=Cy
		W[7,6]=Cy
		W[6,8]=Cz
		W[8,6]=Cz
		W[7,7]=-Cx*G # Correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (55)
		W[8,8]=-Cx*G # Correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (55)
		W*=1e10
		
		# Instrumental part (Cooper-Nathans)
		
		ki,kf=self.calcKiKf()
		thetaM = deg2rad(A1)
		thetaA = deg2rad(A5)
		etaM=self.etaM
		etaA=self.etaA
		alpha0=self.alpha0
		alpha1=self.alpha1
		alpha2=self.alpha2
		alpha3=self.alpha3
		beta0=self.beta0
		beta1=self.beta1
		beta2=self.beta2
		beta3=self.beta3		
		
		# Cooper-Nathans constants
		
		a1=tan(thetaM)/etaM/(ki*1e10)
		a2=1/etaM/(ki*1e10)
		a3=1/alpha1/(ki*1e10)
		a4=1/alpha2/(kf*1e10)
		a5=tan(thetaA)/etaA/(kf*1e10)
		a6=-1/etaA/(kf*1e10)
		a7=2*tan(thetaM)/alpha0/(ki*1e10)
		a8=1/alpha0/(ki*1e10)
		a9=2*tan(thetaA)/alpha1/(kf*1e10)
		a10=-1/alpha3/(kf*1e10)
		a11=sqrt(1/(4*(sin(thetaM))**2*etaM**2*ki**2+beta0**2*ki**2)+1/(beta1**2*ki**2))*1e-10
		a12=sqrt(1/(4*(sin(thetaA))**2*etaA**2*kf**2+beta3**2*kf**2)+1/(beta2**2*kf**2))*1e-10
		
		b0=a1*a2+a7*a8
		b1=a2**2+a3**2+a8**2
		b2=a4**2+a6**2+a10**2
		b3=a5**2+a9**2
		b4=a5*a6+a9*a10
		b5=a1**2+a7**2
		
		# Real part (TAS transmission)		
		
		M11=-0.5*b5*(tan(theta1))**2-0.5*b1+b0*tan(theta1)
		M12=b5/cos(theta1)*tan(theta1)-b0/cos(theta1)
		J11=0.5*hbar/m_n*1/(cos(theta1))**2
		L11=-0.5*b3/(cos(theta2))**2*1/Nf**2
		L12=-0.5*b5/(cos(theta1))**2-0.5*b3/(cos(theta2))**2*Ni**2/Nf**2
		L13=b3/(cos(theta2))**2*Ni/Nf**2
		N11=-0.5*b3*(tan(theta2))**2-0.5*b2+b4*tan(theta2)
		N12=b3/cos(theta2)*tan(theta2)*Ni/Nf-b4/cos(theta2)*Ni/Nf
		N13=-b3/cos(theta2)*tan(theta2)/Nf+b4/cos(theta2)/Nf
		
		# Imaginary part (Larmor precession)
		
		H11=-hbar/2/m_n*1/Nf**2-hbar/2/m_n*(tan(theta2))**2/Nf**2-1/(kf*cos(theta2)*1e10)/Nf
		H12=hbar/2/m_n*1/cos(theta1)**2-hbar/2/m_n*Ni**2/Nf**2/(cos(theta2))**2+ \
			Ni/(ki*cos(theta1)*1e10)-Ni/(kf*cos(theta2)*1e10)*Ni/Nf
		H13=hbar/m_n/Nf**2*Ni/cos(theta2)**2+Ni/(kf*cos(theta2)*1e10)*2/Nf
		J12=-hbar/m_n*tan(theta1)/cos(theta1)
		K11=-hbar/2/m_n/cos(theta2)**2
		K12=hbar/m_n*tan(theta2)/cos(theta2)*Ni/Nf
		K13=-hbar/m_n*tan(theta2)/cos(theta2)/Nf
		
		U11=-L11-1j*H11*tau*1e-12
		U12=-L12-1j*H12*tau*1e-12
		U13=-L13-1j*H13*tau*1e-12		
		V11=-M11-1j*J11*tau*1e-12
		V12=-M12-1j*J12*tau*1e-12
		W11=-N11-1j*K11*tau*1e-12
		W12=-N12-1j*K12*tau*1e-12
		W13=-N13-1j*K13*tau*1e-12
		L44=a11**2-1j*tau*1e-12*hbar/m_n
		L55=a12**2+1j*tau*1e-12*hbar/m_n
		
		len=tau.size
		Linst=zeros([6,6],dtype=complex)	
		Linst_0=zeros([6,6],dtype=complex)		
		Lc=zeros([6,6],dtype=complex)
		Ltot=zeros([9,9],dtype=complex)
		Ltot_0=zeros([6,6],dtype=complex)
		T=array([0,0,0,0,0,Cx,Cy*G,Cz*G])
		y_tau_tot=zeros(len)
		ytot=zeros(len)
		
		for k in xrange(len):
			
			Linst[0,0]=2*U11[k]*1e24
			Linst[0,1]=U13[k]*1e22
			Linst[0,3]=W13[k]*1e22
			Linst[1,0]=U13[k]*1e22
			Linst[1,1]=2*U12[k]*1e20
			Linst[1,2]=V12[k]*1e20
			Linst[1,3]=W12[k]*1e20
			Linst[2,1]=V12[k]*1e20
			Linst[2,2]=2*V11[k]*1e20
			Linst[3,0]=W13[k]*1e22
			Linst[3,1]=W12[k]*1e20
			Linst[3,3]=2*W11[k]*1e20
			Linst[4,4]=L44[k]*1e20
			Linst[5,5]=L55[k]*1e20
			
			Lc=Linst+1j*tau[k]*1e-12*CUM
			
			Ltot=dot(I.T,dot(Lc,I))+N-1j*tau[k]*1e-12*W
			
			linterm=abs(exp(-0.5*(tau[k]**2)*1e-4*dot(T,dot(inv(Ltot[1:,1:]),T.T))))
			
			y_tau_tot[k]=sqrt(abs(1./det(Ltot[1:,1:])))*linterm
		
			
		# Normalize polarization such that P(tau=0)=1
		tau=0
		U11=-L11-1j*H11*tau*1e-12
		U12=-L12-1j*H12*tau*1e-12
		U13=-L13-1j*H13*tau*1e-12		
		V11=-M11-1j*J11*tau*1e-12
		V12=-M12-1j*J12*tau*1e-12
		W11=-N11-1j*K11*tau*1e-12
		W12=-N12-1j*K12*tau*1e-12
		W13=-N13-1j*K13*tau*1e-12
		L44=a11**2-1j*tau*1e-12*hbar/m_n
		L55=a12**2+1j*tau*1e-12*hbar/m_n
		
		Linst_0[0,0]=2*U11*1e24
		Linst_0[0,1]=U13*1e22
		Linst_0[0,3]=W13*1e22
		Linst_0[1,0]=U13*1e22
		Linst_0[1,1]=2*U12*1e20
		Linst_0[1,2]=V12*1e20
		Linst_0[1,3]=W12*1e20
		Linst_0[2,1]=V12*1e20
		Linst_0[2,2]=2*V11*1e20
		Linst_0[3,0]=W13*1e22
		Linst_0[3,1]=W12*1e20
		Linst_0[3,3]=2*W11*1e20
		Linst_0[4,4]=L44*1e20
		Linst_0[5,5]=L55*1e20
		
		Ltot_0=dot(I.T,dot(Linst_0,I))+N		
		ytot=sqrt(abs(det(Ltot_0[1:,1:])))*y_tau_tot
		
		return ytot
	
	def computeInstrumentalResolutionPOP(self,tau): 
# 		Calculates resolution curve taking Popovici approximation
# 		for the TAS transmission function
		
		A1,A2,A3,A4,A5,A6,theta1,theta2,ni,nf,fratio,tmin,tmax,phi_S,sigma_S=self.calcAngles()
		
		Ni=(hbar*1e10/m_n)*ni
		Nf=(hbar*1e10/m_n)*nf
		
		HKLtoRec=self.getReciprocalBasis()	

		IC=zeros([6,6])
		IC[0,1]=1/cos(theta1)
		IC[0,2]=-tan(theta1)
		IC[3,0]=-1/(Nf*cos(theta2))
		IC[3,1]=Ni/(Nf*cos(theta2))
		IC[3,3]=-tan(theta2)
		IC[1,2]=1
		IC[2,4]=1
		IC[4,3]=1
		IC[5,5]=1
		
		# Instrumental part (Popovici)
		
		ki,kf=self.calcKiKf()
		thetaM=deg2rad(A1)
		SM=self.SM
		thetaS=deg2rad(0.5*A4)
		SS=self.SS
		thetaA=deg2rad(A5)
		SA=self.SA
		
		etaM=self.etaM
		etaMv=self.etaM
		etaS=self.etaS
		etaA=self.etaA
		etaAv=self.etaA
		
		if self.guide == 0:
			alpha0=self.guidehdiv*2*pi/ki
			beta0=self.guidevdiv*2*pi/ki
		else:
			alpha0=self.alpha0
			beta0=self.beta0
		
		alpha1=self.alpha1
		alpha2=self.alpha2
		alpha3=self.alpha3

		beta1=self.beta1
		beta2=self.beta2
		beta3=self.beta3
		
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
		
		L0=self.L0
		L1=self.L1
		L2=self.L2
		L3=self.L3
		
		rmh=self.rmh*SM
		rmv=self.rmv
		rah=self.rah*SA
		rav=self.rav
		
		AM=zeros([6,8]) # expressed in AA^-1
		AM[0,0]=0.5*ki/tan(thetaM)
		AM[0,1]=-0.5*ki/tan(thetaM)
		AM[1,1]=ki
		AM[2,3]=ki
		AM[3,4]=0.5*kf/tan(thetaA)
		AM[3,5]=-0.5*kf/tan(thetaA)
		AM[4,4]=kf
		AM[5,6]=kf
		
		FM=zeros([4,4]) # dimensionless
		FM[0,0]=1/etaM**2
		FM[1,1]=1/etaMv**2
		FM[2,2]=1/etaA**2
		FM[3,3]=1/etaAv**2
		
		GM=zeros([8,8]) # dimensionless
		GM[0,0]=1/alpha0**2
		GM[1,1]=1/alpha1**2
		GM[2,2]=1/beta0**2
		GM[3,3]=1/beta1**2		
		GM[4,4]=1/alpha2**2
		GM[5,5]=1/alpha3**2
		GM[6,6]=1/beta2**2
		GM[7,7]=1/beta3**2		
		
		CM=zeros([4,8]) # dimensionless
		CM[0,0]=0.5
		CM[0,1]=0.5
		CM[2,4]=0.5
		CM[2,5]=0.5
		CM[1,2]=0.5/sin(thetaM)
		CM[1,3]=-0.5/sin(thetaM)
		CM[3,6]=0.5/sin(thetaA)
		CM[3,7]=-0.5/sin(thetaA)
		
		TM=zeros([4,13]) # expressed in cm^-1
		TM[0,0]=-0.5/L0
		TM[0,2]=0.5*cos(thetaM)*(1/L1-1/L0)
		TM[0,3]=0.5*sin(thetaM)*(1/L0+1/L1-2*rmh/abs(sin(thetaM)))
		TM[0,5]=0.5*sin(thetaS)/L1
		TM[0,6]=0.5*cos(thetaS)/L1
		TM[1,1]=-0.5/(L0*sin(thetaM))
		TM[1,4]=0.5*(1/L0+1/L1-2*rmv*sin(thetaM))/sin(abs(thetaM))
		TM[1,7]=-0.5/(L1*sin(thetaM))
		TM[2,5]=0.5*sin(thetaS)/L2
		TM[2,6]=-0.5*cos(thetaS)/L2
		TM[2,8]=0.5*cos(thetaA)*(1/L3-1/L2)
		TM[2,9]=0.5*sin(thetaA)*(1/L2+1/L3-2*rah/abs(sin(thetaA)))
		TM[2,11]=0.5/L3
		TM[3,7]=-0.5/(L2*sin(thetaA))
		TM[3,10]=0.5*(1/L2+1/L3-2*rav*sin(thetaA))/abs(sin(thetaA))
		TM[3,12]=-0.5/(L3*sin(thetaA))
		
		DM=zeros([8,13]) # expressed in cm^-1
		DM[0,0]=-1/L0
		DM[0,2]=-cos(thetaM)/L0
		DM[0,3]=sin(thetaM)/L0
		DM[2,1]=-1/L0
		DM[2,4]=1/L0
		DM[1,2]=cos(thetaM)/L1
		DM[1,3]=sin(thetaM)/L1
		DM[1,5]=sin(thetaS)/L1
		DM[1,6]=cos(thetaS)/L1
		DM[3,4]=-1/L1
		DM[3,7]=1/L1
		DM[4,5]=sin(thetaS)/L2
		DM[4,6]=-cos(thetaS)/L2
		DM[4,8]=-cos(thetaA)/L2
		DM[4,9]=sin(thetaA)/L2
		DM[6,7]=-1/L2
		DM[6,10]=1/L2
		DM[5,8]=cos(thetaA)/L3
		DM[5,9]=sin(thetaA)/L3
		DM[5,11]=1/L3
		DM[7,10]=-1/L3
		DM[7,12]=1/L3
		
		# Spatial covariances matrix
		SC=zeros([13,13])
		
		# Source covariance matrix
		if self.sourceshape == 0:
			SC[0,0]=1/16.*sourced**2
			SC[1,1]=1/16.*sourced**2
		else:
			SC[0,0]=1/12.*sourcey**2
			SC[1,1]=1/12.*sourcez**2
		
		# Monochromator covariance matrix
		Smono=zeros([3,3])
		Smono[0,0]=1/12.*monoth**2
		Smono[1,1]=1/12.*monow**2
		Smono[2,2]=1/12.*monoh**2
		SC[2:5,2:5]=Smono
		
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
		SC[5:8,5:8]=Ssam
				
		# Analyzer covariance matrix
		Sana=zeros([3,3])
		Sana[0,0]=1/12.*anath**2
		Sana[1,1]=1/12.*anaw**2
		Sana[2,2]=1/12.*anah**2
		SC[8:11,8:11]=Sana

		# Detector covariance matrix
		if self.detshape == 0:
			SC[11,11]=1/16.*detd**2
			SC[12,12]=1/16.*detd**2
		else:
			SC[11,11]=1/12.*dety**2
			SC[12,12]=1/12.*detz**2		
			
		# Builds full Popovici matrix
		IM=inv(IC)
		SM=inv(SC)			
		LMSpatial=inv(inv(dot(DM,dot(inv(SM+dot(transpose(TM),dot(FM,TM))),transpose(DM))))+GM)
		LMSpatial=dot(AM,dot(LMSpatial,transpose(AM)))
		LMSpatial=inv(dot(IM,dot(LMSpatial,transpose(IM))))
		LMSub=LMSpatial[1:,1:]*1e-20
		
		# Imaginary part of the total resolution matrix (i.e. Larmor precession)	
		H11=-hbar/2/m_n*1/Nf**2-hbar/2/m_n*(tan(theta2))**2/Nf**2-1/(kf*cos(theta2)*1e10)/Nf
		H12=hbar/2/m_n*1/cos(theta1)**2-hbar/2/m_n*Ni**2/Nf**2/(cos(theta2))**2+ \
			Ni/(ki*cos(theta1)*1e10)-Ni/(kf*cos(theta2)*1e10)*Ni/Nf
		J11=0.5*hbar/m_n*1/(cos(theta1))**2
		J12=-hbar/m_n*tan(theta1)/cos(theta1)
		K11=-hbar/2/m_n/cos(theta2)**2
		K12=hbar/m_n*tan(theta2)/cos(theta2)*Ni/Nf
		
		Limag=zeros([5,5])
		Limag[0,0]=2*H12
		Limag[0,1]=J12
		Limag[0,2]=K12
		Limag[1,0]=J12
		Limag[1,1]=2*J11
		Limag[2,0]=K12
		Limag[2,2]=2*K11
		Limag[3,3]=hbar/m_n
		Limag[4,4]=-hbar/m_n
				
		len=tau.size
		Linst=zeros([6,6],dtype=complex)		
		y_tau_instr=zeros(len)
		yinstr=zeros(len)
		
		for k in xrange(len):
			
			Linst=LMSub-1j*tau[k]*1e-12*Limag
			y_tau_instr[k]=sqrt(abs(1./det(Linst)))
		
		# Normalize polarization such that P(tau=0)=1		
		yinstr=sqrt(abs(det(LMSub)))*y_tau_instr
			
		return yinstr

	def computeCurvatureResolutionPOP(self,tau): 
# 		Calculates resolution curve taking Popovici approximation
# 		for the TAS transmission function
		
		A1,A2,A3,A4,A5,A6,theta1,theta2,ni,nf,fratio,tmin,tmax,phi_S,sigma_S=self.calcAngles()
		
		Ni=(hbar*1e10/m_n)*ni
		Nf=(hbar*1e10/m_n)*nf
		
		HKLtoRec=self.getReciprocalBasis()
		
		# Dispersion curvature matrices
		
		Avec=dot(HKLtoRec,self.Avec)
		A=norm(Avec)
		Bvec=dot(HKLtoRec,self.Bvec)
		Cvec=cross(Avec,Bvec)
		RectoAB=array([Avec,Bvec,Cvec]).T
		H=dot(dot((RectoAB.T),self.H),RectoAB) # H is now expressed in the AB basis
		Qvec=dot(HKLtoRec,self.Qvec)
		Q=norm(Qvec)
		phi_AQ=arcsin(norm(cross(Avec,Qvec))/(A*Q))	
		if dot(cross(Avec,Qvec),Cvec)<0:
			phi_AQ*=-1
		ABtoQ=zeros([3,3])
		ABtoQ[0,0]=cos(phi_AQ)
		ABtoQ[0,1]=sin(phi_AQ)
		ABtoQ[1,0]=-sin(phi_AQ)
		ABtoQ[1,1]=cos(phi_AQ)
		ABtoQ[2,2]=1
		H=dot(dot((ABtoQ.T),H),ABtoQ) # H is now expressed in the scattering plane basis
		H*=eV_to_J*1e-3/hbar*1e-20 # meV.AA^2 -> m^2.s^-1
		
		Theta=zeros([3,6])
		Theta[0,0]=cos(phi_S)
		Theta[0,1]=-sin(phi_S)		# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[0,3]=-cos(sigma_S)	# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[0,4]=sin(sigma_S)		# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[1,0]=sin(phi_S)		# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[1,1]=cos(phi_S)		# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[1,3]=-sin(sigma_S)	# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[1,4]=-cos(sigma_S)	# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[2,2]=1
		Theta[2,5]=-1	

		IC=zeros([6,6])
		IC[0,1]=1/cos(theta1)
		IC[0,2]=-tan(theta1)
		IC[3,0]=-1/(Nf*cos(theta2))
		IC[3,1]=Ni/(Nf*cos(theta2))
		IC[3,3]=-tan(theta2)
		IC[1,2]=1
		IC[2,4]=1
		IC[4,3]=1
		IC[5,5]=1
		
		CUM=dot(IC.T,dot(Theta.T,dot(H,dot(Theta,IC))));
		
		# Instrumental part (Popovici)	
		ki,kf=self.calcKiKf()
		thetaM=deg2rad(A1)
		SM=self.SM
		thetaS=deg2rad(0.5*A4)
		SS=self.SS
		thetaA=deg2rad(A5)
		SA=self.SA
		
		etaM=self.etaM
		etaMv=self.etaM
		etaS=self.etaS
		etaA=self.etaA
		etaAv=self.etaA
		
		if self.guide == 0:
			alpha0=self.guidehdiv*2*pi/ki
			beta0=self.guidevdiv*2*pi/ki
		else:
			alpha0=self.alpha0
			beta0=self.beta0
			
		alpha1=self.alpha1
		alpha2=self.alpha2
		alpha3=self.alpha3

		beta1=self.beta1
		beta2=self.beta2
		beta3=self.beta3
		
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
		
		L0=self.L0
		L1=self.L1
		L2=self.L2
		L3=self.L3
		
		rmh=self.rmh*SM
		rmv=self.rmv
		rah=self.rah*SA
		rav=self.rav
		
		AM=zeros([6,8]) # expressed in AA^-1
		AM[0,0]=0.5*ki/tan(thetaM)
		AM[0,1]=-0.5*ki/tan(thetaM)
		AM[1,1]=ki
		AM[2,3]=ki
		AM[3,4]=0.5*kf/tan(thetaA)
		AM[3,5]=-0.5*kf/tan(thetaA)
		AM[4,4]=kf
		AM[5,6]=kf
		
		FM=zeros([4,4]) # dimensionless
		FM[0,0]=1/etaM**2
		FM[1,1]=1/etaMv**2
		FM[2,2]=1/etaA**2
		FM[3,3]=1/etaAv**2

		GM=zeros([8,8]) # dimensionless
		GM[0,0]=1/alpha0**2
		GM[1,1]=1/alpha1**2
		GM[2,2]=1/beta0**2
		GM[3,3]=1/beta1**2		
		GM[4,4]=1/alpha2**2
		GM[5,5]=1/alpha3**2
		GM[6,6]=1/beta2**2
		GM[7,7]=1/beta3**2		
		
		CM=zeros([4,8]) # dimensionless
		CM[0,0]=0.5
		CM[0,1]=0.5
		CM[2,4]=0.5
		CM[2,5]=0.5
		CM[1,2]=0.5/sin(thetaM)
		CM[1,3]=-0.5/sin(thetaM)
		CM[3,6]=0.5/sin(thetaA)
		CM[3,7]=-0.5/sin(thetaA)
		
		TM=zeros([4,13]) # expressed in cm^-1
		TM[0,0]=-0.5/L0
		TM[0,2]=0.5*cos(thetaM)*(1/L1-1/L0)
		TM[0,3]=0.5*sin(thetaM)*(1/L0+1/L1-2*rmh/abs(sin(thetaM)))
		TM[0,5]=0.5*sin(thetaS)/L1
		TM[0,6]=0.5*cos(thetaS)/L1
		TM[1,1]=-0.5/(L0*sin(thetaM))
		TM[1,4]=0.5*(1/L0+1/L1-2*rmv*sin(thetaM))/abs(sin(thetaM))
		TM[1,7]=-0.5/(L1*sin(thetaM))
		TM[2,5]=0.5*sin(thetaS)/L2
		TM[2,6]=-0.5*cos(thetaS)/L2
		TM[2,8]=0.5*cos(thetaA)*(1/L3-1/L2)
		TM[2,9]=0.5*sin(thetaA)*(1/L2+1/L3-2*rah/abs(sin(thetaA)))
		TM[2,11]=0.5/L3
		TM[3,7]=-0.5/(L2*sin(thetaA))
		TM[3,10]=0.5*(1/L2+1/L3-2*rav*sin(thetaA))/abs(sin(thetaA))
		TM[3,12]=-0.5/(L3*sin(thetaA))
		
		DM=zeros([8,13]) # expressed in cm^-1
		DM[0,0]=-1/L0
		DM[0,2]=-cos(thetaM)/L0
		DM[0,3]=sin(thetaM)/L0
		DM[2,1]=-1/L0
		DM[2,4]=1/L0
		DM[1,2]=cos(thetaM)/L1
		DM[1,3]=sin(thetaM)/L1
		DM[1,5]=sin(thetaS)/L1
		DM[1,6]=cos(thetaS)/L1
		DM[3,4]=-1/L1
		DM[3,7]=1/L1
		DM[4,5]=sin(thetaS)/L2
		DM[4,6]=-cos(thetaS)/L2
		DM[4,8]=-cos(thetaA)/L2
		DM[4,9]=sin(thetaA)/L2
		DM[6,7]=-1/L2
		DM[6,10]=1/L2
		DM[5,8]=cos(thetaA)/L3
		DM[5,9]=sin(thetaA)/L3
		DM[5,11]=1/L3
		DM[7,10]=-1/L3
		DM[7,12]=1/L3
		
		# Spatial covariances matrix
		SC=zeros([13,13])
		
		# Source covariance matrix
		if self.sourceshape == 0:
			SC[0,0]=1/16.*sourced**2
			SC[1,1]=1/16.*sourced**2
		else:
			SC[0,0]=1/12.*sourcey**2
			SC[1,1]=1/12.*sourcez**2
		
		# Monochromator covariance matrix
		Smono=zeros([3,3])
		Smono[0,0]=1/12.*monoth**2
		Smono[1,1]=1/12.*monow**2
		Smono[2,2]=1/12.*monoh**2
		SC[2:5,2:5]=Smono
		
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
		SC[5:8,5:8]=Ssam
				
		# Analyzer covariance matrix
		Sana=zeros([3,3])
		Sana[0,0]=1/12.*anath**2
		Sana[1,1]=1/12.*anaw**2
		Sana[2,2]=1/12.*anah**2
		SC[8:11,8:11]=Sana

		# Detector covariance matrix
		if self.detshape == 0:
			SC[11,11]=1/16.*detd**2
			SC[12,12]=1/16.*detd**2
		else:
			SC[11,11]=1/12.*dety**2
			SC[12,12]=1/12.*detz**2		
		
		# Builds full Popovici matrix		
		IM=inv(IC)
		SM=inv(SC)			
		LMSpatial=inv(inv(dot(DM,dot(inv(SM+dot(transpose(TM),dot(FM,TM))),transpose(DM))))+GM)
		LMSpatial=dot(AM,dot(LMSpatial,transpose(AM)))
		LMSpatial=inv(dot(IM,dot(LMSpatial,transpose(IM))))
		LMSub=LMSpatial[1:,1:]
		
		# Imaginary part of the total resolution matrix (i.e. Larmor precession)	
		H12=hbar/2/m_n*1/cos(theta1)**2-hbar/2/m_n*Ni**2/Nf**2/(cos(theta2))**2+ \
			Ni/(ki*cos(theta1)*1e10)-Ni/(kf*cos(theta2)*1e10)*Ni/Nf
		J11=0.5*hbar/m_n*1/(cos(theta1))**2
		J12=-hbar/m_n*tan(theta1)/cos(theta1)
		K11=-hbar/2/m_n/cos(theta2)**2
		K12=hbar/m_n*tan(theta2)/cos(theta2)*Ni/Nf
		
		Limag=zeros([5,5])
		Limag[0,0]=2*H12
		Limag[0,1]=J12
		Limag[0,2]=K12
		Limag[1,0]=J12
		Limag[1,1]=2*J11
		Limag[2,0]=K12
		Limag[2,2]=2*K11
		Limag[3,3]=hbar/m_n
		Limag[4,4]=-hbar/m_n
		
		len=tau.size
		Linst=zeros([5,5],dtype=complex)		
		Lc=zeros([5,5],dtype=complex)
		y_tau_curv=zeros(len)
		ycurv=zeros(len)
		
		for k in xrange(len):
			
			Linst=LMSub-1j*tau[k]*1e8*Limag
			Lc=Linst+1j*tau[k]*1e8*CUM[1:,1:]
			y_tau_curv[k]=sqrt(abs(1./det(Lc)))
		
		# Normalize polarization such that P(tau=0)=1		
		ycurv=sqrt(abs(det(LMSub)))*y_tau_curv
		
		return ycurv	
		
	def computeSampleImperfResolutionPOP(self,tau): 
# 		Calculates resolution curve taking Popovici approximation
# 		for the TAS transmission function
		
		A1,A2,A3,A4,A5,A6,theta1,theta2,ni,nf,fratio,tmin,tmax,phi_S,sigma_S=self.calcAngles()
		
		Ni=(hbar*1e10/m_n)*ni
		Nf=(hbar*1e10/m_n)*nf
		
		HKLtoRec=self.getReciprocalBasis()
		
		IC=zeros([6,6])
		IC[0,1]=1/cos(theta1)
		IC[0,2]=-tan(theta1)
		IC[3,0]=-1/(Nf*cos(theta2))
		IC[3,1]=Ni/(Nf*cos(theta2))
		IC[3,3]=-tan(theta2)
		IC[1,2]=1
		IC[2,4]=1
		IC[4,3]=1
		IC[5,5]=1
		
		Avec=dot(HKLtoRec,self.Avec)
		A=norm(Avec)
		Bvec=dot(HKLtoRec,self.Bvec)
		Cvec=cross(Avec,Bvec)
		
		# Sample imperfections effect matrices
		
		dEdq=self.dEdq*eV_to_J/hbar*1e-13;
		Gmodvec=dot(HKLtoRec,self.Gmodvec)
		Cx=dot(Gmodvec/norm(Gmodvec),Avec/norm(Avec))*dEdq
		Cy=dot(Gmodvec/norm(Gmodvec),Bvec/norm(Bvec))*dEdq
		Cz=dot(Gmodvec/norm(Gmodvec),Cvec/norm(Cvec))*dEdq
		if norm(Cz)>1e-3:
			print '\nSlope is defined out of the scattering plane!!!'
		Gvec=dot(HKLtoRec,self.Gvec)
		G=norm(Gvec)
		
		I=zeros([6,9])
		I[0,0]=1
		I[1,1]=1
		I[2,2]=1
		I[3,3]=1
		I[4,4]=1
		I[5,5]=1
		I[0,6]=-Cx*1e-2
		I[0,7]=-Cy*G*1e-22
		I[0,8]=-Cz*G*1e-22
		
		N=zeros([9,9])
		N[6,6]=1/(self.deltad*G)**2
		N[7,7]=1/self.etaS**2
		N[8,8]=1/self.etaS**2 # Vertical and horizontal mosaic spreads are assumed to be equal
		
		W=zeros([9,9])
		W[6,7]=Cy
		W[7,6]=Cy
		W[6,8]=Cz
		W[8,6]=Cz
		W[7,7]=-Cx*G # Correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (55)
		W[8,8]=-Cx*G # Correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (55)
		W*=1e10
		
		# Instrumental part (Popovici)
		
		ki,kf=self.calcKiKf()
		thetaM=deg2rad(A1)
		SM=self.SM
		thetaS=deg2rad(0.5*A4)
		SS=self.SS
		thetaA=deg2rad(A5)
		SA=self.SA
		
		etaM=self.etaM
		etaMv=self.etaM
		etaS=self.etaS
		etaA=self.etaA
		etaAv=self.etaA
		
		if self.guide == 0:
			alpha0=self.guidehdiv*2*pi/ki
			beta0=self.guidevdiv*2*pi/ki
		else:
			alpha0=self.alpha0
			beta0=self.beta0
		
		alpha1=self.alpha1
		alpha2=self.alpha2
		alpha3=self.alpha3

		beta1=self.beta1
		beta2=self.beta2
		beta3=self.beta3
		
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
		
		L0=self.L0
		L1=self.L1
		L2=self.L2
		L3=self.L3
		
		rmh=self.rmh*SM
		rmv=self.rmv
		rah=self.rah*SA
		rav=self.rav
		
		AM=zeros([6,8]) # expressed in AA^-1
		AM[0,0]=0.5*ki/tan(thetaM)
		AM[0,1]=-0.5*ki/tan(thetaM)
		AM[1,1]=ki
		AM[2,3]=ki
		AM[3,4]=0.5*kf/tan(thetaA)
		AM[3,5]=-0.5*kf/tan(thetaA)
		AM[4,4]=kf
		AM[5,6]=kf
		
		FM=zeros([4,4]) # dimensionless
		FM[0,0]=1/etaM**2
		FM[1,1]=1/etaMv**2
		FM[2,2]=1/etaA**2
		FM[3,3]=1/etaAv**2

		GM=zeros([8,8]) # dimensionless
		GM[0,0]=1/alpha0**2
		GM[1,1]=1/alpha1**2
		GM[2,2]=1/beta0**2
		GM[3,3]=1/beta1**2		
		GM[4,4]=1/alpha2**2
		GM[5,5]=1/alpha3**2
		GM[6,6]=1/beta2**2
		GM[7,7]=1/beta3**2		
		
		CM=zeros([4,8]) # dimensionless
		CM[0,0]=0.5
		CM[0,1]=0.5
		CM[2,4]=0.5
		CM[2,5]=0.5
		CM[1,2]=0.5/sin(thetaM)
		CM[1,3]=-0.5/sin(thetaM)
		CM[3,6]=0.5/sin(thetaA)
		CM[3,7]=-0.5/sin(thetaA)
		
		TM=zeros([4,13]) # expressed in cm^-1
		TM[0,0]=-0.5/L0
		TM[0,2]=0.5*cos(thetaM)*(1/L1-1/L0)
		TM[0,3]=0.5*sin(thetaM)*(1/L0+1/L1-2*rmh/abs(sin(thetaM)))
		TM[0,5]=0.5*sin(thetaS)/L1
		TM[0,6]=0.5*cos(thetaS)/L1
		TM[1,1]=-0.5/(L0*sin(thetaM))
		TM[1,4]=0.5*(1/L0+1/L1-2*rmv*sin(thetaM))/(sin(abs(thetaM)))
		TM[1,7]=-0.5/(L1*sin(thetaM))
		TM[2,5]=0.5*sin(thetaS)/L2
		TM[2,6]=-0.5*cos(thetaS)/L2
		TM[2,8]=0.5*cos(thetaA)*(1/L3-1/L2)
		TM[2,9]=0.5*sin(thetaA)*(1/L2+1/L3-2*rah/abs(sin(thetaA)))
		TM[2,11]=0.5/L3
		TM[3,7]=-0.5/(L2*sin(thetaA))
		TM[3,10]=0.5*(1/L2+1/L3-2*rav*sin(thetaA))/abs(sin(thetaA))
		TM[3,12]=-0.5/(L3*sin(thetaA))
		
		DM=zeros([8,13]) # expressed in cm^-1
		DM[0,0]=-1/L0
		DM[0,2]=-cos(thetaM)/L0
		DM[0,3]=sin(thetaM)/L0
		DM[2,1]=-1/L0
		DM[2,4]=1/L0
		DM[1,2]=cos(thetaM)/L1
		DM[1,3]=sin(thetaM)/L1
		DM[1,5]=sin(thetaS)/L1
		DM[1,6]=cos(thetaS)/L1
		DM[3,4]=-1/L1
		DM[3,7]=1/L1
		DM[4,5]=sin(thetaS)/L2
		DM[4,6]=-cos(thetaS)/L2
		DM[4,8]=-cos(thetaA)/L2
		DM[4,9]=sin(thetaA)/L2
		DM[6,7]=-1/L2
		DM[6,10]=1/L2
		DM[5,8]=cos(thetaA)/L3
		DM[5,9]=sin(thetaA)/L3
		DM[5,11]=1/L3
		DM[7,10]=-1/L3
		DM[7,12]=1/L3
		
		# Spatial covariances matrix
		SC=zeros([13,13])
		
		# Source covariance matrix
		if self.sourceshape == 0:
			SC[0,0]=1/16.*sourced**2
			SC[1,1]=1/16.*sourced**2
		else:
			SC[0,0]=1/12.*sourcey**2
			SC[1,1]=1/12.*sourcez**2
		
		# Monochromator covariance matrix
		Smono=zeros([3,3])
		Smono[0,0]=1/12.*monoth**2
		Smono[1,1]=1/12.*monow**2
		Smono[2,2]=1/12.*monoh**2
		SC[2:5,2:5]=Smono
		
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
		SC[5:8,5:8]=Ssam
				
		# Analyzer covariance matrix
		Sana=zeros([3,3])
		Sana[0,0]=1/12.*anath**2
		Sana[1,1]=1/12.*anaw**2
		Sana[2,2]=1/12.*anah**2
		SC[8:11,8:11]=Sana

		# Detector covariance matrix
		if self.detshape == 0:
			SC[11,11]=1/16.*detd**2
			SC[12,12]=1/16.*detd**2
		else:
			SC[11,11]=1/12.*dety**2
			SC[12,12]=1/12.*detz**2		
			
		# Builds full Popovici matrix
		IM=inv(IC)
		SM=inv(SC)			
		LMSpatial=inv(inv(dot(DM,dot(inv(SM+dot(transpose(TM),dot(FM,TM))),transpose(DM))))+GM)
		LMSpatial=dot(AM,dot(LMSpatial,transpose(AM)))
		LMSpatial=inv(dot(IM,dot(LMSpatial,transpose(IM))))	
		
		# Imaginary part of the total resolution matrix (i.e. Larmor precession)	
		H11=-hbar/2/m_n*1/Nf**2-hbar/2/m_n*(tan(theta2))**2/Nf**2-1/(kf*cos(theta2)*1e10)/Nf
		H12=hbar/2/m_n*1/cos(theta1)**2-hbar/2/m_n*Ni**2/Nf**2/(cos(theta2))**2+ \
			Ni/(ki*cos(theta1)*1e10)-Ni/(kf*cos(theta2)*1e10)*Ni/Nf
		H13=hbar/m_n/Nf**2*Ni/cos(theta2)**2+Ni/(kf*cos(theta2)*1e10)*2/Nf
		J11=0.5*hbar/m_n*1/(cos(theta1))**2
		J12=-hbar/m_n*tan(theta1)/cos(theta1)
		K11=-hbar/2/m_n/cos(theta2)**2
		K12=hbar/m_n*tan(theta2)/cos(theta2)*Ni/Nf
		K13=-hbar/m_n*tan(theta2)/cos(theta2)/Nf
		
		Limag=zeros([6,6])
		Limag[0,0]=2*H11*1e24
		Limag[0,1]=H13*1e22
		Limag[0,3]=K13*1e22	
		Limag[1,0]=H13*1e22
		Limag[1,1]=2*H12*1e20
		Limag[1,2]=J12*1e20
		Limag[1,3]=K12*1e20
		Limag[2,1]=J12*1e20
		Limag[2,2]=2*J11*1e20
		Limag[3,0]=K13*1e22
		Limag[3,1]=K12*1e20
		Limag[3,3]=2*K11*1e20
		Limag[4,4]=hbar/m_n*1e20
		Limag[5,5]=-hbar/m_n*1e20
				
		len=tau.size
		Linst=zeros([6,6],dtype=complex)	
		Lmos=zeros([9,9],dtype=complex)
		Lmostemp=zeros([8,8],dtype=complex)
		Lmos_0=zeros([6,6],dtype=complex)
		T=array([0,0,0,0,0,Cx,Cy*G,Cz*G])
		y_tau_imperf=zeros(len)
		yimperf=zeros(len)
		
		for k in xrange(len):
			
			Linst=LMSpatial-1j*tau[k]*1e-12*Limag
			Lmos=dot(transpose(I),dot(Linst,I))+N-1j*tau[k]*1e-12*W
			Lmostemp=Lmos[1:,1:]
			linterm=abs(exp(-0.5*(tau[k]**2)*(1e-4)*dot(T,dot(inv(Lmostemp),transpose(T)))))
			y_tau_imperf[k]=sqrt(abs(1./det(Lmostemp)))*linterm
				
		# Normalize polarization such that P(tau=0)=1	
		Lmos_0=dot(I.T,dot(LMSpatial,I))+N		
		yimperf=sqrt(abs(det(Lmos_0[1:,1:])))*y_tau_imperf
		
		return yimperf
		
	def computeTotalResolutionPOP(self,tau): 
# 		Calculates resolution curve taking Popovici approximation
# 		for the TAS transmission function
		
		A1,A2,A3,A4,A5,A6,theta1,theta2,ni,nf,fratio,tmin,tmax,phi_S,sigma_S=self.calcAngles()
		
		Ni=(hbar*1e10/m_n)*ni
		Nf=(hbar*1e10/m_n)*nf
		
		HKLtoRec=self.getReciprocalBasis()
		
		# Dispersion curvature matrices
		
		Avec=dot(HKLtoRec,self.Avec)
		A=norm(Avec)
		Bvec=dot(HKLtoRec,self.Bvec)
		Cvec=cross(Avec,Bvec)
		RectoAB=array([Avec,Bvec,Cvec]).T
		H=dot(dot((RectoAB.T),self.H),RectoAB) # H is now expressed in the AB basis
		Qvec=dot(HKLtoRec,self.Qvec)
		Q=norm(Qvec)
		phi_AQ=arcsin(norm(cross(Avec,Qvec))/(A*Q))	
		if dot(cross(Avec,Qvec),Cvec)<0:
			phi_AQ*=-1
		ABtoQ=zeros([3,3])
		ABtoQ[0,0]=cos(phi_AQ)
		ABtoQ[0,1]=sin(phi_AQ)
		ABtoQ[1,0]=-sin(phi_AQ)
		ABtoQ[1,1]=cos(phi_AQ)
		ABtoQ[2,2]=1
		H=dot(dot((ABtoQ.T),H),ABtoQ) # H is now expressed in the scattering plane basis
		H*=eV_to_J*1e-3/hbar*1e-20 # meV.AA^2 -> m^2.s^-1
		
		Theta=zeros([3,6])
		Theta[0,0]=cos(phi_S)
		Theta[0,1]=-sin(phi_S)		# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[0,3]=-cos(sigma_S)	# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[0,4]=sin(sigma_S)		# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[1,0]=sin(phi_S)		# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[1,1]=cos(phi_S)		# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[1,3]=-sin(sigma_S)	# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[1,4]=-cos(sigma_S)	# Sign correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (66)
		Theta[2,2]=1
		Theta[2,5]=-1	

		IC=zeros([6,6])
		IC[0,1]=1/cos(theta1)
		IC[0,2]=-tan(theta1)
		IC[3,0]=-1/(Nf*cos(theta2))
		IC[3,1]=Ni/(Nf*cos(theta2))
		IC[3,3]=-tan(theta2)
		IC[1,2]=1
		IC[2,4]=1
		IC[4,3]=1
		IC[5,5]=1
		
		CUM=dot(IC.T,dot(Theta.T,dot(H,dot(Theta,IC))))
		CUM[1:,1:]*=1e20
		CUM[0,0]*=1e24
		CUM[0,1:]*=1e22
		CUM[1:,0]*=1e22
		
		# Sample imperfections effect matrices
		
		dEdq=self.dEdq*eV_to_J/hbar*1e-13;
		Gmodvec=dot(HKLtoRec,self.Gmodvec)
		Cx=dot(Gmodvec/norm(Gmodvec),Avec/norm(Avec))*dEdq
		Cy=dot(Gmodvec/norm(Gmodvec),Bvec/norm(Bvec))*dEdq
		Cz=dot(Gmodvec/norm(Gmodvec),Cvec/norm(Cvec))*dEdq
		if norm(Cz)>1e-3:
			print '\nSlope is defined out of the scattering plane!!!'
		Gvec=dot(HKLtoRec,self.Gvec)
		G=norm(Gvec)
		
		I=zeros([6,9])
		I[0,0]=1
		I[1,1]=1
		I[2,2]=1
		I[3,3]=1
		I[4,4]=1
		I[5,5]=1
		I[0,6]=-Cx*1e-2
		I[0,7]=-Cy*G*1e-22
		I[0,8]=-Cz*G*1e-22
		
		N=zeros([9,9])
		N[6,6]=1/(self.deltad*G)**2
		N[7,7]=1/self.etaS**2
		N[8,8]=1/self.etaS**2 # Vertical and horizontal mosaic spreads are assumed to be equal
		
		W=zeros([9,9])
		W[6,7]=Cy
		W[7,6]=Cy
		W[6,8]=Cz
		W[8,6]=Cz
		W[7,7]=-Cx*G # Correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (55)
		W[8,8]=-Cx*G # Correction from Groitl et al., J. Appl. Cryst. 36, 1307 (2018), Eq. (55)
		W*=1e10
		
		# Instrumental part (Popovici)
		
		ki,kf=self.calcKiKf()
		thetaM=deg2rad(A1)
		SM=self.SM
		thetaS=deg2rad(0.5*A4)
		SS=self.SS
		thetaA=deg2rad(A5)
		SA=self.SA
		
		etaM=self.etaM
		etaMv=self.etaM
		etaS=self.etaS
		etaA=self.etaA
		etaAv=self.etaA
		
		if self.guide == 0:
			alpha0=self.guidehdiv*2*pi/ki
			beta0=self.guidevdiv*2*pi/ki
		else:
			alpha0=self.alpha0
			beta0=self.beta0
		
		alpha1=self.alpha1
		alpha2=self.alpha2
		alpha3=self.alpha3

		beta1=self.beta1
		beta2=self.beta2
		beta3=self.beta3
		
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
		
		L0=self.L0
		L1=self.L1
		L2=self.L2
		L3=self.L3
		
		rmh=self.rmh*SM
		rmv=self.rmv
		rah=self.rah*SA
		rav=self.rav
		
		AM=zeros([6,8]) # expressed in AA^-1
		AM[0,0]=0.5*ki/tan(thetaM)
		AM[0,1]=-0.5*ki/tan(thetaM)
		AM[1,1]=ki
		AM[2,3]=ki
		AM[3,4]=0.5*kf/tan(thetaA)
		AM[3,5]=-0.5*kf/tan(thetaA)
		AM[4,4]=kf
		AM[5,6]=kf
		
		FM=zeros([4,4]) # dimensionless
		FM[0,0]=1/etaM**2
		FM[1,1]=1/etaMv**2
		FM[2,2]=1/etaA**2
		FM[3,3]=1/etaAv**2

		GM=zeros([8,8]) # dimensionless
		GM[0,0]=1/alpha0**2
		GM[1,1]=1/alpha1**2
		GM[2,2]=1/beta0**2
		GM[3,3]=1/beta1**2		
		GM[4,4]=1/alpha2**2
		GM[5,5]=1/alpha3**2
		GM[6,6]=1/beta2**2
		GM[7,7]=1/beta3**2		
		
		CM=zeros([4,8]) # dimensionless
		CM[0,0]=0.5
		CM[0,1]=0.5
		CM[2,4]=0.5
		CM[2,5]=0.5
		CM[1,2]=0.5/sin(thetaM)
		CM[1,3]=-0.5/sin(thetaM)
		CM[3,6]=0.5/sin(thetaA)
		CM[3,7]=-0.5/sin(thetaA)
		
		TM=zeros([4,13]) # expressed in cm^-1
		TM[0,0]=-0.5/L0
		TM[0,2]=0.5*cos(thetaM)*(1/L1-1/L0)
		TM[0,3]=0.5*sin(thetaM)*(1/L0+1/L1-2*rmh/abs(sin(thetaM)))
		TM[0,5]=0.5*sin(thetaS)/L1
		TM[0,6]=0.5*cos(thetaS)/L1
		TM[1,1]=-0.5/(L0*sin(thetaM))
		TM[1,4]=0.5*(1/L0+1/L1-2*rmv*sin(thetaM))/(sin(abs(thetaM)))
		TM[1,7]=-0.5/(L1*sin(thetaM))
		TM[2,5]=0.5*sin(thetaS)/L2
		TM[2,6]=-0.5*cos(thetaS)/L2
		TM[2,8]=0.5*cos(thetaA)*(1/L3-1/L2)
		TM[2,9]=0.5*sin(thetaA)*(1/L2+1/L3-2*rah/abs(sin(thetaA)))
		TM[2,11]=0.5/L3
		TM[3,7]=-0.5/(L2*sin(thetaA))
		TM[3,10]=0.5*(1/L2+1/L3-2*rav*sin(thetaA))/abs(sin(thetaA))
		TM[3,12]=-0.5/(L3*sin(thetaA))
		
		DM=zeros([8,13]) # expressed in cm^-1
		DM[0,0]=-1/L0
		DM[0,2]=-cos(thetaM)/L0
		DM[0,3]=sin(thetaM)/L0
		DM[2,1]=-1/L0
		DM[2,4]=1/L0
		DM[1,2]=cos(thetaM)/L1
		DM[1,3]=sin(thetaM)/L1
		DM[1,5]=sin(thetaS)/L1
		DM[1,6]=cos(thetaS)/L1
		DM[3,4]=-1/L1
		DM[3,7]=1/L1
		DM[4,5]=sin(thetaS)/L2
		DM[4,6]=-cos(thetaS)/L2
		DM[4,8]=-cos(thetaA)/L2
		DM[4,9]=sin(thetaA)/L2
		DM[6,7]=-1/L2
		DM[6,10]=1/L2
		DM[5,8]=cos(thetaA)/L3
		DM[5,9]=sin(thetaA)/L3
		DM[5,11]=1/L3
		DM[7,10]=-1/L3
		DM[7,12]=1/L3
		
		# Spatial covariances matrix
		SC=zeros([13,13])
		
		# Source covariance matrix
		if self.sourceshape == 0:
			SC[0,0]=1/16.*sourced**2
			SC[1,1]=1/16.*sourced**2
		else:
			SC[0,0]=1/12.*sourcey**2
			SC[1,1]=1/12.*sourcez**2
		
		# Monochromator covariance matrix
		Smono=zeros([3,3])
		Smono[0,0]=1/12.*monoth**2
		Smono[1,1]=1/12.*monow**2
		Smono[2,2]=1/12.*monoh**2
		SC[2:5,2:5]=Smono
		
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
		SC[5:8,5:8]=Ssam
				
		# Analyzer covariance matrix
		Sana=zeros([3,3])
		Sana[0,0]=1/12.*anath**2
		Sana[1,1]=1/12.*anaw**2
		Sana[2,2]=1/12.*anah**2
		SC[8:11,8:11]=Sana

		# Detector covariance matrix
		if self.detshape == 0:
			SC[11,11]=1/16.*detd**2
			SC[12,12]=1/16.*detd**2
		else:
			SC[11,11]=1/12.*dety**2
			SC[12,12]=1/12.*detz**2		
			
		# Builds full Popovici matrix
		IM=inv(IC)
		SM=inv(SC)			
		LMSpatial=inv(inv(dot(DM,dot(inv(SM+dot(transpose(TM),dot(FM,TM))),transpose(DM))))+GM)
		LMSpatial=dot(AM,dot(LMSpatial,transpose(AM)))
		LMSpatial=inv(dot(IM,dot(LMSpatial,transpose(IM))))

		# Imaginary part of the total resolution matrix (i.e. Larmor precession)	
		H11=-hbar/2/m_n*1/Nf**2-hbar/2/m_n*(tan(theta2))**2/Nf**2-1/(kf*cos(theta2)*1e10)/Nf
		H12=hbar/2/m_n*1/cos(theta1)**2-hbar/2/m_n*Ni**2/Nf**2/(cos(theta2))**2+ \
			Ni/(ki*cos(theta1)*1e10)-Ni/(kf*cos(theta2)*1e10)*Ni/Nf
		H13=hbar/m_n/Nf**2*Ni/cos(theta2)**2+Ni/(kf*cos(theta2)*1e10)*2/Nf
		J11=0.5*hbar/m_n*1/(cos(theta1))**2
		J12=-hbar/m_n*tan(theta1)/cos(theta1)
		K11=-hbar/2/m_n/cos(theta2)**2
		K12=hbar/m_n*tan(theta2)/cos(theta2)*Ni/Nf
		K13=-hbar/m_n*tan(theta2)/cos(theta2)/Nf
		
		Limag=zeros([6,6])
		Limag[0,0]=2*H11*1e24
		Limag[0,1]=H13*1e22
		Limag[0,3]=K13*1e22	
		Limag[1,0]=H13*1e22
		Limag[1,1]=2*H12*1e20
		Limag[1,2]=J12*1e20
		Limag[1,3]=K12*1e20
		Limag[2,1]=J12*1e20
		Limag[2,2]=2*J11*1e20
		Limag[3,0]=K13*1e22
		Limag[3,1]=K12*1e20
		Limag[3,3]=2*K11*1e20
		Limag[4,4]=hbar/m_n*1e20
		Limag[5,5]=-hbar/m_n*1e20
				
		len=tau.size
		Linst=zeros([6,6],dtype=complex)	
		Linst_0=zeros([6,6],dtype=complex)		
		Lc=zeros([6,6],dtype=complex)
		Ltot=zeros([9,9],dtype=complex)
		Ltot_0=zeros([6,6],dtype=complex)
		T=array([0,0,0,0,0,Cx,Cy*G,Cz*G])
		y_tau_tot=zeros(len)
		ytot=zeros(len)
		
		for k in xrange(len):
			
			Linst=LMSpatial-1j*tau[k]*1e-12*Limag
			Lc=Linst+1j*tau[k]*1e-12*CUM
			Ltot=dot(I.T,dot(Lc,I))+N-1j*tau[k]*1e-12*W
			linterm=abs(exp(-0.5*(tau[k]**2)*1e-4*dot(T,dot(inv(Ltot[1:,1:]),T.T))))
			y_tau_tot[k]=sqrt(abs(1./det(Ltot[1:,1:])))*linterm
		
		# Normalize polarization such that P(tau=0)=1	
		Ltot_0=dot(I.T,dot(LMSpatial,I))+N		
		ytot=sqrt(abs(det(Ltot_0[1:,1:])))*y_tau_tot
		
		return ytot
	
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
			
		astarvec=array([astar,0,0]) # Tranformation matrix [h,k,l] -> cartesian [a*,b*,c*]
		bstarvec=array([bstar*cos(gammastar),bstar*sin(gammastar),0])
		cstarvec1=cstar*cos(betastar)
		cstarvec2=cstar*(cos(alphastar)-cos(betastar)*cos(gammastar)) \
			/sin(gammastar)
		cstarvec3=sqrt(cstar**2-cstarvec1**2-cstarvec2**2)
		cstarvec=array([cstarvec1,cstarvec2,cstarvec3])
		
		HKLtoRec=array([astarvec,bstarvec,cstarvec]).T
		return HKLtoRec	

def calcPhiNSE(dk,t,ki,dki):	
	return cos(4*6.29595*t*ki*dk)*exp(-4*log(2.)*(dk/dki)**2)
	
if __name__ == "__main__":
	
	fid=open('serescal.par','r')
	params=fid.readlines()
	fid.close()
	parlist=[float(i) for i in params]
	calc=SERESCAL(parlist)	