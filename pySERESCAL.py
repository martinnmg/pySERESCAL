#-*- encoding: ISO-8859-1 -*-
# N. MARTIN 2022/04/28
# (c) Laboratoire Léon Brillouin LLB
# From the 'SERESCAL' MATLAB toolbox by K. Habicht (HZB) 
# References:
# [Habicht2003] K. Habicht et al., J. Appl. Cryst. 36, 1307-1318 (2003)
# [Groitl2018] F. Groitl et al., J. Appl. Cryst. 51, 818-830 (2018)	

import sys
sys.path.append('./@pySERESCAL/')

from Tkinter import *
from tkFileDialog import *
from PIL import Image,ImageTk
import time
from serescal import *
from rescal import *
import matplotlib.pyplot as plt

class pySERESCAL():

	def __init__(self):
		# Creates the app and fill entries with the content of serescal.par
		
		t=time.strftime('%d/%m/%y %H:%M:%S',time.localtime()).split()
		print '\npySERESCAL started on',t[0],'at',t[1],'\n'
		
		self.root=Tk()
		self.root.title('pySERESCAL (Apr. 2022)')
		self.varlist=list()
		self.methodvar=StringVar()
		self.methodvar.set('cn')
		#self.autoNRSE=IntVar()
		#self.autoNRSE.set(1)
		for i in xrange(87):
			self.varlist.append(DoubleVar())
		
		# Menu
		self.menubar=Menu(self.root)
		
		self.filesubmenu=Menu(self.root)
		self.menubar.add_cascade(label='File',menu=self.filesubmenu)
		self.filesubmenu.add_command(label='Import parameters',command=self.importParameters)
		self.filesubmenu.add_command(label='Save parameters',command=self.saveParameters)
		self.filesubmenu.add_command(label='Quit',command=self.quitAll)
				
		self.methodsubmenu=Menu(self.root)
		self.menubar.add_cascade(label='Method',menu=self.methodsubmenu)
		self.methodsubmenu.add_radiobutton(label='Cooper-Nathans', var=self.methodvar, value='cn')
		self.methodsubmenu.add_radiobutton(label='Popovici', var=self.methodvar, value='pop')
		
		self.rescalsubmenu=Menu(self.root)
		self.menubar.add_cascade(label='TAS',menu=self.rescalsubmenu)
		self.rescalsubmenu.add_command(label='Vizualise TAS ellipsoid (test)',command=self.plotTASellipsoids)		
		
		self.serescalsubmenu=Menu(self.root)
		self.menubar.add_cascade(label='NRSE',menu=self.serescalsubmenu)
		self.serescalsubmenu.add_command(label='Plot NRSE resolution',command=self.plotNRSEresolution)		
		self.serescalsubmenu.add_command(label='Save NRSE resolution',command=self.saveNRSEresolution)
		self.serescalsubmenu.add_command(label='Update NRSE settings',command=self.updateNRSEsettings)
		#self.serescalsubmenu.add_checkbutton(label="Auto update?", onvalue=1, offvalue=0, variable=self.autoNRSE)
		
		self.root.config(menu=self.menubar)
		
		# Spectrometer parameters
		Label(self.root,text='TAS parameters',bg='black',fg='red',font=("Helvetica", 12)).grid(row=0,columnspan=10)
		#
		self.dMlabel=Label(self.root,text=u'dM (\u212b)',bg='black',fg='white').grid(row=1, column=0)
		self.dMvalue=Entry(self.root,textvariable=self.varlist[0],bd=0).grid(row=1, column=1)
		self.dAlabel=Label(self.root,text=u'dA (\u212b)',bg='black',fg='white').grid(row=1, column=2)
		self.dAvalue=Entry(self.root,textvariable=self.varlist[1],bd=0).grid(row=1, column=3)
		self.etaMlabel=Label(self.root,text=u'\u03b7M (arcmin)',bg='black',fg='white').grid(row=1, column=4)
		self.etaMvalue=Entry(self.root,textvariable=self.varlist[2],bd=0).grid(row=1, column=5)
		self.etaAlabel=Label(self.root,text=u'\u03b7A (arcmin)',bg='black',fg='white').grid(row=1, column=6)
		self.etaAvalue=Entry(self.root,textvariable=self.varlist[3],bd=0).grid(row=1, column=7)
		#
		self.SMlabel=Label(self.root,text='SM',bg='black',fg='white').grid(row=2, column=0)
		self.SMvalue=Entry(self.root,textvariable=self.varlist[6],bd=0).grid(row=2, column=1)
		self.SSlabel=Label(self.root,text='SS',bg='black',fg='white').grid(row=2, column=2)
		self.SSvalue=Entry(self.root,textvariable=self.varlist[7],bd=0).grid(row=2, column=3)
		self.SAlabel=Label(self.root,text='SA',bg='black',fg='white').grid(row=2, column=4)
		self.SAvalue=Entry(self.root,textvariable=self.varlist[8],bd=0).grid(row=2, column=5)
		#
		self.kfixlabel=Label(self.root,text=u'k (\u212b^-1)',bg='black',fg='white').grid(row=3, column=0)
		self.kfixvalue=Entry(self.root,textvariable=self.varlist[9],bd=0).grid(row=3, column=1)		
		#self.fixlabel=Label(self.root,text='Fixed k (1 = ki; 2 = kf)',bg='black',fg='white').grid(row=3, column=2)
		#self.fixvalue=Entry(self.root,textvariable=self.varlist[10],bd=0).grid(row=3, column=3)
		self.fixlabel=Label(self.root,text='Fixed k',bg='black',fg='white').grid(row=3, column=2)
		self.fixvar=StringVar()
		self.fixvar.set('kf')
		self.fixchoice=OptionMenu(self.root,self.fixvar,'ki','kf',command=self.updateGUIkfix).grid(row=3, column=3, sticky='ew')
		#
		self.guidelabel=Label(self.root,text='In-pile guide?',bg='black',fg='white').grid(row=4, column=0)
		self.guidevar=StringVar()
		self.guidevar.set('Yes')
		self.guidechoice=OptionMenu(self.root,self.guidevar,'Yes','No',command=self.updateGUIguide).grid(row=4, column=1, sticky='ew')
		self.varlist[58].set(0)
		self.guide1label=Label(self.root,text=u'Hor. div. (arcmin/\u212b)',bg='black',fg='white').grid(row=4, column=2)
		self.guide1value=Entry(self.root,textvariable=self.varlist[59],bd=0).grid(row=4, column=3)
		self.guide2label=Label(self.root,text=u'Vert. div. (arcmin/\u212b)',bg='black',fg='white').grid(row=4, column=4)
		self.guide2value=Entry(self.root,textvariable=self.varlist[60],bd=0).grid(row=4, column=5)	
		#
		self.alpha1label=Label(self.root,text=u'\u03b10 (arcmin)',bg='black',fg='white').grid(row=5, column=0)
		self.alpha1value=Entry(self.root,textvariable=self.varlist[11],bd=0).grid(row=5, column=1)
		self.alpha2label=Label(self.root,text=u'\u03b11 (arcmin)',bg='black',fg='white').grid(row=5, column=2)
		self.alpha2value=Entry(self.root,textvariable=self.varlist[12],bd=0).grid(row=5, column=3)
		self.alpha3label=Label(self.root,text=u'\u03b12 (arcmin)',bg='black',fg='white').grid(row=5, column=4)
		self.alpha3value=Entry(self.root,textvariable=self.varlist[13],bd=0).grid(row=5, column=5)
		self.alpha4label=Label(self.root,text=u'\u03b13 (arcmin)',bg='black',fg='white').grid(row=5, column=6)
		self.alpha4value=Entry(self.root,textvariable=self.varlist[14],bd=0).grid(row=5, column=7)
		#
		self.beta1label=Label(self.root,text=u'\u03b20 (arcmin)',bg='black',fg='white').grid(row=6, column=0)
		self.beta1value=Entry(self.root,textvariable=self.varlist[15],bd=0).grid(row=6, column=1)
		self.beta2label=Label(self.root,text=u'\u03b21 (arcmin)',bg='black',fg='white').grid(row=6, column=2)
		self.beta2value=Entry(self.root,textvariable=self.varlist[16],bd=0).grid(row=6, column=3)
		self.beta3label=Label(self.root,text=u'\u03b22 (arcmin)',bg='black',fg='white').grid(row=6, column=4)
		self.beta3value=Entry(self.root,textvariable=self.varlist[17],bd=0).grid(row=6, column=5)
		self.beta4label=Label(self.root,text=u'\u03b23 (arcmin)',bg='black',fg='white').grid(row=6, column=6)
		self.beta4value=Entry(self.root,textvariable=self.varlist[18],bd=0).grid(row=6, column=7)
		
		# Sample parameters
		Label(self.root,text='Sample parameters',bg='black',fg='red',font=("Helvetica", 12)).grid(row=7,columnspan=10)
		#
		self.alabel=Label(self.root,text=u'a (\u212b)',bg='black',fg='white').grid(row=8, column=0)
		self.avalue=Entry(self.root,textvariable=self.varlist[19],bd=0).grid(row=8, column=1)
		self.blabel=Label(self.root,text=u'b (\u212b)',bg='black',fg='white').grid(row=8, column=2)
		self.bvalue=Entry(self.root,textvariable=self.varlist[20],bd=0).grid(row=8, column=3)
		self.clabel=Label(self.root,text=u'c (\u212b)',bg='black',fg='white').grid(row=8, column=4)
		self.cvalue=Entry(self.root,textvariable=self.varlist[21],bd=0).grid(row=8, column=5)
		self.deltaDlabel=Label(self.root,text=u'\u0394d/d',bg='black',fg='white').grid(row=8, column=6)
		self.deltaDvalue=Entry(self.root,textvariable=self.varlist[5],bd=0).grid(row=8, column=7)	
		#
		self.alphalabel=Label(self.root,text=u'\u03b1 (deg)',bg='black',fg='white').grid(row=9, column=0)
		self.alphavalue=Entry(self.root,textvariable=self.varlist[22],bd=0).grid(row=9, column=1)
		self.betalabel=Label(self.root,text=u'\u03b2 (deg)',bg='black',fg='white').grid(row=9, column=2)
		self.betavalue=Entry(self.root,textvariable=self.varlist[23],bd=0).grid(row=9, column=3)
		self.gammalabel=Label(self.root,text=u'\u03b3 (deg)',bg='black',fg='white').grid(row=9, column=4)
		self.gammavalue=Entry(self.root,textvariable=self.varlist[24],bd=0).grid(row=9, column=5)
		self.etaSlabel=Label(self.root,text=u'\u03b7S (arcmin)',bg='black',fg='white').grid(row=9, column=6)
		self.etaSvalue=Entry(self.root,textvariable=self.varlist[4],bd=0).grid(row=9, column=7)
		#
		self.axlabel=Label(self.root,text='Ax',bg='black',fg='white').grid(row=10, column=0)
		self.axvalue=Entry(self.root,textvariable=self.varlist[25],bd=0).grid(row=10, column=1)
		self.aylabel=Label(self.root,text='Ay',bg='black',fg='white').grid(row=10, column=2)
		self.ayvalue=Entry(self.root,textvariable=self.varlist[26],bd=0).grid(row=10, column=3)
		self.azlabel=Label(self.root,text='Az',bg='black',fg='white').grid(row=10, column=4)
		self.azvalue=Entry(self.root,textvariable=self.varlist[27],bd=0).grid(row=10, column=5)
		#
		self.bxlabel=Label(self.root,text='Bx',bg='black',fg='white').grid(row=11, column=0)
		self.bxvalue=Entry(self.root,textvariable=self.varlist[28],bd=0).grid(row=11, column=1)
		self.bylabel=Label(self.root,text='By',bg='black',fg='white').grid(row=11, column=2)
		self.byvalue=Entry(self.root,textvariable=self.varlist[29],bd=0).grid(row=11, column=3)
		self.bzlabel=Label(self.root,text='Bz',bg='black',fg='white').grid(row=11, column=4)
		self.bzvalue=Entry(self.root,textvariable=self.varlist[30],bd=0).grid(row=11, column=5)
		
		# Excitation parameters
		Label(self.root,text='Excitation parameters',bg='black',fg='red',font=("Helvetica", 12)).grid(row=12,columnspan=10)
		#
		self.Qhlabel=Label(self.root,text='Qh (r.l.u.)',bg='black',fg='white').grid(row=13, column=0)
		self.Qhvalue=Entry(self.root,textvariable=self.varlist[31],bd=0).grid(row=13, column=1)
		self.Qklabel=Label(self.root,text='Qk (r.l.u.)',bg='black',fg='white').grid(row=13, column=2)
		self.Qkvalue=Entry(self.root,textvariable=self.varlist[32],bd=0).grid(row=13, column=3)
		self.Qllabel=Label(self.root,text='Ql (r.l.u.)',bg='black',fg='white').grid(row=13, column=4)
		self.Qlvalue=Entry(self.root,textvariable=self.varlist[33],bd=0).grid(row=13, column=5)
		self.Enlabel=Label(self.root,text='E (meV)',bg='black',fg='white').grid(row=13, column=6)
		self.Envalue=Entry(self.root,textvariable=self.varlist[34],bd=0).grid(row=13, column=7)
		#
		self.Ghlabel=Label(self.root,text=u'\u2207h (r.l.u.)',bg='black',fg='white').grid(row=14, column=0)
		self.Ghvalue=Entry(self.root,textvariable=self.varlist[35],bd=0).grid(row=14, column=1)
		self.Gklabel=Label(self.root,text=u'\u2207k (r.l.u.)',bg='black',fg='white').grid(row=14, column=2)
		self.Gkvalue=Entry(self.root,textvariable=self.varlist[36],bd=0).grid(row=14, column=3)
		self.Gllabel=Label(self.root,text=u'\u2207l (r.l.u.)',bg='black',fg='white').grid(row=14, column=4)
		self.Glvalue=Entry(self.root,textvariable=self.varlist[37],bd=0).grid(row=14, column=5)		
		self.dEdQlabel=Label(self.root,text=u'\u2202E/\u2202q (meV.\u212b)',bg='black',fg='white').grid(row=14, column=6)
		self.dEdQvalue=Entry(self.root,textvariable=self.varlist[38],bd=0).grid(row=14, column=7)	
		#
		self.qhlabel=Label(self.root,text='qh (r.l.u.)',bg='black',fg='white').grid(row=15, column=0)
		self.qhvalue=Entry(self.root,textvariable=self.varlist[39],bd=0).grid(row=15, column=1)
		self.qklabel=Label(self.root,text='qk (r.l.u.)',bg='black',fg='white').grid(row=15, column=2)
		self.qkvalue=Entry(self.root,textvariable=self.varlist[40],bd=0).grid(row=15, column=3)
		self.qllabel=Label(self.root,text='ql (r.l.u.)',bg='black',fg='white').grid(row=15, column=4)
		self.qlvalue=Entry(self.root,textvariable=self.varlist[41],bd=0).grid(row=15, column=5)
		self.Egaplabel=Label(self.root,text=u'\u0394E (meV)',bg='black',fg='white').grid(row=15, column=6)
		self.Egapvalue=Entry(self.root,textvariable=self.varlist[85],bd=0).grid(row=15, column=7)		
		#
		self.extypelabel=Label(self.root,text='Type',bg='black',fg='white').grid(row=13, column=8)
		self.exptype=StringVar()
		self.exptype.set('Phonon')
		self.extypechoice=OptionMenu(self.root,self.exptype,'Phonon','Magnon',command=self.updateGUIextype).grid(row=13, column=9, sticky='ew')
		self.varlist[82].set(1)		
		self.anisolabel=Label(self.root,text='Mxy',bg='black',fg='white').grid(row=14, column=8)
		self.extypevalue=Entry(self.root,textvariable=self.varlist[83],bd=0).grid(row=14, column=9)
		self.anisolabel=Label(self.root,text='Mz',bg='black',fg='white').grid(row=15, column=8)
		self.extypevalue=Entry(self.root,textvariable=self.varlist[84],bd=0).grid(row=15, column=9)
		
		# Curvature matrix		
		Label(self.root,text='Curvature matrix',bg='black',fg='red',font=("Helvetica", 12)).grid(row=16,columnspan=10)
		#
		self.Hxxlabel=Label(self.root,text=u'Hxx (meV.\u212b^2)',bg='black',fg='white').grid(row=17, column=0)
		self.Hxxvalue=Entry(self.root,textvariable=self.varlist[42],bd=0).grid(row=17, column=1)
		self.Hxylabel=Label(self.root,text=u'Hxy (meV.\u212b^2)',bg='black',fg='white').grid(row=17, column=2)
		self.Hxyvalue=Entry(self.root,textvariable=self.varlist[43],bd=0).grid(row=17, column=3)
		self.Hxzlabel=Label(self.root,text=u'Hxz (meV.\u212b^2)',bg='black',fg='white').grid(row=17, column=4)
		self.Hxzvalue=Entry(self.root,textvariable=self.varlist[44],bd=0).grid(row=17, column=5)
		#
		self.Hyxlabel=Label(self.root,text=u'Hyx (meV.\u212b^2)',bg='black',fg='white').grid(row=18, column=0)
		self.Hyxvalue=Entry(self.root,textvariable=self.varlist[45],bd=0).grid(row=18, column=1)
		self.Hyylabel=Label(self.root,text=u'Hyy (meV.\u212b^2)',bg='black',fg='white').grid(row=18, column=2)
		self.Hyyvalue=Entry(self.root,textvariable=self.varlist[46],bd=0).grid(row=18, column=3)
		self.Hyzlabel=Label(self.root,text=u'Hyz (meV.\u212b^2)',bg='black',fg='white').grid(row=18, column=4)
		self.Hyzvalue=Entry(self.root,textvariable=self.varlist[47],bd=0).grid(row=18, column=5)
		#
		self.Hzxlabel=Label(self.root,text=u'Hzx (meV.\u212b^2)',bg='black',fg='white').grid(row=19, column=0)
		self.Hzxvalue=Entry(self.root,textvariable=self.varlist[48],bd=0).grid(row=19, column=1)
		self.Hzylabel=Label(self.root,text=u'Hzy (meV.\u212b^2)',bg='black',fg='white').grid(row=19, column=2)
		self.Hzyvalue=Entry(self.root,textvariable=self.varlist[49],bd=0).grid(row=19, column=3)
		self.Hzzlabel=Label(self.root,text=u'Hzz (meV.\u212b^2)',bg='black',fg='white').grid(row=19, column=4)
		self.Hzzvalue=Entry(self.root,textvariable=self.varlist[50],bd=0).grid(row=19, column=5)
		
		# self.cumbasislabel=Label(self.root,text='Basis choice',bg='black',fg='white').grid(row=17, column=6)
		# self.cumbasis=StringVar()
		# self.cumbasis.set('[a*,b*,c*]')		
		# self.cumbasischoice=OptionMenu(self.root,self.cumbasis,'[a*,b*,c*]','[Qx,Qy,Qz]',command=self.updateCurvatureMatrixBasis).grid(row=17, column=7, sticky='ew')
		# self.varlist[86].set(0)
		
		# Spatial parameters
		Label(self.root,text='Spatial parameters',bg='black',fg='red',font=("Helvetica", 12)).grid(row=20,columnspan=10)
		#
		self.sourcelabel=Label(self.root,text='Source shape',bg='black',fg='white').grid(row=21, column=0)
		self.sourcevar=StringVar()
		self.sourcevar.set('Rectangular')
		self.sourcechoice=OptionMenu(self.root,self.sourcevar,'Circular','Rectangular',command=self.updateGUIsource).grid(row=21, column=1, sticky='ew')
		self.varlist[55].set(1)
		self.source1label=Label(self.root,text='Width (cm)',bg='black',fg='white').grid(row=21, column=2)
		self.source1value=Entry(self.root,textvariable=self.varlist[56],bd=0).grid(row=21, column=3)
		self.source2label=Label(self.root,text='Height (cm)',bg='black',fg='white').grid(row=21, column=4)
		self.source2value=Entry(self.root,textvariable=self.varlist[57],bd=0).grid(row=21, column=5)	
		#
		self.samplelabel=Label(self.root,text='Sample shape',bg='black',fg='white').grid(row=22, column=0)
		self.samplevar=StringVar()
		self.samplevar.set('Cylindrical')
		self.samplechoice=OptionMenu(self.root,self.samplevar,'Cylindrical','Cuboidal',command=self.updateGUIsample).grid(row=22, column=1, sticky='ew')
		self.varlist[61].set(0)
		self.sample1label=Label(self.root,text='Diameter (cm)',bg='black',fg='white').grid(row=22, column=2)
		self.sample1value=Entry(self.root,textvariable=self.varlist[62],bd=0).grid(row=22, column=3)
		self.sample2label=Label(self.root,text='Height (cm)',bg='black',fg='white').grid(row=22, column=4)
		self.sample2value=Entry(self.root,textvariable=self.varlist[63],bd=0).grid(row=22, column=5)	
		self.sample3label=Label(self.root,text='--',bg='black',fg='black').grid(row=22, column=6)
		self.sample3value=Label(self.root,text='--',bg='black',fg='black').grid(row=22, column=7)
		#
		self.detlabel=Label(self.root,text='Detector shape',bg='black',fg='white').grid(row=23, column=0)
		self.detvar=StringVar()
		self.detvar.set('Rectangular')
		self.detchoice=OptionMenu(self.root,self.detvar,'Circular','Rectangular',command=self.updateGUIdetector).grid(row=23, column=1, sticky='ew')
		self.varlist[65].set(1)
		self.det1label=Label(self.root,text='Width (cm)',bg='black',fg='white').grid(row=23, column=2)
		self.det1value=Entry(self.root,textvariable=self.varlist[66],bd=0).grid(row=23, column=3)
		self.det2label=Label(self.root,text='Height (cm)',bg='black',fg='white').grid(row=23, column=4)
		self.det2value=Entry(self.root,textvariable=self.varlist[67],bd=0).grid(row=23, column=5)	
		#
		self.monolabel=Label(self.root,text='Monochromator',bg='black',fg='white').grid(row=24, column=0)
		self.monodepthlabel=Label(self.root,text='Depth (cm)',bg='black',fg='white').grid(row=24, column=1)
		self.monodepthvalue=Entry(self.root,textvariable=self.varlist[68],bd=0).grid(row=24, column=2)
		self.monowidthlabel=Label(self.root,text='Width (cm)',bg='black',fg='white').grid(row=24, column=3)
		self.monowidthvalue=Entry(self.root,textvariable=self.varlist[69],bd=0).grid(row=24, column=4)
		self.monoheightlabel=Label(self.root,text='Height (cm)',bg='black',fg='white').grid(row=24, column=5)
		self.monoheightvalue=Entry(self.root,textvariable=self.varlist[70],bd=0).grid(row=24, column=6)		
		#
		self.analabel=Label(self.root,text='Analyzer',bg='black',fg='white').grid(row=25, column=0)
		self.anadepthlabel=Label(self.root,text='Depth (cm)',bg='black',fg='white').grid(row=25, column=1)
		self.anadepthvalue=Entry(self.root,textvariable=self.varlist[71],bd=0).grid(row=25, column=2)
		self.anawidthlabel=Label(self.root,text='Width (cm)',bg='black',fg='white').grid(row=25, column=3)
		self.anawidthvalue=Entry(self.root,textvariable=self.varlist[72],bd=0).grid(row=25, column=4)
		self.anaheightlabel=Label(self.root,text='Height (cm)',bg='black',fg='white').grid(row=25, column=5)
		self.anaheightvalue=Entry(self.root,textvariable=self.varlist[73],bd=0).grid(row=25, column=6)
		#
		self.distlabel=Label(self.root,text='Distances',bg='black',fg='white').grid(row=26, column=0)
		self.L1label=Label(self.root,text='L0 (cm)',bg='black',fg='white').grid(row=26, column=1)
		self.L1value=Entry(self.root,textvariable=self.varlist[74],bd=0).grid(row=26, column=2)
		self.L2label=Label(self.root,text='L1 (cm)',bg='black',fg='white').grid(row=26, column=3)
		self.L2value=Entry(self.root,textvariable=self.varlist[75],bd=0).grid(row=26, column=4)
		self.L3label=Label(self.root,text='L2 (cm)',bg='black',fg='white').grid(row=26, column=5)
		self.L3value=Entry(self.root,textvariable=self.varlist[76],bd=0).grid(row=26, column=6)
		self.L4label=Label(self.root,text='L3 (cm)',bg='black',fg='white').grid(row=26, column=7)
		self.L4value=Entry(self.root,textvariable=self.varlist[77],bd=0).grid(row=26, column=8)
		#
		self.focuslabel=Label(self.root,text='Focus',bg='black',fg='white').grid(row=27, column=0)
		self.rmhlabel=Label(self.root,text='Rmh (1/m)',bg='black',fg='white').grid(row=27, column=1)
		self.rmhvalue=Entry(self.root,textvariable=self.varlist[78],bd=0).grid(row=27, column=2)
		self.rmvlabel=Label(self.root,text='Rmv (1/m)',bg='black',fg='white').grid(row=27, column=3)
		self.rmwvalue=Entry(self.root,textvariable=self.varlist[79],bd=0).grid(row=27, column=4)
		self.rahlabel=Label(self.root,text='Rah (1/m)',bg='black',fg='white').grid(row=27, column=5)
		self.rahvalue=Entry(self.root,textvariable=self.varlist[80],bd=0).grid(row=27, column=6)
		self.ravlabel=Label(self.root,text='Rav (1/m)',bg='black',fg='white').grid(row=27, column=7)
		self.ravvalue=Entry(self.root,textvariable=self.varlist[81],bd=0).grid(row=27, column=8)
		
		# NRSE parameters
		Label(self.root,text='NRSE settings',bg='black',fg='red',font=("Helvetica", 12)).grid(row=28,columnspan=10)
		#
		self.fminlabel=Label(self.root,text='fmin (kHz)',bg='black',fg='white').grid(row=29, column=0)
		self.fminvalue=Entry(self.root,textvariable=self.varlist[51],bd=0).grid(row=29, column=1)
		self.fmaxlabel=Label(self.root,text='fmax (kHz)',bg='black',fg='white').grid(row=29, column=2)
		self.fmaxvalue=Entry(self.root,textvariable=self.varlist[52],bd=0).grid(row=29, column=3)
		self.L4pilabel=Label(self.root,text='L4pi (m)',bg='black',fg='white').grid(row=29, column=4)
		self.L4pivalue=Entry(self.root,textvariable=self.varlist[53],bd=0).grid(row=29, column=5)
		self.L8pilabel=Label(self.root,text='L8pi (m)',bg='black',fg='white').grid(row=29, column=6)
		self.L8pivalue=Entry(self.root,textvariable=self.varlist[54],bd=0).grid(row=29, column=7)
		
		# NRSE settings
		fid=open('pySERESCAL.ini','r')
		params=fid.readlines()
		fid.close()
		params=[float(i) for i in params]
		for i in xrange(params.__len__()):
			self.varlist[i].set(params[i])
		NRSEcalc=SERESCAL(params,self.methodvar.get())
		self.setlist=list()
		for i in xrange(5):
			self.setlist.append(DoubleVar())
		A1,A2,A3,A4,A5,A6,theta1,theta2,ni,nf,fratio,tmin,tmax,phi_S,sigma_S=NRSEcalc.calcAngles()
		#
		self.setlist[0].set(int(1e2*rad2deg(theta1))/1.e2)	
		self.setlist[1].set(int(1e2*rad2deg(theta2))/1.e2)
		self.setlist[2].set(int(1e3*fratio)/1.e3)
		self.setlist[3].set(int(1e2*tmin*1e12)/1.e2)
		self.setlist[4].set(int(1e2*tmax*1e12)/1.e2)
		#
		self.Zetailabel=Label(self.root,text=u'\u03B6i (deg)',bg='black',fg='white').grid(row=30, column=0)
		self.Zetaivalue=Entry(self.root,textvariable=self.setlist[0]).grid(row=30, column=1)
		self.Zetaflabel=Label(self.root,text=u'\u03B6f (deg)',bg='black',fg='white').grid(row=30, column=2)
		self.Zetafvalue=Entry(self.root,textvariable=self.setlist[1]).grid(row=30, column=3)
		self.fratiolabel=Label(self.root,text='Freq. ratio',bg='black',fg='white').grid(row=30, column=4)
		self.fratiovalue=Entry(self.root,textvariable=self.setlist[2]).grid(row=30, column=5)
		#
		self.tminlabel=Label(self.root,text=u'\u03C4min (ps)',bg='black',fg='white').grid(row=31, column=0)
		self.tminvalue=Entry(self.root,textvariable=self.setlist[3]).grid(row=31, column=1)
		self.tmaxlabel=Label(self.root,text=u'\u03C4max (ps)',bg='black',fg='white').grid(row=31, column=2)
		self.tmaxvalue=Entry(self.root,textvariable=self.setlist[4]).grid(row=31, column=3)
		
		# Logo
		img=ImageTk.PhotoImage(Image.open('./@pySERESCAL/logo.png').resize((150, 154), Image.ANTIALIAS))
		self.pySRlogo=Label(self.root, image=img).grid(column=8,row=1,columnspan=3,rowspan=10)	
		
		# Start application
		self.root.configure(background='black')
		self.root.resizable(0,0)
		self.root.mainloop()

	def updateGUIextype(self,extype):
	
		if extype == 'Magnon':
			self.varlist[82].set(0)
		if extype == 'Phonon':
			self.varlist[82].set(1)

		self.root.update()

	def updateGUIkfix(self,fix):
		
		if fix=='ki':
			self.varlist[10].set(1)
		else:
			self.varlist[10].set(2)

		self.root.update()	

	def updateGUIguide(self,isthereguide):
	
		if isthereguide == 'Yes':
			self.varlist[58].set(0)
			self.guidelabel=Label(self.root,text='In-pile guide?',bg='black',fg='white').grid(row=4, column=0)
			self.guidechoice=OptionMenu(self.root,self.guidevar,'Yes','No',command=self.updateGUIguide).grid(row=4, column=1, sticky='ew')
			self.guide1label=Label(self.root,text='Hor. div. (min/AA)',bg='black',fg='white').grid(row=4, column=2)
			self.guide1value=Entry(self.root,textvariable=self.varlist[59],bd=0).grid(row=4, column=3)
			self.guide2label=Label(self.root,text='Vert. div. (min/AA)',bg='black',fg='white').grid(row=4, column=4)
			self.guide2value=Entry(self.root,textvariable=self.varlist[60],bd=0).grid(row=4, column=5)	
		if isthereguide == 'No':
			self.varlist[58].set(1)
			self.guidelabel=Label(self.root,text='In-pile guide?',bg='black',fg='white').grid(row=4, column=0)
			self.guidechoice=OptionMenu(self.root,self.guidevar,'Yes','No',command=self.updateGUIguide).grid(row=4, column=1, sticky='ew')
			self.guide1label=Label(self.root,text='-',bg='black',fg='black').grid(row=4, column=2, sticky='ew')		
			self.guide1value=Label(self.root,text='-',bg='black',fg='black').grid(row=4, column=3, sticky='ew')	
			self.guide2label=Label(self.root,text='-',bg='black',fg='black').grid(row=4, column=4, sticky='ew')	
			self.guide2value=Label(self.root,text='-',bg='black',fg='black').grid(row=4, column=5, sticky='ew')

		self.root.update()	
		
	def updateGUIsource(self,sourcetype):
		
		if sourcetype == 'Rectangular':
			self.varlist[55].set(1)
			self.sourcelabel=Label(self.root,text='Source shape',bg='black',fg='white').grid(row=22, column=0)
			self.sourcechoice=OptionMenu(self.root,self.sourcevar,'Circular','Rectangular',command=self.updateGUIsource).grid(row=22, column=1, sticky='ew')
			Label(self.root,text='Diameter (cm)',bg='black',fg='black').grid(row=22, column=2, sticky='ew')
			self.source1label=Label(self.root,text='Width (cm)',bg='black',fg='white').grid(row=22, column=2)
			self.source1value=Entry(self.root,textvariable=self.varlist[56],bd=0).grid(row=22, column=3)
			self.source2label=Label(self.root,text='Height (cm)',bg='black',fg='white').grid(row=22, column=4)
			self.source2value=Entry(self.root,textvariable=self.varlist[57],bd=0).grid(row=22, column=5)	
		if sourcetype == 'Circular':
			self.varlist[55].set(0)
			self.sourcelabel=Label(self.root,text='Source shape',bg='black',fg='white').grid(row=22, column=0)
			self.sourcechoice=OptionMenu(self.root,self.sourcevar,'Circular','Rectangular',command=self.updateGUIsource).grid(row=22, column=1, sticky='ew')
			self.source1label=Label(self.root,text='Diameter (cm)',bg='black',fg='white').grid(row=22, column=2)
			self.source1value=Entry(self.root,textvariable=self.varlist[56],bd=0).grid(row=22, column=3)
			Label(self.root,text='-',bg='black',fg='black').grid(row=22, column=4, sticky='ew')
			Label(self.root,text='-',bg='black',fg='black').grid(row=22, column=5, sticky='ew')
		
		self.root.update()		
		
	def updateGUIsample(self,sampletype):

		if self.samplevar.get() == 'Cylindrical':
			self.varlist[61].set(0)
			self.samplelabel=Label(self.root,text='Sample shape',bg='black',fg='white').grid(row=23, column=0)
			self.samplechoice=OptionMenu(self.root,self.samplevar,'Cylindrical','Cuboidal',command=self.updateGUIsample).grid(row=23, column=1, sticky='ew')
			Label(self.root,text='-',bg='black',fg='black').grid(row=23, column=2, sticky='ew')
			self.sample1label=Label(self.root,text='Diameter (cm)',bg='black',fg='white').grid(row=23, column=2)
			self.sample1value=Entry(self.root,textvariable=self.varlist[62],bd=0).grid(row=23, column=3)
			self.sample2label=Label(self.root,text='Height (cm)',bg='black',fg='white').grid(row=23, column=4)
			self.sample2value=Entry(self.root,textvariable=self.varlist[63],bd=0).grid(row=23, column=5)
			self.sample3label=Label(self.root,text='-',bg='black',fg='black').grid(row=23, column=6, sticky='ew')
			self.sample3value=Label(self.root,text='-',bg='black',fg='black').grid(row=23, column=7, sticky='ew')
		if sampletype == 'Cuboidal':
			self.varlist[61].set(1)
			self.samplelabel=Label(self.root,text='Sample shape',bg='black',fg='white').grid(row=23, column=0)
			self.samplechoice=OptionMenu(self.root,self.samplevar,'Cylindrical','Cuboidal',command=self.updateGUIsample).grid(row=23, column=1, sticky='ew')
			self.sample1label=Label(self.root,text='Thickness (cm)',bg='black',fg='white').grid(row=23, column=2)
			self.sample1value=Entry(self.root,textvariable=self.varlist[62],bd=0).grid(row=23, column=3)
			Label(self.root,text='-',bg='black',fg='black').grid(row=23, column=4, sticky='ew')
			self.sample2label=Label(self.root,text='Width (cm)',bg='black',fg='white').grid(row=23, column=4)
			self.sample2value=Entry(self.root,textvariable=self.varlist[63],bd=0).grid(row=23, column=5)			
			self.sample3label=Label(self.root,text='Height (cm)',bg='black',fg='white').grid(row=23, column=6)
			self.sample3value=Entry(self.root,textvariable=self.varlist[64],bd=0).grid(row=23, column=7)

		self.root.update()			

	def updateGUIdetector(self,dettype):
	
		if dettype == 'Rectangular':
			self.varlist[65].set(1)
			self.detlabel=Label(self.root,text='Detector shape',bg='black',fg='white').grid(row=24, column=0)
			self.detchoice=OptionMenu(self.root,self.detvar,'Circular','Rectangular',command=self.updateGUIdetector).grid(row=24, column=1, sticky='ew')
			Label(self.root,text='-',bg='black',fg='black').grid(row=24, column=2, sticky='ew')				
			self.det1label=Label(self.root,text='Width (cm)',bg='black',fg='white').grid(row=24, column=2)
			self.det1value=Entry(self.root,textvariable=self.varlist[66],bd=0).grid(row=24, column=3)
			self.det2label=Label(self.root,text='Height (cm)',bg='black',fg='white').grid(row=24, column=4)
			self.det2value=Entry(self.root,textvariable=self.varlist[67],bd=0).grid(row=24, column=5)	
		if dettype == 'Circular':
			self.varlist[65].set(0)
			self.detlabel=Label(self.root,text='Detector shape',bg='black',fg='white').grid(row=24, column=0)
			self.detchoice=OptionMenu(self.root,self.detvar,'Circular','Rectangular',command=self.updateGUIdetector).grid(row=24, column=1, sticky='ew')
			self.det1label=Label(self.root,text='Diameter (cm)',bg='black',fg='white').grid(row=24, column=2)
			self.det1value=Entry(self.root,textvariable=self.varlist[66],bd=0).grid(row=24, column=3)
			self.det2label=Label(self.root,text='-',bg='black',fg='black').grid(row=24, column=4, sticky='ew')	
			self.det2value=Label(self.root,text='-',bg='black',fg='black').grid(row=24, column=5, sticky='ew')

		self.root.update()	

	def updateCurvatureMatrixBasis(self,cumbasis):
	
		if cumbasis == '[a*,b*,c*]':
			self.varlist[86].set(0)
			self.cumbasislabel=Label(self.root,text='Basis choice',bg='black',fg='white').grid(row=17, column=6)
			self.cumbasischoice=OptionMenu(self.root,self.cumbasis,'[a*,b*,c*]','[Qx,Qy,Qz]',command=self.updateCurvatureMatrixBasis).grid(row=17, column=7, sticky='ew')
			self.Hxxlabel=Label(self.root,text=u'Hxx (meV.\u212b^2)',bg='black',fg='white').grid(row=4+13, column=0)
			self.Hxxvalue=Entry(self.root,textvariable=self.varlist[42],bd=0).grid(row=4+13, column=1)
			self.Hxylabel=Label(self.root,text=u'Hxy (meV.\u212b^2)',bg='black',fg='white').grid(row=4+13, column=2)
			self.Hxyvalue=Entry(self.root,textvariable=self.varlist[43],bd=0).grid(row=4+13, column=3)
			self.Hxzlabel=Label(self.root,text=u'Hxz (meV.\u212b^2)',bg='black',fg='white').grid(row=4+13, column=4)
			self.Hxzvalue=Entry(self.root,textvariable=self.varlist[44],bd=0).grid(row=4+13, column=5)
			self.Hyxlabel=Label(self.root,text=u'Hyx (meV.\u212b^2)',bg='black',fg='white').grid(row=4+14, column=0)
			self.Hyxvalue=Entry(self.root,textvariable=self.varlist[45],bd=0).grid(row=4+14, column=1)
			self.Hyylabel=Label(self.root,text=u'Hyy (meV.\u212b^2)',bg='black',fg='white').grid(row=4+14, column=2)
			self.Hyyvalue=Entry(self.root,textvariable=self.varlist[46],bd=0).grid(row=4+14, column=3)
			self.Hyzlabel=Label(self.root,text=u'Hyz (meV.\u212b^2)',bg='black',fg='white').grid(row=4+14, column=4)
			self.Hyzvalue=Entry(self.root,textvariable=self.varlist[47],bd=0).grid(row=4+14, column=5)
			self.Hzxlabel=Label(self.root,text=u'Hzx (meV.\u212b^2)',bg='black',fg='white').grid(row=4+15, column=0)
			self.Hzxvalue=Entry(self.root,textvariable=self.varlist[48],bd=0).grid(row=4+15, column=1)
			self.Hzylabel=Label(self.root,text=u'Hzy (meV.\u212b^2)',bg='black',fg='white').grid(row=4+15, column=2)
			self.Hzyvalue=Entry(self.root,textvariable=self.varlist[49],bd=0).grid(row=4+15, column=3)
			self.Hzzlabel=Label(self.root,text=u'Hzz (meV.\u212b^2)',bg='black',fg='white').grid(row=4+15, column=4)
			self.Hzzvalue=Entry(self.root,textvariable=self.varlist[50],bd=0).grid(row=4+15, column=5)
		if cumbasis == '[Qx,Qy,Qz]':
			self.varlist[86].set(1)
			self.cumbasislabel=Label(self.root,text='Basis choice',bg='black',fg='white').grid(row=17, column=6)
			self.cumbasischoice=OptionMenu(self.root,self.cumbasis,'[a*,b*,c*]','[Qx,Qy,Qz]',command=self.updateCurvatureMatrixBasis).grid(row=17, column=7, sticky='ew')
			self.Hxxlabel=Label(self.root,text=u'Hxx (meV/rlu^2)',bg='black',fg='white').grid(row=4+13, column=0)
			self.Hxxvalue=Entry(self.root,textvariable=self.varlist[42],bd=0).grid(row=4+13, column=1)
			self.Hxylabel=Label(self.root,text=u'Hxy (meV/rlu^2)',bg='black',fg='white').grid(row=4+13, column=2)
			self.Hxyvalue=Entry(self.root,textvariable=self.varlist[43],bd=0).grid(row=4+13, column=3)
			self.Hxzlabel=Label(self.root,text=u'Hxz (meV/rlu^2)',bg='black',fg='white').grid(row=4+13, column=4)
			self.Hxzvalue=Entry(self.root,textvariable=self.varlist[44],bd=0).grid(row=4+13, column=5)
			self.Hyxlabel=Label(self.root,text=u'Hyx (meV/rlu^2)',bg='black',fg='white').grid(row=4+14, column=0)
			self.Hyxvalue=Entry(self.root,textvariable=self.varlist[45],bd=0).grid(row=4+14, column=1)
			self.Hyylabel=Label(self.root,text=u'Hyy (meV/rlu^2)',bg='black',fg='white').grid(row=4+14, column=2)
			self.Hyyvalue=Entry(self.root,textvariable=self.varlist[46],bd=0).grid(row=4+14, column=3)
			self.Hyzlabel=Label(self.root,text=u'Hyz (meV/rlu^2)',bg='black',fg='white').grid(row=4+14, column=4)
			self.Hyzvalue=Entry(self.root,textvariable=self.varlist[47],bd=0).grid(row=4+14, column=5)
			self.Hzxlabel=Label(self.root,text=u'Hzx (meV/rlu^2)',bg='black',fg='white').grid(row=4+15, column=0)
			self.Hzxvalue=Entry(self.root,textvariable=self.varlist[48],bd=0).grid(row=4+15, column=1)
			self.Hzylabel=Label(self.root,text=u'Hzy (meV/rlu^2)',bg='black',fg='white').grid(row=4+15, column=2)
			self.Hzyvalue=Entry(self.root,textvariable=self.varlist[49],bd=0).grid(row=4+15, column=3)
			self.Hzzlabel=Label(self.root,text=u'Hzz (meV/rlu^2)',bg='black',fg='white').grid(row=4+15, column=4)
			self.Hzzvalue=Entry(self.root,textvariable=self.varlist[50],bd=0).grid(row=4+15, column=5)

		self.root.update()			
		
	def updateNRSEsettings(self):
		# Re-calculates NRSE settings when some parameters are changed
	
		parlist=[]
		for i in xrange(self.varlist.__len__()):
			parlist.append(float(self.varlist[i].get()))
		NRSEcalc=SERESCAL(parlist,self.methodvar.get())
		A1,A2,A3,A4,A5,A6,theta1,theta2,ni,nf,fratio,tmin,tmax,phi_S,sigma_S=NRSEcalc.calcAngles()
		self.setlist[0].set(int(1e2*rad2deg(theta1))/1.e2)
		self.setlist[1].set(int(1e2*rad2deg(theta2))/1.e2)
		self.setlist[2].set(int(1e3*fratio)/1.e3)
		self.setlist[3].set(int(1e2*tmin)/1.e2)
		self.setlist[4].set(int(1e2*tmax)/1.e2)
	
	def importParameters(self):
		# Sets GUI values with values read in a parameter file	
		
		fid=askopenfile(mode='r',filetypes=[('pySERESCAL parameters files', '*.par'),('pySERESCAL ini files', '*.ini')])
		if fid is None:
			return
		params=fid.readlines()
		fid.close()
		params=[float(i) for i in params]
		for i in xrange(86):
			self.varlist[i].set(params[i])
		
	def saveParameters(self):
		# Saves GUI values to a file
	
		fid=asksaveasfile(mode='w', filetypes=[('pySERESCAL parameters files', '*.par'),('pySERESCAL ini files', '*.ini')])
		if fid is None:
			return
		for i in xrange(86):
			fid.write(str(self.varlist[i].get())+'\n')
		fid.close()
		
	def plotNRSEresolution(self):
		# Calculates and plots resolution curves in a separate window
		
		self.updateNRSEsettings()
		parlist=[]
		for i in xrange(self.varlist.__len__()):
			parlist.append(float(self.varlist[i].get()))
		NRSEcalc=SERESCAL(parlist,self.methodvar.get())
		NRSEcalc.updateParams(parlist,self.methodvar.get()) 		

		A1,A2,A3,A4,A5,A6,theta1,theta2,ni,nf,fratio,tmin,tmax,phi_S,sigma_S=NRSEcalc.calcAngles()
			
		tau=linspace(0,tmax,1e2+1)		
		
		if self.methodvar.get()=='cn':
			methodstr='Cooper-Nathans'			
		elif self.methodvar.get()=='pop':
			methodstr='Popovici'
		else:
			print '\nError: unknown calculation method!!!'
			return
		
		pinstr,pcurv,pimperf,ptot=NRSEcalc.computeResolutionCurves(tau,self.methodvar.get())
		
		if self.varlist[82].get() == 0:
			extypestr='Magnon'
			pmag_par,pmag_antipar=NRSEcalc.computeMagneticFactor(tau)
			pmag_par*=ptot
			pmag_antipar*=ptot
		else:
			extypestr='Phonon'
		
		print '\nPlotting NRSE resolution curves (',methodstr,') ...\n'
		
		# Rescale Fourier times for further plotting
		tmin*=1e12
		tmax*=1e12
		tau*=1e12
		
		fig = plt.gcf()
		fig.canvas.manager.window.raise_()
		fig.canvas.set_window_title('NRSE resolution curves')
		plt.clf()
		fig.patch.set_facecolor('k')
		plt.rc_context({'axes.edgecolor':'k', 'xtick.color':'w', 'ytick.color':'w', 'figure.facecolor':'w'})
		plt.ion()
		ax=plt.subplot(111)
		plt.plot([tmin,tmin],[-1.1,1.1],'-',color='0.5',linewidth=1.5)
		plt.plot([0,tmax],[1.,1.],'-',color='0.5',linewidth=1.5)
		plt.plot([0,tmax],[.33333,.33333],'--',color='0.5',linewidth=1.5)
		plt.plot([0,tmax],[0,0],'-',color='0.5',linewidth=1.5)
		plt.plot([0,tmax],[-1.,-1.],'-',color='0.5',linewidth=1.5)
		plt.plot([0,tmax],[-.33333,-.33333],'--',color='0.5',linewidth=1.5)
		p0,=plt.plot(tau,pinstr,'r-',linewidth=2.5)
		p1,=plt.plot(tau,pcurv,'b-',linewidth=2.5)
		p2,=plt.plot(tau,pimperf,'g-',linewidth=2.5)
		if self.varlist[82].get() == 0:
			p4,=plt.plot(tau,pmag_par,'y-',linewidth=2.5)			
			p5,=plt.plot(tau,pmag_antipar,'m-',linewidth=2.5)
		else:
			p3,=plt.plot(tau,ptot,'k-',linewidth=2.5)
		plt.xlabel('Fourier time (ps)',fontsize=18,color='w')
		plt.ylabel('Polarisation',fontsize=18,color='w')
		#plt.title('NRSE resolution curves',fontsize=18,color='w')
		plt.xlim(xmin=0,xmax=tau.max())
		if self.varlist[82].get() == 0:
			plt.ylim(ymin=min(minimum(array(pmag_par),array(pmag_antipar)))-0.05,ymax=1.05)
		else:
			plt.ylim(ymin=0,ymax=1.05)		
		if self.varlist[82].get() == 1:		
			plt.legend([p0,p1,p2,p3],['Instrumental','Curvature','Sample imperfections','Total'], \
				loc='center left',bbox_to_anchor=(1,.875))
		else:
			plt.legend([p0,p1,p2,p4,p5],['Instrumental','Curvature','Sample imperfections','Total (++)','Total (+-)'], \
				loc='center left',bbox_to_anchor=(1,.815))		
		if self.varlist[82].get() == 0:
			ptotabs=abs(abs(pmag_par)-.33333)
			tlim_par=tau[int(argmin(ptotabs))]
			if abs(max(pmag_par)) < .33333:
				tlim_par=0
			elif tlim_par > tmax:
				tlim_par=tmax
			ptotabs=abs(abs(pmag_antipar)-.33333)
			tlim_antipar=tau[int(argmin(ptotabs))]
			if max(abs(pmag_antipar)) < .33333:
				tlim_antipar=0
			elif tlim_antipar > tmax:
				tlim_antipar=tmax
			plt.text(1.1*tmax,min(minimum(array(pmag_par),array(pmag_antipar)))+0.2,u'\n'+extypestr+' mode\n'+methodstr+' method\n'\
				'\n$t_{min}$ = '+str(round(tmin,2))+' ps'+\
				'\n$t_{1/3}^{++}$ = '+str(round(tlim_par,2))+' ps'+\
				'\n$t_{1/3}^{+-}$ = '+str(round(tlim_antipar,2))+' ps'+\
				'\n$t_{max}$ = '+str(round(tmax,2))+' ps'
				,fontsize=14,color='w')
		else:
			ptotabs=abs(ptot-.33333)
			tlim=tau[int(argmin(ptotabs))]
			if tlim > tmax:
				tlim=tmax
			plt.text(1.1*tmax,.2,u'\n'+extypestr+' mode\n'+methodstr+' method\n'\
				'\n$t_{min}$ = '+str(round(tmin,2))+' ps'+\
				'\n$t_{1/3}$ = '+str(round(tlim,2))+' ps'+\
				'\n$t_{max}$ = '+str(round(tmax,2))+' ps'
				,fontsize=14,color='w')			
		box=ax.get_position()
		ax.set_position([box.x0, box.y0, box.width*0.6, box.height])
		plt.show()	
		
	def saveNRSEresolution(self):
		# Saves resolution curves to an ASCII file
		
		self.updateNRSEsettings()
		parlist=[]
		for i in xrange(self.varlist.__len__()):
			parlist.append(float(self.varlist[i].get()))
		NRSEcalc=SERESCAL(parlist,self.methodvar.get())
		NRSEcalc.updateParams(parlist,self.methodvar.get()) 		
		A1,A2,A3,A4,A5,A6,theta1,theta2,ni,nf,fratio,tmin,tmax,phi_S,sigma_S=NRSEcalc.calcAngles()	
		
		tau=linspace(0,tmax,1e2+1)	
		
		if self.methodvar.get()=='cn':
			methodstr='Cooper-Nathans'		
		elif self.methodvar.get()=='pop':
			methodstr='Popovici'
		else:
			print '\nError: unknown calculation method!!!'
			return

		pinstr,pcurv,pimperf,ptot=NRSEcalc.computeResolutionCurves(tau,self.methodvar.get())
		
		if self.varlist[82].get() == 0:
			pmag_par,pmag_antipar=NRSEcalc.computeMagneticFactor(tau)
			pmag_par*=ptot
			pmag_antipar*=ptot
			
		fid=asksaveasfile(mode='w', defaultextension='.dat')
		if fid is None:
			return
		if self.varlist[82].get() == 0:			
			fid.write('tau\t P_instr\t P_curv\t P_imperf\t P_total\t P_mag_par\t P_mag_antipar\n')
			for i in xrange(tau.size):
				fid.write(str(tau[i]*1e12)+'\t'+str(pinstr[i])+'\t'+str(pcurv[i])+'\t'+str(pimperf[i])+'\t'+str(ptot[i])+'\t'+str(pmag_par[i])+'\t'+str(pmag_antipar[i])+'\n')
		else:
			fid.write('tau\t P_instr\t P_curv\t P_imperf\t P_total\n')
			for i in xrange(tau.size):
				fid.write(str(tau[i]*1e12)+'\t'+str(pinstr[i])+'\t'+str(pcurv[i])+'\t'+str(pimperf[i])+'\t'+str(ptot[i])+'\n')
		fid.close()

	def plotTASellipsoids(self):
		parlist=[]
		for i in xrange(self.varlist.__len__()):
			parlist.append(float(self.varlist[i].get()))
		TAScalc=RESCAL(parlist,self.methodvar.get())
		
		if self.methodvar.get()=='cn':
			
			methodstr='Cooper-Nathans'
			print '\nPlotting TAS resolution ellipsoids (',methodstr,') ...\n'
			R0,MCN=TAScalc.computeCNmatrix() # get Cooper-Nathans matrix
			th=linspace(0,2*pi,1e4+1)
			HKLtoRec=TAScalc.getReciprocalBasis()
			
			### (Qx,E)-plane
			if norm(TAScalc.qvec) == 0:
				rlutoinvA=norm(dot(HKLtoRec,TAScalc.Qvec))/norm(TAScalc.Qvec)				
			else:
				rlutoinvA=norm(dot(HKLtoRec,TAScalc.qvec))/norm(TAScalc.qvec)					
			# 50 % ellipse coordinates
			r=sqrt(2*log(2)/(MCN[2,2]*(sin(th))**2+MCN[0,0]*(cos(th))**2+2*MCN[0,2]*sin(th)*cos(th))) # 50% ellipsoid		
			xel50=r*cos(th)/rlutoinvA+norm(TAScalc.qvec)*sign(dot(TAScalc.qvec,TAScalc.Avec))			
			yel50=r*sin(th)+TAScalc.En
			# 30 % ellipse coordinates
			r=sqrt(2*log(3)/(MCN[2,2]*(sin(th))**2+MCN[0,0]*(cos(th))**2+2*MCN[0,2]*sin(th)*cos(th))) # 33% ellipsoid		
			xel30=r*cos(th)/rlutoinvA+norm(TAScalc.qvec)*sign(dot(TAScalc.qvec,TAScalc.Avec))			
			yel30=r*sin(th)+TAScalc.En			
			# 20 % ellipse coordinates (to get comfortable space between ellipse and figure borders)
			r=sqrt(2*log(5)/(MCN[2,2]*(sin(th))**2+MCN[0,0]*(cos(th))**2+2*MCN[0,2]*sin(th)*cos(th))) # 20% ellipsoid
			xel20=r*cos(th)/rlutoinvA+norm(TAScalc.qvec)*sign(dot(TAScalc.qvec,TAScalc.Avec))
			yel20=r*sin(th)+TAScalc.En		
			xmin=min(xel20)
			xmax=max(xel20)
			ymin=min(yel20)
			ymax=max(yel20)			
			# Local dispersion curve
			if norm(TAScalc.qvec) == 0:	# Zone center	
				xcen=0
			else:
				xcen=norm(TAScalc.qvec)*sign(dot(TAScalc.qvec,TAScalc.Avec))
			xvel=linspace(xmin,xmax,101)
			#yvel=abs(xvel)*TAScalc.dEdq/rlutoinvA+TAScalc.Egap							
			yvel=sign(xcen)*(xvel-xcen)*TAScalc.dEdq/rlutoinvA+TAScalc.En-TAScalc.dE
			# (Qx,E)-widths
			qxwidth=float(sqrt(8*log(2))/sqrt(MCN[0,0])/rlutoinvA)
			qywidth=float(sqrt(8*log(2))/sqrt(MCN[1,1])/rlutoinvA)	
			ewidth=float(sqrt(8*log(2))/sqrt(MCN[2,2]))
			#
			xwdthline=[TAScalc.En,TAScalc.En]
			xwdthlinex=[xcen-0.5*qxwidth,xcen+0.5*qxwidth]
			ywdthlinex=[xcen,xcen]
			ywdthline=[TAScalc.En-0.5*ewidth,TAScalc.En+0.5*ewidth]			
			# Draw all
			fig = plt.gcf()
			fig.canvas.manager.window.raise_()
			fig.canvas.set_window_title('TAS')
			fig.patch.set_facecolor('k')
			plt.rc_context({'axes.edgecolor':'k', 'xtick.color':'w', 'ytick.color':'w', 'figure.facecolor':'w'})
			plt.ion()
			plt.clf()
			plt.plot(xvel,yvel,'b-',linewidth=20, alpha=0.3)
			plt.plot(xwdthlinex,xwdthline,'r',linewidth=1.5)	
			plt.plot(ywdthlinex,ywdthline,'r',linewidth=1.5)
			plt.plot(xel30,yel30,'r',linewidth=1.5)
			plt.plot(xel50,yel50,'r',linewidth=3.0)			
			plt.xlim(xmin=xmin,xmax=xmax)
			plt.ylim(ymin=ymin,ymax=ymax)
			plt.xlabel('q (r.l.u.)',fontsize=18,color='w')
			plt.ylabel('Energy (meV)',fontsize=18,color='w')
			plt.title('TAS resolution ellipsoid',fontsize=18,color='w')
			plt.text(xmin+0.05*abs(xmax-xmin),ymin+0.05*abs(ymax-ymin), 'Prefactor = %.3f\n$\delta$q = %.3f r.l.u.; $\delta$E = %.3f meV' % (R0,qxwidth,ewidth), fontsize=18)
			plt.grid(True)
			plt.show()
			# Print values 
			print 'Bragg widths:\n-------------'
			print 'Qx: %.3f r.l.u. (FWHM)' % qxwidth
			print 'Qy: %.3f r.l.u. (FWHM)' % qywidth		
			print 'E: %.3f meV (FWHM)'% ewidth
		else:
			print '\nError: unknown calculation method!!!'
			return
	
	def quitAll(self):
		# Kill all existing windows and close app
		
		if plt.get_fignums().__len__() > 0:
			plt.close('all')
			print 'Killing remaining windows ...'
		self.root.destroy()
		t=time.strftime('%d/%m/%y %H:%M:%S',time.localtime()).split()
		print '\npySERSCAL ended on',t[0],'at',t[1],'\n'
		
if __name__ == "__main__":
	
	app=pySERESCAL()