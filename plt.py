#! coding=utf8

import matplotlib.pyplot as plot
plot.ion()

import astropy.table as t
import astropy.units as u
year = 3.15576e+07

File = "Uout.csv"


class Printer:
   def __init__(self,File=File,plot=plot):
      self.T = t.Table.read(File)
      self.plot = plot
      self.colors = ['r','g','b','black','grey','yellow','purple','pink']
      self.linestyles=['-','--','-.']
       
   def draw(self,x='TIME',y=[],clear=True,Nbegin=20,ylim=[]):
      T = self.T[Nbegin:]
      if len(y)==0:
          raise Exception("empty Y");
      xdata = T[x]
      if x=='TIME':
         xdata = xdata/year
      cols  = T[y].columns
      if clear:
         plot.clf()
      for i in xrange(len(cols)):
         self.plot.loglog(xdata,cols[i],label=cols[i].name,\
                   color=self.colors[i%len(self.colors)],\
                   linestyle=self.linestyles[i/len(self.colors)])
      self.plot.legend()
      if ylim==[]:
         ylim=[1E-15,1E-4]
      self.plot.ylim(ylim)
