#!coding=utf8
import logging
log = logging.getLogger()
log.setLevel(logging.INFO)
import random
import matplotlib.pyplot as plt

class Signal:
   typeCode = 0
   typeName = '0'
   def __init__(self):
      pass
  

class AddOn:
   typeCode = 0
   typename = '0'
   def __init__(self):
      self.events = []
      self.loc = None
      pass
   def isMoveableWhenRedistributing(self):
      return True
   def isMoveable(self):
      return True
   def canOverlap(self):
      return False

  
   
   def getEventSeries(self, event):
      pass
   

class Animal(AddOn):
   typeCode = 1
   typeName = 'animal'
   namecode = 0 #(1,2,3,4,5)
   name = '0'
   picture = None

   def __init__(self,namecode,NEventBeforeToBeFired=1):
      self.loc = None
      self.namecode = namecode
      self.NEventToBeFired = NEventBeforeToBeFired
      self.name = chr(64+namecode) # 1 == A
      pass
   
   def setLoc(self,loc):
      self.loc = loc   

   def getString(self):
      return self.name 

   def getEventSeries(self,event):
      if(event.name=='fire'):
         if(self.NEventBeforeToBeFired>0):
            self.NEventBeforeToBeFired = self.NEventBeforeToBeFired-1
            if(self.NEventBeforeToBeFired==0):
               e = event.__class__(event.name,event.mapTable,self)
               self.events.append(e)
               return [e,]
         return []
      else:
          # to be written
          return []

class Generator:
   def __init__(self,seed,N):
      self.random = random.Random()
      self.random.seed(seed)
      self.N = N
   def randint(self):
      return self.random.randint(1,self.N)

class MapCeil:
  
  def __init__(self):
     self.loc = None
     self._empty =None
     self._hasRedirectionTopCeil=None
     self._redirectionTopCeil=None
     self._hasRedirectionBottomCeil=None
     self._redirectionBottomCeil=None

     self.top=None
     self.bottom=None
     self.left=None
     self.right=None
     self.topleft=None
     self.topright=None
     self.bottomleft=None
     self.bottomright=None

     self.topRedirected=None
     self.bottomRedirected=None

     self.generator=None

     self.addOns = None

     self._empty = False
     pass
  def isEmpty(self):
     return self._empty
  def setEmpty(self,TorF):
     self._empty = TorF

  def redirect(self,direction,ceil):
     if(direction == 'top'):
        self.topRedirected = True
     if(direction == 'bottom'):
        self.bottomRedirected = True
     
     setattr(self,direction,ceil)

  def setGenerator(self,generator):
     self.generator = generator

  def getGenerator(self):
     return self.generator

  def setAddOns(self,items):
     self.addOns = items
  
  def getAddOns(self):
     return self.addOns
     

class MapTable:
  def __init__(self,NRows,NCols,Ceil):
    self.NRows = NRows
    self.NCols = NCols
    self.mapTable = []

    for i in xrange(NRows):
       temp = []
       for j in xrange(NCols):
          if (Ceil is MapCeil) or (Ceil in MapCeil.__bases__):
             temp.append(Ceil())
             temp[-1].loc = (i,j)
          else:
             temp.append(None)
       self.mapTable.append(temp)

    self.arrayCallForFilledCeils = []
    self.arrayCanBeFilledWithnOneStepCeils = []
    self.arrayGeneratorCeils = []
      
    pass

  def __getitem__(self,loc):
    (i,j) = loc
    if(i<0 or j<0 or i>=self.NRows or j>= self.NCols):
       return None
    return self.mapTable[i][j]

  def getNoneEmptyCeil(self,loc):
    if self[loc] == None:
       return None
    if self[loc].isEmpty():
       return None
    return self[loc]

  def __setitem__(self,loc,value):
    (i,j) = loc
    self.mapTable[i][j] = value

  def setEmpty(self,location,TorF=True):
    if len(location) == 0:
       return
    if hasattr(location,'__len__'):
       for (i,j) in location:
          self[i,j].setEmpty(TorF)
    else:
       self[location].setEmpty(TorF)    

  def initConnection(self):
    for i in xrange(self.NRows):
       for j in xrange(self.NCols):
          if(self[i,j].isEmpty()):
             break
          self[i,j].top = self.getNoneEmptyCeil((i-1,j))
          self[i,j].bottom = self.getNoneEmptyCeil((i+1,j))
          self[i,j].left = self.getNoneEmptyCeil((i,j-1))
          self[i,j].right = self.getNoneEmptyCeil((i,j+1))
          self[i,j].topleft = self.getNoneEmptyCeil((i-1,j-1))
          self[i,j].topleft = self.getNoneEmptyCeil((i-1,j+1))
          self[i,j].bottomleft = self.getNoneEmptyCeil((i+1,j-1))
          self[i,j].bottomright = self.getNoneEmptyCeil((i+1,j+1))


  def redirect(self,loc1,loc2,direction1='top',direction2='bottom'):
     if (direction1, direction2) in [('top','bottom'),('bottom','top')]:
        if(direction1=='top'): self[loc1].setGenerator(None)
        if(direction2=='top'): self[loc2].setGenerator(None)
        self[loc1].redirect(direction1, self[loc2])
        self[loc2].redirect(direction1, self[loc1])
     else:
        log.info('No such redirection at now!')

  def setGenerator(self,loc,generator):
     self[loc].setGenerator(generator)

  def getGenerator(self,loc,generator):
     self[loc].getGenerator()

  def getNextMove():
     pass

class Viewer:
   def __init__(self, mapTable):
      plt.ion()
      self.ax = plt.axes()
      self.mapTable = mapTable
      pass

   def showMap(self):
      ax = self.ax
      ax.set_xlim(0,self.mapTable.NCols)
      ax.set_ylim(0,self.mapTable.NRows)
      ax.set_xticks([])
      ax.set_yticks([])
      for i in xrange(1,self.mapTable.NRows):
         ax.plot([0,self.mapTable.NCols],[i,i],color='black')
      for i in xrange(1,self.mapTable.NCols):
         ax.plot([i,i],[0,self.mapTable.NRows],color='black')
      

   def showAddOns(self):
      for i in xrange(0, self.mapTable.NRows):
         for j in xrange(0, self.mapTable.NCols):
            loc = (j+0.5, self.mapTable.NRows-0.5-i)
            addOns = self.mapTable[i,j].getAddOns()
            if(addOns != None):
               string = addOns.getString()
               plt.text(loc[0],loc[1],string)
       

class Event:
   def __init__(self,name,mapTable,addOns=None,loc=None,toBeActedTime=0):
      self.name = name
      self.mapTable = mapTable
      self.addOns = addOns # 
      self.loc = loc
      self.toBeActedTime = toBeActedTime

   def getLocation(self):
      return self.loc if not self.loc else self.addOns.loc

   def getLocationsAffected(self):
      pass

   def getEventSeries(self):
      pass

   def operateEvent(self):
      if(self.addOns and self in self.addOns.events):
         self.addOns.events.pop(self.addOns.events.index(self))

class fireEvent(Event):
   def __init__(self,mapTable,addOns=None,loc=None,toBeActedTime=1):
      super(CommonEvent,self).__init__('fire',mapTable,addOns)

   def getLocationsAffected(self):
      pass

   def getEventSeries(self):
      (i,j) = self.getLocation()
      for d in xrange(4):
         ceil = mapTable[i+((d+1)%2)*(d-1),j+(d%2)*(d-2)]
         if(ceil):
            addOns = ceil.getAddOns()
            return addOns.getEventSeries(self)
   def operateEvent(self):
      loc = self.getLocation()
      mapTable[loc].addOns = None
      super(fireEvent,self).operateEvent()

class EventsPool:
   def __init__(self):
      self.eventPool = []

   def isCleanEventPool(self):
      return len(self.eventPool) == 0

   def canbeActedNow(self,e):
      # this function will not check e.toBeActedTime<=0 or not
      # a lot to be written ...
      pass

   def fetchOneToBeActed(self):      
      for i in xrange(len(self.eventPool)):
         e = self.eventPool[i]
         if e.toBeActedTime > 0:
            return None
         if(self.canBeActedNow(e)):
            return self.eventPool.pop(i)
         
   
   def fetchToBeActeds(self):
      result = []
      for i in xrange(len(self.eventPool)):
         e = self.eventPool[i]
         if e.toBeActedTime > 0:
            break
         if(self.canBeActedNow(e)):
            result.append(self.eventPool.pop(i))
      return result

   def insertEvent(self,e):
      self.eventPool.insert(self.findDex(e.toBeActedTime),e)
      
   def findDex(self,t):
      dex = 0
      for  i in xrange(len(self.eventPool)):
         if self.eventPool[i].toBeActedTime <= t:
            dex=dex+1
         else:
            break
      return dex

   def timeSubmitOne(self):
      for e in self.eventPool:
         e.toBeActedTime = e.toBeActedTime - 1

   

class EventHandle:
   def __init__(self,eventPool):
      self.eventPool = eventPool
      pass

   def eventHandle(self,event):
      pass

   def getLocationsAffected(self,event):
      pass

   def getEventsTriggered(self,Event):
      pass

   def operateEvent(self,event):
      pass

   def doEventSeries(self):
      while(not self.isCleanEventPool()):
         while(True):
            event = self.eventPool.fetchOneToBeActed()
            if(event == None):
               break
            es = event.getEventSeries()
            if(es):
               for e in es:
                  self.eventPool.insertEvent(e)
            event.operateEvent()        
         self.eventPool.timeSubmitOne()
    
def starGame():
   NRows = 10
   NCols = 10
   seed = 1200
   NAnimals = 5
   mapTable = MapTable(NRows,NCols,MapCeil)
   mapTable.initConnection()
   generator = Generator(seed,NAnimals)
   for j in xrange(NRows):
      mapTable[0,j].setGenerator = generator
   for i in xrange(NRows):
      for j in xrange(NCols):
         animal = Animal(generator.randint())
         mapTable[i,j].setAddOns(animal)
         animal.setLoc = loc

   #for i in xrange(NRows):
   #   for j in xrange(NCols):
   #      mapTable.arrayCallForFilledCeils.append((i,j))
   #      if(i==0):
   #         mapTable.arrayCanBeFilledWithnOneStepCeils.append((i,j))
   #         mapTable.arrayGeneratorCeils.append((i,j))

   plt.clf()
   plt.ion()
   view = Viewer(mapTable)

   
   view.showMap()
   view.showAddOns()
   
     
  
  
