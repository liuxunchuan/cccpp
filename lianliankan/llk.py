#!coding=utf8
import logging
log = logging.getLogger()
log.setLevel(logging.INFO)
import random
import matplotlib.pyplot as plt
from matplotlib import animation 
import matplotlib.patches as patches 
import time
import threading
import thread

class Signal:
   typeCode = 0
   typeName = '0'
   def __init__(self):
      pass
  

class AddOn(object):
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

   def isFireable(self):
      return False
   
   def getEventSeries(self, event):
      pass
   

class Animal(AddOn):
   typeCode = 1
   typeName = 'animal'
   namecode = 0 #(1,2,3,4,5)
   name = '0'
   picture = None
   colormap = ['r','g','b','black','orange']

   def __init__(self,namecode,NEventBeforeToBeFired=0):
      super(Animal,self).__init__()
      self.loc = None
      self.namecode = namecode
      self.NEventBeforeToBeFired = NEventBeforeToBeFired
      self.name = chr(64+namecode) # 1 == A

      self.color = self.__class__.colormap[namecode-1]
      self.colorcode = namecode
      self._isAlive = True
      self._isFireable = True
      
      pass

   def isAlive(self):
      return self._isAlive

   def isFireable(self):
      return self._isFireable
   
   def setLoc(self,loc):
      self.loc = loc   

   def getString(self):
      return self.name 

   def getColor(self):
      return self.color

   def getEventSeries(self,event):
      
      if(event.name=='fire'):
         if(self.NEventBeforeToBeFired>0):
            self.NEventBeforeToBeFired = self.NEventBeforeToBeFired-1
            if(self.NEventBeforeToBeFired==0):
               e = event.__class__(event.mapTable,self,self.loc,1,view=event.view)
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

  def isWithFireableAddOns(self,loc):
    ceil = self.getNoneEmptyCeil(loc)
    if (not ceil) or (not ceil.addOns):
       return False
    return ceil.addOns.isFireable()
    

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

  def getExplosionFromFocusPoints(self,focusPointLocations=[]):
     """
     [ [loc, [LN,RN,TN,BN],direction,type ], ...  ] 
     direction:
        0: default except for 4-series
        1: hori
        2: vert
     type:
        1: normal
        2: 4-series
        3: bomb-series
        4: bird-series
     """
     result = [ ]
     

     for (i,j) in focusPointLocations:
        
        if not self.isWithFireableAddOns((i,j)):
           continue

        colorcode = self[(i,j)].addOns.colorcode

        LN = 0
        while(True):
           if not self.isWithFireableAddOns((i,j-LN-1)):
              break
           if(self[(i,j-LN-1)].addOns.colorcode != colorcode):
              break
           LN = LN+1

        RN = 0
        while(True):
           if not self.isWithFireableAddOns((i,j+RN+1)):
              break
           if(self[(i,j+RN+1)].addOns.colorcode != colorcode):
              break
           RN = RN+1

        TN = 0
        while(True):
           if not self.isWithFireableAddOns((i-TN-1,j)):
              break
           if(self[(i-TN-1,j)].addOns.colorcode != colorcode):
              break
           TN = TN+1

        BN = 0
        while(True):
           if not self.isWithFireableAddOns((i+BN+1,j)):
              break
           if(self[(i+BN+1,j)].addOns.colorcode != colorcode):
              break
           BN = BN+1

        if(LN+RN<2):
           LN = 0
           BN = 0
 
        if(TN+BN<2):
           TN = 0
           BN = 0

        if( LN+RN<2 and BN+TN<2):
           break

        direction = 0
        typeSeries = 1
        if(LN+RN>=4 or TN+BN>=4):
           typeSeries = 4

        elif(LN+RN>=2 and BN+TN>=2):
           typeSeries = 3
        else:
           if(LN+RN==3):
              direction = 1
              typeSeries = 2
           elif(TN+BN==3):
              direction = 2
              typeSeries = 2

        result.append([(i,j), [LN,RN,TN,BN],direction,typeSeries])

     return result

  def compareSeries(self,a,b):
     """
     return (r1, r2)
     r1:
       0:  series a and b do not interface
       1:  otherwise
     r2:
       -1: a<b when r1==1
        1: a>b when r1==1
        0: a==b when r1==1, or r1==0
        2: just interface
     """
     return (0,0)

     loc1 = a[0]
     loc2 = b[0]

     r1 = 0
     r2 = 0
     if (loc1[0] != loc2[0]) and (loc1[1]!=loc2[1]):
        return (0,0)
     if (loc1[0] == loc2[0]) and (loc1[1]==loc2[1]):
        return (1,0)
     if(loc1[0] == loc2[0]):
        pass            

class Viewer:
   def __init__(self, mapTable):
      plt.ion()
      self.ax = plt.axes()
      self.mapTable = mapTable
      self.do = True
      self.stop = False
      self.flashes = [] #[[flash, (oldpatch1, oldpatch2)]]
      self.closeflash = False
      self.viewEventPool = []
      self.flush = False
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
      self.ax.patches = []
      for i in xrange(0, self.mapTable.NRows):
         for j in xrange(0, self.mapTable.NCols):
            loc = (j+0.5, self.mapTable.NRows-0.5-i)
            addOns = self.mapTable[i,j].getAddOns()
            if(addOns != None):
               patch = patches.Rectangle(
                                (j+0.1,self.mapTable.NRows-i-0.9),   # (x,y)
                                0.8,          # width
                                0.8,          # height
                              )
               patch.set_color(addOns.getColor())
               self.ax.add_patch(patch)
               self.mapTable[(i,j)].addOnsImage = patch

      self.ax.figure.canvas.draw()

   def closeFlash(self):
      self.closeflash = True

   def openFlash(self):
      self.closeflash = False

   def run(self):
      self.flashMovie()

   def flashMovie(self):
    while(True):
      if(self.flush):
         for event in self.viewEventPool:
           if event.name == 'fire':
             loc =  event.getLocation()
             if hasattr(self.mapTable[loc],'addOnsImage'):
               color = self.mapTable[loc].addOnsImage.get_facecolor()
               self.mapTable[loc].addOnsImage.set_visible(False)
               self.flashes.append(
                  [FireRectFlash( (loc[1]+0.5, self.mapTable.NRows-0.5-loc[0]),color=color),()]
               )
         self.viewEventPool = []
         self.flush = False


      if(self.closeflash): return None
      for [flash, oldPatches] in self.flashes:
         dex = self.flashes.index( [flash, oldPatches] )

         for item in oldPatches:
            try:
               idex = self.ax.patches.index(item)
               self.ax.patches.pop(idex)
            except Exception,e:
               pass

         newPatches = []
         try:
           newPatches = flash.next()
           self.flashes[dex][1] = newPatches
         except Exception,e:
           self.flashes.pop(dex)
           
           continue
           #   raise e 

         for item in newPatches:
            self.ax.add_patch(item)

      plt.draw()
      plt.pause(0.0001)
      time.sleep(0.05)

   def handle(self,ev):
     if(ev.name == 'flush'):
       self.flush = True
     else:
       self.viewEventPool.append(ev) 

class Flash(object):
    def __init__(self,*args):
       pass
class FireRectFlash(Flash):
    def __init__(self,loc,length=0.8,fullLength=1,N=5,color='blue',*args):
       super(FireRectFlash,self).__init__(*args)
       self.loc=loc
       self.length = length
       self.fullLength = fullLength
       self.totalN = N
       self.N = N
       self.color = color

    def next(self):
       if(self.N < 1):
          #self.N = self.N - 1
          raise Exception('fireRectEventError')
       Dx = self.length/2.*self.N/self.totalN
       dx = (self.fullLength-self.length)/2.*(self.totalN-self.N)/self.totalN
       (x,y) = self.loc
       self.N = self.N - 1
       out = [
               patches.Rectangle((x-self.length/2.-dx, y-self.length/2.-dx), dx+Dx, dx+Dx  ),
               patches.Rectangle((x+self.length/2.-Dx, y-self.length/2.-dx), dx+Dx, dx+Dx  ),
               patches.Rectangle((x-self.length/2.-dx, y+self.length/2.-Dx), dx+Dx, dx+Dx  ),
               patches.Rectangle((x+self.length/2.-Dx, y+self.length/2.-Dx), dx+Dx, dx+Dx  ),
             ]
       
       for item in out:
          item.set_color(self.color)
       return out

class Event(object):
   def __init__(self,name,mapTable,addOns=None,loc=None,
                toBeActedTime=0,stepsEvent=False,view=None, **kwargs):
      self.name = name
      self.mapTable = mapTable
      self.addOns = addOns # 
      self.loc = loc
      self.toBeActedTime = toBeActedTime
      self.stepsEvent = stepsEvent
      self.view = view
      for key in kwargs:
         self.__dict__[key] = kwargs[key]

   def setAttr(self,key,value):
      self.__dict__[key] = value

   def getLocation(self):
      return self.loc #if not self.loc else self.addOns.loc

   def getLocationsAffected(self):
      pass

   def getEventSeries(self):
      pass

   def operateEvent(self):
      if(self.addOns and self in self.addOns.events):
         self.addOns.events.pop(self.addOns.events.index(self))

class FireEvent(Event):
   def __init__(self,mapTable,addOns=None,loc=None,toBeActedTime=1,
                stepsEvent = False,view=None,**kwargs):
      super(FireEvent,self).__init__(name='fire',mapTable=mapTable,addOns=addOns,
                                     loc=loc,toBeActedTime=toBeActedTime,
                                     stepsEvent=stepsEvent,view=view,
                                     **kwargs)

   def getLocationsAffected(self):
      pass

   def getEventSeries(self):
      #print self.getLocation()
      (i,j) = self.getLocation()
      #print self.getLocation()
      out = []
      for d in xrange(4):
         ceil = self.mapTable[i+((d+1)%2)*(d-1),j+(d%2)*(d-2)]
         if(ceil):
            addOns = ceil.getAddOns()
            if addOns:
               out.extend(addOns.getEventSeries(self))
            
      return out

   def operateEvent(self):
      loc = self.getLocation()
      self.mapTable[loc].addOns = None
      super(FireEvent,self).operateEvent()
      self.view.handle(self)

class EventsPool:
   def __init__(self,mapTable=None,view=None):
      self.mapTable = mapTable
      self.view = view
      self.eventPool = []
      self.stepsEventPool = []
      self.flushEvent = Event(name='flush',mapTable=mapTable)

   def isCleanEventPool(self):
      return len(self.eventPool) == 0

   def canBeActedNow(self,e):
      # this function will not check e.toBeActedTime<=0 or not
      # a lot to be written ...
      return e.toBeActedTime <= 0
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
      if e.stepsEvent:
         self.stepsEventPool.append(e)
         return None
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

   def fetchStepsEvent(self):
      if self.stepsEventPool in [ [], None]:
         return None
      else:
         return self.stepsEventPool.pop(0)

   def run(eventPool):
    while(True):
      while(not eventPool.isCleanEventPool()):
   #for dd in xrange(3):
         
         #print('step')
         while(True):
            event =  eventPool.fetchOneToBeActed()
            if(event == None):
               break
            es = event.getEventSeries()
            #print( 'es:  ', len(es))
            if(es):
               for e in es:
                   eventPool.insertEvent(e)
            event.operateEvent() 
         
         #print('hh')
         eventPool.view.handle(eventPool.flushEvent)
         #print('ha')
         time.sleep(1)
         

         eventPool.timeSubmitOne()
      sevent = eventPool.fetchStepsEvent()
      if sevent == None:
         break

      sevent.toBeActedTime = 0
      sevent.stepsEvent = False
      eventPool.insertEvent(sevent)

    eventPool.view.closeFlash()
    
def starGame():
   NRows = 9
   NCols = 9
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
         animal.setLoc((i,j))

  
   fireEvent = FireEvent(mapTable)

   loc = (4,4)
   addOns = mapTable[loc].addOns   

   fireEvent = FireEvent(mapTable,addOns=addOns,loc=loc,toBeActedTime=0,stepsEvent=True)
   addOns.NEventBeforeToBeFired = 0
   addOns.events.append(fireEvent)
   
   plt.ion()
   plt.figure(1)
   plt.clf()  
   view = Viewer(mapTable)
   
   fireEvent.view = view
 
   eventPool = EventsPool(mapTable,view)
   eventPool.insertEvent(fireEvent)

   #for i in xrange(NRows):
   #   for j in xrange(NCols):
   #      mapTable.arrayCallForFilledCeils.append((i,j))
   #      if(i==0):
   #         mapTable.arrayCanBeFilledWithnOneStepCeils.append((i,j))
   #         mapTable.arrayGeneratorCeils.append((i,j))




   view.showMap()
   view.showAddOns()
   plt.draw()
   plt.pause(0.0001)
   time.sleep(1)

   eventThread = threading.Thread(target=eventPool.run)
   eventThread.setDaemon(True)
   eventThread.start() 


   view.run()
   
   eventThread.join()

   return None

def run():
   while(True):
      starGame()


