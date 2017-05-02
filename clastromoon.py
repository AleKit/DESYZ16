import time,datetime
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib
from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, get_moon, get_sun, Latitude
from astroplan.moon import moon_illumination
from scipy import interpolate
import pytz

#####falta por hacer:
## el tkinter o algo por el estilo para obtener muchos datos en
## un intervalo corto

######-------->>>>documentar

###AstroPy 1.2.1
###astroplan 0.2
###Python 2.7

#this class creates an object with the place of observation and you can introduce a time or an interval of time, 
class MoonPlan(object):

### location of the observer (some locations missing)
    def __init__(self,place='berlin'):

        if place == 'berlin':
            self.__place = EarthLocation(lat=52.5*u.deg, lon=13.4*u.deg, height=34*u.m)
        elif place != 'berlin':
            print('Not defined yet')
            raise SystemExit(0)

### setting time and interval: format 'yyyy-mm-dd hh:mm:ss' (it is UTC)
    def _set_start_time(self, time='now', deltamins = 0, deltadays = 3, points = 180): #npoints = 0):
        if time == 'now':
            self.__time = Time.now()
        else:
            self.__time = Time(time) #e.g. time = '2016-10-21 00:00:00' (UTC)
        self.__deltamins = deltamins
        self.__deltadays = deltadays
        self.__points = points
        #self.__npoints = npoints
        return self.__time #esto es necesario? es el primer elemento de datetimes

### creating the time interval using the data above
    def _set_interval_time(self):
        if self.__deltamins != 0 & self.__deltadays == 0:
            self.__intimes = np.linspace(0,self.__deltadays,self.__points)*u.min
        elif self.__deltadays != 0:
            self.__intimes = np.linspace(0,self.__deltadays,self.__points)*u.d
        else:
          #  print('Something is wrong with deltadays or deltamins')
          #  raise SystemExit(0)
            print('Be careful: deltamins = 0 = deltadays')
        self.__times_1 = self.__time + self.__intimes
        dictimes = {'interval': self.__intimes, 'datetimes': self.__times_1}
        return dictimes

### creating the frame to transform from ra, dec to alt, az; needs pressure, relative humidity, temperature and observed wavelength
    def _set_framing(self,pres = 1015*u.hPa, relhum = 0.7, temp = 15.0*u.deg_C, wavel = 0.5*u.micron):
        self.__frame_1 = AltAz(obstime = self.__times_1, 
                               location = self.__place, 
                               pressure = pres, 
                               relative_humidity = relhum, 
                               temperature = temp, obswl = wavel)


### gives altitude and azimuth of the sun
    def get_sun_altazs(self):
       # self.__set_start_time()##llamar a las funciones privadas aqui y en moon param y ponerles argumentos para que funcione (mirar abajo)
       # self.__set_interval_time()
       # self.__set_framing()
        self.sunaltazs_1 = get_sun(self.__times_1).transform_to(self.__frame_1)
        return self.sunaltazs_1

### gives altitude, azimuth and illumination of the moon
    def get_moon_param(self):
        #self.__set_start_time()
       # self.__set_interval_time()
       # self.__set_framing()
        moon = get_moon(self.__times_1)
        illum = moon_illumination(self.__times_1, self.__place)
        moonaltaz = moon.transform_to(self.__frame_1)

        moon_pars = {'moon illumination': illum, 'moon altitude': moonaltaz.alt, 'moon azimuth': moonaltaz.az} 

        return moon_pars

class MoonInterpolate(object):

    def __init__(self,astronight = -18*u.deg, minalt = 30, maxillum = 0.6): #, mod = 'both'):

        self.__astronight = astronight   #def config parameters: astronight, minalt, maxillum, error if they are not given
        self.__minalt = minalt
        self.__maxillum = maxillum
        self.__mod = None
        self.__moonPlan_1 = None
        self.__moonPlan_2 = None
        self.__moonGraph_1 = None
        self.__moonGraph_2 = None

   # def set_append(self, time, illum, altitude):
   #     deltaobj.append(time)
   #     lumiobj.append(illum)
   #     altiobj.append(altitude)

    def config_mod(self, mod):
        if mod == 'both' or mod == 'graph' or mod == 'text':
            self.__mod = mod
        else:
            print('Error: Wrong word. Only use graph, text or both')
            raise SystemExit(0)
          
    def config_moon_plan(self, place, time, deltamins, deltadays, points, npoints): 
        self.__moonPlan_1 = MoonPlan(place)
        self.__start = self.__moonPlan_1._set_start_time(time, deltamins, deltadays, points)
        print(self.__start)
        dictimes = self.__moonPlan_1._set_interval_time() 
        self.__intervaltime = dictimes['interval']
        self.__moonPlan_1._set_framing()
        self.__sunaltaz_1 = self.__moonPlan_1.get_sun_altazs()
        self.__moon_pars = self.__moonPlan_1.get_moon_param()
        
        self.__moonPlan_2 = MoonPlan(place)
        self.__moonPlan_2._set_start_time(time, deltamins, deltadays, npoints)
        self.__ndictimes = self.__moonPlan_2._set_interval_time()
      #  print(type(self.__ndictimes[interval]))
        self.__moonPlan_2._set_framing()
        self.__nsunaltaz_1 = self.__moonPlan_2.get_sun_altazs()
       # self.__nmoon_pars = self.__moonPlan_2.get_moon_param()

    def config_moon_plan_default(self):
        self.__moonPlan_1 = MoonPlan()

#doc              
    def get_interpolation(self): 

        if self.__moonPlan_1 is None:
            raise TypeError('moonPlan was not created')
        elif self.__mod is None:
            raise TypeError('mod was not chosen')
        else:
            self.__altitude = self.__moon_pars['moon altitude']
            self.__illum = self.__moon_pars['moon illumination']
          
          #  ntime = self.__ndictimes['interval']
            self.__ntime = self.__ndictimes['interval']
            
 

        interpalt = interpolate.interp1d(self.__intervaltime,self.__altitude,kind='cubic')
        interpillum = interpolate.interp1d(self.__intervaltime,self.__illum,kind='cubic')
        delta = []
        lumi = []
        moonaltis = []
        ndelta = []
        nlumi = []
        nmoonaltis = []
        self.__listdelta = []
        self.__listlum = []
        self.__listalt = []
        deltaobj = []
        lumiobj = []
        altiobj = []
        conig = 0
        for c in range(len(self.__ntime)):
#+0.5 little margin
            if self.__nsunaltaz_1.alt[c] <= self.__astronight+0.5*u.deg:
                co = 0
                if interpillum(self.__ntime[c]) < self.__maxillum:
                    ndelta.append(self.__ntime[c])
                    nlumi.append(interpillum(self.__ntime[c]))
                    nmoonaltis.append(interpalt(self.__ntime[c])*u.deg)
            elif co == 0:
                self.__listdelta.append(ndelta)
                self.__listlum.append(nlumi)
                self.__listalt.append(nmoonaltis)
                delta = []
                moonaltis = []
                lumi = []
                co = 1
            if (self.__mod == 'text' or self.__mod == 'both') and c < len(self.__ntime)-1:
     
                if c == 0:
                    dtxt = open ('data.txt', 'w')
                    dtxt.write('Nights we can observe the Moon shadow since ' + str(self.__start) + '(UTC) \n \n') 

                elif c % 500 == 0:
                    print("Step: "+ str(c))

## points with -18 deg (astronight) sun altitude               
                if  (interpillum(self.__ntime[c]) < self.__maxillum and 
                self.__nsunaltaz_1.alt[c] <= self.__astronight and 
                (self.__nsunaltaz_1.alt [c+1] > self.__astronight 
                 or self.__nsunaltaz_1.alt[c-1] > self.__astronight) and 
                interpalt(self.__ntime[c]) >= self.__minalt):
      
                    times_subset = self.__ntime[c-3 : c+4]
                    sunalta_subset = self.__nsunaltaz_1.alt[c-3 : c+4]
                    inter_time = interpolate.interp1d(sunalta_subset, times_subset, kind='cubic')

                    start_time = inter_time(self.__astronight) 
                    inter_moon = interpalt(start_time)
                    inter_ill = interpillum(start_time)
                    dtxt.write('Date and time: ' + str(self.__start+float(start_time)*u.d) + 
                               '\n' + 'Altitude: ' + str(round(inter_moon,3))+ '\n' + 
                               'Illumination: ' + str(round(inter_ill,3)) + '\n \n')    
               
                    if conig == 0:
                        conig = 1
                
## points with 30 deg (minalt) altitude in the night
                elif (interpillum(self.__ntime[c]) < self.__maxillum and 
                      self.__nsunaltaz_1.alt[c] <= self.__astronight and 
                      interpalt(self.__ntime[c]) >= self.__minalt and 
                      (interpalt(self.__ntime[c+1]) < self.__minalt or 
                       interpalt(self.__ntime[c-1]) < self.__minalt)):

                    times_subset = self.__ntime[c-3 : c+4]
            
                    val = []
            
                    for x in range(-3, 4):
                        val.append(interpalt(self.__ntime[c+x]))

                    moon_subset = np.array(val, dtype=float)
                
                    inter_time = interpolate.interp1d(moon_subset, times_subset, kind='cubic')

                    start_time = inter_time(self.__minalt) 
                    inter_moon = interpalt(start_time)
                    inter_ill = interpillum(start_time)
            
                    dtxt.write('Date and time: ' + str(self.__start+float(start_time)*u.d) + '\n' + 
                               'Altitude: ' + str(round(inter_moon,3))+ '\n' + 
                               'Illumination: ' + str(round(inter_ill,3)) + '\n \n') 

                    if conig == 0:
                        conig = 1

## points with 0.6(maxillum) illumination in the night
                elif (interpillum(self.__ntime[c]) < self.__maxillum and 
                      self.__nsunaltaz_1.alt[c] <= self.__astronight and 
                      interpalt(self.__ntime[c]) >= self.__minalt and 
                      (interpillum(self.__ntime[c+1]) > self.__maxillum or 
                       interpillum(self.__ntime[c-1]) > self.__maxillum)):  

                    times_subset = self.__ntime[c-3 : c+4]
                    
                    val = []
                    
                    for x in range(-3, 4):
                        val.append(interpillum(ntimes[c+x]))
                        
                        moon_subset = np.array(val, dtype=float)
                        
                        inter_time = interpolate.interp1d(moon_subset, times_subset, kind='cubic')
                        
                        start_time = inter_time(self.__maxillum)
                        inter_moon = interpalt(start_time)
                        inter_ill = interpillum(start_time)
                        

                        dtxt.write('Date and time: ' + str(self.__start+float(start_time)*u.d) + '\n' + 
                                   'Altitude: ' + str(round(inter_moon,3))+ '\n' + 
                                   'Illumination: ' + str(round(inter_ill,3)) + '\n \n')

                        if conig == 0:
                            conig = 1
        
        if self.__mod == 'text' or self.__mod == 'both':
            if conig == 0:
                dtxt.write('No nights available in this period: '  + str(self.__start) + ' - ' + str(self.__start+ntime[len(self.__ntime)-1]))
            dtxt.close()  


#collects all the useful points to scatter
    def get_scatter(self): #, time, illum, altitude, sunaltaz_1):
        self.__deltaobj = []
        self.__lumiobj = []
        self.__altiobj = []

        for c in range(len(self.__intervaltime)):
            
            if self.__sunaltaz_1.alt[c] <= self.__astronight:
                if self.__illum[c] < self.__maxillum:
                    self.__deltaobj.append(self.__intervaltime[c])
                    self.__lumiobj.append(self.__illum[c])
                    self.__altiobj.append(self.__altitude[c])
      

##this one changes the time list to a u.Quantity class, 
#the illumination list to an array and
#the altitude list to an Astropy Latitude object
#if there is nothing to change, it prints No data
    def change_object_type(self):
        if self.__deltaobj != []:
            self.__timesf = u.Quantity(self.__deltaobj) 
            self.__illumination = np.array(self.__lumiobj, dtype=float)
            self.__moonaltitude = Latitude(self.__altiobj)
        else:
            print('No data')

#this one makes the plots, both the scattering and the interpolation
#it uses the class  MoonGraph  (see below)        
    def config_moon_graph_1(self):
        if self.__mod == 'graph' or self.__mod == 'both':
            self.get_scatter()
            self.change_object_type()
            self.__moonGraph_1 = MoonGraph(self.__astronight, self.__minalt)
            self.__moonGraph_1.filling_graph(self.__ntime,self.__nsunaltaz_1)
            self.__moonGraph_2 = MoonGraph(self.__astronight, self.__minalt)
            self.__moonGraph_1.separate_list(self.__listdelta, self.__listlum, self.__listalt, 1)
            self.__moonGraph_2.point_scatter(self.__timesf,self.__illumination,self.__moonaltitude)

            plt.xlabel('Days from ' + str(self.__start)) 
            plt.hold(True)
            plt.grid()
            plt.show()
        else:
            pass

##printing class: plots, scatters

class MoonGraph(object):

    def __init__(self, astronight, minalt):
        self.__astronight = astronight
        self.__minalt = minalt


    def separate_list(self,listdelta,listlum,listalt,n): 
        for temp,lumi,altit in zip(listdelta,listlum,listalt):
            if temp == []: #to avoid error
                pass
            else:  ###change_object_type
                qtimes = u.Quantity(temp)
                illum = np.array(lumi)
                moonalt = Latitude(altit)
                
            if temp == []: #to avoid error
                pass
            elif n == 0:
                self.blue_plot(qtimes,moonalt) 
            elif n == 1:
                self.red_interp_plot(qtimes,moonalt) 

#astronomical nights in pink shadow, days in white
    def filling_graph(self,ntime,nsunaltazs_1):

        plt.fill_between(ntime.to('d').value, 0, 90,
                         nsunaltazs_1.alt <= self.__astronight, 
                         color='#e7d9e7', zorder=-4)
        
        plt.fill_between(ntime.to('d').value, -90, 90,
                         nsunaltazs_1.alt >= self.__astronight,
                         color='#ffffff', zorder=0)
        

#plt.xlim(-12, 12)
#plt.xticks(np.arange(13)*2 -12)
        plt.ylim(self.__minalt, 70)
       # plt.xlabel('Days from start') #+ str(start) ##
        plt.ylabel('Altitude [deg]')

        font = {'family' : 'normal',
                'weight' : 'bold',
                'size'   : 20}
        
        matplotlib.rc('font', **font)

#blue plotting line
    def blue_plot(self,time,moonalt): #,delta,alt):

        plt.plot(time, moonalt, 'b', zorder = -3)


#red plotting line
    def red_interp_plot(self,time,moonalt): #,delta,alt):

        plt.plot(time, moonalt, 'r', zorder = -3)


   
#this one scatters useful points on the graph
    def point_scatter(self,time,illum,moonalt):
        
        plt.scatter(time, moonalt,
                    c=illum, label='Moon illumination', lw=0, s=120,
                    cmap='viridis', zorder = 1)
        plt.colorbar().set_label('Illumination')
        plt.legend(loc='upper left')
        

#this class is not working properly yet
#file tkmoon3.py
class MoonTkinter(object):
    def __init__(self):
        pass #por si me hace falta

    def config_moon_plan(self, place, time, deltamins, deltadays, points, npoints): 
        self.__moonPlan_1 = MoonPlan(place)
        self.__start = self.__moonPlan_1._set_start_time(time, 
                                                         deltamins, deltadays, points)
        print(self.__start)
        dictimes = self.__moonPlan_1._set_interval_time() 
        self.__intervaltime = dictimes['interval']
        self.__moonPlan_1._set_framing()
        self.__sunaltaz_1 = self.__moonPlan_1.get_sun_altazs()
        self.__moon_pars = self.__moonPlan_1.get_moon_param()

    def guimoon(self):
        app = Tkinter.Tk ()
        app.title("Moon parameters")
        app.geometry('450x100+200+100')
    
 
        Tkinter.Button(app, bg = "#e7d9e7", activebackground = "#e7dbc7", 
                       highlightbackground="cyan", 
                       text='Moonshadow. Click to close!',
                       command=app.destroy).pack(pady=20)
   

        flipper = Tkinter.IntVar()
        o = 0
        f = 0
        def flip_it(o,f):
            
            if flipper.get() == 1:
                
                if o == 0:
                    print("on")
                    o = 1
                    f = 1
                else:
                    print("updating data")
                for c in range(len(self.__intervaltime)):
                    print('Date and time: ' + str(self.__intervaltime[c]) + 
                          '\n Altitude: ' + str(self.__moon_pars['moon altitude'][c])+ 
                          '\n Azimuth: ' + str(self.__moon_pars['moon azimuth'][c]) + 
                          '\n Illumination: ' + str(self.__moon_pars['moon illumination'][c]) + 
                          '\n \n')
                app.after(10000, lambda: flip_it(o,f))
            
            else:
                if f == 1:
                    print("off") #evitar el after en este caso
                
                    f = 0
                    o = 0
                    app.after_cancel(lambda: flip_it(o,f))
        
            
        Tkinter.Checkbutton(app, variable = flipper,
                            command = lambda: flip_it(o,f),
                            text = "Start/Stop").pack()
        app.mainloop()
