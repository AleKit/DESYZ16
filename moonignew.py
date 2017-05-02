#see which nights are the best to observe
#install: astropy v1.2.1, astroplan

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
#requires matplotlib-1.5.2
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)

#Import the packages necessary for finding coordinates and making coordinate transformations

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_moon, solar_system_ephemeris, get_sun, Latitude
#from astropy.coordinates.angle_utilities import angular_separation

from astroplan.moon import moon_phase_angle, moon_illumination
from astroplan import Observer

import datetime
import time

from scipy import interpolate

import pytz

#calculations interval of time, starting time, telescope location
#some empty lists to append
#get the altitude (and azimuth) of the Sun

###input
mod = input('graph, text or both? in quotation marks: ')
if mod != 'graph' and mod != 'text' and mod != 'both':
    print('Error: Wrong word. Write only: graph, text or both')
    raise SystemExit(0)

##falta hacer:
###input numero de horas que se quieren ver
###input numero de puntos en el intervalo

##constants
astronight = -18*u.deg
minalt = 30
maxillum = 0.6
dd = 3
#yr = 2016
#mth = 10
#dy = 21
times = np.linspace(0, dd, 60)*u.d  
ntimes = np.linspace(0, dd, 6000)*u.d 
#start = Time.now()
#para poner un tiempo especifico

#ctime = Time('2017-1-16 00:00:00')
ctime = Time('2016-10-20 00:00:00')
print(ctime)
#start = Time('2016-10-21 23:00:00') - utcoffset
#input write the starting time in the following format y poner el que tengo de ejemplo
berlin = EarthLocation(lat=52.5*u.deg, lon=13.4*u.deg, height=34*u.m)

dstart = time.mktime((2016,10,21,0,0,0,0,0,0))
end = dstart+dd*86400

#if time.localtime(dstart).tm_isdst != time.localtime(end).tm_isdst:
#    print('Hour shifting!!')
#    raise SystemExit(0)
#elif time.localtime(dstart).tm_isdst==0:
#    utcoffset = 1*u.hour  # Berlin Standard Time
#else:
#    utcoffset = 2*u.hour #Berlin Daylight Time
start = ctime #- utcoffset


#illum = []
#delta = []
#moonaltazs = []
#sunaltazs = []
#moonaltaz = Latitude([])
times_1 = start + times
times_2 = start + ntimes
frame_1 = AltAz(obstime=times_1, location=berlin)
frame_2 = AltAz(obstime=times_2, location=berlin)
sunaltazs_1 = get_sun(times_1).transform_to(frame_1)
sunaltazs_2 = get_sun(times_2).transform_to(frame_2)



frame_ori1 = AltAz(obstime = times_1, location = berlin, pressure = 1015*u.hPa, relative_humidity = 0.7, temperature = 15.0*u.deg_C, obswl = 0.5*u.micron)
frame_ori2 = AltAz(obstime = times_2, location = berlin, pressure = 1015*u.hPa, relative_humidity = 0.7, temperature = 15.0*u.deg_C, obswl = 0.5*u.micron)


sunaltazs_ori1 = get_sun(times_1).transform_to(frame_ori1)
sunaltazs_ori2 = get_sun(times_2).transform_to(frame_ori2)


#empty lists of lists
listdelta = []
listlum = []
listalt = []

pic = plt.figure()
pic.suptitle('Moon altitude and illumination', fontsize=38)

#if text
if mod == 'text' or mod == 'both':
    dtxt = open ('data.txt', 'w')
    dtxt.write('Nights we can observe the Moon shadow since ' + str(start) + '\n \n') #he quitado +utcoffset porque voy a usar UTC, lo voy a quitar siempre que lo encuentre

####functions####


#this one calculates illumination and altitude for each time
#if you write graph as input, the function calls chobj
#if you write text as input, the function writes every useful time, illumination and altitude in a plain text file
#if you write both, it does both
#creo que no la estoy usando ya
def separateap():
    for ctime in times:
        t = start + ctime
        frame = AltAz(obstime=t, location=berlin)
        sunaltaz = get_sun(t).transform_to(frame)
        if sunaltaz.alt < astronight:
        #sunaltazs.append(sunaltaz.alt)
        #print(sunaltaz)
        #print(sunaltaz.alt)
        #print('Moon phase angle:', moon_phase_angle(t, berlin))
            lum = moon_illumination(t, berlin)
        #print('Moon illumination:', moon_illumination(t, berlin))

            if lum < maxillum:
                delta.append(ctime)
                illum.append(lum)
                moon = get_moon(t,berlin)
                moonaltaz = moon.transform_to(frame)
           #print("Moon's Altitude = {0.alt:.2}".format(moonaltaz))
                moonaltazs.append(moonaltaz.alt)
               # print(moonaltaz.alt)
                if mod == 'graph':
                    pass
                elif mod == 'text' or mod == 'both':
                    if moonaltaz.alt > Latitude(minalt*u.deg):
                    
                        dtxt.write('Date and time: ' + str(t) + '\n' + 'Altitude: ' + str(moonaltaz.alt)+ '\n' + 'Illumination: ' + str(lum) + '\n \n')


    if mod == 'text':
        pass
    elif mod == 'graph' or mod == 'both':
        chobj(delta,illum,moonaltazs)

#this one calculates parameters without loops
def mooniap(place):
    moon = get_moon(times_1)
    lum = moon_illumination(times_1, place)
    moonaltaz = moon.transform_to(frame_1)
 
    #prinmooniap(times,lum,moonaltaz.alt,ntimes,sunaltazs_1,sunaltazs_2) 


    
    moonaltaz_ori = moon.transform_to(frame_ori1)

    prinmooniap(times,lum,moonaltaz_ori.alt,ntimes,sunaltazs_ori1,sunaltazs_ori2) 
    

#this one transforms lists into the needed classes and calls printogeap
#only useful if you want a graph
def chobj(delta,illum,moonaltazs):
    if delta != []:

#print ('delta tipo', type(delta))
        qtimes = u.Quantity(delta) 
#print ('delta tipo', type(qtimes))
#print('MOONALTAZ TIPO 0:', type(moonaltaz.alt))
#print(moonaltaz.alt)
#print('MOONALTAZS TIPO 1:', type(moonaltazs))
#print(moonaltazs)
#print(type(moonaltazs))
#print(illum)
#print(type(illum))
#print(delta)
#print(type(delta))
        illum = np.array(illum, dtype=float)
#print('ILLUM TIPO:', type(illum))
#print(illum)
#delta = np.array(delta)
#print('DELTA TIPO:', type(delta))
#moonaltazs = np.array(moonaltazs)
        moonalt = Latitude(moonaltazs)
#moonaltazs = np.array(moonaltazs)
#print('MOONALT', moonalt)
#print('MOONALT TIPO 2:', type(moonalt))
#sunaltazs = np.array(sunaltazs)

#print(delta)
#print(moonaltazs)
#print(moonaltazs[11])
#print(sunaltazs)
#print(illum)
        printogeap(qtimes,illum,moonalt)
    else:
        print("No data")

#this one separates the lists of lists to make a nice plot and calls 
#prinseparateap, only useful if you want a graph
def prinlist(listdelta,listlum,listalt,n,sunaltaz):
    #print(n)
    for temp,lumi,altit in zip(listdelta,listlum,listalt):
        if temp == []: #to avoid error
            pass
        else:
            qtimes = u.Quantity(temp)
            illum = np.array(lumi)
            moonalt = Latitude(altit)

        if mod == 'graph' or mod == 'both':
            if temp == []: #to avoid error
                pass
            elif n == 0:
                prinseparateap(qtimes,moonalt) 
            elif n == 1:
                printerpolate(qtimes,moonalt)  


#this one has the parameters for the graph, the astronomical nights are light pink and the days white
##Important: sunalt parameters
def prinparam():
#plt.fill_between(times.to('hr').value, 0, 90,
#                 illum < 0.6, color='0.5', zorder=0)

##############IMPORTANT - SUNALT
    plt.fill_between(ntimes.to('d').value, 0, 90,
                     sunaltazs_ori2.alt <= astronight, color='#e7d9e7', zorder=-4)

    plt.fill_between(ntimes.to('d').value, -90, 90,
                     sunaltazs_ori2.alt >= astronight, color='#ffffff', zorder=0)
 ####IMPORTANT - SUNALT

#plt.fill_between(times.to('hr').value, -50, 70,
#                 sunaltazs_1.alt < 0*u.deg, color='k', zorder=0)
    plt.colorbar().set_label('Illumination')
    plt.legend(loc='upper left')
#plt.xlim(-12, 12)
#plt.xticks(np.arange(13)*2 -12)
    plt.ylim(minalt, 70)
    plt.xlabel('Days from ' + str(start)) #+utco
    plt.ylabel('Altitude [deg]')
   
    font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 28}

    matplotlib.rc('font', **font)
    
 
#this one links points on the graph (not used)
def prinseparateap(qtimes,moonalt): #,delta,alt):

    plt.plot(qtimes, moonalt, 'b', zorder = -3)


#this one draws the interpolation lines
def printerpolate(qtimes,moonalt): #,delta,alt):

    plt.plot(qtimes, moonalt, 'r', zorder = -3)


   
#this one scatters useful points on the graph
def printogeap(qtimes,illum,moonalt):
        
    plt.scatter(qtimes, moonalt,
                c=illum, label='Moon illumination', lw=0, s=120,
                cmap='viridis', zorder = 1)

    prinparam()

#interpolation things
#creo listas de listas para obtener solo los valores que entran dentro de los intervalos correctos
def prinmooniap(qtimes,illum,moonalt,ntimes,sunaltaz,nsunaltaz):

    f = interpolate.interp1d(qtimes,moonalt,kind='cubic')
    g = interpolate.interp1d(qtimes,illum,kind='cubic')
    delta = []
    lumi = []
    moonaltis = []
    ndelta = []
    nlumi = []
    nmoonaltis = []
   # listdelta = []
   # listlum = []
   # listalt = []
    deltaobj = []
    lumiobj = []
    altiobj = []

#este loop es para recolectar todos los puntos utiles para scatter, luego llamo a chobj con estas listas
    for c in range(len(qtimes)):

        if sunaltaz.alt[c] <= astronight:
        #print('Moon phase angle:', moon_phase_angle(t, berlin))
        #print('Moon illumination:', moon_illumination(t, berlin))
   #         co = 0
            if illum[c] < maxillum:
   #             delta.append(qtimes[c])
   #             lumi.append(illum[c])
           #print("Moon's Altitude = {0.alt:.2}".format(moonaltaz))
   #             moonaltis.append(moonalt[c])
                deltaobj.append(qtimes[c])
                lumiobj.append(illum[c])
                altiobj.append(moonalt[c])
   #         elif co == 0:
   #         listdelta.append(delta)
   #         listlum.append(lumi)
   #         listalt.append(moonaltis)
   #         delta = []
   #         moonaltis = []
   #         lumi = []
   #         co = 1
    #prinlist(listdelta,listlum,listalt,0,sunaltaz)

#ahora utilizo lo que he obtenido al interpolar    
    listdelta = []
    listlum = []
    listalt = []
    conig = 0
    for c in range(len(ntimes)):

        if nsunaltaz.alt[c] <= astronight+0.5*u.deg:
        #print('Moon phase angle:', moon_phase_angle(t, berlin))
        #print('Moon illumination:', moon_illumination(t, berlin))
            co = 0
            if g(ntimes[c]) < maxillum:
                ndelta.append(ntimes[c])
                nlumi.append(g(ntimes[c]))
           #print("Moon's Altitude = {0.alt:.2}".format(moonaltaz))
                nmoonaltis.append(f(ntimes[c])*u.deg)
        elif co == 0:
            listdelta.append(ndelta)
            listlum.append(nlumi)
            listalt.append(nmoonaltis)
            delta = []
            moonaltis = []
            lumi = []
            co = 1
    #print(listdelta)

        if (mod == 'text' or mod == 'both') and c < len(ntimes)-1:
            if c % 500 == 0:#Igor
                print("Step: "+ str(c))

            if  g(ntimes[c]) < maxillum and nsunaltaz.alt[c] <= astronight and  (nsunaltaz.alt [c+1] > astronight or nsunaltaz.alt[c-1] > astronight) and f(ntimes[c]) >= minalt:
                #Here I want to interpolate around the critical value wich is exactly -18 deg (astronight).
                
      
                times_subset = ntimes[c-3 : c+4]#Igor: Here I take the 6 values around the closest value
                sunalta_subset = nsunaltaz.alt[c-3 : c+4]
                inter_time = interpolate.interp1d(sunalta_subset, times_subset, kind='cubic')

                start_time = inter_time(astronight) # Igor: It would be better to use a constant instead of hard-coding "-18" 
                inter_moon = f(start_time)
                inter_ill = g(start_time)
            ##prueba
              #  ber = Observer(location=berlin, name="Proto")
              #  print(ber.twilight_morning_astronomical(start+ntimes[c]))
              #  print((start+float(start_time)).format('jd'))
            ##end
                dtxt.write('Date and time: ' + str(start+float(start_time)*u.d) + '\n' + 'Altitude: ' + str(round(inter_moon,3))+ '\n' + 'Illumination: ' + str(round(inter_ill,3)) + '\n \n') #+utcof
                
               # print('\n' + "Mis Valores: " + str(start_time*u.d)+ ' Date and time: ' + str(start+float(start_time)*u.d) + '\n' + 'Altitude: ' + str(inter_moon)+ '\n' + 'Illumination: ' + str(inter_ill) + '\n') #Igor
              
               # print('tipo 1 ' + str(ntimes[c]) + ' Date and time: ' + str(start+ntimes[c]) + '\n' + 'Altitude: ' + str(f(ntimes[c]))+ '\n' + 'Illumination: ' + str(g(ntimes[c])) + '\n \n')#quitar
               
                if conig == 0:
                    conig = 1
                
            elif g(ntimes[c]) < maxillum and nsunaltaz.alt[c] <= astronight and f(ntimes[c]) >= minalt and (f(ntimes[c+1]) < minalt or f(ntimes[c-1]) < minalt):  

                times_subset = ntimes[c-3 : c+4]
                
                val = []
                
                for x in range(-3, 4):
                    val.append(f(ntimes[c+x]))

                moon_subset = np.array(val, dtype=float)

                inter_time = interpolate.interp1d(moon_subset, times_subset, kind='cubic')

                start_time = inter_time(minalt) # Igor: It would be better to use a constant instead of hard-coding "-18" 
                inter_moon = f(start_time)
                inter_ill = g(start_time)
                
               # print('\n' + "Mis Valores: " + str(start_time*u.d)+ ' Date and time: ' + str(start+float(start_time)*u.d) + '\n' + 'Altitude: ' + str(inter_moon)+ '\n' + 'Illumination: ' + str(inter_ill) + '\n') #Igor
                               
               # print('tipo 2 ' + str(ntimes[c]) + ' Date and time: ' + str(start+ntimes[c]) + '\n' + 'Altitude: ' + str(f(ntimes[c]))+ '\n' + 'Illumination: ' + str(g(ntimes[c])) + '\n \n')#quitar
                

                dtxt.write('Date and time: ' + str(start+float(start_time)*u.d) + '\n' + 'Altitude: ' + str(round(inter_moon,3))+ '\n' + 'Illumination: ' + str(round(inter_ill,3)) + '\n \n') #+utcof

                if conig == 0:
                    conig = 1
#######faltaba considerar que la iluminacion cambie de buena a mala
#########
            elif g(ntimes[c]) < maxillum and nsunaltaz.alt[c] <= astronight and f(ntimes[c]) >= minalt and (g(ntimes[c+1]) > maxillum or g(ntimes[c-1]) > maxillum):  

                times_subset = ntimes[c-3 : c+4]
                
                val = []
                
                for x in range(-3, 4):
                    val.append(g(ntimes[c+x]))

                moon_subset = np.array(val, dtype=float)

                inter_time = interpolate.interp1d(moon_subset, times_subset, kind='cubic')

                start_time = inter_time(maxillum) # Igor: It would be better to use a constant instead of hard-coding "-18" 
                inter_moon = f(start_time)
                inter_ill = g(start_time)
                
               # print('\n' + "Mis Valores: " + str(start_time*u.d)+ ' Date and time: ' + str(start+float(start_time)*u.d) + '\n' + 'Altitude: ' + str(inter_moon)+ '\n' + 'Illumination: ' + str(inter_ill) + '\n') #Igor
                               
               # print('tipo 2 ' + str(ntimes[c]) + ' Date and time: ' + str(start+ntimes[c]) + '\n' + 'Altitude: ' + str(f(ntimes[c]))+ '\n' + 'Illumination: ' + str(g(ntimes[c])) + '\n \n')#quitar
                

                dtxt.write('Date and time: ' + str(start+float(start_time)*u.d) + '\n' + 'Altitude: ' + str(round(inter_moon,3))+ '\n' + 'Illumination: ' + str(round(inter_ill,3)) + '\n \n') #+utcof

                if conig == 0:
                    conig = 1
    if (mod == 'text' or mod == 'both') and conig == 0:
        dtxt.write('No nights available in this period: '  + str(start) + ' - ' + str(start+ntimes[len(ntimes)-1])) #+utcof
        
    if mod == 'graph' or mod == 'both':
        prinlist(listdelta,listlum,listalt,1,nsunaltaz)
       # printerpolate(ntimes, f(ntimes))

    #plt.scatter(qtimes, moonalt,
    #            c=illum, label='Moon illumination', lw=0, s=40,
    #            cmap='viridis', zorder = 3)
    #plt.plot(qtimes, moonalt, 'b', zorder = -3)
   # plt.plot(ntimes,f(ntimes),'r',zorder=-2)
    #prinparam()


    #printogeap(qtimes,illum,moonalt)
    chobj(deltaobj,lumiobj,altiobj)
    #if mod == 'text' or mod == 'both':
    #    for j in range(len(qtimes)):   
    #        if moonalt[j] == 30*u.deg or sunaltazs_2.alt[i] == -18*u.deg:


###############################################


#separateap()
mooniap(berlin)

if mod == 'text' or mod == 'both':
    dtxt.close()
    infile = open('data.txt', 'r')
    print(infile.read())

if mod == 'graph' or mod == 'both':
    plt.hold(True)
   # prinlist()


    plt.show() 



