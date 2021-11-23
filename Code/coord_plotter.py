# Used for testing the map plotting functionality.  

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import math


m = Basemap(llcrnrlon=-121,llcrnrlat=20,urcrnrlon=-62,urcrnrlat=51,
    projection='lcc',lat_1=32,lat_2=45,lon_0=-95)
plt.figure(num=1,figsize=(9,7))
m.drawcoastlines()
m.drawstates()
m.drawcountries()
m.drawmapboundary(fill_color='paleturquoise')
m.fillcontinents(color='peachpuff',lake_color='paleturquoise')

test_coords = [(40.44062479999999, -79.9958864), (40.32266845384615, -79.76826980769232), (40.2047121076923, -79.54065321538462), (40.08675576153845, -79.31303662307693), (39.96879941538461, -79.08542003076923), (39.85084306923076, -78.85780343846153), (39.732886723076916, -78.63018684615385), (39.61493037692307, -78.40257025384615), (39.49697403076922, -78.17495366153845), (39.37901768461538, -77.94733706923077), (39.26106133846153, -77.71972047692307), (39.143104992307684, -77.49210388461537), (39.02514864615384, -77.26448729230768), (38.9071923, -77.0368707)]

for (lat, long) in test_coords:
    (x, y) = m(long, lat)
    m.plot(x,y,marker='o',color='Red',markersize=3)
    # plt.annotate(city, xy = (x,y), xytext=(-20,20))

plt.show()