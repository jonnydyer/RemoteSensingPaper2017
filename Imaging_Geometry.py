
import numpy as np
from math import *
from matplotlib import rc
import matplotlib.pyplot as plt

rc('figure', figsize=(3.5,2.5))
rc('legend', fontsize='x-small')
rc('font', family='serif')

mu = 398600.4418  # km^3 / s^2
Re = 6378.1       # km

def earth_angle(ona, alt=600, Re=6371.):
    '''
    Return earth_angle for a given off-nadir-angle.
    Earth angle is the angle between line connecting satellite
    to center of earth and line connecing los / earth intersection
    point and center of earth.
    All units radians
    '''
    alpha = np.arcsin((Re + alt) * np.sin(ona) / Re) - ona
    return alpha

def pass_time(v_gnd, ona0, alt, Re=6378):
    mu = 398600.444
    sat_eca_rate = 1. / sqrt((Re + alt) ** 3 / mu)
    eca = earth_angle(ona0, alt, Re=6371.)
    del_eca = v_gnd / Re
    eca_rate = sat_eca_rate - del_eca
    return eca * 2. / eca_rate

def slew_time(slew_angle):
    '''
    Returns approximate time to complete slew
    for two targets separated by slew_angle
    in the spacecrafts coordinate frame.
    slew_angle in radians
    '''
    angle = [0, 5, 10, 15, 20, 30, 60, 100, 180]
    time = [-13, 8, 13, 16, 18, 21, 26, 29, 37]
    const_time = 13.
    return np.interp(slew_angle * 180. / pi, angle, time, right = 37., left=8.) + const_time

def target_relationships(N=100):
    '''
    Generates N-1 angular separation distances 
    and target lengths using "target of interest" modle
    for N "sequential" targets
    '''
    sep_dist = stats.poisson(100)
    len_dist = stats.poisson(50.)
    stats.poisson()
    return pi/180 * 25. / 100. * sep_dist.rvs(N-1), len_dist.rvs(N)

class Orbit(object):
    def __init__(self, alt):
        super(Orbit, self).__init__()
        self.alt = alt
        
    @property
    def alt(self):
        return self._alt
    
    @alt.setter
    def alt(self, value):
        self._alt = value
        self.a = value + Re                       # km
        self.per = 2 * pi * sqrt(self.a**3 / mu)     # s
        self.V_orb = sqrt(mu / self.a)               # km/s
        self.V_gnd = 2 * pi * Re / self.per          # km/s
        self.rho_horiz = asin(Re / (self.a))
        self.eca_eclipse = asin(Re / self.a) * 2.
        self.eclipse_per = self.eca_eclipse / 2 / pi * self.per
        self.eclipse_frac = (self.per - self.eclipse_per) / self.per
    
    def ona_to_eca(self, ona):
        '''
        Return earth central angle, eca, given spacecraft
        off-nadir angle ona.
        '''
        return np.arcsin(np.sin(ona)*(Re+self.alt)/Re)-ona
    
    def eca_to_ona(self, eca):
        '''
        Return off-nadir-angle, ona, given Earth central
        angle, eca
        '''
        
        return np.arctan(np.sin(self.rho_horiz) * np.sin(eca) /                          (1 - np.sin(self.rho_horiz) * np.cos(eca)))
    
    def strip_len(self, ona1, ona2=None):
        '''
        Returns strip length in km between two off-nadir angles
        '''
        if ona2 is None:
            ona2 = 0.-ona1
        return Re*(self.ona_to_eca(ona1) - self.ona_to_eca(ona2))
    
    def collect_time(self, scanrate, ona1, ona2=None):
        l = self.strip_len(ona1, ona2)
        return l / scanrate
        


# In[6]:

alt = 500       # km
gnd_rate = np.linspace(1., 7., 50)
target_length = 50#np.array([10, 25, 50, 100, 250, 500])
swath = 1.8e-6 * alt * (2560 * 3 - 2*200)
land_min_per_day = 167.7     # From Jim's ISIS analysis, over Antartica


# In[7]:

o = Orbit(alt)
print o.V_orb, o.V_gnd, o.eclipse_per/60.

land_min_per_orb = land_min_per_day / (86400. / o.per)

print 'Max Strip Length: %.1f' % o.strip_len(30.*pi/180.)


# In[8]:

target_spacing = np.array([5., 25., 75, 200]) / alt

targ_per_min = np.zeros([len(gnd_rate), len(target_spacing)])
duty_cycle = np.zeros_like(targ_per_min)
beta = np.zeros_like(targ_per_min)
targ_per_orb = np.zeros_like(targ_per_min)
targ_per_day = np.zeros_like(targ_per_min)
min_img_per_day = np.zeros_like(targ_per_min)
area_per_day = np.zeros_like(targ_per_min)
land_img_min_per_day = np.zeros_like(targ_per_min)
land_area_per_day = np.zeros_like(targ_per_min)

for i,s in enumerate(target_spacing):
    collect_times = target_length / gnd_rate
    eca_limit = collect_times * (o.V_gnd - gnd_rate) / Re / 2.
    ona_limit = o.eca_to_ona(eca_limit)
    slew_times = slew_time(2 * ona_limit + s)
    duty_cycle[:,i] = collect_times / (collect_times + slew_times)#np.tile(slew_times, (len(collect_times), 1)).T)
    targ_per_min[:,i] = 60. / (collect_times + slew_times)
    targ_per_orb[:,i] = o.per * (1.-o.eclipse_frac) / (collect_times + slew_times)
    targ_per_day[:,i] = land_min_per_day * 60. / (collect_times + slew_times)
    beta[:,i] = targ_per_orb[:,i] * collect_times / o.per
    min_img_per_day[:,i] = 24. * 60. * (1.-o.eclipse_frac) * duty_cycle[:,i]
    land_img_min_per_day[:,i] = land_min_per_day * duty_cycle[:,i]
    area_per_day[:,i] = min_img_per_day[:,i] * gnd_rate * 60. * swath
    land_area_per_day[:,i] = land_img_min_per_day[:,i] * gnd_rate * 60. * swath

f = plt.figure(1)
plt.plot(gnd_rate, beta)
plt.xlabel('Scan rate (km/s)')
plt.ylabel('Collection Duty Cycle')
plt.xticks([1,3,5,7])
plt.ylim(0, 0.3)
plt.xlim(1, 7)
plt.legend(['%d km' % t for t in target_spacing * alt], 
           loc='best')
plt.grid(True)
#plt.title(r'Collection Duty Cycle, $\beta$, for an agile ' +\
#         ' \n spacecraft with various target spacings',
#         fontsize='small')
plt.tight_layout()
f.savefig('figures/collection_dc.pgf')

plt.figure(2)
plt.plot(gnd_rate, duty_cycle)
plt.xlabel('Scan rate (km/s)')
plt.ylabel('Imaging Duty Cycle')
plt.legend([str(t) + ' km' for t in target_spacing * alt], 
           loc='best')
plt.grid(True)
#gcf().savefig('duty_cycle.pdf')

plt.figure(3)
plt.plot(gnd_rate, targ_per_day)
plt.xlabel('Scan rate (km/s)')
plt.ylabel('Target per day')
plt.legend([str(t) + ' km' for t in target_spacing * alt], 
           loc='best')
plt.grid(True)
#gcf().savefig('targ_per_day.pdf')

plt.figure(4)
plt.subplot(211)
plt.plot(gnd_rate, area_per_day / 1000.)
plt.grid(True)
plt.ylabel('Area per day (\'1000 sq. km)')
plt.subplot(212, sharex=plt.gca())
plt.plot(gnd_rate, land_area_per_day / 1000.)
plt.xlabel('Scan rate (km/s)')
plt.ylabel('Land area per day (\'1000 sq. km)')
plt.legend([str(t) + ' km' for t in target_spacing * alt], 
           loc='best')
plt.grid(True)
#gcf().savefig('area_per_day.pdf')

plt.figure(5)
plt.subplot(211)
plt.plot(gnd_rate, min_img_per_day)
plt.grid(True)
plt.ylabel('Imaging Minutes per Day')
plt.subplot(212, sharex=plt.gca())
plt.plot(gnd_rate, land_img_min_per_day)
plt.xlabel('Scan rate (km/s)')
plt.ylabel('Imaging Minutes per Day over land')
plt.legend([str(t) + ' km' for t in target_spacing * alt], 
           loc='best',)
plt.grid(True)
#gcf().savefig('min_per_day.pdf')

#print collect_times + np.tile(slew_times, (len(collect_times), 1)).T