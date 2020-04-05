import sys
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 
import numpy as np
import savitzy_golay as sv_g
import copy
from scipy.integrate import quad
from scipy.integrate import trapz

u = open(sys.argv[1],"r") #file with urea information
w = open(sys.argv[2],"r") #file with water information
c = open(sys.argv[3],"r") #file to calculate number of urea molecules (PDB file)
#l=sys.argv[4]

#v = open(sys.argv[4],"r") #calculate box volume
#print c, v

#for fitting the tail of radial distribution function
def func(x, a, b, c, d, e):
     return 1 + a*math.exp(-b*(x - c))*math.sin(d*(x - e))

#for evaluating the function
def KB_func(x, a, b, c, d, e):
     return (a*math.exp(-b*(x - c))*math.sin(d*(x - e)))*4.0*math.pi*x*x*(1.0 - (3*x)/(4.0*15.0) + ((x*x*x)/(16.0*15.0*15.0*15.0)))

#KB integrals for urea
sum_u = 0
xurea = []
yurea = []
for i in u.readlines():
	columns = i.split()
	#print columns
	if (float(columns[0]) > 15.0):
		break	
	yurea.append((float(columns[1])))
	xurea.append(float(columns[0]))
#print yurea
#print xurea
y = np.array(yurea)  # distance upto 15
x = np.array(xurea)
yhat = sv_g.savitzky_golay(y, 11, 3) #probability density

#plt.plot(x,y, color = "black")
urea_rdf = yhat

#print("urea_rdf and dist = ", zip(urea_rdf, x), len(urea_rdf))
#plt.plot(x,urea_rdf)
#plt.show()
counter = 0
indirect_urea_rdf = []
indirect_x = []

# Take only those values in indirect_urea_rdf which occur after the graph crosses 1.0 thrice
for i in range(len(urea_rdf)):
	if(counter == 3):
		indirect_urea_rdf.append(urea_rdf[i])
		indirect_x.append(xurea[i])
		#print indirect_x
	elif(counter % 2 == 0):	
		if(urea_rdf[i] > 1):
			counter += 1
			if counter == 3:
				indirect_urea_rdf.append(urea_rdf[i])
				indirect_x.append(xurea[i])
	elif(counter % 2 == 1):	
		if(urea_rdf[i] < 1):
			counter += 1

indirect_urea_rdf = np.array(indirect_urea_rdf)
indirect_x = np.array(indirect_x)
#print indirect_urea_rdf, indirect_urea_rdf.shape
#print indirect_x
#print np.count_nonzero(urea_rdf)
#initial guesses
einit = indirect_x[0] # This is ru3
#finding second maxima
rmin2 = 0.0
for i in range(0, len(indirect_x) - 2):
	if((indirect_urea_rdf[i + 1] - indirect_urea_rdf[i])*(indirect_urea_rdf[i + 2] - indirect_urea_rdf[i + 1]) < 0.0 and indirect_urea_rdf[i + 1] > 1.0):
		grmax2 = indirect_urea_rdf[i + 1]	
		ainit = grmax2 - 1.0 
		rmax2 = indirect_x[i + 1]
		storeindex_max = i + 2
		break
#finding second minima
for i in range(storeindex_max, len(indirect_x) - 2):
	#print (indirect_urea_rdf[i + 1] - indirect_urea_rdf[i])*(indirect_urea_rdf[i + 2] - indirect_urea_rdf[i + 1]), indirect_x[i+1]
	if((indirect_urea_rdf[i + 1] - indirect_urea_rdf[i])*(indirect_urea_rdf[i + 2] - indirect_urea_rdf[i + 1]) < 0.0 and indirect_urea_rdf[i + 1] < 1.0):
		grmin2 = indirect_urea_rdf[i + 1] 		
		rmin2 = indirect_x[i + 1]
		storeindex = i + 2
		break

#finding fourth unity, i.e. ru4
for i in range(0, len(indirect_x)):
	if(indirect_urea_rdf[i] < 1.0):
		ru4 = indirect_x[i]	
		break

#print("ru3 = ", einit)
#print("rmin2 = ", rmin2)
#print("rmax2 = ", rmax2)
#print("ru4 = ", ru4)

if rmin2 == 0.0: #so we the function g(r) never goes below 1.0 and dies down.
	rmin2 = rmax2 + 1.5
	grmin2 = 0.99
	ru4 = rmax2 + 1.0
#print rmin2, rmax2, grmin2, grmax2, indirect_x[0], ru4
binit = (-1.0/(rmin2 - rmax2))*math.log((1.0 - grmin2)/(grmax2 - 1.0))
cinit = rmax2
dinit = math.pi/(ru4 - indirect_x[0])
binit = abs(binit)
#print binit

while True:
	try:		
		yy = copy.deepcopy(np.array(indirect_urea_rdf))
		xx = copy.deepcopy(np.array(indirect_x))				
		#popt_u, pcov_u = curve_fit(func,xx, yy, p0=[ainit, binit, cinit,dinit,einit])
		popt_u = [ainit, binit, cinit, dinit, einit] 	
		#print popt_u	
		radii = np.linspace(indirect_x[0], 20, num=100)
		fit_curve_urea = np.array([func(i, *popt_u) for i in radii])
		break
	
	except RuntimeError:
		print "oops"

plt.plot(radii,fit_curve_urea,color = "magenta", label = "urea")
plt.plot(x,yhat,"r--", label = "urea")
#plt.show()
#exit(0)

#numerical integration urea
index = np.where(x == indirect_x[0])[0][0] # index in x from where indirect_x starts
#print(np.where(x == indirect_x[0]))
#exit(0)
KB_value_urea = []
for i in range(0,index):
	ww = 1.0 - (3*x[i])/(4.0*15.0) + ((x[i]*x[i]*x[i])/(16.0*15.0*15.0*15.0)) # R = 15.0
        KB_value_urea.append((urea_rdf[i] - 1.0)*4*np.pi*x[i]*x[i]*(ww))

KB_value_urea = np.array(KB_value_urea)
I1 = trapz(KB_value_urea, x[0:index])
I2 = quad(KB_func, indirect_x[0], 15.0, args=(popt_u[0],popt_u[1],popt_u[2],popt_u[3], popt_u[4]))
#print ("I1 = ", I1)
#print ("I2 = ", I2)

#KB integrals for water
xwat = []
ywat = []
for i in w.readlines():
	columns = i.split()
	if (float(columns[0]) > 15.0):
		break	
	ywat.append((float(columns[1])))
	xwat.append(float(columns[0]))

y = np.array(ywat)
x = np.array(xwat)
yhat = sv_g.savitzky_golay(y, 11, 3)
plt.plot(x,y, color = "black")
plt.plot(x,yhat,"g--", label="water")

water_rdf = yhat
counter = 0
indirect_water_rdf = []
indirect_x = []
for i in range(len(water_rdf)):
	if(counter == 3):
		indirect_water_rdf.append(water_rdf[i])
		indirect_x.append(xwat[i])
	elif(counter % 2 == 0):	
		if(water_rdf[i] > 1):
			counter += 1
			if counter == 3:
				indirect_water_rdf.append(water_rdf[i])
				indirect_x.append(xwat[i])
	elif(counter % 2 == 1):	
		if(water_rdf[i] < 1):
			counter += 1

indirect_water_rdf = np.array(indirect_water_rdf)
indirect_x = np.array(indirect_x)

#initial guesses
einit = indirect_x[0]
#finding second maxima
rmin2 = 0.0
for i in range(0, len(indirect_x) - 2):
	if((indirect_water_rdf[i + 1] - indirect_water_rdf[i])*(indirect_water_rdf[i + 2] - indirect_water_rdf[i + 1]) < 0.0 and indirect_water_rdf[i + 1] > 1.0):
		grmax2 = indirect_water_rdf[i + 1]	
		ainit = grmax2 - 1.0 
		rmax2 = indirect_x[i + 1]
		storeindex_max = i + 2
		break

#print len(indirect_x), len(indirect_water_rdf)
#finding second minima
for i in range(storeindex_max, len(indirect_x) - 2):
	if((indirect_water_rdf[i + 1] - indirect_water_rdf[i])*(indirect_water_rdf[i + 2] - indirect_water_rdf[i + 1]) < 0.0 and indirect_water_rdf[i + 1] < 1.0):
		grmin2 = indirect_water_rdf[i + 1] 		
		rmin2 = indirect_x[i + 1]
		storeindex = i + 2	
		break

if rmin2 == 0.0: #so we the function g(r) never goes below 1.0 and dies down.
	rmin2 = rmax2 + 1.5
	grmin2 = 0.99
	ru4 = rmax2 + 1.0

#finding fourth unity
for i in range(0, len(indirect_x)):
	if(indirect_water_rdf[i] < 1.0):
		ru4 = indirect_x[i]
		break
#print rmin2, rmax2, grmin2, grmax2, indirect_x[0], ru4
binit = -1.0/(rmin2 - rmax2)*math.log((1 - grmin2)/(grmax2 - 1))
cinit = rmax2
dinit = math.pi/(ru4 - indirect_x[0])
binit = abs(binit)
#print binit
while True:
	try:
		yy = copy.deepcopy(np.array(indirect_water_rdf))
		xx = copy.deepcopy(np.array(indirect_x))		
		#popt_w, pcov_w = curve_fit(func, indirect_x, indirect_water_rdf,  p0=[ainit, binit, cinit,dinit,einit])
		popt_w = [ainit, binit, cinit, dinit, einit] 
		#print popt_w	
		radii = np.linspace(indirect_x[0], 20, num=100)
		fit_curve_water = np.array([func(i, *popt_w) for i in radii])
		break

	except RuntimeError:
		print "oops"
	

plt.plot(radii,fit_curve_water,color = "cyan", label = "water")
plt.legend()
plt.savefig("rdf_smooth_plot.eps",format="eps",dpi=300)	

#numerical integration water
index = np.where(x == indirect_x[0])[0][0]
print index
KB_value_water = []
for i in range(0,index):
	ww = 1.0 - (3*x[i])/(4.0*15.0) + ((x[i]*x[i]*x[i])/(16.0*15.0*15.0*15.0))
        KB_value_water.append((water_rdf[i] - 1.0)*4*np.pi*x[i]*x[i]*ww)

KB_value_water = np.array(KB_value_water)
I3 = trapz(KB_value_water, x[0:index])
I4 = quad(KB_func, indirect_x[0], 15.0, args=(popt_w[0],popt_w[1],popt_w[2],popt_w[3], popt_w[4]))
#print I3
print I4

prev = -1
count = 0
for i in c.readlines():
	columns = i.split()
	if(len(columns) > 3 and columns[3] == "UREA"):
		if(int(columns[4]) != prev):
			count += 1
			prev = int(columns[4])

#cols = v.readlines()[2].split()
#l = math.ceil(float(cols[1]))
#b = math.ceil(float(cols[5]))
#h = math.ceil(float(cols[9]))


#print ("n_urea = ", count)
print ((I1 + I2[0]) - (I3 + I4[0]))*(count/(33.0*33.0*33.0))
#print (sum_u - sum_w)*(count/(30.0*30.0*30.0))


'''
name = l[0:3]
if(name == "histidine_basic"):
        name ="hisb"
if(l == "asparagine"):
        name = "asn"
if(l == "histidine_neutral"):
        name = "hisn"
if(l == "glutamine"):
        name = "gln"
'''

wrte = open("kb_all.dat", "a")

wrte.write("x" + " " + str(((I3 + I4[0]) - (I1 + I2[0]))*(count/(33.0*33.0*33.0))) + "\n")
 
