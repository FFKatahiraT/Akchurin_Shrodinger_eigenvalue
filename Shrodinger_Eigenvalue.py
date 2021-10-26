import matplotlib.pyplot as plt
import numpy as np

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def plotter(Func, name):
	linestyle=['solid', 'dotted', 'dashed', 'dashdot', (0, (5,10)), (0, (5,1)), (0, (3,10,1,10)), (0, (3,1,1,1))]
	i=0
	for i in range(len(Func)):
		plt.plot(r, Func[i], label=name +' E='+str(lmb*e)+'[eV]', linestyle=linestyle[i])
	plt.title(name)
	plt.legend(loc='best')
	plt.grid()
	plt.ylabel(str(name))
	plt.xlabel(r'$r$ [m]')
	plt.tight_layout()
	plt.savefig('graphs/'+str(name)+str('.eps'))
	plt.close()

def progonka_start_left(A, B, C, F, RB_Condition_alpha, RB_Condition_beta, LB_Condition):
	Func = np.zeros(len(r))
	alpha = np.zeros(len(r))
	beta = np.zeros(len(r))
	#RIGHT BORDER CONDITION 
	alpha[len(r)-1] = RB_Condition_alpha
	beta[len(r)-1] = RB_Condition_beta
	for i in reversed(range(len(r)-1)):	#CALCULATING ALPHA & BETA FROM R TO 0
		alpha[i] = A[i]/(C[i]-B[i]*alpha[i+1])
		beta[i] = (F[i]+B[i]*beta[i+1])/(C[i]-B[i]*alpha[i+1])
	#LEFT BORDER CONDITION
	Func[0] = LB_Condition
	for i in range(len(r)-1):	#CALCULATING A FUNCTION FROM 0 TO R
		Func[i+1] = alpha[i+1]*Func[i]+beta[i+1]		
	return Func

def psi_calc():	
	global lmb
	#SOLVING EQ ON psi
	lmb = 10**-20	#initial approx  sum(V)*h/2
	psi = np.ones(n+1)
	A=np.ones(n+1)
	B=np.ones(n+1)
	C=np.ones(n+1)
	for i in range(1, n):
		A[i] = -h_plank**2/(2*m*h**2)
		B[i] = -h_plank**2/(2*m*h**2)
		C[i] = -h_plank**2/(m*h**2)
	step = 1

	while True:
		F =  (V-lmb)*psi
		# F = np.ones(n+1)
		# for i in range(n+1):
		# 	if V[i]>lmb:
		# 		F[i] = (V[i]-lmb)*psi[i]
		# 	else:
		# 		F[i] = 0 
		#CONSIDER y[i]=alpha[i]*y[i-1]+beta[i]
		#LEFT BORDER CONDITION
		LB_Condition = 10**-60
		#RIGHT BORDER CONDITION 
		RB_Condition_alpha = 10**-60
		RB_Condition_beta = 10**-60	#RB_Condition_beta IS NOT 0 TO AVOID DIVISION BY ZERO ERROR
		psi = progonka_start_left(A, B, C, F, RB_Condition_alpha, RB_Condition_beta, LB_Condition)
		#NORMALIZATION integral{psi^2}=1
		summ = sum(psi**2)*h 	#calc integral
		print(summ, ' summ')
		# print(psi, ' psi')
		# print(F, 'F\n')
		psi = np.sqrt(psi**2/summ)	#bring value to the 1
		if step%10==0:
			print(lmb, ' lmb\n')
			if eps+1 > summ > 1-eps:
				break
			elif summ < 1:
				lmb *= 1.01
			elif summ > 1:
				lmb /= 1.01
		summ_prev = summ
		step+=1
	return psi

##########PARAMETERS##########
m = 1.6726219*10**-27 #kilograms proton mass
h_plank = 6.62607004*10**-34 #[m**2*kg/s]
e = 1.60217662*10**-19	#Electron charge
R = 10**-9	#[m] potential pit size
eps = 0.01	#error value
n = 1000
h = R/(n+0.5)
r = np.zeros(n+1)
for i in range(n+1):
	r[i] = h/2+i*h
V = 10**36*(r-R/2)**2	#particle potential energy
V[:] = 0
V[1:n-n//3] = 0
##############################
psi = psi_calc()

plotter([psi], r'$\psi$')
plotter([psi**2], r'$\psi^2$')
plotter([V], 'potential energy')