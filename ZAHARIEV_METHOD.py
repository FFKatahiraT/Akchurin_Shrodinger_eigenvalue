import numpy as np 
import matplotlib.pyplot as plt

def plotter(Func, E_lambda, name):
	linestyle=['solid', 'dotted', 'dashed', 'dashdot', (0, (5,10)), (0, (5,1)), (0, (3,10,1,10)), (0, (3,1,1,1))]*10
	plt.rcParams.update({'font.size': 8})
	i=0
	for i in range(len(Func)):
		plt.plot(r, Func[i], label='E='+str(round(E_lambda[i]/e,5))+'[eV]', linestyle=linestyle[i])
	plt.ylabel(name)
	plt.legend(loc='upper right')
	plt.grid()
	plt.xlabel(r'$r$ [m]')
	plt.tight_layout()
	plt.savefig('graphs/'+str(name)+str('.png'))
	plt.close()

def normalize(func):
	summ = sum(func)
	for i in range(len(func)):
		func[i] /= summ
	return func

def multiply_diagonal(diagonal, psi, E):
	multiple = 1
	for i in range(len(diagonal)):
		multiple *= diagonal[i]*psi[i]-E
	return multiple

def find_E(V, E, psi):
	#########DIAGONALS#########
	A = E_max*np.ones(n+1)	#A*psi(i-1)+-C*psi(i)+B*psi(i+1) = -F
	A[0] = 0
	C = -E_max*2+V 	
	B = E_max*np.ones(n+1)
	B[n] = 0
	###########################
	det = 1
	while not -eps < det < eps:	#find E
		if det < -eps:
			E /= 1+step
		elif det > eps:
			E *= 1+step
		det = multiply_diagonal(A, psi, E) + multiply_diagonal(C, psi, E) + multiply_diagonal(B, psi, E)
		print(det, 'det')
		print(E, 'E\n')
	return E

def direct_problem(V, E):
	global step
	psi = 1000*np.ones(n+1)
	phi = np.zeros(n+1)
	phi[1] = h 	#border conditions
	counter=0
	psi_prev2=0
	while not -10**5*eps < psi[n] < 10**5*eps:
		print(psi[n], 'psi[n]')
		print(E, 'E\n')
		if psi[n] < -eps:
			E /= 1+step
		elif psi[n] > eps:
			E *= 1+step
		# E = find_E(V, E, psi)
		for i in  range(1,n):	#find phi from left
			phi[i+1] = ((2+(E-V[i])/E_max)*phi[i]-phi[i-1])
		C = 1/(sum(phi**2)*h)**0.5	#Calc C_lambda
		psi = C*phi
		counter += 1
	return E, psi, C

def reversed_problem(E_lambda, C_lambda):
	reversed_V = np.ones(n+1)
	reversed_psi = np.zeros([len(E_lambda), n+1])
	for i in range(len(E_lambda)):	#calc psi&V
		reversed_psi[i][1] = C_lambda[i]*h
	for k in range(2):
		summEPsi1=0
		for j in range(len(E_lambda)):
			summEPsi1 += E_lambda[j]*reversed_psi[j][k]**2
		reversed_V[k] = h*summEPsi1+2*E_max

	for i in range(len(E_lambda)):	#calc psi&V
		for k in range(1,n):
			summEPsi1=0
			reversed_psi[i][k+1] = (-1/E_max*(reversed_V[k]-E_lambda[i])+2)*reversed_psi[i][k]-reversed_psi[i][k-1]
			for j in range(len(E_lambda)):
				summEPsi1 += E_lambda[j]*reversed_psi[j][k+1]**2
			reversed_V[k+1] = h*summEPsi1+2*E_max
	return reversed_V, reversed_psi


##########PARAMETERS##########
m = 1.6726219*10**-27 #kilograms proton mass	1.6726219*10**-27
h_plank = 6.62607*10**-34 #[m**2*kg/s] 6.62607004*10**-34
e = 1.60217662*10**-19	 #Electron charge 1.60217662*10**-19	
R = 1*10**-9	#[m] potential pit size 10**-9
eps = 10**-2	#error value
step = 10**-5
n = 100	#number of mesh splits 
h = R/(n+0.5)
r = np.zeros(n+1)
for i in range(n+1):
	r[i] = h/2+i*h
V = (r-R/2)**2	#particle potential energy
# V = V/max(V)*e
# V[n-30:n+1], V[0:30] = 10*e, 10*e
# V[:] = 0
E_max = -h_plank**2/(2*m*h**2)
E = 5.7*10**-20
psi_list, E_lambda, psi_sqr_list = [],  [], []
C_lambda = []
print(E_max, 'E_max')
##############################

for i in range(3):
	E, psi, C = direct_problem(V, E)
	E_lambda.append(E)
	psi_list.append(psi)
	psi_sqr_list.append(psi**2)
	C_lambda.append(C)
	E *= 3

reversed_V, reversed_psi = reversed_problem(E_lambda, C_lambda)

plotter(psi_list[:3], E_lambda, r'Direct problem $\psi$')
# plotter(psi_sqr_list, E_lambda, r'Direct problem $\psi^2$'))
plotter([V/e], E_lambda, 'potential energy V [eV]')
plotter(reversed_psi[:3], E_lambda, r'Reversed problem $\psi$')
plotter([reversed_V/e], E_lambda, 'Reversed problem potential energy V [eV]')