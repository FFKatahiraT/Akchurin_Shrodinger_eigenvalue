import matplotlib.pyplot as plt
import numpy as np

def plotter(Func, name):
	linestyle=['solid', 'dotted', 'dashed', 'dashdot', (0, (5,10)), (0, (5,1)), (0, (3,10,1,10)), (0, (3,1,1,1))]
	i=0
	for i in range(len(Func)):
		plt.plot(r, Func[i], label=name +'. E='+str(round(E/e,2))+'[eV]', linestyle=linestyle[i])
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

def multiply_diagonal(diagonal, psi, E):
	multiple = 1
	for i in range(len(diagonal)):
		multiple *= diagonal[i]*psi[i]-E
	return multiple

def psi_calc():
	psi = np.ones(n+1)	#psi init value
	summ = sum(psi**2)*h 	#calc integral
	psi = np.sqrt(psi**2/summ)
	E = 2*10**-20	#Energy init value
	summ = 0
	#########DIAGONALS#########
	A = E_max*np.ones(n+1)	#A*psi(i-1)+-C*psi(i)+B*psi(i+1) = -F
	A[0] = 0
	C = E_max*2-V 	
	B = E_max*np.ones(n+1)
	B[n] = 0
	###########################
	while not 1-eps < summ < 1+eps:
		det = 1
		while not -eps < det < eps:	#find E
			if det < -eps or summ > 1+eps:
				E /= 1+10**-3
			elif det > eps or summ < 1-eps:
				E *= 1+10**-3
			det = multiply_diagonal(A, psi, E) - multiply_diagonal(C, psi, E) + multiply_diagonal(B, psi, E)
			print(det, 'det')
			print(E, 'E\n')

		F =  -E*psi
		#CONSIDER y[i]=alpha[i]*y[i-1]+beta[i]
		#LEFT BORDER CONDITION
		LB_Condition = 10**-60
		#RIGHT BORDER CONDITION 
		RB_Condition_alpha = 10**-60
		RB_Condition_beta = 10**-60	#RB_Condition_beta IS NOT 0 TO AVOID DIVISION BY ZERO ERROR
		psi = progonka_start_left(A, B, C, F, RB_Condition_alpha, RB_Condition_beta, LB_Condition)
		summ = sum(psi**2)*h 	#calc integral
		print(summ, ' summ')
		psi = np.sqrt(psi**2/summ)	#bring value to the 1

	return E, psi


##########PARAMETERS##########
m = 1.6726219*10**-27 #kilograms proton mass	1.6726219*10**-27
h_plank = 6.62607*10**-34 #[m**2*kg/s] 6.62607004*10**-34
e = 1.60217662*10**-19	 #Electron charge 1.60217662*10**-19	
R = 1*10**-9	#[m] potential pit size 10**-9
eps = 0.01	#error value
n = 100
h = R/(n+0.5)
r = np.zeros(n+1)
for i in range(n+1):
	r[i] = h/2+i*h
V = (r-R/2)**2	#particle potential energy
V = V/max(V)*e
V[n-30:n+1], V[0:30] = 0, 0
# V[:] = 0
E_max = -h_plank**2/(2*m*h**2)
print(E_max, 'E_max')
##############################

E, psi = psi_calc()

plotter([psi], r'$\psi$')
plotter([psi**2], r'$\psi^2$')
plotter([V/e], 'potential energy V [eV]')