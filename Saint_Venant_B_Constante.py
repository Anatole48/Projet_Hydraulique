import numpy as np
import matplotlib.pyplot as plt
import math
import os, shutil
import moviepy.video.io.ImageSequenceClip

#PARAMETRES

g=9.81

Nx=100
H_L=2
H_R=1
x_L=-1
x_R=1
x_choc=0
dx=(x_R-x_L)/Nx
indice_choc=int((x_choc-x_L)/dx)

#Les cellules 0 et Nx+1 sont des cellules fantomes qui servent de CL.
W=np.zeros((2,Nx+2))
X=np.zeros(Nx)

tmax=50
C=0.9

duration = tmax
dt_video=tmax/100
fps=tmax/(duration*dt_video)
fig, ax = plt.subplots()

#FONCTIONS

def Lambda_L(W,i):
	Lambda_L=min(W[1,i]/W[0,i]-math.sqrt(g*W[0,i]),0)
	return Lambda_L

def Lambda_R(W,i):
	Lambda_R=max(W[1,i]/W[0,i]+math.sqrt(g*W[0,i]),0)
	return Lambda_R

def F(W,i):
	F=np.zeros(2)
	F[0]=W[1,i]
	F[1]=(W[1,i]**2)/W[0,i]+0.5*g*(W[0,i]**2)
	return F


# A Gauche L : i
# A Droite R : i+1

def F_etoile(W,i):
	F_etoile=np.zeros(2)
	F_etoile=(Lambda_R(W,i+1)*F(W,i)-Lambda_L(W,i)*F(W,i+1)+Lambda_L(W,i)*Lambda_R(W,i+1)*(W[:,i+1]-W[:,i]))/(Lambda_R(W,i+1)-Lambda_L(W,i))
	#print('Lambda : ',(Lambda_R(W,i+1)-Lambda_L(W,i)))
	return F_etoile

def Plot_water_high(X,W,temps):
	plt.clf()
	plt.xlabel('Position x')
	plt.ylabel(("Hauteur d'eau h"))
	plt.xlim(-1.1,1.1)
	plt.ylim(-0.1,2.1)
	plt.title("Hauteur d'eau h en fonction de la position x")
	plt.text(0.5,1.8,"Temps : "+str(round(temps,1)))
	plt.plot(X,W[0,1:Nx+1])
	path_fig=image_folder+'/Saint_Venant_t_'+str(round(temps,1))+'.png'
	plt.savefig(path_fig)
	image_path_list.append(path_fig)


#CONDITION INITAL
for i in range(0,Nx):
	X[i]=x_L+i*dx

for i in range(0,Nx+2):
	if i<=indice_choc :
		W[0,i]=H_L
		W[1,i]=0
	else:
		W[0,i]=H_R
		W[1,i]=0
	W[1,i]=0



#RESOLUTION

temps=0
seuil_temps_video=dt_video
image_path_list=[]

# image_folder="/home/anatole/Images/Saint_Venant_Animation" #Attention ce dossier est créé puis détruit à chaque execution
# if os.path.exists(image_folder):
# 	shutil.rmtree(image_folder)
# os.mkdir(image_folder)

Plot_water_high(X,W,temps)
plt.pause(0.1)

dt=C*dx/(2*max(abs(Lambda_R(W,1)),abs(Lambda_L(W,0))))
for i in range(0,Nx+1):
	dt=min(dt,C*dx/(2*max(abs(Lambda_R(W,i+1)),abs(Lambda_L(W,i)))))
temps=dt


while temps<tmax:
	W_inter=W

	dt=C*dx/(2*max(abs(Lambda_R(W,1)),abs(Lambda_L(W,0))))
	for i in range(0,Nx+1):
		W_inter[:,i]=W[:,i]-dt*F_etoile(W,i)
		W_inter[:,i+1]=W[:,i+1]+dt*F_etoile(W,i)
		dt=min(dt,C*dx/(2*max(abs(Lambda_R(W,i+1)),abs(Lambda_L(W,i)))))


	#CONDITION AUX LIMITES
	W_inter[0,0]=W_inter[0,1]
	W_inter[0,Nx+1]=W_inter[0,Nx]

	W_inter[1,0]=W_inter[1,1]
	W_inter[1,Nx+1]=W_inter[1,Nx]


	W=W_inter


	if temps>seuil_temps_video:
		seuil_temps_video=seuil_temps_video+dt_video
		#make_frame(temps)
		Plot_water_high(X,W,temps)
		print('Je plot !')
		plt.pause(0.001)
		print(temps)

	temps=temps+dt

Plot_water_high(X,W,temps)
plt.pause(0.001)


clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_path_list, fps=fps)
clip.write_videofile(r"Saint_Venant_Animation.mp4",fps=30, codec="libx264")