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
F_vect=np.zeros((2,Nx+1))

tmax=5
C=0.9

#Données pour l'enregistrement de la video
duration = tmax*10
dt_video=0.01
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
	return F_etoile

def Plot_water_high(X,W,temps):
	plt.clf()
	plt.xlabel('Position x')
	plt.ylabel(("Hauteur d'eau h"))
	plt.xlim(-1.1,1.1)
	plt.ylim(-0.1,2.1)
	plt.title("Hauteur d'eau h en fonction de la position x")
	plt.text(0.5,1.8,"Temps : "+str(round(temps,2)))
	plt.plot(X,W[0,1:Nx+1])
	path_fig=image_folder+'/Saint_Venant_t_'+str(round(temps,2))+'.png'
	plt.savefig(path_fig)
	image_path_list.append(path_fig)
	plt.pause(0.001) #Pour afficher la hauteur d'eau durant le calcul

def condition_cfl(W):
	cfl=C*dx/(2*max(abs(Lambda_R(W,1)),abs(Lambda_L(W,0))))
	for i in range(0,Nx+1):
		cfl=min(cfl,C*dx/(2*max(abs(Lambda_R(W,i+1)),abs(Lambda_L(W,i)))))
	return(cfl)

def conservation(temps,W): #Ecrit le volume d'eau, l'energie potentiel, cinétique et mécanique dans un fichier pour vérifier leurs bonne conservation
	volume_eau=0
	energie_potentiel=0
	enegie_cinetique=0
	energie_mecanique=0
	for i in range(1,Nx+1):
		volume_eau=volume_eau+W[0,i]*dx
		energie_potentiel=energie_potentiel+g*W[0,i]*dx
		enegie_cinetique=0.5*((W[1,i]/W[0,i])**2)*dx
	energie_mecanique=energie_mecanique+energie_potentiel+enegie_cinetique
	conservation_file.write(str(temps)+" "+str(volume_eau)+" "+str(energie_potentiel)+" "+str(enegie_cinetique)+" "+str(energie_mecanique)+"\n")



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

conservation_file=open("Conservation.txt", "w")
conservation(temps,W)


dt=condition_cfl(W)
temps=temps+dt

while temps<tmax:


	for i in range(0,Nx+1):
		F_vect[:,i]=F_etoile(W,i)

	for i in range(0,Nx+1):
		W[:,i]=W[:,i]-dt*F_vect[:,i]/dx
		W[:,i+1]=W[:,i+1]+dt*F_vect[:,i]/dx


	#CONDITION AUX LIMITES
	W[0,0]=W[0,1]
	W[0,Nx+1]=W[0,Nx]

	W[1,0]=-W[1,1]
	W[1,Nx+1]=-W[1,Nx]


	if temps>seuil_temps_video:
		seuil_temps_video=seuil_temps_video+dt_video
		Plot_water_high(X,W,temps)
		print(temps)
		conservation(temps,W)

	dt=condition_cfl(W)
	temps=temps+dt


Plot_water_high(X,W,temps)


clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_path_list, fps=fps)
clip.write_videofile(r"Saint_Venant_Animation_CL3.mp4", codec="libx264")

conservation_file.close()