import numpy as np
import matplotlib.pyplot as plt
import math
import os, shutil
import moviepy.video.io.ImageSequenceClip


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

def condition_cfl(W):
	cfl=C*dx/(2*max(abs(Lambda_R(W,1)),abs(Lambda_L(W,0))))
	for i in range(0,Nx+1):
		cfl=min(cfl,C*dx/(2*max(abs(Lambda_R(W,i+1)),abs(Lambda_L(W,i)))))
	return(cfl)

def Plot_Water_High(X,W):
	plt.clf()
	plt.xlabel('Position x')
	plt.ylabel(("Hauteur d'eau h"))
	plt.xlim(-1.1,1.1)
	plt.ylim(-0.1,3.1)
	plt.title("Hauteur d'eau h en fonction de la position x")
	plt.text(0,0.5,"Temps : "+str(round(temps,6)))
	plt.text(0,0.2,"dt : "+str(round(dt,16)))
	plt.plot(X,W[0,1:Nx+1])
	plt.pause(0.001)

def Save_Water_High_Plot(X,W,temps,k):
	plt.clf()
	plt.xlabel('Position x')
	plt.ylabel(("Hauteur d'eau h"))
	plt.xlim(-1.1,1.1)
	plt.ylim(-0.1,3.1)
	plt.title("Hauteur d'eau h en fonction de la position x")
	plt.text(0,0.5,"Temps : "+str(round(temps,6)))
	plt.text(0,0.2,"dt : "+str(round(dt,16)))
	plt.plot(X,W[0,1:Nx+1])
	path_fig=image_folder+'/Saint_Venant_t_'+str(k)+'.png'
	plt.savefig(path_fig)
	image_path_list.append(path_fig)
	plt.pause(0.001) #Pour afficher la hauteur d'eau durant le calcul



#PARAMETRES
tmax=5
C=0.4

g=9.81
Nx=100
x_L=-1
x_R=1
H_L=2
H_R=1
dx=(x_R-x_L)/Nx
x_choc=0

Ny=10



#DEFINITION VARIABLES

#Les cellules 0 et Nx+1 sont des cellules fantomes qui servent de CL.
F_vect=np.zeros((2,Nx+1,Ny))	
W=np.zeros((Ny+1,Nx+2))			#W(0,i)=h_xi, W(j,i)=h_xi*u_xi_yj
W_inter=np.zeros((2,Nx+2,Ny))	#W_inter(0,i,j)=h_xi_yj, W_inter(1,xi,yj)=h_xi_yj*u_xi_yj
G=np.zeros((Ny-1,Nx+2))
zeros=np.zeros(Nx+2)

X=np.zeros(Nx)
for i in range(0,Nx):
	X[i]=x_L+i*dx

L=np.zeros(Ny)
Tot=0
for i in range(0,Ny-1):
	L[i]=1/Ny
	Tot=Tot+L[i]
L[Ny-1]=1-Tot



#CONDITIONS INITALES

indice_choc=int((x_choc-x_L)/dx)
for i in range(0,Nx+2):
	if i<=indice_choc :
		W[0,i]=H_L
	else:
		W[0,i]=H_R



#Données pour l'enregistrement de la video
fig, ax = plt.subplots()
image_path_list=[]

# image_folder="/home/anatole/Images/Saint_Venant_Animation_n_couches" #Attention ce dossier est créé puis détruit à chaque execution
# if os.path.exists(image_folder):
# 	shutil.rmtree(image_folder)
# os.mkdir(image_folder)



#RESOLUTION
temps=0
dt=1000
k=1
while temps<tmax and dt>pow(10,-15) and k<500:
	W_inter[:,:,:]=0
	dt=1000

	#CREATION DE W_inter (regroupement des Ny Saint Venant Classiques), CALCULS DES FLUX ET DU PAS DE TEMPS
	for j in range(0,Ny):
		W_inter[0,:,j]=L[j]*W[0,:]
		W_inter[1,:,j]=L[j]*W[j+1,:]
		for i in range(0,Nx+1):
			F_vect[:,i,j]=F_etoile(W_inter[:,:,j],i)

		dt_inter=condition_cfl(W_inter[:,:,j])
		if dt>dt_inter:
			dt=dt_inter
	print(temps,' / ',dt)


	#MISE A JOUR DES Ny Saint Venant Classique dans W_inter
	for j in range(0,Ny):
		for i in range(0,Nx+1):
			W_inter[:,i,j]=W_inter[:,i,j]-dt*F_vect[:,i,j]/dx
			W_inter[:,i+1,j]=W_inter[:,i+1,j]+dt*F_vect[:,i,j]/dx

	#MISE A JOUR DE h^n+1 DANS W
	W[:,:]=0
	for j in range(0,Ny):
		W[0,:]=W[0,:]+W_inter[0,:,j]


	#MISE A JOUR de h^n+1*u_j DANS W
	G[0,:]=L[0]*(W[0,:]-W_inter[0,:,0])/dt
	W[1,:]=W_inter[1,:,0]+dt*(W_inter[1,:,1]*np.maximum(zeros,G[0,:])+W_inter[1,:,0]*np.minimum(zeros,G[0,:]))/L[0]
	for j in range(1,Ny-1):
		G[j,:]=G[j-1,:]+L[j]*(W[0,:]-W_inter[0,:,j])/dt
		W[j+1,:]=W_inter[1,:,j]+dt*(W_inter[1,:,j+1]*np.maximum(zeros,G[j,:])+W_inter[1,:,j]*np.minimum(zeros,G[j,:])-W_inter[1,:,j]*np.maximum(zeros,G[j-1,:])-W_inter[1,:,j-1]*np.minimum(zeros,G[j-1,:]))/L[j]


	#CONDITION AUX LIMITES ANTISYMETRIQUES
	W[0,0]=W[0,1]
	W[0,Nx+1]=W[0,Nx]

	for j in range(1,Ny):
		W[j,0]=-W[j,1]
		W[j,Nx+1]=-W[j,Nx]


	
	#SAUVEGARDE DU PLOT DANS LE DOSSIER "image folder"
	#Plot_water_high(X,W,temps,k)
	k=k+1

	#JUSTE PLOT
	Plot_Water_High(X,W)


	temps=temps+dt


#CREATION DE LA VIDEO
# duration = 10
# fps=len(image_path_list)/duration
# clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_path_list, fps=fps)
# clip.write_videofile(r"Saint_Venant_Animation_"+str(Ny)+"_couches.mp4", codec="libx264")