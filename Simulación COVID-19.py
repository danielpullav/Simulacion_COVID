# =============================================================================
# Importacion de las librerias
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import odeint

# =============================================================================
# Funcion para una mejor visualizacion de la gráfica
# =============================================================================
%matplotlib qt  
pob_azuay = round(712127 + ((712127)/(14483499))*((17268000)-(14483499)), 2) #Poblacion estimada del Azuay

def mod_sir(y=0 , t=0 , N=1 , b=0 , k=0):
    """
    Función que nos permite expresar ecuaciones diferenciales (ecuaiones SIR) para 
    ser utilizadas de una manera más facil en el calculo de las funciones a graficar
    Parameters
    ----------
    y : array
        Valores en y de las 3 funciones a graficar.
    t : array
        Eje x correspondiente al tiempo en dias.
    N : int or float
        Porcentaje poblacional a considerar en la grafica, por defecto el 100%.
    b : float
        Tasa de infeccion o poblacion infectada.
    k : float
        Tasa de recuperacion o poblacion recuperada.

    Returns
    -------
    dsdt : float
        Tasa de cambio de la proporción poblacional susceptible con respecto al tiempo.
    drdt : float
        Tasa de cambio de la proporción poblacional recuperada o fallecida con respecto al tiempo.
    didt : float
        Tasa de cambio de la proporción poblacional infectada con respecto al tiempo.

    """
    s_t, r_t, i_t = y
    dsdt = -b*s_t* i_t/N
    drdt =  k*i_t 
    didt = (b*s_t*i_t/N)-(k*i_t)  
    return dsdt, drdt, didt

t = np.linspace(0, 140, 140,endpoint=True) #Valores del eje x correspondiente a los dias
y_0 = np.array([1 , 0 , 1/pob_azuay]) #Valores iniciales de las 3 funciones

# =============================================================================
# Ingreso y validación del parámetro b por el usuario
# =============================================================================
b=float(input("Ingrese la tasa de infección (b):\t"))
while b<0 or b>1:
    b=float(input("Valor inválido. Reingrese la tasa de infección (b):\t"))
    

# =============================================================================
# Ingreso y validación del parámetro k por el usuario
# =============================================================================
k=float(input("Ingrese la tasa de recuperación (k):\t"))
while k<0 or k>1:
    k=float(input("Valor inválido. Reingrese la tasa de recuperación (k):\t"))
 
N_pob = 1 #Variable que especifica el porcentaje de la poblacion a considerar en la grafica
SIR_data = odeint(mod_sir, y_0, t, args=(N_pob, b , k)) #Uso de la funcion odeint para el calculo de las funciones en 
                                                     #las ecuaciones diferenciale a graficar

# =============================================================================
# Espacio donde se va a ilustrar la grafica
# =============================================================================
fig = plt.figure(figsize=(30,15),facecolor="ivory")
graph_ax = plt.axes((0.1,0.0875,0.825,0.8),facecolor="lightgoldenrodyellow")

def graph(frames):
    """
    Función que permite la graficación de las 3 funciones cuadro por cuadro en un mismo espacio
    para que se pueda apreciar su correspondencia y comparativa 
    Parameters
    ----------
    frames : range object
        Fuente de datos para pasar la función y cada cuadro de animación.

    Returns
    -------
    None.

    """
    graph_ax.clear()
    graph_ax.plot(t[:frames],SIR_data[:,0][:frames],color="r",label="Susceptible",linewidth=9) #plot de proporción susceptible
    graph_ax.plot(t[:frames],SIR_data[:,1][:frames], 'X',color="b",label="Recuperada",linewidth=9) #plot de proporción recuperada
    graph_ax.plot(t[:frames],SIR_data[:,2][:frames],'o--',color="g",label="Infectada",linewidth=9) #plot de proporción infectada
    plt.suptitle("Modelo SIR para la provincia del Azuay\nb={0}   k={1}".format(b,k),fontsize=35,fontweight="bold") #título y presentación de los valores de los parámetros
    graph_ax.set_xlim(0, 140) #límite de la gráfica en el eje x
    graph_ax.set_ylim(-0.01, 1.01)  #límite de la gráfica en el eje y
    graph_ax.set_ylabel("Proporción Poblacional",fontsize=25,fontweight="bold") #etiqueta en el eje y
    graph_ax.set_xlabel("Tiempo en días",fontsize=25,fontweight="bold") #etiqueta en el eje x
    graph_ax.tick_params(labelsize=25) #tamaño de los ticks en los ejes
    graph_ax.legend(loc="upper right",fontsize=25,facecolor="ivory") #leyenda
    plt.grid() #celdas

# =============================================================================
# Animacion de la grafica 
# =============================================================================
animator=animation.FuncAnimation(fig,graph,frames=range(len(t)),repeat=True,interval=0.00001)

# =============================================================================
# Guardado de la grafica como gif
# =============================================================================
#animator.save(r'SimulaciónCOVID.gif')  
    

