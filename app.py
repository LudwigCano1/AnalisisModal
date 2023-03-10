import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px

lista_m = []
lista_k = []

st.title("Análisis Modal")

# Ingresar la cantidad de grados de libertad
st.subheader("Grados de libertad")

n = st.number_input("Número de grados de libertad:",min_value=2,value=2,label_visibility="collapsed")

c1,c2,c3 = st.columns([1,3,3])
with c2: st.subheader("Masas [t]")
with c3: st.subheader("Rigideces [t/m]")

for i in range(n):
    c1,c2,c3 = st.columns([1,3,3])
    with c1: f"Piso {i+1}"
    with c2: lista_m.append(st.number_input(f"Masa {i+1}",min_value=0.0,label_visibility="collapsed"))
    with c3: lista_k.append(st.number_input(f"Rigidez {i+1}",min_value=0.0,label_visibility="collapsed"))

        

# Crear matriz de masas [M]
M = np.zeros(shape=(n,n)) # Matriz de ceros
for i in range(0,n,1):
    M[i,i] = lista_m[i]

# Crear matriz de rigideces [K]
K = np.zeros(shape=(n,n))
for i in range(0,n,1):
    if i == n-1:
        K[i,i] = lista_k[i]
    else:
        K[i,i] = lista_k[i] + lista_k[i+1]
for i in range(1,n,1):
    K[i-1,i] = -lista_k[i]
    K[i,i-1] = -lista_k[i]

# Matriz de rigideces transformada
nK = np.zeros(shape=(n,n))
for i in range(0,n,1):
    nK[i,:] = K[i,:]/M[i,i]

a,x = np.linalg.eig(nK) # a = "Valores propios", b = "Vectores propios"
w = np.sqrt(a)
T = 2 * np.pi / w

# Normalización
for i in range(0,n,1):
    x[:,i] = x[:,i] / max(abs(x[:,i]))

# Normalización con respecto a la masa
Phi = np.zeros(shape=(n,n))
for i in range(n):
    Phi[:,i] = x[:,i] / np.sqrt(x[:,i] @ M @ x[:,i])

# Factores de participación
I = np.ones(shape=(n,1))
r = []
for i in range(n):
    r.append(((x[:,i] @ M @ I)/(x[:,i] @ M @ x[:,i]))[0])

# Masas participativas
mp = []
for i in range(n):
    mp.append(100 * (((I.T @ M @ x[:,i])**2)/(x[:,i] @ M @ x[:,i]))[0] / sum(lista_m))

fig,axes = plt.subplots(nrows=1,ncols=n,figsize=(10,6))

for f in range(1,n+1,1):
    d = x[:,n-f]
    d = np.hstack((np.array([0]),d))
    y = [altura for altura in range(0,n+1,1)]
    axes[f-1].plot(d,y,marker='o',markersize=12,linestyle="--")
    axes[f-1].vlines(0,0,y[n],colors="r")
    axes[f-1].set_title(f"Forma {f} \n T = {T[n-f]:.03} s \n M = {mp[n-f]:.1f} %")
    axes[f-1].set_xlim(-1,1)
    axes[f-1].grid()

tabla = pd.DataFrame()
tabla["Modos"] = [f"Modo {x}" for x in range(1,n+1,1)]
tabla["T (s)"] = T[::-1]
tabla["M_ef (%)"] = mp[::-1]

mpa = [100.0]
for i in range(n-1):
    mpa.append(mpa[-1]-mp[i])

tabla["M_ef_acum (%)"] = mpa[::-1]

# Animación
fig2 = px.line(x=[0],y=[0])
fig2.data = []
fig2.update_layout(xaxis_range=[-1.1,1.1])


with st.sidebar: f = st.number_input("Modo:",min_value=1,max_value=n)

d = x[:,n-f]
d = np.hstack((np.array([0]),d))

for step in np.linspace(0, 20*T[n-f],50):
    fig2.add_scatter(
            visible=False, marker={"size":15},
            line=dict(color="#00CED1", width=4),
            x=d*np.sin(w[n-f]*step/10),
            y=[altura for altura in range(0,n+1,1)])

# Make 10th trace visible
fig2.data[0].visible = True

# Create and add slider
steps = []
for i in range(len(fig2.data)):
    step = dict(
        method="update",
        args=[{"visible": [False] * len(fig2.data)}],  # layout attribute
    )
    step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
    steps.append(step)

sliders = [dict(
    active=0,
    currentvalue={"prefix": "Frequency: "},
    pad={"t": 50},
    steps=steps
)]

fig2.update_layout(
    sliders=sliders,
    width=350,height=150*n,
    margin=dict(l=10,r=10,t=10,b=10)
)

with st.sidebar:
    st.plotly_chart(fig2,theme=None)

p1,p2,p3,p4 = st.tabs(["[ Formas de modo ]","[ Matriz de Masas ]","[ Matriz de Rigideces ]","[ Tabla Resumen ]"])

with p1: st.pyplot(fig)
with p2: M
with p3: K
with p4: tabla
