import numpy as np
import plotly.graph_objects as go
import streamlit as st
import pandas as pd

# 1. Configuracion inicial del simulador
st.set_page_config(layout="wide")
st.title("Simulador Mision Artemis II")
st.markdown("Modelo de trayectoria de retorno libre utilizando Splines Cubicos Parametrizados.")

st.sidebar.header("Panel de Control de Vuelo")
inclinacion = st.sidebar.slider("Inclinacion Orbital (Eje Z)" , min_value=0, max_value=80, value=15)
cruce_x = st.sidebar.slider("Punto de Cruce (Eje X)", min_value=100, max_value=250, value=180)

st.sidebar.header("Tiempo de Mision")
if 'dia_anterior' not in st.session_state:
    st.session_state.dia_anterior = 0.1

dias_objetivo = st.sidebar.slider("Dias de vuelo transcurridos", min_value=0.1, max_value=10.0, value=4.0, step=0.1)
duracion_real_mision = 10.0

# 2. Parametrizacion del tiempo y espacio
t_nodos = np.linspace(0, duracion_real_mision, 5).astype(float)
x_nodos = np.array([-30, cruce_x, 394, cruce_x, -30]).astype(float)
y_nodos = np.array([30, -60, 0, 60, -30]).astype(float)
z_nodos = np.array([0, inclinacion, 0, -inclinacion, 0]).astype(float)

h_t = np.diff(t_nodos)
deltaX = np.diff(x_nodos)
deltaY = np.diff(y_nodos)
deltaZ = np.diff(z_nodos)
n = len(t_nodos)

# 3. Construccion de matrices para el sistema
def construir_Matriz_A():
    A = np.empty((n - 2, n - 2))
    for i in range(n - 2):
        A[i][i] = 2 * (h_t[i] + h_t[i + 1])
        if i < n - 3:
            A[i][i + 1] = h_t[i + 1]
        if i > 0:
            A[i][i - 1] = h_t[i]
    return A

def construir_Vector_B_x():
    B_x = np.empty(n - 2)
    for i in range(n - 2):
        B_x[i] = 6 * ((deltaX[i + 1] / h_t[i + 1]) - (deltaX[i] / h_t[i]))
    return B_x

def construir_Vector_B_y():
    B_y = np.empty(n - 2)
    for i in range(n - 2):
        B_y[i] = 6 * ((deltaY[i + 1] / h_t[i + 1]) - (deltaY[i] / h_t[i]))
    return B_y

def construir_Vector_B_z():
    B_z = np.empty(n - 2)
    for i in range(n - 2):
        B_z[i] = 6 * ((deltaZ[i + 1] / h_t[i + 1]) - (deltaZ[i] / h_t[i]))
    return B_z

A = construir_Matriz_A()
B_x = construir_Vector_B_x()
B_y = construir_Vector_B_y()
B_z = construir_Vector_B_z()

# 4. Funciones de splines y metodos numericos
def pivoteo_parcial(Mat, Vec, k, n_size, fila_mayor):
    for i in range(k + 1, n_size):
        if abs(Mat[i][k]) > abs(Mat[fila_mayor][k]):
            fila_mayor = i
    Mat[[k, fila_mayor]] = Mat[[fila_mayor, k]]
    Vec[k], Vec[fila_mayor] = Vec[fila_mayor], Vec[k]
    return Mat, Vec

def reduccion_cero(Mat, Vec, k, n_size):
    for i in range(k + 1, n_size):
        factor = Mat[i][k] / Mat[k][k]
        for j in range(k, n_size):
            Mat[i][j] -= (factor * Mat[k][j])
        Vec[i] -= (factor * Vec[k])
    return Mat, Vec

def sustitucion_atras(Mat, Vec, n_size):
    x_incognitas = np.zeros(n_size)
    for k in range(n_size - 1, -1, -1):
        suma = Vec[k]
        for i in range(k + 1, n_size):
            suma -= (Mat[k][i] * x_incognitas[i])
        x_incognitas[k] = suma / Mat[k][k]
    return x_incognitas

def eliminacion_gaussiana(Mat, Vec):
    n_size = len(Vec)
    for k in range(n_size):
        fila_mayor = k
        Mat, Vec = pivoteo_parcial(Mat, Vec, k, n_size, fila_mayor)
        Mat, Vec = reduccion_cero(Mat, Vec, k, n_size)
    return sustitucion_atras(Mat, Vec, n_size)

def calcular_coeficientes(ti, val_eje, M, h):
    num_tramos = len(ti) - 1
    coef = np.zeros((num_tramos, 4))
    for i in range(num_tramos):
        coef[i][0] = (M[i + 1] - M[i]) / (6.0 * h[i])
        coef[i][1] = M[i] / 2.0
        coef[i][2] = ((val_eje[i + 1] - val_eje[i]) / h[i]) - ((2.0 * h[i] * M[i] + h[i] * M[i + 1]) / 6.0)
        coef[i][3] = val_eje[i]
    return coef

def evaluar_spline(t_objetivo, ti, coef):
    i = 0
    for j in range(len(ti) - 1):
        if ti[j] <= t_objetivo <= ti[j+1]:
            i = j
            break
    dt = t_objetivo - ti[i]
    a, b, c, d = coef[i]
    return a*dt**3 + b*dt**2 + c*dt + d

# 5. Resolucion fisica por ejes asegurando no destruir los vectores originales
m_internos_x = eliminacion_gaussiana(np.copy(A), np.copy(B_x))
M_x_completo = np.zeros(len(x_nodos))
M_x_completo[1:-1] = m_internos_x
coef_x = calcular_coeficientes(t_nodos, x_nodos, M_x_completo, h_t)

m_internos_y = eliminacion_gaussiana(np.copy(A), np.copy(B_y))
M_y_completo = np.zeros(len(y_nodos))
M_y_completo[1:-1] = m_internos_y
coef_y = calcular_coeficientes(t_nodos, y_nodos, M_y_completo, h_t)

m_internos_z = eliminacion_gaussiana(np.copy(A), np.copy(B_z))
M_z_completo = np.zeros(len(z_nodos))
M_z_completo[1:-1] = m_internos_z
coef_z = calcular_coeficientes(t_nodos, z_nodos, M_z_completo, h_t)

# 6. Renderizado base y configuracion de cursores temporales
resolucion = 1001
t_test = np.linspace(0, duracion_real_mision, resolucion)
x_curva = np.zeros(resolucion)
y_curva = np.zeros(resolucion)
z_curva = np.zeros(resolucion)

for t in range(len(t_test)):
    x_curva[t] = evaluar_spline(t_test[t], t_nodos, coef_x)
    y_curva[t] = evaluar_spline(t_test[t], t_nodos, coef_y)
    z_curva[t] = evaluar_spline(t_test[t], t_nodos, coef_z)

idx_objetivo = int((dias_objetivo / duracion_real_mision) * (resolucion - 1))

# 7. Preparacion de contenedores vacios para interfaz en tiempo real
espacio_alertas = st.empty()
espacio_3d = st.empty()

tab1, tab2 = st.tabs(["Analisis Grafico 2D", "Cinematica y Telemetria"])

with tab1:
    col_izq, col_der = st.columns(2)
    espacio_2d_planta = col_izq.empty()
    espacio_2d_perfil = col_der.empty()

with tab2:
    espacio_metricas = st.empty()
    st.markdown("---")
    espacio_coord_actual = st.empty()
    espacio_tabla = st.empty()
    espacio_descarga = st.empty()

# 8. Funcion de cinematica
def calcular_telemetria_actual(x, y, z, t):
    if len(x) < 2:
        return 0.0, 0.0, 0.0, 0.0, 0.0
        
    v_x = np.gradient(x, t)
    v_y = np.gradient(y, t)
    v_z = np.gradient(z, t)
    
    v_mag_dia = np.sqrt(v_x[-1]**2 + v_y[-1]**2 + v_z[-1]**2)
    velocidad_kmh = v_mag_dia * (1000.0 / 24.0)

    if len(v_x) > 1:
        v_mag_ant = np.sqrt(v_x[-2]**2 + v_y[-2]**2 + v_z[-2]**2)
        delta_v = velocidad_kmh - (v_mag_ant * (1000.0 / 24.0))
    else:
        delta_v = 0.0

    a_x = np.gradient(v_x, t)
    a_y = np.gradient(v_y, t)
    a_z = np.gradient(v_z, t)
    a_mag_dia2 = np.sqrt(a_x[-1]**2 + a_y[-1]**2 + a_z[-1]**2)
    
    aceleracion_kmh2 = a_mag_dia2 * (1000.0 / 576.0)
    fuerzas_g = (a_mag_dia2 * (1000000.0 / (86400.0**2))) / 9.81
    
    distancias_seg = np.sqrt(np.diff(x)**2 + np.diff(y)**2 + np.diff(z)**2)
    distancia_total_km = np.sum(distancias_seg) * 1000.0

    return velocidad_kmh, delta_v, aceleracion_kmh2, fuerzas_g, distancia_total_km

# 9. Motor de renderizado volumetrico fotorrealista
def crear_planeta_renderizado(radio, cx, cy, cz, tipo_cscale):
    # Alta densidad de malla (50x50) para superficie redonda
    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 50)
    x = radio * np.outer(np.cos(u), np.sin(v)) + cx
    y = radio * np.outer(np.sin(u), np.sin(v)) + cy
    z = radio * np.outer(np.ones(np.size(u)), np.cos(v)) + cz
    
    # Generacion de texturas procedimentales para simular accidentes geograficos
    if tipo_cscale == 'tierra':
        surf_color = np.outer(np.sin(u * 4), np.cos(v * 3))
        cscale = [[0.0, '#0d2b45'], [0.3, '#1f77b4'], [0.5, '#2ca02c'], [1.0, '#1b5e20']]
    else:
        surf_color = np.outer(np.cos(u * 6), np.sin(v * 6))
        cscale = [[0.0, '#3a3a3a'], [0.5, '#7f7f7f'], [1.0, '#d3d3d3']]

    return go.Surface(
        x=x, y=y, z=z, 
        surfacecolor=surf_color,
        colorscale=cscale, 
        showscale=False, 
        hoverinfo='skip',
        # Parametros avanzados de iluminacion de hardware WebGL
        lighting=dict(ambient=0.6, diffuse=0.8, roughness=0.4, specular=0.6, fresnel=0.4)
    )

i_safe = max(2, idx_objetivo)

# Grafica 3D
fig3d = go.Figure()
fig3d.add_trace(crear_planeta_renderizado(15, 0, 0, 0, 'tierra'))   # Tierra Estilizada
fig3d.add_trace(crear_planeta_renderizado(5, 384, 0, 0, 'luna'))    # Luna Rugosa
fig3d.add_trace(go.Scatter3d(x=x_curva[:i_safe], y=y_curva[:i_safe], z=z_curva[:i_safe], mode='lines', line=dict(color='#ff4b4b', width=5), name='Artemis II'))
fig3d.add_trace(go.Scatter3d(x=[x_curva[i_safe-1]], y=[y_curva[i_safe-1]], z=[z_curva[i_safe-1]], mode='markers', marker=dict(size=8, color='#ffd700', symbol='diamond'), name='Orion'))

fig3d.update_layout(
    title="Visor Espacial Principal 3D", 
    scene=dict(
        aspectmode="data",
        xaxis=dict(showbackground=False, gridcolor='rgba(128, 128, 128, 0.2)', zerolinecolor='rgba(128, 128, 128, 0.4)'),
        yaxis=dict(showbackground=False, gridcolor='rgba(128, 128, 128, 0.2)', zerolinecolor='rgba(128, 128, 128, 0.4)'),
        zaxis=dict(showbackground=False, gridcolor='rgba(128, 128, 128, 0.2)', zerolinecolor='rgba(128, 128, 128, 0.4)'),
        bgcolor='rgba(0,0,0,0)' # Fondo transparente
    ), 
    paper_bgcolor='rgba(0,0,0,0)',
    plot_bgcolor='rgba(0,0,0,0)',
    height=600, 
    margin=dict(l=0, r=0, b=0, t=40)
)
espacio_3d.plotly_chart(fig3d, use_container_width=True)

# Grafica 2D Planta
fig_planta = go.Figure()
fig_planta.add_trace(go.Scatter(x=[0], y=[0], mode='markers', marker=dict(size=25, color='#1f77b4'), name='Tierra'))
fig_planta.add_trace(go.Scatter(x=[384], y=[0], mode='markers', marker=dict(size=15, color='#7f7f7f'), name='Luna'))
fig_planta.add_trace(go.Scatter(x=x_curva[:i_safe], y=y_curva[:i_safe], mode='lines', line=dict(color='#ff4b4b', width=3), name='Artemis II'))
fig_planta.add_trace(go.Scatter(x=[x_curva[i_safe-1]], y=[y_curva[i_safe-1]], mode='markers', marker=dict(size=10, color='#ffd700', symbol='diamond'), name='Orion'))
fig_planta.update_layout(
    title="Plano X-Y (Planta)", height=400,
    paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)',
    xaxis=dict(gridcolor='rgba(128, 128, 128, 0.2)'), yaxis=dict(gridcolor='rgba(128, 128, 128, 0.2)')
)
espacio_2d_planta.plotly_chart(fig_planta, use_container_width=True)

# Grafica 2D Perfil (Fondo Adaptativo)
fig_perfil = go.Figure()
fig_perfil.add_trace(go.Scatter(x=[0], y=[0], mode='markers', marker=dict(size=25, color='#1f77b4'), name='Tierra'))
fig_perfil.add_trace(go.Scatter(x=[384], y=[0], mode='markers', marker=dict(size=15, color='#7f7f7f'), name='Luna'))
fig_perfil.add_trace(go.Scatter(x=x_curva[:i_safe], y=z_curva[:i_safe], mode='lines', line=dict(color='#ff4b4b', width=3), name='Artemis II'))
fig_perfil.add_trace(go.Scatter(x=[x_curva[i_safe-1]], y=[z_curva[i_safe-1]], mode='markers', marker=dict(size=10, color='#ffd700', symbol='diamond'), name='Orion'))
fig_perfil.update_layout(
    title="Plano X-Z (Perfil Orbital)", height=400,
    paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)',
    xaxis=dict(gridcolor='rgba(128, 128, 128, 0.2)'), yaxis=dict(gridcolor='rgba(128, 128, 128, 0.2)')
)
espacio_2d_perfil.plotly_chart(fig_perfil, use_container_width=True)

# Telemetria Matematica Reactiva
vel, delta, acel, gs, dist = calcular_telemetria_actual(x_curva[:i_safe], y_curva[:i_safe], z_curva[:i_safe], t_test[:i_safe])

with espacio_metricas.container():
    st.markdown("### Detalles de Vuelo en Tiempo Real")
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Velocidad Relativa", f"{np.round(vel):,} km/h", delta=f"{np.round(delta, 2)} km/h")
    c2.metric("Aceleracion", f"{np.round(acel):,} km/h2")
    c3.metric("Fuerza G Aparente", f"{np.round(gs, 4)} G")
    c4.metric("Odometro", f"{np.round(dist):,} km")

# 10. Actualizacion final de estatus y tabla de datos
distancia_luna = np.sqrt((x_curva[idx_objetivo] - 384)**2 + y_curva[idx_objetivo]**2 + z_curva[idx_objetivo]**2)

with espacio_alertas.container():
    st.subheader("Estado de Navegación")
    if dias_objetivo >= 10:
        st.success("Acuatizaje confirmado. La tripulacion ha regresado a salvo a la Tierra.")
    elif distancia_luna < 1.73:
        st.error(f"Impacto lunar. Distancia: {distancia_luna:.2f} unidades.")
    elif distancia_luna <= 15:
        st.warning(f"Maximo acercamiento (Perilunio). Orbita lunar a {distancia_luna:.2f} unidades.")
    else:
        st.info(f"Motores encendidos. Navegando por el espacio profundo. Distancia a la Luna: {distancia_luna:.2f} unidades.")

with espacio_coord_actual.container():
    st.markdown("### Coordenada Actual del Vehículo")
    col_x, col_y, col_z = st.columns(3)
    col_x.metric("Posición X", f"{x_curva[idx_objetivo] * 1000:,.2f} km")
    col_y.metric("Posición Y", f"{y_curva[idx_objetivo] * 1000:,.2f} km")
    col_z.metric("Posición Z", f"{z_curva[idx_objetivo] * 1000:,.2f} km")
    st.markdown("---")

with espacio_tabla.container():
    st.markdown("### Sistema de Coordenadas")
    df_telemetria = pd.DataFrame({
        "Tiempo (Dias)": np.round(t_test[10:idx_objetivo + 1], 2),
        "Coordenada X (km)": np.round(x_curva[10:idx_objetivo + 1] * 1000, 2),
        "Coordenada Y (km)": np.round(y_curva[10:idx_objetivo + 1] * 1000, 2),
        "Coordenada Z (km)": np.round(z_curva[10:idx_objetivo + 1] * 1000, 2)
    })
    st.dataframe(df_telemetria, use_container_width=True)
    
with espacio_descarga.container():
    csv_data = df_telemetria.to_csv(index=False).encode('utf-8')
    st.download_button(
        label="Descargar Telemetria (CSV)",
        data=csv_data,
        file_name="artemisII_telemetria.csv",
        mime="text/csv"
    )

#11. Muestra del motor matemático interno, en específico el sistema de ecuaciones
st.markdown("---")
with st.expander("Ver Motor Matemático Interno"):
    st.markdown("El sistema utiliza Eliminación Gaussiana con Pivoteo Parcial para resolver el sistema de ecuaciones lineales de los Splines Cúbicos. Se asegura continuidad en las segundas derivadas.")
    col_mat1, col_mat2 = st.columns(2)
    with col_mat1:
        st.write("Matriz Tridiagonal A (Coeficientes de tiempo):", A)
    with col_mat2:
        st.write("Vector Resultante B (Eje X):", B_x)
        st.write("Vector Resultante B (Eje Y):", B_y)

st.session_state.dia_anterior = dias_objetivo