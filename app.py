import streamlit as st
import numpy as np
import sympy as sp
from scipy.optimize import root # para sistema Ecuaciones no lineales.

st.set_page_config(page_title="Propagación de Incertidumbre", layout="wide")
st.title("🧮 Análisis de Incertidumbre en Fórmulas")

st.markdown("Calculo de la incertidumbre en la concentración de Co")

# definición de funciones
def Aesp(Cn_i, w_i,lam,tr,td,ti,tv):# Calcula la actividad específica
  C_i = cal_C(lam, tr)
  D_i = cal_D(lam, td)
  H_i = cal_H(tr, tv)
  S_i = cal_S(lam, ti)
  return Cn_i*D_i*C_i*H_i/(S_i*w_i) 
def cal_D(lam, td):# Calcula D del elemento i
  return np.exp(lam*td)
def cal_C(lam, tr):# Calcula C del elemento i
  return lam/(1-np.exp(-lam*tr))
def cal_H(tr,tv):# Calcula H del elemento i
  return tr/tv
def cal_S(lam,ti):# Calcula S del elemento i
  return 1-np.exp(-lam*ti)
def equations(vars, *par):
    alfa = vars[0]
    Aesp_1,k0_1,e_1,Er_1,Q0_1,Aesp_2,k0_2,e_2,Er_2,Q0_2, Aesp_3,k0_3,e_3,Er_3,Q0_3 = par
    eq1 = ((1-(Aesp_2/Aesp_1)*(k0_1/k0_2)*(e_1/e_2))**(-1) - (1-(Aesp_3/Aesp_1)*(k0_1/k0_3)*(e_1/e_3))**(-1))*(Q0_1 - 0.429)/(Er_1**(alfa)) - ((1-(Aesp_2/Aesp_1)*(k0_1/k0_2)*(e_1/e_2))**(-1))*(Q0_2 - 0.429)/(Er_2**(alfa)) + ((1-(Aesp_3/Aesp_1)*(k0_1/k0_3)*(e_1/e_3))**(-1))*(Q0_3 - 0.429)/(Er_3**(alfa))
    return [eq1]
def cal_Q0_alfa_i(Q0,Er,alfa,rho):
  # calcula Q0_alfa del elemento i
  return (Q0-0.429)/(Er*np.exp(-rho*alfa))**alfa + 0.429/((2*alfa+1)*0.55**alfa)
  #return (Q0-0.429)/(Er)**alfa + 0.429/((2*alfa+1)*0.55**alfa)
def cal_f_alfa(Q0_alfa_c,Aesp_c,e_c,k0_c):
  # calcula f
  return ((k0_c[0]/k0_c[1])*(e_c[0]/e_c[1])*Q0_alfa_c[0]  - (Aesp_c[0]/Aesp_c[1])*Q0_alfa_c[1])/(Aesp_c[0]/Aesp_c[1] - (k0_c[0]/k0_c[1])*(e_c[0]/e_c[1]))

# parametros comparadores
# comparadores [Co Au Mo]
k0_c = np.array([1.32, 1.00, 0.0000846])     #
e_c = np.array([0.0113, 0.02272, 0.01605])    #
Q0_c = np.array([2.041, 15.712, 50.365])    #
Cn_c = np.array([36082, 52346, 39272])    # Area del fotopico
w_c = np.array([51.573/1000000, 0.1358/1000000, 2088.591/1000000])     # Masa del radio isotopo o comparador
rho_c = np.array([1.0, 1.0, 1.0])   #
lam_c = np.array([0.000000004167, 0.000002977, 0.00000292])   #
Er_c = np.array([136.0, 5.65, 241])   #

td_c = np.array([306314, 306314, 299161])   # Tiempo de decaimiento del elemento
tr_c = np.array([1500, 1500, 900])          # Tiempo de real del elemento
ti_c = np.array([10800, 10800, 10800])      # Tiempo de irradiación del elemento
tv_c = np.array([1478.0, 1478.0, 866.0])    # Tiempo de vivo

# paramentros de Co
# Co
k0_i = 1.32    #
e_i = 0.0113    #
Q0_i = 2.041    #
Cn_i = 36082/4    # Area del fotopico
w_i = 51.573/1000000     # Masa del radio isotopo o comparador
lam_i = 0.000000004167     # lambda
rho_i = 1.0   #
Er_i = 136.0    #
td_i = 306314.0      # Tiempo de decaimiento del elemento
tr_i = 1500.0      # Tiempo de real del elemento
ti_i = 10800.0      # Tiempo de irradiación del elemento
tv_i = 1478.0      # Tiempo de vivo
td_i = 306314.0      # Tiempo de decaimiento del elemento Co
tr_i = 1500.0      # Tiempo de real del elemento
ti_i = 10800.0      # Tiempo de irradiación del elemento
tv_i = 1478.0      # Tiempo de vivo

# Calculo de Aesp de los comparadores.
Aesp_c = np.zeros(len(k0_c))
Aesp_c[0] = Aesp(Cn_c[0],w_c[0],lam_c[0],tr_c[0],td_c[0],ti_c[0],tv_c[0])
Aesp_c[1] = Aesp(Cn_c[1],w_c[1],lam_c[1],tr_c[1],td_c[1],ti_c[1],tv_c[1])
Aesp_c[2] = Aesp(Cn_c[2],w_c[2],lam_c[2],tr_c[2],td_c[2],ti_c[2],tv_c[2])

st.markdown(f"Valor de Aesp_1: {Aesp_c[0]}")
st.markdown(f"Valor de Aesp_2: {Aesp_c[1]}")
st.markdown(f"Valor de Aesp_3: {Aesp_c[2]}")


# Calculo de alfa 
initial_guesses = [0.0]
par = (Aesp_c[0], k0_c[0], e_c[0], Er_c[0], Q0_c[0], Aesp_c[1], k0_c[1], e_c[1], Er_c[1], Q0_c[1], Aesp_c[2], k0_c[2], e_c[2], Er_c[2], Q0_c[2])
solution = root(equations, x0=initial_guesses, args = par)
alfa = float(solution.x)
st.markdown(f"Valor de alfa: {alfa}")

# calculo de Qo alfa de comparadores
Q0_alfa_c = np.zeros(len(k0_c))
Q0_alfa_c[0] = cal_Q0_alfa_i(Q0_c[0],Er_c[0],alfa,rho_c[0])
Q0_alfa_c[1] = cal_Q0_alfa_i(Q0_c[1],Er_c[1],alfa,rho_c[1])
Q0_alfa_c[2] = cal_Q0_alfa_i(Q0_c[2],Er_c[2],alfa,rho_c[2])

# calculo de f alfa de comparadores
f = cal_f_alfa(Q0_alfa_c,Aesp_c,e_c,k0_c)

# Entrada de fórmulas 
#formula_str = st.text_input("Fórmula", value="A * B / C")
# C = (Cn*exp(-lamb*td)*lamb*tr/((1-exp(-lamb*ti))*(1-exp(-lamb*tr))*w*tv*Aesp_Au))*(1/k0)*(e_Au/e)*((f + Q0_alfa_Au)/(f + (Q0-0.429)/(Er*exp(-rho*alfa))**alfa + 0.429/((2*alfa+1)*0.55**alfa))) #Concentracion
#formula_str = "Cn*exp(-lamb*td)*lamb*tr/((1-exp(-lamb*ti))*(1-exp(-lamb*tr))*w*tv)" # Aesp de Comparador Au
formula_str = "(Cn*exp(-lamb*td)*lamb*tr/((1-exp(-lamb*ti))*(1-exp(-lamb*tr))*w*tv*20562252179.58065))*(1/k0)*(0.00306420020500841/e)*(((-4.549267428275*(1.612/(136*exp(-alfa))**alfa + 0.429/((2*alfa+1)*0.55**alfa))+3.5492674282696*(15.283/(5.65*exp(-alfa))**alfa + 0.429/((2*alfa+1)*0.55**alfa))) + (1.612/(136*exp(-alfa))**alfa + 0.429/((2*alfa+1)*0.55**alfa)))/((-4.549267428275*(1.612/(136*exp(-alfa))**alfa + 0.429/((2*alfa+1)*0.55**alfa))+3.5492674282696*(15.283/(5.65*exp(-alfa))**alfa + 0.429/((2*alfa+1)*0.55**alfa))) + (Q0-0.429)/(Er*exp(-rho*alfa))**alfa + 0.429/((2*alfa+1)*0.55**alfa)))"
# Extraer variables únicas
try:
    variables = sorted(list(sp.sympify(formula_str).free_symbols), key=lambda x: str(x))
except Exception as e:
    st.error(f"Error al interpretar la fórmula: {e}")
    st.stop()
  
# Entrada de valores e incertidumbres
st.subheader("🔢 Valores e Incertidumbres de las Variables ")

valores = {}
incertidumbres = {}

# colocar valores iniciales de las variables.
Val_ini = (Cn_i, Er_i, Q0_i, alfa, e_i, k0_i, lam_i, rho_i, td_i, ti_i, tr_i, tv_i, w_i)
u_v_ini = (Cn_i*0.02, Er_i*0.02, Q0_i*0.05, alfa*0.05, e_i*0.02, k0_i*0.028, lam_i*0.05, rho_i*0.01, td_i*0.01, ti_i*0.01, tr_i*0.01, tv_i*0.01, w_i*0.02) 
i = 0
for var in variables:
    col1, col2 = st.columns(2)
    with col1:
        valor = st.number_input(f"Valor de {var}", value=Val_ini[i], format="%.12f", key=f"{var}_val")
    with col2:
        incertidumbre = st.number_input(f"Incertidumbre de {var}", value=u_v_ini[i], format="%.12f", key=f"{var}_err")
    valores[str(var)] = valor
    incertidumbres[str(var)] = incertidumbre
    i = i + 1

# Cálculo
try:
    # Definir símbolos
    simbolos = {str(v): sp.Symbol(str(v)) for v in variables}
    formula_sym = sp.sympify(formula_str)

    # Calcular valor central
    y_val = float(formula_sym.evalf(subs=valores))

    # Derivadas parciales (sensibilidades)
    st.subheader("📈 Resultados")
    st.write(f"**Valor calculado de la fórmula:** y = {y_val:.6f}")

    u_y_squared = 0
    contribuciones = []

    for v in variables:
        var_name = str(v)
        derivada = formula_sym.diff(v)
        sensibilidad = float(derivada.evalf(subs=valores))
        u_i = incertidumbres[var_name]
        contrib = (sensibilidad * u_i)**2
        u_rel_i = u_i / valores[var_name] if valores[var_name] != 0 else np.nan

        contribuciones.append({
            "Variable": var_name,
            "Sensibilidad ∂y/∂x": sensibilidad,
            "Incertidumbre": u_i,
            "Incertidumbre relativa": u_rel_i,
            "Contribución a u(y)²": contrib,
        })
        u_y_squared += contrib

    u_y = np.sqrt(u_y_squared)
    u_rel_y = u_y / y_val if y_val != 0 else np.nan

    # Calcular porcentaje de contribución
    for c in contribuciones:
        c["% Contribución"] = 100 * c["Contribución a u(y)²"] / u_y_squared if u_y_squared > 0 else np.nan

    # Mostrar resultados
    st.metric("Incertidumbre combinada u(y)", f"{u_y:.6f}")
    st.metric("Incertidumbre relativa de y", f"{u_rel_y:.2%}")

    st.subheader("📊 Contribución de cada variable a la incertidumbre")
    st.dataframe(contribuciones, use_container_width=True)

except Exception as e:
    st.error(f"Ocurrió un error en el cálculo: {e}")
