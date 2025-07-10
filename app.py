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

# Calculo de Aesp de los comparadores.
# parametros comparadores
# comparadores [Co Au Mo]
k0_c = np.array([1.32, 1.00, 0.0000846])     #
e_c = np.array([0.0015999456766857, 0.00306420020500841, 0.0160532014949982 ])    #
Q0_c = np.array([2.041, 15.712, 50.365])    #
Cn_c = np.array([11880, 32416, 40639])    # Area del fotopico
w_c = np.array([40.83662/1000000, 0.1044/1000000, 2090.99/1000000])     # Masa del radio isotopo o comparador
rho_c = np.array([1.0, 1.0, 0.973])   #
lam_c = np.array([0.000000004167, 0.000002977, 0.00000292])   #
Er_c = np.array([136.0, 5.65, 241])   #
td_c = np.array([15962, 15962, 97299])   # Tiempo de decaimiento del elemento
tr_c = np.array([57602, 57602, 10000])   # Tiempo de real del elemento
ti_c = np.array([1800, 1800, 1800])      # Tiempo de irradiación del elemento
tv_c = np.array([57455, 57455, 9953])    # Tiempo de vivo

Aesp_c = np.zeros(len(k0_c))
Aesp_c[0] = Aesp(Cn_c[0],w_c[0],lam_c[0],tr_c[0],td_c[0],ti_c[0],tv_c[0])
Aesp_c[1] = Aesp(Cn_c[1],w_c[1],lam_c[1],tr_c[1],td_c[1],ti_c[1],tv_c[1])
Aesp_c[2] = Aesp(Cn_c[2],w_c[2],lam_c[2],tr_c[2],td_c[2],ti_c[2],tv_c[2])

# Calculo de alfa 
initial_guesses = [0.0]
par = (Aesp_c[0], k0_c[0], e_c[0], Er_c[0], Q0_c[0], Aesp_c[1], k0_c[1], e_c[1], Er_c[1], Q0_c[1], Aesp_c[2], k0_c[2], e_c[2], Er_c[2], Q0_c[2])
solution = root(equations, x0=initial_guesses, args = par)
alfa = solution.x
st.markdown("Valor de alfa: {alfa}")

# Entrada de fórmulas 
#formula_str = st.text_input("Fórmula", value="A * B / C")
formula_str = "Cn*exp(-lamb*td)*lamb*tr/((1-exp(-lamb*ti))*(1-exp(-lamb*tr))*w*tv)" # Aesp de Comparador Au
# Extraer variables únicas
try:
    variables = sorted(list(sp.sympify(formula_str).free_symbols), key=lambda x: str(x))
except Exception as e:
    st.error(f"Error al interpretar la fórmula: {e}")
    st.stop()

formula_str1 = "Cn*exp(-lamb*td)*lamb*tr/((1-exp(-lamb*ti))*(1-exp(-lamb*tr))*w*tv)" # Aesp del elementos Co
# Extraer variables únicas
try:
    variables1 = sorted(list(sp.sympify(formula_str).free_symbols), key=lambda x: str(x))
except Exception as e:
    st.error(f"Error al interpretar la fórmula: {e}")
    st.stop()

# Entrada de valores e incertidumbres
st.subheader("🔢 Valores e Incertidumbres de las Variables ")

valores = {}
incertidumbres = {}

# colocar valores iniciales de las variables.

for var in variables:
    col1, col2 = st.columns(2)
    with col1:
        valor = st.number_input(f"Valor de {var}", value=1.0, format="%.12f", key=f"{var}_val")
    with col2:
        incertidumbre = st.number_input(f"Incertidumbre de {var}", value=0.01, format="%.12f", key=f"{var}_err")
    valores[str(var)] = valor
    incertidumbres[str(var)] = incertidumbre

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
