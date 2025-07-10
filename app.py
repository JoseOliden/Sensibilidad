import streamlit as st
import numpy as np
import sympy as sp

st.set_page_config(page_title="Propagación de Incertidumbre", layout="wide")
st.title("🧮 Análisis de Incertidumbre en Fórmulas")

st.markdown("Ingrese una fórmula con variables (por ejemplo: `A * B / C`)")

# Entrada de fórmula
formula_str = st.text_input("Fórmula", value="A * B / C")

# Extraer variables únicas
try:
    variables = sorted(list(sp.sympify(formula_str).free_symbols), key=lambda x: str(x))
except Exception as e:
    st.error(f"Error al interpretar la fórmula: {e}")
    st.stop()

# Entrada de valores e incertidumbres
st.subheader("🔢 Valores e Incertidumbres de las Variables")

valores = {}
incertidumbres = {}

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
