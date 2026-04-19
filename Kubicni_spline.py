import streamlit as st
import numpy as np
from sympy import symbols, expand, latex, lambdify
import matplotlib.pyplot as plt

# Postavke stranice
st.set_page_config(page_title="Kubični Spline", layout="centered")

# --- NOVO ZAGLAVLJE ---
st.latex(r"\huge \textbf{KUBIČNI SPLINE}")
st.latex(r"\large \text{Izradio: Vedran Bošnjak} \\\\")

# --- INPUT SEKCIJA ---
n = st.number_input(r"$\text{Unesite broj podataka }:$", min_value=3, value=4, step=1)

st.latex(r"\text{Unos koordinata:}")

podatci_input = []
cols = st.columns(2)
for i in range(n):
    with cols[0]:
        x_val = st.number_input(rf"$x_{{{i+1}}}:$", key=f"x{i}", value=float(i))
    with cols[1]:
        y_val = st.number_input(rf"$y_{{{i+1}}}:$", key=f"y{i}", value=float(i**2 if i % 2 == 0 else i+1))
    podatci_input.append((x_val, y_val))

# Gumb za pokretanje
if st.button("IZRAČUNAJ"):
    
    # Sortiranje podataka
    podatci = sorted(podatci_input, key=lambda t: t[0])

    lh = []
    ld = []
    for i in range(1, len(podatci)):
        h = round(podatci[i][0] - podatci[i-1][0], 3)
        lh.append(h)
        d = round((podatci[i][1] - podatci[i-1][1]) / h, 3)
        ld.append(d)

    # MATRICA H
    m_dim = len(podatci) - 2
    H_mat = np.zeros((m_dim, m_dim))
    r_vec = np.zeros(m_dim)

    for i in range(m_dim):
        H_mat[i, i] = 2 * (lh[i] + lh[i+1])
        if i > 0:
            H_mat[i, i-1] = lh[i]
        if i < m_dim - 1:
            H_mat[i, i+1] = lh[i+1]
        r_vec[i] = 6 * (ld[i+1] - ld[i])

    try:
        ls_mid = np.linalg.solve(H_mat, r_vec)
        ls = np.concatenate(([0], ls_mid, [0]))
        ls = [round(val, 3) for val in ls]

        lb = []
        for i in range(len(ld)):
            b = ld[i] - (ls[i+1] + 2 * ls[i]) * lh[i] / 6
            lb.append(round(b, 3))

        x = symbols("x")
        polinomi = []
        for i in range(len(podatci) - 1):
             x1, y1 = podatci[i]
             x2, y2 = podatci[i+1]
             h, s1, s2 = lh[i], ls[i], ls[i+1]

             f = ((s1 / (6 * h)) * (x2 - x)**3 + (s2 / (6 * h)) * (x - x1)**3 + 
                  (y1 / h - s1 * h / 6) * (x2 - x) + (y2 / h - s2 * h / 6) * (x - x1))
             polinomi.append(expand(f))

                # --- ISPIS FUNKCIJE S(x) ---
        st.divider()
        st.latex(r"\text{KUBIČNI SPLINE}")
        
        # Konstrukcija "cases" okruženja za sve polinome
        latex_kod = r"S(x) = \begin{cases} "
        for i in range(len(polinomi)):
            linija = rf"{latex(polinomi[i].evalf(3))} , & {podatci[i][0]} \le x \le {podatci[i+1][0]}"
            # Dodajemo prijelaz u novi red (\\) ako nije zadnji polinom
            if i < len(polinomi) - 1:
                linija += r" \\ "
            latex_kod += linija
        latex_kod += r" \end{cases}"
        
        st.latex(latex_kod)

        # --- GRAF ---
        st.divider()
        st.latex(r"\text{GRAF SPLINEA}")
        fig, ax = plt.subplots(figsize=(10, 6))
        for i in range(len(polinomi)):
            func = lambdify(x, polinomi[i], modules=['numpy'])
            x_vals = np.linspace(podatci[i][0], podatci[i+1][0], 100)
            ax.plot(x_vals, func(x_vals), linewidth=2)
        
        ax.scatter([p[0] for p in podatci], [p[1] for p in podatci], color='red', zorder=5)
        ax.grid(True, linestyle='--', alpha=0.6)
        st.pyplot(fig)

        # --- DETALJAN RAČUN ---
        st.latex(r"\text{RAČUN}")
        c1, c2 = st.columns(2)
        with c1:
            st.latex(r"\textbf{Podatci:}")
            for i, p_val in enumerate(podatci):
                st.latex(rf"T_{{{i+1}}} = ({p_val[0]}, {p_val[1]})")
            
            st.latex(r"\textbf{Koeficijenti } h:")
            for i, h_val in enumerate(lh):
                st.latex(rf"h_{{{i+1}}} = {h_val}")
            
            st.latex(r"\textbf{Koeficijenti } d:")
            for i, d_val in enumerate(ld):
                st.latex(rf"d_{{{i+1}}} = {d_val}")

        with c2:
            st.latex(r"\textbf{Koeficijenti } s:")
            for i, s_val in enumerate(ls):
                st.latex(rf"s_{{{i+1}}} = {s_val}")
            
            st.latex(r"\textbf{Koeficijenti } b:")
            for i, b_val in enumerate(lb):
                st.latex(rf"b_{{{i+1}}} = {b_val}")

        st.latex(r"\textbf{Matrica } H:")
        st.dataframe(H_mat)

        st.latex(r"\textbf{Vektor } r:")
        st.dataframe(r_vec)

    except Exception as e:
        st.error(f"Greška: {e}")
