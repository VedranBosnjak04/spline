import streamlit as st
import numpy as np
from sympy import symbols, expand, latex, lambdify
import matplotlib.pyplot as plt

# Postavke stranice
st.set_page_config(page_title="Kubični Spline", layout="centered")

# --- NOVO ZAGLAVLJE ---
st.latex(r"\Large \textbf{SPLINE KALKULATOR}")
st.latex(r"\normalsize \mathsf{Izradio: Vedran\ Bošnjak} \\\\")

# --- INPUT SEKCIJA ---
n = st.number_input(r"$\textbf{Unesite broj podataka:}$", min_value=3, value=4, step=1)

st.latex(r"\textbf{Unos koordinata:}")

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
        st.latex(r"\textbf{KUBIČNI SPLINE}")
        
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
        st.latex(r"\textbf{GRAF SPLINEA}")
        fig, ax = plt.subplots(figsize=(10, 6))
        for i in range(len(polinomi)):
            func = lambdify(x, polinomi[i], modules=['numpy'])
            x_vals = np.linspace(podatci[i][0], podatci[i+1][0], 100)
            ax.plot(x_vals, func(x_vals), linewidth=2)
        
        ax.scatter([p[0] for p in podatci], [p[1] for p in podatci], color='red', zorder=5)
        ax.grid(True, linestyle='--', alpha=0.6)
        st.pyplot(fig)

        # --- DETALJAN RAČUN ---
                # --- DETALJAN RAČUN (LaTeX stil) ---
        st.divider()
        st.latex(r"\textbf{RAČUN}")
        
        c1, c2 = st.columns(2)
        
        with c1:
            st.markdown('<div style="text-align: left;">', unsafe_allow_html=True)
            
            # Podatci (x i y koordinate)
            p_lat = r"\textbf{Podatci:}\\" + r"\begin{aligned}"
            for i, p_val in enumerate(podatci):
                # Čistimo i x i y vrijednost
                x_cist = "%g" % round(p_val[0], 3)
                y_cist = "%g" % round(p_val[1], 3)
                p_lat += rf"T_{{{i+1}}} &= ({x_cist}, {y_cist}) \\"
            p_lat += r"\end{aligned}"
            st.latex(p_lat)

            # Koeficijenti h
            h_lat = r"\textbf{Koeficijenti } h_k:\\" + r"\begin{aligned}"
            for i, h_val in enumerate(lh):
                h_lat += rf"h_{{{i+1}}} &= {"%g" % round(h_val, 3)} \\"
            h_lat += r"\end{aligned}"
            st.latex(h_lat)

            # Koeficijenti d
            d_lat = r"\textbf{Koeficijenti } d_k:\\" + r"\begin{aligned}"
            for i, d_val in enumerate(ld):
                d_lat += rf"d_{{{i+1}}} &= {"%g" % round(d_val, 3)} \\"
            d_lat += r"\end{aligned}"
            st.latex(d_lat)
            st.markdown('</div>', unsafe_allow_html=True)

        with c2:
            st.markdown('<div style="text-align: left;">', unsafe_allow_html=True)
            
            # Koeficijenti s
            s_lat = r"\textbf{Koeficijenti } s_k:\\" + r"\begin{aligned}"
            for i, s_val in enumerate(ls):
                s_lat += rf"s_{{{i+1}}} &= {"%g" % round(s_val, 3)} \\"
            s_lat += r"\end{aligned}"
            st.latex(s_lat)
            
            # Koeficijenti b
            b_lat = r"\textbf{Koeficijenti } b_k:\\" + r"\begin{aligned}"
            for i, b_val in enumerate(lb):
                b_lat += rf"b_{{{i+1}}} &= {"%g" % round(b_val, 3)} \\"
            b_lat += r"\end{aligned}"
            st.latex(b_lat)
            st.markdown('</div>', unsafe_allow_html=True)

        # Matrica H u LaTeXu
                # --- MATRICA H (LaTeX stil - ČISTI PRIKAZ) ---
        st.markdown('<div style="text-align: left;">', unsafe_allow_html=True)
        
        # Generiranje LaTeX koda za matricu H
        h_matrix_latex = r"\textbf{Matrica } H = \begin{bmatrix} "
        for row in H_mat:
            # Koristimo %g koji miče nepotrebne nule (npr. 2.0 postane 2)
            h_matrix_latex += " & ".join(["%g" % round(val, 3) for val in row]) + r" \\ "
        h_matrix_latex += r"\end{bmatrix}"
        st.latex(h_matrix_latex)

        # --- VEKTOR r (Pravi stupac vektor - ČISTI PRIKAZ) ---
        r_vector_latex = r"\textbf{Vektor } r = \begin{bmatrix} "
        for val in r_vec:
            r_vector_latex += "%g" % round(val, 3) + r" \\ "
        r_vector_latex += r"\end{bmatrix}"
        st.latex(r_vector_latex)
        
        st.markdown('</div>', unsafe_allow_html=True)

        
    except Exception as e:
        st.error(f"Greška: {e}")
