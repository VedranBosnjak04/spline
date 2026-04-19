import streamlit as st
import numpy as np
from sympy import symbols, expand, latex, lambdify
import matplotlib.pyplot as plt

# Postavke stranice
st.set_page_config(page_title="Kubični Spline", layout="centered")

st.title("Izračun Kubičnog Splinea - Vedran Bošnjak")

# --- INPUT SEKCIJA ---
n = st.number_input("Unesite broj podataka:", min_value=3, value=4, step=1)

st.write("### Unos koordinata:")

podatci_input = []
cols = st.columns(2)
for i in range(n):
    with cols[0]:
        x_val = st.number_input(f"x{i+1}:", key=f"x{i}", value=float(i))
    with cols[1]:
        y_val = st.number_input(f"y{i+1}:", key=f"y{i}", value=float(i**2 if i%2==0 else i+2))
    podatci_input.append((x_val, y_val))

# Gumb za pokretanje
if st.button("IZRAČUNAJ"):
    
    # Sortiranje podataka po x (nužno za spline)
    podatci = sorted(podatci_input, key=lambda t: t[0])

    # Računanje h (razmaci) i d (podijeljene razlike)
    lh = []
    ld = []
    for i in range(1, len(podatci)):
        h = podatci[i][0] - podatci[i-1][0]
        lh.append(h)
        d = (podatci[i][1] - podatci[i-1][1]) / h
        ld.append(d)

    # --- POPRAVLJENA LOGIKA MATRICE H ---
    # Dimenzija sustava je n-2 jer su rubni uvjeti s0=0 i sn=0
    m = len(podatci) - 2
    H = np.zeros((m, m))
    r = np.zeros(m)

    for i in range(m):
        # Glavna dijagonala: 2 * (h_i + h_{i+1})
        H[i, i] = 2 * (lh[i] + lh[i+1])
        
        # Susjedni h koeficijenti
        if i > 0:
            H[i, i-1] = lh[i]
        if i < m - 1:
            H[i, i+1] = lh[i+1]
        
        # Desna strana (vektor r)
        r[i] = 6 * (ld[i+1] - ld[i])

    # Rješavanje sustava
    try:
        # ls su naše druge derivacije (momenti)
        ls_middle = np.linalg.solve(H, r)
        # Dodajemo 0 na početak i kraj za prirodni spline
        ls = np.concatenate(([0], ls_middle, [0]))

        # Računanje koeficijenata b
        lb = []
        for i in range(len(ld)):
            # Stabilna formula za linearni dio splinea
            b = ld[i] - (ls[i+1] + 2*ls[i]) * lh[i] / 6
            lb.append(round(b, 4))

        x = symbols("x")
        polinomi = []
        for i in range(len(podatci) - 1):
             x1, y1 = podatci[i]
             x2, y2 = podatci[i+1]
             h = lh[i]
             s1, s2 = ls[i], ls[i+1]

             # Standardna formula za kubični spline segment
             f = ((s1/(6*h)) * (x2 - x)**3 + 
                  (s2/(6*h)) * (x - x1)**3 + 
                  (y1/h - s1*h/6) * (x2 - x) + 
                  (y2/h - s2*h/6) * (x - x1))
             
             polinomi.append(expand(f))

        # --- PRIKAZ REZULTATA ---
        st.divider()
        st.header("KUBIČNI SPLINE")
        
        st.write("$$S(x) = \begin{cases}$$")
        for i in range(len(polinomi)):
            st.latex(rf"{latex(polinomi[i].evalf(4))} , \quad {podatci[i][0]} \le x \le {podatci[i+1][0]}")
        st.write("$$\end{cases}$$")

        # --- GRAF ---
        st.divider()
        st.header("GRAF SPLINEA")
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Crtanje svakog segmenta posebno
        for i in range(len(polinomi)):
            polinom_funkcija = lambdify(x, polinomi[i], modules=['numpy'])
            x_interval = np.linspace(podatci[i][0], podatci[i+1][0], 100)
            y_interval = polinom_funkcija(x_interval)
            
            ax.plot(x_interval, y_interval, linewidth=2.5, label=f"S{i+1}(x)")
        
        # Originalne točke
        x_tocke = [p[0] for p in podatci]
        y_tocke = [p[1] for p in podatci]
        ax.scatter(x_tocke, y_tocke, color='red', s=80, zorder=5, label='Točke')
        
        ax.set_title("Vizualizacija Kubičnog Splinea")
        ax.grid(True, linestyle='--', alpha=0.7)
        ax.legend()
        st.pyplot(fig)
        
        # --- MATEMATIČKI DETALJI ---
        st.divider()
        with st.expander("Vidi detalje proračuna"):
            col1, col2 = st.columns(2)
            with col1:
                st.write("**Koeficijenti h (razmaci):**", lh)
                st.write("**Podijeljene razlike d:**", ld)
                st.write("**Matrica H:**")
                st.dataframe(H)
            with col2:
                st.write("**Vektor r:**", r)
                st.write("**Koeficijenti s (momenti):**", ls)
                st.write("**Koeficijenti b:**", lb)

    except Exception as e:
        st.error(f"Došlo je do greške: {e}")
        st.info("Provjerite da su X vrijednosti različite i unesene pravilno.")

