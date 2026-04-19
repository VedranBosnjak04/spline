import streamlit as st
import numpy as np
from sympy import symbols, expand, Poly, latex, lambdify

# Postavke stranice za ljepši izgled
st.set_page_config(page_title="Kubični Spline", layout="centered")

st.title("Izračun Kubičnog Splinea - Vedran Bošnjak")

# --- INPUT SEKCIJA ---
# Umjesto while True, koristimo brojčani unos
n = st.number_input("Unesite broj podataka:", min_value=4, value=4, step=1)

st.write("### Unos koordinata:")

podatci_input = []
cols = st.columns(2)
for i in range(n):
    with cols[0]:
        x_val = st.number_input(f"x{i+1}:", key=f"x{i}", value=float(i))
    with cols[1]:
        y_val = st.number_input(f"y{i+1}:", key=f"y{i}", value=float(i**2))
    podatci_input.append((round(x_val, 3), round(y_val, 3)))

# Gumb za pokretanje
if st.button("IZRAČUNAJ"):
    
    # --- TVOJ ORIGINALNI KOD (NEPROMIJENJEN) ---
    podatci = sorted(podatci_input, key=lambda t: t[0])

    lh = []
    ld = []
    for i in range(1, len(podatci)):
        h = round(podatci[i][0] - podatci[i-1][0], 3)
        lh.append(h)
        d = round((podatci[i][1] - podatci[i-1][1]) / h, 3)
        ld.append(d)

    H = []
    r = []
    k = 0
    m = len(lh) - 2
    while k <= m:
        red = [0]*(m+1)
        red[k] = 2*(lh[k] + lh[k+1])

        if k == 0:
            if (m+1) > 1: red[k+1] = lh[1]
        elif k == m:
            red[k-1] = lh[1]
        else:
            red[k-1] = lh[1]
            red[k+1] = lh[1]

        H.append(red)
        k += 1

    p = 0
    l = len(ld) - 2
    while p <= l:
        red = []
        red.append(6*(ld[p+1] - ld[p]))
        r.append(red)
        p += 1

    # Rješavanje sustava
    try:
        ls = np.round(np.linalg.solve(H, r).flatten(), 3)
        ls = np.insert(ls, 0, 0)
        ls = np.append(ls, 0)
        ls = ls.tolist()

        lb = []
        for i in range(len(ld)):
            b = ld[i] - (ls[i+1] - ls[i])*lh[i]/6
            lb.append(round(b, 3))

        x = symbols("x")
        polinomi = []
        for i in range(1, len(podatci)):
             x1 = podatci[i-1][0]
             x2 = podatci[i][0]
             y = podatci[i-1][1]
             s1 = ls[i-1]
             s2 = ls[i]
             h = lh[i-1]
             b_coeff = lb[i-1]

             f = (y - s1*(h**2)/6) + b_coeff*(x - x1) + (s1/(6*h))*(x2 - x)**3 + (s2/(6*h))*(x - x1)**3
             f = expand(f)
             pol = Poly(f, x).as_expr().evalf(3)
             polinomi.append(pol)

        # --- UREDAN PRIKAZ REZULTATA ---
        st.divider()
        st.header("------------------ KUBIČNI SPLINE ------------------")
        
        st.write("latex(S(x) = {)")
        for i in range(len(polinomi)):
            # LaTeX prikaz za matematičku ljepotu
            st.latex(rf"{latex(polinomi[i])} , \quad {podatci[i][0]:g} \le x \le {podatci[i+1][0]:g}")
        st.write("latex(})")

        st.divider()

        # --- NOVI DIO: CRTANJE GRAFA ---
        # --- NOVI DIO: CRTANJE GRAFA (Popravljen x_sym) ---
        st.divider()
        st.header("📈 GRAF SPLINEA")
        
        x_puni = []
        y_puni = []
        
        for i in range(len(polinomi)):
            # Ovdje koristimo 'x' jer je tako definirana tvoja varijabla gore
            polinom_funkcija = lambdify(x, polinomi[i], modules=['numpy'])
            
            x_interval = np.linspace(podatci[i][0], podatci[i+1][0], 50)
            y_interval = polinom_funkcija(x_interval)
            
            x_puni.extend(x_interval)
            y_puni.extend(y_interval)
            
        fig, ax = plt.subplots(figsize=(10, 6))
        
        ax.plot(x_puni, y_puni, color='blue', label='Kubični Spline S(x)', linewidth=2)
        
        x_tocke = [p[0] for p in podatci]
        y_tocke = [p[1] for p in podatci]
        ax.scatter(x_tocke, y_tocke, color='red', label='Unesene točke T', s=100, zorder=5)
        
        ax.set_title("Graf Kubičnog Splinea")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.grid(True, linestyle='--', alpha=0.7)
        ax.legend()
        
        st.pyplot(fig)
        
        st.header("--------------------- RAČUN ---------------------")
        
        c1, c2 = st.columns(2)
        c1, c2 = st.columns(2)
        with c1:
            st.write("**Podatci:**")
            for i, p_val in enumerate(podatci):
                st.write(f"T{i+1} = {p_val}")
            
            st.write("**Koeficijenti h:**")
            for i, h_val in enumerate(lh):
                st.write(f"h{i+1} = {h_val}")
            
            st.write("**Koeficijenti d:**")
            for i, d_val in enumerate(ld):
                st.write(f"d{i+1} = {d_val}")

        with c2:
            st.write("**Koeficijenti s:**")
            for i, s_val in enumerate(ls):
                st.write(f"s{i+1} = {s_val}")
            
            st.write("**Koeficijenti b:**")
            for i, b_val in enumerate(lb):
                st.write(f"b{i+1} = {b_val}")

        st.write("**Matrica H:**")
        st.dataframe(np.array(H))

        st.write("**Matrica r:**")
        st.dataframe(np.array(r))

    except Exception as e:
        st.error(f"Došlo je do greške u izračunu: {e}")
