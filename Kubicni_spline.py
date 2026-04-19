import streamlit as st
import numpy as np
from sympy import symbols, expand, Poly, latex

# Postavke za ljepši prikaz brojeva
np.set_printoptions(suppress=True)

st.title("Kubični Spline - Tvoj Kod")

# --- TVOJA LOGIKA (Samo prilagođena za Web sučelje) ---

# Umjesto while True i inputa, koristimo brojčani unos
n = st.number_input("Unesite broj podataka:", min_value=1, step=1, value=4)

if n < 4:
    st.warning("Minimalan broj podataka je 4")
else:
    podatci = []
    
    st.write("### Unesite koordinate:")
    # Koristimo stupce da uštedimo prostor
    for i in range(n):
        col1, col2 = st.columns(2)
        with col1:
            x_val = st.number_input(f"x{i+1} koordinata:", key=f"x{i}", format="%.3f")
        with col2:
            y_val = st.number_input(f"y{i+1} koordinata:", key=f"y{i}", format="%.3f")
        podatci.append((round(x_val, 3), round(y_val, 3)))

    if st.button("POKRENI IZRAČUN"):
        # --- OVDJE POČINJE TVOJA MATEMATIKA (Doslovno kopirana) ---
        podatci.sort(key=lambda t: t[0])

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
                if len(red) > 1: red[k+1] = lh[1]
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

        ls = np.round(np.linalg.solve(H, r).flatten(), 3)
        ls = np.insert(ls, 0, 0)
        ls = np.append(ls, 0)
        ls = ls.tolist()

        lb = []
        for i in range(len(ld)):
            b_coeff = ld[i] - (ls[i+1] - ls[i])*lh[i]/6
            lb.append(b_coeff)

        x = symbols("x")
        polinomi = []
        for i in range(1, len(podatci)):
             x1 = podatci[i-1][0]
             x2 = podatci[i][0]
             y = podatci[i-1][1]
             s1 = ls[i-1]
             s2 = ls[i]
             h = lh[i-1]
             b = lb[i-1]

             f = (y - s1*(h**2)/6) + b*(x - x1) + (s1/(6*h))*(x2 - x)**3 + (s2/(6*h))*(x - x1)**3
             f = expand(f)
             pol = Poly(f, x).as_expr().evalf(3)
             polinomi.append(pol)

        # --- TVOJI PRINTOVI (Prebačeni u st.write/st.latex) ---
        st.header("KUBIČNI SPLINE")
        st.write("S(x) = {")
        for i in range(len(polinomi)):
            st.latex(rf"{latex(polinomi[i])} , \quad {podatci[i][0]:g} \le x \le {podatci[i+1][0]:g}")
        st.write("}")

        st.divider()
        st.subheader("RAČUN")
        st.write("**Podatci:**", podatci)
        st.write("**Koeficijenti h:**", lh)
        st.write("**Koeficijenti d:**", ld)
        
        st.write("**Matrica H:**")
        st.write(np.array(H))
        
        st.write("**Matrica r:**")
        st.write(np.array(r))

        st.write("**Koeficijenti s:**", ls)
        st.write("**Koeficijenti b:**", lb)
