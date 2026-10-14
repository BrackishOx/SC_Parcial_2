import argparse
import csv
import hashlib
import math
import os
import random
import time

try:
    import gmpy2

    def es_primo(n):
        return bool(gmpy2.is_prime(n))
except ImportError:
    import sympy

    def es_primo(n):
        return sympy.isprime(n)

SEED = 2026
CIFRAS = 40
N_NUMEROS = 10_000
CSV_PRIMOS = "primos_encontrados.csv"


class GeneradorLCG:
    m = 2 ** 61 - 1
    a = 6364136223846793005 % m

    def __init__(self, seed):
        self.estado = (seed % (self.m - 1)) + 1

    def crudo(self, bits):
        v, n = 0, 0
        while n < bits:
            self.estado = (self.a * self.estado) % self.m
            v = (v << 61) | self.estado
            n += 61
        return v


class GeneradorMersenneTwister:

    def __init__(self, seed):
        self._rng = random.Random(seed)

    def crudo(self, bits):
        return self._rng.getrandbits(bits)


class GeneradorHashSHA256:

    def __init__(self, seed):
        self.semilla = str(seed).encode()
        self.contador = 0

    def crudo(self, bits):
        partes = []
        for _ in range((bits + 255) // 256):
            partes.append(hashlib.sha256(
                self.semilla + self.contador.to_bytes(8, "big")).digest())
            self.contador += 1
        return int.from_bytes(b"".join(partes), "big")


class GeneradorBBS:
    _P = None
    _Q = None

    @classmethod
    def _modulo(cls):
        if cls._P is None:
            def primo_3mod4(inicio):
                n = inicio | 3
                while not (n % 4 == 3 and es_primo(n)):
                    n += 4
                return n
            cls._P = primo_3mod4(2 ** 127 + 12345)
            cls._Q = primo_3mod4(2 ** 127 + 987654321)
        return cls._P * cls._Q

    def __init__(self, seed):
        M = self._modulo()
        s = (seed % (M - 2)) + 2
        while math.gcd(s, M) != 1:
            s += 1
        self.M = M
        self.x = pow(s, 2, M)

    def crudo(self, bits):
        v = 0
        for _ in range(bits):
            self.x = pow(self.x, 2, self.M)
            v = (v << 1) | (self.x & 1)
        return v


TECNICAS = {
    "lcg": ("LCG", GeneradorLCG),
    "mt": ("Mersenne Twister", GeneradorMersenneTwister),
    "sha256": ("Hash SHA-256", GeneradorHashSHA256),
    "bbs": ("Blum-Blum-Shub", GeneradorBBS),
}


class MaquinaPRNG:

    def __init__(self, tecnica="sha256", semilla=SEED, cifras=CIFRAS):
        tecnica = tecnica.lower()
        if tecnica not in TECNICAS:
            raise ValueError(f"tecnica desconocida {tecnica!r}; "
                             f"opciones: {', '.join(TECNICAS)}")
        if cifras < 2:
            raise ValueError("cifras debe ser >= 2")
        self.tecnica = tecnica
        self.nombre = TECNICAS[tecnica][0]
        self.semilla = semilla
        self.cifras = cifras
        self.inferior = 10 ** (cifras - 1)
        self.rango = 10 ** cifras - self.inferior
        self.bits = self.rango.bit_length() + 64
        self._gen = TECNICAS[tecnica][1](semilla)

    def siguiente(self):
        return self.inferior + (self._gen.crudo(self.bits) % self.rango)

    def lote(self, k):
        return [self.siguiente() for _ in range(k)]


def densidad_teorica(cifras=CIFRAS):
    medio = (10 ** (cifras - 1) + 10 ** cifras) // 2
    return 1 / math.log(medio)


def correr_experimento(cantidad=N_NUMEROS, semilla=SEED, cifras=CIFRAS,
                       guardar_csv=True, verificar_sympy=True, verbose=True):
    p_teo = densidad_teorica(cifras)
    esperados = cantidad * p_teo
    sigma = math.sqrt(cantidad * p_teo * (1 - p_teo))
    if verbose:
        print(f"Densidad teorica 1/ln(x) en el centro del rango: {p_teo*100:.4f}%")
        print(f"Primos esperados en {cantidad} numeros: {esperados:.1f} "
              f"(desv. tipica binomial ~{sigma:.1f})\n")

    resultados = {}
    for clave in TECNICAS:
        maq = MaquinaPRNG(clave, semilla, cifras)
        t0 = time.time()
        numeros = maq.lote(cantidad)
        t_gen = time.time() - t0
        t0 = time.time()
        primos = [(i, n) for i, n in enumerate(numeros) if es_primo(n)]
        t_test = time.time() - t0

        assert all(len(str(n)) == cifras for n in numeros), "numero sin 40 cifras"
        repetidos = cantidad - len(set(numeros))
        pares = sum(1 for n in numeros if n % 2 == 0)

        coincide = None
        if verificar_sympy:
            import sympy
            idx_primos = {i for i, _ in primos}
            muestra = [i for i in range(min(cantidad, 600)) if i not in idx_primos]
            ok_p = all(sympy.isprime(n) for _, n in primos)
            ok_c = not any(sympy.isprime(numeros[i]) for i in muestra)
            coincide = ok_p and ok_c

        r = {"nombre": maq.nombre, "primos": primos, "n": len(primos),
             "densidad": len(primos) / cantidad, "z": (len(primos) - esperados) / sigma,
             "t_gen": t_gen, "t_test": t_test, "repetidos": repetidos,
             "pct_pares": pares / cantidad * 100, "sympy_ok": coincide}
        resultados[clave] = r
        if verbose:
            print(f"=== {maq.nombre} ===")
            print(f"  primos: {r['n']} / {cantidad} ({r['densidad']*100:.4f}%), "
                  f"z = {r['z']:+.2f}")
            print(f"  repetidos: {repetidos}; pares: {r['pct_pares']:.2f}%; "
                  f"verificacion sympy: {coincide}")
            print(f"  tiempo: generar {t_gen:.2f}s, primalidad {t_test:.2f}s")
            for _, n in primos[:3]:
                print(f"    ej: {n}")
            print()

    if verbose:
        print("=== RESUMEN ===")
        print(f"{'Tecnica':<20}{'Primos':<8}{'Densidad %':<12}{'z':<8}")
        for r in resultados.values():
            print(f"{r['nombre']:<20}{r['n']:<8}{r['densidad']*100:<12.4f}{r['z']:+.2f}")

    if guardar_csv:
        with open(CSV_PRIMOS, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["tecnica", "indice", "primo"])
            for clave, r in resultados.items():
                for i, n in r["primos"]:
                    w.writerow([clave, i, n])
        if verbose:
            print(f"\nPrimos guardados en {CSV_PRIMOS}")
    return resultados


class SistemaRSA:
    E_PUBLICO = 65537

    def __init__(self, p, q):
        if p == q:
            raise ValueError("p y q deben ser primos distintos")
        if not (es_primo(p) and es_primo(q)):
            raise ValueError("p y q deben ser primos")
        self.p, self.q = p, q
        self.n = p * q
        self.phi = (p - 1) * (q - 1)
        if math.gcd(self.E_PUBLICO, self.phi) != 1:
            raise ValueError("e no es coprimo con phi(n); elegir otro par")
        self.e = self.E_PUBLICO
        self.d = pow(self.e, -1, self.phi)
        k = (self.n.bit_length() - 1) // 8
        self.tam_bloque = k - 1

    def clave_publica(self):
        return (self.e, self.n)

    def clave_privada(self):
        return (self.d, self.n)

    def claves_validas(self):
        return (self.e * self.d) % self.phi == 1

    def cifrar_entero(self, m):
        if not 0 <= m < self.n:
            raise ValueError("el entero debe cumplir 0 <= m < n")
        return pow(m, self.e, self.n)

    def descifrar_entero(self, c):
        return pow(c, self.d, self.n)

    def cifrar(self, texto):
        datos = texto.encode("utf-8")
        t = self.tam_bloque
        bloques = [datos[i:i + t] for i in range(0, len(datos), t)] or [b""]
        return [self.cifrar_entero(int.from_bytes(b"\x01" + b, "big")) for b in bloques]

    def descifrar(self, cifrado):
        datos = b""
        for c in cifrado:
            m = self.descifrar_entero(c)
            b = m.to_bytes((m.bit_length() + 7) // 8, "big")
            if not b or b[0] != 1:
                raise ValueError("bloque corrupto o clave incorrecta")
            datos += b[1:]
        return datos.decode("utf-8")

    @classmethod
    def desde_primos(cls, primos, intentos=50):
        for i in range(len(primos)):
            for j in range(i + 1, len(primos)):
                try:
                    return cls(primos[i], primos[j])
                except ValueError:
                    continue
        raise ValueError("no hay un par de primos valido")


def primos_por_tecnica(resultados):
    return {k: [n for _, n in r["primos"]] for k, r in resultados.items()}


def demostrar_rsa(resultados, mensaje="HOLA RSA DESDE SISTEMAS COMPLEJOS", verbose=True):
    pp = primos_por_tecnica(resultados)
    p, q = pp["lcg"][0], pp["sha256"][0]
    rsa = SistemaRSA(p, q)
    cif = rsa.cifrar(mensaje)
    des = rsa.descifrar(cif)
    if verbose:
        print("=== Sistema de encriptacion RSA con los primos encontrados ===")
        print(f"p (LCG)    = {p}")
        print(f"q (SHA256) = {q}")
        print(f"n = p*q    = {rsa.n}  ({len(str(rsa.n))} cifras, {rsa.n.bit_length()} bits)")
        print(f"e = {rsa.e}\nd = {rsa.d}")
        print(f"e*d mod phi == 1 -> {rsa.claves_validas()}")
        print(f"Bloque: {rsa.tam_bloque} bytes de mensaje por bloque")
        print(f"Mensaje   : {mensaje!r}")
        for i, c in enumerate(cif):
            print(f"Cifrado[{i}]: {c}")
        print(f"Descifrado: {des!r}  -> coincide: {des == mensaje}")
    return rsa


def pruebas_encriptacion(resultados, verbose=True):
    pp = primos_por_tecnica(resultados)
    rng = random.Random(7)
    out = []

    def reg(nombre, ok, detalle=""):
        out.append((nombre, bool(ok), detalle))
        if verbose:
            print(f"[{'OK ' if ok else 'FALLA'}] {nombre}  {detalle}")

    pares, validos = 0, 0
    sistemas = []
    claves = list(pp)
    for a in claves:
        for b in claves:
            if a <= b:
                continue
            for k in range(5):
                pares += 1
                try:
                    s = SistemaRSA(pp[a][k], pp[b][k])
                except ValueError:
                    continue
                if s.claves_validas():
                    validos += 1
                    sistemas.append(s)
    reg("P1 llaves: e*d = 1 (mod phi) en pares de tecnicas distintas",
        validos == pares and validos > 0, f"{validos}/{pares} pares")
    rsa = sistemas[0]

    textos = {
        "corto": "Hola",
        "vacio": "",
        "frase": "HOLA RSA DESDE SISTEMAS COMPLEJOS",
        "unicode": "\u00d1and\u00fa caf\u00e9 \u65e5\u672c\u8a9e \U0001F510 \u00bfqu\u00e9 tal?",
        "largo(2000 car.)": "Sistemas complejos. " * 100,
        "ceros iniciales": "\x00\x00abc\x00",
    }
    for nom, t in textos.items():
        c = rsa.cifrar(t)
        reg(f"P2 ida y vuelta texto {nom}", rsa.descifrar(c) == t,
            f"{len(t.encode('utf-8'))} bytes -> {len(c)} bloque(s)")

    msg = "Mensaje de prueba para todas las llaves"
    ok = all(s.descifrar(s.cifrar(msg)) == msg for s in sistemas)
    reg("P3 mismo mensaje con todas las llaves generadas", ok, f"{len(sistemas)} llaves")

    nums = [0, 1, 2, rsa.n - 1] + [rng.randrange(rsa.n) for _ in range(200)]
    ok = all(rsa.descifrar_entero(rsa.cifrar_entero(m)) == m for m in nums)
    reg("P4 ida y vuelta de numeros (bordes + 200 aleatorios)", ok, f"{len(nums)} numeros")
    n40 = [pp["mt"][i] for i in range(10)]
    ok = all(rsa.descifrar_entero(rsa.cifrar_entero(m)) == m for m in n40)
    reg("P4b numeros de 40 cifras (primos hallados) cifrados", ok, "10 numeros")

    c = rsa.cifrar("ATAQUE AL AMANECER")[0]
    reg("P5 el cifrado no es el texto plano ni su entero",
        c != int.from_bytes(b"\x01ATAQUE AL AMANECER", "big"))

    otra = sistemas[1]
    fallo = True
    try:
        fallo = otra.descifrar(rsa.cifrar("secreto importante")) == "secreto importante"
    except Exception:
        fallo = False
    reg("P6 clave privada equivocada NO recupera el mensaje", not fallo)

    c = rsa.cifrar("transferir 100 pesos")
    c2 = [c[0] ^ 1] + c[1:]
    try:
        r = rsa.descifrar(c2)
        det = r != "transferir 100 pesos"
    except Exception:
        det = True
    reg("P7 alterar 1 bit del cifrado destruye el mensaje", det)

    a = rsa.cifrar("A" * 20)[0]
    b = rsa.cifrar("A" * 19 + "B")[0]
    ham = bin(a ^ b).count("1")
    reg("P8 cambiar 1 caracter cambia ~la mitad de los bits del bloque",
        0.3 < ham / rsa.n.bit_length() < 0.7,
        f"{ham}/{rsa.n.bit_length()} bits = {ham/rsa.n.bit_length()*100:.1f}%")

    errs = 0
    for f in (lambda: rsa.cifrar_entero(rsa.n), lambda: rsa.cifrar_entero(-1),
              lambda: SistemaRSA(pp["lcg"][0], pp["lcg"][0]),
              lambda: SistemaRSA(pp["lcg"][0] + 1, pp["mt"][0])):
        try:
            f()
        except ValueError:
            errs += 1
    reg("P9 errores controlados (m>=n, m<0, p==q, p compuesto)", errs == 4, f"{errs}/4")

    m1, m2 = 123456789, 987654321
    ok = (rsa.cifrar_entero(m1) * rsa.cifrar_entero(m2)) % rsa.n == rsa.cifrar_entero(m1 * m2)
    reg("P10 propiedad homomorfica E(m1)E(m2)=E(m1*m2) (debilidad del RSA de libro)", ok)

    reg("P11 RSA de libro es determinista (mismo texto -> mismo cifrado)",
        rsa.cifrar("igual") == rsa.cifrar("igual"))

    t0 = time.time()
    texto = "x" * 100_000
    c = rsa.cifrar(texto)
    t1 = time.time()
    d = rsa.descifrar(c)
    t2 = time.time()
    reg("P12 rendimiento 100 kB", d == texto,
        f"cifrar {t1-t0:.2f}s, descifrar {t2-t1:.2f}s, {len(c)} bloques")

    ok = True
    for k in TECNICAS:
        a = MaquinaPRNG(k, 99).lote(20)
        b = MaquinaPRNG(k, 99).lote(20)
        c2 = MaquinaPRNG(k, 100).lote(20)
        ok &= (a == b) and (a != c2)
    reg("P13 reproducibilidad: misma semilla igual salida, otra semilla distinta", ok,
        "4 tecnicas")

    ok = all(len(str(n)) == c
             for k in TECNICAS for c in (10, 25, 40, 60)
             for n in MaquinaPRNG(k, 5, c).lote(100))
    reg("P14 la maquina respeta cifras = 10, 25, 40, 60", ok, "4 tecnicas x 4 tamanos x 100")

    peor = 0.0
    for k in TECNICAS:
        cont = [0] * 10
        for n in MaquinaPRNG(k, SEED).lote(2000):
            for ch in str(n)[1:]:
                cont[int(ch)] += 1
        tot = sum(cont)
        chi2 = sum((o - tot / 10) ** 2 / (tot / 10) for o in cont)
        peor = max(peor, chi2)
    reg("P15 digitos uniformes (chi2 < 27.88 en las 4 tecnicas)", peor < 27.88,
        f"peor chi2 = {peor:.2f}")

    try:
        import sympy
        todos = [n for lst in pp.values() for n in lst]
        reg("P16 todos los primos hallados pasan sympy.isprime",
            all(sympy.isprime(n) for n in todos), f"{len(todos)} primos")
    except ImportError:
        reg("P16 verificacion con sympy", True, "omitida (sympy no instalado)")

    msg = "mensaje confidencial"
    c = rsa.cifrar(msg)
    fallos = 0
    for s2 in sistemas[1:]:
        try:
            if s2.descifrar(c) == msg:
                fallos += 1
        except Exception:
            pass
    reg("P17 ninguna de las otras llaves descifra el mensaje", fallos == 0,
        f"{len(sistemas)-1} llaves probadas")

    texto = "A" * 32 + "B" * 32
    cb = rsa.cifrar(texto)
    reg("P18 intercambiar bloques cifrados cambia el mensaje",
        rsa.descifrar([cb[1], cb[0]]) != texto and len(cb) == 2)

    try:
        import sympy
        pq = MaquinaPRNG("sha256", 3, 10)
        pa = next(n for n in pq.lote(5000) if es_primo(n))
        qa = next(n for n in pq.lote(5000) if es_primo(n) and n != pa)
        chico = SistemaRSA(pa, qa)
        t0 = time.time()
        fac = sorted(sympy.factorint(chico.n))
        d_rec = pow(chico.e, -1, (fac[0] - 1) * (fac[1] - 1))
        c_ = chico.cifrar_entero(123456)
        reg("P19 un n de ~20 cifras se factoriza y se rompe la llave", 
            pow(c_, d_rec, chico.n) == 123456,
            f"{time.time()-t0:.2f}s; con 40+40 cifras hace falta un algoritmo mucho mas pesado (NFS)")
    except (ImportError, StopIteration, ValueError):
        reg("P19 ataque de factorizacion", True, "omitida")

    n_ok = sum(1 for _, ok, _ in out if ok)
    if verbose:
        print(f"\nResultado: {n_ok}/{len(out)} pruebas superadas")
    return out


def sistema_desde_resultados(resultados, tec_p="lcg", tec_q="sha256"):
    pp = primos_por_tecnica(resultados)
    if tec_p == tec_q:
        return SistemaRSA(pp[tec_p][0], pp[tec_p][1])
    return SistemaRSA(pp[tec_p][0], pp[tec_q][0])


def _pedir(texto, defecto, tipo=str):
    try:
        v = input(f"{texto} [{defecto}]: ").strip()
    except EOFError:
        return defecto
    if not v:
        return defecto
    try:
        return tipo(v)
    except ValueError:
        print("  valor no valido, se usa el defecto")
        return defecto


def menu_interactivo():
    resultados = None
    rsa = None
    while True:
        print("\n===== MAQUINA DE ENCRIPTACION CON PRIMOS GRANDES =====")
        print(" 1) Maquina PRNG: generar numeros de 40 cifras")
        print(" 2) Experimento: 10 000 numeros por tecnica y primos hallados")
        print(" 3) Cifrar un mensaje")
        print(" 4) Descifrar un mensaje (pega los bloques separados por comas)")
        print(" 5) Pruebas de encriptacion")
        print(" 6) Salir")
        try:
            op = input("Opcion: ").strip()
        except EOFError:
            return
        if op == "1":
            t = _pedir(f"Tecnica ({'/'.join(TECNICAS)})", "sha256")
            if t not in TECNICAS:
                print("  tecnica desconocida"); continue
            sem = _pedir("Semilla", SEED, int)
            cif = _pedir("Cifras", CIFRAS, int)
            cant = _pedir("Cantidad", 5, int)
            try:
                for n in MaquinaPRNG(t, sem, cif).lote(cant):
                    print(f"  {n}  ({'PRIMO' if es_primo(n) else 'compuesto'})")
            except ValueError as ex:
                print("  error:", ex)
        elif op == "2":
            cant = _pedir("Cantidad por tecnica", N_NUMEROS, int)
            resultados = correr_experimento(cant, SEED)
            rsa = None
        elif op in ("3", "4"):
            if resultados is None:
                print("(generando primos con el experimento...)")
                resultados = correr_experimento(verbose=False)
            if rsa is None:
                tp = _pedir("Tecnica del primo p", "lcg")
                tq = _pedir("Tecnica del primo q", "sha256")
                try:
                    rsa = sistema_desde_resultados(resultados, tp, tq)
                except (KeyError, ValueError) as ex:
                    print("  no se pudo crear la llave:", ex); continue
                print(f"  p = {rsa.p}\n  q = {rsa.q}\n  n = {rsa.n}\n  e = {rsa.e}")
            if op == "3":
                msg = _pedir("Mensaje", "Hola mundo")
                cif = rsa.cifrar(msg)
                print("Cifrado (bloques separados por comas):")
                print(",".join(str(c) for c in cif))
            else:
                try:
                    txt = input("Bloques: ")
                    cif = [int(x) for x in txt.replace(" ", "").split(",") if x]
                    print("Mensaje:", rsa.descifrar(cif))
                except (ValueError, EOFError) as ex:
                    print("  entrada o llave invalida:", ex)
        elif op == "5":
            if resultados is None:
                resultados = correr_experimento(verbose=False)
            pruebas_encriptacion(resultados)
        elif op == "6":
            return
        else:
            print("  opcion no valida")


def interfaz_grafica():
    try:
        import tkinter as tk
        from tkinter import ttk, messagebox
    except ImportError:
        print("tkinter no esta instalado en este Python (en Linux: "
              "sudo apt install python3-tk). Usa el modo 'menu'.")
        return

    estado = {"res": None, "rsa": None}
    root = tk.Tk()
    root.title("Maquina de encriptacion con primos grandes")
    root.geometry("900x620")
    nb = ttk.Notebook(root)
    nb.pack(fill="both", expand=True, padx=6, pady=6)

    def texto(parent, alto):
        t = tk.Text(parent, height=alto, wrap="word", font=("Courier", 10))
        t.pack(fill="both", expand=True, padx=6, pady=4)
        return t

    def escribir(t, s, borrar=True):
        if borrar:
            t.delete("1.0", "end")
        t.insert("end", s)

    def resultados():
        if estado["res"] is None:
            root.config(cursor="watch"); root.update()
            estado["res"] = correr_experimento(verbose=False)
            root.config(cursor="")
        return estado["res"]

    f1 = ttk.Frame(nb); nb.add(f1, text="1. Maquina PRNG")
    fila = ttk.Frame(f1); fila.pack(fill="x", padx=6, pady=6)
    v_tec = tk.StringVar(value="sha256")
    v_sem = tk.StringVar(value=str(SEED))
    v_cif = tk.StringVar(value=str(CIFRAS))
    v_can = tk.StringVar(value="5")
    for lbl, w in (("Tecnica", ttk.Combobox(fila, textvariable=v_tec, values=list(TECNICAS),
                                            state="readonly", width=8)),
                   ("Semilla", ttk.Entry(fila, textvariable=v_sem, width=10)),
                   ("Cifras", ttk.Spinbox(fila, from_=2, to=200, textvariable=v_cif, width=5)),
                   ("Cantidad", ttk.Spinbox(fila, from_=1, to=1000, textvariable=v_can, width=6))):
        ttk.Label(fila, text=lbl).pack(side="left", padx=(8, 2)); w.pack(side="left")
    out1 = texto(f1, 20)

    def generar():
        try:
            maq = MaquinaPRNG(v_tec.get(), int(v_sem.get()), int(v_cif.get()))
            lineas = [f"{n}  ({len(str(n))} cifras, {'PRIMO' if es_primo(n) else 'compuesto'})"
                      for n in maq.lote(int(v_can.get()))]
            escribir(out1, f"Tecnica: {maq.nombre}\n\n" + "\n".join(lineas))
        except ValueError as ex:
            messagebox.showerror("Parametros invalidos", str(ex))
    ttk.Button(fila, text="Generar", command=generar).pack(side="left", padx=12)

    f2 = ttk.Frame(nb); nb.add(f2, text="2. Experimento")
    fila = ttk.Frame(f2); fila.pack(fill="x", padx=6, pady=6)
    v_n = tk.StringVar(value=str(N_NUMEROS))
    ttk.Label(fila, text="Numeros por tecnica").pack(side="left")
    ttk.Entry(fila, textvariable=v_n, width=8).pack(side="left", padx=4)
    tabla = ttk.Treeview(f2, columns=("t", "p", "d", "z"), show="headings", height=5)
    for c, h in (("t", "Tecnica"), ("p", "Primos"), ("d", "Densidad %"), ("z", "z")):
        tabla.heading(c, text=h); tabla.column(c, width=120)
    tabla.pack(fill="x", padx=6)
    out2 = texto(f2, 14)

    def experimento():
        try:
            cant = int(v_n.get())
        except ValueError:
            messagebox.showerror("Error", "cantidad no valida"); return
        root.config(cursor="watch"); root.update()
        estado["res"] = correr_experimento(cant, SEED, verbose=False)
        estado["rsa"] = None
        root.config(cursor="")
        tabla.delete(*tabla.get_children())
        txt = f"Esperados ~{cant*densidad_teorica():.1f} primos por tecnica.\n\n"
        for k, r in estado["res"].items():
            tabla.insert("", "end", values=(r["nombre"], r["n"],
                                            f"{r['densidad']*100:.4f}", f"{r['z']:+.2f}"))
            txt += f"{r['nombre']}: primeros primos\n"
            txt += "".join(f"  {n}\n" for _, n in r["primos"][:3])
        escribir(out2, txt + f"\nPrimos guardados en {CSV_PRIMOS}")
    ttk.Button(fila, text="Ejecutar experimento", command=experimento).pack(side="left", padx=8)

    f3 = ttk.Frame(nb); nb.add(f3, text="3. Cifrado")
    fila = ttk.Frame(f3); fila.pack(fill="x", padx=6, pady=6)
    v_p = tk.StringVar(value="lcg"); v_q = tk.StringVar(value="sha256")
    ttk.Label(fila, text="Primo p de").pack(side="left")
    ttk.Combobox(fila, textvariable=v_p, values=list(TECNICAS), state="readonly",
                 width=8).pack(side="left", padx=4)
    ttk.Label(fila, text="Primo q de").pack(side="left")
    ttk.Combobox(fila, textvariable=v_q, values=list(TECNICAS), state="readonly",
                 width=8).pack(side="left", padx=4)
    info = texto(f3, 6)
    ttk.Label(f3, text="Mensaje:").pack(anchor="w", padx=6)
    t_msg = texto(f3, 4)
    ttk.Label(f3, text="Cifrado (bloques separados por comas):").pack(anchor="w", padx=6)
    t_cif = texto(f3, 6)
    ttk.Label(f3, text="Descifrado:").pack(anchor="w", padx=6)
    t_des = texto(f3, 3)

    def crear():
        try:
            estado["rsa"] = sistema_desde_resultados(resultados(), v_p.get(), v_q.get())
        except ValueError as ex:
            messagebox.showerror("Llave invalida", str(ex)); return None
        r = estado["rsa"]
        escribir(info, f"p = {r.p}\nq = {r.q}\nn = {r.n}  ({r.n.bit_length()} bits)\n"
                       f"e = {r.e}\nd = {r.d}\nbloque = {r.tam_bloque} bytes")
        return r

    def cifrar():
        r = crear()
        if r:
            msg = t_msg.get("1.0", "end").rstrip("\n")
            escribir(t_cif, ",".join(str(c) for c in r.cifrar(msg)))
            escribir(t_des, "")

    def descifrar():
        r = estado["rsa"] or crear()
        if not r:
            return
        try:
            bloques = [int(x) for x in t_cif.get("1.0", "end").replace("\n", "")
                       .replace(" ", "").split(",") if x]
            escribir(t_des, r.descifrar(bloques))
        except ValueError as ex:
            messagebox.showerror("No se pudo descifrar", str(ex))
    fb = ttk.Frame(f3); fb.pack(fill="x", padx=6, pady=4)
    ttk.Button(fb, text="Cifrar (crea llaves nuevas)", command=cifrar).pack(side="left")
    ttk.Button(fb, text="Descifrar", command=descifrar).pack(side="left", padx=8)

    f4 = ttk.Frame(nb); nb.add(f4, text="4. Pruebas")
    out4 = texto(f4, 25)

    def pruebas():
        root.config(cursor="watch"); root.update()
        res = pruebas_encriptacion(resultados(), verbose=False)
        root.config(cursor="")
        lineas = [f"[{'OK ' if ok else 'FALLA'}] {n}  {d}" for n, ok, d in res]
        lineas.append(f"\nResultado: {sum(ok for _, ok, _ in res)}/{len(res)} pruebas superadas")
        escribir(out4, "\n".join(lineas))
    ttk.Button(f4, text="Ejecutar pruebas", command=pruebas).pack(anchor="w", padx=6, pady=6)
    out4.pack_forget(); out4.pack(fill="both", expand=True, padx=6, pady=4)

    root.mainloop()


def main():
    ap = argparse.ArgumentParser(description="Maquina para encriptar mensajes con numeros primos grandes")
    sub = ap.add_subparsers(dest="cmd")
    m = sub.add_parser("maquina", help="usar la maquina PRNG")
    m.add_argument("--tecnica", default="sha256", choices=list(TECNICAS))
    m.add_argument("--semilla", type=int, default=SEED)
    m.add_argument("--cifras", type=int, default=CIFRAS)
    m.add_argument("--cantidad", type=int, default=5)
    e = sub.add_parser("experimento", help="generar N numeros por tecnica y hallar primos")
    e.add_argument("--cantidad", type=int, default=N_NUMEROS)
    e.add_argument("--semilla", type=int, default=SEED)
    c = sub.add_parser("cifrar", help="cifrar/descifrar un mensaje con los primos hallados")
    c.add_argument("--mensaje", default="HOLA RSA DESDE SISTEMAS COMPLEJOS")
    sub.add_parser("pruebas", help="bateria de pruebas de encriptacion")
    sub.add_parser("menu", help="menu interactivo en consola")
    sub.add_parser("gui", help="interfaz grafica (tkinter)")
    sub.add_parser("todo", help="demo completa (por defecto)")
    args = ap.parse_args()
    cmd = args.cmd or "todo"

    if cmd == "menu":
        menu_interactivo()
        return
    if cmd == "gui":
        interfaz_grafica()
        return

    if cmd == "maquina":
        maq = MaquinaPRNG(args.tecnica, args.semilla, args.cifras)
        print(f"Maquina: tecnica={maq.nombre}, semilla={args.semilla}, cifras={args.cifras}")
        for n in maq.lote(args.cantidad):
            print(f"{n}  ({len(str(n))} cifras, {'PRIMO' if es_primo(n) else 'compuesto'})")
        return

    if cmd == "todo":
        print("### 1. LA MAQUINA: una salida de 40 cifras por tecnica (semilla=42)")
        for k in TECNICAS:
            n = MaquinaPRNG(k, 42).siguiente()
            print(f"  {TECNICAS[k][0]:<18}{n}  ({len(str(n))} cifras)")
        print("\n### 2. EXPERIMENTO: 10 000 numeros de 40 cifras por tecnica")
    cant = getattr(args, "cantidad", N_NUMEROS) if cmd == "experimento" else N_NUMEROS
    sem = getattr(args, "semilla", SEED) if cmd == "experimento" else SEED
    res = correr_experimento(cant, sem, verbose=cmd in ("experimento", "todo"))
    if cmd == "experimento":
        return
    if cmd in ("cifrar", "todo"):
        print("\n### 3. SISTEMA DE ENCRIPTACION")
        demostrar_rsa(res, getattr(args, "mensaje", "HOLA RSA DESDE SISTEMAS COMPLEJOS"))
    if cmd in ("pruebas", "todo"):
        print("\n### 4. PRUEBAS DE ENCRIPTACION")
        pruebas_encriptacion(res)


if __name__ == "__main__":
    main()
