# Máquina para encriptar mensajes con números primos grandes

Taller de Sistemas Complejos — Universidad Sergio Arboleda.
Autor: Santiago Alejandro Céspedes Daza.

Programa único (`generadores_primos.py`) que:

1. **Máquina PRNG parametrizable** (técnica, semilla, cifras, cantidad) que genera números pseudoaleatorios de 40 cifras con 4 técnicas: LCG (módulo primo 2^61−1), Mersenne Twister, SHA-256 en modo contador y Blum-Blum-Shub.
2. **Experimento**: 10 000 números de 40 cifras por técnica y búsqueda de los primos (resultado en `primos_encontrados.csv`).
3. **Sistema de encriptación RSA por bloques** construido con los primos encontrados (mensajes UTF-8 de cualquier largo y números).
4. **25 pruebas de encriptación** automáticas.
5. **Menú interactivo** (se abre al ejecutar el programa) con una opción para abrir la **interfaz gráfica** (tkinter).

## Instalación

```bash
pip install -r requirements.txt
```

Si `gmpy2` no se instala, se puede quitar de `requirements.txt`; el programa usa `sympy` como respaldo (más lento).

## Uso

```bash
python generadores_primos.py
```

Abre un **menú interactivo** con estas opciones:

1. Máquina PRNG: generar números de 40 cifras (eliges técnica, semilla, cifras y cantidad).
2. Experimento: 10 000 números por técnica y primos hallados.
3. Cifrar un mensaje.
4. Descifrar un mensaje (pegas los bloques cifrados separados por comas).
5. Pruebas de encriptación (25 pruebas).
6. Demo completa.
7. Abrir la interfaz gráfica.
8. Salir.

También se puede ejecutar cada función directamente:

```bash
python generadores_primos.py gui
python generadores_primos.py pruebas
python generadores_primos.py todo
python generadores_primos.py maquina --tecnica sha256 --semilla 42 --cifras 40 --cantidad 5
python generadores_primos.py experimento --cantidad 10000
python generadores_primos.py cifrar --mensaje "hola mundo"
```

## Resultados (semilla 2026, 10 000 números por técnica)

| Técnica | Primos |
|---|---|
| LCG | 106 |
| Mersenne Twister | 116 |
| SHA-256 | 87 |
| Blum-Blum-Shub | 114 |
| Esperado (TNP) | 109.3 |

## Informe

`informe_prng_primos_ieee.pdf` (formato IEEE, hecho en LaTeX; fuentes en `informe_latex/`).

## Aviso de seguridad

Con primos de 40 cifras, el módulo RSA tiene ~266 bits y se puede factorizar. El sistema es **didáctico**, no apto para proteger datos reales.
