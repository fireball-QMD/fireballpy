# Rama `stationaryQ_GSN`: correcciones GSN al sistema de cargas estacionarias

Rama creada el 2026-06-25 desde `stationaryQ` (`ebd92bc`, "Paso 2"), inicialmente como `D25`.
Referencia teórica: `XC_basic_claude_fable.pdf` (secciones III.D *full system* y III.E
*simplified version*; ecuaciones 117 y 119, respuesta de banda ec. 91).

**Resultado principal**: `stationary_charges` pasa de ser el método menos consistente E↔F a
ser efectivamente variacional. Gap entre el mínimo de energía y el cero de fuerza en H2O:
de ~1.2 mÅ (sin GSN) a **0.01–0.22 mÅ** según la fdata, y funciona también en **base doble**
(donde antes divergía) con gaps ≤ 0.18 mÅ en H2O, dímero de agua y H2O2.

## 1. Punto de partida y diagnóstico

Al crear la rama, `stationaryQ` ya tenía el modelo de Hamiltoniano de cargas estáticas
`H = h⁰ + Σ_α g^α (Q_α − Q⁰_α)` (Coulomb `g_h` + XC `g_xc`, índice de shell global en
primera posición) y el ensamblado único en Kscf=1.

El problema: las fuerzas son la derivada exacta de la energía **a carga fija**; el término
no variacional ∂E/∂Q·∂Q/∂R se omite. `stationary_charges` debería anularlo (∂E/∂Q = 0) pero
su sistema estacionario solo incluía el término local de v^xc y **le faltaba la contribución
GSN/SNXC** (densidad media n̄) — el propio doc lo marcaba con "[REVISAR!: GSN]".

Un experimento previo con Hessiano **numérico** (diferencias finitas sobre las rutinas de
energía, proyecto 16-claude) confirmó el diagnóstico y dejó la lección clave: corregir solo
el vector B no sirve (se lava en la autoconsistencia por la ligadura ΣQ = Ne); lo decisivo
es la **matriz M** (Hessiano del doble-conteo XC) y la respuesta de banda. Esta rama
implementa esa corrección de forma **analítica**.

## 2. Qué se ha implementado

### 2.1 Matrices de densidad `den_or` / `den_sh` (`785d2a5`)

- `M_system.f90` declara y `allocate_system.f90` aloca
  `den_or(γ,µ,ν,ineigh,iatom)` = N^γ_µν = ⟨φµ|(φ^S_γ)²|φν⟩ (ec. 17) y
  `den_sh(γ,α,β,ineigh,iatom)` = M^γ_αβ = ⟨φ^S_α|(φ^S_γ)²|φ^S_β⟩ (ec. 21),
  con γ = índice de shell global (orden de `get_shell_ofatom_issh`).
- Se rellenan en `average_ca_rho.f90` **solo en Kscf=1** (son Q-independientes),
  verificadas contra el ensamblado vivo a ~1e-16 (flag `debug_den`).
- `average_ca_rho.f90` guarda además el solape esférico on-site crudo en el slot self
  de `sm_mat` (Kscf=1).

### 2.2 Correcciones GSN en `stationary_charges.f90` (`load_M`, flag `igsn`)

| igsn | qué añade | gap H2O (superspline) |
|---|---|---|
| 0 | sin GSN | 1.21 mÅ / 0.055 eV/Å |
| 1 | ec. 117 (simplificado, congelado en Q^R) | 3.22 mÅ / 0.35 |
| 2 | ec. 117 + Δα (ec. 119) | 1.83 mÅ / 0.43 |
| 3 | ec. 117 − Δα | 5.16 mÅ / 0.30 |
| **4** | **sistema completo linealizado (III.D) — DEFAULT** | **0.22 mÅ / 0.012** |

`igsn=4`: constantes GSN vivas (coeficientes evaluados en Qin) + respuesta GSN en la matriz
(`gsnM`) + respuesta de banda ec. 91 (`gsnB`, construida con `den_or`/`den_sh`). Al converger
el SCF resuelve la estacionariedad exacta de la energía viva. La matriz queda asimétrica →
`dgesv` tal cual (sin simetrizar).

Con fdata **sin** superspline el gap de H2O baja a **0.00 mÅ / 0.0017 eV/Å** (variacional
exacto): el residuo de 0.22 mÅ es artefacto de la interpolación superspline, no del método.

### 2.3 Shells de polarización y bases con L repetido (`ifix_shells`, `isolver`)

El dímero de agua descubrió que las shells **d** (polarización vacía, Q⁰=0) crean direcciones
casi nulas del sistema estacionario (cond(M) ≈ 3.4·10³): la solución se desliza hasta
|Q_d| ~ 1, donde la linealización 1c de las tablas no vale (Q_d = ±1, no físico, SCF a
15 pasos). En base doble pasa lo mismo con las shells excitadas s*/p* libres (SCF no
converge, Q_p(O) = −3.6).

- **`ifix_shells=2` / `fix_shells='auto'` (DEFAULT)**: fija a carga neutra **todas** las
  shells con Qneutral=0 (d y excitadas), vía `fix_shell_charge`. (`=1`/`'d'` fija solo las d;
  `=0`/`'none'` ninguna; `=3` máscara de usuario.) Desde 2026-07-15 es variable de
  `M_system.f90` controlable desde Python: `Fireball(..., fix_shells='auto'|'d'|'none'|máscara)`,
  donde la máscara es 0/1 por shell global (orden átomo1-shells, átomo2-shells…). **Ojo**:
  `'auto'` no fija nada en fdatas sin shells vacías (todo s,p ocupadas) — ahí considerar
  máscara manual o `isolver=2`. A/B de patrones de fijado en el cuaderno
  (`IMPLEMENTACION_eq117_simplified.md` §5e): cualquier fijado consistente mantiene el gap
  E↔F ~0.1 mÅ, pero solo el criterio Qneutral=0 completo da una PES física.
- **`isolver`**: 0 = `dgesv` directo (DEFAULT); 1 = `dgelsd`/truncamiento SVD (malo para
  scans: saltos de rango → PES discontinua); 2 = Tikhonov suave en dQ = Q − Q⁰ con ligadura
  exacta y µ sin penalizar (ridge ~0.03). Tikhonov acota las cargas pero no sustituye a
  `ifix_dshell`; en la práctica **no ha hecho falta** ni en base doble.

### 2.4 Bugs corregidos

- `buildspline_1d.f90`: leía `z(numz_used+1)` fuera de rango en el último punto del spline
  (crash con `-check bounds` al cargar fdata; `785d2a5`).
- `stationary_charges.f90`: `beta_iatom` sin inicializar en la rama de shells fijas del
  vector B (`785d2a5`).
- `stationary_charges.f90`: `mapindex` alocado con `nssh_tot2` (shells libres) pero escrito
  para todas las shells 1..`nssh_tot` → corrupción de heap → segfault aleatorio en scans
  (afectaba a cualquier uso de shells fijas; los resultados no cambian porque los elementos
  fuera de rango nunca se leían) (`7214b6b`).

## 3. Resultados

Métrica: scan de un átomo a lo largo de una dirección; gap = |λ(E_min) − λ(F=0)| y residuo
max|F + dE/dλ| en la ventana fina (paso 5 mÅ). Todo con `igsn=4` + `ifix_dshell=2` + `dgesv`.

### Base simple (`fdata`: H s · O s,p,d)

| sistema | gap E-min↔F=0 | max\|F+dE/dλ\| |
|---|---|---|
| H2O | 0.01 mÅ | 0.0018 eV/Å |
| dímero de agua | 0.03 mÅ | 0.0041 eV/Å |
| H2O2 | 0.16 mÅ | 0.27 (*) |

### Base doble (`fdata_doble`: H s,s* · O s,p,s*,p*)

| sistema | gap E-min↔F=0 | max\|F+dE/dλ\| |
|---|---|---|
| H2O | 0.11 mÅ | 0.0041 eV/Å |
| dímero de agua | 0.06 mÅ | 0.0011 eV/Å |
| H2O2 | 0.18 mÅ | 0.034 (*) |

(*) Los residuos de H2O2 son dos kinks localizados de ~2 meV en E(λ), comunes a **todos**
los métodos de carga (mulliken_dipole los sufre igual o peor) → artefacto de la PES
(¿cruce de niveles/ocupaciones, nodos de tablas?), no del sistema estacionario.

Referencia `mulliken_dipole`: ~3 mÅ / 0.15 eV/Å en base simple; en base doble **~56 mÅ**
(E_min en λ = +0.040 y F=0 en λ = −0.016 para H2O). El baseline previo a la rama daba
`stationary_charges` en base doble directamente divergente (95.9 mÅ o explosión del SCF).

Cargas convergidas (base doble, en el mínimo): H2O O = −0.298 / H = +0.149; dímero
O = −0.31/−0.28; H2O2 O = −0.120 / H = +0.120. Físicas y con las excitadas a 0 por
construcción.

### Validación en dinámica: NVE del dímero de agua

NVE con Verlet (ASE), 300 K inicial, dt = 0.5 fs, 0.5 ps, misma semilla de velocidades en
todos los métodos. Drift lineal de E_tot / fluctuación RMS respecto a la recta:

| método | base doble: drift / rms | base simple: drift / rms |
|---|---|---|
| **stationary_charges (igsn=4)** | **+2.2 meV/ps / 0.8 meV** | **+20 meV/ps / 4.4 meV** |
| lowdin | +139 / 35 | −55 / 30 |
| weighted_lowdin | +1640 / 118 | −45 / 30 |
| mulliken | +736 / 36 | +488 / 142 |
| mulliken_dipole | +3536 / 657 | −27 / 23 |
| mulliken_dipole_preserving | +3525 / 621 | +159 / 32 |

stationary_charges es el mejor en ambas bases y el único utilizable en base doble: los
métodos con dipolo bombean ~3.5 eV/ps (el término no variacional ∂E/∂Q·∂Q/∂R que sus fuerzas
omiten) y el sistema se autocalienta de 255 K a >500-900 K en 0.5 ps. El drift residual de
stationary_charges no escala con dt (control con dt = 0.25 fs: +3.4 meV/ps) → no es error del
integrador sino la pequeña inconsistencia de fuerzas restante (~0.001-0.004 eV/Å de los
scans). El ranking NVE reproduce el de los gaps E↔F estáticos. Script, datos crudos y figura
E_tot(t): `18_claude/md_nve.py`, `18_claude/nve_data/`, `18_claude/doc/nve_etot_dimer.png`.

## 4. Hallazgos colaterales y pendientes

- Solape esférico on-site tabulado S_αα = 1.0024–1.0028 en vez de 1 (artefacto de tabla en
  y=0 de fdata_superspline) → sesgo ~0.25% en `arho_on`. Revisar en create.
- Kinks ~2 meV en la PES de H2O2, comunes a todos los métodos — investigar ocupaciones o
  interpolación de tablas.
- Residuo restante ~0.004–0.012 eV/Å: candidatos — inconsistencias de tablas (S_αα≠1) y el
  término ∂E/∂Q·∂Q/∂R residual de la linealización 1c que las fuerzas no incluyen.
- Coste: `load_M` reconstruye gsnM/gsnB cada iteración SCF; optimizable si hiciera falta.

## 5. Cómo reproducir

```bash
cd fireballpy && python install.py --intel --fast
# scan E↔F (h2o | h2o2 | dimer), script en 18_claude/scan_generic.py:
python ../scan_generic.py <fdata_path> stationary_charges h2o
```

`scan_generic.py`: malla gruesa ±0.30 Å (paso 0.05) para acotar el mínimo + malla fina
±0.05 (paso 0.005); E_min por parábola local, F=0 por interpolación lineal. Los flags
(`igsn`, `isolver`, `ridge_sc`) son `parameter` en
`src/fireball/SYSTEM/stationary_charges.f90`; `ifix_shells` es variable de `M_system.f90`
y se controla desde Python con el parámetro `fix_shells` de la calculadora (§2.3).

**OJO**: si se cambia `makemunu` hay que regenerar la fdata (los ficheros de interacción
cambian de tamaño; una fdata vieja da `forrtl: severe (59)` al leer `vxc_2c`).

Cuaderno de laboratorio con el detalle de cada sesión: `18_claude/doc/`
(`IMPLEMENTACION_eq117_simplified.md`, `IMPLEMENTACION_fix_gsn.md`, `baseline_EvsF.md`).
