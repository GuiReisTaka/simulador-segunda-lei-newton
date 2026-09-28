"""
fisica.py — fórmulas fechadas dos sistemas simulados (ADO 2 — Leis de Newton).

Este módulo só faz contas com NumPy: não desenha nada e não importa matplotlib.
Assim dá para testar a física separada da interface (ver teste_fisica.py).

Sistemas implementados:
  - 3.2  Força centrípeta (movimento circular uniforme)
  - 3.3  Força de arrasto (linear e quadrático), comparada com a queda livre

Nenhuma integração numérica: todas as funções são fórmulas fechadas, avaliadas
diretamente sobre um vetor de tempo t (np.linspace) ou sobre um escalar.
"""

import numpy as np

G_TERRA = 9.81  # m/s²


# ---------------------------------------------------------------------------
# 3.2 — FORÇA CENTRÍPETA (MCU)
# ---------------------------------------------------------------------------
# Segunda Lei de Newton na direção radial (apontando para o centro):
#     ΣF = m·a_c   com   a_c = v²/R   ⇒   Fc = m·v²/R
# A força centrípeta não é uma força "nova": é a resultante que aponta para o
# centro (aqui, a tração do fio).

def velocidade_angular(v, R):
    """ω = v/R  (definição de velocidade angular no MCU: v = ω·R)."""
    return v / R


def periodo(v, R):
    """T = 2πR/v  (uma volta completa: distância 2πR na velocidade v)."""
    return 2 * np.pi * R / v


def aceleracao_centripeta(v, R):
    """a_c = v²/R  (variação da direção da velocidade; aponta para o centro)."""
    return v**2 / R


def forca_centripeta(m, v, R):
    """Fc = m·v²/R  (2ª Lei de Newton com a_c = v²/R)."""
    return m * v**2 / R


def forca_centripeta_omega(m, w, R):
    """Fc = m·ω²·R  (mesma força, escrita com ω; usa v = ωR)."""
    return m * w**2 * R


def posicao_circular(t, R, w):
    """
    x(t) = R·cos(ωt),  y(t) = R·sin(ωt)
    Posição no MCU: o ângulo cresce linearmente, θ = ωt (ω constante).
    Aceita t escalar ou vetor.
    """
    return R * np.cos(w * t), R * np.sin(w * t)


# ---------------------------------------------------------------------------
# 3.3 — FORÇA DE ARRASTO
# ---------------------------------------------------------------------------
# Queda vertical, eixo y para baixo. 2ª Lei de Newton:
#     m·dv/dt = m·g − F_arrasto
# Velocidade terminal: quando a aceleração zera (F_arrasto = m·g).
# As soluções abaixo são as fórmulas fechadas dadas no enunciado.

# ---- Queda livre (sem arrasto), para comparação ---------------------------

def v_queda_livre(t, g=G_TERRA):
    """v(t) = g·t  (MRUV a partir do repouso)."""
    return g * t


def y_queda_livre(t, g=G_TERRA):
    """y(t) = ½·g·t²  (MRUV a partir do repouso)."""
    return 0.5 * g * t**2


# ---- Arrasto linear: F = −b·v ---------------------------------------------
# m·dv/dt = m·g − b·v  ⇒  v_t = m·g/b,  τ = m/b

def v_terminal_linear(m, b, g=G_TERRA):
    """v_terminal = m·g/b  (aceleração nula: b·v = m·g)."""
    return m * g / b


def tau_linear(m, b):
    """τ = m/b  (constante de tempo da exponencial)."""
    return m / b


def v_linear(t, m, b, g=G_TERRA):
    """v(t) = v_t·(1 − e^(−t/τ))  (solução da EDO linear de 1ª ordem)."""
    vt = v_terminal_linear(m, b, g)
    return vt * (1 - np.exp(-t / tau_linear(m, b)))


def y_linear(t, m, b, g=G_TERRA):
    """y(t) = v_t·t − v_t·τ·(1 − e^(−t/τ))  (integral de v(t) de 0 a t)."""
    vt = v_terminal_linear(m, b, g)
    tau = tau_linear(m, b)
    return vt * t - vt * tau * (1 - np.exp(-t / tau))


# ---- Arrasto quadrático: F = −c·v·|v| --------------------------------------
# m·dv/dt = m·g − c·v²  ⇒  v_t = √(m·g/c)

def v_terminal_quadratico(m, c, g=G_TERRA):
    """v_terminal = √(m·g/c)  (aceleração nula: c·v² = m·g)."""
    return np.sqrt(m * g / c)


def v_quadratico(t, m, c, g=G_TERRA):
    """v(t) = v_t·tanh(g·t/v_t)  (solução da EDO com termo em v²)."""
    vt = v_terminal_quadratico(m, c, g)
    return vt * np.tanh(g * t / vt)


def _log_cosh(x):
    """
    ln(cosh x) numericamente estável: para |x| grande, cosh(x) estoura.
    Identidade: ln cosh x = |x| + ln(1 + e^(−2|x|)) − ln 2.
    """
    ax = np.abs(x)
    return ax + np.log1p(np.exp(-2 * ax)) - np.log(2)


def y_quadratico(t, m, c, g=G_TERRA):
    """y(t) = (v_t²/g)·ln(cosh(g·t/v_t))  (integral de v(t) de 0 a t)."""
    vt = v_terminal_quadratico(m, c, g)
    return (vt**2 / g) * _log_cosh(g * t / vt)


# ---- Tempo para atingir uma fração da velocidade terminal ------------------

def tempo_fracao_linear(m, b, fracao=0.95):
    """
    Resolve 1 − e^(−t/τ) = fração  ⇒  t = −τ·ln(1 − fração).
    Para 95%: t ≈ 3·τ.
    """
    return -tau_linear(m, b) * np.log(1 - fracao)


def tempo_fracao_quadratico(m, c, fracao=0.95, g=G_TERRA):
    """
    Resolve tanh(g·t/v_t) = fração  ⇒  t = (v_t/g)·artanh(fração).
    Para 95%: t ≈ 1,83·v_t/g.
    """
    vt = v_terminal_quadratico(m, c, g)
    return (vt / g) * np.arctanh(fracao)


# ---- Utilitário: coeficiente quadrático equivalente ------------------------

def c_equivalente(m, b, g=G_TERRA):
    """
    Valor de c que dá a MESMA velocidade terminal que o arrasto linear b:
        √(m·g/c) = m·g/b  ⇒  c = b²/(m·g)
    b [kg/s] e c [kg/m] têm unidades diferentes; comparar os dois modelos com o
    "mesmo número" nos dois coeficientes não compara coisas equivalentes.
    """
    return b**2 / (m * g)
