"""
teste_fisica.py — validação das fórmulas de fisica.py, sem abrir nenhuma janela.

A ideia principal: conferir que as fórmulas fechadas realmente satisfazem a
equação de movimento (2ª Lei de Newton). Para isso comparamos a derivada
numérica da solução com o lado direito da equação. (Isso é só um teste: a
simulação em si não usa nenhuma integração numérica.)
"""

import numpy as np

import fisica as f

g = f.G_TERRA
falhas = 0


def confere(nome, condicao):
    global falhas
    print(("OK    " if condicao else "FALHA ") + nome)
    if not condicao:
        falhas += 1


def derivada(func, t, h=1e-6):
    """Derivada central de func em t (usada apenas nos testes)."""
    return (func(t + h) - func(t - h)) / (2 * h)


# --------------------------------------------------------------------------
# 3.2 — Força centrípeta
# --------------------------------------------------------------------------
m, R, v = 1.5, 4.0, 10.0
w = f.velocidade_angular(v, R)

confere("ω = v/R = 2,5 rad/s", np.isclose(w, 2.5))
confere("Fc = m·v²/R = 37,5 N", np.isclose(f.forca_centripeta(m, v, R), 37.5))
confere("m·v²/R == m·ω²·R", np.isclose(f.forca_centripeta(m, v, R),
                                        f.forca_centripeta_omega(m, w, R)))
confere("T = 2π/ω", np.isclose(f.periodo(v, R), 2 * np.pi / w))

t = np.linspace(0, f.periodo(v, R), 1000)
x, y = f.posicao_circular(t, R, w)
confere("trajetória: x² + y² = R²", np.allclose(x**2 + y**2, R**2))
confere("após um período volta ao ponto inicial",
        np.isclose(x[-1], x[0]) and np.isclose(y[-1], y[0], atol=1e-9))

# Aceleração vetorial: a = −ω²·r  (aponta para o centro, módulo v²/R)
t0 = 0.7
ax_ = derivada(lambda s: derivada(lambda u: f.posicao_circular(u, R, w)[0], s, 1e-4), t0, 1e-4)
ay_ = derivada(lambda s: derivada(lambda u: f.posicao_circular(u, R, w)[1], s, 1e-4), t0, 1e-4)
x0, y0 = f.posicao_circular(t0, R, w)
confere("aceleração = −ω²·posição (aponta para o centro)",
        np.isclose(ax_, -w**2 * x0, atol=1e-3) and np.isclose(ay_, -w**2 * y0, atol=1e-3))
confere("|a| = v²/R", np.isclose(np.hypot(ax_, ay_), v**2 / R, atol=1e-3))

# --------------------------------------------------------------------------
# 3.3 — Arrasto linear:  m·dv/dt = m·g − b·v
# --------------------------------------------------------------------------
m, b = 2.0, 1.0
vt = f.v_terminal_linear(m, b)
tau = f.tau_linear(m, b)

confere("linear: v_terminal = m·g/b", np.isclose(vt, m * g / b))
confere("linear: v(0) = 0 e y(0) = 0",
        f.v_linear(0.0, m, b) == 0 and f.y_linear(0.0, m, b) == 0)
confere("linear: v(t → ∞) → v_terminal", np.isclose(f.v_linear(50 * tau, m, b), vt))

ts = np.array([0.1, 0.5, 1.0, 2.5, 5.0, 10.0])
dv = derivada(lambda s: f.v_linear(s, m, b), ts)
confere("linear: m·dv/dt = m·g − b·v (satisfaz a 2ª Lei)",
        np.allclose(m * dv, m * g - b * f.v_linear(ts, m, b), atol=1e-4))
dy = derivada(lambda s: f.y_linear(s, m, b), ts)
confere("linear: dy/dt = v", np.allclose(dy, f.v_linear(ts, m, b), atol=1e-4))

# --------------------------------------------------------------------------
# 3.3 — Arrasto quadrático:  m·dv/dt = m·g − c·v²
# --------------------------------------------------------------------------
c = 0.05
vt_q = f.v_terminal_quadratico(m, c)

confere("quadrático: v_terminal = √(m·g/c)", np.isclose(vt_q, np.sqrt(m * g / c)))
confere("quadrático: v(0) = 0 e y(0) = 0",
        f.v_quadratico(0.0, m, c) == 0 and f.y_quadratico(0.0, m, c) == 0)
confere("quadrático: v(t → ∞) → v_terminal",
        np.isclose(f.v_quadratico(50 * vt_q / g, m, c), vt_q))

dv = derivada(lambda s: f.v_quadratico(s, m, c), ts)
confere("quadrático: m·dv/dt = m·g − c·v² (satisfaz a 2ª Lei)",
        np.allclose(m * dv, m * g - c * f.v_quadratico(ts, m, c)**2, atol=1e-4))
dy = derivada(lambda s: f.y_quadratico(s, m, c), ts)
confere("quadrático: dy/dt = v", np.allclose(dy, f.v_quadratico(ts, m, c), atol=1e-4))
confere("quadrático: y(t) não estoura para t muito grande",
        np.isfinite(f.y_quadratico(1e6, m, c)))

# --------------------------------------------------------------------------
# Comparações com a queda livre
# --------------------------------------------------------------------------
tp = 1e-4
confere("t pequeno: arrasto ≈ queda livre (linear)",
        np.isclose(f.v_linear(tp, m, b), f.v_queda_livre(tp), rtol=1e-3))
confere("t pequeno: arrasto ≈ queda livre (quadrático)",
        np.isclose(f.v_quadratico(tp, m, c), f.v_queda_livre(tp), rtol=1e-3))

tl = np.linspace(0.1, 30, 300)
dif = f.v_queda_livre(tl) - f.v_linear(tl, m, b)
confere("diferença (queda livre − arrasto) cresce com o tempo", np.all(np.diff(dif) > 0))

# --------------------------------------------------------------------------
# Tempo até 95% da velocidade terminal
# --------------------------------------------------------------------------
t95_l = f.tempo_fracao_linear(m, b)
t95_q = f.tempo_fracao_quadratico(m, c)
confere("linear: v(t95) = 0,95·v_terminal",
        np.isclose(f.v_linear(t95_l, m, b), 0.95 * vt))
confere("quadrático: v(t95) = 0,95·v_terminal",
        np.isclose(f.v_quadratico(t95_q, m, c), 0.95 * vt_q))
confere("linear: t95 ≈ 3·τ (2,996·τ)", np.isclose(t95_l / tau, 2.9957, atol=1e-3))

# Com a MESMA velocidade terminal, τ = v_t/g nos dois casos
c_eq = f.c_equivalente(m, b)
confere("c_equivalente dá a mesma v_terminal",
        np.isclose(f.v_terminal_quadratico(m, c_eq), vt))
razao = f.tempo_fracao_quadratico(m, c_eq) / f.tempo_fracao_linear(m, b)
confere("mesma v_t: t95 quadrático / linear ≈ 1,83/3,00 = 0,61",
        np.isclose(razao, 1.8318 / 2.9957, atol=1e-3))

print()
print("Todos os testes passaram." if falhas == 0 else f"{falhas} teste(s) falharam.")
raise SystemExit(0 if falhas == 0 else 1)
