"""
interface.py — simulador interativo (ADO 2 — Leis de Newton).

Sistemas: (3.2) força centrípeta e (3.3) força de arrasto.
Toda a física está em fisica.py; aqui só há a parte gráfica (matplotlib.widgets).

Rodar:  python interface.py
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Button, RadioButtons, Slider

import fisica as f

SISTEMAS = ["Força centrípeta (3.2)", "Força de arrasto (3.3)"]
MODOS = ["Linear:  F = −b·v", "Quadrático:  F = −c·v²"]

INTERVALO_MS = 40                 # tempo entre quadros da animação
DT_ANIM = INTERVALO_MS / 1000     # segundos de "tempo físico" por quadro (tempo real)
N_PONTOS = 500                    # tamanho do vetor de tempo do sistema 3.3


def novo_radio(ax, rotulos):
    """
    RadioButtons sem blit. Com blit (padrão do matplotlib >= 3.7), os círculos são
    redesenhados a cada quadro mesmo quando o eixo está escondido (set_visible(False)),
    e ficavam "soltos" na tela ao trocar de sistema.
    """
    try:
        return RadioButtons(ax, rotulos, active=0, useblit=False)
    except TypeError:  # matplotlib antigo não tem o parâmetro useblit (e não tem o problema)
        return RadioButtons(ax, rotulos, active=0)


def criar_interface():
    fig = plt.figure(figsize=(12, 7.5))
    fig.suptitle("ADO 2 — Leis de Newton: simulador interativo", fontsize=13)
    try:
        fig.canvas.manager.set_window_title("Simulador — Leis de Newton")
    except Exception:
        pass  # alguns backends (ex.: notebook) não têm esse método

    # Estado mínimo compartilhado entre os callbacks
    estado = {"sistema": "circular", "modo": "linear",
              "t": 0.0, "w": None, "R": 5.0}

    # ------------------------------------------------------------------
    # Painel esquerdo: menu, modo do arrasto, botão e texto com resultados
    # ------------------------------------------------------------------
    ax_menu = fig.add_axes([0.02, 0.80, 0.25, 0.12])
    ax_menu.set_title("Sistema", loc="left", fontsize=10)
    radio_sistema = novo_radio(ax_menu, SISTEMAS)

    ax_modo = fig.add_axes([0.02, 0.66, 0.25, 0.11])
    ax_modo.set_title("Modelo de arrasto", loc="left", fontsize=10)
    radio_modo = novo_radio(ax_modo, MODOS)

    ax_botao = fig.add_axes([0.02, 0.59, 0.25, 0.05])
    botao_ceq = Button(ax_botao, "Igualar v_terminal ao linear")

    ax_info = fig.add_axes([0.02, 0.05, 0.25, 0.51])
    ax_info.set_xticks([])
    ax_info.set_yticks([])
    for borda in ax_info.spines.values():
        borda.set_color("0.7")
    ax_info.set_facecolor("#f6f6f6")
    texto_info = ax_info.text(0.05, 0.97, "", va="top", ha="left",
                              family="monospace", fontsize=9,
                              transform=ax_info.transAxes)

    # ------------------------------------------------------------------
    # Sistema 3.2 — eixo da trajetória circular
    # ------------------------------------------------------------------
    ax_circ = fig.add_axes([0.36, 0.26, 0.60, 0.64])
    ax_circ.set_box_aspect(1)        # caixa quadrada + limites iguais = escala igual
    ax_circ.set_xlabel("x (m)")
    ax_circ.set_ylabel("y (m)")
    ax_circ.grid(alpha=0.3)
    linha_orbita, = ax_circ.plot([], [], color="0.6", lw=1)
    ax_circ.plot([0], [0], "k+", ms=10)
    linha_fio, = ax_circ.plot([], [], color="0.3", lw=1.5)
    ponto_massa, = ax_circ.plot([], [], "o", color="C0", ms=14)
    seta_props = dict(arrowstyle="-|>", lw=2, shrinkA=0, shrinkB=0)
    seta_fc = ax_circ.annotate("", xy=(0, 0), xytext=(0, 0),
                               arrowprops=dict(color="red", **seta_props))
    seta_v = ax_circ.annotate("", xy=(0, 0), xytext=(0, 0),
                              arrowprops=dict(color="C2", **seta_props))
    ax_circ.text(0.5, 0.02, "vermelho: Fc (aponta para o centro)  ·  verde: v (tangente)\n"
                 "setas fora de escala",
                 transform=ax_circ.transAxes, ha="center", va="bottom",
                 fontsize=8.5, color="0.3")

    # ------------------------------------------------------------------
    # Sistema 3.3 — dois eixos: v(t) e y(t)
    # ------------------------------------------------------------------
    ax_v = fig.add_axes([0.36, 0.55, 0.60, 0.35])
    ax_y = fig.add_axes([0.36, 0.26, 0.60, 0.22], sharex=ax_v)
    ax_v.set_ylabel("velocidade v (m/s)")
    ax_y.set_ylabel("distância caída (m)")
    ax_y.set_xlabel("tempo t (s)")
    plt.setp(ax_v.get_xticklabels(), visible=False)
    for a in (ax_v, ax_y):
        a.grid(alpha=0.3)

    linha_v, = ax_v.plot([], [], color="C0", lw=2)
    linha_v_livre, = ax_v.plot([], [], "--", color="0.4", lw=1.5,
                               label="queda livre  v = g·t")
    linha_vt = ax_v.axhline(0, color="red", ls=":", lw=1.5, label="v_terminal (fórmula)")
    marca_t95 = ax_v.axvline(0, color="green", ls=":", lw=1.5, label="t95 (95% de v_terminal)")
    ponto_t95, = ax_v.plot([], [], "o", color="green", ms=6)
    linha_y, = ax_y.plot([], [], color="C0", lw=2, label="com arrasto")
    linha_y_livre, = ax_y.plot([], [], "--", color="0.4", lw=1.5, label="queda livre")
    ax_y.legend(loc="upper left", fontsize=8.5)

    # ------------------------------------------------------------------
    # Sliders (3 por sistema). Os do arrasto ficam no mesmo lugar dos do MCU;
    # só os do sistema escolhido ficam visíveis.
    # ------------------------------------------------------------------
    def novo_slider(linha, rotulo, vmin, vmax, vini, fmt):
        y = (0.165, 0.105, 0.045)[linha]
        ax = fig.add_axes([0.46, y, 0.44, 0.03])
        return Slider(ax, rotulo, vmin, vmax, valinit=vini, valfmt=fmt)

    sl_R = novo_slider(0, "Raio R (m)", 1.0, 10.0, 5.0, "%.1f")
    sl_v = novo_slider(1, "Velocidade v (m/s)", 1.0, 20.0, 8.0, "%.1f")
    sl_mc = novo_slider(2, "Massa m (kg)", 0.1, 5.0, 1.0, "%.2f")

    sl_m = novo_slider(0, "Massa m (kg)", 0.5, 10.0, 2.0, "%.2f")
    sl_g = novo_slider(1, "Gravidade g (m/s²)", 1.6, 25.0, f.G_TERRA, "%.2f")
    sl_b = novo_slider(2, "Coef. linear b (kg/s)", 0.1, 5.0, 1.0, "%.2f")
    sl_c = novo_slider(2, "Coef. quadrático c (kg/m)", 0.005, 1.0, 0.05, "%.3f")

    # ------------------------------------------------------------------
    # Atualização do sistema 3.2
    # ------------------------------------------------------------------
    def desenhar_marcador():
        """Posiciona massa, fio e setas no instante estado['t'] (fórmulas do MCU)."""
        R, w = estado["R"], estado["w"]
        x, y = f.posicao_circular(estado["t"], R, w)
        linha_fio.set_data([0, x], [0, y])
        ponto_massa.set_data([x], [y])
        L = 0.35 * R                        # comprimento visual das setas
        # Fc aponta da massa para o centro (direção −r)
        seta_fc.xy = (x - L * x / R, y - L * y / R)
        seta_fc.xyann = (x, y)
        # v é tangente ao círculo: vetor unitário (−y, x)/R
        seta_v.xy = (x - L * y / R, y + L * x / R)
        seta_v.xyann = (x, y)

    def atualizar_circular():
        R, v, m = sl_R.val, sl_v.val, sl_mc.val
        w = f.velocidade_angular(v, R)
        T = f.periodo(v, R)

        # Mantém o ângulo θ = ω·t contínuo quando ω muda (o ponto não "pula")
        if estado["w"] is not None:
            estado["t"] *= estado["w"] / w
        estado["w"], estado["R"] = w, R

        t = np.linspace(0, T, 400)                    # uma volta completa
        x, y = f.posicao_circular(t, R, w)
        linha_orbita.set_data(x, y)
        lim = 1.35 * R
        ax_circ.set_xlim(-lim, lim)
        ax_circ.set_ylim(-lim, lim)
        ax_circ.set_title(f"Movimento circular uniforme — R = {R:.1f} m", fontsize=11)
        desenhar_marcador()

        texto_info.set_text(
            "FORÇA CENTRÍPETA (MCU)\n"
            "──────────────────────────\n\n"
            "Velocidade angular\n"
            "  ω = v/R\n"
            f"  ω = {w:.2f} rad/s\n\n"
            "Período\n"
            "  T = 2πR/v\n"
            f"  T = {T:.2f} s\n\n"
            "Aceleração centrípeta\n"
            "  a_c = v²/R\n"
            f"  a_c = {f.aceleracao_centripeta(v, R):.2f} m/s²\n\n"
            "Força centrípeta\n"
            "  Fc = m·v²/R\n"
            f"  Fc = {f.forca_centripeta(m, v, R):.2f} N"
        )

    # ------------------------------------------------------------------
    # Atualização do sistema 3.3
    # ------------------------------------------------------------------
    def atualizar_arrasto():
        m, g = sl_m.val, sl_g.val
        if estado["modo"] == "linear":
            b = sl_b.val
            vt = f.v_terminal_linear(m, b, g)
            t_max = 6 * vt / g                       # 6 × (v_t/g) = 6τ
            t = np.linspace(0, t_max, N_PONTOS)      # vetor de tempo
            v, y = f.v_linear(t, m, b, g), f.y_linear(t, m, b, g)
            t95 = f.tempo_fracao_linear(m, b)
            rotulo = "com arrasto linear"
            cabecalho = ("ARRASTO LINEAR  (F = −b·v)\n"
                         "──────────────────────────\n\n"
                         "Velocidade terminal\n"
                         "  v_t = m·g/b\n"
                         f"  v_t = {vt:.2f} m/s\n\n"
                         "Constante de tempo\n"
                         "  τ = m/b\n"
                         f"  τ = {f.tau_linear(m, b):.2f} s\n\n")
        else:
            c = sl_c.val
            vt = f.v_terminal_quadratico(m, c, g)
            t_max = 6 * vt / g
            t = np.linspace(0, t_max, N_PONTOS)
            v, y = f.v_quadratico(t, m, c, g), f.y_quadratico(t, m, c, g)
            t95 = f.tempo_fracao_quadratico(m, c, g=g)
            rotulo = "com arrasto quadrático"
            cabecalho = ("ARRASTO QUADRÁTICO (F = −c·v²)\n"
                         "──────────────────────────\n\n"
                         "Velocidade terminal\n"
                         "  v_t = √(m·g/c)\n"
                         f"  v_t = {vt:.2f} m/s\n\n"
                         "Tempo característico\n"
                         "  v_t/g\n"
                         f"  v_t/g = {vt / g:.2f} s\n\n")

        v_livre, y_livre = f.v_queda_livre(t, g), f.y_queda_livre(t, g)

        linha_v.set_data(t, v)
        linha_v.set_label(rotulo)
        linha_v_livre.set_data(t, v_livre)
        linha_vt.set_ydata([vt, vt])
        marca_t95.set_xdata([t95, t95])
        ponto_t95.set_data([t95], [0.95 * vt])
        linha_y.set_data(t, y)
        linha_y_livre.set_data(t, y_livre)

        ax_v.set_xlim(0, t_max)
        ax_v.set_ylim(0, 1.5 * vt)       # a queda livre sai do gráfico: é o esperado
        ax_y.set_ylim(0, 1.05 * y_livre[-1])
        ax_v.legend(loc="lower right", fontsize=8.5)
        ax_v.set_title("Queda com arrasto × queda livre", fontsize=11)

        texto_info.set_text(
            cabecalho +
            f"t95 (95% de v_t) = {t95:.2f} s\n\n"
            f"No fim do gráfico (t = {t_max:.1f} s):\n"
            f"  v com arrasto = {v[-1]:.2f} m/s\n"
            f"    ({100 * v[-1] / vt:.1f}% de v_t)\n"
            f"  v queda livre = {v_livre[-1]:.1f} m/s\n"
            f"  caiu {y[-1]:.0f} m (livre: {y_livre[-1]:.0f} m)"
        )

    # ------------------------------------------------------------------
    # Visibilidade e callbacks
    # ------------------------------------------------------------------
    def atualizar_visibilidade():
        circ = estado["sistema"] == "circular"
        quad = estado["modo"] == "quadratico"
        ax_circ.set_visible(circ)
        for s in (sl_R, sl_v, sl_mc):
            s.ax.set_visible(circ)
        for a in (ax_v, ax_y, ax_modo):
            a.set_visible(not circ)
        for s in (sl_m, sl_g):
            s.ax.set_visible(not circ)
        sl_b.ax.set_visible(not circ and not quad)
        sl_c.ax.set_visible(not circ and quad)
        ax_botao.set_visible(not circ and quad)

    def atualizar(_=None):
        if estado["sistema"] == "circular":
            atualizar_circular()
        else:
            atualizar_arrasto()
        fig.canvas.draw_idle()

    def ao_trocar_sistema(rotulo):
        estado["sistema"] = "circular" if rotulo == SISTEMAS[0] else "arrasto"
        atualizar_visibilidade()
        atualizar()
        # A animação só roda no sistema 3.2 (evita redesenhar à toa no 3.3)
        if estado["sistema"] == "circular":
            animacao.event_source.start()
        else:
            animacao.event_source.stop()

    def ao_trocar_modo(rotulo):
        estado["modo"] = "linear" if rotulo == MODOS[0] else "quadratico"
        atualizar_visibilidade()
        atualizar()

    def ao_clicar_igualar(_):
        # c que dá a mesma v_terminal do arrasto linear (mesmos m, g e b)
        c_eq = f.c_equivalente(sl_m.val, sl_b.val, sl_g.val)
        sl_c.set_val(min(max(c_eq, sl_c.valmin), sl_c.valmax))  # dispara atualizar()

    for s in (sl_R, sl_v, sl_mc, sl_m, sl_g, sl_b, sl_c):
        s.on_changed(atualizar)
    radio_sistema.on_clicked(ao_trocar_sistema)
    radio_modo.on_clicked(ao_trocar_modo)
    botao_ceq.on_clicked(ao_clicar_igualar)

    # ------------------------------------------------------------------
    # Animação do MCU: a cada quadro avança o tempo e relê os valores atuais,
    # então nunca fica presa a uma trajetória antiga depois de mexer nos sliders.
    # ------------------------------------------------------------------
    def quadro(_):
        if estado["sistema"] != "circular":
            return
        T = 2 * np.pi / estado["w"]
        estado["t"] = (estado["t"] + DT_ANIM) % T
        desenhar_marcador()

    animacao = FuncAnimation(fig, quadro, interval=INTERVALO_MS,
                             cache_frame_data=False)

    # Estado inicial: prepara os dois sistemas e mostra o 3.2
    atualizar_arrasto()
    atualizar_circular()
    atualizar_visibilidade()

    # Os objetos precisam ser guardados: se forem coletados, os widgets param de responder
    return {"fig": fig, "animacao": animacao, "estado": estado,
            "radio_sistema": radio_sistema, "radio_modo": radio_modo,
            "botao": botao_ceq,
            "sliders": dict(R=sl_R, v=sl_v, m_circ=sl_mc, m=sl_m, g=sl_g, b=sl_b, c=sl_c)}


if __name__ == "__main__":
    ui = criar_interface()
    plt.show()