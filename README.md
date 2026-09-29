# simulador-segunda-lei-newton

**Atividade de Programação — Leis de Newton e Aplicações**
Mecânica e Física Moderna — Ciência da Computação
Grupo: Guilherme Takahashi, Erick de Paula, Pedro Ulrich

Simulação com interface gráfica (matplotlib.widgets) de dois dos três sistemas propostos
no enunciado (o grupo optou pelos dois abaixo, como o enunciado permite escolher 2 de 3):

- **3.2 Força centrípeta** — movimento circular uniforme, com animação.
- **3.3 Força de arrasto** — queda com arrasto linear ou quadrático, comparada com a queda livre.

Todas as grandezas usam fórmulas fechadas avaliadas sobre um vetor de tempo (NumPy).
Não há integração numérica. A derivação das fórmulas e a discussão dos resultados
estão em `relatorio_ADO2.pdf`.

## Equações principais

**Força centrípeta:** ω = v/R, Fc = m·v²/R, trajetória x(t) = R·cos(ωt), y(t) = R·sin(ωt).

**Força de arrasto** (queda vertical, m·dv/dt = m·g − F_arrasto):
- Linear (F = b·v): v_terminal = m·g/b, v(t) = v_t·(1 − e^(−t/τ)), com τ = m/b.
- Quadrático (F = c·v²): v_terminal = √(m·g/c), v(t) = v_t·tanh(g·t/v_t).

## Abrir a simulação

**No Google Colab (sem instalar nada):**
[![Abrir no Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/GuiReisTaka/simulador-segunda-lei-newton/blob/main/simulacao.ipynb)

Depois de abrir, execute as células em ordem (Ambiente de execução → Executar tudo).

**No computador:**

```bash
pip install -r requirements.txt
python interface.py
```

## Como usar

1. Escolha o sistema no menu *Sistema*.
2. Mexa nos sliders: o gráfico e os valores na caixa de texto atualizam sozinhos.
3. No arrasto, escolha o modelo (linear ou quadrático). O botão *Igualar v_terminal ao linear*
   ajusta `c` para que os dois modelos tenham a mesma velocidade terminal.

## Arquivos

| Arquivo | Conteúdo |
|---|---|
| `fisica.py` | Fórmulas fechadas (sem matplotlib), cada uma comentada com a lei de onde vem |
| `teste_fisica.py` | Validação: as fórmulas satisfazem a 2ª Lei de Newton (`python teste_fisica.py`) |
| `interface.py` | Interface gráfica: menu, sliders, animação |
| `simulacao.ipynb` | Notebook que abre a simulação no Colab |
