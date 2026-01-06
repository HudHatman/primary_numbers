import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

# --- KONFIGURACJA ---
FILE_NAME = 'vortex_scope_data.csv'
KEY_NODES = [105, 115, 135] # Dubn, Moskow, Stop

def normalize_data(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))

def plot_phase_space_overlay():
    try:
        df = pd.read_csv(FILE_NAME)
    except FileNotFoundError:
        print(f"Błąd: Nie znaleziono pliku {FILE_NAME}.")
        return

    # --- PRZYGOTOWANIE DANYCH ---
    # Oś X: Emisja 3D (CH1)
    x = df['CH1'].values
    # Oś Y: Echo Inwersyjne 6D (CH2)
    y = df['CH2'].values
    
    # Kolor: Napięcie Psi (CH3). Odwracamy, aby 0 (Vortex) było "gorące".
    # Używamy logarytmu, aby uwydatnić małe wartości bliskie zeru.
    c_raw = df['CH3'].values + 1e-9 # Unikamy log(0)
    c_log = np.log10(c_raw)
    # Normalizacja do zakresu 0-1 dla colormapy, odwrócona (1 = Vortex hit)
    colors = 1.0 - normalize_data(c_log)

    # Rozmiar: Gęstość Sn (CH4). Normalizujemy do zakresu wielkości punktów.
    sizes = normalize_data(df['CH4'].values) * 150 + 10 # Punkty od 10 do 160 px

    # --- WIZUALIZACJA ---
    plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(12, 12))
    fig.canvas.manager.set_window_title('Baniowski Phase-Space Vortex Integrator')

    # Rysowanie głównej pętli fazowej jako punktów o zmiennym kolorze i rozmiarze
    # Używamy 'inferno' lub 'magma' - kolory od czarnego/fioletowego (wysokie napięcie)
    # do jasnoseledynowego/białego (ZAPAŚĆ VORTEXU - Niskie napięcie)
    sc = ax.scatter(x, y, s=sizes, c=colors, cmap='magma', alpha=0.6, edgecolors='none')

    # Dodanie paska kolorów (legendy)
    cbar = plt.colorbar(sc, label='Proximity to Vortex (CH3 Inverted)')
    cbar.set_ticks([0, 1])
    cbar.set_ticklabels(['High Voltage (No Vortex)', 'Zero Voltage (VORTEX HIT)'])

    # --- OZNACZANIE WĘZŁÓW KLUCZOWYCH ---
    for node in KEY_NODES:
        if node in df['Time'].values:
            row = df.loc[df['Time'] == node]
            nx, ny = row['CH1'].values[0], row['CH2'].values[0]
            
            # Rysowanie celownika na węźle
            ax.scatter(nx, ny, s=300, facecolors='none', edgecolors='cyan', linewidth=2)
            
            # Etykieta
            label_text = f"N={node}\n(Anchor)" if node == 105 else f"N={node}"
            ax.annotate(label_text, xy=(nx, ny), xytext=(nx+0.1, ny+0.1),
                        arrowprops=dict(facecolor='cyan', shrink=0.05, width=1, headwidth=6),
                        color='cyan', fontsize=11, fontweight='bold')

    # --- ESTETYKA WYKRESU ---
    ax.set_xlabel("Oś X: Emisja 3D (CH1: sin(n))", fontsize=12, color='yellow')
    ax.set_ylabel("Oś Y: Inwersja 6D (CH2: asin+acos)", fontsize=12, color='cyan')
    ax.set_title("Zintegrowany Portret Fazowy Węzła Baniowskiego\n"
                 "(Układ Zamknięty: Interakcja Emisja-Zasysanie)", color='white', fontsize=14)

    ax.grid(True, color='#333333', linestyle=':')
    ax.set_aspect('equal', 'box') # Ważne: zachowanie proporcji kwadratu
    
    # Opcjonalnie: usunięcie osi liczbowych dla czystszego obrazu topologii
    # ax.set_xticks([])
    # ax.set_yticks([])

    plt.tight_layout()
    print("Generowanie zintegrowanego obrazu fazowego...")
    plt.show()

if __name__ == "__main__":
    plot_phase_space_overlay()