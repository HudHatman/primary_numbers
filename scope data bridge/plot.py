import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# --- KONFIGURACJA WIDOKU ---
FILE_NAME = 'vortex_scope_data.csv'
THEORY_CONST = 0.3375
KEY_NODES = [105, 115, 135] # Dubn, Moskow, Stop

def plot_vortex_data():
    try:
        # Wczytanie danych
        df = pd.read_csv(FILE_NAME)
    except FileNotFoundError:
        print(f"Błąd: Nie znaleziono pliku {FILE_NAME}. Najpierw uruchom kod C++.")
        return

    # Ustawienie stylu oscyloskopu
    plt.style.use('dark_background')
    fig, axs = plt.subplots(4, 1, figsize=(14, 10), sharex=True)
    fig.canvas.manager.set_window_title('Baniowski Scope Preview - Vortex Analysis')
    
    # Kolory kanałów (standard oscyloskopowy)
    colors = ['#FFFF00', '#00FFFF', '#FF00FF', '#00FF00'] # Yellow, Cyan, Magenta, Green
    labels = [
        'CH1: Emisja 3D (sin(n))', 
        'CH2: Echo Inwersyjne 6D (asin+acos)', 
        'CH3: Napięcie Psi (Vortex Detection)', 
        'CH4: Gęstość Sn (B-Density)'
    ]

    for i in range(4):
        col_name = df.columns[i+1] # CH1, CH2, CH3, CH4
        axs[i].plot(df['Time'], df[col_name], color=colors[i], linewidth=1, label=labels[i])
        axs[i].grid(True, color='#333333', linestyle='--')
        axs[i].legend(loc='upper right', fontsize=9)
        
        # Zaznaczanie punktów węzłowych (Dubn, Moskow, Stop)
        for node in KEY_NODES:
            if node in df['Time'].values:
                val = df.loc[df['Time'] == node, col_name].values[0]
                axs[i].annotate(f'Node {node}', xy=(node, val), xytext=(node+5, val+(0.2 if val >=0 else -0.2)),
                                 arrowprops=dict(facecolor='white', shrink=0.05, width=1, headwidth=5),
                                 fontsize=8, color='white')
                axs[i].axvline(x=node, color='white', alpha=0.2, linestyle=':')

    # Dodatkowe ustawienia osi
    axs[0].set_ylabel("Amplituda [V]")
    axs[3].set_xlabel("Wartość n (Liczba Atomowa / Krok Fazowy)")
    
    plt.suptitle(f"Baniowski Digital Scope Preview\nParametr Rezonansu: {THEORY_CONST}", 
                 color='white', fontsize=16, fontweight='bold')
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    print("Wyświetlanie podglądu... Szukaj zapaści sygnału na CH3 w punktach pierwszych.")
    plt.show()

if __name__ == "__main__":
    plot_vortex_data()