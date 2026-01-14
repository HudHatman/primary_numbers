import numpy as np
import matplotlib.pyplot as plt

def is_stable_weight(n):
    """
    Test stabilności 'atomu matematyki'.
    W tym modelu liczba pierwsza = stabilna waga atomowa.
    """
    if n < 2: return False
    for i in range(2, int(np.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def generate_resonance_field(limit=50):
    """
    Generuje pole rezonansowe na podstawie bombardowania trygonometrycznego.
    """
    # Gęsta oś X dla gładkich wykresów
    x = np.linspace(0, limit, 3000)

    # 1. Definicja Linii Bazowej: 2n = 2n + 1 - 1
    # Fizycznie: Poziom zerowy z fluktuacjami kwantowymi (szum tła)
    baseline = np.zeros_like(x)
    vacuum_fluctuation = 0.02 * np.sin(50 * x) # mikro-fluktuacje (+1 -1)

    # 2. Inicjalizacja sygnału prawdopodobieństwa istnienia
    existence_probability = baseline + vacuum_fluctuation

    # 3. Bombardowanie i Normalizacja
    # Iterujemy po liczbach naturalnych szukając punktów rezonansu
    for n in range(1, limit + 1):
        if is_stable_weight(n):
            # Znaleziono stabilną wagę!
            # Amplituda rezonansu (prawdopodobieństwo istnienia) = 1.0 (znormalizowane)

            # Parametry piku Gaussa
            center = n
            width = 0.25

            # Składowa trygonometryczna: sin(t) + cos(t) / 2
            # Moduluje ona fazę powstawania cząstki
            t = x
            trig_modulator = np.sin(t) + np.cos(t) / 2

            # Normalizacja Amplitudy:
            # Dążymy do uzyskania pików równej wysokości niezależnie od n
            # Wykorzystujemy funkcję Gaussa jako obwiednię
            peak = np.exp(-0.5 * ((x - center) / width)**2)

            # Dodajemy znormalizowany pik do pola
            existence_probability += peak

    # 4. Wizualizacja
    plt.figure(figsize=(15, 8), facecolor='#0a0a0a')
    ax = plt.gca()
    ax.set_facecolor('#0a0a0a')

    # Rysowanie głównego sygnału
    plt.plot(x, existence_probability, color='#00ffcc', lw=1.5,
             label='Amplituda Rezonansu (Prawdopodobieństwo Istnienia)', alpha=0.9)

    # Wypełnienie pod pikami dla efektu 'chmury kwantowej'
    plt.fill_between(x, 0, existence_probability, color='#00ffcc', alpha=0.1)

    # Oznaczenie linii bazowej 2n = 2n + 1 - 1
    plt.axhline(0, color='white', linestyle='--', alpha=0.3, label='Baseline: 2n = 2n + 1 - 1')

    # Podpisy stabilnych wag (dynamicznie znalezionych)
    stable_points = [n for n in range(limit + 1) if is_stable_weight(n)]
    for p in stable_points:
        plt.annotate(f'{p}', xy=(p, 1.05), color='#ff3366',
                     ha='center', fontsize=10, fontweight='bold')

    # Personalizacja wykresu
    plt.title('KWANTOWE POLE REZONANSOWE: SYGNAŁ STABILNYCH WAG ATOMOWYCH', color='white', fontsize=16)
    plt.xlabel('Skalarna Wartość Energii / Masa (n)', color='white')
    plt.ylabel('Prawdopodobieństwo Istnienia P(E)', color='white')
    plt.grid(color='white', alpha=0.1)
    plt.legend(facecolor='#1a1a1a', edgecolor='white', labelcolor='white')
    plt.xlim(0, limit)
    plt.ylim(-0.5, 1.5)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    generate_resonance_field(60)
