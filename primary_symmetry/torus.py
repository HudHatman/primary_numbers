import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# --- KONFIGURACJA MODELU BANIOWSKIEGO ---
# Stała rezonansu definiująca głębokość "zasysania"
B_CONST = 0.3375
# Punkty kluczowe do zaznaczenia na węźle (w stopniach)
KEY_POINTS_DEG = [105, 115, 135]

# --- GENEROWANIE DANYCH PRZESTRZENNYCH ---
# Rozdzielczość symulacji (ilość punktów na jeden obrót)
resolution = 2000
# Kąt theta od 0 do 2*pi (jeden pełny obrót)
theta = np.linspace(0, 2 * np.pi, resolution)

# --- OBLICZANIE 6 FUNKCJI TRYGONOMETRYCZNYCH ---
# Warstwa 1: Emisja 3D (Fale nośne)
sin_t = np.sin(theta)
cos_t = np.cos(theta)
# Tan wymaga przycięcia (clipping) asymptot, aby nie "rozsadzić" wykresu
tan_t = np.clip(np.tan(theta), -4, 4) 

# Warstwa 2: Inwersja 6D (Mechanizm zasysający)
# Używamy złożeń funkcji, aby pokazać "zwijanie" przestrzeni
asin_sin = np.arcsin(sin_t)
acos_cos = np.arccos(cos_t)
atan_tan = np.arctan(tan_t)

# --- RÓWNANIA PARAMETRYCZNE WĘZŁA BANIOWSKIEGO ---
# Definiujemy geometrię torusa, który jest deformowany przez funkcje.
# R - promień główny torusa, r - promień rury
R_major = 4.0
r_minor_base = 1.5

# MODULACJA PROMIENIA (Zasysanie 6D):
# Promień rury "oddycha" w zależności od sumy funkcji inwersyjnych skalowanych przez B_CONST.
# Kiedy asin+acos dąży do PI/2 (rezonans), rura się zwęża.
suction_factor = (np.abs(asin_sin) + np.abs(acos_cos)) / (np.pi/2)
r_modulated = r_minor_base * (1 - B_CONST * 0.5 * (suction_factor - 1))

# Kąt skręcenia torusa (aby stworzyć węzeł)
phi = 3 * theta 

# Współrzędne 3D (X, Y - płaszczyzna fazy, Z - amplituda energii)
# X i Y są modulowane przez "oddech" inwersyjny
X = (R_major + r_modulated * np.cos(phi)) * cos_t
Y = (R_major + r_modulated * np.cos(phi)) * sin_t

# Z jest zdominowane przez emisję (tan), ale tłumione przez inwersję (arctan) i stałą B.
# To tworzy "szpilki" energetyczne.
Z = r_modulated * np.sin(phi) + tan_t * (1 - B_CONST * np.abs(atan_tan)/np.pi)

# --- WIZUALIZACJA (MATPLOTLIB) ---
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')
fig.patch.set_facecolor('black') # Czarne tło dla kontrastu
ax.set_facecolor('black')

# Kolorowanie 4 kwadrantów (po 25% obrotu)
colors = []
for t in theta:
    deg = np.degrees(t)
    if 0 <= deg < 90: colors.append('#FF00FF')   # Q1: Magenta (Inicjacja)
    elif 90 <= deg < 180: colors.append('#00FFFF') # Q2: Cyan (Strefa Dubnu -105)
    elif 180 <= deg < 270: colors.append('#FFFF00')# Q3: Żółty (Strefa Przeskoku 180+)
    else: colors.append('#FF4500')                 # Q4: Czerwony (Rekonstrukcja)

# Rysowanie głównej struktury węzła (jako punkty zlewające się w linię dla gradacji koloru)
# Używamy scatter dla precyzyjnego kolorowania każdego punktu
ax.scatter(X, Y, Z, c=colors, s=2, alpha=0.6, linewidth=0)

# --- ZAZNACZANIE PUNKTÓW KLUCZOWYCH (105, 115, 135) ---
for deg_point in KEY_POINTS_DEG:
    # Znajdź indeks odpowiadający danemu kątowi
    idx = (np.abs(np.degrees(theta) - deg_point)).argmin()
    
    # Współrzędne punktu
    px, py, pz = X[idx], Y[idx], Z[idx]
    
    # Rysowanie "szpilki" i etykiety
    ax.scatter(px, py, pz, color='white', s=100, marker='o', edgecolors='red', linewidth=2)
    
    # Dodanie etykiety z przesunięciem
    label = f"N={deg_point}\n(Sn Node)"
    ax.text(px, py, pz + 0.5, label, color='white', fontsize=10, ha='center')
    
    # Rysowanie linii pionowej wskazującej "kotwiczenie" w punkcie
    ax.plot([px, px], [py, py], [pz - 1, pz + 1], color='white', linestyle='--', alpha=0.5)

# --- USTAWIENIA KOŃCOWE WYKRESU ---
ax.set_title("Matematyczny Węzeł Energii Baniowskiego (Model 6D)\nInterakcja Emisji (tan) i Zasysania (asin/acos) ze stałą B=0.3375", 
             color='white', fontsize=14)
ax.set_xlabel("Oś X (Faza cos)", color='gray')
ax.set_ylabel("Oś Y (Faza sin)", color='gray')
ax.set_zlabel("Oś Z (Amplituda Emisji/Zasysania)", color='gray')

# Ukrycie osi i siatki dla efektu "pustki"
ax.grid(False)
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])
ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False
ax.xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))

# Ustawienie początkowego widoku (aby dobrze widzieć punkty 105-135)
ax.view_init(elev=30, azim=110)

plt.show()
