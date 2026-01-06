import matplotlib.pyplot as plt
import numpy as np
import re

# Stałe modelu BCO
B_CONST = 0.3375
R_BASE = 5.0

# Parsowanie danych z pliku 5_data.txt
def load_bco_data(path):
    with open(path, 'r') as f:
        content = f.read()
    entries = content.split('--------------------')
    parsed = []
    for e in entries:
        a = re.search(r'___angle___\n(-?\d+)', e)
        s = re.search(r'___table_sum___\n(\d+)', e)
        p = re.search(r'___is_sum_prime___\n(true|false)', e)
        if a and s and p:
            parsed.append({'angle': int(a.group(1)), 'sum': int(s.group(1)), 'is_prime': p.group(1) == 'true'})
    return parsed

data = load_bco_data('5_data.txt')

# Generowanie współrzędnych 3D
x, y, z, colors, sizes = [], [], [], [], []

for entry in data:
    angle = entry['angle']
    theta = np.radians(angle)
    r = R_BASE * (1.0 - B_CONST) if entry['is_prime'] else R_BASE
    
    x.append(r * np.cos(theta))
    y.append(r * np.sin(theta))
    z.append(angle / 40.0) # Oś pionowa helisy
    
    if angle == -105: # Kotwica
        colors.append('blue'); sizes.append(150)
    elif angle == -90: # Oś
        colors.append('cyan'); sizes.append(100)
    elif angle == 115: # Akcelerator
        colors.append('red'); sizes.append(180)
    elif entry['is_prime']:
        colors.append('orange'); sizes.append(40)
    else:
        colors.append('gray'); sizes.append(5)

# Rysowanie w 3D
fig = plt.figure(figsize=(12, 10), facecolor='black')
ax = fig.add_subplot(111, projection='3d', facecolor='black')

ax.plot(x, y, z, color='white', alpha=0.3) # Trajektoria
ax.scatter(x, y, z, c=colors, s=sizes, alpha=0.8) # Węzły

# Etykiety punktów kluczowych
for i, a in enumerate([d['angle'] for d in data]):
    if a in [-105, -90, 115]:
        label = "Anchor (-105)" if a == -105 else ("Axis (-90)" if a == -90 else "Accelerator (115)")
        ax.text(x[i], y[i], z[i] + 0.5, label, color='white', fontsize=11, fontweight='bold')

ax.set_title("Struktura 3D Cewki BCO: Oś, Kotwica i Akcelerator", color='white')
ax.set_axis_off()
plt.show()