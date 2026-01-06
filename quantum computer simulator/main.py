import matplotlib.pyplot as plt
import numpy as np
import re
from matplotlib.animation import FuncAnimation

# Konfiguracja modelu
B_CONST = 0.3375
R_BASE = 5.0

def parse_quantum_data(file_path):
    with open(file_path, 'r') as f:
        content = f.read()
    
    entries = content.split('--------------------')
    parsed_data = []
    
    for entry in entries:
        angle_match = re.search(r'___angle___\n(-?\d+)', entry)
        sum_prime_match = re.search(r'___is_sum_prime___\n(true|false)', entry)
        table_match = re.search(r'___table___\n([\s\S]*?)\n___', entry)
        
        if angle_match and sum_prime_match and table_match:
            angle = int(angle_match.group(1))
            is_prime = sum_prime_match.group(1) == 'true'
            # Ekstrakcja wszystkich stanów kwantowych (liczb) z tabeli
            table_vals = [int(x) for x in re.findall(r'\d+', table_match.group(1))]
            
            parsed_data.append({
                'angle': angle,
                'is_prime': is_prime,
                'states': table_vals
            })
    return parsed_data

# Inicjalizacja danych
data = parse_quantum_data('5_data.txt')
data.sort(key=lambda x: x['angle'])

fig = plt.figure(figsize=(10, 10), facecolor='black')
ax = fig.add_subplot(111, projection='3d', facecolor='black')
ax.set_axis_off()

# Obiekty animacji
cloud_scatter = ax.scatter([], [], [], alpha=0.6)
node_point, = ax.plot([], [], [], 'wo', markersize=8)
status_text = ax.text2D(0.05, 0.95, "", transform=ax.transAxes, color='white', fontsize=10)

def update(i):
    d = data[i]
    angle = d['angle']
    theta = np.radians(angle)
    r_core = R_BASE * (1.0 - B_CONST) if d['is_prime'] else R_BASE
    
    # Pozycja jądra
    cx, cy, cz = r_core * np.cos(theta), r_core * np.sin(theta), angle / 40.0
    node_point.set_data([cx], [cy])
    node_point.set_3d_properties([cz])
    node_point.set_color('cyan' if d['is_prime'] else 'gray')
    
    # Symulacja chmury stanów kwantowych wokół jądra
    states = np.array(d['states'])
    if len(states) > 0:
        phis = np.random.uniform(0, 2*np.pi, len(states))
        thetas = np.arccos(np.random.uniform(-1, 1, len(states)))
        # Skalowanie promienia stanu (logarytmiczne dla czytelności skoków energii)
        radii = np.log1p(states) * 0.4 
        
        sx = cx + radii * np.sin(thetas) * np.cos(phis)
        sy = cy + radii * np.sin(thetas) * np.sin(phis)
        sz = cz + radii * np.cos(thetas)
        
        cloud_scatter._offsets3d = (sx, sy, sz)
        cloud_scatter.set_color(plt.cm.magma(np.linspace(0.3, 1, len(states))))
        cloud_scatter.set_sizes(radii * 15)

    status = "VORTEX NODE" if d['is_prime'] else "STABLE POINT"
    status_text.set_text(f"Angle: {angle} | States: {len(states)} | {status}")
    ax.view_init(elev=25, azim=i*0.5)
    return node_point, cloud_scatter, status_text

# Generowanie GIF (co 4-ty krok dla optymalizacji rozmiaru)
ani = FuncAnimation(fig, update, frames=range(0, len(data), 4), interval=80)
ani.save('quantum_simulation.gif', writer='pillow', fps=12)