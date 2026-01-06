import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

# Parametry modelu Baniowskiego
B_CONST = 0.3375
RAD_BASE = 5

# Definicja trajektorii zdarzenia: Oś -> Kotwica -> Akcelerator
path = np.concatenate([
    np.linspace(-90, -105, 40),   # Faza I -> II
    np.linspace(-105, 115, 80)    # Faza II -> III
])

fig = plt.figure(figsize=(10, 8), facecolor='black')
ax = fig.add_subplot(111, projection='3d', facecolor='black')

# Wizualizacja tła fazowego (Przestrzeń 6D)
t = np.linspace(-200, 200, 500)
x_bg = (RAD_BASE + np.cos(3 * np.radians(t))) * np.cos(np.radians(t))
y_bg = (RAD_BASE + np.cos(3 * np.radians(t))) * np.sin(np.radians(t))
z_bg = np.sin(np.radians(t * 2)) * 3
ax.plot(x_bg, y_bg, z_bg, color='gray', alpha=0.2)

# Elementy animowane
point, = ax.plot([], [], [], 'ro', markersize=10)
vortex, = ax.plot([], [], [], 'cyan', alpha=0.5)
txt = ax.text2D(0.05, 0.95, "", transform=ax.transAxes, color='white', fontsize=12)

def update(i):
    angle = path[i]
    theta = np.radians(angle)
    
    # Modulacja promienia stałą B_CONST w punktach Prime
    # -105 (sum 97) i -90 (sum 37) wywołują zapaść (suction)
    is_vortex = abs(angle - (-105)) < 3 or abs(angle - (-90)) < 3
    r = RAD_BASE * (1.0 - (B_CONST if is_vortex else 0))
    
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    z = np.sin(theta * 2) * 3
    
    point.set_data([x], [y])
    point.set_3d_properties([z])
    
    # Pierścień wiru wokół punktu
    p = np.linspace(0, 2*np.pi, 30)
    vx = x + 0.5 * np.cos(p)
    vy = y + 0.5 * np.sin(p)
    vz = np.full_like(p, z)
    vortex.set_data(vx, vy)
    vortex.set_3d_properties(vz)

    # Status operacyjny
    if abs(angle - (-90)) < 2:
        txt.set_text("STATUS: AXIS DETECTED (-90) | SUM: 37 (PRIME)")
        point.set_color('yellow')
    elif abs(angle - (-105)) < 2:
        txt.set_text("STATUS: ANCHOR STABILIZED (-105) | SUM: 97 (PRIME)")
        point.set_color('blue')
    elif abs(angle - 115) < 2:
        txt.set_text("STATUS: ACCELERATOR EMISSION (115) | SUM: 72")
        point.set_color('red')
    else:
        txt.set_text(f"Current Phase Angle: {angle:.1f}")

    return point, vortex, txt

ani = animation.FuncAnimation(fig, update, frames=len(path), interval=50, blit=True)
ax.set_axis_off()
plt.show()