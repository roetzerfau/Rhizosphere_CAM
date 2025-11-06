import matplotlib.pyplot as plt
from matplotlib.offsetbox import (AnnotationBbox, DrawingArea, OffsetImage,
                                  TextArea)
from PIL import Image
import os

pic_path = "./C_N_SVector/"
pic_names = ["cn40_day150.png","cn40_day155.png",
	"cn40_withoutPartMov_day150.png","cn40_withoutPartMov_day155.png"]

positionen = [(0, 0), (1, 0),
	 (0, 1), (1, 1)]  # x, y-Koordinaten für jedes Bild

fig, ax = plt.subplots()

ax.set_xlim(-0.5, 1.5)
ax.set_ylim(-0.5, 1.5)

ax.set_xticks([0,1], ['150', '155'])

ax.set_yticks([0,1], ['1a', '1b'])

# Bilder hinzufügen
len=0.9
for name, (x, y) in zip(pic_names, positionen):
    pic = Image.open(pic_path+name)
    extent = (x-len/2, x +len/2, y-len/2, y+len/2)
    ax.imshow(pic, extent=extent)

ax.set_aspect('equal')
#ax.grid(True)
plt.xlabel("Day")
plt.ylabel("Scenario")
#plt.show()
plt.savefig('pic.png', dpi=500)
