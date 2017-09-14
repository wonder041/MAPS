#!/aptmp/yanzeli/miniconda3/envs/PIL/bin/python
import Image
# im=Image.open("/user1/scl1/yanzeli/Megaviridae/Figures/Tree_1000_170913.png")
im=Image.open("/user1/scl1/yanzeli/Megaviridae/Figures/Tree_1000_unscaled_170913.png")
main=im.crop((150,500,6250,6250))
label=im.crop((3130,6400,5000,7071))
scale=im.crop((4960,6510,6400,6580))
main.paste(scale,(4470,2620))
main.paste(label,(0,5060))
# main.save("/user1/scl1/yanzeli/Megaviridae/Figures/Tree_1000_170913_trimmed.png")
main.save("/user1/scl1/yanzeli/Megaviridae/Figures/Tree_1000_unscaled_170913_trimmed.png")