set -e
rm -f d2.png d2.cool.png
cat input.par | dv-mag -n 30000 -top 120 -precision 1e-6 -balance -o d2
diskvert-pack d2
diskvert-new-plot d2.tgz -o d2.png
eog -n d2.png &
diskvert-cooling2D d2.tgz -o d2.cool.png
eog -n d2.cool.png &
