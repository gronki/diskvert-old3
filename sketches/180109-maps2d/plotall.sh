
plot() {
  python plot2.py $@ &
  python plot3.py $@ &
  python plot4.py $@ &
  python plot6.py $@ &
  wait
}

plot 1 X X X
plot 1 X 0 X
plot 1 X 1 X
plot 1 X 2 X
plot 1 X 3 X
plot 1 X 4 X
plot 1 X 5 X
plot 1 X 6 X
plot 1 0 X X
plot 1 1 X X
plot 1 2 X X
plot 1 3 X X
plot 1 4 X X
plot 1 5 X X
plot 1 6 X X

plot 2 X X X
plot 2 X 0 X
plot 2 X 1 X
plot 2 X 2 X
plot 2 X 3 X
plot 2 X 4 X
plot 2 X 5 X
plot 2 X 6 X
plot 2 0 X X
plot 2 1 X X
plot 2 2 X X
plot 2 3 X X
plot 2 4 X X
plot 2 5 X X
plot 2 6 X X

plot 0 X X X
plot 0 X 0 X
plot 0 X 1 X
plot 0 X 2 X
plot 0 X 3 X
plot 0 X 4 X
plot 0 X 5 X
plot 0 X 6 X
plot 0 0 X X
plot 0 1 X X
plot 0 2 X X
plot 0 3 X X
plot 0 4 X X
plot 0 5 X X
plot 0 6 X X
