
def f(x,y):
    return y

x0 = 0
y0 = 1
xend = 1
step = 0.01

numSteps = (xend - x0)/step

y = y0
x = x0

i=0

print("steps:")
print(numSteps)

while i < numSteps:
    x += step
    y += step * f(x,y)
    i += 1

print("value:")
print(y)
