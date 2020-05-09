def f(x,y):
    return (y-2)*y

x0 = 0
y0 = 1
xend = 5
step = 0.00001

numSteps = (xend - x0)/step

y = y0
x = x0

i=0

print("steps:")
print(numSteps)

while i < numSteps:
    xnew = x + step

    ynew = y + step * f(xnew,y)

    y = y + step/2 * (f(x,y) + f(xnew, ynew))

    i += 1

print("value:")
print(y)
