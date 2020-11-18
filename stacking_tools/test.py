

# compute.py
def compute(inputs, diagnostic_callback):
    intermediate = [x * 2 for x in inputs]
    diagnostic_callback(intermediate)
    return [x + 1 for x in intermediate]

# coordinator.py
def output_diags(values):
    print(min(values), max(values))

compute([1, 2, 3], output_diags)

# option 1
for xarr in matrix:
    for y in xarr:
        intermed = compute(input[x, y])
        result[x, y] = compute2(intermed)

# option 2
for xarr in matrix:
    for y in xarr:
        intermed[x, y] = compute(input[x, y])

for xarr in matrix:
    for y in xarr:
        result[x, y] = compute2(intermed[x, y])
