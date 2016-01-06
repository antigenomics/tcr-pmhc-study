
f = open("../app/src/gromacs/fails/1j8h.log", 'r')
lof = f.read()
for line in log:
	if line.find("Fatal error"):
		break;
