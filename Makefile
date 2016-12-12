examples: spiral_motion lorenz_attractor example

%: examples/%.d src/solver.d 
	dmd src/solver.d $< -of=$@

.PHONY: clean
clean:
	rm -f *.o *.txt
