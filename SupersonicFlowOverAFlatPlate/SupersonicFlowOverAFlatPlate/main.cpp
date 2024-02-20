#include "SupersonicFlowOverAFlatPlate.h"

void run(std::string wall) {
	SupersonicFlowOverAFlatPlate SFOAFP(1.4, 287, 1.789e-5, 288.16, 0.71, 4, 101325, 288.16, 288.16, wall, 0.00001, 69, 69, 0.6, 1e-8);
	Data data;
	SFOAFP.solve(data);
	data.save(wall);
}

int main() {
	run("isothermal");
	run("adiabatic");
	return 0;
}
