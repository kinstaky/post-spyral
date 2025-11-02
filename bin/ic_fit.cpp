#include <TF1.h>
#include <TGraph.h>
#include <TFitResult.h>

int main() {
	constexpr double ic_value[] = {275.0, 776.0, 1147.0, 1590.0, 2063.0};
	for (int i = 0; i < 100; ++i) {
		TF1 *f1 = new TF1("f1", TString::Format("[0]*(x+%d)*(x+%d)+[1]", i, i), 0.0, 2500.0);
		TGraph *g = new TGraph;
		// add points
		for (int i = 0; i < 5; i++) {
			g->AddPoint(i, ic_value[i]);
		}
		auto fit_result = g->Fit(f1, "QRS");
		std::cout << i << ", " << fit_result->Chi2() << "\n";
	}
	return 0;
}