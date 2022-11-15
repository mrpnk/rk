#include <functional>
#include <iomanip>

namespace rk
{
	// General logger for N-component states.
	template<int N>
	struct logger : public sensor<vec<N>> {
		double minInterval = 0.0;
		virtual bool measure(double t, vec<N> const& y, vec<N> const& y_t, vec<N> const& erry, int step, double dt, double& forceDt) override {
			if (t >= lastT + minInterval) {
				
				data.push_back(vec<1+2*N>(t,y,y_t));

				lastT = t;
			}
			return false;
		}

		void print(std::ostream& os) {
			os << std::fixed << std::setprecision(16) << std::scientific;
			for (const auto& v : data) {
				os << v << "\n";
			}
			os.flush();
		}

	private:
		double lastT = -999;
		std::vector<vec<1+2*N>> data;
	};

}
