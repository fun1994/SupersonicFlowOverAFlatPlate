#pragma once
#include <vector>
#include <fstream>

class Data {
public:
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> t;
	std::vector<std::vector<std::vector<double>>> rho;
	std::vector<std::vector<std::vector<double>>> u;
	std::vector<std::vector<std::vector<double>>> v;
	std::vector<std::vector<std::vector<double>>> T;
	std::vector<std::vector<std::vector<double>>> p;
	std::vector<std::vector<std::vector<double>>> Ma;
	void save(std::vector<double>& data, std::string path, std::string filename) {
		std::ofstream file("./data/" + path + "/" + filename + ".txt");
		for (int i = 0; i < data.size(); i++) {
			file << data[i];
			if (i < data.size() - 1) {
				file << " ";
			}
		}
		file.close();
	}
	void save(std::vector<std::vector<double>>& data, std::string path, std::string filename) {
		std::ofstream file("./data/" + path + "/" + filename + ".txt");
		for (int i = 0; i < data.size(); i++) {
			for (int j = 0; j < data[i].size(); j++) {
				file << data[i][j];
				if (j < data[i].size() - 1) {
					file << " ";
				}
			}
			if (i < data.size() - 1) {
				file << "\n";
			}
		}
		file.close();
	}
	void save(std::string path) {
		save(x, path, "x");
		save(y, path, "y");
		save(t, path, "time");
		save(rho[rho.size() - 1], path, "rho");
		save(u[u.size() - 1], path, "u");
		save(v[v.size() - 1], path, "v");
		save(T[T.size() - 1], path, "T");
		save(p[p.size() - 1], path, "p");
		save(Ma[Ma.size() - 1], path, "Ma");
	}
};
